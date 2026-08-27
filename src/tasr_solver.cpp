#include "tasr_solver.h"

#include "solver_backend.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
    constexpr double kTolerance = 1e-7;
    constexpr double kFixingTolerance = 1e-9;

    // Use one consistent tolerance for all constraints/bounds that preserve
    // already-proved lexicographic levels. This is intentionally larger than
    // kFixingTolerance because these values originate from solver-produced
    // floating-point optima and must remain consistent with MIP feasibility
    // tolerances across repeated solves.
    constexpr double kLexTolerance = 1e-6;
    constexpr double kLexFixingTolerance = kLexTolerance;

    // Diagnostics are deliberately separated from the mathematical tolerances.
    // A warning means "worth inspecting"; an error means that the persistent
    // model and the physical loads implied by x disagree too much to trust
    // further bottleneck fixings.
    constexpr double kInvariantWarningTolerance = 1e-7;
    constexpr double kInvariantErrorTolerance = 1e-5;
    constexpr double kBottleneckProofTolerance = 1e-8;
    constexpr double kInf = std::numeric_limits<double>::infinity();

    using opt::BackendOptions;
    using opt::LinearModel;
    using opt::MipStart;
    using opt::ModelBuilder;
    using opt::SolverResult;
    using opt::VariableType;

    struct Layout {
        int demands = 0;
        int times = 0;
        int nodes = 0;
        int arcs = 0;

        std::vector<int> x;       // [d,t,i,j], -1 if unavailable / i==j
        std::vector<int> lambda;  // [t,a], -1 if arc down
        std::vector<int> y;       // [d,t,i,j], only t>=1 and union exists
        std::vector<int> load_row; // [t,a], row for capacity*lambda - traffic
        std::vector<int> max_row;  // [t,a], persistent z >= lambda row
        // In models k >= 2, alpha/beta are indexed directly by the
        // lexicographic iteration number. Indices 0 and 1 remain unused:
        // iteration 1 is represented by the max-load constraints.
        std::vector<int> alpha;
        std::vector<std::vector<int>> beta; // [lexicographic iteration][load position]
        int z = -1;

        std::vector<std::pair<int, int>> load_position; // (t,a)

        std::size_t segmentOffset(int d, int t, int i, int j) const {
            return (((static_cast<std::size_t>(d) * times + t) * nodes + i) * nodes + j);
        }
        std::size_t loadOffset(int t, int a) const {
            return static_cast<std::size_t>(t) * arcs + a;
        }
    };


    // Reduced model used only to construct the first MIP start.
    // It uses one routing variable x[d,i,j] shared by all time periods,
    // while loads are still evaluated separately for every actual (t,a).
    struct CommonRoutingWarmStart {
        std::vector<double> x; // [d,i,j]
        std::vector<double> loads; // all available [t,a] loads
        double max_load = kInf;
    };

    // A load that has been proved unavoidable at one lexicographic level.
    // Once proved, lambda[t,a] can be fixed to this value in all later models.
    struct FixedLoad {
        int time = -1;
        int arc = -1;
        double value = 0.0;
    };

    struct CommonRoutingLayout {
        int demands = 0;
        int times = 0;
        int nodes = 0;
        int arcs = 0;

        std::vector<int> x;       // [d,i,j], -1 if not usable at every t
        std::vector<int> lambda;  // [t,a], -1 if arc down
        int z = -1;

        std::vector<std::pair<int, int>> load_position; // (t,a)

        std::size_t segmentOffset(int d, int i, int j) const {
            return ((static_cast<std::size_t>(d) * nodes + i) * nodes + j);
        }

        std::size_t loadOffset(int t, int a) const {
            return static_cast<std::size_t>(t) * arcs + a;
        }
    };

    bool arcIsAvailable(const InstanceData& instance, int time, int arc) {
        const int external = instance.arc_external_ids[arc];
        return external >= 0 &&
            external < static_cast<int>(instance.link_intervened_at_time[time].size()) &&
            instance.link_intervened_at_time[time][external] == 0;
    }

    int arcNameId(const InstanceData& instance, int arc) {
        const int external = instance.arc_external_ids[arc];
        return external >= 0 ? external : arc;
    }

    bool segmentIsReachable(const SplitCoefficients& splits, int time, int i, int j) {
        return i != j && !splits.nonzero(time, i, j).empty();
    }

    bool segmentIsReachableAtAllTimes(const SplitCoefficients& splits,
        int times,
        int i,
        int j) {
        if (i == j) return false;
        for (int t = 0; t < times; ++t) {
            if (!segmentIsReachable(splits, t, i, j)) return false;
        }
        return true;
    }


    std::optional<double> fixedLoadValue(
        const std::vector<FixedLoad>& fixed_loads,
        int time,
        int arc) {

        for (const FixedLoad& fixed : fixed_loads) {
            if (fixed.time == time && fixed.arc == arc) {
                return fixed.value;
            }
        }
        return std::nullopt;
    }

    bool loadIsFixed(
        const std::vector<FixedLoad>& fixed_loads,
        int time,
        int arc) {

        return fixedLoadValue(fixed_loads, time, arc).has_value();
    }


    void appendTerm(std::vector<int>& indices,
        std::vector<double>& values,
        int column,
        double coefficient) {
        if (column < 0 || std::abs(coefficient) <= 1e-15) return;
        indices.push_back(column);
        values.push_back(coefficient);
    }


    double segmentLoadLowerBound(
        const InstanceData& instance,
        const SplitCoefficients& splits,
        int d,
        int t,
        int i,
        int j) {

        const double volume = instance.demand_matrix[d][t];
        if (volume <= 0.0) return 0.0;

        double lower_bound = 0.0;

        // A selected segment sends r_{i,j,a,t} * v_{d,t} traffic over arc a.
        // Hence this segment alone forces at least this utilization on a.
        for (const auto& entry : splits.nonzero(t, i, j)) {
            const int a = entry.arc_index;
            if (!arcIsAvailable(instance, t, a)) continue;

            const double capacity = instance.arc_capacities[a];
            if (capacity <= 0.0) {
                throw std::runtime_error(
                    "non-positive capacity on available arc " +
                    std::to_string(arcNameId(instance, a)));
            }

            lower_bound = std::max(
                lower_bound,
                entry.coefficient * volume / capacity);
        }

        return lower_bound;
    }

    void addBaseVariables(ModelBuilder& model,
        const InstanceData& instance,
        const SplitCoefficients& splits,
        Layout& layout,
        double segment_pruning_cutoff = kInf,
        double free_load_upper_bound = kInf,
        const std::vector<FixedLoad>& fixed_loads = {}) {

        layout.demands = instance.numDemands();
        layout.times = instance.num_time_slots;
        layout.nodes = instance.numNodes();
        layout.arcs = instance.numArcs();

        const std::size_t segment_count =
            static_cast<std::size_t>(layout.demands) * layout.times *
            layout.nodes * layout.nodes;
        layout.x.assign(segment_count, -1);
        layout.y.assign(segment_count, -1);
        layout.lambda.assign(static_cast<std::size_t>(layout.times) * layout.arcs, -1);
        layout.load_row.assign(static_cast<std::size_t>(layout.times) * layout.arcs, -1);
        layout.max_row.assign(static_cast<std::size_t>(layout.times) * layout.arcs, -1);

        for (int d = 0; d < layout.demands; ++d) {
            for (int t = 0; t < layout.times; ++t) {
                for (int i = 0; i < layout.nodes; ++i) {
                    for (int j = 0; j < layout.nodes; ++j) {
                        if (!segmentIsReachable(splits, t, i, j)) continue;

                        // This pruning must use a valid upper bound for the GLOBAL
                        // maximum load (normally z1), not the progressively smaller
                        // bound for the still-unfixed loads. A segment may legitimately
                        // contribute more than the current free-load boundary to a
                        // previously fixed bottleneck.
                        if (std::isfinite(segment_pruning_cutoff) &&
                            segmentLoadLowerBound(
                                instance, splits, d, t, i, j) >
                            segment_pruning_cutoff + kLexFixingTolerance) {
                            continue;
                        }

                        layout.x[layout.segmentOffset(d, t, i, j)] =
                            model.addColumn(
                                "x_d" + std::to_string(d) +
                                "_t" + std::to_string(t) +
                                "_i" + std::to_string(i) +
                                "_j" + std::to_string(j),
                                0.0, 0.0, 1.0, VariableType::Binary);
                    }
                }
            }
        }

        for (int t = 0; t < layout.times; ++t) {
            for (int a = 0; a < layout.arcs; ++a) {
                if (!arcIsAvailable(instance, t, a)) continue;

                double lambda_lower = 0.0;
                double lambda_upper =
                    std::isfinite(free_load_upper_bound)
                    ? free_load_upper_bound + kLexTolerance
                    : kInf;

                // A proved bottleneck is fixed within the common numerical
                // lexicographic window. It is deliberately exempt from the smaller
                // upper bound imposed on the remaining free loads.
                if (const auto fixed = fixedLoadValue(fixed_loads, t, a)) {
                    lambda_lower =
                        std::max(0.0, *fixed - kLexFixingTolerance);
                    lambda_upper =
                        *fixed + kLexFixingTolerance;
                }

                layout.lambda[layout.loadOffset(t, a)] =
                    model.addColumn(
                        "lambda_a" + std::to_string(arcNameId(instance, a)) +
                        "_t" + std::to_string(t),
                        0.0,
                        lambda_lower,
                        lambda_upper);
                layout.load_position.emplace_back(t, a);
            }
        }

        for (int d = 0; d < layout.demands; ++d) {
            for (int t = 1; t < layout.times; ++t) {
                for (int i = 0; i < layout.nodes; ++i) {
                    for (int j = 0; j < layout.nodes; ++j) {
                        if (i == j) continue;
                        const auto current = layout.x[layout.segmentOffset(d, t, i, j)];
                        const auto previous = layout.x[layout.segmentOffset(d, t - 1, i, j)];
                        if (current < 0 && previous < 0) continue;

                        // y represents |x_t - x_{t-1}| for binary x, hence 0 <= y <= 1.
                        layout.y[layout.segmentOffset(d, t, i, j)] =
                            model.addColumn(
                                "y_d" + std::to_string(d) +
                                "_t" + std::to_string(t) +
                                "_i" + std::to_string(i) +
                                "_j" + std::to_string(j),
                                0.0, 0.0, 1.0);
                    }
                }
            }
        }
    }


    void addBaseConstraints(ModelBuilder& model,
        const InstanceData& instance,
        const SplitCoefficients& splits,
        Layout& layout,
        const std::vector<FixedLoad> & /*fixed_loads*/ = {}) {
        // Arc-load definition:
        //
        //     capacity[a] * lambda[a,t] - routed_traffic[a,t] = 0.
        //
        // lambda is the actual physical utilization, not merely an upper bound.
        // Keeping this as an equality for every arc prevents unused/free lambda
        // variables from carrying artificial slack and corrupting later peeled
        // lexicographic max objectives.
        for (const auto [t, a] : layout.load_position) {
            std::vector<int> index;
            std::vector<double> value;
            appendTerm(index, value, layout.lambda[layout.loadOffset(t, a)],
                instance.arc_capacities[a]);

            for (int d = 0; d < layout.demands; ++d) {
                const double volume = instance.demand_matrix[d][t];
                if (volume <= 0.0) continue;
                for (int i = 0; i < layout.nodes; ++i) {
                    for (int j = 0; j < layout.nodes; ++j) {
                        const int x = layout.x[layout.segmentOffset(d, t, i, j)];
                        if (x < 0) continue;
                        const double r = splits.value(i, j, a, t);
                        appendTerm(index, value, x, -r * volume);
                    }
                }
            }
            layout.load_row[layout.loadOffset(t, a)] =
                model.addRow(
                    "load_a" + std::to_string(arcNameId(instance, a)) +
                    "_t" + std::to_string(t),
                    0.0,
                    0.0,
                    std::move(index),
                    std::move(value));
        }

        // Flow conservation.
        for (int d = 0; d < layout.demands; ++d) {
            const int source = instance.demands[d].source_id;
            const int target = instance.demands[d].target_id;
            for (int t = 0; t < layout.times; ++t) {
                for (int i = 0; i < layout.nodes; ++i) {
                    std::vector<int> index;
                    std::vector<double> value;
                    for (int j = 0; j < layout.nodes; ++j) {
                        if (i == j) continue;
                        appendTerm(index, value,
                            layout.x[layout.segmentOffset(d, t, i, j)], 1.0);
                        appendTerm(index, value,
                            layout.x[layout.segmentOffset(d, t, j, i)], -1.0);
                    }
                    const double rhs = (i == source) ? 1.0 : (i == target ? -1.0 : 0.0);
                    model.addRow(
                        "flow_d" + std::to_string(d) +
                        "_t" + std::to_string(t) +
                        "_i" + std::to_string(i),
                        rhs, rhs, std::move(index), std::move(value));
                }
            }
        }

        // At most one incoming and one outgoing segment per node. These valid
        // inequalities strengthen the path formulation and rule out branching.
        for (int d = 0; d < layout.demands; ++d) {
            for (int t = 0; t < layout.times; ++t) {
                for (int i = 0; i < layout.nodes; ++i) {
                    std::vector<int> out_index, in_index;
                    std::vector<double> out_value, in_value;
                    for (int j = 0; j < layout.nodes; ++j) {
                        if (i == j) continue;
                        appendTerm(out_index, out_value,
                            layout.x[layout.segmentOffset(d, t, i, j)], 1.0);
                        appendTerm(in_index, in_value,
                            layout.x[layout.segmentOffset(d, t, j, i)], 1.0);
                    }
                    model.addRow(
                        "outdegree_d" + std::to_string(d) +
                        "_t" + std::to_string(t) +
                        "_i" + std::to_string(i),
                        -kInf, 1.0, std::move(out_index), std::move(out_value));
                    model.addRow(
                        "indegree_d" + std::to_string(d) +
                        "_t" + std::to_string(t) +
                        "_i" + std::to_string(i),
                        -kInf, 1.0, std::move(in_index), std::move(in_value));
                }
            }
        }

        // Maximum number of segments.
        for (int d = 0; d < layout.demands; ++d) {
            for (int t = 0; t < layout.times; ++t) {
                std::vector<int> index;
                std::vector<double> value;
                for (int i = 0; i < layout.nodes; ++i) {
                    for (int j = 0; j < layout.nodes; ++j) {
                        appendTerm(index, value,
                            layout.x[layout.segmentOffset(d, t, i, j)], 1.0);
                    }
                }
                model.addRow(
                    "maxseg_d" + std::to_string(d) +
                    "_t" + std::to_string(t),
                    -kInf, static_cast<double>(instance.max_segments),
                    std::move(index), std::move(value));
            }
        }

        // Absolute-difference linearization and budget.
        for (int t = 1; t < layout.times; ++t) {
            std::vector<int> budget_index;
            std::vector<double> budget_value;
            for (int d = 0; d < layout.demands; ++d) {
                for (int i = 0; i < layout.nodes; ++i) {
                    for (int j = 0; j < layout.nodes; ++j) {
                        if (i == j) continue;
                        const int y = layout.y[layout.segmentOffset(d, t, i, j)];
                        if (y < 0) continue;
                        const int current = layout.x[layout.segmentOffset(d, t, i, j)];
                        const int previous = layout.x[layout.segmentOffset(d, t - 1, i, j)];

                        std::vector<int> index1, index2;
                        std::vector<double> value1, value2;
                        appendTerm(index1, value1, y, 1.0);
                        appendTerm(index1, value1, current, -1.0);
                        appendTerm(index1, value1, previous, 1.0);
                        model.addRow(
                            "budget_pos_d" + std::to_string(d) +
                            "_t" + std::to_string(t) +
                            "_i" + std::to_string(i) +
                            "_j" + std::to_string(j),
                            0.0, kInf, std::move(index1), std::move(value1));

                        appendTerm(index2, value2, y, 1.0);
                        appendTerm(index2, value2, current, 1.0);
                        appendTerm(index2, value2, previous, -1.0);
                        model.addRow(
                            "budget_neg_d" + std::to_string(d) +
                            "_t" + std::to_string(t) +
                            "_i" + std::to_string(i) +
                            "_j" + std::to_string(j),
                            0.0, kInf, std::move(index2), std::move(value2));

                        appendTerm(budget_index, budget_value, y, 1.0);
                    }
                }
            }
            model.addRow(
                "budget_t" + std::to_string(t),
                -kInf, instance.budget_by_time[t],
                std::move(budget_index), std::move(budget_value));
        }
    }

    void addFirstObjective(ModelBuilder& model,
        const InstanceData& instance,
        Layout& layout,
        double incumbent_upper_bound) {

        layout.z = model.addColumn("z", 1.0, 0.0, kInf);

        for (const auto [t, a] : layout.load_position) {
            layout.max_row[layout.loadOffset(t, a)] =
                model.addRow(
                    "maxload_a" + std::to_string(arcNameId(instance, a)) +
                    "_t" + std::to_string(t),
                    0.0, kInf,
                    { layout.z, layout.lambda[layout.loadOffset(t, a)] },
                    { 1.0, -1.0 });
        }

        // If the auxiliary common-routing model produced a feasible solution,
        // its maximum load U is a valid incumbent upper bound for iteration 1.
        // The warm start is optional, though: if no finite U is available, do
        // not add this cutoff and let the actual T-ASR model solve from scratch.
        if (std::isfinite(incumbent_upper_bound)) {
            model.addRow(
                "incumbent_maxload_bound",
                -kInf,
                incumbent_upper_bound,
                { layout.z },
                { 1.0 });
        }
    }

    // When all previous lexicographic levels have been represented by proved
    // fixed bottlenecks, the next level can be solved much more simply: minimize
    // the maximum load among the remaining (unfixed) arc/time loads.
    void addFreeMaxObjective(
        ModelBuilder& model,
        const InstanceData& instance,
        Layout& layout,
        const std::vector<FixedLoad>& fixed_loads) {

        layout.z = model.addColumn("z_free", 1.0, 0.0, kInf);

        int free_loads = 0;
        for (const auto [t, a] : layout.load_position) {
            if (loadIsFixed(fixed_loads, t, a)) continue;
            ++free_loads;

            model.addRow(
                "free_maxload_a" + std::to_string(arcNameId(instance, a)) +
                "_t" + std::to_string(t),
                0.0,
                kInf,
                { layout.z, layout.lambda[layout.loadOffset(t, a)] },
                { 1.0, -1.0 });
        }

        if (free_loads == 0) {
            throw std::runtime_error(
                "cannot optimize another lexicographic level: all loads are fixed");
        }
    }


    // After bottleneck peeling stops, solve the lexicographic problem only on
    // the still-unfixed residual load vector.
    //
    // Suppose m loads L_1,...,L_m were proved individually unavoidable and
    // fixed. The first residual level R_1 = L_{m+1} was already optimized by
    // the last simple-max solve. Therefore R_1 is preserved directly by the
    // upper bounds lambda <= R_1 + tolerance on every unfixed load.
    //
    // For global iteration k > m+1, define residual_k = k-m. We then need only
    // the residual ordered-sum formulation of order residual_k, not the global
    // order-k formulation. Fixed loads are not included in alpha/beta at all.
    void addResidualKSumObjective(
        ModelBuilder& model,
        const InstanceData& instance,
        Layout& layout,
        int residual_k,
        const std::vector<double>& previous_residual_objectives,
        const std::vector<FixedLoad>& fixed_loads) {

        if (residual_k < 2) {
            throw std::logic_error(
                "addResidualKSumObjective requires residual_k >= 2");
        }
        if (previous_residual_objectives.size() <
            static_cast<std::size_t>(residual_k - 1)) {
            throw std::logic_error(
                "missing previous residual lexicographic objective values");
        }

        const int number_of_loads =
            static_cast<int>(layout.load_position.size());

        std::vector<int> free_positions;
        free_positions.reserve(number_of_loads);
        for (int p = 0; p < number_of_loads; ++p) {
            const auto [t, a] = layout.load_position[p];
            if (!loadIsFixed(fixed_loads, t, a)) {
                free_positions.push_back(p);
            }
        }

        if (static_cast<int>(free_positions.size()) < residual_k) {
            throw std::logic_error(
                "residual k-sum model has fewer unfixed loads than residual_k");
        }

        // Indexed by residual lexicographic level. Entries 0 and 1 are unused:
        // residual level 1 is preserved directly by lambda upper bounds.
        layout.alpha.assign(residual_k + 1, -1);
        layout.beta.assign(
            residual_k + 1,
            std::vector<int>(number_of_loads, -1));

        // Create alpha/beta only for residual levels 2..residual_k.
        // Only the current residual level receives objective coefficients.
        for (int iteration = 2; iteration <= residual_k; ++iteration) {
            const double alpha_cost =
                (iteration == residual_k)
                ? static_cast<double>(residual_k)
                : 0.0;

            // alpha remains unrestricted, as in the canonical ordered-sum
            // formulation.
            layout.alpha[iteration] = model.addColumn(
                "alpha_r" + std::to_string(iteration),
                alpha_cost,
                -kInf,
                kInf);

            for (const int p : free_positions) {
                const double beta_cost =
                    (iteration == residual_k) ? 1.0 : 0.0;
                const auto [t, a] = layout.load_position[p];

                layout.beta[iteration][p] = model.addColumn(
                    "beta_r" + std::to_string(iteration) +
                    "_a" + std::to_string(arcNameId(instance, a)) +
                    "_t" + std::to_string(t),
                    beta_cost,
                    0.0,
                    kInf);
            }
        }

        // beta^i_{a,t} >= lambda_{a,t} - alpha_i for UNFIXED loads only.
        for (int iteration = 2; iteration <= residual_k; ++iteration) {
            for (const int p : free_positions) {
                const auto [t, a] = layout.load_position[p];

                model.addRow(
                    "beta_excess_r" + std::to_string(iteration) +
                    "_a" + std::to_string(arcNameId(instance, a)) +
                    "_t" + std::to_string(t),
                    0.0,
                    kInf,
                    { layout.beta[iteration][p],
                      layout.lambda[layout.loadOffset(t, a)],
                      layout.alpha[iteration] },
                    { 1.0, -1.0, 1.0 });
            }
        }

        // Preserve already optimized residual ordered sums for levels
        // 2..residual_k-1:
        //
        //     i * alpha_i + sum_{free p} beta^i_p <= R_i + tolerance.
        //
        // R_1 needs no alpha/beta representation because all unfixed lambdas
        // already have upper bound R_1 + tolerance in addBaseVariables().
        for (int iteration = 2; iteration < residual_k; ++iteration) {
            std::vector<int> index;
            std::vector<double> value;

            appendTerm(
                index,
                value,
                layout.alpha[iteration],
                static_cast<double>(iteration));

            for (const int p : free_positions) {
                appendTerm(
                    index,
                    value,
                    layout.beta[iteration][p],
                    1.0);
            }

            const double rhs =
                previous_residual_objectives.at(
                    static_cast<std::size_t>(iteration - 1)) +
                kLexTolerance;

            model.addRow(
                "preserve_residual_ksum_r" + std::to_string(iteration),
                -kInf,
                rhs,
                std::move(index),
                std::move(value));
        }
    }

    CommonRoutingWarmStart solveCommonRoutingWarmStart(
        const InstanceData& instance,
        const SplitCoefficients& splits,
        const TasrSolverOptions& options,
        double remaining_time_seconds) {

        ModelBuilder builder;
        CommonRoutingLayout layout;
        layout.demands = instance.numDemands();
        layout.times = instance.num_time_slots;
        layout.nodes = instance.numNodes();
        layout.arcs = instance.numArcs();

        // One routing decision per demand and segment. A segment is admitted
        // only if its forwarding graph exists at every real time period,
        // because this route will be copied unchanged to all t.
        const std::size_t segment_count =
            static_cast<std::size_t>(layout.demands) *
            layout.nodes * layout.nodes;
        layout.x.assign(segment_count, -1);

        for (int d = 0; d < layout.demands; ++d) {
            for (int i = 0; i < layout.nodes; ++i) {
                for (int j = 0; j < layout.nodes; ++j) {
                    if (!segmentIsReachableAtAllTimes(
                        splits, layout.times, i, j)) {
                        continue;
                    }

                    layout.x[layout.segmentOffset(d, i, j)] =
                        builder.addColumn(
                            "warm_x_d" + std::to_string(d) +
                            "_i" + std::to_string(i) +
                            "_j" + std::to_string(j),
                            0.0, 0.0, 1.0, VariableType::Binary);
                }
            }
        }

        // Loads stay time-dependent, so the auxiliary model still sees every
        // actual demand volume, intervention state, and split coefficient.
        layout.lambda.assign(
            static_cast<std::size_t>(layout.times) * layout.arcs, -1);

        for (int t = 0; t < layout.times; ++t) {
            for (int a = 0; a < layout.arcs; ++a) {
                if (!arcIsAvailable(instance, t, a)) continue;

                layout.lambda[layout.loadOffset(t, a)] =
                    builder.addColumn(
                        "warm_lambda_a" +
                        std::to_string(arcNameId(instance, a)) +
                        "_t" + std::to_string(t),
                        0.0, 0.0, kInf);

                layout.load_position.emplace_back(t, a);
            }
        }

        // Same physical load equations as in the full model, but with common
        // x[d,i,j]:
        //
        //     capacity[a] * lambda[a,t] - routed_traffic[a,t] = 0.
        for (const auto [t, a] : layout.load_position) {
            std::vector<int> index;
            std::vector<double> value;

            appendTerm(
                index, value,
                layout.lambda[layout.loadOffset(t, a)],
                instance.arc_capacities[a]);

            for (int d = 0; d < layout.demands; ++d) {
                const double volume = instance.demand_matrix[d][t];
                if (volume <= 0.0) continue;

                for (int i = 0; i < layout.nodes; ++i) {
                    for (int j = 0; j < layout.nodes; ++j) {
                        const int x =
                            layout.x[layout.segmentOffset(d, i, j)];
                        if (x < 0) continue;

                        const double r = splits.value(i, j, a, t);
                        appendTerm(index, value, x, -r * volume);
                    }
                }
            }

            builder.addRow(
                "warm_load_a" +
                std::to_string(arcNameId(instance, a)) +
                "_t" + std::to_string(t),
                0.0, 0.0,
                std::move(index), std::move(value));
        }

        // Because the route is common to all t, path constraints are needed
        // only once per demand.
        for (int d = 0; d < layout.demands; ++d) {
            const int source = instance.demands[d].source_id;
            const int target = instance.demands[d].target_id;

            for (int i = 0; i < layout.nodes; ++i) {
                std::vector<int> flow_index;
                std::vector<double> flow_value;
                std::vector<int> out_index;
                std::vector<double> out_value;
                std::vector<int> in_index;
                std::vector<double> in_value;

                for (int j = 0; j < layout.nodes; ++j) {
                    if (i == j) continue;

                    const int outgoing =
                        layout.x[layout.segmentOffset(d, i, j)];
                    const int incoming =
                        layout.x[layout.segmentOffset(d, j, i)];

                    appendTerm(flow_index, flow_value, outgoing, 1.0);
                    appendTerm(flow_index, flow_value, incoming, -1.0);
                    appendTerm(out_index, out_value, outgoing, 1.0);
                    appendTerm(in_index, in_value, incoming, 1.0);
                }

                const double rhs =
                    (i == source) ? 1.0 :
                    (i == target ? -1.0 : 0.0);

                builder.addRow(
                    "warm_flow_d" + std::to_string(d) +
                    "_i" + std::to_string(i),
                    rhs, rhs,
                    std::move(flow_index), std::move(flow_value));

                builder.addRow(
                    "warm_outdegree_d" + std::to_string(d) +
                    "_i" + std::to_string(i),
                    -kInf, 1.0,
                    std::move(out_index), std::move(out_value));

                builder.addRow(
                    "warm_indegree_d" + std::to_string(d) +
                    "_i" + std::to_string(i),
                    -kInf, 1.0,
                    std::move(in_index), std::move(in_value));
            }

            std::vector<int> segment_index;
            std::vector<double> segment_value;

            for (int i = 0; i < layout.nodes; ++i) {
                for (int j = 0; j < layout.nodes; ++j) {
                    appendTerm(
                        segment_index, segment_value,
                        layout.x[layout.segmentOffset(d, i, j)],
                        1.0);
                }
            }

            builder.addRow(
                "warm_maxseg_d" + std::to_string(d),
                -kInf,
                static_cast<double>(instance.max_segments),
                std::move(segment_index),
                std::move(segment_value));
        }

        // Minimize the maximum load over every real arc/time pair.
        layout.z = builder.addColumn("warm_z", 1.0, 0.0, kInf);

        for (const auto [t, a] : layout.load_position) {
            builder.addRow(
                "warm_maxload_a" +
                std::to_string(arcNameId(instance, a)) +
                "_t" + std::to_string(t),
                0.0, kInf,
                { layout.z, layout.lambda[layout.loadOffset(t, a)] },
                { 1.0, -1.0 });
        }

        MipStart direct_start;
        for (int d = 0; d < layout.demands; ++d) {
            const int source = instance.demands[d].source_id;
            const int target = instance.demands[d].target_id;

            const int direct_col =
                layout.x[layout.segmentOffset(d, source, target)];

            if (direct_col < 0) {
                throw std::runtime_error(
                    "cannot warm-start common-routing model for demand " +
                    std::to_string(d) +
                    ": direct source-target segment is not available "
                    "for all time periods");
            }

            for (int i = 0; i < layout.nodes; ++i) {
                for (int j = 0; j < layout.nodes; ++j) {
                    const int col = layout.x[layout.segmentOffset(d, i, j)];
                    if (col < 0) continue;
                    direct_start.index.push_back(col);
                    direct_start.value.push_back(
                        (i == source && j == target) ? 1.0 : 0.0);
                }
            }
        }

        std::cout
            << "Solving common-routing warm-start model with "
            << opt::solverTypeName(options.solver) << ": "
            << builder.numColumns() << " columns, "
            << builder.numRows() << " rows\n";

        BackendOptions backend_options;
        backend_options.time_limit_seconds = remaining_time_seconds;
        backend_options.mip_relative_gap = options.mip_relative_gap;
        backend_options.solver_output = options.solver_output;

        std::string model_file;
        if (options.write_models) {
            std::ostringstream name;
            name << options.model_prefix << "-"
                << opt::solverTypeName(options.solver)
                << "-common-routing-warm-start.lp";
            model_file = name.str();
        }

        const LinearModel model = builder.build();
        const SolverResult solution = opt::solveModel(
            options.solver, model, backend_options, &direct_start, model_file);

        if (!solution.feasible) {
            throw std::runtime_error(
                std::string(opt::solverTypeName(options.solver)) +
                " common-routing warm-start model has no usable solution; status = " +
                solution.status);
        }

        CommonRoutingWarmStart result;
        result.x.assign(segment_count, 0.0);
        result.loads.reserve(layout.load_position.size());
        result.max_load = 0.0;

        for (int d = 0; d < layout.demands; ++d) {
            for (int i = 0; i < layout.nodes; ++i) {
                for (int j = 0; j < layout.nodes; ++j) {
                    const std::size_t offset =
                        layout.segmentOffset(d, i, j);
                    const int col = layout.x[offset];

                    if (col >= 0) {
                        result.x[offset] = solution.values[col];
                    }
                }
            }
        }

        // Compute U from the returned feasible solution itself rather than
        // relying on the reported objective value.
        for (const auto [t, a] : layout.load_position) {
            const double load = solution.values.at(
                layout.lambda[layout.loadOffset(t, a)]);
            result.loads.push_back(load);
            result.max_load = std::max(result.max_load, load);
        }

        std::cout
            << "  common-routing warm-start maximum load U: "
            << std::setprecision(12)
            << result.max_load << '\n';

        return result;
    }

    MipStart makeCommonRoutingMipStart(
        const Layout& full_layout,
        const std::vector<double>& common_x) {

        const std::size_t expected_size =
            static_cast<std::size_t>(full_layout.demands) *
            full_layout.nodes * full_layout.nodes;

        if (common_x.size() != expected_size) {
            throw std::logic_error(
                "common-routing MIP-start vector has wrong size");
        }

        MipStart start;
        for (int d = 0; d < full_layout.demands; ++d) {
            for (int t = 0; t < full_layout.times; ++t) {
                for (int i = 0; i < full_layout.nodes; ++i) {
                    for (int j = 0; j < full_layout.nodes; ++j) {
                        const int col = full_layout.x[
                            full_layout.segmentOffset(d, t, i, j)];
                        if (col < 0) continue;

                        const std::size_t common_offset =
                            ((static_cast<std::size_t>(d) *
                                full_layout.nodes + i) *
                                full_layout.nodes + j);

                        start.index.push_back(col);
                        start.value.push_back(
                            common_x[common_offset] > 0.5 ? 1.0 : 0.0);
                    }
                }
            }
        }
        return start;
    }

    std::vector<double> loadValues(const Layout& layout,
        const SolverResult& solution) {
        std::vector<double> loads;
        loads.reserve(layout.load_position.size());
        for (const auto [t, a] : layout.load_position) {
            loads.push_back(solution.values.at(layout.lambda[layout.loadOffset(t, a)]));
        }
        return loads;
    }

    // Compute the physical loads implied by x. This is used to select the
    // bottleneck candidate for the fixing test. It avoids accidentally selecting
    // a nonbinding lambda variable that merely contains slack.
    std::vector<double> actualLoadValues(
        const InstanceData& instance,
        const SplitCoefficients& splits,
        const Layout& layout,
        const SolverResult& solution) {

        std::vector<double> loads;
        loads.reserve(layout.load_position.size());

        for (const auto [t, a] : layout.load_position) {
            const double capacity = instance.arc_capacities[a];
            if (capacity <= 0.0) {
                throw std::runtime_error(
                    "non-positive capacity while computing actual load");
            }

            double traffic = 0.0;
            for (int d = 0; d < layout.demands; ++d) {
                const double volume = instance.demand_matrix[d][t];
                if (volume <= 0.0) continue;

                for (int i = 0; i < layout.nodes; ++i) {
                    for (int j = 0; j < layout.nodes; ++j) {
                        const int x =
                            layout.x[layout.segmentOffset(d, t, i, j)];
                        if (x < 0) continue;

                        const double x_value = solution.values.at(x);
                        if (std::abs(x_value) <= 1e-12) continue;

                        const double r = splits.value(i, j, a, t);
                        traffic += r * volume * x_value;
                    }
                }
            }

            loads.push_back(traffic / capacity);
        }

        return loads;
    }

    void reportInvariantDifference(
        const std::string& label,
        double lhs,
        double rhs,
        double warning_tolerance = kInvariantWarningTolerance,
        double error_tolerance = kInvariantErrorTolerance) {

        const double difference = std::abs(lhs - rhs);
        if (difference > error_tolerance) {
            std::ostringstream message;
            message << std::setprecision(16)
                << "MODEL INVARIANT ERROR: " << label
                << ": lhs=" << lhs
                << ", rhs=" << rhs
                << ", |diff|=" << difference
                << " > " << error_tolerance;
            throw std::runtime_error(message.str());
        }

        if (difference > warning_tolerance) {
            std::cerr << std::setprecision(16)
                << "WARNING: " << label
                << ": lhs=" << lhs
                << ", rhs=" << rhs
                << ", |diff|=" << difference
                << " > " << warning_tolerance << '\n';
        }
    }

    void verifyPhysicalLoadConsistency(
        const InstanceData& instance,
        const SplitCoefficients& splits,
        const Layout& layout,
        const SolverResult& solution,
        const std::vector<FixedLoad>& fixed_loads,
        bool simple_max_model,
        double expected_free_max,
        int lexicographic_iteration,
        const char* phase) {

        const std::vector<double> physical =
            actualLoadValues(instance, splits, layout, solution);

        if (physical.size() != layout.load_position.size()) {
            throw std::logic_error(
                "physical-load diagnostic vector has wrong size");
        }

        double largest_free_physical = -kInf;
        double largest_free_lambda = -kInf;

        int largest_free_physical_t = -1;
        int largest_free_physical_a = -1;

        int largest_free_lambda_t = -1;
        int largest_free_lambda_a = -1;
        int largest_free_lambda_position = -1;

        for (int p = 0;
            p < static_cast<int>(layout.load_position.size());
            ++p) {

            const auto [t, a] = layout.load_position[p];

            const int lambda_col =
                layout.lambda[layout.loadOffset(t, a)];

            if (lambda_col < 0 ||
                lambda_col >= static_cast<int>(solution.values.size())) {

                throw std::logic_error(
                    "load diagnostic references invalid lambda column");
            }

            const double lambda_value =
                solution.values[lambda_col];

            const double physical_value =
                physical[p];

            // ---------------------------------------------------------
            // 1. lambda must equal the actual physical utilization.
            // ---------------------------------------------------------
            //
            // Base load row:
            //
            //     capacity * lambda - traffic = 0
            //
            // Thus any material lambda/physical discrepancy is a model or
            // numerical-consistency problem.
            reportInvariantDifference(
                std::string(phase) +
                " lambda vs physical load (t=" +
                std::to_string(t) +
                ",a=" +
                std::to_string(arcNameId(instance, a)) +
                ")",
                lambda_value,
                physical_value);

            // ---------------------------------------------------------
            // 2. Verify already fixed physical bottlenecks.
            // ---------------------------------------------------------
            if (const auto fixed =
                fixedLoadValue(fixed_loads, t, a)) {

                reportInvariantDifference(
                    std::string(phase) +
                    " fixed physical bottleneck (t=" +
                    std::to_string(t) +
                    ",a=" +
                    std::to_string(
                        arcNameId(instance, a)) +
                    ")",
                    physical_value,
                    *fixed);
            }
            else {
                // -----------------------------------------------------
                // 3. Track largest still-free physical load.
                // -----------------------------------------------------
                if (physical_value > largest_free_physical) {
                    largest_free_physical =
                        physical_value;

                    largest_free_physical_t = t;
                    largest_free_physical_a = a;
                }

                // -----------------------------------------------------
                // 4. Track largest still-free lambda value.
                // -----------------------------------------------------
                if (lambda_value > largest_free_lambda) {
                    largest_free_lambda =
                        lambda_value;

                    largest_free_lambda_t = t;
                    largest_free_lambda_a = a;
                    largest_free_lambda_position = p;
                }
            }
        }

        // -------------------------------------------------------------
        // 5. Diagnostics for the persistent simple-max formulation.
        // -------------------------------------------------------------
        if (simple_max_model &&
            std::isfinite(largest_free_physical) &&
            std::isfinite(largest_free_lambda)) {

            const double physical_at_largest_lambda =
                largest_free_lambda_position >= 0
                ? physical[largest_free_lambda_position]
                : kInf;

            std::cout
                << std::setprecision(16)
                << "\n"
                << "  PERSISTENT MAX DIAGNOSTIC"
                << " [iteration "
                << lexicographic_iteration
                << ", phase=" << phase << "]\n"

                << "    solver objective z            = "
                << expected_free_max << "\n"

                << "    largest unfixed lambda        = "
                << largest_free_lambda
                << " at (t="
                << largest_free_lambda_t
                << ", a="
                << arcNameId(
                    instance,
                    largest_free_lambda_a)
                << ")\n"

                << "    physical load at that arc     = "
                << physical_at_largest_lambda
                << "\n"

                << "    largest unfixed physical load = "
                << largest_free_physical
                << " at (t="
                << largest_free_physical_t
                << ", a="
                << arcNameId(
                    instance,
                    largest_free_physical_a)
                << ")\n";

            if (layout.z >= 0 &&
                layout.z <
                static_cast<int>(solution.values.size())) {

                std::cout
                    << "    z variable value              = "
                    << solution.values[layout.z]
                    << "\n";
            }

            std::cout
                << "    z - max(lambda)               = "
                << (expected_free_max -
                    largest_free_lambda)
                << "\n"

                << "    max(lambda) - max(physical)   = "
                << (largest_free_lambda -
                    largest_free_physical)
                << "\n"

                << "    z - max(physical)             = "
                << (expected_free_max -
                    largest_free_physical)
                << "\n\n";

            // ---------------------------------------------------------
            // 6. Existing hard invariant:
            //    objective should equal largest free physical load.
            // ---------------------------------------------------------
            reportInvariantDifference(
                std::string(phase) +
                " solver objective vs largest physical unfixed load",
                expected_free_max,
                largest_free_physical);

            // ---------------------------------------------------------
            // 7. New invariant:
            //    objective should also equal largest free lambda.
            //
            // If this one passes but the physical one fails, then some
            // free lambda contains artificial slack.
            //
            // If this one fails too, then z is constrained by something
            // other than the active free-lambda max rows.
            // ---------------------------------------------------------
            reportInvariantDifference(
                std::string(phase) +
                " solver objective vs largest unfixed lambda",
                expected_free_max,
                largest_free_lambda);

            // ---------------------------------------------------------
            // 8. Existing z-variable consistency check.
            // ---------------------------------------------------------
            if (layout.z >= 0 &&
                layout.z <
                static_cast<int>(solution.values.size())) {

                reportInvariantDifference(
                    std::string(phase) +
                    " solver objective vs z variable",
                    expected_free_max,
                    solution.values[layout.z]);
            }

            // ---------------------------------------------------------
            // 9. Additional diagnostic:
            //    with physical load equalities, the largest free lambda must
            //    correspond to the same physical utilization (up to numerical
            //    feasibility tolerances).
            // ---------------------------------------------------------
            const double slack_at_max_lambda =
                largest_free_lambda -
                physical_at_largest_lambda;

            if (slack_at_max_lambda >
                kInvariantWarningTolerance) {

                std::cerr
                    << std::setprecision(16)
                    << "WARNING in " << phase
                    << " at lexicographic iteration "
                    << lexicographic_iteration
                    << ": largest unfixed lambda differs from physical utilization."
                    << " lambda=" << largest_free_lambda
                    << ", physical="
                    << physical_at_largest_lambda
                    << ", slack="
                    << slack_at_max_lambda
                    << ", at (t="
                    << largest_free_lambda_t
                    << ", a="
                    << arcNameId(
                        instance,
                        largest_free_lambda_a)
                    << ")\n";
            }
        }
    }

    MipStart makeBaseWarmStart(
        const Layout& layout,
        const Layout& previous_layout,
        const SolverResult& previous_solution) {

        MipStart warm;
        auto addStart = [&](int column, double value) {
            if (column < 0) return;
            warm.index.push_back(column);
            warm.value.push_back(value);
            };

        for (std::size_t offset = 0; offset < layout.x.size(); ++offset) {
            const int current_col = layout.x[offset];
            const int previous_col = previous_layout.x[offset];
            if (current_col >= 0 && previous_col >= 0 &&
                previous_col < static_cast<int>(previous_solution.values.size())) {
                addStart(current_col, previous_solution.values[previous_col]);
            }
        }

        for (std::size_t offset = 0; offset < layout.lambda.size(); ++offset) {
            const int current_col = layout.lambda[offset];
            const int previous_col = previous_layout.lambda[offset];
            if (current_col >= 0 && previous_col >= 0 &&
                previous_col < static_cast<int>(previous_solution.values.size())) {
                addStart(current_col, previous_solution.values[previous_col]);
            }
        }

        for (std::size_t offset = 0; offset < layout.y.size(); ++offset) {
            const int current_col = layout.y[offset];
            const int previous_col = previous_layout.y[offset];
            if (current_col >= 0 && previous_col >= 0 &&
                previous_col < static_cast<int>(previous_solution.values.size())) {
                addStart(current_col, previous_solution.values[previous_col]);
            }
        }

        return warm;
    }


    MipStart makeResidualWarmStart(
        const Layout& layout,
        const Layout& previous_layout,
        const SolverResult& previous_solution,
        const std::vector<double>& previous_loads,
        int residual_k,
        const std::vector<FixedLoad>& fixed_loads) {

        MipStart warm;
        auto addStart = [&](int column, double value) {
            if (column < 0) return;
            warm.index.push_back(column);
            warm.value.push_back(value);
            };

        // Cutoff preprocessing can change x-column numbering between
        // lexicographic iterations, so values are copied by semantic offset.
        for (std::size_t offset = 0; offset < layout.x.size(); ++offset) {
            const int current_col = layout.x[offset];
            const int previous_col = previous_layout.x[offset];
            if (current_col >= 0 && previous_col >= 0 &&
                previous_col < static_cast<int>(previous_solution.values.size())) {
                addStart(current_col, previous_solution.values[previous_col]);
            }
        }

        for (std::size_t offset = 0; offset < layout.lambda.size(); ++offset) {
            const int current_col = layout.lambda[offset];
            const int previous_col = previous_layout.lambda[offset];
            if (current_col >= 0 && previous_col >= 0 &&
                previous_col < static_cast<int>(previous_solution.values.size())) {
                addStart(current_col, previous_solution.values[previous_col]);
            }
        }

        for (std::size_t offset = 0; offset < layout.y.size(); ++offset) {
            const int current_col = layout.y[offset];
            const int previous_col = previous_layout.y[offset];
            if (current_col >= 0 && previous_col >= 0 &&
                previous_col < static_cast<int>(previous_solution.values.size())) {
                addStart(current_col, previous_solution.values[previous_col]);
            }
        }

        // Build the previous residual load vector in layout.load_position order.
        std::vector<double> free_previous_loads;
        std::vector<int> free_positions;
        free_previous_loads.reserve(layout.load_position.size());
        free_positions.reserve(layout.load_position.size());

        for (int p = 0;
            p < static_cast<int>(layout.load_position.size());
            ++p) {

            const auto [t, a] = layout.load_position[p];
            if (loadIsFixed(fixed_loads, t, a)) {
                continue;
            }
            if (p >= static_cast<int>(previous_loads.size())) {
                throw std::logic_error(
                    "previous load vector is too short for residual warm start");
            }

            free_positions.push_back(p);
            free_previous_loads.push_back(previous_loads[p]);
        }

        if (residual_k >= 2 && !free_previous_loads.empty()) {
            std::vector<double> sorted = free_previous_loads;
            std::sort(sorted.begin(), sorted.end(), std::greater<double>());

            for (int iteration = 2; iteration <= residual_k; ++iteration) {
                const int sorted_index = std::min<int>(
                    iteration - 1,
                    static_cast<int>(sorted.size()) - 1);
                const double alpha = sorted[sorted_index];

                addStart(layout.alpha[iteration], alpha);

                for (int q = 0;
                    q < static_cast<int>(free_positions.size());
                    ++q) {

                    const int p = free_positions[q];
                    addStart(
                        layout.beta[iteration][p],
                        std::max(0.0, free_previous_loads[q] - alpha));
                }
            }
        }

        return warm;
    }

    struct BottleneckTestResult {
        bool solved_optimally = false;
        bool proved_fixed = false;
        double minimum_candidate_load = kInf;
    };

    // Test whether one particular current bottleneck is unavoidable.
    //
    // Earlier proved bottlenecks are fixed. Every still-free load is bounded by
    // current_boundary. We then minimize the candidate load itself. If its proven
    // optimum cannot go below current_boundary, that arc/time load is unavoidable
    // and can safely be fixed for all subsequent lexicographic levels.
    BottleneckTestResult testBottleneckCandidate(
        const InstanceData& instance,
        const SplitCoefficients& splits,
        const TasrSolverOptions& options,
        double remaining_time_seconds,
        double global_max_cutoff,
        double current_boundary,
        const std::vector<FixedLoad>& fixed_loads,
        int candidate_time,
        int candidate_arc,
        const Layout& incumbent_layout,
        const SolverResult& incumbent_solution,
        int lexicographic_iteration) {

        BottleneckTestResult result;

        ModelBuilder builder;
        Layout layout;

        addBaseVariables(
            builder,
            instance,
            splits,
            layout,
            global_max_cutoff,
            current_boundary,
            fixed_loads);
        addBaseConstraints(builder, instance, splits, layout, fixed_loads);

        const int candidate_col =
            layout.lambda[layout.loadOffset(candidate_time, candidate_arc)];
        if (candidate_col < 0) {
            throw std::logic_error(
                "bottleneck-test candidate load does not exist in the model");
        }

        // ModelBuilder has no objective-coefficient mutator. Introduce a single
        // objective variable equal to the candidate lambda.
        const int test_objective =
            builder.addColumn("bottleneck_test_objective", 1.0, 0.0, kInf);

        builder.addRow(
            "bottleneck_test_link",
            0.0,
            0.0,
            { test_objective, candidate_col },
            { 1.0, -1.0 });

        BackendOptions backend_options;
        backend_options.time_limit_seconds = remaining_time_seconds;
        backend_options.mip_relative_gap = options.mip_relative_gap;
        backend_options.solver_output = options.solver_output;

        MipStart warm_start =
            makeBaseWarmStart(layout, incumbent_layout, incumbent_solution);

        std::string model_file;
        if (options.write_models) {
            std::ostringstream name;
            name << options.model_prefix << "-"
                << opt::solverTypeName(options.solver)
                << "-bottleneck-test-" << lexicographic_iteration << ".lp";
            model_file = name.str();
        }

        std::cout
            << "  Testing whether bottleneck (t=" << candidate_time
            << ", a=" << arcNameId(instance, candidate_arc)
            << ") is unavoidable at load "
            << std::setprecision(12) << current_boundary << "...\n";

        const LinearModel model = builder.build();
        const SolverResult solution = opt::solveModel(
            options.solver,
            model,
            backend_options,
            &warm_start,
            model_file);

        if (!solution.feasible) {
            std::cout
                << "    bottleneck test returned no feasible solution (status "
                << solution.status << "); no fixing is performed.\n";
            return result;
        }

        result.minimum_candidate_load = solution.objective;
        result.solved_optimally = solution.optimal;

        // The temporary objective is supposed to be exactly candidate lambda.
        // Verify that no stale objective coefficient survived a previous
        // persistent solve.
        const double candidate_lambda_value =
            solution.values.at(static_cast<std::size_t>(candidate_col));
        reportInvariantDifference(
            "bottleneck-test ObjVal vs candidate lambda",
            solution.objective,
            candidate_lambda_value);

        // Because candidate lambda is minimized and its load row is
        // capacity*lambda - traffic >= 0, it must be tight at the helper optimum.
        // Recompute the utilization directly from x and verify this explicitly.
        const std::vector<double> helper_physical =
            actualLoadValues(instance, splits, layout, solution);
        int candidate_position = -1;
        for (int p = 0;
            p < static_cast<int>(layout.load_position.size());
            ++p) {
            if (layout.load_position[p].first == candidate_time &&
                layout.load_position[p].second == candidate_arc) {
                candidate_position = p;
                break;
            }
        }
        if (candidate_position < 0) {
            throw std::logic_error(
                "bottleneck-test candidate missing from load-position vector");
        }

        const double candidate_physical_load =
            helper_physical[static_cast<std::size_t>(candidate_position)];

        reportInvariantDifference(
            "bottleneck-test candidate lambda vs physical load from x",
            candidate_lambda_value,
            candidate_physical_load);

        std::cout
            << std::setprecision(16)
            << "    proof diagnostics: boundary=" << current_boundary
            << ", helper ObjVal=" << solution.objective
            << ", helper lambda=" << candidate_lambda_value
            << ", physical load=" << candidate_physical_load
            << ", boundary-helper="
            << (current_boundary - solution.objective)
            << '\n';

        if (!solution.optimal) {
            std::cout
                << "    bottleneck test stopped with status " << solution.status
                << " at candidate load " << std::setprecision(12)
                << solution.objective
                << "; uniqueness is not proved.\n";
            return result;
        }

        const double proof_discrepancy =
            current_boundary - solution.objective;

        if (std::abs(proof_discrepancy) > kBottleneckProofTolerance &&
            solution.objective >= current_boundary - kTolerance) {
            std::cerr << std::setprecision(16)
                << "WARNING: bottleneck would be accepted by the existing "
                << "kTolerance=" << kTolerance
                << " rule, but helper differs from the boundary by "
                << proof_discrepancy
                << ", which exceeds diagnostic proof tolerance "
                << kBottleneckProofTolerance << ".\n";
        }

        if (solution.objective > current_boundary + kInvariantErrorTolerance) {
            throw std::runtime_error(
                "bottleneck helper optimum is materially ABOVE the feasible "
                "current boundary; helper/main models are inconsistent");
        }

        result.proved_fixed =
            solution.objective >= current_boundary - kTolerance;

        if (result.proved_fixed) {
            std::cout
                << "    proved unavoidable: minimum candidate load = "
                << std::setprecision(12) << solution.objective << ".\n";
        }
        else {
            std::cout
                << "    candidate is not individually unavoidable: it can be reduced to "
                << std::setprecision(12) << solution.objective
                << " while all other free loads remain <= "
                << current_boundary << ".\n";
        }

        return result;
    }



    BottleneckTestResult testBottleneckCandidatePersistent(
        opt::SolverSession& session,
        const InstanceData& instance,
        const SplitCoefficients& splits,
        const TasrSolverOptions& options,
        double remaining_time_seconds,
        double current_boundary,
        const std::vector<FixedLoad>& fixed_loads,
        int candidate_time,
        int candidate_arc,
        const Layout& layout,
        int lexicographic_iteration) {

        BottleneckTestResult result;

        // Preserve the just-proved lexicographic level by tightening every
        // still-free load to the current boundary. Already fixed bottlenecks
        // retain their own narrow [L_i-eps, L_i+eps] bounds.
        for (const auto [t, a] : layout.load_position) {
            if (loadIsFixed(fixed_loads, t, a)) continue;
            const int lambda =
                layout.lambda[layout.loadOffset(t, a)];
            session.setColumnBounds(
                lambda,
                0.0,
                current_boundary + kLexFixingTolerance);
        }

        const int candidate_col =
            layout.lambda[layout.loadOffset(candidate_time, candidate_arc)];
        if (candidate_col < 0) {
            throw std::logic_error(
                "persistent bottleneck-test candidate does not exist");
        }

        // Temporarily change the objective from min z_free to min lambda_b.
        // The same native solver model, presolve information and incumbent state
        // are retained; no T-ASR rows or columns are reconstructed.
        session.setObjectiveCoefficient(layout.z, 0.0);
        session.setObjectiveCoefficient(candidate_col, 1.0);

        BackendOptions backend_options;
        backend_options.time_limit_seconds = remaining_time_seconds;
        backend_options.mip_relative_gap = options.mip_relative_gap;
        backend_options.solver_output = options.solver_output;

        // One feasible helper solution strictly below this threshold is already
        // enough to disprove unavoidability.  Ask supported backends (currently
        // Gurobi) to terminate immediately when such an incumbent appears.
        backend_options.stop_when_incumbent_below = true;
        backend_options.incumbent_stop_value =
            current_boundary - kTolerance;

        std::string model_file;
        if (options.write_models) {
            std::ostringstream name;
            name << options.model_prefix << "-"
                << opt::solverTypeName(options.solver)
                << "-bottleneck-test-" << lexicographic_iteration << ".lp";
            model_file = name.str();
        }

        std::cout
            << "  Re-optimizing persistent model to test bottleneck (t="
            << candidate_time
            << ", a=" << arcNameId(instance, candidate_arc)
            << ") at load "
            << std::setprecision(12) << current_boundary << "...\n";

        SolverResult solution;
        try {
            solution = session.optimize(
                backend_options,
                nullptr,
                model_file);
        }
        catch (...) {
            // Always restore the normal peeled min-max objective.
            session.setObjectiveCoefficient(candidate_col, 0.0);
            session.setObjectiveCoefficient(layout.z, 1.0);
            throw;
        }

        // Restore the normal objective for the next lexicographic iteration.
        session.setObjectiveCoefficient(candidate_col, 0.0);
        session.setObjectiveCoefficient(layout.z, 1.0);

        if (!solution.feasible) {
            std::cout
                << "    persistent bottleneck test returned no feasible solution "
                << "(status " << solution.status
                << "); no fixing is performed.\n";
            return result;
        }

        result.minimum_candidate_load = solution.objective;
        result.solved_optimally = solution.optimal;

        // The temporary persistent objective must be exactly candidate lambda.
        // This catches stale objective coefficients left behind by earlier
        // re-optimizations.
        const double candidate_lambda_value =
            solution.values.at(static_cast<std::size_t>(candidate_col));
        reportInvariantDifference(
            "persistent bottleneck-test ObjVal vs candidate lambda",
            solution.objective,
            candidate_lambda_value);

        // Recompute the helper candidate's physical utilization directly from x.
        // Minimizing lambda_b should make its load inequality tight.
        const std::vector<double> helper_physical =
            actualLoadValues(instance, splits, layout, solution);
        int candidate_position = -1;
        for (int p = 0;
            p < static_cast<int>(layout.load_position.size());
            ++p) {
            if (layout.load_position[p].first == candidate_time &&
                layout.load_position[p].second == candidate_arc) {
                candidate_position = p;
                break;
            }
        }
        if (candidate_position < 0) {
            throw std::logic_error(
                "persistent bottleneck-test candidate missing from load-position vector");
        }

        const double candidate_physical_load =
            helper_physical[static_cast<std::size_t>(candidate_position)];

        reportInvariantDifference(
            "persistent bottleneck-test candidate lambda vs physical load from x",
            candidate_lambda_value,
            candidate_physical_load);

        std::cout
            << std::setprecision(16)
            << "    proof diagnostics: boundary=" << current_boundary
            << ", helper ObjVal=" << solution.objective
            << ", helper lambda=" << candidate_lambda_value
            << ", physical load=" << candidate_physical_load
            << ", boundary-helper="
            << (current_boundary - solution.objective)
            << '\n';

        if (solution.incumbent_stop_triggered) {
            std::cout
                << "    early helper stop: found feasible candidate load "
                << std::setprecision(12) << solution.objective
                << " < rejection threshold "
                << backend_options.incumbent_stop_value
                << ". Bottleneck is NOT unavoidable; no optimality proof needed.\n";
            return result;
        }

        if (!solution.optimal) {
            std::cout
                << "    persistent bottleneck test stopped with status "
                << solution.status
                << " at candidate load "
                << std::setprecision(12) << solution.objective
                << "; unavoidability remains unknown.\n";
            return result;
        }

        const double proof_discrepancy =
            current_boundary - solution.objective;

        if (std::abs(proof_discrepancy) > kBottleneckProofTolerance &&
            solution.objective >= current_boundary - kTolerance) {
            std::cerr << std::setprecision(16)
                << "WARNING: bottleneck would be accepted by the existing "
                << "kTolerance=" << kTolerance
                << " rule, but helper differs from the boundary by "
                << proof_discrepancy
                << ", which exceeds diagnostic proof tolerance "
                << kBottleneckProofTolerance << ".\n";
        }

        if (solution.objective > current_boundary + kInvariantErrorTolerance) {
            throw std::runtime_error(
                "bottleneck helper optimum is materially ABOVE the feasible "
                "current boundary; helper/main models are inconsistent");
        }

        result.proved_fixed =
            solution.objective >= current_boundary - kTolerance;

        if (result.proved_fixed) {
            std::cout
                << "    proved unavoidable: minimum candidate load = "
                << std::setprecision(12) << solution.objective << ".\n";
        }
        else {
            std::cout
                << "    candidate is not individually unavoidable: it can be "
                << "reduced to " << std::setprecision(12)
                << solution.objective
                << " while all other free loads remain <= "
                << current_boundary << ".\n";
        }

        return result;
    }


    void prunePersistentSegments(
        opt::SolverSession& session,
        const InstanceData& instance,
        const SplitCoefficients& splits,
        const Layout& layout,
        double global_max_cutoff) {

        if (!std::isfinite(global_max_cutoff)) return;

        std::size_t fixed_to_zero = 0;
        for (int d = 0; d < layout.demands; ++d) {
            for (int t = 0; t < layout.times; ++t) {
                for (int i = 0; i < layout.nodes; ++i) {
                    for (int j = 0; j < layout.nodes; ++j) {
                        const int col =
                            layout.x[layout.segmentOffset(d, t, i, j)];
                        if (col < 0) continue;

                        if (segmentLoadLowerBound(
                            instance, splits, d, t, i, j) >
                            global_max_cutoff + kLexFixingTolerance) {
                            session.setColumnBounds(col, 0.0, 0.0);
                            ++fixed_to_zero;
                        }
                    }
                }
            }
        }

        if (fixed_to_zero > 0) {
            std::cout
                << "  persistent preprocessing fixed "
                << fixed_to_zero
                << " additional x variables to zero using z1.\n";
        }
    }


    void fixPersistentBottleneck(
        opt::SolverSession& session,
        const Layout& layout,
        int time,
        int arc,
        double value) {

        const std::size_t offset = layout.loadOffset(time, arc);
        const int lambda = layout.lambda[offset];
        const int max_row = layout.max_row[offset];

        if (lambda < 0 || max_row < 0) {
            throw std::logic_error(
                "persistent bottleneck references missing model entries");
        }

        // The helper has proved that this physical arc/time load is unavoidable.
        // Keep a small numerical window around the solver-produced value.
        session.setColumnBounds(
            lambda,
            std::max(0.0, value - kLexFixingTolerance),
            value + kLexFixingTolerance);

        // Every load row is already an equality:
        //
        //     capacity * lambda - routed_traffic = 0.
        //
        // Therefore no load-row modification is required here.
        //
        // Peel this bottleneck from the next max objective by changing
        //
        //     z - lambda >= 0
        //
        // into
        //
        //     z >= 0.
        //
        // Removing only the lambda coefficient keeps the native row numbering
        // stable and avoids using GRB_INFINITY as an artificial RHS.
        session.setCoefficient(max_row, lambda, 0.0);
    }


    std::vector<std::vector<std::vector<int>>> extractWaypoints(
        const InstanceData& instance,
        const Layout& layout,
        const SolverResult& solution) {
        std::vector<std::vector<std::vector<int>>> result(
            layout.demands,
            std::vector<std::vector<int>>(layout.times));

        for (int d = 0; d < layout.demands; ++d) {
            for (int t = 0; t < layout.times; ++t) {
                std::vector<int> successor(layout.nodes, -1);
                for (int i = 0; i < layout.nodes; ++i) {
                    for (int j = 0; j < layout.nodes; ++j) {
                        const int col = layout.x[layout.segmentOffset(d, t, i, j)];
                        if (col >= 0 && solution.values[col] > 0.5) {
                            if (successor[i] >= 0 && successor[i] != j) {
                                throw std::runtime_error("solution contains a branching segment path");
                            }
                            successor[i] = j;
                        }
                    }
                }

                const int source = instance.demands[d].source_id;
                const int target = instance.demands[d].target_id;
                int current = source;
                std::vector<unsigned char> visited(layout.nodes, 0);
                for (int step = 0; step <= instance.max_segments; ++step) {
                    if (current == target) break;
                    if (current < 0 || current >= layout.nodes || visited[current]) {
                        throw std::runtime_error("could not extract a simple source-target segment path");
                    }
                    visited[current] = 1;
                    const int next = successor[current];
                    if (next < 0) {
                        throw std::runtime_error("selected segments do not connect demand source to target");
                    }
                    current = next;
                    if (current != target) result[d][t].push_back(current);
                }
                if (current != target) {
                    throw std::runtime_error("segment path exceeds max_segments or contains a cycle");
                }
            }
        }
        return result;
    }

    std::vector<std::vector<std::vector<int>>> extractCommonRoutingWaypoints(
        const InstanceData& instance,
        const std::vector<double>& common_x) {

        const int demands = instance.numDemands();
        const int times = instance.num_time_slots;
        const int nodes = instance.numNodes();
        const std::size_t expected =
            static_cast<std::size_t>(demands) * nodes * nodes;

        if (common_x.size() != expected) {
            throw std::logic_error(
                "common-routing solution has wrong x-vector size");
        }

        std::vector<std::vector<std::vector<int>>> result(
            demands, std::vector<std::vector<int>>(times));

        auto offset = [nodes](int d, int i, int j) {
            return ((static_cast<std::size_t>(d) * nodes + i) * nodes + j);
            };

        for (int d = 0; d < demands; ++d) {
            std::vector<int> successor(nodes, -1);
            for (int i = 0; i < nodes; ++i) {
                for (int j = 0; j < nodes; ++j) {
                    if (common_x[offset(d, i, j)] > 0.5) {
                        if (successor[i] >= 0 && successor[i] != j) {
                            throw std::runtime_error(
                                "common-routing solution contains branching");
                        }
                        successor[i] = j;
                    }
                }
            }

            std::vector<int> waypoints;
            const int source = instance.demands[d].source_id;
            const int target = instance.demands[d].target_id;
            int current = source;
            std::vector<unsigned char> visited(nodes, 0);

            for (int step = 0; step <= instance.max_segments; ++step) {
                if (current == target) break;
                if (current < 0 || current >= nodes || visited[current]) {
                    throw std::runtime_error(
                        "could not extract common-routing segment path");
                }
                visited[current] = 1;
                const int next = successor[current];
                if (next < 0) {
                    throw std::runtime_error(
                        "common-routing segments do not connect source to target");
                }
                current = next;
                if (current != target) waypoints.push_back(current);
            }

            if (current != target) {
                throw std::runtime_error(
                    "common-routing path exceeds max_segments or contains a cycle");
            }

            for (int t = 0; t < times; ++t) {
                result[d][t] = waypoints;
            }
        }

        return result;
    }

    bool lexicographicallyBetterLoads(
        const std::vector<double>& candidate,
        const std::vector<double>& incumbent) {

        if (incumbent.empty()) return true;
        if (candidate.size() != incumbent.size()) {
            throw std::logic_error("load vectors have different sizes");
        }

        std::vector<double> lhs = candidate;
        std::vector<double> rhs = incumbent;
        std::sort(lhs.begin(), lhs.end(), std::greater<double>());
        std::sort(rhs.begin(), rhs.end(), std::greater<double>());

        for (std::size_t i = 0; i < lhs.size(); ++i) {
            if (lhs[i] < rhs[i] - kTolerance) return true;
            if (lhs[i] > rhs[i] + kTolerance) return false;
        }
        return false;
    }

} // namespace

TasrSolver::TasrSolver(const InstanceData& instance,
    const SplitCoefficients& split_coefficients,
    TasrSolverOptions options)
    : instance_(instance), splits_(split_coefficients), options_(std::move(options)) {
}

TasrSolution TasrSolver::solve() {
    using Clock = std::chrono::steady_clock;
    const auto process_start = Clock::now();

    auto elapsedSeconds = [&]() {
        return std::chrono::duration<double>(
            Clock::now() - process_start).count();
        };

    auto remainingSeconds = [&]() {
        if (options_.time_limit_seconds <= 0.0) {
            return 0.0; // backend convention: unlimited
        }
        return std::max(
            0.0,
            options_.time_limit_seconds - elapsedSeconds());
        };

    auto budgetExpired = [&]() {
        return options_.time_limit_seconds > 0.0 &&
            remainingSeconds() <= 0.0;
        };

    const int total_loads = [&]() {
        int count = 0;
        for (int t = 0; t < instance_.num_time_slots; ++t) {
            for (int a = 0; a < instance_.numArcs(); ++a) {
                count += arcIsAvailable(instance_, t, a) ? 1 : 0;
            }
        }
        return count;
        }();

    if (total_loads == 0) {
        throw std::runtime_error(
            "instance has no available arc/time loads");
    }

    int requested_iterations = options_.max_lexicographic_iterations;
    if (requested_iterations < 0) {
        requested_iterations = total_loads;
    }
    requested_iterations = std::max(
        1,
        std::min(requested_iterations, total_loads));

    std::vector<double> previous_objectives;
    std::vector<double> lexicographic_load_values;

    std::vector<double> previous_loads;
    Layout previous_layout;
    SolverResult previous_solution;
    bool have_previous_solution = false;
    int completed_optimal_iterations = 0;

    std::vector<FixedLoad> fixed_bottlenecks;
    bool bottleneck_peeling_active = true;

    // Once individual bottleneck peeling fails, all subsequently optimized
    // levels belong to the residual vector of still-unfixed loads. The first
    // residual maximum was already obtained by the final simple-max solve.
    bool residual_ksum_active = false;
    int residual_fixed_prefix = -1;
    double residual_prefix_objective = 0.0;
    double residual_max_boundary = kInf;

    Layout best_layout;
    SolverResult best_solution;
    std::vector<double> best_loads;
    bool have_best_full_solution = false;

    std::optional<CommonRoutingWarmStart> common_routing_start;
    double initial_upper_bound = kInf;

    // The common-routing model is genuinely a different formulation, so it is
    // still built and solved once as a separate model.
    if (!opt::interruptRequested() && !budgetExpired()) {
        try {
            common_routing_start = solveCommonRoutingWarmStart(
                instance_,
                splits_,
                options_,
                remainingSeconds());
            initial_upper_bound = common_routing_start->max_load;
        }
        catch (const std::exception& error) {
            std::cerr
                << "Warning: common-routing warm start failed: "
                << error.what() << '\n'
                << "Continuing without initial warm start or load cutoff.\n\n";
        }
    }

    // As long as bottleneck peeling succeeds, this session owns one native
    // solver model for ALL peeled lexicographic iterations and all bottleneck
    // tests. No base T-ASR rows/columns are rebuilt in that phase.
    std::unique_ptr<opt::SolverSession> persistent_session;
    Layout persistent_layout;

    for (int k = 1; k <= requested_iterations; ++k) {
        if (opt::interruptRequested()) {
            std::cout
                << "Ctrl+C received; stopping before lexicographic iteration "
                << k << ".\n";
            break;
        }
        if (budgetExpired()) {
            std::cout
                << "Total time limit reached before lexicographic iteration "
                << k << ".\n";
            break;
        }

        const bool use_simple_max =
            (k == 1) ||
            (bottleneck_peeling_active &&
                static_cast<int>(fixed_bottlenecks.size()) == k - 1);

        const double global_max_cutoff =
            (k == 1)
            ? initial_upper_bound
            : previous_objectives.at(0);

        BackendOptions backend_options;
        backend_options.time_limit_seconds = remainingSeconds();
        backend_options.mip_relative_gap = options_.mip_relative_gap;
        backend_options.solver_output = options_.solver_output;

        SolverResult solution;
        Layout layout;
        std::size_t displayed_columns = 0;
        std::size_t displayed_rows = 0;

        if (use_simple_max) {
            if (k == 1) {
                ModelBuilder builder;

                addBaseVariables(
                    builder,
                    instance_,
                    splits_,
                    persistent_layout,
                    initial_upper_bound,
                    initial_upper_bound,
                    {});
                addBaseConstraints(
                    builder,
                    instance_,
                    splits_,
                    persistent_layout,
                    {});
                addFirstObjective(
                    builder,
                    instance_,
                    persistent_layout,
                    initial_upper_bound);

                const LinearModel model = builder.build();
                persistent_session =
                    opt::createSolverSession(options_.solver, model);

                displayed_columns = persistent_session->numColumns();
                displayed_rows = persistent_session->numRows();
            }
            else {
                if (!persistent_session) {
                    throw std::logic_error(
                        "bottleneck peeling requires an active persistent solver session");
                }
                displayed_columns = persistent_session->numColumns();
                displayed_rows = persistent_session->numRows();
            }

            layout = persistent_layout;

            MipStart first_start;
            const MipStart* first_start_ptr = nullptr;

            // Only the first full T-ASR solve needs an explicit MIP start.
            // Afterwards the same native solver model is re-optimized, so the
            // previous native solution state is retained by the backend.
            if (k == 1 && common_routing_start.has_value()) {
                first_start = makeCommonRoutingMipStart(
                    layout,
                    common_routing_start->x);
                first_start_ptr = &first_start;
            }

            std::string model_file;
            if (options_.write_models) {
                std::ostringstream name;
                name << options_.model_prefix << "-"
                    << opt::solverTypeName(options_.solver)
                    << "-lex-" << k
                    << (k >= 2 ? "-peeled-persistent" : "-persistent")
                    << ".lp";
                model_file = name.str();
            }

            std::cout
                << "Lexicographic iteration " << k << '/'
                << requested_iterations << " with "
                << opt::solverTypeName(options_.solver)
                << (k >= 2
                    ? " [persistent simple max over unfixed loads]"
                    : " [persistent first max]")
                << ": " << displayed_columns << " columns, "
                << displayed_rows << " generic rows"
                << ", remaining total time = " << std::fixed
                << std::setprecision(2)
                << (options_.time_limit_seconds > 0.0
                    ? backend_options.time_limit_seconds
                    : 0.0)
                << (options_.time_limit_seconds > 0.0
                    ? " s\n"
                    : " (unlimited)\n");

            solution = persistent_session->optimize(
                backend_options,
                first_start_ptr,
                model_file);
        }
        else {
            // Individual peeling has stopped. Work only on the residual vector
            // of still-unfixed loads. If m loads were fixed before the switch,
            // global iteration k corresponds to residual order q = k-m.
            if (!residual_ksum_active || residual_fixed_prefix < 0) {
                throw std::logic_error(
                    "general residual k-sum requested without initialized "
                    "residual fallback state");
            }

            const int residual_k = k - residual_fixed_prefix;
            if (residual_k < 2) {
                throw std::logic_error(
                    "residual k-sum fallback requires residual order >= 2");
            }

            std::vector<double> previous_residual_objectives;
            previous_residual_objectives.reserve(
                static_cast<std::size_t>(residual_k - 1));

            // Residual cumulative objective R_j equals the global cumulative
            // objective S_{m+j} minus the already-fixed prefix S_m.
            for (int residual_iteration = 1;
                residual_iteration < residual_k;
                ++residual_iteration) {

                const int global_level =
                    residual_fixed_prefix + residual_iteration;
                if (global_level <= 0 ||
                    global_level >
                    static_cast<int>(previous_objectives.size())) {
                    throw std::logic_error(
                        "missing global objective for residual fallback");
                }

                previous_residual_objectives.push_back(
                    previous_objectives[
                        static_cast<std::size_t>(global_level - 1)] -
                    residual_prefix_objective);
            }

            ModelBuilder builder;

            // Segment pruning still uses the GLOBAL maximum L_1. In contrast,
            // every unfixed lambda can now use the much stronger bound R_1,
            // the last simple-max optimum before peeling failed.
            addBaseVariables(
                builder,
                instance_,
                splits_,
                layout,
                global_max_cutoff,
                residual_max_boundary,
                fixed_bottlenecks);
            addBaseConstraints(
                builder,
                instance_,
                splits_,
                layout,
                fixed_bottlenecks);
            addResidualKSumObjective(
                builder,
                instance_,
                layout,
                residual_k,
                previous_residual_objectives,
                fixed_bottlenecks);

            if (opt::interruptRequested() || budgetExpired()) {
                std::cout
                    << "Stopping after model construction because the total "
                    << "time limit or Ctrl+C was reached.\n";
                break;
            }

            MipStart warm_start;
            const MipStart* warm_start_ptr = nullptr;

            if (have_previous_solution) {
                const int residual_k =
                    k - residual_fixed_prefix;

                warm_start = makeResidualWarmStart(
                    layout,
                    previous_layout,
                    previous_solution,
                    previous_loads,
                    residual_k,
                    fixed_bottlenecks);
                warm_start_ptr = &warm_start;
            }

            std::string model_file;
            if (options_.write_models) {
                std::ostringstream name;
                name << options_.model_prefix << "-"
                    << opt::solverTypeName(options_.solver)
                    << "-lex-" << k
                    << "-residual-ksum.lp";
                model_file = name.str();
            }

            displayed_columns = builder.numColumns();
            displayed_rows = builder.numRows();

            const int residual_k_for_log =
                k - residual_fixed_prefix;
            const int residual_free_loads =
                static_cast<int>(layout.load_position.size()) -
                residual_fixed_prefix;

            std::cout
                << "  Residual k-sum setup: fixed prefix m="
                << residual_fixed_prefix
                << ", residual order q="
                << residual_k_for_log
                << ", residual max bound R_1="
                << std::setprecision(12)
                << residual_max_boundary
                << ", unfixed loads="
                << residual_free_loads
                << "\n";

            std::cout
                << "Lexicographic iteration " << k << '/'
                << requested_iterations << " with "
                << opt::solverTypeName(options_.solver)
                << " [residual k-sum fallback, residual order "
                << (k - residual_fixed_prefix) << "]: "
                << displayed_columns << " columns, "
                << displayed_rows << " rows"
                << ", remaining total time = " << std::fixed
                << std::setprecision(2)
                << (options_.time_limit_seconds > 0.0
                    ? backend_options.time_limit_seconds
                    : 0.0)
                << (options_.time_limit_seconds > 0.0
                    ? " s\n"
                    : " (unlimited)\n");

            const LinearModel model = builder.build();
            solution = opt::solveModel(
                options_.solver,
                model,
                backend_options,
                warm_start_ptr,
                model_file);
        }

        if (!solution.feasible) {
            const bool graceful_stop =
                solution.time_limit_reached ||
                solution.interrupted ||
                opt::interruptRequested() ||
                budgetExpired();

            if (graceful_stop) {
                std::cout
                    << "  solve stopped with status "
                    << solution.status
                    << " before finding a new incumbent; retaining the "
                    << "best earlier feasible solution.\n";
                break;
            }

            throw std::runtime_error(
                std::string(opt::solverTypeName(options_.solver)) +
                " did not return a usable feasible solution; status = " +
                solution.status);
        }

        const std::vector<double> candidate_loads =
            loadValues(layout, solution);

        // Cross-check the native solver result against loads recomputed directly
        // from the binary routing decisions. For peeled min-max iterations this
        // also detects stale objective coefficients or an accidentally disabled
        // max row.
        verifyPhysicalLoadConsistency(
            instance_,
            splits_,
            layout,
            solution,
            fixed_bottlenecks,
            use_simple_max,
            solution.objective,
            k,
            "main lexicographic solve");

        if (!have_best_full_solution ||
            lexicographicallyBetterLoads(candidate_loads, best_loads)) {
            best_layout = layout;
            best_solution = solution;
            best_loads = candidate_loads;
            have_best_full_solution = true;
        }

        if (!solution.optimal) {
            std::cout
                << "  solve stopped with status "
                << solution.status
                << "; retaining its best feasible incumbent.\n";
            break;
        }

        double current_component = 0.0;
        double cumulative_objective = 0.0;

        if (use_simple_max) {
            current_component = solution.objective;
            cumulative_objective =
                (previous_objectives.empty()
                    ? 0.0
                    : previous_objectives.back()) +
                current_component;
        }
        else {
            // The fallback objective is the residual ordered sum R_q over only
            // the still-unfixed loads. Convert it back to the global cumulative
            // objective S_k by adding the already-fixed prefix S_m.
            cumulative_objective =
                residual_prefix_objective + solution.objective;

            const double previous_sum =
                previous_objectives.empty()
                ? 0.0
                : previous_objectives.back();

            current_component =
                cumulative_objective - previous_sum;
        }

        // Ordered loads must be non-increasing:
        //     L_1 >= L_2 >= ... >= L_k.
        // A material increase means the current model is no longer representing
        // the intended lexicographic sequence (for example because previously
        // preserved levels and solver feasibility tolerances became inconsistent).
        if (!lexicographic_load_values.empty()) {
            const double previous_component = lexicographic_load_values.back();
            if (current_component >
                previous_component + kLexTolerance) {

                std::ostringstream message;
                message << std::setprecision(16)
                    << "LEXICOGRAPHIC INVARIANT ERROR: ordered load increased at "
                    << "iteration " << k
                    << ": L_" << (k - 1) << "=" << previous_component
                    << ", L_" << k << "=" << current_component
                    << ", increase=" << (current_component - previous_component)
                    << " > lexicographic tolerance " << kLexTolerance;
                throw std::runtime_error(message.str());
            }
        }

        // The sum of the k largest loads cannot decrease as k grows.
        if (!previous_objectives.empty()) {
            const double previous_cumulative = previous_objectives.back();
            if (cumulative_objective + kLexTolerance <
                previous_cumulative) {

                std::ostringstream message;
                message << std::setprecision(16)
                    << "LEXICOGRAPHIC INVARIANT ERROR: cumulative ordered-load "
                    << "objective decreased at iteration " << k
                    << ": previous=" << previous_cumulative
                    << ", current=" << cumulative_objective
                    << ", decrease="
                    << (previous_cumulative - cumulative_objective)
                    << " > lexicographic tolerance " << kLexTolerance;
                throw std::runtime_error(message.str());
            }
        }

        if (current_component < -kLexTolerance) {
            std::ostringstream message;
            message << std::setprecision(16)
                << "LEXICOGRAPHIC INVARIANT ERROR: negative ordered load at "
                << "iteration " << k << ": " << current_component;
            throw std::runtime_error(message.str());
        }

        previous_objectives.push_back(cumulative_objective);
        lexicographic_load_values.push_back(current_component);
        previous_loads = candidate_loads;
        previous_layout = layout;
        previous_solution = solution;
        have_previous_solution = true;
        completed_optimal_iterations = k;

        std::cout
            << "  optimal ordered load L_" << k << " = "
            << std::setprecision(16) << current_component << '\n'
            << "  cumulative objective (sum of " << k
            << " largest loads) = "
            << cumulative_objective << '\n';

        // After z1 is known, exploit the stronger global cutoff without
        // rebuilding the persistent model: variables that would previously
        // have been omitted are simply fixed to zero. Column indices stay stable.
        if (k == 1 && persistent_session) {
            prunePersistentSegments(
                *persistent_session,
                instance_,
                splits_,
                persistent_layout,
                current_component);
        }

        // -------------------------------------------------------------
        // Persistent bottleneck unavoidability test
        // -------------------------------------------------------------
        if (use_simple_max &&
            bottleneck_peeling_active &&
            persistent_session &&
            !opt::interruptRequested() &&
            !budgetExpired()) {

            const std::vector<double> actual_loads =
                actualLoadValues(
                    instance_,
                    splits_,
                    layout,
                    solution);

            // In a simple-max iteration the optimized component must equal the
            // largest physical load among the still-unfixed arc/time pairs.
            double largest_unfixed_physical = -kInf;
            for (int p = 0;
                p < static_cast<int>(layout.load_position.size());
                ++p) {
                const auto [t, a] = layout.load_position[p];
                if (loadIsFixed(fixed_bottlenecks, t, a)) continue;
                largest_unfixed_physical =
                    std::max(largest_unfixed_physical, actual_loads[p]);
            }
            if (std::isfinite(largest_unfixed_physical)) {
                reportInvariantDifference(
                    "pre-proof ordered component vs largest physical unfixed load",
                    current_component,
                    largest_unfixed_physical);
            }

            int candidate_position = -1;
            double candidate_actual_load = -kInf;

            for (int p = 0;
                p < static_cast<int>(layout.load_position.size());
                ++p) {

                const auto [t, a] = layout.load_position[p];
                if (loadIsFixed(fixed_bottlenecks, t, a)) {
                    continue;
                }

                if (actual_loads[p] > candidate_actual_load) {
                    candidate_actual_load = actual_loads[p];
                    candidate_position = p;
                }
            }

            if (candidate_position < 0) {
                std::cout
                    << "  No unfixed bottleneck candidate remains; "
                    << "stopping bottleneck peeling.\n";
                bottleneck_peeling_active = false;
                persistent_session.reset();
            }
            else {
                const auto [candidate_time, candidate_arc] =
                    layout.load_position[candidate_position];

                if (candidate_actual_load <
                    current_component - 10.0 * kTolerance) {

                    std::cout
                        << "  Warning: largest physical unfixed load ("
                        << std::setprecision(12)
                        << candidate_actual_load
                        << ") is below the optimized boundary ("
                        << current_component
                        << "); not attempting a fixing proof.\n";

                    bottleneck_peeling_active = false;
                    persistent_session.reset();
                }
                else {
                    const BottleneckTestResult test =
                        testBottleneckCandidatePersistent(
                            *persistent_session,
                            instance_,
                            splits_,
                            options_,
                            remainingSeconds(),
                            current_component,
                            fixed_bottlenecks,
                            candidate_time,
                            candidate_arc,
                            persistent_layout,
                            k);

                    if (test.proved_fixed) {
                        fixPersistentBottleneck(
                            *persistent_session,
                            persistent_layout,
                            candidate_time,
                            candidate_arc,
                            current_component);

                        fixed_bottlenecks.push_back({
                            candidate_time,
                            candidate_arc,
                            current_component });

                        std::cout
                            << "    FIXED physical load lambda(t="
                            << candidate_time
                            << ", a="
                            << arcNameId(instance_, candidate_arc)
                            << ") ~= "
                            << std::setprecision(12)
                            << current_component
                            << ". The same native solver model will be "
                            << "re-optimized for the next level.\n";
                    }
                    else {
                        // We have fixed exactly the prefix L_1,...,L_m, while
                        // the current simple-max iteration has just optimized
                        // R_1 = L_{m+1} over all remaining loads. Preserve this
                        // residual maximum directly through lambda upper bounds
                        // in every subsequent residual k-sum model.
                        residual_fixed_prefix =
                            static_cast<int>(fixed_bottlenecks.size());

                        residual_prefix_objective =
                            residual_fixed_prefix > 0
                            ? previous_objectives.at(
                                static_cast<std::size_t>(
                                    residual_fixed_prefix - 1))
                            : 0.0;

                        residual_max_boundary = current_component;
                        residual_ksum_active = true;

                        bottleneck_peeling_active = false;
                        persistent_session.reset();

                        std::cout
                            << "    Bottleneck cannot be fixed individually. "
                            << "Persistent peeling stops. The next iteration "
                            << "uses a residual k-sum formulation over only the "
                            << "unfixed loads.\n"
                            << "    Fixed prefix length m = "
                            << residual_fixed_prefix
                            << ", residual R_1 = "
                            << std::setprecision(12)
                            << residual_max_boundary
                            << ".\n";
                    }
                }
            }
        }

        if (opt::interruptRequested() || budgetExpired()) {
            std::cout
                << (opt::interruptRequested()
                    ? "Ctrl+C received; "
                    : "Total time limit reached; ")
                << "returning the best feasible solution found so far.\n";
            break;
        }
    }

    TasrSolution answer;
    answer.lexicographic_objectives = previous_objectives;
    answer.completed_lexicographic_iterations =
        completed_optimal_iterations;

    if (have_best_full_solution) {
        answer.sorted_loads = best_loads;
        std::sort(
            answer.sorted_loads.begin(),
            answer.sorted_loads.end(),
            std::greater<double>());
        answer.waypoints =
            extractWaypoints(
                instance_,
                best_layout,
                best_solution);
    }
    else if (common_routing_start.has_value()) {
        answer.sorted_loads = common_routing_start->loads;
        std::sort(
            answer.sorted_loads.begin(),
            answer.sorted_loads.end(),
            std::greater<double>());
        answer.waypoints =
            extractCommonRoutingWaypoints(
                instance_,
                common_routing_start->x);

        std::cout
            << "No full lexicographic incumbent was available; writing the "
            << "feasible common-routing solution.\n";
    }
    else {
        throw std::runtime_error(
            "stopped before any feasible incumbent was available; "
            "no solution can be written");
    }

    std::cout
        << "Total solver-process time: "
        << std::fixed
        << std::setprecision(2)
        << elapsedSeconds()
        << " seconds\n";

    return answer;
}

void writeSolutionJson(const TasrSolution& solution,
    const std::string& filename) {
    std::ofstream output(filename);
    if (!output) throw std::runtime_error("cannot open solution file: " + filename);

    output << "{\n  \"srpaths\": [";
    bool first = true;
    for (int d = 0; d < static_cast<int>(solution.waypoints.size()); ++d) {
        for (int t = 0; t < static_cast<int>(solution.waypoints[d].size()); ++t) {
            const auto& waypoints = solution.waypoints[d][t];
            if (waypoints.empty()) continue;
            if (!first) output << ',';
            output << "\n    { \"d\": " << d << ", \"t\": " << t << ", \"w\": [";
            for (std::size_t p = 0; p < waypoints.size(); ++p) {
                if (p) output << ", ";
                output << waypoints[p];
            }
            output << "] }";
            first = false;
        }
    }
    if (!first) output << '\n';
    output << "  ]\n}\n";
    if (!output) throw std::runtime_error("failed while writing solution file: " + filename);
}