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
        double max_load_cutoff = kInf) {
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

        for (int d = 0; d < layout.demands; ++d) {
            for (int t = 0; t < layout.times; ++t) {
                for (int i = 0; i < layout.nodes; ++i) {
                    for (int j = 0; j < layout.nodes; ++j) {
                        if (!segmentIsReachable(splits, t, i, j)) continue;

                        // If x[d,t,i,j] = 1 would by itself force an arc load
                        // strictly above the best known first-objective bound,
                        // this variable can safely be fixed to zero.
                        if (std::isfinite(max_load_cutoff) &&
                            segmentLoadLowerBound(
                                instance, splits, d, t, i, j) >
                            max_load_cutoff + kFixingTolerance) {
                            continue;
                        }

                        layout.x[layout.segmentOffset(d, t, i, j)] =
                            model.addColumn(
                                "x_d" + std::to_string(d) +
                                "_t" + std::to_string(t) +
                                "_i" + std::to_string(i) +
                                "_j" + std::to_string(j),
                                0.0, 0.0, 1.0, VariableType::Integer);
                    }
                }
            }
        }

        for (int t = 0; t < layout.times; ++t) {
            for (int a = 0; a < layout.arcs; ++a) {
                if (!arcIsAvailable(instance, t, a)) continue;
                layout.lambda[layout.loadOffset(t, a)] =
                    model.addColumn(
                        "lambda_a" + std::to_string(arcNameId(instance, a)) +
                        "_t" + std::to_string(t),
                        0.0, 0.0, kInf);
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
        const Layout& layout) {
        // Arc load constraints: capacity * lambda - routed traffic >= 0.
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
            model.addRow(
                "load_a" + std::to_string(arcNameId(instance, a)) +
                "_t" + std::to_string(t),
                0.0, kInf, std::move(index), std::move(value));
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

    void addKSumObjective(ModelBuilder& model,
        const InstanceData& instance,
        Layout& layout,
        int k,
        const std::vector<double>& previous_objectives) {

        if (k < 2) {
            throw std::logic_error("addKSumObjective requires k >= 2");
        }
        if (previous_objectives.size() <
            static_cast<std::size_t>(k - 1)) {
            throw std::logic_error(
                "missing previous lexicographic objective values");
        }

        const int number_of_loads =
            static_cast<int>(layout.load_position.size());

        // Index directly by lexicographic iteration. Entries 0 and 1 are
        // intentionally unused because iteration 1 has no alpha/beta model.
        layout.alpha.assign(k + 1, -1);
        layout.beta.assign(
            k + 1,
            std::vector<int>(number_of_loads, -1));

        // Preserve the optimal first objective directly:
        //     lambda_{a,t} <= z_1
        // No alpha_1 or beta_1 variables are introduced.
        const double z1 = previous_objectives.at(0);
        for (const auto [t, a] : layout.load_position) {
            model.addRow(
                "preserve_maxload_a" +
                std::to_string(arcNameId(instance, a)) +
                "_t" + std::to_string(t),
                -kInf,
                z1 + kFixingTolerance,
                { layout.lambda[layout.loadOffset(t, a)] },
                { 1.0 });
        }

        // Create auxiliary variables only for lexicographic levels 2..k.
        // Only level k receives objective coefficients.
        for (int iteration = 2; iteration <= k; ++iteration) {
            const double alpha_cost =
                (iteration == k) ? static_cast<double>(k) : 0.0;

            layout.alpha[iteration] = model.addColumn(
                "alpha_k" + std::to_string(iteration),
                alpha_cost,
                -kInf,
                kInf);

            for (int p = 0; p < number_of_loads; ++p) {
                const double beta_cost =
                    (iteration == k) ? 1.0 : 0.0;
                const auto [t, a] = layout.load_position[p];

                layout.beta[iteration][p] = model.addColumn(
                    "beta_k" + std::to_string(iteration) +
                    "_a" + std::to_string(arcNameId(instance, a)) +
                    "_t" + std::to_string(t),
                    beta_cost,
                    0.0,
                    kInf);
            }
        }

        // beta^i_{a,t} >= lambda_{a,t} - alpha_i, i = 2..k.
        for (int iteration = 2; iteration <= k; ++iteration) {
            for (int p = 0; p < number_of_loads; ++p) {
                const auto [t, a] = layout.load_position[p];

                model.addRow(
                    "beta_excess_k" + std::to_string(iteration) +
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

        // Preserve already optimized k-sum values for levels 2..k-1:
        //     i alpha_i + sum beta_i <= z_i.
        // Level 1 is preserved by lambda <= z_1 above.
        for (int iteration = 2; iteration < k; ++iteration) {
            std::vector<int> index;
            std::vector<double> value;

            appendTerm(
                index,
                value,
                layout.alpha[iteration],
                static_cast<double>(iteration));

            for (const int beta : layout.beta[iteration]) {
                appendTerm(index, value, beta, 1.0);
            }

            const double rhs =
                previous_objectives.at(
                    static_cast<std::size_t>(iteration - 1)) +
                kFixingTolerance;

            model.addRow(
                "preserve_ksum_k" + std::to_string(iteration),
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
                            0.0, 0.0, 1.0, VariableType::Integer);
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

        // Same load equations as in the full model, but with common x[d,i,j].
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
                0.0, kInf,
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

    MipStart makeWarmStart(
        const Layout& layout,
        const Layout& previous_layout,
        const SolverResult& previous_solution,
        const std::vector<double>& previous_loads,
        int k) {

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

        // There are no alpha_1/beta_1 variables. Initialize levels 2..k
        // from the previous load vector.
        if (k >= 2 && !previous_loads.empty()) {
            std::vector<double> sorted = previous_loads;
            std::sort(sorted.begin(), sorted.end(), std::greater<double>());

            for (int iteration = 2; iteration <= k; ++iteration) {
                const int sorted_index = std::min<int>(
                    iteration - 1, static_cast<int>(sorted.size()) - 1);
                const double alpha = sorted[sorted_index];
                addStart(layout.alpha[iteration], alpha);

                for (int p = 0; p < static_cast<int>(previous_loads.size()); ++p) {
                    addStart(layout.beta[iteration][p],
                             std::max(0.0, previous_loads[p] - alpha));
                }
            }
        }

        return warm;
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
        return std::chrono::duration<double>(Clock::now() - process_start).count();
    };

    auto remainingSeconds = [&]() {
        if (options_.time_limit_seconds <= 0.0) return 0.0; // backend convention: unlimited
        return std::max(0.0, options_.time_limit_seconds - elapsedSeconds());
    };

    auto budgetExpired = [&]() {
        return options_.time_limit_seconds > 0.0 && remainingSeconds() <= 0.0;
    };

    const int total_loads = [&]() {
        int count = 0;
        for (int t = 0; t < instance_.num_time_slots; ++t)
            for (int a = 0; a < instance_.numArcs(); ++a)
                count += arcIsAvailable(instance_, t, a) ? 1 : 0;
        return count;
    }();
    if (total_loads == 0) {
        throw std::runtime_error("instance has no available arc/time loads");
    }

    int requested_iterations = options_.max_lexicographic_iterations;
    if (requested_iterations < 0) requested_iterations = total_loads;
    requested_iterations = std::max(1, std::min(requested_iterations, total_loads));

    std::vector<double> previous_objectives;
    std::vector<double> previous_loads;
    Layout previous_layout;
    SolverResult previous_solution;
    bool have_previous_solution = false;
    int completed_optimal_iterations = 0;

    // Best feasible full-model solution seen so far. This can be the solution
    // of the last completed lexicographic iteration or the incumbent returned
    // by an interrupted / time-limited current iteration.
    Layout best_layout;
    SolverResult best_solution;
    std::vector<double> best_loads;
    bool have_best_full_solution = false;

    // The common-routing model is only an acceleration device, but its
    // solution is itself feasible for the original problem and is retained as
    // a last-resort fallback if no full lexicographic model returns an incumbent.
    std::optional<CommonRoutingWarmStart> common_routing_start;
    double initial_upper_bound = kInf;

    if (!opt::interruptRequested() && !budgetExpired()) {
        try {
            const double warm_time = remainingSeconds();
            common_routing_start = solveCommonRoutingWarmStart(
                instance_, splits_, options_, warm_time);
            initial_upper_bound = common_routing_start->max_load;
        }
        catch (const std::exception& error) {
            std::cerr
                << "Warning: common-routing warm start failed: "
                << error.what() << '\n'
                << "Continuing without initial warm start or load cutoff.\n\n";
        }
    }

    for (int k = 1; k <= requested_iterations; ++k) {
        if (opt::interruptRequested()) {
            std::cout << "Ctrl+C received; stopping before lexicographic iteration "
                      << k << ".\n";
            break;
        }
        if (budgetExpired()) {
            std::cout << "Total time limit reached before lexicographic iteration "
                      << k << ".\n";
            break;
        }

        ModelBuilder builder;
        Layout layout;

        const double max_load_cutoff =
            (k == 1) ? initial_upper_bound : previous_objectives.at(0);

        addBaseVariables(builder, instance_, splits_, layout, max_load_cutoff);
        addBaseConstraints(builder, instance_, splits_, layout);

        if (k == 1) {
            addFirstObjective(builder, instance_, layout, initial_upper_bound);
        }
        else {
            addKSumObjective(builder, instance_, layout, k, previous_objectives);
        }

        if (opt::interruptRequested() || budgetExpired()) {
            std::cout << "Stopping after model construction because the total "
                         "time limit or Ctrl+C was reached.\n";
            break;
        }

        BackendOptions backend_options;
        backend_options.time_limit_seconds = remainingSeconds();
        backend_options.mip_relative_gap = options_.mip_relative_gap;
        backend_options.solver_output = options_.solver_output;

        MipStart warm_start;
        const MipStart* warm_start_ptr = nullptr;
        if (k == 1 && common_routing_start.has_value()) {
            warm_start = makeCommonRoutingMipStart(layout, common_routing_start->x);
            warm_start_ptr = &warm_start;
        }
        else if (have_previous_solution) {
            warm_start = makeWarmStart(
                layout, previous_layout, previous_solution, previous_loads, k);
            warm_start_ptr = &warm_start;
        }

        std::string model_file;
        if (options_.write_models) {
            std::ostringstream name;
            name << options_.model_prefix << "-"
                 << opt::solverTypeName(options_.solver)
                 << "-lex-" << k << ".lp";
            model_file = name.str();
        }

        std::cout << "Lexicographic iteration " << k << '/' << requested_iterations
                  << " with " << opt::solverTypeName(options_.solver)
                  << ": " << builder.numColumns() << " columns, "
                  << builder.numRows() << " rows"
                  << ", remaining total time = " << std::fixed
                  << std::setprecision(2)
                  << (options_.time_limit_seconds > 0.0
                          ? backend_options.time_limit_seconds
                          : 0.0)
                  << (options_.time_limit_seconds > 0.0 ? " s\n" : " (unlimited)\n");

        const LinearModel model = builder.build();
        const SolverResult solution = opt::solveModel(
            options_.solver, model, backend_options, warm_start_ptr, model_file);

        if (solution.feasible) {
            const std::vector<double> candidate_loads = loadValues(layout, solution);

            if (!have_best_full_solution ||
                lexicographicallyBetterLoads(candidate_loads, best_loads)) {
                best_layout = layout;
                best_solution = solution;
                best_loads = candidate_loads;
                have_best_full_solution = true;
            }

            // Only an optimal solve may define the fixed objective value used
            // by the next lexicographic iteration. A time-limited/interrupted
            // incumbent is valid output, but must not seed another lex level.
            if (solution.optimal) {
                previous_objectives.push_back(solution.objective);
                previous_loads = candidate_loads;
                previous_layout = layout;
                previous_solution = solution;
                have_previous_solution = true;
                completed_optimal_iterations = k;

                std::cout << "  optimal objective (sum of " << k
                          << " largest loads): " << std::setprecision(12)
                          << solution.objective << '\n';
            }
            else {
                std::cout << "  solve stopped with status " << solution.status
                          << "; retaining its best feasible incumbent.\n";
                break;
            }
        }
        else {
            const bool graceful_stop =
                solution.time_limit_reached || solution.interrupted ||
                opt::interruptRequested() || budgetExpired();

            if (graceful_stop) {
                std::cout << "  solve stopped with status " << solution.status
                          << " before finding a new incumbent; retaining the "
                             "best earlier feasible solution.\n";
                break;
            }

            throw std::runtime_error(
                std::string(opt::solverTypeName(options_.solver)) +
                " did not return a usable feasible solution; status = " +
                solution.status);
        }

        if (opt::interruptRequested() || budgetExpired()) {
            std::cout << (opt::interruptRequested()
                              ? "Ctrl+C received; "
                              : "Total time limit reached; ")
                      << "returning the best feasible solution found so far.\n";
            break;
        }
    }

    TasrSolution answer;
    answer.lexicographic_objectives = previous_objectives;
    answer.completed_lexicographic_iterations = completed_optimal_iterations;

    if (have_best_full_solution) {
        answer.sorted_loads = best_loads;
        std::sort(
            answer.sorted_loads.begin(),
            answer.sorted_loads.end(),
            std::greater<double>());
        answer.waypoints = extractWaypoints(instance_, best_layout, best_solution);
    }
    else if (common_routing_start.has_value()) {
        answer.sorted_loads = common_routing_start->loads;
        std::sort(
            answer.sorted_loads.begin(),
            answer.sorted_loads.end(),
            std::greater<double>());
        answer.waypoints = extractCommonRoutingWaypoints(
            instance_, common_routing_start->x);
        std::cout << "No full lexicographic incumbent was available; writing the "
                     "feasible common-routing solution.\n";
    }
    else {
        throw std::runtime_error(
            "stopped before any feasible incumbent was available; no solution "
            "can be written");
    }

    std::cout << "Total solver-process time: " << std::fixed
              << std::setprecision(2) << elapsedSeconds() << " seconds\n";
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