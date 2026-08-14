#include "solver_backend.h"

#include <algorithm>
#include <cctype>
#include <csignal>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef ROADEF_HAVE_HIGHS
#include "Highs.h"
#endif

#ifdef ROADEF_HAVE_GUROBI
#include "gurobi_c.h"
#endif

namespace opt {

    namespace {
        volatile std::sig_atomic_t g_interrupt_requested = 0;
    }

    void requestInterrupt() noexcept { g_interrupt_requested = 1; }
    void clearInterrupt() noexcept { g_interrupt_requested = 0; }
    bool interruptRequested() noexcept { return g_interrupt_requested != 0; }

    namespace {

        void validateModel(const LinearModel& model) {
            const std::size_t n = model.cost.size();
            if (model.column_names.size() != n ||
                model.lower.size() != n ||
                model.upper.size() != n ||
                model.variable_type.size() != n) {
                throw std::logic_error("inconsistent model column vectors");
            }

            for (const auto& row : model.rows) {
                if (row.index.size() != row.value.size()) {
                    throw std::logic_error("row index/value size mismatch");
                }
                for (const int column : row.index) {
                    if (column < 0 || static_cast<std::size_t>(column) >= n) {
                        throw std::logic_error("constraint contains invalid column index");
                    }
                }
            }
        }

        void validateMipStart(const MipStart& start, std::size_t num_columns) {
            if (start.index.size() != start.value.size()) {
                throw std::logic_error("MIP-start index/value size mismatch");
            }
            for (const int column : start.index) {
                if (column < 0 || static_cast<std::size_t>(column) >= num_columns) {
                    throw std::logic_error("MIP start contains invalid column index");
                }
            }
        }

#ifdef ROADEF_HAVE_HIGHS
        void checkHighsStatus(HighsStatus status, const char* operation) {
            if (status == HighsStatus::kError) {
                throw std::runtime_error(std::string("HiGHS failed while ") + operation);
            }
        }

        bool highsHasUsableSolution(const Highs& highs) {
            switch (highs.getModelStatus()) {
            case HighsModelStatus::kOptimal:
            case HighsModelStatus::kTimeLimit:
            case HighsModelStatus::kIterationLimit:
            case HighsModelStatus::kObjectiveBound:
            case HighsModelStatus::kObjectiveTarget:
            case HighsModelStatus::kInterrupt:
                return true;
            default:
                return false;
            }
        }

        void highsInterruptCallback(const int /*callback_type*/,
            const char* /*message*/,
            const HighsCallbackDataOut* /*data_out*/,
            HighsCallbackDataIn* data_in,
            void* /*user_callback_data*/) {
            if (interruptRequested()) {
                data_in->user_interrupt = 1;
            }
        }

        SolverResult solveHighs(const LinearModel& model,
            const BackendOptions& options,
            const MipStart* mip_start,
            const std::string& model_file) {
            HighsLp lp;
            lp.num_col_ = static_cast<HighsInt>(model.cost.size());
            lp.num_row_ = static_cast<HighsInt>(model.rows.size());
            lp.col_cost_ = model.cost;
            lp.col_lower_ = model.lower;
            lp.col_upper_ = model.upper;
            lp.col_names_ = model.column_names;
            lp.integrality_.reserve(model.variable_type.size());

            for (const VariableType type : model.variable_type) {
                lp.integrality_.push_back(
                    type == VariableType::Integer
                    ? HighsVarType::kInteger
                    : HighsVarType::kContinuous);
            }

            lp.row_names_.reserve(model.rows.size());
            lp.row_lower_.reserve(model.rows.size());
            lp.row_upper_.reserve(model.rows.size());

            // HiGHS wants the sparse matrix in compressed-column format here.
            std::vector<std::vector<std::pair<HighsInt, double>>> columns(model.cost.size());
            for (HighsInt r = 0; r < lp.num_row_; ++r) {
                const SparseRow& row = model.rows.at(static_cast<std::size_t>(r));
                lp.row_names_.push_back(row.name);
                lp.row_lower_.push_back(row.lower);
                lp.row_upper_.push_back(row.upper);

                for (std::size_t p = 0; p < row.index.size(); ++p) {
                    columns.at(static_cast<std::size_t>(row.index[p]))
                        .emplace_back(r, row.value[p]);
                }
            }

            lp.a_matrix_.format_ = MatrixFormat::kColwise;
            lp.a_matrix_.num_col_ = lp.num_col_;
            lp.a_matrix_.num_row_ = lp.num_row_;

            // HighsSparseMatrix starts with internal storage that must not be
            // assumed to be empty. In particular, leaving a pre-existing 0 in
            // start_ and then pushing another 0 shifts every matrix column by one:
            // coefficients intended for column c are written under column c+1.
            // Clear all sparse-matrix arrays explicitly before constructing CSC.
            lp.a_matrix_.start_.clear();
            lp.a_matrix_.index_.clear();
            lp.a_matrix_.value_.clear();

            std::size_t number_of_nonzeros = 0;
            for (const auto& column : columns) {
                number_of_nonzeros += column.size();
            }

            lp.a_matrix_.start_.reserve(model.cost.size() + 1);
            lp.a_matrix_.index_.reserve(number_of_nonzeros);
            lp.a_matrix_.value_.reserve(number_of_nonzeros);
            lp.a_matrix_.start_.push_back(0);

            for (const auto& column : columns) {
                for (const auto& [row, coefficient] : column) {
                    lp.a_matrix_.index_.push_back(row);
                    lp.a_matrix_.value_.push_back(coefficient);
                }
                lp.a_matrix_.start_.push_back(
                    static_cast<HighsInt>(lp.a_matrix_.index_.size()));
            }

            Highs highs;
            checkHighsStatus(
                highs.setCallback(highsInterruptCallback, nullptr),
                "setting interrupt callback");
            checkHighsStatus(
                highs.startCallback(
                    static_cast<int>(HighsCallbackType::kCallbackMipInterrupt)),
                "starting MIP interrupt callback");

            checkHighsStatus(
                highs.setOptionValue("output_flag", options.solver_output),
                "setting output_flag");
            checkHighsStatus(
                highs.setOptionValue("mip_rel_gap", options.mip_relative_gap),
                "setting mip_rel_gap");

            if (options.time_limit_seconds > 0.0) {
                checkHighsStatus(
                    highs.setOptionValue("time_limit", options.time_limit_seconds),
                    "setting time_limit");
            }

            checkHighsStatus(highs.passModel(std::move(lp)), "passing model");

            if (mip_start && !mip_start->empty()) {
                std::vector<HighsInt> indices;
                indices.reserve(mip_start->index.size());
                for (const int index : mip_start->index) {
                    indices.push_back(static_cast<HighsInt>(index));
                }

                checkHighsStatus(
                    highs.setSolution(
                        static_cast<HighsInt>(indices.size()),
                        indices.data(),
                        mip_start->value.data()),
                    "setting MIP start");
            }

            if (!model_file.empty()) {
                checkHighsStatus(highs.writeModel(model_file), "writing model");
            }

            checkHighsStatus(highs.run(), "solving model");

            SolverResult result;
            result.status = highs.modelStatusToString(highs.getModelStatus());
            result.optimal = highs.getModelStatus() == HighsModelStatus::kOptimal;
            result.time_limit_reached =
                highs.getModelStatus() == HighsModelStatus::kTimeLimit;
            result.interrupted =
                highs.getModelStatus() == HighsModelStatus::kInterrupt ||
                interruptRequested();
            result.feasible =
                highsHasUsableSolution(highs) && highs.getSolution().value_valid;

            if (result.feasible) {
                result.objective = highs.getObjectiveValue();
                result.values = highs.getSolution().col_value;
            }

            return result;
        }
#endif

#ifdef ROADEF_HAVE_GUROBI

        struct GurobiEnvDeleter {
            void operator()(GRBenv* env) const noexcept {
                if (env) GRBfreeenv(env);
            }
        };

        struct GurobiModelDeleter {
            void operator()(GRBmodel* model) const noexcept {
                if (model) GRBfreemodel(model);
            }
        };

        using GurobiEnvPtr = std::unique_ptr<GRBenv, GurobiEnvDeleter>;
        using GurobiModelPtr = std::unique_ptr<GRBmodel, GurobiModelDeleter>;

        [[noreturn]] void throwGurobiError(
            GRBenv* env,
            int error,
            const std::string& operation) {

            std::string message = "Gurobi error " + std::to_string(error) +
                " while " + operation;
            if (env) {
                const char* detail = GRBgeterrormsg(env);
                if (detail && *detail) {
                    message += ": ";
                    message += detail;
                }
            }
            throw std::runtime_error(message);
        }

        void checkGurobi(
            GRBenv* env,
            int error,
            const std::string& operation) {

            if (error != 0) {
                throwGurobiError(env, error, operation);
            }
        }

        double gurobiBound(double value) {
            if (std::isinf(value)) {
                return value > 0.0 ? GRB_INFINITY : -GRB_INFINITY;
            }
            return value;
        }

        std::string gurobiStatusString(int status) {
            switch (status) {
            case GRB_LOADED: return "Loaded";
            case GRB_OPTIMAL: return "Optimal";
            case GRB_INFEASIBLE: return "Infeasible";
            case GRB_INF_OR_UNBD: return "InfOrUnbd";
            case GRB_UNBOUNDED: return "Unbounded";
            case GRB_CUTOFF: return "Cutoff";
            case GRB_ITERATION_LIMIT: return "IterationLimit";
            case GRB_NODE_LIMIT: return "NodeLimit";
            case GRB_TIME_LIMIT: return "TimeLimit";
            case GRB_SOLUTION_LIMIT: return "SolutionLimit";
            case GRB_INTERRUPTED: return "Interrupted";
            case GRB_NUMERIC: return "Numeric";
            case GRB_SUBOPTIMAL: return "Suboptimal";
#ifdef GRB_USER_OBJ_LIMIT
            case GRB_USER_OBJ_LIMIT: return "UserObjLimit";
#endif
#ifdef GRB_WORK_LIMIT
            case GRB_WORK_LIMIT: return "WorkLimit";
#endif
#ifdef GRB_MEM_LIMIT
            case GRB_MEM_LIMIT: return "MemLimit";
#endif
            default: return "Status(" + std::to_string(status) + ")";
            }
        }

        void addGurobiConstraint(
            GRBenv* env,
            GRBmodel* grb_model,
            const SparseRow& row) {

            // GRBaddconstr takes non-const pointers although it does not conceptually
            // modify the index/value arrays, so use temporary vectors here.
            std::vector<int> indices = row.index;
            std::vector<double> values = row.value;

            const bool finite_lower = std::isfinite(row.lower);
            const bool finite_upper = std::isfinite(row.upper);

            if (!finite_lower && !finite_upper) {
                return;
            }

            if (finite_lower && finite_upper &&
                std::abs(row.lower - row.upper) <= 1e-12) {

                checkGurobi(
                    env,
                    GRBaddconstr(
                        grb_model,
                        static_cast<int>(indices.size()),
                        indices.empty() ? nullptr : indices.data(),
                        values.empty() ? nullptr : values.data(),
                        GRB_EQUAL,
                        row.lower,
                        row.name.c_str()),
                    "adding equality constraint '" + row.name + "'");
                return;
            }

            if (finite_lower) {
                const std::string name = row.name + "_lb";
                checkGurobi(
                    env,
                    GRBaddconstr(
                        grb_model,
                        static_cast<int>(indices.size()),
                        indices.empty() ? nullptr : indices.data(),
                        values.empty() ? nullptr : values.data(),
                        GRB_GREATER_EQUAL,
                        row.lower,
                        name.c_str()),
                    "adding lower-bound constraint '" + name + "'");
            }

            if (finite_upper) {
                const std::string name = row.name + "_ub";
                checkGurobi(
                    env,
                    GRBaddconstr(
                        grb_model,
                        static_cast<int>(indices.size()),
                        indices.empty() ? nullptr : indices.data(),
                        values.empty() ? nullptr : values.data(),
                        GRB_LESS_EQUAL,
                        row.upper,
                        name.c_str()),
                    "adding upper-bound constraint '" + name + "'");
            }
        }

        GRBenv* sharedGurobiEnvironment() {
            // The solver is sequential, so one process-wide Gurobi environment can
            // safely be reused by the common-routing model and all lexicographic
            // iterations. Gurobi also recommends reusing an environment when possible.
            static GurobiEnvPtr environment = [] {
                GRBenv* raw_env = nullptr;
                const int error = GRBloadenv(&raw_env, nullptr);
                if (error != 0) {
                    throwGurobiError(
                        raw_env, error, "creating the Gurobi environment");
                }
                return GurobiEnvPtr(raw_env);
                }();

            return environment.get();
        }

        int __stdcall gurobiInterruptCallback(GRBmodel* model,
            void* /*cbdata*/,
            int /*where*/,
            void* /*usrdata*/) {
            if (interruptRequested()) {
                GRBterminate(model);
            }
            return 0;
        }

        SolverResult solveGurobi(const LinearModel& model,
            const BackendOptions& options,
            const MipStart* mip_start,
            const std::string& model_file) {
            GRBenv* env = sharedGurobiEnvironment();

            checkGurobi(
                env,
                GRBsetintparam(
                    env,
                    "OutputFlag",
                    options.solver_output ? 1 : 0),
                "setting OutputFlag");

            checkGurobi(
                env,
                GRBsetdblparam(
                    env,
                    "MIPGap",
                    options.mip_relative_gap),
                "setting MIPGap");

            if (options.time_limit_seconds > 0.0) {
                checkGurobi(
                    env,
                    GRBsetdblparam(
                        env,
                        "TimeLimit",
                        options.time_limit_seconds),
                    "setting TimeLimit");
            }

            const int num_vars = static_cast<int>(model.cost.size());

            std::vector<double> lb(num_vars);
            std::vector<double> ub(num_vars);
            std::vector<char> vtype(num_vars);
            std::vector<char*> names(num_vars);

            for (int c = 0; c < num_vars; ++c) {
                lb[c] = gurobiBound(model.lower[static_cast<std::size_t>(c)]);
                ub[c] = gurobiBound(model.upper[static_cast<std::size_t>(c)]);
                vtype[c] = model.variable_type[static_cast<std::size_t>(c)] ==
                    VariableType::Integer
                    ? GRB_INTEGER
                    : GRB_CONTINUOUS;
                names[c] = const_cast<char*>(
                    model.column_names[static_cast<std::size_t>(c)].c_str());
            }

            GRBmodel* raw_model = nullptr;
            checkGurobi(
                env,
                GRBnewmodel(
                    env,
                    &raw_model,
                    "tasr",
                    num_vars,
                    num_vars == 0 ? nullptr : const_cast<double*>(model.cost.data()),
                    num_vars == 0 ? nullptr : lb.data(),
                    num_vars == 0 ? nullptr : ub.data(),
                    num_vars == 0 ? nullptr : vtype.data(),
                    num_vars == 0 ? nullptr : names.data()),
                "creating the Gurobi model");
            GurobiModelPtr grb_model(raw_model);

            for (const SparseRow& row : model.rows) {
                addGurobiConstraint(env, grb_model.get(), row);
            }

            // Materialize all pending additions before setting variable attributes,
            // writing the LP file, or optimizing.
            checkGurobi(
                env,
                GRBupdatemodel(grb_model.get()),
                "updating the Gurobi model");

            if (mip_start && !mip_start->empty()) {
                for (std::size_t p = 0; p < mip_start->index.size(); ++p) {
                    checkGurobi(
                        env,
                        GRBsetdblattrelement(
                            grb_model.get(),
                            "Start",
                            mip_start->index[p],
                            mip_start->value[p]),
                        "setting a Gurobi MIP-start value");
                }
            }

            if (!model_file.empty()) {
                checkGurobi(
                    env,
                    GRBwrite(grb_model.get(), model_file.c_str()),
                    "writing model to '" + model_file + "'");
            }

            checkGurobi(
                env,
                GRBsetcallbackfunc(
                    grb_model.get(), gurobiInterruptCallback, nullptr),
                "setting interrupt callback");

            checkGurobi(
                env,
                GRBoptimize(grb_model.get()),
                "optimizing the model");

            SolverResult result;

            int status = 0;
            checkGurobi(
                env,
                GRBgetintattr(grb_model.get(), "Status", &status),
                "reading model status");

            int solution_count = 0;
            checkGurobi(
                env,
                GRBgetintattr(grb_model.get(), "SolCount", &solution_count),
                "reading solution count");

            result.status = gurobiStatusString(status);
            result.optimal = status == GRB_OPTIMAL;
            result.time_limit_reached = status == GRB_TIME_LIMIT;
            result.interrupted = status == GRB_INTERRUPTED || interruptRequested();
            result.feasible = solution_count > 0;

            if (result.feasible) {
                checkGurobi(
                    env,
                    GRBgetdblattr(grb_model.get(), "ObjVal", &result.objective),
                    "reading objective value");

                result.values.resize(model.cost.size());
                if (!result.values.empty()) {
                    checkGurobi(
                        env,
                        GRBgetdblattrarray(
                            grb_model.get(),
                            "X",
                            0,
                            static_cast<int>(result.values.size()),
                            result.values.data()),
                        "reading variable values");
                }
            }

            return result;
        }
#endif

    } // namespace

    int ModelBuilder::addColumn(std::string name,
        double cost,
        double lower,
        double upper,
        VariableType type) {
        const int column = static_cast<int>(model_.cost.size());
        model_.column_names.push_back(std::move(name));
        model_.cost.push_back(cost);
        model_.lower.push_back(lower);
        model_.upper.push_back(upper);
        model_.variable_type.push_back(type);
        return column;
    }

    void ModelBuilder::addRow(std::string name,
        double lower,
        double upper,
        std::vector<int> index,
        std::vector<double> value) {
        if (index.size() != value.size()) {
            throw std::logic_error("row index/value size mismatch");
        }
        model_.rows.push_back({
            std::move(name),
            lower,
            upper,
            std::move(index),
            std::move(value) });
    }

    LinearModel ModelBuilder::build() const {
        validateModel(model_);
        return model_;
    }

    std::size_t ModelBuilder::numColumns() const noexcept {
        return model_.cost.size();
    }

    std::size_t ModelBuilder::numRows() const noexcept {
        return model_.rows.size();
    }

    SolverType parseSolverType(const std::string& value) {
        std::string normalized = value;
        std::transform(
            normalized.begin(),
            normalized.end(),
            normalized.begin(),
            [](unsigned char c) {
                return static_cast<char>(std::tolower(c));
            });

        if (normalized == "highs") return SolverType::Highs;
        if (normalized == "gurobi" || normalized == "grb") {
            return SolverType::Gurobi;
        }

        throw std::invalid_argument(
            "unknown solver '" + value + "' (use highs or gurobi)");
    }

    const char* solverTypeName(SolverType type) noexcept {
        switch (type) {
        case SolverType::Highs: return "HiGHS";
        case SolverType::Gurobi: return "Gurobi";
        }
        return "unknown";
    }

    bool solverIsCompiled(SolverType type) noexcept {
        switch (type) {
        case SolverType::Highs:
#ifdef ROADEF_HAVE_HIGHS
            return true;
#else
            return false;
#endif

        case SolverType::Gurobi:
#ifdef ROADEF_HAVE_GUROBI
            return true;
#else
            return false;
#endif
        }

        return false;
    }

    SolverResult solveModel(SolverType type,
        const LinearModel& model,
        const BackendOptions& options,
        const MipStart* mip_start,
        const std::string& model_file) {
        validateModel(model);
        if (mip_start) {
            validateMipStart(*mip_start, model.cost.size());
        }

        switch (type) {
        case SolverType::Highs:
#ifdef ROADEF_HAVE_HIGHS
            return solveHighs(model, options, mip_start, model_file);
#else
            throw std::runtime_error(
                "HiGHS support was not compiled into this executable");
#endif

        case SolverType::Gurobi:
#ifdef ROADEF_HAVE_GUROBI
            return solveGurobi(model, options, mip_start, model_file);
#else
            throw std::runtime_error(
                "Gurobi support was not compiled into this executable");
#endif
        }

        throw std::logic_error("invalid solver type");
    }

} // namespace opt
