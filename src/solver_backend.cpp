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
                    (type == VariableType::Integer || type == VariableType::Binary)
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

        class HighsSession final : public SolverSession {
        public:
            explicit HighsSession(const LinearModel& model)
                : num_columns_(model.cost.size()),
                  num_rows_(model.rows.size()) {

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
                        (type == VariableType::Integer || type == VariableType::Binary)
                        ? HighsVarType::kInteger
                        : HighsVarType::kContinuous);
                }

                lp.row_names_.reserve(model.rows.size());
                lp.row_lower_.reserve(model.rows.size());
                lp.row_upper_.reserve(model.rows.size());

                std::vector<std::vector<std::pair<HighsInt, double>>> columns(
                    model.cost.size());

                for (HighsInt r = 0; r < lp.num_row_; ++r) {
                    const SparseRow& row =
                        model.rows.at(static_cast<std::size_t>(r));
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

                checkHighsStatus(
                    highs_.setCallback(highsInterruptCallback, nullptr),
                    "setting interrupt callback");
                checkHighsStatus(
                    highs_.startCallback(
                        static_cast<int>(HighsCallbackType::kCallbackMipInterrupt)),
                    "starting MIP interrupt callback");
                checkHighsStatus(
                    highs_.passModel(std::move(lp)),
                    "passing persistent model");
            }

            SolverResult optimize(
                const BackendOptions& options,
                const MipStart* mip_start,
                const std::string& model_file) override {

                checkHighsStatus(
                    highs_.setOptionValue("output_flag", options.solver_output),
                    "setting output_flag");
                checkHighsStatus(
                    highs_.setOptionValue("mip_rel_gap", options.mip_relative_gap),
                    "setting mip_rel_gap");
                checkHighsStatus(
                    highs_.setOptionValue(
                        "time_limit",
                        options.time_limit_seconds > 0.0
                            ? options.time_limit_seconds
                            : kHighsInf),
                    "setting time_limit");

                if (mip_start && !mip_start->empty()) {
                    validateMipStart(*mip_start, num_columns_);
                    std::vector<HighsInt> indices;
                    indices.reserve(mip_start->index.size());
                    for (const int index : mip_start->index) {
                        indices.push_back(static_cast<HighsInt>(index));
                    }
                    checkHighsStatus(
                        highs_.setSolution(
                            static_cast<HighsInt>(indices.size()),
                            indices.data(),
                            mip_start->value.data()),
                        "setting MIP start");
                }

                if (!model_file.empty()) {
                    checkHighsStatus(
                        highs_.writeModel(model_file),
                        "writing persistent model");
                }

                checkHighsStatus(highs_.run(), "solving persistent model");

                SolverResult result;
                result.status =
                    highs_.modelStatusToString(highs_.getModelStatus());
                result.optimal =
                    highs_.getModelStatus() == HighsModelStatus::kOptimal;
                result.time_limit_reached =
                    highs_.getModelStatus() == HighsModelStatus::kTimeLimit;
                result.interrupted =
                    highs_.getModelStatus() == HighsModelStatus::kInterrupt ||
                    interruptRequested();
                result.feasible =
                    highsHasUsableSolution(highs_) &&
                    highs_.getSolution().value_valid;

                if (result.feasible) {
                    result.objective = highs_.getObjectiveValue();
                    result.values = highs_.getSolution().col_value;
                }
                return result;
            }

            void setColumnBounds(
                int column, double lower, double upper) override {
                checkHighsStatus(
                    highs_.changeColBounds(
                        static_cast<HighsInt>(column), lower, upper),
                    "changing persistent column bounds");
            }

            void setObjectiveCoefficient(
                int column, double value) override {
                checkHighsStatus(
                    highs_.changeColCost(
                        static_cast<HighsInt>(column), value),
                    "changing persistent objective coefficient");
            }

            void setRowBounds(
                int row, double lower, double upper) override {
                checkHighsStatus(
                    highs_.changeRowBounds(
                        static_cast<HighsInt>(row), lower, upper),
                    "changing persistent row bounds");
            }

            void setCoefficient(
                int row, int column, double value) override {
                checkHighsStatus(
                    highs_.changeCoeff(
                        static_cast<HighsInt>(row),
                        static_cast<HighsInt>(column),
                        value),
                    "changing persistent matrix coefficient");
            }

            std::size_t numColumns() const noexcept override {
                return num_columns_;
            }

            std::size_t numRows() const noexcept override {
                return num_rows_;
            }

        private:
            Highs highs_;
            std::size_t num_columns_ = 0;
            std::size_t num_rows_ = 0;
        };
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

        struct GurobiCallbackData {
            bool stop_when_incumbent_below = false;
            double incumbent_stop_value = -kInfinity;
            bool incumbent_stop_triggered = false;
        };

        int __stdcall gurobiInterruptCallback(
            GRBmodel* model,
            void* cbdata,
            int where,
            void* usrdata) {

            auto* data = static_cast<GurobiCallbackData*>(usrdata);

            if (interruptRequested()) {
                GRBterminate(model);
                return 0;
            }

            // For a minimization MIP, one feasible incumbent below the helper's
            // rejection threshold is enough to prove that the tested bottleneck
            // is NOT unavoidable.  There is no need to prove optimality.
            if (data &&
                data->stop_when_incumbent_below &&
                where == GRB_CB_MIPSOL) {

                double incumbent_objective = 0.0;
                const int error = GRBcbget(
                    cbdata,
                    where,
                    GRB_CB_MIPSOL_OBJ,
                    &incumbent_objective);

                if (error == 0 &&
                    incumbent_objective < data->incumbent_stop_value) {
                    data->incumbent_stop_triggered = true;
                    GRBterminate(model);
                }
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

			// Gurobi's MIPFocus=3 setting is recommended for lexicographic optimization.
            checkGurobi(
                env,
                GRBsetintparam(
                    env,
                    "MIPFocus",
                    3),
                "setting MIPFocus");
            
			// Gurobi's Heuristics=0.01 setting is recommended for lexicographic optimization.
            checkGurobi(
                env,
                GRBsetdblparam(
                    env,
                    "Heuristics",
                    0.01),
                "setting Heuristics");
            
            checkGurobi(
                env,
                GRBsetintparam(
                    env,
                    "Cuts",
                    2),
                "setting Cuts");
            

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
                switch (model.variable_type[static_cast<std::size_t>(c)]) {
                case VariableType::Continuous:
                    vtype[c] = GRB_CONTINUOUS;
                    break;

                case VariableType::Integer:
                    vtype[c] = GRB_INTEGER;
                    break;

                case VariableType::Binary:
                    vtype[c] = GRB_BINARY;
                    break;
                }
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

            GurobiCallbackData callback_data;
            callback_data.stop_when_incumbent_below =
                options.stop_when_incumbent_below;
            callback_data.incumbent_stop_value =
                options.incumbent_stop_value;

            checkGurobi(
                env,
                GRBsetcallbackfunc(
                    grb_model.get(), gurobiInterruptCallback, &callback_data),
                "setting interrupt/early-stop callback");

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
            result.incumbent_stop_triggered =
                callback_data.incumbent_stop_triggered;
            result.interrupted =
                (status == GRB_INTERRUPTED || interruptRequested()) &&
                !result.incumbent_stop_triggered;
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

        class GurobiSession final : public SolverSession {
        public:
            explicit GurobiSession(const LinearModel& model)
                : env_(sharedGurobiEnvironment()),
                  num_columns_(model.cost.size()),
                  num_rows_(model.rows.size()),
                  native_rows_(model.rows.size()) {

                const int num_vars = static_cast<int>(model.cost.size());

                std::vector<double> lb(num_vars);
                std::vector<double> ub(num_vars);
                std::vector<char> vtype(num_vars);
                std::vector<char*> names(num_vars);

                for (int c = 0; c < num_vars; ++c) {
                    lb[c] = gurobiBound(model.lower[static_cast<std::size_t>(c)]);
                    ub[c] = gurobiBound(model.upper[static_cast<std::size_t>(c)]);

                    switch (model.variable_type[static_cast<std::size_t>(c)]) {
                    case VariableType::Continuous:
                        vtype[c] = GRB_CONTINUOUS;
                        break;
                    case VariableType::Integer:
                        vtype[c] = GRB_INTEGER;
                        break;
                    case VariableType::Binary:
                        vtype[c] = GRB_BINARY;
                        break;
                    }

                    names[c] = const_cast<char*>(
                        model.column_names[static_cast<std::size_t>(c)].c_str());
                }

                GRBmodel* raw_model = nullptr;
                checkGurobi(
                    env_,
                    GRBnewmodel(
                        env_,
                        &raw_model,
                        "tasr-persistent",
                        num_vars,
                        num_vars == 0
                            ? nullptr
                            : const_cast<double*>(model.cost.data()),
                        num_vars == 0 ? nullptr : lb.data(),
                        num_vars == 0 ? nullptr : ub.data(),
                        num_vars == 0 ? nullptr : vtype.data(),
                        num_vars == 0 ? nullptr : names.data()),
                    "creating persistent Gurobi model");
                model_.reset(raw_model);

                int next_native_row = 0;
                for (std::size_t r = 0; r < model.rows.size(); ++r) {
                    const SparseRow& row = model.rows[r];
                    std::vector<int> indices = row.index;
                    std::vector<double> values = row.value;

                    const bool finite_lower = std::isfinite(row.lower);
                    const bool finite_upper = std::isfinite(row.upper);

                    auto add_one = [&](char sense,
                                       double rhs,
                                       const std::string& name) {
                        checkGurobi(
                            env_,
                            GRBaddconstr(
                                model_.get(),
                                static_cast<int>(indices.size()),
                                indices.empty() ? nullptr : indices.data(),
                                values.empty() ? nullptr : values.data(),
                                sense,
                                rhs,
                                name.c_str()),
                            "adding persistent constraint '" + name + "'");
                        native_rows_[r].push_back(next_native_row++);
                    };

                    if (!finite_lower && !finite_upper) {
                        continue;
                    }

                    if (finite_lower && finite_upper &&
                        std::abs(row.lower - row.upper) <= 1e-12) {
                        add_one(GRB_EQUAL, row.lower, row.name);
                    }
                    else {
                        if (finite_lower) {
                            add_one(
                                GRB_GREATER_EQUAL,
                                row.lower,
                                row.name + "_lb");
                        }
                        if (finite_upper) {
                            add_one(
                                GRB_LESS_EQUAL,
                                row.upper,
                                row.name + "_ub");
                        }
                    }
                }

                checkGurobi(
                    env_,
                    GRBupdatemodel(model_.get()),
                    "updating persistent Gurobi model");

            }

            SolverResult optimize(
                const BackendOptions& options,
                const MipStart* mip_start,
                const std::string& model_file) override {

                checkGurobi(
                    env_,
                    GRBsetintparam(
                        env_,
                        "OutputFlag",
                        options.solver_output ? 1 : 0),
                    "setting OutputFlag");

                checkGurobi(
                    env_,
                    GRBsetdblparam(
                        env_,
                        "MIPGap",
                        options.mip_relative_gap),
                    "setting MIPGap");

                checkGurobi(
                    env_,
                    GRBsetintparam(env_, "MIPFocus", 3),
                    "setting MIPFocus");

                checkGurobi(
                    env_,
                    GRBsetdblparam(env_, "Heuristics", 0.01),
                    "setting Heuristics");

                checkGurobi(
                    env_,
                    GRBsetintparam(env_, "Cuts", 2),
                    "setting Cuts");

                checkGurobi(
                    env_,
                    GRBsetdblparam(
                        env_,
                        "TimeLimit",
                        options.time_limit_seconds > 0.0
                            ? options.time_limit_seconds
                            : GRB_INFINITY),
                    "setting TimeLimit");

                // Explicit user MIP starts persist on a Gurobi model across
                // optimize() calls.  For our persistent T-ASR model that is not
                // what we want: the common-routing start is intended only for the
                // first full solve.  Later solves should reuse Gurobi's incumbent
                // from the previous optimization automatically.
                //
                // Therefore discard every previously stored explicit start before
                // each solve.  This does NOT reset the model or its incumbent.
                checkGurobi(
                    env_,
                    GRBsetintattr(
                        model_.get(),
                        "NumStart",
                        0),
                    "discarding stale persistent Gurobi MIP starts");

                checkGurobi(
                    env_,
                    GRBupdatemodel(model_.get()),
                    "updating persistent Gurobi model after clearing MIP starts");
                GRBenv* model_env = GRBgetenv(model_.get());

                checkGurobi(
                    model_env,
                    GRBsetintparam(
                        model_env,
                        "DualReductions",
                        0),
                    "setting DualReductions");
                // If this particular solve intentionally supplies a fresh start,
                // create exactly one new start and populate only the supplied
                // entries.  Unspecified entries remain undefined.
                if (mip_start && !mip_start->empty()) {
                    validateMipStart(*mip_start, num_columns_);

                    checkGurobi(
                        env_,
                        GRBsetintattr(
                            model_.get(),
                            "NumStart",
                            1),
                        "creating fresh persistent Gurobi MIP start");

                    checkGurobi(
                        env_,
                        GRBsetintparam(
                            env_,
                            "StartNumber",
                            0),
                        "selecting persistent Gurobi MIP start 0");

                    for (std::size_t p = 0;
                         p < mip_start->index.size();
                         ++p) {
                        checkGurobi(
                            env_,
                            GRBsetdblattrelement(
                                model_.get(),
                                "Start",
                                mip_start->index[p],
                                mip_start->value[p]),
                            "setting persistent Gurobi MIP-start value");
                    }
                }

                checkGurobi(
                    env_,
                    GRBupdatemodel(model_.get()),
                    "updating persistent Gurobi model before solve");

                if (!model_file.empty()) {
                    checkGurobi(
                        env_,
                        GRBwrite(model_.get(), model_file.c_str()),
                        "writing persistent model to '" + model_file + "'");
                }

                GurobiCallbackData callback_data;
                callback_data.stop_when_incumbent_below =
                    options.stop_when_incumbent_below;
                callback_data.incumbent_stop_value =
                    options.incumbent_stop_value;

                checkGurobi(
                    env_,
                    GRBsetcallbackfunc(
                        model_.get(),
                        gurobiInterruptCallback,
                        &callback_data),
                    "setting persistent interrupt/early-stop callback");

                checkGurobi(
                    env_,
                    GRBoptimize(model_.get()),
                    "optimizing persistent Gurobi model");

                SolverResult result;

                int status = 0;
                checkGurobi(
                    env_,
                    GRBgetintattr(model_.get(), "Status", &status),
                    "reading persistent model status");

                int solution_count = 0;
                checkGurobi(
                    env_,
                    GRBgetintattr(model_.get(), "SolCount", &solution_count),
                    "reading persistent solution count");

                result.status = gurobiStatusString(status);
                result.optimal = status == GRB_OPTIMAL;
                result.time_limit_reached = status == GRB_TIME_LIMIT;
                result.incumbent_stop_triggered =
                    callback_data.incumbent_stop_triggered;
                result.interrupted =
                    (status == GRB_INTERRUPTED || interruptRequested()) &&
                    !result.incumbent_stop_triggered;
                result.feasible = solution_count > 0;

                if (result.feasible) {
                    checkGurobi(
                        env_,
                        GRBgetdblattr(
                            model_.get(), "ObjVal", &result.objective),
                        "reading persistent objective value");

                    result.values.resize(num_columns_);
                    if (!result.values.empty()) {
                        checkGurobi(
                            env_,
                            GRBgetdblattrarray(
                                model_.get(),
                                "X",
                                0,
                                static_cast<int>(result.values.size()),
                                result.values.data()),
                            "reading persistent variable values");
                    }
                }

                return result;
            }

            void setColumnBounds(
                int column, double lower, double upper) override {
                checkColumn(column);
                checkGurobi(
                    env_,
                    GRBsetdblattrelement(
                        model_.get(), "LB", column, gurobiBound(lower)),
                    "changing persistent Gurobi lower bound");
                checkGurobi(
                    env_,
                    GRBsetdblattrelement(
                        model_.get(), "UB", column, gurobiBound(upper)),
                    "changing persistent Gurobi upper bound");
            }

            void setObjectiveCoefficient(
                int column, double value) override {
                checkColumn(column);
                checkGurobi(
                    env_,
                    GRBsetdblattrelement(
                        model_.get(), "Obj", column, value),
                    "changing persistent Gurobi objective");
            }

            void setRowBounds(
                int row, double lower, double upper) override {
                checkRow(row);
                const auto& native = native_rows_[static_cast<std::size_t>(row)];

                if (native.empty()) {
                    throw std::logic_error(
                        "cannot activate a row that had no finite bounds");
                }

                const bool finite_lower = std::isfinite(lower);
                const bool finite_upper = std::isfinite(upper);
                const bool equality =
                    finite_lower && finite_upper &&
                    std::abs(lower - upper) <= 1e-12;

                if (native.size() == 1) {
                    char sense = 0;
                    double rhs = 0.0;

                    if (equality) {
                        sense = GRB_EQUAL;
                        rhs = lower;
                    }
                    else if (finite_lower && !finite_upper) {
                        sense = GRB_GREATER_EQUAL;
                        rhs = lower;
                    }
                    else if (!finite_lower && finite_upper) {
                        sense = GRB_LESS_EQUAL;
                        rhs = upper;
                    }
                    else if (!finite_lower && !finite_upper) {
                        // A completely free generic row has no direct Gurobi
                        // constraint representation. Keep the native row in place
                        // so that persistent row indices remain stable, but make it
                        // redundant.
                        //
                        //     lhs <= GRB_INFINITY
                        //
                        // This is equivalent to removing the row for our purposes,
                        // without deleting a native constraint and renumbering the
                        // persistent model.
                        sense = GRB_LESS_EQUAL;
                        rhs = GRB_INFINITY;
                    }
                    else {
                        // A single native Gurobi constraint cannot represent a
                        // genuine two-sided range lower <= lhs <= upper.
                        throw std::logic_error(
                            "cannot change single persistent Gurobi row to "
                            "a finite two-sided range");
                    }

                    checkGurobi(
                        env_,
                        GRBsetcharattrelement(
                            model_.get(), "Sense", native[0], sense),
                        "changing persistent Gurobi row sense");

                    checkGurobi(
                        env_,
                        GRBsetdblattrelement(
                            model_.get(), "RHS", native[0], rhs),
                        "changing persistent Gurobi row RHS");

                    return;
                }

                if (native.size() == 2 &&
                    finite_lower && finite_upper && !equality) {
                    checkGurobi(
                        env_,
                        GRBsetcharattrelement(
                            model_.get(),
                            "Sense",
                            native[0],
                            GRB_GREATER_EQUAL),
                        "changing persistent lower-row sense");
                    checkGurobi(
                        env_,
                        GRBsetdblattrelement(
                            model_.get(), "RHS", native[0], lower),
                        "changing persistent lower-row RHS");
                    checkGurobi(
                        env_,
                        GRBsetcharattrelement(
                            model_.get(),
                            "Sense",
                            native[1],
                            GRB_LESS_EQUAL),
                        "changing persistent upper-row sense");
                    checkGurobi(
                        env_,
                        GRBsetdblattrelement(
                            model_.get(), "RHS", native[1], upper),
                        "changing persistent upper-row RHS");
                    return;
                }

                throw std::logic_error(
                    "unsupported persistent Gurobi row-bound change");
            }

            void setCoefficient(
                int row, int column, double value) override {
                checkRow(row);
                checkColumn(column);

                const auto& native =
                    native_rows_[static_cast<std::size_t>(row)];
                for (const int native_row : native) {
                    int constr = native_row;
                    int var = column;
                    double coefficient = value;
                    checkGurobi(
                        env_,
                        GRBchgcoeffs(
                            model_.get(),
                            1,
                            &constr,
                            &var,
                            &coefficient),
                        "changing persistent Gurobi matrix coefficient");
                }
            }

            std::size_t numColumns() const noexcept override {
                return num_columns_;
            }

            std::size_t numRows() const noexcept override {
                return num_rows_;
            }

        private:
            void checkColumn(int column) const {
                if (column < 0 ||
                    static_cast<std::size_t>(column) >= num_columns_) {
                    throw std::out_of_range(
                        "persistent solver column index out of range");
                }
            }

            void checkRow(int row) const {
                if (row < 0 ||
                    static_cast<std::size_t>(row) >= num_rows_) {
                    throw std::out_of_range(
                        "persistent solver row index out of range");
                }
            }

            GRBenv* env_ = nullptr;
            GurobiModelPtr model_;
            std::size_t num_columns_ = 0;
            std::size_t num_rows_ = 0;
            std::vector<std::vector<int>> native_rows_;
        };
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

    int ModelBuilder::addRow(std::string name,
        double lower,
        double upper,
        std::vector<int> index,
        std::vector<double> value) {
        if (index.size() != value.size()) {
            throw std::logic_error("row index/value size mismatch");
        }
        const int row = static_cast<int>(model_.rows.size());
        model_.rows.push_back({
            std::move(name),
            lower,
            upper,
            std::move(index),
            std::move(value) });
        return row;
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


    std::unique_ptr<SolverSession> createSolverSession(
        SolverType type,
        const LinearModel& model) {

        validateModel(model);

        switch (type) {
        case SolverType::Highs:
#ifdef ROADEF_HAVE_HIGHS
            return std::make_unique<HighsSession>(model);
#else
            throw std::runtime_error(
                "HiGHS support was not compiled into this executable");
#endif

        case SolverType::Gurobi:
#ifdef ROADEF_HAVE_GUROBI
            return std::make_unique<GurobiSession>(model);
#else
            throw std::runtime_error(
                "Gurobi support was not compiled into this executable");
#endif
        }

        throw std::logic_error("invalid solver type");
    }

    SolverResult solveModel(SolverType type,
        const LinearModel& model,
        const BackendOptions& options,
        const MipStart* mip_start,
        const std::string& model_file) {

        auto session = createSolverSession(type, model);
        return session->optimize(options, mip_start, model_file);
    }

} // namespace opt
