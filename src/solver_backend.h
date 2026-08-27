#pragma once

#include <cstddef>
#include <csignal>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace opt {

constexpr double kInfinity = std::numeric_limits<double>::infinity();

enum class SolverType {
    Highs,
    Gurobi
};

enum class VariableType {
    Continuous,
    Integer,
    Binary
};

struct SparseRow {
    std::string name;
    double lower = -kInfinity;
    double upper = kInfinity;
    std::vector<int> index;
    std::vector<double> value;
};

struct LinearModel {
    std::vector<std::string> column_names;
    std::vector<double> cost;
    std::vector<double> lower;
    std::vector<double> upper;
    std::vector<VariableType> variable_type;
    std::vector<SparseRow> rows;
};

class ModelBuilder {
public:
    int addColumn(std::string name,
                  double cost,
                  double lower,
                  double upper,
                  VariableType type = VariableType::Continuous);

    int addRow(std::string name,
                double lower,
                double upper,
                std::vector<int> index,
                std::vector<double> value);

    LinearModel build() const;

    std::size_t numColumns() const noexcept;
    std::size_t numRows() const noexcept;

private:
    LinearModel model_;
};

struct BackendOptions {
    double time_limit_seconds = 0.0;
    double mip_relative_gap = 0.0;
    bool solver_output = true;

    // Optional early-stop criterion for minimization MIPs:
    // terminate as soon as a feasible incumbent with objective strictly below
    // incumbent_stop_value is found.  This is used by the bottleneck helper,
    // where one such incumbent is already a certificate that the candidate
    // bottleneck is not unavoidable.
    bool stop_when_incumbent_below = false;
    double incumbent_stop_value = -kInfinity;
};

// Partial MIP start: only listed columns are supplied to the solver.
struct MipStart {
    std::vector<int> index;
    std::vector<double> value;

    bool empty() const noexcept { return index.empty(); }
};

struct SolverResult {
    bool feasible = false;
    bool optimal = false;
    bool time_limit_reached = false;
    bool interrupted = false;

    // True only when the backend deliberately terminated because the requested
    // incumbent threshold was reached.  This is distinct from Ctrl+C/time limit.
    bool incumbent_stop_triggered = false;

    double objective = kInfinity;
    std::vector<double> values;
    std::string status;
};

SolverType parseSolverType(const std::string& value);
const char* solverTypeName(SolverType type) noexcept;
bool solverIsCompiled(SolverType type) noexcept;

// Process-wide Ctrl+C flag. The signal handler only sets this flag; the
// active solver backend observes it through its native callback mechanism.
void requestInterrupt() noexcept;
void clearInterrupt() noexcept;
bool interruptRequested() noexcept;


class SolverSession {
public:
    virtual ~SolverSession() = default;

    virtual SolverResult optimize(
        const BackendOptions& options,
        const MipStart* mip_start = nullptr,
        const std::string& model_file = {}) = 0;

    virtual void setColumnBounds(int column, double lower, double upper) = 0;
    virtual void setObjectiveCoefficient(int column, double value) = 0;
    virtual void setRowBounds(int row, double lower, double upper) = 0;
    virtual void setCoefficient(int row, int column, double value) = 0;

    virtual std::size_t numColumns() const noexcept = 0;
    virtual std::size_t numRows() const noexcept = 0;
};

std::unique_ptr<SolverSession> createSolverSession(
    SolverType type,
    const LinearModel& model);

SolverResult solveModel(SolverType type,
                        const LinearModel& model,
                        const BackendOptions& options,
                        const MipStart* mip_start = nullptr,
                        const std::string& model_file = {});

} // namespace opt
