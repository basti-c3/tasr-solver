#pragma once

#include <cstddef>
#include <csignal>
#include <limits>
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
    Integer
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

    void addRow(std::string name,
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

SolverResult solveModel(SolverType type,
                        const LinearModel& model,
                        const BackendOptions& options,
                        const MipStart* mip_start = nullptr,
                        const std::string& model_file = {});

} // namespace opt
