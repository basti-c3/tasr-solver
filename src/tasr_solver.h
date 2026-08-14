#pragma once

#include "instance_data.h"
#include "split_coefficients.h"
#include "solver_backend.h"

#include <string>
#include <vector>

struct TasrSolverOptions {
    opt::SolverType solver = opt::SolverType::Highs;
    int max_lexicographic_iterations = -1; // -1 means all arc/time loads
    double time_limit_seconds = 590.0;     // total wall-clock budget for the whole solve; 0 = unlimited
    double mip_relative_gap = 0.0;
    bool solver_output = true;
    bool write_models = false;
    std::string model_prefix = "tasr";
};

struct TasrSolution {
    std::vector<std::vector<std::vector<int>>> waypoints; // [d][t]
    std::vector<double> sorted_loads;
    std::vector<double> lexicographic_objectives; // sums of k largest loads
    int completed_lexicographic_iterations = 0;
};

class TasrSolver {
public:
    TasrSolver(const InstanceData& instance,
               const SplitCoefficients& split_coefficients,
               TasrSolverOptions options = {});

    TasrSolution solve();

private:
    const InstanceData& instance_;
    const SplitCoefficients& splits_;
    TasrSolverOptions options_;
};

void writeSolutionJson(const TasrSolution& solution,
                       const std::string& filename);
