#include "instance_data.h"
#include "split_coefficients.h"
#include "tasr_solver.h"

#include <csignal>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

constexpr const char* DATA_DIR = "instances/";

namespace {

    void handleInterrupt(int /*signal*/) {
        opt::requestInterrupt();
    }

    bool parseBooleanArgument(const std::string& value) {
        if (value == "1" || value == "true" || value == "TRUE" ||
            value == "yes" || value == "YES" || value == "on" ||
            value == "ON") {
            return true;
        }

        if (value == "0" || value == "false" || value == "FALSE" ||
            value == "no" || value == "NO" || value == "off" ||
            value == "OFF") {
            return false;
        }

        throw std::invalid_argument(
            "Boolean argument must be one of: "
            "0, 1, false, true, no, yes, off, on"
        );
    }

    void printUsage(const char* executable_name) {
        std::cout
            << "Usage:\n"
            << "  " << executable_name
            << " [--instance NAME] [--output FILE] [--iterations N]"
            << " [--time SECONDS] [--solver highs|gurobi]"
            << " [--write-models 0|1] [--model-prefix PREFIX]\n\n"
            << "  --time is the TOTAL wall-clock budget across warm start and all lexicographic iterations.\n"
            << "  Ctrl+C requests a graceful stop and writes the best feasible incumbent found.\n\n"
            << "Example:\n"
            << "  " << executable_name
            << " --instance setA-01 --iterations 10 --time 590 --solver gurobi\n";
    }

} // namespace

int main(int argc, char** argv) {
    nt::Logging::init(argc, argv);

    opt::clearInterrupt();
    std::signal(SIGINT, handleInterrupt);

    try {
        /*
         * General defaults.
         *
         * These are used in Release mode and also act as fallbacks if
         * no Debug-specific compile definitions are available.
         */
        std::string instance_name = "setA-01";
        std::string output_file;
        TasrSolverOptions options;

        /*
         * Debug-build defaults supplied by CMake.
         */
#ifdef DEBUG_DEFAULT_INSTANCE
        instance_name = DEBUG_DEFAULT_INSTANCE;
#endif

#ifdef DEBUG_DEFAULT_OUTPUT_FILE
        output_file = DEBUG_DEFAULT_OUTPUT_FILE;
#endif

#ifdef DEBUG_DEFAULT_ITERATIONS
        options.max_lexicographic_iterations =
            DEBUG_DEFAULT_ITERATIONS;
#endif

#ifdef DEBUG_DEFAULT_TIME_LIMIT
        options.time_limit_seconds =
            DEBUG_DEFAULT_TIME_LIMIT;
#endif

#ifdef DEBUG_DEFAULT_SOLVER
        options.solver = opt::parseSolverType(DEBUG_DEFAULT_SOLVER);
#endif

#ifdef WRITE_LP_MODELS
        options.write_models = true;
#endif

#ifdef DEBUG_DEFAULT_MODEL_PREFIX
        options.model_prefix =
            DEBUG_DEFAULT_MODEL_PREFIX;
#endif

        /*
         * If no explicit output-file default was supplied, construct it
         * from the instance name.
         */
        if (output_file.empty()) {
            output_file = instance_name + "-srpaths.json";
        }

        /*
  * Named command-line arguments override all defaults.
  */
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];

            if (arg == "--help" || arg == "-h") {
                printUsage(argv[0]);
                return 0;
            }

            if (arg == "--instance") {
                if (i + 1 >= argc) {
                    throw std::invalid_argument(
                        "--instance requires a value"
                    );
                }

                instance_name = argv[++i];

                // If the instance changes and no explicit output is given later,
                // derive a matching output filename.
                output_file =
                    std::string(DATA_DIR)
                    + instance_name
                    + "-srpaths.json";

                continue;
            }

            if (arg == "--output") {
                if (i + 1 >= argc) {
                    throw std::invalid_argument(
                        "--output requires a value"
                    );
                }

                output_file = argv[++i];
                continue;
            }

            if (arg == "--iterations" || arg == "--iter") {
                if (i + 1 >= argc) {
                    throw std::invalid_argument(
                        "--iterations requires a value"
                    );
                }

                options.max_lexicographic_iterations =
                    std::stoi(argv[++i]);

                continue;
            }

            if (arg == "--time" ||
                arg == "--time-limit" ||
                arg == "--timelimit") {

                if (i + 1 >= argc) {
                    throw std::invalid_argument(
                        "--time requires a value"
                    );
                }

                options.time_limit_seconds =
                    std::stod(argv[++i]);

                continue;
            }

            if (arg == "--solver") {
                if (i + 1 >= argc) {
                    throw std::invalid_argument(
                        "--solver requires highs or gurobi");
                }
                options.solver = opt::parseSolverType(argv[++i]);
                continue;
            }

            if (arg == "--write-models") {
                if (i + 1 >= argc) {
                    throw std::invalid_argument(
                        "--write-models requires a value"
                    );
                }

                options.write_models =
                    parseBooleanArgument(argv[++i]);

                continue;
            }

            if (arg == "--model-prefix") {
                if (i + 1 >= argc) {
                    throw std::invalid_argument(
                        "--model-prefix requires a value"
                    );
                }

                options.model_prefix = argv[++i];
                continue;
            }

            throw std::invalid_argument(
                "Unknown command-line argument: " + arg
            );
        }
        const std::string base =
            std::string(DATA_DIR) + instance_name;

        InstanceData instance;
        instance.load(
            base + "-net.json",
            base + "-tm.json",
            base + "-scenario.json"
        );

        SplitCoefficients split_coefficients(instance);
        split_coefficients.compute();

        std::cout
            << "Loaded instance: " << instance_name << '\n'
            << "  nodes:            "
            << instance.numNodes() << '\n'
            << "  links:            "
            << instance.numArcs() << '\n'
            << "  demands:          "
            << instance.numDemands() << '\n'
            << "  time slots:       "
            << instance.num_time_slots << '\n'
            << "  max segments:     "
            << instance.max_segments << '\n'
            << "  nonzero r-values: "
            << split_coefficients.numberOfNonzeroCoefficients()
            << '\n'
            << "  solver:           "
            << opt::solverTypeName(options.solver) << '\n'
            << "  lex iterations:   "
            << options.max_lexicographic_iterations << '\n'
            << "  total time limit: "
            << options.time_limit_seconds << " seconds\n"
            << "  write LP models:  "
            << (options.write_models ? "yes" : "no") << '\n';

        if (options.write_models) {
            std::cout
                << "  model prefix:     "
                << options.model_prefix << '\n';
        }

        std::cout << '\n';

        if (!opt::solverIsCompiled(options.solver)) {
            throw std::runtime_error(
                std::string(opt::solverTypeName(options.solver)) +
                " support is not compiled into this executable");
        }

        TasrSolver solver(
            instance,
            split_coefficients,
            options
        );

        const TasrSolution solution = solver.solve();

        writeSolutionJson(solution, output_file);

        std::cout
            << "\nWrote solution to "
            << output_file << '\n'
            << "Completed lexicographic iterations: "
            << solution.completed_lexicographic_iterations
            << '\n';

        if (!solution.sorted_loads.empty()) {
            std::cout
                << "Maximum link utilization: "
                << solution.sorted_loads.front()
                << '\n';
        }

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr
            << "Failed to run solver: "
            << error.what()
            << '\n';

        printUsage(argv[0]);
        return 1;
    }
}