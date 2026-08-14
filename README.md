# ROADEF T-ASR solver: HiGHS + Gurobi

This version keeps the mathematical T-ASR model solver-independent and can solve the same model with either **HiGHS** or **Gurobi** at runtime.

The Gurobi backend uses **only the Gurobi C API (`gurobi_c.h`)**. It does not use `gurobi_c++.h`. This is intentional for the existing MSYS2 UCRT64 / GCC setup.

## Project structure

```text
Solver/
├── CMakeLists.txt
├── README.md
├── src/
│   ├── main.cpp
│   ├── instance_data.cpp
│   ├── instance_data.h
│   ├── split_coefficients.cpp
│   ├── split_coefficients.h
│   ├── tasr_solver.cpp
│   ├── tasr_solver.h
│   ├── solver_backend.cpp
│   └── solver_backend.h
├── external/
│   └── networktools-main/networktools/   # keep your existing folder
└── instances/                            # keep your existing instances
```

The ZIP contains the files that were changed or supplied in the conversation. Keep your existing `external/networktools-main/networktools` and `instances` directories in the project.

## Environment

In `~/.bashrc`, keep your permanent Gurobi setup, for example:

```bash
export GUROBI_HOME="/c/gurobi1302/win64"
export PATH="$GUROBI_HOME/bin:$PATH"
```

After changing `.bashrc`:

```bash
source ~/.bashrc
```

Check:

```bash
echo "$GUROBI_HOME"
which gurobi_cl
```

## Configure and build both backends

From the MSYS2 **UCRT64** shell in the project root:

```bash
rm -rf build
cmake -B build -G "MinGW Makefiles"
cmake --build build -j
```

Both backends are enabled by default.

The CMake file derives the active UCRT64 prefix from the selected compiler, so an installed HiGHS package under `C:/msys64/ucrt64` should be found automatically.

If your local HiGHS installation still needs an explicit package path, configure as before, for example:

```bash
cmake -B build -G "MinGW Makefiles" \
  -Dhighs_DIR="C:/msys64/ucrt64/lib/cmake/highs"
cmake --build build -j
```

## Run with HiGHS

```bash
./build/solver.exe \
  --instance setA-01 \
  --iterations 10 \
  --time 590 \
  --solver highs
```

## Run with Gurobi

```bash
./build/solver.exe \
  --instance setA-01 \
  --iterations 10 \
  --time 590 \
  --solver gurobi
```

Short form `grb` is also accepted:

```bash
./build/solver.exe --instance setA-01 --iterations 1 --solver grb
```

## Total time limit and Ctrl+C

`--time` is a **single wall-clock budget for the whole solver process**, not a per-iteration limit. For example, `--time 590` gives the common-routing attempt, all lexicographic model construction, and all solver calls together 590 seconds. Before every backend call, only the remaining budget is passed to HiGHS/Gurobi.

If the limit is reached during a lexicographic iteration, the solver keeps the best feasible incumbent returned by that iteration and compares it lexicographically (sorted load vector) against the best solution from the last completed iteration. The better feasible routing is written to the normal `*-srpaths.json` output.

Pressing **Ctrl+C once** requests a graceful stop. HiGHS is interrupted through its MIP interrupt callback and Gurobi through `GRBterminate` from a callback. The same incumbent-selection logic is then used and the best available feasible routing is written before the program exits. If no full-model incumbent exists yet, a feasible common-routing solution is used as fallback when available.

The default total budget is 590 seconds. Set `--time 0` for no time limit.

## Export LP files

HiGHS:

```bash
./build/solver.exe \
  --instance setA-01 \
  --iterations 1 \
  --solver highs \
  --write-models 1 \
  --model-prefix compare
```

Gurobi:

```bash
./build/solver.exe \
  --instance setA-01 \
  --iterations 1 \
  --solver gurobi \
  --write-models 1 \
  --model-prefix compare
```

This creates solver-labelled model files such as:

```text
compare-HiGHS-lex-1.lp
compare-Gurobi-lex-1.lp
```

Gurobi may write a two-sided generic row as two one-sided rows (`_lb` and `_ub`). This is mathematically equivalent to the single ranged row used by the generic model/HiGHS representation.

## Build only one backend

Only HiGHS:

```bash
cmake -B build-highs -G "MinGW Makefiles" \
  -DENABLE_HIGHS=ON \
  -DENABLE_GUROBI=OFF
cmake --build build-highs -j
```

Only Gurobi:

```bash
cmake -B build-gurobi -G "MinGW Makefiles" \
  -DENABLE_HIGHS=OFF \
  -DENABLE_GUROBI=ON
cmake --build build-gurobi -j
```

## Architecture

The mathematical code in `tasr_solver.cpp` builds one solver-neutral `LinearModel` using `ModelBuilder`.

```text
T-ASR formulation
       |
       v
 generic LinearModel
      / \
     /   \
 HiGHS   Gurobi C API
```

`solver_backend.cpp` is the only file that knows about `Highs.h` or `gurobi_c.h`.

Warm starts are also represented generically as column-index/value pairs. The HiGHS backend maps these to `Highs::setSolution`; the Gurobi backend maps them to the variable `Start` attribute.

## First verification I recommend

Before running many lexicographic iterations, solve one iteration with both solvers:

```bash
./build/solver.exe --instance setA-01 --iterations 1 --solver highs
./build/solver.exe --instance setA-01 --iterations 1 --solver gurobi
```

The maximum load / first lexicographic objective should agree within normal numerical tolerances. Then test several iterations.
