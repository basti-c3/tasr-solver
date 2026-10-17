# ROADEF/EURO 2026 – T-ASR Solver

An experimental **C++20 mixed-integer programming approach** to the
[ROADEF/EURO 2026 “Keep the Flow!” challenge](https://roadef.org/challenge/2026/en/index.php),
which asks for robust Segment Routing decisions under time-dependent network
interventions and reconfiguration budgets.

## Approach

The approach was intended to explore whether a direct MILP formulation could
solve T-ASR by:

* computing ECMP split coefficients for each time period;
* selecting segment paths with a sparse MILP containing flow, degree,
maximum-segment, reconfiguration-budget, and exact arc-utilization constraints;
* lexicographically minimizing the sorted arc–time utilization vector, using
persistent bottleneck peeling and a residual ordered-k-sum formulation when peeling fails; and
* supporting both HiGHS and Gurobi through a solver-independent model interface.

## Status

**I got to know the challenge far too late.** Consequently, this approach is
still in its early stages: it currently works only for smaller instances and
does not scale to larger instances yet. It should be considered an experimental
prototype, not a finished or competitive challenge solver.

The mathematical formulation and implementation details are documented in Documentation.pfd.

