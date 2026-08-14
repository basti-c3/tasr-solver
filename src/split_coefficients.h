#pragma once

#include "instance_data.h"

#include <cstddef>
#include <vector>

class SplitCoefficients {
public:
    struct ArcCoefficient {
        int arc_index = -1;
        double coefficient = 0.0;
    };

    using SparseArcVector = std::vector<ArcCoefficient>;

    explicit SplitCoefficients(const InstanceData& instance);

    void compute();

    const SparseArcVector& nonzero(int time, int source, int target) const;

    double value(int source, int target, int arc_index, int time) const;

    int numTimeSlots() const noexcept;
    int numNodes() const noexcept;
    int numArcs() const noexcept;
    std::size_t numberOfNonzeroCoefficients() const noexcept;

private:
    using SegmentTable = std::vector<SparseArcVector>;

    const InstanceData& instance_;
    int time_slots_;
    int nodes_;
    int arcs_;

    // data_[time][source * nodes_ + target]
    std::vector<SegmentTable> data_;

    static constexpr double relative_tolerance_ = 1e-10;
    static constexpr double absolute_tolerance_ = 1e-12;

    std::size_t segmentIndex(int source, int target) const;
    void checkTime(int time) const;
    void checkNode(int node) const;
    void checkArc(int arc_index) const;

    void computeTimeSlot(int time);
    std::vector<double> buildTimeMetrics(int time) const;
    std::vector<double> reverseDijkstra(
        int target,
        const std::vector<double>& time_metrics) const;
    void computeForTarget(
        int time,
        int target,
        const std::vector<double>& time_metrics,
        const std::vector<double>& distance_to_target);

    static bool approximatelyEqual(double lhs, double rhs);
};
