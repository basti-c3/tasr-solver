#include "split_coefficients.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>

namespace {
constexpr double INF = std::numeric_limits<double>::infinity();
using QueueEntry = std::pair<double, int>;
using MinQueue = std::priority_queue<
    QueueEntry,
    std::vector<QueueEntry>,
    std::greater<QueueEntry>>;
}

SplitCoefficients::SplitCoefficients(const InstanceData& instance)
    : instance_(instance),
      time_slots_(instance.num_time_slots),
      nodes_(instance.numNodes()),
      arcs_(instance.numArcs()),
      data_(
          static_cast<std::size_t>(time_slots_),
          SegmentTable(
              static_cast<std::size_t>(nodes_) *
              static_cast<std::size_t>(nodes_))) {
    if (time_slots_ <= 0) {
        throw std::invalid_argument("At least one time slot is required");
    }
    if (nodes_ <= 0) {
        throw std::invalid_argument("At least one node is required");
    }
}

void SplitCoefficients::compute() {
    for (auto& time_table : data_) {
        for (auto& segment : time_table) {
            segment.clear();
        }
    }

    for (int time = 0; time < time_slots_; ++time) {
        computeTimeSlot(time);
    }
}

const SplitCoefficients::SparseArcVector& SplitCoefficients::nonzero(
    int time,
    int source,
    int target) const {
    checkTime(time);
    checkNode(source);
    checkNode(target);
    return data_[time][segmentIndex(source, target)];
}

double SplitCoefficients::value(
    int source,
    int target,
    int arc_index,
    int time) const {
    checkArc(arc_index);
    const auto& entries = nonzero(time, source, target);

    const auto it = std::lower_bound(
        entries.begin(),
        entries.end(),
        arc_index,
        [](const ArcCoefficient& entry, int requested_arc) {
            return entry.arc_index < requested_arc;
        });

    if (it == entries.end() || it->arc_index != arc_index) {
        return 0.0;
    }
    return it->coefficient;
}

int SplitCoefficients::numTimeSlots() const noexcept { return time_slots_; }
int SplitCoefficients::numNodes() const noexcept { return nodes_; }
int SplitCoefficients::numArcs() const noexcept { return arcs_; }

std::size_t SplitCoefficients::numberOfNonzeroCoefficients() const noexcept {
    std::size_t count = 0;
    for (const auto& time_table : data_) {
        for (const auto& segment : time_table) {
            count += segment.size();
        }
    }
    return count;
}

std::size_t SplitCoefficients::segmentIndex(int source, int target) const {
    return static_cast<std::size_t>(source) *
               static_cast<std::size_t>(nodes_) +
           static_cast<std::size_t>(target);
}

void SplitCoefficients::checkTime(int time) const {
    if (time < 0 || time >= time_slots_) {
        throw std::out_of_range("Invalid time index " + std::to_string(time));
    }
}

void SplitCoefficients::checkNode(int node) const {
    if (node < 0 || node >= nodes_) {
        throw std::out_of_range("Invalid node index " + std::to_string(node));
    }
}

void SplitCoefficients::checkArc(int arc_index) const {
    if (arc_index < 0 || arc_index >= arcs_) {
        throw std::out_of_range(
            "Invalid compact arc index " + std::to_string(arc_index));
    }
}

void SplitCoefficients::computeTimeSlot(int time) {
    const std::vector<double> time_metrics = buildTimeMetrics(time);

    for (int target = 0; target < nodes_; ++target) {
        const auto distances = reverseDijkstra(target, time_metrics);
        computeForTarget(time, target, time_metrics, distances);
    }
}

std::vector<double> SplitCoefficients::buildTimeMetrics(int time) const {
    std::vector<double> metrics = instance_.arc_metrics;

    for (const int external_link_id :
         instance_.intervention_links_by_time[time]) {
        const InstanceData::Arc arc =
            instance_.arcByExternalId(external_link_id);

        const int source =
            instance_.node_id[instance_.graph.source(arc)];
        const int target =
            instance_.node_id[instance_.graph.target(arc)];

        const int compact_arc_index =
            instance_.arc_index_matrix[source][target];

        if (compact_arc_index < 0) {
            throw std::runtime_error(
                "Could not map external link id " +
                std::to_string(external_link_id) +
                " to compact arc index");
        }

        // Interventions are represented by an infinite metric.
        metrics[compact_arc_index] = INF;
    }

    return metrics;
}

std::vector<double> SplitCoefficients::reverseDijkstra(
    int target,
    const std::vector<double>& time_metrics) const {
    std::vector<double> distance(nodes_, INF);
    MinQueue queue;

    distance[target] = 0.0;
    queue.emplace(0.0, target);

    while (!queue.empty()) {
        const auto [current_distance, current_node] = queue.top();
        queue.pop();

        if (current_distance > distance[current_node]) {
            continue;
        }

        for (const int arc_index :
             instance_.incoming_arc_indices[current_node]) {
            const double metric = time_metrics[arc_index];
            if (!std::isfinite(metric)) {
                continue;
            }

            const int predecessor = instance_.arc_sources[arc_index];
            const double candidate = current_distance + metric;

            if (candidate < distance[predecessor]) {
                distance[predecessor] = candidate;
                queue.emplace(candidate, predecessor);
            }
        }
    }

    return distance;
}

void SplitCoefficients::computeForTarget(
    int time,
    int target,
    const std::vector<double>& time_metrics,
    const std::vector<double>& distance_to_target) {
    std::vector<std::vector<int>> forwarding_out(nodes_);

    for (int arc_index = 0; arc_index < arcs_; ++arc_index) {
        const double metric = time_metrics[arc_index];
        if (!std::isfinite(metric)) {
            continue;
        }

        const int source = instance_.arc_sources[arc_index];
        const int next = instance_.arc_targets[arc_index];

        if (!std::isfinite(distance_to_target[source]) ||
            !std::isfinite(distance_to_target[next])) {
            continue;
        }

        if (approximatelyEqual(
                distance_to_target[source],
                metric + distance_to_target[next])) {
            forwarding_out[source].push_back(arc_index);
        }
    }

    std::vector<int> node_order(nodes_);
    for (int node = 0; node < nodes_; ++node) {
        node_order[node] = node;
    }

    std::sort(
        node_order.begin(),
        node_order.end(),
        [&](int lhs, int rhs) {
            return distance_to_target[lhs] > distance_to_target[rhs];
        });

    std::vector<double> node_flow(nodes_, 0.0);
    std::vector<double> arc_flow(arcs_, 0.0);
    std::vector<int> touched_arcs;
    touched_arcs.reserve(arcs_);

    for (int source = 0; source < nodes_; ++source) {
        if (source == target ||
            !std::isfinite(distance_to_target[source])) {
            continue;
        }

        std::fill(node_flow.begin(), node_flow.end(), 0.0);
        for (const int arc_index : touched_arcs) {
            arc_flow[arc_index] = 0.0;
        }
        touched_arcs.clear();

        node_flow[source] = 1.0;

        for (const int node : node_order) {
            const double incoming = node_flow[node];
            if (incoming <= absolute_tolerance_ || node == target) {
                continue;
            }

            const auto& eligible_arcs = forwarding_out[node];
            if (eligible_arcs.empty()) {
                throw std::runtime_error(
                    "Reachable node " + std::to_string(node) +
                    " has no forwarding arc toward target " +
                    std::to_string(target));
            }

            const double share =
                incoming / static_cast<double>(eligible_arcs.size());

            for (const int arc_index : eligible_arcs) {
                if (arc_flow[arc_index] == 0.0) {
                    touched_arcs.push_back(arc_index);
                }

                arc_flow[arc_index] += share;
                node_flow[instance_.arc_targets[arc_index]] += share;
            }
        }

        auto& result = data_[time][segmentIndex(source, target)];
        result.reserve(touched_arcs.size());

        for (const int arc_index : touched_arcs) {
            const double coefficient = arc_flow[arc_index];
            if (coefficient > absolute_tolerance_) {
                result.push_back({arc_index, coefficient});
            }
        }

        std::sort(
            result.begin(),
            result.end(),
            [](const ArcCoefficient& lhs, const ArcCoefficient& rhs) {
                return lhs.arc_index < rhs.arc_index;
            });
    }
}

bool SplitCoefficients::approximatelyEqual(double lhs, double rhs) {
    const double scale =
        std::max({1.0, std::abs(lhs), std::abs(rhs)});

    return std::abs(lhs - rhs) <=
           absolute_tolerance_ + relative_tolerance_ * scale;
}
