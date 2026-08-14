#pragma once

#include "networktools.h"

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

class InstanceData {
public:
    using Digraph = nt::graphs::SmartDigraph;
    using Node = Digraph::Node;
    using Arc = Digraph::Arc;
    using DemandGraph = nt::graphs::DemandGraph<Digraph>;
    using DemandArc = DemandGraph::Arc;

    struct BudgetEntry {
        int time = -1;
        double value = 0.0;
    };

    struct InterventionEntry {
        int time = -1;
        std::vector<int> external_link_ids;
    };

    struct DemandRecord {
        int source_id = -1;
        int target_id = -1;
        std::vector<double> values;
        DemandArc arc = nt::INVALID;
    };

    InstanceData();

    // Reads all three files and builds every derived representation.
    void load(const std::string& network_file,
              const std::string& traffic_matrix_file,
              const std::string& scenario_file);

    // Basic sizes.
    int numNodes() const noexcept;
    int numArcs() const noexcept;
    int numDemands() const noexcept;

    // External JSON ID -> Networktools object.
    Node nodeByExternalId(int external_id) const;
    Arc arcByExternalId(int external_id) const;

    // Direct dense lookups. A missing arc is represented by -1 / infinity / 0.
    int arcIndex(int source_id, int target_id) const;
    double metric(int source_id, int target_id) const;
    double capacity(int source_id, int target_id) const;
    double demandValue(int demand_index, int time_slot) const;

    // Networktools graph and property maps.
    Digraph graph;
    Digraph::NodeMap<int> node_id;
    Digraph::NodeMap<nt::String> node_name;
    Digraph::ArcMap<int> link_id;
    Digraph::ArcMap<double> metric_map;
    Digraph::ArcMap<double> capacity_map;

    // Traffic demands as a Networktools DemandGraph.
    DemandGraph demand_graph;
    DemandGraph::ArcMap<int> demand_index_map;
    DemandGraph::ArcMap<nt::DoubleDynamicArray> demand_values_map;

    // Dense representations useful when constructing an optimization model.
    std::vector<Node> nodes_by_external_id;
    std::vector<Arc> arcs_by_external_id;

    std::vector<int> arc_external_ids;
    std::vector<int> arc_sources;
    std::vector<int> arc_targets;
    std::vector<double> arc_metrics;
    std::vector<double> arc_capacities;

    std::vector<std::vector<int>> arc_index_matrix;
    std::vector<std::vector<double>> metric_matrix;
    std::vector<std::vector<double>> capacity_matrix;

    std::vector<std::vector<int>> outgoing_arc_indices;
    std::vector<std::vector<int>> incoming_arc_indices;

    int num_time_slots = 0;
    std::vector<DemandRecord> demands;
    std::vector<std::vector<double>> demand_matrix; // [demand][time]
    std::vector<std::vector<int>> demands_from_node;
    std::vector<std::vector<int>> demands_to_node;

    int max_segments = 0;
    std::vector<BudgetEntry> budgets;
    std::vector<double> budget_by_time;

    std::vector<InterventionEntry> interventions;
    std::vector<std::vector<int>> intervention_links_by_time;
    std::vector<std::vector<unsigned char>> link_intervened_at_time;

private:
    void loadNetwork(const std::string& filename);
    void loadTrafficMatrix(const std::string& filename);
    void loadScenario(const std::string& filename);

    void buildNetworkLookups();
    void validate() const;

    static nt::JSONDocument readJsonDocument(const std::string& filename);
    static void require(bool condition, const std::string& message);
};
