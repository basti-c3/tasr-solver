#include "instance_data.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {
constexpr double INF = std::numeric_limits<double>::infinity();

std::string context(const std::string& file, const std::string& message) {
    return file + ": " + message;
}
} // namespace

InstanceData::InstanceData()
    : graph(),
      node_id(graph),
      node_name(graph),
      link_id(graph),
      metric_map(graph),
      capacity_map(graph),
      demand_graph(graph),
      demand_index_map(demand_graph),
      demand_values_map(demand_graph) {}

void InstanceData::load(const std::string& network_file,
                        const std::string& traffic_matrix_file,
                        const std::string& scenario_file) {
    loadNetwork(network_file);
    buildNetworkLookups();
    loadTrafficMatrix(traffic_matrix_file);
    loadScenario(scenario_file);
    validate();
}

int InstanceData::numNodes() const noexcept {
    return graph.nodeNum();
}

int InstanceData::numArcs() const noexcept {
    return graph.arcNum();
}

int InstanceData::numDemands() const noexcept {
    return static_cast<int>(demands.size());
}

InstanceData::Node InstanceData::nodeByExternalId(int external_id) const {
    if (external_id < 0 ||
        external_id >= static_cast<int>(nodes_by_external_id.size()) ||
        nodes_by_external_id[external_id] == nt::INVALID) {
        throw std::out_of_range("Unknown external node id " +
                                std::to_string(external_id));
    }
    return nodes_by_external_id[external_id];
}

InstanceData::Arc InstanceData::arcByExternalId(int external_id) const {
    if (external_id < 0 ||
        external_id >= static_cast<int>(arcs_by_external_id.size()) ||
        arcs_by_external_id[external_id] == nt::INVALID) {
        throw std::out_of_range("Unknown external link id " +
                                std::to_string(external_id));
    }
    return arcs_by_external_id[external_id];
}

int InstanceData::arcIndex(int source_id, int target_id) const {
    if (source_id < 0 || target_id < 0 ||
        source_id >= static_cast<int>(arc_index_matrix.size()) ||
        target_id >= static_cast<int>(arc_index_matrix.size())) {
        throw std::out_of_range("Node id outside arc-index matrix");
    }
    return arc_index_matrix[source_id][target_id];
}

double InstanceData::metric(int source_id, int target_id) const {
    return metric_matrix.at(source_id).at(target_id);
}

double InstanceData::capacity(int source_id, int target_id) const {
    return capacity_matrix.at(source_id).at(target_id);
}

double InstanceData::demandValue(int demand_index, int time_slot) const {
    return demand_matrix.at(demand_index).at(time_slot);
}

nt::JSONDocument InstanceData::readJsonDocument(const std::string& filename) {
    nt::ByteBuffer file;
    if (!file.load(filename.c_str())) {
        throw std::runtime_error("Cannot open JSON file: " + filename);
    }

    nt::JSONDocument document;
    document.Parse(file.data());

    if (document.HasParseError()) {
        throw std::runtime_error("Invalid JSON document: " + filename);
    }
    if (!document.IsObject()) {
        throw std::runtime_error("JSON root must be an object: " + filename);
    }
    return document;
}

void InstanceData::require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void InstanceData::loadNetwork(const std::string& filename) {
    nt::graphs::GraphProperties<Digraph> properties;
    properties.nodeMap("id", node_id);
    properties.nodeMap("name", node_name);
    properties.arcMap("id", link_id);
    properties.arcMap("metric", metric_map);
    properties.arcMap("capacity", capacity_map);

    nt::graphs::NetworkxIo<Digraph> reader("from", "to");
    if (!reader.readFile(filename.c_str(), graph, &properties)) {
        throw std::runtime_error("Could not read network file: " + filename);
    }

    require(graph.nodeNum() > 0, context(filename, "network has no nodes"));
    require(graph.arcNum() > 0, context(filename, "network has no links"));
}

void InstanceData::buildNetworkLookups() {
    int max_node_id = -1;
    for (Digraph::NodeIt node(graph); node != nt::INVALID; ++node) {
        require(node_id[node] >= 0, "Node IDs must be non-negative");
        max_node_id = std::max(max_node_id, node_id[node]);
    }

    nodes_by_external_id.assign(max_node_id + 1, nt::INVALID);
    for (Digraph::NodeIt node(graph); node != nt::INVALID; ++node) {
        const int id = node_id[node];
        require(nodes_by_external_id[id] == nt::INVALID,
                "Duplicate node id " + std::to_string(id));
        nodes_by_external_id[id] = node;
    }

    const int matrix_size = max_node_id + 1;
    arc_index_matrix.assign(
        matrix_size, std::vector<int>(matrix_size, -1));
    metric_matrix.assign(
        matrix_size, std::vector<double>(matrix_size, INF));
    capacity_matrix.assign(
        matrix_size, std::vector<double>(matrix_size, 0.0));
    outgoing_arc_indices.assign(matrix_size, {});
    incoming_arc_indices.assign(matrix_size, {});

    int max_link_id = -1;
    for (Digraph::ArcIt arc(graph); arc != nt::INVALID; ++arc) {
        require(link_id[arc] >= 0, "Link IDs must be non-negative");
        max_link_id = std::max(max_link_id, link_id[arc]);
    }
    arcs_by_external_id.assign(max_link_id + 1, nt::INVALID);

    arc_external_ids.reserve(graph.arcNum());
    arc_sources.reserve(graph.arcNum());
    arc_targets.reserve(graph.arcNum());
    arc_metrics.reserve(graph.arcNum());
    arc_capacities.reserve(graph.arcNum());

    for (Digraph::ArcIt arc(graph); arc != nt::INVALID; ++arc) {
        const int external_link_id = link_id[arc];
        const int source = node_id[graph.source(arc)];
        const int target = node_id[graph.target(arc)];
        const int compact_index = static_cast<int>(arc_external_ids.size());

        require(arcs_by_external_id[external_link_id] == nt::INVALID,
                "Duplicate link id " + std::to_string(external_link_id));
        require(arc_index_matrix[source][target] == -1,
                "Parallel links are not supported by the dense matrix "
                "representation");

        arcs_by_external_id[external_link_id] = arc;

        arc_external_ids.push_back(external_link_id);
        arc_sources.push_back(source);
        arc_targets.push_back(target);
        arc_metrics.push_back(metric_map[arc]);
        arc_capacities.push_back(capacity_map[arc]);

        arc_index_matrix[source][target] = compact_index;
        metric_matrix[source][target] = metric_map[arc];
        capacity_matrix[source][target] = capacity_map[arc];

        outgoing_arc_indices[source].push_back(compact_index);
        incoming_arc_indices[target].push_back(compact_index);
    }

    for (int i = 0; i < matrix_size; ++i) {
        metric_matrix[i][i] = 0.0;
    }
}

void InstanceData::loadTrafficMatrix(const std::string& filename) {
    nt::JSONDocument document = readJsonDocument(filename);

    require(document.HasMember("num_time_slots") &&
                document["num_time_slots"].IsInt(),
            context(filename, "missing integer num_time_slots"));
    require(document.HasMember("demands") &&
                document["demands"].IsArray(),
            context(filename, "missing demands array"));

    num_time_slots = document["num_time_slots"].GetInt();
    require(num_time_slots > 0,
            context(filename, "num_time_slots must be positive"));

    const auto& json_demands = document["demands"];
    demands.reserve(json_demands.Size());
    demand_matrix.reserve(json_demands.Size());

    demands_from_node.assign(nodes_by_external_id.size(), {});
    demands_to_node.assign(nodes_by_external_id.size(), {});

    for (std::size_t d = 0; d < json_demands.Size(); ++d) {
        const auto& item = json_demands[static_cast<rapidjson::SizeType>(d)];

        require(item.IsObject(),
                context(filename, "each demand must be an object"));
        require(item.HasMember("s") && item["s"].IsInt(),
                context(filename, "demand is missing integer s"));
        require(item.HasMember("t") && item["t"].IsInt(),
                context(filename, "demand is missing integer t"));
        require(item.HasMember("v") && item["v"].IsArray(),
                context(filename, "demand is missing value array v"));

        DemandRecord record;
        record.source_id = item["s"].GetInt();
        record.target_id = item["t"].GetInt();

        const Node source = nodeByExternalId(record.source_id);
        const Node target = nodeByExternalId(record.target_id);

        const auto& values = item["v"];
        require(static_cast<int>(values.Size()) == num_time_slots,
                context(filename,
                        "demand value vector has wrong number of time slots"));

        record.values.resize(num_time_slots);

        for (int time = 0; time < num_time_slots; ++time) {
            require(values[time].IsNumber(),
                context(filename, "demand value must be numeric"));

            const double value = values[time].GetDouble();

            require(value >= 0.0,
                context(filename, "demand values must be non-negative"));

            record.values[time] = value;
        }

        record.arc = demand_graph.addArc(source, target);
        demand_index_map[record.arc] = static_cast<int>(d);

        // Fill the Networktools array directly.
        // DoubleDynamicArray cannot be copy-assigned.
        auto& stored_values = demand_values_map[record.arc];
        stored_values.resize(num_time_slots);

        for (int time = 0; time < num_time_slots; ++time) {
            stored_values[time] = record.values[time];
        }

        demand_matrix.push_back(record.values);
        demands_from_node[record.source_id].push_back(static_cast<int>(d));
        demands_to_node[record.target_id].push_back(static_cast<int>(d));
        demands.push_back(std::move(record));
    }
}

void InstanceData::loadScenario(const std::string& filename) {
    nt::JSONDocument document = readJsonDocument(filename);

    require(document.HasMember("max_segments") &&
                document["max_segments"].IsInt(),
            context(filename, "missing integer max_segments"));
    require(document.HasMember("budget") &&
                document["budget"].IsArray(),
            context(filename, "missing budget array"));
    require(document.HasMember("interventions") &&
                document["interventions"].IsArray(),
            context(filename, "missing interventions array"));

    max_segments = document["max_segments"].GetInt();
    require(max_segments >= 0,
            context(filename, "max_segments must be non-negative"));

    budget_by_time.assign(num_time_slots, 0.0);
    intervention_links_by_time.assign(num_time_slots, {});
    link_intervened_at_time.assign(
        num_time_slots,
        std::vector<unsigned char>(arcs_by_external_id.size(), 0));

    for (const auto& item : document["budget"].GetArray()) {
        require(item.IsObject() &&
                    item.HasMember("t") && item["t"].IsInt() &&
                    item.HasMember("value") && item["value"].IsNumber(),
                context(filename, "invalid budget entry"));

        BudgetEntry entry;
        entry.time = item["t"].GetInt();
        entry.value = item["value"].GetDouble();

        require(entry.time >= 0 && entry.time < num_time_slots,
                context(filename, "budget time index out of range"));
        require(entry.value >= 0.0,
                context(filename, "budget must be non-negative"));

        budgets.push_back(entry);
        budget_by_time[entry.time] = entry.value;
    }

    for (const auto& item : document["interventions"].GetArray()) {
        require(item.IsObject() &&
                    item.HasMember("t") && item["t"].IsInt() &&
                    item.HasMember("links") && item["links"].IsArray(),
                context(filename, "invalid intervention entry"));

        InterventionEntry entry;
        entry.time = item["t"].GetInt();
        require(entry.time >= 0 && entry.time < num_time_slots,
                context(filename, "intervention time index out of range"));

        for (const auto& link_value : item["links"].GetArray()) {
            require(link_value.IsInt(),
                    context(filename, "intervention link id must be integer"));

            const int external_link_id = link_value.GetInt();
            arcByExternalId(external_link_id); // validates existence

            entry.external_link_ids.push_back(external_link_id);
            intervention_links_by_time[entry.time].push_back(external_link_id);
            link_intervened_at_time[entry.time][external_link_id] = 1;
        }

        interventions.push_back(std::move(entry));
    }
}

void InstanceData::validate() const {
    require(static_cast<int>(demands.size()) == demand_graph.arcNum(),
            "Demand vector and DemandGraph size differ");
    require(static_cast<int>(demand_matrix.size()) == numDemands(),
            "Demand matrix has wrong number of rows");
    require(static_cast<int>(budget_by_time.size()) == num_time_slots,
            "Budget vector has wrong size");
    require(static_cast<int>(intervention_links_by_time.size()) ==
                num_time_slots,
            "Intervention vector has wrong size");

    for (int id = 0; id < static_cast<int>(nodes_by_external_id.size()); ++id) {
        require(nodes_by_external_id[id] != nt::INVALID,
                "Node IDs must be contiguous; missing id " +
                    std::to_string(id));
    }
}
