#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <unordered_set>

#include "DataReader.cpp"
using namespace std;

static string pairToString(const pair<double, double>& p) {
    return to_string(p.first) + "," + to_string(p.second);
}

struct Query {
    vector<pair<double, double>> locations;
    double threshold;                      
};

struct Route {
    vector<pair<double, double>> path;
};

struct CRSLTuple {
    int queryId;
    int visitedLocations;
    double similarityUpperBound;
};

class Graph {
public:
    vector<pair<double, double>> vertices;

    // Aristas (vecino y peso)
    map<pair<double, double>, vector<pair<pair<double, double>, double>>> edges;

    void addVertex(const pair<double, double>& vertex) {
        vertices.push_back(vertex);
        edges[vertex] = {};
    }

    void addEdge(const pair<double, double>& v1, const pair<double, double>& v2, double weight) {
        edges[v1].push_back({v2, weight});
    }

    void printGraph() const {
        cout << "Graph Representation:" << endl;
        for (const auto& [vertex, neighbors] : edges) {
            cout << "(" << vertex.first << ", " << vertex.second << ") -> ";
            for (const auto& [neighbor, weight] : neighbors) {
                cout << "[" << "(" << neighbor.first << ", " << neighbor.second << "), " 
                    << "weight: " << weight << "] ";
            }
            cout << endl;
        }
    }
};

class DataManager {
private:
    // Radio de la tierra en metros
    const double R = 6378137.0;

    double degreesToRadians(double degrees) {
        return degrees * M_PI / 180.0;
    }

    // latitud/longitud a cartesianas usando Web Mercator
    pair<double, double> latLongToCartesian(double latitude, double longitude) {
        double x = R * degreesToRadians(longitude);
        double y = R * log(tan(M_PI / 4.0 + degreesToRadians(latitude) / 2.0));
        return {x, y};
    }

    double calculateDistance(const pair<double, double>& p1, const pair<double, double>& p2) {
        return sqrt(pow(p2.first - p1.first, 2) + pow(p2.second - p1.second, 2));
    }

public:
    vector<pair<double, double>> allCartesian;

    DataManager(const vector<pair<double, double>>& points) {
        for (const auto& point : points) {
            allCartesian.push_back(latLongToCartesian(point.first, point.second));
        }
    }

    Graph buildKNNGraph(int k) {
        Graph graph;

        for (const auto& coord : allCartesian) {
            graph.addVertex(coord);
        }

        for (const auto& coord : allCartesian) {
            priority_queue<pair<double, pair<double, double>>, 
                           vector<pair<double, pair<double, double>>>, 
                           greater<>> nearestNeighbors;

            for (const auto& other : allCartesian) {
                if (coord != other) {
                    double distance = calculateDistance(coord, other);
                    nearestNeighbors.push({distance, other});
                }
            }

            int count = 0;
            while (!nearestNeighbors.empty() && count < k) {
                auto [distance, neighbor] = nearestNeighbors.top();
                nearestNeighbors.pop();
                graph.addEdge(coord, neighbor, distance);
                count++;
            }
        }

        return graph;
    }

};

static unordered_map<string, double> dijkstra(const Graph& graph, const pair<double, double>& start) {
    unordered_map<string, double> minCost;
    priority_queue<pair<double, pair<double, double>>, vector<pair<double, pair<double, double>>>, greater<>> pq;

    pq.push({0.0, start});
    minCost[pairToString(start)] = 0.0;

    while (!pq.empty()) {
        auto [currentCost, currentVertex] = pq.top();
        pq.pop();

        if (currentCost > minCost[pairToString(currentVertex)]) {
            continue;
        }

        if (graph.edges.find(currentVertex) == graph.edges.end()) {
            continue; // Si el vértice no está en el grafo, pasamos al siguiente
        }

        for (const auto& [neighbor, weight] : graph.edges.at(currentVertex)) {
            double newCost = currentCost + weight;
            string neighborKey = pairToString(neighbor);
            if (!minCost.count(neighborKey) || newCost < minCost[neighborKey]) {
                minCost[neighborKey] = newCost;
                pq.push({newCost, neighbor});
            }
        }
    }

    return minCost;
}

static vector<int> qloe(const Graph& graph, const vector<Query>& queries, const Route& newRoute) {
    unordered_map<string, unordered_set<int>> queryMap;
    vector<int> matchingQueries;

    for (int i = 0; i < queries.size(); ++i) {
        for (const auto& loc : queries[i].locations) {
            queryMap[pairToString(loc)].insert(i);
        }
    }

    unordered_map<int, CRSLTuple> crslTuples;
    for (int i = 0; i < queries.size(); ++i) {
        crslTuples[i] = {i, 0, static_cast<double>(queries[i].locations.size())};
    }

    unordered_map<string, double> routeCosts;
    for (const auto& vertex : newRoute.path) {
        auto costs = dijkstra(graph, vertex);
        for (const auto& [target, cost] : costs) {
            if (!routeCosts.count(target) || cost < routeCosts[target]) {
                routeCosts[target] = cost;
            }
        }
    }

    for (const auto& [vertex, cost] : routeCosts) {
        if (queryMap.count(vertex)) {
            for (int queryId : queryMap[vertex]) {
                auto& tuple = crslTuples[queryId];
                tuple.visitedLocations++;
                tuple.similarityUpperBound = tuple.visitedLocations + exp(-cost);

                if (tuple.similarityUpperBound >= queries[queryId].threshold) {
                    matchingQueries.push_back(queryId);
                }
            }
        }
    }

    return matchingQueries;
}

static vector<int> qloePlus(const Graph& graph, const vector<Query>& queries, const Route& newRoute) {
    unordered_map<string, unordered_set<int>> queryMap;
    vector<int> matchingQueries;
    
    return matchingQueries;
}