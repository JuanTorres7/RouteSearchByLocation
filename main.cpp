#include "DataManager.cpp"

int main() {
    vector<pair<double, double>> points = {{40.7128, -74.0060}, {34.0522, -118.2437}, {51.5074, -0.1278}, {48.8566, 2.3522}};
    DataManager data(points);

    Graph graph = data.buildKNNGraph(3);
    graph.printGraph();

    vector<Query> queries = {
        {{{40.7128, -74.0060}, {34.0522, -118.2437}}, 1.0},
        {{{51.5074, -0.1278}, {48.8566, 2.3522}}, 1.5}
    };

    Route newRoute = {{{40.7128, -74.0060}, {34.0522, -118.2437}, {51.5074, -0.1278}}};

    vector<int> results = qloe(graph, queries, newRoute);
    cout << "Queries coincidentes: ";
    for (int queryId : results) {
        cout << queryId << " ";
    }
    cout << endl;

    return 0;
}