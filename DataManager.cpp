#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <chrono>
#include <unordered_set>

#include "DataReader.cpp"
using namespace std;

class Query {
public:
    vector<pair<double, double>> points;  // Vector que almacena los puntos
    Query(const vector<pair<double, double>>& inputPoints) : points(inputPoints) {}
    void print() const {
        for (const auto& point : points) {
            cout << "(" << point.first << ", " << point.second << ")" << endl;
        }
    }
};

class RSL {
public:
    Query query;
    double treshold;    

    RSL(const Query& inputQuery, double inputTreshold)
        : query(inputQuery), treshold(inputTreshold) {}
    
    void print() const {
        cout << "Treshold: " << treshold << endl;
        cout << "Query Points:" << endl;
        for (const auto& point : query.points) {
            cout << "(" << point.first << ", " << point.second << ")" << endl;
        }
    }
};


class Graph {
public:
    vector<pair<double, double>> vertices;

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

    const double R = 6378137.0;

    double degreesToRadians(double degrees) {
        return degrees * M_PI / 180.0;
    }

    pair<double, double> latLongToCartesian(double latitude, double longitude) {
        double x = R * degreesToRadians(longitude);
        double y = R * log(tan(M_PI / 4.0 + degreesToRadians(latitude) / 2.0));
        return {x, y};
    }

    double calculateDistance(const pair<double, double>& p1, const pair<double, double>& p2) {
        return sqrt(pow(p2.first - p1.first, 2) + pow(p2.second - p1.second, 2));
    }

    pair<pair<double, double>, double> findFarthestPoint(const pair<double, double>& startPoint, const vector<pair<double, double>>& allPoints) {
        double maxDistance = -1;
        pair<double, double> farthestPoint;

        for (const auto& point : allPoints) {
            double dist = calculateDistance(startPoint, point);
            if (dist > maxDistance) {
                maxDistance = dist;
                farthestPoint = point;
            }
        }
        return {farthestPoint, maxDistance};
    }

    double generateRandomAngle(double min_angle, double max_angle) {
        return min_angle + (rand() % 1000) / 1000.0 * (max_angle - min_angle);
    }

    pair<double, double> generateNewPoint(const pair<double, double>& currentPoint, double angle, double distance) {
        double radianAngle = angle * M_PI / 180.0;
        double newX = currentPoint.first + distance * cos(radianAngle);
        double newY = currentPoint.second + distance * sin(radianAngle);
        return {newX, newY};
    }

    double randomFactor(double min, double max) {
        return min + (rand() / (double)RAND_MAX) * (max - min);
    }


public:
    vector<TaxiTrip> trips;
    vector<pair<double, double>> pickupCartesian;  // Coordenadas cartesianas para pickup
    vector<pair<double, double>> dropoffCartesian; // Coordenadas cartesianas para dropoff
    vector<pair<double, double>> allCartesian;
    vector<int> tripTimes;                         // Tiempos de viaje
    vector<double> tripDistances;                  // Distancias de viaje

    DataManager(const std::string& file_name, int maxRows){
        DataReader data;
        vector<TaxiTrip>taxiTrips = data.readTaxiTrips(file_name, maxRows);
        for (auto trip : taxiTrips){
            trips.push_back(trip);
            tripTimes.push_back(trip.tripTimeSecs);
            tripDistances.push_back(trip.tripDistance);
            pickupCartesian.push_back(latLongToCartesian(trip.pickupLocation.first, trip.pickupLocation.second));
            dropoffCartesian.push_back(latLongToCartesian(trip.dropoffLocation.first, trip.dropoffLocation.second));

            allCartesian.push_back(latLongToCartesian(trip.pickupLocation.first, trip.pickupLocation.second));
            allCartesian.push_back(latLongToCartesian(trip.dropoffLocation.first, trip.dropoffLocation.second));
        }
    }

    void print_result(float tiempo, int cantQue){
        cout << "Para " << cantQue << " se demoró: " << tiempo << "segundos" << endl;
    }

    Query generateRandomQuery() {
        vector<pair<double, double>> randomPoints;
        unordered_set<string> visited;

        if (allCartesian.empty()) {
            cout << "No hay datos de coordenadas disponibles." << endl;
            return Query(randomPoints);
        }

        int numPoints = 3 + rand() % 8;  // Genera un número aleatorio entre 3 y 10

        while (randomPoints.size() < numPoints) {
            pair<double, double> selectedPoint = allCartesian[rand() % allCartesian.size()];
            string pointKey = to_string(selectedPoint.first) + "," + to_string(selectedPoint.second);
            if (visited.find(pointKey) == visited.end()) {
                randomPoints.push_back(selectedPoint);
                visited.insert(pointKey);  // Marcar el punto como visitado
            }
        }

        return Query(randomPoints);
    }

    Graph buildKNNGraph(const vector<pair<double, double>>& coordinates, int k) {
        Graph graph;

        for (const auto& coord : coordinates) {
            graph.addVertex(coord);
        }

        for (const auto& coord : coordinates) {
            priority_queue<pair<double, pair<double, double>>,
                        vector<pair<double, pair<double, double>>>,
                        greater<>> nearestNeighbors;

            for (const auto& other : coordinates) {
                if (coord != other) {
                    double distance = calculateDistance(coord, other);
                    nearestNeighbors.push({distance, other});
                }
            }
            int count = 0;
            while (!nearestNeighbors.empty() && count < k) {
                auto [distance, neighbor] = nearestNeighbors.top();
                nearestNeighbors.pop();

                double adjustedWeight = distance * randomFactor(0.8, 1.2);

                graph.addEdge(coord, neighbor, adjustedWeight);
                count++;
            }
        }

        return graph;
    }

    vector<Query> generateMultipleQueries(int n) {
        vector<Query> queries;

        for (int i = 0; i < n; ++i) {
            // Generar una query aleatoria
            Query query = generateRandomQuery();
            queries.push_back(query);  // Añadir la query al vector
        }

        return queries;
    }

    vector<Graph> generateMultipleSubgraphs(const Graph& grafo, int n) {
        vector<Graph> subgraphs;

        int numVertices = grafo.vertices.size();

        for (int i = 0; i < n; ++i) {
            int size = 2 + rand() % (numVertices - 1);  // Aleatorio entre 2 y numVertices

            Graph subgraph = generateRandomSubgraph(grafo, size);
            subgraphs.push_back(subgraph);  // Añadir el subgrafo al vector
        }

        return subgraphs;
    }


    // Función para mostrar todos los viajes
    void displayTrips() const {
        for (const auto& trip : trips) {
            trip.displayInfo();
            cout << endl;
        }
    }

    // Función para mostrar las coordenadas cartesianas
    void displayCartesianCoordinates() const {
        cout << "Pickup Cartesian Coordinates:" << endl;
        for (const auto& coord : pickupCartesian) {
            cout << "(" << coord.first << ", " << coord.second << ")" << endl;
        }
        cout << "Dropoff Cartesian Coordinates:" << endl;
        for (const auto& coord : dropoffCartesian) {
            cout << "(" << coord.first << ", " << coord.second << ")" << endl;
        }
        cout << endl;
    }

    Graph generateRandomSubgraph(const Graph& graph, int n) {
        Graph subgraph;

        if (graph.vertices.empty()) {
            cout << "Graph is empty!" << endl;
            return subgraph;
        }

        // Obtener un vértice inicial aleatorio
        auto it = graph.vertices.begin();
        advance(it, rand() % graph.vertices.size());
        pair<double, double> currentVertex = *it;

        subgraph.addVertex(currentVertex);

        // Construir la ruta y subgrafo
        for (int i = 1; i < n; ++i) {
            // Verificar si el vértice tiene vecinos
            if (graph.edges.find(currentVertex) == graph.edges.end() || graph.edges.at(currentVertex).empty()) {
                cout << "Vertex has no neighbors. Ending route early." << endl;
                break;
            }

            // Seleccionar un vecino aleatorio
            const auto& neighbors = graph.edges.at(currentVertex);
            int randomIndex = rand() % neighbors.size();
            const auto& [nextVertex, weight] = neighbors[randomIndex];

            // Añadir el siguiente vértice y arista al subgrafo
            subgraph.addVertex(nextVertex);
            subgraph.addEdge(currentVertex, nextVertex, weight);

            currentVertex = nextVertex;
        }

        return subgraph;
    }

    vector<double> calculateLowerCost(const Graph& subgraph, const Query& queryPoints) {
        vector<double> lowerCosts;

        if (subgraph.vertices.empty()) {
            cout << "Subgraph is empty!" << endl;
            return lowerCosts;
        }

        for (const auto& queryPoint : queryPoints.points) {
            double minCost = numeric_limits<double>::infinity();

            for (const auto& subgraphVertex : subgraph.vertices) {
                
                double distance = calculateDistance(queryPoint, subgraphVertex);

                if (distance < minCost) {
                    minCost = distance;
                }
            }

            lowerCosts.push_back(minCost);
        }

        return lowerCosts;
    }

    double get_similarity(const vector<double>& costos) {
        vector<double> costos_normalizados;
        for(auto c : costos){
            costos_normalizados.push_back(c/5000);
        }
        double sum = 0.0;

        for (const double& costo : costos_normalizados) {
            sum += exp(-costo); // Aplicar e^(-costo)
        }

        return sum;
    }

    vector<Query> DRM(vector<RSL> querys, Graph ruta){
        vector<Query> D;
        for(auto rsl: querys){
            double similitud = 0.0;
            vector<double> lower_cost = calculateLowerCost(ruta, rsl.query);
            similitud = get_similarity(lower_cost);
            if(similitud >= rsl.treshold){
                D.push_back(rsl.query);
            }
        }
        


        return D;
    }

    static string pairToString(const pair<double, double>& p) {
        return to_string(p.first) + "," + to_string(p.second);
    }


    static unordered_map<string, double> dijkstra(const Graph& graph, const pair<double, double>& start) {
        unordered_map<string, double> minCost;
        priority_queue<pair<double, pair<double, double>>, vector<pair<double, pair<double, double>>>, greater<>> pq;

        pq.push({0.0, start});
        minCost[pairToString(start)] = 0.0;

        while (pq.empty()) {
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

    vector<int> qloe(const Graph& graph, const vector<RSL>& rsl, const Graph& newRoute) {
        //auto inicio = std::chrono::high_resolution_clock::now();
        unordered_map<string, unordered_set<int>> queryMap;

        for (int i = 0; i < rsl.size()*0.001; ++i) {
            for (const auto& loc : rsl[i].query.points) {
                queryMap[pairToString(loc)].insert(i);
            }
        }

        unordered_map<string, double> routeCosts;
        
        for (const auto& vertex : newRoute.vertices) {
            auto costs = dijkstra(graph, vertex);
            
            for (const auto& [target, cost] : costs) {
                if (!routeCosts.count(target) || cost < routeCosts[target]) {
                    routeCosts[target] = cost;
                }
            }
        }

        vector<int> matchingQueries;

        for (const auto& [vertex, cost] : routeCosts) {
            if (queryMap.count(vertex)) {
                for (int queryId : queryMap[vertex]) {
                    const RSL& currentRSL = rsl[queryId];
                    
                    vector<double> lowerCost = calculateLowerCost(newRoute, currentRSL.query);
                    
                    double similarity = get_similarity(lowerCost);

                    if (similarity >= currentRSL.treshold) {
                        matchingQueries.push_back(queryId);
                    }
                }
            }
        }
        //auto fin = std::chrono::high_resolution_clock::now();
        //auto duracion = std::chrono::duration_cast<std::chrono::milliseconds>(fin - inicio);
        //std::cout << "El tiempo de ejecucion fue: " << duracion.count() << " milisegundos." << std::endl;

        return matchingQueries;
    }

    /*vector<Query> QLOEPlus(pair<double, double> v, vector<pair<double, double>> tau, vector<RSL>& T, unordered_map<pair<double, double>, set<Query>>& M) {
        vector<Query> resultQueries;

        // Lógica para calcular las nuevas queries basadas en el punto v, la ruta tau, y el conjunto T de RSLs
        for (const auto& rsl : T) {
            // Lógica para generar nuevas queries basadas en la consulta y el treshold de cada RSL
            vector<pair<double, double>> newPoints;
            
            // Generamos nuevas consultas a partir de los puntos de la consulta (query) y el treshold
            for (const auto& point : rsl.query.points) {
                double distance = sqrt(pow(v.first - point.first, 2) + pow(v.second - point.second, 2));  // Distancia Euclidiana
                if (distance <= rsl.treshold) {
                    newPoints.push_back(point);
                }
            }

            // Si hay puntos dentro del umbral, creamos una nueva consulta
            if (!newPoints.empty()) {
                Query newQuery(newPoints);
                resultQueries.push_back(newQuery);
            }
        }

        // Lógica para agregar las queries al mapa M
        for (const auto& query : resultQueries) {
            M[v].insert(query);  // Agrega cada consulta al mapa, asociada al punto v
        }

        return resultQueries;
    }*/

};
