#include "DataManager.cpp"


int main(){
    // 1 .Crear querys 3 <= q <= 10 -> QueryGenerate.cpp
    // Crear rutas random 2 < t < n
    // costo
    // probabilidad

    // 1. Crear Querys ---------------------------
    // Obtencion y lectura de querys
    
    int vtx = 20000;
    DataManager data("cleaned_trips_data.csv", 5000);
    Graph grafo;
    //data.displayTrips();
    //data.displayCartesianCoordinates();
    grafo = data.buildKNNGraph(data.allCartesian, 4);
    //grafo.printGraph();

    vector<Query> querys = data.generateMultipleQueries(1000);
    vector<Graph> rutas = data.generateMultipleSubgraphs(grafo, 400);
    vector<RSL> rsls;
    for(auto query : querys){
        rsls.push_back(RSL(query, 0.8));
    }
    cout << "Cantidad de vertices en grafo: " << vtx << endl;
    cout << "Tiempo de ejecucion de DRM" << endl;
    auto inicio = chrono::high_resolution_clock::now();
    for(auto ruta : rutas){
        vector<Query> drm = data.DRM(rsls, ruta);
    }
    auto fin = chrono::high_resolution_clock::now();
    auto duracion = std::chrono::duration_cast<chrono::milliseconds>(fin - inicio);
    cout << "El tiempo de ejecucion fue: " << duracion.count() << " milisegundos." << endl;
    
    cout << "Tiempo de ejecucion de QLOE" << endl;
    auto inicio2 = std::chrono::high_resolution_clock::now();
    for(auto ruta : rutas){
        vector<int> qloe = data.qloe(grafo, rsls, ruta);
    }
    auto fin2 = chrono::high_resolution_clock::now();
    auto duracion2 = chrono::duration_cast<chrono::milliseconds>(fin2 - inicio2);
    cout << "El tiempo de ejecucion fue: " << duracion2.count() << " milisegundos." << endl;

    
    return 0;
}
