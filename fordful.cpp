#include <iostream>
#include <queue>
#include <cstring> 
#include <fstream>
#include <mpi.h>

using namespace std;

#define V 1632405

bool bfs(int** rGraph, int s, int t, int* parent, int rank)
{
    bool* visited = new bool[V];
    memset(visited, 0, sizeof(bool) * V);

    queue<int> q;
    int path_flow=0;
    if (rank == 0) {
        q.push(s);
        visited[s] = true;
        parent[s] = -1;
    }

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (int v = 0; v < V; v++) {
            if (visited[v] == false && rGraph[u][v] > 0) {
                if (v == t) {
                    parent[v] = u;
                    delete[] visited; 
                    return true;
                }
                if (rank == 0) {
                    q.push(v);
                    parent[v] = u;
                    visited[v] = true;
                    path_flow=1;
                    cout << "MPI Send " << rank << endl;
                    MPI_Send(visited, V, MPI_C_BOOL, 1, 0, MPI_COMM_WORLD);
                    cout << "MPI Send " << rank << endl;
                    MPI_Send(&path_flow, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
                }
                else {
                    cout << "MPI " << rank << endl;
                    MPI_Recv(visited, V, MPI_C_BOOL, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    cout << "MPI" << rank << endl;
                    MPI_Recv(&path_flow, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
    }

    delete[] visited; 
    return false;
}

int fordFulkerson(int** graph, int s, int t, int rank)
{
    int u, v;
    int** rGraph = new int*[V];
    for (u = 0; u < V; u++) {
        rGraph[u] = new int[V];
        for (v = 0; v < V; v++) {
            rGraph[u][v] = graph[u][v];
        }
    }

    int* parent = new int[V];
    memset(parent, 0, sizeof(int) * V);

    int max_flow = 0;

    while (bfs(rGraph, s, t, parent, rank)) {
        int path_flow = 1;

        for (v = t; v != s; v = parent[v]) {
            u = parent[v];
            if (rank == 0) {
                rGraph[u][v] -= path_flow;
                rGraph[v][u] += path_flow;
                MPI_Send(&path_flow, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
            }
            else {
                MPI_Recv(&path_flow, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        max_flow += path_flow;
    }

    for (u = 0; u < V; u++) {
        delete[] rGraph[u];
    }
    delete[] rGraph;
    delete[] parent;

    return max_flow;
}

int main(int argc, char* argv[])
{   printf("#Hi hello\n");
     printf("#Debug 2\n");

    MPI_Init(NULL,NULL);
     printf("#DEBUG 1.5\n");

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("DEBUG 2 After RANK AND SIZE\n");

    if (argc != 2) {
        if (rank == 0) {
            cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        }
        MPI_Finalize();
        return 1;
    }

    string inputFile = argv[1];
    ifstream file(inputFile);

    int source = 0;
    int sink = 1632404;

    int** graph = new int*[V];
    for (int i = 0; i < V; i++) {
        graph[i] = new int[V];
        memset(graph[i], 0, sizeof(int) * V);
    }

    int from, to;

    while (file >> from >> to) {
        graph[from][to] = 1;
    }

    if (rank == 0) {
        int max_flow = fordFulkerson(graph, source, sink, rank);
        cout << "The maximum possible flow is " << max_flow << endl;
    }
    else {
        fordFulkerson(graph, source, sink, rank);
    }

    for (int i = 0; i < V; i++) {
        delete[] graph[i];
    }
    delete[] graph;

    MPI_Finalize();

    return 0;
}
