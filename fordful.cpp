#include <iostream>
#include <queue>
#include <string.h>
#include <fstream>
#include <mpi.h>

using namespace std;

#define V 1632405 

bool bfs(int rGraph[V][V], int s, int t, int parent[], int rank)
{
    bool visited[V];
    memset(visited, 0, sizeof(visited));

    queue<int> q;
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
                    return true;
                }
                if (rank == 0) {
                    q.push(v);
                    parent[v] = u;
                    visited[v] = true;
                } else {
                    MPI_Recv(&rGraph[v][u], 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
    }

    return false;
}

int fordFulkerson(int graph[V][V], int s, int t, int rank)
{

        int u, v;
    int rGraph[V][V];
    int parent[V];

    for (u = 0; u < V; u++) {
        for (v = 0; v < V; v++) {
            rGraph[u][v] = graph[u][v];
        }
    }

    memset(parent, 0, sizeof(parent));

    int max_flow = 0;

    while (bfs(rGraph, s, t, parent, rank)) {
        int path_flow = 1;

        for (v = t; v != s; v = parent[v]) {
            u = parent[v];
            if (rank == 0) {
                rGraph[u][v] -= path_flow;
                rGraph[v][u] += path_flow;
            } else {
                MPI_Recv(&path_flow, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        max_flow += path_flow;
    }

    return max_flow;
}


int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

    int graph[V][V] = {0};
    int from, to;

    while (file >> from >> to) {
        graph[from][to] = 1;
    }

    if (rank == 0) {
        int max_flow = fordFulkerson(graph, source, sink, rank);
        cout << "The maximum possible flow is " << max_flow << endl;
    } else {
        fordFulkerson(graph, source, sink, rank);
    }

    MPI_Finalize();

    return 0;
}