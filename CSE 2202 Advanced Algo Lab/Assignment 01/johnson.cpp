#include <bits/stdc++.h>
using namespace std;

const int INF = 1e9 + 7; // Define a constant for infinity
int n = 100;
vector<vector<pair<int, int>>> Graph(n);
vector<int> globalDistance(n, INF); // Global distance array

void bellman_ford(int source, int n) {
    globalDistance[source] = 0;

    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (auto it : Graph[j]) {
                int adjnode = it.first;
                int weight = it.second;
                if (globalDistance[adjnode] > globalDistance[j] + weight) {
                    globalDistance[adjnode] = globalDistance[j] + weight;
                }
            }
        }
    }

    // Check for negative cycles
    for (int j = 0; j < n; ++j) {
        for (auto it : Graph[j]) {
            int adjnode = it.first;
            int weight = it.second;
            if (globalDistance[adjnode] > globalDistance[j] + weight) {
                cout << "Negative cycle detected!" << endl;
                return;
            }
        }
    }
}

void dijkstra(int source, int n) {
    set<pair<int, int>> st;
    vector<int> localDistance(n, INF); // Local distance array for Dijkstra
    localDistance[source] = 0;
    st.insert({0, source});

    while (!st.empty()) {
        auto it = *(st.begin());
        int dis = it.first;
        int node = it.second;
        st.erase(it);

        for (auto neighbor : Graph[node]) {
            int adjnode = neighbor.first;
            int weight = neighbor.second;

            if (localDistance[adjnode] > dis + weight) {
                localDistance[adjnode] = dis + weight;
                st.insert({localDistance[adjnode], adjnode});
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        if (localDistance[i] == INF)
            cout << "Distance from " << source << " to " << i << ": INF" << endl;
        else
            cout << "Distance from " << source << " to " << i << ": " << localDistance[i] << endl;
    }
}

int main() {
    int e, u, v, w;
    cin >> n >> e;

    // Reset graph and distances
    Graph.assign(n, vector<pair<int, int>>());
    globalDistance.assign(n, INF);

    // Read the graph edges
    for (int i = 0; i < e; ++i) {
        cin >> u >> v >> w;
        Graph[u].push_back({v, w});
    }

    // Add a new node and connect to all others with weight 0
    int new_node = n;
    Graph.resize(n + 1); // Increase graph size to account for the new node
    for (int i = 0; i < n; i++) {
        Graph[new_node].push_back({i, 0});
    }

    // Run Bellman-Ford from the new node
    bellman_ford(new_node, n + 1);

    // Reweight the graph
    for (int i = 0; i < n; i++) {
        for (auto& it : Graph[i]) {
            it.second = it.second + globalDistance[i] - globalDistance[it.first];
        }
    }

    // Run Dijkstra from each node
    for (int i = 0; i < n; i++) {
        dijkstra(i, n);
        cout << endl;
    }

    return 0;
}


