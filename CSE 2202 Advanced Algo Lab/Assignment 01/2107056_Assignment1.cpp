#include <iostream>
#include <vector>
#include <unordered_map>
#include <limits>
#include <utility>
#include <algorithm>

using namespace std;

class MinHeap {
private:
    vector<pair<int, int>> heap;       // Vector to store the heap elements as pairs (key, value)
    unordered_map<int, int> key_index; // Map to store the index of each key in the heap

    // Function to swap two elements in the heap
    void swap(int i, int j) {
      std::  swap(key_index[heap[i].first], key_index[heap[j].first]);; // Update indices in the map
       std::swap(heap[i], heap[j]);
    }

    // Function to maintain the heap property by moving the element up
    void heapify_up(int index) {
        while (index > 0) {
            int parent = (index - 1) / 2;
            if (heap[index].second < heap[parent].second) {
                swap(index, parent);
                index = parent;
            } else {
                break;
            }
        }
    }

    // Function to maintain the heap property by moving the element down
    void heapify_down(int index) {
        int left = 2 * index + 1;
        int right = 2 * index + 2;
        int smallest = index;             // Assume the smallest element is at the current index

        // Check if the left child exists and is smaller than the current element
        if (left < heap.size() && heap[left].second < heap[smallest].second) {
            smallest = left;
        }
        // Check if the right child exists and is smaller than the current smallest element
        if (right < heap.size() && heap[right].second < heap[smallest].second) {
            smallest = right;
        }
        // If the smallest element is not at the current index, swap and continue heapifying down
        if (smallest != index) {
            swap(index, smallest);
            heapify_down(smallest);
        }
    }

public:
    // Function to insert a new key-value pair into the heap
    void insert(int key, int value) {
        heap.push_back({key, value});     // Add the new element to the end of the heap
        int index = heap.size() - 1;
        key_index[key] = index;
        heapify_up(index);
    }

    // Function to remove and return the key-value pair with the smallest key
    pair<int, int> extract_min() {
        if (heap.empty()) {
            return {-1, -1};
        }
        auto min = heap.front();          // Get the minimum element from the front of the heap
        swap(0, heap.size() - 1);
        heap.pop_back();
        key_index.erase(min.first);
        if (!heap.empty()) {
            heapify_down(0);
        }
        return min;
    }

    // Function to decrease the value associated with a given key
    void decrease_key(int key, int new_value) {
        int index = key_index[key];
        heap[index].second = new_value;
        heapify_up(index);
    }

    // Function to check if the heap is empty
    bool is_empty() const {
        return heap.empty();
    }
};

// Function to implement Dijkstra's algorithm using the Min-Heap
void dijkstra(const vector<vector<pair<int, int>>>& graph, int source) {
    int n = graph.size();
    vector<int> dist(n, numeric_limits<int>::max()); // Initialize distances from the source
    MinHeap min_heap;                          // Create a Min-Heap

    dist[source] = 0;
    min_heap.insert(source, 0);                // Insert the source vertex into the heap

    while (!min_heap.is_empty()) {
        auto [u, d] = min_heap.extract_min();

        if (d > dist[u]) {
            continue;
        }

        for (const auto& [v, weight] : graph[u]) { // Iterate through the neighboring vertices
            if (dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
                min_heap.insert(v, dist[v]);
                min_heap.decrease_key(v, dist[v]);
            }
        }
    }

    for (int i = 0; i < n; ++i) {                  // Iterate through all vertices
        cout << "Vertex " << i << ": Distance " << dist[i] << endl;
    }
}

int main() {
    int vertices, edges;
    cin >> vertices >> edges;

    vector<vector<pair<int, int>>> graph(vertices); // Create an adjacency list for the graph

    for (int i = 0; i < edges; ++i) {
        int u, v, w;
        cin >> u >> v >> w;                         // Input each edge (start_vertex, end_vertex, weight)
        graph[u].emplace_back(v, w);

    }

    int source = 0;
    dijkstra(graph, source);

    return 0;
}


