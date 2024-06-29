// C++ program for above approach
#include <bits/stdc++.h>
 
using namespace std;
 
// Recursive Function to find the
// Maximal Independent Vertex Set
int n;
vector<int> maior;
vector<int> *graphSets(map<int, vector<int> > graph, int n, vector<int> *I)
{
    if (int(I->size()) == n){
        maior.assign(I->begin(), I->end());
        return I;
    }

    // Branching and bound
    if (I->size() + graph.size() < maior.size())
        return NULL;
    // Base Case - Given Graph has no nodes
    if (graph.size() == 0){
        if (I->size() > maior.size())
            maior.assign(I->begin(), I->end());
        return NULL;
    }
    if (I->size() > maior.size())
        maior.assign(I->begin(), I->end());
 
    // Select a vertex from the graph
    int vCurrent = graph.begin()->first;
 
    // copy graph
    map<int, vector<int> > graph2(graph);
    // remove vertex from graph
    graph2.erase(vCurrent);

    // Loop through its neighbours
    for (auto v : graph.at(vCurrent)) {
        // Delete neighbor from the current subgraph
        if (graph2.count(v)) {
            graph2.erase(v);
        }
    }
    // remove from original
    graph.erase(vCurrent);

    // insert vertex at I
    I->push_back(vCurrent);
    // Call recursive as vexter is in MIS
    vector<int> *res1 = graphSets(graph2, n, I);
    if (res1 != NULL)
        return res1;
    I->pop_back();
 
    // backtracking as vertex is not in MIS
    return graphSets(graph, n, I);
}

int get_row(int n, int k){
    return k/n;
}
int get_col(int n, int k){
    return k%n;
}

void proibe_diagonais(vector<bool> &mat, int n){
    for (int i = 0; i < n; i++){
        mat[i*n + i] = 1;
        mat[i*n + n - i - 1] = 1;
    }
}
 
// Driver Code
int main(int argc, char **argv)
{
    if (argc < 2){
        printf("Usage: %s <board_size>\n", argv[0]);
        return 1;
    }

    vector<vector<int>> E;
    n = atoi(argv[1]);
    vector<bool> mat(n*n, 0);
    proibe_diagonais(mat, n);

    for (int i = 0; i < n*n; i++){
        if (mat[i])
            continue;
        int row_i = get_row(n, i);
        int col_i = get_col(n, i);
        for (int j = i+1; j < n*n; j++){
            int row_j = get_row(n, j);
            int col_j = get_col(n, j);
            if (mat[j])
                continue;
            if (row_i == row_j || col_i == col_j || row_i + col_i == row_j + col_j || row_i - col_i == row_j - col_j)
                E.push_back(vector<int>{i, j});
        }
    }
    map<int, vector<int> > graph;
 
    // Constructs Graph as a dictionary of the following
    // format- graph[VertexNumber V] = list[Neighbors of
    // Vertex V]
    for (int i = 0; i < E.size(); i++) {
        int v1 = E[i][0];
        int v2 = E[i][1];
        if (graph.count(v1) == 0) {
            graph[v1] = vector<int>();
        }
        if (graph.count(v2) == 0) {
            graph[v2] = vector<int>();
        }
        graph[v1].push_back(v2);
        graph[v2].push_back(v1);
    }
 
    // Recursive call considering all vertices in the
    // maximum independent set
    vector<int> I;
    graphSets(graph, n, &I);
    vector<int> *maximalIndependentSet = &maior;
 
    // Prints the Result
    for (auto i : *maximalIndependentSet)
        printf("(%d, %d) ", get_row(n, i), get_col(n, i));
    printf("\n");
    int defer = 0;
    for (int i=0; i<n; i++) {
        for (int j=0; j<2*n; j++) {
            defer = 0;
            if (mat[i*n + j/2]){
                printf("\x1b[1;31m█\x1b[m");
                continue;
            }
            else
                for (int r : *maximalIndependentSet){
                    if ((i*n + j/2) == r){
                        printf("\x1b[1;33m█\x1b[m");
                        defer = 1;
                        break;
                    }
                }
            if (defer);
            else if ((j/2 + i) % 2)
                printf("\x1b[1;30m█\x1b[m");
            else
                printf("█");
        }
        printf("\n");
    }
 
    return 0;
}
