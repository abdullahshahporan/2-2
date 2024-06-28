#include<iostream>
#include<queue>
#include<vector>
using namespace std;
void bfs(vector<vector<int>> &adj, int s_node, vector<bool> &vis)
{
    queue <int> q;
    vis[s_node]=true;
    q.push(s_node);

    while(!q.empty())
    {
        int curr=q.front();
        q.pop();
        cout<<curr<<" ";
        for(int neigh: adj[curr])
        {
            if(!vis[neigh])
            {
                vis[neigh]=true;
                q.push(neigh);
            }
        }
    }
}
void add(vector<vector<int>> & adj, int u, int v)
{
    adj[u].push_back(v);
}
int main()
{
     int ver=5;
     vector<vector<int>> adj(ver);
     add(adj,0,1);
      add(adj,0,2);
       add(adj,1,3);
        add(adj,1,4);
         add(adj,2,4);
         vector<bool>vis(ver,false);
         bfs(adj,0,vis);

}
