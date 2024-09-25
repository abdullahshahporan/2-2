#include<bits/stdc++.h>
using namespace std;
const int MAX_SIZE=100;
void print(double a[MAX_SIZE][MAX_SIZE],result[MAX_SIZE])
int main()
{
    cout<<"Enter No. of equation"<<endl;
    int n;
    cin>>n;
    double a[MAX_SIZE][MAX_SIZE],result[MAX_SIZE];
    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j)
        {
            cin>>a[i][j];
        }
        cin>>result[i];
    }
}
