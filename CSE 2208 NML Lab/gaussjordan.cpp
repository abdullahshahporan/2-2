#include<bits/stdc++.h>

#define SIZE 10

using namespace std;

void gaussJordan(float a[SIZE][SIZE], float x[SIZE], int n) {
    float ratio;

    for(int i = 1; i <= n; i++) {
        if(a[i][i] == 0.0) {
            cout << "Mathematical Error!";
            exit(0);
        }
        for(int j = 1; j <= n; j++) {
            if(i != j) {
                ratio = a[j][i] / a[i][i];
                for(int k = 1; k <= n + 1; k++) {
                    a[j][k] = a[j][k] - ratio * a[i][k];
                }
            }
        }
    }

    for(int i = 1; i <= n; i++) {
        x[i] = a[i][n+1] / a[i][i];
    }
}

int main() {
    float a[SIZE][SIZE], x[SIZE];
    int n;

    cout << setprecision(3) << fixed;

    cout << "Enter number of unknowns: ";
    cin >> n;

    cout << "Enter Coefficients of Augmented Matrix: " << endl;
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n + 1; j++) {
            //cout << "a[" << i << "][" << j << "] = ";
            cin >> a[i][j];
        }
    }

    gaussJordan(a, x, n);

    cout << endl << "Solution: " << endl;
    for(int i = 1; i <= n; i++) {
        cout << "x[" << i << "] = " << x[i] << endl;
    }

    return 0;
}

