#include <bits/stdc++.h>
//#include <cmath>
//This Will Work For Any Equation...
#define f(x, y) sin(x)
using namespace std;

void RK_method(float x0, float y0) {
    float h = 0.1;
    float range = 4 * M_PI;
    int n = (int)((range - x0) / h);

    float k1, k2, k3, k4, yn;

    for (int i = 0; i < n; i++) {
        k1 = h * f(x0, y0);
        k2 = h * f(x0 + (h / 2), y0 + (k1 / 2));
        k3 = h * f(x0 + (h / 2), y0 + (k2 / 2));
        k4 = h * f(x0 + h, y0 + k3);

        yn = y0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        cout <<fixed<<setprecision(5)<< x0 << "\t\t" << y0 << "\t\t" << yn << endl;
       //cout<<yn<<endl;
        x0 = x0 + h;
        y0 = yn;
    }

    cout << "Final Result: " << yn << endl;
}

int main() {
    float x0 = 0, y0 = 0;

    RK_method(x0, y0);

    return 0;
}


