#include<bits/stdc++.h>
using namespace std;

// Function to evaluate the polynomial at value 'a'
double func_value(vector<double>&v1, double a){
    double sum = 0.0;
    for(int i = 0; i < v1.size(); i++){
        sum += v1[i] * pow(a, v1.size()-i-1);
    }
    return sum;
}

// Function to compute the derivative of a polynomial
vector<double> diff(vector<double> &dif){
    int n = dif.size();
    vector<double> v1(n-1);
    for(int i = 0; i < n-1; i++){
        v1[i] = dif[i] * (n-i-1);
    }
    return v1;
}

// Deflate the polynomial by dividing by (x - root)
vector<double> deflate(vector<double> &poly, double root){
    int n = poly.size();
    vector<double> new_poly(n-1);
    new_poly[0] = poly[0];  // leading coefficient stays the same
    for (int i = 1; i < n-1; i++) {
        new_poly[i] = poly[i] + root * new_poly[i-1];
    }
    return new_poly;
}

// Secant method to find one root of the polynomial
double secant_method(vector<double> &v1, double tol, int max_iter){
    double a = 1.0, b = 10.0, e;
    for(int i = 1; i <= max_iter; i++){
        double f_a = func_value(v1, a);
        double f_b = func_value(v1, b);
        e = b - ((f_b * (b - a)) / (f_b - f_a));
        double f_e = func_value(v1, e);
        cout << "Secant Iter No. " << i << " a = " << a << " b = " << b << " c = " << e << endl
             << "Secant Functional value f_c = " << f_e << endl;
        if(fabs(f_e) < tol){
            cout << "Secant Method Found Root " << e << endl;
            return e;
        }
        a = b;
        b = e;
    }
    return e;
}

// Newton method to find one root of the polynomial
double newton_method(vector<double> &v1, vector<double> &dif, double tol, int max_iter){
    double a = 4.0, e, f_a, f_p_a;
    for(int i = 1; i <= max_iter; i++){
        f_a = func_value(v1, a);
        f_p_a = func_value(dif, a);
        if(fabs(f_p_a) < tol){
            cout << "Newton: Derivative too small, stopping iteration." << endl;
            return a;
        }

        e = a - (f_a / f_p_a);
        double f_e = func_value(v1, e);
        cout << "Newton Iter No. " << i << " a = " << a << " c = " << e << endl
             << "Newton Functional value f_c = " << f_e << endl;

        if(fabs(f_e) < tol){
            cout << "Newton Method Found Root " << e << endl;
            return e;
        }
        a = e;
    }
    return a;
}

// Function to find all real roots using Secant Method
void find_all_roots_secant(vector<double> &v1, double tol, int max_iter) {
    int degree = v1.size() - 1;
    vector<double> secant_roots;
    vector<double> v_secant = v1;  // Copy for Secant

    cout << "Using Secant Method to find all roots:" << endl;
    while(degree > 0) {
        // Find one root using Secant Method
        double secant_root = secant_method(v_secant, tol, max_iter);
        secant_roots.push_back(secant_root);

        // Deflate the polynomial
        v_secant = deflate(v_secant, secant_root);
        degree = v_secant.size() - 1;  // Update degree after deflation
    }

    // Output all roots found by Secant method
    cout << "\nAll roots found by Secant Method:" << endl;
    for(double root : secant_roots) {
        cout << root << " ";
    }
    cout << endl;
}

// Function to find all real roots using Newton's Method
void find_all_roots_newton(vector<double> &v1, double tol, int max_iter) {
    int degree = v1.size() - 1;
    vector<double> newton_roots;
    vector<double> v_newton = v1;  // Copy for Newton

    cout << "Using Newton's Method to find all roots:" << endl;
    while(degree > 0) {
        vector<double> dif = diff(v_newton);  // Get the derivative for Newton's method

        // Find one root using Newton's Method
        double newton_root = newton_method(v_newton, dif, tol, max_iter);
        newton_roots.push_back(newton_root);

        // Deflate the polynomial
        v_newton = deflate(v_newton, newton_root);
        degree = v_newton.size() - 1;  // Update degree after deflation
    }

    // Output all roots found by Newton method
    cout << "\nAll roots found by Newton's Method:" << endl;
    for(double root : newton_roots) {
        cout << root << " ";
    }
    cout << endl;
}

int main() {
    int n;
    cout << "Enter the degree of the polynomial: ";
    cin >> n;

    vector<double> v(n+1);  // Polynomial coefficients

    cout << "Enter the coefficients of the polynomial (highest degree first): ";
    for(int i = 0; i < n+1; i++){
        cin >> v[i];
    }

    double tol = 0.0001;
    int max_iter = 100;

    // Call Newton's Method to find all roots
    find_all_roots_newton(v, tol, max_iter);

    cout << endl;

    // Call Secant Method to find all roots
    find_all_roots_secant(v, tol, max_iter);

    return 0;
}
