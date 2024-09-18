#include<bits/stdc++.h>
using namespace std;

double func_value(vector<int>&v1, double a){
    double sum=0.0;
    for(int i = 0;i<v1.size();i++){
        sum = sum + v1[i] * pow(a, v1.size()-i-1);
    }
    return sum;
}

vector<int> diff(vector<int> &dif){
    int n = dif.size(),j=0;
    vector<int> v1(n-1);
    for(int i=n-1;i>0;i--){
        v1[j] = dif[j] * (i);
        j++;
    }
    return v1;
}

void secant_method(vector<int> &v1,double tol, int max_iter){
    double a=0.5,b=1.5,e;
    for(int i=1;i<=max_iter;i++){
        double f_a = func_value(v1,a);
        double f_b = func_value(v1,b);
        e = b - ((f_b * (b-a))/(f_b-f_a));
        double f_e = func_value(v1,e);
        //cout << "Iter No. " << i << " a = " << a << " b = " << b << " c = " << e << endl << "Functional value " << " f_c = " << f_e << endl;
        cout<<f_e<<endl;
        if(fabs(f_e)<tol){
            cout << "Found Root " << e << endl;
            return;
        }
        swap(a,b);
        swap(b,e);
    }
}

void newton(vector<int> &v1, vector<int> &dif, double tol, int max_iter){
    double a = 4.0, e, f_a, f_p_a;
    for(int i = 1; i <= max_iter; i++){
        f_a = func_value(v1, a);
        f_p_a = func_value(dif, a);
        if(fabs(f_p_a) < tol){
            cout << "Derivative too small, stopping iteration." << endl;
            return;
        }

        e = a - (f_a / f_p_a);
        double f_e = func_value(v1, e);
        //cout << "Iter No. " << i << " a = " << a << " c = " << e << endl << "Functional value " << " f_c = " << f_e << endl;;
cout<<f_e<<endl;
        if(fabs(f_e) < tol){
            cout << "Found Root " << e << endl;
            return;
        }
        a = e;
    }
}

int main(){
    int n;
    cout << "Enter the no. of degree : "; cin >> n;
    vector<int> v(n+1),dif;

    for(int i=0;i<n+1;i++){
        cin >> v[i];
    }
    dif = diff(v);
  //  vector<int> v1= {1,-5,6},dif={2,-5};
    double tol = 0.0001;
    cout << "Secant Method : " << endl << endl;
    secant_method(v,tol,10);
    cout << endl << endl;
    cout << "Newton Method : " << endl << endl;
    newton(v,dif,tol,10);
    return 0;
}
