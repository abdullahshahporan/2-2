/*
//Jacobi iteration method
#include<bits/stdc++.h>
using namespace std;
bool is_diag(vector<vector<double>>&coeffs,int n)
{
    for(int i=0;i<n;i++)
    {
        double diag=abs(coeffs[i][i]);
        double sum=0;
        for(int j=0;j<n;j++)
        {
            if(i!=j)
            {
                sum+= abs(coeffs[i][j]);
            }
            
        }
        if(diag<sum)
            return false;
    }
    return true;
}
 bool m_diag(vector<vector<double>>&coeffs,int n)
 {
    for(int i=0;i<n;i++)
    {
        double mdiag= abs(coeffs[i][i]);
        int mR=i;
        for(int j=i+1;j<n;j++)
        {
           if(abs(coeffs[j][i])>mdiag)
           {
            mdiag=abs(coeffs[j][i]);
            mR=j;
           }
        }
        if(mR!=i)
        {
            swap(coeffs[i],coeffs[mR]);
        }
    }
    return is_diag(coeffs,n);
 }
double sol(vector<vector<double>>&coefs, vector<double>&cval,int j,int n)
{
    double cons=coefs[j][n];
    for(int i=0;i<n;i++)
    {
        if(j==i)
        {
         continue;
        }
        cons-=coefs[j][i]*cval[j];
    }
    return cons/coefs[j][j];
}

int main()
{
    int n;
    cin>>n;
    vector<vector<double>>coefs(n,vector<double>(n+1));
    for( int i=0;i<n;i++)
    {
        for(int j=0;j<=n;j++)
        {
           cin>>coefs[i][j];

        }
    }
    if(is_diag(coefs,n)==false)
    {
        if(m_diag(coefs,n)==false)
        {
            cout<<"Not dominant"<<endl;
            return 0;
        }
    }

    double tol=1e-6;
    int mIt=2000;
    vector<double>cval(n,0),nval(n,0);
    for(int i=0;i<n;i++)
    {
        bool f=true;
        for(int j=0;j<n;j++)
        {
            nval[j]=sol(coefs,cval,j,n);
            if(abs(nval[j]-cval[j]>tol))
            f=false;
        }
        cval=nval;
        if(f==true)
        break;
    }
    for(int i=0;i<n;i++)
    {
        cout<<cval[i]<<' ';
    }
    cout<<endl<<endl;

}*/
/*
//Gauss-Seidel Iterative Method
#include<bits/stdc++.h>
using namespace std;
bool is_diag(vector<vector<double>>&coeffs,int n)
{
    for(int i=0;i<n;i++)
    {
        double diag=abs(coeffs[i][i]);
        double sum=0;
        for(int j=0;j<n;j++)
        {
            if(i!=j)
            {
                sum+= abs(coeffs[i][j]);
            }
            
        }
        if(diag<sum)
            return false;
    }
    return true;
}
 bool m_diag(vector<vector<double>>&coeffs,int n)
 {
    for(int i=0;i<n;i++)
    {
        double mdiag= abs(coeffs[i][i]);
        int mR=i;
        for(int j=i+1;j<n;j++)
        {
           if(abs(coeffs[j][i])>mdiag)
           {
            mdiag=abs(coeffs[j][i]);
            mR=j;
           }
        }
        if(mR!=i)
        {
            swap(coeffs[i],coeffs[mR]);
        }
    }
    return is_diag(coeffs,n);
 }
double sol(vector<vector<double>>&coefs, vector<double>&cval,int j,int n)
{
    double cons=coefs[j][n];
    for(int i=0;i<n;i++)
    {
        if(j==i)
        {
         continue;
        }
        cons-=coefs[j][i]*cval[j];
    }
    return cons/coefs[j][j];
}

int main()
{
    int n;
    cin>>n;
    vector<vector<double>>coefs(n,vector<double>(n+1));
    for( int i=0;i<n;i++)
    {
        for(int j=0;j<=n;j++)
        {
           cin>>coefs[i][j];

        }
    }
    if(is_diag(coefs,n)==false)
    {
        if(m_diag(coefs,n)==false)
        {
            cout<<"Not dominant"<<endl;
            return 0;
        }
    }

    double tol=1e-6;
    int mIt=2000;
    vector<double>cval(n,0),nval(n,0);
    for(int i=0;i<n;i++)
    {
        bool f=true;
        for(int j=0;j<n;j++)
        {
            nval[j]=sol(coefs,cval,j,n);
            if(abs(nval[j]-cval[j]>tol))
            f=false;
           cval[j]=nval[j]; 
        }
       
        if(f==true)
        break;
    }
    for(int i=0;i<n;i++)
    {
        cout<<cval[i]<<' ';
    }
    cout<<endl<<endl;

}
*/
/*
//Gauss Elimination Method
#include<bits/stdc++.h>
using namespace std;
bool GaussElimination(vector<vector<double>> &coefficients)
{
    int n=coefficients.size();
    bool flag =false;
    for(int i=0;i<n;i++)
    {
        if(coefficients[i][i] == 0)
        {
            int c =1;
            while((i+c)<n && coefficients[i+c][i] == 0)
            c++;
            if((i+c)==n)
            {
                flag=1; break;
            }
            swap(coefficients[i],coefficients[i+c]); // swapping the entire row
          
        }
        for(int j=i+1; j<n ;j++)
        {
            double frac = coefficients[j][i]/coefficients[i][i];
            for(int k=0;k<=n;k++)
            {
                coefficients[j][k] = coefficients[j][k] - (coefficients[i][k]*frac); // making non diagonal element to zero
            }
        }
    }
    return flag;
}
int checkConsistency(vector<vector<double>> &coefficients)
{
    int n=coefficients.size();
    double sum=0;
    int flag =3; // flag =3 -> assuming no solutions by default
    for(int i=0;i<n;i++)
    {
        sum=0;
        for(int j=0;j<n;j++)
        {
            sum += coefficients[i][j];
        }
        if(sum == coefficients[i][n])
        flag=2; // flag = 2 -> infinite solutions
        
    }
    return flag;
}
void Backsubstitution(vector<vector<double>> &coefficients)
{
    int n=coefficients.size();
    vector<double> ans(n);
    for(int i=n-1;i>=0;i--)
    {
        ans[i]=coefficients[i][n];
        for(int j=i+1;j<n;j++)
        {
            ans[i]  -=coefficients[i][j]*ans[j];
        }
        ans[i] /= coefficients[i][i]; // getting answer
    }
    cout<<"The solutions are(Gauss-Elimination method): ";
    for(int i=0;i<n;i++)
    {
        cout<<ans[i]<<' ';
    }
    cout<<endl;
}
int main()
{
    cout<< setprecision(5)<<fixed;
    cout<<"Enter the number of unknowns:";
    int n;
    cin>>n;
    vector<vector<double>>coefficients(n, vector<double>(n + 1));
  

    cout<<"Enter the coefficients and constants of the equations(Form: ax + by +cz + .. = d):" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            cin >> coefficients[i][j];
        }
    }
    int temp=0;

    bool flag=GaussElimination(coefficients);
    if(flag==true)
    temp=checkConsistency(coefficients);
    if(temp==2){
        cout<<"The system of Linear equations has infinite solutions."<<endl;
        return 0;
    }
    if(flag==3)
    {
        cout<<"The system has no solutions."<<endl;
        return 0;
    }
    cout<<"Matrix after Gaussian Elimination."<<endl;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<=n;j++)
        {
            cout<<coefficients[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<setprecision(5)<<fixed<<endl;

    Backsubstitution(coefficients);
    
}
*/
/*
//Bisection Method
#include <bits/stdc++.h>

using namespace std;

double polynomial_bs(double x, const vector<double> &coefficients)
{
    double result = 0.0;
    int degree = coefficients.size() - 1;
    for (int i = 0; i <= degree; ++i)
    {
        result += coefficients[i] * pow(x, degree - i);
    }
    return result;
}

double find_Root_InInterval(double a, double b, const vector<double> &coefficients, double epsilon, int max_iterations)
{
    double f_a = polynomial_bs(a, coefficients);
    double f_b = polynomial_bs(b, coefficients);

    if (f_a * f_b > 0)
    {
        return NAN;
    }

    int iteration = 0;
    double c;
    while (fabs(b - a) > epsilon && iteration < max_iterations)
    {
        c = (a + b) / 2;
        double f_c = polynomial_bs(c, coefficients);

        if (fabs(f_c) < epsilon)
        {
            return c;
        }

        if (f_a * f_c < 0)
        {
            b = c;
        }
        else
        {
            a = c;
            f_a = f_c;
        }
        iteration++;
    }

    return (a + b) / 2;
}

int main()
{
    int power;
    cout << "Enter the degree (power) of the polynomial: ";
    cin >> power;

    vector<double> coefficients(power + 1);
    cout << "Enter the coefficients (highest to lowest degree): ";
    for (int i = 0; i <= power; ++i)
    {
        cin >> coefficients[i];
    }

    double epsilon = 1e-6;
    int max_iterations = 1000;
    vector<double> roots;

    double a_n = coefficients[0];
    double a_n_minus_1 = (power >= 1) ? coefficients[1] : 0;
    double x_max = abs(-a_n_minus_1 / a_n);

    double interval_limit = x_max;
    if (power >= 2)
    {
        double a_n_minus_2 = coefficients[2];
        interval_limit = sqrt(pow(a_n_minus_1 / a_n, 2) - 2 * (a_n_minus_2 / a_n));
    }

    double lower_bound = -interval_limit;
    double upper_bound = interval_limit;

    int num_subintervals = 200;
    double step = (upper_bound - lower_bound) / num_subintervals;

    for (int i = 0; i < num_subintervals; ++i)
    {
        double a = lower_bound + i * step;
        double b = a + step;

        double root = find_Root_InInterval(a, b, coefficients, epsilon, max_iterations);
        if (!isnan(root))
        {
            bool is_repeated = false;
            for (const double &r : roots)
            {
                if (fabs(r - root) < epsilon)
                {
                    is_repeated = true;
                    break;
                }
            }
            if (!is_repeated)
            {
                roots.push_back(root);
            }
        }
    }

    cout << "Approximate roots of the polynomial are:\n";
    for (const auto &root : roots)
    {
        cout << root << endl;
    }
}*/
/*
//false position method
#include <bits/stdc++.h>

using namespace std;

double TOLERANCE_false_position = 1e-6;
int degree_false_position;
vector<double> coeffs_false_position;

double evaluate_false_position(double x)
{
    double result = 0.0;
    int degree = coeffs_false_position.size() - 1;
    for (int i = 0; i <= degree; ++i)
    {
        result += coeffs_false_position[i] * pow(x, degree - i);
    }
    return result;
}

double findRootFalsePosition(double a, double b, const vector<double> &coefficients, double epsilon, int max_iterations)
{
    double f_a = evaluate_false_position(a);
    double f_b = evaluate_false_position(b);

    if (f_a * f_b > 0)
    {
        return NAN;
    }

    double c;
    int iteration = 0;
    while (fabs(b - a) > epsilon && iteration < max_iterations)
    {
        c = (a * f_b - b * f_a) / (f_b - f_a);
        double f_c = evaluate_false_position(c);

        if (fabs(f_c) < epsilon)
        {
            return c;
        }

        if (f_a * f_c < 0)
        {
            b = c;
            f_b = f_c;
        }
        else
        {
            a = c;
            f_a = f_c;
        }
        iteration++;
    }

    return (a + b) / 2;
}

int main()
{
    cout << "Enter the degree (power) of the polynomial: ";
    cin >> degree_false_position;
    coeffs_false_position.resize(degree_false_position + 1);

    cout << "Enter the coefficients (highest to lowest degree): ";
    for (int i = 0; i <= degree_false_position; i++)
    {
        cin >> coeffs_false_position[i];
    }

    double epsilon = 1e-6;
    int max_iterations = 1000;
    vector<double> roots;

    double a_n = coeffs_false_position[0];
    double a_n_minus_1 = (degree_false_position >= 1) ? coeffs_false_position[1] : 0;
    double x_max = abs(-a_n_minus_1 / a_n);

    double interval_limit = x_max;
    if (degree_false_position >= 2)
    {
        double a_n_minus_2 = coeffs_false_position[2];
        interval_limit = sqrt(pow(a_n_minus_1 / a_n, 2) - 2 * (a_n_minus_2 / a_n));
    }

    double lower_bound = -interval_limit;
    double upper_bound = interval_limit;

    int num_subintervals = 200;
    double step = (upper_bound - lower_bound) / num_subintervals;

    for (int i = 0; i < num_subintervals; ++i)
    {
        double a = lower_bound + i * step;
        double b = a + step;

        double root = findRootFalsePosition(a, b, coeffs_false_position, epsilon, max_iterations);
        if (!isnan(root))
        {
            bool is_repeated = false;
            for (const double &r : roots)
            {
                if (fabs(r - root) < epsilon)
                {
                    is_repeated = true;
                    break;
                }
            }
            if (!is_repeated)
            {
                roots.push_back(root);
            }
        }
    }

    cout << "Roots found using False Position Method:\n";
    for (const auto &root : roots)
    {
        cout << fixed << setprecision(6) << root << endl;
    }
}*/
/*
//secant method
#include <bits/stdc++.h>
using namespace std;

double polynomial(double x, const vector<double> &coefficients)
{
    double result = 0.0;
    int degree = coefficients.size() - 1;
    for (int i = 0; i <= degree; ++i)
    {
        result += coefficients[i] * pow(x, degree - i);
    }
    return result;
}

double findRootInInterval(double x0, double x1, const vector<double> &coefficients, double epsilon, int max_iterations)
{
    double f_x0 = polynomial(x0, coefficients);
    double f_x1 = polynomial(x1, coefficients);

    int iteration = 0;
    double x2;

    while (iteration < max_iterations)
    {
        if (fabs(f_x1 - f_x0) < epsilon)
        {
            return NAN;
        }

        x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0);
        double f_x2 = polynomial(x2, coefficients);

        if (fabs(f_x2) < epsilon)
        {
            return x2;
        }

        x0 = x1;
        f_x0 = f_x1;
        x1 = x2;
        f_x1 = f_x2;

        iteration++;
    }

    return x2;
}

int main()
{
    int power;
    cout << "Enter the degree (power) of the polynomial: ";
    cin >> power;

    vector<double> coefficients(power + 1);
    cout << "Enter the coefficients (highest to lowest degree): ";
    for (int i = 0; i <= power; ++i)
    {
        cin >> coefficients[i];
    }

    double epsilon = 1e-6;
    int max_iterations = 1000;
    vector<double> roots;

    double a_n = coefficients[0];
    double a_n_minus_1 = (power >= 1) ? coefficients[1] : 0;
    double x_max = abs(-a_n_minus_1 / a_n);

    double interval_limit = x_max;
    if (power >= 2)
    {
        double a_n_minus_2 = coefficients[2];
        interval_limit = sqrt(pow(a_n_minus_1 / a_n, 2) - 2 * (a_n_minus_2 / a_n));
    }

    double lower_bound = -interval_limit;
    double upper_bound = interval_limit;

    int num_subintervals = 100;
    double step = (upper_bound - lower_bound) / num_subintervals;

    for (int i = 0; i < num_subintervals; ++i)
    {
        double x0 = lower_bound + i * step;
        double x1 = x0 + step;

        double root = findRootInInterval(x0, x1, coefficients, epsilon, max_iterations);

        if (!isnan(root))
        {
            roots.push_back(root);
        }
    }

    set<double> uniqueRoots;
    for (const auto &root : roots)
    {
        bool foundClose = false;
        for (const auto &uniqueRoot : uniqueRoots)
        {
            if (fabs(root - uniqueRoot) < 0.1)
            {
                foundClose = true;
                break;
            }
        }
        if (!foundClose)
        {
            uniqueRoots.insert(root);
        }
    }

    cout << fixed << setprecision(4);
    cout << "Approximate roots of the polynomial are:\n";
    for (const auto &root : uniqueRoots)
    {
        cout << root << endl;
    }
}*/

//Newton-Raphson Method
#include <bits/stdc++.h>

using namespace std;

double polynomial_newton(double x, const vector<double> &coefficients)
{
    double result = 0.0;
    int degree = coefficients.size() - 1;
    for (int i = 0; i <= degree; ++i)
    {
        result += coefficients[i] * pow(x, degree - i);
    }
    return result;
}

double polynomialDerivative(double x, const vector<double> &coefficients)
{
    double result = 0.0;
    int degree = coefficients.size() - 1;
    for (int i = 0; i < degree; ++i)
    {
        result += coefficients[i] * (degree - i) * pow(x, degree - i - 1);
    }
    return result;
}

double newtonRaphson(double initial_guess, const vector<double> &coefficients, double epsilon, int max_iterations)
{
    double x_n = initial_guess;
    for (int iteration = 0; iteration < max_iterations; ++iteration)
    {
        double f_x_n = polynomial_newton(x_n, coefficients);
        double f_prime_x_n = polynomialDerivative(x_n, coefficients);

        if (fabs(f_prime_x_n) < epsilon)
        {
            return NAN;
        }

        double x_n_plus_1 = x_n - f_x_n / f_prime_x_n;

        if (fabs(x_n_plus_1 - x_n) < epsilon)
        {
            return x_n_plus_1;
        }

        x_n = x_n_plus_1;
    }
    return NAN;
}

int main()
{
    int power;
    cout << "Enter the degree (power) of the polynomial: ";
    cin >> power;

    vector<double> coefficients(power + 1);
    cout << "Enter the coefficients (highest to lowest degree): ";
    for (int i = 0; i <= power; ++i)
    {
        cin >> coefficients[i];
    }

    double epsilon = 1e-6;
    int max_iterations = 1000;
    set<double> uniqueRoots;

    for (double i = -10; i <= 20; i += 1)
    {
        double root = newtonRaphson(i, coefficients, epsilon, max_iterations);
        if (!isnan(root))
        {
            bool isUnique = true;
            for (const auto &uniqueRoot : uniqueRoots)
            {
                if (fabs(root - uniqueRoot) < epsilon)
                {
                    isUnique = false;
                    break;
                }
            }

            if (isUnique)
            {
                uniqueRoots.insert(root);
            }
        }
    }

    cout << fixed << setprecision(4);
    cout << "Approximate roots of the polynomial are:\n";
    for (const auto &root : uniqueRoots)
    {
        cout << root << " ";
    }
    cout << endl;
}