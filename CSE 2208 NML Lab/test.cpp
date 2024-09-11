//jacobi
/*#include<bits/stdc++.h>
#include<iomanip>
#include<math.h>
using namespace std;
int main()
{
    double a1,b1,c1,d1;
   double a2,b2,c2,d2;
   double a3,b3,c3,d3;
   cin>>a1>>b1>>c1>>d1;
   cin>>a2>>b2>>c2>>d2;
   cin>>a3>>b3>>c3>>d3;
    double xp=0,yp=0,zp=0,xn,yn,zn,ex=0,ey=0,ez=0;
    int steps=20;
    cout<<"x\t\ty\t\tz\t\tex\t\tey\t\tez"<<endl;
    do{
        cout<<fixed<<setprecision(4)<<xp<<"\t\t"<<yp<<"\t\t"<<zp<<"\t\t"<<abs(ex)<<"\t\t"<<abs(ey)<<"\t\t"<<abs(ez)<<endl;
        xn=((1/a1)*(d1-b1*yp-c1*zp));
        yn=((1/b2)*(d2-a2*xp-c2*zp));
        zn=((1/c3)*(d3-a3*xp-b3*yp));
        ex=xn-xp;
        ey=yn-yp;
        ez=zn-zp;
        xp=xn;
        yp=yn;
        zp=zn;
    }while(steps--);

}*/
//gauss
/*#include<bits/stdc++.h>
using namespace std;

int main()
{
    double a1,b1,c1,d1;
   double a2,b2,c2,d2;
   double a3,b3,c3,d3;
   cin>>a1>>b1>>c1>>d1;
   cin>>a2>>b2>>c2>>d2;
   cin>>a3>>b3>>c3>>d3;
    double xp=0,yp=0,zp=0,xn,yn,zn,ex=0,ey=0,ez=0;
    int steps=20;
    cout<<"x\t\ty\t\tz\t\tex\t\tey\t\tez"<<endl;
     do{
        cout<<fixed<<setprecision(4)<<xp<<"\t\t"<<yp<<"\t\t"<<zp<<"\t\t"<<abs(ex)<<"\t\t"<<abs(ey)<<"\t\t"<<abs(ez)<<endl;
        xn=((1/a1)*(d1-b1*yp-c1*zp));
        ex=xn-xp;
         xp=xn;
        yn=((1/b2)*(d2-a2*xp-c2*zp));
        ey=yn-yp;
         yp=yn;
        zn=((1/c3)*(d3-a3*xp-b3*yp));
        ez=zn-zp;
        zp=zn;
    }while(steps--);
}
*/
//bisection
/*#include<bits/stdc++.h>
using namespace std;
double f(double x)
{
    return x*x*x-x-1;
}
double bisection(double a, double b, double tolerance)
{
    if(f(a)*f(b)>=0)
    {
        cout<<"Wrong Guess"<<endl;
        return -1;
    }
    double c=a;
    double prev_c=c;
    int iter=1;
    while((b-a)>=tolerance)
    {
        c=(a+b)/2;
        double error=fabs(c-prev_c);
        cout<<fixed<<iter++<<"\t\t\t"<<setprecision(6)<<c<<"\t\t\t"<<error<<endl;
        if(f(c)==0.0||fabs(b-a)<tolerance)
        {
            break;
        }
        if(f(c)*f(a)<0)
        {
            b=c;
        }
        else
        {
            a=c;
        }
        prev_c=c;
        if(iter>=21)
            break;
    }
    return c;
}
int main()
{
    double a,b;
    cin>>a>>b;
    double tolerance=1e-5;
    double root=bisection(a,b,tolerance);
    if(root!=-1)
    {
        cout<<"Root is: "<<setprecision(4)<<root<<endl;
    }
}*/
//False-Positioning
#include<bits/stdc++.h>
using namespace std;
double f(double x)
{
    return x*x*x-x-1;
}
double false_position(double a,double b,double tolerance)
{
    if(f(a)*f(b)>=0)
    {
        cout<<"Wrong think"<<endl;
        return -1;
    }
    double c=a;
    double prev_c=c;
    int iter=1;
    while((b-a)>=tolerance)
    {
        //c=b-(f(b)*(b-a)) /(f(b)-f(a));
        c=(a*f(b)-b*f(a))/(f(b)-f(a));
        double error=fabs(c-prev_c);
        cout<<fixed<<iter++<<"\t\t\t"<<setprecision(4)<<c<<"\t\t\t"<<error<<endl;
        if(f(c)==0.0)
        {
            break;
        }
        else if(f(c)*f(a)<0)
        {
            b=c;
        }
        else
            a=c;

        prev_c=c;
        if(iter>=21)
            break;
    }
    return c;

}
int main()
{
    double a,b;
    cin>>a>>b;
    double tolerance=1e-5;
    double c=false_position(a,b,tolerance);
    if(c!=-1)
    {
        cout<<"Root is "<<setprecision(4)<<c<<endl;
    }
}
