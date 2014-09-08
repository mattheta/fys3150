#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <time.h>


using namespace arma;
using namespace std;

void tridiagonal(mat temp, mat c, mat b, mat a, mat &v, double h, int N);
double f(double x);
void exact(mat &u, double h, int N);
void error(mat v, mat u, mat epsilon, double h, int N);
void LUdecomposition(int N, double h);
int main()
{
    int N;
    double h, duration;
    mat temp,b,c,v,epsilon,a,u;
    clock_t start, finish;
    start = clock();

    N = 100000;                            //matrix size
    cout<<N<<endl;
    h = 1.0/(N+1);                     //step size
    b = ones<mat>(1,N+2);
    v = zeros<mat>(1,N+2);
    u = zeros<mat>(1,N+2);
    c = ones<mat>(1,N+2);
    a = ones<mat>(1,N+2);
    temp = zeros<mat>(1,N+2);
    epsilon = zeros<mat>(1,N+2);

    b = b*2;
    c =c*-1;
    a =a*-1;

    //tridiagonal(temp,c,b,a,v,h, N);
    //exact(u,h,N);
    //error(v,u,epsilon,h,N);
    // start = clock();
    LUdecomposition(N,h);
    //finish = clock();
    //duration=((finish-start)/(double)CLOCKS_PER_SEC);
    //cout<<"duration="<<duration<<endl;

    return 0;
}

double f(double x)
{
    double value;
    value = 100*exp(-10*x);
    return value;
}
void exact(mat &u, double h, int N)
{
    int i;

    for (i=0; i<=N; i++)
    {
        u[i] = 1- (1- exp(-10))*i*h - exp(-10*i*h);
    }

    ofstream myfile;
    myfile.open("exact.txt");

    for (i=0; i<=N+1; i++)
    {
        myfile <<((i)*h)<<"   "<<u[i]<<"   "<<endl;
    }
    myfile.close();

    return;
}
void tridiagonal(mat temp, mat c, mat b, mat a, mat &v, double h, int N)
{
    //Initial parameters
    double btemp;
    int i;
    btemp = b[1];
    v[1] = f(1*h)/btemp;

    for (i=2; i<= N; i++)
    {
        temp[i] = c[i-1]/btemp;
        btemp = b[i]-a[i]*temp[i];
        v[i] = (f(i*h) - a[i]*v[i-1])/btemp;
    }

    for (i=N-1; i>=1; i--)
    {
        v[i] = v[i] - temp[i+1]*v[i+1];
    }

    v = v*h*h;

    // write to file

    ofstream myfile;
    myfile.open("test.txt");

    for (i=0; i<=N+1; i++)
    {
        myfile <<((i)*h)<<"   "<<v[i]<<"   "<<endl;
    }
    myfile.close();

    return;
}
void error(mat v, mat u, mat epsilon,double h ,int N)
{
    int i;

    for (i=0; i<=N+1; i++)
    {
        epsilon[i] = log(abs((v[i]-u[i])/u[i]));
    }

    //plotting
    ofstream myfile;
    myfile.open("error.txt");

    for (i=0; i<=N+1; i++)
    {
        myfile <<((i)*h)<<"   "<<epsilon[i]<<"   "<<endl;
    }
    myfile.close();

    cout<<"max error value: "<<epsilon.max()<<endl;

}
void LUdecomposition(int N, double h)
{
    //LUdecomposition used to solve the linear equation Ax=w

    int i,j;
    mat L,U,P,A,x,y,w;
    clock_t start, finish;
    double duration;

    A = zeros<mat>(N,N);
    U = zeros<mat>(N,N);
    P = zeros<mat>(N,N);
    L = zeros<mat>(N,N);
    x = zeros<mat>(N,1);
    y = zeros<mat>(N,1);
    w = zeros<mat>(N,1);

    start = clock();

    for(i=0; i<N; i++)
    {
        w[i] = h*h*f((i+1)*h);
    }



    //making the correct matrix
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            A(i,j) = 0;

            if(i==j)
            {
                A(i,j)=2;
            }
            else if(j+1==i)
            {
                A(i,j)=-1;
            }
            else if(j-1==i)
            {
                A(i,j)=-1;
            }
        }
    }

    lu(L,U,P,A);
    y = solve(L,w);
    x = solve(U,y);
    finish = clock();
    duration=((finish-start)/(double)CLOCKS_PER_SEC);
    cout<<"LU duration="<<duration<<endl;

    return;
}
