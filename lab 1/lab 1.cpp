#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#define PI 3.14159265
#define RAD PI / 180

using namespace std;


int const n=5, m=5, p=(n+1)*(m+1), N=2*n*m;
double X[p], Y[p];
int T[N][3];

double phi(double x)
{
    return x * x - 3 * x - 2;
};

double psi(double x)
{
    return x * (3 - x);
};


double square(double x1,double x2, double x3, double y1, double y2, double y3)
{
	return fabs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))/2.0;
}

void savetostl(double XX[p], double YY[p], int T[N][3],  int p, int N)
{
    int k=0, ia, ib, ic;
    double v=0;
    fstream TRstl;
    TRstl.open("triangulation.stl",ios::out|ios::app);
    TRstl.close();
    TRstl.open("triangulation.stl",ios::out|ios::in);
    TRstl << "solid <Triangulation>\n";
    
    for(k=0;k<N;k++)
    {
        ia=T[k][0];
        ib=T[k][1];
        ic=T[k][2];
        
        TRstl << "facet normal "<< 0.0 <<" "<< 0.0 <<" "<< 1.0 <<"\n";
        TRstl <<"outer loop\n";
        TRstl <<"vertex ";
        TRstl << XX[ia] <<" "<< YY[ia] <<" "<< 0.0 <<"\n";
        TRstl << "vertex ";
        TRstl << XX[ib] <<" "<< YY[ib] <<" "<< 0.0 <<"\n";
        TRstl << "vertex ";
        TRstl << XX[ic] <<" "<< YY[ic] <<" "<< 0.0 <<"\n";
        TRstl << "endloop\n";
        TRstl << "endfacet\n";
        v+=square(XX[ia],XX[ib],XX[ic],YY[ia],YY[ib],YY[ic]);
        
    }
    cout<<"V="<<v<<endl;
    cout<<"V real=15"<<endl;
    TRstl <<"endsolid";
    TRstl.close();
}


int main()
{
    
    double a = 0, b = 3;
    int k = 0;
    
    for(int i=0; i<=n; i++)
    {
        for(int j=0; j<=m; j++)
        {
            X[k]=a+i*(b-a)/n;
            Y[k]= phi(X[k]) + (psi(X[k]) - phi(X[k])) * j / m;
            k++;
        }
    }
    
    k=0;
    
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<m; j++)
        {
        T[k][0] = (m+1)*i+j;
        T[k][1] = (m+1)*(i+1)+j;
        T[k][2] = (m+1)*i+j+1;
        k++;
            
        T[k][0] = (m+1)*i+j+1;
        T[k][1] = (m+1)*(i+1)+j+1;
        T[k][2] = (m+1)*(i+1)+j;
        k++;
        }
    }
    
    
    savetostl(X, Y, T, p, N);
    system("pause");

    return 0;
}
