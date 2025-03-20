#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#define PI 3.14159265

using namespace std;


int const n=100, m=100, p=(n+1)*(m+1), N=2*n*m;
double X[p], Y[p];
int T[N][3];

double phi (double x)
{
	return 1;
}

double psi (double x)
{
	return 2;
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
        
    }
    TRstl <<"endsolid";
    TRstl.close();
}


int main()
{
    double tethai,tauj,g,alpha=0,betha=2*PI;
    int k = 0;
    
    for(int i=0; i<=n; i++)
    {
    	tethai=alpha+i*(betha-alpha)/n;
        for(int j=0; j<=m; j++)
        {
        	tauj=j/m;
        	g=tauj*psi(tethai)+(1-tauj)*phi(tethai);
        	
            X[k]= g*cos(tethai);
            Y[k]= g*sin(tethai);
            k++;
        }
    }
    
    k=0;
    if (betha!=2*PI)//сектор
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<m; j++)
        {
        T[k][0] = (m+1)*i+j;
        T[k][1] = (m+1)*(i+1)+j;
        T[k][2] = (m+1)*(i+1)+j+1;
        k++;
        
        T[k][0] = (m+1)*i+j;
        T[k][1] = (m+1)*(i+1)+j+1;
        T[k][2] = (m+1)*i+j+1;
        k++;
        }
    }
    else//кольцо
    {
    for(int i=0; i<n-1; i++)
    {
        for(int j=0; j<m; j++)
        {
        T[k][0] = (m+1)*i+j;
        T[k][1] = (m+1)*(i+1)+j;
        T[k][2] = (m+1)*(i+1)+j+1;
        k++;
        
        T[k][0] = (m+1)*i+j;
        T[k][1] = (m+1)*(i+1)+j+1;
        T[k][2] = (m+1)*i+j+1;
        k++;
        }
    }
    	int i=n-1;
    	for(int j=0; j<m; j++)
        {
        T[k][0] = (m+1)*i+j;
        T[k][1] = j;
        T[k][2] = j+1;
        k++;
        
        T[k][0] = (m+1)*i+j;
        T[k][1] = j+1;
        T[k][2] = (m+1)*i+j+1;
        k++;
        }
	}
    
    
    savetostl(X, Y, T, p, N);
    
    return 0;
}
