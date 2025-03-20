#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#define PI 3.14159265

using namespace std;


int const n=10, m=10, p=1+n*m*(m+1)/2, N=n*m*m;
double X[p], Y[p];
int T[N][3];

int NumbPoint(int i,int j,int n)
{
int Nij;
	if(i==0 && j==0)
	{
		Nij=0;
	}
	else
	{
		Nij=n*j*(j-1)/2+i+1;
	}
return Nij;
}

double ro (double x)
{
	double a=1;
	double b=2;
	return a*b/sqrt(a*a*sin(x)*sin(x)+b*b*cos(x)*cos(x));
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
    double phiij;
    int k = 0;
    
    X[0]=0;
	Y[0]=0;
	k++;
    
    for(int j=1; j<=m; j++)
    {
        for(int i=0; i<=n*j-1; i++)
        {
        	phiij=2*PI*i/(j*n);
        	X[k]=j*ro(phiij)*cos(phiij)/m;
        	Y[k]=j*ro(phiij)*sin(phiij)/m;
        	k++;
        }
    }
    
    
    
    
    k=0;
    for(int i=0; i<=n-2; i++)//первые n треугольников
    {
        T[k][0] = NumbPoint(0,0,n);
        T[k][1] = NumbPoint(i,1,n);
        T[k][2] = NumbPoint(i+1,1,n);
        k++;
    }
    
    T[k][0] = NumbPoint(0,0,n);
    T[k][1] = NumbPoint(n-1,1,n);
	T[k][2] = NumbPoint(0,1,n);
    k++;
    
    for(int s=1; s<=n-1 ;s++)//n секторов
    {
    	for(int j=1; j<=m-1 ;j++)
    	{
    		for(int i=(s-1)*j; i<=s*j;i++)
    		{
    			if(i!=s*j)
    			{
    				T[k][0] = NumbPoint(i,j,n);
        			T[k][1] = NumbPoint(i+s-1,j+1,n);
        			T[k][2] = NumbPoint(i+s,j+1,n);
       				k++;
       				
       				T[k][0] = NumbPoint(i,j,n);
        			T[k][1] = NumbPoint(i+s,j+1,n);
        			T[k][2] = NumbPoint(i+1,j,n);
       				k++;
				}
				else
				{
					T[k][0] = NumbPoint(i,j,n);
        			T[k][1] = NumbPoint(i+s-1,j+1,n);
        			T[k][2] = NumbPoint(i+s,j+1,n);
       				k++;
				}
			}
		}
	}
	
	int s=n;
	for(int j=1; j<=m-1 ;j++)
    {
    	for(int i=(s-1)*j; i<=s*j;i++)
    	{
    		if(i<s*j-1)
    		{
    			T[k][0] = NumbPoint(i,j,n);
        		T[k][1] = NumbPoint(i+n-1,j+1,n);
        		T[k][2] = NumbPoint(i+n,j+1,n);
       			k++;
       			
       			T[k][0] = NumbPoint(i,j,n);
        		T[k][1] = NumbPoint(i+n,j+1,n);
        		T[k][2] = NumbPoint(i+1,j,n);
       			k++;
			}
			if(i==s*j-1)
			{
       			T[k][0] = NumbPoint(s*j-1,j,n);
        		T[k][1] = NumbPoint(s*(j+1)-2,j+1,n);
        		T[k][2] = NumbPoint(s*(j+1)-1,j+1,n);
       			k++;
       			
       			T[k][0] = NumbPoint(s*j-1,j,n);
        		T[k][1] = NumbPoint(s*(j+1)-1,j+1,n);
        		T[k][2] = NumbPoint(0,j,n);
       			k++;
			}
			if(i==s*j)
			{
				T[k][0] = NumbPoint(0,j,n);
        		T[k][1] = NumbPoint(s*(j+1)-1,j+1,n);
        		T[k][2] = NumbPoint(0,j+1,n);
       			k++;
			}
		}
	}
    
    
    savetostl(X, Y, T, p, N);
    
    return 0;
}
