#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <windows.h>
#include <stdio.h>
#define PI 3.14159265

using namespace std;

struct dot              //struct of base dot
{
    long double x, y, z;
};

struct vertex : dot       // struct of vertex
{
    int i, j, k, m;
    
    long double RIMT(long double R, int m, int n, int i)    //R with indexes i, m, t
    {
        return i * sqrt((R * R) - (R * (m - 1) / n) * (R * (m - 1) / n)) / n;
    }

    long double Alpha(int i, int j) // angle 
    {
        if (i != 0)
        {
            return 90 * (j - 1) / i;
        }
        else
        {
            return 0;
        }
        
    }   

    vertex()
    {
        i = j = k = m = 0;
    }

    vertex(int i, int j, int k, int m, int n, long double Rad)
    {
        this->i = i;
        this->j = j;
        this->k = k;
        this->m = m;
		
        this->x = Rad*i* cos(Alpha(i, j) * PI / 180)/n;
        this->y = Rad*i * sin(Alpha(i, j) * PI / 180)/n;
        this->z = (m + k - 1) * Rad/n;
    }
};

struct segment          //struct of segment 
{
    vertex Vertexes[3];
};

struct tetrahetron
{
    vertex Vertexes[4];
    segment Segments[7]; // 1:12 2:13 3:14 4:23 5:24 6:34
};

struct Prysm
{
    int i, m, j1, j2;
    long double volume;
    vertex Vertexes[7];
    tetrahetron Tetrahedrons[4];
    
    Prysm()
    {
        i = m = j1 = j2 = 0;
    }

    Prysm(int i, int m, int j1, int j2, int n, long double Rad)
    {
        this->i = i;
        this->m = m;
        this->j1 = j1;
        this->j2 = j2;     
        
        if(i!=n-m+1)
        {
        this->Vertexes[1] = vertex(i, j2, 0, m, n, Rad);
        this->Vertexes[2] = vertex(i-1, j2-1, 0, m, n, Rad);
        this->Vertexes[3] = vertex(i-j1, j2 - 1 + j1, 0, m, n, Rad);
        this->Vertexes[4] = vertex(i, j2, 1, m, n, Rad);
        this->Vertexes[5] = vertex(i - 1, j2 - 1, 1, m, n, Rad);
        this->Vertexes[6] = vertex(i - j1, j2 - 1 + j1, 1, m, n, Rad);
    	}
    	else
    	{
			if(j1==1)//4???????=1
    		{
    			this->Vertexes[1] = vertex(i, j2, 0, m, n, Rad);
       			this->Vertexes[2] = vertex(i-1, j2-1, 0, m, n, Rad);
        		this->Vertexes[3] = vertex(i-j1, j2 - 1 + j1, 0, m, n, Rad);
        		this->Vertexes[4] = vertex(i, j2, 0, m, n, Rad);
        		this->Vertexes[5] = vertex(i - 1, j2 - 1, 1, m, n, Rad);
        		this->Vertexes[6] = vertex(i - j1, j2 - 1 + j1, 1, m, n, Rad);
			}
			if(j1==0)
			{
				this->Vertexes[1] = vertex(i, j2, 0, m, n, Rad);
       			this->Vertexes[2] = vertex(i-1, j2-1, 0, m, n, Rad);
        		this->Vertexes[3] = vertex(i-j1, j2 - 1 + j1, 0, m, n, Rad);
        		this->Vertexes[4] = vertex(i, j2, 0, m, n, Rad);
        		this->Vertexes[5] = vertex(i - 1, j2 - 1, 1, m, n, Rad);
        		this->Vertexes[6] = vertex(i - j1, j2 - 1 + j1, 0, m, n, Rad);
			}
		}

        for (int i = 1; i < 4; i++)
        {
            for (int j = 1; j <= 4; j++)
            {
                this->Tetrahedrons[i].Vertexes[j] = this->Vertexes[j + i - 1];
            }
        }

        for (int i = 1; i < 4; i++)
        {
            for (int j = 1; j <= 3; j++)
            {
                for (int m = 0; m <= 4 - j - 1; m++)
                {
                    this->Tetrahedrons[i].Segments[j * ((int)j / 2 + 1) + m].Vertexes[1] = this->Tetrahedrons[i].Vertexes[j];
                    this->Tetrahedrons[i].Segments[j * ((int)j / 2 + 1) + m].Vertexes[2] = this->Tetrahedrons[i].Vertexes[j + m + 1];
                }
            }
        }
        this->volume = (this->Vertexes[4].z - this->Vertexes[1].z) * abs((this->Vertexes[2].x - this->Vertexes[1].x) * (this->Vertexes[3].y - this->Vertexes[1].y) - (this->Vertexes[3].x - this->Vertexes[1].x) * (this->Vertexes[2].y - this->Vertexes[1].y)) / 2;
    } 
};

int main()
{
    int n, counter = 0, tempn = 0;
    double Rad;
    long double volume = 0;
 
    cout << "n: ";
    cin >> n;

    while (n > 600)
    {
        system("cls");
        cout << "n is too lagre!" << endl;
        cout << "n: ";
        cin >> n;
    }

    cout << "Rad: ";
    cin >> Rad;

    tempn = n * n * n;
    Prysm* prysms = new Prysm[tempn];

    for (int m = 1; m <= n; m++)
    {
        for (int i = 1 ; i <= n-m+1; i++)
        {
            for (int j1 = 0; j1 <= 1; j1++)
            {
                for (int j2 = 2; j2 <= i - j1 + 1; j2++)
                {
                    prysms[counter] = Prysm(i, m, j1, j2, n, Rad);
                    volume += prysms[counter].volume;
                    counter++;
                }
            }
        }
    }

    volume *= 8;
    cout.precision(20);
    cout << "sum=" << volume << endl << "V=" << 4 * PI * Rad * Rad * Rad / 3 << endl << "% =" << (volume / (4 * PI * Rad * Rad * Rad / 3) - 1) * 100 << " runtime = " << clock() / 1000.0 << " sec" << endl;

    ofstream file("approx.obj");

    if (!file)
    {
        cout << "FILE DOESN'T CREATED!" << endl;
    }
    else
    {
        cout << "FILE SUCCESFULLY CREATED" << endl;

        for (int i = 0; i < counter; i++)
        {
        	//if(prysms[i].Vertexes[4].z!=prysms[i].Vertexes[1].z)
        	///
            	for (int j = 1; j < 7; j++)
            	{
              		//  file.setf(ios_base::fixed);
                	file << "v " << prysms[i].Vertexes[j].x << " " << prysms[i].Vertexes[j].y << " " << prysms[i].Vertexes[j].z << endl;              
            	}
            	file << "f " << 6 * i + 1 << " " << 6 * i + 2 << " " << 6 * i + 3 << endl << "f " << 6 * i + 4 << " " << 6 * i + 5 << " " << 6 * i + 6 << endl << "f " << 6 * i + 1 << " " << 6 * i + 2 << " " << 6 * i + 5 << " " << 6 * i + 4 << endl << "f " << 6 * i + 2 << " " << 6 * i + 3 << " " << 6 * i + 6 << " " << 6 * i + 5 << endl << "f " << 6 * i + 1 << " " << 6 * i + 3 << " " << 6 * i + 6 << " " << 6 * i + 4 << endl;
        	//}
        	/*else
        	{
        		if (prysms[i].Vertexes[6].z!=prysms[i].Vertexes[3].z)//??????????? ???????? ?? ??????
        		{
        			for(int j = 1; j < 7; j++)
        			{
        				file << "v " << prysms[i].Vertexes[j].x << " " << prysms[i].Vertexes[j].y << " " << prysms[i].Vertexes[j].z << endl;
					}
					file << "f " << 6 * i + 1 << " " << 6 * i + 2 << " " << 6 * i + 3 << endl << "f " << 6 * i + 1 << " " << 6 * i + 5 << " " << 6 * i + 6 << endl << "f " << 6 * i + 2 << " " << 6 * i + 3 << " " << 6 * i + 6 << " " << 6 * i + 5 << endl << "f " << 6 * i + 1 << " " << 6 * i + 2 << " " << 6 * i +5 << endl << "f " << 6 * i + 1 << " " << 6 * i +3 << " " << 6 * i +6 << endl ;
				}
				else//??????????? ???????? ? ??????
				{
					for(int j = 1; j < 7; j++)
        			{
        				file << "v " << prysms[i].Vertexes[j].x << " " << prysms[i].Vertexes[j].y << " " << prysms[i].Vertexes[j].z << endl;
					}
					file << "f " << 6 * i + 1 << " " << 6 * i + 2 << " " << 6 * i + 3 << endl << "f " << 6 * i + 1 << " " << 6 * i + 2 << " " << 6 * i + 5 << endl << "f " << 6 * i + 1 << " " << 6 * i + 3 << " " << 6 * i + 5 << endl << "f " << 6 * i + 3 << " " << 6 * i + 2 << " " << 6 * i + 5 << endl;
				}
			}*/
		}       
        cout << "FILE SUCCESFULLY WRITTEN" << endl;
        cout << clock()/1000.0 << " sec write-time" << endl;

        file.close();
        
        system("explorer approx.obj");
    }

    delete[] prysms;
    return 0;
}