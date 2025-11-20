#include <iostream>
#include <thread>
#include <random>
#include <numeric>
#include <sstream>
#include <fstream>
#include <vector>
#include <ctime>
#include <cmath>
#include <unordered_map>
#include <utility>
#include <algorithm>

double J = 1.0;

double Magnetization(const std::vector<int>&spinsm)
{
	int n=spinsm.size();
	double m = 0.0;
	for(int i=0;i<n;i++)
	{
		m = m+spinsm[i];
	}
	return (m/double(n));
}

double Energy(const std::vector<int>&spinsf,const std::vector<int>&oper_list, const double &h, const std::vector<int>&up)
{
        int l = spinsf.size();
        int p = oper_list.size();
        std::vector<int> spinsm(l);
        double E =0.0;
        for(int i=0; i<l;i++)
        {
                spinsm[i] = spinsf[i];
        }
        for(int i =0; i<p/2;i++)
        {
                int site = oper_list[i]/3;
                if(oper_list[i]%3==2)
                {
                        spinsm[site] = -spinsm[site];
                }
        }
        for(int i=0; i<l; i++)
        {
                if(spinsm[i] ==spinsm[up[i]])
                {
                        E= E+2*J;
                }
                E=E+h;
        }

        return (E/double(l)-h-J);
}



void diagonal_update(std::vector<int>&spinsf,std::vector<int>&spinsm, std::vector<int>&spinsl, std::vector<int>&oper_list, std::uniform_real_distribution<double> &dis2,std::uniform_int_distribution<int>&dis3, const int &Ns, const int &Nb, std::mt19937 &gen, const std::vector<int>&up, const double&h)
{
	int L = spinsl.size();
	int p = oper_list.size();
	for(int i= 0;i<L;i++)
	{
		spinsl[i] =spinsf[i];
	}
	for(int i = 0; i<p; i++)
	{
		int site = oper_list[i]/3;
		if(oper_list[i]%3 ==2)
		{
			spinsl[site] = -spinsl[site];
		}
		else
		{
			double prob_accept = h*Ns/(h*Ns+2*J*Nb);
			while(1){
				int site2 = dis3(gen);
				if(dis2(gen)<prob_accept)
				{
					oper_list[i] = site2*3+1;
					break;
				}
				else
				{
					if(spinsl[site2]==spinsl[up[site2]])
					{
						oper_list[i] = site2*3+3;
						break;
					}
				}
			}
		}
		if(i == p/2-1)
		{
			for(int j = 0; j< L;j++)
			{
				spinsm[j] = spinsl[j];
			}
		}
	}
}

void link_update(const std::vector<int>&oper_list, std::vector<int>&X_list, std::vector<int>&V_first, std::vector<int>&V_last, const std::vector<int>&up)
{
	int N = up.size();
	int p = oper_list.size();
	for(int i=0;i<p;i++)
	{
		if(oper_list[i]%3==0)
		{
			int v0 = 4*i+N;
			int site = oper_list[i]/3-1;
			int i1 = site;
			int i2 = up[site];
			int v1 = V_last[i1];
			int v2 = V_last[i2];
			if(v1!=-1)
			{
				X_list[v1] = v0;
				X_list[v0] = v1;
			}
			else
			{
				V_first[i1]=v0;
			}
			if(v2!=-1)
			{
				X_list[v2] = v0+1;
				X_list[v0+1]=v2;

			}
			else
			{
				V_first[i2] = v0+1;
			}
			V_last[i1] = v0+2;
			V_last[i2] = v0+3;
		}
		else
		{
			int v0 = 4*i+N;
			int site =oper_list[i]/3;
			int v = V_last[site];
			if(v!=-1)
			{
				X_list[v] = v0;
				X_list[v0] = v;
			}
			else
			{
				V_first[site] = v0;
			}
			V_last[site] = v0+1;
		}
	}
	for(int i =0; i< N ;i++)
	{
		int vf = V_first[i];
		int vl = V_last[i];
		if(V_first[i]!=-1)
		{
			X_list[i] = vf;
			X_list[vf] = i;
		}
		if(V_last[i]!=-1)
		{
			X_list[i+4*p+N] = vl;
			X_list[vl] =i+4*p+N;
		}
		if(vf == -1 && vl == -1)
		{
			X_list[i] = 4*p+N+i;
			X_list[4*p+N+i] = i;
		}

	}

}

int find(int x,std::vector<int>&clusters)
{
	if(clusters[x]==x)
	{
		return x;
	}
	return clusters[x] = find(clusters[x],clusters);
}

void unite(int x, int y, std::vector<int>&clusters)
{
	int rx = find(x, clusters);
	int ry = find(y, clusters);
	if(rx !=ry)
	{
		clusters[ry] =rx;
	}
}

void off_diagonal_update(std::vector<int>&X_list, std::vector<int>&V_first,std::vector<int>&V_last,  std::vector<int>&oper_list, std::uniform_real_distribution<double> &dis2, std::mt19937 &gen,std::vector<int>&spinsf, std::vector<int>&spinsl)
{
	int p = oper_list.size();
	int n = V_first.size();
	
	std::vector<int> clusters(p*4+2*n,-1);
	for(int i=0;i<4*p+2*n;i++)
	{
		if(X_list[i]!=-1)
		{
			clusters[i]=i;
		}
	}
	for(int i=0;i<n;i++)
	{
		if(X_list[i]!=-1)
		{
			unite(i,X_list[i],clusters);
		}
		if(X_list[i+4*p+n]!=-1)
		{
			unite(i+4*p+n, X_list[i+4*p+n],clusters);
		}
	}
	for(int i=0;i<p;i++)
	{
		if(oper_list[i]%3==0)
		{
			int v0 = 4*i+n;
			int v1 = v0+1;
			int v2 = v0+2;
			int v3 = v0+3;
			unite(v0,v1,clusters);
			unite(v1,v2,clusters);
			unite(v2,v3,clusters);
			unite(v0,X_list[v0],clusters);
			unite(v1,X_list[v1],clusters);
			unite(v2,X_list[v2],clusters);
			unite(v3,X_list[v3],clusters);
		}
		else
		{
			int v0 = 4*i+n;
			int v1 = v0+1;
			unite(v0,X_list[v0],clusters);
			unite(v1,X_list[v1],clusters);
		}
	}

	for(int i = 0; i<4*p+2*n; i++)
	{
		if(clusters[i]!=-1)
		{
			clusters[i] = find(i, clusters);
		}
	}


	std::unordered_map<int,bool>flip_clusters;

	for(int i=0; i<4*p+2*n; i++)
	{
		if(clusters[i] == -1)
		{
			continue;
		}
		int cluster_index =clusters[i];
		if(flip_clusters.find(cluster_index)==flip_clusters.end())
		{
			flip_clusters.insert({cluster_index,(dis2(gen)<0.5)});
		}
		if(flip_clusters[cluster_index])
		{
			clusters[i] = -2;
		}
	}
	
	for(int i=0; i<p; i++)
	{
		if(oper_list[i]%3==1)
		{	
			if(clusters[4*i+n]==-2)
				oper_list[i]++;
			if(clusters[4*i+1+n]==-2)
				oper_list[i]++;
			if(oper_list[i]%3==0)
				oper_list[i]=oper_list[i]-2;
		}
		else if(oper_list[i]%3==2)
		{
			if(clusters[4*i+n]==-2)
				oper_list[i]--;
			if(clusters[4*i+1+n]==-2)
				oper_list[i]--;
			if(oper_list[i]%3==0)
				oper_list[i]=oper_list[i]+2;
		}
	}

	for(int i=0;i<n;i++)
	{
		if(clusters[i]==-2)
		{
			spinsf[i]=-spinsf[i];
		}
	}
	for(int i=4*p+n;i<4*p+2*n;i++)
	{
		if(clusters[i]==-2)
                {
                        spinsl[i-4*p-n]=-spinsl[i-4*p-n];
                }

	}
}

std::pair<double,double> Measure2(const std::vector<double>&Magnet, const int &binning)
{
        int n = Magnet.size();
        int l = n/binning;
        std::vector<double> Re_M(l);
        for(int i=0; i<n; i++)
        {
                Re_M[i/binning] = Re_M[i/binning] + Magnet[i]/(double)binning;
        }
        double M_aver = 0.0;
        double M_aver2 = 0.0;
        for(int i=0; i<l; i++)
        {
                M_aver = M_aver + Re_M[i]/l;
                M_aver2 = M_aver2 + Re_M[i]*Re_M[i]/l;
        }
        double V = M_aver2-M_aver*M_aver;
        V = V/(double)(l-1);

        return {M_aver, std::sqrt(V)};
}


int main()
{
	//variable
	int L = 8;
	int sweep = 1048576;
	int tsweep = 1000;
	int binning = 8;
	int p = 256;
	double h = 1.0;

	clock_t start = std::clock();

	std::random_device rd;
	std::mt19937 gen(rd());

	std::vector<int> spinsf(L);
	std::vector<int> spinsl(L);
	std::vector<int> spinsm(L);
	std::uniform_int_distribution<int> dis(0,1);
	std::uniform_real_distribution<double> dis2(0.0,1.0);
	std::uniform_int_distribution<int> dis3(0,L-1);

	for(int i=0; i<L; i++)
	{
		int number = dis(gen) == 0 ? -1:1;
		spinsf[i] = number;
		spinsm[i] = number;
		spinsl[i] = number;
	}
	//Ns is the number of sites and Nb is the number of bonds 
	int Ns = L;
	int Nb = L;

	std::vector<int> up(L), down(L);
	for(int i=0; i<L;i++)
	{
		up[i] = (i-1+L)%L;
		down[i] = (i+1+L)%L;
	}

	std::ostringstream name;
	name << "RESULT_SSE_cluster_p_ground_:"<<int(h) << "_" <<L;
	std::ofstream out(name.str());
	for(int j= 0; j<1; j++)
	{
		p = p*2;
		std::vector<int> oper_list(p,4);
		std::vector<double> Magnet2(sweep);
		std::vector<double> logw(sweep);
		std::vector<double> energy(sweep);
		for(int i=0; i<tsweep; i++)
		{
			diagonal_update(spinsf,spinsm, spinsl, oper_list,dis2,dis3,Ns,Nb,gen,up,h);
			std::vector<int> X_list(4*p+2*L,-1);
			std::vector<int> V_first(L,-1);
			std::vector<int> V_last(L,-1);
			link_update(oper_list, X_list, V_first, V_last, up);
			off_diagonal_update(X_list, V_first,V_last, oper_list,dis2,gen,spinsf, spinsl);
		}

		for(int i=0; i<sweep; i++)
		{
			energy[i] = Energy(spinsf,oper_list,h,up);
			diagonal_update(spinsf,spinsm, spinsl, oper_list, dis2, dis3, Ns,Nb, gen,up,h);
                        Magnet2[i]=std::abs(Magnetization(spinsm));
			std::vector<int> X_list(4*p+2*L,-1);
			std::vector<int> V_first(L,-1);
			std::vector<int> V_last(L,-1);
			link_update(oper_list, X_list, V_first, V_last, up);
			off_diagonal_update(X_list, V_first,V_last, oper_list,dis2,gen,spinsf, spinsl);
			//std::this_thread::sleep_for(std::chrono::seconds(1));
		}
		double EM=0.0;
		double MM=0.0;
		for (int i=0; i<sweep; i++)
	       	{
			EM+=energy[i]/double(sweep);
			MM+=Magnet2[i]/double(sweep);
		}

		/*for(int i= 0;i <sweep; i++)
		{
			out << i << "\t" << Energy[i] << "\n";
		}*/

		out << p << "\t" << MM  <<"\t"<< EM<< "\n";
	}
	

	out.close();
		

	clock_t end = std::clock();

	double time_taken = double(end-start)/CLOCKS_PER_SEC;
	std::cout<< time_taken << std::endl;

	return 0;
}
