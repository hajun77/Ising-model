#include <iostream>
#include <random>
#include <sstream>
#include <fstream>
#include <vector>
#include <ctime>
#include <cmath>
#include <unordered_map>
#include <utility>

double J = 1;

double Magnetization(const std::vector<int>&spins)
{
	int n=spins.size();
	double m = 0.0;
	for(int i=0;i<n;i++)
	{
		m = m+spins[i];
	}
	return m/n;
}

double diagonal_update(std::vector<int>&spins, std::vector<int>&oper_list, const double &beta, std::uniform_real_distribution<double> &dis2,std::uniform_int_distribution<int>&dis3, const int &Ns, const int &Nb, std::mt19937 &gen, const std::vector<int>&up, const double&h)
{
	int L = spins.size();
	int p = oper_list.size();
	int n_inv = 0;
	for(int i = 0; i<p; i++)
	{
		if(oper_list[i]==0)
		{
			n_inv++;
		}
	}
	for(int i = 0; i<p; i++)
	{
		double prob_accept = 0.0;
		if(oper_list[i]==0)
		{
			if(n_inv==0)
			{
				prob_accept = 1;
			}
			else
			{
				prob_accept = beta*(h*Ns+2*J*Nb)/double(n_inv);
			}
			if(dis2(gen)<prob_accept)
			{
				double prob_hj = (h*Ns)/double(h*Ns+2*J*Nb);
				int site = dis3(gen);
				if(dis2(gen) < prob_hj)
				{
					oper_list[i] = site*3+1;
					n_inv--;
				}
				else
				{
					if(spins[site]==spins[up[site]])
					{
						oper_list[i] = site*3+3;
						n_inv--;
					}
				}
			}
		}
		else if(oper_list[i]%3==0||oper_list[i]%3==1)
		{
			double prob_reject = (n_inv+1)/(beta*(h*Ns+2*J*Nb));
			if(dis2(gen)<prob_reject)
			{
				oper_list[i] = 0;
				n_inv++;
			}
		}
		else if(oper_list[i]%3==2)
		{
			int site = oper_list[i]/3;
			spins[site] = -spins[site];
		}
	}
	return (double(p-n_inv)/(beta*L));
}

void link_update(const std::vector<int>&oper_list, std::vector<int>&X_list, std::vector<int>&V_first, std::vector<int>&V_last, const std::vector<int>&up)
{
	int N = up.size();
	int p = oper_list.size();
	for(int i=0;i<p;i++)
	{
		if(oper_list[i]==0)
		{
			continue;
		}
		else if(oper_list[i]%3==0)
		{
			int v0 = 4*i;
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
			int v0 = 4*i;
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

	for(int i =0; i<N; i++)
	{
		int f = V_first[i];
		if(f!=-1)
		{
			int l =V_last[i];
			X_list[f] = l;
			X_list[l] = f;
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

void off_diagonal_update(std::vector<int>&X_list, std::vector<int>&V_first, std::vector<int>&oper_list, std::uniform_real_distribution<double> &dis2, std::mt19937 &gen, std::vector<int>&spins)
{
	int p = oper_list.size();
	int n = V_first.size();
	
	std::vector<int> clusters(p*4,-1);
	for(int i=0;i<4*p;i++)
	{
		if(X_list[i]!=-1)
		{
			clusters[i]=i;
		}
	}
	for(int i=0;i<p;i++)
	{
		if(oper_list[i]==0)
		{
			continue;
		}
		else if(oper_list[i]%3==0)
		{
			int v0 = 4*i;
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
			int v0 = 4*i;
			int v1 = v0+1;
			unite(v0,X_list[v0],clusters);
			unite(v1,X_list[v1],clusters);
		}
	}

	for(int i = 0; i<4*p; i++)
	{
		if(X_list[i]!=-1)
		{
			clusters[i] = find(i, clusters);
		}
	}


	std::unordered_map<int,bool>flip_clusters;

	for(int i=0; i<4*p; i++)
	{
		if(X_list[i] == -1)
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
		if(oper_list[i]==0)
		{
			continue;
		}
		else if(oper_list[i]%3==1)
		{	
			if(clusters[4*i]==-2)
				oper_list[i]++;
			if(clusters[4*i+1]==-2)
				oper_list[i]++;
			if(oper_list[i]%3==0)
				oper_list[i]=oper_list[i]-2;
		}
		else if(oper_list[i]%3==2)
		{
			if(clusters[4*i]==-2)
				oper_list[i]--;
			if(clusters[4*i+1]==-2)
				oper_list[i]--;
			if(oper_list[i]%3==0)
				oper_list[i]=oper_list[i]+2;
		}
	}

	for(int i=0;i<n;i++)
	{
		if(V_first[i]==-1)
		{
			if(dis2(gen)<0.5)
			{
				spins[i]=-spins[i];
			}
		}
		else if(clusters[V_first[i]]==-2)
		{
			spins[i]=-spins[i];
		}
	}

}

std::pair<double,double> Measure(const std::vector<double>&Energy, const int &binning)
{
	int n = Energy.size();
	int l = n/binning;
	std::vector<double> Re_E(l);
	for(int i=0; i<n; i++)
	{
		Re_E[i/binning] = Re_E[i/binning] + Energy[i]/(double)binning;
	}
	double E_aver = 0.0;
	double E_aver2 = 0.0;
	for(int i=0; i<l; i++)
	{
		E_aver = E_aver + Re_E[i]/l;
		E_aver2 = E_aver2 + Re_E[i]*Re_E[i]/l;
	}
	double V = E_aver2-E_aver*E_aver;
	V = V/(double)(l-1);
	
	return {E_aver, std::sqrt(V)};
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
	int L = 32;
	int sweep = 262144;
	int tsweep = 100;
	int binning = 8;
	int p = 32;
	double h = 1.0;

	clock_t start = std::clock();

	double T = 0.05;
	double beta = 1.0/T;

	std::random_device rd;
	std::mt19937 gen(rd());

	std::vector<int> spins(L);
	std::uniform_int_distribution<int> dis(0,1);
	std::uniform_real_distribution<double> dis2(0.0,1.0);
	std::uniform_int_distribution<int> dis3(0,L-1);

	for(int i=0; i<L; i++)
	{
		spins[i] = dis(gen) == 0 ? -1:1;
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
	name << "RESULT_SSE_cluster_p_beta:" << int(beta) << "_" <<L;
	std::ofstream out(name.str());
	//out <<"h" << "\t" << "Energy" <<"\t"<< "error" <<"\t"<< "Magnet" << "\t" << "error" << "\n";
	for(int j= 0; j<14; j++)
	{
		std::vector<int>vertex_link(4*p);
		std::vector<int> oper_list(p);
		std::vector<double> Energy(sweep);
		std::vector<double> Magnet(sweep);
		//h = h+0.01;
		p = p * 2
		for(int i=0; i<tsweep; i++)
		{
			double T_E = diagonal_update(spins, oper_list,beta,dis2,dis3,Ns,Nb,gen,up,h);
			std::vector<int> X_list(4*p,-1);
			std::vector<int> V_first(L,-1);
			std::vector<int> V_last(L,-1);
			link_update(oper_list, X_list, V_first, V_last, up);
			off_diagonal_update(X_list, V_first, oper_list,dis2,gen,spins);
		}

		for(int i=0; i<sweep; i++)
		{
			Energy[i]=diagonal_update(spins, oper_list, beta, dis2, dis3, Ns,Nb, gen,up,h);
			std::vector<int> X_list(4*p,-1);
			std::vector<int> V_first(L,-1);
			std::vector<int> V_last(L,-1);
			link_update(oper_list, X_list, V_first, V_last, up);
			off_diagonal_update(X_list, V_first, oper_list,dis2,gen,spins);
			Magnet[i]=std::abs(Magnetization(spins));
		}
		
		auto [E_M,E_E] = Measure(Energy, binning);
		auto [M_M,M_ME] = Measure2(Magnet, binning);
		/*for(int i= 0;i <sweep; i++)
		{
			out << i << "\t" << Energy[i] << "\n";
		}*/

		out << h << "\t" << E_M << "\t" << E_E <<"\t" << M_M << "\t" << M_ME << "\n";
	}
	

	out.close();
		

	clock_t end = std::clock();

	double time_taken = double(end-start)/CLOCKS_PER_SEC;
	std::cout<< time_taken << std::endl;

	return 0;
}
