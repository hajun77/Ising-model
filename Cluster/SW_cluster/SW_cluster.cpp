#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <numeric>
#include <fstream>
#include <ctime>
#include <unordered_map>
#include <sstream>

void Save_spin(const std::vector<int>&spins, const std::string&filename)
{
	std::ofstream out(filename);
	for (int s : spins){
		out << s << " ";
	}
	out.close();

}

std::vector<int> Load_spin(const std::string&filename, int N)
{
	std::vector<int> spins(N);
	std::ifstream in(filename);
	for(int i=0;i<N;i++){
		in >> spins[i];
	}
	in.close();
	return spins;
}


double Energy(const std::vector<int>&spins, const std::vector<int>&up, const std::vector<int>&down, const std::vector<int>&left, const std::vector<int>&right)
{
	double energy=0.0;
	int N=spins.size();

	for(int i=0; i<N; i++)
	{
		int s_i=spins[i];
		int neighbor_sum=spins[up[i]]+spins[down[i]]+spins[left[i]]+spins[right[i]];
		energy-=s_i*neighbor_sum;
	}

	return energy/(2*N);
}

double Magnetization(const std::vector<int>&spins)
{
	double M=0;
	int N=spins.size();
	for(int i=0;i<N;i++){
		M+=spins[i];
	}
	return M/double(N);
}

int find(int x, std::vector<int>&clusters) {
	if (clusters[x] == x) return x;
	return clusters[x] = find(clusters[x], clusters);
}

void unite(int x, int y, std::vector<int>&clusters) {
	int rx = find(x, clusters);
	int ry = find(y, clusters);
	if (rx != ry) {
		clusters[ry] = rx;
	}
}

void cluster_update(std::vector<int>&spins,  const std::vector<int>&down, const std::vector<int>&right, double prob, std::mt19937 &gen, std::uniform_real_distribution<double> &dis2)
{
	int N=spins.size();

	std::vector<int> clusters(N);
	for(int i=0;i<N;i++) {
		clusters[i]=i;
	}

	for(int i=0;i<N;i++) {
		if (spins[i]==spins[right[i]]&& dis2(gen) < prob) {
			unite(i, right[i], clusters);
		}
		if (spins[i]==spins[down[i]] && dis2(gen) < prob) {
			unite(i, down[i], clusters);
		}

	}
	for (int i = 0; i < N; i++)
		clusters[i] = find(i, clusters);
	
	std::unordered_map <int, bool>flip_cluster;

	for( int i = 0; i < N;i++)
	{
		int cluster_index = clusters[i];
		if (flip_cluster.find(cluster_index)==flip_cluster.end())
			flip_cluster.insert({cluster_index, (dis2(gen)<0.5)});
		if(flip_cluster[cluster_index])
			spins[i]=-spins[i];
	}
	
}

double binning(std::vector<double>data)
{
	int N = data.size();
	double squ_ava0 = 0;
	double ava0 = 0;
	for(int i=0;i<N;i++){
		squ_ava0 = squ_ava0+data[i]*data[i];
		ava0 = ava0+data[i];
	}
	squ_ava0 = squ_ava0/N;
	ava0 = ava0/N;
	ava0 = ava0*ava0;
	double delta0 = (squ_ava0-ava0)/(N-1);
	while(data.size() > 128){
		std::vector<double>current(data.size()/2);
		for(int j=0; j<current.size(); j++)
		{
			current[j] = (data[2*j]+data[2*j+1])/2;
		}

		data = current;
	}
	N = data.size();
	double squ_aval = 0;
	double aval = 0;
	for(int i=0;i<N;i++)
	{
		squ_aval = squ_aval+data[i]*data[i];
		aval = aval + data[i];
	}
	squ_aval = squ_aval/N;
	aval = aval/N;
	aval =aval * aval;
	double deltal = (squ_aval-aval)/(N-1);
	return 0.5*(deltal/delta0-1);
}

int main()
{
	int L=100;
	int sweeps=524288;
	int tsweeps=100;
	int interval=0;
	int N=L*L;

	clock_t start = std::clock();

	double T=2.27;
	double beta= 1.0/T;

	std::random_device rd;
	std::mt19937 gen(rd());

	std::vector<int> spins(N);
	std::uniform_int_distribution dis(0,1);
	std::uniform_real_distribution<double> dis2(0.0,1.0);
	//std::ifstream test("last_spin.txt");
	/*if (test)
	{
		spins = Load_spin("last_spin.txt",N);
	}
	else
	{
		for(int i=0; i<N; i++){
			spins[i]=dis(gen) == 0 ? -1:1 ;
		}
	}*/
	for(int i=0; i<N;i++)
	{
		spins[i] = dis(gen) == 0 ? -1:1;
	}
	std::vector<int> up(N), down(N), left(N), right(N);
	for(int i=0;i<N;i++){
		up[i] = (i-L+N) % N;
		down[i] = (i+L) % N;
		left[i] = (i % L == 0) ? i + (L-1) : i-1;
		right[i] = ((i+1) % L == 0) ? i - (L-1) : i +1;	
	}
	
	double prob = 1-exp(-2*beta);

	std::vector<double> tau(100);
	for(int j=0;j<tsweeps;j++)
	{
		cluster_update(spins, down, right, prob, gen,dis2);
	}
	for(int j=0;j<100;j++){
		//std::vector<double> E(sweeps);
		std::vector<double> M(sweeps);
		for(int i=0;i<sweeps;i++) {
			cluster_update(spins, down, right, prob, gen, dis2);
			//E[i]=Energy(spins, up, down, left, right);
			M[i]=std::abs(Magnetization(spins));
		}
		tau[j] = binning(M);
	}

	std::ostringstream name;
	name<<"RESULT_SW_Cluster_"<< int(T*100) << "_" << L;
	std::ofstream out(name.str());
	for(int i=0;i< 100;i++)
	{
		out << "tau:" << "\t" << tau[i] << "\n"; 
	}

	out.close();

	//Save_spin(spins,"last_spin.txt");

	std::cout << "# escape time : " << (double)(std::clock()-start)/CLOCKS_PER_SEC << std::endl;
	
	return 0;
}	
