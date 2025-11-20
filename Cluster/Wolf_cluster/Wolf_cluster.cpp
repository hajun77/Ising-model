#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <ctime>

using namespace std;

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
	return M/N;
}

void grow_cluster(std::vector<int>&spins, const std::vector<int>&up, const std::vector<int>&down, const std::vector<int>&left, const std::vector<int>&right, int loc, int state, double &prob, std::mt19937&gen, std::uniform_real_distribution<double> &dis2, double &C_number)
{

	int loc_up = up[loc];
	int loc_down = down[loc];
	int loc_left = left[loc];
        int loc_right = right[loc];

	spins[loc] = -spins[loc];

	C_number++;

	if (state == spins[loc_up] && dis2(gen) < prob)
	{
		grow_cluster(spins, up, down, left, right, loc_up, state, prob,gen, dis2, C_number);
	}

	if (state == spins[loc_down] && dis2(gen) < prob)
	{
		grow_cluster(spins, up, down, left, right, loc_down, state, prob,gen, dis2, C_number);
	}

	if (state == spins[loc_left] && dis2(gen)<prob)
	{
		grow_cluster(spins, up, down, left, right, loc_left, state, prob,gen, dis2, C_number);
	}

	if (state == spins[loc_right] && dis2(gen)<prob)
	{
		grow_cluster(spins, up, down, left, right, loc_right, state, prob,gen,dis2, C_number);
	}

}

void update_spin(std::vector<int>&spins, const std::vector<int>&up, const std::vector<int>&down, const std::vector<int>&left, const std::vector<int>&right, double &prob, std::mt19937&gen,std::uniform_real_distribution<double> &dis2, double&C_number)
{
	int N=spins.size();
	

	std::uniform_int_distribution dis(0,N-1);
	int chosen_spin = dis(gen);
	int state = spins[chosen_spin];

	grow_cluster(spins, up, down, left, right, chosen_spin, state, prob,gen,dis2,C_number);

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
	aval = aval * aval;
        double deltal = (squ_aval-aval)/(N-1);
        return (0.5)*(deltal/delta0-1);
}

	
	

int main()
{
	int L=100;
	int sweeps=1048576;
	int tsweeps=1000;
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
	for(int i=0; i<N; i++){
		spins[i]=dis(gen) == 0 ? -1:1 ;
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
	for(int i=0;i<tsweeps;i++){
		double C_number =0;
		update_spin(spins, up, down, left, right, prob,gen,dis2,C_number);
	}
	for(int j=0;j<100;j++){
		//std::vector<double> E(sweeps);
		std::vector<double> M(sweeps);
		double C_N=0;

		for(int i=0; i<sweeps; i++){
			double C_number=0;
			update_spin(spins, up, down, left, right, prob,gen,dis2,C_number);
			C_N =C_N+C_number/sweeps;
			//E[i] = Energy(spins, up, down, left, right);
			M[i] = std::abs(Magnetization(spins));
		}
		C_N = C_N/N;
		tau[j] = binning(M)*C_N;
	}

	for(int i=0;i<100;i++)
	{
		std::cout <<"tau" << "\t" << tau[i] << "\n";
	}

	std::cout << "# escape time :" << (double)(std::clock()-start)/CLOCKS_PER_SEC << std::endl;
}	
