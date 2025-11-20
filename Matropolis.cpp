#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <ctime>
double Magnetization(const std::vector<int> spins) {
    int M=0;
    int N = spins.size();
    for (int i = 0; i < N; i++) {
        M+=spins[i];
    }
    return static_cast<double>(M) / N;
}

double Energy(const std::vector<int>& spins,
                      const std::vector<int>& up,
                      const std::vector<int>& down,
                      const std::vector<int>& left,
                      const std::vector<int>& right
                      ) {
    double energy = 0.0;
    int N = spins.size();

    for (int i = 0; i < N; ++i) {
        int s_i = spins[i];
        int neighbor_sum = spins[up[i]] + spins[down[i]] + spins[left[i]] + spins[right[i]];
        energy += -s_i * neighbor_sum;
    }
    return energy /(2*N);
}
void Metropolis( std::vector<int>& spins,
                    const std::vector<int>& up,
                    const std::vector<int>& down,
                    const std::vector<int>& left,
                    const std::vector<int>& right,
                    const std::vector<double>& exp_table,
                    std::mt19937& gen) {
    int N = spins.size();
    std::uniform_real_distribution<> dist(0.0, 1.0);

    for (int i = 0; i < N; ++i) {
        int s = spins[i];
        int neighbor_sum = spins[up[i]] + spins[down[i]] + spins[left[i]] + spins[right[i]];
        int deltaE = 2*s * neighbor_sum;

        if (deltaE <= 0) {
            spins[i] *= -1;
        }
        else if ( dist(gen) < exp_table[deltaE + 8])  {
            spins[i] *= -1;
        }
    }
}

int main() {
    int NX = 100;
    int NY = 100;
    int steps = 1048576;
    int tsteps = 0;
    int output_interval = 0;
    int N = NX * NY;

    clock_t start = std::clock();

    std::random_device rd;
    std::mt19937 gen(rd());

    for (int ti = 227; ti <= 227; ++ti) {
        double T = ti / 100.0;
        double beta = 1.0 / T;

        std::vector<int> spins(N);
        std::uniform_int_distribution<> dis(0, 1);
        for (int i = 0; i < N; ++i) {
            spins[i] = dis(gen) == 0 ? -1 : 1;
        }

        std::vector<int> up(N), down(N), left(N), right(N);
        for (int idx = 0; idx < N; ++idx) {
            up[idx] = (idx - NY + N) % N;
            down[idx] = (idx + NY) % N;
            left[idx] = (idx % NY == 0) ? idx + (NY - 1) : idx - 1;
            right[idx] = (idx % NY == NY - 1) ? idx - (NY - 1) : idx + 1;
        }

        std::vector<double> exp_table(17, 0.0);
        for (int deltaE = -8; deltaE <= 8; deltaE += 4) {
            exp_table[deltaE + 8] = std::exp(-beta * deltaE);
        }

        std::ostringstream filename;
        filename << "ising_output_T" << static_cast<int>(ti) << "_100.txt";
        std::ofstream outfile(filename.str());

	std::vector<double> E(steps-tsteps);
	std::vector<double> M(steps-tsteps);

        for (int i = 0; i < steps; ++i) {
            Metropolis(spins, up, down, left, right, exp_table, gen);
            M[i] = Magnetization(spins);
            E[i] = Energy(spins, up, down, left, right);
	}
	for (int i=0;i<steps; i++){
		outfile << i << "\t" << E[i] <<"\t"<< M[i] << "\n";
	}
    }

    std::cout << "# escape time : " << (double)(std::clock()-start)/CLOCKS_PER_SEC << std::endl;

    return 0;
}
