// Importing packages
#include <iostream>
#include <iomanip>

#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <functional>
#include <cassert>

// Initializing standard normal distribution
std::random_device rd;
std::mt19937 gen(rd());
std::normal_distribution<double> dist(0.0, 1.0);

// Loop blocking size
const int BLOCK_SIZE = 64;

// There is no general constructor in this file because they all have unique inputs.
// Consider putting each model into unique cpp file and putting testing environment
// inside.
class BlackKarasinski {
    private:
        // Theta depends on time
        double r_0, kappa, sigma;

        // Uncomment based on desired theta model
        double get_theta(const double& time, const double& r_init) {
            // return r_init; // Constant
            return r_init + 0.001 * time; // Linear with time
            // return r_init * exp(0.001 * time); // Exponential
        }

    public:
        BlackKarasinski(double r_0_con, double kappa_con, double sigma_con) {
            r_0 = r_0_con;
            kappa = kappa_con;
            sigma = sigma_con;
        }

        std::vector<double> simulated_value(const double& num_time_steps, const double& T) {
            double d_t = T / num_time_steps;

            std::vector<double> rates(num_time_steps, 0);
            rates[0] = r_0;

            for (int i = 1; i < num_time_steps; i += BLOCK_SIZE) {
                for (int j = i; j < std::min<int>(i + BLOCK_SIZE, num_time_steps); ++j) {
                    double drift = rates[j - 1] * (get_theta(j * d_t, r_0) + std::pow(sigma, 2) / 2 - kappa * std::log(rates[j - 1]));
                    double diffusion = sigma * rates[j - 1] * sqrt(d_t) * dist(gen);

                    rates[j] = rates[j - 1] + drift * d_t + diffusion;
                }
            }

            return rates;
        }
};

int main() {
    // BlackKarasinski bk_testing_class(0.05, 0.2, 0.02);
    // std::vector<double> bk_temp = bk_testing_class.simulated_value(10000, 1);
    // std::cout << "Black-Karasinski simulated rate: " << bk_temp[bk_temp.size() - 1] << std::endl;

    return 0;
}