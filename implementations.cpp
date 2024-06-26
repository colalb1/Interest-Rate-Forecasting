// Importing packages
#include <iostream>
#include <iomanip>

#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <random>

// Declaring namespaces
using namespace std;

// Initializing standard normal distribution
normal_distribution<double> distribution(0., 1.);

// Model variables
double r_0 = 0, theta, kappa, sigma;    // mean level (the price the process reverts to), reversion rate, and volatility
double d_t = 0.01;                      // time increment

// Interest rate models

class Vasicek {
    private: 
        double B(const double& t, const double& T) {
            return (1 - exp(-kappa * (T - t))) / kappa;
        }

        double A(const double& t, const double& T) {
            double temp_B = B(t, T);
            return exp((theta - pow(sigma / kappa, 2) / 2) * (temp_B - T + t) - pow(sigma * temp_B, 2) / (4 * kappa));
        }

    public:
        double exact_value(const double& t, const double& T) {
            return exp(A(t, T) - r_0 * B(t, T));
        }

        vector<double> simulated_value(const int num_sims, const int num_time_steps, const double& T) {
            double d_t = T / num_time_steps;
            vector<double> rates(d_t, 0);

            for (int i = 1; i < num_time_steps + 1; i++) {
                rates[i] = rates[i - 1] + kappa * (theta - rates[i - 1]) * d_t + sigma * sqrt(d_t) * distribution(num_sims);
            }

            return rates;
        }
};


int main() {
    double r_0 = 0, theta, kappa, sigma;
    std::cin >> theta >> kappa >> sigma;
}