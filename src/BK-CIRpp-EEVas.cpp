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
        auto r_0, kappa, sigma;
        std::vector<auto> theta;

        // Uncomment based on desired theta model
        double get_theta(auto time, auto r_init) {
            return r_init; // Constant
            // return r_init + 0.001 * time; // Linear with time
            // return r_init * exp(0.001 * time); // Exponential
        }

    public:
        BlackKarasinski(auto r_0_con, auto kappa_con, auto sigma_con, std::vector<auto> theta_con) {
            r_0 = r_0_con;
            kappa = kappa_con;
            sigma = sigma_con;
            theta = theta_con;
        }

        std::vector<auto> simulated_value(long int num_time_steps, double T) {
            auto d_t = T / num_time_steps;

            std::vector<auto> rates(num_time_steps, 0);
            rates[0] = r_0;

            for (int i = 1; i < num_time_steps; ++i) {
                auto drift = rates[i - 1] * (get_theta(i * d_t, r_0) + std::pow(sigma, 2) / 2 - kappa * std::log(rates[i - 1]));
                auto diffusion = sigma * rates[i - 1] * sqrt(d_t) * dist(gen);

                rates[i] = rates[i - 1] + drift * d_t + diffusion;
            }

            return rates;
        }
};

// pg 102
class CIRpp {
    private:
        auto x_0, kappa, theta, sigma;

        auto f_cir(const auto time, const auto kappa, const auto theta, const auto sigma) {
            auto h = sqrt(std::pow(kappa, 2) + 2 * std::pow(sigma, 2));

            auto first_term = 2 * kappa * theta * (exp(t * h) - 1) / (2 * h + (kappa + h) * (exp(t * h) - 1));
            auto second_term = x_0 * 4 * std::pow(h, 2) * exp(t * h) / std::pow(2 * h + (kappa + h) * (exp(t * h) - 1), 2);

            return first_term + second_term;
        }

        auto phi_cir(const auto time, const auto alpha) {
            return f_M(0, time) - f_cir(time, kappa, theta, sigma);
        }

    public:

};


class ExtendedExponentialVasicek {
    private:

    public:

};


int main() {
    return 0;
}