#ifndef BK_CIRPP_EEVAS_HEADER
#define BK_CIRPP_EEVAS_HEADER

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
        double r_0, kappa, sigma;
        std::vector<double> theta;

        // Uncomment based on desired theta model
        double get_theta(double time, double r_init) {
            return r_init; // Constant
            // return r_init + 0.001 * time; // Linear with time
            // return r_init * exp(0.001 * time); // Exponential
        }

    public:
        BlackKarasinski(double r_0_con, double kappa_con, double sigma_con, std::vector<double> theta_con) {
            r_0 = r_0_con;
            kappa = kappa_con;
            sigma = sigma_con;
            theta = theta_con;
        }

        std::vector<double> simulated_value(long int num_time_steps, double T) {
            double d_t = T / num_time_steps;

            std::vector<double> rates(num_time_steps, 0);
            rates[0] = r_0;

            for (int i = 1; i < num_time_steps; ++i) {
                double drift = rates[i - 1] * (get_theta(i * d_t, r_0) + std::pow(sigma, 2) / 2 - kappa * std::log(rates[i - 1]));
                double diffusion = sigma * rates[i - 1] * sqrt(d_t) * dist(gen);

                rates[i] = rates[i - 1] + drift * d_t + diffusion;
            }

            return rates;
        }
};

// pg 102
class CIRpp {
    private:
        double x_0, kappa, theta, sigma;

    public:

};


class ExtendedExponentialVasicek {
    private:

    public:

};


# endif