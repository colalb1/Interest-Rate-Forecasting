#ifndef MY_HEADER_H
#define MY_HEADER_H

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

// Initializing standard normal distribution
std::random_device rd;
std::mt19937 gen(rd());
std::normal_distribution<double> dist(0, 1);
auto distribution = bind(dist, gen);

// Loop blocking size
// Hopefully this will speed up performance. This needs to be checked with Intel VTune once the program is more developed and diagnostics are more insightful.
const int BLOCK_SIZE = 64;

// Initializing global variables that the models will inherit
class GeneralModel {
    protected:
        // Model variables
        double r_0, theta, kappa, sigma, d_t;    // intial rate, mean level (the price the process reverts to), reversion rate, volatility, and time increment

        // Class constructor
        GeneralModel(double r_0_con, double theta_con, double kappa_con, double sigma_con, double d_t_con) {
            r_0 = r_0_con;
            theta = theta_con;
            kappa = kappa_con;
            sigma = sigma_con;
            d_t = d_t_con;
        }
};


// Interest rate models
class Vasicek: private GeneralModel {
    private: 
        double B(const double& t, const double& T) {
            return (1 - exp(-kappa * (T - t))) / kappa;
        }

        double A(const double& t, const double& T) {
            double temp_B = B(t, T);
            return (theta - pow(sigma / kappa, 2) / 2) * (temp_B - T + t) - pow(sigma * temp_B, 2) / (4 * kappa);
        }

    public:
        // Constructor, uses same construction as parent.
        Vasicek(double r_0_con, double theta_con, double kappa_con, double sigma_con, double d_t_con):GeneralModel(r_0_con, theta_con, kappa_con, sigma_con, d_t_con) {};

        double exact_value(const double& t, const double& T) {
            return exp(A(t, T) - r_0 * B(t, T));
        }

        double expected_rate(const double& t, const double& T) {
            return r_0 * exp(-kappa * (T - t)) + theta * (1 - exp(-kappa * (T - t)));
        }

        double expected_variance(const double& t, const double& T) {
            return pow(sigma, 2) * (1 - exp(-2 * kappa * (T - t))) / (kappa * 2);
        }

        std::vector<double> simulated_value(const int& num_sims, const int& num_time_steps, const double& T) {
            const double d_t = T / num_time_steps;
            std::vector<double> rates(num_time_steps, 0);
            rates[0] = r_0;

            for (int i = 1; i < num_time_steps + 1; i += BLOCK_SIZE) {
                for (int j = i; j < std::min(i + BLOCK_SIZE, num_time_steps + 1); ++j) {
                    rates[j] = rates[j - 1] + kappa * (theta - rates[j - 1]) * d_t + sigma * sqrt(d_t) * distribution(num_sims);
                }
            }

            return rates;
        }
};

class CIR: private GeneralModel {
    private:
        double B(const double& t, const double& T) {
            const double tau = T - t;
            const double h = sqrt(pow(kappa, 2) + 2 * pow(sigma, 2));

            return 2 * (exp(tau * h) - 1) / (2 * h + (kappa + h) * (exp(tau * h) - 1));
        }

        double A(const double& t, const double& T) {
            const double tau = T - t;
            const double h = sqrt(pow(kappa, 2)) + 2 * pow(sigma, 2);

            return pow(2 * h * (exp((kappa + h) * tau / 2)) / (2 * h + (kappa + h) * (exp(tau * h) - 1)), 2 * kappa * theta / pow(sigma, 2));
        }
    
    public:
        // Constructor
        CIR(double r_0_con, double theta_con, double kappa_con, double sigma_con, double d_t_con):GeneralModel(r_0_con, theta_con, kappa_con, sigma_con, d_t_con) {};

        double exact_value(const double& t, const double& T) {
            return A(t, T) * exp(-r_0 * B(t, T));
        }

        double expected_rate(const double& t, const double& T) {
            return r_0 * exp(-kappa * (T - t)) + theta * (1 - exp(-kappa * (T - t)));
        }

        double expected_variance(const double& t, const double& T) {
            return r_0 * pow(sigma, 2) / kappa * (exp(-kappa * (T - t)) - exp(-2 * kappa * (T - t))) + theta * pow(sigma, 2) * pow(1 - exp(-kappa * (T - t)), 2) / (2 * kappa);
        }

        std::vector<double> simulated_value(const int& num_sims, const int& num_time_steps, const double& T) {
            const double d_t = T / num_time_steps;
            std::vector<double> rates(num_time_steps, 0);
            rates[0] = r_0;

            for (int i = 1; i < num_time_steps + 1; i += BLOCK_SIZE) {
                for (int j = i; j < std::min(i + BLOCK_SIZE, num_time_steps + 1); ++j) {
                    rates[j] = rates[j - 1] + kappa * (theta - rates[j - 1]) * d_t + sigma * sqrt(std::max<double>(rates[j - 1], 0.0) * d_t) * distribution(num_sims);
                }
            }

            return rates;
        }
};

#endif