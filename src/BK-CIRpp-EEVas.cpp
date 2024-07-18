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

        std::vector<double> simulated_value(double num_time_steps, double T) {
            auto d_t = T / num_time_steps;

            std::vector<double> rates(num_time_steps, 0);
            rates[0] = r_0;

            for (int i = 1; i < num_time_steps; ++i) {
                auto drift = rates[i - 1] * (get_theta(i * d_t, r_0) + std::pow(sigma, 2) / 2 - kappa * std::log(rates[i - 1]));
                auto diffusion = sigma * rates[i - 1] * sqrt(d_t) * dist(gen);

                rates[i] = rates[i - 1] + drift * d_t + diffusion;
            }

            return rates;
        }
};

// I took a fair number of creative liberties on this model for simplicity purposes.
// If you desire a greater amount of realism, I encourage you to fork and add details.
class CIRpp {
    private:
        double x_0, kappa, theta, sigma;

        std::vector<double> maturities, yields;

        // Simple zero-coupon bond price. This is not the focus of the project, and rather on the methods.
        // Add realistic dynamics as you please.
        auto ZCBP(const double rate, const double expiration_time) {
            std::cout << "ZCB value: " << exp(-rate * expiration_time) << std::endl;
            return exp(-rate * expiration_time);
        }

        auto spot_rate(const double rate, const double expiration_time) {
            std::cout << "Spot rate: " << -log(ZCBP(rate, expiration_time)) / expiration_time << std::endl;
            return -log(ZCBP(rate, expiration_time)) / expiration_time;
        }

        // Simple implementation. Will make more realistic, time permitting.
        auto instant_forward_rate(const double rate, const double expiration_time, const double d_R_d_T) {
            std::cout << "instant forward rate: " << spot_rate(rate, expiration_time) + expiration_time * d_R_d_T << std::endl;
            return spot_rate(rate, expiration_time) + expiration_time * d_R_d_T;
        }

        auto f_cir(const double time, const double kappa, const double theta, const double sigma) {
            auto h = sqrt(std::pow(kappa, 2) + 2 * std::pow(sigma, 2));

            auto first_term = 2 * kappa * theta * (exp(time * h) - 1) / (2 * h + (kappa + h) * (exp(time * h) - 1));
            auto second_term = x_0 * 4 * std::pow(h, 2) * exp(time * h) / std::pow(2 * h + (kappa + h) * (exp(time * h) - 1), 2);

            return first_term + second_term;
        }

        auto phi_cir(const double time, const double last_rate) {
            return instant_forward_rate(last_rate, time, 0.001) - f_cir(time, kappa, theta, sigma); // 0.001 is placeholder for now
        }

    public:
        CIRpp(double x_0_con, 
              double theta_con, 
              double kappa_con, 
              double sigma_con, 
              std::vector<double> maturities_con, 
              std::vector<double> yields_con) {
            x_0 = x_0_con;
            kappa = kappa_con;
            sigma = sigma_con;
            theta = theta_con;
            maturities = maturities_con;
            yields = yields_con;
        }

        // "xates" represents the solution to the CIR, "rates" represents xates + (shift imposed by \phi)
        // See the documentation for more mathematical context
        std::vector<double> simulated_value(double num_time_steps, double T) {
            auto d_t = T / num_time_steps;

            std::vector<double> xates(num_time_steps, 0), rates(num_time_steps, 0);
            xates[0] = x_0; rates[0] = xates[0] + phi_cir(0, 0);

            for (int i = 1; i < num_time_steps; ++i) {
                xates[i] = std::max(0.0, rates[i - 1] + kappa * (theta - rates[i - 1]) * d_t + sigma * sqrt(std::max<double>(rates[i - 1], 0.0) * d_t) * dist(gen));
                rates[i] = xates[i - 1] + phi_cir(i * T, rates[i - 1]);
            }

            return rates;
        }

};


class ExtendedExponentialVasicek {
    private:

    public:

};


int main() {
    // BlackKarasinski bk_testing_class(0.05, 0.2, 0.02);
    // std::vector<double> bk_temp = bk_testing_class.simulated_value(10000, 1);
    // std::cout << "Black-Karasinski simulated rate: " << bk_temp[bk_temp.size() - 1] << std::endl;


// double x_0_con, 
//               double theta_con, 
//               double kappa_con, 
//               double sigma_con, 
//               std::vector<double> maturities_con, 
//               std::vector<double> yields_con

    CIRpp cirpp_testing_class(0.05, 0.1, 0.2, 0.02, {1}, {0.02});
    std::vector<double> cirpp_temp = cirpp_testing_class.simulated_value(10000, 1);
    std::cout << "CIR++ simulated rate: " << cirpp_temp[cirpp_temp.size() - 1] << std::endl;


    // std::vector<double> maturities = {1, 2, 3, 4, 5}; // Maturities in years
    // std::vector<double> yields = {0.02, 0.025, 0.03, 0.035, 0.04}; // Spot rates
    return 0;
}