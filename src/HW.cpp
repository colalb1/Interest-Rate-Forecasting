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

class General {
    protected:
        void check_same_length(const std::vector<std::vector<double>>& vecs) {
            if (vecs.empty()) {
                throw std::invalid_argument("The list of vectors is empty.");
            }

            auto first_length = vecs[0].size();
            
            for (const auto& vec: vecs) {
                if (vec.size() != first_length) {
                    throw std::invalid_argument("(At least) two input vectors do not have the same length.");
                }
            }
        }

        int get_parameter_index(int step_num, int num_steps, int num_intervals) {
            int interval_size = num_steps / num_intervals;
            return step_num / interval_size;
        }
};

class HullWhite: private General {
    private:
        double r_0, kappa;
        std::vector<double> theta, sigma; // One could also write function of time that represent each of these; having the vector as an input is the choice I made for simplicity.
                                          // Aka, I'm letting the user bear the cross of calculating the function. It is more modular this way. If you disagree, send me an email: dormantemail22@gmail.com

    public:
        HullWhite(double r_0_con, 
                  std::vector<double> theta_con, 
                  double kappa_con, 
                  std::vector<double> sigma_con) {
            r_0 = r_0_con;
            theta = theta_con;
            kappa = kappa_con;
            sigma = sigma_con;
        }

        std::vector<double> simulated_value(const int& num_time_steps, const double& T) {
            check_same_length({theta, sigma});

            double d_t = T / num_time_steps;

            std::vector<double> rates(num_time_steps, 0);
            rates[0] = r_0;

            for (int i = 1; i < num_time_steps; i += BLOCK_SIZE) {
                for (int j = i; j < std::min<int>(i + BLOCK_SIZE, num_time_steps); ++j) {
                    int parameter_index = get_parameter_index(j, num_time_steps, theta.size());
                    double theta_curr = theta[parameter_index], sigma_curr = sigma[parameter_index];

                    rates[j] = rates[j - 1] + (theta_curr - kappa * rates[j - 1]) * d_t + sigma_curr * sqrt(d_t) * dist(gen);
                }
            }

            return rates;
        }
};

// This implementation is technically G2++.
// The Hull-White Two-Factor model is mathematically equivalent to G2++: 
// https://www.math.kth.se/matstat/seminarier/reports/M-exjobb12/120220b.pdf proof starts at botton of page 19
class HullWhiteTwoFactor: private General {
    private:
        double r_0, kappa_1, kappa_2, rho; // rho = instantaneous correlation of two-dimensional Brownian motion
        std::vector<double> theta, sigma_1, sigma_2;

    public:
        HullWhiteTwoFactor(double r_0_con, 
                           std::vector<double> theta_con, 
                           std::vector<double> kappas_con, 
                           std::vector<std::vector<double>> sigmas_con,
                           double rho_con) {
            r_0 = r_0_con;
            theta = theta_con;

            kappa_1 = kappas_con[0];
            kappa_2 = kappas_con[1];

            sigma_1 = sigmas_con[0];
            sigma_2 = sigmas_con[1];

            rho = rho_con;
        }


        std::vector<double> simulated_value(const int& num_time_steps, const double& T) {
            check_same_length({theta, sigma_1, sigma_2});
            assert(-1 <= rho && rho <= 1);

            auto d_t = T / num_time_steps;

            std::vector<double> rates(num_time_steps, 0);
            double x_rates = 0, y_rates = 0;
            rates[0] = r_0;

            for (int i = 0; i < num_time_steps; i += BLOCK_SIZE) {
                for (int j = i; j < std::min<int>(i + BLOCK_SIZE, num_time_steps); ++j) {
                    int parameter_index = get_parameter_index(j, num_time_steps, theta.size());
                    auto W_1 = sqrt(d_t) * dist(gen), W_2 = rho * W_1 + sqrt(1 - pow(rho, 2)) * sqrt(d_t) * dist(gen);

                    auto theta_curr = theta[parameter_index], sigma_1_curr = sigma_1[parameter_index], sigma_2_curr = sigma_2[parameter_index];

                    x_rates += -kappa_1 * x_rates * d_t + sigma_1_curr * sqrt(d_t) * W_1;
                    y_rates += -kappa_2 * y_rates * d_t + sigma_2_curr * sqrt(d_t) * W_2;

                    rates[j] = std::max(r_0 + x_rates + y_rates, 0.0);
                }
            }

            return rates;
        }
};


int main() {
    // HullWhite hw_testing_class(0.05, {0.03, 0.02, 0.025}, 0.2, {0.01, 0.015, 0.02});
    // std::vector<double> hw_temp = hw_testing_class.simulated_value(10000, 1);
    // std::cout << "Hull-White simulated rate: " << hw_temp[hw_temp.size() - 1] << std::endl;


    // HullWhiteTwoFactor hw2_testing_class(0.08, {0.03, 0.02, 0.025}, {0.2, 0.99}, {{0.65, 0.6, 0.7}, {0.01, 0.015, 0.02}}, 0.5);
    // std::vector<double> hw2_temp = hw2_testing_class.simulated_value(10000, 1);
    // std::cout << "Two-Factor Hull-White simulated rate: " << hw2_temp[hw2_temp.size() - 1] << std::endl;

    return 0;
}