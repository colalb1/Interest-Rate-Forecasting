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

class HullWhite {
    private:
        double r_0, kappa;
        std::vector<double> theta, sigma; // One could also write function of time that represent each of these; having the vector as an input is the choice I made for simplicity.
                                          // Aka, I'm letting the user bear the cross of calculating the function. It is more modular this way. If you disagree, send me an email: dormantemail22@gmail.com

        void check_same_length(const std::vector<std::vector<double>>& vecs) {
            if (vecs.empty()) {
                throw std::invalid_argument("The list of vectors is empty.");
            }

            auto first_length = vecs[0].size();
            
            for (const auto& vec: vecs) {
                if (vec.size() != first_length) {
                    throw std::invalid_argument("(At least) two vectors do not have the same length");
                }
            }
        }

        int get_parameter_index(int step_num, int num_steps, int num_intervals) {
            int interval_size = num_steps / num_intervals;
            return step_num / interval_size;
        }

    public:
        HullWhite(double r_0_con, std::vector<double> theta_con, double kappa_con, std::vector<double> sigma_con) {
            r_0 = r_0_con;
            theta = theta_con;
            kappa = kappa_con;
            sigma = sigma_con;
        }

        std::vector<double> simulated_value(const int& num_time_steps, const double& T) {
            check_same_length({theta, sigma});

            auto d_t = T / num_time_steps;

            std::vector<double> rates(num_time_steps, 0);
            rates[0] = r_0;

            for (int i = 1; i < num_time_steps; ++i) {
                int parameter_index = get_parameter_index(i, num_time_steps, theta.size());
                auto theta_curr = theta[parameter_index], sigma_curr = sigma[parameter_index];

                rates[i] = rates[i - 1] + (theta_curr - kappa * rates[i - 1]) * d_t + sigma_curr * sqrt(d_t) * dist(gen);
            }

            return rates;
        }
};


int main() {
    // HullWhite hw_testing_class(0.05, {0.03, 0.02, 0.025}, 0.2, {0.01, 0.015, 0.02});
    // std::vector<double> hw_temp = hw_testing_class.simulated_value(10000, 1);
    // std::cout << "Hull-White simulated rate: " << hw_temp[hw_temp.size() - 1] << std::endl;

    return 0;
}