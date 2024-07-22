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

// This is technically the Extended Hull-White model since kappa is time dependent as well
class HullWhite {
    private:
        double r_0;
        std::vector<double> theta, kappa, sigma; // One could also write function of time that represent each of these; having the vector as an input is the choice I made for simplicity.
                                                 // Aka, I'm letting the user bear the cross of calculating the function. It is more modular this way. If you disagree, send me an email: dormantemail22@gmail.com
    public:
        HullWhite(double r_0_con, std::vector<double> theta_con, std::vector<double> kappa_con, std::vector<double> sigma_con) {
            r_0 = r_0_con;
            theta = theta_con;
            kappa = kappa_con;
            sigma = sigma_con;
        }

        std::vector<double> simulated_value(const int& num_time_steps, const double& T) {
            auto d_t = T / num_time_steps;

            std::vector<double> rates(num_time_steps, 0);
            rates[0] = r_0;

            for (int i = 1; i < num_time_steps; ++i) {
                rates[i] = rates[i - 1] + (theta[i - 1] - kappa[i - 1] * rates[i - 1]) * d_t + sigma[i - 1] * sqrt(d_t) * dist(gen);
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