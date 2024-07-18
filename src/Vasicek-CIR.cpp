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
std::normal_distribution<double> dist(0.0, 1.0);

// Loop blocking size. Should be checked with Intel VTune for each unique machine.
const int BLOCK_SIZE = 64;

// Initializing global variables
class GeneralModel {
    protected:
        // Model variables
        double r_0, theta, kappa, sigma;    // initial rate, mean level (the price the process reverts to), reversion rate, volatility

        // Class constructor
        GeneralModel(double r_0_con, double theta_con, double kappa_con, double sigma_con) {
            r_0 = r_0_con;
            theta = theta_con;
            kappa = kappa_con;
            sigma = sigma_con;
        }
};


// Basic interest rate models
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
        // Constructor, uses parent contructor
        Vasicek(double r_0_con, 
                double theta_con, 
                double kappa_con, 
                double sigma_con): GeneralModel(r_0_con, 
                                                theta_con, 
                                                kappa_con, 
                                                sigma_con) {};

        double exact_value(const double& t, const double& T) {
            return exp(A(t, T) - r_0 * B(t, T));
        }

        double expected_rate(const double& t, const double& T) {
            return r_0 * exp(-kappa * (T - t)) + theta * (1 - exp(-kappa * (T - t)));
        }

        double expected_variance(const double& t, const double& T) {
            return pow(sigma, 2) * (1 - exp(-2 * kappa * (T - t))) / (kappa * 2);
        }

        std::vector<double> simulated_value(const int& num_time_steps, const double& T) {
            const double d_t = T / num_time_steps;
            std::vector<double> rates(num_time_steps, 0);
            rates[0] = r_0;

            for (int i = 1; i < num_time_steps; i += BLOCK_SIZE) {
                for (int j = i; j < std::min(i + BLOCK_SIZE, num_time_steps); ++j) {
                    rates[j] = rates[j - 1] + (theta - kappa * rates[j - 1]) * d_t + sigma * sqrt(d_t) * dist(gen);
                }
            }

            return rates;
        }
};


class ExponentialVasicek: private GeneralModel {
    public:
        // Constructor, uses parent contructor
        ExponentialVasicek(double r_0_con, 
                           double theta_con, 
                           double kappa_con, 
                           double sigma_con): GeneralModel(r_0_con, 
                                                           theta_con, 
                                                           kappa_con, 
                                                           sigma_con) {};


        // Expected rate and expected variance assume that T\to\infty
        // See page 71 of this book: Interest Rate Models - Theory and Practice: With Smile, Inflation and Credit by Brigo and Mercurio
        double expected_rate(const double& t, const double& T) {
            return exp(theta / kappa + std::pow(sigma, 2) / (4 * kappa));
        }

        double expected_variance(const double& t, const double& T) {
            return exp(2 * theta / kappa + std::pow(sigma, 2) / (2 * kappa)) * (exp(std::pow(sigma, 2) / (2 * kappa)) - 1);
        }

        std::vector<double> simulated_value(const int& num_time_steps, const double& T) {
            const auto d_t = T / num_time_steps;
            std::vector<double> rates(num_time_steps, 0);
            rates[0] = r_0;

            for (int i = 1; i < num_time_steps; i += BLOCK_SIZE) {
                for (int j = i; j < std::min(i + BLOCK_SIZE, num_time_steps); ++j) {
                    auto drift = rates[j - 1] * (theta + std::pow(sigma, 2) - kappa * log(rates[j - 1])) * d_t;
                    auto diffusion = sigma * rates[j - 1] * sqrt(d_t) * dist(gen);

                    rates[j] = rates[j - 1] + drift + diffusion;
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
        CIR(double r_0_con, 
            double theta_con, 
            double kappa_con, 
            double sigma_con): GeneralModel(r_0_con, 
                                            theta_con, 
                                            kappa_con, 
                                            sigma_con) {};

        double exact_value(const double& t, const double& T) {
            return A(t, T) * exp(-r_0 * B(t, T));
        }

        double expected_rate(const double& t, const double& T) {
            return r_0 * exp(-kappa * (T - t)) + theta * (1 - exp(-kappa * (T - t)));
        }

        double expected_variance(const double& t, const double& T) {
            return r_0 * pow(sigma, 2) / kappa * (exp(-kappa * (T - t)) - exp(-2 * kappa * (T - t))) + theta * pow(sigma, 2) * pow(1 - exp(-kappa * (T - t)), 2) / (2 * kappa);
        }

        std::vector<double> simulated_value(const int& num_time_steps, const double& T) {
            const double d_t = T / num_time_steps;
            std::vector<double> rates(num_time_steps, 0);
            rates[0] = r_0;

            for (int i = 1; i < num_time_steps; i += BLOCK_SIZE) {
                for (int j = i; j < std::min(i + BLOCK_SIZE, num_time_steps); ++j) {
                    rates[j] = std::max(0.0, rates[j - 1] + kappa * (theta - rates[j - 1]) * d_t + sigma * sqrt(std::max<double>(rates[j - 1], 0.0) * d_t) * dist(gen));
                }
            }

            return rates;
        }
};


int main() {
    // Uncomment models for testing

    Vasicek vas_testing_class(0.1, 0.02, 10, 0.06);
    std::cout << "Vasicek exact rate: " << vas_testing_class.exact_value(0, 1) << std::endl;
    std::cout << "Vasicek expected rate: " << vas_testing_class.expected_rate(0, 1) << std::endl;
    std::cout << "Vasicek expected variance: " << vas_testing_class.expected_variance(0, 1) << std::endl;
    std::vector<double> vas_temp = vas_testing_class.simulated_value(10000, 1);
    std::cout << "Vasicek simulated rate: " << vas_temp[vas_temp.size() - 1] << std::endl;


    // CIR cir_testing_class(0.1, 0.02, 10, 0.06);
    // std::cout << "CIR exact rate: " << cir_testing_class.exact_value(0, 1) << std::endl;
    // std::cout << "CIR expected rate: " << cir_testing_class.expected_rate(0, 1) << std::endl;
    // std::cout << "CIR expected variance: " << cir_testing_class.expected_variance(0, 1) << std::endl;
    // std::vector<double> cir_temp = cir_testing_class.simulated_value(10000, 1);
    // std::cout << "CIR simulated rate: " << cir_temp[cir_temp.size() - 1] << std::endl;


    ExponentialVasicek exp_vas_testing_class(0.1, 0.02, 10, 0.06);
    std::cout << "Exp Vas expected rate: " << exp_vas_testing_class.expected_rate(0, 1) << std::endl;
    std::cout << "Exp Vas expected variance: " << exp_vas_testing_class.expected_variance(0, 1) << std::endl;
    std::vector<double> exp_vas_temp = exp_vas_testing_class.simulated_value(10000, 1);
    std::cout << "Exp Vas simulated rate: " << exp_vas_temp[exp_vas_temp.size() - 1] << std::endl;

    return 0;
}