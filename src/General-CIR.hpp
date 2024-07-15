#ifndef GENERAL_CIR_HEADER
#define GENERAL_CIR_HEADER

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

// Loop blocking size
const int BLOCK_SIZE = 64;

// Defining Expects
void Expects(bool condition, const char* message = "Precondition failed") {
    if (!condition) {
        throw std::logic_error(message);
    }
}

// Initializing global variables that the models will inherit
class GeneralModelAlpha {
    protected:
        // Model variables
        double r_0, theta, kappa, sigma;    // initial rate, mean level (the price the process reverts to), reversion rate, volatility

        // Constructor
        GeneralModelAlpha(double r_0_con, double theta_con, double kappa_con, double sigma_con) {
            r_0 = r_0_con;
            theta = theta_con;
            kappa = kappa_con;
            sigma = sigma_con;
        }

        std::vector<double> alphas, etas;

        // Overloaded constructor for AlphaCIR
        GeneralModelAlpha(double r_0_con, double theta_con, double kappa_con, std::vector<double> alphas_con, std::vector<double> etas_con) {
            r_0 = r_0_con;
            theta = theta_con;
            kappa = kappa_con;
            alphas = alphas_con;
            etas = etas_con;
        }

        // Helper function for alpha-stable distribution generator
        double kappa_sign(double& val) {
            if (val < 1) {
                return val;
            } else if (val > 1) {
                return val - 2;
            }
            return 0;
        }

        // Distribution for models using non-Brownian motion models
        // Using Chamber-Mallows-Stuck method: https://www.sciencedirect.com/science/article/pii/0167715295001131
        // alpha = shape (stability) parameter in (0, 2], beta = skewness param in [-1, 1], mu = location aka r_0, sigma = dispersion > 0
        double alpha_stable_pdf(double alpha, double beta, double mu, double sigma) {
            Expects(0 < alpha && alpha <= 2 && -1 < beta && beta < 1 && sigma > 0, "Alpha-stable pdf conditions violated");

            // Defining generator and distributions for path simulations
            std::default_random_engine generator;
            std::uniform_real_distribution<double> unif_distr(-M_PI_2, M_2_PI);
            std::exponential_distribution<double> exp_distr(1.0);
            const double gamma_0 = -M_PI_2 * beta * kappa_sign(alpha) / alpha;

            double rand_unif_value = unif_distr(generator);
            double rand_exp_value = exp_distr(generator);

            double solution = 0;

            if (alpha == 1) {
                const double first_term = (M_PI_2 + beta * rand_unif_value) * std::tan(rand_unif_value);
                const double second_term = beta * std::log(rand_exp_value * std::cos(rand_unif_value) / (M_PI_2 + beta * rand_unif_value));

                double X = M_2_PI * (first_term - second_term);

                solution = sigma * X + M_2_PI * beta * sigma * std::log(sigma) + mu;
            } else {
                const double S_term = std::pow(1 + std::pow(beta * std::tan(M_PI_2 * alpha), 2), 1 / (2 * alpha));
                const double B_term = std::atan(beta * std::tan(M_PI_2 * alpha)) / alpha;

                const double sin_cos_term = std::sin(alpha * (rand_unif_value + B_term)) / std::pow(std::cos(rand_unif_value), 1 / alpha);
                const double cos_pow_term = std::pow(std::cos(rand_unif_value - alpha * (rand_unif_value + B_term)) / rand_exp_value, (1 - alpha) / alpha);

                const double X = S_term * sin_cos_term * cos_pow_term;


                solution = sigma * X + mu;
            }

            return solution;
        }
};

// Model from this article: https://arxiv.org/abs/2402.07503
class StableCIR: private GeneralModelAlpha {
    private:
        // Shape and skew
        double alpha = 2, beta = 0;

    public:
        StableCIR(double alpha_con, 
                  double beta_con,
                  double r_0_con,
                  double theta_con,
                  double kappa_con, 
                  double sigma_con): alpha(alpha_con), 
                                     beta(beta_con), 
                                     GeneralModelAlpha(r_0_con, theta_con, kappa_con, sigma_con) {}
        
        // Create the alpha-stable path
        std::vector<double> simulated_value(long int num_time_steps, double T) {
            double d_t = T / num_time_steps;

            std::vector<double> rates(num_time_steps, 0);
            rates[0] = r_0;

            for (int i = 1; i < num_time_steps; i += BLOCK_SIZE) {
                for (int j = i; j < std::min(static_cast<long int> (i + BLOCK_SIZE), num_time_steps); ++j) {
                    double dZ_alpha = alpha_stable_pdf(alpha, beta, r_0, sigma) * std::pow(d_t, 1 / alpha);

                    // diffeq from 2.23 of "Affine term structure models driven by independent Levy process"
                    rates[j] = rates[j - 1] + kappa * (theta - rates[j - 1]) * d_t + sigma * std::pow(rates[j - 1], 1 / alpha) * dZ_alpha;
                }
            }
            
            return rates;
        }
};

// The alpha-CIR model is better than CIR for real-world modeling since it allows large fluctuations since it has a tail-fatness parameter
// and reduces overestimation when interests rates are low, a common issue with the CIR model.
class AlphaCIR: private GeneralModelAlpha {
    private:
        // first is an indicator for first iteration
        long double d_calc(const double& alpha, const double& eta, bool first = false) {
            if (first && alpha == 2) {
                return 2 * eta;
            }

            // Checks \alpha\in (1, 2)
            Expects(1 < alpha && alpha < 2, "Variance calculation violated");
            return eta * alpha * (alpha - 1) / tgamma(2 - alpha);
        }

    public:
        AlphaCIR(double r_0_con,
                 double theta_con,
                 double kappa_con,
                 std::vector<double> alphas_con, // Input should be reverse sorted; READ THE DOCUMENTATION
                 std::vector<double> etas_con): GeneralModelAlpha(r_0_con, 
                                                                  theta_con, 
                                                                  kappa_con, 
                                                                  alphas_con, 
                                                                  etas_con) {}
        
        // Create the alpha-stable path
        std::vector<double> simulated_value(long int num_time_steps, double T) {
            double d_t = T / num_time_steps;

            std::vector<double> d_variance(etas.size(), 0);
            d_variance[0] = d_calc(alphas[0], etas[0], true);

            for (int i = 1; i < d_variance.size(); ++i) {
                d_variance[i] = d_calc(alphas[i], etas[i]);
            }

            std::vector<double> rates(num_time_steps, 0);
            rates[0] = r_0;

            for (int i = 1; i < num_time_steps; i += BLOCK_SIZE) {
                for (int j = i; j < std::min(static_cast<long int> (i + BLOCK_SIZE), num_time_steps); ++j) {
                    double dZ_alpha_sum = 0;

                    for (int k = 0; k < d_variance.size(); ++k) {
                        dZ_alpha_sum += alpha_stable_pdf(alphas[k], etas[k], r_0, d_variance[k]) * std::pow(d_variance[k] * d_t * rates[j - 1], 1 / alphas[k]);
                    }

                    rates[j] = rates[j - 1] + kappa * (theta - rates[j - 1]) * d_t + dZ_alpha_sum;
                }
            }
            
            return rates;
        }
};


# endif