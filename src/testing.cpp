// Importing header files
#include "Vasicek-CIR.hpp"
#include "General-CIR.hpp"

int main() {
    // Uncomment models for testing

    // Vasicek vas_testing_class(0.1, 0.02, 10, 0.06);
    // std::cout << "Vasicek exact rate: " << vas_testing_class.exact_value(0, 1) << std::endl;
    // std::cout << "Vasicek expected rate: " << vas_testing_class.expected_rate(0, 1) << std::endl;
    // std::cout << "Vasicek expected variance: " << vas_testing_class.expected_variance(0, 1) << std::endl;
    // std::vector<double> vas_temp = vas_testing_class.simulated_value(10000, 1);
    // std::cout << "Vasicek simulated rate: " << vas_temp[vas_temp.size() - 1] << std::endl;


    // CIR cir_testing_class(0.1, 0.02, 10, 0.06);
    // std::cout << "CIR exact rate: " << cir_testing_class.exact_value(0, 1) << std::endl;
    // std::cout << "CIR expected rate: " << cir_testing_class.expected_rate(0, 1) << std::endl;
    // std::cout << "CIR expected variance: " << cir_testing_class.expected_variance(0, 1) << std::endl;
    // std::vector<double> cir_temp = cir_testing_class.simulated_value(10000, 1);
    // std::cout << "CIR simulated rate: " << cir_temp[cir_temp.size() - 1] << std::endl;


    // StableCIR stable_testing_class(1, 0, 0.12, 0.02, 3, 0.4);
    // std::vector<double> stb_temp = stable_testing_class.simulated_value(1000000, 1);
    // std::cout << "Stable CIR simulated rate: " << stb_temp[stb_temp.size() - 1] << std::endl;
    

    // Defining alpha and variance values for alpha-CIR

    // std::vector<double> temp_alphas = {2, 1.5};
    // std::vector<double> temp_variances = {0.15, 0.3 * tgamma(0.5) / 0.75};

    // AlphaCIR alph_testing_class(0.12, 0.02, 500000, temp_alphas, temp_variances);
    // std::vector<double> alph_temp = alph_testing_class.simulated_value(1000000, 1);
    // std::cout << "Alpha-CIR simulated rate: " << alph_temp[alph_temp.size() - 1] << std::endl;
}