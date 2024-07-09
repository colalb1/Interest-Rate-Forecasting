// Importing header files
// #include "Vasicek-CIR.hpp"
#include "General-CIR.hpp"

int main() {
    // CIR testing_class(0.1, 0.02, 10, 0.06);
    // std::cout << "CIR exact rate: " << testing_class.exact_value(0, 1) << std::endl;
    // std::cout << "CIR expected rate: " << testing_class.expected_rate(0, 1) << std::endl;
    // std::cout << "CIR expected variance: " << testing_class.expected_variance(0, 1) << std::endl;
    // std::vector<double> temp = testing_class.simulated_value(10000, 1);
    // std::cout << "CIR simulated rate: " << temp[temp.size() - 1] << std::endl;

    // StableCIR testing_class(1, 0, 0.12, 0.02, 3, 0.4);
    

    std::vector<double> temp_alphas = {2, 1.5};
    std::vector<double> temp_variances = {0.15, 0.3 * tgamma(0.5) / 0.75};


    AlphaCIR testing_class(0.12, 0.02, 500000, temp_alphas, temp_variances);
    std::vector<double> temp = testing_class.simulated_value(1000000, 1);
    std::cout << "Alpha-CIR simulated rate: " << temp[temp.size() - 1] << std::endl;
}