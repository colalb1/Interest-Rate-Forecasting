// Importing header files
// #include "Vasicek-CIR.hpp"
#include "General-CIR.hpp"

int main() {
    // Vasicek testing_class(0.05, 0.02, 10, 0.06);
    // std::cout << "Vasicek exact rate: " << testing_class.exact_value(0, 1) << std::endl;
    // std::cout << "Vasicek expected rate: " << testing_class.expected_rate(0, 1) << std::endl;
    // std::cout << "Vasicek expected variance: " << testing_class.expected_variance(0, 1) << std::endl;

    StableCIR testing_class(1, 0, 0.12, 0.02, 3, 0.4);
    std::vector<double> temp = testing_class.simulated_value(1000, 1);
    std::cout << "Stable-CIR simulated rate: " << temp[temp.size() - 1] << std::endl;
}