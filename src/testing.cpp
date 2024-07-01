// Importing header files
#include "Vasicek-CIR.hpp"

int main() {
    Vasicek testing_class(0.05, 0.02, 3, 0.15, 0.001);
    std::cout << "Vasicek exact rate: " << testing_class.exact_value(0, 1) << std::endl;
    std::cout << "Vasicek expected rate: " << testing_class.expected_rate(0, 1) << std::endl;
    std::cout << "Vasicek expected variance: " << testing_class.expected_variance(0, 1) << std::endl;

    std::vector<double> temp = testing_class.simulated_value(100000, 1000000, 1);
    std::cout << "Vasicek simulated rate: " << temp[0] << std::endl; 
}