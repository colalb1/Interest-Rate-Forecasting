// Importing header files
#include "Vasicek-CIR.hpp"
// #include "General-CIR.hpp"

int main() {
    Vasicek testing_class(0.05, 0.02, 10, 0.06);
    std::cout << "Vasicek exact rate: " << testing_class.exact_value(0, 1) << std::endl;
    std::cout << "Vasicek expected rate: " << testing_class.expected_rate(0, 1) << std::endl;
    std::cout << "Vasicek expected variance: " << testing_class.expected_variance(0, 1) << std::endl;

    std::vector<double> temp = testing_class.simulated_value(1000000, 1);
    std::cout << "Vasicek simulated rate: " << temp[temp.size() - 1] << std::endl; 

    CIR testing_class1(0.05, 0.02, 10, 0.06);
    std::cout << "CIR exact rate: " << testing_class1.exact_value(0, 1) << std::endl;
    std::cout << "CIR expected rate: " << testing_class1.expected_rate(0, 1) << std::endl;
    std::cout << "CIR expected variance: " << testing_class1.expected_variance(0, 1) << std::endl;

    std::vector<double> temp1 = testing_class1.simulated_value(1000000, 1);
    std::cout << "CIR simulated rate: " << temp1[temp1.size() - 1] << std::endl; 
}