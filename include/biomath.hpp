#ifndef BIOMATH_HPP
#define BIOMATH_HPP 1

#include <vector>
#include <string>

namespace bioinfo {
    double factorial(int n);
    unsigned long int fibonacci(unsigned int n);
    double binomialDistribution(unsigned int n, unsigned int x, double p);

    class TotalPermutations {
        private:
            unsigned int n;
            std::vector<std::vector<unsigned int>> p;
        public:
            TotalPermutations(unsigned int n);
            std::string getPermutationSummary();
    };
}

#endif