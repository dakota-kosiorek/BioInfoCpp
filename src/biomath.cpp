#include <biomath.hpp>
#include <stdexcept>
#include <limits.h>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <algorithm>

namespace bioinfo {
    // Calculate the factorial of a number `n`.
    double factorial(int n) {
        double fn = 1;
        unsigned int i;

        if (n < 0) {
            throw std::invalid_argument("ERROR: Factorial input cannot be negative!");
        } else {
            for (i = 1; i <= n ; ++i) {
                /*if (i > 0 && fn > ULONG_MAX / i) {
                    throw std::overflow_error("ERROR: Factorial result overflowed!");
                }*/

                fn *= i;
            }
        }

        return fn;
    }

    // Calculate the fibonacci of a numer `n`.
    unsigned long int fibonacci(unsigned int n) {
        unsigned long int fn = 0;
        unsigned int i;
        unsigned long int currNum;

        std::vector<unsigned long int> nums = {1, 1};

        for(i = 2; i < n; ++i) {
            if (nums.at(i - 2) > 0 && nums.at(i - 1) > ULONG_MAX - nums.at(i - 2)) {
                throw std::overflow_error("ERROR: Fibonacci result overflowed!");
            }

            currNum = nums.at(i - 1) + nums.at(i - 2);

            nums.push_back(currNum);
        }

        fn = nums.back();
        return fn;
    }

    // Calculate the binomial distribution with `n` number of trials, `x` number of times for a specific 
    // outcome within `n` trials, and `p` probability of success on a single trial.
    double binomialDistribution(unsigned int n, unsigned int x, double p) {
        double prob = 0.0;
        double q = 1 - p;

        prob = ( factorial(n) / (factorial(n-x) * factorial(x)) ) * pow(p, x) * pow(q, n - x);

        return prob;
    }

    // Create a new `TotalPermutations` object that has all possible permuations of the positive integers with a 
    // permutations length of `n`.
    TotalPermutations::TotalPermutations(unsigned int n) {
        (*this).n = n;
        std::vector<std::vector<unsigned int>> p;
        std::vector<unsigned int> v;

        int i;

        if (n > 0) {
            for (i = 1; i < n + 1; ++i) {
                v.push_back(i);
            }

            do {
                p.push_back(v);
            } while (std::next_permutation(v.begin(), v.end()));
        }

        (*this).p = p;
    }

    // Returns a string summary with the total number of permutaions on the first line with a new permutations on each of the
    // following lines.
    std::string TotalPermutations::getPermutationSummary() {
        std::stringstream ss;
        std::vector<std::vector<unsigned int>>::iterator row;
        std::vector<unsigned int>::iterator col;

        ss << factorial(n) << "\n";

        for (row = (*this).p.begin(); row != (*this).p.end(); row++) {
            for (col = row->begin(); col != row->end(); col++) {
                ss << *col << " ";
            }

            if (row != (*this).p.end() - 1) {
                ss << "\n";
            }
        }

        return ss.str();
    }
}