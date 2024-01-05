// Standard libs
#include <iostream>
#include <iomanip>
#include <vector>
// Bioinformatics libs
#include <fundamentals.hpp>
#include <genetics.hpp>
#include <analysis.hpp>
#include <biomath.hpp>

int main() {
    std::cout << std::fixed << std::setprecision(3);

    bioinfo::TotalPermutations a = bioinfo::TotalPermutations(5);

    std::cout << a.getPermutationSummary() << std::endl;

    return 0;
}