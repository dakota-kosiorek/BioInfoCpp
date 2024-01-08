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

    std::vector<bioinfo::DNAString> vec = bioinfo::readDNAStringFile("rosalind_grph.txt");

    std::vector<bioinfo::DNAString>::iterator it;

    for (it = vec.begin(); it != vec.end(); it++) {
        std::cout << it->getHeader() << std::endl;
        std::cout << it->getSequence() << std::endl;
    }

    return 0;
}