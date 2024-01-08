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
    bioinfo::AdjacencyList dsal = bioinfo::AdjacencyList(vec, 3);

    std::cout << dsal.toString() << std::endl;

    return 0;
}