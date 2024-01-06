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

    bioinfo::RNAString a = bioinfo::RNAString(
        std::string("Rosalind_10"),
        std::string("ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG")
    );

    bioinfo::RNAString b = bioinfo::RNAString(
        std::string("Rosalind_12"),
        std::string("ATCGGTCGAA")
    );

    bioinfo::RNAString c = bioinfo::RNAString(
        std::string("Rosalind_15"),
        std::string("ATCGGTCGAGCGTGT")
    );

    std::vector<bioinfo::RNAString> vec;

    vec.push_back(b);
    vec.push_back(c);

    bioinfo::RNAString mrna = bioinfo::spliceRNA(a, vec);
    bioinfo::AAString protein = bioinfo::AAString(mrna, bioinfo::GeneticCode::STANDARD_GENETIC_CODE);

    std::cout << "pre-mRNA: " << a.getSequence() << std::endl;
    std::cout << "mRNA:     " << mrna.getSequence() << std::endl;
    std::cout << "Protein:  " << protein.getSequence() << std::endl;

    return 0;
}