#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP 1

#include "fundamentals.hpp"
#include <vector>
#include <string>

namespace bioinfo {
    typedef std::unordered_map<char, double> MassTable;

    namespace MassTables {
        const MassTable MONOISOTOPIC_MASS_TABLE = {
            { 'A', 71.03711 },  { 'C', 103.00919 }, { 'D', 115.02694 }, { 'E', 129.04259 },
            { 'F', 147.06841 }, { 'G', 57.02146 },  { 'H', 137.05891 }, { 'I', 113.08406 },
            { 'K', 128.09496 }, { 'L', 113.08406 }, { 'M', 131.04049 }, { 'N', 114.04293 },
            { 'P', 97.05276 },  { 'Q', 128.05858 }, { 'R', 156.10111 }, { 'S', 87.03203 },
            { 'T', 101.04768 }, { 'V', 99.06841 },  { 'W', 186.07931 }, { 'Y', 163.06333 }
        };
    }

    struct DirectedEdge {
        std::string tail;
        std::string head;
    };

    class AdjacencyList {
        private:
            std::vector<DirectedEdge> dsde;
            template <typename T> void internalConstructor(std::vector<T> &vec, unsigned int ok);
        public:
            AdjacencyList(std::vector<DNAString> &vec, unsigned int ok);
            std::string toString();

    };

    std::vector<unsigned int> exactDNAStringMotif(DNAString &ds, DNAString &motif, bool overlap);
    unsigned int hammingDistance(DNAString &s, DNAString &t);
    double proteinMass(AAString &as, const MassTable &mt);
    unsigned int inferredRNACount(AAString &as, const AATranscribableUnitTable &ut, unsigned int m);
    RNAString spliceRNA(RNAString &s, std::vector<RNAString> &introns);
}

#endif