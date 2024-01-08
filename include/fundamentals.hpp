#ifndef FUNDAMENTALS_HPP
#define FUNDAMENTALS_HPP 1

#include <string>
#include <iostream>
#include <unordered_map>
#include <vector>

namespace bioinfo {
    typedef std::unordered_map<std::string, char> AATable;
    typedef std::unordered_map<char, unsigned char> AATranscribableUnitTable;

    namespace GeneticCode {
        // Hashtable for a RNA codon to a single letter amino acid character 
        const AATable STANDARD_GENETIC_CODE = {
            { "UUU", 'F' }, { "UUC", 'F' }, { "UUA", 'L' }, { "UUG", 'L' },
            { "UCU", 'S' }, { "UCC", 'S' }, { "UCA", 'S' }, { "UCG", 'S' },
            { "UAU", 'Y' }, { "UAC", 'Y' }, { "UAA", '*' }, { "UAG", '*' },
            { "UGU", 'C' }, { "UGC", 'C' }, { "UGA", '*' }, { "UGG", 'W' },

            { "CUU", 'L' }, { "CUC", 'L' }, { "CUA", 'L' }, { "CUG", 'L' },
            { "CCU", 'P' }, { "CCC", 'P' }, { "CCA", 'P' }, { "CCG", 'P' },
            { "CAU", 'H' }, { "CAC", 'H' }, { "CAA", 'Q' }, { "CAG", 'Q' },
            { "CGU", 'R' }, { "CGC", 'R' }, { "CGA", 'R' }, { "CGG", 'R' },

            { "AUU", 'I' }, { "AUC", 'I' }, { "AUA", 'I' }, { "AUG", 'M' },
            { "ACU", 'T' }, { "ACC", 'T' }, { "ACA", 'T' }, { "ACG", 'T' },
            { "AAU", 'N' }, { "AAC", 'N' }, { "AAA", 'K' }, { "AAG", 'K' },
            { "AGU", 'S' }, { "AGC", 'S' }, { "AGA", 'R' }, { "AGG", 'R' },

            { "GUU", 'V' }, { "GUC", 'V' }, { "GUA", 'V' }, { "GUG", 'V' },
            { "GCU", 'A' }, { "GCC", 'A' }, { "GCA", 'A' }, { "GCG", 'A' },
            { "GAU", 'D' }, { "GAC", 'D' }, { "GAA", 'E' }, { "GAG", 'E' },
            { "GGU", 'G' }, { "GGC", 'G' }, { "GGA", 'G' }, { "GGG", 'G' }
        };
    };

    namespace GeneticCodeTranscribableUnits {
        const AATranscribableUnitTable STANDARD_GENETIC_CODE = {
            {'A', 4}, {'R', 6}, {'N', 2}, {'D', 2}, {'C', 2}, 
            {'Q', 2}, {'E', 2}, {'G', 4}, {'H', 2}, {'I', 3}, 
            {'L', 6}, {'K', 2}, {'M', 1}, {'F', 2}, {'P', 4}, 
            {'S', 6}, {'T', 4}, {'W', 1}, {'Y', 2}, {'V', 4}, 
            {'*', 3}
        };
    };

    class DNAString {
        private:
            std::string header;
            std::string sequence;
            unsigned int sequenceLength;
        public:
            DNAString(std::string h, std::string s);
            std::string getHeader();
            std::string getSequence();            
            void setHeader(std::string h);
            void setSequence(std::string s);

            unsigned int getSequenceLength();
    };

    class RNAString {
        private:
            std::string header;
            std::string sequence;
            unsigned int sequenceLength;
        public:
            RNAString();
            RNAString(std::string h, std::string s);
            RNAString(DNAString &ds);
            std::string getHeader();
            std::string getSequence();            
            void setHeader(std::string h);
            void setSequence(std::string s);

            unsigned int getSequenceLength();
    };

    class AAString {
        private:
            std::string header;
            std::string sequence;
            unsigned int sequenceLength;
        public:
            AAString(std::string h, std::string s, const AATable &code);
            AAString(DNAString &ds, const AATable &code);
            AAString(RNAString &rs, const AATable &code);

            std::string getHeader();
            std::string getSequence();            
            void setHeader(std::string h);
            void setSequence(std::string s, const AATable &code);

            unsigned int getSequenceLength();
    };

    std::string toUpper(std::string s);
    std::string transcribe(std::string s);
    std::string translate(std::string s, const AATable &code);
    DNAString reverseComplement(DNAString &ds);

    std::vector<DNAString> readDNAStringFile(std::string &fn);
    std::vector<DNAString> readDNAStringFile(const char *fnp);
};

#endif // FUNDAMENTALS_HPP