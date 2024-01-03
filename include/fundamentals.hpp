#ifndef FUNDAMENTALS_HPP
#define FUNDAMENTALS_HPP 1

#include <string>
#include <iostream>
#include <unordered_map>

namespace bioinfo {
    typedef std::unordered_map<std::string, char> AATable;

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
};

#endif // FUNDAMENTALS_HPP