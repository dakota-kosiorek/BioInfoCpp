#ifndef GENETICS_HPP
#define GENETICS_HPP 1

#include <unordered_map>
#include <limits.h>
#include <stdexcept>

namespace bioinfo {
    struct MendelianInheritanceStatistics {
        double recessivePhenotypeFrequency = 0.0;
        double dominantPhenotypeFrequency = 0.0;
    } typedef MendelianInheritanceStatistics;

    struct GenotypesCount {
        unsigned int dd = 0; // AA-AA mating
        unsigned int dh = 0; // AA-Aa mating
        unsigned int dr = 0; // AA-aa mating
        unsigned int hh = 0; // Aa-Aa mating
        unsigned int hr = 0; // Aa-aa mating
        unsigned int rr = 0; // aa-aa mating
    } typedef GenotypesCount;
    
    MendelianInheritanceStatistics mendelianInheritance(unsigned int k, unsigned int m, unsigned int n);
    double calculateExpectedOffspring(GenotypesCount mg, unsigned int n);
    double tomsIndependentAlleles(unsigned int k, unsigned long int n);

    class WascallyWabbits {
        private:
            unsigned int rabbitPairs = 0;
            unsigned int lifespan = 0;
            std::unordered_map<unsigned int, unsigned long int> timeline;
            unsigned long int fibonacci(unsigned int n);
            unsigned long int MortalFibonacci(unsigned int n, unsigned int m);
        public:
            WascallyWabbits(unsigned int k);
            unsigned long int simulate(unsigned int n);
            unsigned long int simulateMortal(unsigned int n, unsigned int m);
    };
}

#endif