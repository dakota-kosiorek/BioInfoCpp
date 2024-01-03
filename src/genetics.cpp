#include <genetics.hpp>
#include <biomath.hpp>
#include <unordered_map>
#include <limits.h>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <iostream>

namespace bioinfo {
    /*  
        Predict the frequency of individuals in a population possesing dominant and recessive phenotypes for the next generation 
        of a population of organisms given the current number of individuals of the population who are homozygous dominant (`k`), 
        heterozygous (`m`), and homozygous recessive (`n`) following a probability based mendelian inheritance approach.
    */
    MendelianInheritanceStatistics mendelianInheritance(unsigned int k, unsigned int m, unsigned int n) {
        unsigned int totalPop = k + m + n;
        unsigned int totalDomAlleles = (k * 2) + m;
        unsigned int totalRecAlleles = (n * 2) + m;
        unsigned int totalAlleles = totalPop * 2;

        // Two recessive organisms mating
        double rr = ((double) n / totalPop) * ((double) (n - 1) / (totalPop - 1));
        // Two heterozygous organisms mating
        double hh = ((double) m / totalPop) * ((double) (m - 1) / (totalPop - 1));
        // Heterozygous and recessive organisms mating
        double hr = ((double) m / totalPop) * ((double) n / (totalPop - 1)) + ((double) n / totalPop) * ((double) m / (totalPop - 1));

        double recessivePhenotypeFrequency = rr + hh * 0.25 + hr * 0.5;
        double dominantPhenotypeFrequency = 1 - recessivePhenotypeFrequency;

        MendelianInheritanceStatistics stats;
        
        stats.recessivePhenotypeFrequency = recessivePhenotypeFrequency;
        stats.dominantPhenotypeFrequency = dominantPhenotypeFrequency;

        return stats;
    }

    // Create a new WascallyWabbits object, with `k` numbers of rabbit pairs reproduced every generations.
    WascallyWabbits::WascallyWabbits(unsigned int k) {
        (*this).rabbitPairs = k;
    }

    // Calculate how many rabbits exists after `n` generations.
    unsigned long int WascallyWabbits::simulate(unsigned int n) {
        return fibonacci(n);
    }

    // Calculate how many rabbits exists after `n` generations with a rabbit lifespan of `m`.
    unsigned long int WascallyWabbits::simulateMortal(unsigned int n, unsigned int m) {
        return MortalFibonacci(n, m);
    }

    unsigned long int WascallyWabbits::fibonacci(unsigned int n) {
        unsigned long int fn = 0;

        if (n == 0) {
            fn = 0;
        } else if (n == 1) {
            fn = 1;
        }  else {
            unsigned long int fnMinusOne;
            unsigned long int fnMinusTwo;

            // See if Fn-1 has not yet been added to the timeline hashtable
            // Calculate and add it if it hasn't
            if ((*this).timeline.find(n - 1) == (*this).timeline.end()) {
                fnMinusOne = fibonacci(n - 1);
                (*this).timeline[n - 1] = fnMinusOne;
            } else {
                fnMinusOne = (*this).timeline[n - 1];
            }

            fnMinusTwo = (*this).timeline[n - 2];

            // See if there will be an overflow when adding together Fn-1 and Fn-2
            if (fnMinusTwo > 0 && fnMinusOne > ULONG_MAX - fnMinusTwo) {
                throw std::overflow_error("ERROR: WascallyWabbits calculated too many wabbits!");
            }

            fn = fnMinusOne + fnMinusTwo * (*this).rabbitPairs;
        }

        return fn;
    }

    unsigned long int WascallyWabbits::MortalFibonacci(unsigned int n, unsigned int m) {
        unsigned long int fn = 0;
        unsigned int i;
        unsigned long int currPop;

        std::vector<unsigned long int> populations = {1, 1};

        for(i = 2; i < n; ++i) {
            if (populations.at(i - 2) > 0 && populations.at(i - 1) > ULONG_MAX - populations.at(i - 2)) {
                throw std::overflow_error("ERROR: WascallyWabbits calculated too many wabbits!");
            }

            currPop = populations.at(i - 1) + populations.at(i - 2) * (*this).rabbitPairs;

            if (i == m) {
                currPop = currPop - 1;
            } else if (i > m) {
                currPop = currPop - populations.at(i - m - 1);
            }

            populations.push_back(currPop);
        }

        fn = populations.back();
        return fn;
    }

    // Calculate the expected number of offspring displaying the dominant phenotype in the next generation with `n` number of
    // offspring per mating gruop
    double calculateExpectedOffspring(GenotypesCount mg, unsigned int n) {
        double expected = 0.0;
        unsigned int i = 0;

        unsigned int genotypes[6] = {mg.dd, mg.dh, mg.dr, mg.hh, mg.hr, mg.rr};
        double dominantPhenotypeCoefficients[6] = {1.0, 1.0, 1.0, 0.75, 0.5, 0.0};

        for (i = 0; i <= 6; ++i) {
            expected += genotypes[i] * dominantPhenotypeCoefficients[i] * n;
        }


        return expected;
    }

    // Calculate the probability that at least `n` Aa Bb organisms will belong to the `k`-th generation of "Tom's" family tree and 
    // that Mendel's Second law  holds for the factors.
    // - Tom is the 0th geneation organism and has the Aa Bb genotype
    // - Toms has two children, who each have two children, and so on
    // - Each organism always mates with another Aa Bb organism
    double tomsIndependentAlleles(unsigned int k, unsigned long int n) {
        double prob = 0.0;
        unsigned int i;
        unsigned int pop = pow(2, k);

        for (i = n; i < pop + 1; ++i) {
            prob += binomialDistribution(pop, i, 0.25);
        }

        return prob;
    }
}