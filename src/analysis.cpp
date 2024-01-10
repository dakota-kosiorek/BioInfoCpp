#include <analysis.hpp>
#include <fundamentals.hpp>
#include <string>
#include <iostream>
#include <stdexcept>
#include <sstream>

namespace bioinfo {
    // Create a new `AdjacencyList` object that generates the overlap graph of a vector of DNAStrings and 
    // stores the headers of each directed edge that make up the graph.
    AdjacencyList::AdjacencyList(std::vector<DNAString> &vec, unsigned int ok) {
        (*this).internalConstructor(vec, ok);
    }

    template <typename T> void AdjacencyList::internalConstructor(std::vector<T> &vec, unsigned int ok) {
        std::vector<DNAString>::iterator outer;
        std::vector<DNAString>::iterator inner;

        DirectedEdge dsde;
        std::string prefix;
        std::string suffix;

        for (outer = vec.begin(); outer != vec.end(); outer++) {
            for (inner = vec.begin(); inner != vec.end(); inner++) {
                if (outer != inner && outer->getSequenceLength() > ok && inner->getSequenceLength() > ok) {

                    suffix = outer->getSequence().substr(outer->getSequenceLength() - ok);
                    prefix = inner->getSequence().substr(0, ok);

                    if (prefix == suffix) {
                        dsde.tail = outer->getHeader();
                        dsde.head = inner->getHeader();

                        (*this).dsde.push_back(dsde);
                    }
                }
            }
        }
    }

    // Return the `AdjacencyList` directed edge vector as a string.
    std::string AdjacencyList::toString() {
        std::stringstream ss;

        std::vector<DirectedEdge>::iterator it;

        for (it = (*this).dsde.begin(); it != (*this).dsde.end(); it++) {
            ss << it->tail << " " << it->head;

            if(it + 1 != (*this).dsde.end()) {
                ss << "\n";
            }
        }

        return ss.str();
    }

    // Get the indices in the DNAString `ds` where a DNA motif (`motif`) is found .
    std::vector<unsigned int> exactDNAStringMotif(DNAString &ds, DNAString &motif, bool overlap) {
        std::vector<unsigned int> positions;
        unsigned int pos = 0;

        if (motif.getSequenceLength() > ds.getSequenceLength()) {
            // throw error that motif is greater than ds
        } else if (motif.getSequenceLength() == 0) {
            // throw error that motif is of size 0
        } else {
            while ((pos = ds.getSequence().find(motif.getSequence(), pos)) < ds.getSequenceLength()) {
                positions.push_back(pos);

                if (overlap) {
                    // With overlap
                    pos += 1;
                } else {
                    // Without overlap
                    pos += motif.getSequenceLength();
                }
            }
        }

        return positions;
    }

    // Get the hamming distance between the sequences of two DNAStrings (`s` and `t`).
    unsigned int hammingDistance(DNAString &s, DNAString &t) {
        unsigned int hd = 0;
        unsigned int i = 0;

        if (s.getSequenceLength() != t.getSequenceLength()) {
            throw std::invalid_argument("ERROR: DNAString sequences are not the same length!");
        }

        for (i = 0; i < s.getSequenceLength(); ++i) {
            if (s.getSequence().at(i) != t.getSequence().at(i)) {
                ++hd;
            }
        }

        return hd;
    }

    // Calculate the total mass of a protein `as` in daltons based on a mass table `mt`.
    double proteinMass(bioinfo::AAString &as, const MassTable &mt) {
        double pm = 0.0;
        unsigned int i;
        unsigned int asLength = as.getSequenceLength();
        MassTable::const_iterator it;
        char currAA;

        for (i = 0; i < asLength; ++i) {
            currAA = as.getSequence().at(i);
            it = mt.find(currAA);

            if (it != mt.end()) {
                pm += mt.at(currAA);
            }
        }

        return pm;
    }

    // Calculate how many possible mRNA strands an inputted protein sequence `as` could of come from applied with the modulus 
    // operator at a value of `m`.
    unsigned int inferredRNACount(AAString &as, const AATranscribableUnitTable &ut, unsigned int m) {
        unsigned int i;
        unsigned int asLength = as.getSequenceLength();
        unsigned int result = 0;

        AATranscribableUnitTable::const_iterator it;
        unsigned char tUnits;
        char currAA;

        if (asLength > 0) {
            for (i = 1; i < asLength; ++i) {
                currAA = as.getSequence().at(i);
                it = ut.find(currAA);

                if (it != ut.end() && result != 0) {
                    result = (result * ut.at(currAA)) % m;
                } else if (it != ut.end() && result == 0) {
                    result = ut.at(currAA);
                }
            }
        }

        return result;
    }

    // Remove the introns of a RNAString `s` given by a vector of RNAString objects `introns`. The sequence of each item in the 
    // vector is used as the intron sequence and the first occurrence of each sequence is removed from the original RNA sequence
    // given by `s`.
    RNAString spliceRNA(RNAString &s, std::vector<RNAString> &introns) {
        RNAString mrna = RNAString(s.getHeader(), std::string(""));
        std::string seq = s.getSequence();

        std::vector<RNAString>::iterator it;
        unsigned int loc;

        // For every intron loop through pre-mrna sequence and remove the intron
        for (it = introns.begin(); it != introns.end(); it++) {
            loc = seq.find(it->getSequence());

            if (loc != std::string::npos) {
                seq.erase(loc, it->getSequenceLength());
            }
        }

        mrna.setSequence(seq);
        return mrna;
    }
}