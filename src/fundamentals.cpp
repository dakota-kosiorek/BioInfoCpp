#include <fundamentals.hpp>
#include <string>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <fstream>

namespace bioinfo {
    // Create a new empty DNAString
    DNAString::DNAString() {
        (*this).header = "";
        (*this).sequence = "";
        (*this).sequenceLength = 0;
    }

    // Create a new DNAString with a header and sequence
    DNAString::DNAString(std::string h, std::string s) {
        (*this).header = h;
        (*this).sequence = toUpper(s);
        (*this).sequenceLength = s.length();
    }

    // Get the header of the DNAString
    std::string DNAString::getHeader() {
        return (*this).header;
    }

    // Get the sequence of the DNAString
    std::string DNAString::getSequence() {
        return (*this).sequence;
    }

    // Get the how many nucleotides are in the DNAString
    unsigned int DNAString::getSequenceLength() {
        return (*this).sequenceLength;
    }

    // Change the header of the DNAString
    void DNAString::setHeader(std::string h) {
        (*this).sequence = h;
    }

    // Change the seuqence of the DNAString
    void DNAString::setSequence(std::string s) {
        (*this).sequence = toUpper(s);
        (*this).sequenceLength = s.length();
    }

    // --------------------------------------------------------------------------

    // Create a new empty RNAString
    RNAString::RNAString() {
        (*this).header = "";
        (*this).sequence = "";
        (*this).sequenceLength = 0;
    }

    // Create a new RNAString with a header and sequence
    RNAString::RNAString(std::string h, std::string s) {
        (*this).header = h;
        (*this).sequence = transcribe(s);
        (*this).sequenceLength = s.length();
    }

    RNAString::RNAString(DNAString &ds) {
        (*this).header = ds.getHeader();
        (*this).sequence = transcribe(ds.getSequence());
        (*this).sequenceLength = ds.getSequenceLength();
    }

    // Get the header of the RNAString
    std::string RNAString::getHeader() {
        return (*this).header;
    }

    // Get the sequence of the RNAString
    std::string RNAString::getSequence() {
        return (*this).sequence;
    }

    // Get the how many nucleotides are in the RNAString
    unsigned int RNAString::getSequenceLength() {
        return (*this).sequenceLength;
    }

    // Change the header of the RNAString
    void RNAString::setHeader(std::string h) {
        (*this).sequence = h;
    }

    // Change the seuqence of the RNAString
    void RNAString::setSequence(std::string s) {
        (*this).sequence = transcribe(s);
        (*this).sequenceLength = s.length();
    }

    // --------------------------------------------------------------------------

    // Create a new AAString with header `h`, sequence `s`, and using a certain genetic code table (`code`)
    AAString::AAString(std::string h, std::string s, const AATable &code) {
        (*this).header = h;
        (*this).sequence = toUpper(s);
        (*this).sequenceLength = (*this).sequence.length();
    }

    // Create a new AAString from a DNAstring object (`ds`) using a certain genetic code table (`code`)
    AAString::AAString(DNAString &ds, const AATable &code) {
        (*this).header = ds.getHeader();
        (*this).sequence = translate(transcribe(ds.getSequence()), code);
        (*this).sequenceLength = (*this).sequence.length();
    }

    // Create a new AAString from a RNAstring object (`ds`) using a certain genetic code table (`code`)
    AAString::AAString(RNAString &rs, const AATable &code) {
        (*this).header = rs.getHeader();
        (*this).sequence = translate(rs.getSequence(), code);
        (*this).sequenceLength = (*this).sequence.length();
    }

    // Get the header of the AAString
    std::string AAString::getHeader() {
        return (*this).header;
    }

    // Get the sequence of the AAString
    std::string AAString::getSequence() {
        return (*this).sequence;
    }

    // Get how many amino acids are in the sequence (including stop codons)
    unsigned int AAString::getSequenceLength() {
        return (*this).sequenceLength;
    }

    // Change the header of the AAString
    void AAString::setHeader(std::string h) {
        (*this).sequence = h;
    }

    // Change the seuqence of the AAString
    void AAString::setSequence(std::string s, const AATable &code) {
        (*this).sequence = translate(transcribe(s), code);
        (*this).sequenceLength = (*this).sequence.length();
    }

    // --------------------------------------------------------------------------

    // Return a version of the string `s` with all uppercase letters
    std::string toUpper(std::string s) {
        unsigned int i;
        unsigned int sLength = s.length();

        for (i = 0; i < sLength; ++i) {
            s.at(i) = std::toupper(s.at(i));
        }

        return s;
    }

    // Transcribe a string of DNA to RNA
    std::string transcribe(std::string s) {
        unsigned int i;
        unsigned int sLength = s.length();

        for (i = 0; i < sLength; ++i) {
            s.at(i) = std::toupper(s.at(i));

            if (s.at(i) == 'T') {
                s.at(i) = 'U';
            }
        }

        return s;
    }

    // Translate a string of RNA to AA
    std::string translate(std::string s, const AATable &code) {
        unsigned int i;
        unsigned int sLength = s.length();
        AATable::const_iterator it;
        std::string window = "";
        std::string newSeq = "";
        
        if (sLength > 2) {
            for (i = 0; i < sLength-2; i += 3) {
                s.at(i) = std::toupper(s.at(i));
                s.at(i+1) = std::toupper(s.at(i+1));
                s.at(i+2) = std::toupper(s.at(i+2));

                window = s.substr(i, 3);
                it = code.find(window);

                if (it != code.end()) {
                    newSeq += code.at(window);
                } else {
                    newSeq += "X";
                }
            }
        }

        return newSeq;
    }

    // Reverse complement the sequence of a DNAString
    DNAString reverseComplement(DNAString &ds) {
        DNAString rc = DNAString(ds.getHeader(), "");
        std::string s = "";

        unsigned int i;
        char oldC;
        char newC;

        for (i = 0; i < ds.getSequenceLength(); ++i) {
            oldC = ds.getSequence().at(ds.getSequenceLength() - i - 1);

            switch (oldC) {
                case 'A':
                    newC = 'T';
                    break;
                case 'C':
                    newC = 'G';
                    break;
                case 'G':
                    newC = 'C';
                    break;
                case 'T':
                    newC = 'A';
                    break;
                default:
                    newC = 'N';
            }

            s += newC;
        }

        rc.setSequence(s);
        return rc;
    }

    // Read a FASTA file with name `fn` and return a vector of DNAString objects.
    std::vector<DNAString> readDNAStringFile(std::string &fn) {
        std::vector<DNAString> vec;
        std::ifstream file(fn);

        std::string txt;
        std::string header = "";
        std::string seq = "";

        if (!file.good()) {
            throw std::invalid_argument("ERROR: readDNAStringFile could not open file!");
        } else {
            while (std::getline(file, txt)) {
                if (txt.empty()) {
                    continue;
                } else if (txt.at(0) == '>') {
                    if (!header.empty()) {
                        DNAString ds = DNAString(header, seq);
                        vec.push_back(ds);
                    }

                    header = txt.substr(1);
                    while (header.back() == '\n' || header.back() == '\r') {
                        header.pop_back();
                    }
                    seq = "";
                } else {
                    seq += txt;
                    while (seq.back() == '\n' || seq.back() == '\r') {
                        seq.pop_back();
                    }
                }
            }

            if (!header.empty()) {
                DNAString ds = DNAString(header, seq);
                vec.push_back(ds);
            }
        }

        file.close();

        return vec;
    }

    // Read a FASTA file with name `fn` and return a vector of DNAString objects.
    std::vector<DNAString> readDNAStringFile(const char *fnp) {
        std::string fn = std::string(fnp);
        return readDNAStringFile(fn);
    }
}