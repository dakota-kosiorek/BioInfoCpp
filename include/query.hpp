#ifndef QUERY_HPP
#define QUERY_HPP 1

#include <string>
#include <fundamentals.hpp>

#define UNIPROT_URL "http://www.uniprot.org/uniprot/"

namespace bioinfo {
    std::vector<AAString> queryUniProt(std::vector<std::string> &vec);
}

#endif