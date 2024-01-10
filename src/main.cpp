// Standard libs
#include <iostream>
#include <iomanip>
#include <vector>
// Bioinformatics libs
#include <fundamentals.hpp>
#include <genetics.hpp>
#include <analysis.hpp>
#include <biomath.hpp>
#include <query.hpp>

int main() {
    std::cout << std::fixed << std::setprecision(3);

    const std::string host = "uniprot.org";
    const std::string path = "/uniprot/A2Z669";

    std::string response = bioinfo::httpGet(host, path);

    if (!response.empty()) {
        std::cout << "Response:" << std::endl;
        std::cout << response << std::endl;
    }

    return 0;
}