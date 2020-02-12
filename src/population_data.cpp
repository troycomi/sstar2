#include "sstar2/population_data.hpp"

void population_data::read_data(std::istream &pop_file,
        std::set<std::string> target,
        std::set<std::string> reference,
        std::set<std::string> exclude){
    // read in 3 column data from pop_file and place entries into appropriate
    // public set based on provided arguments.  Target indivs also go into
    // target to population map
    std::string line;
    std::string indiv, pop, superpop;
    while(std::getline(pop_file, line)){
        std::istringstream iss(line);
        if (!(iss >> indiv >> pop >> superpop))
            break;
        
        if(target.find(indiv) != target.end() ||
                target.find(pop) != target.end() ||
                target.find(superpop) != target.end())
            std::cout << "";
    }
}
