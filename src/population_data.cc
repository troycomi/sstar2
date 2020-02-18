#include "sstar2/population_data.h"

void PopulationData::read_data(std::istream &pop_file,
        std::set<std::string> &target,
        std::set<std::string> &reference,
        std::set<std::string> &exclude){
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
                target.find(superpop) != target.end()){
            targets.insert(indiv);
            target_to_population.insert({indiv, pop});
        }

        if(reference.find(indiv) != reference.end() ||
                reference.find(pop) != reference.end() ||
                reference.find(superpop) != reference.end())
            references.insert(indiv);

        if(exclude.find(indiv) != exclude.end() ||
                exclude.find(pop) != exclude.end() ||
                exclude.find(superpop) != exclude.end())
            excluded.insert(indiv);
    }
}
