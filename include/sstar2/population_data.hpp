#pragma once
#include <iostream>
#include <sstream>
#include <set>
#include <map>
#include <string>

class population_data{
    public:
        std::set<std::string> targets, reference, exclude;
        std::map<std::string, std::string> target_to_population;

        void read_data(std::istream &pop_file,
                std::set<std::string> target,
                std::set<std::string> reference,
                std::set<std::string> exclude);
};
