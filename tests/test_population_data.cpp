#include <iostream>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "sstar2/population_data.hpp"

TEST(PopulationFile, CanInitialize){
    std::istringstream infile(
            "samp\tpop\tsuper_pop\n"
            "msp_0\tNeand1\tNeand1\n"
            "msp_1\tNeand2\tNeand2\n"
            "msp_2\tAFR\tAFR\n"
            "msp_3\tAFR\tAFR\n"
            "msp_4\tAFR\tAFR\n"
        );
    population_data pop;
    std::set<std::string> target, reference, exclude;
    pop.read_data(infile, target, reference, exclude);
}
