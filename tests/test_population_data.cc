#include <iostream>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "sstar2/population_data.h"

using testing::Pair;
using testing::ElementsAre;

TEST(PopulationFile, CanInitialize){
    PopulationData pop;

    ASSERT_TRUE(pop.targets.empty());
    ASSERT_TRUE(pop.references.empty());
    ASSERT_TRUE(pop.excluded.empty());
    ASSERT_TRUE(pop.target_to_population.empty());

}

TEST(PopulationFile, ReadDataNoSets){
    std::istringstream infile(
            "samp\tpop\tsuper_pop\n"
            "msp_0\tNeand1\tNeand1\n"
            "msp_1\tNeand2\tNeand2\n"
            "msp_2\tAFR\tAFR2\n"
            "msp_3\tAFR1\tAFR\n"
            "msp_4\tEUR\tEUR\n"
        );
    PopulationData pop;
    std::set<std::string> target, reference, exclude;
    pop.read_data(infile, target, reference, exclude);

    ASSERT_TRUE(pop.targets.empty());
    ASSERT_TRUE(pop.references.empty());
    ASSERT_TRUE(pop.excluded.empty());
    ASSERT_TRUE(pop.target_to_population.empty());
}

TEST(PopulationFile, ReadDataOneMatch){
    std::istringstream infile(
            "samp\tpop\tsuper_pop\n"
            "msp_0\tNeand1\tNeand\n"
            "msp_1\tNeand2\tNeand\n"
            "msp_2\tAFR\tAFR2\n"
            "msp_3\tAFR1\tAFR\n"
            "msp_4\tEUR\tEUR\n"
        );
    PopulationData pop;
    std::set<std::string> target, reference, exclude;
    target.insert("msp_1");  // one indiv
    reference.insert("Neand1");  // one pop
    exclude.insert("AFR2");  // one superpop
    pop.read_data(infile, target, reference, exclude);

    ASSERT_THAT(pop.targets, ElementsAre("msp_1"));
    ASSERT_THAT(pop.references, ElementsAre("msp_0"));
    ASSERT_THAT(pop.excluded, ElementsAre("msp_2"));
    ASSERT_THAT(pop.target_to_population,
            ElementsAre(Pair("msp_1", "Neand2")));
}

TEST(PopulationFile, ReadDataMatchMany){
    std::istringstream infile(
            "samp\tpop\tsuper_pop\n"
            "msp_0\tNeand1\tNeand\n"
            "msp_1\tNeand2\tNeand\n"
            "msp_2\tAFR\tAFR2\n"
            "msp_3\tAFR1\tAFR\n"
            "msp_4\tEUR\tEUR\n"
        );
    PopulationData pop;
    std::set<std::string> target, reference, exclude;
    target.insert("Neand");  // two super pop
    reference.insert("arf");  // nothing
    exclude.insert("AFR");   // one pop, one superpop
    exclude.insert("msp_4"); 
    exclude.insert("Neand2"); 
    exclude.insert("Neand");  // overlap target
    pop.read_data(infile, target, reference, exclude);

    ASSERT_THAT(pop.targets, ElementsAre("msp_0", "msp_1"));
    ASSERT_TRUE(pop.references.empty());
    ASSERT_THAT(pop.excluded, ElementsAre("msp_0", "msp_1",
                "msp_2", "msp_3", "msp_4"));
    ASSERT_THAT(pop.target_to_population,
            ElementsAre(Pair("msp_0", "Neand1"), Pair("msp_1", "Neand2")));
}
