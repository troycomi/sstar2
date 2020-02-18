#include <iostream>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "sstar2/window_generator.h"

using testing::ElementsAre;
using testing::Pair;

class Generator_Input : public ::testing::Test{
    protected:
        void SetUp(){
            vcf_str = (
                    "##fileformat=VCFv4.2\n"
                    "##source=msprime 0.6.1\n"
                    "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
                    "##contig=<ID=1,length=10000000>\n"
                    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
                    "\tmsp_0\tmsp_1\tmsp_2\tmsp_3\tmsp_4\tmsp_5\n"
                    "1\t1\t. \tA\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0\t0|0\n"
                    "1\t2\t. \tA\tT\t.\tPASS\t.\tGT\t1|1\t1|1\t1|1\t1|0\t1|1\t1|1\n"
                    "1\t5\t. \tA\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|1\t0|0\n"
                    "1\t6\t. \tA\tT\t.\tPASS\t.\tGT\t1|1\t1|1\t1|0\t1|1\t1|1\t0|0\n"
                    "1\t9\t. \tA\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t1|1\t1|1\t1|1\t0|0\n"
                    "1\t15\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|1\t0|0\n"
                    "1\t16\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|1\t0|1\n"
                    //                               ^no use^  ^  ref ^  ^tar ^ ex
                    );
            pop_str = (
                    "samp\tpop\tsuper_pop\n"
                    "msp_0\tNeand1\tNeand\n"
                    "msp_1\tNeand2\tNeand\n"
                    "msp_2\tAFR\tAFR\n"
                    "msp_3\tAFR\tAFR\n"
                    "msp_4\tEUR\tEUR\n"
                    "msp_5\tASN\tASN\n"
                   );
        }
        std::string vcf_str, pop_str;
};

TEST_F(Generator_Input, CanInitializeValidate){
    WindowGenerator gen(10, 10);  // same size ok
    ASSERT_THROW(WindowGenerator gen_bad(10, 11),
            std::invalid_argument);
}

TEST_F(Generator_Input, CanInitialize){
    WindowGenerator gen(10, 5);
    std::set<std::string> target, reference, exclude;
    target.insert("EUR");
    reference.insert("AFR");
    exclude.insert("ASN");
    std::istringstream vcf(vcf_str);
    std::istringstream pop(pop_str);
    gen.initialize(vcf, pop, target, reference, exclude);
    ASSERT_THAT(gen.population.targets, ElementsAre("msp_4"));
    ASSERT_THAT(gen.population.references, ElementsAre("msp_2", "msp_3"));
    ASSERT_THAT(gen.population.excluded, ElementsAre("msp_5"));
    ASSERT_THAT(gen.population.target_to_population,
            ElementsAre(Pair("msp_4", "EUR")));
    ASSERT_STREQ(gen.vcf_line.chromosome.c_str(), "1");  // first line
    ASSERT_EQ(gen.vcf_line.haplotypes.size(), 4*2);  // 4 indivs, 2 haplotypes
    // two buckets for 10 / 5
    ASSERT_EQ(gen.window.buckets.size(), 2);
    // one target
    for (auto b : gen.window.buckets)
        ASSERT_EQ(b.genotypes.size(), 1);
} 

TEST_F(Generator_Input, CanInitializeOddWindow){
    WindowGenerator gen(10, 3);
    std::set<std::string> target, reference, exclude;
    target.insert("EUR");
    target.insert("msp_0");
    target.insert("msp_2");
    reference.insert("AFR");
    exclude.insert("ASN");
    std::istringstream vcf(vcf_str);
    std::istringstream pop(pop_str);
    gen.initialize(vcf, pop, target, reference, exclude);
    ASSERT_THAT(gen.population.targets, ElementsAre("msp_0", "msp_2", "msp_4"));
    ASSERT_THAT(gen.population.references, ElementsAre("msp_2", "msp_3"));
    ASSERT_THAT(gen.population.excluded, ElementsAre("msp_5"));
    ASSERT_THAT(gen.population.target_to_population,
            ElementsAre(Pair("msp_0", "Neand1"),
                Pair("msp_2", "AFR"),
                Pair("msp_4", "EUR")
                ));
    ASSERT_STREQ(gen.vcf_line.chromosome.c_str(), "1");
    ASSERT_EQ(gen.vcf_line.haplotypes.size(), 5*2);  // 4 indivs, 2 haplotypes
    // four buckets for 10 / 3
    ASSERT_EQ(gen.window.buckets.size(), 4);
    // three target
    for (auto b : gen.window.buckets)
        ASSERT_EQ(b.genotypes.size(), 3);
} 

TEST_F(Generator_Input, CanYieldWindow){
    // repeated calls to yield window update chrom, start and end
    WindowGenerator gen(10, 5);
    std::set<std::string> target, reference, exclude;
    target.insert("EUR");
    reference.insert("AFR");
    exclude.insert("ASN");
    std::istringstream vcf(vcf_str);
    std::istringstream pop(pop_str);
    gen.initialize(vcf, pop, target, reference, exclude);

    ASSERT_TRUE(gen.next_window());
    ASSERT_STREQ(gen.window.chromosome.c_str(), "1");
    ASSERT_EQ(gen.window.start, 0);
    ASSERT_EQ(gen.window.end, 10);
    std::cout << gen.window;

    ASSERT_TRUE(gen.next_window());
    ASSERT_STREQ(gen.window.chromosome.c_str(), "1");
    ASSERT_EQ(gen.window.start, 5);
    ASSERT_EQ(gen.window.end, 15);
    std::cout << gen.window;

    ASSERT_FALSE(gen.next_window());
}
