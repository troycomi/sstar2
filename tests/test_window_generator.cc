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
            vcf_str2 = (
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t"
                    "FORMAT\tmsp_0\tmsp_1\tmsp_2\tmsp_3\tmsp_4\tmsp_5\n"
                    "1\t6\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n"
                    "1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n"
                    "1\t9\t.\tAA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n"
                    "1\t10\t.\tA\tTT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n"
                    "1\t11\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n"
                    "1\t12\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n"
                    "1\t13\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|0\t0|0\t0|0\n"
                    "1\t14\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t1|1\t1|1\t./.\n"
                    "1\t15\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t0|0\t1|1\t./.\n"
                    "1\t16\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t0|0\t1|1\t./.\n"
                    "1\t17\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\t0|1\t0|0\t0|0\t0|0\t./.\n"
                    "1\t19\t.\tA\tT\t.\tPASS\t.\tGT:X\t"
                    "1|1:X\t1|1:X\t1|1:X\t0|0:X\t1|1:X\t1|0:X\n"
                    "1\t20\t.\tA\tT\t.\tPASS\t.\tGT:X:X\t"
                    "1|1:X:X\t1|1:X:X\t1|1:X:X\t0|0:X:X\t1|1:X:X\t1|0:X:X\n"
                    );
        }
        std::string vcf_str, vcf_str2, pop_str;
};

TEST_F(Generator_Input, StepCanInitialize){
    WindowGenerator gen(std::unique_ptr<Window>(new StepWindow(5, 10)));
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
    ASSERT_EQ(gen.vcf_line.genotypes.size(), 4);  // 4 indivs
} 

TEST_F(Generator_Input, StepCanInitializeOddWindow){
    WindowGenerator gen(std::unique_ptr<Window>(new StepWindow(3, 10)));
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
    ASSERT_EQ(gen.vcf_line.genotypes.size(), 5);  // 4 indivs
} 

TEST_F(Generator_Input, StepCanYieldWindow){
    // repeated calls to yield window update chrom, start and end
    WindowGenerator gen(std::unique_ptr<Window>(new StepWindow(5, 10)));
    std::set<std::string> target, reference, exclude;
    target.insert("EUR");
    reference.insert("AFR");
    exclude.insert("ASN");
    std::istringstream vcf(vcf_str);
    std::istringstream pop(pop_str);
    gen.initialize(vcf, pop, target, reference, exclude);

    ASSERT_TRUE(gen.next_window());
    ASSERT_STREQ(gen.window->chromosome.c_str(), "1");
    ASSERT_EQ(gen.window->start, 0);
    ASSERT_EQ(gen.window->end, 10);
    ASSERT_EQ(gen.window->total_snps(), 2);
    ASSERT_EQ(gen.window->reference_snps(), 1);
    ASSERT_EQ(gen.window->individual_snps(0), 1);

    ASSERT_TRUE(gen.next_window());
    ASSERT_STREQ(gen.window->chromosome.c_str(), "1");
    ASSERT_EQ(gen.window->start, 5);
    ASSERT_EQ(gen.window->end, 15);
    ASSERT_EQ(gen.window->total_snps(), 2);
    ASSERT_EQ(gen.window->reference_snps(), 1);
    ASSERT_EQ(gen.window->individual_snps(0), 1);

    ASSERT_FALSE(gen.next_window());
}

TEST_F(Generator_Input, StepCanYieldWindow2){
    // repeated calls to yield window update chrom, start and end
    WindowGenerator gen(std::unique_ptr<Window>(new StepWindow(2, 5)));
    std::set<std::string> target, reference, exclude;
    target.insert("msp_0");
    target.insert("msp_4");
    target.insert("msp_5");
    reference.insert("msp_1");
    reference.insert("msp_2");
    exclude.insert("msp_3");
    std::istringstream vcf(vcf_str2);
    std::istringstream pop(pop_str);
    gen.initialize(vcf, pop, target, reference, exclude);

    std::vector<int> total{0, 2, 2, 2, 3, 4, 4, 4, 3};
    std::vector<int> refs{0, 0, 0, 0, 0, 1, 3, 4, 3};
    std::vector<int> ind0{0, 2, 2, 2, 3, 3, 1, 0, 0};
    std::vector<int> ind1{0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<int> ind2{0, 0, 0, 0, 0, 0, 0, 0, 0};
    for(unsigned int i = 0; i < total.size(); ++i){
        ASSERT_TRUE(gen.next_window());
        ASSERT_STREQ(gen.window->chromosome.c_str(), "1");
        ASSERT_EQ(gen.window->start, 0 + i*2);
        ASSERT_EQ(gen.window->end, 5 + i*2);
        ASSERT_EQ(gen.window->total_snps(), total[i]);
        ASSERT_EQ(gen.window->reference_snps(), refs[i]);
        ASSERT_EQ(gen.window->individual_snps(0), ind0[i]);
        ASSERT_EQ(gen.window->individual_snps(1), ind1[i]);
        ASSERT_EQ(gen.window->individual_snps(2), ind2[i]);
    }

    ASSERT_FALSE(gen.next_window());
}

TEST_F(Generator_Input, StepCanYieldOddWindow){
    // repeated calls to yield window update chrom, start and end
    WindowGenerator gen(std::unique_ptr<Window>(new StepWindow(3, 10)));
    std::set<std::string> target, reference, exclude;
    target.insert("EUR");
    reference.insert("AFR");
    exclude.insert("ASN");
    std::istringstream vcf(vcf_str);
    std::istringstream pop(pop_str);
    gen.initialize(vcf, pop, target, reference, exclude);

    ASSERT_TRUE(gen.next_window());
    ASSERT_STREQ(gen.window->chromosome.c_str(), "1");
    ASSERT_EQ(gen.window->start, 0);
    ASSERT_EQ(gen.window->end, 10);
    ASSERT_EQ(gen.window->total_snps(), 2);
    ASSERT_EQ(gen.window->reference_snps(), 1);
    ASSERT_EQ(gen.window->individual_snps(0), 1);

    ASSERT_TRUE(gen.next_window());
    ASSERT_STREQ(gen.window->chromosome.c_str(), "1");
    ASSERT_EQ(gen.window->start, 3);
    ASSERT_EQ(gen.window->end, 13);
    ASSERT_EQ(gen.window->total_snps(), 2);
    ASSERT_EQ(gen.window->reference_snps(), 1);
    ASSERT_EQ(gen.window->individual_snps(0), 1);

    ASSERT_TRUE(gen.next_window());
    ASSERT_STREQ(gen.window->chromosome.c_str(), "1");
    ASSERT_EQ(gen.window->start, 6);
    ASSERT_EQ(gen.window->end, 16);
    ASSERT_EQ(gen.window->total_snps(), 1);
    ASSERT_EQ(gen.window->reference_snps(), 0);
    ASSERT_EQ(gen.window->individual_snps(0), 1);

    ASSERT_FALSE(gen.next_window());
}
