#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <iostream>

#include "sstar2/sstar.h"

using testing::ElementsAre;

TEST(SSTAR, CanWriteHeader){
    std::ostringstream outfile;
    SStarCaller sstar;
    sstar.write_header(outfile);
    ASSERT_STREQ(outfile.str().c_str(),
            "chrom\twinstart\twinend\tn_snps\tn_ind_snps\tn_region_ind_snps\t"
            "ind_id\tpop\ts_star\tnum_s_star_snps\ts_star_snps\thap_1_s_start\t"
            "hap_1_s_end\thap_2_s_start\thap_2_s_end\ts_start\ts_end\t"
            "n_s_star_snps_hap1\tn_s_star_snps_hap2\ts_star_haps\t"
            "callable_bases\n");
}

class SStarFixtureEmpty : public ::testing::Test{
    protected:
        void SetUp(){
            std::set<std::string> target, reference, exclude;
            target.insert("targ");
            reference.insert("ref");
            std::string vcf_str(
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t"
                    "FORMAT\tmsp_0\tmsp_1\tmsp_2\tmsp_3\tmsp_4\tmsp_5\tref\n"
                    "1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t1|0\t0|0\t1|1\t1|1\t0|0\t0|0\n"
                    );
            std::string pop_str(
                    "samp\tpop\tsuper_pop\n"
                    "msp_0\tpop0\ttarg\n"
                    "msp_1\tpop1\ttarg\n"
                    "msp_2\tpop2\ttarg\n"
                    "msp_3\tpop3\ttarg\n"
                    "msp_4\tpop4\ttarg\n"
                    "msp_5\tpop5\ttarg\n"
                    "ref\t.\tref\n"
                    );
            std::istringstream vcf(vcf_str);
            std::istringstream pop(pop_str);
            generator.initialize(vcf, pop, target, reference, exclude);
            generator.next_window();
        }
        WindowGenerator generator{5, 2};
        SStarCaller sstar;
};

class SStarFixtureNormal : public ::testing::Test{
    protected:
        void SetUp(){
            std::set<std::string> target, reference, exclude;
            target.insert("targ");
            reference.insert("ref");
            std::string vcf_str(
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t"
                    "FORMAT\tmsp_0\tmsp_1\tmsp_2\tmsp_3\tmsp_4\tmsp_5\tref\n"
                    "1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t1|0\t0|0\t1|1\t1|1\t0|0\t0|0\n"
                    "1\t15\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t1|1\t0|0\t0|0\t0|0\t0|0\t0|0\n"
                    "1\t30\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|1\t0|0\t0|0\t0|0\t0|0\t0|0\n"
                    "1\t45\t.\tA\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|0\t0|1\t0|0\t0|0\t0|0\n"
                    "1\t5382\t.\tA\tT\t.\tPASS\t. \tGT\t0|0\t0|0\t1|0\t0|0\t0|0\t0|0\t1|0\n"
                    "1\t28610\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|1\t0|0\t0|0\t0|0\t0|0\n"
                    "1\t32662\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t1|0\t0|0\t0|0\t0|0\t1|0\n"
                    "1\t35985\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t1|1\t0|0\t0|0\t0|0\t0|0\n"
                    );
            std::string pop_str(
                    "samp\tpop\tsuper_pop\n"
                    "msp_0\tpop0\ttarg\n"
                    "msp_1\tpop1\tnothing\n"
                    "msp_2\tpop2\ttarg\n"
                    "msp_3\tpop3\ttarg\n"
                    "msp_4\tpop4\ttarg\n"
                    "msp_5\tpop5\ttarg\n"
                    "ref\t.\tref\n"
                    );
            std::istringstream vcf(vcf_str);
            std::istringstream pop(pop_str);
            generator.initialize(vcf, pop, target, reference, exclude);
            generator.next_window();
        }
        WindowGenerator generator{50000, 10000};
        SStarCaller sstar;
};

TEST_F(SStarFixtureEmpty, CanWriteEmptyWindow){
    std::ostringstream outfile;
    sstar.write_window(outfile, generator);
    ASSERT_STREQ(outfile.str().c_str(),
            // chrom, start end, snps, indiv and pop
            "1\t0\t5\t0\t0\t0\tmsp_0\tpop0\t"
            // sstar related stuff, callable bases
            "0\t0\t.\t0\t0\t0\t0\t0\t0\t0\t0\t.\t5\n"
            "1\t0\t5\t0\t0\t0\tmsp_1\tpop1\t"
            "0\t0\t.\t0\t0\t0\t0\t0\t0\t0\t0\t.\t5\n"
            "1\t0\t5\t0\t0\t0\tmsp_2\tpop2\t"
            "0\t0\t.\t0\t0\t0\t0\t0\t0\t0\t0\t.\t5\n"
            "1\t0\t5\t0\t0\t0\tmsp_3\tpop3\t"
            "0\t0\t.\t0\t0\t0\t0\t0\t0\t0\t0\t.\t5\n"
            "1\t0\t5\t0\t0\t0\tmsp_4\tpop4\t"
            "0\t0\t.\t0\t0\t0\t0\t0\t0\t0\t0\t.\t5\n"
            "1\t0\t5\t0\t0\t0\tmsp_5\tpop5\t"
            "0\t0\t.\t0\t0\t0\t0\t0\t0\t0\t0\t.\t5\n"
            );
}

TEST_F(SStarFixtureNormal, CanWriteWindow){
    std::ostringstream outfile;
    sstar.write_window(outfile, generator);
    ASSERT_STREQ(outfile.str().c_str(),
            // chrom, start end, snps, indiv and pop
            "1\t0\t50000\t8\t4\t6\tmsp_0\tpop0\t"
            // sstar related stuff, callable bases
            "10035\t3\t10,30,45\t10\t45\t0\t0\t10\t45\t3\t0\t1,1,1\t50000\n"
            "1\t0\t50000\t8\t2\t4\tmsp_2\tpop2\t"
            "0\t0\t.\t0\t0\t0\t0\t0\t0\t0\t0\t.\t50000\n"
            "1\t0\t50000\t8\t2\t4\tmsp_3\tpop3\t"
            "0\t0\t.\t0\t0\t0\t0\t0\t0\t0\t0\t.\t50000\n"
            "1\t0\t50000\t8\t1\t3\tmsp_4\tpop4\t"
            "0\t0\t.\t0\t0\t0\t0\t0\t0\t0\t0\t.\t50000\n"
            "1\t0\t50000\t8\t0\t2\tmsp_5\tpop5\t"
            "0\t0\t.\t0\t0\t0\t0\t0\t0\t0\t0\t.\t50000\n"
            );
    ASSERT_FALSE(generator.next_window());
}

TEST(SStarMismatches, TiedScoreUseLonger){
    std::vector<WindowGT>genotypes{
        {7462931, 3},
        {7467931, 3},
        {7468819, 1},
        {7492334, 2},
        {7493191, 2}};
    SStarCaller sstar;
    ASSERT_EQ(sstar.sstar(genotypes), 34372);
    for (auto &gt : genotypes)
        std::cout << gt.position << ',' << gt.genotype << ' ';
    std::cout << '\n';
    ASSERT_THAT(genotypes, ElementsAre(
                WindowGT{7462931, 3},
                WindowGT{7467931, 3},
                WindowGT{7468819, 1},
                WindowGT{7492334, 2},
                WindowGT{7493191, 2}
        ));
}
