#include <iostream>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "sstar2/vcf_file.h"

TEST(VCF_File, CanInitialize){
    VcfFile f;
    std::set<std::string> indivs = {"msp_1", "msp_3", "msp_5", "msp_00"};
    // 3 individuals match
    ASSERT_EQ(f.initialize_individuals(
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                "msp_0\tmsp_1\tmsp_2\tmsp_3\tmsp_4\tmsp_5",
                indivs), 3);

    ASSERT_EQ(f.individual_map.size(), 3);
    ASSERT_STREQ(f.individual_map[10].c_str(), "msp_1");
    ASSERT_STREQ(f.individual_map[12].c_str(), "msp_3");
    ASSERT_STREQ(f.individual_map[14].c_str(), "msp_5");
}

TEST(VCF_File, CanInitializeNoOverlap){
    VcfFile f;
    std::set<std::string> indivs = {"ms_1", "ms_3", "ms_5", "ms_00"};
    // 3 individuals match
    ASSERT_EQ(f.initialize_individuals(
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                "msp_0\tmsp_1\tmsp_2\tmsp_3\tmsp_4\tmsp_5",
                indivs), 0);

    ASSERT_EQ(f.individual_map.size(), 0);
}

TEST(VCF_File, CanInitializeEmptyIndivs){
    VcfFile f;
    std::set<std::string> indivs = {};
    // 3 individuals match
    ASSERT_EQ(f.initialize_individuals(
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                "msp_0\tmsp_1\tmsp_2\tmsp_3\tmsp_4\tmsp_5",
                indivs), 6);

    ASSERT_EQ(f.individual_map.size(), 6);
}

class VCF_File_F : public ::testing::Test{
    protected:
        void SetUp() override {
            std::set<std::string> indivs = {"msp_1", "msp_3", "msp_5", "msp_00"};
            vcf.initialize_individuals(
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                    "msp_0\tmsp_1\tmsp_2\tmsp_3\tmsp_4\tmsp_5",
                    indivs);
        }

        VcfFile vcf;
};

TEST_F(VCF_File_F, CanInitializeEntry){
    VcfFile empty;
    VcfEntry e = empty.initialize_entry();
    ASSERT_STREQ(e.chromosome.c_str(), "");
    ASSERT_EQ(e.genotypes.size(), 0);

    e = vcf.initialize_entry();
    ASSERT_STREQ(e.chromosome.c_str(), "");
    ASSERT_EQ(e.genotypes.size(), 3);  // 3 match
}

TEST_F(VCF_File_F, CanParseLine){
    VcfEntry entry = vcf.initialize_entry();
    // no execption parse
    ASSERT_TRUE(vcf.parse_line(
            "1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0\t0|0",
            entry));
    ASSERT_STREQ(entry.chromosome.c_str(), "1");
    ASSERT_EQ(entry.position, 7);
    ASSERT_EQ(entry.reference, 'A');
    ASSERT_EQ(entry.alternative, 'T');
    ASSERT_THAT(entry.genotypes, ::testing::ElementsAre(0, 0, 0));

    // with haplotypes
    ASSERT_TRUE(vcf.parse_line(
            "2\t8\t.\tC\tG\t.\tPASS\t.\tGT\t1|0\t1|0\t0|1\t0|1\t0|0\t1|1",
            entry));
    ASSERT_STREQ(entry.chromosome.c_str(), "2");
    ASSERT_EQ(entry.position, 8);
    ASSERT_EQ(entry.reference, 'C');
    ASSERT_EQ(entry.alternative, 'G');
    ASSERT_THAT(entry.genotypes, ::testing::ElementsAre(1, 2, 3));

    // with other vcf entries, and dots
    ASSERT_TRUE(vcf.parse_line(
            "2\t8\t.\tC\tG\t.\tPASS\t.\tGT:and:others\t1|0\t1|0:asdf:asdf"
            "\t0|1\t.\t1|0\t./.",
            entry));
    ASSERT_STREQ(entry.chromosome.c_str(), "2");
    ASSERT_EQ(entry.position, 8);
    ASSERT_EQ(entry.reference, 'C');
    ASSERT_EQ(entry.alternative, 'G');
    ASSERT_THAT(entry.genotypes, ::testing::ElementsAre(1, 0, 0));

    // failing format
    ASSERT_THROW(vcf.parse_line(
            "2\t8\t.\tC\tG\t.\tPASS\t.\tnot:GT\t1|0\t1|0:asdf:asdf"
            "\t0|1\t.\t1|0\t./.",
            entry), std::invalid_argument);

    // skip for biallelic
    ASSERT_FALSE(vcf.parse_line(
            "2\t8\t.\tCA\tG\t.\tPASS\t.\tGT:and:others\t1|0\t1|0:asdf:asdf"
            "\t0|1\t.\t1|0\t./.",
            entry));
    ASSERT_FALSE(vcf.parse_line(
            "2\t8\t.\tC\tAG\t.\tPASS\t.\tGT:and:others\t1|0\t1|0:asdf:asdf"
            "\t0|1\t.\t1|0\t./.",
            entry));
    ASSERT_FALSE(vcf.parse_line(
            "2\t8\t.\tCT\tAG\t.\tPASS\t.\tGT:and:others\t1|0\t1|0:asdf:asdf"
            "\t0|1\t.\t1|0\t./.",
            entry));
}

TEST_F(VCF_File_F, ParseLineUnphasedMakesWarning){
    VcfEntry entry = vcf.initialize_entry();
    // warn for unphased haplotypes
    testing::internal::CaptureStderr();
    ASSERT_TRUE(vcf.parse_line(
            "3\t10\t.\tC\tG\t.\tPASS\t.\tGT:and:others\t1/0\t1/0:asdf:asdf"
            "\t0/1\t.\t1/0\t./.",
            entry));
    ASSERT_STREQ(entry.chromosome.c_str(), "3");
    ASSERT_EQ(entry.position, 10);
    ASSERT_THAT(entry.genotypes, ::testing::ElementsAre(1, 0, 0));
    std::string output = testing::internal::GetCapturedStderr();
    ASSERT_STREQ(output.c_str(),
            "WARNING: Detected unphased haplotype at chrom 3 and pos 10!\n");

    // but only once
    testing::internal::CaptureStderr();
    ASSERT_TRUE(vcf.parse_line(
            "2\t8\t.\tC\tG\t.\tPASS\t.\tGT:and:others\t1/0\t1|0:asdf:asdf"
            "\t0|1\t.\t1|0\t./.",
            entry));
    ASSERT_THAT(entry.genotypes, ::testing::ElementsAre(1, 0, 0));
    ASSERT_STREQ(entry.chromosome.c_str(), "2");
    ASSERT_EQ(entry.position, 8);
    output = testing::internal::GetCapturedStderr();
    ASSERT_STREQ(output.c_str(), "");
}
