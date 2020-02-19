#include <iostream>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "sstar2/vcf_file.h"

class VCF_Entry_Test : public ::testing::Test{
    protected:
        void SetUp() override {
            entry_5.chromosome = "1";
            entry_5.position = 10;
            entry_5.reference = 'A';
            entry_5.alternative = 'T';
            entry_5.haplotypes[3] = true;
            entry_5.haplotypes[4] = true;
            entry_5.haplotypes[5] = true;
            entry_5.haplotypes[6] = true;

            entry_0.chromosome = "2";
            entry_0.position = 11;
            entry_0.reference = 'G';
            entry_0.alternative = 'C';
        }

        VcfEntry entry_5{"1", 5};
        VcfEntry entry_0{"1", 0};

};

TEST_F(VCF_Entry_Test, CanInitialize){
    ASSERT_STREQ(entry_5.to_str().c_str(),
            "1\t10\tA\tT\t0\t0\t0\t1\t1\t1\t1\t0\t0\t0\n");
    ASSERT_STREQ(entry_0.to_str().c_str(),
            "2\t11\tG\tC\n");
}

TEST_F(VCF_Entry_Test, CanGetGenotype){
    ASSERT_EQ(entry_5.genotype(0), 0);
    ASSERT_EQ(entry_5.genotype(1), 1);
    ASSERT_EQ(entry_5.genotype(2), 2);
    ASSERT_EQ(entry_5.genotype(3), 1);
    ASSERT_EQ(entry_5.genotype(4), 0);
    ASSERT_THROW(entry_5.genotype(5), std::out_of_range);
    ASSERT_THROW(entry_0.genotype(0), std::out_of_range);
}

TEST_F(VCF_Entry_Test, CanGetGTCode){
    ASSERT_EQ(entry_5.gt_code(0), 0);
    ASSERT_EQ(entry_5.gt_code(1), 2);
    ASSERT_EQ(entry_5.gt_code(2), 3);
    ASSERT_EQ(entry_5.gt_code(3), 1);
    ASSERT_EQ(entry_5.gt_code(4), 0);
}

TEST_F(VCF_Entry_Test, CanFindAnyHaplotype){
    ASSERT_FALSE(entry_5.any_haplotype(std::vector<unsigned int> {0}));
    ASSERT_TRUE(entry_5.any_haplotype(std::vector<unsigned int> {0, 1}));
    ASSERT_TRUE(entry_5.any_haplotype(std::vector<unsigned int> {0, 2}));
    ASSERT_TRUE(entry_5.any_haplotype(std::vector<unsigned int> {0, 3, 4}));
    ASSERT_FALSE(entry_5.any_haplotype(std::vector<unsigned int> {0, 4}));
    ASSERT_FALSE(entry_5.any_haplotype(std::vector<unsigned int> {}));
}

TEST_F(VCF_Entry_Test, CanCountHaplotype){
    ASSERT_EQ(entry_5.count_haplotypes(std::vector<unsigned int> {0}), 0);
    ASSERT_EQ(entry_5.count_haplotypes(std::vector<unsigned int> {0, 1}), 1);
    ASSERT_EQ(entry_5.count_haplotypes(std::vector<unsigned int> {0, 2}), 2);
    ASSERT_EQ(entry_5.count_haplotypes(std::vector<unsigned int> {0, 3, 4}), 1);
    ASSERT_EQ(entry_5.count_haplotypes(std::vector<unsigned int> {1, 2, 3, 4}), 4);
    ASSERT_EQ(entry_5.count_haplotypes(std::vector<unsigned int> {0, 4}), 0);
    ASSERT_EQ(entry_5.count_haplotypes(std::vector<unsigned int> {}), 0);
}
