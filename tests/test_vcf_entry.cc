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
            entry_5.genotypes[1] = 2;
            entry_5.genotypes[2] = 3;
            entry_5.genotypes[3] = 1;
        }

        VcfEntry entry_5{"1", 5};
};

TEST_F(VCF_Entry_Test, CanGetGenotype){
    ASSERT_EQ(entry_5.genotypes[0], 0);
    ASSERT_EQ(entry_5.genotypes[1], 2);
    ASSERT_EQ(entry_5.genotypes[2], 3);
    ASSERT_EQ(entry_5.genotypes[3], 1);
    ASSERT_EQ(entry_5.genotypes[4], 0);
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
