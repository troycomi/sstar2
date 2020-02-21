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
    ASSERT_EQ(gen.vcf_line.genotypes.size(), 4);  // 4 indivs
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
    ASSERT_EQ(gen.vcf_line.genotypes.size(), 5);  // 4 indivs
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
    std::cout << gen.window;
    ASSERT_STREQ(gen.window.chromosome.c_str(), "1");
    ASSERT_EQ(gen.window.start, 0);
    ASSERT_EQ(gen.window.end, 10);
    ASSERT_EQ(gen.window.total_snps(), 2);
    ASSERT_EQ(gen.window.reference_snps(), 1);
    ASSERT_EQ(gen.window.individual_snps(0), 1);

    ASSERT_TRUE(gen.next_window());
    ASSERT_STREQ(gen.window.chromosome.c_str(), "1");
    ASSERT_EQ(gen.window.start, 5);
    ASSERT_EQ(gen.window.end, 15);
    ASSERT_EQ(gen.window.total_snps(), 2);
    ASSERT_EQ(gen.window.reference_snps(), 1);
    ASSERT_EQ(gen.window.individual_snps(0), 1);

    ASSERT_FALSE(gen.next_window());
}

TEST_F(Generator_Input, CanYieldWindow2){
    // repeated calls to yield window update chrom, start and end
    WindowGenerator gen(5, 2);
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
        ASSERT_STREQ(gen.window.chromosome.c_str(), "1");
        ASSERT_EQ(gen.window.start, 0 + i*2);
        ASSERT_EQ(gen.window.end, 5 + i*2);
        ASSERT_EQ(gen.window.total_snps(), total[i]);
        ASSERT_EQ(gen.window.reference_snps(), refs[i]);
        ASSERT_EQ(gen.window.individual_snps(0), ind0[i]);
        ASSERT_EQ(gen.window.individual_snps(1), ind1[i]);
        ASSERT_EQ(gen.window.individual_snps(2), ind2[i]);
    }

    ASSERT_FALSE(gen.next_window());
}

TEST_F(Generator_Input, CanYieldOddWindow){
    // repeated calls to yield window update chrom, start and end
    WindowGenerator gen(10, 3);
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
    ASSERT_EQ(gen.window.total_snps(), 2);
    ASSERT_EQ(gen.window.reference_snps(), 1);
    ASSERT_EQ(gen.window.individual_snps(0), 1);

    ASSERT_TRUE(gen.next_window());
    ASSERT_STREQ(gen.window.chromosome.c_str(), "1");
    ASSERT_EQ(gen.window.start, 3);
    ASSERT_EQ(gen.window.end, 13);
    ASSERT_EQ(gen.window.total_snps(), 2);
    ASSERT_EQ(gen.window.reference_snps(), 1);
    ASSERT_EQ(gen.window.individual_snps(0), 1);

    ASSERT_TRUE(gen.next_window());
    ASSERT_STREQ(gen.window.chromosome.c_str(), "1");
    ASSERT_EQ(gen.window.start, 6);
    ASSERT_EQ(gen.window.end, 16);
    ASSERT_EQ(gen.window.total_snps(), 1);
    ASSERT_EQ(gen.window.reference_snps(), 0);
    ASSERT_EQ(gen.window.individual_snps(0), 1);

    ASSERT_FALSE(gen.next_window());
}

TEST(WindowBucket, CanResetBucket){

    WindowBucket bucket;
    ASSERT_EQ(bucket.start, 0);
    ASSERT_EQ(bucket.end, 0);
    ASSERT_EQ(bucket.site_count, 0);
    ASSERT_EQ(bucket.reference_count, 0);
    ASSERT_EQ(bucket.genotypes.size(), 0);

    // empty bucket
    bucket.reset_bucket(0, 10);
    ASSERT_EQ(bucket.start, 0);
    ASSERT_EQ(bucket.end, 10);
    ASSERT_EQ(bucket.site_count, 0);
    ASSERT_EQ(bucket.reference_count, 0);
    ASSERT_EQ(bucket.genotypes.size(), 0);

    // modify some values
    bucket.site_count += 5;
    bucket.reference_count += 3;
    bucket.genotypes.resize(2);
    bucket.genotypes[0].push_back({1, 1});
    bucket.genotypes[0].push_back({2, 2});
    bucket.genotypes[1].push_back({2, 3});
    ASSERT_EQ(bucket.start, 0);
    ASSERT_EQ(bucket.end, 10);
    ASSERT_EQ(bucket.site_count, 5);
    ASSERT_EQ(bucket.reference_count, 3);
    ASSERT_EQ(bucket.genotypes.size(), 2);
    ASSERT_THAT(bucket.genotypes[0],
            ElementsAre(WindowGT{1, 1}, WindowGT{2, 2}));
    ASSERT_THAT(bucket.genotypes[1],
            ElementsAre(WindowGT{2, 3}));

    // empty bucket
    bucket.reset_bucket(5, 15);
    ASSERT_EQ(bucket.start, 5);
    ASSERT_EQ(bucket.end, 15);
    ASSERT_EQ(bucket.site_count, 0);
    ASSERT_EQ(bucket.reference_count, 0);
    ASSERT_EQ(bucket.genotypes.size(), 2);
    ASSERT_EQ(bucket.genotypes[0].size(), 0);
    ASSERT_EQ(bucket.genotypes[1].size(), 0);
}

TEST(WindowGenerator, CanResetWindow){
    Window window;
    ASSERT_EQ(window.start, 0);
    ASSERT_EQ(window.end, 0);
    ASSERT_EQ(window.buckets.size(), 0);

    // reset empty
    window.reset_window("chrom", 10, 5);  // length step
    ASSERT_STREQ(window.chromosome.c_str(), "chrom");
    ASSERT_EQ(window.start, 0);
    ASSERT_EQ(window.end, 10);
    ASSERT_EQ(window.callable_bases, 10);
    ASSERT_EQ(window.buckets.size(), 0);

    // add some buckets
    window.buckets.resize(2);
    window.buckets[0].genotypes.resize(1);
    window.buckets[1].genotypes.resize(1);
    window.buckets[0].site_count += 5;
    window.buckets[0].reference_count += 3;
    window.buckets[0].genotypes[0].push_back({0, 1});
    window.buckets[0].genotypes[0].push_back({1, 2});
    ASSERT_EQ(window.buckets[0].start, 0);  //defaults
    ASSERT_EQ(window.buckets[0].end, 0);
    ASSERT_EQ(window.buckets[0].site_count, 5);  //set above
    ASSERT_EQ(window.buckets[0].reference_count, 3);
    ASSERT_EQ(window.buckets[0].genotypes.size(), 1);  //one indiv, 2 GTs
    ASSERT_EQ(window.buckets[0].genotypes[0].size(), 2);
    ASSERT_THAT(window.buckets[0].genotypes[0],
            ElementsAre(WindowGT{0, 1}, WindowGT{1, 2}));

    // then reset
    window.reset_window("new chrom", 10, 5);
    ASSERT_STREQ(window.chromosome.c_str(), "new chrom");
    ASSERT_EQ(window.start, 0);
    ASSERT_EQ(window.end, 10);
    ASSERT_EQ(window.buckets.size(), 2);
    ASSERT_EQ(window.buckets[0].start, 0);  // set by window
    ASSERT_EQ(window.buckets[0].end, 5);
    ASSERT_EQ(window.buckets[0].site_count, 0);  // clear these
    ASSERT_EQ(window.buckets[0].reference_count, 0);
    ASSERT_EQ(window.buckets[0].genotypes.size(), 1);
    ASSERT_EQ(window.buckets[0].genotypes[0].size(), 0);
    ASSERT_EQ(window.buckets[1].start, 5);
    ASSERT_EQ(window.buckets[1].end, 10);
    ASSERT_EQ(window.buckets[1].genotypes.size(), 1);
}

TEST(WindowGenerator, CanStepWindow){
    Window window;

    // add some buckets
    window.buckets.resize(2);
    window.reset_window("chrom", 10, 5);  // set bucket start/end
    window.buckets[0].genotypes.resize(1);
    window.buckets[1].genotypes.resize(1);
    window.buckets[0].site_count += 5;
    window.buckets[0].reference_count += 3;
    window.buckets[0].genotypes[0].push_back({1, 1});
    window.buckets[0].genotypes[0].push_back({2, 2});

    window.buckets[1].site_count += 2;
    window.buckets[1].reference_count += 4;
    window.buckets[1].genotypes[0].push_back({2, 3});
    window.buckets[1].genotypes[0].push_back({3, 2});

    // then step
    window.step(5);
    ASSERT_STREQ(window.chromosome.c_str(), "chrom");
    ASSERT_EQ(window.start, 5);
    ASSERT_EQ(window.end, 15);
    ASSERT_EQ(window.buckets.size(), 2);
    ASSERT_EQ(window.buckets[0].start, 5);  // set by window
    ASSERT_EQ(window.buckets[0].end, 10);
    ASSERT_EQ(window.buckets[0].site_count, 2);  // from previous bucket '1'
    ASSERT_EQ(window.buckets[0].reference_count, 4);
    ASSERT_EQ(window.buckets[0].genotypes.size(), 1);
    ASSERT_THAT(window.buckets[0].genotypes[0],
            ElementsAre(WindowGT{2, 3}, WindowGT{3, 2}));
    // these are cleared
    ASSERT_EQ(window.buckets[1].start, 10);
    ASSERT_EQ(window.buckets[1].end, 15);
    ASSERT_EQ(window.buckets[1].genotypes.size(), 1);
}

TEST(WindowGenerator, CanGetValues){
    Window window;

    // add some buckets
    window.buckets.resize(3);
    window.reset_window("chrom", 10, 5);  // set bucket start/end

    window.buckets[0].genotypes.resize(2);
    window.buckets[1].genotypes.resize(2);
    window.buckets[2].genotypes.resize(2);

    ASSERT_EQ(window.total_snps(), 0);
    ASSERT_EQ(window.reference_snps(), 0);
    ASSERT_EQ(window.individual_snps(0), 0);
    ASSERT_EQ(window.individual_snps(1), 0);

    window.buckets[0].site_count += 5;
    window.buckets[0].reference_count += 3;
    window.buckets[0].genotypes[0].push_back({0, 1});
    window.buckets[0].genotypes[0].push_back({1, 2});

    window.buckets[1].site_count += 2;
    window.buckets[1].reference_count += 4;
    window.buckets[1].genotypes[0].push_back({0, 3});
    window.buckets[1].genotypes[1].push_back({1, 2});

    window.buckets[2].site_count += 3;
    window.buckets[2].reference_count += 0;
    window.buckets[2].genotypes[0].push_back({0, 3});
    window.buckets[2].genotypes[1].push_back({0, 3});
    window.buckets[2].genotypes[1].push_back({1, 2});

    ASSERT_EQ(window.total_snps(), 10);  // sum site count
    ASSERT_EQ(window.reference_snps(), 7);  // sum reference count
    ASSERT_EQ(window.individual_snps(0), 4);
    ASSERT_EQ(window.individual_snps(1), 3);
}

TEST(WindowGenerator, CanRecordValues){
    Window window;

    // add some buckets
    window.buckets.resize(2);
    window.reset_window("chrom", 10, 5);  // set bucket start/end

    window.buckets[0].genotypes.resize(3);
    window.buckets[1].genotypes.resize(3);

    ASSERT_EQ(window.total_snps(), 0);
    ASSERT_EQ(window.reference_snps(), 0);
    ASSERT_EQ(window.individual_snps(0), 0);
    ASSERT_EQ(window.individual_snps(1), 0);

    VcfEntry line{"chrom", 5};
    std::vector<unsigned int> targets{0, 1, 4};

    line.reference = 'A';
    line.alternative = 'T';
    // genotypes 2-3 are ref and not directly checked

    line.position = 1;
    line.genotypes[0] = 0;
    line.genotypes[1] = 0;
    line.genotypes[4] = 0;

    window.record(line, targets, 2);

    ASSERT_EQ(window.total_snps(), 1);
    ASSERT_EQ(window.reference_snps(), 1);
    ASSERT_EQ(window.individual_snps(0), 0);
    ASSERT_EQ(window.individual_snps(1), 0);
    ASSERT_EQ(window.individual_snps(2), 0);

    line.position = 2;
    window.record(line, targets, 0);  // just update total

    ASSERT_EQ(window.total_snps(), 2);
    ASSERT_EQ(window.reference_snps(), 1);
    ASSERT_EQ(window.individual_snps(0), 0);
    ASSERT_EQ(window.individual_snps(1), 0);
    ASSERT_EQ(window.individual_snps(2), 0);

    line.position = 3;
    line.genotypes[0] = 2;
    line.genotypes[1] = 3;
    line.genotypes[4] = 1;

    window.record(line, targets, 0);

    ASSERT_EQ(window.total_snps(), 3);
    ASSERT_EQ(window.reference_snps(), 1);
    ASSERT_EQ(window.individual_snps(0), 1);
    ASSERT_EQ(window.individual_snps(1), 1);
    ASSERT_EQ(window.individual_snps(2), 1);

    ASSERT_THAT(window.buckets[0].genotypes[0], 
            ElementsAre(WindowGT{3, 2}));
    ASSERT_THAT(window.buckets[0].genotypes[1], 
            ElementsAre(WindowGT{3, 3}));
    ASSERT_THAT(window.buckets[0].genotypes[2], 
            ElementsAre(WindowGT{3, 1}));

    line.position = 4;
    line.genotypes[0] = 2;
    line.genotypes[1] = 2;
    line.genotypes[4] = 0;

    window.record(line, targets, 0);

    ASSERT_EQ(window.total_snps(), 4);
    ASSERT_EQ(window.reference_snps(), 1);
    ASSERT_EQ(window.individual_snps(0), 2);
    ASSERT_EQ(window.individual_snps(1), 2);
    ASSERT_EQ(window.individual_snps(2), 1);

    ASSERT_THAT(window.buckets[0].genotypes[0], 
            ElementsAre(WindowGT{3, 2}, WindowGT{4, 2}));
    ASSERT_THAT(window.buckets[0].genotypes[1], 
            ElementsAre(WindowGT{3, 3}, WindowGT{4, 2}));
    ASSERT_THAT(window.buckets[0].genotypes[2], 
            ElementsAre(WindowGT{3, 1}));

    line.position = 5;
    line.genotypes[0] = 0;
    line.genotypes[1] = 1;
    line.genotypes[4] = 1;

    window.record(line, targets, 0);

    ASSERT_EQ(window.total_snps(), 5);
    ASSERT_EQ(window.reference_snps(), 1);
    ASSERT_EQ(window.individual_snps(0), 2);
    ASSERT_EQ(window.individual_snps(1), 3);
    ASSERT_EQ(window.individual_snps(2), 2);

    ASSERT_THAT(window.buckets[0].genotypes[0], 
            ElementsAre(WindowGT{3, 2}, WindowGT{4, 2}));
    ASSERT_THAT(window.buckets[0].genotypes[1], 
            ElementsAre(WindowGT{3, 3}, WindowGT{4, 2}, WindowGT{5, 1}));
    ASSERT_THAT(window.buckets[0].genotypes[2], 
            ElementsAre(WindowGT{3, 1}, WindowGT{5, 1}));

    line.position = 6;
    line.genotypes[0] = 2;
    line.genotypes[1] = 1;
    line.genotypes[4] = 3;

    window.record(line, targets, 0);

    ASSERT_EQ(window.callable_bases, 10);
    ASSERT_EQ(window.total_snps(), 6);
    ASSERT_EQ(window.reference_snps(), 1);
    ASSERT_EQ(window.individual_snps(0), 3);
    ASSERT_EQ(window.individual_snps(1), 4);
    ASSERT_EQ(window.individual_snps(2), 3);

    ASSERT_THAT(window.buckets[0].genotypes[0], 
            ElementsAre(WindowGT{3, 2}, WindowGT{4, 2}));
    ASSERT_THAT(window.buckets[0].genotypes[1], 
            ElementsAre(WindowGT{3, 3}, WindowGT{4, 2}, WindowGT{5, 1}));
    ASSERT_THAT(window.buckets[0].genotypes[2], 
            ElementsAre(WindowGT{3, 1}, WindowGT{5, 1}));

    ASSERT_THAT(window.buckets[1].genotypes[0], 
            ElementsAre(WindowGT{6, 2}));
    ASSERT_THAT(window.buckets[1].genotypes[1], 
            ElementsAre(WindowGT{6, 1}));
    ASSERT_THAT(window.buckets[1].genotypes[2], 
            ElementsAre(WindowGT{6, 3}));

    line.position = 7;

    window.record(line, targets, 1);  // only reference is recorded

    ASSERT_EQ(window.callable_bases, 10);
    ASSERT_EQ(window.total_snps(), 7);
    ASSERT_EQ(window.reference_snps(), 2);
    ASSERT_EQ(window.individual_snps(0), 3);
    ASSERT_EQ(window.individual_snps(1), 4);
    ASSERT_EQ(window.individual_snps(2), 3);

    ASSERT_THAT(window.buckets[0].genotypes[0], 
            ElementsAre(WindowGT{3, 2}, WindowGT{4, 2}));
    ASSERT_THAT(window.buckets[0].genotypes[1], 
            ElementsAre(WindowGT{3, 3}, WindowGT{4, 2}, WindowGT{5, 1}));
    ASSERT_THAT(window.buckets[0].genotypes[2], 
            ElementsAre(WindowGT{3, 1}, WindowGT{5, 1}));

    ASSERT_THAT(window.buckets[1].genotypes[0], 
            ElementsAre(WindowGT{6, 2}));
    ASSERT_THAT(window.buckets[1].genotypes[1], 
            ElementsAre(WindowGT{6, 1}));
    ASSERT_THAT(window.buckets[1].genotypes[2], 
            ElementsAre(WindowGT{6, 3}));
}
