#include <iostream>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "sstar2/window.h"

using testing::ElementsAre;
using testing::Pair;

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

TEST(StepWindow, CanInitializeValidate){
    StepWindow window(10, 10);  // same size ok
    ASSERT_THROW(StepWindow win_bad(11, 10),
            std::invalid_argument);
}

TEST(StepWindow, CanStartWindow){
    StepWindow window(5, 10);
    ASSERT_EQ(window.start, 0);
    ASSERT_EQ(window.end, 0);

    // add some buckets
    window.initialize(1);
    VcfEntry line{"chrom", 1};
    window.start_window(line);
    ASSERT_STREQ(window.chromosome.c_str(), "chrom");
    ASSERT_EQ(window.start, 0);
    ASSERT_EQ(window.end, 10);
    ASSERT_EQ(window.callable_bases.totalLength(), 10);

    // add stuff
    std::vector<unsigned int> targets{0};

    line.reference = 'A';
    line.alternative = 'T';
    line.position = 1;
    line.genotypes[0] = 0;
    window.record(line, targets, 2);

    line.position = 2;
    line.genotypes[0] = 1;
    window.record(line, targets, 0);

    ASSERT_EQ(window.total_snps(), 2);
    ASSERT_EQ(window.reference_snps(), 1);
    ASSERT_EQ(window.individual_snps(0), 1);

    // then reset
    line.chromosome = "new chrom";
    window.start_window(line);
    ASSERT_STREQ(window.chromosome.c_str(), "new chrom");
    ASSERT_EQ(window.start, 0);
    ASSERT_EQ(window.end, 10);
    ASSERT_EQ(window.total_snps(), 0);
    ASSERT_EQ(window.reference_snps(), 0);
    ASSERT_EQ(window.individual_snps(0), 0);
}

TEST(StepWindow, CanStartWindowStep){
    StepWindow window(5, 10);

    // add some buckets
    window.initialize(1);
    VcfEntry line{"chrom", 1};
    window.start_window(line);

    // add stuff
    std::vector<unsigned int> targets{0};

    line.reference = 'A';
    line.alternative = 'T';
    line.position = 1;
    line.genotypes[0] = 0;

    window.record(line, targets, 2);
    line.position = 6;
    line.genotypes[0] = 1;
    window.record(line, targets, 0);

    // then step
    window.start_window(line);
    ASSERT_STREQ(window.chromosome.c_str(), "chrom");
    ASSERT_EQ(window.start, 5);
    ASSERT_EQ(window.end, 15);
    ASSERT_EQ(window.total_snps(), 1);
    ASSERT_EQ(window.reference_snps(), 0); // cleared in first bucket
    ASSERT_EQ(window.individual_snps(0), 1);
}

TEST(StepWindow, CanGetValues){
    StepWindow window(5, 10);

    // add some buckets
    window.initialize(2);
    VcfEntry line{"chrom", 2};
    window.start_window(line);  // set bucket start/end

    ASSERT_EQ(window.total_snps(), 0);
    ASSERT_EQ(window.reference_snps(), 0);
    ASSERT_EQ(window.individual_snps(0), 0);
    ASSERT_EQ(window.individual_snps(1), 0);

    // add stuff
    std::vector<unsigned int> targets{0, 1};

    line.reference = 'A';
    line.alternative = 'T';
    line.position = 1;
    line.genotypes[0] = 1;
    line.genotypes[1] = 0;
    window.record(line, targets, 0);

    line.position = 2;
    line.genotypes[0] = 2;
    line.genotypes[1] = 1;
    window.record(line, targets, 0);
    window.record(line, targets, 3);

    line.position = 7;
    line.genotypes[0] = 2;
    line.genotypes[1] = 1;
    window.record(line, targets, 0);
    window.record(line, targets, 3);

    ASSERT_EQ(window.total_snps(), 5);  // sum site count
    ASSERT_EQ(window.reference_snps(), 2);  // sum reference count
    ASSERT_EQ(window.individual_snps(0), 3);
    ASSERT_EQ(window.individual_snps(1), 2);
}

TEST(StepWindow, CanRecordValues){
    StepWindow window(5, 10);

    // add some buckets
    window.initialize(3);
    VcfEntry line{"chrom", 5};
    window.start_window(line);  // set bucket start/end

    ASSERT_EQ(window.total_snps(), 0);
    ASSERT_EQ(window.reference_snps(), 0);
    ASSERT_EQ(window.individual_snps(0), 0);
    ASSERT_EQ(window.individual_snps(1), 0);

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

    line.position = 6;
    line.genotypes[0] = 2;
    line.genotypes[1] = 1;
    line.genotypes[4] = 3;

    window.record(line, targets, 0);

    ASSERT_EQ(window.callable_bases.totalLength(), 10);
    ASSERT_EQ(window.total_snps(), 6);
    ASSERT_EQ(window.reference_snps(), 1);
    ASSERT_EQ(window.individual_snps(0), 3);
    ASSERT_EQ(window.individual_snps(1), 4);
    ASSERT_EQ(window.individual_snps(2), 3);

    line.position = 7;

    window.record(line, targets, 1);  // only reference is recorded

    ASSERT_EQ(window.callable_bases.totalLength(), 10);
    ASSERT_EQ(window.total_snps(), 7);
    ASSERT_EQ(window.reference_snps(), 2);
    ASSERT_EQ(window.individual_snps(0), 3);
    ASSERT_EQ(window.individual_snps(1), 4);
    ASSERT_EQ(window.individual_snps(2), 3);
}

TEST(RangedWindow, CanInitializeValidateString){
    RangedWindow window(10, 10, "1-2");  // same size ok
    RangedWindow window2(10, 10, "chrom:1-2");
    ASSERT_THROW(RangedWindow win_bad1(11, 10, "1-2"),
            std::invalid_argument);
    ASSERT_THROW(RangedWindow win_bad2(10, 10, "12"),
            std::invalid_argument);
    ASSERT_THROW(RangedWindow win_bad3(10, 10, "a-2"),
            std::invalid_argument);
    ASSERT_THROW(RangedWindow win_bad4(10, 10, "1-a"),
            std::invalid_argument);
    ASSERT_THROW(RangedWindow win_bad5(10, 10, "1-2:chrom"),
            std::invalid_argument);
    ASSERT_THROW(RangedWindow win_bad5(10, 10, ":1-2"),
            std::invalid_argument);
}

TEST(RangedWindow, CanInitializeValidateFile){
    std::istringstream input("chrom1 1 10\nchrom2 2 20");
    RangedWindow window(10, 10, input);  // same size ok
    ASSERT_THROW(RangedWindow win_bad1(11, 10, input),
            std::invalid_argument);
    input.str("chrom1 1 10\nchrom22 20");
    input.clear();
    testing::internal::CaptureStderr();
    RangedWindow window1(10, 10, input);  // should complain about input
    std::string errout = testing::internal::GetCapturedStderr();
    ASSERT_STREQ(errout.c_str(), "Unable to parse line from region file\n" 
            "chrom22 20\n");
}

TEST(RangedWindow, CanStartWindow){
    RangedWindow window(5, 10, "3-15");
    ASSERT_EQ(window.start, 0);
    ASSERT_EQ(window.end, 0);

    // add some buckets
    window.initialize(1);
    VcfEntry line{"chrom", 1};
    window.start_window(line);
    ASSERT_STREQ(window.chromosome.c_str(), "chrom");
    ASSERT_EQ(window.start, 3);
    ASSERT_EQ(window.end, 13);
    ASSERT_EQ(window.callable_bases.totalLength(), 10);

    // add stuff
    std::vector<unsigned int> targets{0};

    line.reference = 'A';
    line.alternative = 'T';
    line.position = 1;
    line.genotypes[0] = 0;
    window.record(line, targets, 2);

    line.position = 8;
    line.genotypes[0] = 1;
    window.record(line, targets, 0);

    ASSERT_EQ(window.total_snps(), 1);
    ASSERT_EQ(window.reference_snps(), 0);
    ASSERT_EQ(window.individual_snps(0), 1);

    // then reset
    line.chromosome = "new chrom";
    window.start_window(line);
    ASSERT_STREQ(window.chromosome.c_str(), "new chrom");
    ASSERT_EQ(window.start, 3);
    ASSERT_EQ(window.end, 13);
    ASSERT_EQ(window.total_snps(), 0);
    ASSERT_EQ(window.reference_snps(), 0);
    ASSERT_EQ(window.individual_snps(0), 0);

    RangedWindow window2(5, 10, "chrom:3-15");
    ASSERT_EQ(window2.start, 0);
    ASSERT_EQ(window2.end, 0);

    // add some buckets
    window2.initialize(1);
    line.chromosome = "chrom";
    window2.start_window(line);
    ASSERT_STREQ(window2.chromosome.c_str(), "chrom");
    ASSERT_EQ(window2.start, 3);
    ASSERT_EQ(window2.end, 13);
    ASSERT_EQ(window2.callable_bases.totalLength(), 10);

    // add stuff
    line.reference = 'A';
    line.alternative = 'T';
    line.position = 1;
    line.genotypes[0] = 0;
    window2.record(line, targets, 2);

    line.position = 8;
    line.genotypes[0] = 1;
    window2.record(line, targets, 0);

    ASSERT_EQ(window2.total_snps(), 1);
    ASSERT_EQ(window2.reference_snps(), 0);
    ASSERT_EQ(window2.individual_snps(0), 1);

    // then reset
    line.chromosome = "new chrom";
    window2.start_window(line);
    // will start by reseting start/end and counts, but new chrom isn't in regions
    ASSERT_STREQ(window2.chromosome.c_str(), "");
    ASSERT_EQ(window2.start, 3);
    ASSERT_EQ(window2.end, 13);
    ASSERT_EQ(window2.total_snps(), 0);
    ASSERT_EQ(window2.reference_snps(), 0);
    ASSERT_EQ(window2.individual_snps(0), 0);
}

