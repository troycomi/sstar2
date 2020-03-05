#include <iostream>
#include <sstream>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <vector>

#include "sstar2/validator.h"

TEST(BedValidator, CanTestBed){
    std::istringstream infile(
            "1\t0\t10\n"
            "1\t20\t30\n"
            "2\t20\t30\n"
            "4\t20\t30\n"
        );
    BedFile b(&infile);

    ASSERT_FALSE(b.inBed("1", 0));
    ASSERT_TRUE(b.inBed("1", 1));
    ASSERT_TRUE(b.inBed("1", 10));
    ASSERT_FALSE(b.inBed("1", 11));

    ASSERT_FALSE(b.inBed("1", 20));
    ASSERT_TRUE(b.inBed("1", 21));
    ASSERT_TRUE(b.inBed("1", 30));
    ASSERT_FALSE(b.inBed("1", 31));

    ASSERT_FALSE(b.inBed("2", 20));
    ASSERT_TRUE(b.inBed("2", 21));
    ASSERT_TRUE(b.inBed("2", 30));
    ASSERT_FALSE(b.inBed("2", 31));

    ASSERT_FALSE(b.inBed("3", 20));
    ASSERT_FALSE(b.inBed("3", 21));

    ASSERT_FALSE(b.inBed("4", 20));
    ASSERT_TRUE(b.inBed("4", 21));
    ASSERT_TRUE(b.inBed("4", 30));
    ASSERT_FALSE(b.inBed("4", 31));
    // bed is empty
    ASSERT_FALSE(b.inBed("asdf", 5));
    ASSERT_FALSE(b.inBed("asdf", 6));
    ASSERT_FALSE(b.inBed("asdf", 7));
    ASSERT_TRUE(b.inBed("4", 30));
}

TEST(BedValidator, PositiveBedValidator){
    std::istringstream infile(
            "1\t0\t10\n"
            "1\t20\t30\n"
            "2\t20\t30\n"
            "4\t20\t30\n"
        );
    PositiveBedValidator validator(&infile);
    VcfEntry entry("", 0);

    entry.chromosome = "1";
    entry.position = 0;
    ASSERT_FALSE(validator.isValid(entry));
    entry.position = 1;
    ASSERT_TRUE(validator.isValid(entry));
    entry.position = 10;
    ASSERT_TRUE(validator.isValid(entry));
    entry.position = 11;
    ASSERT_FALSE(validator.isValid(entry));

    entry.chromosome = "4";
    entry.position = 20;
    ASSERT_FALSE(validator.isValid(entry));
    entry.position = 21;
    ASSERT_TRUE(validator.isValid(entry));
    entry.position = 30;
    ASSERT_TRUE(validator.isValid(entry));
    entry.position = 31;
    ASSERT_FALSE(validator.isValid(entry));

    // bed is empty
    entry.chromosome = "asdf";
    ASSERT_FALSE(validator.isValid(entry));
    ASSERT_FALSE(validator.isValid(entry));
    ASSERT_FALSE(validator.isValid(entry));
}

TEST(BedValidator, PositiveBedCallable){
    std::istringstream infile(
            "1\t0\t10\n"
            "1\t20\t30\n"
            "2\t20\t30\n"
            "4\t20\t30\n"
        );
    PositiveBedValidator validator(&infile);
    BaseRegions region;
    VcfEntry entry("", 0);

    region.set("1", 0, 100);
    validator.updateCallable(region);  // this reads to end
    ASSERT_EQ(region.totalLength(), 20);

    entry.chromosome = "1";
    entry.position = 5;
    ASSERT_TRUE(validator.isValid(entry));

    entry.chromosome = "4";
    entry.position = 20;
    ASSERT_FALSE(validator.isValid(entry));
    entry.position = 21;
    ASSERT_TRUE(validator.isValid(entry));
    entry.position = 30;
    ASSERT_TRUE(validator.isValid(entry));
    entry.position = 31;
    ASSERT_FALSE(validator.isValid(entry));

    region.set("4", 0, 100);
    validator.updateCallable(region);
    ASSERT_EQ(region.totalLength(), 10);

    // bed is empty
    entry.chromosome = "asdf";
    ASSERT_FALSE(validator.isValid(entry));
    ASSERT_FALSE(validator.isValid(entry));
    ASSERT_FALSE(validator.isValid(entry));

    region.set("1", 0, 100);
    validator.updateCallable(region);
    ASSERT_EQ(region.totalLength(), 20);
    region.set("2", 0, 100);
    validator.updateCallable(region);
    ASSERT_EQ(region.totalLength(), 10);
    region.set("4", 0, 100);
    validator.updateCallable(region);
    ASSERT_EQ(region.totalLength(), 10);
}

TEST(BedValidator, NegativeBedValidator){
    std::istringstream infile(
            "1\t0\t10\n"
            "1\t20\t30\n"
            "2\t20\t30\n"
            "4\t20\t30\n"
        );
    NegativeBedValidator validator(&infile);
    VcfEntry entry("", 0);

    entry.chromosome = "1";
    entry.position = 0;
    ASSERT_TRUE(validator.isValid(entry));
    entry.position = 1;
    ASSERT_FALSE(validator.isValid(entry));
    entry.position = 10;
    ASSERT_FALSE(validator.isValid(entry));
    entry.position = 11;
    ASSERT_TRUE(validator.isValid(entry));

    entry.chromosome = "4";
    entry.position = 20;
    ASSERT_TRUE(validator.isValid(entry));
    entry.position = 21;
    ASSERT_FALSE(validator.isValid(entry));
    entry.position = 30;
    ASSERT_FALSE(validator.isValid(entry));
    entry.position = 31;
    ASSERT_TRUE(validator.isValid(entry));

    // bed is empty
    entry.chromosome = "asdf";
    ASSERT_TRUE(validator.isValid(entry));
    ASSERT_TRUE(validator.isValid(entry));
    ASSERT_TRUE(validator.isValid(entry));
}

TEST(BedValidator, NegativeBedCallable){
    std::istringstream infile(
            "1\t0\t10\n"
            "1\t20\t30\n"
            "2\t20\t30\n"
            "4\t20\t30\n"
        );
    NegativeBedValidator validator(&infile);
    BaseRegions region;
    region.set("1", 0, 100);
    VcfEntry entry("", 0);

    validator.updateCallable(region);
    ASSERT_EQ(region.totalLength(), 80);

    entry.chromosome = "1";
    entry.position = 21;
    ASSERT_FALSE(validator.isValid(entry));
    region.set("1", 0, 100);
    validator.updateCallable(region);
    ASSERT_EQ(region.totalLength(), 80);

    entry.chromosome = "4";
    entry.position = 20;
    ASSERT_TRUE(validator.isValid(entry));

    region.set("1", 0, 100);
    validator.updateCallable(region);
    ASSERT_EQ(region.totalLength(), 80);

    // bed is empty
    region.set("4", 0, 100);
    validator.updateCallable(region);
    ASSERT_EQ(region.totalLength(), 90);
    entry.chromosome = "asdf";
    ASSERT_TRUE(validator.isValid(entry));
    ASSERT_TRUE(validator.isValid(entry));
    ASSERT_TRUE(validator.isValid(entry));
    region.set("1", 0, 100);
    validator.updateCallable(region);
    ASSERT_EQ(region.totalLength(), 80);
    region.set("2", 0, 100);
    validator.updateCallable(region);
    ASSERT_EQ(region.totalLength(), 90);
    region.set("4", 0, 100);
    validator.updateCallable(region);
    ASSERT_EQ(region.totalLength(), 90);
}

TEST(FixValidator, CanValidateFixation){
    VcfEntry entry("", 5);
    FixationValidator valid(
            std::vector<unsigned int>{0, 2},
            std::vector<unsigned int>{3});
    ASSERT_FALSE(valid.isValid(entry));
    // positions 1 and 4 don't matter
    for(uint8_t i = 0; i < 4; ++i){
        entry.genotypes[1] = i;
        for(uint8_t j = 0; j < 4; ++j){
            entry.genotypes[4] = j;
            ASSERT_FALSE(valid.isValid(entry));
        }
    }
    entry.genotypes[0] = 1;
    // positions 1 and 4 don't matter
    for(uint8_t i = 0; i < 4; ++i){
        entry.genotypes[1] = i;
        for(uint8_t j = 0; j < 4; ++j){
            entry.genotypes[4] = j;
            ASSERT_TRUE(valid.isValid(entry));
        }
    }
    entry.genotypes[0] = 3;
    ASSERT_TRUE(valid.isValid(entry));
    entry.genotypes[2] = 3;
    ASSERT_TRUE(valid.isValid(entry));
    entry.genotypes[3] = 2;
    ASSERT_TRUE(valid.isValid(entry));
    entry.genotypes[3] = 3;
    // fixed!
    ASSERT_FALSE(valid.isValid(entry));
}

TEST(BaseRegions, CanSet){
    BaseRegions region;
    std::ostringstream output;
    output << region;
    ASSERT_STREQ(output.str().c_str(), "No region\n");
    ASSERT_EQ(region.totalLength(), 0);
    output.str("");

    std::string chrom = "chr1";
    region.set(chrom, 4, 7);
    output << region;
    ASSERT_STREQ(output.str().c_str(), "chr1:4,7,\n");
    output.str("");

    chrom = "chr2";
    region.set(chrom, 6, 9);
    output << region;
    ASSERT_STREQ(output.str().c_str(), "chr2:6,9,\n");
    ASSERT_EQ(region.totalLength(), 3);
}

TEST(BaseRegions, CanAdd){
    BaseRegions region;
    std::ostringstream output;

    region.add("chr1", 4, 7);
    output << region;
    ASSERT_STREQ(output.str().c_str(), "chr1:4,7,\n");
    ASSERT_EQ(region.totalLength(), 3);

    output.str("");
    region.add("chr1", 16, 19);
    output << region;
    ASSERT_STREQ(output.str().c_str(), "chr1:4,7,16,19,\n");
    ASSERT_EQ(region.totalLength(), 6);

    output.str("");
    region.add("chr2", 16, 19);
    output << region;
    ASSERT_STREQ(output.str().c_str(), "chr1:4,7,16,19,\nchr2:16,19,\n");
    ASSERT_EQ(region.totalLength(), 9);

    output.str("");
    region.add("chr2", 19, 21);
    output << region;
    ASSERT_STREQ(output.str().c_str(), "chr1:4,7,16,19,\nchr2:16,19,19,21,\n");
    ASSERT_EQ(region.totalLength(), 11);
}

void setTargetSingle(BaseRegions &target){
    std::string chrom = "1";
    target.set(chrom, 50, 100);
}

TEST(BaseRegions, CanIntersect){
    // empty cases
    BaseRegions target;
    BaseRegions other;
    ASSERT_EQ(target.totalLength(), 0);

    target.intersect(other);
    ASSERT_EQ(target.totalLength(), 0);

    std::string chrom = "1";
    target.add(chrom, 50, 100);
    target.intersect(other);
    ASSERT_EQ(target.totalLength(), 0);

    setTargetSingle(target);
    std::string chrom2 = "2";
    other.add(chrom2, 50, 100);
    target.intersect(other);
    ASSERT_EQ(target.totalLength(), 0);

    // no overlap
    setTargetSingle(target);
    other.set(chrom, 25, 35);
    target.intersect(other);
    ASSERT_EQ(target.totalLength(), 0);

    // partial overlap
    setTargetSingle(target);
    other.add(chrom, 45, 55);
    target.intersect(other);
    ASSERT_EQ(target.totalLength(), 5);

    // full overlap
    setTargetSingle(target);
    other.add(chrom, 65, 75);
    target.intersect(other);
    ASSERT_EQ(target.totalLength(), 15);

    // partial overlap end
    setTargetSingle(target);
    other.add(chrom, 95, 105);
    target.intersect(other);
    ASSERT_EQ(target.totalLength(), 20);

    // no overlap end
    setTargetSingle(target);
    other.add(chrom, 115, 125);
    target.intersect(other);
    ASSERT_EQ(target.totalLength(), 20);

    std::ostringstream output;
    output << target;
    ASSERT_STREQ(output.str().c_str(), "1:50,55,65,75,95,100,\n");
}

TEST(BaseRegions, CanIntersectSingletons){
    BaseRegions target;
    BaseRegions other;

    std::string chrom = "1";
    std::vector<unsigned long> starts{25, 25, 25, 25, 25,  // start before
        50, 50, 50, 50,  // start on other
        75, 75, 75,  // start between
        100, 100, 105  // start on or after end
    };
    std::vector<unsigned long> ends{35, 50, 55, 100, 105,
        50, 55, 100, 105,
        85, 100, 105,
        100, 105, 110
    };
    std::vector<unsigned long> lengths{0, 0, 5, 50, 50,
        0, 5, 50, 50,
        10, 25, 25,
        0, 0, 0
    };

    for(size_t i = 0; i < starts.size(); ++i){
        setTargetSingle(target);
        other.set(chrom, starts[i], ends[i]);
        target.intersect(other);
        ASSERT_EQ(target.totalLength(), lengths[i]);
        // and reverse
        setTargetSingle(target);
        other.intersect(target);
        ASSERT_EQ(other.totalLength(), lengths[i]);
    }
}

TEST(BaseRegions, CanSubtract){
    // empty cases
    BaseRegions target;
    BaseRegions other;
    ASSERT_EQ(target.totalLength(), 0);

    target.subtract(other);
    ASSERT_EQ(target.totalLength(), 0);

    std::string chrom = "1";
    target.add(chrom, 50, 100);
    target.subtract(other);
    ASSERT_EQ(target.totalLength(), 50);

    setTargetSingle(target);
    std::string chrom2 = "2";
    other.add(chrom2, 50, 100);
    target.subtract(other);
    ASSERT_EQ(target.totalLength(), 50);

    // no overlap
    setTargetSingle(target);
    other.set(chrom, 25, 35);
    target.subtract(other);
    ASSERT_EQ(target.totalLength(), 50);

    // partial overlap
    setTargetSingle(target);
    other.add(chrom, 45, 55);
    target.subtract(other);
    ASSERT_EQ(target.totalLength(), 45);

    // full overlap
    setTargetSingle(target);
    other.add(chrom, 65, 75);
    target.subtract(other);
    ASSERT_EQ(target.totalLength(), 35);

    // partial overlap end
    setTargetSingle(target);
    other.add(chrom, 95, 105);
    target.subtract(other);
    ASSERT_EQ(target.totalLength(), 30);

    // no overlap end
    setTargetSingle(target);
    other.add(chrom, 115, 125);
    target.subtract(other);
    ASSERT_EQ(target.totalLength(), 30);

    std::ostringstream output;
    output << target;
    ASSERT_STREQ(output.str().c_str(), "1:55,65,75,95,\n");
}

TEST(BaseRegions, CanSubtractSingletons){
    BaseRegions target;
    BaseRegions other;

    std::string chrom = "1";
    std::vector<unsigned long> starts{25, 25, 25, 25, 25,  // start before
        50, 50, 50, 50,  // start on other
        75, 75, 75,  // start between
        100, 100, 105  // start on or after end
    };
    std::vector<unsigned long> ends{35, 50, 55, 100, 105,
        50, 55, 100, 105,
        85, 100, 105,
        100, 105, 110
    };
    std::vector<unsigned long> lengths{50, 50, 45, 0, 0,
        50, 45, 0, 0,
        40, 25, 25,
        50, 50, 50
    };
    std::vector<unsigned long> lengths2{10, 25, 25, 25, 30,
        0, 0, 0, 5,
        0, 0, 5,
        0, 5, 5
    };

    for(size_t i = 0; i < starts.size(); ++i){
        setTargetSingle(target);
        other.set(chrom, starts[i], ends[i]);
        target.subtract(other);
        ASSERT_EQ(target.totalLength(), lengths[i]);
        // and reverse
        setTargetSingle(target);
        other.subtract(target);
        ASSERT_EQ(other.totalLength(), lengths2[i]);
    }
}

void setTargetTriplet(BaseRegions &target){
    std::string chrom = "1";
    target.set(chrom, 50, 100);
    target.add(chrom, 150, 200);
    target.add(chrom, 250, 300);
}

TEST(BaseRegions, CanIntersectTriplets){
    BaseRegions target;
    BaseRegions other;
    std::string chrom = "1";

    std::vector<unsigned long> starts{
        25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325
    };
    // match region endpoints, add 5 otherwise
    std::vector<unsigned long> ends{
        30, 50, 80, 100, 130, 150, 180, 200, 230, 250, 280, 300, 330
    };

    for(auto & start : starts){
        for(auto & end : ends){
            if (end < start)
                continue;
            unsigned long length = 150;
            if(start < 50)
                ;
            else if(start >= 50 && start <= 100)
                length -= start - 50;
            else if(start >= 100 && start <= 150)
                length -= 50;
            else if(start >= 150 && start <= 200)
                length -= 50 + (start - 150);
            else if(start >= 200 && start <= 250)
                length -= 100;
            else if(start >= 250 && start <= 300)
                length -= 100 + (start - 250);
            else if(start > 300)
                length -= 150;

            if(end < 50)
                length -= 150;
            else if(end >= 50 && end <= 100)
                length -= 100 + (100 - end);
            else if(end >= 100 && end <= 150)
                length -= 100;
            else if(end >= 150 && end <= 200)
                length -= 50 + (200 - end);
            else if(end >= 200 && end <= 250)
                length -= 50;
            else if(end >= 250 && end <= 300)
                length -= 300 - end;

            setTargetTriplet(target);
            other.set(chrom, start, end);
            target.intersect(other);
            ASSERT_EQ(target.totalLength(), length);
            // and reverse
            setTargetTriplet(target);
            other.intersect(target);
            ASSERT_EQ(other.totalLength(), length);
        }
    }
}

TEST(BaseRegions, CanSubtractTriplets){
    BaseRegions target;
    BaseRegions other;
    std::string chrom = "1";

    std::vector<unsigned long> starts{
        25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325
    };
    // match region endpoints, add 5 otherwise
    std::vector<unsigned long> ends{
        30, 50, 80, 100, 130, 150, 180, 200, 230, 250, 280, 300, 330
    };
    // enumerate because the code was not coming out easily
    std::vector<unsigned long> lengths{
        150, 150, 120, 100, 100, 100, 70, 50, 50, 50, 20, 0, 0,
        1, 150, 120, 100, 100, 100, 70, 50, 50, 50, 20, 0, 0,
        1, 1, 145, 125, 125, 125, 95, 75, 75, 75, 45, 25, 25,
        1, 1, 1, 150, 150, 150, 120, 100, 100, 100, 70, 50, 50,
        1, 1, 1, 1, 150, 150, 120, 100, 100, 100, 70, 50, 50,

        1, 1, 1, 1, 1, 150, 120, 100, 100, 100, 70, 50, 50,
        1, 1, 1, 1, 1, 1, 145, 125, 125, 125, 95, 75, 75,
        1, 1, 1, 1, 1, 1, 1, 150, 150, 150, 120, 100, 100,
        1, 1, 1, 1, 1, 1, 1, 1, 150, 150, 120, 100, 100,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 150, 120, 100, 100,

        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 145, 125, 125,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 150, 150,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 150,
    };
    std::vector<unsigned long> lengths2{
        5, 25, 25, 25, 55, 75, 75, 75, 105, 125, 125, 125, 155,
        1, 0, 0, 0, 30, 50, 50, 50, 80, 100, 100, 100, 130,
        1, 1, 0, 0, 30, 50, 50, 50, 80, 100, 100, 100, 130,
        1, 1, 1, 0, 30, 50, 50, 50, 80, 100, 100, 100, 130,
        1, 1, 1, 1, 5, 25, 25, 25, 55, 75, 75, 75, 105,

        1, 1, 1, 1, 1, 0, 0, 0, 30, 50, 50, 50, 80,
        1, 1, 1, 1, 1, 1, 0, 0, 30, 50, 50, 50, 80,
        1, 1, 1, 1, 1, 1, 1, 0, 30, 50, 50, 50, 80,
        1, 1, 1, 1, 1, 1, 1, 1, 5, 25, 25, 25, 55,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 30,

        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 30,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 30,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5,
    };

    lengths.insert(lengths.end(), 13*13, 500);
    lengths2.insert(lengths2.end(), 13*13, 500);

    for(auto i = 0; i < starts.size(); ++i){
        for(auto j = i; j < ends.size(); ++j){
            setTargetTriplet(target);
            other.set(chrom, starts[i], ends[j]);
            target.subtract(other);
            ASSERT_EQ(target.totalLength(), lengths[i*13+j]);
            // and reverse
            setTargetTriplet(target);
            other.subtract(target);
            ASSERT_EQ(other.totalLength(), lengths2[i*13+j]);
        }
    }
}
