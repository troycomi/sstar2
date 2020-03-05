#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include "sstar2/window_generator.h"

class SStarCaller{
    long match_bonus;
    long mismatch_penalty;

    const char* emptyline = (
            "0\t0\t.\t" // sstar, num snps, snps
            "0\t0\t0\t0\t"  // hap1 and 2 start end
            "0\t0\t"  // sstar start and end
            "0\t0\t.\t");  // n snps hap 1 and 2, sstar haps

    public:
        SStarCaller() :
            match_bonus(5000), mismatch_penalty(-10000) {};
        SStarCaller(long bonus, long penalty) :
            match_bonus(bonus), mismatch_penalty(penalty) {};

        void write_header(std::ostream &output);
        // write the current window in generator
        void write_window(std::ostream &output,
                WindowGenerator &generator);
        // calculates sstar and updates the windowGT to include just snps
        long sstar(std::vector<WindowGT> &genotypes);
};
