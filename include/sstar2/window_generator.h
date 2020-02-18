// class for yielding windows of a vcf file
// interacts with a VcfFile to read lines
// and PopulationData to determine which files are targets
// the window contains only targets and reference statistics are summarized

#pragma once
#include <iostream>
#include <set>
#include <deque>
#include <vector>
#include <algorithm>
#include "sstar2/vcf_file.h"
#include "sstar2/population_data.h"

struct WindowGT {
    int position_index;
    short genotype;  // encoded haplotypes, 3 -> homozygous alt
};

struct WindowBucket {
    long start, end;
    int site_count=0,  // number of snps after removing excluded, fixed, and homozygous ref
        reference_count=0;  // number in any ref
    std::vector<int> positions;
    std::vector<std::vector<WindowGT>> genotypes;
    void reset_bucket(int start, int end);
};

struct Window {
    std::string chromosome;
    // contains positions in (start, end]
    long start=0, end=0;
    std::deque<WindowBucket> buckets;
    void reset_window(std::string chrom, int length, int step);
    void step(int step);
};

std::ostream& operator<<(std::ostream &strm, const Window &window);
std::ostream& operator<<(std::ostream &strm, const WindowBucket &bucket);

class WindowGenerator{

    int length, step;
    bool terminated = false;

    std::istream *vcf = nullptr;
    std::string vcf_string;
    std::vector<int> targets;
    std::vector<std::string> target_names;
    std::vector<int> references;
    std::vector<int> excluded;
    void initialize_vcf();
    void initialize_window();
    bool next_line();

    public:
        VcfFile vcf_file;
        VcfEntry vcf_line{"", 0};
        PopulationData population;
        Window window;

        WindowGenerator(int window_length, int step_size);

        void initialize(
                std::istream &vcf_input,
                std::istream &pop_file,
                std::set<std::string> &target,
                std::set<std::string> &reference,
                std::set<std::string> &exclude);

        bool next_window();
};
