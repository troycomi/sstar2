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
    unsigned long position;
    short unsigned int genotype;  // encoded haplotypes, 3 -> homozygous alt
    bool operator==(const WindowGT& rhs) const;

    WindowGT(unsigned long pos, short unsigned int gt):
        position(pos), genotype(gt){};
};

struct WindowBucket {
    unsigned long start=0, end=0;
    unsigned int site_count=0,  // number of snps after removing excluded, fixed, and homozygous ref
        reference_count=0;  // number in any ref
    std::vector<std::vector<WindowGT>> genotypes;
    void reset_bucket(unsigned int start, unsigned int end);
};

struct Window {
    std::string chromosome;
    // contains positions in (start, end]
    unsigned long start=0, end=0, callable_bases=0;
    std::deque<WindowBucket> buckets;
    void reset_window(std::string chrom, unsigned int length, unsigned int step);
    void step(unsigned int step);
    void record(const VcfEntry &entry, const std::vector<unsigned int> &targets,
            unsigned int reference_haplotypes);
    unsigned int total_snps() const;
    unsigned int reference_snps() const;
    unsigned int individual_snps(unsigned int individual) const;
};

std::ostream& operator<<(std::ostream &strm, const Window &window);
std::ostream& operator<<(std::ostream &strm, const WindowBucket &bucket);

class WindowGenerator{

    unsigned int length, step;
    bool terminated = false;

    std::istream *vcf = nullptr;
    std::string vcf_string;
    std::vector<unsigned int> references;
    std::vector<unsigned int> excluded;
    void initialize_vcf();
    void initialize_window();
    bool next_line();

    public:
        VcfFile vcf_file;
        VcfEntry vcf_line{"", 0};
        PopulationData population;
        Window window;
        std::vector<unsigned int> targets;
        std::vector<std::string> target_names;

        WindowGenerator(unsigned int window_length, unsigned int step_size);

        void initialize(
                std::istream &vcf_input,
                std::istream &pop_file,
                std::set<std::string> &target,
                std::set<std::string> &reference,
                std::set<std::string> &exclude);

        bool next_window();
};
