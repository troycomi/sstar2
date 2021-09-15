#pragma once
#include <deque>
#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>
#include <sstream>
#include "sstar2/validator.h"

struct WindowGT {
    unsigned long position;
    short unsigned int genotype;  // encoded haplotypes, 3 -> homozygous alt
    bool operator==(const WindowGT& rhs) const;

    WindowGT(unsigned long pos, short unsigned int gt):
        position(pos), genotype(gt){};
};

// TODO need to make a subclass flat window to handle non-step motions
// TODO support for arbitrary regions can use this window with modifications
// TODO the next_window code is pretty tightly coupled to the window class...

class Window {
    public:
        std::string chromosome;
        // contains positions in (start, end]
        unsigned long start=0, end=0;
        BaseRegions callable_bases;

        // called once on start of next window so can prepare
        virtual void start_window(VcfEntry &entry) = 0;
        // called once per line, return true if this line is
        // beyond the window and loop should break
        virtual bool should_break(VcfEntry &entry) const = 0;
        // called after should break, return true if the line should
        // be recorded after other validators run
        virtual bool should_record(VcfEntry &entry) const = 0;
        virtual void record(const VcfEntry &entry,
                const std::vector<unsigned int> &targets,
                unsigned int reference_haplotypes) = 0;

        virtual void fill_genotypes(std::vector<WindowGT> &genotypes, int individual) const = 0;
        virtual unsigned int total_snps() const = 0;
        virtual unsigned int reference_snps() const = 0;
        virtual unsigned int individual_snps(unsigned int individual) const = 0;
        virtual void initialize(unsigned int num_targets) = 0;
        virtual ~Window() = default;
};

struct WindowBucket {
    unsigned long start=0, end=0;
    unsigned int site_count=0,  // number of snps after removing excluded, fixed, and homozygous ref
        reference_count=0;  // number in any ref
    std::vector<std::vector<WindowGT>> genotypes;
    void reset_bucket(unsigned int start, unsigned int end);
};

// A simple, concrete window that yields a given length and step over all
// positions.
class StepWindow : public Window {
    virtual void reset(std::string &chrom);
    virtual void next();

    protected:
        std::deque<WindowBucket> buckets;
        unsigned int step, length;

    public:
        StepWindow(unsigned int window_step, unsigned int window_length);
        void initialize(unsigned int num_targets);
        void start_window(VcfEntry &entry);
        bool should_break(VcfEntry &entry) const;
        bool should_record(VcfEntry &entry) const;
        void record(const VcfEntry &entry, const std::vector<unsigned int> &targets,
                unsigned int reference_haplotypes);
        void fill_genotypes(std::vector<WindowGT> &genotypes, int individual) const;
        unsigned int total_snps() const;
        unsigned int reference_snps() const;
        unsigned int individual_snps(unsigned int individual) const;
};

// a stepping window that has a set of start/end positions
// similar to step windows but the next and reset code has additional checks
class RangedWindow : public StepWindow {
    std::map<std::string, std::pair<unsigned long, unsigned long>> regions;
    unsigned long window_start = 0, window_end = 0;
    void reset(std::string &chrom);
    void next();

    public:
        RangedWindow(unsigned int window_step, unsigned int window_length,
                std::istream &regions);
        RangedWindow(unsigned int window_step, unsigned int window_length,
                const std::string &region);
        bool should_break(VcfEntry &entry) const;
        bool should_record(VcfEntry &entry) const;
};

std::ostream& operator<<(std::ostream &strm, const Window &window);
std::ostream& operator<<(std::ostream &strm, const WindowBucket &bucket);
