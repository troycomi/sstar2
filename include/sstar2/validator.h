#pragma once
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <map>

#include "sstar2/vcf_file.h"

class BaseRegions{
    // contains a chromosome pointing to a list of start/end positions
    // works to store a bed file or callable bases
    std::map<std::string, std::list<unsigned long>> positions;

    public:
        // add region, clearing if needed.  Handles overlaps, assumes sorted input
        void add(const std::string &chrom, unsigned long start, unsigned long end);
        // set region to chromosome with a single entry
        void set(const std::string &chrom, unsigned long start, unsigned long end);
        // get total callable bases
        unsigned long totalLength();
        void intersect(const BaseRegions &other);
        void subtract(const BaseRegions &other);
        void write(std::ostream &strm) const;
        // get the first chromosome
        const std::string getChromosome() const;
        unsigned long getEnd(const std::string chromosome);
        bool inRegion(const std::string chromosome, unsigned long position);
};
std::ostream& operator<<(std::ostream &strm, const BaseRegions region);

class Validator{
    public:
        virtual bool isValid(const VcfEntry &entry) = 0;
        virtual void updateCallable(BaseRegions &callable) = 0;
};

class FixationValidator : public Validator{
    std::vector<unsigned int> targets;
    std::vector<unsigned int> references;

    public:
        FixationValidator(
                std::vector<unsigned int> target_inds,
                std::vector<unsigned int> reference_inds) :
            targets(target_inds), references(reference_inds) {};

        bool isValid(const VcfEntry &entry);
        void updateCallable(BaseRegions &callable){}
};

class BedFile{
    std::istream *bedfile;
    std::string chromosome;
    unsigned long start, end;
    void readline();
    BaseRegions regions;

    public:
        BedFile(std::istream *file) : bedfile(file), chromosome("*") {};
        bool inBed(std::string chrom, unsigned long position);
        void intersect(BaseRegions &callable);
        void subtract(BaseRegions &callable);
};

class PositiveBedValidator : public Validator{
    BedFile bedfile;

    public:
        PositiveBedValidator(std::istream *file) :
            bedfile(file) {};

        bool isValid(const VcfEntry &entry);
        void updateCallable(BaseRegions &callable);
};

class NegativeBedValidator : public Validator{
    BedFile bedfile;

    public:
        NegativeBedValidator(std::istream *file) :
            bedfile(file) {};

        bool isValid(const VcfEntry &entry);
        void updateCallable(BaseRegions &callable);
};
