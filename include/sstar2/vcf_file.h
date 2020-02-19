// class for representing lines of a vcf for chosen individuals
// includes chrom, pos, ref and alt
// only keeps snps and only parses format == GT
// Genotypes are assumed phased, .'s become 0's stored as haplotypes

#pragma once
#include <vector>
#include <string>
#include <map>
#include <set>
#include <iostream>
#include <string.h>

struct VcfEntry{
    std::string chromosome;
    unsigned long int position;
    char reference;
    char alternative;
    std::vector<uint8_t> haplotypes;

    VcfEntry(std::string chrom, size_t individuals) :
        chromosome(chrom), haplotypes(individuals * 2) {}
    std::string to_str() const;
    short unsigned int genotype(unsigned int individual) const;
    short unsigned int gt_code(unsigned int individual) const;
    bool any_haplotype(const std::vector<unsigned int> &individuals) const;
    unsigned int count_haplotypes(const std::vector<unsigned int> &individuals) const;
};

class VcfFile{
    // if user has been warned about unphased data by this
    bool warned_unphased = false;

    public:
        // map of individual to position in vcf file
        std::map <unsigned int, std::string> individual_map;

        unsigned int initialize_individuals(const std::string line,
                const std::set<std::string> individuals);
        VcfEntry initialize_entry();
        bool parse_line(const char* line, VcfEntry &entry);
};
