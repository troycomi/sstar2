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
    std::vector<uint8_t> genotypes;

    VcfEntry(std::string chrom, size_t individuals) :
        chromosome(chrom), genotypes(individuals) {}
    std::string to_str() const;
    bool any_haplotype(const std::vector<unsigned int> &individuals) const;
    unsigned int count_haplotypes(const std::vector<unsigned int> &individuals) const;
};

class VcfFile{
    // if user has been warned about unphased data
    bool warned_unphased = false;
    std::vector<unsigned int> individual_indices;

    public:
        // map of individual to position in vcf file
        std::map <unsigned int, std::string> individual_map;

        unsigned int initialize_individuals(const std::string line,
                const std::set<std::string> &individuals);
        VcfEntry initialize_entry();
        bool parse_line(const char* line, VcfEntry &entry);
};
