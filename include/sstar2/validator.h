#pragma once
#include "sstar2/vcf_file.h"

class Validator{
    public:
        virtual bool isValid(const VcfEntry &entry) = 0;
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
};
