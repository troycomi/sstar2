// class for yielding windows of a vcf file
// interacts with a VcfFile to read lines
// and PopulationData to determine which files are targets

#pragma once
#include <set>
#include <vector>
#include "sstar2/vcf_file.h"
#include "sstar2/population_data.h"
#include "sstar2/validator.h"
#include "sstar2/window.h"

class WindowGenerator{

    bool terminated = false;

    std::istream *vcf = nullptr;
    std::string vcf_string;
    std::vector<unsigned int> references;
    std::vector<unsigned int> excluded;
    std::vector<std::unique_ptr<Validator>> validators;

    void initialize_vcf();
    bool next_line();
    bool entry_is_valid();

    public:
        VcfFile vcf_file;
        VcfEntry vcf_line{"", 0};
        PopulationData population;
        std::unique_ptr<Window> window;
        std::vector<unsigned int> targets;
        std::vector<std::string> target_names;
        std::vector<std::string> population_names;

        WindowGenerator(std::unique_ptr<Window> wind) : window(std::move(wind)){};

        void initialize(
                std::istream &vcf_input,
                std::istream &pop_file,
                std::set<std::string> &target,
                std::set<std::string> &reference,
                std::set<std::string> &exclude);
        void add_validator(std::unique_ptr<Validator> validator);
        unsigned int callable_length();
        bool next_window();
};
