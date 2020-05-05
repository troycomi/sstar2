#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <set>

#include <CLI/CLI.hpp>

#include "sstar2/sstar.h"
#include "sstar2/window_generator.h"
#include "sstar2/validator.h"

int main(int argc, char** argv)
{
    CLI::App app{"Fast, lean sstar rewrite"};

    std::string vcf_file;
    app.add_option("-v,--vcf", vcf_file,
            "Input vcf file; can accept input redirection")
        ->required()->check(CLI::ExistingFile);
    std::string popfile;
    app.add_option("-p,--popfile", popfile,
            "Population file; tsv with indiv, pop, superpop")
        ->required()->check(CLI::ExistingFile);

    std::vector<std::string> targets, references, excluded;
    app.add_option("-t,--targets", targets,
            "Comma separated list of target populations or individuals")
        ->required()->delimiter(',');

    app.add_option("-r,--references", references,
            "Comma separated list of reference populations or individuals")
        ->required()->delimiter(',');

    app.add_option("-e,--excluded", excluded,
            "Comma separated list of excluded populations or individuals; default to none excluded")
        ->delimiter(',');

    unsigned int length = 50000, step = 10000;
    app.add_option("-l,--length", length, "Window length; default 50,000");
    app.add_option("-s,--step", step, "Window step; default 10,000");

    int bonus = 5000, penalty = -10000;
    app.add_option("--match-bonus", bonus, "Match bonus for sstar; default 5000");
    app.add_option("--mismatch-penalty", penalty, "Mismatch penalty for sstar; default -10000");

    std::string positiveBed = "";
    app.add_option("--include-bed", positiveBed,
            "Bed file with regions to include")
        ->check(CLI::ExistingFile);

    std::string negativeBed = "";
    app.add_option("--exclude-bed", negativeBed,
            "Bed file with regions to exclude")
        ->check(CLI::ExistingFile);

    std::string outfile = "-";
    app.add_option("-o,--output", outfile,
            "Output file; can accept input redirection; default stdout");

    CLI11_PARSE(app, argc, argv);

    std::ifstream vcf, popdata, posBed, negBed;
    vcf.open(vcf_file);
    popdata.open(popfile);

    std::set<std::string> target_set, reference_set, excluded_set;
    for (const auto &indiv : targets)
        target_set.insert(indiv);
            
    for (const auto &indiv : references)
        reference_set.insert(indiv);
            
    for (const auto &indiv : excluded)
        excluded_set.insert(indiv);

    std::ofstream of;
    std::streambuf * buf;
    if (outfile.compare("-") == 0){
        buf = std::cout.rdbuf();
    }
    else{
        of.open(outfile);
        buf = of.rdbuf();
    }
    std::ostream output(buf);

    WindowGenerator generator(std::unique_ptr<Window>(new StepWindow(step, length)));
    generator.initialize(vcf, popdata, target_set, reference_set, excluded_set);

    // add validators
    if(positiveBed != ""){
        posBed.open(positiveBed);
        generator.add_validator(std::unique_ptr<Validator>(
                    new PositiveBedValidator(&posBed)));
    }

    if(negativeBed != ""){
        negBed.open(negativeBed);
        generator.add_validator(std::unique_ptr<Validator>(
                    new NegativeBedValidator(&negBed)));
    }

    SStarCaller sstar{bonus, penalty};
    sstar.write_header(output);

    while (generator.next_window())
        sstar.write_window(output, generator);

    vcf.close();
    popdata.close();
    if(of.is_open())
        of.close();
    return 0;
}
