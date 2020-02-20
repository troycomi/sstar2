#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <set>

#include "cxxopts/cxxopts.hpp"
#include "sstar2/sstar.h"
#include "sstar2/window_generator.h"

int main(int argc, char** argv)
{
    cxxopts::Options options("test", "A brief description");

    options.add_options()
        ("v,vcf", "Input vcf file; can accept input redirection",
         cxxopts::value<std::string>())
        ("p,popfile", "Population file; tsv with indiv, pop, superpop",
         cxxopts::value<std::string>())
        ("t,targets", "Comma separated list of target populations or individuals",
         cxxopts::value<std::vector<std::string>>())
        ("r,references", "Comma separated list of reference populations or individuals",
         cxxopts::value<std::vector<std::string>>())
        ("l,length", "Window length", cxxopts::value<unsigned int>()->default_value("50000"))
        ("s,step", "Window step", cxxopts::value<unsigned int>()->default_value("10000"))
        ("o,output", "Output file; can accept input redirection, default stdout",
         cxxopts::value<std::string>())
        // TODO add bonus and penalty options, excluded individuals
        ("h,help", "Print usage")
        ;

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    if (result.count("vcf") == 0){
        std::cerr << "Must specify an input vcf!\n";
        exit(EXIT_FAILURE);
    }

    if (result.count("popfile") == 0){
        std::cerr << "Must specify a population file!\n";
        exit(EXIT_FAILURE);
    }

    if (result.count("targets") == 0){
        std::cerr << "Must specify a target!\n";
        exit(EXIT_FAILURE);
    }

    if (result.count("references") == 0){
        std::cerr << "Must specify a reference!\n";
        exit(EXIT_FAILURE);
    }

    std::ifstream vcf, popfile;
    vcf.open(result["vcf"].as<std::string>());
    popfile.open(result["popfile"].as<std::string>());

    std::set<std::string> targets, references, excluded;
    for (auto t : result["targets"].as<std::vector<std::string>>())
        targets.insert(t);
            
    for (auto t : result["references"].as<std::vector<std::string>>())
        references.insert(t);

    unsigned int length, step;
    length = result["length"].as<unsigned int>();
    step = result["step"].as<unsigned int>();

    std::ofstream of;
    std::streambuf * buf;
    if (result.count("output")){
        of.open(result["output"].as<std::string>());
        buf = of.rdbuf();
    }
    else{
        buf = std::cout.rdbuf();
    }
    std::ostream output(buf);

    WindowGenerator generator(length, step);
    generator.initialize(vcf, popfile, targets, references, excluded);
    SStarCaller sstar;
    sstar.write_header(output);

    while (generator.next_window())
        sstar.write_window(output, generator);

    vcf.close();
    popfile.close();
    if (result.count("output"))
        of.close();
    return 0;
}
