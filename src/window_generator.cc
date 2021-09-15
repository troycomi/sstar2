#include "sstar2/window_generator.h"

void WindowGenerator::initialize(
                std::istream &vcf_input,
                std::istream &pop_file,
                std::set<std::string> &target,
                std::set<std::string> &reference,
                std::set<std::string> &exclude){
    // sets up vcf file and population data to prepare for returning
    population.read_data(pop_file, target, reference, exclude);
    vcf = &vcf_input;
    initialize_vcf();
    window->initialize(targets.size());
}

void WindowGenerator::initialize_vcf(){
    // initialize vcf_file and vcf_line
    std::set<std::string> individuals;
    individuals.insert(population.targets.begin(), population.targets.end());
    individuals.insert(population.references.begin(), population.references.end());
    individuals.insert(population.excluded.begin(), population.excluded.end());
    while(std::getline(*vcf, vcf_string)){
        if(vcf_string[1] != '#'){
            // setup indiviual mapping
            vcf_file.initialize_individuals(vcf_string, individuals);
            // setup entry
            vcf_line = vcf_file.initialize_entry();
            break;
        }
    }
    // setup index into haplotype for each input
    // store target names because those will be needed later
    int ind = 0;
    bool found;
    for (const auto &el : vcf_file.individual_map){
        found = false;
        if(population.targets.find(el.second) != population.targets.end()){
            targets.push_back(ind);
            target_names.push_back(el.second);
            population_names.push_back(population.target_to_population.at(el.second));
            found = true;
        }
        if(population.references.find(el.second) != population.references.end()){
            references.push_back(ind);
            found = true;
        }
        if(population.excluded.find(el.second) != population.excluded.end()){
            excluded.push_back(ind);
            found = true;
        }
        if (!found)
            throw std::invalid_argument("VCF file is yeilding extra individuals");
        ++ind;
    }
    // include minimal validators
    validators.push_back(std::unique_ptr<Validator>(
                new FixationValidator(targets, references)));
    // read first line
    next_line();
}

void WindowGenerator::add_validator(std::unique_ptr<Validator> validator){
    validators.push_back(std::move(validator));
}

bool WindowGenerator::next_window(){
    // updates public facing window to next value
    // at the start, vcf_line is either empty or contains the last,
    // unprocessed line

    if(terminated){
        // free validators
        for(auto & validator : validators)
            validator.reset();
        window.reset();
        return false;  // all done
    }

    window->start_window(vcf_line);

    do {
        if(window->should_break(vcf_line)){
            return true;
        }

        // validate other properties
        if(window->should_record(vcf_line) && entry_is_valid()){
            unsigned int ref_haps = vcf_line.count_haplotypes(references);
            window->record(vcf_line, targets, ref_haps);
        }

    }while(next_line());
    // at this point, no more lines are available, but the window is valid
    // need to record no lines are left and return false the next time...
    terminated = true;  // for next time
    return true;
}

bool WindowGenerator::next_line(){
    // updates vcf_line to new value, returns true when the line is valid
    // false when the end of file was reached
    while(std::getline(*vcf, vcf_string)){
        if(! vcf_file.parse_line(vcf_string.c_str(), vcf_line))
            continue;
        if(vcf_line.any_haplotype(excluded))
            continue;
        return true;
    }
    return false;
}

bool WindowGenerator::entry_is_valid(){
    for (auto &validator : validators)
        if(! validator->isValid(vcf_line))
            return false;
    return true;
}

unsigned int WindowGenerator::callable_length(){
    for(auto const &validator : validators){
        validator->updateCallable(window->callable_bases);
    }
    return window->callable_bases.totalLength();
}
