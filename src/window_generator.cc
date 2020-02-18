#include "sstar2/window_generator.h"

WindowGenerator::WindowGenerator(int window_length, int step_size) :
    length(window_length), step(step_size) {
        if(step > length)
            throw std::invalid_argument("Window step can not be larger than length");
    };

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
    initialize_window();
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
    for (auto el : vcf_file.individual_map){
        found = false;
        if(population.targets.find(el.second) != population.targets.end()){
            targets.push_back(ind);
            target_names.push_back(el.second);
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
    // read first line
    next_line();
}

void WindowGenerator::initialize_window(){
    // sets up the window structures for holding incoming vcf lines
    int num_buckets = length / step;
    num_buckets += (length % step == 0) ? 0 : 1;  //ceiling operation
    window.buckets.resize(num_buckets);
    // initilize size of genotype vector
    int num_targets = targets.size();
    for (int i = 0; i < num_buckets; ++i){
        window.buckets[i].genotypes.resize(num_targets);
    }
}

bool WindowGenerator::next_window(){
    // updates public facing window to next value
    // at the start, vcf_line is either empty or contains the last,
    // unprocessed line

    if(terminated)
        return false;  // all done

    // last line was a new chromosome
    if(vcf_line.chromosome.compare(window.chromosome) != 0)
        window.reset_window(vcf_line.chromosome, length, step);

    else  // starting new window
        window.step(step);

    do {
        // new chromosome
        if(vcf_line.chromosome.compare(window.chromosome) != 0)
            return true;

        // position is past end
        if (vcf_line.position > window.end)
            return true;

        // validate other properties
        int targ_haps = vcf_line.count_haplotypes(targets);
        int ref_haps = vcf_line.count_haplotypes(references);
        if(targ_haps == 0 && ref_haps == 0)  // not found
            continue;
        if(targ_haps == targets.size()*2 &&
                ref_haps == references.size()*2)  // fixed
            continue;

        std::cout << "VALID: " << vcf_line.to_str();
        // TODO make a record function
        // test private methods (probably ok to be public)
        // test more of the struct methods
    }while(next_line());
    // at this point, no more lines are available, but the window is valid
    // need to record no lines are left and return false the next time...
    terminated = true;  // for next time
    return true;
}

bool WindowGenerator::next_line(){
    while(std::getline(*vcf, vcf_string)){
        if(! vcf_file.parse_line(vcf_string.c_str(), vcf_line))
            continue;
        if(vcf_line.any_haplotype(excluded))
            continue;
        return true;
    }
    return false;
}

void Window::reset_window(std::string chrom, int length, int step){
    chromosome = chrom;
    start = 0;
    end = length;
    for(int i = 0; i < buckets.size(); i++)
        buckets[i].reset_bucket(i*step, (i+1)*step);
}

void Window::step(int step){
    // increment start and end, prepare bucket
    start += step;
    buckets[0].reset_bucket(end, end + step);
    end += step;
    // cycle buckets
    std::rotate(buckets.begin(),
            buckets.begin()+1,
            buckets.end());
}

void WindowBucket::reset_bucket(int start, int end){
    this->start = start;
    this->end = end;
    site_count = 0;
    reference_count = 0;
    positions.clear();
    for (int i = 0; i < genotypes.size(); i++)
        genotypes[i].clear();
}

std::ostream& operator<<(std::ostream &strm, const Window &window){
    strm << "chrom: " << window.chromosome << "\tstart: "
        << window.start << "\tend: " << window.end << "\nBuckets:\n";
    for (int i = 0; i < window.buckets.size(); ++i)
        strm << "\t#" << i << ":\n" << window.buckets[i] << '\n';
    return strm;
}

std::ostream& operator<<(std::ostream &strm, const WindowBucket &bucket){
    strm << "\tstart: " << bucket.start << "\tend: " << bucket.end << "\nwith "
        << bucket.site_count << " sites\t " << bucket.reference_count
        << " references\n";
    // for each position
    for (int i = 0; i < bucket.positions.size(); ++i){
        strm << '\t' << bucket.positions[i];
        // for each individual
        for (int j = 0; j < bucket.genotypes.size(); j++){
            bool found = false;
            // for each GT entry
            for (auto gt : bucket.genotypes[j])
                if (gt.position_index == i){
                    found = true;
                    strm << gt.genotype;
                    break;
                }

            if (! found)
                strm << '0';
            strm << ' ';
        }
        strm << "\n";
    }
    return strm;
}
