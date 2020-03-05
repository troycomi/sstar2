#include "sstar2/window_generator.h"

WindowGenerator::WindowGenerator(unsigned int window_length,
        unsigned int step_size) :
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
    for (const auto &el : vcf_file.individual_map){
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
    // include minimal validators
    validators.push_back(std::unique_ptr<Validator>(
                new FixationValidator(targets, references)));
    // read first line
    next_line();
}

void WindowGenerator::initialize_window(){
    // sets up the window structures for holding incoming vcf lines
    int num_buckets = length / step;
    num_buckets += (length % step == 0) ? 0 : 1;  //ceiling operation
    window.buckets.resize(num_buckets);
    // initilize size of genotype vector
    unsigned int num_targets = targets.size();
    for (auto &bucket : window.buckets){
        bucket.genotypes.resize(num_targets);
    }
}

void WindowGenerator::add_validator(std::unique_ptr<Validator> validator){
    validators.push_back(std::move(validator));
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
        if(entry_is_valid()){
            unsigned int ref_haps = vcf_line.count_haplotypes(references);
            window.record(vcf_line, targets, ref_haps);
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

void Window::record(const VcfEntry &entry,
        const std::vector<unsigned int> &targets,
        unsigned int reference_haplotypes){
    // record the provided entry to the window
    // assumes in the range of (start, end], not fixed and present in ref or target

    // iterate in reverse since will mostly be inserting into last bucket
    for(auto bucket = buckets.rbegin(); bucket != buckets.rend(); ++bucket){
        // skip bucket outside of entry
        if(bucket->start >= entry.position || entry.position > bucket->end)
            continue;
        ++bucket->site_count;
        if(reference_haplotypes != 0){
            ++bucket->reference_count;
            return;  // only care about non-ref snps
        }
        short unsigned int gt;
        unsigned int indiv = 0;
        // for each target, if gt != 0, recorde position and gt
        for(unsigned int target : targets){
            if((gt = entry.genotypes[target]) != 0)
                bucket->genotypes[indiv].emplace_back(entry.position, gt);
            ++indiv;
        }
        return;
    }
}

unsigned int Window::total_snps() const{
    unsigned int result = 0;
    for (const auto & bucket : buckets){
        result += bucket.site_count;
    }
    return result;
}

unsigned int Window::reference_snps() const{
    unsigned int result = 0;
    for (const auto & bucket : buckets){
        result += bucket.reference_count;
    }
    return result;
}

unsigned int Window::individual_snps(unsigned int individual) const{
    unsigned int result = 0;
    for (const auto & bucket : buckets){
        result += bucket.genotypes[individual].size();
    }
    return result;
}

void Window::reset_window(std::string chrom,
        unsigned int length, unsigned int step){
    chromosome = chrom;
    start = 0;
    end = length;
    callable_bases.set(chrom, start, end);
    for(unsigned int i = 0; i < buckets.size(); ++i)
        buckets[i].reset_bucket(i*step, (i+1)*step);
}

void Window::step(unsigned int step){
    // increment start and end, prepare bucket
    start += step;
    buckets[0].reset_bucket(buckets.back().end,
            buckets.back().end + step);
    end += step;
    // update callable_bases
    callable_bases.set(chromosome, start, end);
    // cycle buckets
    std::rotate(buckets.begin(),
            buckets.begin()+1,
            buckets.end());
}

unsigned int WindowGenerator::callable_length(){
    for(auto const &validator : validators){
        validator->updateCallable(window.callable_bases);
    }
    return window.callable_bases.totalLength();
}

void WindowBucket::reset_bucket(unsigned int start, unsigned int end){
    this->start = start;
    this->end = end;
    site_count = 0;
    reference_count = 0;
    for (auto &gt : genotypes)
        gt.clear();
}

std::ostream& operator<<(std::ostream &strm, const Window &window){
    strm << "chrom: " << window.chromosome << "\tstart: "
        << window.start << "\tend: " << window.end << "\nBuckets:\n";
    for (unsigned int i = 0; i < window.buckets.size(); ++i)
        strm << "\t#" << i << ":\n" << window.buckets[i] << '\n';
    return strm;
}

std::ostream& operator<<(std::ostream &strm, const WindowBucket &bucket){
    strm << "\tstart: " << bucket.start << "\tend: " << bucket.end << "\n\twith "
        << bucket.site_count << " sites\t " << bucket.reference_count
        << " references\n";
    int indiv = 0;
    for (const auto &genotype : bucket.genotypes){
        strm << "\t\tindiv " << indiv << "\t";
        for (const auto &gt : genotype)
            strm << "(" << gt.position << ", " << gt.genotype << ")\t";
        strm << "\n";
        ++indiv;
    }
    return strm;
}

bool WindowGT::operator==(const WindowGT& rhs) const{
    return position == rhs.position &&
        genotype == rhs.genotype;
}
