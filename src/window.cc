#include "sstar2/window.h"

bool WindowGT::operator==(const WindowGT& rhs) const{
    return position == rhs.position &&
        genotype == rhs.genotype;
}

void WindowBucket::reset_bucket(unsigned int start, unsigned int end){
    this->start = start;
    this->end = end;
    site_count = 0;
    reference_count = 0;
    for (auto &gt : genotypes)
        gt.clear();
}

StepWindow::StepWindow(unsigned int window_step, unsigned int window_length) :
    length(window_length), step(window_step) {
        if(step > length)
            throw std::invalid_argument("Window step can not be larger than length");
}

void StepWindow::initialize(unsigned int num_targets){
    // sets up the window structures for holding incoming vcf lines
    int num_buckets = length / step;
    num_buckets += (length % step == 0) ? 0 : 1;  //ceiling operation
    buckets.resize(num_buckets);
    // initilize size of genotype vector
    for (auto &bucket : buckets){
        bucket.genotypes.resize(num_targets);
    }
}

void StepWindow::start_window(VcfEntry &entry){
    // last line was a new chromosome
    if(entry.chromosome.compare(chromosome) != 0)
        reset(entry.chromosome);

    else  // starting new window
        next();
}

bool StepWindow::should_break(VcfEntry &entry) const {
    // new chromosome
    if(entry.chromosome.compare(chromosome) != 0)
        return true;

    // position is past end
    if (entry.position > end)
        return true;
    return false;
}

bool StepWindow::should_record(VcfEntry &entry) const {
    // for step windows, should always record once valid, otherwise would break
    return true;
}

void StepWindow::record(const VcfEntry &entry,
        const std::vector<unsigned int> &targets,
        unsigned int reference_haplotypes){
    // record the provided entry to the window

    // skip outside of start/end
    if(start >= entry.position || entry.position > end)
        return;

    // iterate in reverse since will mostly be inserting into last bucket
    for(auto bucket = buckets.rbegin(); bucket != buckets.rend(); ++bucket){
        // skip bucket outside of entry, will also stop recording outside region
        if(bucket->start >= entry.position || entry.position > bucket->end)
            continue;
        ++bucket->site_count;
        if(reference_haplotypes != 0){
            ++bucket->reference_count;
            return;  // only care about non-ref snps
        }
        short unsigned int gt;
        unsigned int indiv = 0;
        // for each target, if gt != 0, record position and gt
        for(unsigned int target : targets){
            if((gt = entry.genotypes[target]) != 0)
                bucket->genotypes[indiv].emplace_back(entry.position, gt);
            ++indiv;
        }
        return;
    }
}

unsigned int StepWindow::total_snps() const{
    unsigned int result = 0;
    for (const auto & bucket : buckets){
        result += bucket.site_count;
    }
    return result;
}

unsigned int StepWindow::reference_snps() const{
    unsigned int result = 0;
    for (const auto & bucket : buckets){
        result += bucket.reference_count;
    }
    return result;
}

unsigned int StepWindow::individual_snps(unsigned int individual) const{
    unsigned int result = 0;
    for (const auto & bucket : buckets){
        result += bucket.genotypes[individual].size();
    }
    return result;
}

void StepWindow::fill_genotypes(std::vector<WindowGT> &genotypes, int individual) const{
    for(const auto &bucket : buckets)
        genotypes.insert(genotypes.end(),
                bucket.genotypes[individual].begin(),
                bucket.genotypes[individual].end());
}

void StepWindow::reset(std::string &chrom){
    chromosome = chrom;
    start = 0;
    end = length;
    callable_bases.set(chrom, start, end);
    for(unsigned int i = 0; i < buckets.size(); ++i)
        buckets[i].reset_bucket(i*step, (i+1)*step);
}

void StepWindow::next(){
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

std::ostream& operator<<(std::ostream &strm, const Window &window){
    strm << "chrom: " << window.chromosome << "\tstart: "
        << window.start << "\tend: " << window.end << "\n";
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

RangedWindow::RangedWindow(unsigned int window_step, unsigned int window_length,
        std::istream &region_file) :
    StepWindow(window_step, window_length) {
        // regions is a tab delimited file with chrom, start, end regions
        std::string chrom, line;
        unsigned long start, end;
        while(std::getline(region_file, line)){
            std::stringstream lineStream(line);
            if(lineStream >> chrom >> start >> end)
                regions[chrom] =
                    std::pair<unsigned long, unsigned long>(start, end);
            else
                std::cerr << "Unable to parse line from region file\n" << line << '\n';
        }
    }

RangedWindow::RangedWindow(unsigned int window_step, unsigned int window_length,
        const std::string &region) :
    StepWindow(window_step, window_length) {
        // add constructor to take a string with [chrom:]start-end
        auto colon = region.find(':');
        colon = colon == -1 ? 0 : colon + 1;
        auto hyphen = region.find('-');
        if(hyphen == -1 || colon == 1)  // no hyphen or colon is first character
            throw std::invalid_argument("Region must be given as `[{chrom}:]{start}-{end}`");
        unsigned long start = std::stoul(region.substr(colon, hyphen-colon));
        unsigned long end = std::stoul(region.substr(hyphen+1));
        if(colon == 0){  // no chromosome provided
            window_start = start;
            window_end = end;
        }
        else{
            regions[region.substr(0, colon-1)] =
                std::pair<unsigned long, unsigned long>(start, end);
        }
    }

void RangedWindow::reset(std::string &chrom){
    if(regions.empty())
        chromosome = chrom;  // keep last start/end

    else{
        auto region = regions.find(chrom);
        if(region == regions.end())
            chromosome = "";  // keep last start/end
        else{
            chromosome = chrom;
            window_start = region->second.first;
            window_end = region->second.second;
        }
    }
    start = window_start;
    end = start + length;
    callable_bases.set(chrom, start, end);
    for(unsigned int i = 0; i < buckets.size(); ++i)
        buckets[i].reset_bucket(start + i*step, start + (i+1)*step);
}

void RangedWindow::next(){
    // TODO how to deal with end and buckets?
}

bool RangedWindow::should_break(VcfEntry &entry) const{
    return false;
}

bool RangedWindow::should_record(VcfEntry &entry) const{
    return false;
}
