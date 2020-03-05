#include "sstar2/validator.h"

void BaseRegions::add(const std::string &chrom, unsigned long start, unsigned long end){
    // add region, clearing if needed.  assumes sorted input with no overlap
    if(chrom != chromosome){
        set(chrom, start, end);
        return;
    }
    positions.push_back(start);
    positions.push_back(end);
}

void BaseRegions::set(const std::string &chrom, unsigned long start, unsigned long end){
    // set region to chromosome with a single entry
    chromosome = chrom;
    positions.clear();
    positions.push_back(start);
    positions.push_back(end);
}

unsigned long BaseRegions::totalLength(){
    // get total callable bases
    if(positions.empty())
        return 0;
    unsigned long result = 0;
    // this keeps start and end in step
    for(auto start = positions.begin(), end = std::next(start);
            start != positions.end();
            std::advance(start, 2), std::advance(end, 2))
        result += *end - *start;
    return result;
}

void BaseRegions::intersect(const BaseRegions &other){
    // assume both are strictly increasing with no overlap within segments
    if(other.positions.empty() || chromosome != other.chromosome){
        positions.clear();
        return;
    }
    auto mystart = positions.begin(), myend = std::next(mystart);
    auto otherstart = other.positions.begin(), otherend = std::next(otherstart);
    while(mystart != positions.end() && otherstart != other.positions.end()){
        if(*mystart >= *otherend){  // other too far behind
            std::advance(otherstart, 2);
            std::advance(otherend, 2);
            continue;
        }
        if(*otherstart >= *myend){  // mine too far behind
            // remove this as it doesn't match
            mystart = positions.erase(mystart, std::next(myend));
            myend = std::next(mystart);
            continue;
        }
        // now we have some overlap
        if(*otherstart > *mystart)  // start must be cut
            *mystart = *otherstart;
        // end can be used again...
        if(*otherend < *myend){
            // set this end to otherend, new start there as well
            positions.insert(myend, 2, *otherend);
            std::advance(mystart, 2);
            std::advance(otherstart, 2);
            std::advance(otherend, 2);
        }
        else{  // otherend > mine, go to next
            std::advance(mystart, 2);
            std::advance(myend, 2);
        }
    }
    if(mystart != positions.end())  // remove any remaining
        positions.erase(mystart, positions.end());
}

void BaseRegions::subtract(const BaseRegions &other){
    // assume both are strictly increasing with no overlap within segments
    // nothing to subtract
    if(other.positions.empty() || chromosome != other.chromosome){
        return;
    }
    auto mystart = positions.begin(), myend = std::next(mystart);
    auto otherstart = other.positions.begin(), otherend = std::next(otherstart);
    while(mystart != positions.end() && otherstart != other.positions.end()){
        if(*mystart >= *otherend){  // other too far behind
            std::advance(otherstart, 2);
            std::advance(otherend, 2);
            continue;
        }
        if(*otherstart >= *myend){  // mine too far behind
            std::advance(mystart, 2);
            std::advance(myend, 2);
            continue;
        }
        // now we have some overlap
        if(*otherstart > *mystart){  // start after mine
            // keep first part
            positions.insert(myend, 2, *otherstart);
            std::advance(mystart, 2);
        }

        if(*otherend < *myend){  // other end in region
            // move end
            *mystart = *otherend;
            // advance other, mine can be used
            std::advance(otherstart, 2);
            std::advance(otherend, 2);
        }
        else{  // end beyond, remove this consider other again
            mystart = positions.erase(mystart, std::next(myend));
            myend = std::next(mystart);
        }
    }
}

void BaseRegions::write(std::ostream &strm) const{
    if (chromosome == ""){
        strm << "No region\n";
    }
    else{
        strm << chromosome << ':';
        for (auto pos : positions)
            strm << pos << ',';
        strm << '\n';
    }
}

unsigned long BaseRegions::getEnd() const{
    if(positions.size() == 0)
        return 0;
    return positions.back();
}

std::ostream& operator<<(std::ostream &strm, const BaseRegions region){
    region.write(strm);
    return strm;
}

bool FixationValidator::isValid(const VcfEntry &entry){
    unsigned int targ_haps = entry.count_haplotypes(targets);
    unsigned int ref_haps = entry.count_haplotypes(references);
    if(targ_haps == 0 && ref_haps == 0)  // not found
        return false;
    if(targ_haps == targets.size()*2 &&
            ref_haps == references.size()*2)  // fixed
        return false;
    return true;
}

bool BedFile::inBed(std::string chrom, unsigned long position){
    // test if chrom/position is in bed file
    // assume queries are sorted in same order as bed file!
    // true if position is in (start, end]
    for(;;){
        // current location is beyond query
        if(chrom < chromosome)
            return false;
        if(chrom == chromosome){
            if(position <= start)
                return false;
            else if(end < position)
                readline();
            else
                return true;  // start < position <= end
        }
        if(chrom > chromosome)
            readline();
        if(chromosome == "")
            return false;
    }
}

void BedFile::readline(){
    std::string line;
    if(std::getline(*bedfile, line)){
        std::istringstream iss(line);
        if (!(iss >> chromosome >> start >> end))
            chromosome = "";
        else
            regions.add(chromosome, start, end);
    }
    else
        chromosome = "";
}

void BedFile::intersect(BaseRegions &callable){
    callable.intersect(regions);
}

void BedFile::subtract(BaseRegions &callable){
    callable.subtract(regions);
}

bool PositiveBedValidator::isValid(const VcfEntry &entry){
    return bedfile.inBed(entry.chromosome, entry.position);
}

void PositiveBedValidator::updateCallable(BaseRegions &callable){
    // need to check if end is in bedfile to force it to read through end
    bedfile.inBed(callable.getChromosome(), callable.getEnd());
    bedfile.intersect(callable);
}

bool NegativeBedValidator::isValid(const VcfEntry &entry){
    return ! bedfile.inBed(entry.chromosome, entry.position);
}

void NegativeBedValidator::updateCallable(BaseRegions &callable){
    // need to check if end is in bedfile to force it to read through end
    bedfile.inBed(callable.getChromosome(), callable.getEnd());
    bedfile.subtract(callable);
}

