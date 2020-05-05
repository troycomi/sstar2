#include "sstar2/validator.h"

void BaseRegions::add(const std::string &chrom, unsigned long start, unsigned long end){
    // add region.  assumes sorted input with no overlap
    if(positions.count(chrom) == 0)
        positions.insert(std::pair<std::string, std::list<unsigned long>>(
                    chrom, {start, end}));
    else{
        positions[chrom].push_back(start);
        positions[chrom].push_back(end);
    }
}

void BaseRegions::set(const std::string &chrom, unsigned long start, unsigned long end){
    // set region to chromosome with a single entry
    positions.clear();
    add(chrom, start, end);
}

unsigned long BaseRegions::totalLength(){
    // get total callable bases
    if(positions.empty())
        return 0;
    unsigned long result = 0;
    // this keeps start and end in step
    for(auto &pair : positions){
        for(auto start = pair.second.begin(), end = std::next(start);
                start != pair.second.end();
                std::advance(start, 2), std::advance(end, 2))
            result += *end - *start;
    }
    return result;
}

void BaseRegions::intersect(const BaseRegions &other){
    // assume both are strictly increasing with no overlap within segments
    for (auto &pos : positions){
        auto otherpos = other.positions.find(pos.first);
        // not in other or is empty
        if(otherpos == other.positions.end() || otherpos->second.empty()){
            pos.second.clear();
            return;
        }
        auto mystart = pos.second.begin(), myend = std::next(mystart);
        auto otherstart = otherpos->second.begin(), otherend = std::next(otherstart);
        while(mystart != pos.second.end()
                && otherstart != otherpos->second.end()){
            if(*mystart >= *otherend){  // other too far behind
                std::advance(otherstart, 2);
                std::advance(otherend, 2);
                continue;
            }
            if(*otherstart >= *myend){  // mine too far behind
                // remove this as it doesn't match
                mystart = pos.second.erase(mystart, std::next(myend));
                myend = std::next(mystart);
                continue;
            }
            // now we have some overlap
            if(*otherstart > *mystart)  // start must be cut
                *mystart = *otherstart;
            // end can be used again...
            if(*otherend < *myend){
                // set this end to otherend, new start there as well
                pos.second.insert(myend, 2, *otherend);
                std::advance(mystart, 2);
                std::advance(otherstart, 2);
                std::advance(otherend, 2);
            }
            else{  // otherend > mine, go to next
                std::advance(mystart, 2);
                std::advance(myend, 2);
            }
        }
        if(mystart != pos.second.end())  // remove any remaining
            pos.second.erase(mystart, pos.second.end());
    }
}

void BaseRegions::subtract(const BaseRegions &other){
    // assume both are strictly increasing with no overlap within segments
    for (auto &pos : positions){
        auto otherpos = other.positions.find(pos.first);
        // not in other or is empty
        // nothing to subtract
        if(otherpos == other.positions.end() || otherpos->second.empty())
            continue;
    
        auto mystart = pos.second.begin(), myend = std::next(mystart);
        auto otherstart = otherpos->second.begin(), otherend = std::next(otherstart);
        while(mystart != pos.second.end() && otherstart != otherpos->second.end()){
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
                pos.second.insert(myend, 2, *otherstart);
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
                mystart = pos.second.erase(mystart, std::next(myend));
                myend = std::next(mystart);
            }
        }
    }
}

const std::string BaseRegions::getChromosome() const{
    return positions.begin()->first;
}

unsigned long BaseRegions::getEnd(const std::string &chromosome){
    if(positions.count(chromosome) == 0 || positions[chromosome].empty())
        return 0;
    return positions[chromosome].back();
}

bool BaseRegions::inRegion(const std::string &chromosome, unsigned long position){
    if(positions.count(chromosome) == 0)
        return false;
    // TODO search backward, exit if past possible
    // search through positions for one in (start, end]
    for(auto end = positions[chromosome].rbegin(), start = std::next(end);
            end != positions[chromosome].rend();
            std::advance(start, 2), std::advance(end, 2)){
        if(*start < position && position <= *end)
            return true;
        else if(position > *end)
            return false;
    }
        
    return false;

}

void BaseRegions::write(std::ostream &strm) const{
    if (positions.empty()){
        strm << "No region\n";
    }
    else{
        for(auto &position : positions){
            strm << position.first << ':';
            for (auto &pos : position.second)
                strm << pos << ',';
            strm << '\n';
        }
    }
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
        if(chrom == chromosome){
            if(position <= start)
                return regions.inRegion(chrom, position);
            else if(end < position)
                readline();
            else
                return true;  // start < position <= end
        }
        else if(chromosome == "")
            return regions.inRegion(chrom, position);
        else
            readline();
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
    auto chrom = callable.getChromosome();
    bedfile.inBed(chrom, callable.getEnd(chrom));
    bedfile.intersect(callable);
}

bool NegativeBedValidator::isValid(const VcfEntry &entry){
    return ! bedfile.inBed(entry.chromosome, entry.position);
}

void NegativeBedValidator::updateCallable(BaseRegions &callable){
    // need to check if end is in bedfile to force it to read through end
    auto chrom = callable.getChromosome();
    bedfile.inBed(chrom, callable.getEnd(chrom));
    bedfile.subtract(callable);
}

