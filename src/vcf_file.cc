#include "sstar2/vcf_file.h"
#include <sstream>

std::string VcfEntry::to_str(void) const{
    std::ostringstream sstr;
    sstr << chromosome << '\t'
        << position << '\t'
        << reference << '\t'
        << alternative;
    for (const auto &b : genotypes)
        sstr << '\t' << b;
    sstr << '\n';
    return sstr.str();
}

bool VcfEntry::any_haplotype(const std::vector<unsigned int> &individuals) const{
    for (const auto &indiv : individuals){
        if (genotypes[indiv])
            return true;
    }
    return false;
}

unsigned int VcfEntry::count_haplotypes(const std::vector<unsigned int> &individuals) const{
    unsigned int result = 0;
    for (const auto &indiv : individuals){
        result += genotypes[indiv] - (genotypes[indiv] >> 1);
        // this maps 0 to 0, 1 to 1, 2 to 1 and 3 to 1
    }
    return result;
}

unsigned int VcfFile::initialize_individuals(const std::string &line,
        const std::set<std::string> &individuals){
    // matches individuals to the location in the vcf file
    // line is the header line of vcf with individual names
    std::stringstream stream(line);
    std::string token;
    unsigned int index = 0;
    while( std::getline(stream, token, '\t') ){
        if(index > 8 &&
                (individuals.empty() ||
                 individuals.find(token) != individuals.end())){
            individual_map.insert(std::pair<unsigned int, std::string>(index, token));
            individual_indices.push_back(index);
        }
        ++index;
    }
    return individual_map.size();
}

VcfEntry VcfFile::initialize_entry(){
    return VcfEntry("", individual_map.size());
}

bool VcfFile::parse_line(const char* line, VcfEntry &entry){
    // returns true if line is updated
    const char *start, *end;
    unsigned int token = 0, geno_count = 0;
    start = end = line;
    auto current_indiv = individual_indices.begin();
    for(;;){
        // move end to next tab
        end = start + strcspn(start, "\t");

        switch (token){
            case 0:  // chromosome 
                if(entry.chromosome.compare(0, end-start, start) != 0)
                    entry.chromosome = std::string(start, end-start);
                break;

            case 1:  // position 
                // faster version of? entry.position = std::stoi(start, nullptr);
                entry.position = 0;
                for( ; start != end; ++start)
                    entry.position = (*start - '0') + entry.position * 10;
                break;

            case 3:  // ref
                if ((end - start) > 1)
                    return false;
                entry.reference = *start;
                break;

            case 4:  // alt
                if ((end - start) > 1)
                    return false;
                entry.alternative = *start;
                break;

            case 2:  // ID 
            case 5:  // QUAL 
            case 6:  // FILTER
            case 7:  // INFO
                break;

            case 8:  // FORMAT
                if(strncmp("GT", start, 2) != 0)
                    throw std::invalid_argument("FORMAT must start with GT");
                break;

            default:  // individuals
                if (current_indiv != individual_indices.end()
                        && token == *current_indiv){
                    if (*start == '.'){  // unknown (either . or ./.)
                        entry.genotypes[geno_count] = 0;
                    }
                    else{
                        if (*(start+1) != '|' && !warned_unphased){
                            std::cerr << "WARNING: Detected unphased "
                                "haplotype at chrom " << entry.chromosome <<
                                " and pos " << entry.position << "!\n";
                            warned_unphased = true;
                        }
                        entry.genotypes[geno_count] = (*start == '1') + 
                            ((*(start + 2) == '1') << 1);
                    }
                    ++current_indiv;
                    ++geno_count;
                }
                break;
        }
        // break if at end, else increment to start next token
        if(*end == '\0')
            break;
        start = ++end;
        ++token;
    }
    return true;
}
