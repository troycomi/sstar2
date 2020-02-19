#include "sstar2/vcf_file.h"
#include <sstream>

std::string VcfEntry::to_str(void) const{
    std::ostringstream sstr;
    sstr << chromosome << '\t'
        << position << '\t'
        << reference << '\t'
        << alternative;
    for (auto b : haplotypes)
        sstr << '\t' << (b ? '1' : '0');
    sstr << '\n';
    return sstr.str();
}

short unsigned int VcfEntry::genotype(unsigned int individual) const{
    individual <<= 1;  // *= 2

    return haplotypes.at(individual) + haplotypes.at(individual + 1);
}

short unsigned int VcfEntry::gt_code(unsigned int individual) const{
    individual <<= 1;  // *= 2

    return haplotypes.at(individual) + (haplotypes.at(individual + 1) << 1);
}

bool VcfEntry::any_haplotype(const std::vector<unsigned int> &individuals) const{
    for (auto indiv : individuals){
        indiv <<= 1;  // *=2
        if (haplotypes[indiv] | haplotypes[indiv + 1])
            return true;
    }
    return false;
}

unsigned int VcfEntry::count_haplotypes(const std::vector<unsigned int> &individuals) const{
    unsigned int result = 0;
    for (auto indiv : individuals){
        indiv <<= 1;  // *=2
        result += haplotypes[indiv] + haplotypes[indiv + 1];
    }
    return result;
}

unsigned int VcfFile::initialize_individuals(const std::string line,
        const std::set<std::string> individuals){
    // matches individuals to the location in the vcf file
    // line is the header line of vcf with individual names
    std::stringstream stream(line);
    std::string token;
    unsigned int index = 0;
    while( std::getline(stream, token, '\t') ){
        if(index > 8 &&
                (individuals.empty() ||
                 individuals.find(token) != individuals.end()))
            individual_map.insert(std::pair<unsigned int, std::string>(index, token));
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
    unsigned int token = 0, haplo_count = 0;
    start = end = line;
    auto current_indiv = individual_map.begin();
    for(;;){
        // move end to next tab
        for(; *end != '\t' && *end != '\0'; ++end)  ;

        switch (token){
            case 0:  // chromosome 
                if(entry.chromosome.compare(0, end-start, start) != 0)
                    entry.chromosome = std::string(start, end-start);
                break;

            case 1:  // position 
                entry.position = std::stoi(start, nullptr);
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
                if (token == current_indiv->first){
                    if (*start == '.'){  // unknown (either . or ./.)
                        entry.haplotypes[haplo_count] = 0;
                        entry.haplotypes[haplo_count+1] = 0;
                    }
                    else{
                        if (*(start+1) != '|' && !warned_unphased){
                            std::cerr << "WARNING: Detected unphased "
                                "haplotype at chrom " << entry.chromosome <<
                                " and pos " << entry.position << "!\n";
                            warned_unphased = true;
                        }
                        entry.haplotypes[haplo_count] = (*start == '1');
                        entry.haplotypes[haplo_count+1] = (*(start + 2) == '1');
                    }
                    ++current_indiv;
                    haplo_count += 2;
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
