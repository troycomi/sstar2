#include "sstar2/validator.h"

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
