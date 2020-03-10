#include "sstar2/sstar.h"

void SStarCaller::write_header(std::ostream &output){
    output << 
        "chrom\twinstart\twinend\t"
        "n_snps\tn_ind_snps\tn_region_ind_snps\t"
        "ind_id\tpop\t"
        "s_star\tnum_s_star_snps\ts_star_snps\t"
        "hap_1_s_start\thap_1_s_end\t"
        "hap_2_s_start\thap_2_s_end\t"
        "s_start\ts_end\t"
        "n_s_star_snps_hap1\tn_s_star_snps_hap2\t"
        "s_star_haps\t"
        "callable_bases\n";
}

void SStarCaller::write_window(std::ostream &output,
                WindowGenerator &generator){
    unsigned int total_snps = generator.window.total_snps();
    unsigned int ref_snps = generator.window.reference_snps();
    unsigned int callable = generator.callable_length();
    for(unsigned int i = 0; i < generator.targets.size(); ++i){
        unsigned int indiv_snps = generator.window.individual_snps(i);
        output << generator.window.chromosome << '\t'
            << generator.window.start << '\t'
            << generator.window.end << '\t'
            << total_snps << '\t'
            << indiv_snps << '\t'
            << indiv_snps + ref_snps << '\t'
            << generator.target_names[i] << '\t'
            << generator.population_names[i] << '\t';
        if (indiv_snps <= 2)
            output << emptyline;
        else{
            // build genotypes vector
            std::vector<WindowGT> genotypes;
            genotypes.reserve(indiv_snps);
            generator.window.fill_genotypes(genotypes, i);
            long s_score = sstar(genotypes);
            output << s_score << '\t'
                << genotypes.size() << '\t';

            unsigned long hap1start=0, hap2start=0, hap1end=0, hap2end=0;
            int hap1count = 0, hap2count = 0;
            bool first = true;
            for (const auto &gt : genotypes){
                // print positions joined on a comma
                if (first)
                    first = !first;
                else
                    output << ',';
                output << gt.position;

                // update haplo start and end
                if(gt.genotype == 1 || gt.genotype == 3){
                    ++hap1count;
                    hap1end = gt.position;
                    if(hap1start == 0)
                        hap1start = gt.position;
                }
                if(gt.genotype == 2 || gt.genotype == 3){
                    ++hap2count;
                    hap2end = gt.position;
                    if(hap2start == 0)
                        hap2start = gt.position;
                }
            }
            // if only contains one snp, set to 0
            if (hap1start == hap1end)
                hap1start = hap1end = 0;
            if (hap2start == hap2end)
                hap2start = hap2end = 0;

            output << '\t' << hap1start << '\t'
                << hap1end << '\t'
                << hap2start << '\t'
                << hap2end << '\t'
                << genotypes.front().position << '\t'  // s start
                << genotypes.back().position << '\t'  // s end
                << hap1count << '\t'
                << hap2count << '\t';
            // write comma joined haplotypes
            first = true;
            for (const auto &gt : genotypes){
                if (first)
                    first = !first;
                else
                    output << ',';
                output << gt.genotype;
            }
            output << '\t';
        }
        output << callable << '\n';
    }
}

long SStarCaller::sstar(std::vector<WindowGT> &genotypes){
    size_t nsnps = genotypes.size();
    // start with 10 mismatches as no-score without worring about overflow
    std::vector<int> scores(nsnps, mismatch_penalty*10);
    long new_score, append_score, bp_dist;
    // using snps as a 2d vector below
    std::vector<uint8_t> snps(nsnps*nsnps, false);  // true if snps is used
    short int gt;
    for (size_t k = 0; k < nsnps; ++k){
        for (size_t j = 0; j < k; ++j){
            bp_dist = genotypes[k].position - genotypes[j].position;
            if(bp_dist < 10)
                continue;

            // see NOTE XOR below
            gt = genotypes[k].genotype ^ genotypes[j].genotype;
            new_score = ((gt == 0) | (gt == 3))?
                match_bonus + bp_dist :
                mismatch_penalty;

            append_score = scores[j] + new_score;

            if (scores[k] < append_score){
                scores[k] = append_score;
                std::copy(snps.begin() + j*nsnps,
                        snps.begin() + (j+1)*nsnps,
                        snps.begin() + k*nsnps);
                snps[k*nsnps + k] = true;
            }
            if (scores[k] < new_score){
                scores[k] = new_score;
                std::fill(snps.begin() + k*nsnps,
                        snps.begin() + (k+1) * nsnps,
                        false);
                snps[k*nsnps + j] = true;
                snps[k*nsnps + k] = true;
            }
        }
    }
    auto maxScore = std::max_element(scores.begin(), scores.end());

    // update genotypes based on maxScore
    std::vector<WindowGT> gts;
    gts.reserve(nsnps);
    auto gt_input = genotypes.begin();
    // copy over snps which are used in max score
    auto maxIndex = maxScore - scores.begin();
    auto maxSnp = snps.begin() + maxIndex * nsnps;
    for(int i = 0; i < nsnps; ++i){
        if(*maxSnp)
            gts.push_back(*gt_input);
        ++gt_input;
        ++maxSnp;
    }

    genotypes.assign(gts.begin(), gts.end());
    return *maxScore;
}

// NOTE XOR
// valid genotypes are 1, 2, 3 at this point
// bitwise xor maps as follows:
// k  j  gt
// 1  1  0
// 1  2  3
// 1  3  2
// 2  2  0
// 2  3  1
// 3  3  0
// at this point gt is equal to 0 if j, k match or 3 if 
// j/k = 1/2
