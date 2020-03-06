![](https://github.com/troycomi/sstar2/workflows/UnitTests/badge.svg)
![](https://github.com/troycomi/sstar2/workflows/AcceptanceTests/badge.svg)

# sstar2
C++ implementation of sstar analysis

## About
This project aims to improve [freezing-archer](https://github.com/bvernot/freezing-archer)
for sstar analyses.  Since that repository contains many modules and is very
flexible, sstar2 will not fully replace all features and currently provides
equivalent output for a subset of situations.  Additional uses will be included
as they are requested; please submit an issue to request support!

sstar2 depends on c++ 11 standard libraries and [CLI11](https://github.com/CLIUtils/CLI11/),
which is included in this repository.  sstar2 is **leaner and faster**, around 40x
faster on simulated data.  Currently, sstar2 handles the case of `--no-pvalues`
without masking.  Input vcf files may contain multiple chromosomes as long as
they are sorted.  The entire region is processed.  Output column names match
freezing-archer, but unsupported and deprecated columns related to p values are
removed.

## Installation
The executable is produced using
[cmake>=3.14](https://cliutils.gitlab.io/modern-cmake/chapters/intro/installing.html).
Download or clone the repository.  From the root directory,
run the following to make the sstar2 executable:
```bash
mkdir build
cd build
cmake ..
cmake --build .
mv src/sstar2 ${FINAL_LOCATION}
```
Where the last step moves the standalone executable to wherever you want to
run sstar2 from and is optional.

To run unit tests, perform the following from the build directory:
```bash
cmake .. -DBUILD_TESTING=True
cmake --build .
ctest
```

## Usage
The following options can be retrieved with `sstar2 --help`:
```bash
Options:
-h,--help                   Print this help message and exit
-v,--vcf TEXT:FILE REQUIRED Input vcf file; can accept input redirection
-p,--popfile TEXT:FILE REQUIRED
Population file; tsv with indiv, pop, superpop
-t,--targets TEXT ... REQUIRED
Comma separated list of target populations or individuals
-r,--references TEXT ... REQUIRED
Comma separated list of reference populations or individuals
-e,--excluded TEXT ...      Comma separated list of excluded populations or individuals;
                            default to none excluded
-l,--length UINT            Window length; default 50,000
-s,--step UINT              Window step; default 10,000
--match-bonus INT           Match bonus for sstar; default 5000
--mismatch-penalty INT      Mismatch penalty for sstar; default -10000
--include-bed TEXT:FILE     Bed file with regions to include
--exclude-bed TEXT:FILE     Bed file with regions to exclude
-o,--output TEXT            Output file; can accept input redirection; default stdout
```

To convert from freezing-archer:
```bash
-vcf file.vcf                -> --vcf file.vcf.gz
-vcfz file.vcf.gz            -> --vcf <(zcat file.vcf.gz)
-ref-pops AFR -ref-inds ind1 -> --references AFR,ind1
-winlen 50000                -> --length 50000
-winstep 10000               -> --step 10000
--regions myfile.bbg         -> --include-bed myfile.bed
--exclude-region myfile.bbg  -> --exclude-bed myfile.bed
```
Other options are similar.  Note that only a single bed file is accepted for
include or exclude (though both may be specified in a run).  If multiple bed
files need to be merged that has to be done in a separate step.  However,
sstar2 accepts standards bed files instead of binary bed files.

```bash
./sstar2 \
    --vcf <(zcat 1.mod.vcf.gz) \
    --popfile base.popfile \
    --targets EUR,ASN \
    --references AFR \
    -l 50000 \
    -s 10000 \
    | head -n 5

chrom   winstart   winend  n_snps  n_ind_snps   n_region_ind_snps       ind_id    pop     s_star  num_s_star_snps   s_star_snps             hap_1_s_start   hap_1_s_end     hap_2_s_start   hap_2_s_end     s_start s_end   n_s_star_snps_hap1    n_s_star_snps_hap2       s_star_haps     callable_bases
1       0          50000   451     4            194                     msp_110   EUR     45603   4                 5382,28610,32662,35985  0               0               5382            32662           5382    35985   1                     3                        2,2,2,1         50000
1       0          50000   451     0            190                     msp_111   EUR     0       0                 .                       0               0               0               0               0       0       0                     0                        .               50000
1       0          50000   451     3            193                     msp_112   EUR     28616   3                 17506,26875,36122       0               0               17506           36122           17506   36122   1                     2                        2,1,2           50000
1       0          50000   451     2            192                     msp_113   EUR     0       0                 .                       0               0               0               0               0       0       0                     0                        .               50000
```
