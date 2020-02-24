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
The executable is produced using cmake.  Download or clone the repository.
From the root directory, run the following to make the sstar2 executable:
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
-o,--output TEXT            Output file; can accept input redirection; default stdout
```

To convert from freezing-archer:
```bash
-vcfz file.vcf               -> --vcf file.vcf.gz
-vcfz file.vcf.gz            -> --vcf <(zcat file.vcf.gz)
-ref-pops AFR -ref-inds ind1 -> --references AFR,ind1
-winlen 50000                -> --length 50000
-winstep 10000               -> --step 10000
```
Other options are similar.
