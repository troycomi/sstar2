#!/bin/bash

set -euo pipefail

indiv=$1
test_type=$2

if [[ $indiv == "all" ]]; then
    indiv={1..10}
fi
actual="test.windowcalc.gz"
popfile="http://tigress-web.princeton.edu/~tcomi/sstar_tests/base.popfile"
# eval echo expands the variable {1..10}
indivfile=$(eval echo \
    "http://tigress-web.princeton.edu/~tcomi/sstar_tests/vcfs/${indiv}.mod.vcf.gz")
resultfile=$(eval echo \
    "http://tigress-web.princeton.edu/~tcomi/sstar_tests/outputs/${indiv}.windowcalc.gz")

if [[ $test_type == "no_exclude" ]]; then

    # cut removes deprecated calls
    # awk changes start and end to correct values
    cmp \
        <(wget -qO - $resultfile | \
            zcat | cut -f1-11,42-47,49- | \
            awk 'NR < 10 {print $0; next} !/^chrom/ {print $0}' |
            awk 'BEGIN {OFS="\t"} NR != 1 {if ( $11 == "." ){$16 = 0; $17 = 0} else {n=split($11, a, ","); $16=a[1]; $17=a[n]}} { print $0}') \
        <(\
        src/sstar2 \
            --vcf <(wget -qO - $indivfile | zcat \
            | awk 'NR < 10 {print $0; next} !/^#/{print $0}' \
            ) \
            --popfile <(wget -qO - $popfile) \
            --targets EUR ASN \
            --references AFR \
            -l 50000 -s 10000 )

else
    echo "Unknown test type $test_type"
fi

# expected=/tigress/tcomi/abwolf_abc/SplitPop/testing/RegionFiles/{1..10}.windowcalc.gz
# actual=/tigress/tcomi/abwolf_abc/SplitPop/testing/RegionFiles/test_merge.windowcalc.gz

# # cut removes deprecated calls
# # awk changes start and end to correct values
# # line=$(cmp <(zcat $expected | cut -f1-11,42-47,49- | awk 'BEGIN {OFS="\t"} NR != 1 {if ( $11 == "." ){$16 = 0; $17 = 0} else {n=split($11, a, ","); $16=a[1]; $17=a[n]}} { print $0}') $actual | awk '{print $NF}')
# line=$(cmp <(zcat /tigress/tcomi/abwolf_abc/SplitPop/testing/RegionFiles/{1..10}.windowcalc.gz | \
    #     cut -f1-11,42-47,49- | \
    #     awk 'NR < 10 {print $0; next} !/^chrom/ {print $0}' |
#     awk 'BEGIN {OFS="\t"} NR != 1 {if ( $11 == "." ){$16 = 0; $17 = 0} else {n=split($11, a, ","); $16=a[1]; $17=a[n]}} { print $0}') <(zcat $actual) | awk '{print $NF}')
# if [[ ! -z $line ]]; then
#     echo "MISMATCH! Line: $line"
#     awk -v line=$line 'NR==line{print "Actual:   "$0; exit}' <(zcat $actual)
#     exit
# fi
# exit
