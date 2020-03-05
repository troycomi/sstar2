#!/bin/bash

set -euo pipefail

indiv=$1
test_type=$2
sstar="src/sstar2"
# to work for local tests
if [[ ! -f $sstar ]]; then
    sstar="../../../build/src/sstar2"
fi

if [[ $indiv == "all" ]]; then
    indiv={1..10}
fi
popfile="http://tigress-web.princeton.edu/~tcomi/sstar_tests/base.popfile"
# eval echo expands the variable {1..10}
indivfile=$(eval echo \
    "http://tigress-web.princeton.edu/~tcomi/sstar_tests/vcfs/${indiv}.mod.vcf.gz")
beds=(space_100 space_1000 rand_0 rand_1 rand_2)

if [[ $test_type == "no_exclude" ]]; then
    resultfile=$(eval echo \
        "http://tigress-web.princeton.edu/~tcomi/sstar_tests/outputs/${indiv}.windowcalc.gz")

    # cut removes deprecated calls
    # awk changes start and end to correct values
    cmp \
        <(wget -qO - $resultfile | \
            zcat | cut -f1-11,42-47,49- | \
            awk 'NR < 10 {print $0; next} !/^chrom/ {print $0}' |
            awk 'BEGIN {OFS="\t"} NR != 1 {if ( $11 == "." ){$16 = 0; $17 = 0} else {n=split($11, a, ","); $16=a[1]; $17=a[n]}} { print $0}') \
        <(\
        $sstar \
            --vcf <(wget -qO - $indivfile | zcat \
            | awk 'NR < 10 {print $0; next} !/^#/{print $0}' \
            ) \
            --popfile <(wget -qO - $popfile) \
            --targets EUR ASN \
            --references AFR \
            -l 50000 -s 10000 )

elif [[ $test_type == "exclude" ]]; then
    resultfile=$(eval echo \
        "http://tigress-web.princeton.edu/~tcomi/sstar_tests/outputs_exclude/${indiv}.windowcalc.gz")

    # cut removes deprecated calls
    # awk changes start and end to correct values
    cmp \
        <(wget -qO - $resultfile | \
            zcat | cut -f1-11,42-47,49- | \
            awk 'NR < 10 {print $0; next} !/^chrom/ {print $0}' |
            awk 'BEGIN {OFS="\t"} NR != 1 {if ( $11 == "." ){$16 = 0; $17 = 0} else {n=split($11, a, ","); $16=a[1]; $17=a[n]}} { print $0}') \
        <(\
        $sstar \
            --vcf <(wget -qO - $indivfile | zcat \
            | awk 'NR < 10 {print $0; next} !/^#/{print $0}' \
            ) \
            --popfile <(wget -qO - $popfile) \
            --targets EUR \
            --excluded ASN \
            --references AFR \
            -l 50000 -s 10000 )

elif [[ $test_type == "bed_single" ]]; then
    # cut removes deprecated calls
    # awk changes start and end to correct values
    for b in ${beds[@]}; do
        bedfile="http://tigress-web.princeton.edu/~tcomi/sstar_tests/bed_files/${b}.bed.gz"
        # include
        resultfile=$(eval echo \
            "http://tigress-web.princeton.edu/~tcomi/sstar_tests/bed_outputs/${indiv}_${b}_.windowcalc.gz")

        cmp \
            <(wget -qO - $resultfile | \
                zcat | cut -f1-11,42-47,49- | \
                awk 'NR < 10 {print $0; next} !/^chrom/ {print $0}' |
                awk 'BEGIN {OFS="\t"} NR != 1 {if ( $11 == "." ){$16 = 0; $17 = 0} else {n=split($11, a, ","); $16=a[1]; $17=a[n]}} { print $0}') \
            <(\
            $sstar \
                --vcf <(wget -qO - $indivfile | zcat \
                | awk 'NR < 10 {print $0; next} !/^#/{print $0}' \
                ) \
                --popfile <(wget -qO - $popfile) \
                --targets EUR ASN\
                --references AFR \
                --include-bed <(wget -qO - $bedfile | zcat) \
                -l 50000 -s 10000 )

        # exclude
        resultfile=$(eval echo \
            "http://tigress-web.princeton.edu/~tcomi/sstar_tests/bed_outputs/${indiv}__${b}.windowcalc.gz")

        cmp \
            <(wget -qO - $resultfile | \
                zcat | cut -f1-11,42-47,49- | \
                awk 'NR < 10 {print $0; next} !/^chrom/ {print $0}' |
                awk 'BEGIN {OFS="\t"} NR != 1 {if ( $11 == "." ){$16 = 0; $17 = 0} else {n=split($11, a, ","); $16=a[1]; $17=a[n]}} { print $0}') \
            <(\
            $sstar \
                --vcf <(wget -qO - $indivfile | zcat \
                | awk 'NR < 10 {print $0; next} !/^#/{print $0}' \
                ) \
                --popfile <(wget -qO - $popfile) \
                --targets EUR ASN \
                --references AFR \
                --exclude-bed <(wget -qO - $bedfile | zcat) \
                -l 50000 -s 10000 )
    done

else
    echo "Unknown test type $test_type"
fi
