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
url_base="http://tigress-web.princeton.edu/~tcomi/sstar_tests"
popfile="${url_base}/base.popfile"

# eval echo expands the variable {1..10}
indivfile=$(eval echo "${url_base}/vcfs/${indiv}.mod.vcf.gz")
beds=(space_100 space_1000 rand_0 rand_1 rand_2)

read_result() {
    # cut removes deprecated calls
    # awk changes start and end to correct values
    wget -qO - $1 |
        zcat | cut -f1-11,42-47,49- |
        awk 'NR < 10 {print $0; next} !/^chrom/ {print $0}' |
        awk 'BEGIN {OFS="\t"}
    NR != 1 {
        if ( $11 == "." )
            {$16 = 0; $17 = 0}
        else
            {n=split($11, a, ","); $16=a[1]; $17=a[n]}}
    {print $0}'
}

run_sstar() {
    $sstar \
        --vcf <(wget -qO - $indivfile | zcat \
        | awk 'NR < 10 {print $0; next} !/^#/{print $0}' \
        ) \
        --popfile <(wget -qO - $popfile) \
        -l 50000 -s 10000 $1
}

# called on include only beds
run_sstar_i() {
    $sstar \
        --vcf <(wget -qO - $indivfile | zcat \
        | awk 'NR < 10 {print $0; next} !/^#/{print $0}' \
        ) \
        --popfile <(wget -qO - $popfile) \
        -l 50000 -s 10000 \
        --targets EUR ASN --references AFR \
        --include-bed <(wget -qO - $1 | zcat)
}

# called on exclude only beds
run_sstar_e() {
    $sstar \
        --vcf <(wget -qO - $indivfile | zcat \
        | awk 'NR < 10 {print $0; next} !/^#/{print $0}' \
        ) \
        --popfile <(wget -qO - $popfile) \
        -l 50000 -s 10000 \
        --targets EUR ASN --references AFR \
        --exclude-bed <(wget -qO - $1 | zcat)
}

# called on both beds
run_sstar_ie() {
    $sstar \
        --vcf <(wget -qO - $indivfile | zcat \
        | awk 'NR < 10 {print $0; next} !/^#/{print $0}' \
        ) \
        --popfile <(wget -qO - $popfile) \
        -l 50000 -s 10000 \
        --targets EUR ASN --references AFR \
        --include-bed <(wget -qO - $1 | zcat) \
        --exclude-bed <(wget -qO - $2 | zcat)
}

if [[ $test_type == "no_exclude" ]]; then
    echo "$indiv no exclude"
    resultfile=$(eval echo "${url_base}/outputs/${indiv}.windowcalc.gz")

    cmp \
        <(read_result "$resultfile") \
        <(run_sstar "--targets EUR ASN --references AFR")

elif [[ $test_type == "exclude" ]]; then
    echo "$indiv exclude"
    resultfile=$(eval echo "${url_base}/outputs_exclude/${indiv}.windowcalc.gz")

    cmp \
        <(read_result "$resultfile") \
        <(run_sstar "--targets EUR --excluded ASN --references AFR")

elif [[ $test_type == "bed_single" ]]; then
    for b in ${beds[@]}; do
        echo $indiv $b
        bedfile="${url_base}/bed_files/${b}.bed.gz"

        # include
        resultfile=$(eval echo "${url_base}/bed_outputs/${indiv}_${b}_.windowcalc.gz")
        cmp \
            <(read_result "$resultfile") \
            <(run_sstar_i $bedfile)

        # exclude
        resultfile=$(eval echo  "${url_base}/bed_outputs/${indiv}__${b}.windowcalc.gz")
        cmp \
            <(read_result "$resultfile") \
            <(run_sstar_e $bedfile)
    done

elif [[ $test_type == "bed_pair" ]]; then
    for ind in ${!beds[@]}; do
        for ind2 in $(seq $(($ind + 1)) $((${#beds[@]} - 1))); do
            # bed names
            b1=${beds[$ind]}
            b2=${beds[$ind2]}
            echo "$indiv $b1 | $b2"

            # urls
            bedfile1="${url_base}/bed_files/${beds[$ind]}.bed.gz"
            bedfile2="${url_base}/bed_files/${beds[$ind2]}.bed.gz"

            resultfile=$(eval echo "${url_base}/bed_outputs/${indiv}_${b1}_${b2}.windowcalc.gz")
            # include 1
            cmp \
                <(read_result "$resultfile") \
                <(run_sstar_ie $bedfile1 $bedfile2)

            # exclude 1
            resultfile=$(eval echo  "${url_base}/bed_outputs/${indiv}_${b2}_${b1}.windowcalc.gz")

            cmp \
                <(read_result "$resultfile") \
                <(run_sstar_ie $bedfile2 $bedfile1)

        done
    done

else
    echo "Unknown test type $test_type"
    exit 1
fi
