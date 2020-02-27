#!/bin/bash -l

set -ex

mkdir -p cmake_dir
mkdir -p build_tmp
mkdir -p cmake_sources
rm -rf cmake_dir/* build_tmp/*

test_type=$1
chrom_type=$2
v="3.16.2"
fn=cmake-$v-Linux-x86_64.tar.gz

if [ ! -f cmake_sources/$fn ]; then
    wget -qO cmake_sources/$fn "https://cmake.org/files/v${v%.*}/$fn"
fi

tar -xzf cmake_sources/$fn --strip-components=1 -C $PWD/cmake_dir

export PATH=$PWD/cmake_dir/bin:$PATH

cmake --version

cd build_tmp 
cmake ..
cmake --build . 

if [[ $chrom_type == "all" ]]; then
    /acceptance_test.sh all $test_type
else
    for i in {1..10}; do
        /acceptance_test.sh $i $test_type
    done
fi
