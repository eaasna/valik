#!/bin/bash
cd prepare
set -Eeuo pipefail

ref_input="../dream/ref.fasta"
p=50
for shape in "11111111" "1001"
do
    echo "Creating IBF for shape $shape"
    seg_meta="s${shape}.bin"
    index="s${shape}.ibf"
    mkdir -p "s$shape"

    if [ ! -s s$shape ]; then
        rm "s$shape/*"
    fi
    
    valik split "$ref_input" --shape $shape --pattern "$p" --out "$seg_meta" --write-out
    valik build --output "$index" --fast --ref-meta $seg_meta --kmer-count-min 1
    
    mv ref.*.header s$shape/
    mv ref.*.minimiser s$shape/
done

mkdir -p bins
rm bins/*
mv ../dream/ref_*.fasta bins/
