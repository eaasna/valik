#!/bin/bash
cd split
set -Eeuo pipefail

#----------- Simulate chromosomes of various lengths -----------

SEED=${1}

i=1
for length in 210 490 420 280 4
do
    echo "Simulating chromosome with length $length"
    chr_out="chr"$i".fasta"
    mason_genome -l $length -o $chr_out -s $SEED &> /dev/null

    sed -i "s/^>.*$/>chr$i/g" $chr_out
    let i=i+1
done

cat chr*.fasta > database.fasta
rm chr*.fasta

for pattern in 20
do
    for bins in 4 16
    do
        out_dir=write_out_${pattern}_${bins}
        mkdir -p ${out_dir}

        valik split database.fasta --pattern ${pattern} --seg-count ${bins} --out ${out_dir}/reference_metadata.bin --without-parameter-tuning
    done
done
