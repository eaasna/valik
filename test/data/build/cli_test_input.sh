#!/bin/bash
cd build
set -Eeuo pipefail

# raptor_data_simulation has to be built from source
# git@github.com:eaasna/raptor_data_simulation.git
LENGTH=8192 # 8*2^10 = 8KiB
SEED=${1}
BIN_NUMBER=${2}
HAPLOTYPE_COUNT=${3}    # each haplotype is 1024bp

bin_length=$((LENGTH / BIN_NUMBER))
echo "Simulating $BIN_NUMBER bins with reference length of $LENGTH and bin_length of $bin_length"

# Simulate reference
echo "Simulating genome"
mason_genome -l $LENGTH -o ref.fasta -s $SEED &>/dev/null

# Evenly distribute it over bins
echo "Splitting genome into bins"
split_sequence --input ref.fasta --length $bin_length --parts $BIN_NUMBER
# We do not need the reference anymore
rm ref.fasta
# Rename the bins to .fa
for i in *.fasta; do mv $i $(basename $i .fasta).fa; done
# Simulate haplotypes for each bin
echo "Generating haplotypes"
for i in *.fa
do
   mason_variator \
       -ir $i \
       -n $HAPLOTYPE_COUNT \
       -of $(basename $i .fa).fasta \
       -ov $(basename $i .fa).vcf &>/dev/null
   rm $i
   rm $i.fai
   rm $(basename $i .fa).vcf
done

# Create metadata list of bin file names
seq -f "bin_%0g.fasta" 0 $((BIN_NUMBER-1)) > bin_paths.txt

valik split bin_paths.txt --metagenome --out bin_meta.bin
