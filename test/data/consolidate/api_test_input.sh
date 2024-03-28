#!/bin/bash
cd consolidate
set -Eeuo pipefail

rm -f *_full.gff
rm -f *bins*overlap_segment_metadata.bin
rm -f *bins*overlap_reference_metadata.bin

ref_file=multi_seq_ref.fasta
read_length=150
read_dir=reads_$read_length
mkdir -p $read_dir
query_file=multi_seq_query.fasta
read_count=2
error_rate=0.06

i=1
# make multiple copies of identical chromosomes so that each read has multiple matches across the genome
for length in 560 560 560 464 464 464 464 210 
do
    echo "Simulating chromosome with length $length"
    chr_out="chr"$i".fasta"
    if (( $length == 560 )) ; then
        SEED=102
    elif (( $length == 464 )) ; then
        SEED=19
    else
        SEED=42
    fi
    mason_genome -l $length -o $chr_out -s $SEED &>/dev/null

    sed -i "s/^>.*$/>chr$i/g" $chr_out
    
    #----------- Sample reads from reference sequence -----------
    echo "Generating $read_count reads of length $read_length with error rate $error_rate"
    generate_local_matches \
        --matches-out $read_dir/chr$i.fasta \
        --genome-out $read_dir/genome.fasta \
        --max-error-rate $error_rate \
        --num-matches $read_count \
        --min-match-length $read_length \
        --max-match-length $read_length \
        --reverse \
        --ref-len $length \
        --seed $SEED \
        $chr_out
    
    sed -i "s/>.*/&_$i/" $read_dir/chr$i.fasta
    let i=i+1
done

# duplicate the sequence of the last chromosome so that reads have multiple matches per chromosome
cat chr$(expr $i - 1).fasta > chr${i}.fasta
tail -n +2 chr$(expr $i - 1).fasta >> chr${i}.fasta
tail -n +2 chr$(expr $i - 1).fasta >> chr${i}.fasta
rm chr$(expr $i - 1).fasta

cat chr*.fasta > $ref_file
rm chr*.fasta
cat $read_dir/chr*.fasta > query_tmp.fasta
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < query_tmp.fasta > $query_file
rm -r $read_dir
rm query_tmp.fasta

min_len=50

for bin in 8 16
do
        valik split $ref_file --out ${bin}bins${min_len}overlap_reference_metadata.bin --seg-count $bin --pattern $min_len --without-parameter-tuning --write-out        
        cat multi_seq_ref_*.fasta > multi_seq_ref.segments.fasta
        rm multi_seq_ref_*.fasta
        rm seg_files.txt
        grep ">" multi_seq_ref.segments.fasta | cut -c 2- | awk -F'_' '{print $1 "\t" $2 "\t" $3}' > segments.tsv
        while read -r id start len;
        do
                end=$(echo $start + $len | bc)
                stellar -e $error_rate -l $min_len -o multi_seq_ref_${id}_${start}_${len}.gff \
                        --sequenceOfInterest $id --segmentBegin $start \
                        --segmentEnd $end $ref_file $query_file > /dev/null

        done < segments.tsv

        rm segments.tsv multi_seq_ref.segments.fasta

        cat multi_seq_ref_*.gff > ${bin}bins${min_len}overlap_dream_all.gff
        rm multi_seq_ref_*
done
