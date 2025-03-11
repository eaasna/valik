#!/bin/bash
cd split
set -Eeuo pipefail

if [ -z "${VALIK_TMP}" ]; then
    echo "no VALIK_TMP folder given"
    exit 127
fi

mkdir -p $VALIK_TMP

#----------- Split multiple sequences of various lengths -----------

seg_input="various_chromosome_lengths.fasta"
for p in 20
do
    for b in 4 16
    do
        echo "Splitting the genome into $b segments that overlap by $p"
        seg_meta="multi/"$p"overlap"$b"bins.bin"
        valik split "$seg_input" --pattern "$p" --seg-count "$b" --out "$seg_meta" --without-parameter-tuning
    done
done

#----------- Split and search a single reference sequence -----------

# Split parameters
seg_overlap="150"     # how much adjacent segments overlap

# Build parameters
k=13
ibf_size="32k"

# Search parameters
errors=1              # max allowed errors
pattern=50            # min local match length

ref_input="single_reference.fasta"
query="single_query.fasta"
for b in 4 16
do
    echo "Splitting the genome into $b segments that overlap by $seg_overlap"
    seg_meta="single/"$seg_overlap"overlap"$b"bins.bin"

    valik split "$ref_input" --pattern "$seg_overlap" --seg-count "$b" --out "$seg_meta" --without-parameter-tuning

    for w in 13 15
    do
        echo "Creating IBF for w=$w and k=$k where segments overlap by $seg_overlap"
        index="single/"$seg_overlap"overlap"$b"bins"$w"window.ibf"
        valik build --kmer "$k" --window "$w" --size "$ibf_size" --output "$index" --ref-meta "$seg_meta" --without-parameter-tuning

        echo "Searching IBF with $errors errors"
        search_out="single/"$seg_overlap"overlap"$b"bins"$w"window"$errors"errors.gff"
        error_rate=$(echo $errors/$pattern| bc -l )
        valik search --distribute --index "$index" --query "$query" --output "$search_out" --error-rate "$error_rate" \
                     --pattern "$pattern" --query-every 1 --ref-meta "$seg_meta" --threads 1 --very-verbose \
                     --without-parameter-tuning --cart-max-capacity 3 --max-queued-carts 10
    rm "$search_out"    # only look at .gff.out
    done
done

rm -r $VALIK_TMP stellar.disabled.fasta
