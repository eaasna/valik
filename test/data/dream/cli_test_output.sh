#!/bin/bash
cd dream
set -Eeuo pipefail

if [ -z "${VALIK_TMP}" ]; then
    echo "no VALIK_TMP folder given"
    exit 127
fi

mkdir -p $VALIK_TMP

#----------- Index and search the reference genome -----------

# Split parameters
seg_overlap="150"     # how much adjacent segments overlap

# Build parameters
k=13
ibf_size="32k"

# Search parameters
pattern=50            # min local match length
pat_overlap=49        # how much adjacent patterns overlap

ref_input="ref.fasta"
query="query.fasta"

e=1
er=$(echo $e/$pattern | bc -l)
for b in 4 16
do
    echo "Splitting the genome into $b segments that overlap by $seg_overlap"
    seg_meta="seg_meta"$seg_overlap"overlap"$b"bins.bin"
    valik split "$ref_input" --pattern "$seg_overlap" --seg-count "$b" --out "$seg_meta" --without-parameter-tuning

    for w in 13 15
    do
        echo "Creating IBF for w=$w and k=$k where segments overlap by $seg_overlap"
        index=$b"bins"$w"window.ibf"
        valik build --kmer "$k" --window "$w" --size "$ibf_size" --output "$index" --ref-meta "$seg_meta"

        echo "Searching IBF with $e errors"
        dist_out=$b"bins"$w"window"$e"error.gff"
        #local_out="local"$b"bins"$w"window"$e"error.gff"
        valik search --distribute --index "$index" --query "$query" --output "$dist_out" --error-rate "$er" --pattern "$pattern" --overlap "$pat_overlap" \
                     --ref-meta "$seg_meta" --repeatPeriod 1 --repeatLength 10 --numMatches 2
        #valik search --index "$index" --query "$query" --output "$local_out" --error-rate "$er" --pattern "$pattern" --overlap "$pat_overlap" \
        #             --ref-meta "$seg_meta" --repeatPeriod 1 --repeatLength 10 --numMatches 1

        rm $VALIK_TMP/*
	rm $index
	done
    rm "$seg_meta"
done

stellar_out="stellar.gff"
stellar $ref_input $query -e $er -l $pattern -o $stellar_out --repeatPeriod 1 --repeatLength 10 --numMatches 2 &> /dev/null
rm stellar.disabled.fasta

rm -r $VALIK_TMP
