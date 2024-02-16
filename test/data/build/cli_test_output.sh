#!/bin/bash
cd build
set -Eeuo pipefail

for w in 19 23
do
    echo "Creating IBF for w=$w and k=19"
    output="8bins"$w"window.ibf"
    valik build --window "$w" --kmer 19 --output "$output" --size 32k --ref-meta bin_meta.bin
done

rm bin_meta.bin
