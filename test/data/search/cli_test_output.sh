#!/bin/bash

cd search

tau=1.0
p_max=1.0

for w in 19 23
do
    for e in 0 1
    do
        for p in 50 100
        do
            for o in 1 40
            do
                echo "Searching the IBF (w=$w, k=19) for an approximate local match of length $p with $e error(s)"
                echo "Potential matches overlap by $o bp"
                output="8bins"$w"window"$e"error"$p"pattern"$o"overlap.out"
                valik search --index ../build/8bins${w}window.ibf --query query.fq --output "$output" --error "$e" --pattern "$p" --overlap "$o" --tau "$tau" --p_max "$p_max"
	        done
        done
    done
done
