#!/bin/bash

seeds=`cat ../seeds_div2.txt`

rm -r previous_div2; mkdir previous_div2
mv *.eps *.jld *.pdf *.log previous_div2

for seed in ${seeds[*]}
do
    echo "Performing with seed $seed..."
    ./perform_one.jl --seed $seed
done
