#!/bin/bash

seeds=`cat ../seeds_div4.txt`

rm -r previous_div4; mkdir previous_div4
mv *.eps *.jld *.pdf *.log previous_div4

for seed in ${seeds[*]}
do
    echo "Performing with seed $seed..."
    ./perform_one.jl --seed $seed
done
