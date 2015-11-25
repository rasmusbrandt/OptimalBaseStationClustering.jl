#!/bin/bash

seeds=`cat ../seeds_div3.txt`

rm -r previous_div3; mkdir previous_div3
mv *.eps *.jld *.pdf *.log previous_div3

for seed in ${seeds[*]}
do
    echo "Performing with seed $seed..."
    ./perform_one.jl --seed $seed
done
