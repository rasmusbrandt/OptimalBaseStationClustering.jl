#!/bin/bash

seeds=`cat ../seeds_div1.txt`

rm -r previous_div1; mkdir previous_div1
mv *.eps *.jld *.pdf *.log previous_div1

for seed in ${seeds[*]}
do
    echo "Performing with seed $seed..."
    ./perform_one.jl --seed $seed
done
