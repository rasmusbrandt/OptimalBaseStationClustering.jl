#!/bin/bash

seeds=`cat ../seeds.txt`

rm -r previous; mkdir previous
mv *.eps *.jld *.pdf *.log previous

for seed in ${seeds[*]}
do
    echo "Performing with seed $seed..."
    ./perform_one.jl --seed $seed
done

echo "Merging..."
./merge.jl

echo "Plotting..."
./plot.jl
