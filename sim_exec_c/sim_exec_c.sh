#!/bin/bash

path="../simulations_c/"
sample=100000
# Array to loop over
ratio=(4 6 8)

for ratioi in "${ratio[@]}"
do 
	julia th_free_scalar_main.jl "$ratioi" "$sample"
done


