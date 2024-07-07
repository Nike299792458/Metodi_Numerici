#!/bin/bash

path="../simulations_c/"
sample=5000000
# Array to loop over
Nt=(4 6 8 7 10)
ratio=4


for Nti in "${Nt[@]}"
do 
	julia th_free_scalar.jl  "$ratio" "$Nti" "$sample"
done

