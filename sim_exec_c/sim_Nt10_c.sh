#!/bin/bash

path="../simulations_c/Nt=10"
sample=5000000
# Array to loop over
#questa va il Ntequals10.log
Nt=(10)
ratio=5


for Nti in "${Nt[@]}"
do 
	julia th_free_scalar.jl  "$ratio" "$Nti" "$sample"
done


