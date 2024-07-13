#!/bin/bash

path="../simulations_c/"
sample=5000000
# Array to loop over
Nt=( 5 6 8 10)
ratio=5


for Nti in "${Nt[@]}"
do 
	julia doublers_main.jl  "$ratio" "$Nti" "$sample"
done