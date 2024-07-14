#!/bin/bash

#questa va in spectrum.log
path="../simulations_c/spectrum_c"
sample=10000000

Nt=(60)
for Nti in "${Nt[@]}"
do 
	julia spectrum_main.jl   "$Nti" "$sample" "-v"
done
