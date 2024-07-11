#!/bin/bash

#questa va in spectrum.log
path="../simulations_c/spectrum_c"
sample=100000000

Nt=60

do 
	julia spectrum_analysis.jl   "$Nt" "$sample"
done
