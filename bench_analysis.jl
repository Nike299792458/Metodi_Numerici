using ArgParse, DelimitedFiles, LinearAlgebra, Printf, Statistics
include("free_scalar.jl")

Nt=2
Ns=2
ratio=2
stvol=Nt*Ns
obs1=[]
therm=10000
blocksize=2000


path = joinpath(["..", "simulations_c"])
fname="free_scalar_th_sample=1.0e+07ratio=2Nt=02Ns=02.txt"
fr = joinpath([path, fname])
w = open(fr, "r") do io
    readdlm(io, header = true)
end
O1_j = JackKnife(w[1][therm:end,1], blocksize)
push!(obs1, mean(O1_j ))
println(obs1)