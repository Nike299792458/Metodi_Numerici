using Combinatorics, LinearAlgebra, Random

#= Cosa trovo in questo file:
Osservabili, heathbath, metropolis, microcanonico
notazione= stvol= spacetime volume
=#
STDIM = 2 #spacetime dimensionality

#O1 = 1/(Nt*Ns^{STDIM-1}) Σ_n (hatm^2 ϕ_n^2)
function O1(stvol::Int, mhat::Int, lattice::Array{Float64})
   ris= sum(mhat*mhat*(lattice .^2))/stvol
    return ris
end

#O2 = 1/(Nt*Ns^{STDIM-1}) Σ_r  Σ_{mu>0} (ϕ(n+μ)-ϕ_n)^2
function O2(stvol::Int, mhat::Int, Nt::Int, lattice::Array{Float64})
    ris= mhat*mhat*dot(lattice,circshift(lattice,-Nt))/stvol
    return ris
end
#// O3 = 1/(Nt*Ns^{STDIM-1}) Σ_n (ϕ_(n+0)-ϕ_n)^2
function O3(stvol::Int, mhat::Int, lattice::Array{Float64})
    rows = [collect(row) for row in eachcol(lattice)]
    a=0
    for i = 1:size(lattice,1)
    a+=dot(rows[i], circshift(rows[i],-1))
    end
    ris= mhat*mhat*a/stvol
    return ris
end
#=
Ns è il numero di colonne (lunghezza riga ) e Nt il numero di righe (lunghezza colonne)
per passare dal singolo indice r ai due indici (i,j):
j= celing-division tra 1 e Nt e i=modulo in base Nt di r
=#

function Σ_n(lattice::Array{Float64}, r::Int, Nt::Int, Ns::Int)
    i=mod1(r,Nt)
    j=cld(r,Nt)
    lattice[mod1(i+1, Nt), j] + lattice[mod1(i-1, Nt), j] +
                   lattice[i, mod1(j+1, Ns)] + lattice[i, mod1(j-1, Ns)]
end
    

function heathbath!(lattice::Array{Float64}, r::Int, mhat::Int, Nt::Int, Ns::Int)
    std = 1.0/sqrt(mhat*mhat+2.0*STDIM)
    avg = Σ_n(lattice, r, Nt, Ns)/(mhat*mhat+2.0*STDIM)
    lattice[r] = avg+std*randn()

    return 1
end

function overrelax!(lattice::Array{Float64}, r::Int, mhat::Int, Nt::Int, Ns::Int)
    avg =  Σ_n(lattice, r, Nt, Ns)/(mhat*mhat+2.0*STDIM)
    new = 2.0*avg-lattice[r]
    lattice[r] = new

    return 1
end


function Blocking(x::Array, blocksize::Int)
    length = size(x,1) ÷ blocksize 
    x_cut = x[1:length*blocksize]
    x_b = reshape(x_cut, (blocksize, length))

    return mean(x_b, dims = 1), length
end

function JackKnife(x::Array, blocksize::Int)
    x_blocked, length = Blocking(x, blocksize)
    jk(a) = let s = sum(a); [s-v for v in a] end
    x_j = vec(jk(x_blocked)/(length-1))
    return x_j
end
