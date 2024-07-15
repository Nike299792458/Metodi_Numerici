using Combinatorics, LinearAlgebra, Random

#= Cosa trovo in questo file:
Osservabili con discretizzzazione sbagliata, heathbath, microcanonico, JackKnife, Blocking
notazione: stvol= spacetime volume, STDIM spacetime dimensionality
le funzioni contrassegnate come _v sono una prova di implementazione vettoriale, ma non hanno portato a risulati compatibili 
=#
STDIM = 2

#O1 = 1/(Nt*Ns^{STDIM-1}) Σ_n (mhat^2 ϕ_n^2)
function O1d(stvol::Int, mhat::Float64, lattice::Array{Float64})
   ris= mhat*mhat*dot(lattice,lattice)/stvol
    return ris
end

#O2 = 1/(4*Nt*Ns^{STDIM-1}) Σ_n Σ_{mu>0} (ϕ(n+μ)-ϕ_(n-μ))^2
function O2d_v(stvol::Int, Nt::Int, lattice::Array{Float64})
    lattice=reshape(lattice, :)
    ris= dot(circshift(lattice,-Nt)-circshift(lattice,+Nt),circshift(lattice,-Nt)-circshift(lattice,+Nt))/(stvol*4)
    return ris
end

function O2d(stvol::Int, Nt::Int ,Ns::Int, lattice::Array{Float64})
    ris=0
    for r in LinearIndices(lattice)
        i=mod1(r,Nt)
        j=cld(r,Nt)
        aux=(lattice[i, mod1(j+1, Ns)] - lattice[i,mod1(j-1, Ns)])/2
        ris = ris + aux*aux
    end
    ris=ris/stvol
    return ris
end

#// O3 = 1/(4*Nt*Ns^{STDIM-1}) Σ_n (ϕ_(n+0)-ϕ_(n-0))^2
function O3d_v(stvol::Int, lattice::Array{Float64})
    rows = [collect(row) for row in eachcol(lattice)]
    a=0
    for i in 1:size(lattice,1)
    a+=dot(circshift(rows[i],-1)-circshift(rows[i],+1),circshift(rows[i],-1)-circshift(rows[i],+1))
    end
    ris= a/(stvol*4)
    return ris
end

function O3d(stvol::Int, Nt::Int, lattice::Array{Float64})
    ris=0
    for r in LinearIndices(lattice)
        i=mod1(r,Nt)
        j=cld(r,Nt)
        aux=(lattice[mod1(i+1, Nt), j]-lattice[mod1(i-1, Nt), j])/2
        ris = ris + aux*aux
    end
    ris=ris/stvol
    return ris
end
#=
Ns è il numero di colonne (lunghezza riga ) e Nt il numero di righe (lunghezza colonne)
per passare dal singolo indice r ai due indici (i,j):
j= celing-division tra 1 e Nt e i=modulo in base Nt di r
=#

function Σ_nd(lattice::Array{Float64}, r::Int, Nt::Int, Ns::Int)
    i=mod1(r,Nt)
    j=cld(r,Nt)
    lattice[mod1(i+2, Nt), j] + lattice[mod1(i-2, Nt), j] +
                   lattice[i, mod1(j+2, Ns)] + lattice[i, mod1(j-2, Ns)]
end
    

function heathbathd!(lattice::Array{Float64}, r::Int, mhat::Float64, Nt::Int, Ns::Int)
    std = 1.0/sqrt(mhat*mhat+0.5*STDIM)
    avg = Σ_nd(lattice, r, Nt, Ns)/(mhat*mhat+0.5*STDIM)/4
    lattice[r] = avg+std*randn()

    return 1
end

function overrelaxd!(lattice::Array{Float64}, r::Int, mhat::Float64, Nt::Int, Ns::Int)
    avg =  Σ_nd(lattice, r, Nt, Ns)/(mhat*mhat+0.5*STDIM)/4
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
