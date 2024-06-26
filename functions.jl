using Combinatorics, LinearAlgebra, Random

#= Cosa trovo in questo file:
Osservabili, heathbath, metropolis, microcanonico
notazione= stvol= spacetime volume
=#


function init_lattice(L::Int)
    lattice = 
    return lattice
end

#O1 = (1/stvol) Σ_n (hatm^2 ϕ_n^2)
function O1(stvol::Int, mhat::Float64, lattice::Array{Float64})
   ris= sum(mhat*mhat*(lattice .^2))
    return ris
end

#O2 = (1/stvol) Σ_r  Σ_{mu>0} (ϕ(n+μ)-ϕ_n)^2
function O2(stvol::Int, mhat::Int, Nt::Int, lattice::Array{Int})
    ris= mhat*mhat*dot(lattice,circshift(lattice,-Nt))/stvol
    return ris
end
#// O3 = (1/stvol) Σ_n (ϕ_(n+0)-ϕ_n)^2
function O3(stvol::Int, mhat::Int, lattice::Array{Int})
    ris= mhat*mhat*dot(lattice,circshift(lattice,-1))/stvol
    return ris
end

function metropolis!(lattice::Array{Float64}, r::Int, Δ::Int, η::Float64)
    ΔS = Δ*(2*rand(Float64)-1)
    trial = (lattice[r]+ΔS)
    Eold = 
    Enew = 

    if Enew < Eold
        lattice[r] = trial
        return 1
    elseif rand(Float64) < exp(-(Enew-Eold))
        lattice[r] = trial
        return 1
    end

    return 0
end

function heathbath!(lattice::Array{Float64}, r::Int, eta::Float64)
    std = 1.0/sqrt(eta + 2.0/eta)
    avg = (circshift(lattice, 1)[r]+circshift(lattice, -1)[r])/(eta*(eta + 2.0/eta))
    
    lattice[r] = avg+std*randn()

    return 1
end

function overrelax!(lattice::Array{Float64}, r::Int, eta::Float64)
    avg = (circshift(lattice, 1)[r]+circshift(lattice, -1)[r])/(eta*(eta + 2.0/eta))
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
