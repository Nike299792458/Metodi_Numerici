using LinearAlgebra, Random

STDIM = 2
function ϕ_fixed_p(lattice, p, Nt, Ns)
    stvol = Nt * Ns^(STDIM-1)
    O = Complex{Float64}[0.0 + 0.0im for _ in 1:Nt]

    for r in 1:stvol
        i=mod1(r,Nt)
        j=cld(r,Nt)
        sum = 0.0
        for j in 1:Ns
            for k in 1:(STDIM-1)
             sum += p[k]* j
            end
        end
        O[i] += lattice[r] * exp(-1im * sum)  
    end

    for i in 1:Nt
        O[i] /= sqrt(stvol)
    end
    return O
end

function t_corr(lattice, Nt, δt)
    p = zeros(Float64, STDIM-1)
    corr_p0 = zeros(Float64, δt )

    O = ϕ_fixed_p(lattice, p, Nt, δt)

    corr_p0= dot(O, circshift(O , δt))
    return corr_p0/Nt
end

function Σ_n(lattice::Array{Float64}, r::Int, Nt::Int, Ns::Int)
    i=mod1(r,Nt)
    j=cld(r,Nt)
    lattice[mod1(i+1, Nt), j] + lattice[mod1(i-1, Nt), j] + lattice[i, mod1(j+1, Ns)] + lattice[i, mod1(j-1, Ns)]
end
    

function heathbath!(lattice::Array{Float64}, r::Int, mhat::Float64, Nt::Int, Ns::Int)
    std = 1.0/sqrt(mhat*mhat+2.0*STDIM)
    avg = Σ_n(lattice, r, Nt, Ns)/(mhat*mhat+2.0*STDIM)
    lattice[r] = avg+std*randn()

    return 1
end 

function overrelax!(lattice::Array{Float64}, r::Int, mhat::Float64, Nt::Int, Ns::Int)
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
