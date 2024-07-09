using LinearAlgebra
using Random
using Plots
#=
# Funzione per convertire indice lineare in coordinate cartesiane
function lex_to_cart_st(r, Nt, Ns, STDIM)
    coord = zeros(Int, STDIM)
    coord[1] = r % Nt
    r = div(r, Nt)
    for i in 2:STDIM
        coord[i] = r % Ns
        r = div(r, Ns)
    end
    return coord
end

# Funzione per calcolare il campo scalare a momento fissato
function phi_at_fixed_momentum(lattice, p, Nt, Ns, STDIM)
    stvol = Nt * Ns^(STDIM-1)
    O = Complex{Float64}[0.0 + 0.0im for _ in 1:Nt]

    for r in 1:stvol
        coord = lex_to_cart_st(r-1, Nt, Ns, STDIM)
        sum = 0.0
        for i in 2:STDIM
            sum += p[i-1] * coord[i]
        end
        O[coord[1]+1] += lattice[r] * exp(-1im * sum)  # Nota: uso di -1im per il fattore di fase corretto
    end

    for i in 1:Nt
        O[i] /= sqrt(stvol)
    end
    return O
end
=#
function phi_at_fixed_momentum(lattice, p, Nt, Ns, STDIM)
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
#=
# Funzione per calcolare i correlatori temporali a momento zero
function time_corr_meas(lattice, Nt, Ns, STDIM)
    p = zeros(Float64, STDIM-1)
    corr_p0 = zeros(Float64, div(Nt, 4))

    O = phi_at_fixed_momentum(lattice, p, Nt, Ns, STDIM)

    for i in 1:div(Nt, 4)
        for j in 1:Nt
            k = (j + i - 1) % Nt + 1
            corr_p0[i] += real(conj(O[j]) * O[k])
        end
        corr_p0[i] /= Nt
    end
    return corr_p0
end
=#
function time_corr_meas(lattice, Nt, δt::Int, STDIM)
    p = zeros(Float64, STDIM-1)
    corr_p0 = zeros(Float64, δt )

    O = phi_at_fixed_momentum(lattice, p, Nt, Ns, STDIM)

    for i in 1:δt
    for j in 1:Nt
        #k = (j + i - 1) % Nt + 1
        k = mod1(j+i, Nt)
        corr_p0[i] += real(conj(O[j]) * O[k])
    end
    corr_p0[i] /= Nt
end
return corr_p0
end
δt=div(Nt, 4)-1
# Parametri del reticolo
Nt = 16  # Numero di punti nel reticolo temporale
Ns = 16  # Numero di punti nel reticolo spaziale
STDIM = 4  # Dimensione dello spaziotempo

# Creiamo un reticolo di esempio con valori casuali
lattice = randn(Float64, Nt * Ns^(STDIM-1))

# Calcoliamo la correlazione temporale a momento zero
corr_p0 = time_corr_meas(lattice, Nt, 1, STDIM)
println(corr_p0)

# Plot dei correlatori temporali a impulso zero
time_steps = 0:div(Nt, 4)-1
plot(time_steps, corr_p0, marker=:o, label="Correlatori temporali a impulso zero", xlabel="Tempo", ylabel="Correlatore", title="Correlatori Temporali a Impulso Zero")
