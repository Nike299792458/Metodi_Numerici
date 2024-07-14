using ArgParse, DelimitedFiles, LinearAlgebra, Printf, Statistics
include("free_scalar.jl")

doublers= false
blocksize = 3000
therm = 10000
ratio = 5
sample = 5000000
Nt = 8
Nt_b=Nt*ratio
path= "/Users/nicoletognetti/uni/Magistrale/MetodiNumerici/simulations_c/Tequalsm/"
cd(path)

if doublers == false
    dfname = @sprintf("Tequals_ratio=%.iNt=%i_sample=%.1e.txt",ratio, Nt, sample)
else 
    dfname = @sprintf("doublers_Nt=%i_sample=%.1e.txt", Nt, sample)
end


ϵ_norm_j=[]
ϵ_norm_b_j=[]
ϵ_norm = []
ϵ_normv = []
ϵ_norm_b =[]
ϵ_normv_b =[]
ϵ_norm_r = []
ϵ_normv_r = []
obs1 = []
obs1v = []
obs1_b = []
obs1v_b = []
obs1_r = []
obs1v_r = []


fname= @sprintf("fs_th_sample=%.1eratio=%.iNt=%2.2iTonm=1.00.txt", sample, ratio, Nt)
w = open(fname, "r") do io
    readdlm(io, header = true)
end

O1_j = Nt*Nt*JackKnife(w[1][therm:end,1], blocksize)
O2_j = Nt*Nt*JackKnife(w[1][therm:end,2], blocksize)
O3_j = Nt*Nt*JackKnife(w[1][therm:end,3], blocksize)
ϵ_norm_j=(O1_j+O2_j-O3_j)/2
        

push!(ϵ_norm, mean(ϵ_norm_j))
push!(ϵ_normv, std(ϵ_norm_j, corrected = false).*sqrt(length(ϵ_norm_j)-1))
push!(obs1 , mean(O1_j))
push!(obs1v , std(O1_j, corrected = false).*sqrt(length(O1_j)-1))


       
fname= @sprintf("fs_th_sample=%.1eratio=%.iNt_b=%2.2iNt=%2.2iTonm=1.00.txt",  sample, ratio, Nt_b, Nt)
w = open(joinpath([path, fname]), "r") do io
    readdlm(io, header = true)
end

O1_b_j = Nt*Nt*JackKnife(w[1][therm:end,1], blocksize)
O2_b_j = Nt*Nt*JackKnife(w[1][therm:end,2], blocksize)
O3_b_j = Nt*Nt*JackKnife(w[1][therm:end,3], blocksize)
ϵ_norm_b_j=(O1_b_j+O2_b_j-O3_b_j)/2
        

        
push!(ϵ_norm_b, mean(ϵ_norm_b_j))
push!(ϵ_normv_b, std(ϵ_norm_b_j, corrected = false).*sqrt(length(ϵ_norm_b_j)-1))
push!(obs1_b , mean(O1_b_j))
push!(obs1v_b , std(O1_b_j, corrected = false).*sqrt(length(O1_b_j)-1))


ϵ_norm_r = ϵ_norm - ϵ_norm_b
ϵ_normv_r = sqrt.((ϵ_normv ./ ϵ_norm).^2 .+ (ϵ_normv_b ./ ϵ_norm_b).^2) .* ϵ_norm_r
obs1_r= obs1 - obs1_b
obs1v_r = sqrt.((obs1v ./ obs1).^2 .+ (obs1v_b ./ obs1_b).^2) .* obs1_r


w = open(joinpath([path, dfname]), "w") do io
    writedlm(io, ["Nt"  "ϵ_norm" "ϵ_normv" "obs1" "obs1v" ], ",")
    writedlm(io, [Nt ϵ_norm_r ϵ_normv_r obs1_r obs1v_r], ",")
end
println("Done! Data stored in $(joinpath([path, dfname]))")




