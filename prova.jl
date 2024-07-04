using CSV, DataFrames, Dates, Printf
include("free_scalar.jl")



matrix = ([1.0 2.0 3.0 4.0; 5.0 6.0 7.0 8.0; 9.0 10.0 11.0 12.0; 13.0 14.0 15.0 16.0])
Nt= size(matrix,1)
matrix=reshape(matrix,:)
circshift(matrix,-3)
Ns= size(matrix, 2)
rows = [collect(row) for row in eachcol(matrix)]
circshift(rows[1],-1)-circshift(rows[1],+1)


show(matrix)
rows = [collect(row) for row in eachcol(matrix)]
a=0
for i in 1:size(matrix,1)
    println(i)
end
for i = 1:size(matrix,1)
a=circshift(rows[i],-1)
println(a)
end


size(matrix, 1)#1 sta per righe
somma=0
for i = 1:3
    somma+=dot(array_righe[i], circshift(array_righe[i],-1))
end
println(somma)
STDIM= 3
Nt=3
Ns=4
stvolume=Nt
    for i in 1:STDIM
        stvolume=stvolume*Ns
        println(stvolume)
    end
println(stvolume)

timestamps=[4,5,6,7,8,10]

for i in timestamps
    # initializing...
    Nt=i
    println(Nt)
end


ratio=4
sample= 1000
doublers=false
Ns=40
Nt=4
Nt_b=10
if !isdir(path)
    mkpath(path)
end
#Apertura lettura e chiusura files
dfname = @sprintf("data_ratio=%i_sample=%.1e_doublers=%i.txt", ratio, sample, doublers)
startp1 = @sprintf "free_scalar_th_sample=%.1eratio=%.iNs=%iNt=%i" sample ratio Ns Nt
paths1 = filter(startswith(startp1), readdir(path))
startp2 = @sprintf "free_scalar_th_sample=%.1eratio=%.iNt_bar=%iNt=%i" sample ratio Nt_b Nt
paths2 = filter(startswith(startp2), readdir(path))
temporal_dim = []
ϵ_norm = []
ϵ_norm_b = []
ϵ_norm_r = []
for (i,fname) in enumerate(paths1)
    Nt = parse(Int, fname[end-5:end-4])

    local w = open(joinpath([path, fname]), "r") do io
        readdlm(io, header = true)
    end

    O1_j = JackKnife(w[1][therm:end,1], blocksize)
    O2_j = JackKnife(w[1][therm:end,2], blocksize)
    O3_j = JackKnife(w[1][therm:end,3], blocksize)
    ϵ_norm_j =(O1_j+O2_j-O3_j)/2

    push!(temporal_dim, Nt)
    push!(ϵ_norm, mean(ϵ_norm_j))

end

for (i,fname) in enumerate(paths2)
    Nt = parse(Int, fname[end-5:end-4])  

    local k = open(joinpath([path, fname]), "r") do io
        readdlm(io, header = true)
    end

    O1_b_j = JackKnife(k[1][therm:end,1], blocksize)
    O2_b_j = JackKnife(k[1][therm:end,2], blocksize)
    O3_b_j = JackKnife(k[1][therm:end,3], blocksize)
    ϵ_norm_b_j =(O1_b_j+O2_b_j-O3_b_j)/2

    push!(ϵ_norm_b, mean(ϵ_norm_b_j))

end

ϵ_norm_r= ϵ_norm-ϵ_norm_b

w = open(joinpath([path, dfname]), "w") do io
    writedlm(io, ["temporal_dim"  "ϵ_norm" ], ",")
    writedlm(io, [temporal_dim ϵ_norm_r ], ",")
end
println("Done! Data stored in $(joinpath([path, dfname]))")


