using CSV, DataFrames, Dates
include("free_scalar.jl")



matrix = ([1.0 2.0 3.0 4.0; 5.0 6.0 7.0 8.0; 9.0 10.0 11.0 12.0])
Ns= size(matrix, 2)

Nt= size(matrix,1)


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

Σ_n(matrix, 1)
