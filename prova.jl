using CSV, DataFrames, Dates
include("functions.jl")


matrix = ([1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0])
#matrix= reshape(matrix, 9, 1)
show(matrix)
#nuova=show(circshift(matrix ,-1 ))
#O3(1,1,matrix)
array_righe = [collect(row) for row in eachcol(matrix)]
#circshift(array_righe[1],1)
size(matrix, 1)#1 sta per righe
somma=0
for i = 1:3
    somma+=dot(array_righe[i], circshift(array_righe[i],-1))
end
println(somma)
STDIM= 3
Nt=3
Ns=2
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