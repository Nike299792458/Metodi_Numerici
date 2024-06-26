using CSV, DataFrames, Dates
include("functions.jl")


matrix = ([1 2 3; 4 5 6; 7 8 9])
#matrix= reshape(matrix, 9, 1)
show(matrix)
nuova=show(circshift(matrix ,-1 ))
O3(1,1,matrix)
array_righe = [collect(row) for row in eachcol(matrix)]
array_righe[1]
circshift(array_righe[1],1)
for i in (1,2,3)
    ris=dot(array_righe[i], circshift(array_righe[i],-1))
    println(ris)
end
