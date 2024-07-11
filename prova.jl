using CSV, DataFrames, Dates, Printf
include("free_scalar.jl")



matrix = ([1.0 2.0 3.0 4.0 5.0; 1.0 1.0 1.0 1.0 1.0; 2.0 2.0 1.0 1.0 1.0; 3.0 2.0 1.0 2.0 3.0])

Nt= size(matrix,1) #ho 4 righe (4 è la lunghezza della dimensione temporale)
Ns= size(matrix, 2)
matrix=reshape(matrix, :)
circshift(matrix,-Nt)-matrix

matrix = ([1.0 2.0 3.0 4.0 5.0; 1.0 1.0 1.0 1.0 1.0; 2.0 2.0 1.0 1.0 1.0; 3.0 2.0 1.0 2.0 3.0])
rows = [collect(row) for row in eachcol(matrix)]
for i in 1:4
    println(rows[i])
    println(circshift(rows[i],+1)-rows[i])
end
mod1(7,4)
matrix[]
(mod1(10,4))
cld(10,4)
mod1(1,60)