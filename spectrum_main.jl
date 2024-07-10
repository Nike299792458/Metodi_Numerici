using ArgParse, Dates, DelimitedFiles, Random, Printf
include("spectrum_c.jl")

function parse_cmd()
    s = ArgParseSettings()

    @add_arg_table s begin
        "Nt"
            help = "the number of time steps"
            required = true
            arg_type = Int
        "sample"
            help = "Number of steps in the simulation"
            required = true
            arg_type = Int
        "--path", "-p"
            help = "the path where files are stored"
            default = joinpath(["..", "simulations_c"])
            required = false
            arg_type = String
        "--verbose", "-v"
            help = "if given, percentage of execution is printed"
            action = :store_true
    end
    return parse_args(s)
end

function main()
    parsed_args = parse_cmd()
    Nt = parsed_args["Nt"]
    sample = parsed_args["sample"]
    path = parsed_args["path"]
    verbose = parsed_args["verbose"]
    #Tonm fa le veci di β*h*ω=5
    Tonm=5
    Nt=60
    Ns=10
    mhat=1/(Nt*Tonm) #mhat che va a 0 è l'analogo di η a 0
    # initializing...
    lattice = zeros(Float64, Nt, Ns)
    acc = 0.

    # simulation parameters
    orsteps = 5
    measevery = 10

    # files management
    if !isdir(path)
        mkpath(path)
    end
    fname = @sprintf "spectrum_sample=%.1eNt=%2.2i.txt" sample Nt
    fr = joinpath([path, fname])
    if !isfile(fr)
        touch(fr)
    end
    # writing header
    ct = permutedims(["ct_$i" for i in 0:Nt÷4-1 ])
    open(fr, "w") do infile
        writedlm(infile, [ct], " ")
    end
    
    start = now()
    datafile = open(fr, "a")
    for iter in 1:sample
        for r in LinearIndices(lattice)
                    acc+=heathbath!(lattice, r, mhat, Nt, Ns)
                for _ in 1:orsteps
                    acc+=overrelax!(lattice, r, mhat, Nt, Ns)
                end
        end
        #t_corr in δt è un array con l'unica entrata diversa da 0 in posizione δt e si sfrutta questo
        if iter % measevery == 0
            corr_matrix = t_corr(lattice, Nt, Ns,1)[1]#Julia inizia a contare da uno corrisponde a δt=0 
            for δt in 1:Nt÷4-1
                corr_matrix = hcat(corr_matrix, t_corr(lattice, Nt, Ns,δt)[δt] ) 
            end
             writedlm(datafile, corr_matrix, " ")
        end
        
        if verbose && iter % (sample÷100) == 0
            print("$((100*iter÷sample))% \r")
        end
    end
    close(datafile)
    elapsed = Dates.canonicalize(Dates.round((now()-start), Dates.Second))
    println("\n$(round(now(), Dates.Second));\nNₜ = $Nt, elapsed time $(elapsed)\n")
end
main()
