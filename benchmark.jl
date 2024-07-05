using ArgParse, Dates, DelimitedFiles, Random, Printf, LinearAlgebra
include("free_scalar.jl")
function main()
    ratio = 2
    sample = 10000000
    path = joinpath(["..", "simulations_c"])
   

    
    println(@sprintf "Starting simulation: sample=%.1e ratio=%.i " sample ratio )

   
        Nt=2
        Ns=2
        # initializing...
        lattice = zeros(Float64, Nt,Ns)
        acc=0
        
        # simulation parameters
        orsteps = 8
        measevery = 10


        #mhat::Float64=(1.0)/Nt
        mhat=1.0
      
      
        stvol=Nt*Ns #stvol=Nt*Ns^{STDIM-1}
       

        
        # files management
        if !isdir(path)
            mkpath(path)
        end
        fname = @sprintf("free_scalar_th_sample=%.1eratio=%.iNt=%2.2iNs=%2.2i.txt" , sample, ratio, Nt, Ns)
        fr = joinpath([path, fname])
        if !isfile(fr)
            touch(fr)
        end
        # writing header
        open(fr, "w") do infile
            writedlm(infile, ["obs1" "obs2" "obs3"], " ")
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

            if iter%measevery == 0
                obs1 = O1(stvol, mhat,lattice)
                obs2 = O2(stvol, Nt, lattice)
                obs3 = O3(stvol, lattice)
                writedlm(datafile, [obs1 obs2 obs3 ], " ")
            end
            
        end
        close(datafile)

    
       #= Nt_b = ratio*Nt
        Ns_b = Ns
        stvol_b = Nt_b # stvol=Nt*Ns^{STDIM-1}
        for i in 1:(STDIM-1)
            stvol_b = stvol_b * Ns_b
        end
        lattice_b = zeros(Float64, Nt_b, Ns_b)
        acc = 0

        # files management
        if !isdir(path)
            mkpath(path)
        end
        fname = @sprintf("free_scalar_th_sample=%.1eratio=%.iNt_b=%2.2iNt=%2.2i.txt" , sample, ratio, Nt_b, Nt)
        fr = joinpath([path, fname])
        if !isfile(fr)
            touch(fr)
        end
        # writing header
        open(fr, "w") do infile
            writedlm(infile, ["obs1_b" "obs2_b" "obs3_b"], " ")
        end

        datafile = open(fr, "a")
        for iter in 1:sample
            for r in LinearIndices(lattice_b)
                acc += heathbath!(lattice_b, r, mhat, Nt_b, Ns_b)
                for _ in 1:orsteps
                    acc += overrelax!(lattice_b, r, mhat, Nt_b, Ns_b)
                end
            end

            if iter % measevery == 0
                obs1_b = O1(stvol_b, mhat, lattice_b)
                obs2_b = O2(stvol_b, Nt_b, lattice_b)
                obs3_b = O3(stvol_b, lattice_b)
                writedlm(datafile, [obs1_b obs2_b obs3_b ], " ")
            end
        end
        close(datafile)

        =#
        elapsed = Dates.canonicalize(Dates.round((now() - start), Dates.Second))
        println("\n$(round(now(), Dates.Second));\nNₜ = $Nt,Ns = $Ns, elapsed time $(elapsed)\n")
    
end
main()