using FermiPairing 

#========================================================
HF-CID ground state with scan.
The following relation would be helpful:

Number of points 
= 1 + (final - initial)/distance. 
========================================================#

let
    
    # Print
    println("Start the HF-CID scan")

    # System inputs
    m = 16
    n = 8
    ham = :RBCS
    bc = :PBC

    # Scan inputs
    hpinit = 0.02 
    dhp = 0.02
    nP = 2

    # Initialize 
    ehf = zeros(Float64, nP)
    eci = zeros(Float64, nP)
    Gval = zeros(Float64, nP)
    Gval[1] = hpinit 
    
    # Start the loop
    for j = 1:nP
      
        # Print
        println("GorD value: ", Gval[j])
 
        # HF 
        ehf[j] = pairinghf(m, n, Gval[j], ham, bc)

        # CID
        eci[j], evec = pairingcid(m, n, Gval[j], ham, bc)

        # Update G
        j == nP && continue 
        Gval[j + 1] = Gval[j] + dhp

    end

    # Prepare outputs 
    f = open("results/ScanOutput.txt", "w")

        println(f, "-------------------------")
        println(f, "m: $m, n: $n")

        println(f, "-------------------------")
        if ham == :RBCS
            println(f, "RBCS G values:")
        else
            println(f, "1D XXZ $bc D values:")
        end
        for i in Gval 
            println(f, i) 
        end

        println(f, "-------------------------")
        println(f, "HF energies:")
        for i in ehf
            println(f, i) 
        end

        println(f, "-------------------------")
        println(f, "HF-CID energies:")
        for i in eci 
            println(f, i) 
        end

    close(f)

    # Print
    println("Check output file")

    # Exit    
    nothing    
end

