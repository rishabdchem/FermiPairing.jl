using FermiPairing 

#========================================================
AGP-based methods with scan.
The following relation would be helpful:

Number of points 
= 1 + (final - initial)/distance. 
========================================================#

let
    
    # Print
    println("Start the AGP methods scan")

    # System inputs
    m = 16
    n = 8
    ham = :RBCS
    bc = :PBC

    # Methods inputs
    iteragp = 5000 
    pval = 0.0
    jk = 2

    # Scan inputs
    hpinit = 0.02 
    dhp = 0.02
    nP = 2

    # Initialize 
    eagp = zeros(Float64, nP)
    ejci = zeros(Float64, nP)
    zeromd = zeros(Int, nP)
    Gval = zeros(Float64, nP)
    Gval[1] = hpinit 
    
    # Start the loop
    for j = 1:nP
      
        # Print
        println("GorD value: ", Gval[j])

        # AGP guess 
        etaguess = pagpfrompcid(m, n, Gval[j], ham, bc)
    
        # AGP 
        @time eagp[j], eta = pairingagp(m, n, Gval[j], etaguess, iteragp, ham, bc)
        
        # J-CI-AGP  
        @time ejci[j], evec, zeromd[j] = pairinglcagp_pivot(normalizeagp!(eta, n), n, pval, jk, Gval[j], ham, bc)
 
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
        println(f, "AGP energies:")
        for i in eagp
            println(f, i) 
        end

        println(f, "-------------------------")
        println(f, "LC-AGP energies:")
        for i in ejci 
            println(f, i) 
        end

        println(f, "-------------------------")
        println(f, "LC-AGP near-zero modes:")
        for i in zeromd 
            println(f, i) 
        end

    close(f)

    # Print
    println("Check output file")

    # Exit    
    nothing    
end

