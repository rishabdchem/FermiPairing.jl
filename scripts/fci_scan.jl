using FermiPairing 

#========================================================
Exact ground state with scan.
Applicable to small systems only.
The following relation would be helpful:

Number of points 
= 1 + (final - initial)/distance. 
========================================================#

let
    
    # Print
    println("Start the DOCI scan")

    # System inputs
    m = 8
    n = 4
    ham = :XXZ
    bc = :OBC

    # Scan inputs
    hpinit = -1.0 
    dhp = 0.01
    nP = 5 

    # Initialize 
    efci = zeros(Float64, nP)
    Gval = zeros(Float64, nP)
    Gval[1] = copy(hpinit) 

    # Start the loop
    for j = 1:nP
      
        # Print
        println("GorD value: ", Gval[j])

        # FCI
        efci[j], evec = pairingfci(m, n, Gval[j], ham, bc)
 
        # Update G
        j == nP && continue 
        Gval[j + 1] = Gval[j] + dhp

    end

    # Print
    f = open("results/ScanOutput.txt", "w")
    println(f, "-------------------------")
    println(f, "m: $m, n: $n")

    println(f, "-------------------------")
    if ham == :RBCS
        println(f, "RBCS G values:")
    else
        println(f, "1D XXZ $bc D values:")
    end
    for i = 1:nP
        #println(f, Gval[nP-i+1])
        println(f, Gval[i])
    end
 
    println(f, "-------------------------")
    println(f, "FCI energies:")
    for i = 1:nP
        #println(f, efci[nP-i+1])
        println(f, efci[i])
    end

    # Close output
    close(f)

    # Print
    println("Check output file")

    # Exit    
    nothing    
end

