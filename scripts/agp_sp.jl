using FermiPairing 

#=====================================
AGP-based methods.
=====================================#

let
    
    # Print
    println("Start the AGP methods")

    # System inputs
    m = 8
    n = 4
    GorD = 0.3 
    ham = :RBCS
    bc = :OBC

    # Methods inputs
    iteragp = 5000 
    pval = 0.0
    jk = 2

    # AGP guess 
    etaguess = pagpfrompcid(m, n, GorD, ham, bc)

    # AGP 
    @time eagp, eta = pairingagp(m, n, GorD, etaguess, iteragp, ham, bc)
    
    # J-CI-AGP  
    @time ejci, evec, zeromd = pairinglcagp_pivot(normalizeagp!(eta, n), n, pval, jk, GorD, ham, bc)

    # Print
    println("-------------------------")
    println("m: $m, n: $n")
    if ham == :RBCS
        println("RBCS with G: $GorD")
    else
        println("1D XXZ with $bc and D: $GorD")
    end
    println("-------------------------")
    println("eagp: $eagp")
    println("ejci: $ejci")
    println("-------------------------")
    println("Normalized AGP coefficients:")
    println()
    display(normalizeagp!(eta, n)) 
    println("-------------------------")
    println("LC-AGP has $zeromd number of near-zero modes")
    println("NOCI coefficients:")
    println()
    display(evec) 
    println()

    # Exit    
    nothing    
end    

