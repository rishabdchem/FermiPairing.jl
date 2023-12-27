using FermiPairing 

#==============================
BTS ground state.
==============================#

let
    
    # Print
    println("Start the BTS method")

    # System inputs
    m = 8
    n = 4
    G = 0.5 
    ham = :RBCS
    bc = :OBC

    # Guess
    etavec = pagpfrompcid(m, n, G, ham, bc) 
    eagp, eta = pairingagp(m, n, G, etavec, 5000, ham, bc)
    eagp = agphamoverlap(eta, eta, n, G, ham, bc) / agpoverlap(eta, eta, n)
    btsguess = btsfromagp(eta, n)

    # Optimize
    en, btsmat = btsmin(m, n, G, btsguess, 5000, ham, bc)

    # Print
    println("-------------------------")
    println("m: $m, n: $n")
    if ham == :RBCS
        println("RBCS with G: $G")
    else
        println("1D XXZ with $bc and D: $G")
    end
    println("-------------------------")
    println("eagp: $eagp")
    println("ebts: $en")
    println()

    # Exit    
    nothing    
end    

