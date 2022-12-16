using FermiPairing 

#=====================================
HF-based methods.
=====================================#

let
    
    # Print
    println("HF and CID methods start")

    # User inputs
    m = 10
    n = 5
    GorD = 0.3 
    ham = :RBCS
    bc = :OBC

    # HF 
    @time ehf = pairinghf(m, n, GorD, ham, bc)

    # CID
    @time eci, evec = pairingcid(m, n, GorD, ham, bc)

    # Print
    println("-------------------------")
    println("m: $m, n: $n")
    if ham == :RBCS
        println("RBCS with G: $GorD")
    else
        println("1D XXZ with $bc and D: $GorD")
    end
    println("-------------------------")
    println("ehf: $ehf")
    println("eci: $eci")
    println("-------------------------")
    println("CI coefficients:")
    println()
    display(evec) 
    println()

    # Exit    
    nothing    
end    

