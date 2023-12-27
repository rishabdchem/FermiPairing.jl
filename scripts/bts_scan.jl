using FermiPairing 
const fp = FermiPairing

#========================================================
BTS ground state with scan.
The following relation would be helpful:

Number of points 
= 1 + (final - initial)/distance. 
========================================================#

let
    
    # Print
    println("Start the BTS scan")

    # System inputs
    m = 12 
    n = 6
    ham = :RBCS
    bc = :PBC

    # Scan inputs
    hpinit = 0.02 
    dhp = 0.02
    nP = 2

    # Initialize 
    etaguess = zeros(Float64, m) 
    btsguess = zeros(Float64, btpdim(m, n))
    Gval = copy(hpinit)

    # Open output files
    f1 = open("results/ScanOutput.txt", "w")
    f2 = open("results/EtaScan.txt", "w")
    if ham == :RBCS
        println(f1, "RBCS with m: $m, n: $n")
        println(f2, "RBCS with m: $m, n: $n")
    else
        println(f1, "1D XXZ with m: $m, n: $n, bc: $bc")
        println(f2, "1D XXZ with m: $m, n: $n, bc: $bc")
    end
    println(f1, "-------------------------")
 
    # Start the loop
    for j = 1:nP
      
        # Guess
        if j == 1
            etaguess = pagpfrompcid(m, n, Gval, ham, bc)
            eagp, eta = pairingagp(m, n, Gval, etaguess, 5000, ham, bc)
            btsguess = btsfromagp(eta, n)
        end
    
        # Optimize
        @time ebts, btsmat = btsmin(m, n, Gval, btsguess, 10000, ham, bc)

        # Print
        println("-------------")
        println("GorD: $Gval, EBTS: $ebts")

        # Copy result 
        rGval = round(Gval, digits = 3)
        println(f1, "GorD: $rGval, EBTS:   $ebts")
        println(f2, "-------------------------")
        println(f2, "GorD: $Gval")
        println(f2, "-------------------------")
        for p in fp.btpmattovec(btsmat)
            println(f2, p) 
        end

        # Update 
        j == nP && continue 
        Gval += dhp
        btsguess = fp.btpmattovec(btsmat) 

    end

    # Close output
    close(f1)
    close(f2)

    # Print
    println("Check output files")

    # Exit    
    nothing    
end

