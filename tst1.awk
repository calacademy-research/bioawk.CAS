function fillruns(runs) {
    if(ORD==78 || ORD==110) {
        if(Ns) {
            runs[++r] = Ns "N"
            Ns = 0; ACTGs = 0
        }
        ACTGs++; totACTGs++
    } else {
        if(ACTGs) {
            runs[++r] = ACTGs
            ACTGs = 0; Ns = 0
        }
        Ns++; totNs++
    }
}
function donothing(dummy){i+=ORD}
BEGIN{start=systime()}
0{
    delete runs
    print $name,length($seq)
    applytochars($seq, donothing(runs))
    exit
}

1{  ctg_start = systime()
    Ns = 0; ACTGs = 0; r= 0; delete runs
    applytochars($seq, fillruns(runs))
    runs[++r] = Ns ? Ns "N" : ACTGs
    nrun = r

    # compact into contigs
    ctgbreak = 10 # combine runs when Ns are < break
    for(c=i=1; i <= r; i++) {
        if(match(runs[i], "N$") && int(runs[i])>=ctgbreak) {
            contigs[++c] = runs[i]
            c++
        }
        else
            contigs[c] += int(runs[i])
    }
    ctgs=c
    num_contigs = (ctgs==1) ? ctgs : (ctgs-1)/2
    if(num_contigs==1) next

    # print "len "length($seq), "Ns "totNs, "other " totACTGs, "tNs+other = "totNs+totACTGs
    printf "%s: %d contigs, length %d [%ds]\n", $name, num_contigs, length($seq), (systime()-ctg_start)

    if(0) {
        printf "%s: ", $name
        for(r=1; r<=nrun; r++)
            printf "%s ", runs[r]
        printf "\n"
    }
    if (0) {
        for(c=1; c<=ctgs; c++)
            printf "%s ", contigs[c]
        printf "\n"
    }

     # exit
}
END {
    fin=systime()
    print "\ntotal time " (fin-start) "s"
}