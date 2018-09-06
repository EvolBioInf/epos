# fold.awk folds a site frequency spectrum
# Example: awk -f scripts/fold.awk data/testNewtonU.dat
# Author: Bernhard Haubold, haubold@evolbio.mpg.de
# Date: September 6, 2018
!/^#/ {
    f[++n] = $2
}
END {
    n++
    for(i=1; i<n/2; i++)
	f[i] += f[n-i]
    print "#r\tf(r)"
    for(i=1; i<=n/2; i++)
	printf "%d\t%d\n", i, f[i]
}
