set terminal epslatex color solid
set size 5/5., 4/3.
set format xy "\\large $10^{%T}$"
set output "kap.tex"
set logscale xy
set xlabel "\\large Time (Generations)
set ylabel "\\large $N_{\\rm e}$"
plot "epos4.dat" using 1:2 title "\\large 5\\% quantile"  with lines lw 3,\
     "epos4.dat" using 1:3 title "\\large median"       with lines lw 3,\
     "epos4.dat" using 1:4 title "\\large 95\\% quantile" with lines lw 3
