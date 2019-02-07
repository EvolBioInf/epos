set terminal epslatex color solid
set size 5/5., 4/3.
set format xy "\\large $10^{%T}$"
set output "fig2b.tex"
set logscale xy
set xlabel "\\large Time (Generations)
set ylabel "\\large $N_{\\rm e}$"
plot[10**2:5*10**4] "epos2.dat" using 1:2 title "\\large 5\\% quantile"  with lines lw 3,\
                    "epos2.dat" using 1:3 title "\\large median"       with lines lw 3,\
                    "fig2b_e.dat"         title "\\large expected"     with lines lw 3,\
                    "epos2.dat" using 1:4 title "\\large 95\\% quantile" with lines lw 3
