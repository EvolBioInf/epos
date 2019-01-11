set logscale xy
set xlabel "Time (Generations)
set ylabel "N_e"
plot[10**2:5*10**4] "epos2.dat" using 1:2 title "5% quantile"  with lines,\
                    "epos2.dat" using 1:3 title "median"       with lines,\
                    "fig2b_e.dat"         title "expected"     with lines,\
                    "epos2.dat" using 1:4 title "95% quantile" with lines
