set logscale xy
set xlabel "Time (Generations)
set ylabel "N_e"
plot "epos4.dat" using 1:2 title "5% quantile"  with lines,\
     "epos4.dat" using 1:3 title "median"       with lines,\
     "epos4.dat" using 1:4 title "95% quantile" with lines
