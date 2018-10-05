set logscale xy
f(x)=1e6/2
set xlabel "Time (Generations)
set ylabel "N_e"
plot "epos1.dat" using 1:2 title "5% quantile"  with lines,\
     "epos1.dat" using 1:3 title "median"       with lines,\
     f(x)                  title "expected"     with lines,\
     "epos1.dat" using 1:4 title "95% quantile" with lines