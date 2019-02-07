set terminal epslatex color solid
set output "con.tex"
set size 5/5., 4/3.
set format xy "\\large $10^{%T}$"
#set format xy "\\large$%g$"
set logscale xy
f(x)=1e6/2
set xlabel "\\large Time (Generations)
set ylabel "\\large $N_{\\rm e}$"
#set arrow nohead from 5*10**5,10**5 to 5*10**5,10**7
plot[10**3:][] "epos1.dat" using 1:2 title "\\large 5\\% quantile"  with lines lw 3,\
     "epos1.dat" using 1:3 title "\\large median"       with lines lw 3,\
     f(x)                 title "\\large expected"      with lines lw 3,\
     "epos1.dat" using 1:4 title "\\large 95\\% quantile" with lines lw 3
