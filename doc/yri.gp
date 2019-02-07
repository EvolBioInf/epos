set terminal epslatex color solid
set size 5/5., 4/3.
set logscale xy

set output "yriSi.tex"
# set format xy "\\Large $%.3t\\times 10^{%T}$"
# set format xy "\\large$%g$"
# set logscale x
set xlabel "\\large Years ($\\times 10^3$)
set ylabel "\\large $N_{\\rm e} (\\times 10^3)$"
plot[1:2000][1:10000] "yriSi.dat" using ($1*24/1000):($2/1000) title "\\large 2.5\\% quantile"  with lines lw 3,\
         "yriSi.dat" using ($1*24/1000):($3/1000) title "\\large median"           with lines lw 3,\
         "yriSi.dat" using ($1*24/1000):($4/1000) title "\\large 97.5\\% quantile" with lines lw 3

set output "yriNoSi.tex"
# set format xy "\\Large $%.3t\\times 10^{%T}$"
# set format xy "\\large$%g$"
# set logscale x
set xlabel "\\large Years ($\\times 10^3$)
set ylabel "\\large $N_{\\rm e} (\\times 10^3)$"
plot[1:2000][1:10000] "yriNoSi.dat" using ($1*24/1000):($2/1000) title "\\large 2.5\\% quantile"  with lines lw 3,\
         "yriNoSi.dat" using ($1*24/1000):($3/1000) title "\\large median"           with lines lw 3,\
         "yriNoSi.dat" using ($1*24/1000):($4/1000) title "\\large 97.5\\% quantile" with lines lw 3