./epos -l 1e7 -u 1.2e-8 -U ../data/testU.dat > tmp.out
DIFF=$(diff tmp.out ../data/testU.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(epos, unfolded, greedy)\t\tpass\n"
else
    printf "Test(epos, unfolded, greedy)\t\tfail\n"
    echo ${DIFF}
fi

./epos -l 1e7 -u 1.2e-8 ../data/testF.dat > tmp.out
DIFF=$(diff tmp.out ../data/testFg.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(epos, folded,   greedy)\t\tpass\n"
else
    printf "Test(epos, folded,   greedy)\t\tfail\n"
    echo ${DIFF}
fi

./epos -E 10 -l 1e7 -u 1.2e-8 ../data/testF.dat > tmp.out
DIFF=$(diff tmp.out ../data/testFe.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(epos, folded,   exhaustive)\tpass\n"
else
    printf "Test(epos, folded,   exhaustive)\tfail\n"
    echo ${DIFF}
fi

./epos ../data/kap144i.dat > tmp.out
DIFF=$(diff tmp.out ../data/kap144i.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(epos, kap144i,  greedy)\t\tpass\n"
else
    printf "Test(epos, kap144i,  greedy)\t\tfail\n"
    echo ${DIFF}
fi

./epos -E 3 ../data/kap144i.dat > tmp.out
DIFF=$(diff tmp.out ../data/kap144ie.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(epos, kap144i,  exhaustive)\tpass\n"
else
    printf "Test(epos, kap144i,  exhaustive)\tfail\n"
    echo ${DIFF}
fi

./epos -s 13 -k 5 ../data/kap144i.dat > tmp.out
DIFF=$(diff tmp.out ../data/kap144ix.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(epos, kap144i,  greedy, x-val)\tpass\n"
else
    printf "Test(epos, kap144i,  greedy, x-val)\tfail\n"
    echo ${DIFF}
fi

./epos -s 13 -E 3 -k 5 ../data/kap144i.dat > tmp.out
DIFF=$(diff tmp.out ../data/kap144iex.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(epos, kap144i,  exhaustive, x-val)\tpass\n"
else
    printf "Test(epos, kap144i,  exhaustive, x-val)\tfail\n"
    echo ${DIFF}
fi

rm tmp.out
