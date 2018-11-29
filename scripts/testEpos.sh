../build/epos -l 10000000 -u 1.2e-8 -U ../data/testU.dat > tmp.out
DIFF=$(diff tmp.out ../data/testU.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(unfolded, greedy)\t\tpass\n"
else
    printf "Test(unfolded, greedy)\t\tfail\n"
    echo ${DIFF}
fi

../build/epos -l 10000000 -u 1.2e-8 ../data/testF.dat > tmp.out
DIFF=$(diff tmp.out ../data/testFg.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(folded,   greedy)\t\tpass\n"
else
    printf "Test(folded,   greedy)\t\tfail\n"
    echo ${DIFF}
fi

../build/epos -E 10 -l 10000000 -u 1.2e-8 ../data/testF.dat > tmp.out
DIFF=$(diff tmp.out ../data/testFe.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(folded,   exhaustive)\tpass\n"
else
    printf "Test(folded,   exhaustive)\tfail\n"
    echo ${DIFF}
fi

../build/epos ../data/kap144i.dat > tmp.out
DIFF=$(diff tmp.out ../data/kap144i.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(kap144i,  greedy)\t\tpass\n"
else
    printf "Test(kap144i,  greedy)\t\tfail\n"
    echo ${DIFF}
fi

../build/epos -E 3 ../data/kap144i.dat > tmp.out
DIFF=$(diff tmp.out ../data/kap144ie.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(kap144i,  exhaustive)\tpass\n"
else
    printf "Test(kap144i,  exhaustive)\tfail\n"
    echo ${DIFF}
fi

rm tmp.out
