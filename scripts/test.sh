./epos -l 10000000 -u 1.2e-8 -U data/testNewtonU.dat > tmp.out
DIFF=$(diff tmp.out data/testNewtonU.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(Newton, unfolded)\t\t\tpass\n"
else
    printf "Test(Newton, unfolded)\t\t\tfail: %s\n" ${DIFF}
fi

./epos -l 10000000 -u 1.2e-8 data/testNewtonF.dat > tmp.out
DIFF=$(diff tmp.out data/testNewtonF.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(Newton, folded)\t\t\tpass\n"
else
    printf "Test(Newton, folded)\t\t\tfail: %s\n" ${DIFF}
fi

./epos data/kap144i.dat > tmp.out
DIFF=$(diff tmp.out data/kap144i.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(Newton, kap144i, greedy)\t\tpass\n"
else
    printf "Test(Newton, kap144i, greedy)\t\tfail: %s\n" ${DIFF}
fi

./epos -E 3 data/kap144i.dat > tmp.out
DIFF=$(diff tmp.out data/kap144ie.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(Newton, kap144i, exhaustive)\tpass\n"
else
    printf "Test(Newton, kap144i, exhaustive)\tfail: %s\n" ${DIFF}
fi

rm tmp.out
