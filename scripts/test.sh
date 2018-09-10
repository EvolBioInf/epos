./epos -l 10000000 -u 1.2e-8 -U data/testNewtonU.dat > tmp.out
DIFF=$(diff tmp.out data/testNewtonU.out)
if [ "$DIFF" == "" ] 
then
    echo "Test(Newton, unfolded): pass"
else
    echo "Test(Newton, unfolded): fail" ${DIFF}
fi

./epos -l 10000000 -u 1.2e-8 data/testNewtonF.dat > tmp.out
DIFF=$(diff tmp.out data/testNewtonF.out)
if [ "$DIFF" == "" ] 
then
    echo "Test(Newton, folded): pass"
else
    echo "Test(Newton, folded): fail" ${DIFF}
fi

./epos data/kap144i.dat > tmp.out
DIFF=$(diff tmp.out data/kap144i.out)
if [ "$DIFF" == "" ] 
then
    echo "Test(Newton, kap144i): pass"
else
    echo "Test(Newton, kap144i): fail" ${DIFF}
fi
