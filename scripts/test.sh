./epos -l 10000000 -u 1.2e-8 -U data/testNewtonU.dat > tmp.out
DIFF=$(diff tmp.out data/testNewtonU.out)
if [ "$DIFF" == "" ] 
then
    echo "Test(Newton, unfolded): pass"
else
    echo "Test(Newton, unfolded): fail" ${DIFF}
fi
rm tmp.out

./epos -l 10000000 -u 1.2e-8 data/testNewtonF.dat > tmp.out
DIFF=$(diff tmp.out data/testNewtonF.out)
if [ "$DIFF" == "" ] 
then
    echo "Test(Newton, folded): pass"
else
    echo "Test(Newton, folded): fail" ${DIFF}
fi
rm tmp.out
