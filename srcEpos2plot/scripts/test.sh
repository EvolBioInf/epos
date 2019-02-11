bzcat ../data/example.epos.bz2 | epos2plot -r > tmp.out
DIFF=$(diff tmp.out data/exampleR.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(epos2plot, raw)\t\tpass\n"
else
    printf "Test(epos2plot, raw)\t\tfail\n"
    echo ${DIFF}
fi

bzcat ../data/example.epos.bz2 | epos2plot > tmp.out
DIFF=$(diff tmp.out data/exampleQ.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(epos2plot, quantile)\tpass\n"
else
    printf "Test(epos2plot, quantile)\tfail\n"
    echo ${DIFF}
fi

bzcat ../data/example.epos.bz2 | epos2plot -t 1 > tmp.out
DIFF=$(diff tmp.out data/exampleQ1.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(epos2plot, quantile, step)\tpass\n"
else
    printf "Test(epos2plot, quantile, step)\tfail\n"
    echo ${DIFF}
fi

rm tmp.out
