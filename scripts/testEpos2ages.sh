../build/epos2ages -n 144 ../data/kap144i.out > tmp.out
DIFF=$(diff tmp.out ../data/kap144ia.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(epos2ages)\tpass\n"
else
    printf "Test(epos2ages)\tfail\n"
    echo ${DIFF}
fi

