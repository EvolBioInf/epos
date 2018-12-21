# valgrind --track-origins=yes --leak-check=full --show-leak-kinds=all --log-file="epos.val" --dsymutil=yes ./epos ../Data/kap144*.dat
# valgrind --track-origins=yes --leak-check=full --show-leak-kinds=all --log-file="epos.val" --dsymutil=yes ./epos sim10.dat
# valgrind --track-origins=yes --leak-check=full --show-leak-kinds=all --log-file="epos.val" --dsymutil=yes ./epos -l 10000000 -u 1.2e-8 -b 5 -U -b 2 data/testNewtonU.dat
# valgrind --track-origins=yes --leak-check=full --show-leak-kinds=all --log-file="epos.val" --dsymutil=yes ./epos -l 10000000 -u 1.2e-8 -b 5 -b 2 data/testNewtonF.dat

# valgrind --track-origins=yes --leak-check=full --show-leak-kinds=all --log-file="testNewton.val" --dsymutil=yes ./testNewton -l 10000000 -U data/testNewtonU.dat
# valgrind --track-origins=yes --leak-check=full --show-leak-kinds=all --log-file="test.val" --dsymutil=yes ./epos -t
valgrind --track-origins=yes --leak-check=full --show-leak-kinds=all --log-file="epos.val" --dsymutil=yes ./epos -l 10000000 ../data/testF.dat


