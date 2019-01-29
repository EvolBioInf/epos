# valgrind --track-origins=yes --leak-check=full --show-leak-kinds=all --log-file="epos.val" --dsymutil=yes ./epos -X 1 -x 1 -u 1.2e-8 -l 3e9 ../data/FoldSFS_YRI.txt

valgrind --track-origins=yes --leak-check=full --show-leak-kinds=all --log-file="epos.val" --dsymutil=yes ./epos -X 1 -x 1 -s 13 -u 1.2e-8 -l 1e7 ../data/testF.dat
