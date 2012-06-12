pkill mtb
cd tbDst/
dmake -f makerelease
cd ../
mpirun -np 8 tbDst/mtb -c 3 ../data/tseq/seq10 ../data/tseq/seq10_2 ../data/tseq/seq10_3 -d 0 -o d3_1
