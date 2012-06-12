mpicc -o mmDst -lpthread  -g  master.c mcheckp.c mpiMaster.c lq.c partitioning.c  moa.c utils.c moaDst.c slave.c scheckp.c  scoring.c  scores.c  mpiSlave.c
echo "Finished Compiling & Building Score Computation"
mpicc -o mtb -lpthread  -g  mtb.c mcheckp.c stb.c scheckp.c moaDst.c utils.c moa.c
echo "Finished Compiling & Building Trace Back"



