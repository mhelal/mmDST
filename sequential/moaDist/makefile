moaDist : moaDst.o moamsa.o moa.o scores.o utils.o lq.o
	mpicc -o moaDist -lpthread ../moaDst.c ../moamsa.c ../moa.c \
	../scores.c ../utils.c ../lq.c

moaDst.o : ../moaDst.c ../moa.h ../scores.h ../utils.h ../moamsa.h
	 mpicc -c ../moaDst.c

moamsa.o : ../moamsa.c ../moa.h ../scores.h ../utils.h ../moamsa.h
	cc -c ../moamsa.c

moa.o :  ../moa.c ../moa.h ../utils.h ../moamsa.h ../lq.h
	cc -c ../moa.c

scores.o : ../scores.c ../scores.h
	cc -c ../scores.c

utils.o : ../utils.c ../moa.h ../utils.h ../moamsa.h
	cc -c ../utils.c

lq.o : ../lq.c ../lq.h
	cc -c ../lq.c

clean : 
	rm moaDist  moadist.o moamsa.o moa.o scores.o utils.o lq.o