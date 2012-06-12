moaLib : moaLb.o moa.o utils.o lq.o
	gcc ../moaLb.c ../moa.c ../utils.c ../lq.c -o moaLib

moaLb.o : ../moaLb.c
	cc -c ../moaLb.c

moa.o :  ../moa.c ../moa.h ../utils.h ../moamsa.h ../lq.h
	cc -c ../moa.c

utils.o : ../utils.c ../moa.h ../utils.h ../moamsa.h
	cc -c ../utils.c

lq.o : ../lq.c ../lq.h
	cc -c ../lq.c

clean : 
	rm moaLib moaLb.o moa.o utils.o lq.o