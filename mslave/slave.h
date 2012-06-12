#ifndef SLAVEH_
#define SLAVEH_

#include <semaphore.h>
#include "moa.h"

typedef struct MOAPart_rec {
  int processed; /*0 if not yet, 1 otherwise*/
  MOA_rec * msaAlgn;
  char * * sequences;
  long waveNo;
} MOAPartition;

typedef struct Slave_rec {
  int stype;
  long seqNum;
  long * seqLen;
  MOAPartition * MOAPart;
  long partitionsCount;
  long processedPartitions;
  int ComputationPhase;
	int compThreadFinished;
  pthread_mutex_t *CPmut;
  sem_t cpSem; /* to sync bet receive part & compute scores*/
  sem_t retSem; /* to sync bet compute scores & return to master*/
  int OCempty;
  pthread_mutex_t *OCmut;
  sem_t dpSem;
  int DPempty;
  pthread_mutex_t *DPmut;
  sem_t ocSem;
} Slave;


typedef struct OCI_Tag {
	long waveNo;
  long cellIndex; /* Overlapping Cell Index */
  long cellScore; /* Overlapping Cell Score */
  int fromProc; /* Overlapping Cell Score */
} OCIType;

/*Overlapping Cells incoming*/
long OCin_ub;
OCIType * OCin;

typedef struct OCO_Tag {
	long waveNo;
  long cellIndex; /* Overlapping Cell Index */
  long cellScore; /* Overlapping Cell Score */
  int toProc; /* Overlapping Cell Score */
} OCOType;


/*Overlapping Cells outgoing*/
long OCout_ub;
OCOType * OCout;
#ifndef depProc_tag_
#define depProc_tag_

typedef struct depProc_tag {
  unsigned long cellIndex; /* Overlapping Cell Index */
  unsigned long partIndex; /* Overlapping Partition Index */
  long waveNo; /* Overlapping Cell Wave No */
  int depProc; /* dep Proc */
  int sent; /* sent Flag */
} depProc_rec;
#endif

depProc_rec * depProclist;
long depProclistCount;

void distributedSlaveTraceBack (void * arg, char * * sequences); 
void sendMaxScore (void * arg);
 
int initSlaveMemory (Slave * * slave, long seqNum, long * seqLen, int stype);
void freeSlaveMemory (Slave * slave);

int checkPointSlave (void * threadarg);
int restoreSlaveCheckPoint (void * threadarg);

int RecvCompPhase(int * ComputationPhase);
int ReceivePartition (Slave  * slave);

int receiveCellDepProc (Slave * slave, long * depProclistCount, depProc_rec * * depProclist);

int SendSlavetoProcGate (int ComputationPhase );

int receiveOC (Slave * slave, long * OCin_ub, OCIType * * OCin);
int trySendOC (unsigned long cellIndex, long cellScore);
void sendOC (int toProc, unsigned long cellIndex, long cellScore);

int checkPrevPartitions (void * threadarg, unsigned long cellIndex, long * cellScore, long waveNo);
int checkRecvOC (unsigned long cellIndex, long * cellScore);

int receiveMOAIndices (MOAPartition  * MOAPart);

void * SSendThread (void *args);
void * SCompThread (void *args);
void * SRecvThread (void *args);
void SlaveProcess (long seqNum, long * seqLen, int stype);

#endif
