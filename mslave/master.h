#ifndef MASTERH_
#define MASTERH_

#include <semaphore.h>
#include "globals.h"
#include "moa.h"

typedef struct Master_rec {
	long waveNo;
  MOA_rec * msaAlgn;
  int * ComputationPhase;
  pthread_mutex_t *  mut;
  sem_t qSem;
  long seqNum;
  char * * sequences;
  long * seqLen;
  int stype;
} Master;


typedef struct DS_struct {
  int proc; /* The processor ID*/
  unsigned long partIndex; /* The partition ID is the flat index of the first cell*/
  long waveNo;
  int sent;
  long * lbCells; /* lower Border cells */
  long * hbCells; /* higher Border cells*/
  long lbCellsCount;
  long hbCellsCount;
} Distribution;

Distribution * DS;
long DS_ubound;

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

int initMasterMemory (Master * * master, long seqNum, long * seqLen, int stype, char * * sequences);
void freeMasterMemory (Master * master);
int checkPointMaster (void * threadarg);
int restoreMasterCheckPoint (void * threadarg);
long LaterWavesDistributeDiagonals (long nextWavePartStart, Master * master, long waveNo, long waveSize);

int FirstWaveDistributeDiagonals (Master * master, long waveNo, long waveSize);

void getNextWave (Master * master, long waveNo, MOA_rec * MOA_nghb);

void manageDependancies (Master *  master, long partitionSize, long waveNo, long startPartNo, long EndPartNo);

void * MSendThread (void *args);
void * MPartThread (void *args);
void * MRecvThread (void *args);
void MasterProcess (long seqNum, char * * sequences, long * seqLen, int stype);

void SendPartition (Master * master, int toProcessor, MOA_rec * MOA_partition, long waveNo);
int sendCellDepProc (long waveNo, int toProc, unsigned long cellIndex, unsigned long partIndex, int depProc);
int SendEndPartitioingSignal();
int ReceiveSlavetoProcGate (int * ComputationPhase);
int sendMOAIndices (Master * master, int ProcDest, MOA_rec * MOA_in, long waveNo);
#endif
