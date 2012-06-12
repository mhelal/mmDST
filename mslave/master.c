#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <mpi.h>
#include <semaphore.h>
#include <errno.h>
#include "master.h"
#include "slave.h"
#include "lq.h"
#include "moaDst.h"
#include "utils.h"


void * MRecvThread(void *m) {
	char msg[MID_MESSAGE_SIZE];
  int procToEnq, errCode;
  Master * master;
  int notdone, i, ComputationPhase;
  master = (Master *)m;

  notdone = 0;
  /* first 2 Waves need to be done by only one slave proc 
	procToEnq = 1;
  master->ComputationPhase[1] = 1;
	while ((master->waveNo <= 3) && (notdone == 0)) {
    errCode = pthread_mutex_lock (master->mut);
    checkMutexErrorCode (1, errCode);
    enqueue(&processors, procToEnq);
		master->ComputationPhase[procToEnq] = ComputationPhase;
    errCode = pthread_mutex_unlock (master->mut);
    checkMutexErrorCode (1, errCode);
    sem_post (&master->qSem);
		// Then Receive the proc again to re-assign to it 
		mprintf (2, "Master Receiving Thread waiting for First Proc in Wave <=2\n", 2);
    procToEnq = ReceiveSlavetoProcGate(&ComputationPhase);
		master->ComputationPhase[procToEnq] = ComputationPhase;
    // Check if done
		sprintf (msg, "Master Received First Proc in Wave <=2 with ComputationPhase %d procToEnq %d comp %d\n", ComputationPhase, procToEnq, master->ComputationPhase[procToEnq]);
		mprintf (2, msg, 2);
    if (master->ComputationPhase[procToEnq] > 1) 
			notdone = 1;	// then done
	}
  /* After Wave 2 Add all slaves to queue */
  for (procToEnq=1;procToEnq<ClusterSize;procToEnq++) {
    errCode = pthread_mutex_lock (master->mut);
    checkMutexErrorCode (1, errCode);
    enqueue(&processors, procToEnq);
    master->ComputationPhase[procToEnq] = 1;
    errCode = pthread_mutex_unlock (master->mut);
    checkMutexErrorCode (1, errCode);
    //sem_post (&master->qSem);
  }
  sprintf (msg, "Added %d processors to the queue\n", ClusterSize);
  mprintf (1, msg, 2);
  sem_post (&master->qSem);
  //notdone = 0;
  /*
  while (notdone == 0) {
    procToEnq = ReceiveSlavetoProcGate(&ComputationPhase);
    //errCode = pthread_mutex_lock (master->mut);
    //master->ComputationPhase[procToEnq] = ComputationPhase;
    //enqueue(&processors, procToEnq);
    //errCode = pthread_mutex_unlock (master->mut);
    //sem_post (&master->qSem);
    // Check if done
    notdone = 1; // assume we are done
    for (i=1;i<ClusterSize;i++) {
      if (master->ComputationPhase[i] < 2) // if one processor not yet done
				notdone = 0;	// then not done
    }    
  }  
	*/
  mprintf(1, "Master Proc Gate Thread Leaving\n", 2);  
  return NULL;
}

void * MPartThread (void *m) {
  Master * master;
  int newProcessor, notdone;
  int errCode;
  int ComputationRounds, ComputationPhase;
  struct timeval timeout;  
  long oldprevWavePartStart, prevWavePartStart, currWavePartStart, nextWavePartStart;
  int MPI_return;
	char msg[MID_MESSAGE_SIZE];
	long i, j, seqNum;
int  stype;

  master = (Master *)m;

  ComputationPhase = 1;
  // 1. Synchronize Computation Phase 1 and Send Partitions
  DS_ubound = 0;
  master->waveNo = 1;
  oldprevWavePartStart = prevWavePartStart = 0;
  currWavePartStart = FirstWaveDistributeDiagonals (master, master->waveNo, partitionSize);
  while (currWavePartStart  < DS_ubound) {
    currNow = getTime();
    if (((currNow->tm_yday * 1440) + (currNow->tm_hour * 60) + currNow->tm_min) > ((prevNow->tm_yday * 1440) + (prevNow->tm_hour * 60) + prevNow->tm_min + 10)) {
      
      prevNow->tm_hour = currNow->tm_hour;
      prevNow->tm_isdst = currNow->tm_isdst;
      prevNow->tm_mday = currNow->tm_mday;
      prevNow->tm_min = currNow->tm_min;
      prevNow->tm_mon = currNow->tm_mon;
      prevNow->tm_sec = currNow->tm_sec;
      prevNow->tm_wday = currNow->tm_wday;
      prevNow->tm_yday = currNow->tm_yday;
      prevNow->tm_year = currNow->tm_year;
    }
    master->waveNo ++;
    sprintf(msg,"Master Thread 2 going to wave %ld \n", master->waveNo);
    mprintf(8, msg, 1);
    nextWavePartStart = LaterWavesDistributeDiagonals (currWavePartStart, master, master->waveNo, partitionSize);
    sprintf(msg,"Master Thread 2 after wave %ld \n", master->waveNo);
    mprintf(2, msg, 1);
    if (master->waveNo > 2) {
    sprintf(msg, "Master Thread 2 before manageDependancies at wave %ld \n", master->waveNo);
    mprintf(2, msg, 1);
    manageDependancies (master, partitionSize, master->waveNo - 2, oldprevWavePartStart, prevWavePartStart);
    sprintf(msg,"Master Thread 2 after manageDependancies at wave %ld \n", master->waveNo);
    mprintf(2, msg, 1);
    }
    oldprevWavePartStart = prevWavePartStart;
    prevWavePartStart = currWavePartStart; 
    currWavePartStart = nextWavePartStart;
    sprintf (msg, "old %ld prev %ld curr %ld next %ld", oldprevWavePartStart, prevWavePartStart, currWavePartStart,  nextWavePartStart);
    mprintf(8, msg, 1);
  }
  
    mprintf(2, "Master Thread 2 before before last  manageDependancies  \n", 1);
    manageDependancies (master, partitionSize, master->waveNo - 1, oldprevWavePartStart, prevWavePartStart);
    mprintf(2, "Master Thread 2 before last  manageDependancies  \n", 1);
    manageDependancies (master, partitionSize, master->waveNo, prevWavePartStart, nextWavePartStart);
    sprintf(msg, "Master Thread 2 after last  manageDependancies  last wave %ld\n", master->waveNo);
  	mprintf(2, msg, 1);

  
  // 2. Synchronize Computation Phase 2
  
  ComputationPhase = 2;
  SendEndPartitioingSignal();
  mprintf(1, "Master Partitioning Thread sending the end signal & Leaving\n", 1);  
	
	for (i = 0; i < DS_ubound; i++) {
		sprintf (msg, "Part %ld sent to %d \n", i, DS[i].proc);
		mprintf (10, msg, 1);
		sprintf (msg, "Part %ld has LB %ld: ", i, DS[i].lbCellsCount);
		mprintf (20, msg, 1);
		for (j=0;j<DS[i].lbCellsCount ; j++) {   
			sprintf (msg, "%ld ", DS[i].lbCells[j]);
			mprintf (20, msg, 1);
		}
		sprintf (msg, " and HB %ld: ", DS[i].hbCellsCount);
		mprintf (20, msg, 1);
		for (j=0;j<DS[i].hbCellsCount ; j++) {   
			sprintf (msg, " %ld ", DS[i].hbCells[j]);
			mprintf (20, msg, 1);
		}
		mprintf (20, "\n", 1);
	}
	
  return NULL;
}

void * MSendThread (void *m) {
  return NULL;
}




void MasterProcess (long seqNum, char * * sequences, long * seqLen, int stype) {
  pthread_t SendProc, PartProc, RecvProc;
  int MPI_return, ret;
  Master * master;
  long i, j;

	ret = initMasterMemory(&master, seqNum, seqLen, stype, sequences);  
	if (ret != 0) {
		mprintf (1, "Error Initializing Master Data, Exiting\n", 1);
		fflush (stdout);
		return;
	}
  //pthread_create (&SendProc, NULL, MSendThread, master);

  pthread_create (&PartProc, NULL, MPartThread, master);
  pthread_create (&RecvProc, NULL, MRecvThread, master);

  // Join Master Threads
  //pthread_join (SendProc, NULL);
  //mprintf(1, "Master Joined Sending Thread 1\n");


  pthread_join (PartProc, NULL);
  mprintf(1, "Master Joined Partitiong Thread 2\n", 1);
  pthread_join (RecvProc, NULL);
  mprintf(1, "Master Joined Receiving Thread 1\n", 1);

  checkPointMaster (master);
	
  freeMasterMemory (master);
	/*
	ret = initMasterMemory(&master, seqNum, seqLen, stype, sequences);  
	restoreMasterCheckPoint (master);
  freeMasterMemory (master);
*/
}

int main (int argc, char * argv[]) {
  MPI_Group orig_group;
	int ret;
  int stype = 0, MPI_return;
  char * * sequences = NULL, msg[500];
  long * seqLen = NULL;
  long seqNum, i, j;
	char ufilename[SHORT_MESSAGE_SIZE];
	
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &myProcid);
  MPI_Comm_size (MPI_COMM_WORLD, &ClusterSize);
  MPI_Comm_group (MPI_COMM_WORLD, &orig_group);
  MPI_Comm_create(MPI_COMM_WORLD, orig_group, &MOAMSA_COMM_WORLD);

  // 1. Process Arguments
  processArguments(argc, argv, &seqNum, &sequences, &seqLen, &stype);
	
	strcpy (ufilename, outputfilename);
	sprintf (outputfilename, "mmsa_%s", ufilename);
  init_output();
	sprintf (msg, "Program Arguments: debuglevel = %d maxAlignmentsNumber = %d Epsilons= %d Alignment Type = %d stype %d outputfilename = %s partitionSize = %d\n", pdebug, maxAlignmentsNumber, Epsilons, AlignmentType, stype, outputfilename, partitionSize);
	mprintf (1, msg, 1);

  Mode = Distributed;
  threadnum = 1;

  prevNow = NULL;
  prevNow = (struct tm *) mmalloc (sizeof(struct tm));
  currNow = getTime();
  prevNow->tm_hour = currNow->tm_hour;
  prevNow->tm_isdst = currNow->tm_isdst;
  prevNow->tm_mday = currNow->tm_mday;
  prevNow->tm_min = currNow->tm_min;
  prevNow->tm_mon = currNow->tm_mon;
  prevNow->tm_sec = currNow->tm_sec;
  prevNow->tm_wday = currNow->tm_wday;
  prevNow->tm_yday = currNow->tm_yday;
  prevNow->tm_year = currNow->tm_year;

	/*
  for (i=1;i<ClusterSize;i++) {
  /// 2. Broadcast DataSize
 	  MPI_return = MPI_Send ( &stype, 1, MPI_TYPE_DSSType, i, MPI_TAG_DSSType, MOAMSA_COMM_WORLD);

  	checkMPIErrorCode (MPI_return);
	  MPI_return = MPI_Send (&seqNum, 1, MPI_TYPE_DSDimn, i, MPI_TAG_DSDimn, MOAMSA_COMM_WORLD);
  	checkMPIErrorCode (MPI_return);
	  for (j=0;j<seqNum;j++) {
			MPI_return = MPI_Send (&seqLen[j], 1, MPI_TYPE_DSShape, i, MPI_TAG_DSShape, MOAMSA_COMM_WORLD);
  		checkMPIErrorCode (MPI_return);
	 	}
	}
	*/
  if (myProcid == 0)
    MasterProcess(seqNum, sequences, seqLen, stype);
  else
    SlaveProcess (seqNum, seqLen, stype);

  
  mprintf(1, "Master before finalize \n", 1);
  fflush(stdout);
  MPI_Finalize ();
  mprintf(1, "Master after finalize \n", 1);
  fflush(stdout);
  return 0;
}

