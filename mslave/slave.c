#include <stdio.h>
#include <pthread.h>
#include <mpi.h>
#include <string.h>
#include <semaphore.h>
#include "slave.h"
#include "globals.h"
#include "moaDst.h"
#include "utils.h"
#include "scoring.h"

int TestEndofComputation(int ComputationPhase) {
return 0;
}

void * SSendThread (void * threadarg) {
  Slave * slave;
  int errCode;
	char msg[MID_MESSAGE_SIZE];

  slave = (Slave *) threadarg;
	sem_wait(&slave->retSem);
  while (slave->ComputationPhase == 1) {
    //SendSlavetoProcGate (slave->ComputationPhase);
    sprintf(msg, "received MOA Part %ld, starting %ld ComputationPhase %d , and returned \n", slave->partitionsCount, slave->MOAPart[slave->partitionsCount-1].msaAlgn->indexes[0], slave->ComputationPhase);
		mprintf(8, msg, 3);
	  sem_wait(&slave->retSem);
  	mprintf (3, "returned to Master\n", 3);
 }
  //SendSlavetoProcGate (slave->ComputationPhase);
  mprintf (3, "returned finally to Master\n", 3);
return NULL;
}

void * SCompThread (void * threadarg) {
  Slave * slave;
  int errCode;
	char msg[MID_MESSAGE_SIZE];
	slave = (Slave *) threadarg;
	
// first time outside the loop, to receive first
	sem_wait(&slave->cpSem);
	//errCode = pthread_mutex_lock (slave->CPmut);
  checkMutexErrorCode (1, errCode);
  if (slave->processedPartitions < slave->partitionsCount) {
  	slave->processedPartitions ++;
	  ComputePartitionScores (slave);
		slave->MOAPart[slave->processedPartitions-1].processed = 1;	  
		//printMOA(slave->MOAPart[slave->processedPartitions-1].msaAlgn);
		//printMOAIndices(slave->MOAPart[slave->processedPartitions-1].msaAlgn);
  }
	sprintf (msg, "computed %ld out of %ld\n", slave->processedPartitions, slave->partitionsCount);
	mprintf(3, msg, 1);
	//errCode = pthread_mutex_unlock (slave->CPmut);
  checkMutexErrorCode (1, errCode);
  sem_post (&slave->retSem);
	sprintf (msg, "Will Continue computing %ld out of %ld\n", slave->processedPartitions, slave->partitionsCount);
	mprintf(3, msg, 1);

	while ((slave->partitionsCount > 0) && (slave->processedPartitions <= slave->partitionsCount)) {
	//while (slave->ComputationPhase == 1) {
		sem_wait(&slave->cpSem);
	  if (slave->processedPartitions < slave->partitionsCount) {
  		slave->processedPartitions ++;
	    ComputePartitionScores (slave);
    	PrintPrevChains(slave->MOAPart[slave->processedPartitions-1].msaAlgn);
			errCode = pthread_mutex_lock (slave->CPmut);
  		checkMutexErrorCode (1, errCode);
			slave->MOAPart[slave->processedPartitions-1].processed = 1;
			errCode = pthread_mutex_unlock (slave->CPmut);
  		checkMutexErrorCode (1, errCode);
		  printMOA(slave->MOAPart[slave->processedPartitions-1].msaAlgn);
			printMOAIndices(slave->MOAPart[slave->processedPartitions-1].msaAlgn);
		}
		else /* if the signal is received and the partitions counter didn't increase, it means no more, and increment here to get out of the loop*/
  		slave->processedPartitions ++; 
	  sem_post (&slave->retSem);
  }
  
	slave->compThreadFinished = 1;
	sprintf (msg, "left computation Thread with %ld out of %ld\n", slave->processedPartitions, slave->partitionsCount);
	mprintf(3, msg, 1);
return NULL;

}

int checkPrevPartitions (void * threadarg, unsigned long cellIndex, long * cellScore, long waveNo) {
	long i, j;
	int found = -1;

	for (i=0;((i<OCout_ub) && (found == -1));i++) {
		if ((OCout[i].waveNo == waveNo) || (OCout[i].waveNo == waveNo - 1)) {
			if (cellIndex == OCout[i].cellIndex) {
				(*cellScore) = OCout[i].cellScore;
				found = 0;
			}
		}
	}
	return found;
}

int checkRecvOC (unsigned long cellIndex, long * cellScore) {
	long i, j;
	int received = -1;

	for (i=0;((i<OCin_ub) && (received == -1));i++) {
		if (cellIndex == OCin[i].cellIndex) {
			(*cellScore) = OCin[i].cellScore;
			received = 0;
		}
	}
	return received;
}
/*once a dependancy is received, this is called to check if it is computed to send it to waiting slave processors, otherwise, it will remain in the queue to be sent when it is computed, or before the slave exits*/

int trySendDPOC (depProc_rec *  depProc) {
	long i, j;
	//char msg[MID_MESSAGE_SIZE];

	for (i=0;((i<OCout_ub) && (depProc->sent == 0));i++) {
			//sprintf (msg, "trySendDPOC i %d cell %d score %d\n ", i, OCout[i].cellIndex, OCout[i].cellScore);
			//mprintf(3, msg, 3);
		if (depProc->cellIndex == OCout[i].cellIndex) {
			sendOC (depProc->depProc, OCout[i].cellIndex, OCout[i].cellScore);
			depProc->sent = 1;
		}	
	}
	return depProc->sent;
}

/*once HB cell score is computed, this is called to check if there are dependand processors to send it to them, otherwise, it will remain in the OC queue to be sent when a depProc is received, or before the slave exits*/

int trySendOC (unsigned long cellIndex, long cellScore) {
	long i, j;

	int sent = 0;
	for (i=0;i<depProclistCount;i++) {
		if ((depProclist[i].cellIndex == cellIndex) && (depProclist[i].sent == 0)) {
			sendOC (depProclist[i].depProc, cellIndex,cellScore);
			sent = 1;
			depProclist[i].sent = 1;
		}
	}
	return sent;
}

int	checkunSendOC () {
	long i, j;
	char msg[MID_MESSAGE_SIZE];

	int sent = 0;
	for (i=0;i<depProclistCount;i++) {
		if (depProclist[i].sent == 0) {
			sprintf (msg, "Need to send cell %ld in part %ld to %d\n ", depProclist[i].cellIndex, depProclist[i].partIndex, depProclist[i].depProc);
			mprintf(3, msg, 3);
		
			for (i=0;i<OCout_ub;i++) {
				if (depProclist[i].cellIndex == OCout[j].cellIndex) {
					sprintf (msg, "found cell %ld in part %ldto %d\n ", depProclist[i].cellIndex, depProclist[i].partIndex, depProclist[i].depProc);
					mprintf (3, msg, 3);
					sendOC (depProclist[i].depProc,depProclist[i].cellIndex,OCout[j].cellScore);
					sent = 1;
					depProclist[i].sent = 1;
				}
			}
		}
	}
	return sent;
}

void * SRecvThread (void * threadarg) {
	int i, ibuf, flag;
	float fbuf;
	char cbuf;
	MPI_Status status;
  Slave * slave;
  int errCode;
	int DepProcflag, OCFlag;
	int done = 0;
	char msg[MID_MESSAGE_SIZE];

  slave = (Slave *) threadarg;
  depProclistCount = 0;
  DepProcflag = OCFlag = 1;
  slave->ComputationPhase = 1;
	while ((slave->ComputationPhase == 1) || (DepProcflag == 1) || (slave->compThreadFinished == 0)) {
		flag = 0; 
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
		while (flag == 0) {
			MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, MOAMSA_COMM_WORLD, &flag, &status);
			if (slave->compThreadFinished != 0) {
				break;
			}
		}
		if (slave->compThreadFinished != 0) {
			break;
		}
		sprintf (msg, "Received MPI_ANY_TAG %d CompPhase %d DepProcFlag %d CompThreadFinished %d\n", status.MPI_TAG, slave->ComputationPhase, DepProcflag, slave->compThreadFinished);
    mprintf(3, msg, 2);
		switch  (status.MPI_TAG) {
			case MPI_TAG_CompPhase:
			  RecvCompPhase(&slave->ComputationPhase);
				sprintf (msg, "Received CompPhase %d\n", slave->ComputationPhase);
    		mprintf(3, msg, 2);
  			if (slave->ComputationPhase == 1)  {
		  		ReceivePartition (slave);  				
				  mprintf (10, "partitioned sequences: \n", 2);
				  for (i=0; i<slave->MOAPart[slave->partitionsCount-1].msaAlgn->dimn; i++) {
				  	sprintf (msg, "%s\n", slave->MOAPart[slave->partitionsCount-1].sequences[i]);
						mprintf (10, msg, 2);
					}
					fflush(stdout);
				}
				else {/* to return to master the last time*/
					sem_post (&slave->cpSem);
					//sem_post (&slave->retSem);
				}
			break;
			case MPI_TAG_DPFlag:  
		    DepProcflag = receiveCellDepProc (slave, &depProclistCount, &depProclist);  
			  if (DepProcflag == 1) {
	    		sprintf (msg, "received cell %ld in part %ld to %d \n", depProclist[depProclistCount - 1].cellIndex, depProclist[depProclistCount - 1].partIndex, depProclist[depProclistCount - 1].depProc);
	    		mprintf(3, msg, 2);
	    		trySendDPOC (&depProclist[depProclistCount - 1]);
					if (depProclist[depProclistCount - 1].sent == 1)
						sprintf (msg, "sent %ld \n", depProclist[depProclistCount - 1].cellIndex);
					else
						sprintf (msg, "did not send %ld \n", depProclist[depProclistCount - 1].cellIndex);
		    		mprintf(3, msg, 2);
 	  		}
 	  		else {
 	  			mprintf (3, "received the end for DP\n", 2);
 	  		}
 	  			
			break;

			case MPI_TAG_OCFlag :    
		    OCFlag = receiveOC (slave, &OCin_ub, &OCin);
			  if (OCFlag == 1) {
	    		sprintf (msg, "received IC  %ld with score %ld \n", OCin[OCin_ub - 1].cellIndex, OCin[OCin_ub - 1].cellScore);
					mprintf (3, msg, 2);
					sem_post (&slave->ocSem);
 	  		}
			break;
		}
	}
return NULL;
}

void SlaveProcess (long seqNum, long * seqLen, int stype) {
  pthread_t SendProc, SCompProc, RecvProc;

  myMasterid = 0;
  Slave * slave;
	
  initSlaveMemory (&slave, seqNum, seqLen, stype);
  // 1. Create Slave threads
	fflush(stdout);
  pthread_create (&SendProc, NULL, SSendThread, slave);
  pthread_create (&SCompProc, NULL, SCompThread, slave);
  pthread_create (&RecvProc, NULL, SRecvThread, slave);

  // Join Master Threads
  pthread_join (SendProc, NULL);
  mprintf(1, "Joined Sending Thread 1\n", 1);
  pthread_join (SCompProc, NULL);
  mprintf (1, "Joined Score Computation Thread\n", 1);
  pthread_join (RecvProc, NULL);
  mprintf (1, "Joined Receiving Thread\n", 1);

	checkunSendOC ();

	checkPointSlave (slave);
	freeSlaveMemory (slave);
  initSlaveMemory (&slave, seqNum, seqLen, stype);
  restoreSlaveCheckPoint (slave);

	mprintf (1, "read Slave data from checkpoint file\n", 1);

	freeSlaveMemory (slave);

  mprintf(1, "joined threads & exiting\n", 1);
}

/*
int main (int argc, char ** argv) {
  MPI_Status status;
  int stype = 0, MPI_return;
  char * * sequences = NULL, msg[500];
  long * seqLen = NULL;
  long seqNum, i, j;
	char ufilename[SHORT_MESSAGE_SIZE];


  MPI_Group orig_group;
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
  MPI_return = MPI_Recv ( &stype, 1, MPI_TYPE_DSSType, 0, MPI_TAG_DSSType, MOAMSA_COMM_WORLD, &status);

  checkMPIErrorCode (MPI_return);
	sprintf(msg, "stype %d\n", stype);
	mprintf(3, msg, 1);

  MPI_return = MPI_Recv (&seqNum, 1, MPI_TYPE_DSDimn, 0, MPI_TAG_DSDimn, MOAMSA_COMM_WORLD, &status);
  checkMPIErrorCode (MPI_return);
	sprintf(msg, "seqNum %d\n", seqNum);
	mprintf(3, msg, 1);
  seqLen =  (long *) mmalloc (sizeof(long) * seqNum);
  for (j=0;j<seqNum;j++) {
		MPI_return = MPI_Recv (&seqLen[j], 1, MPI_TYPE_DSShape, 0, MPI_TAG_DSShape, MOAMSA_COMM_WORLD, &status);
  	checkMPIErrorCode (MPI_return);
		sprintf(msg, "seqLen %d \n", seqLen[j]);
		mprintf(3, msg, 1);
	}
  
	*
		sprintf(msg, " ClusterSize%d \n", ClusterSize);
		mprintf(1, msg, 1);
  //  SlaveProcess (seqNum, seqLen, stype);
	
	mprintf(1, "before finalize \n", 1);
  MPI_Finalize ();
	mprintf(1, "after finalize \n", 1);
  return 0;
}

*/
