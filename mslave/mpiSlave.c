#include <stdio.h>
#include<string.h>
#include <stdlib.h>
#include <mpi.h>
#include "moa.h"
#include "moaDst.h"
#include "globals.h"
#include "utils.h"
#include "slave.h"
#include <pthread.h>

int ReceivePartition (Slave  * slave) {
	int errCode;
	char msg[SHORT_MESSAGE_SIZE];

  mprintf (20, "in ReceivePartition \n", threadnum);
	//errCode = pthread_mutex_lock (slave->CPmut);
  checkMutexErrorCode (1, errCode);
  slave->partitionsCount ++;

  if (slave->partitionsCount == 1) 
  	slave->MOAPart = (MOAPartition *) mmalloc (sizeof(MOAPartition));    
  else 
   	slave->MOAPart = (MOAPartition *) realloc (slave->MOAPart, sizeof(MOAPartition) * (slave->partitionsCount));
    
  /* receive MOA partition indices*/
  mprintf (20, "will receive MOA partition \n", threadnum);
  receiveMOAIndices  (&slave->MOAPart[slave->partitionsCount-1]);
	//errCode = pthread_mutex_unlock (slave->CPmut);
 	checkMutexErrorCode (1, errCode);

	sprintf (msg, "Received Partition %ld starting @ %ld\n", slave->partitionsCount, slave->MOAPart[slave->partitionsCount-1].msaAlgn->indexes[0]);
  mprintf(1, msg, threadnum);
	slave->MOAPart[slave->partitionsCount-1].processed = 0;
  sem_post (&slave->cpSem);
return 0;
}

int receiveMOAIndices (MOAPartition  * MOAPart) {
  int src;
  int msgsize;
  MPI_Status status;
  MPI_Request request;
  char * buffer, cmessage, msg[SHORT_MESSAGE_SIZE];
  unsigned int	membersize, maxsize;
  int dest, tag;
  int position;
  long i, j, lmessage;
  unsigned long ulmessage;
  long * ind;
  int MPI_return;
  
  createMOAStruct (&MOAPart->msaAlgn);
  /*3) Wave No */
  MPI_return = NBReceive(&MOAPart->waveNo, 1, MPI_TYPE_PWave, 0, MPI_TAG_PWave, &request, &status);
  checkMPIErrorCode (2, MPI_return);

  /*4) dimn */
  MPI_return = NBReceive(&MOAPart->msaAlgn->dimn, 1, MPI_TYPE_PDimn, 0, MPI_TAG_PDimn, &request, &status);
  checkMPIErrorCode (2, MPI_return);
 
		
    //printf("Slave %d received dimn %d\nSlave %d received shape = ", myProcid, MOAPart->msaAlgn->dimn, myProcid);
		//fflush(stdout);
  /*5) shape */
  MOAPart->msaAlgn->shape = (long *) mcalloc(MOAPart->msaAlgn->dimn, sizeof(long));
  for (i=0;i<MOAPart->msaAlgn->dimn;i++) {
    MPI_return = NBReceive(&lmessage, 1, MPI_TYPE_PShape, 0, MPI_TAG_PShape, &request, &status);
  checkMPIErrorCode (2, MPI_return);
    MOAPart->msaAlgn->shape[i] = lmessage;
  }
  //6) elements_ub  
  MPI_return = NBReceive(&MOAPart->msaAlgn->elements_ub, 1, MPI_TYPE_PElmUB, 0, MPI_TAG_PElmUB, &request, &status);
  checkMPIErrorCode (2, MPI_return);
  sprintf(msg, "received elements_ub %ld\n", MOAPart->msaAlgn->elements_ub);
  mprintf(12, msg, threadnum);
  //7) indexes 
  MOAPart->msaAlgn->elements = (MOA_elm * ) mcalloc (MOAPart->msaAlgn->elements_ub, sizeof(MOA_elm));
  MOAPart->msaAlgn->indexes = (unsigned long * ) mcalloc (MOAPart->msaAlgn->elements_ub, sizeof(unsigned long));
  ind =  (long *) mmalloc (sizeof(long) * MOAPart->msaAlgn->dimn);
  
  
  for (i=0;i<MOAPart->msaAlgn->elements_ub;i++) {
    MPI_return = NBReceive(&ulmessage, 1, MPI_TYPE_PIndices, 0, MPI_TAG_PIndices, &request, &status);
    checkMPIErrorCode (2, MPI_return);
    MOAPart->msaAlgn->indexes[i] = ulmessage;
    MOAPart->msaAlgn->elements[i].val = 0;
  }
  sprintf(msg, "received indices[%ld] = %ld\n", i-1, MOAPart->msaAlgn->indexes[i-1]);
	mprintf(12, msg, threadnum);
    // 8) Receive Corresponding Residues from the Sequences read
  MOAPart->sequences = (char * * ) mmalloc (sizeof (char *) * MOAPart->msaAlgn->dimn);
  
  for (j=0;j<MOAPart->msaAlgn->dimn;j++) {
  }
  for (i=0;i< MOAPart->msaAlgn->dimn;i++) {
	    MOAPart->sequences[i] = (char *) mmalloc (sizeof (char) * (MOAPart->msaAlgn->shape[i] + 1));
      MPI_return = NBReceive(&cmessage, 1, MPI_TYPE_PResidue, 0, MPI_TAG_PResidue, &request, &status);
      checkMPIErrorCode (2, MPI_return);
      MOAPart->sequences[i][0] = cmessage;
      sprintf(msg, "received partition characters MOAPart->sequences[%ld][0] %c\n", i, MOAPart->sequences[i][0]);
			mprintf(12, msg, threadnum);
    }
      
  for (i=0;i< MOAPart->msaAlgn->dimn;i++) {
    for (j=1;j<MOAPart->msaAlgn->shape[i];j++) {
      MPI_return = NBReceive(&cmessage, 1, MPI_TYPE_PResidue, 0, MPI_TAG_PResidue, &request, &status);
      checkMPIErrorCode (2, MPI_return);
      MOAPart->sequences[i][j] = cmessage;
      sprintf(msg, "received partition characters MOAPart->sequences[%ld][%ld] %c\n", i, j, MOAPart->sequences[i][j]);
			mprintf(12, msg, threadnum);
    }
  }
  for (j=0;j<MOAPart->msaAlgn->dimn;j++) 
	  MOAPart->sequences[j][MOAPart->msaAlgn->shape[j]] = '\0';
  free (ind);
  return 0;
}

int SendSlavetoProcGate (int ComputationPhase ) {
  MPI_Request request; /*Non-Blocking send and receive status*/
	char msg[SHORT_MESSAGE_SIZE];
  int MPI_return;
  sprintf (msg, "Will Return to ProcGate with CompPhase %d\n", ComputationPhase);
	mprintf (9, msg, threadnum);
  MPI_return = NBSend (&ComputationPhase, 1, MPI_TYPE_ProcEnq, 0, MPI_TAG_ProcEnq, &request);
  checkMPIErrorCode (3, MPI_return);
  mprintf (9, "Returned to ProcGate\n", threadnum);
  return 0;
}

int receiveCellDepProc (Slave * slave, long * depProclistCount, depProc_rec * * depProclist) {
  int MPI_return, flag, DepProcflag;
  MPI_Request request;/*Non-blocking send or receive status*/
  MPI_Status status; /*Message receive strcucture*/
  struct timeval timeout;
	char msg[SHORT_MESSAGE_SIZE];
  long waveNo;
  unsigned long cellIndex, partIndex; /* Intermediate Cell Index */
  int depProc; /* dep Proc */

  MPI_return = NBReceive (&DepProcflag, 1, MPI_TYPE_DPFlag, 0, MPI_TAG_DPFlag, &request, &status);
  sprintf (msg, "received DepProcflag = %d\n", DepProcflag);
	mprintf (10, msg, threadnum);
  if (DepProcflag == 1) { // Still sending, otherwise stop thread

    MPI_return = NBReceive (&cellIndex, 1, MPI_TYPE_DPCellIndex, 0, MPI_TAG_DPCellIndex, &request, &status);
    
    MPI_return = NBReceive (&partIndex, 1, MPI_TYPE_DPPartIndex, 0, MPI_TAG_DPPartIndex, &request, &status);

    MPI_return = NBReceive (&waveNo, 1, MPI_TYPE_DPWave, 0, MPI_TAG_DPWave, &request, &status);

    MPI_return = NBReceive (&depProc, 1, MPI_TYPE_DPProcID, 0, MPI_TAG_DPProcID, &request, &status);
		pthread_mutex_lock (slave->DPmut);
		(*depProclistCount) ++;
    if ((*depProclistCount) == 1) 
			(*depProclist) = (depProc_rec *) mmalloc (sizeof(depProc_rec));	
   	else 
			(*depProclist) = (depProc_rec *) realloc ((*depProclist), sizeof(depProc_rec) * ((*depProclistCount)));
    (*depProclist)[(*depProclistCount) - 1].depProc = depProc;
    (*depProclist)[(*depProclistCount) - 1].cellIndex = cellIndex;
    (*depProclist)[(*depProclistCount) - 1].partIndex = partIndex;
    (*depProclist)[(*depProclistCount) - 1].waveNo = waveNo;
    (*depProclist)[(*depProclistCount) - 1].sent = 0;
    pthread_mutex_unlock (slave->DPmut);
		sem_post (&slave->dpSem);
  }
  
  return DepProcflag;
}

void sendOC (int toProc, unsigned long cellIndex, long cellScore) {
  MPI_Request request;/*Non-blocking send or receive status*/
  int MPI_return, ICflag;
	char msg[SHORT_MESSAGE_SIZE];

	  sprintf (msg, "in sendOC cellIndex %ld cellScore %ld toProc %d\n", cellIndex, cellScore, toProc);
		mprintf (12, msg, threadnum);
	if (toProc != myProcid) {

		mprintf (12, "after if\n", threadnum);
		ICflag = 1;
	  MPI_return = NBSend (&ICflag, 1, MPI_TYPE_OCFlag, toProc, MPI_TAG_OCFlag, &request);
  	checkMPIErrorCode (3, MPI_return);
		
		sprintf (msg, "send flag and return %d\n", MPI_return);
		mprintf (12, msg, threadnum);

	  MPI_return = NBSend (&cellIndex, 1, MPI_TYPE_OCCellIndex, toProc, MPI_TAG_OCCellIndex, &request);
	  checkMPIErrorCode (3, MPI_return);

		sprintf (msg, "send cellIndex and return %d\n", MPI_return);
		mprintf (12, msg, threadnum);
  
  	MPI_return = NBSend (&cellScore, 1, MPI_TYPE_OCCellScore, toProc, MPI_TAG_OCCellScore, &request);
	  checkMPIErrorCode (3, MPI_return);

	  sprintf (msg, "sent cell %ld score %ld to %d MPI_return %d\n", cellIndex, cellScore, toProc, MPI_return);
		mprintf (12, msg, threadnum);
	}
}

int receiveOC (Slave * slave, long * OCin_ub, OCIType * * OCin){
  MPI_Status status; /*Message receive strcucture*/
  MPI_Request request;/*Non-blocking send or receive status*/
	char msg[SHORT_MESSAGE_SIZE];
  int MPI_return;
	unsigned long cellIndex;
	long cellScore;
  int flag, OCFlag; /* dep Proc */
  struct timeval timeout;  

  MPI_return = NBReceive (&OCFlag, 1, MPI_TYPE_OCFlag, MPI_ANY_SOURCE, MPI_TAG_OCFlag, &request, &status);

  sprintf (msg, "received OCFlag = %d\n", OCFlag);
	mprintf (10, msg, threadnum);
  if (OCFlag == 1) { // Still sending, otherwise stop thread

  	MPI_return = NBReceive (&cellIndex, 1, MPI_TYPE_OCCellIndex, status.MPI_SOURCE, MPI_TAG_OCCellIndex, &request, &status); 
  
  	MPI_return = NBReceive (&cellScore, 1, MPI_TYPE_OCCellScore, status.MPI_SOURCE, MPI_TAG_OCCellScore, &request, &status);

		(*OCin_ub) ++;
    if ((*OCin_ub) == 1) 
			(*OCin) = (OCIType *) mmalloc (sizeof(OCIType));	
   	else 
			(*OCin) = (OCIType *) realloc ((*OCin), sizeof(OCIType) * ((*OCin_ub)));
    (*OCin)[(*OCin_ub) - 1].cellIndex = cellIndex;
    (*OCin)[(*OCin_ub) - 1].cellScore = cellScore;
		(*OCin)[(*OCin_ub) - 1].fromProc = status.MPI_SOURCE;
	}
  return OCFlag;
}

int RecvCompPhase(int * ComputationPhase) {
  int i, MPI_return;
	char msg[SHORT_MESSAGE_SIZE];
	MPI_Request request;
  MPI_Status status;


  MPI_return = NBReceive (&(*ComputationPhase), 1, MPI_TYPE_CompPhase, myMasterid, MPI_TAG_CompPhase, &request, &status);
  checkMPIErrorCode (2, MPI_return);
  sprintf (msg, "Received Computationphase %d\n", (*ComputationPhase));
	mprintf (12, msg, threadnum);
  return MPI_return;
}
