#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include "slave.h"
#include "globals.h"
#include "utils.h"
#include "moaDst.h"
#define _GNU_SOURCE

void assemblePathParts (long pathParts, char * * * algnseq, long * aSeqLen, long seqNum) {
	long i, j, k, comseqLen;
	char  * * completePath, msg[MID_MESSAGE_SIZE];
	completePath = (char * *) mmalloc (seqNum * sizeof(char *));    
	comseqLen = 0;
	for (i=0;i<seqNum;i++) {
		comseqLen = 0;
		completePath[i] = NULL;
		for (j = pathParts-1;j>=0;j--) {
			completePath[i] = (char *) realloc (completePath[i], sizeof(char) * (comseqLen + aSeqLen[j]));
			sprintf(msg, "completePath[%ld] has %ld and will realloc %ld \n", i, comseqLen, (comseqLen + aSeqLen[j]));  
			mprintf(20, msg, threadnum);
			for (k = 0;k<aSeqLen[j];k++) {
				completePath[i][comseqLen] = algnseq[j][i][k];
				comseqLen ++;
			}
		}
	}
	PrintOptimalPath (seqNum, &completePath, comseqLen);
}

void distributedMasterTraceBack (long startCellScore, long startCellIndex, int startProc, long seqNum) {
	int MPI_return, currProc, done = 0;
	long maxCellScore, maxCellIndex, * aSeqLen, endCellIndex, localCellIndex, maxlocalCellIndex, i, j, pathParts;
	char * * * algnseq, msg[MID_MESSAGE_SIZE];
	MPI_Request request;
	MPI_Status status;
	maxCellScore = startCellScore;
	maxCellIndex = startCellIndex;
	currProc = startProc;
	pathParts = 0;
	algnseq = NULL;
	aSeqLen = NULL;
	while (done == 0) {
		/*send to processor containing the maxScore to trace back*/
		MPI_return = NBSend (&done, 1, MPI_LONG, currProc, 2, &request);
		sprintf(msg, "DMTB sent flag %d to proc %d first\n", done, currProc);  
		mprintf(3, msg, threadnum);
		MPI_return = NBSend (&maxCellIndex, 1, MPI_LONG, currProc, 3, &request);
		sprintf(msg, "DMTB sent maxCellIndex %ld to proc %d\n", maxCellIndex, currProc);  
		mprintf(3, msg, threadnum);
		/*Receive aligned sequences*/
		pathParts ++;
		algnseq = (char * * *) realloc (algnseq, pathParts * sizeof(char * *));
		aSeqLen = (long *) realloc (aSeqLen, pathParts * sizeof(long));
		algnseq[pathParts-1] = (char * *) mmalloc (seqNum * sizeof(char *));    
		MPI_return = NBReceive (&aSeqLen[pathParts-1], 1, MPI_LONG, currProc, 4, &request, &status);
		sprintf(msg, "DMTB received path of length %ld :", aSeqLen[pathParts-1]);  
		mprintf(3, msg, threadnum);
		for (i=0;i<seqNum;i++) {
			algnseq[pathParts-1][i] = (char *) mmalloc (aSeqLen[pathParts-1] * sizeof(char));    
			for (j=0;j<aSeqLen[pathParts-1];j++) {
				MPI_return = NBReceive (&algnseq[pathParts-1][i][j], 1, MPI_CHAR, currProc, 5, &request, &status);
				sprintf(msg, " %c ", algnseq[pathParts-1][i][j]);  
				mprintf(3, msg, threadnum);
			}
		}
		mprintf(3, "\n", threadnum);
		/*receive the last global index in this alignment*/
		MPI_return = NBReceive (&maxCellIndex, 1, MPI_LONG, currProc, 6, &request, &status);
		MPI_return = NBReceive (&currProc, 1, MPI_INT, currProc, 7, &request, &status);
		sprintf(msg, "DMTB received maxCellIndex %ld & currProc %d\n", maxCellIndex, currProc);  
		mprintf(3, msg, threadnum);
		/*test if end of aligned sequence is zero to exit */
		if (maxCellIndex == 0)
			done = 2;
		else if (currProc <= 0) { /* determine the next tracing Processor*/
			done = 1;
			maxlocalCellIndex = 0;
			for (i=1;i<ClusterSize;i++) {
				MPI_return = NBSend (&done, 1, MPI_INT, i, 2, &request);
				MPI_return = NBSend (&maxCellIndex, 1, MPI_LONG, i, 8, &request);
				sprintf(msg, "DMTB sent flag %d to proc %ld second\n", done, i);  
				mprintf(3, msg, threadnum);
				MPI_return = NBReceive (&localCellIndex, 1, MPI_INT, i, 9, &request, &status);
				sprintf(msg, "received local index %ld from %ld\n", localCellIndex, i);  
				mprintf(3, msg, threadnum);
				if (localCellIndex > maxlocalCellIndex) {
					maxlocalCellIndex = localCellIndex;
					currProc = i;
				}
			}
			sprintf(msg, "maxlocalIndex %ld \n", maxlocalCellIndex);  
			mprintf(3, msg, threadnum);
			done = 0;			
		}
	}
	/*send to all processes to exit tracing*/
	done = 2;
	for (i=1;i<ClusterSize;i++) {
		MPI_return = NBSend (&done, 1, MPI_LONG, i, 2, &request);
		sprintf(msg, "DMTB sent flag %d to proc %ld\n", done, i);  
		mprintf(3, msg, threadnum);
	}

	assemblePathParts (pathParts, algnseq, aSeqLen, seqNum);
}

void getMaxScore (long * maxCellScore, long * maxCellIndex, int * maxProc) {
	int i;
	long maxScore, maxIndex;
	char msg[MID_MESSAGE_SIZE];

	MPI_Status status;
	MPI_Request request;

	for (i=1;i<ClusterSize;i++) {
		NBReceive (&maxScore, 1, MPI_LONG, i, 0, &request, &status);
		NBReceive (&maxIndex, 1, MPI_LONG, i, 1, &request, &status);
		sprintf (msg, "Master received maxScore = %ld, maxIndex = %ld from %d\n", maxScore, maxIndex, i);
		mprintf (3, msg, threadnum);
		if (i == 1) {
			(*maxCellScore) = maxScore;
			(*maxCellIndex) = maxIndex;
			(*maxProc) = i;
		}
		else if (maxScore > (*maxCellScore)) {
			(*maxCellScore) = maxScore;
			(*maxCellIndex) = maxIndex;
			(*maxProc) = i;
		}
	}
	sprintf (msg, "Master maximum  = %ld, maximum Index = %ld from %d\n", (*maxCellScore), (*maxCellIndex), (*maxProc));
	mprintf (3, msg, threadnum);
	
}

int main (int argc, char * argv[]) {
  MPI_Group orig_group;
  char * * sequences = NULL, msg[MID_MESSAGE_SIZE];
  long * seqLen = NULL;
  long seqNum, i, j, maxCellIndex, maxCellScore;
  int stype, MPI_return, maxProc;
  char ufilename[SHORT_MESSAGE_SIZE];
  int ret;
  Slave * slave;
	
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &myProcid);
  MPI_Comm_size (MPI_COMM_WORLD, &ClusterSize);
  MPI_Comm_group (MPI_COMM_WORLD, &orig_group);
  MPI_Comm_create(MPI_COMM_WORLD, orig_group, &MOAMSA_COMM_WORLD);
  threadnum = 1;
  // 1. Process Arguments
  processArguments(argc, argv, &seqNum, &sequences, &seqLen, &stype);
  strcpy (ufilename, outputfilename);
  sprintf (outputfilename, "mmtb_%s", ufilename);
  init_output();
  sprintf (msg, "Program Arguments: debuglevel = %d maxAlignmentsNumber = %d Epsilons= %d Alignment Type = %d stype %d outputfilename = %s partitionSize = %d\n", pdebug, maxAlignmentsNumber, Epsilons, AlignmentType, stype, outputfilename, partitionSize);
  mprintf (1, msg, threadnum);
  
  //2. Synchronize with Slaves to trace back
  // one path for now
  
  if (myProcid == 0) {
    getMaxScore (&maxCellScore, &maxCellIndex, &maxProc);
    distributedMasterTraceBack (maxCellScore, maxCellIndex, maxProc, seqNum);
  }
  else {
    //2. Load Tensor Partitions Computed
    initSlaveMemory (&slave, 0, NULL, 0);
    mprintf (3, "Initialized Slave data and will read checkpoint file\n", threadnum);
    ret = restoreSlaveCheckPoint (slave);
    if (ret != 0)
      mprintf (3, " Could not read Slave data from checkpoint file\n", threadnum);
    else {
      sprintf (msg, "read Slave data from checkpoint file partitions %ld last score in last partition is %ld sqm %ld sqlen0 %ld %ld %ld\n", slave->partitionsCount, slave->MOAPart[slave->partitionsCount - 1].msaAlgn->elements[slave->MOAPart[slave->partitionsCount - 1].msaAlgn->elements_ub - 1].val, slave->seqNum, slave->seqLen[0], slave->seqLen[1], slave->seqLen[2]);
      mprintf(3, msg, threadnum);
    }	
    //3. Synchronize with Master to trace back
    sendMaxScore (slave);
    distributedSlaveTraceBack (slave, sequences);
    freeSlaveMemory (slave);
  }
/*
  aSeqLen = (long *)   mmalloc (sizeof(long));
  aSeqLen[0] = 0;
  prevCell = -1;
  algnseq = (char * * *) mmalloc(sizeof(char * *));
  algnseq[0] = (char * *) mmalloc (seqNum * sizeof(char *));    
  if (AlignmentType == Global) { // if Global Alignment 
    //PrintPrevChains(msaAlgn);
    // Get Max Cell on Last Border as Current Cell 
    alignmentsNo = 0;
    currentScore = getMaxOnLastBorder (msaAlgn, &currentCell);
    traceBack (seqNum, sequences, seqLen, msaAlgn, stype, &algnseq, &aSeqLen, &alignmentsNo, &currentCell, &currentScore, 0);
  }
  else { // if Local Alignment 
    alignmentsNo = -1;
    for (k = 0;k<maxAlignmentsNumber;k++) {
      if (k == 0)
	currentScore = MOA_max(msaAlgn, 0, 0, &currentCell);
      else {
	currentScore = MOA_max(msaAlgn, 1, currentScore, &currentCell);
      }
      if ((prevCell != currentCell) && (currentScore > 0)){
	alignmentsNo ++;
	algnseq = (char * * *) realloc (algnseq, (alignmentsNo + 1) * sizeof(char * *));
     if (algnseq == NULL) {
		sprintf(msg, "Could not reallocate memory for Aligned Sequence Set %d!\n", alignmentsNo + 1);
		mprintf(1, msg, threadnum);
		return -1;
      }
	algnseq[alignmentsNo] = (char * *) mmalloc (seqNum * sizeof(char *));
	aSeqLen = (long *)  realloc (aSeqLen, (alignmentsNo + 1) * sizeof(long));
     if (aSeqLen == NULL) {
		sprintf(msg, "Could not reallocate memory for Aligned Sequence Length %d!\n", alignmentsNo + 1);
		mprintf(1, msg, threadnum);
		return -1;
      }
	traceBack_loc (seqNum, sequences, seqLen, msaAlgn, stype, &algnseq, &aSeqLen, &alignmentsNo, &currentCell, &currentScore, 0);
	}
	prevCell = currentCell;
	prevScore = currentScore;
    }
  }

  // D. Print the resulting Alignemnts 
  PrintASeq (seqNum, sequences, seqLen, &algnseq, aSeqLen, alignmentsNo+1) ;
*/
	mprintf(1, "before finalize \n", threadnum);
  MPI_Finalize ();
	mprintf(1, "after finalize \n", threadnum);
	return 0;
}
