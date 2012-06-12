#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "slave.h"
#include "globals.h"
#include "moa.h"
#include "utils.h"
#include "moaDst.h"

int myProcid, ClusterSize;
MPI_Comm MOAMSA_COMM_WORLD;


long getPartIndex (void * arg, long cellIndex, long * localCellIndex, long startPart) {
	long i, j;
	char msg[MID_MESSAGE_SIZE];
	Slave * slave;
	slave = (Slave *) arg;
	(*localCellIndex) = -1;
	for (i=startPart;i<slave->partitionsCount;i++) {
		for (j=0;j<slave->MOAPart[i].msaAlgn->elements_ub;j++) {
			if (cellIndex == slave->MOAPart[i].msaAlgn->indexes[j]) {
				(*localCellIndex) = j;
				return i;
			}
		}		
	}	
	return -1;
}

long getTraceBackPartition (void * arg, long cellIndex, long * localCellIndex) {
	long partIndex, maxpartIndex, maxlocalCellIndex, locCellIndex;
	char msg[MID_MESSAGE_SIZE];
	Slave * slave;
	slave = (Slave *) arg;
	
	partIndex = getPartIndex (slave, cellIndex, &locCellIndex, 0);
	maxlocalCellIndex = locCellIndex;
	maxpartIndex = partIndex;
	sprintf(msg, "in getTraceBackPartition locCellIndex %ld partIndex %ld\n", locCellIndex, partIndex);  
	mprintf(20, msg);
	while ((partIndex >= 0) && (partIndex < slave->partitionsCount-1)) {
		partIndex = getPartIndex (slave, cellIndex, &locCellIndex, partIndex + 1);
		sprintf(msg, "in getTraceBackPartition locCellIndex %ld partIndex %ld \n", locCellIndex, partIndex);  
		mprintf(20, msg);
		if ((locCellIndex > maxlocalCellIndex) && (partIndex != -1) && (partIndex < slave->partitionsCount-1)) {
			maxlocalCellIndex = locCellIndex;
			maxpartIndex = partIndex;
		}
	}
	(*localCellIndex) = maxlocalCellIndex;
	sprintf(msg, "leaving getTraceBackPartition with (*localCellIndex) %ld maxpartIndex %ld\n", (*localCellIndex), maxpartIndex);  
	mprintf(20, msg);
	return maxpartIndex;
}

long getMaxScoreinLowerNeighbors (long * ind, MOA_rec * MOA_in) {
	MOA_rec * rslt;
	int ret;
	long i, maxIndex, maxScore;
  char msg[MID_MESSAGE_SIZE];

	mprintf(20, "in getMaxScoreinLowerNeighbors\n");
	ret = MOAGetLowerNeighbors (ind, MOA_in, rslt);
	if (ret < 0)
		return -1;
	sprintf(msg, "rslt->elements_ub %ld\n", rslt->elements_ub);  
	mprintf(20, "in getMaxScoreinLowerNeighbors\n");
	maxScore = 	rslt->elements[0].val;
	maxIndex = 0;
	for (i=1;i<rslt->elements_ub;i++) {
		if (rslt->elements[i].val > maxScore) {
			maxScore = rslt->elements[i].val;
			maxIndex = i;
		}
	}
	return maxIndex;
}

int checkPrevInLocalPartitions (void * arg, long cellIndex, long * locCellIndex, long * partIndex) {
	long i, j, currPartIndex = (*partIndex);
	Slave * slave;
  char msg[MID_MESSAGE_SIZE];

	slave = (Slave *) arg;
	mprintf(20, "in checkPrevInLocalPartitions\n");

	for (i=0;i<currPartIndex;i++) {
		for (j=0;j<slave->MOAPart[i].msaAlgn->elements_ub;j++) {
			if ((cellIndex == slave->MOAPart[i].msaAlgn->indexes[j]) && (slave->MOAPart[i].msaAlgn->elements[j].prev_ub > 0)) {
				(*partIndex) = i;
				(*locCellIndex) = j;
				return i;
			}
		}
	}
	return -1;
}
	/* check in Received Cells from Other Processors */
long checkPrevInRemotePartitions (long cellIndex, long * partIndex) {
	long i;
  char msg[MID_MESSAGE_SIZE];

	mprintf(20, "in checkPrevInRemotePartitions\n");
	for (i=0;i<OCin_ub;i++) {
		if (cellIndex == OCin[i].cellIndex) {
				(*partIndex) = i;
				sprintf(msg, "Found remote partition in proc %d OC[%ld] with global index %ld\n", OCin[i].fromProc, i, OCin[i].cellIndex);  
				mprintf(3, msg);
				return i;
		}		
	}	
	return -1;
}

long traceBack_one (void * arg, char * * sequences, long startIndex, char * * * algnseq, long * endCellIndex, int * nextProc) {
  long currentScore, currentCell, i, j, aSeqLen = 0;
  long partIndex, hLNCellIndex, prevScore, prevCell = 0;
  long * old_index = NULL; /* multidimensional index*/
  long * new_index = NULL; /* multidimensional index of new cell*/
  char temp, msg[MID_MESSAGE_SIZE];
	Slave * slave;
	slave = (Slave *) arg;
	MOA_rec * msaAlgn;

  aSeqLen = 0;
	(*nextProc) = -1;
  /* Get Max Cell on Last Border as Current Cell*/
  //currentScore = getMaxOnLastBorder (msaAlgn, &currentCell);
	partIndex = getTraceBackPartition (slave, startIndex, &currentCell);
	msaAlgn = slave->MOAPart[partIndex].msaAlgn;
	currentScore = msaAlgn->elements[currentCell].val;
  (*endCellIndex) = msaAlgn->indexes[currentCell];

  old_index = (long *)  mcalloc ((msaAlgn->dimn), sizeof(long));
  new_index = (long *)  mcalloc ((msaAlgn->dimn), sizeof(long));
		/*decide whether local or global address need to be used here*/
	  //Gamma_Inverse(currentCell, msaAlgn->shape, msaAlgn->dimn, new_index);
  Gamma_Inverse(msaAlgn->indexes[currentCell], slave->seqLen, slave->seqNum, new_index);

  /*
  mprintf(6, "All Sequences:\n");
  for (i=0;i<slave->seqNum;i++) {
  	old_index[i] = new_index[i];
  	for (j=0;j<slave->seqLen[i];j++) {
			sprintf(msg, " %c ", sequences[i][j]);
			mprintf(6, msg);
		}
			mprintf(6, "\n");
	}
  */
  
  /* Iterate making the Neighbor Cell that Computed the Current Cell, the new current cell till the origin is reached*/
  while (currentCell > 0) {
    prevScore = currentScore; 
    prevCell = currentCell; 

    if ((msaAlgn->elements[prevCell].prev_ub > 0) && (msaAlgn->elements[prevCell].prev != NULL)) {
      currentCell = msaAlgn->elements[prevCell].prev[0];
    	currentScore = msaAlgn->elements[currentCell].val;
	  	(*endCellIndex) = msaAlgn->indexes[currentCell];
			/*decide whether local or global address need to be used here*/
		  //Gamma_Inverse(currentCell, msaAlgn->shape, msaAlgn->dimn, new_index);
	  	Gamma_Inverse(msaAlgn->indexes[currentCell], slave->seqLen, slave->seqNum, new_index);
		}
		/*else ifLowerBordercell and not lowerin whole tensor, continue getting prev_ub in previous partitions - means don't send prev data in OC, otherwise, send computing processor id*/
	/* check in prevPartition in this Processor*/

		else if (checkPrevInLocalPartitions (slave, msaAlgn->indexes[currentCell], &currentCell, &partIndex) > 0) {
			msaAlgn = slave->MOAPart[partIndex].msaAlgn;
			sprintf (msg, "Found in local partition %ld currentCell %ld with prev_ub %ld\n", partIndex, currentCell, msaAlgn->elements[currentCell].prev_ub);
			mprintf(3, msg);
			currentCell = msaAlgn->elements[currentCell].prev[0];
			sprintf (msg, "CurrentCell %ld \n", currentCell);
			mprintf(3, msg);
    	currentScore = msaAlgn->elements[currentCell].val;
	  	(*endCellIndex) = msaAlgn->indexes[currentCell];
			/* decide whether local or global address need to be used here */
		  //Gamma_Inverse(currentCell, msaAlgn->shape, msaAlgn->dimn, new_index);
	  	Gamma_Inverse(msaAlgn->indexes[currentCell], slave->seqLen, slave->seqNum, new_index);
		}
		/* check in Received Cells from Other Processors */
		else if (checkPrevInRemotePartitions (msaAlgn->indexes[currentCell], &partIndex) >= 0) {
			(*nextProc) = OCin[partIndex].fromProc;
			sprintf(msg, "Found in remote partition in proc %d\n", (*nextProc));
			mprintf(3, msg);
		    	currentCell = OCin[partIndex].cellIndex;
		    	currentScore = OCin[partIndex].cellScore;
		  	(*endCellIndex) = OCin[partIndex].cellIndex;
			/*decide whether local or global address need to be used here*/
			//Gamma_Inverse(currentCell, msaAlgn->shape, msaAlgn->dimn, new_index);
		  	Gamma_Inverse(OCin[partIndex].cellIndex, slave->seqLen, slave->seqNum, new_index);
		}
	/*
		else if (hLNCellIndex = getMaxScoreinLowerNeighbors (new_index, msaAlgn) > 0)
			currentCell = hLNCellIndex;
*/
    else {
	currentCell = 0;
    	currentScore = msaAlgn->elements[currentCell].val;
	(*endCellIndex) = msaAlgn->indexes[currentCell];
	/*decide whether local or global address need to be used here*/
	//Gamma_Inverse(currentCell, msaAlgn->shape, msaAlgn->dimn, new_index);
	Gamma_Inverse(msaAlgn->indexes[currentCell], slave->seqLen, slave->seqNum, new_index);
	}
	/* only get aligned residue, in case there was a move in the current partition */
	sprintf (msg, "proc %d\n", (*nextProc));
	mprintf (20, msg);
	if ((*nextProc) == -1) {
		sprintf(msg, "prevCell = %ld prevScore = %ld currentCell = %ld currentScore = %ld prev_ub %ld\n", prevCell, prevScore, currentCell, currentScore, msaAlgn->elements[prevCell].prev_ub);
		mprintf(3, msg);
	 	for (i=0;i<slave->seqNum;i++) {
  	    		if (aSeqLen == 0)
				(*algnseq)[i] = (char *) mmalloc (sizeof(char) * (aSeqLen + 1));
	      		else
				(*algnseq)[i] = (char *) realloc ((*algnseq)[i], sizeof(char) * (aSeqLen+1));
     
			if ((*algnseq)[i] == NULL ) {
				sprintf(msg, "Can not allocate memory to Alligned Sequence %ld\n", i);
				mprintf(1, msg);
				return -1;
		    	  }
	      		if ((old_index[i] > new_index[i])) /* && (isalpha( sequences[i][old_index[i] - 1]))) */
				(*algnseq)[i][aSeqLen] = sequences[i][old_index[i]];
		    	else
				(*algnseq)[i][aSeqLen] = GAPCHAR;
			sprintf(msg, "old_index[%ld] = %ld new_index = %ld is %c\n", i,  old_index[i], new_index[i],(*algnseq)[i][aSeqLen]);
			mprintf(20, msg);
	  	  	old_index[i] = new_index[i];
    		}
	    	aSeqLen ++;
	}
	else
	/* Otherwise, go to continue tracing in the other partition*/
	currentCell = 0;
  }

  sprintf(msg, "\n aSeqLen %ld ", aSeqLen);  
	mprintf(20, msg);
     
  /* Reverse the contents of the Aligned Sequences*/
  for (i=0;i<slave->seqNum;i++) {
    for (j=0;j<aSeqLen/2;j++) {
      temp = (*algnseq)[i][j];
      (*algnseq)[i][j] = (*algnseq)[i][aSeqLen - j - 1];
      (*algnseq)[i][aSeqLen - j - 1] = temp;
    }
  }
  if (new_index != NULL) 
    free(new_index);
  if (old_index != NULL) 
    free(old_index);
  return aSeqLen;
}
void traceBack (long seqNum, char * * sequences, long * seqLen, MOA_rec * msaAlgn, int stype, char * * * * algnseq, long * * aSeqLen, int * alignmentsNo, long * currentCell, long * currentScore, int Prev_index) {
  
  long i, j, k=0;
  long prevScore, prevCell=0;
  long * old_index = NULL; /* multidimensional index*/
  long * new_index = NULL; /* multidimensional index of new cell*/
  char temp, msg[100];
  int myAlign = (*alignmentsNo);
  int myPrevInd = Prev_index;

  old_index = (long *)  mcalloc ((msaAlgn->dimn), sizeof(long));
  new_index = (long *)  mcalloc ((msaAlgn->dimn), sizeof(long));

  Gamma_Inverse((*currentCell), msaAlgn->shape, msaAlgn->dimn, new_index);
  for (i=0;i<seqNum;i++) 
    old_index[i] = new_index[i];
  
  /* Iterate making the Neighbor Cell that Computed the Current Cell, the new current cell till the origin is reached*/
  while ((*currentCell) > 0) {
    prevScore = (*currentScore); 
    prevCell = (*currentCell); 
    if (myPrevInd != 0) {
      (*currentCell) = msaAlgn->elements[prevCell].prev[myPrevInd];
      (*currentScore) = msaAlgn->elements[(*currentCell)].val;
      myPrevInd = 0;
    }
    else if ((((*alignmentsNo) + 1) < maxAlignmentsNumber) && (msaAlgn->elements[prevCell].prev_ub > 1) && (msaAlgn->elements[prevCell].prev != NULL)) {
      for (i=1;i<msaAlgn->elements[prevCell].prev_ub;i++) {
				(*alignmentsNo) ++; 
				(*algnseq) = (char * * *) realloc ((*algnseq), ((*alignmentsNo) + 1) * sizeof(char * *));
		    if ((*algnseq) == NULL) {
					sprintf(msg, "Could not reallocate memory for New Aligned Sequences Set %ld!\n", alignmentsNo+1);
					mprintf(1, msg);
					return;
      	}
		    (*algnseq)[(*alignmentsNo)] = (char * *) mmalloc (seqNum * sizeof(char *));
				(*aSeqLen) = (long *)  realloc ((*aSeqLen), ((*alignmentsNo) + 1) * sizeof(long));
		    if ((*aSeqLen) == NULL) {
					sprintf(msg, "Could not reallocate memory for New Aligned Sequences Length %ld!\n", alignmentsNo+1);
					mprintf(1, msg);
					return;
      	}
				(*aSeqLen)[(*alignmentsNo)] = (*aSeqLen)[myAlign];
			  sprintf(msg, "\n copying %ld: ", (*aSeqLen)[myAlign]);
				mprintf(6, msg);
				for (j=0;j<seqNum;j++) {
				  (*algnseq)[(*alignmentsNo)][j] = (char *) mmalloc (sizeof(char) * ((*aSeqLen)[myAlign] + 1));
				  for (k=0;k<(*aSeqLen)[myAlign];k++) {
				    (*algnseq)[(*alignmentsNo)][j][k] = (*algnseq)[myAlign][j][k];
	    		  sprintf(msg, "%c ", (*algnseq)[myAlign][j][k]);
						mprintf(6, msg);
				  }
	  			mprintf (6, "\n");
				}
				sprintf(msg, "\n myAlign = %d : (*alignmentsNo) = %d, prev_ub: %ld prevCell = %ld prevScore = %ld currentCell = %ld currentScore = %ld ", myAlign, (*alignmentsNo), msaAlgn->elements[prevCell].prev_ub, prevCell, prevScore, (*currentCell), (*currentScore)  );
				mprintf(1, msg);

				traceBack (seqNum, sequences, seqLen, msaAlgn, stype, algnseq, aSeqLen, alignmentsNo, currentCell, currentScore, i);
      } 
      (*currentCell) = msaAlgn->elements[prevCell].prev[0];
      (*currentScore) = msaAlgn->elements[(*currentCell)].val;
			sprintf(msg, "\n myAlign = %d : (*alignmentsNo) = %d, prev_ub: %ld prevCell = %ld prevScore = %ld currentCell = %ld currentScore = %ld ", myAlign, (*alignmentsNo), msaAlgn->elements[prevCell].prev_ub, prevCell, prevScore, (*currentCell), (*currentScore)  );
			mprintf(6, msg);
    }
    else if ((msaAlgn->elements[prevCell].prev_ub > 0) && (msaAlgn->elements[prevCell].prev != NULL)) {
      (*currentCell) = msaAlgn->elements[prevCell].prev[0];
      (*currentScore) = msaAlgn->elements[(*currentCell)].val;
    }
    else {
      (*currentCell) = 0;
      (*currentScore) = msaAlgn->elements[(*currentCell)].val;
    }
    sprintf(msg, "\n myAlign = %d : prev_ub: %ld prevCell = %ld prevScore = %ld currentCell = %ld currentScore = %ld ", myAlign, msaAlgn->elements[prevCell].prev_ub, prevCell, prevScore, (*currentCell), (*currentScore)  );
		mprintf(6, msg);
    fflush(stdout);
    Gamma_Inverse((*currentCell), msaAlgn->shape, msaAlgn->dimn, new_index);
    for (i=0;i<seqNum;i++) {
      if ((*aSeqLen)[myAlign] == 0) {
				(*algnseq)[myAlign][i] = (char *) mmalloc (sizeof(char) * ((*aSeqLen)[myAlign] + 1));
      }
      else {
				(*algnseq)[myAlign][i] = (char *) realloc ((*algnseq)[myAlign][i], sizeof(char) * (((*aSeqLen)[myAlign])+1));
		    if ((*algnseq)[myAlign][i] == NULL) {
					sprintf(msg, "Could not reallocate memory for %d Aligned %ld Sequence at Character %ld!\n", myAlign, i+1, (((*aSeqLen)[myAlign])+1));
					mprintf(1, msg);
					return;
      	}
			}
      if ((*algnseq)[myAlign][i] == NULL ) {
				sprintf(msg, "Can not allocate memory to Alligned Sequence %ld\n", i);
				mprintf(1, msg);
				return;
      }
      if ((old_index[i] > new_index[i])) /* && (isalpha( sequences[i][old_index[i] - 1]))) */
				(*algnseq)[myAlign][i][(*aSeqLen)[myAlign]] = sequences[i][old_index[i] - 1];
      else
				(*algnseq)[myAlign][i][(*aSeqLen)[myAlign]] = GAPCHAR;
			sprintf(msg, "\n aSeqLen: %ld:  old_index[%ld] = %ld new_index = %ld is %c  ", (*aSeqLen)[myAlign], i,  old_index[i], new_index[i],(*algnseq)[myAlign][i][(*aSeqLen)[myAlign]] );
			mprintf(6, msg);
      old_index[i] = new_index[i];
    }
    (*aSeqLen)[myAlign] ++;
  }   
  /* Reverse the contents of the Aligned Sequences*/
  for (i=0;i<seqNum;i++) {
    for (j=0;j<(*aSeqLen)[myAlign]/2;j++) {
      temp = (*algnseq)[myAlign][i][j];
      (*algnseq)[myAlign][i][j] = (*algnseq)[myAlign][i][(*aSeqLen)[myAlign] - j - 1];
      (*algnseq)[myAlign][i][(*aSeqLen)[myAlign] - j - 1] = temp;
    }
  }
  if (new_index != NULL) 
    free(new_index);
  if (old_index != NULL) 
    free(old_index);
}
void traceBack_loc (long seqNum, char * * sequences, long * seqLen, MOA_rec * msaAlgn, int stype, char * * * * algnseq, long * * aSeqLen, int * alignmentsNo, long * currentCell, long * currentScore, int Prev_index) {
  
  long i, j, k=0;
  long prevScore, prevCell=0;
  long * old_index = NULL; /* multidimensional index*/
  long * new_index = NULL; /* multidimensional index of new cell*/
  char temp;
  int myAlign = (*alignmentsNo);
  int myPrevInd = Prev_index;
	char msg[100];

  old_index = (long *)  mcalloc ((msaAlgn->dimn), sizeof(long));
  new_index = (long *)  mcalloc ((msaAlgn->dimn), sizeof(long));

  Gamma_Inverse((*currentCell), msaAlgn->shape, msaAlgn->dimn, new_index);
  for (i=0;i<seqNum;i++) 
    old_index[i] = new_index[i];
  
  /* Iterate making the Neighbor Cell that Computed the Current Cell, the new current cell till the origin is reached*/
  while ((*currentScore) > 0) {
    prevScore = (*currentScore); 
    prevCell = (*currentCell); 
    if (myPrevInd != 0) {
      (*currentCell) = msaAlgn->elements[prevCell].prev[myPrevInd];
      (*currentScore) = msaAlgn->elements[(*currentCell)].val;
      myPrevInd = 0;
    }
    else if ((((*alignmentsNo) + 1) < maxAlignmentsNumber) && (msaAlgn->elements[prevCell].prev_ub > 1) && (msaAlgn->elements[prevCell].prev != NULL)) {
      for (i=1;i<msaAlgn->elements[prevCell].prev_ub;i++) {
				(*alignmentsNo) ++; 
				(*algnseq) = (char * * *) realloc ((*algnseq), ((*alignmentsNo) + 1) * sizeof(char * *));
     		if ((*algnseq) == NULL) {
					sprintf(msg, "Could not reallocate memory for Aligned Sequence Set Pointer %d!\n", (*alignmentsNo) + 1);
					mprintf(1, msg);
					return;
      	}
				(*algnseq)[(*alignmentsNo)] = (char * *) mmalloc (seqNum * sizeof(char *));
				(*aSeqLen) = (long *)  realloc ((*aSeqLen), ((*alignmentsNo) + 1) * sizeof(long));
     		if ((*aSeqLen) == NULL) {
					sprintf(msg, "Could not reallocate memory for Aligned Sequence Length %d!\n", (*alignmentsNo) + 1);
					mprintf(1, msg);
					return;
      	}
				(*aSeqLen)[(*alignmentsNo)] = (*aSeqLen)[myAlign];
			  sprintf(msg, "\n copying %ld: ", (*aSeqLen)[myAlign]);
				mprintf(6, msg);
				for (j=0;j<seqNum;j++) {
				  (*algnseq)[(*alignmentsNo)][j] = (char *) mmalloc (sizeof(char) * ((*aSeqLen)[myAlign] + 1));
				  for (k=0;k<(*aSeqLen)[myAlign];k++) {
				    (*algnseq)[(*alignmentsNo)][j][k] = (*algnseq)[myAlign][j][k];
	    		  sprintf(msg, "%c ", (*algnseq)[myAlign][j][k]);
						mprintf(1, msg);
				  }
				  mprintf (6, "\n");
				}
			  sprintf(msg, "\n myAlign = %d : (*alignmentsNo) = %d, prev_ub: %ld prevCell = %ld prevScore = %ld currentCell = %ld currentScore = %ld ", myAlign, (*alignmentsNo), msaAlgn->elements[prevCell].prev_ub, prevCell, prevScore, (*currentCell), (*currentScore)  );
				mprintf(6, msg);

				traceBack (seqNum, sequences, seqLen, msaAlgn, stype, algnseq, aSeqLen, alignmentsNo, currentCell, currentScore, i);
      } 
      (*currentCell) = msaAlgn->elements[prevCell].prev[0];
      (*currentScore) = msaAlgn->elements[(*currentCell)].val;
			sprintf(msg, "\n myAlign = %d : (*alignmentsNo) = %d, prev_ub: %ld prevCell = %ld prevScore = %ld currentCell = %ld currentScore = %ld ", myAlign, (*alignmentsNo), msaAlgn->elements[prevCell].prev_ub, prevCell, prevScore, (*currentCell), (*currentScore)  );
			mprintf(6, msg);
    }
    else if ((msaAlgn->elements[prevCell].prev_ub > 0) && (msaAlgn->elements[prevCell].prev != NULL)) {
      (*currentCell) = msaAlgn->elements[prevCell].prev[0];
      (*currentScore) = msaAlgn->elements[(*currentCell)].val;
    }
    else {
      (*currentCell) = 0;
      (*currentScore) = msaAlgn->elements[(*currentCell)].val;
    }
    sprintf(msg, "\n myAlign = %d : prev_ub: %ld prevCell = %ld prevScore = %ld currentCell = %ld currentScore = %ld ", myAlign, msaAlgn->elements[prevCell].prev_ub, prevCell, prevScore, (*currentCell), (*currentScore)  );
		mprintf(6, msg);
    Gamma_Inverse((*currentCell), msaAlgn->shape, msaAlgn->dimn, new_index);
    for (i=0;i<seqNum;i++) {
      if ((*aSeqLen)[myAlign] == 0) {
				(*algnseq)[myAlign][i] = (char *) mmalloc (sizeof(char) * ((*aSeqLen)[myAlign] + 1));
      }
      else {
				(*algnseq)[myAlign][i] = (char *) realloc ((*algnseq)[myAlign][i], sizeof(char) * (((*aSeqLen)[myAlign])+1));
     		if ((*algnseq)[myAlign][i] == NULL) {
					sprintf(msg, "Could not reallocate memory for %d Aligned %ld Sequence at Character %ld!\n", myAlign, i+1, (((*aSeqLen)[myAlign])+1));
					mprintf(1, msg);
					return;
      	}
    	}
    	if ((*algnseq)[myAlign][i] == NULL ) {
				sprintf(msg, "Can not allocate memory to Alligned Sequence %ld\n", i);
				mprintf(1, msg);
				return;
    	}
    	if ((old_index[i] > new_index[i])) /* && (isalpha( sequences[i][old_index[i] - 1]))) */ 
				(*algnseq)[myAlign][i][(*aSeqLen)[myAlign]] = sequences[i][old_index[i] - 1];
    	else
				(*algnseq)[myAlign][i][(*aSeqLen)[myAlign]] = GAPCHAR;
			sprintf(msg, "\n aSeqLen: %ld:  old_index[%ld] = %ld new_index = %ld is %c  ", (*aSeqLen)[myAlign], i,  old_index[i], new_index[i],(*algnseq)[myAlign][i][(*aSeqLen)[myAlign]] );
			mprintf(6, msg);
    	old_index[i] = new_index[i];
  	}
  	(*aSeqLen)[myAlign] ++;
  } /*End of While Loop*/  

  /* Reverse the contents of the Aligned Sequences */
  for (i=0;i<seqNum;i++) {
    for (j=0;j<(*aSeqLen)[myAlign]/2;j++) {
      temp = (*algnseq)[myAlign][i][j];
      (*algnseq)[myAlign][i][j] = (*algnseq)[myAlign][i][(*aSeqLen)[myAlign] - j - 1];
      (*algnseq)[myAlign][i][(*aSeqLen)[myAlign] - j - 1] = temp;
    }
  }
  if (new_index != NULL) 
    free(new_index);
  if (old_index != NULL) 
    free(old_index);
}

void distributedSlaveTraceBack (void * arg, char * * sequences) {
	Slave * slave;
	slave = (Slave *) arg;
	long maxScore, maxIndex, aSeqLen, endCellIndex, localCellIndex, maxlocalCellIndex, partIndex, i, j, * old_index = NULL;
	//long partPaths;
	int MPI_return, isLower, nextProc, done = 0;
	MPI_Request request;
	MPI_Status status;
	char * * /***/ algnseq, msg[MID_MESSAGE_SIZE];

	//partPaths = 1;
	while (done < 2) {
		/*receive tracing flag, 0: trace back, Otherwise: finish & exit*/
		MPI_return = NBReceive (&done, 1, MPI_INT, 0, 2, &request, &status);
		sprintf(msg, "DSTB received flag %d \n", done);  
		mprintf(3, msg);
		if (done == 0) {
			/*receive starting index*/
			MPI_return = NBReceive (&maxIndex, 1, MPI_LONG, 0, 3, &request, &status);
			sprintf(msg, "DSTB received maxIndex %ld \n", maxIndex);  
			mprintf(3, msg);
			/*Perform Partition trace back*/
			//if (partPaths == 1)
			//  	algnseq = (char * * *) mmalloc (sizeof(char * *));    
			//else
			//  	algnseq = (char * * *) mmalloc (partPaths * sizeof(char * *));    
		  	algnseq/*[partPaths-1]*/ = (char * *) mmalloc (slave->seqNum * sizeof(char *));    
			aSeqLen = traceBack_one (slave, sequences, maxIndex, &algnseq/*[partPaths-1]*/, &endCellIndex, &nextProc);
		  	sprintf(msg, "DSTB returned path of length %ld : ", aSeqLen);  
			mprintf(3, msg);
			/*send aligned sequences to Master*/		
			MPI_return = NBSend (&aSeqLen, 1, MPI_LONG, 0, 4, &request);
		  	for (i=0;i<slave->seqNum;i++) {
		    		for (j=0;j<aSeqLen;j++) {
					MPI_return = NBSend (&algnseq/*[partPaths-1]*/[i][j], 1, MPI_CHAR, 0, 5, &request);
				  	sprintf(msg, " %c ", algnseq/*[partPaths-1]*/[i][j]);  
					mprintf(6, msg);
				}
			}
			//partPaths ++;
			mprintf(3, "\n");
			/*send the last global index in this alignment to Master*/
			MPI_return = NBSend (&endCellIndex, 1, MPI_LONG, 0, 6, &request);
			MPI_return = NBSend (&nextProc, 1, MPI_INT, 0, 7, &request);
		}
		else if (done == 1) {
			MPI_return = NBReceive (&maxIndex, 1, MPI_LONG, 0, 8, &request, &status);
			sprintf(msg, "distributedSlaveTraceBack received maxIndex %ld \n", maxIndex);  
			mprintf(3, msg);
			partIndex = getTraceBackPartition (slave, maxIndex, &maxlocalCellIndex);
			sprintf(msg, "distributedSlaveTraceBack found %ld cell in %ld\n", localCellIndex, partIndex);  
			mprintf(3, msg);
			MPI_return = NBSend (&maxlocalCellIndex, 1, MPI_LONG, 0, 9, &request);
			sprintf(msg, "distributedSlaveTraceBack found %ld cell & informed Master\n", maxlocalCellIndex);  
			mprintf(3, msg);
			if (old_index != NULL)
				free (old_index);
			old_index = NULL;
		}
	}
}

void sendMaxScore (void * arg) {
	long i, j, maxScore, maxIndex, * wholeIndex;
	char msg[MID_MESSAGE_SIZE];
	MPI_Request request;
	int MPI_return;
	Slave * slave;
	slave = (Slave *) arg;
	maxScore = 0;
	maxIndex = 0;
	
	wholeIndex = NULL;
	wholeIndex = (long * ) mmalloc (sizeof(long) * slave->seqNum);
	for (i=0;i<slave->partitionsCount;i++) {
		sprintf (msg, "i %ld elm_ub %ld\n", i, slave->MOAPart[i].msaAlgn->elements_ub);
		mprintf(20, msg);
		for (j=0;j<slave->MOAPart[i].msaAlgn->elements_ub;j++) {
			sprintf (msg, "j %ld index %ld sqnum %ld sqlen %ld %ld %ld\n", j, slave->MOAPart[i].msaAlgn->indexes[j], slave->seqNum, slave->seqLen[0], slave->seqLen[1], slave->seqLen[2]);
			mprintf(20, msg);
			Gamma_Inverse (slave->MOAPart[i].msaAlgn->indexes[j], slave->seqLen, slave->seqNum, wholeIndex);
			if (isHigherBorderCell(wholeIndex, slave->seqNum, slave->seqLen) == 1)
				if (slave->MOAPart[i].msaAlgn->elements[j].val > maxScore) {
					maxScore = slave->MOAPart[i].msaAlgn->elements[j].val;
					maxIndex = slave->MOAPart[i].msaAlgn->indexes[j];
				}
			sprintf (msg, "testing index %ld score %ld\n", slave->MOAPart[i].msaAlgn->indexes[j], maxScore = slave->MOAPart[i].msaAlgn->elements[j].val);
			mprintf (20, msg);
		}	
	}
	MPI_return = NBSend (&maxScore, 1, MPI_LONG, 0, 0, &request);
	MPI_return = NBSend (&maxIndex, 1, MPI_LONG, 0, 1, &request);
	sprintf (msg, "DSTB maxScore%ld, maxIndex = %ld sent to Master\n", maxScore, maxIndex);
	mprintf (3, msg);
	if (wholeIndex != NULL)
		free (wholeIndex);
	wholeIndex = NULL;
}
/*
int main (int argc, char * argv[]) {
  MPI_Group orig_group;
  char * * sequences = NULL;
  long * seqLen = NULL;
  long seqNum, i, j;
	int stype, MPI_return, ret;
	char msg[MID_MESSAGE_SIZE];
	MPI_Status status;
	Slave * slave;
	char ufilename[SHORT_MESSAGE_SIZE];

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &myProcid);
  MPI_Comm_size (MPI_COMM_WORLD, &ClusterSize);
  MPI_Comm_group (MPI_COMM_WORLD, &orig_group);
  MPI_Comm_create(MPI_COMM_WORLD, orig_group, &MOAMSA_COMM_WORLD);

	/*1. Receive initial Data Size Info

  MPI_return = MPI_Recv ( &stype, 1, MPI_TYPE_DSSType, 0, MPI_TAG_DSSType, MOAMSA_COMM_WORLD, &status);

  checkMPIErrorCode (MPI_return);
	printf("Slave [%d] stype %d\n", myProcid, stype);
	fflush(stdout);

  MPI_return = MPI_Recv (&seqNum, 1, MPI_TYPE_DSDimn, 0, MPI_TAG_DSDimn, MOAMSA_COMM_WORLD, &status);
  checkMPIErrorCode (MPI_return);
	printf("Slave [%d] seqNum %d\n", myProcid, seqNum);
	fflush(stdout);
  seqLen =  (long *) mmalloc (sizeof(long) * seqNum);
  for (j=0;j<seqNum;j++) {
		MPI_return = MPI_Recv (&seqLen[j], 1, MPI_TYPE_DSShape, 0, MPI_TAG_DSShape, MOAMSA_COMM_WORLD, &status);
  	checkMPIErrorCode (MPI_return);
		printf("Slave [%d] seqLen %d \n", myProcid, seqLen[j]);
		fflush(stdout);
	}
	*
  // 1. Process Arguments
  processArguments(argc, argv, &seqNum, &sequences, &seqLen, &stype);
	strcpy (ufilename, outputfilename);
	sprintf (outputfilename, "mmtb_%s", ufilename);
  init_output();
	sprintf (msg, "Program Arguments: debuglevel = %d maxAlignmentsNumber = %d Epsilons= %d Alignment Type = %d stype %d outputfilename = %s partitionSize = %d\n", pdebug, maxAlignmentsNumber, Epsilons, AlignmentType, stype, outputfilename, partitionSize);
	mprintf (1, msg);

	//2. Load Tensor Partitions Computed
  initSlaveMemory (&slave, 0, NULL, 0);
	mprintf (1, "Initialized Slave data and will read checkpoint file\n");
  ret = restoreSlaveCheckPoint (slave);
	if (ret != 0)
		mprintf (1, "read Slave data from checkpoint file\n");
	else {
		sprintf (msg, "read Slave data from checkpoint file partitions %d last score in last partition is %d sqm %d sqlen0 %d %d %d\n", slave->partitionsCount, slave->MOAPart[slave->partitionsCount - 1].msaAlgn->elements[slave->MOAPart[slave->partitionsCount - 1].msaAlgn->elements_ub - 1].val, slave->seqNum, slave->seqLen[0], slave->seqLen[1], slave->seqLen[2]);
		mprintf(1, msg);
	}	
	//3. Synchronize with Master to trace back
	sendMaxScore (slave);
	distributedSlaveTraceBack (slave, sequences);
	freeSlaveMemory (slave);

	mprintf(1, "before finalize \n");
  MPI_Finalize ();
	mprintf(1, "after finalize \n");
}
*/
