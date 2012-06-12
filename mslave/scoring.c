#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>      /* include the character library */
#include <pthread.h>
#include <mpi.h>
#include <errno.h>
#include "globals.h"
#include "utils.h"
#include "moaDst.h"
#include "moa.h"
#include "slave.h"
#include "scores.h"
#include "scoring.h"

int subScore (char char1, char char2, int stype) {
  int score = 0;
  if (stype == 1) { /* linear score */
    if ((char1 == GAPCHAR) || (char2 == GAPCHAR))
      score = -6;
    else if (char1 == char2)
      score = 3;
    else
      score = -2;
  }
  else if (stype == 2) /* PAM250 if protein */
    score  = PAM250 (char1, char2);
  else if (stype == 3) /* BLOSUM if protein */
    score  = BLOSUM (char1, char2);
  return score;
}

int gapScore(int stype) {
  int score = 0;
  score = subScore ('A', GAPCHAR, stype);
  return score;
}

/* Function to return the neighbors that got decremented , and those that didn't*/

long getRelation (long * cell, long * neighbor, long dimn, long * * decremented) {
  long i, cnt, ndecr = 0;
	char msg[SHORT_MESSAGE_SIZE];
  cnt = 0;
  ndecr = 0;
  for (i = 0; i< dimn; i++) {
    if (neighbor[i] < cell[i]) {
      (*decremented)[cnt] = i;
	    sprintf (msg, "decremented[%ld] = %ld\n", cnt, (*decremented)[cnt]);
  	  mprintf(12, msg, threadnum);
      cnt ++;
    }
    else 
      ndecr ++;
    sprintf (msg, "ngb %ld cell %ld cnt %ld ndcr %ld\n", neighbor[i], cell[i], cnt, ndecr);
    mprintf(12, msg, threadnum);
  }
  return ndecr;
}

/* Function to examine all temporary neighbors to determine the current score of the new cell*/
long getNeghbScore (long * m_index, long * * neighbor, MOA_rec * msaAlgn, long * * decremented, long * * * pwScores, char * * sequences, MOA_rec * NghbMOA, int stype, long totpwiseScore) {
	long score;
  long j, k, ndecr = 0;
	char msg[MID_MESSAGE_SIZE];
  
    /* get relation between this neighbor and the current cell */
  ndecr = getRelation (m_index, (*neighbor), msaAlgn->dimn, decremented);
    /* if we have 2 or more decremented neighbors*/
    /* sum pairwise scores for each pair of dimensions that got decremented*/
    totpwiseScore = 0;
    if ((msaAlgn->dimn  - ndecr) > 1) {
			sprintf(msg, "\ndecr %ld ndecr %ld  decremented Neighb ", (msaAlgn->dimn  - ndecr), ndecr);
			mprintf(12, msg, threadnum);
			for (j = 0; j < (msaAlgn->dimn  - ndecr);j++) {
	  		sprintf(msg, " i %ld = %ld ", j, (*decremented)[j]);
				mprintf(12, msg, threadnum);
			}
    }
      
    for (j=0;j<(msaAlgn->dimn  - ndecr) - 1; j++) {
			for (k=j+1;k<(msaAlgn->dimn  - ndecr); k++) {
	  		totpwiseScore += (*pwScores)[(*decremented)[j]][(*decremented)[k]];
		    sprintf(msg, "\n pw[%ld][%ld] = %ld seq1 %c seq2 %c ", (*decremented)[j], (*decremented)[k], (*pwScores)[(*decremented)[j]][(*decremented)[k]], sequences[j][m_index[j]-1], sequences[k][m_index[k]-1]);
				mprintf(12, msg, threadnum);
			}
    }
      
    /* multiple the number of dimensions that did not get decremented by the gap score */
    /* add both to the previoues score in the neighbor's cell in the alignment matrix to get the temporary score for this neighbor*/
	score = totpwiseScore + (ndecr * gapScore(stype));
	
	return score;
}

long initLBCell (long * m_index, long dimn, int stype) {
	long iscore, indexTau, i;

	indexTau = m_index[0]+1;
	//sprintf (msg, "%d %d\n",indexTau,  m_index[0] + 1);
	for (i = 1; i < dimn;i++) {
		indexTau = indexTau * (m_index[i] + 1);
		//sprintf (msg, "%d %d\n",indexTau,  m_index[i]);
	}      
	iscore = (gapScore(stype) * indexTau);
	return iscore;
}

int getNeighborScores (char * * sequences, long * seqLen, long * m_index, MOA_rec * msaAlgn, int stype, long * LNCount, long * * lnScores, long * * lnIndices, long * * * pwScores, long * * neighbor, long * * decremented, long * * wm_index, MOA_rec * NghbMOA) {

  long i, j, k, l, ndecr = 0;
  long totpwiseScore = 0;
  char char_a, char_b, msg[SHORT_MESSAGE_SIZE];
  /*Fill in the pair wise score matrix for the current cell*/
	mprintf(12, "Pairwise scores matrix: \n", threadnum);
  totpwiseScore = 0;
  for (l=0;l<msaAlgn->dimn  - 1; l++) {
    for (k=l+1;k<msaAlgn->dimn; k++) {
      if (Mode == Distributed) {
	      char_a = sequences[l][m_index[l]];
  	    char_b = sequences[k][m_index[k]];
      }
      else {
	      char_a = sequences[l][m_index[l]-1];
  	    char_b = sequences[k][m_index[k]-1];
      }

      (*pwScores)[l][k] = subScore(char_a, char_b, stype);
      totpwiseScore = totpwiseScore + (*pwScores)[l][k];
			sprintf(msg, "%ld, %ld, %c %c %ld ",l, k, char_a, char_b, (*pwScores)[l][k]);
			mprintf(12, msg, threadnum);
    }
    sprintf(msg, "\n tpw %ld ", totpwiseScore);
		mprintf(12, msg, threadnum);
  }
  
  /*get neighbors of the current cell */
  (*LNCount) = MOAGetLowerNeighbors (m_index, msaAlgn, NghbMOA);
  if (((*LNCount) == -1) && (Mode != Distributed)){
  	mprintf (1, "\nNo Lower neighbors while expecting, returning\n", threadnum);
    fflush(stdout);
    return -1;
  }
  /* it is distributed and it is a global lower border cell*/
  else if ((*LNCount) == -1) {
  	return -2; /* to let the caller initialize the border cell score*/
  }
  /*printMOA(NghbMOA); */
  
  /* loop throught neighbors */
  //mprintf(12, "with neighbors:", threadnum);
  for (i=0;i<NghbMOA->elements_ub - 1;i++) {
    /*sprintf(msg, "nghb %d containts %d and index %d\n", i, NghbMOA->elements[i].val, NghbMOA->indexes[i]);*/
		
    Gamma_Inverse(NghbMOA->indexes[i], msaAlgn->shape, msaAlgn->dimn, (*neighbor));
    //Gamma_Inverse(NghbMOA->indexes[i], msaAlgn->shape, msaAlgn->dimn, (*wm_index));    
		(*lnScores)[i] = msaAlgn->elements[NghbMOA->indexes[i]].val + getNeghbScore (m_index, neighbor, msaAlgn, decremented, pwScores, sequences, NghbMOA, stype, totpwiseScore);
		(*lnIndices)[i] = NghbMOA->indexes[i];
		sprintf(msg, " lnScores[%ld] %ld totpwiseScore %ld gapscore %ld ", i, (*lnScores)[i], totpwiseScore, ndecr * gapScore(stype));
		mprintf(12, msg, threadnum);
  }
  return 0;
}

int PrintPrevChains (MOA_rec *  msaAlgn) {
  long i, j = 0;
  char msg[MID_MESSAGE_SIZE];

  for (i = 0; i< msaAlgn->elements_ub; i++)  {
    sprintf (msg, "\nelm %ld with score %ld has %ld previous: ", i, msaAlgn->elements[i].val, msaAlgn->elements[i].prev_ub);
		mprintf(15, msg, threadnum);
   for (j = 0; j< msaAlgn->elements[i].prev_ub; j++)  {
     sprintf (msg, " %ld ", msaAlgn->elements[i].prev[j]);
		 mprintf(15, msg, threadnum);
   }
  }
  return 0;
}

int getPrevCells(long findex, long score, long LNCount, long * lnScores, long * lnIndices, MOA_rec *  msaAlgn) {
  int i = 0;
	char msg[MID_MESSAGE_SIZE];

  msaAlgn->elements[findex].prev_ub = 0;
  for (i = 0; i < LNCount; i++) {
		if (lnScores[i] == score) {
    	if ( msaAlgn->elements[findex].prev_ub == 0) 
	    	msaAlgn->elements[findex].prev = (unsigned long *) mmalloc (sizeof(unsigned long));
     	else 
      	msaAlgn->elements[findex].prev  = (unsigned long *) realloc (msaAlgn->elements[findex].prev, sizeof(unsigned long) * (msaAlgn->elements[findex].prev_ub + 1));

			if (msaAlgn->elements[findex].prev == NULL) {
				mprintf(1, "Could not reallocate memory for MSA Align Previous Cells pointer!\n", threadnum);
		   	return -1;
       }
     
     	msaAlgn->elements[findex].prev[msaAlgn->elements[findex].prev_ub] = lnIndices[i];
     	msaAlgn->elements[findex].prev_ub ++;      
	   	sprintf (msg, "\nelm = %ld prev_ub = %ld lnIndices[%d] = %ld\n", findex, msaAlgn->elements[findex].prev_ub, i, lnIndices[i] );
		 	mprintf(12, msg, threadnum);
  	}
	}		
  return 0;
}

long getScore (void * threadarg, long findex, long * m_index, long * LNCount, long * * lnScores, long * * lnIndices, long * * * pwScores, long * * neighbor, long * * decremented, long * * wm_index, MOA_rec * NghbMOA ) {

  long score, i = 0;
  int ret;
  char * * sequences;
  long seqNum, * seqLen;
  MOA_rec * msaAlgn; 
  int stype;
  Slave * slave;
	char msg[MID_MESSAGE_SIZE];

  slave = (Slave *) threadarg;

  sequences = slave->MOAPart[slave->processedPartitions-1].sequences;
  seqNum = slave->seqNum;
  seqLen = slave->seqLen;
  msaAlgn = slave->MOAPart[slave->processedPartitions-1].msaAlgn;
  stype = slave->stype;

  if (Mode == Distributed) {

		Gamma_Inverse(msaAlgn->indexes[findex], seqLen, seqNum, (*wm_index));    
     /*if (GLB) Global lower border cell from the whole tensor, then initialize with gabscores multiplied by indices if global, or zero if local alignment*/
    if (isLowerBorderCell((*wm_index), seqNum) == 1) {
	  	//mprintf (12, "Global Lower Border, initializing ", threadnum);
    		     
			score = initLBCell ((*wm_index), seqNum,stype);
				
	  	sprintf (msg, "Scoring GLB %ld init %ld \n", msaAlgn->indexes[findex], score);
			mprintf (10, msg, threadnum);
   	}
    		// is (LLB) Local lower Border cell from local partition, then block till it is received from the adjacent processor 
    else if (isLowerBorderCell(m_index, msaAlgn->dimn) == 1) {
	  	//mprintf (12, "Local Lower Border, receiving \n ", threadnum);

 	  	/*search for the overlapping cell in previously calculated partitions in local processor first*/

			ret = checkPrevPartitions (slave, msaAlgn->indexes[findex], &score, slave->MOAPart[slave->processedPartitions-1].waveNo);
			if (ret == -1) {
				/*if not found, search the Overlapping Cells scores received already from other slave processors*/
	
  			ret = checkRecvOC (msaAlgn->indexes[findex], &score);
				while (ret == -1) {
					sprintf (msg, "Waiting for an OC %ld in Local Partition %ld\n", msaAlgn->indexes[findex], slave->processedPartitions-1);
					mprintf (3, msg, threadnum);
	  			sem_wait(&slave->ocSem);
	  			ret = checkRecvOC(msaAlgn->indexes[findex], &score);
				/*if not found, block to receive it later from another slave after calculating its score*/
					//recvIC (msaAlgn->indexes[findex], &score);
//				score = -9999; /* for now, later block to receive*/
   			}
			}
   		sprintf (msg, "Scoring LLB %ld scores %ld \n", msaAlgn->indexes[findex], score);
	   	mprintf (10, msg, threadnum);
   	}
   	else {
	  	ret = getNeighborScores (sequences, seqLen, m_index, msaAlgn, stype, LNCount, lnScores, lnIndices, pwScores, neighbor, decremented, wm_index, NghbMOA);
			if (ret == -2) {
		/*I should never come here, if everything is alright, because this case should be trapped in the GLB case above*/
   			Gamma_Inverse(msaAlgn->indexes[findex], seqLen, seqNum, (*wm_index));    
      
   			score = initLBCell ((*wm_index), msaAlgn->dimn,stype);
   			sprintf (msg, "Scoring non-IC %ld, init %ld \n", msaAlgn->indexes[findex], score);
   			mprintf(10, msg, threadnum);
   		}
			else {
  /* take the maximum temporary score as the score for the current cell*/
  			score = a_max((*lnScores), (*LNCount));
  			i = getPrevCells(findex, score, (*LNCount), (*lnScores), 	(*lnIndices), msaAlgn);
  	 		sprintf (msg, "Scoring IC %ld scores %ld \n", msaAlgn->indexes[findex], score);
	   		mprintf(10, msg, threadnum);
  		}
   	}
 	}
	else {
	/*this is a repeatition for the last case above, Internal Cell for Distributed score calculation, repeated for the sequential calculation, however, I should have a separate function for sequential, should eventually be removed form here*/
  	ret = getNeighborScores (sequences, seqLen, m_index, msaAlgn, stype, LNCount, lnScores, lnIndices, pwScores, neighbor, decremented, wm_index, NghbMOA);
		if (ret == -2) {
		/*I should never come here, if everything is alright, because this case should be trapped in the GLB case above*/
   		Gamma_Inverse(msaAlgn->indexes[findex], seqLen, seqNum, (*wm_index));    
   		score = initLBCell ((*wm_index), msaAlgn->dimn,stype);
   		sprintf (msg, "Scoring non-IC %ld, init %ld \n", msaAlgn->indexes[findex], score);
   		mprintf(10, msg, threadnum);
   	}
		else {
  /* take the maximum temporary score as the score for the current cell*/
  		score = a_max((*lnScores), (*LNCount));
  		i = getPrevCells(findex, score, (*LNCount), (*lnScores), 	(*lnIndices), msaAlgn);
   		sprintf (msg, "Scoring IC %ld scores %ld \n", msaAlgn->indexes[findex], score);
   		mprintf(10, msg, threadnum);
  	}
  }
	if (AlignmentType == Local) {
		if (score <0)
			score = 0;
	}

  if (isHigherBorderCell(m_index, msaAlgn->dimn, msaAlgn->shape) == 1) {
  	OCout_ub ++;
  	if (OCout_ub == 1)
	  	OCout = (OCOType *) mmalloc(sizeof(OCOType));
  	else
	  	OCout = (OCOType *) realloc(OCout, sizeof(OCOType) * OCout_ub);
	  	
	  OCout[OCout_ub - 1].cellIndex = msaAlgn->indexes[findex];
	  OCout[OCout_ub - 1].cellScore = score;
	  OCout[OCout_ub - 1].waveNo = slave->MOAPart[slave->processedPartitions-1].waveNo;
	  OCout[OCout_ub - 1].toProc = 0;
    if (trySendOC (msaAlgn->indexes[findex], score) == 0) 
    	mprintf (10, "OC was not sent\n", threadnum);
	}
	msaAlgn->elements[findex].val = score;
	for (i = 0; i < (*LNCount); i++) {
		sprintf(msg, "lnScores[%ld] = %ld ", i, (*lnScores)[i]);
 		mprintf(12, msg, threadnum);		
	}
	for (i = 0; i < (*LNCount); i++) {
		sprintf(msg, "lnIndices[%ld] = %ld ", i, (*lnIndices)[i]);
 		mprintf(12, msg, threadnum);
	}
	sprintf(msg, "score %ld  \n", score);
	mprintf(12, msg, threadnum);
  return score;
}

void initTensor ( MOA_rec * msaAlgn, int stype) {
  long i, indexTau, j = 0;
  long * m_index = NULL;
	char msg[MID_MESSAGE_SIZE];

  m_index = (long *) mcalloc  (msaAlgn->dimn,  sizeof(long));
  for (i=0;i<msaAlgn->elements_ub;i++) {
    Gamma_Inverse(i, msaAlgn->shape, msaAlgn->dimn, m_index);
    if (isLowerBorderCell(m_index, msaAlgn->dimn) == 1) {
      indexTau = m_index[0];
      for (j=1;j<msaAlgn->dimn;j++) {
	indexTau = indexTau * (m_index[j] + 1);
      }
      msaAlgn->elements[i].val = (gapScore(stype) * indexTau);
    }
  }
  msaAlgn->elements[0].val = 0;
  if (m_index != NULL)
    free(m_index);
}

void * ComputePartitionScores (void * threadarg) {
  long i, j, startIndex, score, borderCell = 0;
  long * m_index = NULL; /* multidimensional index*/
  long LNCount, CalLnCount;
  long * lnScores = NULL; /* Lower neighbor temporary scores*/
  long * lnIndices = NULL; /* Lower neighbor Indices */
  long * * pwScores = NULL; /* holds the scores of pairwise match / mistmatch scores*/
  MOA_rec * NghbMOA = NULL;
  long * neighbor = NULL; /* multidimensional index of the current neighbor*/
  long * decremented = NULL; /* list of decremented indices in the multidimensional index of the current neighbor*/
  long * wm_index = NULL; /*multidimensional index from the whole tensor*/
  
  char * * sequences;
  long * seqLen;
  MOA_rec * msaAlgn; 
  int stype;
  Slave * slave;
  char msg[MID_MESSAGE_SIZE];

  slave = (Slave *) threadarg;

  sequences = slave->MOAPart[slave->processedPartitions-1].sequences;
  seqLen = slave->seqLen;
  msaAlgn = slave->MOAPart[slave->processedPartitions-1].msaAlgn;
  stype = slave->stype;
  CalLnCount = (long) mpow(2, msaAlgn->dimn) - 1;
  neighbor = (long *)  mcalloc ((msaAlgn->dimn), sizeof(long));
  decremented = (long *)  mcalloc ((msaAlgn->dimn), sizeof(long));
  wm_index = (long *)  mcalloc ((msaAlgn->dimn), sizeof(long));
  /*NghbMOA = (MOA_rec *) mmalloc(sizeof(MOA_rec));*/
  createMOAStruct (&NghbMOA);

  lnScores = (long *) mcalloc (CalLnCount, sizeof(long));
  lnIndices = (long *) mcalloc (CalLnCount, sizeof(long));

  /*get pairwise
 scores of the corresponding residues of the current cell in the alignment tensor*/
  /* Compute pairwise scores for the current index*/
  pwScores = (long * *) mcalloc (msaAlgn->dimn, sizeof(long *));
  m_index = (long *)  mcalloc ((msaAlgn->dimn), sizeof(long));
	
  for (i=0; i<msaAlgn->dimn; i++) {
    pwScores[i] = (long *) mcalloc (msaAlgn->dimn, sizeof(long));
	  if (Mode == Distributed) { 
  	  m_index [i] = 0; /* to retrieve the first non-border cell*/
		}    
		else {
			m_index [i] = 1; /* to retrieve the first non-border cell to save iterations on the loop below*/
		}
  }
  if (Mode == Distributed) 
		startIndex = 0;
	else
  	startIndex = Gamma(m_index, msaAlgn->dimn, msaAlgn->shape, msaAlgn->dimn, 1);
  
  
  /* loop the MOAAlign tensor for scores*/
  for (i = startIndex; i< msaAlgn->elements_ub; i++)  {
    sprintf(msg, "scoring cell %ld\n", i);
   	mprintf(12, msg, threadnum);
    currNow = getTime();
    if (((currNow->tm_yday * 1440) + (currNow->tm_hour * 60) + currNow->tm_min) > ((prevNow->tm_yday * 1440) + (prevNow->tm_hour * 60) + prevNow->tm_min + 10)) {

      sprintf(msg, "10 minutes elapsed and doing cell = %ld from total %ld and part %ld from parts %ld \n", i, msaAlgn->elements_ub, slave->processedPartitions-1, slave->partitionsCount);
      mprintf	(1, msg, threadnum);
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
    if (i != startIndex) /*because we already have it from above*/ {
      Gamma_Inverse(i, msaAlgn->shape, msaAlgn->dimn, m_index);
		}
    LNCount = 0;
    
    if ((Mode == Distributed) || ((Mode != Distributed) && (isLowerBorderCell(m_index, msaAlgn->dimn) == 0))) {
      score = getScore(slave, i, m_index, &LNCount,  &lnScores, &lnIndices, &pwScores, &neighbor, &decremented, &wm_index, NghbMOA);
    //if (LNCount > 0)
      msaAlgn->elements[i].val = score;
    }
    sprintf(msg, " cell %ld index %ld score %ld \n", i, msaAlgn->indexes[i], msaAlgn->elements[i].val);
		mprintf(10, msg, threadnum);
  }

  if (m_index != NULL)
    free(m_index);
  if (lnScores != NULL)
    free(lnScores);
  if (lnIndices != NULL)
    free(lnIndices);
 
  deleteMOA (NghbMOA);
  if (neighbor != NULL)
    free(neighbor);
  if (decremented != NULL)
    free(decremented);
  if (wm_index != NULL)
    free(wm_index);
  if (pwScores != NULL) {
    for (i=0; i<msaAlgn->dimn; i++) {
      if (pwScores[i] != NULL) 
        free( pwScores[i]);
    }
    free(pwScores);
  }
  
  slave->MOAPart[slave->processedPartitions-1].msaAlgn = msaAlgn;
  return NULL;
}
