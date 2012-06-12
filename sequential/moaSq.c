#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>      /* include the character library */
#include "moa.h"
#include "utils.h"
#include <errno.h>
#include "moamsa.h"


int main (int argc, char ** argv) {
  long seqNum, i, k  = 0;
 
  char * * sequences = NULL;
  long * seqLen = NULL;
  char * * * algnseq = NULL;
  long * aSeqLen = NULL;
  int alignmentsNo=0;
  MOA_rec * msaAlgn = NULL;
 
  int stype = 0;
  long currentScore, currentCell = 0;
  long prevScore, prevCell = 0;
  Mode = Sequential;
  processArguments(argc, argv, &seqNum, &sequences, &seqLen, &stype);
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
    

  
  /* A. Crteate MOA Alignment Tensor */
  /*msaAlgn = (MOA_rec *) mmalloc(sizeof(MOA_rec));*/
  createMOAStruct (&msaAlgn);
  
  createMOA(seqLen /* shape*/, seqNum /* dimension*/, msaAlgn /* MOA structure*/,0,0);
  if (pdebug == 1)
    mprintf(outputfilename, "MOA dimn %d, elm ub %d\n", 2, msaAlgn->dimn, msaAlgn->elements_ub); 
  /* B. Fill the Tensor */
 
  if (AlignmentType == Global) 
    initTensor (msaAlgn, stype);
  fillTensor (sequences, msaAlgn, stype);
  if (pdebug == 1) 
    printMOA(msaAlgn);

  /* C. trace back */
  aSeqLen = (long *)   mmalloc (sizeof(long));
  aSeqLen[0] = 0;
  prevCell = -1;
  algnseq = (char * * *) mmalloc(sizeof(char * *));
  algnseq[0] = (char * *) mmalloc (seqNum * sizeof(char *));    
  if (AlignmentType == Global) { /* if Global Alignment */
    /*PrintPrevChains(msaAlgn);
    // Get Max Cell on Last Border as Current Cell */
    alignmentsNo = 0;
    currentScore = getMaxOnLastBorder (msaAlgn, &currentCell);
    traceBack (seqNum, sequences, seqLen, msaAlgn, stype, &algnseq, &aSeqLen, &alignmentsNo, &currentCell, &currentScore, 0);
  }
  else { /* if Local Alignment */
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
		mprintf(outputfilename, "Could not reallocate memory for Aligned Sequence Set %d!\n", 1, alignmentsNo + 1);
		return -1;
      }
	algnseq[alignmentsNo] = (char * *) mmalloc (seqNum * sizeof(char *));
	aSeqLen = (long *)  realloc (aSeqLen, (alignmentsNo + 1) * sizeof(long));
     if (aSeqLen == NULL) {
		mprintf(outputfilename, "Could not reallocate memory for Aligned Sequence Length %d!\n", 1, alignmentsNo + 1);
		return -1;
      }
	traceBack_loc (seqNum, sequences, seqLen, msaAlgn, stype, &algnseq, &aSeqLen, &alignmentsNo, &currentCell, &currentScore, 0);
	}
	prevCell = currentCell;
	prevScore = currentScore;
    }
  }

  /* D. Print the resulting Alignemnts */
  PrintASeq (seqNum, sequences, seqLen, &algnseq, aSeqLen, alignmentsNo+1) ;

  /* Free all Memory Allocations & Exit. */
  deleteMOA (msaAlgn);
  if (sequences != NULL) {
    for (i=0;i<seqNum;i++) {
      if (sequences[i] != NULL) 
        free(sequences[i]);
	}
    free(sequences);
  }
  if (algnseq != NULL) {
    for (k=0;k<=alignmentsNo;k++) {
	  if (algnseq[k] != NULL) {
        for (i=0;i<seqNum;i++) {
          if (algnseq[k][i] != NULL)
	        free(algnseq[k][i]);
        }
        free(algnseq[k]);
	  }
    }
    free(algnseq);
  }
  if (seqLen != NULL)
    free(seqLen);
  if (aSeqLen != NULL)
    free(aSeqLen);
  if (prevNow != NULL)
    free(prevNow);

  return 0;

}


