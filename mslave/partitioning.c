#include <stdio.h>
#include <stdlib.h>
#include <semaphore.h>
#include "globals.h"
#include "moaDst.h"
#include "utils.h"
#include "moa.h"
#include "master.h"
#include "lq.h"
#include <pthread.h>

/* ***************************************************************
// **********************   DistributeDiagonals  *****************
// ***************************************************************/

void manageDependancies (Master *  master, long partitionSize, long waveNo, long startPartNo, long EndPartNo) {
	char msg[SHORT_MESSAGE_SIZE];
  long i, j, k, l, m, * ind;
  int notsent;
  
  ind = (long *) mmalloc (sizeof(long) * master->msaAlgn->dimn);
  /*search the partitions of this wave*/
  for (i = startPartNo; (i < EndPartNo) && ((DS[i].waveNo == waveNo) || (DS[i].waveNo == waveNo-1)); i++) {
  	/*for all HB Cells*/
    for (j = 0; j < DS[i].hbCellsCount ; j++) {   
      Gamma_Inverse(DS[i].hbCells[j], master->msaAlgn->shape, master->msaAlgn->dimn, ind);
      /*that are not LB in the whole tensor*/
      if (isLowerBorderCell(ind,master->msaAlgn->dimn) == 0) {
     	 /*and all partitions*/
				for (k = startPartNo; (k < DS_ubound) && (DS[k].waveNo <= (waveNo+1)); k++) {	
	     	 /*that has an LB equal to the previous HB*/
	  			for (l = 0; l < DS[k].lbCellsCount ; l++) { 
	    			if ( DS[i].hbCells[j] ==  DS[k].lbCells[l]) {
	      			/* if not sent before, send it now */
				      notsent = 0;
							if (depProclistCount > 0) {
					    	for (m=0; (m<depProclistCount) && (notsent == 0);m++) {
									if ((depProclist[m].cellIndex == DS[i].hbCells[j]) && (depProclist[m].depProc == DS[k].proc))
									  notsent = 1;
						    }
							}
							sprintf (msg, "Index %ld notsent %d DS[i].proc %d in wave %ld DS[k].proc %d in wave %ld\n", DS[i].partIndex, notsent, DS[i].proc,  DS[i].waveNo, DS[k].proc, DS[k].waveNo);
							mprintf (15, msg, threadnum);
							if ((notsent == 0) && (DS[i].proc != DS[k].proc)) {
								sendCellDepProc (waveNo, DS[i].proc, DS[i].hbCells[j], DS[i].partIndex, DS[k].proc);
								sprintf (msg, " sent to %d Cell %ld wave %ld to %d in wave %ld\n", DS[i].proc, DS[i].hbCells[j], DS[i].waveNo, DS[k].proc, DS[k].waveNo);
							  mprintf (10, msg, threadnum);
								depProclistCount ++; /* update list to avoid sending again */
								if (depProclistCount == 1)
								  depProclist = (depProc_rec *) mmalloc (sizeof(depProc_rec));
								else
								  depProclist = (depProc_rec *) realloc (depProclist, sizeof(depProc_rec) * depProclistCount);
								depProclist[depProclistCount - 1].cellIndex = DS[i].hbCells[j];
								depProclist[depProclistCount - 1].depProc = DS[k].proc;
     					}

	    			}		
	  			}
				}
      }
    }
  }
  free (ind);  
}

void getNextProc (Master *  master, int * newProcessor) {
  int errCode;
  
  (*newProcessor) = -1;
	
  //sem_wait(&master->qSem);
  errCode = pthread_mutex_lock (master->mut);
  checkMutexErrorCode (2, errCode);
  (*newProcessor) = dequeue(&processors); 
  enqueue(&processors, (*newProcessor)); /* to make it a circular queue that doesn't stop for the proc to return*/
  errCode = pthread_mutex_unlock (master->mut);
  checkMutexErrorCode (2, errCode);
  
}

int IsAssignedtoprocessor (long partIndex) {
  long i;
  for (i = 0; i < DS_ubound; i++) {
    if (DS[i].partIndex == partIndex) 
      return 1;
  }  
  return 0;
}

long LaterWavesDistributeDiagonals (long nextWavePartStart, Master * master, long waveNo, long waveSize) {
  long prevWaveEnd;
  MOA_rec * MOA_nghb;
  int moreneighb;
  int newProcessor;
  long i, j;
  long * ind;
	char msg[MID_MESSAGE_SIZE];

  ind =  (long *) mmalloc (sizeof(long) * master->msaAlgn->dimn);
  prevWaveEnd = DS_ubound;  
  /*loop for this wave partitions*/
  for (i = nextWavePartStart; i<prevWaveEnd; i++) {
	
    createMOAStruct (&MOA_nghb);
    
    Gamma_Inverse(DS[i].partIndex, master->msaAlgn->shape, master->msaAlgn->dimn, ind);    
    
 	 /*get the partition*/
    moreneighb = MOAGetHigherNeighbors (waveSize, ind, master->msaAlgn, MOA_nghb);
    
    if (moreneighb > 0) {
      
	 	 /*send it to the next available processor*/
      getNextProc(master, &newProcessor);
      sprintf(msg, "Master returned From GetProc i %ld waveNo %ld\n", i, waveNo);
      mprintf (15, msg, threadnum); 
      SendPartition (master, newProcessor, MOA_nghb, waveNo);
      sprintf(msg, "Part %ld waveNo %ld to proc %d starting %ld\n", i, waveNo, newProcessor, DS[i].partIndex);
      mprintf (1, msg, threadnum); 
      DS[i].proc = newProcessor;
      DS[i].sent = 1;
      DS[i].lbCellsCount = 0;
      DS[i].hbCellsCount = 0;

	 	 /*identify the partitions HB and LB Cells*/
      
      sprintf(msg, "identify the partitions HB and LB Cells\n");
      mprintf (15, msg, threadnum); 
      for (j = 0; j < MOA_nghb->elements_ub;j++) {
	Gamma_Inverse(j, MOA_nghb->shape, MOA_nghb->dimn, ind);
	sprintf(msg, "After Gam_Inv\n");
	mprintf (15, msg, threadnum); 
	if (isLowerBorderCell(ind, MOA_nghb->dimn) == 1) {
	  DS[i].lbCellsCount ++;
	  if (DS[i].lbCellsCount == 1)
	    DS[i].lbCells = (long * ) mmalloc (sizeof(long));
	  else
	    DS[i].lbCells = (long * ) realloc (DS[i].lbCells, sizeof(long) * DS[i].lbCellsCount);
	  
	  DS[i].lbCells[DS[i].lbCellsCount - 1] = MOA_nghb->indexes[j];
	}
	mprintf(15, "Finished LB Cells\n", threadnum);
	if ( isHigherBorderCell (ind, MOA_nghb->dimn, MOA_nghb->shape) == 1) {
	  DS[i].hbCellsCount ++;
	  if (DS[i].hbCellsCount == 1)
	    DS[i].hbCells = (long * ) mmalloc (sizeof(long));
	  else
	    DS[i].hbCells = (long * ) realloc (DS[i].hbCells, sizeof(long) * DS[i].hbCellsCount);
	  
	  DS[i].hbCells[DS[i].hbCellsCount - 1] = MOA_nghb->indexes[j];
	}
	mprintf(15, "Finished HB \n", threadnum);
      }
      /*read next corner cells to form the basis for the next wave, according to the partitions in this wave to form a breadth first search, then send all next wave partitions at once, in the later call to "LaterWavesDistributeDiagonals"*/
      sprintf(msg, "Master getting Next Wave from i %ld in wave %ld\n", i, waveNo);
      mprintf (15, msg, threadnum); 
      getNextWave (master, waveNo, MOA_nghb);
      deleteMOA (MOA_nghb);  
    }
  }
  free(ind);
  //if (prevWaveEnd < DS_ubound)
    return  prevWaveEnd;
  //else
  //  return 0;
}

void getNextWave (Master * master, long waveNo, MOA_rec * MOA_nghb) {
  long i, j, k, l, border, flatIndex; 
  long * nghb_ind;
  long * ind;
  long combnum; /* number of combinations for decrementing i number of indices*/
  long * * combin; /* combinations of lower neighbor indices*/
	char msg[SHORT_MESSAGE_SIZE];

  nghb_ind =  (long *) mmalloc (sizeof(long) * MOA_nghb->dimn);
  ind =  (long *) mmalloc (sizeof(long) * MOA_nghb->dimn);
  for (j = 1;j<MOA_nghb->dimn; j++)  {
    
    /* number of the different combinations of which indices in the multidimensional index to give the shape bound so that all edge points can be processed */
    combnum = (Factorial(MOA_nghb->dimn)/(Factorial(MOA_nghb->dimn-j) * Factorial(j)));
    /* create memory for combinations matrix */
    combin = (long * *) calloc (combnum, sizeof(long *));
    for (k=0; k<combnum; k++)
      combin[k] = (long *) calloc (j, sizeof(long));
    if (combin == NULL) {
      mprintf (1, "getNextWave Error: Can not allocate memory to combinations matrix\n", threadnum);
      return;
    }
    /* get matrix of all possible combinations */
    Combinations(MOA_nghb->dimn, j, &combin);
    /* loop to decrement the selected indices */
    for (k=0; k<combnum; k++) {
      for (l = 0; l <MOA_nghb->dimn; l++) {
				nghb_ind[l] = 0;
      }
      for (l=0; l<j; l++) {
				nghb_ind[combin[k][l]-1] = MOA_nghb->shape[combin[k][l]-1] - 1;
      }
      for (l = 0; l <MOA_nghb->dimn; l++) {
				sprintf(msg, "nghb_ind[%ld] = %ld\n", l, nghb_ind[l]);
				mprintf (15, msg, threadnum);
      }
      /* get the index of the edge cell in the current partition */
      flatIndex = Gamma( nghb_ind, MOA_nghb->dimn,  MOA_nghb->shape, MOA_nghb->dimn, 1);
			sprintf (msg, "Local Index = %ld Global Index %ld\n", flatIndex, MOA_nghb->indexes[flatIndex]);
			mprintf (15, msg, threadnum);
      /* get the index of the edge cell in the whole tensor */      
      Gamma_Inverse(MOA_nghb->indexes[flatIndex], master->msaAlgn->shape, master->msaAlgn->dimn, ind);
      /*check if it is an edge cell in the whole tensor*/
      border = 1;
      for (l=0; l<master->msaAlgn->dimn; l++) {
				if (ind[l] >= master->msaAlgn->shape[l] - 1) 
				  border = 0;
      }
      
      /*check if it is previously assigned to a processor*/
      if (IsAssignedtoprocessor (MOA_nghb->indexes[flatIndex]) == 1) {
				border = 0;
      }
      /*if not edge, and not assigned before, then it is a new partition*/
      if (border == 1) {
				sprintf (msg, "Assigning to proc index %ld\n", MOA_nghb->indexes[flatIndex]);
				mprintf (15, msg, threadnum);
				DS_ubound ++;
				if (DS_ubound == 1)
				  DS = (Distribution * ) mmalloc (sizeof(Distribution));
				else
				  DS = (Distribution * ) realloc (DS, sizeof(Distribution) * DS_ubound);
				DS[DS_ubound - 1].partIndex = MOA_nghb->indexes[flatIndex];
				DS[DS_ubound - 1].waveNo = waveNo + 1;
				DS[DS_ubound - 1].sent = 0;

			}
    }
  }
  /* do the final edge where all at the shape limit*/
  for (l = 0; l <MOA_nghb->dimn; l++) {
    nghb_ind[l] = MOA_nghb->shape[l] - 1;
		sprintf(msg, "nghb_ind[%ld] = %ld\n", l, nghb_ind[l]);
		mprintf (15, msg, threadnum);
  }
  flatIndex = Gamma( nghb_ind, MOA_nghb->dimn,  MOA_nghb->shape, MOA_nghb->dimn, 1);
	sprintf (msg, "Local Index = %ld Global Index %ld\n", flatIndex, MOA_nghb->indexes[flatIndex]);
	mprintf (15, msg, threadnum);
  
  Gamma_Inverse(MOA_nghb->indexes[flatIndex], master->msaAlgn->shape, master->msaAlgn->dimn, ind);    
  /*check if it is an edge cell in the whole tensor*/
  border = 1;
  for (l=0; l<master->msaAlgn->dimn; l++) {
    if (ind[l] >= master->msaAlgn->shape[l] - 1) 
      border = 0;
  }
  /*check if it is previously assigned to a processor*/
  if (IsAssignedtoprocessor (MOA_nghb->indexes[flatIndex]) == 1) {
		border = 0;
  }
  /*if not edge, and not assigned before, then it is a new partition*/
  if (border == 1) {
				sprintf (msg, "Assigning to proc index %ld\n", MOA_nghb->indexes[flatIndex]);
				mprintf (15, msg, threadnum);
    DS_ubound ++;
    if (DS_ubound == 1)
      DS = (Distribution * ) mmalloc (sizeof(Distribution));
    else
      DS = (Distribution * ) realloc (DS, sizeof(Distribution) * DS_ubound);
    DS[DS_ubound - 1].partIndex = MOA_nghb->indexes[flatIndex];
    DS[DS_ubound - 1].waveNo = waveNo + 1;
    DS[DS_ubound - 1].sent = 0;
  }
  free(nghb_ind);
  free(ind);
 }

int FirstWaveDistributeDiagonals (Master * master, long waveNo, long waveSize) {
  MOA_rec * MOA_nghb = NULL;
  long * m_ind = NULL;
  long nextWavePartStart, i; 
  int moreneighb;
  int newProcessor;
  long * ind = NULL;
	char msg[MID_MESSAGE_SIZE];

  createMOAStruct (&MOA_nghb);
 
  ind = (long *) mmalloc ( master->msaAlgn->dimn * sizeof(long));

  for (i = 0; i <  master->msaAlgn->dimn; i++) {
    ind[i] = 0;
  }
  moreneighb = MOAGetHigherNeighbors (waveSize, ind, master->msaAlgn, MOA_nghb);  
  sem_wait(&master->qSem);
  getNextProc(master, &newProcessor);
  if (moreneighb <= 0)
    return 0;
  SendPartition (master, newProcessor, MOA_nghb, waveNo);
  sprintf(msg, "Part 0 waveNo %ld to proc %d starting %ld\n", waveNo, newProcessor, MOA_nghb->indexes[0]);
  mprintf (1, msg, threadnum); 
  
  DS_ubound ++;
  if (DS_ubound == 1)
    DS = (Distribution * ) mmalloc (sizeof(Distribution));
  else
    DS = (Distribution * ) realloc (DS, sizeof(Distribution) * DS_ubound);
  DS[DS_ubound - 1].proc = newProcessor;
  DS[DS_ubound - 1].partIndex = MOA_nghb->indexes[0];
  DS[DS_ubound - 1].waveNo = waveNo;
  DS[DS_ubound - 1].sent = 1;
  DS[DS_ubound - 1].lbCellsCount = 0;
  DS[DS_ubound - 1].hbCellsCount = 0;

  sprintf(msg, "identify the partitions HB and LB Cells\n");
  mprintf (15, msg, threadnum); 
  /*identify the partitions HB and LB Cells*/
  
  m_ind =  (long *) mmalloc (sizeof(long) * MOA_nghb->dimn);
  for (i = 0; i < MOA_nghb->elements_ub;i++) {
    Gamma_Inverse(i, MOA_nghb->shape, MOA_nghb->dimn, m_ind);
    sprintf(msg, "After Gam_Inv\n");
    mprintf (15, msg, threadnum); 
    
    if (isLowerBorderCell(m_ind, MOA_nghb->dimn) == 1) {
      DS[DS_ubound - 1].lbCellsCount ++;
      if (DS[DS_ubound - 1].lbCellsCount == 1)
	DS[DS_ubound - 1].lbCells = ( long* ) mmalloc (sizeof(long));
      else
	DS[DS_ubound - 1].lbCells = (long * ) realloc (DS[DS_ubound - 1].lbCells, sizeof(long) * DS[DS_ubound - 1].lbCellsCount);
      
      DS[DS_ubound - 1].lbCells[DS[DS_ubound - 1].lbCellsCount - 1] = MOA_nghb->indexes[i];
    }
    mprintf(15, "Finished LB Cells\n", threadnum);
    if ( isHigherBorderCell (m_ind, MOA_nghb->dimn, MOA_nghb->shape) == 1) {
      DS[DS_ubound - 1].hbCellsCount ++;
      if (DS[DS_ubound - 1].hbCellsCount == 1)
	DS[DS_ubound - 1].hbCells = ( long* ) mmalloc (sizeof(long));
      else
	DS[DS_ubound - 1].hbCells = (long * ) realloc (DS[DS_ubound - 1].hbCells, sizeof(long) * DS[DS_ubound - 1].hbCellsCount);
      
      DS[DS_ubound - 1].hbCells[DS[DS_ubound - 1].hbCellsCount - 1] = MOA_nghb->indexes[i];
    }
  }
  mprintf(15, "Finished HB Cells\n", threadnum);
  nextWavePartStart = DS_ubound;  
  getNextWave (master, waveNo, MOA_nghb);
  mprintf(15, "After getting First Next Wave \n", threadnum);
  if (m_ind != NULL)
    free(m_ind);
  mprintf(15, "After 1 \n", threadnum);
  //if (ind != NULL)
  //  free(ind);
  mprintf(15, "After 2 \n", threadnum);
  deleteMOA (MOA_nghb);  
  mprintf(15, "After 3 \n", threadnum);
  if (nextWavePartStart < DS_ubound)
    return nextWavePartStart;
  else
    return 0;
}
