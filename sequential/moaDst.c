#include <stdio.h>
#include<string.h>
#include <stdlib.h>
#include <ctype.h>      /* include the character library */
#include <pthread.h>
#include "moa.h"
#include "utils.h"
#include <mpi.h>
#include "lq.h"
#include "moamsa.h"

const char GABCHAR = '-';


void * ProcessorGateProcess (void *args);
void * PartitioningProcess (void *args);

void SendPartition (int toProcessor, MOA_rec * MOA_partition) {
  MPI_Request request;/*Non-blocking send or receive status*/
  int ComputationPhase = 1;
  MPI_Isend (&ComputationPhase, 1, MPI_INT, toProcessor, 98, MPI_COMM_WORLD, &request);
  sendMOAIndices (toProcessor, MOA_partition); 
}
/* ***************************************************************
// **********************   DistributeDiagonals  *****************
// ***************************************************************/
void DistributeDiagonals (Master * master, long waveNo, long waveSize, long * ind, MOA_rec * MOA_in, int currProcessor) {
  MOA_rec * MOA_nghb;
  long i, j, k, l, border, flatIndex; 
  int moreneighb;
  long * nghb_ind;
  long combnum; /* number of combinations for decrementing i number of indices*/
  long * * combin; /* combinations of lower neighbor indices*/
  int newProcessor;

  createMOAStruct (&MOA_nghb);
 
  moreneighb = MOAGetHigherNeighbors (waveSize, ind, MOA_in, MOA_nghb);
  if (moreneighb <= 0)
    return;
  else {
    mprintf(outputfilename, "\nWave %d moreneighb %d MOA_nghb->elements_ub %d from point : ", 3, waveNo , moreneighb, MOA_nghb->elements_ub);
    for (j = 0; j <  MOA_in->dimn; j++) {
      mprintf(outputfilename, "%ld ", 1, ind[j]);     
    }
    
    mprintf(outputfilename, " has %d points with Indices: ", 1, moreneighb);
    for (i=0;i<MOA_nghb->elements_ub; i++) {
      mprintf(outputfilename, "%ld ", 1, MOA_nghb->indexes[i]);    
      /*assign point to the next available processor*/
      if (MOA_in->proc[MOA_nghb->indexes[i]] != 0) {
	printf ("\n partition starting at %d was assigned from before.", MOA_in->proc[MOA_nghb->indexes[i]]);
	return;
      }
      MOA_in->proc[MOA_nghb->indexes[i]] = currProcessor;
      /* send to that processor*/
    }

    SendPartition (currProcessor, MOA_nghb);


    
    nghb_ind =  (long *) mmalloc (sizeof(long) * MOA_nghb->dimn);
      
    for (j = 1; j< MOA_nghb->dimn; j++)  {
      
      /* number of the different combinations of which indeces in the multidimensional index to give the shape bound so that all edge points can be processed */
      combnum = (Factorial(MOA_nghb->dimn)/(Factorial(MOA_nghb->dimn-j) * Factorial(j)));
      /* create memory for combinations matrix */
      combin = (long * *) calloc (combnum, sizeof(long *));
      for (k=0; k<combnum; k++)
	combin[k] = (long *) calloc (j, sizeof(long));
      if (combin == NULL) {
	printf ("Can not allocate memory to combinations matrix\n");
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
	flatIndex = Gamma( nghb_ind, MOA_nghb->dimn,  MOA_nghb->shape, MOA_nghb->dimn, 1);
	
	Gamma_Inverse(MOA_nghb->indexes[flatIndex], MOA_in->shape, MOA_in->dimn, ind);
	
	border = 1;
	for (l=0; l<MOA_in->dimn; l++) {
	  if (ind[l] >= MOA_in->shape[l] - 1) {
	    border = 0;
	  }
	}
	if (border == 1) {
	  pthread_mutex_lock (master->mut);
	  newProcessor = dequeue(&processors);
	  pthread_mutex_unlock (master->mut);
	  while (newProcessor == -1) {
	    pthread_cond_wait (master->notEmpty, master->mut);    
	    pthread_mutex_lock (master->mut);
	    newProcessor = dequeue(&processors);
	    pthread_mutex_unlock (master->mut);
	  }
	printf ("\n getting neighbors for proc %d, %d for partition starting at: %d", currProcessor, newProcessor, MOA_nghb->indexes[flatIndex]);
	  MOA_in->dproc_ub[MOA_nghb->indexes[flatIndex]] ++;
	printf ("\n neighbors for proc %d are %d", currProcessor, MOA_in->dproc_ub[MOA_nghb->indexes[flatIndex]]);
	  if (MOA_in->dproc_ub[MOA_nghb->indexes[flatIndex]] == 1) 

	      MOA_in->dproc[MOA_nghb->indexes[flatIndex]] = (int *) mcalloc (1, sizeof(int));

	  else

	      
	    MOA_in->dproc[MOA_nghb->indexes[flatIndex]] = (int *) realloc (MOA_in->dproc[MOA_nghb->indexes[flatIndex]],  MOA_in->dproc_ub[MOA_nghb->indexes[flatIndex]]);


	  MOA_in->dproc[MOA_nghb->indexes[flatIndex]][MOA_in->dproc_ub[MOA_nghb->indexes[flatIndex]]] = newProcessor;

	  DistributeDiagonals (master, waveNo + 1, waveSize, ind, MOA_in, newProcessor);

	}
      }
    }
    /* do the final edge where all at the shape limit*/
    for (l = 0; l <MOA_nghb->dimn; l++) {
      nghb_ind[l] = MOA_nghb->shape[l] - 1;
    }
    flatIndex = Gamma( nghb_ind, MOA_nghb->dimn,  MOA_nghb->shape, MOA_nghb->dimn, 1);
    
    Gamma_Inverse(MOA_nghb->indexes[flatIndex], MOA_in->shape, MOA_in->dimn, ind);
    
    border = 1;
    for (l=0; l<MOA_in->dimn; l++) {
      if (ind[l] >= MOA_in->shape[l] - 1) {
	border = 0;
      }
    }
    if (border == 1) {
      getDiagonals (waveNo + 1, waveSize, ind, MOA_in);
    }
    free(nghb_ind);
  }
  deleteMOA (MOA_nghb);  
}
void * ProcessorGateProcess(void *m) {
  int i, procToEnq;
  Master * master;
  MPI_Status status; /*Message receive strcucture*/
	
  master = (Master *)m;
  pthread_mutex_lock (master->mut);
  for (i=0;i<ClsuterSize;i++) {
    if (myProcid != i) {
      enqueue(&processors, i);
      pthread_cond_signal (master->notEmpty);
    }
  }
  pthread_mutex_unlock (master->mut);
  while (1) {
    MPI_Recv (&procToEnq, 1, MPI_INT, MPI_ANY_SOURCE, 99,MPI_COMM_WORLD, &status);
	printf ("\n[%d] returned to the master scheduler", procToEnq);
    pthread_mutex_lock (master->mut);
    enqueue(&processors, procToEnq);
    pthread_mutex_unlock (master->mut);
  }
  return NULL;
}

void * PartitioningProcess(void *m) {
  long i, WaveSize;
  long * ind = NULL;
  MPI_Request request;/*Non-blocking send or receive status*/
  Master * master;
  int ComputationPhase; /* 1 = Get Scores for more partitions - 2 = Trace Back*/
  int newProcessor;
  ComputationPhase = 1;								
  master = (Master *)m;
  ind = (long *) mmalloc ( master->msaAlgn->dimn * sizeof(long));
  
  /*printMOA(MOA1);*/
  
  for (i = 0; i <  master->msaAlgn->dimn; i++) {
    ind[i] = 1;
  }
  pthread_mutex_lock (master->mut);
  newProcessor = dequeue(&processors);
  pthread_mutex_unlock (master->mut);
  while (newProcessor == -1) {
    pthread_cond_wait (master->notEmpty, master->mut);    
    pthread_mutex_lock (master->mut);
    newProcessor = dequeue(&processors);
    pthread_mutex_unlock (master->mut);
  }
  DistributeDiagonals (master, 1, partitionSize, ind, master->msaAlgn, newProcessor);
  ComputationPhase = 2;
  for (i=0;i<ClsuterSize;i++) {
    if (i != myProcid)
      MPI_Isend (&ComputationPhase, 1, MPI_INT, i, 98, MPI_COMM_WORLD, &request);
  }
  
  free(ind);
  return NULL;
}

void SlaveProcess () {
  MOA_rec * MOA_in;
  long i;
  char buf;
  int ComputationPhase; /* 1 = Get Scores for more partitions - 2 = Trace Back*/
  MPI_Request request;/*Non-blocking send or receive status*/
  MPI_Status status; /*Message receive structure*/
  ComputationPhase = 1;
  while (ComputationPhase == 1) {
    MPI_Recv (&ComputationPhase, 1, MPI_INT, MPI_ANY_SOURCE, 98,MPI_COMM_WORLD, &status);
    if (ComputationPhase == 1) {
      /* receive MOA partition indices*/
      createMOAStruct (&MOA_in);
      recivMOAIndices  (&MOA_in);
      /*start computing the score*/
      /*wait for dependancy*/
      /*finish and notify scheduler of finishing to recive further jobs*/
      printf("\n[%d] doing the following partition:", myProcid);
      printMOA (MOA_in);
      deleteMOA(MOA_in);
      MPI_Isend (&myProcid, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, &request);
    }
    else {
      printf("\n[%d] going to phase 2", myProcid);
      fflush(stdout);
    }
  }
}
int recivMOAIndices (MOA_rec * * MOA_in) {
  int src;
  int msgsize;
  MPI_Status status;
  MPI_Request request;
  char * buffer;
  unsigned int	membersize, maxsize;
  int dest, tag;
  int position;
  long i, j, lmessage;

  /*1) dimn */
  MPI_Recv(&(*MOA_in)->dimn, 1, MPI_LONG,
	   MPI_ANY_SOURCE , 1 , MPI_COMM_WORLD, &status);
  /*2) shape */
  (*MOA_in)->shape = (long *) mcalloc((*MOA_in)->dimn, sizeof(long));
  for (i=0;i<(*MOA_in)->dimn;i++) {
    MPI_Recv(&lmessage, 1, MPI_LONG,
	     MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);
    (*MOA_in)->shape[i] = lmessage;
  }
  /*3) elements_ub  */
  MPI_Recv(&(*MOA_in)->elements_ub, 1, MPI_LONG,
	   MPI_ANY_SOURCE , 3 , MPI_COMM_WORLD, &status);
  
  /*4) indexes */
  (*MOA_in)->elements = (MOA_elm * ) mcalloc ((*MOA_in)->elements_ub, sizeof(MOA_elm));
  (*MOA_in)->indexes = (signed long * ) mcalloc ((*MOA_in)->elements_ub, sizeof(signed long));
  for (i=0;i<(*MOA_in)->elements_ub;i++) {
    MPI_Recv(&lmessage, 1, MPI_LONG,
	     MPI_ANY_SOURCE, 4, MPI_COMM_WORLD, &status);
    (*MOA_in)->indexes[i] = lmessage;
    (*MOA_in)->elements[i].val = i;
  }
  
}


void Diagonalize (MOA_rec * MOA_in) {

}

int sendMOAIndices (int ProcDest, MOA_rec * MOA_in) {
   unsigned int	membersize, maxsize;
    int			position;
    int			dest, tag;
    char		*buffer;
    long i, j, PrevCellsCount, lmessage;
    MPI_Request request; 
/* Send the MOA Buffer*/
  /*1) dimn */
    MPI_Isend(&MOA_in->dimn, 1, MPI_LONG,
	    ProcDest, 1, MPI_COMM_WORLD, &request );
  /*2) shape */
  for (i=0;i<MOA_in->dimn;i++) {
    lmessage = MOA_in->shape[i];
    MPI_Isend(&MOA_in->shape[i], 1, MPI_LONG,
	    ProcDest, 2, MPI_COMM_WORLD, &request );
  }
    /*3) elements_ub  */
    MPI_Isend(&MOA_in->elements_ub, 1, MPI_LONG,
	    ProcDest, 3, MPI_COMM_WORLD, &request );
/*4) indexes */
  for (i=0;i<MOA_in->elements_ub;i++) {
    lmessage = MOA_in->indexes[i];
    MPI_Isend(&lmessage, 1, MPI_LONG,
	    ProcDest, 4, MPI_COMM_WORLD, &request );
  }

}


void sendMOA (int ProcDest, MOA_rec * MOA_in) {
   unsigned int	membersize, maxsize;
    int			position;
    int			dest, tag;
    char		*buffer;
    long i, j, PrevCellsCount, lmessage;
    MPI_Request request; 
/* Send the MOA Buffer*/
  /*1) dimn */
    MPI_Isend(&MOA_in->dimn, 1, MPI_LONG,
	    1, 1, MPI_COMM_WORLD, &request );
  /*2) shape */
  for (i=0;i<MOA_in->dimn;i++) {
    lmessage = MOA_in->shape[i];
    MPI_Isend(&MOA_in->shape[i], 1, MPI_LONG,
	    1, 2, MPI_COMM_WORLD, &request );
  }
    /*3) elements_ub  */
    MPI_Isend(&MOA_in->elements_ub, 1, MPI_LONG,
	    1, 3, MPI_COMM_WORLD, &request );
/*4) indexes */
  for (i=0;i<MOA_in->elements_ub;i++) {
    lmessage = MOA_in->indexes[i];
    MPI_Isend(&lmessage, 1, MPI_LONG,
	    1, 4, MPI_COMM_WORLD, &request );
  }
/*5) elements */
  for (i=0;i<MOA_in->elements_ub;i++) {
/*5.1) values*/
    lmessage = MOA_in->elements[i].val;
    MPI_Isend(&lmessage, 1, MPI_LONG,
	    1, 5, MPI_COMM_WORLD, &request );

/*5.2) flags*/
    lmessage = MOA_in->elements[i].flag;
    MPI_Isend(&lmessage, 1, MPI_LONG,
	    1, 6, MPI_COMM_WORLD, &request );
/*5.3) Parent Cells Count*/
    lmessage = MOA_in->elements[i].prev_ub;
    MPI_Isend(&lmessage, 1, MPI_LONG,
	    1, 9, MPI_COMM_WORLD, &request );
/*5.4) Parent Cells Indices */
    if (MOA_in->elements[i].prev_ub > 0 ) {
      for (j=0;j<MOA_in->elements[i].prev_ub;j++) {
	lmessage = MOA_in->elements[i].prev[j];
	MPI_Isend(&lmessage, 1, MPI_LONG,
		  1, 10, MPI_COMM_WORLD, &request );
      }
    }
  }
}

int recivMOA (MOA_rec * MOA_in) {
    int             src;
    int             msgsize;
    MPI_Status      status;
    MPI_Request      request;
    char		*buffer;
   unsigned int	membersize, maxsize;
    int			dest, tag;
    int			position;
    long i, j, lmessage;

  /*1) dimn */
    MPI_Recv(&MOA_in->dimn, 1, MPI_LONG,
	    MPI_ANY_SOURCE , 1 , MPI_COMM_WORLD, &status);
  printf("\nI am [%d] received MOA_in->dimn = %d\n", myProcid, MOA_in->dimn);
  /*2) shape */
    MOA_in->shape = (long *) mcalloc(MOA_in->dimn, sizeof(long));
  for (i=0;i<MOA_in->dimn;i++) {
    MPI_Recv(&lmessage, 1, MPI_LONG,
	    MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);
    MOA_in->shape[i] = lmessage;
    printf("I am [%d] received MOA_in->shape[%d] = %d\n", myProcid, i, MOA_in->shape[i]);
  }
    /*3) elements_ub  */
    MPI_Recv(&MOA_in->elements_ub, 1, MPI_LONG,
	    MPI_ANY_SOURCE , 3 , MPI_COMM_WORLD, &status);

    printf("I am [%d] received elements_ub = %d\n", myProcid, MOA_in->elements_ub);
/*4) indexes */
    MOA_in->indexes = (signed long * ) mcalloc (MOA_in->elements_ub, sizeof(signed long));
  for (i=0;i<MOA_in->elements_ub;i++) {
    MPI_Recv(&lmessage, 1, MPI_LONG,
	    MPI_ANY_SOURCE, 4, MPI_COMM_WORLD, &status);
    MOA_in->indexes[i] = lmessage;
  }
/*5) elements */
    MOA_in->elements = (MOA_elm * ) mcalloc (MOA_in->elements_ub, sizeof(MOA_elm));
  for (i=0;i<MOA_in->elements_ub;i++) {
/*5.1) values*/
    MPI_Recv(&lmessage, 1, MPI_LONG,
	    MPI_ANY_SOURCE, 5, MPI_COMM_WORLD, &status);
    MOA_in->elements[i].val = lmessage;
    printf("I am [%d] received MOA_in->elements[%d].val = %d\n", myProcid, i, MOA_in->elements[i].val);

/*5.2) flags*/
    MPI_Recv(&lmessage, 1, MPI_LONG,
	    MPI_ANY_SOURCE, 6, MPI_COMM_WORLD, &status);
    MOA_in->elements[i].flag = lmessage;
    /*5.3) Parent Cells Count */
    MPI_Recv(&lmessage, 1, MPI_LONG,
	    MPI_ANY_SOURCE, 9, MPI_COMM_WORLD, &status);
    MOA_in->elements[i].prev_ub = lmessage;
    if (MOA_in->elements[i].prev_ub > 0 ) {
      /*5.4) Parent Cells Indices */
      MOA_in->elements[i].prev = (unsigned long * ) mcalloc (MOA_in->elements[i].prev_ub, sizeof(unsigned long));
      for (j=0;j<MOA_in->elements[i].prev_ub;j++) {
	MPI_Recv(&lmessage, 1, MPI_LONG,
		 MPI_ANY_SOURCE, 10, MPI_COMM_WORLD, &status);
	MOA_in->elements[i].prev[j] = lmessage;
      }
    }
    else
      MOA_in->elements[i].prev = NULL;
  }
}


void MasterProcess (int argc, char ** argv) {
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

  long flatIndex;
  long rslt_ub, strides; 
  int finished = 1;
  long * rslt = NULL;
  long * rsltInd = NULL;
  pthread_t PGProc, PartProc;
  Master * master;
  printf("In Master Process \n");
  processArguments(argc, argv, &seqNum, &sequences, &seqLen, &stype);
  createMOAStruct (&msaAlgn);
  
  createMOA(seqLen /* shape*/, seqNum /* dimension*/, msaAlgn /* MOA structure*/,-1,0);
  master = ( Master*)malloc (sizeof (Master));
  master->empty = 1;
  master->mut = (pthread_mutex_t *) malloc (sizeof (pthread_mutex_t));
  pthread_mutex_init (master->mut, NULL);
  master->notEmpty = (pthread_cond_t *) malloc (sizeof (pthread_cond_t));
  pthread_cond_init (master->notEmpty, NULL);

  if (master == NULL) return;
  master->msaAlgn = msaAlgn;
  pthread_create (&PGProc, NULL, ProcessorGateProcess, master);
  pthread_create (&PartProc, NULL, PartitioningProcess, master);
  /* Send MOA structure*/
  //sendMOAIndices (1, msaAlgn);
  pthread_join (PartProc, NULL);
  pthread_cancel (PGProc);
  deleteMOA (msaAlgn);  
}
void   Distribute (int argc, char ** argv) {
  int i;
 
  if (myProcid == 0) {
    MasterProcess(argc, argv);
  }
  else  {
    SlaveProcess ();
  }
}

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

  /* initialize MPI and get own id (rank) */
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &myProcid);
  MPI_Comm_size (MPI_COMM_WORLD, &ClsuterSize);

  Mode = Distributed;
  Distribute (argc, argv);
  
  MPI_Finalize ();
  return EXIT_SUCCESS;

}

void testcomm () {

  int msgsrc, i;/*Message Sender, and counters*/
  int msgval_r;/*Single Message value to be exchanged*/
  MPI_Status status; /*Message receive strcucture*/
  MPI_Request request;/*Non-blocking send or receive status*/
  int tag = 1;/*Message Exchange tag*/
  
  /*  
  for (i=0;i<ClsuterSize;i++) {
    MPI_Isend (&myProcid, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &request);
    printf("Iam %d sent my id to %d\n", myProcid, i);
  }
  for (i=0;i<ClsuterSize;i++) {
    MPI_Recv (&msgval_r, 1, MPI_INT, MPI_ANY_SOURCE, tag,MPI_COMM_WORLD, &status);
    msgsrc = status.MPI_SOURCE;   
    printf("Iam %d received %d from %d\n", myProcid, msgval_r, msgsrc);
  }
  */
}
