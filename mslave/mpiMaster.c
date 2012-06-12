#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>      /* include the character library */
#include <pthread.h>
#include <mpi.h>
#include <sys/time.h>
#include "master.h"
#include "moaDst.h"
#include "utils.h"

void SendPartition (Master * master, int toProcessor, MOA_rec * MOA_partition, long waveNo) {
  MPI_Request request;/*Non-blocking send or receive status*/
  int ComputationPhase = 1;
  int MPI_return;
	char msg[MID_MESSAGE_SIZE];
	
  MPI_return = NBSend (&ComputationPhase, 1, MPI_TYPE_CompPhase, toProcessor, MPI_TAG_CompPhase, &request);
  checkMPIErrorCode (3, MPI_return);
  sprintf (msg,"Master sent ComputationPhase %d\n", ComputationPhase);
	mprintf (9, msg, 3);	
  sendMOAIndices (master, toProcessor, MOA_partition, waveNo); 
  mprintf(9,"Master sent MOA Partition\n", 3);

} 


int sendMOAIndices (Master * master, int ProcDest, MOA_rec * MOA_in, long waveNo) {
  unsigned int	membersize, maxsize;
  int	position;
  int	dest, tag;
  char *buffer, cmessage, msg[SHORT_MESSAGE_SIZE];
  long i, j, PrevCellsCount, lmessage;
  unsigned long ulmessage;
  long * ind = NULL;
  long * ind2 = NULL;
  int MPI_return;
  MPI_Request request; 
  
  /* Send the MOA Buffer*/
  /*3) Wave No */
  MPI_return = NBSend(&waveNo, 1, MPI_TYPE_PWave, ProcDest, MPI_TAG_PWave, &request );
  checkMPIErrorCode (3, MPI_return);

  /*3) dimn */
  MPI_return = NBSend(&MOA_in->dimn, 1, MPI_TYPE_PDimn, ProcDest, MPI_TAG_PDimn, &request );
  checkMPIErrorCode (3, MPI_return);

  /*4) shape */
  for (i=0;i<MOA_in->dimn;i++) {
    lmessage = MOA_in->shape[i];
    MPI_return = NBSend(&MOA_in->shape[i], 1, MPI_TYPE_PShape, ProcDest, MPI_TAG_PShape, &request );
  checkMPIErrorCode (3, MPI_return);
  }
  //5) elements_ub  
  ulmessage = MOA_in->elements_ub;
  MPI_return = NBSend(&ulmessage, 1, MPI_TYPE_PElmUB, ProcDest, MPI_TAG_PElmUB, &request );
  checkMPIErrorCode (3, MPI_return);
  //6) indexes 
  sprintf (msg,"Master sending %ld elements_ub\n", MOA_in->elements_ub);
	mprintf(12, msg, 3);
  ind = (long *) mmalloc (sizeof(long) * master->msaAlgn->dimn);
  for (i=0;i<MOA_in->elements_ub;i++) {
    ulmessage = MOA_in->indexes[i];
    sprintf (msg,"Master sending index[%ld] %ld\n", i, MOA_in->indexes[i]);
		mprintf(12, msg, 3);
    MPI_return = NBSend(&ulmessage, 1, MPI_TYPE_PIndices, ProcDest, MPI_TAG_PIndices, &request );
	  checkMPIErrorCode (3, MPI_return);
  
    sprintf (msg, "Master finished sending indices with i = %ld\n", i);
		mprintf(12, msg, 3);
   }
    // 7) Send Corresponding Residues from the Sequences read
    //Gamma_Inverse(i, MOA_in->shape, MOA_in->dimn, ind);
    //if (isLowerBorderCell(ind, MOA_in->dimn) == 0) {
      //Gamma_Inverse(MOA_in->indexes[i], master->msaAlgn->shape, master->msaAlgn->dimn, ind);
  ind2 = (long *) mmalloc (sizeof(long) * master->msaAlgn->dimn);

    for (i=0;i<MOA_in->dimn;i++) 
			ind [i] = 0;
    Gamma_Inverse(MOA_in->indexes[0], master->msaAlgn->shape, master->msaAlgn->dimn, ind2);
    for (i=0;i<MOA_in->dimn;i++) {
			cmessage = master->sequences[i][ind2[i]];
			MPI_return = NBSend(&cmessage, 1, MPI_TYPE_PResidue, ProcDest, MPI_TAG_PResidue, &request );
			checkMPIErrorCode (3, MPI_return);
			sprintf(msg, "sent sequences[%ld][%ld] = %c  to Processor = %d\n", i, ind2[i], cmessage, ProcDest);
			mprintf(12, msg, 3);
     }
				
  for (i=0;i< MOA_in->dimn;i++) {
	  for (j=1;j< MOA_in->shape[i];j++) {
	  	ind[i] = j;
  		lmessage = Gamma (ind, MOA_in->dimn, MOA_in->shape, MOA_in->dimn, 1);
    Gamma_Inverse(MOA_in->indexes[lmessage], master->msaAlgn->shape, master->msaAlgn->dimn, ind2);

			cmessage = master->sequences[i][ind2[i]];
			MPI_return = NBSend(&cmessage, 1, MPI_TYPE_PResidue, ProcDest, MPI_TAG_PResidue, &request );
			checkMPIErrorCode (3, MPI_return);
			sprintf(msg, "sent sequences[%ld][%ld] = %c  to Processor = %d\n", i, ind2[j], cmessage, ProcDest);
			mprintf(12, msg, 3);
			ind[i]=0;
     }
	}
  
  if (ind != NULL)
	  free (ind);    
  if (ind2 != NULL)
	  free (ind2);    
return 0;
}



int ReceiveSlavetoProcGate (int * ComputationPhase) {
  MPI_Status status; /*Message receive strcucture*/
  MPI_Request request;
  int MPI_return, flag, procToEnq;
  struct timeval timeout;  
	char msg[SHORT_MESSAGE_SIZE];

  MPI_return = MPI_Irecv (ComputationPhase, 1,  MPI_TYPE_ProcEnq, MPI_ANY_SOURCE, MPI_TAG_ProcEnq, MOAMSA_COMM_WORLD, &request);
  checkMPIErrorCode (2, MPI_return);
  MPI_return = MPI_Test (&request, &flag, &status);
  checkMPIErrorCode (2, MPI_return);
    
  while(flag == 0) {
    timeout.tv_sec = 1;
    timeout.tv_usec = 0;
    select(10, NULL, NULL, NULL, &timeout);
    MPI_return = MPI_Test (&request, &flag, &status);
    checkMPIErrorCode (2, MPI_return);
  }
  procToEnq = status.MPI_SOURCE;
  sprintf (msg,"Master ProcGate receieved Slave %d with CompPhase %d\n", procToEnq, (*ComputationPhase));  	
	mprintf(12, msg, 2);
	return procToEnq;
}

int sendCellDepProc (long waveNo, int toProc, unsigned long cellIndex, unsigned long partIndex, int depProc) {
  int MPI_return, DepProcflag;
	char msg[SHORT_MESSAGE_SIZE];
  MPI_Request request;/*Non-blocking send or receive status*/
  DepProcflag = 1;

  MPI_return = NBSend (&DepProcflag, 1, MPI_TYPE_DPFlag, toProc, MPI_TAG_DPFlag, &request);
  checkMPIErrorCode (3, MPI_return);

  MPI_return = NBSend (&cellIndex, 1, MPI_TYPE_DPCellIndex, toProc, MPI_TAG_DPCellIndex, &request);
  checkMPIErrorCode (3, MPI_return);

  MPI_return = NBSend (&partIndex, 1, MPI_TYPE_DPPartIndex, toProc, MPI_TAG_DPPartIndex, &request);
  checkMPIErrorCode (3, MPI_return);

  MPI_return = NBSend (&waveNo, 1, MPI_TYPE_DPWave, toProc, MPI_TAG_DPWave, &request);
  checkMPIErrorCode (3, MPI_return);

  MPI_return = NBSend (&depProc, 1, MPI_TYPE_DPProcID, toProc, MPI_TAG_DPProcID, &request);
  checkMPIErrorCode (3, MPI_return);

  sprintf (msg,"Master sent cellIndex %ld in part %ld & depProc %d in wave %ld\n", cellIndex, partIndex, depProc, waveNo);
	mprintf(12, msg, 3);
  return MPI_return;
}

int SendEndPartitioingSignal() {
  int i, MPI_return;
  MPI_Request request;
	int ComputationPhase = 2;
  for (i=1;i<ClusterSize;i++) {
    MPI_return = NBSend (&ComputationPhase, 1, MPI_TYPE_DPFlag, i, MPI_TAG_DPFlag, &request); /*to stop the DepProc Thread in the Slaves*/
    checkMPIErrorCode (3, MPI_return);
    MPI_return = NBSend (&ComputationPhase, 1, MPI_TYPE_CompPhase, i, MPI_TAG_CompPhase, &request); /* to stop the Compute Partition Score Thread in the slave*/
  checkMPIErrorCode (3, MPI_return);
  }
return 0;
}


