#ifndef MOADSTH_
#define MOADSTH_

#include <mpi.h>

/* Communications Channels & Types*/
/* Data Size Dimn or Sequences number */
#define MPI_TAG_DSSType 0 
#define MPI_TYPE_DSSType MPI_INT 
/* Data Size Dimn or Sequences number */
#define MPI_TAG_DSDimn 1 
#define MPI_TYPE_DSDimn MPI_LONG 
/* Data Size Shape */
#define MPI_TAG_DSShape 2 
#define MPI_TYPE_DSShape MPI_LONG 
/* Partition Wave No */
#define MPI_TAG_PWave 3 
#define MPI_TYPE_PWave MPI_LONG 
/* Partition Dimension */
#define MPI_TAG_PDimn 4 
#define MPI_TYPE_PDimn MPI_LONG 
/* Partition Shape */
#define MPI_TAG_PShape 5 
#define MPI_TYPE_PShape MPI_LONG 
/* Partition Elements Upper Bound */
#define MPI_TAG_PElmUB 6 
#define MPI_TYPE_PElmUB MPI_UNSIGNED_LONG 
/* Partition Indices */
#define MPI_TAG_PIndices 7 
#define MPI_TYPE_PIndices MPI_UNSIGNED_LONG 
/* Partition Sequence Residues */
#define MPI_TAG_PResidue 8 
#define MPI_TYPE_PResidue MPI_CHAR 

/* Overlapping Cells Flag */
#define MPI_TAG_OCFlag 91 
#define MPI_TYPE_OCFlag MPI_INT 
/* Overlapping Cells Index */
#define MPI_TAG_OCCellIndex 92
#define MPI_TYPE_OCCellIndex MPI_UNSIGNED_LONG
/* Overlapping Cells Score */
#define MPI_TAG_OCCellScore 93
#define MPI_TYPE_OCCellScore MPI_LONG
/* Dependant Processors Flag */
#define MPI_TAG_DPFlag 94
#define MPI_TYPE_DPFlag MPI_INT 
/* Dependant Processors Cell */
#define MPI_TAG_DPCellIndex 95 
#define MPI_TYPE_DPCellIndex MPI_UNSIGNED_LONG 
/* Dependant Processors Partition Index */
#define MPI_TAG_DPPartIndex 96
#define MPI_TYPE_DPPartIndex MPI_UNSIGNED_LONG 
/* Dependant Processors Wave No */
#define MPI_TAG_DPWave 97
#define MPI_TYPE_DPWave MPI_LONG 
/* Dependant Processors Proc ID */
#define MPI_TAG_DPProcID 98
#define MPI_TYPE_DPProcID MPI_INT 
/* ComputationPhase */
#define MPI_TAG_CompPhase 99
#define MPI_TYPE_CompPhase MPI_INT 
/* processorsToEnq */
#define MPI_TAG_ProcEnq 100
#define MPI_TYPE_ProcEnq MPI_INT

MPI_Comm MOAMSA_COMM_WORLD;

void checkMutexErrorCode (int caller, int errCode);
void checkMPIErrorCode (int caller, int errmsg);

int NBReceive (void * buffer, int count, MPI_Datatype type, int source, int tag, MPI_Request * request, MPI_Status * status);
int NBSend (void * buffer, int count, MPI_Datatype type, int destination, int tag, MPI_Request * request);

#endif

