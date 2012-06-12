/**********************************************************
* Author: Manal Helal																			*
* Last Modification Date: Fri 12 Jan 2007 03:39:51 AM EST *
* Project : MMSA - Multiple Sequence Alignment Based on 	*
* 					Mathematics of Arrays - PhD Experimentation		*
* File: master.h, header for master process functions			*
***********************************************************/
#ifndef MAINH_
#define MAINH_


#include <mpi.h>
#include "globals.h"
#include "moa.h"

/* New Overlapping Cells categorized by waves, to reduce search space at every wave*/
typedef struct WaveOCO_Tag {
    int  sent; /* Once I send OCO to all dependant processors, I delete them, I don't need to keep them in the buffer*/
    MOATypeShape * cellIndex;  /* Overlapping Cell Index */
    MOATypeElmVal cellScore;  /* Overlapping Cell Score */
    int  depProc_ub; /* dep Proc count - increase while I receive dependency from master*/
    int  * depProc;  /* dep Proc added here*/
} WaveOCOType;

typedef struct OCO_Tag {
	long wavesOC; /*Total OCout in current indexed wave*/ 
	WaveOCOType * WOCO; /*Waves OC list - as before*/
} OCOType;

typedef struct WaveOCI_Tag {
	MOATypeShape * cellIndex; /* Overlapping Cell Index */
	MOATypeElmVal cellScore; /* Overlapping Cell Score */
} WaveOCIType;

typedef struct OCI_Tag {
	long wavesOC; /*Total ICin in current indexed wave*/ 
	WaveOCIType * WOCI;
} OCIType;

typedef struct Proc_rec {
    int stype;
    volatile int compFinished;
    MOATypeDimn seqNum;
    char * * seqName;
    MOATypeShape * seqLen;
    char * * sequences;
    MOA_rec * msaAlgn;
    //char * * currPartSeq;
    long globalWaveNo, waveNo, partitionsCount, partNo, computedPartitions, sendOC, sendOCPart, partitionSize;
    int * validPartitions;
    /*Communication Buffer*/
    MOATypeDimn commBufSize;
    MPI_Request * mpi_requests;
    int sendRequests;
    long  * buffer_size;
    MOATypeElmVal * * buffer;
    int depProcCount, * proc_index;
    /*Overlapping Cells outgoing*/
    OCOType * OCout;
    /*Overlapping Cells incoming*/
    OCIType * OCin;
} ProcessData;

typedef struct ScoringDataType_rec { /* ScoringData */
    MOATypeDimn seqNum;
    MOATypeShape * seqLen;
    char * * sequences;
    MOA_rec * msaAlgn; /* The scoring tensor or a partition of it*/ 
    MOA_rec * NghbMOA; /* The lower neighbors window*/
    MOATypeElmVal score;
    MOATypeShape * p_index; /* global multidimensional partition index-First element in the partition*/
    MOATypeShape * lm_index; /* local multidimensional index in the current partition*/
    MOATypeShape * gm_index; /*global multidimensional index from the whole tensor*/
    MOATypeShape * posDimn; /* dimensions where the cell is overlapping on to get its enclosing partition*/
    MOATypeShape * neighbor; /*current neighbor global Index, that will read the neighbor local Index if exists locally*/
    MOATypeDimn LNCount, CalLnCount;
    long waveNo, partNo;
    MOATypeInd findex; /* Current Local Flat Index*/
    int * lnInSearchSpace; /* Lower neighbor In Search Space 1 = true, 0 = false */
    MOATypeElmVal * * pwScores; /* holds the scores of pairwise match / mistmatch scores*/
    MOATypeShape * decremented; /* list of decremented indices in the multidimensional index of the current neighbor*/
    MOATypeShape * depPart; /*global index of partition index of local lower border cell that need to be received*/
    int stype, depPartOutsideEpsilon;
} ScoringData;

typedef struct reg_rec {
    double maxColScore; /* Maximum Alignment Columns SP Score up to the number of seeds*/
    double minColScore; /* Minimum Alignment Columns SP Score up to the number of seeds*/
    MOATypeInd maxStartColScore; /* Maximum Alignment Columns SP Score up to the number of seeds*/
    MOATypeInd maxEndColScore; /* Maximum Alignment Columns SP Score up to the number of seeds*/
    MOATypeInd minStartColScore; /* Minimum Alignment Columns SP Score up to the number of seeds*/
    MOATypeInd minEndColScore; /* Minimum Alignment Columns SP Score up to the number of seeds*/
} RegionData;

typedef struct tb_rec {
    MOATypeDimn seqNum;
    int currProc; /* The current trace back processor*/
    MOATypeElmVal maxCellScore; /* The current maximum score*/
    /*allocated in seqNum in initTBData in utils.c*/
    MOATypeShape * maxCellIndex; /*Global Index pf the current Maximum Score tp resume tracing from*/
    MOATypeShape * partIndex; /*The current trace back partition Index*/
    /*initialized to zeros and NULLS in initTBData in utils.c, then allocated as partial paths are found*/
    long pathParts; /* No of partial paths found in processors so far (in master, total all over processors, otherwise, total of local partial paths, which is not needed as yet, since only the last partial part will be processed at a time, and sent to master for assembly, and no need to keep a record for it, so in slaves, it is always 1, and overwritten if more is found.)*/
    char * * * algnseq; /*all partial paths, again in slave only one*/
    MOATypeShape * aSeqLen; /*length of each partial path*/
    
    /*Allocated only after end of distrbuted trace back to concatenate partial paths*/
    MOATypeShape comseqLen; /*The complete sequence alignement Length*/
    char  * * completePath; /* The complete sequence Alignment - concatenation of all partial alignments*/
    MOATypeElmVal * alignColScore; /* Alignment Columns SP Score*/
    MOATypeElmVal alignSPScore; /*Total SP score for all columns - Alignment Score*/
    long double a, c, g, t, indel;
    long double r, n, d, q, e, h, i, l, k, m, f, p, s, w, y, v, b, z, x; /*remaining amino acids*/
    long double alignShannonEntropy; /*Shannon Entropy of the alignment produced*/
    long double meanSPColScore;
    long double * * alignDistaneMeasures; /*Pair wise distance measures matrix*/
    long double * alignGapScore; /*the Gap count in each sequence divided by the square root of each length*/
    RegionData regD; /*maximum and minimum regions up to the specified number of seeds*/
}TracebackData;

typedef struct w_rec {
    /* Temporary variables to print all parts in waves whether valid or not starting with All*/
    MOATypeDimn seqNum;
    MOATypeShape * seqLen;
    MOATypeShape * * waveMiddle; /*current wave Middle partition Index, fair distribution of wave length over dimensions starting from the middle index*/
    long double * waveLength;
    long * AllpartsInWave, AllPartOrder;
    long partitionSize, calcParts, partsTotal, wavesTotal;
    long * partsInWave; /*Number of Parts in each wave*/
    MOATypeInd duplicatesTotal; 
    MOATypeShape * * * AllpartsInWaveIndices, * * * partsInWaveIndices; /*multidimensional Index of part index of each part in each wave. 
     * first dimension is waves, second dimensions is part order in this wave, third dimension is 
     * composite of the k dimensionality of the dataset*/
}WavesData; /*Waves Partitioning Data*/


/* Trace back functions*/
/* Master functions*/
void assemblePathParts (TracebackData * tbData);

void distributedMasterTraceBack (ProcessData * pData, WavesData * wData, TracebackData * tbData);

void getmaxCellScore (ProcessData * pData, TracebackData * tbData);
void getLocalMaxScore (ProcessData * pData, WavesData * wData, TracebackData * tbData) ;

MOATypeElmVal getLocalMaxCellScore (ProcessData * pData, WavesData * wData, TracebackData * tbData, int searchAllPartitions);

void * tbMaster (ProcessData * pData, WavesData * wData);
void ExitProcess (ProcessData * pData, ScoringData * sData, WavesData * wData);

/* Slave functions*/

void * tbSlave (ProcessData * pData, WavesData * wData);

MOATypeShape traceBack (ProcessData * pData, WavesData * wData, TracebackData * tbData);

void distributedSlaveTraceBack (ProcessData * pData, WavesData * wData, TracebackData * tbData);
void sendmaxCellScore (ProcessData * pData, WavesData * wData, TracebackData * tbData);

/* Checkpointing functions*/
int restoreWavesCalculations (ProcessData * pData, WavesData * wData);
int checkPointWavesCalculations (ProcessData * pData, WavesData * wData);

int restorePartitionCheckPoint (ProcessData * pData, WavesData * wData, long waveNo, long partNo);
int checkPointPartition (ProcessData * pData, ScoringData * sData);

int restoreCheckPoint (ProcessData * pData, WavesData * wData);
int checkPoint (ProcessData * pData, ScoringData * sData);
int checkPointOC (ProcessData * pData, long waveNo);
int restoreOCCheckPoint (ProcessData * pData, long waveNo);



/*Scoring Functions*/
void MainProcess (MOATypeDimn seqNum, char * * sequences, char * * seqName, MOATypeShape * seqLen, int stype, long partitionSize);

int initProcessMemory (ProcessData * * pData, ScoringData * * sData, WavesData * * wData, MOATypeDimn seqNum, MOATypeShape * seqLen, char * * sequences, char * * seqName, int stype, long partitionSize);

void freeProcessMemory (ProcessData * * pData, ScoringData * * sData, WavesData * * wData);
int PrintPrevChains (MOA_rec *  msaAlgn);

/* Dynamic Programming Scoring Function*/
void * DPComputeScores (ProcessData * pData, ScoringData * sData, WavesData * wData);
/* Sum of Pairs Scoring Function*/
void * SPComputeScores (ProcessData * pData, ScoringData * sData, WavesData * wData);

void * ScoreCompThread (ProcessData * pData, ScoringData * ScoringData, WavesData * wData);

/*HyperPlane Functions*/

long double  getPythagoreanOriginDistance (MOATypeShape * point, MOATypeDimn dimn);
long double  getPythagoreanDistance (MOATypeShape * x_point, MOATypeShape * y_point, MOATypeDimn dimn);
long double  getPythagoreanOriginDirection (MOATypeShape * point, MOATypeDimn dimn);
long double  getPythagoreanDirection (MOATypeShape * x_point, MOATypeShape * y_point, MOATypeDimn dimn);
long double  getMoADistance (MOATypeShape * x_point, MOATypeShape * y_point, MOATypeDimn dimn);
long double  getVectorsAngle (MOATypeShape * x_point, MOATypeShape * y_point, MOATypeDimn dimn);
long double  getVectorsRankedAngle (MOATypeShape * x_point, MOATypeShape * y_point, MOATypeDimn dimn);
MOATypeShape isVectorInHyperPlane (MOATypeShape * x_point, MOATypeShape * y_point, MOATypeDimn dimn, MOATypeShape * z_point);
double getWaveLength (MOATypeShape * shape, MOATypeDimn dimn, long partitionSize, long waveNo);
void getFirstPartitionInWave (MOATypeShape * shape, MOATypeDimn dimn, long partitionSize, long waveNo, MOATypeShape * * partIndex);
void getSecondPartitionInWave (MOATypeShape * shape, MOATypeDimn dimn, long partitionSize, long waveNo, MOATypeShape * firstPartition, MOATypeShape * * partIndex);
void getLastPartitionInWave (MOATypeShape * shape, MOATypeDimn dimn, long partitionSize, long waveNo, MOATypeShape * * partIndex);
void getMiddlePartitionInWave (MOATypeShape * shape, MOATypeDimn dimn, long partitionSize, long waveNo, MOATypeShape * * partIndex);

/* Partitioning Functions*/
int getWavesPartsDuplicates(MOATypeDimn dimn, MOATypeShape * shape, long pSize, long * wTotal, long * pTotal, MOATypeInd * dTotal);

int getProcID (WavesData * wData, long waveNo, long partNo);

void getNextPartition (WavesData * wData, long * waveNo, long * partNo);
void getPrevPartition (WavesData * wData, long * waveNo, long * partNo);
long getWavePartsTotal (long WaveNo, MOATypeDimn dimn);

long double distanceFromMiddle(WavesData * wData, MOATypeShape * partIndex, long waveNo);
//void getNextPIndex(int * more, MOATypeDimn dimn, MOATypeDimn act_dimn, long waveNo, MOATypeShape * * PIndex);

void distributeRemDistance (MOATypeShape * * PIndex, MOATypeDimn * middleDimn, MOATypeDimn dimn, MOATypeDimn startDimn, MOATypeInd startDist, MOATypeInd remdist, MOATypeDimn remDimn);

void getNextPIndex_middle(int * more, long waveNo, MOATypeDimn dimn, MOATypeDimn * middleDimn, MOATypeShape * * PIndex);

void getNextPIndex_ordered(int * more, MOATypeDimn dimn, MOATypeDimn act_dimn, long waveNo, MOATypeShape * * PIndex);

int notPreviouslyVisited (WavesData * wData, long partsinCurrentWave, long waveNo, MOATypeShape * ind);

int notPreviouslyVisited_dbg (WavesData * wData, long waveNo, MOATypeShape * ind); 

int addPartitionIndex (ProcessData * pData, WavesData * wData, long * PartOrder, MOATypeShape * dist, long waveNo, int * morePartitions);
long getPIndicesinWave (ProcessData * pData, WavesData * wData, int * more, long waveNo);

long calcWaves (ProcessData * pData, WavesData * wData);

int getPartition (MOATypeShape * startIndex, MOATypeDimn dimn, MOATypeShape * shape, MOA_rec * * part, long partSize);

int IsCellInPart (MOATypeShape * cellIndex, MOATypeShape * partIndex, MOATypeDimn dimn, MOATypeShape * shape, long partSize);
int getLocalIndex (MOATypeShape * cellIndex, MOATypeShape * partIndex, MOATypeDimn dimn, MOATypeShape * shape, long partSize, MOATypeShape * * localIndex);
MOATypeInd getLocalFlatIndex (MOATypeShape * globalIndex, MOA_rec * msaAlgn, MOATypeShape * shape, long partSize);

long  getCellPartIndex (MOATypeInd cellIndex, MOATypeDimn dimn, MOATypeShape * shape, long partSize, MOATypeInd * *  partIndex);

long * getCellAllPartIndex (MOATypeInd cellIndex, MOATypeDimn dimn, MOATypeShape * shape, long partSize, MOATypeInd * * * partIndex, long * parts_ub);

int  getPositionalPartitionIndex (int flag, MOATypeShape * posDimn, MOATypeShape * cellIndex, MOATypeDimn dimn, MOATypeShape * shape, long partitionSize, MOATypeShape * * partIndex);
int  getPartitionIndex (MOATypeShape * cellIndex, MOATypeDimn dimn, MOATypeShape * shape, long partitionSize, MOATypeShape * * partIndex);
int getPartitionDetails (ProcessData * pData, WavesData * wData, MOATypeShape * * partIndex, MOATypeShape * cellIndex, long * waveNo, long * partNo, int * procID);

int isPartInSearchSpace(MOATypeShape * mPartIndex, WavesData * wData);
long isCellInSearchSpace(MOATypeShape * * cellIndex, WavesData * wData, long * waveNo);

int isPartitionInWave (WavesData * wData, long waveNo, MOATypeShape * partIndex);
long getPartIndex (ProcessData * pData, MOATypeInd globalIndex, MOATypeInd * localCellIndex, MOATypeInd startPart);
long getWaveNo (MOATypeDimn seqNum, long partitionSize, MOATypeShape * partIndex);
long getPartNo (WavesData * wData, MOATypeShape * partIndex, long waveNo);
long getPartitionPosition (WavesData * wData, MOATypeShape * partIndex, long * waveNo, long * partNo);

/* Dependency Functions*/

int receiveOC (ProcessData * pData, ScoringData * sData);

int prepareWaveOC (long waveNo, ProcessData * pData);
int sendOCtoHigherProcessors (ProcessData * pData);
int sendOCtoLowerProcessors (ProcessData * pData);

int addOC (ProcessData * pData, ScoringData * sData, WavesData * wData);

int checkPrevPartitions (ProcessData * pData, MOATypeShape * cellIndex, MOATypeElmVal * score);
int checkRecvOC (ProcessData * pData, WavesData * wData, MOATypeShape * cellIndex, MOATypeElmVal * score, MOATypeInd * startIndex);


#endif
