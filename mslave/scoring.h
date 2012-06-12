#ifndef SCORINGH_
#define SCORINGH_

#include "moa.h"

int subScore (char char1, char char2, int stype);
int gapScore(int stype);
long getRelation (long * cell, long * neighbor, long dimn, long * * decremented);
int getNeighborScores (char * * sequences, long * seqLen, long * m_index, MOA_rec * msaAlgn, int stype, long * LNCount, long * * lnScores, long * * lnIndices, long * * * pwScores, long * * neighbor, long * * decremented, long * * wm_index, MOA_rec * NghbMOA);
int PrintPrevChains (MOA_rec *  msaAlgn);
int getPrevCells(long findex, long score, long LNCount, long * lnScores, long * lnIndices, MOA_rec *  msaAlgn);
long getScore (void * threadarg, long findex, long * m_index, long * LNCount, long * * lnScores, long * * lnIndices, long * * * pwScores, long * * neighbor, long * * decremented, long * * wm_index, MOA_rec * NghbMOA );
void initTensor ( MOA_rec * msaAlgn, int stype);
void * IntermCellsThread (void * threadarg);
void *  DepProcThread (void * threadarg);
void * neighborsCache (void * threadarg);
void * computeScoreThread (void * threadarg);
void * ComputePartitionScores (void * threadarg);
void fillTensor (char * * sequences, long * seqLen, MOA_rec * msaAlgn, int stype);
int traceBack_one (long seqNum, char * * sequences, long * seqLen, MOA_rec * msaAlgn, int stype, char * * * algnseq);
void traceBack (long seqNum, char * * sequences, long * seqLen, MOA_rec * msaAlgn, int stype, char * * * * algnseq, long * * aSeqLen, int * alignmentsNo, long * currentCell, long * currentScore, int Prev_index);
void traceBack_loc (long seqNum, char * * sequences, long * seqLen, MOA_rec * msaAlgn, int stype, char * * * * algnseq, long * * aSeqLen, int * alignmentsNo, long * currentCell, long * currentScore, int Prev_index);
#endif
