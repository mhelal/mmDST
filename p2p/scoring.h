#ifndef SCORINGH_
#define SCORINGH_

#include "moa.h"

long getRelation (ScoringData * sData, MOATypeInd neighbIndex);
int getNeighborScores (ProcessData * pData, ScoringData * sData, WavesData * wData);
int PrintPrevChains (MOA_rec *  msaAlgn);
MOATypeElmVal initLBCell (MOATypeShape * lm_index, MOATypeDimn dimn, int stype);
int  getCellScore(ProcessData * pData, ScoringData * sData, WavesData * wData, MOATypeShape * cellIndex, MOATypeElmVal * score, int * inSearchSpace, int NeighborSearch, MOATypeInd NeighbIndex);
MOATypeElmVal getCellDPScore (ProcessData * pData, ScoringData * sData, WavesData * wData);
MOATypeElmVal getCellSPScore (ProcessData * pData, ScoringData * sData, WavesData * wData);
MOATypeElmVal getNeghbScore (ScoringData * sData, MOATypeInd neighbIndex);
#endif
