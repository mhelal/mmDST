#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "../moaDst.h"   
#include "../main.h"
#include "../globals.h"
#include "../moa.h"
#include "../utils.h"
#include "../scores.h"
 

int myProcid, ClusterSize;
MPI_Comm MOAMSA_COMM_WORLD;


MOATypeShape processTBStep (ProcessData * pData, TracebackData * tbData, MOATypeShape * * old_index, MOATypeShape * * new_index) {
    MOATypeDimn i;
    int allGaps;
#ifndef NDEBUG
    char msg[MID_MESSAGE_SIZE];
#endif

    allGaps = 1;
    for (i=0;i<tbData->seqNum;i++) {
        if (tbData->aSeqLen[tbData->pathParts-1] == 0)
            tbData->algnseq[tbData->pathParts-1][i] = mmalloc (2 * sizeof *tbData->algnseq[tbData->pathParts-1][i]);
        else 
            tbData->algnseq[tbData->pathParts-1][i] = realloc (tbData->algnseq[tbData->pathParts-1][i], (tbData->aSeqLen[tbData->pathParts-1]+2) * sizeof *(tbData->algnseq[tbData->pathParts-1][i]));

        if (tbData->algnseq[tbData->pathParts-1][i] == NULL ) {
            printf("[%d]Can not allocate memory to Alligned Sequence %ld\n", myProcid, i);
            return -1;
        }

        if (((*old_index)[i] > (*new_index)[i]))  {
            tbData->algnseq[tbData->pathParts-1][i][tbData->aSeqLen[tbData->pathParts-1]] = pData->sequences[i][(*old_index)[i]];
            allGaps = 0;
        }
        else 
            tbData->algnseq[tbData->pathParts-1][i][tbData->aSeqLen[tbData->pathParts-1]] = GAPCHAR;
#ifndef NDEBUG
        sprintf(msg, "old_index[%ld] = %ld new_index = %ld is %c\n", i,  (*old_index)[i], (*new_index)[i],tbData->algnseq[tbData->pathParts-1][i][tbData->aSeqLen[tbData->pathParts-1]]);
        mprintf(20, msg, 1);
#endif
        (*old_index)[i] = (*new_index)[i];
    }
    if (allGaps == 0)
        tbData->aSeqLen[tbData->pathParts-1] ++;
    return tbData->aSeqLen[tbData->pathParts-1];
}

MOATypeInd locateCell (ProcessData * pData, WavesData * wData, TracebackData * tbData, MOATypeShape * * cellIndex, MOATypeShape * * partIndex, int * currProc) {
    long currWaveNo = pData->waveNo, currPartNo = pData->partNo;
    MOATypeInd currLocalIndex = -1; /*  initialize in case not found, to get out of the local trace back*/                                        
    MOATypeDimn k;
    
    if (Mode == Distributed) {
        if  (getPartitionIndex ((*cellIndex), tbData->seqNum, pData->seqLen, wData->partitionSize, partIndex) == 0) {
           if (isPartInSearchSpace((*partIndex), wData) == 0) {
                getPartitionPosition (wData, (*partIndex), &pData->waveNo, &pData->partNo);
                if ((currWaveNo != pData->waveNo) || (currPartNo != pData->partNo)) {
                    (*currProc) = getProcID (wData, pData->waveNo, pData->partNo);
                    if (myProcid ==(*currProc)) {                                                 
                        if (restorePartitionCheckPoint(pData, wData, pData->waveNo, pData->partNo) != 0) {
                            printf ("[%d] Error Retrieving partition file. Exitiing.\n", myProcid);
                           return -1;
                        }        
                        currLocalIndex = getLocalFlatIndex ((*cellIndex), pData->msaAlgn, pData->seqLen, pData->partitionSize);
                        if (currLocalIndex == -1)
                            currLocalIndex = 0; /* to get out of the local trace back*/
                    }
                } /*End if not new partition, cell is in same partition*/
                else {
                    currLocalIndex = getLocalFlatIndex ((*cellIndex), pData->msaAlgn, pData->seqLen, pData->partitionSize);                    
                }
            }/*End looking for main partition Index in search space*/
            else {
                /*if part is not in search space, then cellIndex is used from another side-partition. search all partitions in the waves were the cell belong*/
                if ((pData->partNo = isCellInSearchSpace(cellIndex, wData, &pData->waveNo)) >= 0) {
                    for (k=0;k<tbData->seqNum;k++) 
                        (*partIndex)[k] = wData->partsInWaveIndices[pData->waveNo][pData->partNo][k];
                    (*currProc) = getProcID (wData, pData->waveNo, pData->partNo);
                    if (myProcid ==(*currProc)) {                                                 
                        if (restorePartitionCheckPoint(pData, wData, pData->waveNo, pData->partNo) != 0) {
                            printf ("[%d] Error Retrieving partition file. Exitiing.\n", myProcid);
                            return -1;
                        }           
                        currLocalIndex = getLocalFlatIndex ((*cellIndex), pData->msaAlgn, pData->seqLen, pData->partitionSize);
                    }                        
                }                
            } /*End looking for cell in other partitions in search space*/                    
        } /*End getting main partition index*/        
    }
    /* Sequential Mode*/
    else {
        currLocalIndex = Gamma((*cellIndex), pData->msaAlgn->dimn, pData->msaAlgn->shape, pData->msaAlgn->dimn, 1);
        (*currProc) = 0;
    }
    return currLocalIndex;   
}

MOATypeInd checkCellNeighbors (ProcessData * pData, WavesData * wData, TracebackData * tbData, MOATypeShape * * cellIndex, MOATypeElmVal * currentScore) {
    MOATypeInd j, currLocalIndex, maxNghbIndex;
    MOATypeDimn k;
    int maxNghbProc, currProc;
    long waveNo, partNo;
    MOATypeElmVal maxNghbScore;
    MOA_rec * NghbMOA = NULL;
    MOATypeShape * currentIndex = NULL; /* multidimensional current global index*/
    MOATypeShape * partIndex = NULL;
#ifndef NDEBUG
    char msg[MID_MESSAGE_SIZE];
#endif
    currProc = tbData->currProc;
    partIndex = mcalloc ((MOATypeInd) tbData->seqNum, (MOATypeInd) sizeof *partIndex);
    if (partIndex == NULL) {
        printf ("[%d] Could not allocate memopry for part Index in checkCellNeighbors. Exiting\n", myProcid);
        return -1;
    }
    currentIndex = mcalloc ((MOATypeInd) pData->seqNum, (MOATypeInd) sizeof *currentIndex);
    if (currentIndex == NULL) {
        printf ("[%d] Could not allocate memopry for currentIndex in checkCellNeighbors. Exiting\n", myProcid);
        return -1;
    }
    for (k=0;k<pData->seqNum;k++) 
        currentIndex[k] = (*cellIndex)[k];
    createMOAStruct (&NghbMOA);
    maxNghbIndex = currLocalIndex = -1;/* initialize in case not found, to get out of the local trace back*/
    if (getLowerNeighbors (2, currentIndex, pData->seqNum, pData->seqLen, &NghbMOA) == 0) {
        /* loop throught neighbors */                
        for (j=0;j<NghbMOA->elements_ub - 1;j++) {                    
#ifndef NDEBUG
            printf ("checkCellNeighbors: before 1 locateCell currProc = %d \n", currProc);
            fflush(stdout);
#endif
            currLocalIndex = locateCell (pData, wData, tbData, &NghbMOA->indexes[j], &partIndex, &currProc);
#ifndef NDEBUG
            printf ("checkCellNeighbors: after 1 locateCell currProc = %d currLocalIndex %lld maxNghbIndex %lld \n", currProc, currLocalIndex, maxNghbIndex);
            fflush(stdout);
#endif
            if (currLocalIndex >= 0) {
                (*currentScore) = pData->msaAlgn->elements[currLocalIndex].val;         
                if (maxNghbIndex == -1) {
                    maxNghbScore = (*currentScore);
                    maxNghbIndex = j;
                    maxNghbProc = currProc;
                }
                else if ((*currentScore) > maxNghbScore) {
                    maxNghbScore = (*currentScore);
                    maxNghbIndex = j;
                    maxNghbProc = currProc;
                }
            }
        }
        if (maxNghbIndex >= 0) {
            (*currentScore) = maxNghbScore;
            tbData->currProc = maxNghbProc;
            for (k=0;k<pData->seqNum;k++) 
                currentIndex[k] = NghbMOA->indexes[maxNghbIndex][k];
            currLocalIndex = locateCell (pData, wData, tbData, &currentIndex, &partIndex, &currProc);
            if (currLocalIndex >= 0) 
                (*currentScore) = pData->msaAlgn->elements[currLocalIndex].val;
        }
    }
    for (k=0;k<pData->seqNum;k++) 
        (*cellIndex)[k] = currentIndex[k];

    if (NghbMOA != NULL)
        deleteMOA (NghbMOA);           
    if (currentIndex != NULL)
        free (currentIndex);
    currentIndex = NULL;    
    if (partIndex != NULL) 
        free (partIndex);
    partIndex = NULL;
    if (maxNghbIndex == -1)
        tbData->currProc = currProc;
    return currLocalIndex;   
}
MOATypeShape traceBack (ProcessData * pData, WavesData * wData, TracebackData * tbData) {
    MOATypeInd currLocalIndex, j;
    MOATypeDimn k;
    MOATypeElmVal prevScore, currentScore;
    long OCIwaveNo, OCIndex;
    MOATypeShape * localIndex = NULL; /* multidimensional current global index*/
    MOATypeShape * globalIndex = NULL; /* multidimensional current global index*/
    MOATypeShape * old_globalIndex = NULL; /* multidimensional old global index*/
    char temp;
    int cellcase, ret, movehappened;
#ifndef NDEBUG
    char msg[MID_MESSAGE_SIZE];
#endif

    globalIndex = mcalloc ((MOATypeInd) tbData->seqNum, (MOATypeInd) sizeof *globalIndex);
    old_globalIndex = mcalloc ((MOATypeInd) tbData->seqNum, (MOATypeInd) sizeof *old_globalIndex);
    tbData->aSeqLen[tbData->pathParts-1] = 0;
    tbData->currProc = myProcid;
    currLocalIndex = locateCell (pData, wData, tbData, &tbData->maxCellIndex, &tbData->partIndex, &tbData->currProc);
    if (currLocalIndex >= 0) {
        for (k=0;k<tbData->seqNum;k++) {
            globalIndex[k] = tbData->maxCellIndex[k];
            old_globalIndex[k] = globalIndex[k];
        }
        prevScore = currentScore = pData->msaAlgn->elements[currLocalIndex].val;
    }
    else
        goto cleanExit;

    /* Iterate making the Neighbor Cell that Computed the Current Cell, the new current cell till the origin is reached*/
#ifndef NDEBUG
    mprintf (3, "globalIndex partOrder currLocalIndex score prev_ub prev[0] prevcell prevscore\n", 1);
#endif

    while ((tbData->currProc == myProcid) && (currLocalIndex >= 0))  {
      
        movehappened = 0;
        cellcase = 2;
        if ((pData->msaAlgn->elements[currLocalIndex].prev_ub > 0) && (pData->msaAlgn->elements[currLocalIndex].prev != NULL)) {
            cellcase = 1;
            //printf ("[%d] current cell has previous chains\n", myProcid);
            for (k=0;k<tbData->seqNum;k++) 
                globalIndex[k] = pData->msaAlgn->elements[currLocalIndex].prev[0][k];
            movehappened = 1;
            currLocalIndex = locateCell (pData, wData, tbData, &globalIndex, &tbData->partIndex, &tbData->currProc);
            if (currLocalIndex >= 0) {
                currentScore = pData->msaAlgn->elements[currLocalIndex].val;
            }
        }
        else  
            currLocalIndex = checkCellNeighbors (pData, wData, tbData, &globalIndex, &currentScore);
#ifndef NDEBUG
        sprintf (msg, "[%d] -> %d case %d Score %lld->%lld Index {%lld", myProcid, tbData->currProc, cellcase, prevScore, currentScore, old_globalIndex[0]);
        for (k=1;k<tbData->seqNum;k++) 
            sprintf (msg, "%s, %lld", msg, old_globalIndex[k]);
        sprintf (msg, "%s}->{%lld", msg, globalIndex[0]);
        for (k=1;k<tbData->seqNum;k++) 
            sprintf (msg, "%s, %lld", msg, globalIndex[k]);
        sprintf (msg, "%s} locIdx %lld pIdx {%lld", msg, currLocalIndex, tbData->partIndex[0]);
        for (k=1;k<tbData->seqNum;k++) 
            sprintf (msg, "%s, %lld", msg, tbData->partIndex[k]);
        sprintf (msg, "%s} \n", msg);
        mprintf (0, msg, 1);
        printf (msg);
#endif
        movehappened = 0;
        for (k=0;k<tbData->seqNum;k++) 
            if (old_globalIndex[k] > globalIndex[k])
                movehappened = 1;
        /* only get aligned residue, in case there was a move in the current partition */
#ifndef NDEBUG
        sprintf (msg, "proc %d\n", tbData->currProc);
        mprintf (20, msg, 1);
#endif
        if (movehappened == 1) 
            tbData->aSeqLen[tbData->pathParts-1] = processTBStep (pData, tbData, &old_globalIndex, &globalIndex);
        prevScore = currentScore; 
    } /* end while loop*/
#ifndef NDEBUG
    sprintf(msg, "\n aSeqLen %ld ", tbData->aSeqLen[tbData->pathParts-1]);  
    mprintf(20, msg, 1);
#endif 
    /* Reverse the contents of the Aligned Sequences*/
    for (k=0;k<tbData->seqNum && tbData->aSeqLen[tbData->pathParts-1] > 0;k++) {
        for (j=0;j<tbData->aSeqLen[tbData->pathParts-1]/2;j++) {
            temp = tbData->algnseq[tbData->pathParts-1][k][j];
            tbData->algnseq[tbData->pathParts-1][k][j] = tbData->algnseq[tbData->pathParts-1][k][tbData->aSeqLen[tbData->pathParts-1] - j - 1];
            tbData->algnseq[tbData->pathParts-1][k][tbData->aSeqLen[tbData->pathParts-1] - j - 1] = temp;
        }
    }
    for (k=0;k<tbData->seqNum;k++) 
        tbData->maxCellIndex[k] = globalIndex[k];
#ifndef NDEBUG
    printf ("[%d] Trace Back Exiting with tbData->maxCellIndex {%lld", myProcid, tbData->maxCellIndex[0]);
    for (k=1;k<tbData->seqNum;k++)
        printf (", %lld", tbData->maxCellIndex[k]);
    printf ("}\n");
    fflush (stdout);
    sprintf (msg, "[%d] end of tb currLocalIndex %lld globalIndex {%lld", myProcid, currLocalIndex, tbData->maxCellIndex[0]);
    for (k=1;k<tbData->seqNum;k++)
        sprintf (msg, "%s, %lld", msg,  tbData->maxCellIndex[k]);
    sprintf (msg, "%s}proc %d\n", msg, tbData->currProc);
    mprintf(3, msg, 1);
    printf (msg);
#endif
cleanExit:    
    if (old_globalIndex != NULL)
        free (old_globalIndex);
    old_globalIndex = NULL;
    if (globalIndex != NULL)
        free (globalIndex);
    globalIndex = NULL;
    return tbData->aSeqLen[tbData->pathParts-1];
}

MOATypeElmVal getLocalMaxCellScore (ProcessData * pData, WavesData * wData, TracebackData * tbData, int searchAllPartitions) {
    char msg[MID_MESSAGE_SIZE];
    MOATypeInd i, j, startIndex;
    MOATypeDimn k;

	
#ifndef NDEBUG
    sprintf (msg, "prtcount %ld\n", pData->partitionsCount);
    mprintf(20, msg, 1);
#endif
    if (Mode == Distributed) {
        pData->waveNo = wData->wavesTotal - 1;
        pData->partNo = wData->partsInWave[pData->waveNo]-1;
        tbData->maxCellScore = 0;
        if (getProcID(wData, pData->waveNo, pData->partNo) != myProcid) 
            getPrevPartition (wData, &pData->waveNo, &pData->partNo);
    }
    i = 0;
    while (pData->partNo != -1) {
        if (Mode == Distributed) {
            if (restorePartitionCheckPoint(pData, wData, pData->waveNo, pData->partNo) != 0) {
                printf ("Error Retrieving partition file. Exitiing.\n");
                return -1;
            }
        }
#ifndef NDEBUG
        sprintf (msg, "i %ld elm_ub %ld\n", i, pData->msaAlgn->elements_ub);
        mprintf(20, msg, 1);
#endif
        if (i==0) {
            tbData->maxCellScore = pData->msaAlgn->elements[pData->msaAlgn->elements_ub-1].val;
            if (Mode == Distributed) {
                for (k=0;k<pData->seqNum;k++)
                    tbData->maxCellIndex[k] = pData->msaAlgn->indexes[pData->msaAlgn->elements_ub-1][k];
            }
            else 
                Gamma_Inverse(pData->msaAlgn->elements_ub-1, pData->msaAlgn->shape,  pData->msaAlgn->dimn, &tbData->maxCellIndex, 1);
            startIndex = pData->msaAlgn->elements_ub-2;
        }
        else
            startIndex = pData->msaAlgn->elements_ub-1;
        for (j=startIndex;j>=0;j--) {
#ifndef NDEBUG
            sprintf (msg, "j %ld index %ld sqnum %ld sqlen %ld %ld %ld\n", j, pData->msaAlgn->indexes[j], pData->seqNum, pData->seqLen[0], pData->seqLen[1], pData->seqLen[2]);
            mprintf(20, msg, 1);
#endif

            if (Mode != Distributed) {
                Gamma_Inverse(j, pData->msaAlgn->shape,  pData->msaAlgn->dimn, &pData->msaAlgn->indexes[j], 1);
            }
            if (isHigherBorderCellandNotLower(pData->msaAlgn->indexes[j], pData->seqNum, pData->seqLen) == 1)
                if (pData->msaAlgn->elements[j].val > tbData->maxCellScore) {
                    tbData->maxCellScore = pData->msaAlgn->elements[j].val;
                    if (Mode == Distributed) {
                        for (k=0;k<pData->seqNum;k++)
                            tbData->maxCellIndex[k] = pData->msaAlgn->indexes[j][k];
                    }
                    else 
                        Gamma_Inverse(j, pData->msaAlgn->shape,  pData->msaAlgn->dimn, &tbData->maxCellIndex, 1);
                }
#ifndef NDEBUG
            sprintf (msg, "testing index %ld score %ld\n", pData->msaAlgn->indexes[j], pData->msaAlgn->elements[j].val);
            mprintf (20, msg, 1);
#endif
        }
        if ((Mode == Distributed) && (searchAllPartitions == 1))             
            getPrevPartition (wData, &pData->waveNo, &pData->partNo);
        else
            return tbData->maxCellScore;
            
        i++;
    }
#ifndef NDEBUG
    sprintf (msg, "[%d] maxCellScore%lld, maxCellIndex = {%lld", myProcid, tbData->maxCellScore, tbData->maxCellIndex[0]);
    for (k=1;k<pData->seqNum;k++)
        sprintf (msg, "%s, %lld", msg, tbData->maxCellIndex[k]);
    sprintf (msg, "%s} sent to Master\n", msg);
    mprintf (20, msg, 1);
    printf (msg);
#endif
    return tbData->maxCellScore;
}

void sendmaxCellScore (ProcessData * pData, WavesData * wData, TracebackData * tbData) {
	char msg[MID_MESSAGE_SIZE];
	MPI_Request request;
	int MPI_return;
        
	tbData->maxCellScore = getLocalMaxCellScore(pData, wData, tbData, 0);
	MPI_return = MPI_Send (&tbData->maxCellScore, 1, MPI_LONG_LONG, 0, 0, MOAMSA_COMM_WORLD);
	MPI_return = MPI_Send (tbData->maxCellIndex, pData->seqNum, MPI_LONG_LONG, 0, 1, MOAMSA_COMM_WORLD);
}

