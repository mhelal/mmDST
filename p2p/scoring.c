/****************************************************************
* Author: Manal Helal                                           *
* Last Modification Date: Fri 12 Jan 2007 03:39:51 AM EST       *
* Project : mmDst - Distributed Multiple sequence slignment based on 	       *
* 					Mathematics of arrays - PhD Experimentation		 *
* File: scoring.c
* Description: contain the scoring functions used in 	   	    *
*       this project                                            *
* Function:
*		ComputePartitionScores
*     subScore
*     gapScore
*     initLBCell
*     getScore
*     getPrevCells
*     getNeighborScores
*     getNeghbScore
*     getRelation
*		PrintPrevChains
****************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>      /* include the character library */
#include <pthread.h>
#include <mpi.h>
#include <errno.h>
#include "globals.h"
#include "utils.h"
#include "moa.h"
#include "main.h"
#include "scores.h"
#include "scoring.h"

/* Sum of Pairs Scoring Function*/
/***********************************************************************
	Function: getCellSPScore
		Computes Sum of Pairs score for cell of index findex
***********************************************************************/
MOATypeElmVal getCellSPScore (ProcessData * pData, ScoringData * sData, WavesData * wData) {
    MOATypeDimn l, k;
 #ifndef NDEBUG
    char msg[MID_MESSAGE_SIZE];
    int dbglevel = 1;
#endif  
    sData->score = 0;
    for (l=0;l<pData->seqNum  - 1; l++) {
        for (k=l+1;k<pData->seqNum; k++) {
            sData->score += subScore(sData->sequences[l][sData->gm_index[l]-1], sData->sequences[k][sData->gm_index[k]-1], sData->stype);
#ifndef NDEBUG
            sprintf(msg, "[%d]>getCellSPScore: (%lld, %lld), (%c-%c) => score %lld\n", myProcid,l, k, sData->sequences[l][sData->gm_index[l]-1], sData->sequences[k][sData->gm_index[k]-1], sData->score);
            mprintf(dbglevel+2, msg, 1);
#endif
        }
#ifndef NDEBUG
        sprintf(msg, "[%d]>getCellSPScore:  pairwise scores %lld\n", myProcid, sData->score);
        mprintf(dbglevel+2, msg, 1);
#endif
    }    
    return sData->score;    
}
void * SPComputeScores (ProcessData * pData, ScoringData * sData, WavesData * wData) {
    MOATypeInd i; /*tensor Cells iterator*/
    MOATypeDimn k; /*Multidimensional Index Iterator*/
    char command[2000];

#ifndef NDEBUG
    char msg[MID_MESSAGE_SIZE];
    int dbglevel = 1;
#endif

    sData->msaAlgn = pData->msaAlgn;
    sData->waveNo = pData->waveNo;
    sData->partNo = pData->partNo;
    sData->p_index = pData->msaAlgn->indexes[0];
    
    /* loop the MOAAlign tensor for scores*/
    for (i = 0; i< sData->msaAlgn->elements_ub; i++)  {
        sData->findex = i;
        Gamma_Inverse(i, sData->msaAlgn->shape, sData->msaAlgn->dimn, &sData->lm_index, 1);
        //Gamma_Inverse(sData->msaAlgn->indexes[i], pData->seqLen, pData->seqNum, &sData->gm_index, 1);
        if (Mode != Distributed) {
            for (k=0;k<pData->seqNum;k++) {
                sData->gm_index[k] = sData->lm_index[k];
                sData->msaAlgn->indexes[i][k] = sData->lm_index[k];
            }
        }
        else {
            for (k=0;k<pData->seqNum;k++) 
                sData->gm_index[k] = sData->msaAlgn->indexes[i][k];                
        }
        analyzeCellPosition (sData->lm_index, i, &sData->msaAlgn, pData->seqLen);
        sData->LNCount = 0;
        sData->msaAlgn->elements[i].val = getCellSPScore(pData, sData, wData);
        currNow = getTime();
	if (currNow == NULL) {
		printf ("Could not read current time. Exiting.\n");
		return NULL;
	}
        if (isTimeDiffEquals(currNow, prevNow, 'm', 10) == 1) {
            printf("[%d]time (%d, %d, %d, %d) Scored cell %lld/%lld index {%lld", myProcid, currNow->tm_yday, currNow->tm_hour, currNow->tm_min, currNow->tm_sec, i, sData->msaAlgn->elements_ub, sData->msaAlgn->indexes[i][0]);
            for (k=1;k<pData->seqNum;k++) 
                printf (", %lld", sData->msaAlgn->indexes[i][k]);
            printf ("} score %lld part Order %ld/%ld Part %ld/%ld(T:%ld) in wave %ld/%ld\n", sData->msaAlgn->elements[i].val, pData->partNo, wData->partsInWave[pData->waveNo], pData->computedPartitions+1, pData->partitionsCount, wData->partsTotal, pData->waveNo, wData->wavesTotal);
            cpTime (currNow, &prevNow);
            sprintf (command, "prstat 1 1 > /export/home/mhelal1/thesis/exp/run/prstatus/prst_%s_%d_%d_%d_%d", outputfilename, currNow->tm_yday, currNow->tm_hour, currNow->tm_min, currNow->tm_sec);
            i = system (command);
        }
#ifndef NDEBUG
        sprintf(msg, "####[%lld] index %lld score %lld llOC %d, lhOC %d, glOC %d, ghOC %d\n", i, sData->msaAlgn->indexes[i], sData->msaAlgn->elements[i].val, sData->msaAlgn->elements[i].cp.llOC, sData->msaAlgn->elements[i].cp.lhOC, sData->msaAlgn->elements[i].cp.glOC, sData->msaAlgn->elements[i].cp.ghOC);
        mprintf(dbglevel, msg, 1);
#endif
    } /*End MOAAlign tensor's Cell Scoring Loop*/
#ifndef NDEBUG
    mprintf(dbglevel, "END LOOP ***************************\n", 1);
#endif
    pData->msaAlgn = sData->msaAlgn;
    return NULL;
}
/* Dynamic Programming Scoring Function*/
void * DPComputeScores (ProcessData * pData, ScoringData * sData, WavesData * wData) {
    MOATypeInd i; /*tensor Cells iterator*/
    MOATypeDimn k; /*Multidimensional Index Iterator*/
    char command[2000];

#ifndef NDEBUG
    char msg[MID_MESSAGE_SIZE];
    int dbglevel = 1;
#endif

    sData->msaAlgn = pData->msaAlgn;
    sData->waveNo = pData->waveNo;
    sData->partNo = pData->partNo;
    sData->p_index = pData->msaAlgn->indexes[0];
    
    /* loop the MOAAlign tensor for scores*/
    for (i = 0; i< sData->msaAlgn->elements_ub; i++)  {
        sData->findex = i;
        Gamma_Inverse(i, sData->msaAlgn->shape, sData->msaAlgn->dimn, &sData->lm_index, 1);
        //Gamma_Inverse(sData->msaAlgn->indexes[i], pData->seqLen, pData->seqNum, &sData->gm_index, 1);
        if (Mode != Distributed) {
            for (k=0;k<pData->seqNum;k++) {
                sData->gm_index[k] = sData->lm_index[k];
                sData->msaAlgn->indexes[i][k] = sData->lm_index[k];
            }
        }
        else {
            for (k=0;k<pData->seqNum;k++) 
                sData->gm_index[k] = sData->msaAlgn->indexes[i][k];                
        }
        analyzeCellPosition (sData->lm_index, i, &sData->msaAlgn, pData->seqLen);
        sData->LNCount = 0;
        sData->msaAlgn->elements[i].val = getCellDPScore(pData, sData, wData);
        currNow = getTime();
	if (currNow == NULL) {
		printf ("Could not read current time. Exiting.\n");
		return NULL;
	}
        if (isTimeDiffEquals(currNow, prevNow, 'm', 10) == 1) {
            printf("[%d]time (%d, %d, %d, %d) Scored cell %lld/%lld index {%lld", myProcid, currNow->tm_yday, currNow->tm_hour, currNow->tm_min, currNow->tm_sec, i, sData->msaAlgn->elements_ub, sData->msaAlgn->indexes[i][0]);
            for (k=1;k<pData->seqNum;k++) 
                printf (", %lld", sData->msaAlgn->indexes[i][k]);
            printf ("} score %lld part Order %ld/%ld Part %ld/%ld(T:%ld) in wave %ld/%ld\n", sData->msaAlgn->elements[i].val, pData->partNo, wData->partsInWave[pData->waveNo], pData->computedPartitions+1, pData->partitionsCount, wData->partsTotal, pData->waveNo, wData->wavesTotal);
            cpTime (currNow, &prevNow);
            sprintf (command, "prstat 1 1 > /export/home/mhelal1/thesis/exp/run/prstatus/prst_%s_%d_%d_%d_%d", outputfilename, currNow->tm_yday, currNow->tm_hour, currNow->tm_min, currNow->tm_sec);
            i = system (command);
        }
#ifndef NDEBUG
        sprintf(msg, "####[%lld] index %lld score %lld llOC %d, lhOC %d, glOC %d, ghOC %d\n", i, sData->msaAlgn->indexes[i], sData->msaAlgn->elements[i].val, sData->msaAlgn->elements[i].cp.llOC, sData->msaAlgn->elements[i].cp.lhOC, sData->msaAlgn->elements[i].cp.glOC, sData->msaAlgn->elements[i].cp.ghOC);
        mprintf(dbglevel, msg, 1);
#endif
    } /*End MOAAlign tensor's Cell Scoring Loop*/
#ifndef NDEBUG
    mprintf(dbglevel, "END LOOP ***************************\n", 1);
#endif
    pData->msaAlgn = sData->msaAlgn;
    return NULL;
}


/*******************************************************
	Function: initLBCell
*******************************************************/
MOATypeElmVal initLBCell (MOATypeShape * lm_index, MOATypeDimn dimn, int stype) {
	MOATypeElmVal iscore, indexTau;
	MOATypeDimn i;
  	if (initializationMethod == IndexProductInit)
		indexTau = lm_index[0]+1;
	else if (initializationMethod == IndexSumInit)
		indexTau = lm_index[0];
	else {
		return gapScore(stype);
		
	}
	for (i = 1; i < dimn;i++) {
	  	if (initializationMethod == IndexProductInit)
			indexTau = indexTau * (lm_index[i] + 1);
		else if (initializationMethod == IndexSumInit)
			indexTau = indexTau + lm_index[i];
	}      
	iscore = (gapScore(stype) * indexTau);
	return iscore;
}

/***********************************************************************
	Function: getCellDPScore
		Computes Dynamic Programming score for cell of index findex, 
	1. by initializing : If a Global lower Overlapping Cell in distributed mode, or lower border cell generally in sequential mode
	2. or reading from a local other partition, if a local lower Overlapping Cell 
	3. or receiving from a remote processor, if a rempte local lower Overlapping Cell 
	4. or retrieve its lower neighbors scores, if in sequential mode, or an inner local cell. 
***********************************************************************/
MOATypeElmVal getCellDPScore (ProcessData * pData, ScoringData * sData, WavesData * wData) {
    MOATypeInd startOC;
    int ret, inSearchSpace;
    MOATypeDimn i;
#ifndef NDEBUG
    MOATypeDimn k;
    char msg[MID_MESSAGE_SIZE];
    int dbglevel = 3;
#endif

	
#ifndef NDEBUG
    sprintf (msg, "[%d] gloc %d lloc %d lhoc %d Cell: {%lld", myProcid, sData->msaAlgn->elements[sData->findex].cp.glOC, sData->msaAlgn->elements[sData->findex].cp.llOC, sData->msaAlgn->elements[sData->findex].cp.lhOC, sData->gm_index[0]);
    mprintf (dbglevel, msg, 1); 
    for (i = 1; i < sData->seqNum; i++)  {
    	sprintf (msg, ", %lld", sData->gm_index[i]);
        mprintf (dbglevel, msg, 1); 
    }
    mprintf (dbglevel, "} ", 1);
#endif
    if (Mode == Distributed) {    
        if (sData->msaAlgn->elements[sData->findex].cp.glOC == 1) {
            /*if (GLB) Global lower border cell from the whole tensor, then initialize with gabscores multiplied by indices if global, or zero if local alignment*/
            sData->score = initLBCell (sData->gm_index, pData->seqNum,pData->stype);
#ifndef NDEBUG
            sprintf (msg, ">Initialized: 1. score %lld\n", sData->score);
            mprintf (dbglevel, msg, 1);
#endif
        } else if (sData->msaAlgn->elements[sData->findex].cp.llOC == 1) { 
            /*search for the local lower Overlapping Cell in previously calculated partitions in local processor first*/
            /* if not found, then block till it is received from the adjacent processor  */
            ret = getCellScore (pData, sData, wData, sData->gm_index, &sData->score, &inSearchSpace, 1, -1);
            
            while ((inSearchSpace == 0) && (ret < 0))
                ret = receiveOC(pData, sData); 
        } else {
            /*Inner Cell not gloc or lloc, score*/
            ret = getNeighborScores (pData, sData, wData);
            if (ret == -2) {
                /* Lower Border Cell - Initialize ===========================*/
                sData->score = initLBCell (sData->gm_index, pData->seqNum,sData->stype);
#ifndef NDEBUG
                sprintf (msg, ">initLBCell: 6. score %lld\n", sData->score);
                mprintf(dbglevel, msg, 1);
#endif
            } 
            else {
#ifndef NDEBUG
                sprintf (msg, ">getNeighborScores: 5. score %lld\n", sData->score);
                mprintf(dbglevel, msg, 1);
#endif
            }
        }
    } else {
	/* If Sequential Scoring  Mode NOT Distributed ============================*/
        if ((sData->msaAlgn->elements[sData->findex].cp.glOC == 1) || (sData->msaAlgn->elements[sData->findex].cp.llOC == 1)) {
            /* Lower Border Cell - Initialize ===========================*/
            sData->score = initLBCell (sData->gm_index, pData->seqNum, sData->stype);
#ifndef NDEBUG
            sprintf (msg, "[%d]>getScore: Seq Lower Bound => score %lld\n", myProcid, sData->score);
            mprintf(dbglevel, msg, 1);
#endif
        }
        else {
            ret = getNeighborScores (pData, sData, wData);
            if (ret == -2) {
                /* Can not get neighbors, Initialize, and check design :) ==*/
                sData->score = initLBCell (sData->gm_index, pData->seqNum,sData->stype);
#ifndef NDEBUG
                sprintf (msg, "[%d]>getScore: Seq Inner Cell Initialized, error to check=> score %lld\n", myProcid, sData->score);
                mprintf(dbglevel, msg, 1);
#endif
            } 
        }
    }
    
    if (AlignmentType == Local) {
        if (sData->score < 0)
            sData->score = 0;
    }
#ifndef NDEBUG
    sprintf (msg, "[%d]>getScore: element %lld score %lld\n", myProcid, sData->msaAlgn->indexes[sData->findex], sData->msaAlgn->elements[sData->findex].val);
    mprintf (dbglevel, msg, 1);
    sprintf (msg, "will check lhOC\n");
    mprintf(dbglevel, msg, 1);
#endif
    if ((Mode == Distributed) && (sData->msaAlgn->elements[sData->findex].cp.lhOC == 1) && (sData->msaAlgn->elements[sData->findex].cp.glOC == 0)) {
#ifndef NDEBUG
        mprintf (dbglevel, "Will Add to OCout Buffer", 1);
#endif		
        addOC (pData, sData, wData);
#ifndef NDEBUG
        sprintf (msg, "after addOC wavesOC %ld\n", pData->OCout[pData->waveNo].wavesOC);
        mprintf(dbglevel, msg, 1);
#endif
    }
    return sData->score;
}


/*****************************************************************************
	Function: getRelation 
		Returns the neighbors that got decremented , and those that didn't
*****************************************************************************/


long getRelation (ScoringData * sData, MOATypeInd neighbIndex) {
    long i, cnt, ndecr = 0;
#ifndef NDEBUG
    char msg[SHORT_MESSAGE_SIZE];
    int dbglevel = 6;
#endif	
    cnt = 0;
    ndecr = 0;
    for (i = 0; i< sData->msaAlgn->dimn; i++) {
        if (sData->NghbMOA->indexes[neighbIndex][i] < sData->gm_index[i]) {
            sData->decremented[cnt] = i;
#ifndef NDEBUG
            sprintf (msg, "[%d]>getRelation: decr[%ld] = %lld\n", myProcid, cnt, sData->decremented[cnt]);
            mprintf(dbglevel, msg, 1);
#endif
            cnt ++;
        }
        else 
            ndecr ++;
#ifndef NDEBUG
        sprintf (msg, "[%d]>getRelation: ngb %lld cell %lld cnt %ld ndcr %ld\n", myProcid, sData->NghbMOA->indexes[neighbIndex][i], sData->gm_index[i], cnt, ndecr);
        mprintf(dbglevel, msg, 1);
#endif
    } 
    return ndecr;
}
/********************************************************************************
	Function: getNeghbScore 
		Examine a specific neighbor's score based on its distance from the current 
		cell being scored
********************************************************************************/
MOATypeElmVal getNeghbScore (ScoringData * sData, MOATypeInd neighbIndex) {
    MOATypeElmVal totpwiseScore, score;
    MOATypeDimn j, k, ndecr = 0;
#ifndef NDEBUG
    char msg[MID_MESSAGE_SIZE];
    int dbglevel = 5;
#endif	
    /* get relation between this neighbor and the current cell */
    //ndecr = getRelation (sData, neighbIndex);
    /* if we have 2 or more decremented neighbors*/
    /* sum pairwise scores for each pair of dimensions that got decremented*/
    totpwiseScore = 0;
    ndecr = 0;

/*old code before merging the getRelation loop*/
//    for (j=0;j<(sData->msaAlgn->dimn  - ndecr) - 1; j++) {
//        for (k=j+1;k<(sData->msaAlgn->dimn  - ndecr); k++) {
//                totpwiseScore += sData->pwScores[sData->decremented[j]][sData->decremented[k]];
//                sprintf(msg, "[%d]>getNeghbScore: pw[%lld][%lld] = %lld (%c-%c)\n", myProcid, sData->decremented[j], sData->decremented[k], sData->pwScores[sData->decremented[j]][sData->decremented[k]], sData->sequences[j][sData->gm_index[j]-1], sData->sequences[k][sData->gm_index[k]-1]);
    
    // to merge the getRelation loop
    int gapPenalty = 0;
    for (j=0;j<sData->msaAlgn->dimn-1; j++) {
        if (sData->NghbMOA->indexes[neighbIndex][j] < sData->gm_index[j]) {
            for (k=j+1;k<sData->msaAlgn->dimn; k++) {            
                if (sData->NghbMOA->indexes[neighbIndex][k] < sData->gm_index[k]) {
                    totpwiseScore += sData->pwScores[j][k];
                    
#ifndef NDEBUG
                    sprintf(msg, "[%d]>getNeghbScore: pw[%lld][%lld] = %lld (%c-%c)\n", myProcid, j, k, sData->pwScores[j][k], sData->sequences[j][sData->gm_index[j]-1], sData->sequences[k][sData->gm_index[k]-1]);
                    mprintf(dbglevel, msg, 1);
#endif
                }
            }
        }
        else {
            ndecr++;        
            if ((sData->NghbMOA->elements[neighbIndex].prev != NULL) && (sData->NghbMOA->elements[neighbIndex].prev_ub > 0)) { 
                if (sData->NghbMOA->elements[neighbIndex].prev[0][j] < sData->NghbMOA->indexes[neighbIndex][j]) 
                    gapPenalty += gapExtension;
                else 
                    gapPenalty += gapOpenning;
            }
            else 
                gapPenalty += gapOpenning;
        }
    }
    j = sData->msaAlgn->dimn-1;
    if (sData->NghbMOA->indexes[neighbIndex][j] == sData->gm_index[j]) {
        ndecr++;        
        if ((sData->NghbMOA->elements[neighbIndex].prev != NULL) && (sData->NghbMOA->elements[neighbIndex].prev_ub > 0)) { 
            if (sData->NghbMOA->elements[neighbIndex].prev[0][j] < sData->NghbMOA->indexes[neighbIndex][j]) 
                gapPenalty += gapExtension;
            else 
                gapPenalty += gapOpenning;
        }
        else 
            gapPenalty += gapOpenning;
    }
    
    /* multiply the number of dimensions that did not get decremented by the gap score */
    /* add both to the previoues score in the neighbor's cell in the alignment matrix to get the temporary score for this neighbor*/

    /*inserting one residue, is like deleting one residue, so ndecr of 2 is like ndec of 1 in a 3D, and ndecr of 3 is like ndecr of 1 in 4D and */
    div_t div_val;
    div_val = div(sData->msaAlgn->dimn, 2);
    
    if (ndecr<=div_val.quot)
        score = totpwiseScore + gapPenalty;
    else {
        /*because we already added up the gapPenalty for all decremented indices based on whether it should be gap openning or gap extension
         Now, if the number of decremented indices is more than half the dimension, we will divide what we added up and multiply this average
         on the remaining decremented indices*/
        div_val = div(gapPenalty, ndecr);
        score = totpwiseScore + ((sData->msaAlgn->dimn-ndecr) * div_val.quot);
    }
  
    return score;
}

MOATypeElmVal averageNeighborsScore (ProcessData * pData, ScoringData * sData, WavesData * wData, MOATypeShape * cellIndex) {
    /*Otherwise, average the score value based on neighbor's lower neighbors that are inside the search space*/
    MOATypeElmVal averageScore = 0;
    MOATypeInd j, n, NeighbFlatIndex;
    MOA_rec * NnghbMOA;
    int ret, inSearchSpace;
    if (getNeighbors (2, cellIndex, sData->msaAlgn->dimn, sData->seqLen, &NnghbMOA) == 0) {
        MOATypeElmVal * NneighborScores = mmalloc ((NnghbMOA->elements_ub - 1) * sizeof *NneighborScores);
        n = 0; /*Number of Neighbors scored*/
        for (j=0;j<NnghbMOA->elements_ub - 1;j++) { 
            ret = getCellScore (pData, sData, wData, NnghbMOA->indexes[j], &NneighborScores[n], &inSearchSpace, 1, -1);
            /*if (ret < 0)
                NneighborScores[j] =  a_average (NneighborScores, j-1); /*Last resort - average the values read so far*/
            if (ret == 0) 
                n ++;            
        } /*End loop for neighbor's neighbors*/
        averageScore = a_average (NneighborScores, n);        
        free (NneighborScores);
        deleteMOA(NnghbMOA);
    }
    return averageScore;
}

int getCellScore (ProcessData * pData, ScoringData * sData, WavesData * wData, MOATypeShape * cellIndex, MOATypeElmVal * score, int * inSearchSpace, int NeighborSearch, MOATypeInd NeighbIndex) {
    int ret = 0;
    MOATypeDimn k;
    MOATypeInd NeighbFlatIndex;
    /*Check if cellIndex is found in the current scoring partition*/
    if ((NeighborSearch == 1) && (IsCellInPart(cellIndex, sData->p_index, sData->seqNum, sData->seqLen, pData->partitionSize) == 0) && 
        (getLocalIndex (cellIndex, sData->p_index, sData->seqNum, sData->seqLen, pData->partitionSize, &sData->neighbor) == 0)) {
        NeighbFlatIndex = Gamma(sData->neighbor, sData->msaAlgn->dimn, sData->msaAlgn->shape,  sData->msaAlgn->dimn, 1);
       (*score) = sData->msaAlgn->elements[NeighbFlatIndex].val;
       if (sData->msaAlgn->elements[NeighbFlatIndex].prev != NULL && sData->msaAlgn->elements[NeighbFlatIndex].prev_ub > 0 && 
           sData->NghbMOA != NULL && NeighbIndex >= 0 && NeighbIndex < sData->NghbMOA->elements_ub) {
           sData->NghbMOA->elements[NeighbIndex].prev = mmalloc(sizeof *sData->NghbMOA->elements[NeighbIndex].prev);
           sData->NghbMOA->elements[NeighbIndex].prev_ub = 1;
           sData->NghbMOA->elements[NeighbIndex].prev[0] = mmalloc(sData->seqNum * sizeof *sData->NghbMOA->elements[NeighbIndex].prev[0]);
           for (k=0;k<sData->seqNum;k++)
                sData->NghbMOA->elements[NeighbIndex].prev[0][k] = sData->msaAlgn->elements[NeighbFlatIndex].prev[0][k];
       }
    }
    else {
        /*check if neighbor's partition is included in search space*/
        MOATypeShape * partIndex = mmalloc (pData->seqNum * sizeof *partIndex);
        if  (getPartitionIndex (cellIndex, pData->seqNum, pData->seqLen, wData->partitionSize, &partIndex) == 0) {
            if ((*inSearchSpace) = isPartInSearchSpace(partIndex, wData) == 0) {
                long waveNo, partNo;
                getPartitionPosition (wData, partIndex, &waveNo, &partNo);
                if (partNo >= 0) {
                    if (myProcid == getProcID (wData, waveNo, partNo)) {                        
                        /*Check if Neighbor is found in other local partitions OCout Buffer*/
                        if(checkPrevPartitions(pData, cellIndex, score) != 0) {
                            /*average the neighboring (up to 2 strides) cell scores*/
                            (*score) = averageNeighborsScore(pData, sData, wData, cellIndex);
                        }
                    }
                    /*Check if Neighbor is already received from other processors in OCin Buffer*/
                    else if (checkRecvOC(pData, wData, cellIndex, score, 0) != 0)    
                        ret = -1;
                }
            }
        }
        free (partIndex);
    }
    return ret;
}

/****************************************************************************
	Function: getNeighborScores - based on Global partition neighbors (the right one I hope)
		Examine all neighbors scores to determine the current score of the new 
		cell being scored.
****************************************************************************/
int getNeighborScores (ProcessData * pData, ScoringData * sData, WavesData * wData) {
    MOATypeInd i, locCellIndex, partIndex, NeighbFlatIndex;
    MOATypeDimn k, l;
    MOATypeElmVal maxScore = 0, nghbScore = 0;
    char msg[SHORT_MESSAGE_SIZE];
    int inSearchSpace;
#ifndef NDEBUG
    int dbglevel = 4;

    sprintf (msg, "[%d]>getNeighborScores: cell {%lld", myProcid, sData->gm_index[0]);
    mprintf (dbglevel, msg, 1);
    for (k=1;k<sData->seqNum;k++) {
        sprintf (msg, ", %lld", sData->gm_index[k]);
        mprintf (dbglevel, msg, 1);
    }
    mprintf (dbglevel, "}\n", 1);
#endif
    /*Fill in the pair wise score matrix for the current cell*/
    for (l=0;l<pData->seqNum  - 1; l++) {
        for (k=l+1;k<pData->seqNum; k++) {
            sData->pwScores[l][k] = subScore(sData->sequences[l][sData->gm_index[l]-1], sData->sequences[k][sData->gm_index[k]-1], sData->stype);
#ifndef NDEBUG
            sprintf(msg, "[%d]>getNeighborScores: (%lld, %lld), (%c-%c) => score %lld\n", myProcid,l, k, sData->sequences[l][sData->gm_index[l]-1], sData->sequences[k][sData->gm_index[k]-1], sData->pwScores[l][k]);
            mprintf(dbglevel+2, msg, 1);
#endif
        }
    }
  
    /*get neighbors of the current cell */
    if (getLowerNeighbors (2, sData->gm_index, sData->msaAlgn->dimn, sData->seqLen, &sData->NghbMOA) != 0) {
        sprintf(msg, "[%d]>getNeighborScores: No Lower neighbors while expecting, returning\n", myProcid);
        mprintf(0, msg, 1);
        fflush(stdout);
        return -1; /* Error - need debugging*/
    }
    /* loop throught neighbors */        
    for (i=0;i<sData->NghbMOA->elements_ub - 1;i++) { 
        /*Check if Sequential Processing, read neighbor score locally*/
        if (Mode != Distributed)  {
            MOATypeInd NeighbFlatIndex = Gamma(sData->NghbMOA->indexes[i], sData->msaAlgn->dimn, sData->msaAlgn->shape,  sData->msaAlgn->dimn, 1);
            nghbScore = sData->NghbMOA->elements[i].val = sData->msaAlgn->elements[NeighbFlatIndex].val;
            if (sData->msaAlgn->elements[NeighbFlatIndex].prev != NULL && sData->msaAlgn->elements[NeighbFlatIndex].prev_ub > 0) {
                sData->NghbMOA->elements[i].prev = mmalloc(sizeof *sData->NghbMOA->elements[i].prev);
                sData->NghbMOA->elements[i].prev_ub = 1;
                sData->NghbMOA->elements[i].prev[0] = mmalloc(sData->seqNum * sizeof *sData->NghbMOA->elements[i].prev[0]);
                for (k=0;k<sData->seqNum;k++)
                        sData->NghbMOA->elements[i].prev[0][k] = sData->msaAlgn->elements[NeighbFlatIndex].prev[0][k];
            }
        }
        else 
            getCellScore (pData, sData, wData, sData->NghbMOA->indexes[i], &nghbScore, &inSearchSpace, 1, i);
        nghbScore += getNeghbScore (sData, i);
        if (i == 0) {
            maxScore = nghbScore;            
            sData->msaAlgn->elements[sData->findex].prev = mmalloc ((MOATypeInd) sizeof *(sData->msaAlgn->elements[sData->findex].prev));
            if (sData->msaAlgn->elements[sData->findex].prev == NULL) {
                printf ("[%d]Error creating memory for cells Previous cells chains.\n", myProcid);
                return -1;
            }
            sData->msaAlgn->elements[sData->findex].prev[0] = NULL;
            sData->msaAlgn->elements[sData->findex].prev[0] = mcalloc ((MOATypeInd) sData->seqNum, (MOATypeInd) sizeof *(sData->msaAlgn->elements[sData->findex].prev[0]));
            if (sData->msaAlgn->elements[sData->findex].prev[0] == NULL) {
                printf ("[%d]Error creating memory for cells Previous cells chains.\n", myProcid);
                return -1;
            }

            for (k=0;k<sData->seqNum;k++) 
                sData->msaAlgn->elements[sData->findex].prev[0][k] = sData->NghbMOA->indexes[i][k];

            sData->msaAlgn->elements[sData->findex].prev_ub = 1;
        }
        else if (nghbScore > maxScore) {
            if (sData->msaAlgn->elements[sData->findex].prev != NULL) {
                for (k=0;k<sData->msaAlgn->elements[sData->findex].prev_ub;k++) {
                    if (sData->msaAlgn->elements[sData->findex].prev[k] != NULL)
                        free (sData->msaAlgn->elements[sData->findex].prev[k]);
                    sData->msaAlgn->elements[sData->findex].prev[k] = NULL;
                }
                free (sData->msaAlgn->elements[sData->findex].prev);
            }
            sData->msaAlgn->elements[sData->findex].prev = NULL;            
            sData->msaAlgn->elements[sData->findex].prev  = mmalloc ((MOATypeInd) sizeof *(sData->msaAlgn->elements[sData->findex].prev));
            if (sData->msaAlgn->elements[sData->findex].prev == NULL) {
                printf ("[%d]Error creating memory for cells Previous cells chains.\n", myProcid);
                return -1;
            }
            maxScore = nghbScore;
            sData->msaAlgn->elements[sData->findex].prev[0] = NULL;
            sData->msaAlgn->elements[sData->findex].prev[0] = mcalloc ((MOATypeInd) sData->seqNum, (MOATypeInd) sizeof *(sData->msaAlgn->elements[sData->findex].prev[0]));
            if (sData->msaAlgn->elements[sData->findex].prev[0] == NULL) {
                printf ("[%d]Error creating memory for cells Previous cells chains.\n", myProcid);
                return -1;
            }

            for (k=0;k<sData->seqNum;k++) 
                sData->msaAlgn->elements[sData->findex].prev[0][k] = sData->NghbMOA->indexes[i][k];        
            sData->msaAlgn->elements[sData->findex].prev_ub = 1;
        }
        else if (nghbScore == maxScore) {
            sData->msaAlgn->elements[sData->findex].prev  = realloc (sData->msaAlgn->elements[sData->findex].prev, ((MOATypeInd) (sData->msaAlgn->elements[sData->findex].prev_ub + 1)) * ((MOATypeInd) sizeof *(sData->msaAlgn->elements[sData->findex].prev)));
            if (sData->msaAlgn->elements[sData->findex].prev == NULL) {
                printf ("[%d]Error reallocating memory for Previous cells chains.\n", myProcid);
                return -1;
            }
            sData->msaAlgn->elements[sData->findex].prev[sData->msaAlgn->elements[sData->findex].prev_ub] = NULL;
            sData->msaAlgn->elements[sData->findex].prev[sData->msaAlgn->elements[sData->findex].prev_ub] = mcalloc ((MOATypeInd) sData->seqNum, ((MOATypeInd) sizeof *(sData->msaAlgn->elements[sData->findex].prev[sData->msaAlgn->elements[sData->findex].prev_ub])));
            if (sData->msaAlgn->elements[sData->findex].prev[sData->msaAlgn->elements[sData->findex].prev_ub] == NULL) {
                printf ("[%d]Error reallocating memory for Previous cells chains.\n", myProcid);
                return -1;
            }	

            for (k=0;k<sData->seqNum;k++)
                sData->msaAlgn->elements[sData->findex].prev[sData->msaAlgn->elements[sData->findex].prev_ub][k] = sData->NghbMOA->indexes[i][k];

            sData->msaAlgn->elements[sData->findex].prev_ub ++;
        }
    }

    sData->score = maxScore;
    return 0;  
}
int PrintPrevChains (MOA_rec *  msaAlgn) {
    MOATypeInd i;
    MOATypeDimn j, k;

    printf ("Index\t\t\tprev_ub\t\t\tprev[0]\n");
    for (i=0;i<msaAlgn->elements_ub;i++) {
        printf ("{%lld", msaAlgn->indexes[i][0]);
        for (k=1;k<msaAlgn->dimn;k++)
            printf (", %lld", msaAlgn->indexes[i][k]);
        printf ("}\t\t\t%lld\t\t\t{", msaAlgn->elements[i].prev_ub);
        for (j=0;j<msaAlgn->elements[i].prev_ub;j++) {
            printf ("{%lld", msaAlgn->elements[i].prev[0][0]);
            for (k=1;k<msaAlgn->dimn;k++)
                printf (", %lld", msaAlgn->elements[i].prev[0][k]);
            printf ("} ");
        }
        printf ("}\n");
    }

    return 0;
}

