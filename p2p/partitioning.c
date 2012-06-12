
/***************************************************************
* Author: Manal Helal														                             *
* Last Modification Date: Fri 12 Jan 2007 03:39:51 AM EST                            *
* Project : MMSA - Multiple Sequence Alignment Based on 	                         *
* 					Mathematics of Arrays - PhD Experimentation		                     *
* File: partitioning.c, contain the MOA socring tensor	                                     * 
* partitioning & scheduling functions                                                             *
* Function:
*		initProcessMemory
*     int checkPoint (ProcessData * pData, long startPart);
*     int checkPart (ProcessData * pData, long PartNo);
*     int restoreCheckPoint (ProcessData * * pData, WavesData * * wData);
*		calcWaves
*     getPartitionsNum
*     getWavePartsTotalp
*     getPIndicesinWave
*     getNextPIndex
*     addPartitionIndex
*     notPreviouslyVisited
*     getProcID
*     freeProcessMemory
*     getPartition
*     getNextPartition
*     IsCellInPart
***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <mpi.h>
#include "moaDst.h"
#include "globals.h"
#include "utils.h"
#include "moa.h"
#include "main.h"

/***************************************************************************
	Function: initProcessMemory
	Description:
		Initialize Process Data (pData) with initial values.
		
***************************************************************************/
int initProcessMemory (ProcessData * * pData, ScoringData * * sData, WavesData * * wData, MOATypeDimn seqNum, MOATypeShape * seqLen, char * * sequences, char * * seqName, int stype, long partitionSize) {

    MOATypeInd i;

    (*pData) = NULL;
    (*pData) =  mmalloc ((MOATypeInd) sizeof *(*pData));
    if ((*pData) == NULL)
            return -1;

    (*sData) = NULL;
    (*sData) =  mmalloc ((MOATypeInd) sizeof *(*sData));
    if ((*sData) == NULL)
            return -1;

    (*wData) = NULL;
    (*wData) =  mmalloc ((MOATypeInd) sizeof *(*wData));
    if ((*wData) == NULL)
            return -1;
    (*pData)->stype = stype;
    (*sData)->stype = stype;
    (*pData)->seqNum = seqNum;
    (*sData)->seqNum = seqNum;
    (*pData)->seqLen = seqLen;
    (*wData)->seqLen = seqLen;
    (*sData)->seqLen = seqLen;
    (*pData)->mpi_requests = NULL;
    createMOAStruct(&(*pData)->msaAlgn);
    (*pData)->msaAlgn->dimn = seqNum;
    (*pData)->msaAlgn->shape = mcalloc ((MOATypeInd) seqNum, (MOATypeInd) sizeof *(*pData)->msaAlgn->shape);
    (*pData)->commBufSize = seqNum+2;
    (*pData)->buffer_size = NULL;
    (*pData)->buffer = NULL;
    (*pData)->proc_index = NULL;
    (*pData)->depProcCount = 0;
    (*pData)->globalWaveNo = 0;
    (*pData)->seqName = seqName;
    (*pData)->sequences = sequences;
    
    (*sData)->sequences = sequences;
    (*sData)->msaAlgn = (*pData)->msaAlgn;
    /*Multidimensional Indices required during scoring*/
    createMOAStruct (&(*sData)->NghbMOA);
    (*sData)->NghbMOA->dimn = seqNum;
    (*sData)->NghbMOA->shape = mcalloc ((MOATypeInd) seqNum, (MOATypeInd) sizeof *(*pData)->msaAlgn->shape);
    
    /*Multidimensional Indices required during scoring*/
    (*sData)->decremented = NULL; 
    (*sData)->decremented = mcalloc ((MOATypeInd) seqNum, (MOATypeInd) sizeof *(*sData)->decremented);
    if ((*sData)->decremented == NULL) {
        printf ("Error Allocating Memory for Scoring Decremented Indices.\n");
        return -1; 
    }

    (*sData)->p_index = NULL; 
    (*sData)->p_index = mcalloc ((MOATypeInd) seqNum, (MOATypeInd) sizeof *(*sData)->p_index);
    if ((*sData)->p_index == NULL) {
        printf ("Error Allocating Memory for Scoring part Indices.\n");
        return -1; 
    }

    (*sData)->lm_index = NULL; 
    (*sData)->lm_index = mcalloc ((MOATypeInd) seqNum, (MOATypeInd) sizeof *(*sData)->lm_index);
    if ((*sData)->lm_index == NULL) {
        printf ("Error Allocating Memory for Scoring Local Indices.\n");
        return -1; 
    }


    (*sData)->gm_index = NULL; 
    (*sData)->gm_index = mcalloc ((MOATypeInd) seqNum, (MOATypeInd) sizeof *(*sData)->gm_index);
    if ((*sData)->gm_index == NULL) {
        printf ("Error Allocating Memory for Scoring Global Indices.\n");
        return -1; 
    }

    (*sData)->depPart = NULL; 
    (*sData)->depPart = mcalloc ((MOATypeInd) seqNum, (MOATypeInd) sizeof *(*sData)->depPart);
    if ((*sData)->depPart == NULL) {
        printf ("Error Allocating Memory for Scoring dependency part Indices.\n");
        return -1; 
    }

    (*sData)->posDimn = NULL; 
    (*sData)->posDimn = mcalloc ((MOATypeInd) seqNum, (MOATypeInd) sizeof *(*sData)->posDimn);
    if ((*sData)->posDimn == NULL) {
        printf ("Error Allocating Memory for Scoring dependency neighboring position Indices.\n");
        return -1; 
    }

    (*sData)->neighbor = NULL; 
    (*sData)->neighbor = mcalloc ((MOATypeInd) seqNum, (MOATypeInd) sizeof *(*sData)->neighbor);
    if ((*sData)->neighbor == NULL) {
        printf ("Error Allocating Memory for Scoring dependency neighboring position Indices.\n");
        return -1; 
    }

    /* Dependent Neighboring Indices and Scores*/
    (*sData)->CalLnCount = (MOATypeDimn) mpow(2, seqNum) - 1;
    if ((*sData)->CalLnCount < seqNum) {
        printf ("Error calculating Number of neighbors [%lld] for %lld sequences. Exiting.\n", (*sData)->CalLnCount, seqNum);
        return -1; 
    }

    (*pData)->validPartitions = NULL;
    (*pData)->validPartitions = mcalloc ((MOATypeInd) (*sData)->CalLnCount, ((MOATypeInd) sizeof *(*pData)->validPartitions));
    if ((*pData)->validPartitions == NULL) {
        printf ("Error Allocating memory for validPartitions. Exiting\n");
        return -1;
    }
    
    /*(*sData)->lnScores = NULL; 
    (*sData)->lnScores = mcalloc ((MOATypeInd) (*sData)->CalLnCount, (MOATypeInd) sizeof *(*sData)->lnScores);
    if ((*sData)->lnScores == NULL) {
        printf ("Error Allocating Memory for Neighbor's Scores.\n");
        return -1; 
    }

    (*sData)->lnIndices = NULL;
    (*sData)->lnIndices = mcalloc ((MOATypeInd) (*sData)->CalLnCount, (MOATypeInd) sizeof *(*sData)->lnIndices);
    if ((*sData)->lnIndices == NULL) {
        printf ("Error Allocating Memory for Neighbor's Indices.\n");
        return -1; 
    }
    for (i=0; i<(*sData)->CalLnCount; i++)  {
        (*sData)->lnIndices[i] = NULL;
        (*sData)->lnIndices[i] = mcalloc ((MOATypeInd) seqNum, (MOATypeInd) sizeof *(*sData)->lnIndices[i]);
        if ((*sData)->lnIndices[i] == NULL) {
            printf ("Error Allocating Memory for Neighbor's Indices Details.\n");
            return -1; 
        }
    }
    
    (*sData)->lnInSearchSpace = NULL; 
    (*sData)->lnInSearchSpace = mcalloc ((MOATypeInd) (*sData)->CalLnCount, (MOATypeInd) sizeof *(*sData)->lnInSearchSpace);
    if ((*sData)->lnInSearchSpace == NULL) {
        printf ("Error Allocating Memory for Neighbor's in Search Space.\n");
        return -1; 
    }*/

    /* Create Memory for the pairwise scores for each current index*/
    (*sData)->pwScores = NULL; 
    (*sData)->pwScores = mcalloc ((MOATypeInd) seqNum, (MOATypeInd) sizeof *(*sData)->pwScores); 
    if ((*sData)->pwScores == NULL) {
        printf ("Error Allocating Memory for Pairwise Scores.\n");
        return -1; 
    }
    for (i=0; i<seqNum; i++)  {
        (*sData)->pwScores[i] = NULL; 
        (*sData)->pwScores[i] = mcalloc ((MOATypeInd) seqNum, ((MOATypeInd) sizeof *((*sData)->pwScores[i])));
        if ((*sData)->pwScores[i] == NULL) {
            printf ("Error Allocating Memory for Pairwise 	Scores[%lld].\n", i);
            return -1; 
        }
    }


    (*pData)->OCin = NULL;
    (*pData)->OCout = NULL;
    (*pData)->waveNo = (*pData)->partNo = (*pData)->partitionsCount = (*pData)->computedPartitions = (*pData)->sendOC = (*pData)->sendOCPart = (*pData)->compFinished = 0;

    (*pData)->partitionSize = partitionSize;
    (*wData)->partitionSize = partitionSize;
    (*wData)->waveMiddle = NULL;
    (*wData)->waveLength = NULL;
    (*wData)->seqNum = seqNum;
    (*wData)->partsInWave = NULL;	
    (*wData)->partsInWaveIndices = NULL;
    (*wData)->wavesTotal = (*wData)->partsTotal = 0;
    if (Mode != Distributed) {
        (*wData)->partsInWave = mcalloc ((MOATypeInd) 1, (MOATypeInd) sizeof *((*wData)->partsInWave));	
        if ((*wData)->partsInWave == NULL) {
            printf ("Error Allocating Memory for partsInWave. Exiting.\n");
            return -1; 
        }
        (*wData)->partsInWaveIndices = mcalloc ((MOATypeInd) 1, (MOATypeInd) sizeof *((*wData)->partsInWaveIndices));	 
        if ((*wData)->partsInWaveIndices == NULL) {
            printf ("Error Allocating Memory for partsInWaveIndices. Exiting.\n");
            return -1; 
        }
        (*wData)->partsInWaveIndices[0] = NULL;
        (*wData)->partsInWaveIndices[0] = mcalloc ((MOATypeInd) 1, (MOATypeInd) sizeof *((*wData)->partsInWaveIndices[0]));	 
        if ((*wData)->partsInWaveIndices[0] == NULL) {
            printf ("Error Allocating Memory for partsInWaveIndices[0]. Exiting.\n");
            return -1; 
        }
        (*wData)->wavesTotal = (*wData)->partsTotal = (*wData)->partsInWave[0] = 1;
        (*wData)->partsInWaveIndices[0][0] = mcalloc ((MOATypeInd) seqNum, (MOATypeInd) sizeof *((*wData)->partsInWaveIndices[0][0]));	 
        for (i=0; i<seqNum; i++)  
            (*wData)->partsInWaveIndices[0][0][i] = 0;
    }
    
    return 0;
}
/***************************************************************************
	Function: checkPointWavesCalculations
	Description:
		write waves Calculations done by process of rank 0 to be read by other peers, instead of redoing the calculationg in all processes.
***************************************************************************/

int checkPointWavesCalculations (ProcessData * pData, WavesData * wData) {
    FILE * sfile;
    int ret;
    long i, j;
    MOATypeDimn k;
    char sfilename[FILENAME_MAX_LENGTH], bkfilename[FILENAME_MAX_LENGTH];
	
    /* file chpM... stores waves info ============================================
            total partitions
            total waves
            parts in each wave followed by wave indexes. 
    ============================================================================= */
    /* Save a backup to the previous saved check point file */
    sfilename[0] = '\0'; 
    sprintf(sfilename, "../out/cwaves%s", outputfilename);
    sprintf(bkfilename, "../out/cwaves%s_bk", outputfilename);
    file_copy( sfilename, bkfilename, get_bytes(0.0, sfilename), get_bytes(100.0, bkfilename) );

    if (( sfile= fopen (sfilename, "w")) == NULL) {
        mprintf(0, "Can not write checkpoint file[Wave Info].\n", 3);
        return -1;
    }
    /* In the main file, store once, when startPart = 0, the general information, waves info, that doesn't change during computation,  */
    /*1.Write Parts Total */
    ret = fprintf (sfile, "%ld\n", wData->partsTotal);
    if (ret < 0) {
        printf ("[%d]Error (1)writing waves calculations file, Exiting!\n", myProcid);
        return -1;
    }
    /*2.Write Waves Total */
    ret = fprintf (sfile, "%ld\n", wData->wavesTotal);
    if (ret < 0) {
        printf ("[%d]Error (2)writing waves calculations file, Exiting!\n", myProcid);
        return -1;
    }
    for (i=0;i<wData->wavesTotal;i++) {
        /*3.Write Parts Total in this wave*/
        ret = fprintf (sfile, "%ld\n", wData->partsInWave[i]);
        if (ret < 0) {
            printf ("[%d]Error (3)writing waves calculations file, Exiting!\n", myProcid);
            return -1;
        }
        for (j=0;j<wData->partsInWave[i];j++) {	      	
            /*4.Write this Part Index in this wave*/
            for (k=0;k<pData->seqNum;k++) {
                ret = fprintf(sfile, "%lld\n", wData->partsInWaveIndices[i][j][k]);
                if (ret < 0) {
                    printf ("[%d]Error (4)writing waves calculations file, Exiting!	\n", myProcid);
                    sprintf(bkfilename, "../out/cwaves%s_bk", outputfilename);
file_copy( sfilename, bkfilename, get_bytes(0.0, sfilename), get_bytes(100.0, bkfilename) );
                    return -1;
                }
            }
        }
    }
    fclose(sfile);
    return 0;
}

/***************************************************************************
	Function: restoreWavesCalculations
	Description:
		read waves Calculations done by process of rank 0  instead of redoing the calculationg in all processes.
***************************************************************************/
int restoreWavesCalculations (ProcessData * pData, WavesData * wData) {
    FILE * sfile;
    char sfilename[FILENAME_MAX_LENGTH];
    char line[LINE_MAX];
    long i, j;
    MOATypeDimn k;
    int which_proc, firstPartFound = 0;
    char msg[SHORT_MESSAGE_SIZE];
#ifndef NDEBUG
    int dbglevel = 10;
#endif
    /* file cwaves... stores waves info ============================================
            total partitions
            total waves
            parts in each wave followed by wave indexes. 
    ============================================================================= */
    sfilename[0] = '\0'; 
    sprintf(sfilename, "../out/cwaves%s", outputfilename);
    pData->partitionsCount = 0;
    if (( sfile= fopen (sfilename, "r")) == NULL) {
        sprintf(msg, "[%d]>restoreCheckPoint: Invalid checkpoint file, exiting [Wave Info - (%s)]\n", myProcid, sfilename);
        mprintf(0, msg, 3);
        return -1;
    }
    /*1.Read Parts Total */
    if (fgets (line, LINE_MAX, sfile) == NULL) {
        mprintf(0, "Invalid checkpoint file, exiting [reading Wave Info line 1].\n", 3);
        return -1;
    }
    wData->partsTotal = atol (line);
#ifndef NDEBUG
    sprintf(msg, "[%d]>restoreCheckPoint: read Total Parts %ld\n", myProcid, wData->partsTotal);
    mprintf(dbglevel, msg, 3);
#endif
    /*2.Read Waves Total */
    if (fgets (line, LINE_MAX, sfile) == NULL) {
        mprintf(0, "Invalid checkpoint file, exiting [reading Wave Info line 2].\n", 3);
        return -1;
    }
    wData->wavesTotal = atol (line);
#ifndef NDEBUG
    sprintf(msg, "[%d]>restoreCheckPoint: read Total Waves %ld\n", myProcid, wData->wavesTotal);
    mprintf(dbglevel, msg, 3);
#endif
    /* if there are waves, allocate memory for the parts totals array and indices 2D array*/
    if (wData->wavesTotal > 0) {
        wData->partsInWave = mcalloc ((MOATypeInd) wData->wavesTotal, (MOATypeInd) sizeof *(wData->partsInWave));	
        wData->partsInWaveIndices = mmalloc ((MOATypeInd) wData->wavesTotal * (MOATypeInd) sizeof *(wData->partsInWaveIndices));	
        for (i=0;i<wData->wavesTotal;i++) {
            /*3.Read Parts Total in this wave*/
            if (fgets (line, LINE_MAX, sfile) == NULL) {
                mprintf(1, "Invalid checkpoint file, exiting [reading Wave Info line 2+].\n", 3);
                return -1;
            }
            wData->partsInWave[i] = atol(line);
#ifndef NDEBUG
            sprintf(msg, "[%d]>restoreCheckPoint: read partsInWave[%ld] %ld\n", myProcid, i, wData->partsInWave[i]);

    mprintf(dbglevel, msg, 3);
#endif
            if (wData->partsInWave[i] > 0) {
                wData->partsInWaveIndices[i] = mmalloc ((MOATypeInd) wData->partsInWave[i] *  (MOATypeInd) sizeof *(wData->partsInWaveIndices[i]));	
                for (j=0;j<wData->partsInWave[i];j++) {	      	
                    wData->partsInWaveIndices[i][j]= mcalloc ((MOATypeInd) pData->seqNum,  (MOATypeInd) sizeof *(wData->partsInWaveIndices[i][j]));	
                    for (k=0;k<pData->seqNum;k++) {
                        /*4.Read this Part Index in this wave*/
                        if (fgets (line, LINE_MAX, sfile) == NULL) {
                            mprintf(1, "Invalid checkpoint file, exiting [reading Wave Info wave indexes].\n", 3);
                            return -1;
                        }
                        wData->partsInWaveIndices[i][j][k] = atoll(line);
                    }
                    which_proc = getProcID (wData, i, j);
                    if (which_proc == myProcid) {
                        pData->partitionsCount  ++;
                        if (firstPartFound == 0)  {
                            pData->waveNo = i;
                            pData->partNo = j;
                            firstPartFound = 1;
                        }
                    }
#ifndef NDEBUG
                    sprintf(msg, "[%d]>restoreCheckPoint: read waves: %ld, parts: %ld\n", myProcid, i, j);
                    mprintf(dbglevel, msg, 3);
#endif
                } /*End Loop through Parts*/
            } /*End if wave contains parts*/
        } /*End Loop through waves*/
    }/*End Check if waves>0*/
    /*Create & initialize the memory for the OCout for all waves, and delete in freeProcessMemory*/
    pData->OCout = mmalloc(((MOATypeInd) wData->wavesTotal) * ((MOATypeInd) sizeof *(pData->OCout)));
    if ( pData->OCout == NULL) {
        mprintf(1, "Couldn't create memory for OCout waves Buffers. Exiting.\n", 3);
        printf("Couldn't create memory for OCout waves Buffers. Exiting.\n");
        return -1;
    }
    pData->OCin = mmalloc(((MOATypeInd) wData->wavesTotal) * ((MOATypeInd) sizeof *(pData->OCin)));
    if ( pData->OCin == NULL) {
        mprintf(1, "Couldn't create memory for OCin waves Buffers. Exiting.\n", 3);
        printf("Couldn't create memory for OCin waves Buffers. Exiting.\n");
        return -1;
    }
    
    for (i=0;i<wData->wavesTotal;i++) {
        pData->OCout[i].wavesOC = 0;
        pData->OCout[i].WOCO = NULL;
        pData->OCin[i].wavesOC = 0;
        pData->OCin[i].WOCI = NULL;
    }
    fclose(sfile);
#ifndef NDEBUG
    printf ("[%d] read %ld local partitions in wave %ld / %ld, starting at part order %ld index {%lld", myProcid, pData->partitionsCount, pData->waveNo, wData->wavesTotal, pData->partNo, wData->partsInWaveIndices[pData->waveNo][pData->partNo][0]);
    for (k=1;k<pData->seqNum;k++) 
        printf (", %lld", wData->partsInWaveIndices[pData->waveNo][pData->partNo][k]);
    printf ("}\n");
    printf ("[%d] For %ld parts and partitionSize %ld, k = %lld, Distributed Estimated Scoring Total Memory %lld\n", myProcid, pData->partitionsCount, wData->partitionSize,  pData->seqNum, (MOATypeInd) (pData->partitionsCount * ((((MOATypeInd) sizeof (pData->seqNum)) + ((MOATypeInd) pData->seqNum * sizeof (pData->seqLen)) + ((MOATypeInd) sizeof (MOATypeInd)) + ((MOATypeInd) sizeof (MOATypeInd) * ((MOATypeInd) powl(wData->partitionSize, pData->seqNum))) + ((MOATypeInd) sizeof (MOA_elm) * ( (MOATypeInd) powl(wData->partitionSize, pData->seqNum)))))));
    fflush (stdout);
#endif
    return 0;
}

/***************************************************************************
	Function: checkPointPartition
	Description:
		CheckPoint Partition Data into its own file.
***************************************************************************/

int checkPointPartition (ProcessData * pData, ScoringData * sData) {
    FILE * sfile;
    MOATypeInd i, j, k;
    char sfilename[FILENAME_MAX_LENGTH], bkfilename[FILENAME_MAX_LENGTH];
    sprintf(sfilename, "../out/Prt%s%dw%ldp%ld", outputfilename, myProcid, sData->waveNo, sData->partNo);
    sprintf(bkfilename, "../out/Prt%s%dw%ldp%ld_bk", outputfilename, myProcid, sData->waveNo, sData->partNo);
    file_copy( sfilename, bkfilename, get_bytes(0.0, sfilename), get_bytes(100.0, bkfilename) );
    if (( sfile= fopen (sfilename, "w")) == NULL) {
        mprintf(1, "Can not write checkpoint file [Partition Data].\n", 3);
        return -1;
    }
    for (i = 0;i < pData->msaAlgn->elements_ub;i++) {
        /*content 1 ... score value*/
        fprintf (sfile, "%lld\n", pData->msaAlgn->elements[i].val);	
        /*content 2 ... , total previous cells (led to its score) recorded, list of all previous cells global indecies*/
        fprintf (sfile, "%lld\n", pData->msaAlgn->elements[i].prev_ub);	
        for (j = 0; j < pData->msaAlgn->elements[i].prev_ub; j++)  {
             for (k=0;k<pData->seqNum;k++) 
                fprintf (sfile, "%lld\n", pData->msaAlgn->elements[i].prev[j][k]);
        }
    }
    fclose(sfile);
    return 0;
}
/***************************************************************************
	Function: restorePartitionCheckPoint
	Description:
		read specific partition data from its checkpoint file
***************************************************************************/

int restorePartitionCheckPoint (ProcessData * pData, WavesData * wData, long waveNo, long partNo) {
    FILE * sfile;
    char sfilename[FILENAME_MAX_LENGTH];
    char line[LINE_MAX];
    MOATypeInd i, j, k;
#ifndef NDEBUG
    char msg[SHORT_MESSAGE_SIZE];
    int dbglevel = 10;
#endif		
    pData->waveNo = waveNo;
    pData->partNo = partNo;
    sprintf(sfilename, "../out/Prt%s%dw%ldp%ld", outputfilename, myProcid, waveNo, partNo);
    if (( sfile= fopen (sfilename, "r")) == NULL) {
        printf("[%d] Can not read checkpoint file %s.\n", myProcid, sfilename);
        return -1;
    }
    getPartition (wData->partsInWaveIndices[waveNo][partNo], pData->seqNum, pData->seqLen, &pData->msaAlgn, wData->partitionSize);
    if (pData->msaAlgn->elements_ub > 0) {
        /* line ..: Partition All elements score values, total previous cells (led to its score) recorded, list of all previous cells global indecies*/
        for (i = 0;i < pData->msaAlgn->elements_ub;i++) {
            /* msaAlgn->element score value */
            if (fgets (line, LINE_MAX, sfile) == NULL) {
                mprintf(1, "Invalid checkpoint file, exiting 13.\n", 3);
                return -1;
            }	
            /*1. score*/
            pData->msaAlgn->elements[i].val  = atoll(line);
#ifndef NDEBUG
            sprintf(msg, "[%d]>restoreCheckPoint: read MOAPart[%ld]{.msaAlgn->indexes[%lld] %lld, elements[%lld].val %lld}\n", myProcid, partNo, i, pData->msaAlgn->indexes[i], i, pData->msaAlgn->elements[i].val);
            mprintf(dbglevel, msg, 3);
#endif
            /* msaAlgn->element number of prev elements */
            if (fgets (line, LINE_MAX, sfile) == NULL) {
                mprintf(1, "Invalid checkpoint file, number of prev elements.\n", 3);
                return -1;
            }	
            /*2. prev_ub*/
            pData->msaAlgn->elements[i].prev_ub  = atoll(line);
            pData->msaAlgn->elements[i].prev = NULL;

#ifndef NDEBUG
            sprintf(msg, "[%d]>restoreCheckPoint: read MOAPart[%ld].msaAlgn->elements[%lld].prev_ub %lld{", myProcid, partNo, i, pData->msaAlgn->elements[i].prev_ub);
#endif
            if (pData->msaAlgn->elements[i].prev_ub > 0 ) {
                pData->msaAlgn->elements[i].prev = mcalloc ((MOATypeInd) pData->msaAlgn->elements[i].prev_ub, (MOATypeInd) sizeof *pData->msaAlgn->elements[i].prev);
                if (pData->msaAlgn->elements[i].prev == NULL) {
                    mprintf(1, "Can not allocate memory while restoring partition check point. Exiting. 2.\n", 3);
                    return -1;
                }
                for (j = 0; j < pData->msaAlgn->elements[i].prev_ub; j++)  {
                    pData->msaAlgn->elements[i].prev[j] = NULL;
                    pData->msaAlgn->elements[i].prev[j] = mcalloc ((MOATypeInd) pData->seqNum, (MOATypeInd) sizeof *pData->msaAlgn->elements[i].prev[j]);
                    if (pData->msaAlgn->elements[i].prev[j] == NULL) {
                       mprintf(1, "Can not allocate memory while restoring partition check point. Exiting. 2.\n", 3);
                       return -1;
                    }
                    for (k=0;k<pData->seqNum;k++) {
                        if (fgets (line, LINE_MAX, sfile) == NULL) {
                            mprintf(1, "Invalid checkpoint file, prev elements.\n", 3);
                            return -1;
                        }	
                        /*3. prev*/
                        pData->msaAlgn->elements[i].prev[j][k]  = atoll(line);
#ifndef NDEBUG
                        sprintf(msg, "%s%lld", msg, pData->msaAlgn->elements[i].prev[j][k]);
#endif
                    }
                } /* end of prev elements loop */
#ifndef NDEBUG
                sprintf(msg, "%s}\n", msg);
                mprintf(dbglevel, msg, 3);
#endif
            } /* end of prev elements condition */
        } /* end of elements loop */
    } /* end of elements_ub >0 if condition */
    fclose(sfile);
    if (TBFlag == 0) {
        if (pData->waveNo > pData->seqNum+1)
            j = pData->waveNo - pData->seqNum;
        else
            j = 0;
        for (i=j;i<=pData->waveNo;i++) {
            if ((pData->OCout[i].wavesOC == 0) && (pData->OCout[i].WOCO == NULL))
                restoreOCCheckPoint(pData,  i);
        }
    }
    return 0;
}
/***************************************************************************
	Function: checkPoint Current Status File
	Description:
		checkPoint Process Data into files.
***************************************************************************/
int checkPoint (ProcessData * pData, ScoringData * sData) {
    FILE * sfile;
    MOATypeInd i, j, k;
    char sfilename[FILENAME_MAX_LENGTH], bkfilename[FILENAME_MAX_LENGTH];

    /* file chpS... stores computation status =======================================
    line 1: computed partitions
    line 2: current wave number
    line 3: current partition Order
    line 4 - line 4+ 4 * n: total incomming cells (n), (wave no, cell index, cell score, from process)
    line 4 + 4*n + 1 to ...: total outgoing cells, 
    =============================================================================== */
    /* Save a backup to the previous saved check point file */
    sprintf(sfilename, "../out/chpS%s%d", outputfilename, myProcid);
    sprintf(bkfilename, "../out/chpS%s%d_bk", outputfilename, myProcid);
    file_copy( sfilename, bkfilename, get_bytes(0.0, sfilename), get_bytes(100.0, bkfilename) );

    if (( sfile= fopen (sfilename, "w")) == NULL) {
        mprintf(0, "Can not write checkpoint file [computation status].\n", 3);
        return -1;
    }
    /*store the current computation status information, and OCin, and OCout*/	
    /* line 1: Process partitionsCount */
    //fprintf (sfile, "%ld\n", pData->partitionsCount);
    /* line 2: computed partitions */
    fprintf (sfile, "%ld\n", pData->computedPartitions);
    /* line 3: current wave number */
    fprintf (sfile, "%ld\n", pData->waveNo);
    /* line 4: current partition Order */ 
    //fprintf (sfile, "%ld\n", pData->partNo);
    fclose(sfile);
    if (pData->waveNo > sData->waveNo)
        checkPointOC(pData, sData->waveNo);
    if (checkPointPartition (pData, sData) < 0 )
        return -1;

    return 0;
}

int checkPointOC (ProcessData * pData, long waveNo) {
    FILE * sfile;
    MOATypeInd i, j, k;
    char sfilename[FILENAME_MAX_LENGTH], bkfilename[FILENAME_MAX_LENGTH];
    sprintf(sfilename, "../out/oc%sp%dw%ld", outputfilename, myProcid, waveNo);
    if (( sfile= fopen (sfilename, "w")) == NULL) {
        mprintf(0, "Can not write checkpoint file [computation status].\n", 3);
        return -1;
    }
    i = waveNo;
    //for (i=0;i<pData->OCin_ub;i++) {
        /*lines ...: Total Overlapping Cells for this wave*/
        fprintf (sfile, "%ld\n",pData->OCin[i].wavesOC);
        /*lines ...: Overlapping Incoming Cells cellIndex, cellScore, fromProc*/
        for (j=0;j<pData->OCin[i].wavesOC;j++) {
            for (k=0;k<pData->seqNum;k++) 
                fprintf (sfile, "%lld\n", pData->OCin[i].WOCI[j].cellIndex[k]);
            fprintf (sfile, "%lld\n", pData->OCin[i].WOCI[j].cellScore);	
        }
    //}
    /* if still computing scores then save extra info to resume*/
    /*if (slave->ComputationPhase < 2) {*/
    /*line x: Overlapping Outgoing Cells Scores*/
    //fprintf (sfile, "%ld\n", pData->OCout_ub);
    //for (i=0;i<pData->OCout_ub;i++) {
        /*lines ...: Total Overlapping Outgoing for this wave*/
        fprintf (sfile, "%ld\n",pData->OCout[i].wavesOC); 
        for (j=0;j<pData->OCout[i].wavesOC;j++) {
                /*lines ...: Overlapping Outgoing Cells partIndex, cellIndex, cellScore, sent flag, totaldependent Proccessors, and list of Processors*/
            for (k=0;k<pData->seqNum;k++) 
                fprintf (sfile, "%lld\n", pData->OCout[i].WOCO[j].cellIndex[k]);
            fprintf (sfile, "%lld\n", pData->OCout[i].WOCO[j].cellScore);	
            fprintf (sfile, "%d\n", pData->OCout[i].WOCO[j].sent);	
            fprintf (sfile, "%d\n", pData->OCout[i].WOCO[j].depProc_ub);	
            for (k=0;k<pData->OCout[i].WOCO[j].depProc_ub;k++) {
                fprintf (sfile, "%d\n", pData->OCout[i].WOCO[j].depProc[k]);	
            }
        }
    //}
    fclose(sfile);
    return 0;    
}

/***************************************************************************
	Function: restoreCheckPoint Current Status CheckPoint File
	Description:
		read process data from checkpoint files
***************************************************************************/
int restoreOCCheckPoint (ProcessData * pData, long waveNo) {
    FILE * sfile;
    char sfilename[FILENAME_MAX_LENGTH];
    char line[LINE_MAX];
    long i, j;
    MOATypeDimn k;
#ifndef NDEBUG
    char msg[SHORT_MESSAGE_SIZE];
    int dbglevel = 10;
#endif
    sprintf(sfilename, "../out/oc%sp%dw%ld", outputfilename, myProcid, waveNo);

    if (( sfile= fopen (sfilename, "r")) == NULL) {
        mprintf(1, "Invalid checkpoint file, exiting [reading computation status].\n", 3);
        return -1;
    }

    i = waveNo;
    if (pData->OCin != NULL) {
        if ((pData->OCin[i].WOCI != NULL) && (pData->OCin[i].wavesOC > 0)) {
            for (j=0;j<pData->OCin[i].wavesOC;j++) {
                if (pData->OCin[i].WOCI[j].cellIndex != NULL) {
                    free (pData->OCin[i].WOCI[j].cellIndex);
                    pData->OCin[i].WOCI[j].cellIndex = NULL;
                }
            }
            free (pData->OCin[i].WOCI);                                  
            pData->OCin[i].WOCI = NULL;
        }
        pData->OCin[i].wavesOC = 0;
    }
    if (pData->OCout != NULL) {
        if ((pData->OCout[i].WOCO != NULL) && (pData->OCout[i].wavesOC > 0)) {
            for (j=0;j<pData->OCout[i].wavesOC;j++) {
                if (pData->OCout[i].WOCO[j].cellIndex != NULL) {
                    free (pData->OCout[i].WOCO[j].cellIndex);
                    pData->OCout[i].WOCO[j].cellIndex = NULL;
                }
                if ((pData->OCout[i].WOCO[j].depProc_ub > 0) && (pData->OCout[i].WOCO[j].depProc  != NULL)) {
                    free(pData->OCout[i].WOCO[j].depProc);	
                    pData->OCout[i].WOCO[j].depProc = NULL;
                }
            }
            free (pData->OCout[i].WOCO);                                  
            pData->OCout[i].WOCO = NULL;
        }
        pData->OCout[i].wavesOC = 0;
    }
            
//for (i=0;i<pData->OCin_ub;i++) {
    /*lines ...: Total Overlapping Cells for this wave*/
    if (fgets (line, LINE_MAX, sfile) == NULL) {
        mprintf(0, "Invalid checkpoint file, reading wave number.\n", 3);
        return -1;
    }
    pData->OCin[i].wavesOC = atol(line);
    pData->OCin[i].WOCI = mmalloc (((MOATypeInd) pData->OCin[i].wavesOC) * ((MOATypeInd) sizeof *(pData->OCin[i].WOCI)));
    if (pData->OCin[i].WOCI == NULL) {
        printf("Failed to allocate memory for OCin at wave %ld while restoring files.\n", i);
        return -1;
    }
    for (j=0;j<pData->OCin[i].wavesOC;j++) {
        pData->OCin[i].WOCI[j].cellIndex = mcalloc ((MOATypeInd) pData->seqNum, (MOATypeInd) sizeof * pData->OCin[i].WOCI[j].cellIndex);
        /*lines ...: Overlapping Incoming Cells cellIndex, cellScore, fromProc*/
        /* cell index */
        for (k=0;k<pData->seqNum;k++) {
            if (fgets (line, LINE_MAX, sfile) == NULL) {
                mprintf(0, "Invalid checkpoint file, reading cell index.\n", 3);
                return -1;
            }
            pData->OCin[i].WOCI[j].cellIndex[k] = atoll(line);
        }
        /* cell score */
        if (fgets (line, LINE_MAX, sfile) == NULL) {
            mprintf(0, "Invalid checkpoint file, reading cell score.\n", 3);
            return -1;
        }
        pData->OCin[i].WOCI[j].cellScore = atoll(line) ;
    }
// }


// for (i=0;i<pData->OCout_ub;i++) {
    /*lines ...: Total Overlapping Outgoing for this wave*/
    if (fgets (line, LINE_MAX, sfile) == NULL) {
            mprintf(0, "Invalid checkpoint file, reading wave number.\n", 3);
            return -1;
    }
    pData->OCout[i].wavesOC = atol(line);
    pData->OCout[i].WOCO = mmalloc (((MOATypeInd) pData->OCout[i].wavesOC) * ((MOATypeInd) sizeof *(pData->OCout[i].WOCO)));
    if (pData->OCout[i].WOCO == NULL) {
        printf("Failed to allocate memory for OCout at wave %ld while restoring files.\n", i);
        return -1;
    }
    for (j=0;j<pData->OCout[i].wavesOC;j++) {
        pData->OCout[i].WOCO[j].cellIndex = mcalloc ((MOATypeInd) pData->seqNum, ((MOATypeInd) sizeof *pData->OCout[i].WOCO[j].cellIndex));
        /*lines ...: Overlapping Outgoing Cells cellIndex, cellScore, sent flag, totaldependent Proccessors, and list of Processors*/

        // cellIndex 
        for (k=0;k<pData->seqNum;k++) {
            if (fgets (line, LINE_MAX, sfile) == NULL) {
                    mprintf(1, "Invalid checkpoint file, reading cell index.\n", 3);
                    return -1;
            }
            pData->OCout[i].WOCO[j].cellIndex[k] = atoll(line);
        }
        // cellScore 
        if (fgets (line, LINE_MAX, sfile) == NULL) {
                mprintf(1, "Invalid checkpoint file, reading cell score.\n", 3);
                return -1;
        }
        pData->OCout[i].WOCO[j].cellScore = atoll(line) ;
        // sent 
        if (fgets (line, LINE_MAX, sfile) == NULL) {
                mprintf(1, "Invalid checkpoint file, reading sent switch.\n", 3);
                return -1;
        }
        pData->OCout[i].WOCO[j].sent = atoi(line);
        if (RestoreFlag) 
                pData->OCout[i].WOCO[j].sent = 0;
        // depProc_ub 
        if (fgets (line, LINE_MAX, sfile) == NULL) {
                mprintf(1, "Invalid checkpoint file, reading number of dependent processes.\n", 3);
                return -1;
        }
        pData->OCout[i].WOCO[j].depProc_ub = atoi(line) ;
#ifndef NDEBUG
        sprintf(msg, "[%d]>restoreCheckPoint: read OCout[%ld](%ld, {%lld", myProcid, i, j, pData->OCout[i].WOCO[j].cellIndex[0]);
        mprintf(dbglevel, msg, 3);
        for (k=1;k<pData->seqNum;k++) {
            sprintf(msg, ", %lld", pData->OCout[i].WOCO[j].cellIndex[k]);
            mprintf(dbglevel, msg, 3);
        }
        sprintf(msg, "}, %lld, %d, %d}\n", pData->OCout[i].WOCO[j].cellScore, pData->OCout[i].WOCO[j].sent, pData->OCout[i].WOCO[j].depProc_ub);
        mprintf(dbglevel, msg, 3);
#endif
        if (pData->OCout[i].WOCO[j].depProc_ub > 0) {
            pData->OCout[i].WOCO[j].depProc = NULL;
            pData->OCout[i].WOCO[j].depProc = mmalloc (((MOATypeInd) pData->OCout[i].WOCO[j].depProc_ub) * ((MOATypeInd) sizeof *(pData->OCout[i].WOCO[j].depProc)));
            if (pData->OCout[i].WOCO[j].depProc == NULL) {
                printf("Failed to allocate memory for OCout %d depProc at wave %ld while restoring files.\n", pData->OCout[i].WOCO[j].depProc_ub, i);
                return -1;
            }
            for (k=0;k<pData->OCout[i].WOCO[j].depProc_ub;k++) {
                if (fgets (line, LINE_MAX, sfile) == NULL) {
                    mprintf(1, "Invalid checkpoint file, reading dependent process.\n", 3);
                    return -1;
                }
                pData->OCout[i].WOCO[j].depProc[k] = atoi(line);
#ifndef NDEBUG
                sprintf(msg, "%s%d", msg, pData->OCout[i].WOCO[j].depProc[k]);
#endif
            }
        }		
#ifndef NDEBUG	
        sprintf(msg, "%s}\n", msg);
        mprintf(dbglevel, msg, 3);
#endif
    }// end loop OCout in waves
//    } // end loop OCout waves
  
    fclose(sfile);
    return 0;
}
int restoreCheckPoint (ProcessData * pData, WavesData * wData) {
    FILE * sfile;
    char sfilename[FILENAME_MAX_LENGTH];
    char line[LINE_MAX];
    long i, j;
    MOATypeDimn k;
    ldiv_t partsInCluster;
    long break_point;   
    char msg[SHORT_MESSAGE_SIZE];
#ifndef NDEBUG
    int dbglevel = 10;
#endif

	 
    if (((Mode == Distributed) || (ClusterSize > 1)) &&  (restoreWavesCalculations (pData, wData) < 0)) {
        printf ("Failed to restore Wave Calculations. Exiting Restore check points\n");
        fflush (stdout);        
        return -1;
    }

  
    /* file chpS... read computation status =======================================
    line 1: computed partitions
    line 2: current wave number
    line 3: current partition Order
    line 4 - line 4+ 4 * n: total incommeng cells (n), (wave no, cell index, cell score, from process)
    line 4 + 4*n + 1 to ...: total outgoing cells, 
    =============================================================================== */
    sprintf(sfilename, "../out/chpS%s%d", outputfilename, myProcid);

    if (( sfile= fopen (sfilename, "r")) == NULL) {
        mprintf(1, "Invalid checkpoint file, exiting [reading computation status].\n", 3);
        return -1;
    }
    /* line 1: Process partitionsCount 
    if (fgets (line, LINE_MAX, sfile) == NULL) {
        sprintf(msg, "[%d]>restoreCheckPoint: exiting [reading Wave Info line 1] - file (%s)\n", myProcid, sfilename);
        mprintf(0, msg, 3);
        return -1;
    }
    pData->partitionsCount = atol (line);
#ifndef NDEBUG
    sprintf(msg, "[%d]>restoreCheckPoint: read partitionsCount %ld\n", myProcid, pData->partitionsCount);
    mprintf(dbglevel, msg, 3);
#endif*/
    /* line 2: computed partitions */
    if (fgets (line, LINE_MAX, sfile) == NULL) {
        mprintf(1, "Invalid checkpoint file, exiting [reading Computed Partitions].\n", 3);
        return -1;
    }
    pData->computedPartitions = atol(line);
     
#ifndef NDEBUG
    sprintf(msg, "[%d]>restoreCheckPoint: read computedPartitions %ld\n", myProcid, pData->computedPartitions);
    mprintf(dbglevel, msg, 3);
#endif
    /* line 3: current wave number */
    if (fgets (line, LINE_MAX, sfile) == NULL) {
        mprintf(1, "Invalid checkpoint file, exiting [reading current wave number].\n", 3);
        return -1;
    }
    pData->globalWaveNo =  pData->waveNo = atol(line);
#ifndef NDEBUG	
    sprintf(msg, "[%d]>restoreCheckPoint: read currWaveNo %ld\n", myProcid, pData->waveNo);
    mprintf(dbglevel, msg, 3);
#endif
   
    fclose(sfile);

    if ((RestoreFlag == 1) && (TBFlag != 1)) { 
        if (myProcid == 0)
            pData->partNo = 0;
        else {
            partsInCluster = ldiv(wData->partsInWave[pData->waveNo], ClusterSize);
            if ((wData->partsInWave[pData->waveNo] < ClusterSize) && (partsInCluster.quot == 0) && (myProcid >= partsInCluster.rem))
                pData->partNo = -1;         
            else {
                if (partsInCluster.rem > 0) {
                    if (myProcid < partsInCluster.rem) {
                        pData->partNo = (partsInCluster.quot + 1) * myProcid;
                    } else {
                        break_point = partsInCluster.rem * (partsInCluster.quot + 1);
                        pData->partNo = partsInCluster.quot * (myProcid - partsInCluster.rem) + break_point;
                    }
                } else {
                    pData->partNo = partsInCluster.quot * myProcid;
                }
            }            
        }
    } 
    /*else if (TBFlag == 1) {
        pData->globalWaveNo = pData->waveNo = 0;
    }*/
    else { //if (TBFlag == 1) {
        pData->globalWaveNo = pData->waveNo = 0;
        if (myProcid == 0)
            pData->partNo = 0;
        else {
            getNextPartition (wData, &pData->waveNo, &pData->partNo);
        }
    }

    if (RestoreFlag == 1) {
        if (pData->waveNo > pData->seqNum+1)
            j = pData->waveNo - pData->seqNum;
        else
            j = 0;
        for (i=j;i<=pData->waveNo;i++) {
            if ((pData->OCout[i].wavesOC == 0) && (pData->OCout[i].WOCO == NULL))
                restoreOCCheckPoint(pData,  i);
        }
    }
    else  { //if (TBFlag != 1)
        if (restorePartitionCheckPoint(pData, wData, pData->waveNo, pData->partNo) < 0) {
            return -1;
        }
    }
    return 0;
}

/* ======================================================================
	function calcWaves by MoA Dependency:
	Input:
		dimn: dimension of array
		shape: array of lengths of each dimension.
		partSize: size of partition.
	Output:
		myPartsCount: number of partitions.
		myCurrWave: current wave to be processed.
		myCurrPart: current partition to be processed.
======================================================================= */
long double  getPythagoreanDistance (MOATypeShape * x_point, MOATypeShape * y_point, MOATypeDimn dimn) {
    long double distance;
    MOATypeDimn k; /*Multidimensional Index Iterator*/
    if (dimn <= 0)
        return -999;
    distance = powl((long double) y_point[0] - (long double) x_point[0], (long double)  2);
    for (k=1;k<dimn;k++)
        distance += powl((long double) y_point[k] - (long double) x_point[k], (long double)  2);
    
    distance = sqrtl (distance);
    return distance;
}

long double  getPythagoreanOriginDistance (MOATypeShape * point, MOATypeDimn dimn) {
    long double distance;
    MOATypeDimn k; /*Multidimensional Index Iterator*/
    if (dimn <= 0)
        return -999;
    distance = powl((long double) point[0], (long double)  2);
    for (k=1;k<dimn;k++) {
        distance += powl((long double) point[k], (long double)  2);
    }
    distance = sqrtl (distance);
    return distance;
}

long double  getPythagoreanOriginDirection (MOATypeShape * point, MOATypeDimn dimn) {
    long double distance, direction;
    MOATypeDimn k; /*Multidimensional Index Iterator*/
    if (dimn <= 0)
        return -999;
    distance = getPythagoreanOriginDistance(point, dimn);
    direction = point[0];
    for (k=1;k<dimn;k++) 
        direction += point[k];
    direction = direction / distance;
    return direction;
}

long double  getPythagoreanDirection (MOATypeShape * x_point, MOATypeShape * y_point, MOATypeDimn dimn) {
    long double x_distance, y_distance, direction;
    MOATypeDimn k; /*Multidimensional Index Iterator*/
    if (dimn <= 0)
        return -999;
    x_distance = getPythagoreanOriginDistance(x_point, dimn);
    y_distance = getPythagoreanOriginDistance(y_point, dimn);
    direction = x_point[0] * y_point[0];
    for (k=1;k<dimn;k++) 
        direction += x_point[k] * y_point[k];
    direction = direction / (x_distance * y_distance);
    return direction;
}

long double  getVectorsAngle (MOATypeShape * x_point, MOATypeShape * y_point, MOATypeDimn dimn) {
    long double x_distance, y_distance, angle;
    int sign = -2; /*initialized as not yet read*/
    MOATypeDimn k; /*Multidimensional Index Iterator*/
    if (dimn <= 0)
        return -999;
    x_distance = getPythagoreanOriginDistance(x_point, dimn);
    y_distance = getPythagoreanOriginDistance(y_point, dimn);
    angle = x_point[0] * y_point[0];
    if (x_point[0] < y_point[0])
        sign = 1;
    else if (x_point[0] > y_point[0])
        sign = -1;
    
    for (k=1;k<dimn;k++) {
        angle += x_point[k] * y_point[k];
        if (sign == -2) {
            if (x_point[k] < y_point[k])
                sign = 1;
            else if (x_point[k] > y_point[k])
                sign = -1;
        }
    }
    if (sign == -2)  /*if all elements in both vectors are equal, and no sign defined so far*/
        sign = 1;    
    angle = acos (angle / (x_distance * y_distance)) * sign;
    return angle;
}

long double  getVectorsRankedAngle (MOATypeShape * x_point, MOATypeShape * y_point, MOATypeDimn dimn) {
    long double x_distance, y_distance, angle;
    int sign = -2; /*initialized as not yet read*/
    MOATypeDimn k; /*Multidimensional Index Iterator*/
    if (dimn <= 0)
        return -999;
    x_distance = getPythagoreanOriginDistance(x_point, dimn);
    y_distance = getPythagoreanOriginDistance(y_point, dimn);
    //angle = ((long double) (x_point[0]+1) * (y_point[0]+1)) * (powl(dimn, 2) * ((long double)(y_point[0]+1)));
    angle = ((long double) (x_point[0]+1) * (y_point[0]+1)) * (powl(dimn, 2) + dimn);
    if (x_point[0] < y_point[0])
        sign = 1;
    else if (x_point[0] > y_point[0])
        sign = -1;
    
    for (k=1;k<dimn;k++) {
        //angle += ((long double)(x_point[k]+1) * (y_point[k]+1)) * (powl(dimn-k, 2) * ((long double) (y_point[k]+1)));
        angle += ((long double)(x_point[k]+1) * (y_point[k]+1)) * (powl(dimn-k, 2) + (dimn-k));
        if (sign == -2) {
            if (x_point[k] < y_point[k])
                sign = 1;
            else if (x_point[k] > y_point[k])
                sign = -1;
        }
    }
    if (sign == -2)  /*if all elements in both vectors are equal, and no sign defined so far*/
        sign = 1;    
    angle = (angle / (x_distance * y_distance)) * sign;
    return angle;
}



int  getVectorsProjection (MOATypeShape * x_point, MOATypeShape * y_point, MOATypeDimn dimn, MOATypeShape * * projection) {
    long double x_Dot_u = 0, y_distance;
    MOATypeDimn k; /*Multidimensional Index Iterator*/
    if (dimn <= 0)
        return -1;
    y_distance = getPythagoreanOriginDistance(y_point, dimn);
    for (k=0;k<dimn;k++) {
        (*projection)[k] = y_point[k] / y_distance;
        x_Dot_u += x_point[k] * (*projection)[k];
    }
    for (k=0;k<dimn;k++) 
        (*projection)[k] = x_Dot_u * (*projection)[k];
    return 0;
}

int  getVectorsNormal (MOATypeShape * * points, MOATypeDimn dimn, MOATypeShape * * projection) {
    MOATypeDimn k, l; /*Multidimensional Index Iterator - Verify this function check cross product in MoA*/
    if (dimn <= 0)
        return -1;
    for (k=0;k<dimn;k++) {
        (*projection)[k] = 1;
        for (l=1;l<dimn;l++) {
            (*projection)[k] = (*projection)[k] * (points[k][l] - points[k][0]);
        }
    }
    return 0;
}

MOATypeShape isVectorInHyperPlane (MOATypeShape * x_point, MOATypeShape * y_point, MOATypeDimn dimn, MOATypeShape * z_point) {
    int ret = -1;
    MOATypeDimn k; /*Multidimensional Index Iterator*/
    /*Testing if hytperplane defined in points x and y points, includes also z point, when the x.(z-y) = 0, 
     * then z is a point in the same hyperplane
     instead of testing with the zero perpendicular, which doesn't say anything, fix the getNormal (cross Product)
     *above, and test with it      
     */
    MOATypeShape x_dot_z_minus_y;
    x_dot_z_minus_y = 0;
    for (k=0;k<dimn;k++) {
        x_dot_z_minus_y += x_point[k] * (z_point[k] - y_point[k]);
    }
    return x_dot_z_minus_y;
}

void getFirstPartitionInWave (MOATypeShape * shape, MOATypeDimn dimn, long partitionSize, long waveNo, MOATypeShape * * partIndex) {
    MOATypeDimn k; /*Multidimensional Index Iterator*/
    long remDist = waveNo * (partitionSize-1);
    for (k=0;k<dimn;k++) {
        if (shape[k] - partitionSize > remDist)
            (*partIndex)[k] = remDist;
        else
            (*partIndex)[k] = shape[k] - partitionSize;
        remDist -= (*partIndex)[k];
    }
}

void getSecondPartitionInWave (MOATypeShape * shape, MOATypeDimn dimn, long partitionSize, long waveNo, MOATypeShape * firstPartition, MOATypeShape * * partIndex) {
    MOATypeDimn k; /*Multidimensional Index Iterator*/
    int done = -1;
    for (k=0;k<dimn;k++)
        (*partIndex)[k] = firstPartition[k];
    for (k=dimn-2;(k>=0) && (done == -1);k--) {
        if ((firstPartition[k] > partitionSize - 1) && (firstPartition[k+1] < shape[k+1] - partitionSize)) {
            (*partIndex)[k] = firstPartition[k] -  partitionSize + 1;
            (*partIndex)[k+1] = firstPartition[k+1] +  partitionSize - 1;
            done = 0;
        }       
    }
}

void getLastPartitionInWave (MOATypeShape * shape, MOATypeDimn dimn, long partitionSize, long waveNo, MOATypeShape * * partIndex) {
    MOATypeDimn k; /*Multidimensional Index Iterator*/
    long remDist = waveNo * (partitionSize-1);
    for (k=dimn-1;k>=0;k--) {
        if (shape[k] - partitionSize > remDist)
            (*partIndex)[k] = remDist;
        else
            (*partIndex)[k] = shape[k] - partitionSize;
        remDist -= (*partIndex)[k];
    }
}

void getMiddlePartitionInWave (MOATypeShape * shape, MOATypeDimn dimn, long partitionSize, long waveNo, MOATypeShape * * partIndex) {
    MOATypeDimn remDimn = dimn; 
    MOATypeShape remDist = waveNo; 
    MOATypeDimn startDimn = (MOATypeDimn) floorl((long double) dimn/2);
    MOATypeShape startDist = (MOATypeShape) ceill((long double) remDist/remDimn);
    MOATypeDimn i;

    if (shape[startDimn] > startDist * (partitionSize-1))
        (*partIndex)[startDimn] = startDist;
    else
        (*partIndex)[startDimn] = (shape[startDimn] - partitionSize) / (partitionSize-1);
    remDist = remDist-(*partIndex)[startDimn];
    remDimn --;		
    startDist = (MOATypeShape) ceill((long double) remDist/remDimn);
    if (isEven(dimn) == 1) {
            if (startDimn == 0)
                startDist = remDist;
            
            if (shape[startDimn] > startDist * (partitionSize-1))
                (*partIndex)[dimn-startDimn-1] = startDist;
            else
                (*partIndex)[dimn-startDimn-1] = (shape[startDimn] - partitionSize) / (partitionSize-1);
            
            remDist = remDist-(*partIndex)[dimn-startDimn-1];
            remDimn --;		
            startDist = (MOATypeShape) ceill((long double) remDist/remDimn);
            startDimn --;
    }
    startDimn--;
    for (i=startDimn;i>=0;i--) { 
        if (shape[i] > startDist * (partitionSize-1))
            (*partIndex)[i] = startDist;        
        else
            (*partIndex)[i] = (shape[i] - partitionSize) / (partitionSize-1);        
        remDist = remDist-(*partIndex)[i];
        remDimn --;		
        startDist = (MOATypeShape) ceill((long double) remDist/remDimn);
        if (i == 0)
            startDist = remDist;
        if (shape[dimn-i-1] > startDist * (partitionSize-1))
            (*partIndex)[dimn-i-1] = startDist;
        else
            (*partIndex)[dimn-i-1] = (shape[dimn-i-1] - partitionSize) / (partitionSize-1);
        remDist = remDist-(*partIndex)[dimn-i-1];
        remDimn --;
        if (remDimn > 0)		 
            startDist = (MOATypeShape) ceill((long double) remDist/remDimn);
        else
            startDist = 0;
    }
    
    for (i=0;i<dimn;i++) 
        (*partIndex)[i] *= (partitionSize-1);
}


double getWaveLength (MOATypeShape * shape, MOATypeDimn dimn, long partitionSize, long waveNo) {
    double waveLength;
    MOATypeShape * firstPartition = NULL;
        
    firstPartition = mcalloc (dimn, sizeof *firstPartition);
    getFirstPartitionInWave (shape, dimn, partitionSize, waveNo, &firstPartition);
    waveLength = (double) getVectorsRankedAngle (firstPartition, firstPartition, dimn);
    if (firstPartition != NULL)
        free(firstPartition);
    
    return waveLength;
}

int getPartitionProcID (MOATypeShape * partIndex, MOATypeShape * shape, MOATypeDimn dimn, long waveNo, long partitionSize) {
    MOATypeDimn k; /*Multidimensional Index Iterator*/
    long break_point;
    int procID = -1;
    MOATypeShape * firstPartition = NULL;
    MOATypeShape * lastPartition = NULL;
                
    if (waveNo == 0) 
        procID = 0;
    else {
        firstPartition = mcalloc (dimn, sizeof *firstPartition);
        lastPartition = mcalloc (dimn, sizeof *lastPartition);
        getFirstPartitionInWave (shape, dimn, partitionSize, waveNo, &firstPartition);
        getLastPartitionInWave (shape, dimn, partitionSize, waveNo, &lastPartition);
    
        double partitionAngle = (double) getVectorsRankedAngle (partIndex, firstPartition, dimn);
        double smallestAngle = (double) getVectorsRankedAngle (lastPartition, firstPartition, dimn);
        procID = (int)  ceil( (((double) ClusterSize) * partitionAngle)/ smallestAngle);
    
        if (firstPartition != NULL)
            free(firstPartition);
        if (lastPartition != NULL)
            free(lastPartition);
    }
    procID = procID % ClusterSize;
    return procID;
}

long double  getMoADistance (MOATypeShape * x_point, MOATypeShape * y_point, MOATypeDimn dimn) {
    long double distance = 0;
    MOATypeDimn k; /*Multidimensional Index Iterator*/
    int sign = -2;
    if (dimn <= 0)
        return -999;
    for (k=0;k<dimn;k++) {
        distance +=((long double) llabs((long long) powl((long double) y_point[k]- x_point[k], (long double) 2))) * (powl(dimn-k, 2));
         if (sign == -2) {
            if (y_point[k] < x_point[k])
                sign = 1;
            else if (y_point[k] > x_point[k])
                sign = -1;
        }
    }
    if (sign == -2)  /*if all elements in both vectors are equal, and no sign defined so far*/
        sign = 1;    
    distance = distance * sign;
    return distance;
}

long double getMoAWaveLength (WavesData * wData, long waveNo) {
    long double waveLength;
    MOATypeShape * firstPartition = NULL;
        
    firstPartition = mcalloc (wData->seqNum, sizeof *firstPartition);
    getFirstPartitionInWave (wData->seqLen, wData->seqNum, wData->partitionSize, waveNo, &firstPartition);
    waveLength = getMoADistance (firstPartition, wData->waveMiddle[waveNo], wData->seqNum);
    if (firstPartition != NULL)
        free(firstPartition);
    
    return waveLength;
}

long double  distanceFromMiddle(WavesData * wData, MOATypeShape * partIndex, long waveNo) {
    long double  distance = getMoADistance (partIndex, wData->waveMiddle[waveNo], wData->seqNum);    
    return distance;
}

int isPartWithinEpsilonPercentage(MOATypeShape * partIndex, WavesData * wData, long waveNo) {
    long double partDistance = distanceFromMiddle(wData, partIndex, waveNo);
    if ((waveNo == 0) || (waveNo == wData->wavesTotal-1))
        return 0;
    float ps = (float) partDistance / (float) wData->waveLength[waveNo] * 100;
    if ( labs((long) ps) < Epsilons)
        return 0;
    else
        return -1;
}
int isPartWithinEpsilon(MOATypeShape * partIndex, WavesData * wData, long waveNo) {
    MOATypeDimn k; /*Multidimensional Index Iterator*/
    int validPartition = 0;
    long indicesSpan = Epsilons * (wData->partitionSize - 1);
    if ((waveNo == 0) || (waveNo == wData->wavesTotal-1))
        return 0;
    for (k=0;k<wData->seqNum && validPartition == 0;k++) {
        if ((Epsilons != 0) && ((partIndex[k] < wData->waveMiddle[waveNo][k] - indicesSpan) || (partIndex[k] > wData->waveMiddle[waveNo][k] + indicesSpan)))
            validPartition = -1;
    }
    return validPartition;
}


int getWavesParts (ProcessData * pData, WavesData * wData, MOATypeShape * dist, long * waveNo) {
    MOA_rec * HParts = NULL;
    long j, HPWaveNo;
    int which_proc, * validPartitions = NULL;
    MOATypeInd p, i, nextParts;

    MOATypeDimn k; /*Multidimensional Index Iterator*/
    MOATypeDimn CalLnCount = (MOATypeDimn) mpow(2, pData->seqNum) - 1;
    validPartitions = mmalloc (CalLnCount * sizeof *validPartitions);
    nextParts = MOAGetHigherPartitions (wData->partitionSize, 2, dist, pData->seqNum, pData->seqLen, &HParts, &validPartitions);
    /* Breadth search at each node will span k waves */
    //for (i=1;i<nextParts;i++) { /*Loop through Higher Neighboring partitions*/
    for (p =  floorl((long double) nextParts/2);p>=0;p--) { /*Loop through Higher Neighboring partitions starting from the middle*/
        /* Check left side of the middle*/
        i = p;
        if (validPartitions[i] == 1) { /*check valid partitions*/
            /*retrieve partition order in the wave to define its processor*/
            HPWaveNo = getWaveNo (wData->seqNum, wData->partitionSize, HParts->indexes[i]);

            if (HPWaveNo == -1) 
                continue;
            if ((*waveNo) < HPWaveNo) {
                wData->partsInWave = realloc (wData->partsInWave, ((MOATypeInd) (HPWaveNo+1)) * ((MOATypeInd) sizeof *wData->partsInWave));
                if (wData->partsInWave == NULL) {
                    mprintf(1, "Couldn't re-allocate memory for Number of Partitions per wave vector to Calculate waves.\n", 3);
                    printf("Couldn't re-allocate memory for Number of Partitions per wave vector to Calculate waves.\n");		
                    return -1;
                }
                wData->partsInWaveIndices  = realloc (wData->partsInWaveIndices, ((MOATypeInd) (HPWaveNo+1)) * ((MOATypeInd) sizeof *wData->partsInWaveIndices));
                if (wData->partsInWaveIndices == NULL) {
                    mprintf(1, "Couldn't re-allocate memory for partitions Indices matrix to Calculate waves.\n", 3);
                    printf("Couldn't re-allocate memory for partitions Indices matrix to Calculate waves.\n");		
                    return -1;
                }
                wData->waveLength = realloc (wData->waveLength, ((MOATypeInd) (HPWaveNo+1)) * ((MOATypeInd) sizeof *wData->waveLength));
                wData->waveMiddle  = realloc (wData->waveMiddle, ((MOATypeInd) (HPWaveNo+1)) * ((MOATypeInd) sizeof *wData->waveMiddle));
                if (wData->waveMiddle == NULL) {
                    mprintf(1, "Couldn't re-allocate memory for waveMiddle.\n", 3);
                    printf("Couldn't re-allocate memory for waveMiddle.\n");		
                    return -1;
                }                
                for (k=(*waveNo)+1;k<=HPWaveNo;k++) {
                    wData->waveMiddle[k] = mmalloc (((MOATypeInd) pData->seqNum) * ((MOATypeInd) sizeof *wData->waveMiddle[k]));
                    if (wData->waveMiddle[k] == NULL) {
                        mprintf(1, "Couldn't re-allocate memory for waveMiddle Indices.\n", 3);
                        printf("Couldn't re-allocate memory for waveMiddle Indices.\n");		
                        return -1;
                    }                
                    getMiddlePartitionInWave (wData->seqLen, wData->seqNum, wData->partitionSize, k, &wData->waveMiddle[k]) ;
                    wData->waveLength[k] = getMoAWaveLength(wData, k);
                    wData->partsInWave[k] = 0;
                    wData->partsInWaveIndices[k] = NULL;
                }
                (*waveNo) = HPWaveNo;
            }
            //if ((Epsilons == 0) || ((Epsilons != 0) && (wData->partsInWave[HPWaveNo] < Epsilons))) {  
            if ((Epsilons == 0) || ((Epsilons != 0) && (isPartWithinEpsilon(HParts->indexes[i], wData, HPWaveNo) == 0))) {  
                if (isPartitionInWave (wData, HPWaveNo, HParts->indexes[i]) == 0) {
                    if (wData->partsInWaveIndices[HPWaveNo] == NULL) 
                        wData->partsInWaveIndices[HPWaveNo] = mmalloc (((MOATypeInd) (wData->partsInWave[HPWaveNo]+1)) * ((MOATypeInd) sizeof *(wData->partsInWaveIndices[HPWaveNo])));
                    else
                        wData->partsInWaveIndices[HPWaveNo] = realloc (wData->partsInWaveIndices[HPWaveNo], ((MOATypeInd) (wData->partsInWave[HPWaveNo]+1)) * ((MOATypeInd) sizeof *(wData->partsInWaveIndices[HPWaveNo])));
                    if (wData->partsInWaveIndices[HPWaveNo] == NULL) {
                        mprintf(1, "Couldn't re-allocate memory for partitions Indices matrix to Calculate waves 2.\n", 3);
                        printf("Couldn't re-allocate memory for partitions Indices matrix to Calculate waves 2.\n");		
                        return -1;
                    }
                    wData->partsInWaveIndices[HPWaveNo][wData->partsInWave[HPWaveNo]] = NULL;
                    wData->partsInWaveIndices[HPWaveNo][wData->partsInWave[HPWaveNo]] = mmalloc (((MOATypeInd) pData->seqNum) * ((MOATypeInd) sizeof *(wData->partsInWaveIndices[HPWaveNo][wData->partsInWave[HPWaveNo]])));
                    if (wData->partsInWaveIndices[HPWaveNo][wData->partsInWave[HPWaveNo]] == NULL) {
                        mprintf(1, "Couldn't re-allocate memory for partitions Indices matrix to Calculate waves 3.\n", 3);
                        printf("Couldn't re-allocate memory for partitions Indices matrix to Calculate waves 3.\n");		
                        return -1;
                    }
                    for (k=0;k<pData->seqNum;k++) 
                            wData->partsInWaveIndices[HPWaveNo][wData->partsInWave[HPWaveNo]][k] = HParts->indexes[i][k];
                    wData->partsInWave[HPWaveNo]++;  		
                    wData->partsTotal ++;


                    getWavesParts (pData, wData, HParts->indexes[i], waveNo);
                } /*End Check if partition is not added before*/
            } /* End Search Space check*/
        } /* End if Valid Partitions*/
        /* Check right side of the middle*/
        i = nextParts-p-1;
        if (validPartitions[i] == 1) { /*check valid partitions*/
            /*retrieve partition order in the wave to define its processor*/
            HPWaveNo = getWaveNo (wData->seqNum, wData->partitionSize, HParts->indexes[i]);

            if (HPWaveNo == -1) 
                continue;
            if ((*waveNo) < HPWaveNo) {
                wData->partsInWave = realloc (wData->partsInWave, ((MOATypeInd) (HPWaveNo+1)) * ((MOATypeInd) sizeof *wData->partsInWave));
                if (wData->partsInWave == NULL) {
                    mprintf(1, "Couldn't re-allocate memory for Number of Partitions per wave vector to Calculate waves.\n", 3);
                    printf("Couldn't re-allocate memory for Number of Partitions per wave vector to Calculate waves.\n");		
                    return -1;
                }
                wData->partsInWaveIndices  = realloc (wData->partsInWaveIndices, ((MOATypeInd) (HPWaveNo+1)) * ((MOATypeInd) sizeof *wData->partsInWaveIndices));
                if (wData->partsInWaveIndices == NULL) {
                    mprintf(1, "Couldn't re-allocate memory for partitions Indices matrix to Calculate waves.\n", 3);
                    printf("Couldn't re-allocate memory for partitions Indices matrix to Calculate waves.\n");		
                    return -1;
                }
                wData->waveLength = realloc (wData->waveLength, ((MOATypeInd) (HPWaveNo+1)) * ((MOATypeInd) sizeof *wData->waveLength));
                wData->waveMiddle  = realloc (wData->waveMiddle, ((MOATypeInd) (HPWaveNo+1)) * ((MOATypeInd) sizeof *wData->waveMiddle));
                if (wData->waveMiddle == NULL) {
                    mprintf(1, "Couldn't re-allocate memory for waveMiddle.\n", 3);
                    printf("Couldn't re-allocate memory for waveMiddle.\n");		
                    return -1;
                }                

                for (k=(*waveNo)+1;k<=HPWaveNo;k++) {
                    wData->waveMiddle[k] = mmalloc (((MOATypeInd) pData->seqNum) * ((MOATypeInd) sizeof *wData->waveMiddle[k]));
                    if (wData->waveMiddle[k] == NULL) {
                        mprintf(1, "Couldn't re-allocate memory for waveMiddle Indices.\n", 3);
                        printf("Couldn't re-allocate memory for waveMiddle Indices.\n");		
                        return -1;
                    }                
                    getMiddlePartitionInWave (wData->seqLen, wData->seqNum, wData->partitionSize, k, &wData->waveMiddle[k]) ;
                    wData->waveLength[k] = getMoAWaveLength(wData, k);
                    wData->partsInWave[k] = 0;
                    wData->partsInWaveIndices[k] = NULL;
                }
                (*waveNo) = HPWaveNo;
            }
            //if ((Epsilons == 0) || ((Epsilons != 0) && (distanceFromMiddle(wData, HParts->indexes[i], HPWaveNo) < Epsilons))) {  
            if ((Epsilons == 0) || ((Epsilons != 0) && (isPartWithinEpsilon(HParts->indexes[i], wData, HPWaveNo) == 0))) {  
                if (isPartitionInWave (wData, HPWaveNo, HParts->indexes[i]) == 0) {
                    if (wData->partsInWaveIndices[HPWaveNo] == NULL) 
                        wData->partsInWaveIndices[HPWaveNo] = mmalloc (((MOATypeInd) (wData->partsInWave[HPWaveNo]+1)) * ((MOATypeInd) sizeof *(wData->partsInWaveIndices[HPWaveNo])));
                    else
                        wData->partsInWaveIndices[HPWaveNo] = realloc (wData->partsInWaveIndices[HPWaveNo], ((MOATypeInd) (wData->partsInWave[HPWaveNo]+1)) * ((MOATypeInd) sizeof *(wData->partsInWaveIndices[HPWaveNo])));
                    if (wData->partsInWaveIndices[HPWaveNo] == NULL) {
                        mprintf(1, "Couldn't re-allocate memory for partitions Indices matrix to Calculate waves 2.\n", 3);
                        printf("Couldn't re-allocate memory for partitions Indices matrix to Calculate waves 2.\n");		
                        return -1;
                    }
                    wData->partsInWaveIndices[HPWaveNo][wData->partsInWave[HPWaveNo]] = NULL;
                    wData->partsInWaveIndices[HPWaveNo][wData->partsInWave[HPWaveNo]] = mmalloc (((MOATypeInd) pData->seqNum) * ((MOATypeInd) sizeof *(wData->partsInWaveIndices[HPWaveNo][wData->partsInWave[HPWaveNo]])));
                    if (wData->partsInWaveIndices[HPWaveNo][wData->partsInWave[HPWaveNo]] == NULL) {
                        mprintf(1, "Couldn't re-allocate memory for partitions Indices matrix to Calculate waves 3.\n", 3);
                        printf("Couldn't re-allocate memory for partitions Indices matrix to Calculate waves 3.\n");		
                        return -1;
                    }
                    for (k=0;k<pData->seqNum;k++) 
                            wData->partsInWaveIndices[HPWaveNo][wData->partsInWave[HPWaveNo]][k] = HParts->indexes[i][k];
                    wData->partsInWave[HPWaveNo]++;  		
                    wData->partsTotal ++;
                    getWavesParts (pData, wData, HParts->indexes[i], waveNo);
                } /*End Check if partition is not added before*/
            } /* End Search Space check*/
        } /* End if Valid Partitions*/
    } /* End For Loop for higher partitions*/   
    if (HParts != NULL)
        deleteMOA (HParts);  
    HParts = NULL;
    if (validPartitions != NULL)
        free (validPartitions);
    validPartitions = NULL;
    return 0;
}

long calcWaves (ProcessData * pData, WavesData * wData) {
    long i, j, waveNo;
    MOATypeDimn k; /*Multidimensional Index Iterator*/
    
    int which_proc;
    char msg[SHORT_MESSAGE_SIZE];
    int dbglevel = 0;
	
    /* free space of partitions in wave in case it was calculated before*/
    if (wData->partsInWave != NULL) {

        if (wData->partsInWaveIndices != NULL) {
            for (i=0;i<wData->wavesTotal;i++) {
                if (wData->partsInWaveIndices[i] != NULL) {
                    for (j=0;j<wData->partsInWave[i];j++) {
                        if (wData->partsInWaveIndices[i][j] != NULL)
                            free (wData->partsInWaveIndices[i][j]);
                        wData->partsInWaveIndices[i][j] = NULL;
                    }
                    free (wData->partsInWaveIndices[i]);
                    wData->partsInWaveIndices[i] = NULL;
                }
            }
        }
        free (wData->partsInWaveIndices);
        wData->partsInWaveIndices = NULL;

        free (wData->partsInWave);
        wData->partsInWave = NULL;
    }
    
    /* compute total partitions and duplicates */
    getWavesPartsDuplicates(pData->seqNum, pData->seqLen, wData->partitionSize, &wData->wavesTotal, &wData->partsTotal, &wData->duplicatesTotal);

    /* print debug information (debug level 1) =================*/

    printf ("Dimn\tEpsilon\tPart Size\tLength\t\t\t\tWaveNo\tParts Total\tDuplicates\tComp Parts\n");
    printf ("K\tE\tS\tL\t\t\t\tT\tP\tD\tCP\n");
    printf ("%lld\t%ld\t%ld\t{%lld", pData->seqNum, Epsilons, wData->partitionSize, pData->seqLen[0]);
    for (k=1;k<pData->seqNum;k++) 
        printf (", %lld", pData->seqLen[k]);
    printf ("\t\t\t\t%ld\t%ld\t%lld\t", wData->wavesTotal, wData->partsTotal, wData->duplicatesTotal);

#ifndef NDEBUG
    sprintf (msg, ">>>>For part Size %ld, Complete Waves Total is %ld\n>>>>Partitions Total is %ld\n>>>>Duplicates Total is %lld\n", wData->partitionSize, wData->wavesTotal, wData->partsTotal, wData->duplicatesTotal);
    mprintf (dbglevel, msg, 1);
#endif
    /*First Partition at wave 0 is always the origin : Index Vector is all zeros, and done by processor 0*/
    wData->calcParts = waveNo = 0;
    wData->partsInWave = NULL;
    wData->partsInWave = mmalloc ((MOATypeInd) sizeof *wData->partsInWave);
    if (wData->partsInWave == NULL) {
        mprintf(1, "Couldn't create memory for Number of Partitions per wave vector to Calculate waves.\n", 3);
        printf("Couldn't create memory for Number of Partitions per wave vector to Calculate waves.\n");		
        return -1;
    }
    wData->partsInWave[0] = 1;
    wData->partsInWaveIndices = NULL;
    wData->partsInWaveIndices = mmalloc ((MOATypeInd) sizeof *wData->partsInWaveIndices);

    if (wData->partsInWaveIndices == NULL) {
        mprintf(1, "Couldn't create memory for partitions Indices matrix to Calculate waves.\n", 3);
        printf("Couldn't create memory for partitions Indices matrix to Calculate waves.\n");		
        return -1;
    }
    wData->partsInWaveIndices[0] = NULL;
    wData->partsInWaveIndices[0] = mmalloc ((MOATypeInd) sizeof *(wData->partsInWaveIndices[0]));
    if (wData->partsInWaveIndices[0] == NULL) {
        mprintf(1, "Couldn't re-allocate memory for partitions Indices matrix to Calculate waves 2.\n", 3);
        printf("Couldn't re-allocate memory for partitions Indices matrix to Calculate waves 2.\n");		
        return -1;
    }
    wData->partsInWaveIndices[0][0] = NULL;
    wData->partsInWaveIndices[0][0] = mmalloc (((MOATypeInd) pData->seqNum) * ((MOATypeInd) sizeof *(wData->partsInWaveIndices[0][0])));
    if (wData->partsInWaveIndices[0][0] == NULL) {
        mprintf(1, "Couldn't re-allocate memory for partitions Indices matrix to Calculate waves 3.\n", 3);
        printf("Couldn't re-allocate memory for partitions Indices matrix to Calculate waves 3.\n");		
        return -1;
    }
    wData->waveLength = mmalloc ((MOATypeInd) sizeof *wData->waveLength);
    wData->waveMiddle  = mmalloc ((MOATypeInd) sizeof *wData->waveMiddle);
    if (wData->waveMiddle == NULL) {
        mprintf(1, "Couldn't re-allocate memory for waveMiddle.\n", 3);
        printf("Couldn't re-allocate memory for waveMiddle.\n");		
        return -1;
    }
    wData->waveMiddle[0] = mmalloc (((MOATypeInd) pData->seqNum) * ((MOATypeInd) sizeof *wData->waveMiddle[0]));
    for (k=0;k<pData->seqNum;k++) {
        wData->waveMiddle[0][k] = 0;
        wData->partsInWaveIndices[0][0][k] = 0;
    }
    wData->waveLength[0] = 0;
    pData->partNo = 0;
    pData->waveNo = 0;
    pData->partitionsCount = wData->partsTotal = 1;
    getWavesParts (pData, wData, wData->partsInWaveIndices[0][0], &waveNo);
    wData->wavesTotal = waveNo+1;
    /*Create & initialize the memory for the OCout & OCin for all waves, and delete in freeProcessMemory*/
    pData->OCout = mmalloc(((MOATypeInd) wData->wavesTotal) * ((MOATypeInd) sizeof *(pData->OCout)));
    if (pData->OCout == NULL) {
        mprintf(1, "Couldn't create memory for OCout. Exiting.\n", 3);
        printf("Couldn't create memory for OCout. Exiting.\n");		
        return -1;
    }
    pData->OCin = mmalloc(((MOATypeInd) wData->wavesTotal) * ((MOATypeInd) sizeof *(pData->OCin)));
    if (pData->OCin == NULL) {
        mprintf(1, "Couldn't create memory for OCin. Exiting.\n", 3);
        printf("Couldn't create memory for OCin. Exiting.\n");		
        return -1;
    }
    /*Initialize the OCout per waves*/
    for (waveNo=0;waveNo<wData->wavesTotal;waveNo++) {
        pData->OCout[waveNo].wavesOC = 0;
        pData->OCout[waveNo].WOCO = NULL;
        pData->OCin[waveNo].wavesOC = 0;
        pData->OCin[waveNo].WOCI = NULL;
    }
    printf ("%ld\n", wData->partsTotal);
    printf ("Calculated %ld parts in %ld waves.\n[%d]read %ld local partitions in wave %ld / %ld, starting at part order %ld index {%lld", wData->partsTotal, wData->wavesTotal, myProcid, pData->partitionsCount, pData->waveNo, wData->wavesTotal, pData->partNo, wData->partsInWaveIndices[pData->waveNo][pData->partNo][0]);
    for (k=1;k<pData->seqNum;k++) 
        printf (", %lld", wData->partsInWaveIndices[pData->waveNo][pData->partNo][k]);
    printf ("}\n");
    fflush(stdout);
    /* debug output =========================================================*/
    sprintf (msg, ">>>>Total %ld Parts in  %ld Waves\n", wData->calcParts, waveNo);
    //if (myProcid == 0) 	
    //	printf (">>>>Total %ld Parts in  %ld Waves\nParts in Wave \n", calcParts, waveNo);
    mprintf (0, msg, 1);
    i = 0;
    pData->partitionsCount = 0;
    printf ("Wave No\t\tCalculated\t\tActual\n");
    pData->partitionsCount = 0;
    for (waveNo = 0;waveNo<wData->wavesTotal;waveNo++) {
        if (waveNo >= ceill((double) wData->wavesTotal/2))
            wData->calcParts = getWavePartsTotal ((wData->wavesTotal - waveNo), pData->seqNum);
        else
            wData->calcParts = getWavePartsTotal ((waveNo+1), pData->seqNum);        
        printf ("%ld\t\t%ld\t\t%ld\n", waveNo, wData->calcParts, wData->partsInWave[waveNo]);
        for (j=0;j<wData->partsInWave[waveNo];j++) {
            which_proc = getProcID (wData, waveNo, j);
            if (which_proc == myProcid)
                pData->partitionsCount  ++;
        }
    }
#ifndef NDEBUG
    for (waveNo = 0;waveNo<wData->wavesTotal;waveNo++) {			
        if (waveNo >= ceill((double) wData->wavesTotal/2))
            wData->calcParts = getWavePartsTotal ((wData->wavesTotal - waveNo), pData->seqNum);
        else
            wData->calcParts = getWavePartsTotal ((waveNo+1), pData->seqNum);        
        printf ("Wave %ld has %ld parts Calculated %ld\nSerial\tWave\tPart No\tProc\tPart Index\t Part Distance\n", waveNo, wData->partsInWave[waveNo], wData->calcParts);
       sprintf (msg, "Wave %ld has %ld parts\nSerial\tWave\tPart No\tProc\tPart Index\t Part Angle\n", waveNo, wData->partsInWave[waveNo]);
        mprintf (0, msg, 1);
        for (j=0;j<wData->partsInWave[waveNo];j++) {
            i++;
            //which_proc = getPartitionProcID (wData->partsInWaveIndices[waveNo][j], pData->seqLen, pData->seqNum, waveNo, pData->partitionSize);
            which_proc = getProcID (wData, waveNo, j);
            if (which_proc == myProcid)
                pData->partitionsCount  ++;
            printf ("%ld\t%ld\t%ld\t%d\t{%lld", i, waveNo, j, which_proc, wData->partsInWaveIndices[waveNo][j][0]);
            sprintf (msg, "%ld\t%ld\t%ld\t%d\t{%lld", i, waveNo, j, which_proc, wData->partsInWaveIndices[waveNo][j][0]);
            mprintf (0, msg, 1);
            for (k=1;k<pData->seqNum;k++) {
                printf (", %lld", wData->partsInWaveIndices[waveNo][j][k]);
                sprintf (msg, ", %lld", wData->partsInWaveIndices[waveNo][j][k]);
                mprintf (0, msg, 1);
            }
            double partDistance = (double)  distanceFromMiddle(wData,wData->partsInWaveIndices[waveNo][j], waveNo);
            float ps = (float) partDistance / (float) wData->waveLength[waveNo] * 100;
            printf ("}\t%f\t%f\t%f\n", partDistance, ps, (float) wData->waveLength[waveNo] );
            sprintf (msg, "}\t%f\t%f\t%f\n", partDistance, ps, (float) wData->waveLength[waveNo] );
            mprintf (0, msg, 1);
        }
    }
#endif
    /* ================================ end debug output ==============*/

    mprintf(dbglevel,"out of calcwaves \n", 1);
    printf ("[%d] For %ld parts Distributed Estimated Scoring Total Memory %lld\n", myProcid, pData->partitionsCount, ((MOATypeInd) pData->partitionsCount) * (((MOATypeInd) sizeof (pData->seqNum)) + ((MOATypeInd) pData->seqNum * (MOATypeInd) sizeof *pData->seqLen) + ((MOATypeInd) sizeof (MOATypeInd)) + ((MOATypeInd) sizeof (MOATypeInd) * ((MOATypeInd) powl(wData->partitionSize, pData->seqNum))) + ((MOATypeInd) sizeof (MOA_elm) * ( (MOATypeInd) powl(wData->partitionSize, pData->seqNum)))));

    return wData->wavesTotal;
}
/* ======================================================================
	function calcWaves_ip by Integer Partitions:
	Input:
		dimn: dimension of array
		shape: array of lengths of each dimension.
		partSize: size of partition.
	Output:
		myPartsCount: number of partitions.
		myCurrWave: current wave to be processed.
		myCurrPart: current partition to be processed.
======================================================================= */
long calcWaves_ip (ProcessData * pData, WavesData * wData) {
    long i, j, waveNo;
    MOATypeDimn k; /*Multidimensional Index Iterator*/
    MOATypeInd elements_ub;
    int firstPartFound = 0, which_proc;
    int more = 1;
#ifndef NDEBUG
    char msg[SHORT_MESSAGE_SIZE];
    int dbglevel = 1;
#endif

    /* free space of partitions in wave in case it was calculated before*/
    if (wData->partsInWave != NULL) {

        if (wData->partsInWaveIndices != NULL) {
            for (i=0;i<wData->wavesTotal;i++) {
                if (wData->partsInWaveIndices[i] != NULL) {
                    for (j=0;j<wData->partsInWave[i];j++) {
                        if (wData->partsInWaveIndices[i][j] != NULL)
                            free (wData->partsInWaveIndices[i][j]);
                        wData->partsInWaveIndices[i][j] = NULL;
                    }
                    free (wData->partsInWaveIndices[i]);
                    wData->partsInWaveIndices[i] = NULL;
                }
            }
        }
        free (wData->partsInWaveIndices);
        wData->partsInWaveIndices = NULL;

        free (wData->partsInWave);
        wData->partsInWave = NULL;
    }
    
    /* compute total partitions and duplicates */
    getWavesPartsDuplicates(pData->seqNum, pData->seqLen, wData->partitionSize, &wData->wavesTotal, &wData->partsTotal, &wData->duplicatesTotal);

    /* print debug information (debug level 1) =================*/
    printf ("[%d]For part Size %ld, Estimated wavesTotal %ld, partsTotal %ld, duplicates Total %lld\n", myProcid, wData->partitionSize, wData->wavesTotal, wData->partsTotal, wData->duplicatesTotal);

#ifndef NDEBUG
    sprintf (msg, ">>>>Complete Waves Total is %ld\n>>>>Partitions Total is %ld\n>>>>Duplicates Total is %lld\n", wData->wavesTotal, wData->partsTotal, wData->duplicatesTotal);
    mprintf (dbglevel, msg, 1);
#endif
    wData->calcParts = waveNo = 0;
#ifndef NDEBUG
    /*Temporary Printing Variables*/
    wData->AllpartsInWave = mmalloc ((MOATypeInd) sizeof *wData->AllpartsInWave);
    wData->AllpartsInWaveIndices  = mmalloc ((MOATypeInd) sizeof *wData->AllpartsInWaveIndices);
#endif
    wData->partsInWave = mmalloc ((MOATypeInd) sizeof *wData->partsInWave);
    wData->partsInWaveIndices = mmalloc ((MOATypeInd) sizeof *wData->partsInWaveIndices);
    pData->partitionsCount = 0;
    while (waveNo < wData->wavesTotal) {		
#ifndef NDEBUG
        sprintf (msg, "in loop at part %ld at wave %ld\n", wData->calcParts, waveNo);
        mprintf (1, msg, 1);
        sprintf(msg, "calcWaves: loop parts [%ld] of [%ld]: My Total Partitions: %ld\n", wData->calcParts, wData->partsTotal, pData->partitionsCount);	
        mprintf(dbglevel+1, msg, 1);
#endif
#ifndef NDEBUG
        /*Temporary Printing Variables*/
        wData->AllPartOrder = 0;
        wData->AllpartsInWave = realloc (wData->AllpartsInWave, ((MOATypeInd) (waveNo+1)) * ((MOATypeInd) sizeof *wData->AllpartsInWave));
        wData->AllpartsInWaveIndices  = realloc (wData->AllpartsInWaveIndices, ((MOATypeInd) (waveNo+1)) * ((MOATypeInd) sizeof *wData->AllpartsInWaveIndices));
        wData->AllpartsInWave[waveNo] = 0;
#endif
        wData->partsInWave = realloc (wData->partsInWave, ((MOATypeInd) (waveNo+1)) * ((MOATypeInd) sizeof *wData->partsInWave));
        wData->partsInWaveIndices  = realloc (wData->partsInWaveIndices, ((MOATypeInd) (waveNo+1)) * ((MOATypeInd) sizeof *wData->partsInWaveIndices));
        if (waveNo > ceill((double) wData->wavesTotal/2))
            wData->partsInWave[waveNo] = getWavePartsTotal ((wData->wavesTotal - waveNo), pData->seqNum);
        else
            wData->partsInWave[waveNo] = getWavePartsTotal ((waveNo+1), pData->seqNum);
            
#ifndef NDEBUG
        printf ("calcParts[%ld]=%ld\n", waveNo, wData->partsInWave[waveNo]);
#endif
        wData->calcParts = wData->partsInWave[waveNo];
        wData->partsInWave[waveNo] = getPIndicesinWave(pData, wData, &more, waveNo);
#ifndef NDEBUG
        printf ("actualParts[%ld]=%ld/%ld\n", waveNo, wData->partsInWave[waveNo], wData->calcParts);
#endif
        if (wData->partsInWave[waveNo] == 0)
            break;
        for (i=0;i<wData->partsInWave[waveNo];i++) {
            which_proc = getProcID (wData, waveNo, i);
#ifndef NDEBUG
            sprintf(msg, "[%d]>calcWaves: wave [%ld]- partition [%ld] - Process %d - more %d\n", myProcid, waveNo, i, which_proc, more);	
            mprintf(dbglevel+2, msg, 1);
#endif
            if (which_proc == myProcid) {
                pData->partitionsCount  ++;
                if (firstPartFound == 0)  {
                    pData->partNo = i;
                    pData->waveNo = waveNo;
                    firstPartFound = 1;
                }
            }
        }
		
        waveNo++;
#ifndef NDEBUG
        printf ("Calculating Waves: partitions %ld/%ld for wave %ld/%ld\n", wData->calcParts, wData->partsTotal, waveNo, wData->wavesTotal);
        fflush(stdout);
#endif
    }

    wData->wavesTotal = waveNo;
    wData->partsTotal = wData->calcParts;
    /*Create & initialize the memory for the OCout & OCin for all waves, and delete in freeProcessMemory*/
    pData->OCout = mmalloc(((MOATypeInd) wData->wavesTotal) * ((MOATypeInd) sizeof *(pData->OCout)));
    if (pData->OCout == NULL) {
        mprintf(1, "Couldn't create memory for OCout. Exiting.\n", 3);
        printf("Couldn't create memory for OCout. Exiting.\n");		
        return -1;
    }
    pData->OCin = mmalloc(((MOATypeInd) wData->wavesTotal) * ((MOATypeInd) sizeof *(pData->OCin)));
    if (pData->OCin == NULL) {
        mprintf(1, "Couldn't create memory for OCin. Exiting.\n", 3);
        printf("Couldn't create memory for OCin. Exiting.\n");		
        return -1;
    }
    /*Initialize the OCout per waves*/
    for (waveNo=0;waveNo<wData->wavesTotal;waveNo++) {
        pData->OCout[waveNo].wavesOC = 0;
        pData->OCout[waveNo].WOCO = NULL;
        pData->OCin[waveNo].wavesOC = 0;
        pData->OCin[waveNo].WOCI = NULL;
    }
#ifndef NDEBUG
    printf ("[%d] read %ld local partitions in wave %ld / %ld, starting at part order %ld index {%lld", myProcid, pData->partitionsCount, pData->waveNo, wData->wavesTotal, pData->partNo, wData->partsInWaveIndices[pData->waveNo][pData->partNo][0]);
    for (k=1;k<pData->seqNum;k++) 
        printf (", %lld", wData->partsInWaveIndices[pData->waveNo][pData->partNo][k]);
    printf ("}\n");
    /* debug output =========================================================*/
    sprintf (msg, ">>>>Total %ld Parts in  %ld Waves\n", wData->calcParts, waveNo);
    //if (myProcid == 0) 	
    //	printf (">>>>Total %ld Parts in  %ld Waves\nParts in Wave \n", calcParts, waveNo);
    mprintf (0, msg, 1);
    i = 0;
    for (waveNo = 0;waveNo<wData->wavesTotal;waveNo++) {			
        printf ("Wave %ld has %ld parts\nSerial\tWave\tPart No\tProc\tPart Index\n", waveNo, wData->partsInWave[waveNo]);
        sprintf (msg, "Wave %ld has %ld parts\nSerial\tWave\tPart No\tProc\tPart Index\n", waveNo, wData->partsInWave[waveNo]);
        mprintf (dbglevel+1, msg, 1);
        for (j=0;j<wData->partsInWave[waveNo];j++) {
            i++;
            which_proc = getProcID (wData, waveNo, j);
            printf ("%ld\t%ld\t%ld\t%d\t{%lld", i, waveNo, j, which_proc, wData->partsInWaveIndices[waveNo][j][0]);
            sprintf (msg, "%ld\t%ld\t%ld\t%d\t{%lld", i, waveNo, j, which_proc, wData->partsInWaveIndices[waveNo][j][0]);
            mprintf (dbglevel+1, msg, 1);
            for (k=1;k<pData->seqNum;k++) {
                printf (", %lld", wData->partsInWaveIndices[waveNo][j][k]);
                sprintf (msg, ", %lld", wData->partsInWaveIndices[waveNo][j][k]);
                mprintf (dbglevel+1, msg, 1);
            }
            printf ("}\n");
            mprintf (dbglevel+1, "}\n", 1);
        }
    }
    /* ================================ end debug output ==============*/

    mprintf(dbglevel,"out of calcwaves \n", 1);
    printf ("[%d] For %ld parts Distributed Estimated Scoring Total Memory %lld\n", myProcid, pData->partitionsCount, ((MOATypeInd) pData->partitionsCount) * (((MOATypeInd) sizeof (pData->seqNum)) + ((MOATypeInd) pData->seqNum * (MOATypeInd) sizeof *pData->seqLen) + ((MOATypeInd) sizeof (MOATypeInd)) + ((MOATypeInd) sizeof (MOATypeInd) * ((MOATypeInd) powl(wData->partitionSize, pData->seqNum))) + ((MOATypeInd) sizeof (MOA_elm) * ( (MOATypeInd) powl(wData->partitionSize, pData->seqNum)))));
#endif

    return wData->wavesTotal;
}
/* ***************************************************************
*  function name:  getWavesPartsNumbers
*  Description:    Calculated Total Number of Partitions & duplicates
*  Input:
*		dimn: dimension of array
*		shape: array of lengths of each dimension.
*		pSize: size of partition.
*  Output:
*		wTotal: Total number of waves.
*		pTotal: Total number of partitions.
*		dTotal: Total number of duplicates
 ***************************************************************/
int getWavesPartsDuplicates(MOATypeDimn dimn, MOATypeShape * shape, long pSize, long * wTotal, long * pTotal, MOATypeInd * dTotal) {
	MOATypeDimn i;
	MOATypeInd PrevOC,  * dupInDim = NULL;
	long * partInDim = NULL;
	ldiv_t wavesTotalDiv;

	partInDim = mmalloc (((MOATypeInd) dimn) * ((MOATypeInd) sizeof *partInDim));
	if (partInDim == NULL)
		return -1;
	dupInDim = mmalloc (((MOATypeInd) dimn) * ((MOATypeInd) sizeof *dupInDim));
	if (dupInDim == NULL)
		return -1;
	(*pTotal) = (long) ceill((double) (shape[0] - 1) / (double) (pSize - 1));
	PrevOC = (*pTotal) - 1;
        dupInDim[0]= (MOATypeInd) (*pTotal) - 1;
	(*dTotal) = dupInDim[0];
	(*wTotal) = (long) ceill((double) (shape[0] - 1) / (double) (pSize - 1));
	for (i = 1; i < dimn; i++) {
                (*wTotal) +=  (long) ceill((double) (shape[i] - 1) /(double)  (pSize - 1)) - 1;
		/* to accomodate for variable partition size, this equation need to be in a loop where the staing partition size is the same from both ends of the shape of the current dimension, then increase till the optimal towards the middle of the shape of the dimension*/
		partInDim[i] = (long) ceill((double) (shape[i] - 1) /(double) (pSize - 1));
		(*pTotal) *= partInDim[i];
		
                dupInDim[i] = dupInDim[i-1] * ((MOATypeInd) shape[i]) + PrevOC * ((MOATypeInd) powl(2, i-1));
		(*dTotal) += dupInDim[i];		
		PrevOC *= partInDim[i] - 1;
	}     
	if (partInDim != NULL)
		free (partInDim);
	if (dupInDim != NULL)
		free (dupInDim);
	return 0;
}
int getWavesPartsDuplicates_old(MOATypeDimn dimn, MOATypeShape * shape, long pSize, long * wTotal, long * pTotal, MOATypeInd * dTotal) {
	MOATypeDimn i;
	MOATypeInd PrevOC,  * dupInDim = NULL;
	long * partInDim = NULL;
	ldiv_t wavesTotalDiv;

	partInDim = mmalloc (((MOATypeInd) dimn) * ((MOATypeInd) sizeof *partInDim));
	if (partInDim == NULL)
		return -1;
	dupInDim = mmalloc (((MOATypeInd) dimn) * ((MOATypeInd) sizeof *dupInDim));
	if (dupInDim == NULL)
		return -1;
	(*pTotal) = 1;
	(*dTotal) = 0;
	(*wTotal) = 0;
	PrevOC = 1;
	for (i = 0; i < dimn; i++) {
                (*wTotal) +=  shape[i] - pSize;
		/* to accomodate for variable partition size, this equation need to be in a loop where the staing partition size is the same from both ends of the shape of the current dimension, then increase till the optimal towards the middle of the shape of the dimension*/
		partInDim[i] = (long) ceill((double) (shape[i] - 1) / (pSize - 1));
		(*pTotal) *= partInDim[i];
		
		if (i==0)
			dupInDim[0]= (MOATypeInd) partInDim[0] - 1;
		else
			dupInDim[i] = dupInDim[i-1] * ((MOATypeInd) shape[i]) + PrevOC * ((MOATypeInd) powl(2, i-1));
		(*dTotal) += dupInDim[i];		
		PrevOC *= partInDim[i] - 1;
	}     
	//(*wTotal) = (long) ceill((double) (*wTotal)/(pSize-1)) + dimn;
        wavesTotalDiv = ldiv((*wTotal), pSize-1);
	(*wTotal) = (long) wavesTotalDiv.rem  + wavesTotalDiv.quot + 1;
	if (partInDim != NULL)
		free (partInDim);
	if (dupInDim != NULL)
		free (dupInDim);
	return 0;
}

/**************************************************************
	Function: getWavePartsTotal
	Description:
		Compute maximum number of total partions that can be done in a wave (WaveNo) for Dimension (dimn), the specific seqNum and seqLen, might make the real total partitions less.
	Reurns: Number of partition in the specified wave 
**************************************************************/
long getWavePartsTotal (long WaveNo, MOATypeDimn dimn) {
  long i, WaveParts;

	WaveParts = 0;
	if ((WaveNo < 1) || (dimn < 2))
		WaveParts = -1;
	else if (WaveNo == 1)
		WaveParts = 1;
	else if (dimn == 2)
		WaveParts = WaveNo;
	else {
		for (i=0;i<WaveNo;i++) {
			WaveParts += getWavePartsTotal(WaveNo-i, dimn-1);
		}
	}
	return WaveParts;
}

/***************************************************************************
	Function: getPIndicesinWave
	Description:
		Get partition indices in a wave (waveNo)
	Input
		partsInWave: array contains number of partitions in each wave.
		dimn: dimension of MOA array.
		shape: array of lengths of each dimension.
		waveNo: Wave Number to compute its partition indices.
		partSize: Partition Size.
	Output:
		partsInWaveIndices: Array of partition indices of the specied wave, flag of morePartitions if last indexed partition is not reached yet
	Reurn: -1 incase of error.
***************************************************************************/
void swapDist(MOATypeDimn i, MOATypeDimn j, MOATypeShape * * dist) {
    MOATypeShape temp;
    temp = (*dist)[i];
    (*dist)[i] = (*dist)[j];
    (*dist)[j] = temp;
}
int mpermute(MOATypeDimn dimn, MOATypeShape * * dist, int * more, MOATypeDimn * edge, MOATypeDimn * first, MOATypeDimn * last, MOATypeDimn middle) {
    MOATypeDimn i, j;
    int swapped = 0;
    
    if (((*edge) <= middle) && ((*first)<=middle) && ((*last) >= middle) && ((*dist)[(*first)] != (*dist)[(*last)])) {
        swapDist((*first), (*last), dist);
        swapped = 1;
    }
    (*last) --;
    if ((*last) < middle) {
        (*first) ++;
        (*last) = dimn - 1;
    }
    if ((*first) > middle) {
        (*edge) ++;
        (*first) = 0;
        (*last) = dimn - 1;
    }
    if ((*edge) > middle) {
        (*more) = 0;
        swapped = 1;
    }
    else
        (*more) = 1;
    
    return swapped;
}


long getPIndicesinWave (ProcessData * pData, WavesData * wData, int * morePartitions, long waveNo) {
    MOATypeDimn middleDimn, j;
    MOATypeShape * dist, * dist_orig;
    long mypartsInWave, i;
    int more = 0, pmore = 0;
    MOATypeDimn edge, first, last, middle;
    char msg[SHORT_MESSAGE_SIZE];
	
#ifndef NDEBUG
    int dbglevel = 2;
#endif

    mypartsInWave = wData->partsInWave[waveNo];
    if (mypartsInWave<=0)
            return -1;

    dist = mcalloc ((MOATypeInd) pData->seqNum, (MOATypeInd) sizeof *dist);
    dist_orig = mcalloc ((MOATypeInd) pData->seqNum, (MOATypeInd) sizeof *dist_orig);
    i = 0;
    middleDimn = -1;
    wData->partsInWaveIndices[waveNo] = mmalloc ((MOATypeInd) sizeof *(wData->partsInWaveIndices[waveNo]));
#ifndef NDEBUG
    /* Temporary variable for printing all Paritiong indices whether valid or not*/
    wData->AllpartsInWaveIndices[waveNo] = mmalloc ((MOATypeInd) sizeof *wData->AllpartsInWaveIndices[waveNo]);
#endif
    do  {
        getNextPIndex_middle (&more, waveNo, pData->seqNum, &middleDimn, &dist);
        addPartitionIndex (pData, wData, &i, dist, waveNo, morePartitions);
        if (((Epsilons != 0) && (i >= Epsilons))  || (i >= wData->partsInWave[waveNo]))
            more = 0;
        else {
            for (j=0;j<pData->seqNum;j++)
                dist_orig[j] = dist[j];
            /* Ordered Permutations */
            pmore = 1;
            asort (pData->seqNum, dist);
            addPartitionIndex (pData, wData, &i, dist, waveNo, morePartitions);
            if (((Epsilons != 0) && (i >= Epsilons))  || (i >= wData->partsInWave[waveNo])) {
                pmore = 0;
                more = 0;
            }
            while (pmore ==1) {
                permute(pData->seqNum, dist, &pmore);
                if (dist[(MOATypeShape) floorl((double) pData->seqNum/2)] == waveNo)
                    more = 0;
                addPartitionIndex (pData, wData, &i, dist, waveNo, morePartitions);
                if (((Epsilons != 0) && (i >= Epsilons))  || (i >= wData->partsInWave[waveNo])) {
                    pmore = 0;
                    more = 0;
                }
            }
            /* Middle Permutations
            pmore = 1;
            edge = first = 0;
            last = dimn-1;
            middle = (MOATypeDimn) floorl((long double) dimn/2);
            while (pmore ==1) {
                while (mpermute(dimn, &dist, &pmore, &edge, &first, &last, middle) == 0);
                addPartitionIndex (&i, dist, partsInWave, dimn, shape, waveNo, partSize, partsInWaveIndices, mpartsInWaveIndices, morePartitions);
                /*sprintf (msg, "(%lld", dist[0]);
                for (j=1;j<dimn;j++)
                    sprintf (msg, "%s, %lld", msg, dist[j]);
                sprintf (msg, "%s)\n", msg);
                mprintf (0, msg, 1);
                if ((Epsilons != 0) && (i >= Epsilons))
                    pmore = 0;
            }
             */
            for (j=0;j<pData->seqNum;j++)
                dist[j] = dist_orig[j];
        }
    } while (more == 1);
    free (dist);
    free (dist_orig);
#ifndef NDEBUG
    sprintf (msg, ">>>>Actual Partitions in Wave %ld is %ld\n", waveNo, i);
    mprintf (dbglevel, msg, 1);
#endif
    return i;
}

/******************************************************************************
	Function: getNextPIndex
	Description:
		recursive function to compute Partition Indexes of a wave.
	Input/output:
		more: to inicate first call and next calls.
	Input:
		dimn: dimension of MOA array. originally equal act_dimn and recursively decreased by 1.
		act_dimn: actual dimension of MOA array.
		waveNo: Wave Number
	Output:
		PIndex: Array of partition Indices.
******************************************************************************/

void distributeRemDistance (MOATypeShape * * PIndex, MOATypeDimn * middleDimn, MOATypeDimn dimn, MOATypeDimn startDimn, MOATypeShape startDist, MOATypeShape remDist, MOATypeDimn remDimn) {
	MOATypeDimn i;
	for (i=startDimn;i>=0;i--) { 
		(*PIndex)[i] = startDist;
		remDist = remDist-(*PIndex)[i];
		remDimn --;		
		startDist = (MOATypeShape) ceill((long double) remDist/remDimn);
		if (i == 0)
			(*PIndex)[dimn-i-1] = remDist;
		else
			(*PIndex)[dimn-i-1] = startDist;
		remDist = remDist-(*PIndex)[dimn-i-1];
		remDimn --;
		if (remDimn > 0)		 
			startDist = (MOATypeShape) ceill((long double) remDist/remDimn);
		else
			startDist = 0;
		if ((*PIndex)[i] > 0)
			(*middleDimn) = i;
	}

}

void getNextPIndex_middle(int * more, long waveNo, MOATypeDimn dimn, MOATypeDimn * middleDimn, MOATypeShape * * PIndex) {
	MOATypeDimn startDimn, remDimn;
	MOATypeShape startDist, remDist;


	if ((*more) == 0) {  /* first call ===========================*/
#ifndef NDEBUG
		mprintf (10, " 000 ", 1);
#endif
		remDist = waveNo;
		remDimn = dimn;
		startDist = (MOATypeShape) ceill((long double) remDist/remDimn);
		(*middleDimn) = (MOATypeDimn) floorl((long double) dimn/2);
		startDimn = (*middleDimn);
		(*PIndex)[(*middleDimn)] = startDist;
		remDist = remDist-startDist;
		remDimn --;		
		startDist = (MOATypeShape) ceill((long double) remDist/remDimn);
		if (isEven(dimn) == 1) {
			if ((*middleDimn) == 0)
				(*PIndex)[dimn-(*middleDimn)-1] = remDist;
			else
				(*PIndex)[dimn-(*middleDimn)-1] = startDist;
			remDist = remDist-(*PIndex)[dimn-(*middleDimn)-1];
			remDimn --;		
			startDist = (MOATypeShape) ceill((long double) remDist/remDimn);
			if ((*PIndex)[(*middleDimn)] > 1)
				(*middleDimn) --;
			startDimn --;
		}
		if ((*PIndex)[(*middleDimn)] > 1)
			(*middleDimn) --;
		startDimn--;
		distributeRemDistance (PIndex, middleDimn, dimn, startDimn, startDist, remDist, remDimn);
	}
	else {
                if (((*middleDimn) != (dimn-(*middleDimn)-1)) && ((*PIndex)[dimn-(*middleDimn)-1] > 0)) {
                    (*PIndex)[(*middleDimn)] ++;
                    (*PIndex)[dimn-(*middleDimn)-1] --;
                }
                else {
                    (*PIndex)[(*middleDimn)] --;
                    (*PIndex)[(*middleDimn)+1] ++;
                    if ((*PIndex)[(*middleDimn)] == 0)
                        (*middleDimn) ++;
                    remDist = waveNo;
                    remDimn = dimn;
                    for (startDimn = (MOATypeDimn) floorl((long double) dimn/2);startDimn>=0;startDimn--) {
                        if (startDimn >= (*middleDimn)) {
                            remDist -= (*PIndex)[startDimn]+(*PIndex)[dimn-startDimn-1];
                            remDimn -= 2;
                            if ((isEven(dimn) == 1) && (startDimn == (MOATypeDimn) floorl((long double) dimn/2))) 
                                startDimn --;
                       }
                        else 
                            (*PIndex)[startDimn] = (*PIndex)[dimn-startDimn-1] = 0;      
                    }

                    if ((remDimn > 0) && (remDist > 0)) {
                        startDist = (MOATypeShape) ceill((long double) remDist/remDimn);
                        distributeRemDistance (PIndex, middleDimn, dimn, (*middleDimn)-1, startDist, remDist, remDimn);
                    }
                }
	}	
	if ((*PIndex)[(MOATypeDimn) floorl((long double) dimn/2)] == waveNo)  
		(*more) = 0;
	else
		(*more) = 1;

}


/********************************************************************************************
Function  getNextPIndex_ordered
*********************************************************************************************/
void getNextPIndex_ordered(int * more, MOATypeDimn dimn, MOATypeDimn act_dimn, long waveNo, MOATypeShape * * PIndex) {
	MOATypeInd startDist;
	MOATypeDimn i;

	/* Method 2: Change the kth index in each k dimension, to result in asending order index to avoid sorting later
	 get starting Distance from origin: if called to the first time, distribute the distance from origin over dimensions to get the middle point of the current wave (wave Number divided by the dimension), else, it is the last point assigned to the first dimension, to continue from there*/
	if ((*more) == 0) {  /* first call ===========================*/
		if (waveNo > 0)
			startDist = (MOATypeInd) ceill((long double) waveNo/act_dimn);
		else
			startDist = 0;
		/*Now assign the middle point in the wave, to the middle dimension*/
		(*PIndex)[dimn-1] = startDist;
		for (i=dimn-2;i>=0;i--) 
			(*PIndex)[i] = 0;
	}
	else  /* next calls ===========================================*/
		startDist = (*PIndex)[dimn-1];
		
#ifndef NDEBUG
	/*sprintf (msg, "in d %ld, %ld, %ld w %ld d %ld s%ld m %d\n", (*PIndex)[0], (*PIndex)[1], (*PIndex)[2], waveNo, dimn, startDist, (*more));
	mprintf (10, msg, 1);*/
#endif
	/* All Terminal Cases first, then the 2D terminal cases, then the rest N-D cases
	 test for wave 0 for all dimensions is the same, all zeros*/
	if (waveNo == 0) {
#ifndef NDEBUG
		mprintf (10, " 000 ", 1);
#endif
		for (i=0;i<act_dimn;i++) 
			(*PIndex)[i] = 0;
		(*more) = 0;
	}
	else if	((*PIndex)[dimn-1] == waveNo) {
	/* test for last case where the starting Distance is already equal the wave number, then the rest of dimensions' indices are zeros*/
#ifndef NDEBUG
		mprintf (10, " 111 ", 1);
#endif
		for (i=dimn-2;i>=0;i--) 
			(*PIndex)[i] = 0;
		(*more) = 0;
	}
	else if (( (*PIndex)[dimn-1] * dimn == waveNo) && ((*PIndex)[dimn-1] != (*PIndex)[dimn-2])) {
	/* All Equal Indices, starting middle point	  */
#ifndef NDEBUG
		mprintf (10, " 222 ", 1);
#endif
		for (i=dimn-2;i>=0;i--) 
			(*PIndex)[i] = (*PIndex)[dimn-1];
		(*more) = 1;
	}
	else if (dimn == 2) { 
	/* Cases for 2D*/
		if (((*PIndex)[dimn-2] != waveNo - startDist)) {
		/* Firs Case, where first dimension = starting middle Distance, and second dimension is not yet assigned,*/
#ifndef NDEBUG
			mprintf (10, " 333 ", 1);
#endif
			(*PIndex)[dimn-2] = waveNo - startDist;
			(*more) = 1;
		}
		else if (startDist+1<waveNo) {
		/* Increase one and decrease the other*/
#ifndef NDEBUG
			mprintf (10, " 444 ", 1);
#endif
			startDist++;
			(*PIndex)[dimn-1] = startDist;
			(*PIndex)[dimn-2] = waveNo - startDist;
			(*more) = 1;
		} 
		else if (startDist+1==waveNo) {
		/* no more and return, we should never reach this case, since it is covered in case 111*/
#ifndef NDEBUG
			mprintf (10, " 555 ", 1);
#endif
			(*PIndex)[dimn-1] = waveNo;
			(*PIndex)[dimn-2] = 0;
			(*more) = 0;
 		}
	}
	else {
	/* more the 2D*/
		if  ((( (*PIndex)[dimn-1] * dimn == waveNo) && ((*PIndex)[dimn-1] == (*PIndex)[dimn-2])) || 
		(((*PIndex)[dimn-1] + (*PIndex)[dimn-2] == waveNo) && ((*PIndex)[dimn-1]+1<waveNo)))  {
		/* next case after all are equal
		 First and second dimension are equal waveNo, i.e. all permutations are finished for the previous first dimension,  then increase the first dimension,*/
#ifndef NDEBUG
			mprintf (10, " 666  ", 1);
#endif
			startDist++;
			(*PIndex)[dimn-1] = startDist;
			for (i=dimn-2;i>=0;i--) 
				(*PIndex)[i] = 0;
			if (waveNo - startDist != 0) {
				(*more) = 0;
				getNextPIndex_ordered (more, dimn-1, act_dimn, waveNo - startDist, PIndex);
				(*more) = 1;
			}
			else
				(*more) = 0;
		}
		else if ((*PIndex)[dimn-1] + (*PIndex)[dimn-2] != waveNo) {
		/* First and second dimension are not yet equal waveNo, i.e. there more permutations for the current first dimension value*/
#ifndef NDEBUG
			mprintf (10, " 777 ", 1);
#endif
			(*PIndex)[dimn-1] = startDist;
			getNextPIndex_ordered (more, dimn-1, act_dimn, waveNo - startDist, PIndex);
			(*more) = 1;
		}	
		else if ((*PIndex)[dimn-1]+1==waveNo) {
		/* last case, where the final increase will be the waveNo, and the rest need to be zeros			mprintf (10, " 888 ", 1);*/
			(*PIndex)[dimn-1] = waveNo;
			for (i=dimn-2;i>=0;i--) 
				(*PIndex)[i] = 0;
			(*more) = 0;
		} 
	}
#ifndef NDEBUG
	/*sprintf (msg, "out d %ld, %ld, %ld w %ld d %ld s%ld m %d\n", (*PIndex)[0], (*PIndex)[1], (*PIndex)[2], waveNo, dimn, startDist, (*more));
	mprintf (10, msg, 1);*/
#endif
}
/*********************************************************************************
	Function: addPartitionIndex
	Description:
	Input:
		dimn: dimension of MOA array.
		waveNo: Wave Number
		partSize: Partition Size
	Output:
	
*********************************************************************************/
int isPartitionInWave (WavesData * wData, long waveNo, MOATypeShape * partIndex) {
    MOATypeDimn k;
    int ret = 0;
    long i;
    if (waveNo < wData->wavesTotal) {
        for (i=0;i<wData->partsInWave[waveNo] && ret == 0;i++) {
            ret = 1;
            for (k=0;k<wData->seqNum && ret == 1;k++) {
                if (partIndex[k] != wData->partsInWaveIndices[waveNo][i][k])
                    ret = 0;
            }
            
        }
    }
    return ret;
}
int addPartitionIndex (ProcessData * pData, WavesData * wData, long * PartOrder, MOATypeShape * dist, long waveNo, int * morePartitions) {
    int localMorePartitions = 0;
    int ret, validIndex;
    MOATypeDimn k;
    MOATypeShape * ind;

    ind = mcalloc ((MOATypeInd) pData->seqNum, (MOATypeInd) sizeof *ind);
    validIndex = 1;
    ret = 0;
    for (k=0;k<pData->seqNum;k++) {
        ind[k] = dist[k] * (wData->partitionSize-1);
        /* It is invalid if it is more than the shape (lengths of sequences), or if it at the last one, this means that it is only the overlapping cells, with nothing internal to compute on its own, which would be redundancy and wouldn't happen actually*/
        if ((ind[k] < 0) || (ind[k] >= (pData->seqLen[k]-1)))
            validIndex = 0;

        /*If a reduced search space is limited with Epsilons value, then don;t take more partitions in this wave*/
        if ((Epsilons != 0) && ((*PartOrder) >= Epsilons)) 
            validIndex = 0;
        /* if we add the partition size to the index, we can know that value of
       the last cell in this parition, and check if it is the last cell in the 
       whole tensor or more, to decide that it is the last partition or there is still more*/
        if (((ind[k] + (wData->partitionSize-1)) < (pData->seqLen[k]-1)) && ((*morePartitions) == 1))
            localMorePartitions = 1;
    }
    /*if (*morePartitions) is already zero, it should never be changed locally*/
    if ((*morePartitions) == 1)
        (*morePartitions) = localMorePartitions;
#ifndef NDEBUG
    /*Temporary Printing Variables*/	
    if (notPreviouslyVisited_dbg(wData, waveNo, ind) == 0) { 
        if (myProcid == 0) {
            printf ("%ld/%ld\t%ld/%ld\t{%lld", waveNo, wData->wavesTotal, wData->calcParts, wData->partsTotal, ind[0]);	
            for (k=1;k<pData->seqNum;k++) 
                printf (", %lld", ind[k]);

            printf ("}\t\t%d\t\t%ld\n", validIndex, (*PartOrder));	
        }
        wData->AllpartsInWave[waveNo] ++;
        wData->AllpartsInWaveIndices[waveNo] = realloc (wData->AllpartsInWaveIndices[waveNo], ((MOATypeInd) (wData->AllPartOrder+1)) * ((MOATypeInd) sizeof *wData->AllpartsInWaveIndices[waveNo]));
        wData->AllpartsInWaveIndices[waveNo][wData->AllPartOrder] = mmalloc (((MOATypeInd) pData->seqNum) * ((MOATypeInd) sizeof *(wData->AllpartsInWaveIndices[waveNo][(*PartOrder)])));
        for (k=0;k<pData->seqNum;k++) 
            wData->AllpartsInWaveIndices[waveNo][wData->AllPartOrder][k] = ind[k];
        wData->AllPartOrder++;  		
        wData->AllpartsInWave[waveNo] = wData->AllPartOrder;
    }
#endif
    if ((notPreviouslyVisited(wData, (*PartOrder), waveNo, ind) == 0)  && (validIndex == 1) ) {
            wData->partsInWaveIndices[waveNo] = realloc (wData->partsInWaveIndices[waveNo], ((MOATypeInd) ((*PartOrder)+1)) * ((MOATypeInd) sizeof *(wData->partsInWaveIndices[waveNo])));
            wData->partsInWaveIndices[waveNo][(*PartOrder)] = mmalloc (((MOATypeInd) pData->seqNum) * ((MOATypeInd) sizeof *(wData->partsInWaveIndices[waveNo][(*PartOrder)])));
            for (k=0;k<pData->seqNum;k++) 
                    wData->partsInWaveIndices[waveNo][(*PartOrder)][k] = ind[k];
            (*PartOrder)++;  		
            ret = 1;
    }
    free (ind);
    return ret;
}

/********************************************************************
	Function: notPreviouslyVisited_atAll - for printing of all partition indices whether valid or not
********************************************************************/
int notPreviouslyVisited_dbg (WavesData * wData, long waveNo, MOATypeShape * ind) {
    long i, j, startPart;
     MOATypeDimn k;
   int found = 0;
    /* Check if it is not a valid Index as well*/
    if (ind == NULL)
        return  1;
    for (i=0;(i<=waveNo) && (found == 0); i++) {
        /* Check if this is the last wave being computed, take the current total parts found so far in the argument, otherwise, use the saved in the array*/
        if (i==waveNo)
            startPart = wData->AllPartOrder;
        else
            startPart = wData->AllpartsInWave[i];
            for (j=0;(j<startPart) && (found == 0); j++) {
                found = 1;
                for (k=0;(k<wData->seqNum) && (found == 1);k++) 
                    if (wData->AllpartsInWaveIndices[i][j][k] != ind[k]) {
                        found = 0;
                        break;
                    }
            }
    }
    return found;
}

/********************************************************************
	Function: notPreviouslyVisited
********************************************************************/
int notPreviouslyVisited (WavesData * wData, long partsinCurrentWave, long waveNo, MOATypeShape * ind) {
    long i, j, startPart;
    MOATypeDimn k;
    int found = 0;
    /* Check if it is not a valid Index as well*/
    if (ind == NULL)
        return  1;
    //for (i=0;(i<=waveNo) && (found == 0); i++) {
    i=waveNo;
    if (i==waveNo)
        startPart = partsinCurrentWave;
    else
        startPart = wData->partsInWave[i];
    for (j=0;(j<startPart) && (found == 0); j++) {
        found = 1;
        for (k=0;(k<wData->seqNum) && (found ==1);k++) 
            if (wData->partsInWaveIndices[i][j][k] != ind[k]) {
                found = 0;
                break;
            }
    }
    //}
    return found;
}

/*********************************************************************
	Function: freeProcessMemory
		Free all allocated data to Process Data
*********************************************************************/
void freeProcessMemory (ProcessData * * pData, ScoringData * * sData, WavesData * * wData) {
    MOATypeInd i, j;
    MOATypeDimn k;

    /* Process Data ====================================*/
    if ((*pData) != NULL) {
        /* Incomming Ovelaping Cells (Received Cells) ==================  */
        if ((*pData)->OCin != NULL) {
            for (i=0;i<(*wData)->wavesTotal;i++) {
                if (((*pData)->OCin[i].wavesOC > 0) && ((*pData)->OCin[i].WOCI != NULL)) {
                    for (j=0;j<(*pData)->OCin[i].wavesOC;j++) {
                        if ((*pData)->OCin[i].WOCI[j].cellIndex != NULL) {
                            free ((*pData)->OCin[i].WOCI[j].cellIndex);
                            (*pData)->OCin[i].WOCI[j].cellIndex = NULL;
                        }
                    }
                    free ((*pData)->OCin[i].WOCI);
                    (*pData)->OCin[i].WOCI = NULL;
                }
            }
            free((*pData)->OCin);
            (*pData)->OCin = NULL;
        }

        /* Outgoing Ovelapping Cells (Sent Cells) ====================*/
        if ((*pData)->OCout  != NULL) {
            for (i=0;i<(*wData)->wavesTotal;i++) {
                if (((*pData)->OCout[i].wavesOC > 0) && ((*pData)->OCout[i].WOCO != NULL)) {
                   for (j=0;j<(*pData)->OCout[i].wavesOC;j++) {
                        if ((*pData)->OCout[i].WOCO[j].cellIndex != NULL) {
                            free ((*pData)->OCout[i].WOCO[j].cellIndex);
                            (*pData)->OCout[i].WOCO[j].cellIndex = NULL;
                        }
                        if (((*pData)->OCout[i].WOCO[j].depProc_ub > 0) && ((*pData)->OCout[i].WOCO[j].depProc  != NULL)) {
                            free((*pData)->OCout[i].WOCO[j].depProc);	
                            (*pData)->OCout[i].WOCO[j].depProc = NULL;
                        }
                    }
                    free ((*pData)->OCout[i].WOCO);
                    (*pData)->OCout[i].WOCO = NULL;
                }
            }
            free((*pData)->OCout);
            (*pData)->OCout = NULL;
        }
        if ((*pData)->validPartitions != NULL)
            free((*pData)->validPartitions);
        (*pData)->validPartitions = NULL;
        if ((*pData)->msaAlgn != NULL) {
            /* Delete MOA ==================================== */
            deleteMOA ((*pData)->msaAlgn);
        }
        /* sequences =======================================*/
        if (((*pData)->seqNum > 0) && ((*pData)->seqLen != NULL)) {
            if ((*pData)->sequences != NULL) {
                for (i=0;i<(*pData)->seqNum;i++) {
                    if (((*pData)->seqLen[i] > 0) && ((*pData)->sequences[i] != NULL)) {
                        free ((*pData)->sequences[i]);
                        (*pData)->sequences[i] = NULL;
                        if ((*pData)->seqName != NULL) {
                            if ((*pData)->seqName[i] != NULL) {
                                free ((*pData)->seqName[i]);
                                (*pData)->seqName[i] = NULL;
                            }
                        }
                    }
                }
                free ((*pData)->sequences);
                (*pData)->sequences = NULL;
            }
            free ((*pData)->seqLen);
            (*pData)->seqLen = NULL;
        }
        if ((*pData)->seqName != NULL) 
            free ((*pData)->seqName);
        (*pData)->seqName = NULL;

        (*pData)->seqLen = NULL;
        (*pData)->seqNum = 0;
        (*pData)->partitionsCount = 0;
        (*pData)->computedPartitions = 0;
        free ((*pData));
        (*pData) = NULL;
    }

    if ((*wData) != NULL) {
        if ((*wData)->partsInWaveIndices != NULL) {
            for (i=0;i<(*wData)->wavesTotal;i++) {
                if ((*wData)->waveMiddle != NULL) {
                    if ((*wData)->waveMiddle[i] != NULL) 
                        free ((*wData)->waveMiddle[i]);
                    (*wData)->waveMiddle[i] = NULL;
                }
                if ((*wData)->partsInWaveIndices[i] != NULL) {
                    for (j=0;j<(*wData)->partsInWave[i];j++) {
                        free ((*wData)->partsInWaveIndices[i][j]);
                        (*wData)->partsInWaveIndices[i][j] = NULL;
                    }
                    free ((*wData)->partsInWaveIndices[i]);
                    (*wData)->partsInWaveIndices[i] = NULL;
                }
            }
            free ((*wData)->partsInWaveIndices); 
            (*wData)->partsInWaveIndices = NULL;
        }
        if ((*wData)->partsInWave != NULL)
            free((*wData)->partsInWave);
        (*wData)->partsInWave = NULL;
        if ((*wData)->waveMiddle != NULL) {
            free((*wData)->waveMiddle);
        }
        (*wData)->waveMiddle = NULL;
        if ((*wData)->waveLength != NULL) {
            free((*wData)->waveLength);
        }
        (*wData)->waveLength = NULL;
        free ((*wData));
        (*wData) = NULL;
    }


    /* free allocated work memory ==================================*/
    if ((*sData) != NULL) {
        if ((*sData)->p_index != NULL)
            free((*sData)->p_index);
        (*sData)->p_index = NULL;
        if ((*sData)->gm_index != NULL)
            free((*sData)->gm_index);
        (*sData)->gm_index = NULL;
        if ((*sData)->lm_index != NULL)
            free((*sData)->lm_index);
        (*sData)->lm_index = NULL;
        if ((*sData)->depPart != NULL)
            free((*sData)->depPart);
        (*sData)->depPart = NULL;
        if ((*sData)->neighbor != NULL)
            free((*sData)->neighbor);
        (*sData)->neighbor = NULL;
        /*if ((*sData)->lnScores != NULL)
            free((*sData)->lnScores);
        (*sData)->lnScores = NULL;
        if ((*sData)->lnIndices != NULL){
            for (i=0; i<(*sData)->CalLnCount; i++)  {
                if ((*sData)->lnIndices[i] != NULL) {
                    free ((*sData)->lnIndices[i]);
                    (*sData)->lnIndices[i] = NULL;
                }
            }
            free((*sData)->lnIndices);
            (*sData)->lnIndices = NULL;
        }
        if ((*sData)->lnInSearchSpace != NULL)
            free((*sData)->lnInSearchSpace);
        (*sData)->lnInSearchSpace = NULL;*/

        if ((*sData)->decremented != NULL)
            free((*sData)->decremented);
        (*sData)->decremented = NULL;

        if ((*sData)->NghbMOA != NULL)
            deleteMOA ((*sData)->NghbMOA);

        if ((*sData)->posDimn != NULL)
            free ((*sData)->posDimn);
        (*sData)->posDimn = NULL;
        if ((*sData)->pwScores != NULL) {
            for (i=0; i<(*sData)->seqNum; i++) {
                if ((*sData)->pwScores[i] != NULL) {
                    free((*sData)->pwScores[i]);
                    (*sData)->pwScores[i] = NULL;
                }
            }
            free((*sData)->pwScores);
            (*sData)->pwScores = NULL;
        }
        free((*sData));
        (*sData) = NULL;
    }
}
/***************************************************************************
	Function: isPartInSearchSpace
	Description:
		returns true if partition index is included in the reduced search space, false otherwise
***************************************************************************/
int isPartInSearchSpace(MOATypeShape * partIndex, WavesData * wData) {
    int found = -1;
    long i, startingWave;
    MOATypeDimn k;
    
    if (partIndex == NULL) {
        printf ("You need to provide the multidimensional partition index to test if it is included in search space. Exiting\n");
        return -1;
    }    
    long waveNo = getWaveNo (wData->seqNum, wData->partitionSize, partIndex);
    if (waveNo > wData->wavesTotal-1)
        return -1;

    for (i=0;(i<wData->partsInWave[waveNo]) && (found == -1);i++) { 
        found = 0; /* if it remains zero, then partition found*/
        for (k=0;k<wData->seqNum;k++) {
            if (partIndex[k] != wData->partsInWaveIndices[waveNo][i][k]) {
                found = -1;
                break;
            }
        }
    }
    return found;
}
/***************************************************************************
	Function: isCellInSearchSpace
	Description:
		returns a partition No (j)  and waveNo, is cellIndex is in partition that is included in the reduced search space, -1 otherwise
***************************************************************************/
long isCellInSearchSpace(MOATypeShape * * cellIndex, WavesData * wData, long * waveNo) {
    long partNo;
    MOATypeInd i;
    MOATypeDimn k;
    MOATypeShape temp, * partIndex = NULL;
    MOA_rec * part;
    
    if (cellIndex == NULL) {
        printf ("You need to provide the multidimensional partition index to test if it is included in search space. Exiting\n");
        return -1;
    }    
    partIndex = mmalloc ((MOATypeInd) wData->seqNum * ((MOATypeInd) sizeof *partIndex));
    if (getNeighbors (2, (*cellIndex), wData->seqNum, wData->seqLen, &part) == 0) {
        for (i=0;i< part->elements_ub;i++) {
            if  (getPartitionIndex (part->indexes[i], wData->seqNum, wData->seqLen, wData->partitionSize, &partIndex) == 0) {
                if (isPartInSearchSpace(partIndex, wData) == 0) {
                    for (k=0;k< wData->seqNum;k++) 
                        (*cellIndex)[k] = part->indexes[i][k];
                    getPartitionPosition (wData, partIndex, waveNo, &partNo);
                    if (partIndex != NULL)
                        free (partIndex);
                    partIndex = NULL;
                    deleteMOA(part);
                    return partNo;
                }
            }
        }
    }
    if (partIndex != NULL)
        free (partIndex);
    partIndex = NULL;
    (*waveNo) = -1;
    return -1;
}
/***************************************************************************
	Function: getPositionalPartitionIndex
	Description:

Get the Partition Index of a specific cell, based on the given indicators
Input:
Flag: 0 == return the Partition where the cell is a lower border cell in all dimensions, which is the partition index itself, if there is one, otherwise return -1.
1 == return the Partition where the cell is a lower border cell on the specific dimensions specified in argument posDimn == 1,
2 == return the Partition where the cell is a higher border cell in all dimensions, which is the last cell in the parition
3 == return the Partition where the cell is a higher border cell on the specific dimension specified in argument posDimn  == 1,
4 == return the Partition where the cell is a perfect internal cell

Should pass either the flat index (cellIndex), or the m_index, otherwise initialise, with -ve value for first, and null for second.
dimn and shape to manipulate the index
the partition size to get the partition size

Output:
Will return the partition multidimensional index in the pointer to the partIndex array, and the flat index in the return long value.

***************************************************************************/
int  getPositionalPartitionIndex (int flag, MOATypeShape * posDimn, MOATypeShape * cellIndex, MOATypeDimn dimn, MOATypeShape * shape, long partitionSize, MOATypeShape * * partIndex) {
    MOATypeDimn i;
    int  valid = 0;

    if (cellIndex == NULL) {
        printf ("You need to provide either the flat Index, or the multidimensional index of a cell, to calculate its partition Index. Exiting\n");
        return -1;
    }
    for (i=0;(i<dimn) && (valid == 0);i++) {
        switch (flag) {
            case 0 :
                (*partIndex)[i] = (cellIndex[i] + 1) % (partitionSize - 1);
                (*partIndex)[i] = (cellIndex[i] + 1) - (*partIndex)[i];
                break;
            case 1 :
                if (posDimn[i] == 1) {
                    (*partIndex)[i] = (cellIndex[i] + 1) % (partitionSize - 1);
                    (*partIndex)[i] = (cellIndex[i] + 1) - (*partIndex)[i];
                }
                else {
                    (*partIndex)[i] = cellIndex[i] % (partitionSize - 1);
                    (*partIndex)[i] = cellIndex[i] - (*partIndex)[i];
                }
                break;
            case 2 :
                (*partIndex)[i] = (cellIndex[i] - 1) % (partitionSize - 1);
                (*partIndex)[i] = (cellIndex[i] - 1) - (*partIndex)[i];
                break;
            case 3 :
                if (posDimn[i] == 1) {
                    (*partIndex)[i] = (cellIndex[i] - 1) % (partitionSize - 1);
                    (*partIndex)[i] = (cellIndex[i] - 1) - (*partIndex)[i];
                }
                else {
                    (*partIndex)[i] = cellIndex[i] % (partitionSize - 1);
                    (*partIndex)[i] = cellIndex[i] - (*partIndex)[i];
                }
                break;
            case 4 :
                break;
            default  :
                valid = -1;
                break;
        }
        if ((*partIndex)[i] < 0)
            (*partIndex)[i] = 0;
        if ((*partIndex)[i] >= shape[i]-1)
            (*partIndex)[i] = shape[i] - partitionSize;
    }
    return valid;
}

/*********************************************************************************
 *  Function : getPartitionIndex
 *  takes input of hyper-tensor description (dimn and shape), and a cell index, and return the 
 *  minimum partition Index that encapsulate the cell index. based on the partition size.
*********************************************************************************/
int  getPartitionIndex (MOATypeShape * cellIndex, MOATypeDimn dimn, MOATypeShape * shape, long partitionSize, MOATypeShape * * partIndex) {
    MOATypeDimn i;

    if (cellIndex == NULL) {
        printf ("You need to provide either the multidimensional index of a cell, to calculate its partition Index. Exiting\n");
        return -1;
    }
    for (i=0;i<dimn;i++) {
        (*partIndex)[i] = cellIndex[i] % ((MOATypeShape) (partitionSize - 1));
        if ((*partIndex)[i] == 0)
            (*partIndex)[i] = cellIndex[i] - partitionSize + 1;
        else
            (*partIndex)[i] = cellIndex[i] - (*partIndex)[i];
        if ((*partIndex)[i] < 0)
            (*partIndex)[i] = 0;
        if ((*partIndex)[i] >= shape[i]-1)
            (*partIndex)[i] = shape[i] - partitionSize;
    }
    return 0;
}
/*********************************************************************************
 *  Function : getWaveNo 
*********************************************************************************/

long getWaveNo (MOATypeDimn seqNum, long partitionSize, MOATypeShape * partIndex) {
    long  waveNo = -1;
    MOATypeDimn k;
    MOATypeShape * distance = NULL;

    distance = mmalloc ((MOATypeInd) seqNum * ((MOATypeInd) sizeof *distance));
    if (distance == NULL) {
        printf ("You need to provide either the flat Index, or the multidimensional distance of a partition Index. Exiting\n");
        return -1;
    }
    for (k=0;k<seqNum;k++) 
        distance[k] = partIndex[k] / ((MOATypeShape)(partitionSize-1));
    waveNo = Sum (distance, seqNum);
    if (distance != NULL)
        free (distance);
    return waveNo;
}

/*********************************************************************************
 *  Function : getPartNo 
*********************************************************************************/
long getPartNo (WavesData * wData, MOATypeShape * partIndex, long waveNo) {
    long j, partNo = -1;
    MOATypeDimn k;
    int found = -1;

    if (partIndex == NULL) {
        printf ("You need to provide the multidimensional partition index to calculate its wave Number and partition Number. Exiting\n");
        return -1;
    }
    for (j=0;(j<wData->partsInWave[waveNo]) && (found == -1);j++) { 
        found = 0;
        for (k=0;k<wData->seqNum;k++) {
            if (wData->partsInWaveIndices[waveNo][j][k] != partIndex[k]) {
                found = -1;
                break;
            }
        } 
        if (found == 0)  {
            partNo = j;
            break;
        }
    }
    
    return partNo;
}

/*********************************************************************************
 *  Function : getPartitionPosition 
*********************************************************************************/
long getPartitionPosition (WavesData * wData, MOATypeShape * partIndex, long * waveNo, long * partNo) {
    long j;
    MOATypeDimn k;

    if (partIndex == NULL) {
        printf ("You need to provide either the multidimensional partition index to calculate its wave Number and partition Number. Exiting\n");
        return -1;
    }
    (*waveNo) = getWaveNo (wData->seqNum, wData->partitionSize, partIndex);
    (*partNo) = getPartNo (wData, partIndex, (*waveNo));
    return (*partNo);
}


/*********************************************************************************
	Function: getPartitionDetails
		identify minimum partition index from a global cell index, its waveNo and partNo and its processor
 *              ( hence whether it is local or not)
	Input:
		pData: process data to search all local partitions
		globalIndex: globalIndex to be found in local partitions
	Output:
		partIndex: as a return value, NULL if not found
		localIndex: local index in the local partition found, corresponding to the globalIndex, NULL if not local
 *              procID: will be equal myProcid if local, otherwise if remote.
 *              0 in return if local partition, -1 if remote or invalid.
*********************************************************************************/
int getPartitionDetails (ProcessData * pData, WavesData * wData, MOATypeShape * * partIndex, MOATypeShape * cellIndex, long * waveNo, long * partNo, int * procID) {
    (*procID) = (*waveNo) = (*partNo) = -1;
    if  (getPartitionIndex (cellIndex, pData->seqNum, pData->seqLen, wData->partitionSize, partIndex) == 0) {
        getPartitionPosition (wData, (*partIndex), waveNo, partNo);
        if ((*partNo) >= 0) {
            (*procID) = getProcID (wData, (*waveNo), (*partNo));
            return 0;
        }
        return -1;
    }
    return -1;
}
MOATypeInd getLocalFlatIndex (MOATypeShape * globalIndex, MOA_rec * msaAlgn, MOATypeShape * shape, long partSize) {
    MOATypeInd localFlatIndex;
    MOATypeDimn k;
    int found;
    MOATypeShape * localIndex = NULL;
    localIndex = malloc (((MOATypeInd) msaAlgn->dimn) * ((MOATypeInd) sizeof *localIndex));
    if (localIndex == NULL) {
        printf ("Error creating memory for temporary index. Exiting\n");
        return -1;
    }
    getLocalIndex (globalIndex, msaAlgn->indexes[0], msaAlgn->dimn, shape, partSize, &localIndex);
    localFlatIndex = Gamma(localIndex, msaAlgn->dimn, msaAlgn->shape, msaAlgn->dimn, 1);
    if (localIndex != NULL) 
        free (localIndex);
    localIndex = NULL;
    if ((localFlatIndex >= 0) && (localFlatIndex < msaAlgn->elements_ub))
        return localFlatIndex;
    return -1;
}

/*********************************************************************************
	Function: getPartition
		Construct Partition MOA record
	Input:
		flatIndex:
		dimn: dimension of MOA array.
		shape: array of lengths of each dimention.
		waveNo: Wave Number to compute its partition indices.
		partSize: Partition Size.
	Output:
		part: MOA Partition record.
*********************************************************************************/
int getPartition (MOATypeShape * startingIndex, MOATypeDimn dimn, MOATypeShape * shape, MOA_rec * * part, long partSize) {
    return getHigherNeighbors (partSize, startingIndex, dimn, shape, part);
}

/***************************************************************************
	Function: getNextPartition - used in scoring
	Input/Output:
		waveNo: Identify Current Wave and contains next partition wave
		PartOrder: Identify current partition and returns with next partition
***************************************************************************/
void getNextPartition (WavesData * wData, long * waveNo, long * partNo) {
    ldiv_t partsInCluster;
    long break_point, currWave, currPart;
    int valid = 1;
    currWave = (*waveNo);
    currPart = (*partNo);
    
    if (((*partNo) == wData->partsInWave[(*waveNo)] - 1) && ((*waveNo) == wData->wavesTotal - 1)) {
        (*partNo) = -1;
        (*waveNo) = -1;
        return ;
    }
    (*partNo) ++;
    if ((*partNo) > wData->partsInWave[(*waveNo)]) {
        valid = 0;
    }
    else if(getProcID(wData, (*waveNo), (*partNo)) != myProcid) {
        valid = 0;
    }

    while (((*waveNo) < wData->wavesTotal - 1)  && (valid == 0)) {
        (*waveNo) ++;  /* Need to debug this to fix for initial case of restore*/
        partsInCluster = ldiv(wData->partsInWave[(*waveNo)], ClusterSize);
        if ((wData->partsInWave[(*waveNo)] < ClusterSize) && (partsInCluster.quot == 0) && (myProcid >= partsInCluster.rem)) {
            (*partNo) = -1;
            //valid = 0;
        }
        else {
            valid = 1;
            if (partsInCluster.rem > 0) {
                if (myProcid < partsInCluster.rem) {
                    (*partNo) = (partsInCluster.quot + 1) * myProcid;
                } else {
                    break_point = partsInCluster.rem * (partsInCluster.quot + 1);
                    (*partNo) = partsInCluster.quot * (myProcid - partsInCluster.rem) + break_point;
                }
            } else {
                (*partNo) = partsInCluster.quot * myProcid;
            }
        }
        /*if (((*partNo) < currPart) && ((*waveNo) == currWave) || (valid == 0)) {
            (*partNo) = -1;
             valid = 0;
            (*waveNo) ++;
        }*/
    }	
}
/***************************************************************************
	Function: getPrevPartition - used in trace back
	Input/Output:
		waveNo: Identify Current Wave and contains previous partition in wave
		PartNo: Identify current partition and returns with previous partition
***************************************************************************/
void getPrevPartition (WavesData * wData, long * waveNo, long * partNo) {
    ldiv_t partsInCluster;
    long break_point;
    int valid = 0;
    
    if (((*partNo) == 0) && ((*waveNo) == 0)) {
        (*partNo) = -1;
        (*waveNo) = -1;
        return ;
    }
    if (((*partNo) > 0)  && (valid == 0)) {
        (*partNo) --;
        if (getProcID(wData, (*waveNo), (*partNo)) == myProcid) 
            valid = 1;
    }
    while (((*waveNo) > 0)  && (valid == 0)) {        
        (*waveNo) --;
        partsInCluster = ldiv(wData->partsInWave[(*waveNo)], ClusterSize);
        if ((wData->partsInWave[(*waveNo)] < ClusterSize) && (partsInCluster.quot == 0) && (myProcid >= partsInCluster.rem)) 
            (*partNo) = -1;
        else {
            valid = 1;
            if (partsInCluster.rem > 0) {
                if (myProcid < partsInCluster.rem) {
                    (*partNo) = ((partsInCluster.quot + 1) * myProcid) + partsInCluster.quot;
                } else {
                    break_point = partsInCluster.rem * (partsInCluster.quot + 1);
                    (*partNo) = (partsInCluster.quot * (myProcid - partsInCluster.rem)) + partsInCluster.quot + break_point - 1;
                }
            } else {
                (*partNo) = (partsInCluster.quot * myProcid) + partsInCluster.quot - 1;
            }
            if ((*partNo) >= wData->partsInWave[(*waveNo)]) {
                (*partNo) = -1;
                valid = 0;
            }
        }
    }	
}
/************************************************************
	Function: getProcID
	Description:
		Retutns Process ID (task) that process Wave (waveNo) 
		and Partition (PartOrder)
***********************************************************/
int getProcID (WavesData * wData, long waveNo, long partNo) {
	ldiv_t partsInCluster, procValue;
	long break_point;
	int procID = -1;
        if ((waveNo < 0) || (partNo < 0))
            return -1;
	/*printf (" w %ld o %ld ap %ld ", waveNo, PartOrder, partsInWave[waveNo]);*/
	/*****************************************************************
	partsInWave[waveNo]: total partitions in wave [wp]
	PartOrder: order (sequence) of partition in wave
	ClusterSize: totals tasks (processes) running
	*****************************************************************/
	if (wData->partsInWave != NULL) {
		if (partNo == 0)
			procID = 0;
		else if (partNo < wData->partsInWave[waveNo]) {
			partsInCluster = ldiv(wData->partsInWave[waveNo], ClusterSize);
			if (partsInCluster.rem > 0) {
				/* set break point */
				break_point = partsInCluster.rem * (partsInCluster.quot + 1);
				if (partNo < break_point) {
					procValue = ldiv(partNo, (partsInCluster.quot + 1));
					procID = procValue.quot;
				} else {
					procValue = ldiv((partNo-break_point), partsInCluster.quot);
					procID = procValue.quot + partsInCluster.rem;
				}
			} else {
				procValue = ldiv(partNo, partsInCluster.quot);
				procID = procValue.quot;
			}
		}
	}
	return procID;
}

/**********************************************************************************************
 ************   Given a CellIndex, PartitionsSize and a PartIndex, returns 1 if  **************
 ************   the cellindex falls inside the partition, and zero otherwise ******************
 **********************************************************************************************/

int IsCellInPart (MOATypeShape * cellIndex, MOATypeShape * partIndex, MOATypeDimn dimn, MOATypeShape * shape, long partSize) {
    MOATypeDimn i;
    MOATypeShape remShape;
    for (i=0;i<dimn;i++) {
        if (partIndex[i] + partSize < shape[i])
            remShape = partIndex[i] + partSize - 1;
        else
            remShape = shape[i]-1;
        if ((cellIndex[i] < partIndex[i]) || (cellIndex[i] > remShape)) {
            return -1;
        }
    }
    return 0;
}

/**********************************************************************************************
 ************   Given a global CellIndex, PartitionsSize and a PartIndex, returns *************
 **************** the local address in a partition, given that IsCellInPart is  ***************
 ***************** called first , otherwise will return -1 ****************************** *****
 **********************************************************************************************/

int getLocalIndex (MOATypeShape * cellIndex, MOATypeShape * partIndex, MOATypeDimn dimn, MOATypeShape * shape, long partSize, MOATypeShape * * localIndex) {
    MOATypeDimn i;
    MOATypeShape remShape;
    for (i=0;i<dimn;i++) {
        if (partIndex[i] + partSize < shape[i])
            remShape = partIndex[i] + partSize - 1;
        else
            remShape = shape[i]-1;
        if ((cellIndex[i] < partIndex[i]) || (cellIndex[i] > remShape)) {
           return -1;
        }
       else
           (*localIndex)[i] = cellIndex[i] - partIndex[i];
    }
    return 0;
}

/*****************************************************************
	Function: getDepProcs
*****************************************************************/
void getDepProcs (ProcessData * pData, ScoringData * sData, WavesData * wData, long OCIndex) {
    MOATypeDimn k;
    long i, partNo, waveNo;
    MOATypeShape * mIndex = NULL, * NghbInd = NULL;
    int pfound, nproc;
    MOA_rec * HParts = NULL;
#ifndef NDEBUG
    char msg[SHORT_MESSAGE_SIZE];
    int dbglevel = 5;
#endif
	
#ifndef NDEBUG
    sprintf (msg, "[%d]>getDepProcs:Wave %ld /%ld: cell {%lld\n", myProcid, pData->waveNo, wData->wavesTotal, pData->OCout[pData->waveNo].WOCO[OCIndex].cellIndex[0]);
    for (k=1;k<pData->seqNum;k++) 
        sprintf(msg, "%s, %lld", msg, pData->OCout[pData->waveNo].WOCO[OCIndex].cellIndex[k]);
    sprintf (msg, "%s}\n", msg);
    mprintf (dbglevel, msg, 2);
#endif
    if (pData->waveNo+1 < wData->wavesTotal)  {  // check if there are remaining higher waves
    /* OC in the current wave, expected to send to the next wave only, if there is any*/
        mIndex = mcalloc ((MOATypeInd) pData->seqNum, ((MOATypeInd) sizeof *mIndex));
        if (mIndex == NULL) {
            printf ("Couldn't create memory for m_index for neighbors in calcDep. Exiting\n");
            return;
        }
        MOAGetHigherPartitions (wData->partitionSize, 2, sData->p_index, pData->seqNum, pData->seqLen, &HParts, &pData->validPartitions); /* get all higher neighbors of the current Cell Index*/
#ifndef NDEBUG
        sprintf (msg, "[%d]>getDepProcs: OC Neighbors %lld ", myProcid, HParts->elements_ub);
        mprintf (dbglevel, msg, 2);
#endif

#ifndef NDEBUG
        sprintf (msg, "Cell {%lld ", pData->OCout[pData->waveNo].WOCO[OCIndex].cellIndex[0]);
        for (k=1;k<pData->seqNum;k++) 
            sprintf (msg, "%s, %lld ", msg, pData->OCout[pData->waveNo].WOCO[OCIndex].cellIndex[k]);
        sprintf(msg, "%s}\n", msg);
        mprintf (dbglevel, msg, 2);
#endif
        for (i=1;i<HParts->elements_ub;i++) { /*Loop through Higher Neighboring partitions*/
            if (pData->validPartitions[i] == 1) { /*check valid partitions*/
                /*retrieve partition order in the wave to define its processor*/
                partNo = -1;
                waveNo = getWaveNo (wData->seqNum, wData->partitionSize, HParts->indexes[i]);
                if (waveNo == -1) /*if the above loop ended and no partition was found, then break to next neighboring partition*/
                    continue;
                partNo = getPartNo (wData, HParts->indexes[i], waveNo);
                if (partNo == -1) /*if the above loop ended and no partition was found, then break to next neighboring partition*/
                    continue;
                nproc = getProcID (wData, waveNo, partNo);
#ifndef NDEBUG
                sprintf (msg, " proc %d ", nproc);
                mprintf (dbglevel, msg, 2);							
#endif
                if (nproc != myProcid) { /*if not a local partition, and need to be sent to another processor*/
#ifndef NDEBUG
                    sprintf (msg, "[%d] w=%ld: sending to proc %d CI {%lld ", myProcid, pData->waveNo, nproc, pData->OCout[pData->waveNo].WOCO[OCIndex].cellIndex[0]);
                    for (k=1;k<pData->seqNum;k++) 
                        sprintf (msg, "%s, %lld", msg, pData->OCout[pData->waveNo].WOCO[OCIndex].cellIndex[k]);
                    sprintf (msg, "%s} \n", msg);
                    mprintf (dbglevel, msg, 2);
#endif
                    /*Check if no processor list created, create one*/
                    if (pData->OCout[pData->waveNo].WOCO[OCIndex].depProc_ub == 0) {
                        pData->OCout[pData->waveNo].WOCO[OCIndex].depProc_ub ++;
                        pData->OCout[pData->waveNo].WOCO[OCIndex].depProc = NULL;
                        pData->OCout[pData->waveNo].WOCO[OCIndex].depProc = mmalloc ((MOATypeInd) sizeof *pData->OCout[pData->waveNo].WOCO[OCIndex].depProc);
                        if (pData->OCout[pData->waveNo].WOCO[OCIndex].depProc == NULL) {
                            printf ("Couldn't allocate memory for dep Proc in OCout. Exiting\n");
                            return;
                        }
                        pData->OCout[pData->waveNo].WOCO[OCIndex].depProc[0] = nproc;
                    }
                    /*else, check if processor already included, or else, re-allocate and append the new processor to the list*/
                    else {
                        pfound = 0;
                        for (k = 0; k < pData->OCout[pData->waveNo].WOCO[OCIndex].depProc_ub;k++) {
                            if (pData->OCout[pData->waveNo].WOCO[OCIndex].depProc[k] == nproc)
                                pfound = 1;
                        }
                        if (pfound == 0) {
                            pData->OCout[pData->waveNo].WOCO[OCIndex].depProc_ub ++;
                            pData->OCout[pData->waveNo].WOCO[OCIndex].depProc = realloc (pData->OCout[pData->waveNo].WOCO[OCIndex].depProc, ((MOATypeInd) pData->OCout[pData->waveNo].WOCO[OCIndex].depProc_ub) * ((MOATypeInd) sizeof *pData->OCout[pData->waveNo].WOCO[OCIndex].depProc));
                            if (pData->OCout[pData->waveNo].WOCO[OCIndex].depProc == NULL) {
                                printf ("Couldn't reallocate memory for deo Proc in OCout. Exiting\n");
                                return;
                            }
                            pData->OCout[pData->waveNo].WOCO[OCIndex].depProc[pData->OCout[pData->waveNo].WOCO[OCIndex].depProc_ub-1] = nproc;
                        }
                    }
                } /*end check non-local partition*/
#ifndef NDEBUG
                sprintf (msg,  " dep_ub %d ", pData->OCout[pData->waveNo].WOCO[OCIndex].depProc_ub);
                mprintf (dbglevel, msg, 2);
#endif
            } /* End check valid Partition*/
        } /* end loop for Higher partitions*/
#ifndef NDEBUG
        mprintf (dbglevel, "\n", 2);			
#endif
    } // end check remaining waves

    if (mIndex != NULL)
        free (mIndex); 
    if (NghbInd != NULL)
        free (NghbInd);
    if (HParts != NULL)
        deleteMOA (HParts);
}



/******************************************************************************
	Function: addOC
		add overlapping cell to Outgoing cells array (OCout)
******************************************************************************/
int addOC (ProcessData * pData, ScoringData * sData, WavesData * wData) {
    MOATypeDimn k;
    long i;
    int found = -1;
#ifndef NDEBUG
    int dbglevel = 4;
    char msg[LONG_MESSAGE_SIZE];
#endif
	
#ifndef NDEBUG
    /*printf ("[%d]>addOC: count[%ld] w[%ld], pi{%lld", myProcid, pData->OCout_ub, pData->waveNo, sData->p_index[0]);
    for (k=1;k<pData->seqNum;k++) 
        printf (", %lld", sData->p_index[k]);
    printf ("}, ci{%lld", sData->gm_index[0]);
    for (k=1;k<pData->seqNum;k++) 
        printf (", %lld", sData->gm_index[k]);
    printf("}, cs[%lld]\n", sData->score);*/
#endif
    for (i=0;(i<pData->OCout[pData->waveNo].wavesOC) && (found == -1);i++) {
        found = 0;
        for (k=0;k<pData->seqNum;k++) {
            if (pData->OCout[pData->waveNo].WOCO[i].cellIndex[k] != sData->gm_index[k]) {
                found = -1;
                break;
            }
        }
        if (found == 0) {
#ifndef NDEBUG
            sprintf (msg, "[%d]>addOC: found %ld ", myProcid, i);
            mprintf (dbglevel, msg, 1);		
#endif
            pData->OCout[pData->waveNo].WOCO[i].cellScore = sData->score;
        }
    }
    
    if (found == -1) {
#ifndef NDEBUG
        sprintf (msg, "[%d]>addOC: Not found added\n", myProcid);
        mprintf (dbglevel, msg, 1);		
#endif
        pData->OCout[pData->waveNo].wavesOC ++;
        if (pData->OCout[pData->waveNo].WOCO == NULL)
            pData->OCout[pData->waveNo].WOCO = mmalloc(((MOATypeInd) pData->OCout[pData->waveNo].wavesOC) * ((MOATypeInd) sizeof *(pData->OCout[pData->waveNo].WOCO)));
        else
            pData->OCout[pData->waveNo].WOCO = realloc(pData->OCout[pData->waveNo].WOCO, ((MOATypeInd) pData->OCout[pData->waveNo].wavesOC) * ((MOATypeInd) sizeof *(pData->OCout[pData->waveNo].WOCO)));
        if (pData->OCout[pData->waveNo].WOCO == NULL) {
            mprintf(1, "Couldn't create memory for OCout[waveNo].WOCO while adding an OC. Exiting.\n", 3);
            printf("Couldn't create memory for OCout[waveNo = %ld].WOCO while adding an OC in wave %ld. Exiting.\n", pData->waveNo,  pData->waveNo);
            return -1;
        }
        pData->OCout[pData->waveNo].WOCO[pData->OCout[pData->waveNo].wavesOC-1].cellIndex = NULL;
        pData->OCout[pData->waveNo].WOCO[pData->OCout[pData->waveNo].wavesOC-1].cellIndex = mmalloc (((MOATypeInd) pData->seqNum) * ((MOATypeInd) sizeof *(pData->OCout[pData->waveNo].WOCO[pData->OCout[pData->waveNo].wavesOC-1].cellIndex)));
        if (pData->OCout[pData->waveNo].WOCO[pData->OCout[pData->waveNo].wavesOC-1].cellIndex == NULL) {
            mprintf(1, "Couldn't create memory for OCout[waveNo].WOCO cellIndex while adding an OC. Exiting.\n", 3);
            printf("Couldn't create memory for OCout[waveNo = %ld].WOCO cellIndex while adding an OC in wave %ld. Exiting.\n", pData->waveNo,  pData->waveNo);
            return -1;
        }

        for (k=0;k<pData->seqNum;k++) {
            pData->OCout[pData->waveNo].WOCO[pData->OCout[pData->waveNo].wavesOC-1].cellIndex[k] = sData->gm_index[k];
        }
        pData->OCout[pData->waveNo].WOCO[pData->OCout[pData->waveNo].wavesOC-1].cellScore = sData->score;
        pData->OCout[pData->waveNo].WOCO[pData->OCout[pData->waveNo].wavesOC-1].depProc_ub = 0;
        pData->OCout[pData->waveNo].WOCO[pData->OCout[pData->waveNo].wavesOC-1].depProc = NULL;
        pData->OCout[pData->waveNo].WOCO[pData->OCout[pData->waveNo].wavesOC-1].sent = 0;
        getDepProcs (pData, sData, wData, pData->OCout[pData->waveNo].wavesOC-1);
    } 
    return 0;
}
#ifndef NDEBUG
void checkMPIErrorCode (char * module, int dbglevel, int caller, int errmsg) {
	char msg[MID_MESSAGE_SIZE];
	
	switch(errmsg) {
		case MPI_SUCCESS:
			//sprintf (msg, "$$$$[%s]: Success Proc[%d]\n", module, myProcid);
			//mprintf(dbglevel, msg, caller);
			break;
		case MPI_ERR_COMM:
			sprintf (msg, "$$$$[%s]: Invalid communicator[%d]\n", module, myProcid);
			mprintf(dbglevel, msg, caller);
			break;
		case MPI_ERR_TYPE:
			sprintf (msg, "$$$$[%s]: Invalid datatype argument[%d]\n", module, myProcid);
			mprintf(dbglevel, msg, caller);
			break;
		case MPI_ERR_COUNT:
			sprintf (msg, "$$$$[%s]: Invalid count argument Proc[%d]\n", module, myProcid);
			mprintf(dbglevel, msg, caller);
			break;
		case MPI_ERR_TAG:
			sprintf (msg, "$$$$[%s]: Invalid  tag  argument Proc[%d]\n", module, myProcid);
			mprintf(dbglevel, msg, caller);
			break;
		case MPI_ERR_RANK:
			sprintf (msg, "$$$$[%s]: Invalid source or destination rank Proc[%d]\n", module, myProcid);
			mprintf(dbglevel, msg, caller);
			break;
	}
}
#endif
/*******************************************************************
	Function: prepareWaveOC
		Prepares Ovelapping Cells in Buffers to dependent processors
*******************************************************************/
int prepareWaveOC (long waveNo, ProcessData * pData) {
    int temp, p, dest, found;
    long i, j, l;
    MOATypeDimn k, bIterator;
    MPI_Status * status;
#ifndef NDEBUG
    int dbglevel = 1;
    char msg[SHORT_MESSAGE_SIZE];
#endif	
    
    /*printf ("[%d] in prepareWaveOC @ wave %ld part {%lld", myProcid, waveNo, pData->msaAlgn->indexes[0][0]);
    for (k=1;k<pData->seqNum;k++) 
        printf (", %lld", pData->msaAlgn->indexes[0][k]);
    printf ("}\n");
    fflush (stdout);*/
    if (pData->sendRequests > 0) {
        if (pData->mpi_requests != NULL) {
            status = mmalloc (pData->sendRequests * sizeof *status);
            MPI_Waitall(pData->sendRequests, pData->mpi_requests, status);
            free (status);
            free(pData->mpi_requests);
        }
    }
    pData->mpi_requests = NULL;
    pData->sendRequests = 0;
    /* free communication buffer allocated memory =====================*/
    if (pData->buffer != NULL) {
        for (p=0;p<pData->depProcCount;p++) {
            if (pData->buffer[p] != NULL)
                free (pData->buffer[p]);
            pData->buffer[p] = NULL;
        }
        free (pData->buffer);
    }
    pData->buffer = NULL;
    if (pData->buffer_size != NULL)
        free ( pData->buffer_size);
    pData->buffer_size = NULL;
    if (pData->proc_index != NULL)
        free (pData->proc_index);
    pData->proc_index = NULL;
    pData->depProcCount = 0;

    pData->proc_index = mmalloc (((MOATypeInd) ClusterSize) * ((MOATypeInd) sizeof *pData->proc_index));
    if (pData->proc_index == NULL) {
        return -1;
    }
    pData->buffer_size = mmalloc (((MOATypeInd) ClusterSize) * ((MOATypeInd) sizeof *pData->buffer_size));
    if (pData->buffer_size == NULL) {
        return -1;
    }
    pData->buffer = mmalloc (((MOATypeInd) ClusterSize) * ((MOATypeInd) sizeof * pData->buffer));
    if ( pData->buffer == NULL) {
        return -1;
    }
    for (p=0;p<ClusterSize;p++) {
        pData->proc_index[p] = -1;
        pData->buffer_size[p] = 0;
    }
    pData->depProcCount = 0;		
    for (i=0; i<pData->OCout[waveNo].wavesOC;i++) {
        if (pData->OCout[waveNo].WOCO[i].sent == 0) {
            if (pData->OCout[waveNo].WOCO[i].depProc_ub > 0) {		
                for (j=0; j<pData->OCout[waveNo].WOCO[i].depProc_ub;j++) {
                    bIterator = 0;
                    found = 0;
                    for (p=0; p<pData->depProcCount; p++) {
                        if (pData->proc_index[p] == pData->OCout[waveNo].WOCO[i].depProc[j]) {
                             pData->buffer_size[p] += pData->commBufSize;
                            pData->buffer[p] = realloc (pData->buffer[p], ((MOATypeInd)  pData->buffer_size[p]) * ((MOATypeInd) sizeof *pData->buffer[p]));
                            found = 1;
                            break;
                        }
                    }                    
                    if (found == 0) {
                        /* add the destination process ================== */
                        if (pData->depProcCount >= ClusterSize) {
                            /* Error exceeds number of processes =====*/
#ifndef NDEBUG
                            sprintf(msg, "[%d]>sendOC: Process Count [%d] Exceeds Cluster Size [%d], exiting!\n", myProcid, pData->depProcCount, ClusterSize);
                            mprintf(dbglevel, msg, 1);
#endif	
                            fflush(stdout);
                            return 0;
                        } else { 
                            pData->buffer[pData->depProcCount] = mmalloc (((MOATypeInd) pData->commBufSize) * ((MOATypeInd) sizeof *pData->buffer[pData->depProcCount]));
                            pData->proc_index[pData->depProcCount] = pData->OCout[waveNo].WOCO[i].depProc[j];
                            pData->buffer_size[pData->depProcCount] = pData->commBufSize;
                            p = pData->depProcCount;
                            pData->depProcCount++;
                        }
                    }
                    /* fill buffer with OCout elements ================ */

                    l =  pData->buffer_size[p] - pData->commBufSize;
                    pData->buffer[p][l++] = (MOATypeElmVal) waveNo; /*wave No*/
                    for (k=0;k<pData->seqNum;k++) {
                        pData->buffer[p][l++] = (MOATypeElmVal) pData->OCout[waveNo].WOCO[i].cellIndex[k];
                    }
                    pData->buffer[p][l++] = (MOATypeElmVal) pData->OCout[waveNo].WOCO[i].cellScore;

#ifndef NDEBUG
                    /*printf ("[%d] Sending cellIndex {%lld ", myProcid, pData->OCout[waveNo].WOCO[i].cellIndex[0]);
                    for (k=1;k<pData->seqNum;k++) 
                        printf (", %lld", pData->OCout[waveNo].WOCO[i].cellIndex[k]);
                    printf ("} score %lld  in wave %ld to proc %d\n", pData->OCout[waveNo].WOCO[i].cellScore, waveNo, pData->proc_index[p]);
                    */
                    sprintf (msg, "[%d]sendOC: i %ld ci  {%lld", myProcid, waveNo, pData->OCout[waveNo].WOCO[i].cellIndex[0]);
                    for (k=1;k<pData->seqNum;k++) 
                        sprintf (msg, "%s, %lld", msg, pData->OCout[waveNo].WOCO[i].cellIndex[k]);
                    sprintf (msg, "%s} cs %lld to proc %d\n", msg, pData->OCout[waveNo].WOCO[i].cellScore, pData->proc_index[p]);
                    mprintf (dbglevel, msg, 2);
                    bIterator = pData->commBufSize;
                    sprintf (msg, "[%d]sendOC: buf (%d,%ld) i %lld ci {%lld", myProcid , p, pData->buffer_size[p], i, pData->buffer[p][l-bIterator]);
                    for (k=1;k<pData->seqNum;k++)  {
                        bIterator --;
                        sprintf (msg, "%s, %lld", msg, pData->buffer[p][l-bIterator]);
                    }
                    bIterator --;
                    sprintf (msg, "%s} cs %lld\n", msg, pData->buffer[p][l-bIterator]);
                    mprintf (dbglevel, msg, 2);
#endif
                } /* end loop for dependent processes */
            }
            pData->OCout[waveNo].WOCO[i].sent = 1;
        } /* end selection of particular OCout (partIndex and not sent)*/
    }
    return 0;
}

int sendOCtoHigherProcessors (ProcessData * pData) {
    int p, MPI_return;
    long l;
    MOATypeDimn k, bIterator;
#ifndef NDEBUG
    int dbglevel = 1;
    char msg[SHORT_MESSAGE_SIZE];
#endif	
    if (pData->depProcCount > 0) {
        pData->mpi_requests = mmalloc(2 * pData->depProcCount * sizeof *pData->mpi_requests);
        /* Send Overalpping Cells =================*/
        for (p=0; p<pData->depProcCount;p++) {
            if (pData->proc_index[p] > myProcid) {
            MPI_return = MPI_Isend(&pData->buffer_size[p], 1, MPI_LONG, pData->proc_index[p], MOAMSA_SEND_RECEIVE_TAG, MOAMSA_COMM_WORLD, &pData->mpi_requests[pData->sendRequests]);
            pData->sendRequests ++;
            MPI_return = MPI_Isend(pData->buffer[p],  pData->buffer_size[p], MPI_LONG_LONG, pData->proc_index[p], MOAMSA_SEND_RECEIVE_TAG, MOAMSA_COMM_WORLD, &pData->mpi_requests[pData->sendRequests]);
            pData->sendRequests ++;
#ifndef NDEBUG
            //printf ("[%d]>sendOC: Element Count[%ld] To Process [%d]\n", myProcid,  pData->buffer_size[p], pData->proc_index[p]);
            //fflush(stdout);
            sprintf (msg, "[%d]>sendOC: Element Count[%ld] To Process [%d]\n", myProcid,  pData->buffer_size[p], pData->proc_index[p]);
            mprintf (dbglevel, msg, 2);
            //printf ("(  wn,   pi,   ci,     cs)\n");
            //fflush(stdout);
            mprintf (dbglevel, "(  wn,   pi,   ci,     cs)\n", 2);
            for (l=0; l< pData->buffer_size[p]; l+=pData->commBufSize) {
                //printf ("(%4lld, {%lld", pData->buffer[p][l], pData->buffer[p][l+1]);
                //fflush(stdout);
                sprintf (msg, "(%4lld, {%lld", pData->buffer[p][l], pData->buffer[p][l+1]);
                mprintf (dbglevel, msg, 2);
                bIterator = 2;
                for (k=1;k<pData->seqNum;k++) {
                    //printf (", %lld", pData->buffer[p][bIterator]);
                    //fflush(stdout);
                    sprintf (msg, ", %lld", pData->buffer[p][l+bIterator]);
                    mprintf (dbglevel, msg, 2);
                    bIterator ++;
                }
                //printf ("}, {%4lld", pData->buffer[p][l+2]);
                //fflush(stdout);
                sprintf (msg, "}, %6lld)\n", pData->buffer[p][l+pData->commBufSize-1]);
                mprintf (dbglevel, msg, 2);
            }
#endif
        } /*End check if Higher Processor*/
        } /*End Loop for Processors*/
#ifndef NDEBUG
        sprintf (msg, "[%d]>sendOC:sent OC to total processes [%d]\n", myProcid,pData->depProcCount);
        mprintf (dbglevel, msg, 2);
#endif
    }
    return 0;
}
int sendOCtoLowerProcessors (ProcessData * pData){
    int p, MPI_return;
    long l;
    MOATypeDimn k, bIterator;
#ifndef NDEBUG
    int dbglevel = 1;
    char msg[SHORT_MESSAGE_SIZE];
#endif	
    if (pData->depProcCount > 0) {
        /* Send Overalpping Cells =================*/
        for (p=0; p<pData->depProcCount;p++) {
            if (pData->proc_index[p] < myProcid) {
            //printf ("[%d] sending to %d buffer size %ld\n", myProcid, pData->proc_index[p],  pData->buffer_size[p]);
            MPI_return = MPI_Isend(&pData->buffer_size[p], 1, MPI_LONG, pData->proc_index[p], MOAMSA_SEND_RECEIVE_TAG, MOAMSA_COMM_WORLD, &pData->mpi_requests[pData->sendRequests]);
            pData->sendRequests ++;
            //printf ("[%d] MPI_return %d \n", myProcid, MPI_return);
            MPI_return = MPI_Isend(pData->buffer[p],  pData->buffer_size[p], MPI_LONG_LONG, pData->proc_index[p], MOAMSA_SEND_RECEIVE_TAG, MOAMSA_COMM_WORLD, &pData->mpi_requests[pData->sendRequests]);
            pData->sendRequests ++;
            //printf ("[%d] buffer sent MPI_return %d \n", myProcid, MPI_return);
#ifndef NDEBUG
            //printf ("[%d]>sendOC: Element Count[%ld] To Process [%d]\n", myProcid,  pData->buffer_size[p], pData->proc_index[p]);
            sprintf (msg, "[%d]>sendOCtoLowerProcessors: Element Count[%ld] To Process [%d]\n", myProcid,  pData->buffer_size[p], pData->proc_index[p]);
            mprintf (dbglevel, msg, 2);
            //printf ("(  wn,   pi,   ci,     cs)\n");
            mprintf (dbglevel, "(  wn,   pi,   ci,     cs)\n", 2);
            for (l=0; l< pData->buffer_size[p]; l+=pData->commBufSize) {
                //printf ("(%4lld, {%lld", pData->buffer[p][l], pData->buffer[p][l+1]);
                sprintf (msg, "(%4lld, {%lld", pData->buffer[p][l], pData->buffer[p][l+1]);
                mprintf (dbglevel, msg, 2);
                bIterator = 2;
                for (k=1;k<pData->seqNum;k++) {
                    //printf (", %lld", pData->buffer[p][bIterator]);
                    sprintf (msg, ", %lld", pData->buffer[p][l+bIterator]);
                    mprintf (dbglevel, msg, 2);
                    bIterator ++;
                }
                //printf ("}, {%4lld", pData->buffer[p][l+2]);
                sprintf (msg, "}, %6lld)\n", pData->buffer[p][l+pData->commBufSize-1]);
                mprintf (dbglevel, msg, 2);
            }
#endif
        } /*End check if Lower Processor*/
        } /*End Loop for Processors*/
#ifndef NDEBUG
        sprintf (msg, "[%d]>sendOCtoLowerProcessors:sent OC to total processes [%d] sendRequests %d\n", myProcid,pData->depProcCount, pData->sendRequests);
        mprintf (dbglevel, msg, 2);
#endif
    }
    return 0;
}


/*******************************************************************
	Function: receiveOC
	called from scoring.c:getScore()
	int source; // read from MPI_Probe return status.MPI_SOURCE
	long cellIndex; // Cell Index to be received, return 1 if found, 0 otherwise
	long * cellScore; // output score received
*******************************************************************/
int receiveOC (ProcessData * pData, ScoringData * sData) {
    MOATypeDimn k, bIterator;
    long RwaveNo;
    MOATypeInd j;
    MOATypeElmVal * buffer = NULL;
    int source;
    long i, buffer_size;
    MPI_Status status;
#ifndef NDEBUG
    char msg[SHORT_MESSAGE_SIZE];
    int dbglevel = 1;
#endif
    int received = -1;
    int MPI_return = 0;
#ifndef NDEBUG
    sprintf (msg, "[%d] Recv in wave %ld  PI {%lld", myProcid,  pData->waveNo, sData->p_index[0]);
    for (k=1;k<pData->seqNum;k++)
        sprintf (msg, "%s, %lld", msg, sData->p_index[k]);        
    sprintf (msg, "%s} Receiving: { %lld", msg, sData->gm_index[0]);
    for (k=1;k<pData->seqNum;k++)
        sprintf (msg, "%s, %lld", msg, sData->gm_index[k]);
    sprintf (msg, "%s}\n", msg);
    mprintf (dbglevel, msg, 3);
#endif	
    sData->score = 0;
    buffer_size = 0;	
    /*printf ("[%d] Recv in wave %ld  PI {%lld", myProcid,  pData->waveNo, sData->p_index[0]);
    for (k=1;k<pData->seqNum;k++)
        printf (", %lld", sData->p_index[k]);
        
    printf ("} Receiving: { %lld", sData->gm_index[0]);
    for (k=1;k<pData->seqNum;k++)
        printf (", %lld", sData->gm_index[k]);
     
    printf ("[%d] Will Probe in receiveOC \n", myProcid);    
    fflush(stdout);*/
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MOAMSA_COMM_WORLD, &MPI_return, &status);
    //while (received == -1) {
    //printf ("} MPI ret %d condition %d\n", MPI_return, MPI_return == 1);
        //fflush (stdout);
    while (MPI_return == 1) {
        source = status.MPI_SOURCE;
        MPI_return = MPI_Recv(&buffer_size, 1, MPI_LONG, source, MOAMSA_SEND_RECEIVE_TAG, MOAMSA_COMM_WORLD, &status);
#ifndef NDEBUG
        checkMPIErrorCode ("mmDst", 1, 2, MPI_return);
#endif
        if (buffer_size > 0) {
            buffer = mmalloc (((MOATypeInd) buffer_size) * ((MOATypeInd) sizeof *buffer));
            if (buffer == NULL) {
#ifndef NDEBUG
                sprintf(msg, "[%d]>receiveOC: Couldn't allocate memory for OCin receiving buffer, exiting!\n", myProcid);
                mprintf(0, msg, 3);
    #endif
                return -1;
            } 
            MPI_return = MPI_Recv(buffer, buffer_size, MPI_LONG_LONG, source, MOAMSA_SEND_RECEIVE_TAG, MOAMSA_COMM_WORLD, &status);

#ifndef NDEBUG
            checkMPIErrorCode ("mmDst", 1, 2, MPI_return);
            sprintf  (msg, "[%d]>receiveOC: Received buffer: size: %ld from [%d]\n", myProcid, buffer_size, status.MPI_SOURCE);
            mprintf (dbglevel, msg, 3);		  
            sprintf (msg, "Recieved idx(  wn,  pi,  ci,    cs,  fp)\n");
            mprintf (dbglevel, msg, 3);
#endif
            for (i=0;i<buffer_size;i+=pData->commBufSize) {
                bIterator = 0;
                RwaveNo = (long) buffer[i];
                bIterator++;
#ifndef NDEBUG
                sprintf (msg, "[%d]>receiveOC: loop (%ld / %ld)\n", myProcid, i, buffer_size);
                mprintf (dbglevel + 1, msg, 3);
#endif
                if (pData->OCin[RwaveNo].wavesOC == 0)
                    pData->OCin[RwaveNo].WOCI = mmalloc ((MOATypeInd) sizeof *pData->OCin[RwaveNo].WOCI);
                else
                    pData->OCin[RwaveNo].WOCI = realloc (pData->OCin[RwaveNo].WOCI, ((MOATypeInd) (pData->OCin[RwaveNo].wavesOC + 1)) * ((MOATypeInd) sizeof *pData->OCin[RwaveNo].WOCI));

                if (pData->OCin[RwaveNo].WOCI == NULL) {
#ifndef NDEBUG
                    sprintf(msg, "[%d]>receiveOC: Couldn't allocate memory for OCin[].WOCI, exiting!\n", myProcid);
                    mprintf(0, msg, 1);
#endif
                    return -1;
                }
                pData->OCin[RwaveNo].WOCI[pData->OCin[RwaveNo].wavesOC].cellIndex = mmalloc  (((MOATypeInd) pData->seqNum)* ((MOATypeInd) sizeof *(pData->OCin[RwaveNo].WOCI[pData->OCin[RwaveNo].wavesOC].cellIndex)));
                if (pData->OCin[RwaveNo].WOCI[pData->OCin[RwaveNo].wavesOC].cellIndex == NULL) {
#ifndef NDEBUG
                    sprintf(msg, "[%d]>receiveOC: Couldn't allocate memory for OCin[].WOCI[].partIndex, exiting!\n", myProcid);
                    mprintf(0, msg, 1);
#endif
                    return -1;
                }

                if (received < 0)
                    received = -1;
                 for (k=0;k<pData->seqNum;k++) {
                    pData->OCin[RwaveNo].WOCI[pData->OCin[RwaveNo].wavesOC].cellIndex[k] = (MOATypeInd) buffer[i+bIterator];       
                    if (received == -1)
                        if (pData->OCin[RwaveNo].WOCI[pData->OCin[RwaveNo].wavesOC].cellIndex[k] != sData->gm_index[k])
                            received = -2;                                
                    bIterator++;
                }

                pData->OCin[RwaveNo].WOCI[pData->OCin[RwaveNo].wavesOC].cellScore = buffer[i+bIterator];

                if (received == -1) {
                    received = 0;
                    sData->score = pData->OCin[RwaveNo].WOCI[pData->OCin[RwaveNo].wavesOC].cellScore;
                }
#ifndef NDEBUG
                sprintf (msg, " wave %4ld(%4ld {%4lld", RwaveNo, RwaveNo, pData->OCin[RwaveNo].WOCI[pData->OCin[RwaveNo].wavesOC].cellIndex[0]);
                for (k=1;k<pData->seqNum;k++)                    
                    sprintf (msg, "%s, %4lld", msg, pData->OCin[RwaveNo].WOCI[pData->OCin[RwaveNo].wavesOC].cellIndex[k]);
                sprintf (msg, "%s} %6lld %4d)\n ", msg, pData->OCin[RwaveNo].WOCI[pData->OCin[RwaveNo].wavesOC].cellScore, status.MPI_SOURCE);
                mprintf (dbglevel, msg, 3);
#endif

                pData->OCin[RwaveNo].wavesOC ++;

            } /*End Buffer loop*/
            if (buffer != NULL)
                free (buffer);
            buffer = NULL;
        }/* End Buffer size check*/
    //} /*End Probing success condition*/
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MOAMSA_COMM_WORLD, &MPI_return, &status);
    } /*End while receive required OC*/
#ifndef NDEBUG
    if (received == 0)
        sprintf  (msg, "[%d]>receiveOC: FOUND cell  score [%lld]\n", myProcid, sData->score);
    else
        sprintf  (msg, "[%d]>receiveOC: NOT FOUND\n", myProcid);
    mprintf (dbglevel, msg, 3);		  
#endif
    //printf ("[%d] receiveOC flag %d score %lld\n", myProcid, received, sData->score);
    return received;
}

/***********************************************************************
	Function: checkRecvOC
		returns score of cell (cellIndex) if received in previous receive command
	long cellIndex; // searching for this cellIndex
	long * cellScore; // if cellIndex is found, cellScore is update
	long * startIndex; // Partition start index, to start searching from within the current wave
	long waveNo; // current waveNo to search only the 2 waves before it for dependency, actually one wave, but I think I am paranoid
***********************************************************************/
int checkRecvOC (ProcessData * pData, WavesData * wData, MOATypeShape * cellIndex, MOATypeElmVal * score, MOATypeInd * startIndex) {
    MOATypeDimn k;
    MOATypeInd i, j;
    MOATypeInd startingWave, startingIndex;
    int received = -1;
#ifndef NDEBUG
    char msg[MID_MESSAGE_SIZE];
    int dbglevel = 1;
#endif

#ifndef NDEBUG
    sprintf (msg, "[%d]checkRecvOC: Current Wave: %ld - Current Partition: %ld \n", myProcid, pData->waveNo, pData->partNo);
    mprintf(dbglevel, msg, 3);
    if (startIndex != NULL)
        sprintf(msg, "[%d]>checkRecvOC: Incoming Overlapping cells start=[%lld] waveNo=[%ld] looking for cell Index: {%lld", myProcid, (*startIndex), pData->waveNo, cellIndex[0]);
    else
        sprintf(msg, "[%d]>checkRecvOC: Incoming Overlapping cells waveNo=[%ld] looking for cell Index: {%lld", myProcid, pData->waveNo, cellIndex[0]);
    for (k=1;k<pData->seqNum;k++)
        sprintf (msg, "%s, %lld", msg, cellIndex[k]);
    sprintf (msg, "%s}\n", msg);
    mprintf(dbglevel, msg, 3);
#endif
    i = 0;
    if (pData->waveNo > pData->seqNum)
            startingWave = pData->waveNo - pData->seqNum;
    else 
            startingWave = 0;
    if (startIndex != NULL)
            startingIndex = (*startIndex);
    else 
            startingIndex = 0;
    for (j=startingWave;j<pData->waveNo && j <wData->wavesTotal && received==-1;j++) {
#ifndef NDEBUG
        sprintf (msg, "[%d] OCin[%lld].wavesOC %ld{\n#### ", myProcid, j, pData->OCin[j].wavesOC);
        mprintf(dbglevel, msg, 3);
#endif
        for (i=startingIndex;((i<pData->OCin[j].wavesOC) && (received == -1));i++) {
#ifndef NDEBUG
            sprintf(msg, "%lld ", pData->OCin[j].WOCI[i].cellIndex);
            mprintf(dbglevel, msg, 3);
#endif
            received = 0;
            for (k=0;k<pData->seqNum;k++) {
                if (cellIndex[k] != pData->OCin[j].WOCI[i].cellIndex[k]) {
                    received = -1;
                    break;
                }
            }
            if (received == 0)
                (*score) = pData->OCin[j].WOCI[i].cellScore;
        }
#ifndef NDEBUG
        mprintf(dbglevel, "\n}\n", 3);
#endif
    }
    if ((startIndex != NULL) && (i > 0))
        (*startIndex) = i-1;
    else if (startIndex != NULL)
        (*startIndex) = 0;
#ifndef NDEBUG
    if (received == 0)
        sprintf  (msg, "[%d]>checkRecvOC: FOUND\n", myProcid);
    else
        sprintf  (msg, "[%d]>receiveOC: NOT FOUND\n", myProcid);
    mprintf(dbglevel, msg, 3);
#endif
    //printf ("checkRecvOC flag %d score %lld\n", received, (*score));
    return received;
}

/************************************************************************
	Function: checkPrevPartitions
		returns score of cell (cellIndex) if computed in previous partitions
************************************************************************/
int checkPrevPartitions (ProcessData * pData, MOATypeShape * cellIndex, MOATypeElmVal * score) {
#ifndef NDEBUG
    char msg[MID_MESSAGE_SIZE];
    int dbglevel = 5;
#endif
    long i, j, startingWave, endingWave;
    MOATypeDimn k;
    int found = -1;
#ifndef NDEBUG
    sprintf (msg, "in checkPrevPartitions w %ld wavesOC %ld ", pData->waveNo-1, pData->OCout[pData->waveNo-1].wavesOC);
    mprintf (dbglevel, msg, 1);
#endif
    if (pData->waveNo - pData->seqNum >= 0)
        startingWave = pData->waveNo - pData->seqNum;
    else
        startingWave = 0;
    if (Mode != Distributed)
        endingWave = 1;
    else 
        endingWave = pData->waveNo;
    
    // search only k (seqNum) waves before the current one, for the required partition (dependency go only that far).
    for (i=startingWave;(i<endingWave) && (found == -1);i++) { 
        for (j=0;((j<pData->OCout[i].wavesOC) && (found == -1));j++) {
            found = 0; /* if it remains zero, then cell found*/
            for (k=0;k<pData->seqNum;k++) {
                if (cellIndex[k] != pData->OCout[i].WOCO[j].cellIndex[k]) {
                    found = -1;
                    break;
                }
            }
            if (found == 0)
                (*score) = pData->OCout[i].WOCO[j].cellScore;
        }
    }
#ifndef NDEBUG
    sprintf (msg, " found %d\n", found);
    mprintf (dbglevel, msg, 1);
#endif
    //printf ("checkPrevPartitions flag %d score %lld\n", found, (*score));
    return found;
}
