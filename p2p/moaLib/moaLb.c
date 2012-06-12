   
/**********************************************************
* Author: Manal Helal																			*
* Last Modification Date: Fri 12 Jan 2007 03:39:51 AM EST *
* Project : MMSA - Multiple Sequence Alignment Based on 	*
* 					Mathematics of Arrays - PhD Experimentation		*
* File: moaLb.c, a library for testing the basic MOA 			*
* functions & some MSA functions in isolation.				 		*
***********************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>      /* include the character library */
#include <math.h>
#include <time.h>
#include <limits.h>
#include "../main.h"
#include "../moa.h"
#include "../globals.h"
#include "../utils.h"
#include "../scores.h"

//const char GABCHAR = '-';


int checkCells (MOA_rec * MOA_in, MOATypeInd * sendCells, MOATypeInd * recvCells, MOATypeInd * CompCells) {
  MOATypeShape * ind = NULL;
  MOATypeInd i;
  (*sendCells) = 0;
  (*recvCells) = 0;
  (*CompCells) = 0;
  ind = mmalloc ( MOA_in->dimn * sizeof *ind);
  for (i = 0; i <  MOA_in->elements_ub; i++) {
    Gamma_Inverse(i, MOA_in->shape,  MOA_in->dimn, &ind, 1);
    if (isHigherBorderCellandNotLower (ind, MOA_in->dimn, MOA_in->shape) == 1) {
      (*sendCells)++;
      (*CompCells)++;
    }
    else if (isLowerBorderCell (ind, MOA_in->dimn) == 1)
      (*recvCells)++;
    else
      (*CompCells)++;
  }
  free(ind);
  return 0;
}

int generateRandomSequence (MOATypeShape seqLen, char * * sequence) {
    MOATypeShape i;
    double r;
    int x, M;
    M = 4;
    for (i=0;i<seqLen;i++) {
        r = (double) rand() / ((double) (RAND_MAX) + (double) (1));
        x = (int) (r * M) +1;
        if (x > 15) {
            printf ("Error generating Random number in DNA range.\n");
            return -1;
        }
        (*sequence)[i] = encodeDNA (x);
    }
    return 0;
}

int IsHighBorderCellinPart (MOATypeDimn dimn, MOATypeShape * shape, long pSize, MOATypeShape * cellIndex, MOATypeShape * partIndex) {
    int found = 0, ret = 0;
    MOATypeInd i;
    MOATypeDimn k;
    MOATypeShape * cellmInd = NULL;
    MOA_rec * part;
    getPartition (partIndex, dimn, shape, &part, pSize);

    for (i=0;(i<part->elements_ub) && (found == 0);i++) {            
        found = 0;
        for (k=0;(k<dimn) && (found == 0);k++) {
            if (part->indexes[i][k] != cellIndex[k]) {
                found = -1;
                break;
            }
        }
        if (found == 0){
            cellmInd = mcalloc (dimn, sizeof *cellmInd);
            Gamma_Inverse(i, part->shape, part->dimn, &cellmInd, 1);
            for (k=0;(k<dimn) && (ret == 0);k++) {
                if (cellmInd[k] == part->shape[k]-1)
                    ret = 1;
            }	
            free (cellmInd);
            cellmInd = NULL;
        }	
    }
    deleteMOA (part);
    return ret;
}

int IsLowBorderCellinPart (MOATypeDimn dimn, MOATypeShape * shape, long pSize,  MOATypeShape * cellIndex, MOATypeShape * partIndex) {
    int found = 0, ret = 0;
    MOATypeInd i; 
    MOATypeDimn k;
    MOATypeShape * cellmInd = NULL;
    MOA_rec * part;
    getPartition (partIndex, dimn, shape, &part, pSize);

    for (i=0;(i<part->elements_ub) && (found == 0);i++) {
        found = 0;
        for (k=0;(k<dimn) && (found == 0);k++) {
            if (part->indexes[i][k] != cellIndex[k]) {
                found = -1;
                break;
            }
        }
        if (found == 0){
            found = 1;
            cellmInd = mcalloc (dimn, sizeof *cellmInd);
            Gamma_Inverse(i, part->shape, part->dimn, &cellmInd, 1);
            for (k=0;(k<dimn) && (ret == 0);k++) {
                    if (cellmInd[k] == 0)
                            ret = 1;
            }	
            free (cellmInd);
            cellmInd = NULL;
        }	
    }
    deleteMOA (part);
    return ret;
}


void Menu (int *choice) {

    char    local;

    printf ("\n1: Re-Define MOA\n");
    printf ("2: Print MOA\n");
    printf ("3: Gamma\n");
    printf ("4: Gamma_inverse\n");
    printf ("5: Take\n");
    printf ("6: Drop\n");
    printf ("7: MOAGetLowerNeighbors\n");
    printf ("8: MOAGetHigherNeighbors\n");
    printf ("10: Get Send/Recv/Comp Borders Number for a Range of Partitioning Sizes\n");
    printf ("11: Get Send/Recv/Comp Borders Number for a Range of Dimensions, Shapes and Partitioning Sizes\n");
    printf ("12: Get Partitions Number\n");
    printf ("13: Partition & Distribute\n");
    printf ("14: get Sizes\n");
    //printf ("15: New Pattern Parallelism\n");
    printf ("16: Get Parts Total in a Wave\n");
    printf ("17: Parts Total Growth in Waves as Dimension grows\n");
    printf ("18: Analyze Cell Index in a Part\n");
    printf ("19: get Part's Proc\n");
    printf ("20: get Proc's Next Part\n");
    printf ("21: get Proc's Previous Part\n");
    printf ("22: get Partition Index of a Cell Index\n");
    printf ("23: get Higher Partitions\n");
    printf ("24: get Partition Details\n");
    printf ("25: get Dependand Processors\n");
    printf ("26: Compute Scores\n");
    printf ("27: Trace Back\n");
    printf ("28: Compute Scores for one Partition and prepapre OC\n");
    printf ("29: Calculate SP Score from a Fasta Alignment File\n");
    printf ("30: Vector Pythagorean Distance From Origin\n");
    printf ("31: 2 Vectors Pythagorean Distance\n");
    printf ("32: Vector Pythagorean Direction\n");
    printf ("33: 2 Vectors Pythagorean Direction Difference\n");
    printf ("34: 2 Vectors MoA Difference\n");
    printf ("35: 2 Vectors Angle\n");
    printf ("36: 2 Vectors Ranked Angle\n");
    printf ("37: Is Point included in Hyperplane\n");
    printf ("38: Get the Wave Length\n");
    printf ("39: Get Partition's Distance from the wave middle\n");
    printf ("40: Read Sequences and convert to Fasta\n");
    printf ("41: get Middle Partition for a wave and a dimension\n");
    printf ("42: Generate Random Sequence\n");
    printf ("43: Convert Fasta Alignment Format to MSF Format\n");
    printf ("0: Quit\n");
}


void initGlobalVariables () {  
int x =1;
    pdebug = 1; 
    threadnum = 1;
    myProcid = 0;
    ClusterSize = 1;
    /*open once to erase previous contents then close */
    strcpy(outPrefix , "lb");
    strcpy(outputfilename, "m");
    if (init_output() != 0) {
        printf("Error inititliazing output files. Exiting\n");
        return ;
    }

    Epsilons = 0;
    AlignmentType = Global;
    SchedMethod = RR;
    maxAlignmentsNumber = 20;
    strcpy(outputfilename, "mm");
    RestoreFlag = 0;
    Mode = Sequential;
    gapOpenning = -4;
    gapExtension = -2;
    mismatchscore = -1;
    matchscore = 1;
    initializationMethod = IndexProductInit;
    alignResolution = 100;
    prevNow = NULL;
    prevNow = (struct tm *) mmalloc (sizeof(struct tm));
    currNow = NULL;
    currNow = getTime();
    if (currNow == NULL) {
        printf ("Could nor read current time. Exiting.\n");
	return;
    }
    prevNow->tm_hour = currNow->tm_hour;
    prevNow->tm_isdst = currNow->tm_isdst;
    prevNow->tm_mday = currNow->tm_mday;
    prevNow->tm_min = currNow->tm_min;
    prevNow->tm_mon = currNow->tm_mon;
    prevNow->tm_sec = currNow->tm_sec;
    prevNow->tm_wday = currNow->tm_wday;
    prevNow->tm_yday = currNow->tm_yday;
    prevNow->tm_year = currNow->tm_year;

} 
 
int main (int argc, char ** argv) {
    char sfilename[SHORT_MESSAGE_SIZE];
    char sfilename_2[SHORT_MESSAGE_SIZE];
    ProcessData * pData = NULL; 
    ScoringData * sData = NULL; 
    WavesData * wData = NULL; 
    TracebackData * tbData = NULL;  
    MOATypeShape * shape = NULL;
    MOATypeShape * ind = NULL;
    MOATypeShape * ind2 = NULL;
    MOATypeShape * posDimn = NULL, fromshape, toshape;
    MOATypeInd flatIndex, partIndex, duplicatesTotal;
    long fromSize, toSize, MypartsTotal, wave, PartOrder, strides;
    MOATypeDimn dimn, fromdimn, todimn, k;
    MOATypeInd sendCells, recvCells, CompCells;
    MOATypeInd i, j, l, rslt_ub; 
    int choice, ret, * validPartitions = NULL;
    MOATypeShape * rslt = NULL;
    MOATypeShape * rsltInd = NULL;
    MOA_rec * MOA1 = NULL; 
    MOA_rec * MOA2 = NULL;
    char * * sequences, * * seqName;
    tbData = mmalloc ((MOATypeInd) sizeof *tbData);
    if (tbData == NULL) {
        printf ("Could not create trace back data memory. Exiting\n");
        return -1;
    }
  
    initGlobalVariables (); 
    createMOAStruct(&MOA1);
    ind =  mmalloc ( ((MOATypeInd) MOA1->dimn) * ((MOATypeInd) sizeof *ind));
    ind2 = mmalloc ( ((MOATypeInd) MOA1->dimn) * ((MOATypeInd) sizeof *ind2));
    posDimn = mmalloc ( ((MOATypeInd) MOA1->dimn) * ((MOATypeInd) sizeof *posDimn));

   Menu (&choice);
    do {
    printf ("\nEnter op code:" );    
    scanf ("%d", &choice);
    printf ("choice %d \n", choice);
    
    switch (choice) {
        default:
            printf ("Invalid menu choice - try again\n");
            break;
      
        case 0:
            printf ("\n good Bye\n");
            break;
        case 1:     
            deleteMOA (MOA1); 
            MOA1 = NULL;
            if (ind != NULL)
                free(ind);
            ind = NULL;
            if (ind2 != NULL)
                free(ind2);
            ind2 = NULL;
            if (posDimn != NULL)
                free(posDimn);
            posDimn = NULL;
            if (shape != NULL)
                free(shape);
            shape = NULL;
            createMOAStruct(&MOA1);
            printf ("\nEnter dimensionality:");
            scanf("%lld", &MOA1->dimn);
            ind = mmalloc ( ((MOATypeInd) MOA1->dimn) * ((MOATypeInd) sizeof *ind));
            ind2 = mmalloc ( ((MOATypeInd)MOA1->dimn) * ((MOATypeInd) sizeof *ind2));
            posDimn = mmalloc ( ((MOATypeInd) MOA1->dimn) * ((MOATypeInd) sizeof *posDimn));
            shape = mmalloc ( ((MOATypeInd) MOA1->dimn) * ((MOATypeInd) sizeof *shape));
            for (i = 0; i < MOA1->dimn; i++) {
                printf ("\nEnter shape at dimension %lld:", i);
                scanf("%lld", &shape[i]);  
            }

            if (createMOA(shape /* shape*/, MOA1->dimn /* dimension*/, MOA1 /* MOA structure*/, -1, 0) != 0)
                break;
            if (MOA1->shape == NULL) {
                printf ("Error allocating memory for MoA structure.\n");
            }
            freeProcessMemory (&pData, &sData, &wData);
            printf ("\nEnter partition size:");
            scanf("%lld", &flatIndex);
            ret = initProcessMemory(&pData, &sData, &wData, MOA1->dimn, MOA1->shape, NULL, NULL, 1, (long) flatIndex);  
            if (ret != 0) {
                mprintf (0, "Error Initializing Process Data, Exiting\n", 1);
                fflush (stdout);
                return -1;
            }
            
            /*printMOA(MOA1); // we didn't create the tensor to print it*/
            if (shape != NULL)
                free(shape);
            shape = NULL;
            printf("created the MOA\n");
            break;
    case 2:  
        if (MOA1 != NULL) {
            printMOA1(MOA1, 0);
            printMOA1(MOA1, 1);
        }
        else
            printf ("Please Define an MOA Structure first.\n");
        break;
    case 3:
        if (MOA1 != NULL) {
            for (i = 0; i < MOA1->dimn; i++) {
                printf ("\nEnter index at dimn %ld:", i);
                scanf("%lld", &ind[i]);
            }      
            flatIndex = Gamma(ind, MOA1->dimn, MOA1->shape,  MOA1->dimn, 1);
            printf("\n the flatIndex = %lld \n", flatIndex);
        }
        else
            printf ("Please Define an MOA Structure first.\n");
        break;
    case 4:
        if (MOA1 != NULL) {
            printf ("\nEnter flat index :" );
            scanf("%lld", &flatIndex);
            Gamma_Inverse(flatIndex, MOA1->shape,  MOA1->dimn, &ind, 1);
            printf("%lld", ind[0]);
            for (i = 1; i < MOA1->dimn; i++) {
                printf(", %lld", ind[i]);
            }
        }
        else
            printf ("Please Define an MOA Structure first.\n");
        break;
    case 5:
        if (MOA1 != NULL) {
            createMOAStruct (&MOA2);
            for (i = 0; i <  MOA1->dimn; i++) {  
                printf ("\nEnter take index at dimn %ld:", i);
                scanf("%lld", &ind[i]);
            }
            Take (ind, MOA1->dimn, MOA1, &MOA2);
            printMOA_scr(MOA2, 0);
            printMOA_scr(MOA2, 1);
            deleteMOA (MOA2);  
            MOA2 = NULL;
        }
        else
            printf ("Please Define an MOA Structure first.\n");
        break;
    case 6:
        if (MOA1 != NULL) {
            createMOAStruct (&MOA2);
            for (i = 0; i <  MOA1->dimn; i++) {
                printf ("\nEnter Drop index at dimn %lld:", i);
                scanf("%lld", &ind[i]);
            }
            DropInd(ind, MOA1->dimn, MOA1->shape, &MOA2); 
            printMOA_scr(MOA2, 0);
            printMOA_scr(MOA2, 1);
            deleteMOA (MOA2);  
            MOA2 = NULL;
        }
        else
            printf ("Please Define an MOA Structure first.\n");
        break;
    case 7:
	if (MOA1 != NULL) {
            createMOAStruct (&MOA2);
            for (i = 0; i <  MOA1->dimn; i++) {
                printf ("\nEnter Cell index at dimn %lld:", i);
                scanf("%lld", &ind[i]);
            }
            getLowerNeighbors (2, ind, MOA1->dimn, MOA1->shape, &MOA2);
            printf("With Indices: ");
            for (i=0;i<MOA2->elements_ub; i++) {                
                printf ("{%lld", MOA2->indexes[i][0]);
                for (k=1;k<MOA2->dimn;k++) {
                    printf(", %lld", MOA2->indexes[i][k]);
                }
                printf ("} ");
            }            
            deleteMOA (MOA2); 
            MOA2 = NULL;
        }
	else
            printf ("Please Define an MOA Structure first.\n");
        break;
    case 8:
        if (MOA1 != NULL) {
            createMOAStruct (&MOA2);
            for (i = 0; i <  MOA1->dimn; i++) {
                printf ("\nEnter Cell index at dimn %lld:", i);
                scanf("%lld", &ind[i]);
            }
            printf ("\nEnter The Stride Size:");
            scanf("%lld", &flatIndex);
            if (getHigherNeighbors (flatIndex, ind, MOA1->dimn, MOA1->shape, &MOA2) == 0) {
                printf("With Indices: ");
                for (i=1;i<MOA2->elements_ub; i++) {                
                    printf ("{%lld", MOA2->indexes[i][0]);
                    for (k=1;k<MOA2->dimn;k++) {
                        printf(", %lld", MOA2->indexes[i][k]);
                    }
                    printf ("} ");
                }            
                printf("\n");
                deleteMOA (MOA2);
                MOA2 = NULL;
            }
        }
        else
            printf ("Please Define an MOA Structure first.\n");
        break;
    case 10:
        printf ("\nEnter The tensor Dimensionality:");
        scanf("%lld", &dimn);
        shape = NULL;
        shape = mmalloc (dimn * sizeof *shape);
        printf ("\nEnter The Starting Size of the Partition:");
        scanf("%ld", &fromSize);
        printf ("\nEnter The Ending Size of the Partition:");
        scanf("%ld", &toSize);
        printf ("Cells Analysis for:\n");
        printf ("Partition \t\tTot_Comm \t\tSend  \t\tRecv  \t\tCompute \t\tTotal\n");
        for (j = fromSize; j <  toSize; j++) {
            createMOAStruct (&MOA2);
            MOA2->dimn = dimn;
            for (i = 0; i <  dimn; i++) {
                shape[i] = j;
            }
            createMOA(shape /* shape*/, dimn /* dimension*/, MOA2 /* MOA structure*/, -1, 0);
            checkCells (MOA2, &sendCells, &recvCells, &CompCells);
            printf ("%lld \t\t%lld \t\t%lld \t\t%lld \t\t%lld \t\t%lld\n", j,  sendCells+recvCells, sendCells, recvCells, CompCells, MOA2->elements_ub);
            deleteMOA (MOA2);
        }
        if (shape != NULL)
            free (shape);
        break;
    case 11:
        printf ("\nEnter The tensor Starting Dimensionality:");
        fflush(stdout);
        scanf("%lld", &fromdimn);
        printf ("\nEnter The tensor Ending Dimensionality:");
        scanf("%lld", &todimn);
        printf ("\nEnter The tensor Starting Shape:");
        scanf("%lld", &fromshape);		
        printf ("\nEnter The tensor Starting Shape:");
        scanf("%lld", &toshape);		
        printf ("\nEnter The Starting Size of the Partition:");
        scanf("%ld", &fromSize);
        printf ("\nEnter The Ending Size of the Partition:");
        scanf("%ld", &toSize);
        printf ("Comm Border Cells for:\n");
        printf ("Dimn\tShape\tPartition \t\tTot_Comm \t\tSend  \t\tRecv  \t\tCompute \t\tTotal\n");
        for (k = fromdimn; k <  todimn; k++) {
            for (l = fromshape; l <  toshape; l++) {
                for (j = fromSize; j <  toSize; j++) {
                    createMOAStruct (&MOA2);
                    //printf ("\na\n");
                    MOA2->dimn = k;
                    shape = NULL;
                    shape = mmalloc (k * sizeof *shape);
                    for (i = 0; i <  k; i++) {
                        shape[i] = l;
                    }
                    createMOA(shape /* shape*/, dimn /* dimension*/, MOA2 /* MOA structure*/, -1, 0);
                    //printf ("\n1\n");
                    checkCells (MOA2, &sendCells, &recvCells, &CompCells);
                    //printf ("\n2\n");
                    printf ("%lld \t%lld \t%lld \t\t%lld \t\t%lld \t\t%lld \t\t%lld \t\t%lld\n", k, l, j,  sendCells+recvCells, sendCells, recvCells, CompCells, MOA2->elements_ub);
                    //printf ("\n3\n");
                    deleteMOA (MOA2);
                    //printf ("\n4\n");
                    if (shape != NULL)
                        free (shape);
                    shape = NULL;
                }
            }
         }
        break;
    case 12:
        //nk = mk(p-2) + (mk-1) +2
        // nk = mk(p-1) +1
        //mk=(nk-1)/(p-1) mk is no pf partitions at each dimension k and p is the partitioning size, and nk is the length at dimn k
        // equation for duplicates at dimension k
        // D(k) = D(k-1) * nk + (Pi= 0 to k (mi-1)) * 2^(k-1)
        printf ("\nEnter The tensor Dimensionality:");
        scanf("%lld", &dimn);
        shape = NULL;
        shape = mmalloc (dimn * sizeof *shape);
        for (i = 0; i < dimn; i++) {
            printf ("\nEnter shape at dimension %lld:", i);
            scanf("%lld", &shape[i]);
        }
        createMOAStruct (&MOA2);
        createMOA(shape /* shape*/, dimn /* dimension*/, MOA2 /* MOA structure*/, -1, 0);
        printf ("\nEnter The Size of the Partition:");
        scanf("%ld", &fromSize);
        getWavesPartsDuplicates(MOA2->dimn, MOA2->shape, fromSize, &wData->wavesTotal, &wData->partsTotal, &wData->duplicatesTotal);
        printf ("Waves Total is %ld, Partitions Total is %ld Duplicates Total is %lld\n", wData->wavesTotal, wData->partsTotal, wData->duplicatesTotal);
        deleteMOA (MOA2);
        if (shape != NULL)
            free (shape);
        break;
    case 13:
        if (MOA1 != NULL) {  
            printf ("\nEnter The Size of the Cluster:");
            scanf("%d", &ClusterSize);
            printf ("\nEnter myProcid:");
            scanf("%d", &myProcid); 
            printf ("\nEnter Epsilons:");
            scanf("%d", &Epsilons); 
            //Epsilons = 90;
            wData->wavesTotal = calcWaves (pData, wData);
            printf ("proc %d has %ld/%ld parts in %ld waves starts in wave %ld, part %ld Part Index {%lld", myProcid, pData->partitionsCount, wData->partsTotal, wData->wavesTotal, pData->waveNo, pData->partNo, wData->partsInWaveIndices[pData->waveNo][pData->partNo][0]);
            for (k=1;k<dimn; k++)
                printf (", %lld", wData->partsInWaveIndices[pData->waveNo][pData->partNo][k]);
            printf ("}\n");
            
	}
	else
            printf ("Please Define an MOA Structure first.\n");
        break;
    case 14:
      getSizes ();
	/*
    case 15:
      for (i=2;i<21;i++) {
	rslt = calcWaves (i);
	printf ("parallelism clusters in one wave for dimn i %ld is: %ld", i, rslt[0]);
	for (j=1;j<i;j++) {
	  printf (", %ld ", rslt[j]);
	}
	printf ("\n");
	free (rslt);
      }
      break;
*/
        break;
    case 16:
        printf ("Enter the Dimension\n");
        scanf("%lld", &dimn);
        printf ("Enter the Wave No\n");
        scanf("%ld", &wave);
        printf ("There are %ld Parts in %lld dimesnion in wave Number %ld\n", getWavePartsTotal (wave, dimn), dimn, wave);
	break;
    case 17:
        printf ("Parts Total Growth with Waves and Dimensions\nDimn\tWave\tParts\t\n");
        for (i=2;i<10;i++) { /*Dimns*/
            for (j=1;j<30;j++)  /*Waves*/
                    printf ("%lld\t%lld\t%ld\n", i, j,  getWavePartsTotal (j, i));
        }
	break;
    case 18:
        printf ("Enter the Dimension\n");
        scanf("%lld", &dimn);
        shape = mmalloc (dimn * sizeof *shape);
        for (i = 0; i < dimn; i++) {
            printf ("\nEnter shape at dimension %ld:", i);
            scanf("%lld", &shape[i]);
        }
        printf ("Enter the Partition Size\n");
        scanf("%ld", &fromSize);
        for (i = 0; i < dimn; i++) {
            printf ("Enter the Cell Index[%ld]\n", i);
            scanf("%lld", &ind[i]);
        }
        for (i = 0; i < dimn; i++) {
            printf ("Enter the Cell Index[%ld]\n", i);
            scanf("%lld", &ind2[i]);
        }
        if (IsHighBorderCellinPart (dimn, shape, fromSize, ind, ind2) == 1) 
            printf ("It is a High Border Cell\n");
        else
            printf ("It is NOT a High Border Cell\n");

        if (IsLowBorderCellinPart (dimn, shape, fromSize, ind, ind2) == 1) 
            printf ("It is a Low Border Cell\n");
        else
            printf ("It is NOT a Low Border Cell\n");
        if (shape != NULL)
            free (shape);
        shape = NULL;
        break;
    case 19:
        if (wData->partsInWaveIndices == NULL) {
            printf ("You have to Calculate Waves first\n");			
        }
        else {
            printf ("Enter wave No\n");
            scanf("%ld", &pData->waveNo);
            printf ("Enter part Order\n");
            scanf("%ld", &pData->partNo);
            myProcid = getProcID (wData, pData->waveNo, pData->partNo);
            printf ("part %ld in Wave %ld with Part Index {%lld", pData->partNo, pData->waveNo, wData->partsInWaveIndices[pData->waveNo][pData->partNo][0]);
            for (i=1;i<pData->seqNum;i++)
                printf (", %lld", wData->partsInWaveIndices[pData->waveNo][pData->partNo][i]);
            printf ("} is processed by %d proc\n", myProcid);
        }	
	break;
    case 20:
        if (wData->partsInWaveIndices == NULL) {
            printf ("You have to Calculate Waves first\n");			
        }
        else {
            printf ("Looping through getNextPartition shows:\nProc ID\t\tWave No\t\tPart No\t\tPart Index\n");
            for (k=0;k<ClusterSize;k++) {
                pData->waveNo = 0;
                pData->partNo = 0;
                myProcid = k;                
                if (getProcID(wData, pData->waveNo, pData->partNo) == myProcid) {
                    printf ("%d\t\t%ld\t\t%ld\t\t{%lld", myProcid,  pData->waveNo, pData->partNo, wData->partsInWaveIndices[pData->waveNo][pData->partNo][0]);
                    for (i=1;i<pData->seqNum;i++)
                        printf (", %lld", wData->partsInWaveIndices[pData->waveNo][pData->partNo][i]);
                    printf ("}\n");
                }
                getNextPartition (wData, &pData->waveNo, &pData->partNo);
                while ((pData->partNo != -1) && (pData->waveNo < wData->wavesTotal)){
                    printf ("%d\t\t%ld\t\t%ld\t\t{%lld", myProcid,  pData->waveNo, pData->partNo, wData->partsInWaveIndices[pData->waveNo][pData->partNo][0]);
                    for (i=1;i<pData->seqNum;i++)
                        printf (", %lld", wData->partsInWaveIndices[pData->waveNo][pData->partNo][i]);
                    printf ("}\n");
                    getNextPartition (wData, &pData->waveNo, &pData->partNo);
                }               
            }
        }
        break;
     case 21:
        if (wData->partsInWaveIndices == NULL) {
            printf ("You have to Calculate Waves first\n");			
        }
        else {
            printf ("Looping through getPrevPartition shows:\nProc ID\t\tWave No\t\tPart No\t\tPart Index\n");
            /*myProcid = 7;
            pData->waveNo = wData->wavesTotal - 1;
            pData->partNo = wData->partsInWave[pData->waveNo]-1;
            tbData->maxCellScore = 0;
            if (getProcID(wData, pData->waveNo, pData->partNo) != myProcid) 
                getPrevPartition (wData, &pData->waveNo, &pData->partNo);
            for (i=0;i<(pData->partitionsCount) && (pData->partNo != -1);i++) {
                printf ("[%d] doing part (%ld, %ld)\n", myProcid, pData->waveNo, pData->partNo);
                getPrevPartition (wData, &pData->waveNo, &pData->partNo);
            }*/
                       
            for (k=0;k<ClusterSize;k++) {
                pData->waveNo = wData->wavesTotal-1;
                pData->partNo = wData->partsInWave[pData->waveNo]-1;
                myProcid = k;       
                if (getProcID(wData, pData->waveNo, pData->partNo) == myProcid) {
                    printf ("%d\t\t%ld\t\t%ld\t\t{%lld", myProcid,  pData->waveNo, pData->partNo, wData->partsInWaveIndices[pData->waveNo][pData->partNo][0]);
                    for (i=1;i<pData->seqNum;i++)
                        printf (", %lld", wData->partsInWaveIndices[pData->waveNo][pData->partNo][i]);
                    printf ("}\n");
                }                
                getPrevPartition (wData, &pData->waveNo, &pData->partNo);
                while ((pData->partNo != -1) && (pData->waveNo >= 0)) {
                    printf ("%d\t\t%ld\t\t%ld\t\t{%lld", myProcid,  pData->waveNo, pData->partNo, wData->partsInWaveIndices[pData->waveNo][pData->partNo][0]);
                    for (i=1;i<pData->seqNum;i++)
                        printf (", %lld", wData->partsInWaveIndices[pData->waveNo][pData->partNo][i]);
                    printf ("}\n");
                    getPrevPartition (wData, &pData->waveNo, &pData->partNo);
                }
            }
             
        }
 
        break;
   case 22:
	if (wData->partsInWaveIndices == NULL) 
            printf ("You have to Calculate Waves first\n");			
	else {		
            for (i=0;i<pData->seqNum;i++) {
                printf ("Enter the Cell Index[%lld]\n", i);
                scanf("%lld", &ind[i]);
                posDimn[i] = 1;
            }
            if (getPartitionIndex (ind, pData->seqNum, pData->seqLen, wData->partitionSize, &ind2) == 0) {
                printf ("The partIndex = {%lld", ind2[0]);
                for (i=1;i<pData->seqNum;i++) 
                    printf (", %lld", ind2[i]);
                printf ("}\n");
            }
            else
                printf ("Invalid cell Index");
        }	
	break;
    case 23:
	if (wData->partsInWaveIndices == NULL) 
            printf ("You have to Calculate Waves first\n");			
	else {		
            for (i=0;i<pData->seqNum;i++) {
                printf ("Enter the Cell Index[%lld]\n", i);
                scanf("%lld", &ind[i]);
                posDimn[i] = 1;
            }
            if (validPartitions != NULL)
                free(validPartitions);
            validPartitions = NULL;
            validPartitions = mcalloc ((MOATypeInd) sData->CalLnCount, ((MOATypeInd) sizeof *validPartitions));
            if (ind2 == NULL) {
                printf ("Couldn't create memory for validPartitions in calcDep. Exiting\n");
                return -1;
            }
            if (MOAGetHigherPartitions (wData->partitionSize, 2, ind, pData->seqNum, pData->seqLen, &MOA2, &validPartitions) > 0) {
                printf ("The Higher Neighboring Partitions are:\n PartIndex\t\tWave No.\tPart No\t\t Processor\n");
                for (i=1;i<MOA2->elements_ub;i++) {
                    if (validPartitions[i] == 1) {
                        printf ("{");
                        for (j=0;j<pData->seqNum;j++) {
                            if (j>0)
                                printf (", ");
                            printf ("%lld", MOA2->indexes[i][j]);
                        }
                        getPartitionPosition (wData, MOA2->indexes[i], &pData->waveNo, &pData->partNo);
                        if (pData->partNo >= 0) 
                            myProcid = getProcID (wData, pData->waveNo, pData->partNo);
                        printf ("}\t\t%ld\t\t%ld/%ld\t\t%d", pData->waveNo, pData->partNo,wData->partsInWave[pData->waveNo], myProcid);                            
                    }
                    printf ("\n");
                }
            }
            else
                printf ("Invalid cell Index");
        }	
	break;
    case 24:
	if (wData->partsInWaveIndices == NULL) 
            printf ("You have to Calculate Waves first\n");			
	else {		
            for (i=0;i<pData->seqNum;i++) {
                printf ("Enter the Cell Index[%lld]\n", i);
                scanf("%lld", &ind[i]);
                posDimn[i] = 1;
            }
            if(getPartitionDetails (pData, wData, &ind2, ind, &pData->waveNo, &pData->partNo, &myProcid) == 0) {
                printf ("Cell's Partition Index = {%lld", ind2[0]);
                for (i=1;i<pData->seqNum;i++) {
                    printf (", %lld", ind2[i]);
                }
                printf ("} in waveNo %ld partNo %ld processed by processor %d\n", pData->waveNo, pData->partNo, myProcid);
            }
            else
                printf ("Invalid cell Index");
        }
	break;
    case 25:
	if (wData->partsInWaveIndices == NULL) 
            printf ("You have to Calculate Waves first\n");			
	else {		
            for (i=0;i<pData->seqNum;i++) {
                printf ("Enter the Cell Index[%lld]\n", i);
                scanf("%lld", &sData->gm_index[i]);
                posDimn[i] = 1;
            }
            sData->score = -1;
            
            if(getPartitionDetails (pData, wData, &sData->p_index, sData->gm_index, &pData->waveNo, &pData->partNo, &myProcid) == 0) {
                printf ("Cell's Partition Index = {%lld", sData->p_index[0]);
                for (i=1;i<pData->seqNum;i++) {
                    printf (", %lld", sData->p_index[i]);
                }
                printf ("} in waveNo %ld partNo %ld processed by processor %d\n", pData->waveNo, pData->partNo, myProcid);
                addOC (pData, sData, wData);
                printf ("Dependent Processors are: {");
                for (i=0;i<pData->OCout[pData->waveNo].WOCO[pData->OCout[pData->waveNo].wavesOC-1].depProc_ub;i++) {
                    if (i>0)
                        printf (", ");
                    printf ("%d", pData->OCout[pData->waveNo].WOCO[pData->OCout[pData->waveNo].wavesOC-1].depProc[i]);
                }
                printf ("}\n");
            }
            else
                printf ("Invalid cell Index");
        }
	break;
    case 26:
        /*Read Sequences Data*/
        printf ("Enter Scoring Type\n");
        scanf("%d", &ret);

        
        printf ("Enter Number of Sequences\n");
        scanf("%lld", &dimn);
        seqName = mmalloc (((MOATypeInd) dimn) * ((MOATypeInd) sizeof *seqName));
        shape = mcalloc ((MOATypeInd) dimn, (MOATypeInd) sizeof *shape);
        sequences = mmalloc (((MOATypeInd) dimn) * (MOATypeInd) sizeof *sequences);
        for (i=0;i<dimn;i++) {
            printf ("Enter %lld Sequence File:\n", i);
            seqName[i] = NULL;
            seqName[i] = mmalloc (((MOATypeInd) (LINE_MAX)) * ((MOATypeInd) sizeof *seqName[i]));
            if (seqName[i] == NULL) {
                printf ("Error Reading File. Exiting\n");
                break;
            }
            scanf("%s", &sfilename);
            strcpy(seqName[i], sfilename);
            shape[i] = readInputSequence(&seqName[i], &sequences[i]);
        }
       
        /*freeProcessMemory (&pData, &sData, &wData);
        printf ("Enter Fasta Sequences File:\n", i);
        scanf("%s", &sfilename);
        dimn = readSequencesFastaFile(sfilename, &sequences, &seqName, &shape);
         */
        ret = initProcessMemory(&pData, &sData, &wData, dimn, shape, sequences, seqName, ret, 3);  
            if (ret != 0) {
                mprintf (0, "Error Initializing Process Data, Exiting\n", 1);
                fflush (stdout);
                return -1;
            }            
        pData->msaAlgn =  NULL;
        createMOAStruct(&pData->msaAlgn);
        if (createMOA(pData->seqLen, pData->seqNum, pData->msaAlgn, 0, 0) < 0)
            break;		
        wData->wavesTotal = 1;
        wData->AllpartsInWave =  mmalloc((MOATypeInd) sizeof *wData->AllpartsInWave);
        wData->AllpartsInWave[0] = 1;
        pData->waveNo =  0;
        pData->partNo = 0;
        pData->partitionsCount = 1;
        /*pData->OCout = mmalloc((MOATypeInd) sizeof *(pData->OCout));
        if (pData->OCout == NULL) {
            mprintf(1, "Couldn't create memory for OCout. Exiting.\n", 3);
            printf("Couldn't create memory for OCout. Exiting.\n");			
            return -1;
        }
        pData->OCout[0].wavesOC = 0;
        pData->OCout[0].WOCO = NULL;*/
        /* Compute Scored for Current Partition*/
        DPComputeScores (pData, sData, wData);
        //printMOA_scr(pData->msaAlgn, 0);
        //pData->computedPartitions ++;
        //checkPoint (pData, sData);		
        break;
    case 27:
        if (tbData != NULL) {
            freetbData(&tbData);
         }
        tbData = NULL;
        if (initTBData(&tbData, pData->seqNum, pData->seqLen) != 0) {
            printf ("Failed to allocate memory for trace back structure. Exiting\n");
            return -1;
        }
        
        tbData->aSeqLen = mcalloc ((MOATypeInd) 1, (MOATypeInd) sizeof *(tbData->aSeqLen));
        tbData->algnseq = mmalloc ((MOATypeInd) sizeof *(tbData->algnseq));
        tbData->pathParts = 1;
        tbData->algnseq[0] = NULL;
        getLocalMaxCellScore (pData, wData, tbData, 0);
        tbData->aSeqLen = mcalloc ((MOATypeInd) 1, (MOATypeInd) sizeof *(tbData->aSeqLen));
        tbData->aSeqLen[0] = 0;
        tbData->algnseq = mmalloc ((MOATypeInd) sizeof *(tbData->algnseq));
        tbData->algnseq[0] = NULL;
        tbData->algnseq[0] = mmalloc (pData->seqNum * sizeof *(tbData->algnseq[0]));    
        tbData->pathParts = 1;
        traceBack (pData, wData, tbData);
        tbData->comseqLen = tbData->aSeqLen[0];
        tbData->completePath = tbData->algnseq[0];
        calcAlignmentSPScore (pData, tbData);
        outputAlignment(pData, tbData, 0);
        editAlignment (pData, tbData);
        calcAlignmentSPScore (pData, tbData);
	outputAlignment(pData, tbData, 0);
        calcAlignmentSPScore (pData, tbData);
        break;
    case 28:
        /*Read Sequences Data*/
        printf ("Enter Scoring Type\n");
        scanf("%d", &ret);
        printf ("Enter Number of Sequences\n");
        scanf("%lld", &dimn);
        pData->seqName = mmalloc (((MOATypeInd) pData->seqNum) * ((MOATypeInd) sizeof *pData->seqName));
        pData->seqLen = mcalloc ((MOATypeInd) pData->seqNum, (MOATypeInd) sizeof *pData->seqLen);
        pData->sequences = mmalloc (((MOATypeInd) pData->seqNum) * (MOATypeInd) sizeof *pData->sequences);

        for (i=0;i<pData->seqNum;i++) {
            printf ("Enter %lld Sequence File:\n", i);
            pData->seqName[i] = NULL;
            pData->seqName[i] = mmalloc (((MOATypeInd) (LINE_MAX)) * ((MOATypeInd) sizeof *pData->seqName[i]));
            if (pData->seqName[i] == NULL) {
                printf ("Error Reading File. Exiting\n");
                break;
            }
            scanf("%s", &sfilename);
            strcpy(pData->seqName[i], sfilename);
            pData->seqLen[i] = readInputSequence(&pData->seqName[i], &pData->sequences[i]);
        }
        printf ("\nEnter partition size:");
        scanf("%lld", &flatIndex);
        ret = initProcessMemory(&pData, &sData, &wData, dimn, pData->seqLen, pData->sequences, pData->seqName, pData->stype, (long) flatIndex);  
        /*Partition Waves*/
        Mode = Distributed;
        printf ("\nEnter The Size of the Cluster:");
        scanf("%d", &ClusterSize);
        printf ("\nEnter myProcid:");
        scanf("%d", &myProcid); 
        wData->wavesTotal = calcWaves (pData, wData);
        /*Score one Partition*/
        pData->msaAlgn =  NULL;
        pData->waveNo = 0;
        pData->partNo = 0;
        getPartition (wData->partsInWaveIndices[pData->waveNo][pData->partNo], pData->seqNum, pData->seqLen, &pData->msaAlgn, wData->partitionSize);
        
        
        /* Compute Scored for Current Partition*/
        DPComputeScores (pData, sData, wData);
        printMOA_scr(pData->msaAlgn, 0);
        prepareWaveOC (0, pData);   
        pData->computedPartitions ++;
        break;
    case 29:
        freeProcessMemory (&pData, &sData, &wData);
        if (tbData != NULL)
            freetbData(&tbData);
        pData =  mmalloc ((MOATypeInd) sizeof *pData);
        if (pData == NULL)
            return -1;

        pData =  mmalloc ((MOATypeInd) sizeof *pData);
        if (pData == NULL)
            return -1;

        tbData = NULL;
        tbData =  mmalloc ((MOATypeInd) sizeof *tbData);
        if (tbData == NULL)
            return -1;
	printf ("Enter Fasta Alignemt File:\n");
	scanf("%s", &sfilename);
        pData->stype = 0;
	readFastaAlignment (sfilename, pData, tbData);	
        calcAlignmentSPScore (pData, tbData);
	outputAlignment(pData, tbData, 0);
        break;
    case 30:
        printf ("Enter the Dimension\n");
        scanf("%lld", &dimn);
        ind = mmalloc (dimn * sizeof *shape);
        for (i=0;i<dimn;i++) {
            printf ("Enter the vectors Index[%lld]\n", i);
            scanf("%lld", &ind[i]);
        }
        printf ("Pythagorean Distance from Origin = %f\n", (double) getPythagoreanOriginDistance (ind, dimn));
	if (ind != NULL)
            free(ind);
        break;
    case 31:
        printf ("Enter the Dimension\n");
        scanf("%lld", &dimn);
        ind = mmalloc (dimn * sizeof *shape);
        ind2 = mmalloc (dimn * sizeof *shape);
        for (i=0;i<dimn;i++) {
            printf ("Enter the X vectors Index[%lld]\n", i);
            scanf("%lld", &ind[i]);
        }
        for (i=0;i<dimn;i++) {
            printf ("Enter the Y vectors Index[%lld]\n", i);
            scanf("%lld", &ind2[i]);
        }
        printf ("Pythagorean Distance between the 2 vectors = %f\n", (double) getPythagoreanDistance (ind, ind2, dimn));
	if (ind != NULL)
            free(ind);
	if (ind2 != NULL)
            free(ind2);
        break;
    case 32:
        printf ("Enter the Dimension\n");
        scanf("%lld", &dimn);
        ind = mmalloc (dimn * sizeof *shape);
        for (i=0;i<dimn;i++) {
            printf ("Enter the vectors Index[%lld]\n", i);
            scanf("%lld", &ind[i]);
        }
        printf ("Pythagorean Origin Direction of the vector = %f\n", (double) getPythagoreanOriginDirection (ind, dimn));
	if (ind != NULL)
            free(ind);
        break;
    case 33:
        printf ("Enter the Dimension\n");
        scanf("%lld", &dimn);
        ind = mmalloc (dimn * sizeof *shape);
        ind2 = mmalloc (dimn * sizeof *shape);
        for (i=0;i<dimn;i++) {
            printf ("Enter the X vectors Index[%lld]\n", i);
            scanf("%lld", &ind[i]);
        }
        for (i=0;i<dimn;i++) {
            printf ("Enter the Y vectors Index[%lld]\n", i);
            scanf("%lld", &ind2[i]);
        }
        printf ("Pythagorean Direction Difference between the 2 vectors = %f\n", (double) getPythagoreanDirection (ind, ind2, dimn));
	if (ind != NULL)
            free(ind);
	if (ind2 != NULL)
            free(ind2);
        break;
    case 34:
        printf ("Enter the Dimension\n");
        scanf("%lld", &dimn);
        ind = mmalloc (dimn * sizeof *shape);
        ind2 = mmalloc (dimn * sizeof *shape);
        for (i=0;i<dimn;i++) {
            printf ("Enter the X vectors Index[%lld]\n", i);
            scanf("%lld", &ind[i]);
        }
        for (i=0;i<dimn;i++) {
            printf ("Enter the Y vectors Index[%lld]\n", i);
            scanf("%lld", &ind2[i]);
        }
        printf ("MoA Direction Difference between the 2 vectors = %f\n", (double) getMoADistance (ind, ind2, dimn));
	if (ind != NULL)
            free(ind);
	if (ind2 != NULL)
            free(ind2);
        break;
    case 35:
        printf ("Enter the Dimension\n");
        scanf("%lld", &dimn);
        ind = mmalloc (dimn * sizeof *shape);
        ind2 = mmalloc (dimn * sizeof *shape);
        for (i=0;i<dimn;i++) {
            printf ("Enter the X vectors Index[%lld]\n", i);
            scanf("%lld", &ind[i]);
        }
        for (i=0;i<dimn;i++) {
            printf ("Enter the Y vectors Index[%lld]\n", i);
            scanf("%lld", &ind2[i]);
        }
        printf ("Angle Difference between the 2 vectors = %f\n", (double) getVectorsAngle (ind, ind2, dimn));
	if (ind != NULL)
            free(ind);
	if (ind2 != NULL)
            free(ind2);
        break;
    case 36:
        //printf ("Enter the Dimension\n");
        //scanf("%lld", &dimn);
        dimn = 6;
        strides = 3;
        //printf ("Enter the Wave No\n");
        //scanf("%lld", &wave);
        wave = 8;
        ind = mmalloc (dimn * sizeof *shape);
        ind2 = mmalloc (dimn * sizeof *shape);
        shape = mmalloc (dimn * sizeof *shape);
        /*for (i=0;i<dimn;i++) {
            printf ("Enter the X vectors Index[%lld]\n", i);
            scanf("%lld", &ind[i]);
        }*/
        ind[0] = 1;
        ind[1] = 1;
        ind[2] = 1;
        ind[3] = 1;
        ind[4] = 1;
        ind[5] = 1;
        shape[0] = 11;
        shape[1] = 11;
        shape[2] = 11;
        shape[3] = 11;
        shape[4] = 11;
        shape[5] = 11;
        //getFirstPartitionInWave (shape, dimn, strides, wave, &ind);
        /*for (i=0;i<dimn;i++) {
            printf ("Enter the Y vectors Index[%lld]\n", i);
            scanf("%lld", &ind2[i]);
        }*/
         FILE * f = fopen ("/export/home/mhelal1/thesis/exp/SunStudioProjects/moaLib/list.txt", "r");
         if (f == NULL) {
            printf("File could not be opened!\n");
            return -1;
        }
        char elm = fgetc(f);
        printf ("Ranked Angles From Index {%lld", ind[0]);
        for (i=1;i<dimn;i++) 
            printf (", %lld", ind[i]);
        printf ("} to rest of indices in Wave\nAngle\n");
        while (elm != EOF) {
            i = 0;
            while (i<dimn) {
                if ((elm != ' ') && (elm != '\n') ){
                    ind2[i] = atoll (&elm);
                    i ++;
                }
                elm = fgetc(f);
            } 
            printf ("%f\n", (double) getVectorsRankedAngle (ind2, ind, dimn));
        }
	if (ind != NULL)
            free(ind);
	if (ind2 != NULL)
            free(ind2);
	if (shape != NULL)
            free(shape);
        fclose(f);
        break;
    case 37:
        printf ("Enter the Dimension\n");
        scanf("%lld", &dimn);
        ind = mmalloc (dimn * sizeof *shape);
        ind2 = mmalloc (dimn * sizeof *shape);
        shape = mmalloc (dimn * sizeof *shape);
        for (i=0;i<dimn;i++) {
            ind[i] = 0;
        }
        for (i=0;i<dimn;i++) {
            printf ("Enter the Y vectors Index[%lld]\n", i);
            scanf("%lld", &ind2[i]);
        }
        for (i=0;i<dimn;i++) {
            printf ("Enter the Z vectors Index[%lld]\n", i);
            scanf("%lld", &shape[i]);
        }
        printf ("Is Vector in HyperPlane = %lld\n",  isVectorInHyperPlane (ind, ind2, dimn, shape));
	if (ind != NULL)
            free(ind);
	if (ind2 != NULL)
            free(ind2);
	if (shape != NULL)
            free(shape);
        break;
    case 38:
        printf ("Enter the Dimension\n");
        scanf("%lld", &dimn);
        shape = mmalloc (dimn * sizeof *shape);
        for (i=0;i<dimn;i++) {
            printf ("Enter the shape[%lld]\n", i);
            scanf("%lld", &shape[i]);
        }
        printf ("Enter the partitionSize\n");
        scanf("%lld", &strides);
        printf ("Enter the Wave No\n");
        scanf("%lld", &wave);
        printf ("The Wave Length = %f\n", (double) getWaveLength (shape, dimn, strides, wave));
	if (shape != NULL)
            free(shape);
        break;
    case 39:
        printf ("Enter the Dimension\n");
        scanf("%lld", &dimn);
        shape = mmalloc (dimn * sizeof *shape);
        ind = mmalloc (dimn * sizeof *shape);
        ind2 = mmalloc (dimn * sizeof *shape);
        for (i=0;i<dimn;i++) {
            printf ("Enter the Shape [%lld]\n", i);
            scanf("%lld", &shape[i]);
        }
        for (i=0;i<dimn;i++) {
            printf ("Enter the Partition Index [%lld]\n", i);
            scanf("%lld", &ind2[i]);
        }
        printf ("Enter the partitionSize\n");
        scanf("%lld", &strides);
        printf ("Enter the Wave No\n");
        scanf("%lld", &wave);
        getMiddlePartitionInWave (shape, dimn, strides, wave, &ind) ;
        
        printf ("The Wave Length = %f\n", (double) getMoADistance (ind2, ind, dimn));
	if (shape != NULL)
            free(shape);
	if (ind != NULL)
            free(ind);
	if (ind2 != NULL)
            free(ind2);
        break;
    case 40:
	printf ("Enter Path to files:\n");
	scanf("%s", &sfilename);
        readAndConvertToFasta (sfilename);
        break;
    case 41:
	printf ("Enter wave No:\n");
	scanf("%ld", &pData->waveNo);
        ind = mmalloc (MOA1->dimn * sizeof *ind);
        getMiddlePartitionInWave (MOA1->shape, MOA1->dimn, pData->partitionSize, pData->waveNo, &ind);
        printf ("The middle partition is: {%lld", ind[0]);
        for (i=1;i<MOA1->dimn;i++) 
            printf (", %lld", ind[i]);
        printf ("}\n");

        if (ind != NULL)
            free(ind);
        break;
        case 42:
            printf ("Enter wave No:\n");
            scanf("%ld", &i);
            sequences =  mmalloc (sizeof* sequences);
            sequences[0] = mmalloc ((i+1) * sizeof* sequences[0]);
            generateRandomSequence ((MOATypeShape) i, &sequences[0]);
            printf ("Randomly Generate Sequence is %s", sequences[0]);
            free (sequences[0]);
            free (sequences);
            break;
    case 43:
        freeProcessMemory (&pData, &sData, &wData);
        if (tbData != NULL)
            freetbData(&tbData);
        pData =  mmalloc ((MOATypeInd) sizeof *pData);
        if (pData == NULL)
            return -1;

        pData =  mmalloc ((MOATypeInd) sizeof *pData);
        if (pData == NULL)
            return -1;

        tbData = NULL;
        tbData =  mmalloc ((MOATypeInd) sizeof *tbData);
        if (tbData == NULL)
            return -1;
	printf ("Enter Fasta Alignemt input File:\n");
	scanf("%s", &sfilename);
        printf ("Enter MSF Alignemt output File:\n");
	scanf("%s", &sfilename_2);
        pData->stype = 0;
	readFastaAlignment (sfilename, pData, tbData);	
        outputMSFAlignment (sfilename_2, pData, tbData);
        outputMSFAlignment_2 (sfilename_2, pData, tbData);
        break;
    }
    } while (choice != 0);
    if (posDimn != NULL)
        free(posDimn);
    if (ind != NULL)
      free(ind);
    if (ind2 != NULL)
        free(ind2);
    if (shape != NULL)
        free (shape);
    freeProcessMemory (&pData, &sData, &wData);
    if (tbData != NULL)
        freetbData(&tbData);

    deleteMOA (MOA1);  
}
