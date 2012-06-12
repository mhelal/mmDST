#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <mpi.h>
#include <malloc.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <errno.h>
#include "../main.h"
#include "../globals.h"
#include "../utils.h"
#include "../moaDst.h"
#define _GNU_SOURCE



void assemblePathParts (TracebackData * tbData) {
    int ret;
    MOATypeInd i, j, k;
    char  msg[MID_MESSAGE_SIZE];

    tbData->completePath = mmalloc (tbData->seqNum * sizeof *(tbData->completePath));    
    tbData->comseqLen = 0;

    for (i=0;i<tbData->seqNum;i++) {
        tbData->comseqLen = 0;
        tbData->completePath[i] = NULL;
        for (j = tbData->pathParts-1;j>=0;j--) {
            tbData->completePath[i] = realloc (tbData->completePath[i], (tbData->comseqLen + tbData->aSeqLen[j]) * sizeof *(tbData->completePath[i]) );
            sprintf(msg, "completePath[%ld] has %ld and will realloc %ld \n", i, tbData->comseqLen, (tbData->comseqLen + tbData->aSeqLen[j]));  
            mprintf(20, msg, 1);
            for (k = 0;k<tbData->aSeqLen[j];k++) {
                tbData->completePath[i][tbData->comseqLen] = tbData->algnseq[j][i][k];
                tbData->comseqLen ++;
            }
        }
    }
}

void * tbMaster (ProcessData * pData, WavesData * wData) {
    int MPI_return, foundproc, done = 0; /*currProc*/
    MOATypeDimn i, k;
    long long receviedLongLong;
    
    char msg[MID_MESSAGE_SIZE];
    MPI_Request request;
    MPI_Status status;
    TracebackData * tbData;
    tbData = NULL;
    if (initTBData(&tbData, pData->seqNum, pData->seqLen) != 0) {
        printf ("Failed to allocate memory for trace back structure. Exiting\n");
        return NULL;
    }
    if (AlignmentType == Local) 
        getmaxCellScore (pData, tbData);    
    else {
        //tbData->maxCellScore = getLocalMaxCellScore(pData, wData, tbData, 0);
        pData->waveNo = wData->wavesTotal - 1;
        pData->partNo = wData->partsInWave[pData->waveNo]-1;
        tbData->maxCellScore = 0;
        if (getProcID(wData, pData->waveNo, pData->partNo) != myProcid) 
            getPrevPartition (wData, &pData->waveNo, &pData->partNo);
        if (pData->partNo != -1) {
            if (restorePartitionCheckPoint(pData, wData, pData->waveNo, pData->partNo) != 0) {
                printf ("Error Retrieving partition file. Exitiing.\n");
                return NULL;
            }
            for (k=0;k<pData->seqNum;k++)
                tbData->maxCellIndex[k] = pData->msaAlgn->indexes[pData->msaAlgn->elements_ub-1][k];
        }
        else {
            printf ("Error Retrieving part No. Exitiing.\n");
            return NULL;
        }
        
        tbData->currProc = 0;
    }
    
    //currProc = tbData->currProc;
    while (done == 0) {
        printf ("currProc = %d done = %d maxCellIndex { %lld", tbData->currProc, done, tbData->maxCellIndex[0]);
        for (k=1;k<tbData->seqNum;k++)
            printf (", %lld", tbData->maxCellIndex[k]);
        printf ("} in proc %d\n", tbData->currProc);
        /*send to processor containing the maxCellScore to trace back*/
        if (tbData->currProc == 0) {
                /*Perform Partition trace back*/
                tbData->pathParts ++;
                if (tbData->pathParts == 1) {
                    tbData->aSeqLen = mmalloc (((MOATypeInd) sizeof *(tbData->aSeqLen)));
                    if (tbData->aSeqLen == NULL) {
                        printf ("Failed to reallocate memory for %ld partial alignments Lengths. Exiting.\n", tbData->pathParts);
                        fflush (stdout);
                        return NULL;
                    }
                    tbData->algnseq = mmalloc (((MOATypeInd) sizeof *(tbData->algnseq)));
                    if (tbData->aSeqLen == NULL) {
                        printf ("Failed to reallocate memory for %ld partial alignments. Exiting.\n", tbData->pathParts);
                        fflush (stdout);
                        return NULL;
                    }
                }
                else {
                    tbData->aSeqLen = realloc (tbData->aSeqLen, ((MOATypeInd) tbData->pathParts) * ((MOATypeInd) sizeof *(tbData->aSeqLen)));
                    if (tbData->aSeqLen == NULL) {
                        printf ("Failed to reallocate memory for %ld partial alignments Lengths. Exiting.\n", tbData->pathParts);
                        fflush (stdout);
                        return NULL;
                    }
                    tbData->algnseq = realloc (tbData->algnseq, ((MOATypeInd) tbData->pathParts) * ((MOATypeInd) sizeof *(tbData->algnseq)));
                    if (tbData->algnseq == NULL) {
                        printf ("Failed to reallocate memory for %ld partial alignments. Exiting.\n", tbData->pathParts);
                        fflush (stdout);
                        return NULL;
                    }
                }
                tbData->algnseq[tbData->pathParts-1] = NULL;
                tbData->algnseq[tbData->pathParts-1] = mmalloc (((MOATypeInd) tbData->seqNum) * ((MOATypeInd) sizeof *(tbData->algnseq[tbData->pathParts-1])));    
                if (tbData->algnseq[tbData->pathParts-1] == NULL) {
                    printf ("Failed to reallocate memory for %ld partial alignments sequences. Exiting.\n", tbData->pathParts);
                    fflush (stdout);
                    return NULL;
                }
                tbData->aSeqLen[tbData->pathParts-1] = traceBack(pData, wData, tbData);
                sprintf(msg, "DMTB returned Local path of length %ld\n", tbData->aSeqLen[tbData->pathParts-1]);  
                mprintf(3, msg, 1);
        }
        else {
                /*1. Send Tracing flag (done = 0) to the processor where the maximum score was found*/
                MPI_return = MPI_Send (&done, 1, MPI_INT, tbData->currProc, 2, MOAMSA_COMM_WORLD);
#ifndef NDEBUG
                sprintf(msg, "DMTB sent flag %d to proc %d first\n", done, tbData->currProc);  
                mprintf(3, msg, 1);
#endif
                /*2. send starting global index*/
                MPI_return = MPI_Send (tbData->maxCellIndex, tbData->seqNum, MPI_LONG_LONG, tbData->currProc, 3, MOAMSA_COMM_WORLD);
#ifndef NDEBUG
                sprintf(msg, "DMTB sent maxCellIndex to proc %d\n", tbData->currProc);  
                mprintf(3, msg, 1);
#endif


                tbData->pathParts ++;
                if (tbData->pathParts == 1) {
                    tbData->aSeqLen = mmalloc (((MOATypeInd) sizeof *tbData->aSeqLen));
                    if (tbData->aSeqLen == NULL) {
                        printf ("Failed to reallocate memory for %ld partial alignments Lengths. Exiting.\n", tbData->pathParts);
                        fflush (stdout);
                        return NULL;
                    }
                    tbData->algnseq = mmalloc (((MOATypeInd) sizeof *tbData->algnseq));
                    if (tbData->aSeqLen == NULL) {
                        printf ("Failed to reallocate memory for %ld partial alignments. Exiting.\n", tbData->pathParts);
                        fflush (stdout);
                        return NULL;
                    }
                }
                else {
                    tbData->aSeqLen = realloc (tbData->aSeqLen, ((MOATypeInd) tbData->pathParts) * ((MOATypeInd) sizeof *tbData->aSeqLen));
                    if (tbData->aSeqLen == NULL) {
                        printf ("Failed to reallocate memory for %ld partial alignments Lengths. Exiting.\n", tbData->pathParts);
                        fflush (stdout);
                        return NULL;
                    }
                    tbData->algnseq = realloc (tbData->algnseq, ((MOATypeInd) tbData->pathParts) * ((MOATypeInd) sizeof *tbData->algnseq));
                    if (tbData->aSeqLen == NULL) {
                        printf ("Failed to reallocate memory for %ld partial alignments. Exiting.\n", tbData->pathParts);
                        fflush (stdout);
                        return NULL;
                    }
                }
                tbData->algnseq[tbData->pathParts-1] = NULL;
                tbData->algnseq[tbData->pathParts-1] = mmalloc (((MOATypeInd) tbData->seqNum) * ((MOATypeInd) sizeof *(tbData->algnseq[tbData->pathParts-1])));  
                if (tbData->algnseq[tbData->pathParts-1] == NULL) {
                    printf ("Failed to reallocate memory for %ld partial alignments sequences. Exiting.\n", tbData->pathParts);
                    fflush (stdout);
                    return NULL;
                }
                tbData->aSeqLen[tbData->pathParts-1] = 2;
                /*3. receive partital alignment length*/
                MPI_return = MPI_Recv (&receviedLongLong, 1, MPI_LONG_LONG, tbData->currProc, 4, MOAMSA_COMM_WORLD, &status);
                tbData->aSeqLen[tbData->pathParts-1] = receviedLongLong;
#ifndef NDEBUG
                printf("Master received length %ld  MPI_return %d \n", tbData->aSeqLen[tbData->pathParts-1], MPI_return);  
                fflush(stdout);
                sprintf(msg, "DMTB received Remote path of length %ld :", tbData->aSeqLen[tbData->pathParts-1]);  
                mprintf(3, msg, 1);
#endif
                /*4. receive the partital alignment itself*/
                for (i=0;i<tbData->seqNum && tbData->aSeqLen[tbData->pathParts-1] > 0;i++) {
                    tbData->algnseq[tbData->pathParts-1][i] = NULL;
                    tbData->algnseq[tbData->pathParts-1][i] = mmalloc (((MOATypeInd) (tbData->aSeqLen[tbData->pathParts-1] +1)) * ((MOATypeInd) sizeof *(tbData->algnseq[tbData->pathParts-1][i] )));    
                    if (tbData->algnseq[tbData->pathParts-1][i] == NULL) {
                        printf ("Failed to reallocate memory for %ld partial alignments sequences %lld residues. Exiting.\n", tbData->pathParts, tbData->aSeqLen[tbData->pathParts-1]);
                        fflush (stdout);
                        return NULL;
                    }
                    MPI_return = MPI_Recv (tbData->algnseq[tbData->pathParts-1][i], tbData->aSeqLen[tbData->pathParts-1], MPI_CHAR, tbData->currProc, 5, MOAMSA_COMM_WORLD, &status);
                    tbData->algnseq[tbData->pathParts-1][i][tbData->aSeqLen[tbData->pathParts-1]] = '\0';
#ifndef NDEBUG
                    printf("Master received aligned seq for seq %lld MPI_return %d = %s\n", i, MPI_return, tbData->algnseq[tbData->pathParts-1][i]);  
                    fflush(stdout);
                    sprintf(msg, " %s ", tbData->algnseq[tbData->pathParts-1][i]);  
                    mprintf(3, msg, 1);
#endif
                }
                /*5. receive the last global index in this partial alignment*/
                MPI_return = MPI_Recv (tbData->maxCellIndex, tbData->seqNum, MPI_LONG_LONG, tbData->currProc, 7, MOAMSA_COMM_WORLD, &status);
#ifndef NDEBUG
                printf ("Master received MPI_return %d maxCellIndex { %lld", MPI_return, tbData->maxCellIndex[0]);
                for (i=1;i<tbData->seqNum;i++)
                    printf (", %lld", tbData->maxCellIndex[i]);
                printf ("}\n");
                fflush(stdout);
#endif
                /*6. receive the next Processor where a next remote score was found from this partial alignment*/
                MPI_return = MPI_Recv (&tbData->currProc, 1, MPI_INT, tbData->currProc, 8, MOAMSA_COMM_WORLD, &status);
                //tbData->currProc = currProc;
#ifndef NDEBUG
                sprintf(msg, "DMTB received currProc %d, maxCellIndex {%lld", tbData->currProc, tbData->maxCellIndex[0]);  
                for (k=1;k<tbData->seqNum;k++)
                    sprintf(msg, "%s, %lld", msg, tbData->maxCellIndex[k]);
                sprintf(msg, "%s}\n", msg);
                mprintf(3, msg, 1);
#endif
            }
            /*test if end of aligned sequence is zero to exit */
            done = 2;
            for (k=0;k<tbData->seqNum;k++)
                if (tbData->maxCellIndex[k] != 0)
                    done = 0;
            /* else, determine the next tracing Processor*/
            if ((done == 0) && (tbData->currProc < 0)) { 
                /*if other slave processes don't have this index, check the master*/
                if (getPartitionDetails (pData, wData, &tbData->partIndex, tbData->maxCellIndex, &pData->waveNo, &pData->partNo, &tbData->currProc) == 0) { 
#ifndef NDEBUG
                    sprintf(msg, "DMTB max cell {%lld", tbData->maxCellIndex[0]);  
                    for (k=1;k<tbData->seqNum;k++)
                        sprintf (msg, "%s, %lld", msg, tbData->maxCellIndex[k]);
                    sprintf (msg, "%s} in proc %d\n", msg, tbData->currProc);
                    mprintf(3, msg, 1);
#endif
                    done = 0;
                }
                else /*other wise, end of tracing*/
                    done = 2;
            }
    } /* End While*/
    /*send to all processes to exit tracing*/
    done = 2;
    for (i=1;i<ClusterSize;i++) {
        MPI_return = MPI_Send (&done, 1, MPI_INT, i, 2, MOAMSA_COMM_WORLD);
        sprintf(msg, "DMTB sent flag %d to proc %ld\n", done, i);  
        mprintf(3, msg, 1);
    }
    assemblePathParts (tbData);
    calcAlignmentSPScore(pData, tbData);
    outputAlignment(pData, tbData, 1);
    editAlignment (pData, tbData);
    calcAlignmentSPScore (pData, tbData);
    outputAlignment(pData, tbData, 1);
    outputFastaAlignment(pData, tbData);    
    outputMSFAlignment(outputfilename, pData, tbData);    
    if (tbData != NULL)
        freetbData(&tbData);
    return NULL;
}

void * tbSlave (ProcessData * pData, WavesData * wData) {
    MOATypeInd localCellIndex, maxlocalCellIndex;
    MOATypeDimn i;
    int MPI_return, isLower, done = 0;
    MPI_Request request;
    MPI_Status status;
    char msg[MID_MESSAGE_SIZE];
    TracebackData * tbData;

    tbData = NULL;
    if (initTBData(&tbData, pData->seqNum, pData->seqLen) != 0) {
        printf ("Failed to allocate memory for trace back structure. Exiting\n");
        return NULL;
    }

    if (AlignmentType == Local) 
        sendmaxCellScore (pData, wData, tbData);
    

    tbData->aSeqLen = mcalloc ((MOATypeInd) 1, (MOATypeInd) sizeof *(tbData->aSeqLen));
    tbData->algnseq = mmalloc ((MOATypeInd) sizeof *(tbData->algnseq));
    tbData->pathParts = 1;
    tbData->algnseq[0] = NULL;
    tbData->algnseq[0] = mmalloc (tbData->seqNum * sizeof *(tbData->algnseq[0]));    
    while (done < 2) {
        /*receive tracing flag, 0: trace back, Otherwise: finish & exit*/
#ifndef NDEBUG
        mprintf (3, "before receive new flag from master", 1);
#endif
        /*1. Receive flag */
        MPI_return = MPI_Recv (&done, 1, MPI_INT, 0, 2, MOAMSA_COMM_WORLD, &status);
#ifndef NDEBUG
        sprintf(msg, "DSTB received flag %d \n", done);  
        mprintf(3, msg, 1);
#endif

        /*Check tracing flag (done = 0)*/
        if (done == 0) {
            /*2. receive starting global index*/
            MPI_return = MPI_Recv (tbData->maxCellIndex, tbData->seqNum, MPI_LONG_LONG, 0, 3, MOAMSA_COMM_WORLD, &status);
#ifndef NDEBUG
            printf ("[%d] Received maxCellIndex { %lld", myProcid, tbData->maxCellIndex[0]);
            for (i=1;i<tbData->seqNum;i++)
                printf (", %lld", tbData->maxCellIndex[i]);
            printf ("}\n");
            fflush(stdout);
            sprintf(msg, "DSTB received maxCellIndex {%lld", tbData->maxCellIndex[0]);  
            for (i=1;i<tbData->seqNum;i++) 
                sprintf (msg, "%s, %lld", msg, tbData->maxCellIndex[i]);
            sprintf (msg, "}\n", msg);
            mprintf(3, msg, 1);
#endif
            /*Perform local Partitions trace back*/
            tbData->aSeqLen[0] = traceBack (pData, wData, tbData);
#ifndef NDEBUG
            sprintf(msg, "DSTB returned path of length %ld \n", tbData->aSeqLen[0]);  
            mprintf(3, msg, 1);
#endif

            /*3. send partital alignment length*/	
            MPI_return = MPI_Send (&tbData->aSeqLen[0], 1, MPI_LONG_LONG, 0, 4, MOAMSA_COMM_WORLD);
#ifndef NDEBUG
            printf("[%d]  sent length %ld  MPI_return %d \n", myProcid, tbData->aSeqLen[0], MPI_return);  
            fflush(stdout);
#endif
            /*4. send the partital alignment itself*/
            for (i=0;i<tbData->seqNum;i++) {
                tbData->algnseq[0][i][tbData->aSeqLen[0]] = '\0';
                MPI_return = MPI_Send (tbData->algnseq[0][i], tbData->aSeqLen[0], MPI_CHAR, 0, 5, MOAMSA_COMM_WORLD);
#ifndef NDEBUG
                sprintf(msg, " sent aligned seq for seq %lld MPI_return %d = %s ", i, MPI_return, tbData->algnseq[0][i]);  
                mprintf(3, msg, 1);
                printf("[%d]  sent aligned seq for seq %lld MPI_return %d = %s\n", myProcid, i, MPI_return, tbData->algnseq[0][i]);  
                fflush(stdout);
#endif
            }
#ifndef NDEBUG
            mprintf(3, "\n", 1);
#endif
            /*5. send the last global index in this partial alignment to Master*/		
            MPI_return = MPI_Send (tbData->maxCellIndex, tbData->seqNum, MPI_LONG_LONG, 0, 7, MOAMSA_COMM_WORLD);
#ifndef NDEBUG
            printf ("[%d] sent MPI_return %d maxCellIndex { %lld", myProcid, MPI_return, tbData->maxCellIndex[0]);
            for (i=1;i<tbData->seqNum;i++)
                printf (", %lld", tbData->maxCellIndex[i]);
            printf ("}\n");
            fflush(stdout);
            sprintf(msg, "after send maxCellIndex {%ld", tbData->maxCellIndex[0]);
            mprintf(3, msg, 1);
            for (i=1;i<tbData->seqNum;i++) {
                sprintf(msg, ", %ld", tbData->maxCellIndex[i]);
                mprintf(3, msg, 1);
            }
            mprintf(3, "}\n", 1);
#endif
            /*6. send the next Processor where a next remote score was found from this partial alignment to Master*/
            MPI_return = MPI_Send (&tbData->currProc, 1, MPI_INT, 0, 8, MOAMSA_COMM_WORLD);
#ifndef NDEBUG
            printf("[%d]  sent new  tbData->currProc = %d MPI_return %d \n", myProcid, tbData->currProc, MPI_return);  
            fflush(stdout);
            sprintf(msg, "after send next proc %d\n", tbData->currProc);
            mprintf(3, msg, 1);
#endif
        }
    }
    if (tbData != NULL)
        freetbData(&tbData);
    return NULL;
}

void getmaxCellScore (ProcessData * pData, TracebackData * tbData) {
    int i;
    MOATypeElmVal recvMaxCellScore;
    MOATypeShape * recvMaxCellIndex = NULL;
    char msg[MID_MESSAGE_SIZE];
    MOATypeDimn k;

    MPI_Status status;
    MPI_Request request;
    recvMaxCellIndex =  mcalloc ((MOATypeInd) tbData->seqNum, (MOATypeInd) sizeof *recvMaxCellIndex);
    if (recvMaxCellIndex == NULL) {
        printf ("Failed to allocate memory for the multidimensional Index. Exiting.\n");
        return;
    }
    for (i=1;i<ClusterSize;i++) {
        MPI_Recv (&recvMaxCellScore, 1, MPI_LONG_LONG, i, 0, MOAMSA_COMM_WORLD, &status);
        MPI_Recv (recvMaxCellIndex, pData->seqNum, MPI_LONG_LONG, i, 1, MOAMSA_COMM_WORLD, &status);
        sprintf (msg, "Master received recvMaxCellScore = %lld, recvMaxCellIndex = {%lld ", recvMaxCellScore, recvMaxCellIndex[0]);
        mprintf (3, msg, 1);
        for (k=1;k<pData->seqNum;k++) {
            sprintf (msg, ", %lld", recvMaxCellIndex[k]);
            mprintf (3, msg, 1);
        }
        sprintf  (msg, "} from %d\n", i);
        mprintf (3, msg, 1);
        if (i == 1) {
            tbData->maxCellScore = recvMaxCellScore;
            for (k=0;k<pData->seqNum;k++)
                tbData->maxCellIndex[k] = recvMaxCellIndex[k];
            tbData->currProc = i;
        }
        else if (recvMaxCellScore > tbData->maxCellScore) {
            tbData->maxCellScore = recvMaxCellScore;
            for (k=0;k<pData->seqNum;k++)
                tbData->maxCellIndex[k] = recvMaxCellIndex[k];
            tbData->currProc = i;
        }
    }
    sprintf (msg, "Master received max proc %d score  = %lld, Index = {%lld", tbData->currProc, tbData->maxCellScore, tbData->maxCellIndex[0]);
    for (k=1;k<pData->seqNum;k++)
        sprintf (msg, "%s, %lld", msg, tbData->maxCellIndex[k]);
    sprintf (msg, "%s}\n", msg);
    mprintf (3, msg, 1);
    if (recvMaxCellIndex != NULL)
        free (recvMaxCellIndex);
    recvMaxCellIndex = NULL;
}

void ExitProcess (ProcessData * pData, ScoringData * sData, WavesData * wData) {  

    freeProcessMemory (&pData, &sData, &wData);

    mprintf(1, "before finalize \n", 1);
    MPI_Finalize ();
    mprintf(1, "after finalize \n", 1);
}

int main (int argc, char * argv[]) {
    MPI_Group orig_group; 
    pthread_t MasterThread, SlaveThread;
    MOATypeDimn seqNum;
    MOATypeShape * seqLen = NULL;
    char * * sequences = NULL,  * * seqName = NULL, msg[MID_MESSAGE_SIZE];
    long partitionSize;
    long ID1;
    int stype, MPI_return;
    /*char ufilename[SHORT_MESSAGE_SIZE];*/
    int ret;
    ProcessData * pData;
    ScoringData * sData;
    WavesData * wData;
	
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myProcid);
    MPI_Comm_size (MPI_COMM_WORLD, &ClusterSize);
    MPI_Comm_group (MPI_COMM_WORLD, &orig_group);
    MPI_Comm_create(MPI_COMM_WORLD, orig_group, &MOAMSA_COMM_WORLD);
    TBFlag = 1;
    prevNow = NULL;
    prevNow = mmalloc ((MOATypeInd) sizeof *prevNow);
    currNow = NULL;
    //currNow = mmalloc ((MOATypeInd) sizeof *currNow);
    currNow = getTime();
    if (currNow == NULL) {
        printf ("Could not read current time. Exiting.\n");
        return -1;
    }
    ID1 = getpid();

    printf("[%d]>PID = %ld Started at time (%d, %d, %d, %d)\n", myProcid, ID1, currNow->tm_yday, currNow->tm_hour, currNow->tm_min, currNow->tm_sec);
    /* 1. Process Arguments*/
    processArguments(argc, argv, &seqNum, &sequences, &seqName, &seqLen, &stype, &partitionSize);
    /*strcpy (ufilename, outputfilename);*/
    strcpy (outPrefix, "b");
    /*sprintf (outputfilename, "mmtb%s", ufilename);*/
    init_output();
    sprintf (msg, "Program Arguments: debuglevel = %d maxAlignmentsNumber = %d Epsilons= %ld Alignment Type = %d stype %d outputfilename = %s partitionSize = %ld\n", pdebug, maxAlignmentsNumber, Epsilons, AlignmentType, stype, outputfilename, partitionSize);
    mprintf (1, msg, 1);
  
    /*2. Do the Master Tasks: Synchronize with Slaves to trace back*/
    /* one path for now*/
  
    ret = initProcessMemory(&pData, &sData, &wData, seqNum, seqLen, sequences, seqName, stype, partitionSize);  
    if (ret != 0) { 
        mprintf (1, " Could not initialize process memory\n", 1);
        ExitProcess (pData, sData, wData);
        return -1;
    }
    mprintf (1, "Initialized Slave data and will read checkpoint file\n", 1);
    /*3. Load Tensor Partitions Computed*/
    ret = restoreCheckPoint (pData, wData);
    if (ret != 0) {
        printf (" Could not read Slave data from checkpoint file\n");
        ExitProcess (pData, sData, wData);
        return -1;
    }
    if (pData->computedPartitions <= 0) {
        mprintf (1, " No Partitions read in this process\n", 1);
        /*ExitProcess (pData);*/
        /*return -1;*/
    }
    else {
        sprintf (msg, "read Slave data from checkpoint file partitions %ld last score in last partition is %ld sqm %ld sqlen0 %ld %ld %ld\n", pData->partitionsCount, pData->msaAlgn->elements[pData->msaAlgn->elements_ub - 1].val, pData->seqNum, pData->seqLen[0], pData->seqLen[1], pData->seqLen[2]);
        mprintf(3, msg, 1);
    }

    if (myProcid == 0) 
        tbMaster(pData, wData);
    else
    /* do the slave tasks*/
        tbSlave(pData, wData);
    if (myProcid == 0) {
        /* Getting Process Resources Usage ===================== */
        struct rusage usageRec;
        double utime, stime;

        ret = getrusage(RUSAGE_SELF, &usageRec);
        if (ret == 0) {
            //printf ("[%d]Resources Usage: UTime %ld, STime %ld, Mem %ld, Virt %ld\n", myProcid, usageRec.ru_utime.tv_sec, usageRec.ru_stime.tv_sec, usageRec.ru_maxrss, usageRec.ru_ixrss);
            utime = (double) usageRec.ru_utime.tv_sec + 1.e-6 * (double) usageRec.ru_utime.tv_usec;
            stime = (double) usageRec.ru_stime.tv_sec + 1.e-6 * (double) usageRec.ru_stime.tv_usec;	
            //printf ("[%d]Resources Usage: UTime %f, STime %f\n", myProcid, utime, stime);
        }
        else
            printf ("[%d]Failed to retrieve Process Resources Usage, errno %d\n", myProcid, errno);

        //struct mallinfo info;
        //info = mallinfo();

        //printf("[%d] STime\tUTime\theap\tMemory\t\n",myProcid);
        //printf("[%d] %f\t%f\t%d\t%d\n", myProcid, stime, utime, info.arena, info.usmblks + info.uordblks);

        printf("[%d] STime\tUTime\n",myProcid);
        printf("[%d] %f\t%f\n", myProcid, stime, utime);
        currNow = getTime();
        printf("[%d]>Finalized at time (%d, %d, %d, %d)\n", myProcid, currNow->tm_yday, currNow->tm_hour, currNow->tm_min, currNow->tm_sec);
    }
    ExitProcess (pData, sData, wData);
    return 0;
}
