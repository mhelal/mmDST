#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
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

void doPartitionTraceBack (ProcessData * pData, WavesData * wData, TracebackData * tbData) {
    char msg[MID_MESSAGE_SIZE];
    /*Perform Partition trace back*/
    tbData->pathParts ++;
    if (tbData->pathParts == 1) {
        tbData->aSeqLen = mmalloc (((MOATypeInd) sizeof *(tbData->aSeqLen)));
        if (tbData->aSeqLen == NULL) {
            printf ("Failed to reallocate memory for %ld partial alignments Lengths. Exiting.\n", tbData->pathParts);
            fflush (stdout);
            return;
        }
        tbData->algnseq = mmalloc (((MOATypeInd) sizeof *(tbData->algnseq)));
        if (tbData->aSeqLen == NULL) {
            printf ("Failed to reallocate memory for %ld partial alignments. Exiting.\n", tbData->pathParts);
            fflush (stdout);
            return;
        }
    }
    else {
        tbData->aSeqLen = realloc (tbData->aSeqLen, ((MOATypeInd) tbData->pathParts) * ((MOATypeInd) sizeof *(tbData->aSeqLen)));
        if (tbData->aSeqLen == NULL) {
            printf ("Failed to reallocate memory for %ld partial alignments Lengths. Exiting.\n", tbData->pathParts);
            fflush (stdout);
            return;
        }
        tbData->algnseq = realloc (tbData->algnseq, ((MOATypeInd) tbData->pathParts) * ((MOATypeInd) sizeof *(tbData->algnseq)));
        if (tbData->algnseq == NULL) {
            printf ("Failed to reallocate memory for %ld partial alignments. Exiting.\n", tbData->pathParts);
            fflush (stdout);
            return;
        }
    }
    tbData->algnseq[tbData->pathParts-1] = NULL;
    tbData->algnseq[tbData->pathParts-1] = mmalloc (((MOATypeInd) tbData->seqNum) * ((MOATypeInd) sizeof *(tbData->algnseq[tbData->pathParts-1])));    
    if (tbData->algnseq[tbData->pathParts-1] == NULL) {
        printf ("Failed to reallocate memory for %ld partial alignments sequences. Exiting.\n", tbData->pathParts);
        fflush (stdout);
        return;
    }
    tbData->aSeqLen[tbData->pathParts-1] = traceBack(pData, wData, tbData);
    sprintf(msg, "DMTB returned Local path of length %ld\n", tbData->aSeqLen[tbData->pathParts-1]);  
    mprintf(3, msg, 1);    
}
void * tbMaster (ProcessData * pData, WavesData * wData) {
    int MPI_return, foundproc, done = 0; /*currProc*/
    MOATypeDimn i, k;
    long long receviedLongLong;
    
    char msg[MID_MESSAGE_SIZE];
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
            printf ("Error Retrieving part No. Exiting.\n");
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
        myProcid = tbData->currProc;
        doPartitionTraceBack(pData, wData, tbData);

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
    
    assemblePathParts (tbData);
    calcAlignmentSPScore(pData, tbData);
    outputAlignment(pData, tbData, 1);
    editAlignment (pData, tbData);
    calcAlignmentSPScore (pData, tbData);
    outputAlignment(pData, tbData, 1);
    outputFastaAlignment(pData, tbData);    
    outputMSFAlignment(pData, tbData);    
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


int main (int argc, char * argv[]) {
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
        freeProcessMemory (&pData, &sData, &wData);
        return -1;
    }
    /*3. Load Tensor Partitions Computed*/
    ret = restoreCheckPoint (pData, wData);
    if (ret != 0) {
        printf (" Could not read Slave data from checkpoint file\n");
        freeProcessMemory (&pData, &sData, &wData);
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

    tbMaster(pData, wData);
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
            printf ("Failed to retrieve Process Resources Usage, errno %d\n", errno);

        //struct mallinfo info;
        //info = mallinfo();

        //printf("[%d] STime\tUTime\theap\tMemory\t\n",myProcid);
        //printf("[%d] %f\t%f\t%d\t%d\n", myProcid, stime, utime, info.arena, info.usmblks + info.uordblks);

        printf("STime\tUTime\n");
        printf("%f\t%f\n", stime, utime);
        currNow = getTime();
        printf("Finalized at time (%d, %d, %d, %d)\n", currNow->tm_yday, currNow->tm_hour, currNow->tm_min, currNow->tm_sec);
   freeProcessMemory (&pData, &sData, &wData);

    return 0;
}
