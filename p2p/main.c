/******************************************************************************
* FILE: main.c
* DESCRIPTION:  
* AUTHOR: Manal E. Helal
* LAST REVISED:
* Function:
*		main
*		MainProcess
*     ScoreCompThread
*     addOC
*     getDepProcs
*     sendOC
*     checkRecvOC
*     checkPrevPartitions
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <malloc.h>
#include <pthread.h>
#include <mpi.h>
#include "globals.h"
#include "utils.h"
#include "moaDst.h"
#include "main.h"

#include <errno.h>
//#include <linux/unistd.h>       /* for _syscallX macros/related stuff */
//#include <linux/kernel.h>       /* for struct sysinfo */

 	

//_syscall1(int, sysinfo, struct sysinfo *, info);




/*******************************************************************
	Function: ScoreCompThread
		Score Computation
	Input/Output:
		m: Process Data Structure
*******************************************************************/
void * ScoreCompThread (ProcessData * pData, ScoringData * sData, WavesData * wData) {
    MOATypeDimn k, j;
    long ocWave;
#ifndef NDEBUG
    int dbglevel = 1;
    char msg[MID_MESSAGE_SIZE];
#endif
	
#ifndef NDEBUG
    sprintf (msg, "[%d]>ScoreCompThread: Loop To Compute Scores for total [%ld] Partitions in this process\nEnter loop ******************************************\n", myProcid, pData->partitionsCount);
    mprintf (dbglevel, msg, 1);
#endif
    //while ((pData->waveNo < wData->wavesTotal) && (pData->partitionsCount > 0)) {
    sData->waveNo = pData->waveNo;
    pData->mpi_requests = NULL;
    pData->sendRequests = 0;
    while (pData->globalWaveNo < wData->wavesTotal) {
        while (pData->waveNo == pData->globalWaveNo) {
            getPartition (wData->partsInWaveIndices[pData->waveNo][pData->partNo], pData->seqNum, pData->seqLen, &pData->msaAlgn, wData->partitionSize);
            /* Compute Scored for Current Partition*/
#ifndef NDEBUG
            sprintf (msg, "[%d]>ScoreCompThread[%ld]: Will call ComputePartitionScores\n", myProcid, pData->computedPartitions);
            mprintf (dbglevel, msg, 1);
#endif
            if (Algorithm == DP) 
                DPComputeScores (pData, sData, wData);
            else if (Algorithm == SP)
                SPComputeScores (pData, sData, wData);
            //printMOA_scr(pData->msaAlgn, 0);
#ifndef NDEBUG
            sprintf (msg, "[%d]>ScoreCompThread[%ld]: Will call printMOA\n", myProcid, pData->computedPartitions);
            mprintf (dbglevel, msg, 1);		
            /* Print elements ======================================================= */
            printMOA(2, pData->msaAlgn, pData->sequences, 0);
            /* Print Indexes ========================================================*/
            printMOA(2, pData->msaAlgn, pData->sequences, 1);
            sprintf (msg, "[%d]>ScoreCompThread[%ld]: - wave: %ld order: %ld has %lld elm\n", myProcid, pData->computedPartitions-1, pData->waveNo, pData->partNo, pData->msaAlgn->elements_ub);
            mprintf (dbglevel, msg, 1);
            sprintf (msg, "[%d]>ScoreCompThread[%ld]: Will call getNextPartition\n", myProcid, pData->computedPartitions-1);
            mprintf (dbglevel, msg, 1);
#endif
            /****** For testing checkpoint resume ****************************
            if (!RestoreFlag && pData->computedPartitions == force_exit) {
                sprintf (msg, "[%d]>ScoreCompThread[%ld]: Forced to exit\n", myProcid, pData->computedPartitions);
                mprintf (dbglevel, msg, 1);
                break;
            }
            *****************************************************************/ 
            printf ("[%d]W %ld/%ld PO %ld/%ld Part %ld/%ld(T:%ld) PI {%lld ", myProcid, pData->waveNo, wData->wavesTotal, pData->partNo, wData->partsInWave[pData->waveNo], pData->computedPartitions+1, pData->partitionsCount, wData->partsTotal, pData->msaAlgn->indexes[0][0]);
            for (k=1;k<pData->seqNum;k++) 
                printf (", %lld", pData->msaAlgn->indexes[0][k]);
            printf ("}\n");
            fflush(stdout);
            getNextPartition (wData, &pData->waveNo, &pData->partNo);
            checkPoint (pData, sData);
            pData->computedPartitions ++;
        }
        if ((pData->waveNo > sData->waveNo) || (pData->waveNo != pData->globalWaveNo)) {
#ifndef NDEBUG
            sprintf (msg, "[%d]>ScoreCompThread[%ld] w%ld: Will communicate\n", myProcid, pData->computedPartitions, pData->waveNo);
            mprintf (dbglevel, msg, 1);
#endif
            if (pData->waveNo > sData->waveNo) {
                prepareWaveOC (sData->waveNo, pData);   
                sendOCtoHigherProcessors(pData);
            }
            MPI_Barrier(MOAMSA_COMM_WORLD);
            receiveOC(pData, sData); 
            if (pData->waveNo > sData->waveNo) 
                sendOCtoLowerProcessors(pData);            
            MPI_Barrier(MOAMSA_COMM_WORLD);
            receiveOC(pData, sData); 
            pData->globalWaveNo ++;
            /*Delete OC in waves before k (dimn or seqNum) waves from the current wave.*/
            if (pData->waveNo > pData->seqNum + 1) {
                ocWave = pData->waveNo - pData->seqNum - 1;
                if (pData->OCin != NULL) {
                    if ((pData->OCin[ocWave].WOCI != NULL) && (pData->OCin[ocWave].wavesOC > 0)) {
                        for (j=0;j<pData->OCin[ocWave].wavesOC;j++) {
                            if (pData->OCin[ocWave].WOCI[j].cellIndex != NULL) {
                                free (pData->OCin[ocWave].WOCI[j].cellIndex);
                                pData->OCin[ocWave].WOCI[j].cellIndex = NULL;
                            }
                        }
                        free (pData->OCin[ocWave].WOCI);                                  
                        pData->OCin[ocWave].WOCI = NULL;
                    }
                    pData->OCin[ocWave].wavesOC = 0;
                }
                if (pData->OCout != NULL) {
                    if ((pData->OCout[ocWave].WOCO != NULL) && (pData->OCout[ocWave].wavesOC > 0)) {
                        for (j=0;j<pData->OCout[ocWave].wavesOC;j++) {
                            if (pData->OCout[ocWave].WOCO[j].cellIndex != NULL) {
                                free (pData->OCout[ocWave].WOCO[j].cellIndex);
                                pData->OCout[ocWave].WOCO[j].cellIndex = NULL;
                            }
                            if ((pData->OCout[ocWave].WOCO[j].depProc_ub > 0) && (pData->OCout[ocWave].WOCO[j].depProc  != NULL)) {
                                free(pData->OCout[ocWave].WOCO[j].depProc);	
                                pData->OCout[ocWave].WOCO[j].depProc = NULL;
                            }
                        }
                        free (pData->OCout[ocWave].WOCO);                                  
                        pData->OCout[ocWave].WOCO = NULL;
                    }
                    pData->OCout[ocWave].wavesOC = 0;
                }
            }
        }

        //if ((pData->computedPartitions >= pData->partitionsCount) ||  (pData->partNo < 0)){
        if (pData->partNo < 0) {
            pData->waveNo = wData->wavesTotal; /*to get out of the main loop*/
        }
    }
    //PrintPrevChains (pData->msaAlgn);
    pData->compFinished = 1;		
    return NULL;
}

/* ============================================================================
	function MainProcess:
		seqNum: Number of sequences in string sequences (ex. 3)
		sequences: holds the sequences. (ex. [GTGCAACGTACT])
		seqLen: Array of the length of each sequence in sequences. (ex. 5,4,3)
		======= which means that we have three sequences GTGCA, ACGT, and ACT.
		stype: Scoring type (1: linear score, 2: PAM250 if protein, 3: BLOSUM if protein)
		partitionSize: Partion Size
============================================================================== */
void MainProcess (MOATypeDimn seqNum, char * * sequences, char * * seqName, MOATypeShape * seqLen, int stype, long partitionSize) {
    ProcessData * pData = NULL;
    ScoringData * sData = NULL;
    WavesData * wData = NULL;
    MOATypeDimn k;
    MOATypeInd i;
    int ret, startflag;
    struct rusage usageRec;
    double utime, stime;
    MPI_Status status;
    double t_start, t_finish;
    char command[2000];

#ifndef NDEBUG
    char msg[SHORT_MESSAGE_SIZE];
    int dbglevel = 0;
    MOATypeInd j;
#endif

    t_start = MPI_Wtime();
    /* print the input arguments ============================================*/
    //PrintSequencies (0, seqNum, sequences, seqLen);
#ifndef NDEBUG
    sprintf(msg, ">>>>MainProcess: Scoring Type: %d\n>>>>Partition Size: %ld\n", stype, partitionSize);	
    mprintf(dbglevel, msg, 1);
#endif

    /* Initialize Process Memory pData (function located in partitioning.c*/

    ret = initProcessMemory(&pData, &sData, &wData, seqNum, seqLen, sequences, seqName, stype, partitionSize);  
    if (ret != 0) {
        mprintf (0, ">>>>MainProcess: Error Initializing Process Data, Exiting\n", 1);
        fflush (stdout);
        return;
    }
    /* if restore previouse run read check point data here, do not calculate waves */
    pData->OCout = NULL;
    pData->OCin = NULL;
    if (Mode != Distributed) {

        /* Construct MOA record */
        pData->msaAlgn =  NULL;
        createMOAStruct(&pData->msaAlgn);
        if (createMOA(seqLen, seqNum, pData->msaAlgn, 0, 0) < 0)
            return;		
        wData->wavesTotal = 1;
        wData->AllpartsInWave =  mmalloc((MOATypeInd) sizeof *wData->AllpartsInWave);
        wData->AllpartsInWave[0] = 1;
        pData->waveNo =  0;
        pData->partNo = 0;
        pData->partitionsCount = 1;
        /*pData->OCout = mmalloc((MOATypeInd) sizeof *(pData->OCout));
        if (pData->OCout == NULL) {
            mprintf(1, "Couldn't create memory for OCout while adding an OC. Exiting.\n", 3);
            printf("Couldn't create memory for OCout while . Exiting.\n",);
            return;
        }
        pData->OCout[0].wavesOC = 0;
        pData->OCout[0].WOCO = NULL;
        */
#ifndef NDEBUG
        sprintf (msg, "[%d]>ScoreCompThread[%ld]: Will call ComputePartitionScores\n", myProcid, pData->computedPartitions);
        mprintf (dbglevel, msg, 1);
#endif
            /* Compute Scored for Current Partition*/
        if (Algorithm == DP) 
            DPComputeScores (pData, sData, wData);
        else if (Algorithm == SP)
            DPComputeScores (pData, sData, wData);
        /* Print elements ======================================================= */
        //printMOA_scr(pData->msaAlgn, 0);
        /* Print Indexes ========================================================*/
        //printMOA_scr(pData->msaAlgn, 1);

        pData->computedPartitions ++;
        checkPoint (pData, sData);		
    }
    else {
         if (RestoreFlag == 1) {
            /* restore data */
            restoreCheckPoint (pData, wData);
            pData->globalWaveNo =  pData->waveNo;
        } else {
            pData->partNo = 0;
            pData->waveNo = 0;
            if (myProcid == 0) {
                printf ("[%d] Calculating waves and partitions .... ", myProcid);
                calcWaves (pData, wData);
                currNow = getTime();
                printf("[%d] Done Calculating waves and partitions time (%d, %d, %d, %d)\n", myProcid, currNow->tm_yday, currNow->tm_hour, currNow->tm_min, currNow->tm_sec);
                fflush(stdout);
                if( checkPointWavesCalculations (pData, wData) == 0) {
                    startflag = 1;
                    for (i=1; i<ClusterSize; i++)
                        MPI_Send(&startflag, 1, MPI_INT, i, 0, MOAMSA_COMM_WORLD);
                }
                else {
                    printf ("[%d]Couldn't write Waves calculations, Exiting\n", myProcid);
                    return;
                }
            }
            else {
                //printf ("[%d] waiting for start flag\n", myProcid);
                MPI_Recv(&startflag, 1, MPI_INT, 0, 0, MOAMSA_COMM_WORLD, &status);
                //printf ("[%d] received start flag = %d\n", myProcid, startflag);
                printf ("[%d] Reading waves and partitions .... ", myProcid);
                if (restoreWavesCalculations(pData, wData) != 0) {
                    printf ("[%d]Couldn't read Waves calculations, Exiting\n", myProcid);
                    return;
                }
                printf ("done.\n");
                fflush(stdout);
            }
        }	

#ifndef NDEBUG
        sprintf(msg, "[%d]>MainProcess: Current Wave: %ld - Current Partition: %ld - Total Partitions in Process: %ld\n", myProcid, pData->waveNo, pData->partNo, pData->partitionsCount);	
        mprintf(dbglevel, msg, 1);
#endif
        ScoreCompThread (pData, sData, wData);
    }
    if (myProcid == 0) {
        t_finish = MPI_Wtime();
        /* Getting Process Resources Usage ===================== */
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

        printf("STime\tUTime\n");
        printf("%f\t%f\n", stime, utime);
        printf ("Elsp-time: %f\n", t_finish - t_start);
        fflush(stdout);
        sprintf (command, "prstat 1 1 > /export/home/mhelal1/thesis/exp/run/prstatus/prst_%s", outputfilename);
        i = system (command);
    }
    /* Free allocated memory and exit routine ===================== */

    freeProcessMemory (&pData, &sData, &wData);
}

/* =============================================================================
	The main function:
		Initialize MPI, 
		Process arguments (processArguments)
		Initialize output debugging files (init_output)
		call MainProcess function, 
		and finalizes MPI.

============================================================================= */
int main(int argc, char **argv) {
    MPI_Group orig_group;	
    char * * sequences = NULL;/*, ufilename[SHORT_MESSAGE_SIZE];*/
    char * * seqName = NULL;
    long partitionSize;
    MOATypeDimn seqNum;
    MOATypeShape * seqLen = NULL;
    int stype;
    long ID1;
#ifndef NDEBUG
    char msg[SHORT_MESSAGE_SIZE];
#endif			
    /* MPI Initiliaztion ==================================================*/
    MPI_Init(&argc, &argv);
    /* Get my Rank in myProcid*/
    MPI_Comm_rank(MPI_COMM_WORLD, &myProcid);
    /* Get the number of processes running in the ClusterSize*/
    MPI_Comm_size(MPI_COMM_WORLD, &ClusterSize);
    /* Get the group associated with the communicator ===============*/
    MPI_Comm_group (MPI_COMM_WORLD, &orig_group);
    MPI_Comm_create(MPI_COMM_WORLD, orig_group, &MOAMSA_COMM_WORLD);
    /* ============================================ end MPI initialization*/
    TBFlag = 0;
    /* 1. Process Arguments, Sequences, lengths, scoring type, partition Size, output prefix, ... etc*/
    processArguments(argc, argv, &seqNum, &sequences, &seqName, &seqLen, &stype, &partitionSize);
    /*  Thought I can run the same program in Distributed and Sequential Mode, but only now the Distributed one is being tested, so ignore this one*/

    /* Initialize output debugging files*/
    strcpy (outPrefix, "c");
    if (init_output() == 0) {
        /* Initialize timing variables, previous Time (prevNow) and Current Time  (currNow)*/
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
#ifndef NDEBUG
        sprintf(msg, "[%d]>Started at time (%d, %d, %d, %d)\n", myProcid, currNow->tm_yday, currNow->tm_hour, currNow->tm_min, currNow->tm_sec);	
        mprintf(0, msg, 1);
#endif
        cpTime (currNow, &prevNow);

        /* Calling the main process ==========================================*/
        MainProcess(seqNum, sequences, seqName, seqLen, stype, partitionSize);
        currNow = getTime();
if (myProcid == 0) {
        printf("[%d]>Finalized at time (%d, %d, %d, %d)\n", myProcid, currNow->tm_yday, currNow->tm_hour, currNow->tm_min, currNow->tm_sec);
}
#ifndef NDEBUG
        sprintf(msg, "[%d]>Finalized at time (%d, %d, %d, %d)\n", myProcid, currNow->tm_yday, currNow->tm_hour, currNow->tm_min, currNow->tm_sec);	
        mprintf(0, msg, 1);
#endif
        if (prevNow != NULL)
            free (prevNow);
    }
    if (close_output () != 0) 
        printf ("[%d] Error closing output files\n", myProcid);
    /* Finalize MPI ==================================================*/
    MPI_Finalize();

    return 0;
} /* of main */
