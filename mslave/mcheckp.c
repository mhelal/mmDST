#include <stdio.h>
#include <stdlib.h>
#include <limits.h> 
#include <string.h>
#include <pthread.h>
#include <semaphore.h>
#include "master.h"
#include "moaDst.h"
#include "utils.h"

#define _GNU_SOURCE

int initMasterMemory (Master * * master, long seqNum, long * seqLen, int stype, char * * sequences) {
	char msg[MID_MESSAGE_SIZE];

  (*master) = ( Master*)mmalloc (sizeof (Master));
  if ((*master) == NULL) 
		return -1;
  (*master)->seqNum = seqNum;
  (*master)->seqLen = seqLen;
  (*master)->stype = stype;
  (*master)->sequences = sequences;
  createMOAStruct (&(*master)->msaAlgn);
  createMOA(seqLen /*shape*/, seqNum /*dimension*/, (*master)->msaAlgn /*MOA structure*/,-1,0);
  sprintf (msg, " MSAAlign dimn %ld shape[%ld] %ld elements_up %ld\n", (*master)->msaAlgn->dimn, (*master)->msaAlgn->dimn-1, (*master)->msaAlgn->shape[(*master)->msaAlgn->dimn-1], (*master)->msaAlgn->elements_ub);
  mprintf(1, msg, threadnum);
  (*master)->ComputationPhase = NULL;
  (*master)->ComputationPhase = (int *) mmalloc (sizeof (int) * ClusterSize);
  if ((*master)->ComputationPhase == NULL) 
		return -1;
  (*master)->mut = NULL;
  (*master)->mut = (pthread_mutex_t *) mmalloc (sizeof (pthread_mutex_t));
  if ((*master)->mut == NULL) 
		return -1;
  pthread_mutex_init ((*master)->mut, NULL);
  sem_init(&(*master)->qSem, 0, 0);
	return 0;
}

void freeMasterMemory (Master * master) {
	long i;
	char msg[SHORT_MESSAGE_SIZE];

  if (DS != NULL) {
		sprintf (msg, "Master DS_ubound = %ld\n", DS_ubound);
		mprintf(12, msg, threadnum);
    for (i=0;i<DS_ubound;i++) {
			sprintf (msg, "Master i = %ld lb %p hb %p\n", i, DS[i].lbCells, DS[i].hbCells);
			mprintf(12, msg, threadnum);
      if (DS[i].lbCells != NULL) 
				free (DS[i].lbCells);
			DS[i].lbCells = NULL;
      if (DS[i].hbCells != NULL) 
				free (DS[i].hbCells);
			DS[i].hbCells = NULL;
			sprintf (msg, "Master freed DS[%ld]\n", i);
			mprintf(12, msg, threadnum);
    }
    free (DS);
		DS = NULL;
		sprintf (msg, "Master freed all DS master = %p\n", master);
		mprintf(12, msg, threadnum);
  }
  if (master != NULL) {
		if (master->msaAlgn != NULL)
		  deleteMOA (master->msaAlgn);
		master->msaAlgn = NULL;
		mprintf (12, "Master freed msaAlgn\n", threadnum);
	  if (master->mut != NULL) {
			pthread_mutex_destroy(master->mut);
			free (master->mut);
		}
	  master->mut = NULL;
		mprintf (12, "Master freed mut\n", threadnum);
	  /*if (master->ComputationPhase != NULL)
	    free (master->ComputationPhase);
		master->ComputationPhase = NULL;
			mprintf (12, "Master freed ComputationPhase\n", threadnum);
			*/
    free(master);
		master = NULL;
		mprintf (12, "Master freed all its data\n", threadnum);
  }
}
int checkPointMaster (void * threadarg) {
	FILE * sfile;
  Master * master;
	long i, j, k;

  master = (Master *) threadarg;
	char sfilename[20], bkfilename[20];
	sprintf(sfilename, "master_%d_chp", myProcid);
	sprintf(bkfilename, "master_%d_chp_bk", myProcid);
	file_copy( sfilename, bkfilename, get_bytes(0.0, sfilename), get_bytes(100.0, bkfilename) );
		
	if (( sfile= fopen (sfilename, "w")) == NULL) {
		mprintf(1, "Can not write checkpoint file\n", threadnum);
		return -1;
	}
	fprintf (sfile, "pid=%d\n", myProcid);
	fprintf (sfile, "cls=%d\n", ClusterSize);
	for (i=0;i<ClusterSize;i++)
	  fprintf (sfile, "%d\n", master->ComputationPhase[i]);
	fprintf (sfile, "stp=%d\n", master->stype);
	fprintf (sfile, "sqm=%ld\n", master->seqNum);
	for (i = 0;i < master->seqNum;i++) {
		fprintf (sfile, "%ld\n", master->seqLen[i]);
		fprintf (sfile, "%s\n", master->sequences[i]);
	}
	fprintf (sfile, "wve=%ld\n", master->waveNo);
	fprintf (sfile, "prt=%ld\n", DS_ubound);
	for (i=0;i<DS_ubound;i++) {
		fprintf (sfile, "%ld\n", DS[i].partIndex);
		fprintf (sfile, "%ld\n", DS[i].waveNo);
		fprintf (sfile, "%d\n", DS[i].sent);
		fprintf (sfile, "%ld\n", DS[i].lbCellsCount);
		fprintf (sfile, "%ld\n", DS[i].hbCellsCount);
		for (j=0;j<DS[i].lbCellsCount;j++) {
		  fprintf (sfile, "%ld\n", DS[i].lbCells[j]);
		}
		for (j=0;j<DS[i].hbCellsCount;j++) {
		  fprintf (sfile, "%ld\n", DS[i].hbCells[j]);
		}
	}
	fprintf (sfile, "dpc=%ld\n", depProclistCount);
	for (i=0;i<depProclistCount;i++) {
		fprintf (sfile, "%ld\n", depProclist[i].cellIndex);
		fprintf (sfile, "%ld\n", depProclist[i].partIndex);
		fprintf (sfile, "%ld\n", depProclist[i].waveNo);
		fprintf (sfile, "%d\n", depProclist[i].depProc);
		fprintf (sfile, "%d\n", depProclist[i].sent);
	}
	fclose (sfile);
	return 0;
}

int restoreMasterCheckPoint (void * threadarg) {
	FILE * sfile;
	char sfilename[20];
	char msg[SHORT_MESSAGE_SIZE];
	char * t1 = NULL;
	int sProcid;
	long i, j, k;
  Master * master;
	char * read;
char line[LINE_MAX];

  master = (Master *) threadarg;

	sprintf(sfilename, "master_%d_chp", myProcid);

	if (( sfile= fopen (sfilename, "r")) == NULL) {
		mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
		return -1;
  }

	while ((read = fgets (line, LINE_MAX, sfile))  != NULL) {
		if ((line[0] == 'd') && (line[1] == 'p') && (line[2] == 'c')) {
			t1 = strtok(line,"=");
			t1 = strtok(NULL,"=");
			depProclistCount = atol (t1);
			sprintf(msg, "read depProclistCount %ld\n", depProclistCount);
			mprintf(6, msg, threadnum);
			depProclist = (depProc_rec * ) mmalloc (depProclistCount *sizeof(depProc_rec));
			for (i=0;i<depProclistCount;i++) {
				read = fgets (line, LINE_MAX, sfile);
				if (read == NULL) {
					mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
					return -1;
				}					
				depProclist[i].cellIndex = atol (line);
				sprintf(msg, "read depProclist[%ld].cellIndex %ld\n", i, depProclist[i].cellIndex);
				mprintf(6, msg, threadnum);
				read = fgets (line, LINE_MAX, sfile);
				if (read == NULL) {
					mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
					return -1;
				}					
				depProclist[i].partIndex = atol (line);
				sprintf(msg, "read depProclist[%ld].partIndex %ld\n", i,depProclist[i].partIndex);
				mprintf(6, msg, threadnum);
				read = fgets (line, LINE_MAX, sfile);
				if (read == NULL) {
					mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
					return -1;
				}					
				depProclist[i].waveNo = atol (line);
				sprintf(msg, "read depProclist[%ld].waveNo %ld\n", i,depProclist[i].waveNo);
				mprintf(6, msg, threadnum);
				read = fgets (line, LINE_MAX, sfile);
				if (read == NULL) {
					mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
					return -1;
				}					
				depProclist[i].depProc = atol (line);
				sprintf(msg, "read depProclist[%ld].depProc %d\n", i,depProclist[i].depProc);
				mprintf(6, msg, threadnum);
				read = fgets (line, LINE_MAX, sfile);
				if (read == NULL) {
					mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
					return -1;
				}					
				depProclist[i].sent = atol (line);
				sprintf(msg, "read depProclist[%ld].sent %d\n", i,depProclist[i].sent);
				mprintf(6, msg, threadnum);
			}
		}
		else if ((line[0] == 'p') && (line[1] == 'r') && (line[2] == 't')) {
			t1 = strtok(line,"=");
			t1 = strtok(NULL,"=");
			DS_ubound = atol (t1);
			sprintf(msg, "read DS_ubound %ld\n", DS_ubound);
			mprintf(6, msg, threadnum);
			DS = (Distribution * ) mmalloc (DS_ubound * sizeof (Distribution));
			for (i=0;i<DS_ubound;i++) {
				read = fgets (line, LINE_MAX, sfile);
				if (read == NULL) {
					mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
					return -1;
				}					
				DS[i].partIndex = atol (line);
				sprintf(msg, "read DS[%ld].partIndex %ld\n", i,DS[i].partIndex );
				mprintf(6, msg, threadnum);
				read = fgets (line, LINE_MAX, sfile);
				if (read == NULL) {
					mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
					return -1;
				}					
				DS[i].waveNo = atol (line);
				sprintf(msg, "read DS[%ld].waveNo %ld\n", i,DS[i].waveNo);
				mprintf(6, msg, threadnum);
				read = fgets (line, LINE_MAX, sfile);
				if (read == NULL) {
					mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
					return -1;
				}					
				DS[i].sent = atol (line);
				sprintf(msg, "read DS[%ld].sent %d\n", i,DS[i].sent);
				mprintf(6, msg, threadnum);
				read = fgets (line, LINE_MAX, sfile);
				if (read == NULL) {
					mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
					return -1;
				}					
				DS[i].lbCellsCount = atol (line);
				sprintf(msg, "read DS[%ld].lbCellsCount %ld\n", i,DS[i].lbCellsCount);
				mprintf(6, msg, threadnum);
				read = fgets (line, LINE_MAX, sfile);
				if (read == NULL) {
					mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
					return -1;
				}					
				DS[i].hbCellsCount = atol (line);
				sprintf(msg, "read DS[%ld].hbCellsCount %ld\n", i,DS[i].hbCellsCount);
				mprintf(6, msg, threadnum);
				DS[i].lbCells = (long *) mmalloc (DS[i].lbCellsCount * sizeof(long));
				DS[i].hbCells = (long *) mmalloc (DS[i].hbCellsCount * sizeof(long));
				for (j=0;j<DS[i].lbCellsCount;j++) {
				  read = fgets (line, LINE_MAX, sfile);
				  if (read == NULL) {
				    mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
				    return -1;
				  }					
				  DS[i].lbCells[j] = atol (line);
				  sprintf(msg, "read DS[%ld].lbCells[%ld] %ld\n", i,j, DS[i].lbCells[j]);
				  mprintf(6, msg, threadnum);
				}
				for (j=0;j<DS[i].hbCellsCount;j++) {
				  read = fgets (line, LINE_MAX, sfile);
				  if (read == NULL) {
				    mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
				    return -1;
				  }					
				  DS[i].hbCells[j] = atol (line);
				  sprintf(msg, "read DS[%ld].hbCells[%ld] %ld\n", i,j, DS[i].hbCells[j]);
				  mprintf(6, msg, threadnum);
				}
			}
		}
		else if ((line[0] == 'w') && (line[1] == 'v') && (line[2] == 'e')) {
			t1 = strtok(line,"=");
			t1 = strtok(NULL,"=");
		        master->waveNo = atol (t1);
			sprintf(msg, "read master->waveNo %ld\n", master->waveNo);
			mprintf(6, msg, threadnum);
		}
		else if ((line[0] == 's') && (line[1] == 'q') && (line[2] == 'm')) {
			t1 = strtok(line,"=");
			t1 = strtok(NULL,"=");
			master->seqNum = atol (t1);
			sprintf(msg, "read master->seqNum %ld\n", master->seqNum);
			mprintf(6, msg, threadnum);
			master->seqLen = (long * ) mcalloc (master->seqNum, sizeof(long));
			master->sequences = (char * * ) mmalloc (master->seqNum * sizeof(char *));
			for (i=0;i<master->seqNum;i++) {
				read = fgets (line, LINE_MAX, sfile);
				if (read == NULL) {
					mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
					return -1;
				}					
				master->seqLen[i] = atol (line);
				sprintf(msg, "read master->seqLen[%ld] %ld\n", i, master->seqLen[i]);
				mprintf(6, msg, threadnum);
				read = fgets (line, LINE_MAX, sfile);
				if (read == NULL) {
					mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
					return -1;
				}			
				master->sequences[i] = (char * ) mmalloc (strlen(line) * sizeof(char));
		
				master->sequences[i] = line;
				sprintf(msg, "read master->sequences[%ld] %ld\n", i, master->sequences[i]);
				mprintf(6, msg, threadnum);
			}
		}
		else if ((line[0] == 's') && (line[1] == 't') && (line[2] == 'p')) {
			t1 = strtok(line,"=");
			t1 = strtok(NULL,"=");
			master->stype = atol (t1);
			sprintf(msg, "read stype %d\n", master->stype);
			mprintf(6, msg, threadnum);
		}
		else if ((line[0] == 'c') && (line[1] == 'l') && (line[2] == 's')) {
			t1 = strtok(line,"=");
			t1 = strtok(NULL,"=");
			ClusterSize = atol (t1);
			sprintf(msg, "read ClusterSize %d\n", ClusterSize);
			mprintf(6, msg, threadnum);
			master->ComputationPhase = (int *) mmalloc (sizeof(int) * ClusterSize);
			for (i=0;i<ClusterSize;i++){
				read = fgets (line, LINE_MAX, sfile);
				if (read == NULL) {
					mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
					return -1;
				}					
				master->ComputationPhase[i] = atol (line);
				sprintf(msg, "read master->seqLen[%ld] %ld\n", i, master->seqLen[i]);
				mprintf(6, msg, threadnum);
			  			  
			}
		}
		else if ((line[0] == 'p') && (line[1] == 'i') && (line[2] == 'd')) {
			t1 = strtok(line,"=");
			t1 = strtok(NULL,"=");
			sProcid = atol (t1);
			for (i=0;i<ClusterSize;i++) {
				read = fgets (line, LINE_MAX, sfile);
				if (read == NULL) {
					mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
					return -1;
				}					
				master->ComputationPhase[i] = atol (line);
			}
		}

	}
  fclose (sfile);

		//if (t1 != NULL)
			//free(t1);
	return 0;
}
