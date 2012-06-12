#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <semaphore.h>
#include "slave.h"
#include "globals.h"
#include "moaDst.h"
#include "utils.h"
#include <pthread.h>

int initSlaveMemory (Slave * * slave, long seqNum, long * seqLen, int stype) {
  (*slave) = (Slave *) mmalloc (sizeof(Slave));
  (*slave)->stype = stype;
  (*slave)->seqNum = seqNum;
  (*slave)->seqLen = seqLen;

	(*slave)->compThreadFinished = 0;
  OCin = NULL;
  OCin_ub = 0;
  OCout = NULL;
  OCout_ub = 0;
  (*slave)->DPmut = (pthread_mutex_t *) mmalloc (sizeof (pthread_mutex_t));
  pthread_mutex_init ((*slave)->DPmut, NULL);
  sem_init(&(*slave)->dpSem, 0, 0);

  (*slave)->OCmut = (pthread_mutex_t *) mmalloc (sizeof (pthread_mutex_t));
  pthread_mutex_init ((*slave)->OCmut, NULL);
  sem_init(&(*slave)->ocSem, 0, 0);

  (*slave)->CPmut = (pthread_mutex_t *) mmalloc (sizeof (pthread_mutex_t));
  pthread_mutex_init ((*slave)->CPmut, NULL);
  sem_init(&(*slave)->cpSem, 0, 0);
  sem_init(&(*slave)->retSem, 0, 0);

  (*slave)->partitionsCount = 0;
  (*slave)->processedPartitions = 0;
  return 0;
}

void freeSlaveMemory (Slave * slave) {
  long i, j;
	if ((depProclistCount > 0) && (depProclist != NULL)) 
		free(depProclist);
	depProclistCount = 0;
	depProclist = NULL;

	if ((OCin_ub > 0) && (OCin != NULL))
		free(OCin);
	OCin_ub = 0;
	OCin = NULL;
	

  if ((OCout_ub > 0) && (OCout  != NULL))
      free(OCout);
	OCout_ub = 0;
	OCout = NULL;    
  
	if (slave != NULL) {
		if ((slave->seqNum > 0) && (slave->seqLen != NULL))
	    free (slave->seqLen);
		slave->seqLen = NULL;
		slave->seqNum = 0;
		if ((slave->partitionsCount > 0) && (slave->MOAPart!= NULL)) {
			for (i=0;i<slave->partitionsCount;i++) {
				
				for (j=0;j<slave->MOAPart[i].msaAlgn->dimn;j++) {
					if (slave->MOAPart[i].sequences[j] != NULL)
						free (slave->MOAPart[i].sequences[j]);
					slave->MOAPart[i].sequences[j] = NULL;
				}
				free(slave->MOAPart[i].sequences);
				slave->MOAPart[i].sequences = NULL;
			  deleteMOA (slave->MOAPart[i].msaAlgn);
			}
		}
		free (slave->MOAPart);
		slave->partitionsCount = 0;
		slave->processedPartitions = 0;
		slave->ComputationPhase = 0;
		slave->MOAPart = NULL;
		if (slave->DPmut != NULL)
  	  free (slave->DPmut);
		slave->DPmut = NULL;
  	sem_destroy(&slave->dpSem);
		if (slave->OCmut != NULL)
	    free (slave->OCmut);
		slave->OCmut = NULL;
	  sem_destroy(&slave->ocSem);
		if (slave->CPmut != NULL)
	    free (slave->CPmut);
		slave->CPmut = NULL;
	  sem_destroy(&slave->cpSem);
  	sem_destroy(&slave->retSem);
    free(slave);
  }
		slave = NULL;
}

int checkPointSlave (void * threadarg) {
	FILE * sfile;
  Slave * slave;
	long i, j, k;

  slave = (Slave *) threadarg;
	char sfilename[20], bkfilename[20];
	sprintf(sfilename, "slave_%d_chp", myProcid);
	sprintf(bkfilename, "slave_%d_chp_bk", myProcid);
	file_copy( sfilename, bkfilename, get_bytes(0.0, sfilename), get_bytes(100.0, bkfilename) );
		
	if (( sfile= fopen (sfilename, "w")) == NULL) {
		mprintf(1, "Can not write checkpoint file.\n", threadnum);
		return -1;
	}
	fprintf (sfile, "pid=%d\n", myProcid);
	fprintf (sfile, "cmp=%d\n", slave->ComputationPhase);
	fprintf (sfile, "stp=%d\n", slave->stype);
	fprintf (sfile, "sqm=%ld\n", slave->seqNum);
	for (i = 0;i < slave->seqNum;i++) {
		fprintf (sfile, "%ld\n", slave->seqLen[i]);
	}
	fprintf (sfile, "prt=%ld\n", slave->partitionsCount);
	for (i = 0;i < slave->partitionsCount;i++) {
		fprintf (sfile, "ppc=%d\n", slave->MOAPart[i].processed);
		fprintf (sfile, "pdm=%ld\n", slave->MOAPart[i].msaAlgn->dimn);
		for (j = 0;j < slave->MOAPart[i].msaAlgn->dimn;j++) {
			fprintf (sfile, "%ld\n", slave->MOAPart[i].msaAlgn->shape[j]);
			for (k = 0;k < slave->MOAPart[i].msaAlgn->shape[j];k++) {
				fprintf (sfile, "%c\n", slave->MOAPart[i].sequences[j][k]);
			}
		}
		fprintf (sfile, "pel=%ld\n", slave->MOAPart[i].msaAlgn->elements_ub);
		for (j = 0;j < slave->MOAPart[i].msaAlgn->elements_ub;j++) {
			fprintf (sfile, "%ld\n", slave->MOAPart[i].msaAlgn->elements[j].val);			
			fprintf (sfile, "%ld\n", slave->MOAPart[i].msaAlgn->elements[j].prev_ub);			
			for (k = 0; k < slave->MOAPart[i].msaAlgn->elements[j].prev_ub; k++)  {
				fprintf (sfile, "%ld\n", slave->MOAPart[i].msaAlgn->elements[j].prev[k]);			
			}
		}
		for (j = 0;j < slave->MOAPart[i].msaAlgn->elements_ub;j++) {
			fprintf (sfile, "%ld\n", slave->MOAPart[i].msaAlgn->indexes[j]);			
		}
	}
	/*Saving Overlapping Incoming Cells Scores Received*/
	fprintf (sfile, "icn=%ld\n", OCin_ub);
	for (i=0;i<OCin_ub;i++) {
		fprintf (sfile, "%ld\n", OCin[i].waveNo);
		fprintf (sfile, "%ld\n", OCin[i].cellIndex);
		fprintf (sfile, "%ld\n", OCin[i].cellScore);	
		fprintf (sfile, "%d\n", OCin[i].fromProc);	
	}
	/* if still computing scores then save extra info to resume*/
	if (slave->ComputationPhase < 2) {
		/*Saving Dependancy Information*/
		fprintf (sfile, "dpn=%ld\n", depProclistCount);
		for (i=0;i<depProclistCount;i++) {
			fprintf (sfile, "%ld\n", depProclist[i].cellIndex);
			fprintf (sfile, "%d\n", depProclist[i].depProc);
		}
		/*Saving Overlapping Outgoing Cells Scores*/
		fprintf (sfile, "ocn=%ld\n", OCout_ub);
		for (i=0;i<OCout_ub;i++) {
			fprintf (sfile, "%ld\n", OCout[i].waveNo);
			fprintf (sfile, "%ld\n", OCout[i].cellIndex);
			fprintf (sfile, "%ld\n", OCout[i].cellScore);	
			fprintf (sfile, "%d\n", OCout[i].toProc);	
		}
	}
	fclose(sfile);
	return 0;
}

int restoreSlaveCheckPoint (void * threadarg) {
	FILE * sfile;
	char sfilename[20], msg[SHORT_MESSAGE_SIZE];
	char line[LINE_MAX];
	char * t1 = NULL;
	int sProcid;
	long i, j, k;
  Slave * slave;
	char * read;

  slave = (Slave *) threadarg;

	sprintf(sfilename, "slave_%d_chp", myProcid);

	if (( sfile= fopen (sfilename, "r")) == NULL) {
		mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
		return -1;
  }

	sprintf(sfilename, "slave_%d_chp", myProcid);
	while ((read = fgets (line, LINE_MAX, sfile)) != NULL) {
		//if (slave->ComputationPhase < 2) {
			if ((line[0] == 'o') && (line[1] == 'c') && (line[2] == 'n')) {
				t1 = strtok(line,"=");
				t1 = strtok(NULL,"=");
				OCout_ub = atol (t1);
				sprintf(msg, "read OCout_ub %ld\n", OCout_ub);
				mprintf(20, msg, threadnum);
				OCout = (OCOType *) mmalloc (sizeof(OCOType) * OCout_ub);	
				for (i=0;i<OCout_ub;i++) {
					read = fgets (line, LINE_MAX, sfile);
					if (read == NULL) {
						mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
						return -1;
					}
					OCout[i].waveNo = atol(line);
					sprintf(msg, "read OCout[%ld].waveNo %ld\n", i, OCout[i].waveNo);
					mprintf(20, msg, threadnum);
					read = fgets (line, LINE_MAX, sfile);
					if (read == NULL) {
						mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
						return -1;
					}
					OCout[i].cellIndex = atol(line);
					sprintf(msg, "read OCout[%ld].cellIndex %ld\n", i, OCout[i].cellIndex);
					mprintf(20, msg, threadnum);
					read = fgets (line, LINE_MAX, sfile);
					if (read == NULL) {
						mprintf(6, "Can not read checkpoint file, exiting.\n", threadnum);
						return -1;
					}
					OCout[i].cellScore = atol(line) ;
					sprintf(msg, "read OCout[%ld].cellScore %ld\n", i, OCout[i].cellScore);
					mprintf(20, msg, threadnum);
					read = fgets (line, LINE_MAX, sfile);
					if (read == NULL) {
						mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
						return -1;
					}
					OCout[i].toProc = atol(line) ;
					sprintf(msg, "read OCout[%ld].toProc %d\n", i, OCout[i].toProc);
					mprintf(20, msg, threadnum);
				}
			}		
			else if ((line[0] == 'd') && (line[1] == 'p') && (line[2] == 'n')) {
				t1 = strtok(line,"=");
				t1 = strtok(NULL,"=");
				depProclistCount = atol (t1);
				sprintf(msg, "read depProclistCount %ld\n", depProclistCount);
				mprintf(20, msg, threadnum);
				depProclist = (depProc_rec *) mmalloc (sizeof(depProc_rec) * depProclistCount);	
				for (i=0;i<depProclistCount;i++) {
					read = fgets (line, LINE_MAX, sfile);
					if (read == NULL) {
						mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
						return -1;
					}
					depProclist[i].cellIndex = atol(line);
					sprintf(msg, "read depProclist[%ld].cellIndex %ld\n", i, depProclist[i].cellIndex);
					mprintf(20, msg, threadnum);
					read = fgets (line, LINE_MAX, sfile);
					if (read == NULL) {
						mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
						return -1;
					}
					depProclist[i].depProc = atol(line) ;
					sprintf(msg, "read depProclist[%ld].depProc %d\n", i, depProclist[i].depProc);
					mprintf(20, msg, threadnum);
				}					
			}
			else if ((line[0] == 'i') && (line[1] == 'c') && (line[2] == 'n')) {
				t1 = strtok(line,"=");
				t1 = strtok(NULL,"=");
				OCin_ub = atol (t1);
				sprintf(msg, "read OCin_ub %ld\n", OCin_ub);
				mprintf(20, msg, threadnum);
				OCin = (OCIType *) mmalloc (sizeof(OCIType) * OCin_ub);	
				for (i=0;i<OCin_ub;i++) {
					read = fgets (line, LINE_MAX, sfile);
					if (read == NULL) {
						mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
						return -1;
					}
					OCin[i].waveNo = atol(line);
					sprintf(msg, "read OCin[%ld].waveNo %ld\n", i, OCin[i].waveNo);
					mprintf(20, msg, threadnum);
					read = fgets (line, LINE_MAX, sfile);
					if (read == NULL) {
						mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
						return -1;
					}
					OCin[i].cellIndex = atol(line);
					sprintf(msg, "read OCin[%ld].cellIndex %ld\n", i, OCin[i].cellIndex);
					mprintf(20, msg, threadnum);
					read = fgets (line, LINE_MAX, sfile);
					if (read == NULL) {
						mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
						return -1;
					}
					OCin[i].cellScore = atol(line) ;
					sprintf(msg, "read OCin[%ld].cellScore %ld\n", i, OCin[i].cellScore);
					mprintf(20, msg, threadnum);
					read = fgets (line, LINE_MAX, sfile);
					if (read == NULL) {
						mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
						return -1;
					}
					OCin[i].fromProc = atol(line) ;
					sprintf(msg, "read OCin[%ld].fromProc %d\n", i, OCin[i].fromProc);
					mprintf(20, msg, threadnum);
				}
			}		
		//}
		else if ((line[0] == 'p') && (line[1] == 'r') && (line[2] == 't')) {
			t1 = strtok(line,"=");
			t1 = strtok(NULL,"=");
			slave->partitionsCount = atol (t1);
			sprintf(msg, "read slave->partitionsCount %ld\n", slave->partitionsCount);
			mprintf(20, msg, threadnum);
	  	slave->MOAPart = (MOAPartition *) mmalloc (sizeof(MOAPartition) * slave->partitionsCount);    
			mprintf(20, "created Partitions count\n", threadnum);
			for (i = 0;i < slave->partitionsCount;i++) {
				read = fgets (line, LINE_MAX, sfile);
				if (read == NULL) {
					mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
					return -1;
				}	
				if ((line[0] == 'p') && (line[1] == 'p') && (line[2] == 'c')) {				
					mprintf(20, "created Partitions count\n", threadnum);
					t1 = strtok(line,"=");
					t1 = strtok(NULL,"=");
					slave->MOAPart[i].processed = atol (t1);
					sprintf(msg, "read slave->MOAPart[%ld].processed %d\n", i, slave->MOAPart[i].processed);
					mprintf(20, msg, threadnum);
				}
			  createMOAStruct (&slave->MOAPart[i].msaAlgn);
				read = fgets (line, LINE_MAX, sfile);
				if (read == NULL) {
					mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
					return -1;
				}	
				if ((line[0] == 'p') && (line[1] == 'd') && (line[2] == 'm')) {				
					t1 = strtok(line,"=");
					t1 = strtok(NULL,"=");
					slave->MOAPart[i].msaAlgn->dimn = atol (t1);
					sprintf(msg, "read slave->MOAPart[%ld].msaAlgn->dimn %ld\n", i, slave->MOAPart[i].msaAlgn->dimn);
					mprintf(20, msg, threadnum);
				  slave->MOAPart[i].msaAlgn->shape = (long * ) mcalloc (slave->MOAPart[i].msaAlgn->dimn, sizeof(long));
					slave->MOAPart[i].sequences = (char * *) mmalloc(slave->MOAPart[i].msaAlgn->dimn * sizeof(char *));
					for (j = 0;j < slave->MOAPart[i].msaAlgn->dimn;j++) {
						read = fgets (line, LINE_MAX, sfile);
						if (read == NULL) {
							mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
							return -1;
						}					
						slave->MOAPart[i].msaAlgn->shape[j] = atol (line);
						sprintf(msg, "read slave->MOAPart[%ld].msaAlgn->shape[%ld] %ld\n", i, j, slave->MOAPart[i].msaAlgn->shape[j]);
						mprintf(20, msg, threadnum);
						slave->MOAPart[i].sequences[j] = (char *) mmalloc(slave->MOAPart[i].msaAlgn->shape[j] * sizeof(char));
						for (k = 0;k < slave->MOAPart[i].msaAlgn->shape[j];k++) {
							read = fgets (line, LINE_MAX, sfile);
							if (read == NULL) {
								mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
								return -1;
							}					
							slave->MOAPart[i].sequences[j][k] = line[0];
							sprintf(msg, "read slave->MOAPart[%ld].sequences[%ld][%ld] %c\n", i, j, k, slave->MOAPart[i].sequences[j][k]);
							mprintf(20, msg, threadnum);
						}
					}
				}

				read = fgets (line, LINE_MAX, sfile);
				if (read == NULL) {
					mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
					return -1;
				}	
				if ((line[0] == 'p') && (line[1] == 'e') && (line[2] == 'l')) {				
					t1 = strtok(line,"=");
					t1 = strtok(NULL,"=");
					slave->MOAPart[i].msaAlgn->elements_ub = atol (t1);
					sprintf(msg, "read slave->MOAPart[%ld].msaAlgn->elements_ub %ld\n", i, slave->MOAPart[i].msaAlgn->elements_ub);
					mprintf(20, msg, threadnum);
				  slave->MOAPart[i].msaAlgn->elements = (MOA_elm * ) mcalloc (slave->MOAPart[i].msaAlgn->elements_ub, sizeof(MOA_elm));
				  slave->MOAPart[i].msaAlgn->indexes = (unsigned long * ) mcalloc (slave->MOAPart[i].msaAlgn->elements_ub, sizeof(unsigned long));
					for (j = 0;j < slave->MOAPart[i].msaAlgn->elements_ub;j++) {
						read = fgets (line, LINE_MAX, sfile);
						if (read == NULL) {
							mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
							return -1;
						}					
						slave->MOAPart[i].msaAlgn->elements[j].val = atol (line);
						sprintf(msg, "read slave->MOAPart[%ld].msaAlgn->elements[%ld].val %ld\n", i, j, slave->MOAPart[i].msaAlgn->elements[j].val);
						mprintf(20, msg, threadnum);

						read = fgets (line, LINE_MAX, sfile);
						if (read == NULL) {
							mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
							return -1;
						}					
						slave->MOAPart[i].msaAlgn->elements[j].prev_ub = atol (line);
						sprintf(msg, "read slave->MOAPart[%ld].msaAlgn->elements[%ld].prev_ub %ld\n", i, j, slave->MOAPart[i].msaAlgn->elements[j].prev_ub);
						mprintf(20, msg, threadnum);

						if (slave->MOAPart[i].msaAlgn->elements[j].prev_ub > 0) {
						slave->MOAPart[i].msaAlgn->elements[j].prev = (unsigned long *) mmalloc (slave->MOAPart[i].msaAlgn->elements[j].prev_ub * sizeof (unsigned long));

							for (k = 0; k < slave->MOAPart[i].msaAlgn->elements[j].prev_ub; k++)  {
								read = fgets (line, LINE_MAX, sfile);
								if (read == NULL) {
									mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
									return -1;
								}					
								slave->MOAPart[i].msaAlgn->elements[j].prev[k] = atol (line);
								sprintf(msg, "read slave->MOAPart[%ld].msaAlgn->elements[%ld].prev[%ld] %ld\n", i, j, k, slave->MOAPart[i].msaAlgn->elements[j].prev[k]);
								mprintf(20, msg, threadnum);
							}
						}
					}
					for (j = 0;j < slave->MOAPart[i].msaAlgn->elements_ub;j++) {
						read = fgets (line, LINE_MAX, sfile);
						if (read == NULL) {
							mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
							return -1;
						}					
						slave->MOAPart[i].msaAlgn->indexes[j] = atol (line);
						sprintf(msg, "read slave->MOAPart[%ld].msaAlgn->indexes[%ld] %ld\n", i, j, slave->MOAPart[i].msaAlgn->indexes[j]);
						mprintf(20, msg, threadnum);
					}
				}
			}
		}
		else if ((line[0] == 's') && (line[1] == 'q') && (line[2] == 'm')) {
			t1 = strtok(line,"=");
			t1 = strtok(NULL,"=");
			slave->seqNum = atol (t1);
			sprintf(msg, "read slave->seqNum %ld\n", slave->seqNum);
			mprintf(20, msg, threadnum);
			slave->seqLen = (long * ) mcalloc (slave->seqNum, sizeof(long));
			for (i=0;i<slave->seqNum;i++) {
				read = fgets (line, LINE_MAX, sfile);
				if (read == NULL) {
					mprintf(1, "Can not read checkpoint file, exiting.\n", threadnum);
					return -1;
				}					
				slave->seqLen[i] = atol (line);
				sprintf(msg, "read slave->seqLen[%ld] %ld\n", i, slave->seqLen[i]);
				mprintf(20, msg, threadnum);
			}
		}
		else if ((line[0] == 's') && (line[1] == 't') && (line[2] == 'p')) {
			t1 = strtok(line,"=");
			t1 = strtok(NULL,"=");
			slave->stype = atol (t1);
			sprintf(msg, "read stype %d\n", slave->stype);
			mprintf(20, msg, threadnum);
		}
		else if ((line[0] == 'c') && (line[1] == 'm') && (line[2] == 'p')) {
			t1 = strtok(line,"=");
			t1 = strtok(NULL,"=");
			slave->ComputationPhase = atol (t1);
			sprintf(msg, "read ComputationPhase %d\n", slave->ComputationPhase);
			mprintf(20, msg, threadnum);
		}
		else if ((line[0] == 'p') && (line[1] == 'i') && (line[2] == 'd')) {
			t1 = strtok(line,"=");
			t1 = strtok(NULL,"=");
			sProcid = atol (t1);
			sprintf(msg, "read sProcid %d\n", sProcid);
			mprintf(20, msg, threadnum);
		}
	}
  fclose(sfile);
		//if (t1 != NULL)
			//free(t1);
	return 0;
}
