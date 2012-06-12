/***********************************************************
* Author: Manal Helal                                      *
* Last Modification Date: Fri 12 Jan 2007 03:39:51 AM EST  *
* Project : MMSA - Multiple Sequence Alignment Based on 	  *
* 					Mathematics of Arrays - PhD Experimentation *
* File: utils.c, a library for global utility functions    *
* Function:
*		processArguments
*		readInputSequence
*     printSeq
*     print_OCout
*     to_proc_cells
*     print_outging_cells
*     print_OCin
*     PrintSequencies
*     PrintASeq
*     PrintOptimalPath
*     getTime
*     cpTime
*     isTimeDiffEquals
*     init_output
*     isEven
*     mprintf
*     getSizes
*     mmalloc
*     mcalloc
*     strrev
*     mpow
*     a_max
*     a_min
*     Factorial
*     l2Comb
*     maxVal
*     Combinations
*     file_copy
*     get_bytes
***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include <dirent.h>
#include "globals.h"
#include "moaDst.h"
#include "moa.h"
#include "utils.h"
#include "main.h"
#include "scores.h"

int processArguments (int argc, char ** argv, MOATypeDimn * seqNum, char * * * sequences, char * * * seqName, MOATypeShape * * seqLen, int * stype, long * partitionSize) {
  /*
    -a = Algorithm used, the default is 0) dynamic programming, otherwise 1) Sum-of-Pairs Scores, 2) Rubber Band Heuristic 
    -b = source file for all sequences in fasta format
    -c = number of sequences
    -s = Scoring Type, 0 means either default values or supplied values, other as implemented in scores.h and scores.c
    -l = local alignment, other wise the default is global alignment
    -m = maximum number of optimal paths to be retrieved
    -e = Epsilons Value for search space reduction. percentage of partitions to be done from all available in a wave
    -n = the type of Epsilon used, 0: default is fixed distance from middle on each dimension, 1: Percentage of partition from all partitions available in a wave
    -o = the output log filename (will be prefixed with the process name and suffixed with the processor number)
    -d = the debug level (what to be written in the og files)
    -p = the partition size
    -r = Restore from checkpoint files & resume from there
    -g = gapOpening to be used in scoring if scoring type = 1
    -j = gapExtension to be used in scoring if scoring type = 1
    -t = matchscore to be used in scoring if scoring type = 1
    -f = mismatchscore to be used in scoring if scoring type = 1
    -i = Border Cell Scores Initialization Method : default 0= IndexProduct, 1= IndexSum, 2= GapPenalty
    -x = Execution Method, 0=Sequential or 1=Distributed (default)    -X (debrecated Master/Slave Design) = Scheduling Method used, 0=Round Robin 1=Bag of Tasks.2= Dependency Based Scheduling
    -z =  Number of Processors if needed to passed as argument to sequential processing simulations
 - available k, q, u, v, w, y
   */
	int argRead = 1;
	long i;
	/* Initialize Arguments default Values*/
	(*seqNum) = 0;
	/* Epsilons Value is used to reduce the search space by scoring only the
	 * partitions indexed +/- Epsilons from m_i where m is the middle of the 		 * sequence i length (dimension i shape value), if zero, then all search 	  * space is used */
	Epsilons = 0;
        EpsilonsType = 0;
	(*stype) = 0;
	AlignmentType = Global;
        Algorithm = DP;
	SchedMethod = RR;
	maxAlignmentsNumber = 20;
	strcpy(outputfilename, "mm");
	pdebug = 0;
	RestoreFlag = 0;
	(*partitionSize) = 3;
	Mode = Distributed;
	gapOpenning = -4;
        gapExtension = -2;
	mismatchscore = -1;
	matchscore = 1;
        alignResolution = 100;
	initializationMethod = IndexProductInit;
	/* read & process MOAMSA Arguments */
	while (argRead < argc) {
		/*0. Read Scoring Algorithm */
		if (strcmp(argv[argRead],"-a") == 0) {
                    if (argc  <= argRead)  {
                        printf("More Arguments expected\n");
                        return -1;
                    }
                    if (argRead <argc) 
                        Algorithm = atoi(argv[++argRead]);
                    else {
                        printf("Expected Numeric Argument after -a for Scoring Algorithm, 1) Dynami Prigramming - Default, 2) Sum of Pairs. 3) Rubber Band Method.\n");
                        return -1;
                    }                     
                }
		/*0. Read SeqNum and Sequences File Names */
                else if (strcmp(argv[argRead],"-b") == 0) {
                        argRead++;
                        char * FastaFileName = NULL;
                        FastaFileName = mmalloc (((MOATypeInd) (LINE_MAX)) * ((MOATypeInd) sizeof *FastaFileName));
                        if (FastaFileName == NULL) {
                            printf("Couldn't allocate memory for FastaFileName. Exiting.\n");
                            return -1;
                        }
                        strcpy(FastaFileName, argv[argRead]);
                        (*seqNum) = readSequencesFastaFile(FastaFileName, sequences, seqName, seqLen);
                        if (FastaFileName != NULL)
                            free (FastaFileName);
                        FastaFileName = NULL;
		}
		/*1. Read SeqNum and Sequences File Names */
		else if (strcmp(argv[argRead],"-c") == 0) {
			(*seqNum) = atoi(argv[++argRead]);
			if ((*seqNum) <= 1) {
				printf("Sequence Num must be more than 1\n");
				return -1;
			}
			if (argc  <= ((*seqNum) + 3))  {
				printf("More Arguments expected\n");
				return -1;
			}
			(*sequences) = mmalloc (((MOATypeInd)(*seqNum)) * (MOATypeInd) sizeof *(*sequences));
			(*seqLen) = mcalloc ((MOATypeInd) (*seqNum), (MOATypeInd) sizeof *(*seqLen));
			(*seqName) = mmalloc (((MOATypeInd) (*seqNum)) * ((MOATypeInd) sizeof *(*seqName)));
			for (i=0;i<(*seqNum);i++) {	
				argRead++;
				(*seqName)[i] = mmalloc (((MOATypeInd) (LINE_MAX)) * ((MOATypeInd) sizeof *(*seqName)[i]));
				strcpy((*seqName)[i], argv[argRead]);
				(*seqLen)[i] = readInputSequence(&(*seqName)[i], &(*sequences)[i]);
				if ((*seqLen)[i] == -1) {
					printf("Sequence %ld File Could not open. Exiting.\n", i);
					return -1;
				}
				// printSeq (sequences[i], seqLen[i]);
			}
		}
		/*2. Read Scoring Type */
		else if (strcmp(argv[argRead], "-s") == 0) {      
			if (argc  <= argRead)  {
				printf("More Arguments expected\n");
				return -1;
			}
                        if (argRead <argc) {
                            (*stype) = atoi(argv[++argRead]);
                            SetScoreMatrix(NULL, (*stype));
                        }
			else {
				printf("Expected Numeric Argument after -s for Scoring Type required.\n");
				return -1;
			} 
		}
		/* 3. Decide if Local or Global Alignment */
		else if (strcmp(argv[argRead], "-l") == 0) {
			AlignmentType = Local;
		}
		else if (strcmp(argv[argRead],"-r") == 0) {
                    RestoreFlag = 1;
		}    
		/* 4. Know Maximum Alignments Required to be returned */
		else if (strcmp(argv[argRead],"-m") == 0) {
			if (argc  < argRead)  {
				printf("More Arguments expected\n");
				return -1;
			}
			if (argRead <argc) 
				maxAlignmentsNumber = atoi(argv[++argRead]);
			else {
				printf("Expected Numeric Argument after -m for maximum Alignemnts required.\n");
				return -1;
			} 
		}
		/* 5. Decide the Epsilons Value for search space reduction */
		else if (strcmp(argv[argRead],"-e") == 0) {
			if (argc  < argRead)  {
				printf("More Arguments expected\n");
				return -1;
			}
			if (argRead <argc) 
				Epsilons= atoi(argv[++argRead]);
			else {
				printf("Expected Numeric Argument after -e for Epsilons Value for percentage of partitions to be done in a wave.\n");
				return -1;
			} 
		}
		/* 6. Decide the number of partitions to be done in each wave */
		else if (strcmp(argv[argRead],"-n") == 0) {
			if (argc  < argRead)  {
				printf("More Arguments expected\n");
				return -1;
			}
			if (argRead <argc) 
				EpsilonsType= atoi(argv[++argRead]);
			else {
				printf("Expected Numeric Argument after -n for type of Epsilon, 0: Default means fixed distance, 1: Percentage of paritions to be scored from all available in a wave of partitions.\n");
				return -1;
			} 
		}
		/* 7. Decide the output filename */
		else if (strcmp(argv[argRead], "-o") == 0) {
			if (argc  <= argRead)  {
				printf("More Arguments expected\n");
				return -1;
			}
			if (argRead <argc) 
				strcpy (outputfilename, argv[++argRead]);
			else {
				printf("Expected Numeric Argument after -o for outputfilename.\n");
				return -1;
			} 
		}
		/* 8. Decide the debug Level */
		else if (strcmp(argv[argRead], "-d") == 0) {
			if (argc  < argRead)  {
				printf("More Arguments expected\n");
				return -1;
			}
			if (argRead <argc) 
				pdebug = atoi(argv[++argRead]);
			else {
				printf("Expected Numeric Argument after -d for debug level required.\n");
				return -1;
			} 
		}
		/* 9. Decide the partition size - for the distributed version */
		else if (strcmp(argv[argRead], "-p") == 0) {
			if (argc  < argRead)  {
				printf("More Arguments expected\n");
				return -1;
			}
			if (argRead <argc) 
				(*partitionSize)  = atoi(argv[++argRead]);
			else {
				printf("Expected Numeric Argument after -p for partitions Size.\n");
				return -1;
			} 
		}
		/* 10. Decide the gapPenalty */
		else if (strcmp(argv[argRead], "-g") == 0) {
			if (argc  < argRead)  {
				printf("More Arguments expected\n");
				return -1;
			}
			if (argRead <argc) 
				gapOpenning  = atoi(argv[++argRead]);
			else {
				printf("Expected Numeric Argument after -g for gap Openning Penalty.\n");
				return -1;
			} 
		}
		/* 11. Decide the gapPenalty */
		else if (strcmp(argv[argRead], "-j") == 0) {
			if (argc  < argRead)  {
				printf("More Arguments expected\n");
				return -1;
			}
			if (argRead <argc) 
				gapExtension  = atoi(argv[++argRead]);
			else {
				printf("Expected Numeric Argument after -j for gap Extension Penalty.\n");
				return -1;
			} 
		}
		/* 11. Decide the match score */
		else if (strcmp(argv[argRead], "-t") == 0) {
			if (argc  < argRead)  {
				printf("More Arguments expected\n");
				return -1;
			}
			if (argRead <argc) 
				matchscore  = atoi(argv[++argRead]);
			else {
				printf("Expected Numeric Argument after -t for match score.\n");
				return -1;
			} 
		}
		/* 12. Decide the mismatch score */
		else if (strcmp(argv[argRead], "-f") == 0) {
			if (argc  < argRead)  {
				printf("More Arguments expected\n");
				return -1;
			}
			if (argRead <argc) 
				mismatchscore  = atoi(argv[++argRead]);
			else {
				printf("Expected Numeric Argument after -i for mis-match score.\n");
				return -1;
			} 
		}
		/* 13. Decide the initialization Method */
		else if (strcmp(argv[argRead], "-i") == 0) {
			if (argc  < argRead)  {
				printf("More Arguments expected\n");
				return -1;
			}
			if (argRead <argc) {
				initializationMethod  = atoi(argv[++argRead]);
				if (initializationMethod > 2) {
					printf("Error reading initialization Method value. initialization Method should be 0 for Index Product (* Gap Penalty), 1 for Index Sum (* Gap Penalty), or 2 for Gap Penalty. Exiting.\n");
					return -1;
				}
			}
			else {
				printf("Expected Numeric Argument after -i for mis-match score.\n");
				return -1;
			} 
		}
		/* 14. Decide the Execution Method (Distributed / Sequential) */
		else if (strcasecmp(argv[argRead], "-x") == 0) {
			if (argc  < argRead)  {
				printf("More Arguments expected\n");
				return -1;
			}
			if (argRead <argc)  {
				Mode = atoi(argv[++argRead]); /*should 0 or 1 only*/
				if (Mode > 1) {
					printf("Error reading execution mode value. Mode should either be 0 for Distributed or 1 for Sequential. Exiting.\n");
					return -1;
				}
			}
				
			else {
				printf("Expected Numeric Argument after -x for Execution Mode.\n");
				return -1;
			} 
		}		
		/* 15. Decide the Sceduling Method */
		else if (strcasecmp(argv[argRead], "-X") == 0) {
			if (argc  < argRead)  {
				printf("More Arguments expected\n");
				return -1;
			}
			if (argRead <argc) 
				SchedMethod = atoi(argv[++argRead]);
			else {
				printf("Expected Numeric Argument after -X for scheduling method required.\n");
				return -1;
			} 
		}
		/* 16.  number of processors for simulations*/
		else if (strcasecmp(argv[argRead], "-z") == 0) {
			if (argc  < argRead)  {
				printf("More Arguments expected\n");
				return -1;
			}
			if (argRead <argc) 
				ClusterSize = atoi(argv[++argRead]);
			else {
				printf("Expected Numeric Argument after -z for number of processors for simulations.\n");
				return -1;
			} 
		}
		/* 17. Display Help - available k, q, u, v, w, y*/
		else if (strcasecmp(argv[argRead], "-h") == 0) {
                    printf("mpirun -np [Number of Processors] [mpirun arguments] [mmDst executable Name] [mmDst Executable arguments from the following].\n");
                    printf("-a = Algorithm used, the default is 0) dynamic programming, otherwise 1) Sum-of-Pairs Scores, 2) Rubber Band Heuristic.\n");
                    printf("-b = [Mandatory if -c is not used] source file for all sequences in fasta format.\n");
                    printf("-c = [Mandatory if -b is not used] number of sequences.\n");
                    printf("-s = Scoring Type, Default is 0 meaning either default values or supplied values in -g -j -t -f, other as implemented in scores.h and scores.c.\n");
                    printf("-l = local alignment, other wise the default is global alignment.\n");
                    printf("-m = maximum number of optimal paths to be retrieved.\n");
                    printf("-e = Epsilons Value for search space reduction. percentage of partitions to be done from all available in a wave.\n");
                    printf("-n = the type of Epsilon used, 0: default is fixed distance from middle on each dimension, 1: Percentage of partition from all partitions available in a wave.\n");
                    printf("-o = the output log filename (will be prefixed with the process name and suffixed with the processor number).\n");
                    printf("-d = the debug level (what to be written in the og files).\n");
                    printf("-p = the partition size.\n");
                    printf("-r = Restore from checkpoint files & resume from there.\n");
                    printf("-g = gapOpening to be used in scoring if scoring type = 1.\n");
                    printf("-j = gapExtension to be used in scoring if scoring type = 1.\n");
                    printf("-t = matchscore to be used in scoring if scoring type = 1.\n");
                    printf("-f = mismatchscore to be used in scoring if scoring type = 1.\n");
                    printf("-i = Border Cell Scores Initialization Method : default 0= IndexProduct, 1= IndexSum, 2= GapPenalty.\n");
                    printf("-x = Execution Method, 0=Sequential or 1=Distributed (default)    -X (debrecated Master/Slave Design) = Scheduling Method used, 0=Round Robin 1=Bag of Tasks.2= Dependency Based Scheduling.\n");
                    printf("-z =  Number of Processors if needed to passed as argument to sequential processing simulations.\n");
                    printf("-h =  Displays this help screen.\n");
                    return -1;
		}
		argRead ++;
	}
	return 0;
}

MOATypeShape readInputSequence (char * * fileName, char * * sequence) {
  char line[LINE_MAX];
  char elm, msg[SHORT_MESSAGE_SIZE];
  MOATypeShape i = 0;
  FILE * f = fopen ((*fileName), "r");
  if (f == NULL) {
    sprintf(msg, "File %s could not be opened!\n", (*fileName));
    mprintf (1, msg, 1);
    return -1;
  }
  (*sequence) = NULL; 
  (*sequence) = mmalloc (((MOATypeInd) 2) * ((MOATypeInd) sizeof *(*sequence)));      
  if ((*sequence) == NULL) {
    sprintf(msg, "Could not allocate memory for sequence %s!\n", (*fileName));
    mprintf (1, msg, 1);
    fclose(f);
    return -1;
  }
  //read the sequence Name first line
  fgets (line, LINE_MAX, f);
  strcpy((*fileName), line);
  (*fileName)[strlen((*fileName))-1] = '\0';
  // read character by character.
  i++;
  (*sequence)[0]= GAPCHAR;
  while (!feof(f)) {
    elm = fgetc(f);
    if (elm != EOF)  {
      if (((elm >= 'A') && (elm <= 'Z')) ||((elm >= 'a') && (elm <= 'z'))) {
	i++;
	(*sequence) = realloc ((*sequence), ((MOATypeInd) sizeof(char)) * ((MOATypeInd) (i+1)));
	if ((*sequence) == NULL) {
	  sprintf(msg, "Could not allocate memory for sequence %s!\n", (*fileName));
	  mprintf (1, msg, 1);
	  fclose(f);
	  return -1;
	}
	else {
	  (*sequence)[i-1] = toupper(elm);
	}
      }
    }
  }  
  (*sequence)[i] = '\0';
  fclose(f);
  return i;
}

MOATypeDimn readSequencesFastaFile (char * sfilename, char * * * sequences, char * * * seqName, MOATypeShape * * seqLen) {
  char line[LINE_MAX];
  char elm, msg[SHORT_MESSAGE_SIZE];
  MOATypeShape i, a, s = 0;
  MOATypeDimn seqNum;
  FILE * f = fopen (sfilename, "r");
  if (f == NULL) {
    sprintf(msg, "File %s could not be opened!\n", sfilename);
    mprintf (1, msg, 1);
    return -1;
  }
    (*seqLen) = NULL;
    (*sequences) = NULL;    
    (*seqName) = NULL;
    seqNum = 0;
    while (!feof(f)) {
    //read the sequence Name first line
    if (fgets (line, LINE_MAX, f) == NULL) {
        fclose(f);
        return seqNum;
    }         
    if (strcmp (line, "\n") != 0) {
        seqNum++;
        (*seqName) = realloc ((*seqName), (MOATypeInd) (seqNum) *((MOATypeInd) sizeof *(*seqName))); 
        if ((*seqName) == NULL) {
                sprintf(msg, "Could not allocate memory for aligned sequence %s!\n", sfilename);
                mprintf (1, msg, 1);
                fclose(f);
                return -1;
        }
        (*seqName)[seqNum-1] = NULL;
        (*seqName)[seqNum-1] = mmalloc (((MOATypeInd) LINE_MAX) * ((MOATypeInd) sizeof *(*seqName)[seqNum-1]));
        if ((*seqName)[seqNum-1] == NULL) {
                sprintf(msg, "Could not allocate memory for aligned sequence %s!\n", sfilename);
                mprintf (1, msg, 1);
                fclose(f);
                return -1;
        }
        //read the sequence Name first line
        strcpy((*seqName)[seqNum-1], line);
        (*seqName)[seqNum-1][strlen((*seqName)[seqNum-1])-1] = '\0';
        (*seqLen)  = realloc ((*seqLen), (MOATypeInd) (seqNum) * (MOATypeInd) sizeof *(*seqLen));
        if ((*seqLen) == NULL) {
                sprintf(msg, "Could not allocate memory for aligned sequence %s!\n", sfilename);
                mprintf (1, msg, 1);
                fclose(f);
                return -1;
        }
        (*seqLen)[seqNum-1] = 0;
        (*sequences) = realloc ((*sequences), ((MOATypeInd) seqNum) * ((MOATypeInd) sizeof *(*sequences)));
        if ((*sequences) == NULL) {
                sprintf(msg, "Could not allocate memory for aligned sequence %s!\n", sfilename);
                mprintf (1, msg, 1);
                fclose(f);
                return -1;
        }
        (*sequences)[seqNum-1] = NULL;
        // read character by character from the end of the line.
        (*sequences)[seqNum-1] = mmalloc ((MOATypeInd) sizeof(char));
        if ((*sequences)[seqNum-1] == NULL) {
            sprintf(msg, "Could not allocate memory for sequence %s!\n", sfilename);
            mprintf (1, msg, 1);
            fclose(f);
            return -1;
        }
        (*sequences)[seqNum-1][0]= GAPCHAR;
        s = 1;
      elm = fgetc(f);
      while ((elm != '\0' ) && (elm != '>') && (elm != ' ') && (elm != EOF)) {
        if(((elm >= 'A') && (elm <= 'Z')) ||((elm >= 'a') && (elm <= 'z'))) {
            s++;
            (*sequences)[seqNum-1] = realloc ((*sequences)[seqNum-1], ((MOATypeInd) (s+1)) * ((MOATypeInd) sizeof *((*sequences)[seqNum-1])));
            if ((*sequences)[seqNum-1] == NULL) {
              sprintf(msg, "Could not allocate memory for sequence %lld in file %s!\n", seqNum, sfilename);
              mprintf (1, msg, 1);
              fclose(f);
              return -1;
            }
            else {
              (*sequences)[seqNum-1][s-1] = toupper(elm);
            }
        }
        elm = fgetc(f);
      } // end of Sequence line
        (*sequences)[seqNum-1][s] = '\0';
        (*seqLen)[seqNum-1] = s;
    }
  } // end of file
  fclose(f);
 
  return seqNum;
}

int outputFastaFormat (char * fileName, char * sequence, MOATypeShape seqLen) {
    MOATypeShape i;
    chdir("out");
    FILE * outfile = fopen (fileName, "w");
    if (outfile == NULL) {
        printf("File %s could not be opened!\n", fileName);
        return -1;
    }
    fprintf(outfile,">%s\n", fileName);
    
    for (i=1;i<seqLen;i++) {
        fprintf(outfile,"%c", sequence[i]);
    }
    fclose(outfile);
    chdir("../");
    return 0;
}

int readAndConvertToFasta (char * folderName) {
    FILE *fp;
    DIR *dirp;
    struct dirent *dp;
    struct stat buf;
    char * seqName = NULL;
    char * seqFileName = NULL;
    char * tempToken = NULL;
    char * sequence = NULL;
    char * dot;
    char *token;
    size_t nameLen, tempLen;
    MOATypeShape seqLen;
    
    
    struct dirent **namelist;
    int n;

    chdir(folderName);
    n = scandir(folderName, &namelist, 0, alphasort);
    if (n < 0)
        perror("scandir");
    else {
        seqName = mmalloc (((MOATypeInd) LINE_MAX) * ((MOATypeInd) sizeof *seqName));
        if (seqName == NULL) {
            printf ("Couldn't Allocate memory to sequence Namse in readAndConvertToFasta function. Exiting.\n");
            return -1;
        }
        seqFileName = mmalloc (((MOATypeInd) LINE_MAX) * ((MOATypeInd) sizeof *seqFileName));
        if (seqFileName == NULL) {
            printf ("Couldn't Allocate memory to sequence Namse in readAndConvertToFasta function. Exiting.\n");
            return -1;
        }
        tempToken = mmalloc (((MOATypeInd) LINE_MAX) * ((MOATypeInd) sizeof *tempToken));
        if (tempToken == NULL) {
            printf ("Couldn't Allocate memory to sequence Namse in readAndConvertToFasta function. Exiting.\n");
            return -1;
        }
        while(n--) {
            /*if (stat(namelist[n]->d_name, &buf) == -1)
                perror("stat\n");

            if (S_ISDIR(buf.st_mode)) {
                printf("%s is a directory\n", namelist[n]->d_name);
                if ((strcmp(namelist[n]->d_name,".") != 0) && (strcmp(namelist[n]->d_name,"..") != 0))
                    readAndConvertToFasta(namelist[n]->d_name);
            }
            else if  (S_ISREG(buf.st_mode)) {*/
                nameLen = strlen(namelist[n]->d_name);
                seqName[0] = '\0';
                seqFileName[0] = '\0';
                
                strncpy(seqName, namelist[n]->d_name, nameLen);
                seqName[nameLen]  = '\0';

                dot = strrchr(namelist[n]->d_name, '.');
                if (dot != NULL)  {
                  if ((token = strtok(namelist[n]->d_name, ".")) != NULL) {
                      strcpy(tempToken, token);                          
                      seqFileName = "";
                      tempLen = strlen(token);
                      nameLen = tempLen;
                      while ((token = strtok(NULL, ".")) != NULL) {
                          strcat(seqFileName, tempToken);                              
                          strcpy(tempToken, token);                          
                          tempLen = strlen(token);
                          nameLen += tempLen;
                          if ((token = strtok(NULL, ".")) != NULL) {
                              strcat(seqFileName, tempToken);                              
                              strcpy(tempToken, token);                          
                              tempLen = strlen(token);
                              nameLen += tempLen;
                          }
                      }
                      nameLen -= tempLen;
                }
                
                    //seqFileName = strtok(namelist[n]->d_name, ".");                
                    //nameLen -= 4;
                }
                else 
                    strncpy(seqFileName, namelist[n]->d_name, nameLen);                    
                
                seqFileName[nameLen]  = '\0';
                if (strcmp(seqFileName, "Nocardia_carnea_AF430036") == 0)
                    printf ("Here is the bug\n");
                seqLen = readInputSequence (&seqName, &sequence);
                outputFastaFormat (seqFileName, sequence, seqLen);
                printf("%s\t\t%lld\n", seqFileName, seqLen);
                if (sequence != NULL)
                    free (sequence);
                sequence = NULL;
            //}
            free(namelist[n]);            
        }
        free(namelist);
        if (seqName != NULL)
            free (seqName);
        seqName = NULL;
        if (seqFileName != NULL)
            free (seqFileName);
        seqFileName = NULL;
        if (tempToken != NULL)
            free (tempToken);
        tempToken = NULL;
    }    
    /*
    dirp = opendir(folderName);
    chdir(folderName);

     /* Look at each entry in turn *
     while ((dp = readdir(dirp)) != NULL) {

        /* Now stat the file to get more information *
        if (stat(dp->d_name, &buf) == -1)
            perror("stat\n");

        if (S_ISDIR(buf.st_mode)) {
            printf("%s is a directory\n", dp->d_name);
            if ((strcmp(dp->d_name,".") != 0) && (strcmp(dp->d_name,"..") != 0))
                readAndConvertToFasta(dp->d_name);
        }
        else if  (S_ISREG(buf.st_mode)) {
            printf("Processing Sequence %s ..... ", dp->d_name);
            seqName = mmalloc (((MOATypeInd) (LINE_MAX)) * ((MOATypeInd) sizeof *seqName));
            seqFileName = mmalloc (((MOATypeInd) (LINE_MAX)) * ((MOATypeInd) sizeof *seqFileName));
            if (seqName == NULL) {
                printf ("Couldn't Allocate memory to sequence Namse in readAndConvertToFasta function. Exiting.\n");
                return -1;
            }
            nameLen = strlen(dp->d_name);
            strcpy(seqName, dp->d_name);
            strncpy(seqFileName, dp->d_name, nameLen);
            seqLen = readInputSequence (&seqName, &sequence);
            outputFastaFormat (seqFileName, sequence, seqLen);
            printf("length: %lld .... done.\n", seqLen);
            if (seqName != NULL)
                free (seqName);
            seqName = NULL;
            if (seqFileName != NULL)
                free (seqFileName);
            seqFileName = NULL;
            if (sequence != NULL)
                free (sequence);
            sequence = NULL;
        }
     }

     (void) closedir(dirp);  
     */ 
     return 0;
}

int readFastaAlignment (char * sfilename, ProcessData * pData, TracebackData * tbData) {
  char line[LINE_MAX];
  char elm, msg[SHORT_MESSAGE_SIZE];
  MOATypeShape i, a, s = 0;
  FILE * f = fopen (sfilename, "r");
  if (f == NULL) {
    sprintf(msg, "File %s could not be opened!\n", sfilename);
    mprintf (1, msg, 1);
    return -1;
  }
    pData->seqLen = NULL;
    pData->sequences = NULL;
    tbData->completePath = NULL;
    pData->seqName = NULL;
    tbData->comseqLen = 0;
    pData->seqNum = 0;
    while (!feof(f)) {
    //read the sequence Name first line
    if (fgets (line, LINE_MAX, f) == NULL) {
        fclose(f);
        return pData->seqNum;
    }
    
        
    if (strcmp (line, "\n") != 0) {
        pData->seqNum++;
        tbData->seqNum = pData->seqNum;
        pData->seqName = realloc (pData->seqName, (MOATypeInd) (pData->seqNum) *((MOATypeInd) sizeof *(pData->seqName))); 
        if (pData->seqName == NULL) {
                sprintf(msg, "Could not allocate memory for aligned sequence %s!\n", sfilename);
                mprintf (1, msg, 1);
                fclose(f);
                return -1;
        }
        pData->seqName[pData->seqNum-1] = NULL;
        pData->seqName[pData->seqNum-1] = mmalloc (((MOATypeInd) (LINE_MAX)) * ((MOATypeInd) sizeof *(pData->seqName[pData->seqNum-1])));
        if (pData->seqName[pData->seqNum-1] == NULL) {
                sprintf(msg, "Could not allocate memory for aligned sequence %s!\n", sfilename);
                mprintf (1, msg, 1);
                fclose(f);
                return -1;
        }
        //read the sequence Name first line
        strcpy(pData->seqName[pData->seqNum-1], line);
        pData->seqName[pData->seqNum-1][strlen(pData->seqName[pData->seqNum-1])-1] = '\0';
        pData->seqLen  = realloc (pData->seqLen, (MOATypeInd) (pData->seqNum) * (MOATypeInd) sizeof *(pData->seqLen));
        if (pData->seqLen == NULL) {
                sprintf(msg, "Could not allocate memory for aligned sequence %s!\n", sfilename);
                mprintf (1, msg, 1);
                fclose(f);
                return -1;
        }
        pData->seqLen[pData->seqNum-1] = 0;
        pData->sequences = realloc (pData->sequences, ((MOATypeInd) (pData->seqNum)) * ((MOATypeInd) sizeof *(pData->sequences)));
        if (pData->sequences == NULL) {
                sprintf(msg, "Could not allocate memory for aligned sequence %s!\n", sfilename);
                mprintf (1, msg, 1);
                fclose(f);
                return -1;
        }
        pData->sequences[pData->seqNum-1] = NULL;
        tbData->completePath = realloc (tbData->completePath, (MOATypeInd) (pData->seqNum) * (MOATypeInd) sizeof *(tbData->completePath));
        if (tbData->completePath == NULL) {
                sprintf(msg, "Could not allocate memory for aligned sequence %s!\n", sfilename);
                mprintf (1, msg, 1);
                fclose(f);
                return -1;
        }
        tbData->completePath[pData->seqNum-1] = NULL;
        // read character by character from the end of the line.
        pData->sequences[pData->seqNum-1] = mmalloc ((MOATypeInd) sizeof(char));
        if (pData->sequences[pData->seqNum-1] == NULL) {
            sprintf(msg, "Could not allocate memory for sequence %s!\n", sfilename);
            mprintf (1, msg, 1);
            fclose(f);
            return -1;
        }
        pData->sequences[pData->seqNum-1][0]= GAPCHAR;
        a = 0;
        s = 1;
      elm = fgetc(f);
      while ((elm != '\0' ) && (elm != '\n') && (elm != ' ') && (elm != EOF)) {

          if (((elm >= 'A') && (elm <= 'Z')) ||((elm >= 'a') && (elm <= 'z')) || (elm == GAPCHAR)) {
          a++;


                if(((elm >= 'A') && (elm <= 'Z')) ||((elm >= 'a') && (elm <= 'z'))) {
                    s++;
                    pData->sequences[pData->seqNum-1] = realloc (pData->sequences[pData->seqNum-1], ((MOATypeInd) (s+1)) * ((MOATypeInd) sizeof *(pData->sequences[pData->seqNum-1])));
                    if (pData->sequences[pData->seqNum-1] == NULL) {
                      sprintf(msg, "Could not allocate memory for sequence %s!\n", sfilename);
                      mprintf (1, msg, 1);
                      fclose(f);
                      return -1;
                    }
                    else {
                      pData->sequences[pData->seqNum-1][s-1] = toupper(elm);
                    }
                }
                tbData->completePath[pData->seqNum-1] = realloc (tbData->completePath[pData->seqNum-1], ((MOATypeInd) (a+1)) * ((MOATypeInd)  sizeof *(tbData->completePath[pData->seqNum-1])));
                if (tbData->completePath[pData->seqNum-1] == NULL) {
                        sprintf(msg, "Could not allocate memory for aligned sequence %s!\n", sfilename);
                        mprintf (1, msg, 1);
                        fclose(f);
                        return -1;
                }
                else {
                        tbData->completePath[pData->seqNum-1][a-1] = toupper(elm);
                }

         }
        elm = fgetc(f);
      } // end of Sequence line
        pData->sequences[pData->seqNum-1][s] = '\0';
        tbData->completePath[pData->seqNum-1][a] = '\0';
        pData->seqLen[pData->seqNum-1] = s;
        if ((pData->seqNum > 1) && (tbData->comseqLen != a))
            printf ("Faulty Alignment File, aligned sequences are not of equal length\n.");
        else
            tbData->comseqLen = a;
    }
  } // end of file
  fclose(f);
 
  return 0;
}

MOATypeElmVal calcAlignmentColumnScore (ProcessData * pData, TracebackData * tbData, MOATypeInd k) {
    int newGap = 1;
    MOATypeDimn i, j;
    tbData->alignColScore[k] = 0;
    for (i=0;i<tbData->seqNum-1;i++) {
            for (j=i+1;j<tbData->seqNum;j++) {
                if (((tbData->completePath[i][k] != GAPCHAR) && (tbData->completePath[j][k] == GAPCHAR)) || ((tbData->completePath[i][k] == GAPCHAR) && (tbData->completePath[j][k] != GAPCHAR))) {
                  tbData->alignColScore[k] -= 1;
                }
                else if ((tbData->completePath[i][k] == tbData->completePath[j][k]) && (tbData->completePath[i][k] != GAPCHAR)) 
                  tbData->alignColScore[k] += 1;
            }
    }

    return tbData->alignColScore[k];
}
long double getMatchesCount (TracebackData * tbData, MOATypeDimn seq1, MOATypeDimn seq2) {
    long double matches = 0;
    MOATypeInd i;
    for (i=0;i<tbData->comseqLen;i++) {
        if (tbData->completePath[seq1][i] == tbData->completePath[seq2][i]) 
            matches ++;
    }
    return matches;
}

long double getGapCount (TracebackData * tbData, MOATypeDimn seq) {
    long double gapCount = 0;
    MOATypeInd i;
    for (i=0;i<tbData->comseqLen;i++) {
        if (tbData->completePath[seq][i] == GAPCHAR) 
            gapCount ++;
    }
    return gapCount;
}

long double calcAlignmentScores (ProcessData * pData, TracebackData * tbData) {
    MOATypeInd iter;
    MOATypeDimn kiter, liter;

    tbData->a = tbData->c = tbData->g = tbData->t = tbData->indel = tbData->r = tbData->n = tbData->d = tbData->q = tbData->e = tbData->h = tbData->i = tbData->l = tbData->k = tbData->m = tbData->f = tbData->p = tbData->s = tbData->w = tbData->y = tbData->v = tbData->b = tbData->z = tbData->x = 0;
    tbData->alignShannonEntropy = 0;
    if (tbData->alignDistaneMeasures != NULL) {
        for (kiter=0;kiter<tbData->seqNum;kiter++) {
            if (tbData->alignDistaneMeasures[kiter] != NULL) 
                free (tbData->alignDistaneMeasures[kiter]);
            tbData->alignDistaneMeasures[kiter] = NULL;
        }
        free (tbData->alignDistaneMeasures);
    }
    tbData->alignDistaneMeasures = NULL;
    tbData->alignDistaneMeasures = mmalloc (((MOATypeInd) tbData->seqNum) * ((MOATypeInd) sizeof *(tbData->alignDistaneMeasures)));
    if (tbData->alignGapScore != NULL) 
        free (tbData->alignGapScore);
    tbData->alignGapScore = NULL;
    tbData->alignGapScore = mmalloc (((MOATypeInd) tbData->seqNum) * ((MOATypeInd) sizeof *(tbData->alignGapScore)));
    for (kiter=0;kiter<tbData->seqNum;kiter++) {
        for (iter=0;iter<tbData->comseqLen;iter++) {
            if (((tbData->completePath[kiter][iter] >= 'A') && (tbData->completePath[kiter][iter] <= 'Z')) ||((tbData->completePath[kiter][iter] >= 'a') && (tbData->completePath[kiter][iter] <= 'z')) || (tbData->completePath[kiter][iter] == GAPCHAR)) {
                if (tbData->completePath[kiter][iter] == 'a' || tbData->completePath[kiter][iter] == 'A') 
                    tbData->a ++;
                if (tbData->completePath[kiter][iter] == 'c' || tbData->completePath[kiter][iter] == 'C') 
                    tbData->c ++;
                if (tbData->completePath[kiter][iter] == 'g' || tbData->completePath[kiter][iter] == 'G') 
                    tbData->g ++;
                if (tbData->completePath[kiter][iter] == 't' || tbData->completePath[kiter][iter] == 'T') tbData->t ++;
                if (tbData->completePath[kiter][iter] == GAPCHAR) tbData->indel ++;
                if (tbData->completePath[kiter][iter] == 'r' || tbData->completePath[kiter][iter] == 'R')	tbData->r++;
                if (tbData->completePath[kiter][iter] == 'n' || tbData->completePath[kiter][iter] == 'N')	tbData->n++;
                if (tbData->completePath[kiter][iter] == 'd' || tbData->completePath[kiter][iter] == 'D')	tbData->d++;
                if (tbData->completePath[kiter][iter] == 'q' || tbData->completePath[kiter][iter] == 'Q')	tbData->q++;
                if (tbData->completePath[kiter][iter] == 'e' || tbData->completePath[kiter][iter] == 'E')	tbData->e++;
                if (tbData->completePath[kiter][iter] == 'i' || tbData->completePath[kiter][iter] == 'I')	tbData->i++;
                if (tbData->completePath[kiter][iter] == 'h' || tbData->completePath[kiter][iter] == 'H')	tbData->h++;
                if (tbData->completePath[kiter][iter] == 'l' || tbData->completePath[kiter][iter] == 'L')	tbData->l++;
                if (tbData->completePath[kiter][iter] == 'k' || tbData->completePath[kiter][iter] == 'K')	tbData->k++;
                if (tbData->completePath[kiter][iter] == 'm' || tbData->completePath[kiter][iter] == 'M')	tbData->m++;
                if (tbData->completePath[kiter][iter] == 'f' || tbData->completePath[kiter][iter] == 'F')	tbData->f++;
                if (tbData->completePath[kiter][iter] == 'p' || tbData->completePath[kiter][iter] == 'P')	tbData->p++;
                if (tbData->completePath[kiter][iter] == 's' || tbData->completePath[kiter][iter] == 'S')	tbData->s++;
                if (tbData->completePath[kiter][iter] == 'w' || tbData->completePath[kiter][iter] == 'W')	tbData->w++;
                if (tbData->completePath[kiter][iter] == 'y' || tbData->completePath[kiter][iter] == 'Y')	tbData->y++;
                if (tbData->completePath[kiter][iter] == 'v' || tbData->completePath[kiter][iter] == 'V')	tbData->v++;
                if (tbData->completePath[kiter][iter] == 'b' || tbData->completePath[kiter][iter] == 'B')	tbData->b++;
                if (tbData->completePath[kiter][iter] == 'z' || tbData->completePath[kiter][iter] == 'Z')	tbData->z++;
                if (tbData->completePath[kiter][iter] == 'x' || tbData->completePath[kiter][iter] == 'X')	tbData->x++;
            }
        }
        tbData->alignDistaneMeasures[kiter] = mcalloc (((MOATypeInd) tbData->seqNum), ((MOATypeInd) sizeof *tbData->alignDistaneMeasures[kiter]));
        for (liter=kiter+1;liter<tbData->seqNum;liter++) {
            long double matches = getMatchesCount(tbData, kiter, liter);
            tbData->alignDistaneMeasures[kiter][liter] = matches / sqrtl ((long double) (pData->seqLen[kiter] - 1) * (pData->seqLen[liter] - 1));
        }
        long double gapCount = getGapCount(tbData, kiter);
        tbData->alignGapScore[kiter] = gapCount / sqrtl ((long double) (pData->seqLen[kiter] - 1));
    }
    if (tbData->a > 0) tbData->alignShannonEntropy += tbData->a*logl(tbData->a);
    if (tbData->c > 0) tbData->alignShannonEntropy += tbData->c*logl(tbData->c);
    if (tbData->g > 0) tbData->alignShannonEntropy += tbData->g*logl(tbData->g);
    if (tbData->t > 0) tbData->alignShannonEntropy += tbData->t*logl(tbData->t);
    if (tbData->r > 0) tbData->alignShannonEntropy += tbData->r*logl(tbData->r);
    if (tbData->n > 0) tbData->alignShannonEntropy += tbData->n*logl(tbData->n);
    if (tbData->d > 0) tbData->alignShannonEntropy += tbData->d*logl(tbData->d);
    if (tbData->q > 0) tbData->alignShannonEntropy += tbData->q*logl(tbData->q);
    if (tbData->e > 0) tbData->alignShannonEntropy += tbData->e*logl(tbData->e);
    if (tbData->i > 0) tbData->alignShannonEntropy += tbData->i*logl(tbData->i);
    if (tbData->h > 0) tbData->alignShannonEntropy += tbData->h*logl(tbData->h);
    if (tbData->l > 0) tbData->alignShannonEntropy += tbData->l*logl(tbData->l);
    if (tbData->k > 0) tbData->alignShannonEntropy += tbData->k*logl(tbData->k);
    if (tbData->m > 0) tbData->alignShannonEntropy += tbData->m*logl(tbData->m);
    if (tbData->f > 0) tbData->alignShannonEntropy += tbData->f*logl(tbData->f);
    if (tbData->p > 0) tbData->alignShannonEntropy += tbData->p*logl(tbData->p);
    if (tbData->s > 0) tbData->alignShannonEntropy += tbData->s*logl(tbData->s);
    if (tbData->w > 0) tbData->alignShannonEntropy += tbData->w*logl(tbData->w);
    if (tbData->y > 0) tbData->alignShannonEntropy += tbData->y*logl(tbData->y);
    if (tbData->v > 0) tbData->alignShannonEntropy += tbData->v*logl(tbData->v);
    if (tbData->b > 0) tbData->alignShannonEntropy += tbData->b*logl(tbData->b);
    if (tbData->z > 0) tbData->alignShannonEntropy += tbData->z*logl(tbData->z);
    if (tbData->x > 0) tbData->alignShannonEntropy += tbData->x*logl(tbData->x);
    if (tbData->indel > 0) tbData->alignShannonEntropy += tbData->indel*logl(tbData->indel);
    
    return tbData->alignShannonEntropy;
}
MOATypeElmVal calcAlignmentSPScore (ProcessData * pData, TracebackData * tbData) {
    int k, l, startSearch, endSearch;
    MOATypeInd i, j;
    calcAlignmentScores (pData, tbData);
    tbData->alignSPScore = 0;
    tbData->meanSPColScore = 0;
    if (tbData->alignColScore != NULL)
        free (tbData->alignColScore);
    tbData->alignColScore = NULL;
    tbData->alignColScore = mmalloc (((MOATypeInd) tbData->comseqLen) * ((MOATypeInd) sizeof *(tbData->alignColScore)));

    for (i=0;i<tbData->comseqLen;i++) {
        tbData->alignColScore[i] = calcAlignmentColumnScore(pData, tbData, i);
        tbData->alignSPScore += tbData->alignColScore[i];
    }

    double accColScore;
    int alignResolution2 = alignResolution * 2;
    tbData->regD.maxColScore = 0;
    for (i=0;i<tbData->comseqLen-alignResolution2;i++) { /*searching all alignment length*/
        accColScore = tbData->alignColScore[i];
        for (j=i+1;j<i+alignResolution2;j++) { /* for the required region length*/
            accColScore += tbData->alignColScore[j];  /*should multiply by sequence weight here*/   
        }
        accColScore = (double) accColScore / alignResolution;
        if (i==0) {
            tbData->regD.maxColScore = accColScore;
            tbData->regD.maxStartColScore = i;
            tbData->regD.maxEndColScore = i+alignResolution2;
            tbData->regD.minColScore = accColScore;
            tbData->regD.minStartColScore = i;
            tbData->regD.minEndColScore = i+alignResolution2;
        }
        else {
            if(accColScore > tbData->regD.maxColScore) {
                tbData->regD.maxColScore = accColScore;
                tbData->regD.maxStartColScore = i;
                tbData->regD.maxEndColScore = i+alignResolution2;
            }
            if(accColScore < tbData->regD.minColScore) {
                tbData->regD.minColScore = accColScore;
                tbData->regD.minStartColScore = i;
                tbData->regD.minEndColScore = i+alignResolution2;
            }
        }
    }



    return tbData->alignSPScore;
}

void swapResidues (TracebackData * tbData, MOATypeInd i, MOATypeInd j, MOATypeInd k) {
    char temp;
    temp = tbData->completePath[i][j];
    tbData->completePath[i][j] = tbData->completePath[i][k];
    tbData->completePath[i][k] = temp;
}

void editAlignment (ProcessData * pData, TracebackData * tbData) {
	MOATypeInd i, j, k, maxSwapPosition;
        MOATypeElmVal maxColumnScore;

	for (j=0;j<tbData->comseqLen;j++) { // traverse columns
		for (i=0;i<pData->seqNum;i++) { // traverse rows                    
                    if (tbData->completePath[i][j] == GAPCHAR) {
                        maxColumnScore = tbData->alignColScore[j];
                        maxSwapPosition = j;
        		for (k=j+1;(k<tbData->comseqLen) && (tbData->completePath[i][k] != GAPCHAR);k++) { // traverse remaining columns till first residue
                            if (tbData->completePath[i][k] != GAPCHAR) {
                                swapResidues (tbData, i, j, k);
                                tbData->alignColScore[j] = calcAlignmentColumnScore(pData, tbData, j);
                                if (tbData->alignColScore[j] >= maxColumnScore) {
                                    maxColumnScore = tbData->alignColScore[j];
                                    maxSwapPosition = k;
                                }
                                swapResidues (tbData, i, k, j);
                                break;
                            }
                        }
        		for (k=j-1;(k>=0) && (tbData->completePath[i][k] != GAPCHAR);k--) { // traverse remaining columns till first residue
                            if (tbData->completePath[i][k] != GAPCHAR) {
                                swapResidues (tbData, i, j, k);
                                tbData->alignColScore[j] = calcAlignmentColumnScore(pData, tbData, j);
                                if (tbData->alignColScore[j] >= maxColumnScore) {
                                    maxColumnScore = tbData->alignColScore[j];
                                    maxSwapPosition = k;
                                }
                                swapResidues (tbData, i, k, j);
                                break;
                            }
                        }
                        swapResidues (tbData, i, maxSwapPosition, j);
                        tbData->alignColScore[j] = maxColumnScore;                        
                    }
		}
	}
}


int outputAlignment (ProcessData * pData, TracebackData * tbData, int flag) {
	FILE * outfile;
	char sfilename[SHORT_MESSAGE_SIZE];
	MOATypeInd i, j;

        if (flag == 1) {
            if (strlen(outputfilename) <= 0)
                    strcpy(outputfilename,  outPrefix);
            sprintf(sfilename, "../out/Align%s.aln", outputfilename);
            if (( outfile= fopen (sfilename, "w")) == NULL) {    
            //if (( outfile= open (sfilename, O_WRONLY|O_APPEND)) == -1) {
                    printf("Can not Open Alignment output file %s, exiting.\n", sfilename);
                    fflush(stdout);
                    return -1;
            }

            fprintf(outfile,"\n Aligned Sequences: ");
        }
        if (tbData->a > 0) { printf ("a %Lf\n", tbData->a); if (flag == 1) fprintf (outfile, "a %Lf\n", tbData->a);}
        if (tbData->c > 0) { printf ("c %Lf\n", tbData->c); if (flag == 1) fprintf (outfile, "c %Lf\n", tbData->c);}
        if (tbData->g > 0) { printf ("g %Lf\n", tbData->g); if (flag == 1) fprintf (outfile, "g %Lf\n", tbData->g);}
        if (tbData->t > 0) { printf ("t %Lf\n", tbData->t); if (flag == 1) fprintf (outfile, "t %Lf\n", tbData->t);}
        if (tbData->r > 0) { printf ("r %Lf\n", tbData->r); if (flag == 1) fprintf (outfile, "r %Lf\n", tbData->r);}
        if (tbData->n > 0) { printf ("n %Lf\n", tbData->n); if (flag == 1) fprintf (outfile, "n %Lf\n", tbData->n);}
        if (tbData->d > 0) { printf ("d %Lf\n", tbData->d); if (flag == 1) fprintf (outfile, "d %Lf\n", tbData->d);}
        if (tbData->q > 0) { printf ("q %Lf\n", tbData->q); if (flag == 1) fprintf (outfile, "q %Lf\n", tbData->q);}
        if (tbData->e > 0) { printf ("e %Lf\n", tbData->e); if (flag == 1) fprintf (outfile, "e %Lf\n", tbData->e);}
        if (tbData->i > 0) { printf ("i %Lf\n", tbData->i); if (flag == 1) fprintf (outfile, "i %Lf\n", tbData->i);}
        if (tbData->h > 0) { printf ("h %Lf\n", tbData->h); if (flag == 1) fprintf (outfile, "h %Lf\n", tbData->h);}
        if (tbData->l > 0) { printf ("l %Lf\n", tbData->l); if (flag == 1) fprintf (outfile, "l %Lf\n", tbData->l);}
        if (tbData->k > 0) { printf ("k %Lf\n", tbData->k); if (flag == 1) fprintf (outfile, "k %Lf\n", tbData->k);}
        if (tbData->m > 0) { printf ("m %Lf\n", tbData->m); if (flag == 1) fprintf (outfile, "m %Lf\n", tbData->m);}
        if (tbData->f > 0) { printf ("f %Lf\n", tbData->f); if (flag == 1) fprintf (outfile, "f %Lf\n", tbData->f);}
        if (tbData->p > 0) { printf ("p %Lf\n", tbData->p); if (flag == 1) fprintf (outfile, "p %Lf\n", tbData->p);}
        if (tbData->s > 0) { printf ("s %Lf\n", tbData->s); if (flag == 1) fprintf (outfile, "s %Lf\n", tbData->s);}
        if (tbData->w > 0) { printf ("w %Lf\n", tbData->w); if (flag == 1) fprintf (outfile, "w %Lf\n", tbData->w);}
        if (tbData->y > 0) { printf ("y %Lf\n", tbData->y); if (flag == 1) fprintf (outfile, "y %Lf\n", tbData->y);}
        if (tbData->v > 0) { printf ("v %Lf\n", tbData->v); if (flag == 1) fprintf (outfile, "v %Lf\n", tbData->v);}
        if (tbData->b > 0) { printf ("b %Lf\n", tbData->b); if (flag == 1) fprintf (outfile, "b %Lf\n", tbData->b);}
        if (tbData->z > 0) { printf ("z %Lf\n", tbData->z); if (flag == 1) fprintf (outfile, "z %Lf\n", tbData->z);}
        if (tbData->x > 0) { printf ("x %Lf\n", tbData->x); if (flag == 1) fprintf (outfile, "x %Lf\n", tbData->x);}
        if (tbData->indel > 0) { printf ("indel %Lf\n", tbData->indel); if (flag == 1) fprintf (outfile, "indel %Lf\n", tbData->indel);}
	printf("\n Aligned Sequences: ");
	for (i=0;i<pData->seqNum;i++) {
            if (flag == 1) 
		fprintf (outfile, "\n%s\t\t", pData->seqName[i]);
            printf ("\n%s\t\t", pData->seqName[i]);
            for (j=0;j<tbData->comseqLen;j++) {
                if (((tbData->completePath[i][j] >= 'A') && (tbData->completePath[i][j] <= 'Z')) ||((tbData->completePath[i][j] >= 'a') && (tbData->completePath[i][j] <= 'z')) || (tbData->completePath[i][j] == GAPCHAR)) {
                    if (flag == 1) 
                        fprintf(outfile, "%c", tbData->completePath[i][j]);
                    printf("%c", tbData->completePath[i][j]);
                }
            }
            if (flag == 1) 
                fprintf (outfile, " %lld", pData->seqLen[i]-1);
            printf (" %lld", pData->seqLen[i]-1);
	}
        if (flag == 1) 
            fprintf(outfile, "\n\t\t");	
        printf ("\n\t\t");
	for (j=0;j<tbData->comseqLen;j++) {
            if (flag == 1) 
		fprintf(outfile, "%lld ", tbData->alignColScore[j]);
            printf("%lld ", tbData->alignColScore[j]);
	}
        if (flag == 1) {
            fprintf (outfile, "\nThe Alignment SP Score = %lld\n", tbData->alignSPScore);
            fprintf (outfile, "\nThe Alignment Shannon Entropy = %Lf\n", tbData->alignShannonEntropy);
        }
	printf ("\nThe Alignment SP Score = %lld\n", tbData->alignSPScore);
	printf ("\nThe Alignment Shannon Entropy = %Lf\n Distance Measures are:\n", tbData->alignShannonEntropy);
        
        printf ("\t\t\t\t%15s", pData->seqName[0]);
        if (flag == 1) 
            fprintf (outfile, "\t\t\t\t%15s", pData->seqName[0]);
	for (i=1;i<pData->seqNum;i++) {
            printf ("\t\t%15s", pData->seqName[i]);
            if (flag == 1) 
                fprintf (outfile, "\t\t%15s", pData->seqName[i]);
        }
	for (i=0;i<pData->seqNum;i++) {
            printf ("\n%s\t\t%15Lf", pData->seqName[i], tbData->alignDistaneMeasures[i][0]);
            if (flag == 1) 
                fprintf (outfile, "\n%s\t\t%15Lf", pData->seqName[i], tbData->alignDistaneMeasures[i][0]);
            for (j=1;j<pData->seqNum;j++) {
                printf ("\t\t%15Lf", tbData->alignDistaneMeasures[i][j]);
                if (flag == 1) 
                    fprintf (outfile, "\t\t%15Lf", tbData->alignDistaneMeasures[i][j]);
            }            
        }
        printf ("\n\nGap Scores:\n\t\t\t\t%Lf", tbData->alignGapScore[0]);
        if (flag == 1) 
            fprintf (outfile, "\n\nGap Scores:\n\t\t\t\t%Lf", tbData->alignGapScore[0]);
	for (i=1;i<pData->seqNum;i++) {
            printf ("\t\t%Lf", tbData->alignGapScore[i]);
            if (flag == 1) 
                fprintf (outfile, "\t\t%Lf", tbData->alignGapScore[i]);
        }
        printf ("\nHighest Regions of Similarity\nScore\t\tFrom\t\tto\n");
        printf ("%fd\t\t%lld\t\t%lld\n", tbData->regD.maxColScore, tbData->regD.maxStartColScore, tbData->regD.maxEndColScore);
        printf ("\nLowest Regions of Similarity\nScore\t\tFrom\t\tto\n");
        printf ("%fd\t\t%lld\t\t%lld\n", tbData->regD.minColScore, tbData->regD.minStartColScore, tbData->regD.minEndColScore);
        if (flag == 1)  {
            fprintf (outfile, "\nHighest Regions of Similarity\nScore\t\tFrom\t\tto\n");
            fprintf (outfile, "%fd\t\t%lld\t\t%lld\n", tbData->regD.maxColScore, tbData->regD.maxStartColScore, tbData->regD.maxEndColScore);
            fprintf (outfile, "\nLowest Regions of Similarity\nScore\t\tFrom\t\tto\n");
            fprintf (outfile, "%fd\t\t%lld\t\t%lld\n", tbData->regD.minColScore, tbData->regD.minStartColScore, tbData->regD.minEndColScore);
        }
        if (flag == 1)  {
            if (fclose(outfile) != 0) {
            	printf ("ERROR Closing the Alignment output file %s!", sfilename);
            	return -1;   
            }
        }
	return 0;
}
int outputFastaAlignment (ProcessData * pData, TracebackData * tbData) {
	FILE * outfile;
	char sfilename[FILENAME_MAX_LENGTH];
	MOATypeInd i, j;

        if (strlen(outputfilename) <= 0)
                strcpy(outputfilename,  outPrefix);
        sprintf(sfilename, "../out/FASTA_%s.aln", outputfilename);
        if (( outfile= fopen (sfilename, "w")) == NULL) {    
        //if (( outfile= open (sfilename, O_WRONLY|O_APPEND)) == -1) {
                printf("Can not Open Alignment output file %s, exiting.\n", sfilename);
                fflush(stdout);
                return -1;
        }
	for (i=0;i<pData->seqNum;i++) {
            fprintf (outfile, ">%s\n", pData->seqName[i]);
            for (j=0;j<tbData->comseqLen;j++) 
                if (((tbData->completePath[i][j] >= 'A') && (tbData->completePath[i][j] <= 'Z')) ||((tbData->completePath[i][j] >= 'a') && (tbData->completePath[i][j] <= 'z')) || (tbData->completePath[i][j] == GAPCHAR)) 
                   fprintf(outfile, "%c", tbData->completePath[i][j]);
            fprintf(outfile, "\n");	
	}

        if (fclose(outfile) != 0) {
            printf ("ERROR Closing the Alignment output file %s!", sfilename);
            return -1;   
        }
	return 0;
}

size_t getLongestSeqName (MOATypeDimn seqNum, char * * seqNames) {
    MOATypeDimn k;
    size_t maxSize = strlen(seqNames[0]);
    for (k=1;k<seqNum;k++) {
        if (strlen(seqNames[k]) > maxSize)
            maxSize = strlen(seqNames[k]);
    }
    return maxSize;
}

int outputMSFAlignment_2 (char * sfilename, ProcessData * pData, TracebackData * tbData) {
	FILE * outfile;
	MOATypeInd i, j;
        char msffilename[FILENAME_MAX_LENGTH];
       
        if (strlen(sfilename) <= 0)
                strcpy(sfilename,  outPrefix);
 
        if (strstr(sfilename, "../out/") == NULL)
            sprintf(msffilename, "../out/%s.msf_2", sfilename);
        else
            sprintf(msffilename, "%s.msf_2", sfilename);

        if (( outfile= fopen (sfilename, "w")) == NULL) {    
        //if (( outfile= open (sfilename, O_WRONLY|O_APPEND)) == -1) 
                printf("Can not Open Alignment output file %s, exiting.\n", sfilename);
                fflush(stdout);
                return -1;
        }
        fprintf (outfile, "mmDst\n\n\n   MSF:   %lld  Type: P    Check:  xxxx   ..\n\n", tbData->comseqLen);
	for (i=0;i<pData->seqNum;i++) 
            fprintf (outfile, " Name: %s oo  Len:   %lld  Check:  xxxx  Weight:  xxx\n", pData->seqName[i], pData->seqLen[i]);
        fprintf (outfile, "\n//\n\n\n");
        int pos = 0;
        while (pos < tbData->comseqLen) {
        	for (i=0;i<pData->seqNum;i++) { 
                fprintf (outfile, "%20s", pData->seqName[i]);
                for (j=pos;(j<tbData->comseqLen) || (pos < 50);j++) {
                    if (((tbData->completePath[i][j] >= 'A') && (tbData->completePath[i][j] <= 'Z')) ||((tbData->completePath[i][j] >= 'a') && (tbData->completePath[i][j] <= 'z')) || (tbData->completePath[i][j] == GAPCHAR)) {
                       fprintf(outfile, "%c", tbData->completePath[i][j]);
                       pos ++;
                       if ((pos % 50) == 0)
                        fprintf(outfile, "\n");
                       if ((pos % 10) == 0)
                        fprintf(outfile, " ");
                    }
                }
            }
            fprintf(outfile, "\n\n");
        }
        if (fclose(outfile) != 0) {
            printf ("ERROR Closing the Alignment output file %s!", sfilename);
            return -1;   
        }
	return 0;
}

int outputMSFAlignment (char * sfilename, ProcessData * pData, TracebackData * tbData) {
	FILE * outfile;
	MOATypeInd i, j, cnt;
        char msffilename[FILENAME_MAX_LENGTH];
        
        if (strlen(sfilename) <= 0)
                strcpy(sfilename,  outPrefix);
        if (strstr(sfilename, "../out/") == NULL)
            sprintf(msffilename, "../out/%s.msf", sfilename);
        else
            sprintf(msffilename, "%s.msf", sfilename);
        if (( outfile= fopen (sfilename, "w")) == NULL) {    
        //if (( outfile= open (sfilename, O_WRONLY|O_APPEND)) == -1) {
                printf("Can not Open Alignment output file %s, exiting.\n", sfilename);
                fflush(stdout);
                return -1;
        }
        cnt = 0;
        while (cnt < tbData->comseqLen) {
            for (j=cnt;j<tbData->comseqLen;j++)  {
                for (i=0;i<pData->seqNum;i++) {
                    fprintf (outfile, "%12s", pData->seqName[i]);
                    if (((tbData->completePath[i][j] >= 'A') && (tbData->completePath[i][j] <= 'Z')) ||((tbData->completePath[i][j] >= 'a') && (tbData->completePath[i][j] <= 'z')) || (tbData->completePath[i][j] == GAPCHAR))  {
                        fprintf(outfile, "%c", tbData->completePath[i][j]);
                        cnt ++;
                        div_t div_val;
                        div_val = div(cnt, 50);
                        if (div_val.rem == 0) {
                            fprintf(outfile, "\n");
                            break;
                        }
                        div_val = div(cnt, 10);
                        if (div_val.rem == 0)
                           fprintf(outfile, " ");
                    }
                }
                fprintf(outfile, "\n");
            }
        }

        if (fclose(outfile) != 0) {
            printf ("ERROR Closing the Alignment output file %s!", sfilename);
            return -1;   
        }
	return 0;
}


#ifndef NDEBUG
void printSeq (int dbglevel, char * sequence , int sq_sz){
	int i;
	char msg[MID_MESSAGE_SIZE];
  
	sprintf(msg, "\nsequence[%d] > |", sq_sz);
	mprintf (dbglevel, msg, 1);
	for (i=0; i < sq_sz;i++) {
		sprintf(msg, "%c", sequence[i]);
		mprintf (dbglevel, msg, 1);
	}
	mprintf (dbglevel, "|", 1);
}
#endif
#ifndef NDEBUG
// Print Overlapping Outgoing Cells
void print_OCout(ProcessData * pData, WavesData * wData, int db_level) {
    long i, j;
    MOATypeDimn k;
    char msg[SHORT_MESSAGE_SIZE];

    sprintf (msg, "\nOverlapping Outgoing Cells: count[%ld] from proc [%d]\n", wData->wavesTotal, myProcid);
    mprintf (db_level, msg, 1);		
    for (i=0; i<wData->wavesTotal;i++) {
        for (j=0; j<pData->OCout[i].wavesOC;i++) {
            sprintf (msg, "cellIndex{%lld", pData->OCout[i].WOCO[j].cellIndex[0]);
            for (k=1;k<pData->seqNum;k++) 
                sprintf (msg, "%s, %lld", msg, pData->OCout[i].WOCO[j].cellIndex[k]);
            sprintf (msg, "%s}, cellScore[%lld], sent[%d]\n", msg, pData->OCout[i].WOCO[j].cellScore, pData->OCout[i].WOCO[j].sent);
            mprintf (db_level, msg, 1);
            if (pData->OCout[i].WOCO[j].depProc_ub > 0) {		
                sprintf (msg, "dependent proc: count[%d] {", pData->OCout[i].WOCO[j].depProc_ub);
                mprintf (db_level, msg, 1);		
                for (k=0; k<pData->OCout[i].WOCO[j].depProc_ub;k++) {
                    sprintf (msg, "%d ", pData->OCout[i].WOCO[j].depProc[k]);
                    mprintf (db_level, msg, 1);		
                }
                mprintf (db_level, "}\n", 1);		
            }
        }
    }
}
#endif
#ifndef NDEBUG
void to_proc_cells(ProcessData * pData, WavesData * wData, int proc, int db_level) {
    long i, j;
    MOATypeDimn k;
    char msg[SHORT_MESSAGE_SIZE];
    int cell_count;
	
    mprintf (db_level, "(  wn,   pi,   ci,     cs,   sent)\n", 1);
    cell_count = 0;
    for (i=0; i<wData->wavesTotal;i++) {
        for (j=0; j<pData->OCout[i].wavesOC;j++) {
            if (pData->OCout[i].WOCO[j].depProc_ub > 0) {		
                for (k=0; k<pData->OCout[i].WOCO[j].depProc_ub;k++) {
                    if (proc == pData->OCout[i].WOCO[j].depProc[k]) {
                        sprintf (msg, "(%4ld, {%4lld", msg, pData->OCout[i].WOCO[j].cellIndex[0]);
                        for (k=1;k<pData->seqNum;k++) 
                            sprintf (msg, "%s, %4lld", msg, pData->OCout[i].WOCO[j].cellIndex[k]);
                        sprintf (msg, "%s}, %6lld, %6d)\n", msg, pData->OCout[i].WOCO[j].cellScore, pData->OCout[i].WOCO[j].sent);
                        mprintf (db_level, msg, 1);
                        cell_count++;
                    }
                }
            }
        }
    }
    sprintf (msg, "Overlapping cells count: %d\n", cell_count);
    mprintf (db_level, msg, 1);
}
#endif

#ifndef NDEBUG
void print_outging_cells(ProcessData * pData, WavesData * wData, int db_level) {
    long i, j, l;
    char msg[SHORT_MESSAGE_SIZE];
    int * to_proc, proc_count, found, k;

    sprintf (msg, "\nOverlapping Outgoing Cells: from proc [%d]\n",  myProcid);
    mprintf (db_level, msg, 1);
    proc_count = 0;		
    for (i=0; i<wData->wavesTotal;i++) {
        for (j=0; j<pData->OCout[i].wavesOC;j++) {
            if (pData->OCout[i].WOCO[j].depProc_ub > 0) {		
                for (l=0; l<pData->OCout[i].WOCO[j].depProc_ub;l++) {
                    found = 0;
                    for (k=0; k<proc_count; k++) {
                        if (to_proc[k] == pData->OCout[i].WOCO[j].depProc[l]) {
                            found = 1;
                            break;
                        }
                    }
                    if (found == 0) {
                        proc_count++;
                        if (proc_count == 1) 
                            to_proc = mmalloc ((MOATypeInd) sizeof *to_proc);
                        else 
                            to_proc = realloc (to_proc, ((MOATypeInd) proc_count) * ((MOATypeInd) sizeof *to_proc));
                        to_proc[proc_count - 1] = pData->OCout[i].WOCO[j].depProc[l];
                    }
                }
            }
        }
    }
    if (proc_count > 0) {
        sprintf (msg, "send to proc. count[%d]\n", proc_count);
        mprintf (db_level, msg, 1);		
        for (i=0; i<proc_count;i++) {
            sprintf (msg, "to proc[%d]\n", to_proc[i]);
            mprintf (db_level, msg, 1);		
            to_proc_cells(pData, wData, to_proc[i], db_level);
        }
    }
}
#endif

#ifndef NDEBUG
// Print Overlapping Incoming Cells
void print_OCin(ProcessData * pData, WavesData * wData) {
    long i, j;
    MOATypeDimn k;
    char msg[SHORT_MESSAGE_SIZE];

    sprintf (msg, "\nOverlapping Incoming Cells:\n");
    mprintf (2, msg, 1);		
    for (i=0; i<wData->wavesTotal;i++) {
        for (j=0; j<pData->OCin[i].wavesOC;j++) {
            sprintf (msg, "(%4ld, {%4lld", msg, pData->OCin[i].WOCI[j].cellIndex[0]);
            for (k=1;k<pData->seqNum;k++) 
                sprintf (msg, "%s, %4lld", msg, pData->OCin[i].WOCI[j].cellIndex[k]);
            sprintf (msg, "%s}, %6lld)\n", msg, pData->OCin[i].WOCI[j].cellScore);
            mprintf (2, msg, 1);
        }
    }
}
#endif
#ifndef NDEBUG
void PrintSequencies (int dbglevel, MOATypeDimn seqNum, char * * sequences, MOATypeShape * seqLen) {
	MOATypeDimn i;
	char msg[MID_MESSAGE_SIZE];
  
	sprintf(msg, "\n[%d]>Sequences [%lld]: {", myProcid, seqNum);	
	mprintf(dbglevel, msg, 1);
	for (i=0;i<seqNum;i++) {
		printSeq (dbglevel, sequences[i], seqLen[i]);
	}
	mprintf(dbglevel, "\n}\n", 1);
}
#endif
#ifndef NDEBUG
void PrintASeq (MOATypeDimn seqNum, char * * sequences, MOATypeShape * seqLen, char * * * * algnseq, MOATypeShape * aseqLen, int alignmentsNo) {
	MOATypeShape i, j, k = 0;
	char msg[MID_MESSAGE_SIZE];
	int dbglevel = 1;
	
	sprintf(msg, "\n Sequences [%lld]: ", seqNum);
	mprintf(dbglevel, msg, 1);
	for (i=0;i<seqNum;i++) {
		printSeq (dbglevel, sequences[i], seqLen[i]);
	}
	sprintf(msg, "\n Found %d Optimal Alignments: \n", alignmentsNo);
	mprintf (dbglevel, msg, 1);
	for (k=0;k<alignmentsNo;k++) {
		sprintf(msg, "\n Aligned Sequences %lld: ", k+1);
		mprintf (dbglevel, msg, 1);
		for (i=0;i<seqNum;i++) {
			mprintf (dbglevel, "\n> ", 1);
			for (j=0;j<aseqLen[k];j++) {
				sprintf(msg, "%c ", (*algnseq)[k][i][j]);
				mprintf (dbglevel, msg, 1);
			}
		}
	}
}
#endif

struct tm * getTime () {
  time_t ltime;
  struct tm *m_now;
  
  time( &ltime );
  
  /*printf( "UNIX time and date:\t\t\t%s", _ctime64( &ltime ) );*/
  
  /* Convert to time structure and adjust for PM if necessary. */
  m_now = localtime( &ltime );
  
  if( m_now->tm_hour == 0 )  /* Adjust if midnight hour. */
    m_now->tm_hour = 12;
  return m_now;
}

void cpTime (struct tm  * currNow, struct tm  * * prevNow) {
	(*prevNow)->tm_hour = currNow->tm_hour;
	(*prevNow)->tm_isdst = currNow->tm_isdst;
	(*prevNow)->tm_mday = currNow->tm_mday;
	(*prevNow)->tm_min = currNow->tm_min;
	(*prevNow)->tm_mon = currNow->tm_mon;
	(*prevNow)->tm_sec = currNow->tm_sec;
	(*prevNow)->tm_wday = currNow->tm_wday;
	(*prevNow)->tm_yday = currNow->tm_yday;
	(*prevNow)->tm_year = currNow->tm_year;
}

int isTimeDiffEquals (struct tm  * currNow, struct tm  * prevNow, char unit, int value) {
	int result = 0;
  switch (unit) {
    default:
      printf ("Invalid menu choice - try again\n");
    break;
      
    case 'm':
			if (((currNow->tm_yday * 1440) + (currNow->tm_hour * 60) + currNow->tm_min) > ((prevNow->tm_yday * 1440) + (prevNow->tm_hour * 60) + prevNow->tm_min + value)) {
				result = 1;
			}
    break;
	}
	return result;

}

int init_output() {
	char sfilename[FILENAME_MAX_LENGTH];
   
  
	if (strlen(outputfilename) <= 0)
		strcpy(outputfilename,  outPrefix);
	sprintf(sfilename, "../out/%s%sp%dt%d", outPrefix, outputfilename, myProcid, 1);
	if (( outfile1 = fopen (sfilename, "w")) == NULL) {
		printf("Proc [%d] Can not Initialize output file %s, exiting.\n", myProcid, sfilename);
		fflush(stdout);
		return -1;
	}
    
	sprintf(sfilename, "../out/%s%sp%dt%d", outPrefix, outputfilename, myProcid, 2);
	if (( outfile2 = fopen (sfilename, "w")) == NULL) {
		printf("Proc [%d] Can not Initialize output file %s, exiting.\n", myProcid, sfilename);
		fflush(stdout);
		return -1;
	}
	sprintf(sfilename, "../out/%s%sp%dt%d", outPrefix, outputfilename, myProcid, 3);
	if (( outfile3 = fopen (sfilename, "w")) == NULL) {
		printf("Proc [%d] Can not Initialize output file %s, exiting.\n", myProcid, sfilename);
		fflush(stdout);
		return -1;
	}
  
	return 0;
}

int close_output () {
	//if (close(outfile) == -1)
	if (fclose(outfile1) != 0) {
		printf ("ERROR Closing the debugging file!");   
		fflush(stdout);
		return -1;
    	}
	if (fclose(outfile2) != 0) {
		printf ("ERROR Closing the debugging file!");   
		fflush(stdout);
		return -1;
    	}
	if (fclose(outfile3) != 0) {
		printf ("ERROR Closing the debugging file!");   
		fflush(stdout);
		return -1;
    	}
  	return 0;
}

int isEven(long i) 
{
  return !(i%2);
}

int initTBData (TracebackData * * tbData, MOATypeDimn seqNum, MOATypeShape * seqLen){
    MOATypeDimn i;
    
    (*tbData) = mmalloc (sizeof *(*tbData));
    if ((*tbData) == NULL) {
        printf ("Can not allocate memory for trace back data. Exiting\n");
        return -1;
    }
    (*tbData)->seqNum = seqNum;
    (*tbData)->maxCellIndex = NULL;
    (*tbData)->maxCellIndex = mcalloc ((MOATypeInd) seqNum, (MOATypeInd) sizeof *(*tbData)->maxCellIndex);
    if ((*tbData)->maxCellIndex == NULL)
        return -1;
    (*tbData)->partIndex = NULL;
    (*tbData)->partIndex = mcalloc ((MOATypeInd) seqNum, (MOATypeInd) sizeof *(*tbData)->partIndex);
    if ((*tbData)->partIndex == NULL)
        return -1;
    (*tbData)->pathParts = 0; /*partial alignments, will be increemented as found*/
    (*tbData)->aSeqLen = NULL; /*partial alignments lengths, will be created as found*/
    (*tbData)->algnseq = NULL; /*partial aligned sequences, will be created as found*/
    /*Allocated after end of distrbuted trace back*/
    (*tbData)->completePath = NULL;
    (*tbData)->alignColScore = NULL;
    (*tbData)->alignDistaneMeasures = NULL;
    (*tbData)->alignGapScore = NULL;
    return 0;
}

int freetbData (TracebackData * * tbData) {
    MOATypeDimn k;
    long i;
    if ((*tbData) != NULL) {
        if ((*tbData)->maxCellIndex != NULL) 
            free ((*tbData)->maxCellIndex);
        (*tbData)->maxCellIndex = NULL;
        if ((*tbData)->partIndex != NULL) 
            free ((*tbData)->partIndex);
        (*tbData)->partIndex = NULL;
        if ((*tbData)->algnseq != NULL) {
            for (i=0;i<(*tbData)->pathParts;i++) {
                if ((*tbData)->algnseq[i] != NULL) {
                    /*for (k=0;k<(*tbData)->seqNum;k++) {
                        if ((*tbData)->algnseq[i][k] != NULL) 
                                free ((*tbData)->algnseq[i][k]);
                        (*tbData)->algnseq[i][k] = NULL;
                    }*/
                    free ((*tbData)->algnseq[i]);
                    (*tbData)->algnseq[i] = NULL;
                }
            }
            free ((*tbData)->algnseq);
            (*tbData)->algnseq = NULL;
	}
	if ((*tbData)->aSeqLen != NULL) 
            free ((*tbData)->aSeqLen);
	(*tbData)->aSeqLen = NULL;
	if ((*tbData)->completePath != NULL) {
            for (k=0;k<(*tbData)->seqNum;k++) {
                if ((*tbData)->completePath[k] != NULL) 
                    free ((*tbData)->completePath[k]);
                (*tbData)->completePath[k] = NULL;
            }		
            free ((*tbData)->completePath);
            (*tbData)->completePath = NULL;
	}
	if ((*tbData)->alignColScore != NULL) 
            free ((*tbData)->alignColScore);
	(*tbData)->alignColScore = NULL;
	if ((*tbData)->alignGapScore != NULL) 
            free ((*tbData)->alignGapScore);
	(*tbData)->alignGapScore = NULL;
        if ((*tbData)->alignDistaneMeasures != NULL) {
            for (k=0;k<(*tbData)->seqNum;k++) {
                if ((*tbData)->alignDistaneMeasures[k] != NULL)
                    free ((*tbData)->alignDistaneMeasures[k]);
                (*tbData)->alignDistaneMeasures[k] = NULL;
            }
            free ((*tbData)->alignDistaneMeasures);                
        }
	(*tbData)->alignDistaneMeasures = NULL;
        free ((*tbData));
	(*tbData) = NULL;
    }
    return 0;
}
 
int mprintf (int dbglevel, const char *msg, int thread_num) {
	int ret = 0;

	if (dbglevel <= pdebug) {
		switch (thread_num) {
			case 1:
				ret = fprintf(outfile1, msg);
				break;
			case 2:
				ret = fprintf(outfile2, msg);
				break;
			case 3:
				ret = fprintf(outfile3, msg);
				break;
			default:
				ret = fprintf(outfile1, msg);
				break;
		}
		
		 if (ret < 0)
			printf ("[%d] Error writing to output file in thread %d\n", myProcid, thread_num);
		
	}
	return ret;
}

size_t getSizes () {
  printf ("size_t = %d \n", (int) sizeof(size_t));
  //printf ("size_t(0) = %d \n", size_t(0));
  //printf ("size_t(0)/sizeof(double) = %fd \n", (size_t(0)) / sizeof(double));
  printf ("int = %d \n", (int) sizeof(int));
  printf ("long = %d \n", (int) sizeof(long));
  printf ("long long = %d \n", (int) sizeof(long long));
  printf ("sizeof(long) %d \n", (int) sizeof(long));
  printf ("sizeof(unsigned long) %d \n", (int) sizeof(unsigned long));
  printf ("sizeof(unsigned long * ) %d \n", (int) sizeof(unsigned long * ));  printf ("Size of:\n Int %d unsigned Int %d\n Long %d unsigned Long %d\n Long Long %d unsigned Long Long %d\n", (int) sizeof(int), (int) sizeof(unsigned int), (int) sizeof(long), (int) sizeof(unsigned long), (int) sizeof(long long), (int) sizeof(unsigned long long));
  printf ("Data Values Ranges:\nShort Min %d, Short Max %d, Unsig Short Max %d\nInt Min %d, Int Max %d, Unsig Int Max %ud\nLong Min %ld, Long Max %ld, Unsig Long Max %uld\n", SHRT_MIN, SHRT_MAX, USHRT_MAX, INT_MIN, INT_MAX, UINT_MAX, LONG_MIN, LONG_MAX, ULONG_MAX);

  printf ("sizeof(MOA_elm) %d\n", (int) sizeof(MOA_elm));
  printf ("sizeof(MOA_rec) %d \n", (int) sizeof(MOA_rec));
  printf ("sizeof(ProcessData) %d \n", (int) sizeof(ProcessData));

  return sizeof(ProcessData);
}

void * mmalloc(MOATypeInd size) {
  void * alloc_mem = NULL;
  char msg[MID_MESSAGE_SIZE];
  alloc_mem = malloc (size);
  if (alloc_mem == NULL ) {
    sprintf(msg, "mmalloc: Can not allocate memory  with size %lld \n",  size);
    mprintf (0, msg, 1);
    return NULL;
  }
  return alloc_mem;
}

void * mcalloc(MOATypeInd num, MOATypeInd size) {
  void * alloc_mem = NULL;
  char msg[MID_MESSAGE_SIZE];
  
  alloc_mem = calloc (num, size);
  if (alloc_mem == NULL ) {
    sprintf(msg, "mcalloc: Can not allocate memory with num %lld and size %lld \n",  num, size);
    mprintf (0, msg, 1);
    return NULL;
  }
  return alloc_mem;
}


/*
  reinventing the wheel? probably.. :)
  this code is in the public domain
  application of shell-swap algorithm avoids inefficiency of
  multiple strlen invocations. ;)
*/

char* strrev(char *pB)
{
  char *pE = pB;
  char *s = pB;

  /* avoid trouble */
  if(!pB || !*pB) return s;

  /* start looking for the null, then backtrack to last char */
  while (*++pE);
  pE--;

  /* repeat until begin/end pointers meet */
  while (pB < pE)
  {
    /* swap begin+n and end-n chars */
    char hold = *pB;
    /* gratuitous use of postfix increment with dereference */ 
    *pB++ = *pE;
    *pE-- = hold;
  }
  return s;
}



long mpow (long x, long y) {
  long powerV, i;
  powerV = 1;
  for (i=1;i<=y;i++) {
    powerV=powerV*x;
  }

  return powerV;
}
MOATypeElmVal  a_max(MOATypeElmVal * values, MOATypeInd ubound) {
  MOATypeElmVal maxVal;
  MOATypeInd i;
  
  maxVal = values[0];
  for (i=1;i<ubound;i++) {
    if (maxVal < values[i])
      maxVal = values[i];
  }
  return maxVal;
}
MOATypeElmVal  a_min(MOATypeElmVal * values, MOATypeInd ubound) {
  MOATypeElmVal minVal;
  MOATypeInd i;
  
  minVal = values[0];
  for (i=1;i<ubound;i++) {
    if (minVal > values[i])
      minVal = values[i];
  }
  return minVal;
}
MOATypeElmVal  a_average(MOATypeElmVal * values, MOATypeInd ubound) {
  MOATypeElmVal avgVal = 0;
  MOATypeInd i;
  for (i=0;i<ubound;i++) {
      avgVal +=  values[i];
  }    
  if (ubound > 0)
    avgVal = (MOATypeElmVal) ceill ((long double ) avgVal / (long double ) ubound);
  return avgVal;
}

long Factorial (long n) {
  if (n==0)
    return 1;
  else
    return (n * Factorial(n-1));
}

void l2Comb(long * n, long k, long * * * Combin)
{
  int i, j, row, column;
  
  row = 0;
  column = 0;
  for (i=0;i<k;i++){
    for (j=i+1;j<k;j++){
      (*Combin)[row][column] = n[i];
      column++;
      (*Combin)[row][column] = n[j];
      row++;
      column = 0;
    }
  }
}
/* max val in 3 values 
int max (int val1, int val2, int val3) {
  int maxVal = 0;
  if (val1  > maxVal)
    maxVal  = val1;
  if (val2  > maxVal)
    maxVal  = val2;
  if (val3  > maxVal)
    maxVal  = val3;
  
  return maxVal;
}
*/
/*max value from a linear array */
int maxVal (int * values, int valLen) {
  int i, maxVal = 0;
  for (i=0; i < valLen; i++)
    if (values[i] > maxVal)
      maxVal  = values[i];
  
  return maxVal;
}
/* Algorithm by Donald Knuth. */
/* C implementation by Glenn C. Rhoads */
void Combinations(long n, long k, long * * * Combin)
{
  long i, j=1, *c, x, row, column;
  c = mmalloc(((MOATypeInd)  (k+3)) * ((MOATypeInd) sizeof *c));
  
  for (i=1; i <= k; i++) c[i] = i;
  c[k+1] = n+1;
  c[k+2] = 0;
  j = k;
  row = 0;
  column = 0;
  
 visit:
  for (i=k; i >= 1; i--) {
    (*Combin)[row][column] = c[i];
    column++;
  }
  column = 0;
  row++;
  
  if (j > 0) {x = j+1; goto incr;}
  
  if (c[1] + 1 < c[2])
    {
      c[1] += 1;
      goto visit;
    }
  
  j = 2;
  
 do_more:
  c[j-1] = j-1;
  x = c[j] + 1;
  if (x == c[j+1]) {j++; goto do_more;}
  
  if (j > k) {
    return;
  }
  
 incr:
  c[j] = x;
  j--;
  goto visit;
}

void reverse(char s[])
   {
   int c, i, j;
   
   for (i=0, j=strlen(s)-1; i < j;i++, j--)
      {

      c = s[i];
      s[i] = s [j];
      s[j] = c;
      }
   }
/*Function for actual copying of file, return 1 on success, 0 otherwise*/
int file_copy( char *oldname, char *newname, long bytes1, long bytes2)
     /*bytes1 denotes the start position in bytes, bytes2 denotes the end position*/
{
  FILE *fold, *fnew;
  int c;
  
  if ( ( fold = fopen( oldname, "rb" ) ) == NULL )
    return -1;
  
  if ( ( fnew = fopen( newname, "wb" ) ) == NULL  )
    {
      fclose ( fold );
      return -1;
    }
  
  /*Set file pointer to proper start location, as specified by user*/
  if ( ( fseek(fold, bytes1, SEEK_SET) != 0 ) )
    {
      fprintf(stderr, "\nError using fseek().");
      return -1;
    }
  
  while (1)
    {
      c=fgetc(fold);
      
      /*Continue copying until end of file or until the requested limit has been reached*/
      if ( !feof( fold ) && ftell(fold) <= bytes2)
	fputc( c, fnew );
      else
	break;
    }
  
  fclose ( fnew );
  fclose ( fold );
  
  return 1;
}


/*This function finds how many bytes need to copied, calculating it from the percentage inputted by the user*/
long get_bytes(float percent, char *source)
{
  long bytes;
  FILE *fold;
  
  if (percent<=100)
    {
      if ( ( fold = fopen( source, "rb" ) ) == NULL )
	{
	  //puts("Error opening source file");
	  return -1;
	}
      if ( (fseek(fold, 0, SEEK_END))!=0)
	{
	  fprintf(stderr, "\nError using fseek().");
	  return -1;
	}
      
      bytes=ftell(fold);
      bytes*=(percent/100);
    }
  else 
    {
      printf("Error in input\n");
      return -1;
    }
  
  fclose(fold);
  return bytes;
}


