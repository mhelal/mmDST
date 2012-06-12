/**********************************************************
* Author: Manal Helal																			*
* Last Modification Date: Fri 12 Jan 2007 03:39:51 AM EST *
* Project : MMSA - Multiple Sequence Alignment Based on 	*
* 					Mathematics of Arrays - PhD Experimentation		*
* File: globals.h, header file for all global variables		*
***********************************************************/
#ifndef GLOBALSH_
#define GLOBALSH_

#include "moa.h"


#define LONG_MESSAGE_SIZE 1000 
#define MID_MESSAGE_SIZE 500 
#define SHORT_MESSAGE_SIZE 200 
#define FILENAME_MAX_LENGTH 100
#define LINE_MAX 80

FILE * outfile1;
FILE * outfile2;
FILE * outfile3;
struct tm * prevNow, * currNow;

int pdebug, threadnum;
int RestoreFlag, TBFlag; /*Scoring Restoration and Trace Back Flags*/
int maxAlignmentsNumber;
long Epsilons;
int EpsilonsType;
int alignResolution; /*For region of similarity plots and identification of highest and lowest*/


char outputfilename[SHORT_MESSAGE_SIZE];
char outPrefix[5];
enum initializationEnum {IndexProductInit, IndexSumInit, GapPenaltyInit};
enum initializationEnum initializationMethod;
enum ModeEnum {Sequential, Distributed};
enum ModeEnum Mode;
enum AlignmentTypeEnum {Global, Local};
enum AlignmentTypeEnum AlignmentType;

enum SchedulingMethod {RR, BT, DB}; /* Bag of Tasks, Round Robin, Dependency Based */
enum ScAlgorithm {DP, SP, RB}; /* Tensor Scoring Algorithm, Dynamic Programming, Sum-of-Pairs, RubberBand */
enum ScAlgorithm Algorithm;
enum SchedulingMethod SchedMethod;

int myProcid;/*My Process ID*/
int ClusterSize;/*size of computing cluster*/
#define MOAMSA_SEND_RECEIVE_TAG 1
/* MOAMSA_SEND_RECEIVE_TAG is message type (tag) for MPI_Send and MPI_Recv*/

#endif
