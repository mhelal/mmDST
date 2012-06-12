/**********************************************************
* Author: Manal Helal																			*
* Last Modification Date: Fri 12 Jan 2007 03:39:51 AM EST *
* Project : MMSA - Multiple Sequence Alignment Based on 	*
* 					Mathematics of Arrays - PhD Experimentation		*
* File: utils.c, a header for global utility functions		*
***********************************************************/
#ifndef UTILSH_
#define UTILSH_

#include "moa.h"
#include "main.h"

#define Error( Str )        FatalError( Str )
#define FatalError( Str )   fprintf( stderr, "%s\n", Str ), exit( 1 )

int readAndConvertToFasta (char * sfilename);
int readFastaAlignment (char * sfilename, ProcessData * pData, TracebackData * tbData);
MOATypeDimn readSequencesFastaFile (char * sfilename, char * * * sequences, char * * * seqName, MOATypeShape * * seqLen);
MOATypeElmVal calcAlignmentSPScore (ProcessData * pData, TracebackData * tbData);
int outputAlignment (ProcessData * pData, TracebackData * tbData, int flag);
int outputFastaAlignment (ProcessData * pData, TracebackData * tbData);
int outputMSFAlignment (char * sfilename, ProcessData * pData, TracebackData * tbData);
int outputMSFAlignment_2 (char * sfilename, ProcessData * pData, TracebackData * tbData);
void editAlignment (ProcessData * pData, TracebackData * tbData);
int initTBData (TracebackData * * tbData, MOATypeDimn seqNum, MOATypeShape * seqLen);
int freetbData (TracebackData * * tbData);
size_t getSizes ();
void reverse(char s[]);
void printSeq (int dbglevel, char * sequence , int sq_sz);
void print_OCout(ProcessData * pData, WavesData * wData, int db_level);
void to_proc_cells(ProcessData * pData, WavesData * wData, int proc, int db_level);
void print_outging_cells(ProcessData * pData, WavesData * wData, int db_level);
void print_OCin(ProcessData * pData, WavesData * wData);
void PrintSequencies (int dbglevel, MOATypeDimn seqNum, char * * sequences, MOATypeShape * seqLen);
void PrintASeq (MOATypeDimn seqNum, char * * sequences, MOATypeShape * seqLen, char * * * * algnseq, MOATypeShape * aseqLen, int alignmentsNo);
int isEven(long i) ;
int processArguments (int argc, char ** argv, MOATypeDimn * seqNum, char * * * sequences, char * * * seqName, MOATypeShape * * seqLen, int * stype, long * partitionSize);
MOATypeShape readInputSequence (char * * fileName, char * * sequence);
void PrintSequencies (int dbglevel, MOATypeDimn seqNum, char * * sequences, MOATypeShape * seqLen);
void printSeq (int dbglevel, char * sequence , int sq_sz);
void PrintASeq (MOATypeDimn seqNum, char * * sequences, MOATypeShape * seqLen, char * * * * algnseq, MOATypeShape * aseqLen, int alignmentsNo);
void PrintOptimalPath (MOATypeDimn seqNum, char * * * algnseq, MOATypeShape aseqLen);
struct tm * getTime ();
int isTimeDiffEquals (struct tm  * currNow, struct tm  * prevNow, char unit, int value);
void cpTime (struct tm * currNow, struct tm  * * prevNow);
int init_output() ;
int close_output ();
int mprintf(int dbglevel, const char *msg, int thread_num);
char* strrev(char *pB);
long mpow (long x, long y);
MOATypeElmVal  a_max(MOATypeElmVal * values, MOATypeInd ubound);
MOATypeElmVal  a_min(MOATypeElmVal * values, MOATypeInd ubound);
MOATypeElmVal  a_average(MOATypeElmVal * values, MOATypeInd ubound);
long Factorial (long n);
void l2Comb(long * n, long k, long * * * Combin);
/*int max (int val1, int val2, int val3);*/
int maxVal (int * values, int valLen);
void Combinations(long n, long k, long * * * Combin);
void * mmalloc(MOATypeInd size);
void * mcalloc(MOATypeInd num, MOATypeInd size);
int file_copy(char *oldname, char *newname, long bytes1, long bytes2);
long get_bytes(float percent, char *source);

void dsort ( long n, MOATypeShape a[] );
void asort ( long n, MOATypeShape a[] );
void permute ( long n, MOATypeShape a[], int *more );

#endif
