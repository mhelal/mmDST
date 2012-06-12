#ifndef UTILSH_
#define UTILSH_


int processArguments (int argc, char ** argv, long * seqNum, char * * * sequences, long * * seqLen, int * stype);
long readInput (char * fileName, char * * sequence);
void printSeq (char * sequence , int sq_sz);
void PrintASeq (long seqNum, char * * sequences, long * seqLen, char * * * * algnseq, long * aseqLen, int alignmentsNo);
void PrintOptimalPath (long seqNum, char * * * algnseq, long aseqLen);
struct tm * getTime ();
int init_output() ;
int mprintf(int dbglevel, const char *msg, int thread_num);
char* strrev(char *pB);
long mpow (long x, long y);
long  a_max(long * values, long ubound);
long  a_min(long * values, long ubound);
long Factorial (long n);
void l2Comb(long * n, long k, long * * * Combin);
/*int max (int val1, int val2, int val3);*/
int maxVal (int * values, int valLen);
void Combinations(long n, long k, long * * * Combin);
void * mmalloc(size_t size);
void * mcalloc(size_t num, size_t size);
int file_copy(char *oldname, char *newname, long bytes1, long bytes2);
long get_bytes(float percent, char *source);
#endif
