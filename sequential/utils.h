#ifndef UTILSH_
#define UTILSH_
struct tm * getTime ();
int mprintf(const char * stream, const char *template, int count, ...);
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

#endif
