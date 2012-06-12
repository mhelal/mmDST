#ifndef GLOBALSH_
#define GLOBALSH_

#include <time.h>

#define GAPCHAR '-' 
#define LONG_MESSAGE_SIZE 1000 
#define MID_MESSAGE_SIZE 500 
#define SHORT_MESSAGE_SIZE 200 

struct tm * prevNow, *currNow;

int threadnum;
int pdebug;
int maxAlignmentsNumber;
int Epsilons;
char outputfilename[SHORT_MESSAGE_SIZE];
int partitionSize;
enum ModeEnum {Sequential, Distributed};
enum ModeEnum Mode;
enum AlignmentTypeEnum {Global, Local};
enum AlignmentTypeEnum AlignmentType;

int myMasterid;/*My Master Process ID - each slave has only one at any given time*/
int myProcid;/*My Process ID*/
int ClusterSize;/*size of computing cluster*/

#endif
