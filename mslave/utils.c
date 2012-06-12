#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "globals.h"
#include "utils.h"


int processArguments (int argc, char ** argv, long * seqNum, char * * * sequences, long * * seqLen, int * stype) {
/*
-c = number of sequences
-s = Scoring Type
-l = local alignment
-g = global alignment - default
-m = maximum number of optimal paths to be retrieved
-e = Epsilons Value for search space reduction
-o = the output log filename (will be prefixed with the process name and suffixed with the processor number)
-d = the debug level (what to be written in the og files)
-p = the partition size
-n = number of slave nodes to be created
*/
  int argRead = 1;
	char msg[MID_MESSAGE_SIZE];
  long i;
  /* Initialize Arguments default Values*/
  (*seqNum) = 0;
  Epsilons = 0;
  (*stype) = 1;
  AlignmentType = Global;
  maxAlignmentsNumber = 20;
  strcpy(outputfilename, "mmsaOut");
  pdebug = 0;
  partitionSize = 1;
  /* read & process MOAMSA Arguments */
  while (argRead < argc) {
    /*1. Read SeqNum and Sequences File Names */
    if (strcmp(argv[argRead],"-c") == 0) {
      /*printf ("\nread -c");       */
      (*seqNum) = atoi(argv[++argRead]);
      if ((*seqNum) <= 1) {
				printf("Sequence Num must be more than 1\n");
				return -1;
      }
      if (argc  <= ((*seqNum) + 3))  {
				printf("More Arguments expected\n");
				return -1;
      }
      (*sequences) = (char * *) mcalloc ((*seqNum), sizeof(char *));
      (*seqLen) = (long *)  mcalloc ((*seqNum), sizeof(long));
      for (i=0;i<(*seqNum);i++) {	
				(*seqLen)[i] =readInput(argv[++argRead], &(*sequences)[i]);
				if ((*seqLen)[i] == -1) {
	  			printf("Sequence %ld File Could not open. Exiting.\n", i);
	  			return -1;
				}
				//(*seqLen)[i] ++;

	/* printSeq (sequences[i], seqLen[i]); */
      }
  		if ((*seqNum) <= 1) {
    		printf("Sequence Num must be more than 1\n");
    		return -1;
  		}
    }
    /*2. Read Scoring Type */
    else if (strcmp(argv[argRead], "-s") == 0) {      
      if (argc  <= argRead)  {
				printf("More Arguments expected\n");
				return -1;
      }
      if (argRead <argc) 
      	(*stype) = atoi(argv[++argRead]);
			else {
      	printf("Expected Numeric Argument after -s for Scoring Type required.\n");
      	return -1;
      } 
    }
    /* 3. Decide if Local or Global Alignment */
    else if (strcmp(argv[argRead], "-l") == 0) {
      AlignmentType = Local;
    }
      else if (strcmp(argv[argRead],"-g") == 0) {
	AlignmentType = Global;
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
      	printf("Expected Numeric Argument after -e for Epsilons Value.\n");
      	return -1;
      } 
   }

    /* 6. Decide the output filename */
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

    /* 7. Decide the debug Level */
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
    /* 8. Decide the partition size - for the distributed version */
    else if (strcmp(argv[argRead], "-p") == 0) {
      if (argc  < argRead)  {
				printf("More Arguments expected\n");
				return -1;
      }
      if (argRead <argc) 
				partitionSize = atoi(argv[++argRead]);
      else {
      	printf("Expected Numeric Argument after -p for partitions Size.\n");
      	return -1;
      } 
    }
    argRead ++;
  }
return 0;
}

long readInput (char * fileName, char * * sequence) {
  char elm, msg[SHORT_MESSAGE_SIZE];
  long i = 0;
  FILE * f = fopen (fileName, "r");
  if (f == NULL) {
    sprintf(msg, "File %s could not be opened!\n", fileName);
		mprintf (1, msg, threadnum);
    return -1;
  }
  (*sequence) = NULL; 
	(*sequence) = (char *) mmalloc (sizeof(char) * 2);      
  if ((*sequence) == NULL) {
		sprintf(msg, "Could not allocate memory for sequence %s!\n", fileName);
		mprintf (1, msg, threadnum);
		fclose(f);
		return -1;
  }
  i++;
	(*sequence)[0]= GAPCHAR;
  while (!feof(f)) {
    elm = fgetc(f);
    if (elm != EOF)  {
    if (((elm >= 'A') && (elm <= 'Z')) ||((elm >= 'a') && (elm <= 'z'))) {
      i++;
	  	(*sequence) = (char *) realloc ((*sequence), sizeof(char) * (i+1));
      if ((*sequence) == NULL) {
				sprintf(msg, "Could not allocate memory for sequence %s!\n", fileName);
				mprintf (1, msg, threadnum);
				fclose(f);
				return -1;
      }
      else {
				(*sequence)[i-1] = elm;
      }
      }
    }
  }
  
	(*sequence)[i] = '\0';
  fclose(f);
  return i;
}

void printSeq (char * sequence , int sq_sz){
  int i;
	char msg[MID_MESSAGE_SIZE];

  mprintf (2, "\n> ", threadnum);
  for (i=0; i< sq_sz-1;i++) {
    sprintf(msg, "%c ", sequence[i]);
		mprintf (1, msg, threadnum);
	}
}

void PrintASeq (long seqNum, char * * sequences, long * seqLen, char * * * * algnseq, long * aseqLen, int alignmentsNo) {
  long i, j, k = 0;
	char msg[MID_MESSAGE_SIZE];

  mprintf(1, "\n Sequences Read: ", threadnum);
  for (i=0;i<seqNum;i++) {
    printSeq (sequences[i], seqLen[i]);
  }
  sprintf(msg, "\n Found %d Optimal Alignments: \n", alignmentsNo);
	mprintf (1, msg, threadnum);
  for (k=0;k<alignmentsNo;k++) {
    sprintf(msg, "\n Aligned Sequences %ld: ", k+1);
		mprintf (1, msg, threadnum);
    for (i=0;i<seqNum;i++) {
      mprintf (1, "\n> ", threadnum);
      for (j=0;j<aseqLen[k];j++) {
				sprintf(msg, "%c ", (*algnseq)[k][i][j]);
				mprintf (1, msg, threadnum);
			}
    }
  }
}

void PrintOptimalPath (long seqNum, char * * * algnseq, long aseqLen) {
  long i, j, k = 0;
	char msg[MID_MESSAGE_SIZE];

    sprintf(msg, "\n Aligned Sequences %ld: ", k+1);
		mprintf (1, msg, threadnum);
    for (i=0;i<seqNum;i++) {
      mprintf (1, "\n> ", threadnum);
      for (j=0;j<aseqLen;j++) {
				sprintf(msg, "%c ", (*algnseq)[i][j]);
				mprintf (1, msg, threadnum);
			}
    }
}

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

int init_output() {
	char sfilename[SHORT_MESSAGE_SIZE];
  FILE * outfile;
  int i;

  
	if (strlen(outputfilename) <= 0)
		strcpy(outputfilename,  "mmsa");
for (i=1;i<=3;i++) {
	sprintf(sfilename, "%s_proc%d_thrd%d_log", outputfilename, myProcid, i);

  if (( outfile= fopen (sfilename, "w")) == NULL) {
  	printf("Proc [%d] Can not Initialize output file %s, exiting.\n", myProcid, sfilename);
    fflush(stdout);
    return -1;
  }
if (fclose(outfile) != 0) {
   printf ("ERROR Closing the debugging file!");   
   fflush(stdout);
   return -1;}
}
  return 0;
}

int mprintf (int dbglevel, const char *msg, int thread_num) {

  FILE * outfile;
	int ret;
	char sfilename[SHORT_MESSAGE_SIZE], umsg[LONG_MESSAGE_SIZE];
	if (dbglevel <= pdebug) {
		if (strlen(outputfilename) <= 0)
			strcpy(outputfilename,  "mmsa");
		sprintf(sfilename, "%s_proc%d_thrd%d_log", outputfilename, myProcid, thread_num);
  	if (( outfile= fopen (sfilename, "a")) == NULL) {  	//if (( outfile= open (sfilename, O_WRONLY|O_APPEND)) == -1) {  		printf("Can not Open output file, exiting.\n");
  	 	fflush(stdout);
   	 	return -1;
		}
  	//printf("Proc [%d] Writing %s to output file %s.\n", myProcid, msg, sfilename);
    //fflush(stdout);
		//if (dbglevel > 1)
		//	sprintf(umsg, "Proc [%d] %s", myProcid, msg);
		ret = fprintf(outfile, msg); 
		//ret = write (outfile, msg, (int) sizeof(msg));
		
	//if (close(outfile) == -1)
	if (fclose(outfile) != 0)
	   printf ("ERROR Closing the debugging file!");   
	}
	return 0;
}
void * mmalloc(size_t size) {
	void * alloc_mem = NULL;
	char msg[MID_MESSAGE_SIZE];

	alloc_mem = (long *) malloc (size);
    if (alloc_mem == NULL ) {
		sprintf(msg, "mmalloc: Can not allocate memory  with size %ld \n", size);
		mprintf (1, msg, threadnum);
	  return NULL;
    }
	return alloc_mem;
}

void * mcalloc(size_t num, size_t size) {
	void * alloc_mem = NULL;
	char msg[MID_MESSAGE_SIZE];

	alloc_mem = (long *) calloc (num, size);
    if (alloc_mem == NULL ) {
			sprintf(msg, "mcalloc: Can not allocate memory with num %ld and size %ld \n", num, size);
			mprintf (1, msg, threadnum);
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
long  a_max(long * values, long ubound) {
  long maxVal, i;
  
  maxVal = values[0];
  for (i=1;i<ubound;i++) {
    if (maxVal < values[i])
      maxVal = values[i];
  }
  return maxVal;
}
long  a_min(long * values, long ubound) {
  long minVal, i;
  
  minVal = values[0];
  for (i=1;i<ubound;i++) {
    if (minVal > values[i])
      minVal = values[i];
  }
  return minVal;
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
    c = mmalloc( (k+3) * sizeof(long));
    
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
