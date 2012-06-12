#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "globals.h"

int mprintf(int dbglevel, const char *fmt, int count, ...) {
  FILE * outfile;
	va_list ap;
	char sfilename[20];

	printf ("Proc [%d] 1 fmt %s pd %d/%d count %d\n", myProcid, fmt, dbglevel, pdebug, count);
	fflush(stdout);
	if (dbglevel <= pdebug) {
		printf ("Proc [%d] 2\n", myProcid);
		fflush(stdout);
		if (strlen(outputfilename) <= 0)
			strcpy(outputfilename,  "mmsa");
		sprintf(sfilename, "%s_%d_log", outputfilename, myProcid);
		printf ("Proc [%d] 4\n", myProcid);
		fflush(stdout);

  	if (( outfile= fopen (sfilename, "a")) == NULL) {
  		printf("Can not Open output file, exiting.\n");
   	 fflush(stdout);
   	 return -1;
		}
		printf ("Proc [%d] 5\n", myProcid);
		fflush(stdout);
		va_start(ap, count);
		vfprintf(outfile, fmt, ap); 
		va_end(ap);
		//fprintf (outfile, "Proc[%d] from printf count %d\n", myProcid, count);
  	printf ("Proc [%d] 7\n", myProcid);
		fflush(stdout);

		fclose(outfile);
	}
		printf ("Proc [%d] 8	\n", myProcid);
		fflush(stdout);
	return 0;
}

int mprintfn(int dbglevel, const char *msg) {
  FILE * outfile;
	char sfilename[20];
	if (dbglevel <= pdebug) {
		if (strlen(outputfilename) <= 0)
			strcpy(outputfilename,  "mmsa");
		sprintf(sfilename, "%s_%d_log", outputfilename, myProcid);
		printf ("Proc [%d] 4\n", myProcid);
		fflush(stdout);

  	if (( outfile= fopen (sfilename, "a")) == NULL) {
  		printf("Can not Open output file, exiting.\n");
   	 fflush(stdout);
   	 return -1;
		}
		fprintf(outfile, msg); 
		fclose(outfile);
	}
	return 0;
}
