#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main (int argc, char * argv[]) {
	int i, starti, slaveNodes;
	char arguments[1000];
	char command[2000];
	char timeCommand[1000];

	slaveNodes = 2;
	if (argc < 3) {
	  printf ("Please enter enough arguments, at least 2 sequences files to align.\n");
	  fflush(stdout);
	  return -1;
	}
	fflush(stdout);
	if (strcmp(argv[1],"-n") == 0) {
  	slaveNodes = atoi(argv[2]);	
		strcpy(arguments, "");
		starti = 3;
	}
	else {
		strcpy(arguments, argv[1]);
		starti = 2;
	}
	for (i=starti;i<argc;i++) {
		/*-n = number of slave nodes to be created*/
		if (strcmp(argv[i],"-n") == 0) {
      slaveNodes = atoi(argv[++i]);	
		}
		else {
			strcat(arguments, " ");
			strcat(arguments, argv[i]);
		}
	}
	sprintf (timeCommand, "/usr/bin/time -f '\\t%%E \\t%%U \\t%%S \\t%%M \\t%%t \\t%%K \\t%%D \\t%%p \\t%%X \\t%%Z \\t%%F \\t%%R \\t%%W \\t%%c \\t%%w \\t%%I \\t%%O \\t%%r \\t%%s \\t%%k \\t%%C \\t%%x\n' -a -o ");

	sprintf (command, "%s res_mmDst_all.out mpirun -np %d %s res_mmDst.out ./mmDst %s ", timeCommand, slaveNodes, timeCommand, arguments);
	//printf ("%s \n", command);
	i = system (command);
	printf ("Returning From Score Calculation with return Value %d \n", i);
	if (i >= 0) {
		sprintf (command, "%s res_mtb_all.out mpirun -np %d %s res_mtb.out ./mtb %s ", timeCommand, slaveNodes, timeCommand, arguments);
		//printf ("%s \n", command);
		i = system (command);
	}
	printf ("Returning From Trace Back with return Value %d\n", i);
}
