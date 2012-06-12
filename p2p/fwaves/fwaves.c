/******************************************************************************
* FILE: fwaves.c
* DESCRIPTION:  
* AUTHOR: Manal E. Helal
* LAST REVISED:
* Function:
*		main
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LINE_MAX 80

/* =============================================================================
	The main function:

============================================================================= */
int main(int argc, char **argv) {
    int ClusterSize, i, uFlag = 0;
    char outputfilename[20];
    char sfilename[20];
    FILE * sfile;
    char line[LINE_MAX];
    long computedPartitions, waveNo, newWaveNo;
    if (argc < 5) {
        printf ("Few Arguments Read.\nFirst argument need to be Number of processors.\nSecond argument need to be output filename.\nThird argument need to be the new Wave Number to be set.\nFourth arg is a flag: 1 to update. 0 to read only. Exiting.\n");
        return -1;
    }
    else {
        printf ("%s: arg1 %s arg2 %s\n", argv[0], argv[1], argv[2]);
        ClusterSize = atoi(argv[1]);
        strcpy (outputfilename, argv[2]);
        newWaveNo = atoi(argv[3]);
        uFlag = atoi(argv[4]);
        
        printf ("%s: read ClusterSize %d outputfilename %s\nProc\t\tParts\t\tWave No\n", argv[0], ClusterSize, outputfilename);
        for (i=0;i<ClusterSize;i++) {
            strcpy (sfilename, "");
            sfilename[0] = '\0';
            strcpy (outputfilename, argv[2]);
            sprintf(sfilename, "../out/chpS%s%d", outputfilename, i);

            if (( sfile = fopen (sfilename, "r")) == NULL) {
                printf("%d: Invalid checkpoint file: %s, exiting.\n", i, sfilename);
                return -1;
            }    
            if (fgets (line, LINE_MAX, sfile) == NULL) {
                printf("%d: Invalid checkpoint file: %s, exiting [reading Computed Partitions].\n", i, sfilename);
                return -1;
            }
            computedPartitions = atol(line);
            if (fgets (line, LINE_MAX, sfile) == NULL) {
                printf("%d: Invalid checkpoint file: %s, exiting [reading current wave number].\n", i, sfilename);
                return -1;
            }
            waveNo = atol(line);
            fclose(sfile);
            printf ("%d\t\t%ld\t\t%ld\n", i, computedPartitions, waveNo);
            fflush(stdout);
            if (uFlag == 1) {
                if (( sfile= fopen (sfilename, "w")) == NULL) {
                    printf("%d: Can not write  checkpoint file: %s, exiting.\n", i, sfilename);
                    return -1;
                }
                fprintf (sfile, "%ld\n", computedPartitions);
                fprintf (sfile, "%ld\n", newWaveNo);
                fclose(sfile);
                printf ("%d: Updated computedPartitions %ld waveNo %ld\n", i, computedPartitions, newWaveNo);
                fflush(stdout);
            }
        }
    }
    return 0;
} /* of main */
