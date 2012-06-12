#include <stdio.h>
#include<string.h>
#include <stdlib.h>
#include <ctype.h>      /* include the character library */
#include "moa.h"
#include "moamsa.h"
#include "utils.h"
#include "lq.h"

const char GABCHAR = '-';



int main (int argc, char ** argv) {
  MOA_rec * MOA1 = NULL; 
  MOA_rec * MOA2 = NULL;

  long * shape = NULL;
  long * ind = NULL;
  long * ind2 = NULL;
  long flatIndex;
  long i, rslt_ub, strides; 
  int finished = 1;
  long * rslt = NULL;
  long * rsltInd = NULL;

  /*open once to erase previous contents then close */
  outputfilename = "moamsa.out";
  if (( outfile= fopen (outputfilename, "w")) == NULL) {
    printf("Can not Open output file, exiting.\n");
    return;
  }
  fclose(outfile);
  createMOAStruct(&MOA1);
  ind = (long *) mmalloc ( MOA1->dimn * sizeof(long));
  ind2 = (long *) mmalloc ( MOA1->dimn * sizeof(long));
  finished = 1;
  while (finished != 0) {
    switch (finished) {
    case 0:
    default:
      printf("\n Good Bye\n");
      break;
    case 1:     
      deleteMOA (MOA1); 
	  MOA1 = NULL;
      free(ind);
	  ind = NULL;
      free(ind2);
	  ind2 = NULL;
      createMOAStruct(&MOA1);
      printf ("\nEnter dimensionality:");
      scanf("%ld", & MOA1->dimn);
      ind = (long *) mmalloc ( MOA1->dimn * sizeof(long));
      ind2 = (long *) mmalloc ( MOA1->dimn * sizeof(long));
      shape = (long *) mmalloc ( MOA1->dimn * sizeof(long));
      for (i = 0; i < MOA1->dimn; i++) {
	    printf ("\nEnter shape at dimension %ld:", i);
	    scanf("%ld", &shape[i]);
      }

      createMOA(shape /* shape*/, MOA1->dimn /* dimension*/, MOA1 /* MOA structure*/, -1, 0);
      /*printMOA(MOA1); // we didn't create the tensor to print it*/
      printf("created the MOA\n");
      break;
    case 2:     
      printMOA(MOA1);
      break;
    case 3:
      for (i = 0; i < MOA1->dimn; i++) {
	printf ("\nEnter index at dimn %ld:", i);
	scanf("%ld", &ind[i]);
      }
      
      flatIndex = Gamma(ind, MOA1->dimn, shape,  MOA1->dimn, 1);
      printf("\n the flatIndex = %ld \n", flatIndex);
      break;
    case 4:
      printf ("\nEnter flat index :" );
      scanf("%ld", &flatIndex);
      Gamma_Inverse(flatIndex, shape,  MOA1->dimn, ind);
      for (i = 0; i < MOA1->dimn; i++) {
	printf("%ld", ind[i]);
      }
      break;
    case 5:
      MOA2 = (MOA_rec *) mmalloc(sizeof(MOA_rec));
      for (i = 0; i <  MOA1->dimn; i++) {
	    printf ("\nEnter take index at dimn %ld:", i);
	    scanf("%ld", &ind[i]);
      }
      Take (ind, MOA1->dimn, MOA1, MOA2);
      printMOA(MOA2);
      
      deleteMOA (MOA2);  
	  MOA2 = NULL;
      break;
    case 6:
      MOA2 = (MOA_rec *) mmalloc(sizeof(MOA_rec));
      for (i = 0; i <  MOA1->dimn; i++) {
	    printf ("\nEnter Drop index at dimn %ld:", i);
	    scanf("%ld", &ind[i]);
      }
      Drop (ind, MOA1->dimn, MOA1, MOA2);
      printMOA(MOA2);
      deleteMOA (MOA2);  
	  MOA2 = NULL;
     break;
      /* case 5:
      printf ("\nEnter flat index :" );
      scanf("%d", &flatIndex);
      Gamma_Inverse(flatIndex, shape,  MOA1->dimn, ind);
      LNghbCount  = GetLowerNeighbors ( MOA1->dimn, shape, ind, &lNeighbors);
      printf("\n Has %d Neighbors :", LNghbCount);
      for (i = 0; i < LNghbCount; i++) {
	printf(" %d ", lNeighbors[i]);
      }
      break;*/
    case 7:
      MOA2 = (MOA_rec *) mmalloc(sizeof(MOA_rec));
      for (i = 0; i <  MOA1->dimn; i++) {
	    printf ("\nEnter Cell index at dimn %ld:", i);
	    scanf("%ld", &ind[i]);
      }
      MOAGetLowerNeighbors (ind, MOA1, MOA2);
      printMOA(MOA2);
      
      printf("With Indices ");
      for (i=0;i<MOA2->elements_ub; i++)
	    printf("%ld ", MOA2->indexes[i]);
      printf("\n");
      deleteMOA (MOA2); 
	  MOA2 = NULL;
     break;
    case 8:
      MOA2 = (MOA_rec *) mmalloc(sizeof(MOA_rec));
      for (i = 0; i <  MOA1->dimn; i++) {
	    printf ("\nEnter Cell index at dimn %ld:", i);
	    scanf("%ld", &ind[i]);
      }
      printf ("\nEnter The Stride Size:");
      scanf("%ld", &flatIndex);
      MOAGetHigherNeighbors (flatIndex, ind, MOA1, MOA2);
      printMOA(MOA2);
      
      printf("With Indices ");
      for (i=0;i<MOA2->elements_ub; i++)
	     printf("%ld ", MOA2->indexes[i]);
      printf("\n");
      deleteMOA (MOA2);
	  MOA2 = NULL;
     break;
    case 9:
      printf ("\nEnter The Size of the Wave:");
      scanf("%ld", &flatIndex);
      for (i = 0; i <  MOA1->dimn; i++) {
	    ind[i] = 1;
      }
      getDiagonals (1, flatIndex, ind, MOA1);
     break;
    case 99:
      for (i = 0; i <  MOA1->dimn; i++) {
	    printf ("\nEnter Navigation Director at dimn %ld:", i);
	    scanf("%ld", &ind[i]);
      }
      for (i = 0; i <  MOA1->dimn; i++) {
	    printf ("\nEnter Starting Cell index at dimn %ld:", i);
	    scanf("%ld", &ind2[i]);
      }
      printf ("\nEnter Strides: ");
      scanf("%ld", &strides);
      
      Navigate (ind, strides, ind2, MOA1, &rslt, &rsltInd, &rslt_ub);
      printf("rslt_ub = %ld \n", rslt_ub);
      
      for (i=0;i<rslt_ub; i++)
      	printf("%ld at %ld \n", rslt[i],rsltInd [i]);
      printf("\n");
      free(rslt);
	  rslt = NULL;
      free(rsltInd);
	  rsltInd = NULL;
      break;

    }
    printf ("\n1: Re-Define MOA\n");
    printf ("2: Print MOA\n");
    printf ("3: Gamma\n");
    printf ("4: Gamma_inverse\n");
    printf ("5: Take\n");
    printf ("6: Drop\n");
    printf ("7: MOAGetLowerNeighbors\n");
    printf ("8: MOAGetHigherNeighbors\n");
    printf ("9: GetDiagonals\n");
    printf ("0: Quit\n");
    printf ("\nEnter op code:" );
    scanf("%d", &finished);
    if (finished == 0) 
      printf("\n Good Bye\n");
  }
  free(ind);
  free(ind2);
  deleteMOA (MOA1);  
}
