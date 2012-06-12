/****************************************************************
* Author: Manal Helal                                           *
* Last Modification Date: Fri 12 Jan 2007 03:39:51 AM EST       *
* Project : MMSA - Multiple Sequence Alignment Based on 	       *
* 					Mathematics of Arrays - PhD Experimentation		 *
* File: moa.c
* Description: contain the basic MOA functions used in 		    *
*       this project                                            *
* Function:
*		Tau
*		Gamma
*     deleteMOA
*     Gamma_Inverse
*     createMOAStruct
*     printMOA
*     isLowerBorderCell
*     isHigherBorderCell
*     isHigherBorderCellandNotLower
*     MOAGetGlobalLowerNeighbors - MOAGetLowerNeighbors
*     NextIndex
*     copyIndicesElm
*     copyIndices
*     Take
*     TakeInd
*     Drop
*     DropInd
*     VecIsEqual
*     MOAGetHigherNeighbors
****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "globals.h"
#include "moa.h"
#include "utils.h"


/*****************************************************
	Function: createMOAStruct
		Allocate memory for MOA record and initialize its elements.
*****************************************************/
void createMOAStruct(MOA_rec * * MOA_val) {

  (*MOA_val) = mmalloc((MOATypeInd) sizeof *(*MOA_val));
  (*MOA_val)->dimn = 1;
  (*MOA_val)->shape = NULL;
  (*MOA_val)->indexes = NULL;
  (*MOA_val)->elements = NULL;
}

int createMOA(MOATypeShape * shape, MOATypeDimn dimn, MOA_rec * MOA_val, int callFlag, int cellValue)
{
    MOATypeInd i;
	
    MOA_val->elements_ub = Tau(shape, dimn);
    if (MOA_val->elements_ub == 0) {
        mprintf(1, "The smallest type of MOA array is a scalar that requires at least one element! NO MOA type supports zero element list!\n", threadnum);
        return -1;
    }
    MOA_val->dimn = dimn;
    if (dimn == 0) {
        MOA_val->shape = mcalloc ((MOATypeInd) 1, (MOATypeInd) sizeof *MOA_val->shape);
        MOA_val->shape[0] = 1;
    }
    else {
        MOA_val->shape = mcalloc((MOATypeInd) MOA_val->dimn, (MOATypeInd) sizeof *MOA_val->shape);
            for (i = 0; i< dimn; i++) 
                MOA_val->shape[i] = shape[i];            
   }
    if (MOA_val->shape == NULL) {
        printf("Error Allocating memory for array shape. Exiting Function.");
        return -1;
    }
    /* if callFlag is -1 means, don't allocate memory for the tensor, this is used with sparse arrays implementation. */

    if (callFlag != -1) {
        /*MOA_val->elements_ub =  Tau(shape, dimn); */
        MOA_val->elements = mmalloc (MOA_val->elements_ub * (MOATypeInd) sizeof *MOA_val->elements);
        if (MOA_val->elements == NULL) {
            printf("Error Allocating memory for array elements. Exiting Function.");
            return -1;
        }
        MOA_val->indexes = mmalloc (MOA_val->elements_ub * (MOATypeInd) sizeof *MOA_val->indexes);    
        if (MOA_val->indexes == NULL) {
            printf("Error Allocating memory for array indexes. Exiting Function.");
            return -1;
        } 
        for (i = 0; i< MOA_val->elements_ub; i++) {	
            MOA_val->elements[i].prev_ub = 0;	
            MOA_val->elements[i].prev_ub = 0;
            MOA_val->elements[i].prev = NULL;
            MOA_val->indexes[i] = NULL;
            MOA_val->indexes[i] = mcalloc (MOA_val->dimn, (MOATypeInd) sizeof *MOA_val->indexes[i]);    
            if (MOA_val->indexes[i] == NULL) {
                printf("Error Allocating memory for array indexes 2. Exiting Function.");
                return -1;
            }
            /* if callFlag is 0, means initialize cell values with 0, which is done with mcalloc */	
            if (callFlag == 1) {
                /* if callFlag is 1 meaning, initialize the array values with the cellValue argument value */
                MOA_val->elements[i].val = cellValue;
            }
        }
    }
    return 0;
}
/* ***************************************************************
* Function Name:  Tau
* Description: Returns Number of elements in an MOA structure 
*       based on dimension (array_ub) and shape (array_in) - regular MOA
*****************************************************************/

MOATypeInd Tau (MOATypeShape * array_in, MOATypeDimn array_ub) {
	MOATypeDimn k;
	MOATypeInd rslt;
        if (array_in[0] == 0)
            rslt = 1;
        else
            rslt = array_in[0];
	for ( k = 1; k< array_ub; k++) {
            if (array_in[k] != 0)
		rslt = rslt * array_in[k];	
	}
	return (rslt);
} /* end Tau */

/* ***************************************************************
* Function Name:  Sum
* Description: Returns total summation of elements in a linear array 
*****************************************************************/

MOATypeInd Sum (MOATypeShape * array_in, MOATypeDimn array_ub) {
    MOATypeDimn k;
    MOATypeInd rslt;
    rslt = array_in[0];
    for ( k = 1; k< array_ub; k++) {
            rslt = rslt + array_in[k];	
    }
    return (rslt);
} /* end Tau */


/*********************************************************************
	Function: deleteMOA
		free allocated memory for MOA record
*********************************************************************/
void deleteMOA(MOA_rec * MOA) {
    MOATypeInd i, j;
#ifndef NDEBUG
    int dbglevel = 30;
#endif

    if (MOA != NULL) {	
#ifndef NDEBUG
        mprintf (dbglevel, "deleteMOA 1\n", 1);
#endif
        if (MOA->elements != NULL) {
#ifndef NDEBUG
            mprintf (dbglevel, "deleteMOA 2\n", 1);
#endif
            for (i = 0; i< MOA->elements_ub; i++) {
#ifndef NDEBUG
                mprintf (dbglevel, "deleteMOA 3\n", 1);
#endif
                if ((MOA->elements[i].prev_ub > 0) && (MOA->elements[i].prev != NULL)) {
                    for (j = 0; j<MOA->elements[i].prev_ub; j++) {
                        if (MOA->elements[i].prev[j] != NULL) {
                            free(MOA->elements[i].prev[j]);
                            MOA->elements[i].prev[j] = NULL;            
                        }
                    }
                    free(MOA->elements[i].prev);
                    MOA->elements[i].prev = NULL;
                }
	        if (MOA->indexes != NULL) {
                    if (MOA->indexes[i] != NULL) {
                        free (MOA->indexes[i]);
                        MOA->indexes[i] = NULL;
                    }
		}
            }
            free(MOA->elements);
            MOA->elements = NULL;
#ifndef NDEBUG
            mprintf (dbglevel, "deleteMOA 4\n", 1);
#endif
        }
        if (MOA->indexes != NULL)
            free(MOA->indexes);
        MOA->indexes = NULL;
#ifndef NDEBUG
        mprintf (dbglevel, "deleteMOA 5\n", 1);
#endif
        if (MOA->shape != NULL)
            free(MOA->shape);
        MOA->shape = NULL;
#ifndef NDEBUG
        mprintf (dbglevel, "deleteMOA 6\n", 1);
#endif    
        free(MOA);
        MOA = NULL;
#ifndef NDEBUG
        mprintf (dbglevel, "deleteMOA 7\n", 1);
#endif
    }
}

/* ***************************************************************
* Function: Gamma
* Description: Get flay index of MOA array
*****************************************************************/
MOATypeInd Gamma (MOATypeShape * m_index, MOATypeDimn ind_ub, MOATypeShape * arr_shape, MOATypeDimn shape_ub, int Front) {
	MOATypeInd i;
	MOATypeInd  rslt;
	rslt = m_index[shape_ub - 1];
	for (i = shape_ub - 2; i>=0;i--) {
		rslt = (rslt * arr_shape[i]) + m_index[i];
   	}
	return rslt;

}

/***************************************************************
 **********************  Gamma Inverse************************
 *************************************************************/

int Gamma_Inverse (MOATypeInd ind, MOATypeShape * arr_shape, MOATypeDimn shape_ub, MOATypeShape * * rslt, int thrd) {
	MOATypeDimn i;
	MOATypeInd ind_temp;
	ldiv_t  div_result;
	char msg[MID_MESSAGE_SIZE];

	if (ind >= Tau(arr_shape, shape_ub)) {
		sprintf (msg, "Gamma_Inverse: The index is beyond the upper bound of the array of the defined shape! (%lld/%lld)\n", ind, Tau(arr_shape, shape_ub));
		mprintf(0, msg, thrd);
		return -1;
	}
	
	ind_temp = ind;
	for (i = 0; i< shape_ub; i ++) {
                (*rslt)[i] = 0;
                if (arr_shape[i] > 0) {
                    div_result = ldiv((long) ind_temp, (long) arr_shape[i]);
                    (*rslt)[i] = (MOATypeShape) div_result.rem; /* ind_temp Mod arr_shape->GetAt(i) */
                    ind_temp = (MOATypeInd) div_result.quot; /* ind_temp \ arr_shape->GetAt(i); */
                }
	}
	return 0;
} /* end Gamma_Inverse */



#ifndef NDEBUG
void printMOA_Matrix(int dbglevel, MOATypeDimn dimn, int ident, void * elements, int elm_type, MOATypeInd total_elements, MOATypeDimn first_dimn, char * * sequences) {
	MOATypeInd i;
	MOATypeDimn k, j, dim2cnt, elm_value;
	MOATypeShape * * MOA_indexes;
	char msg[MID_MESSAGE_SIZE];
	MOA_elm * MOA_elements;
	
	if (elm_type == 0) MOA_elements = elements;
	else MOA_indexes = elements;
	for (j = 0; j < ident; j++) 
		mprintf(dbglevel, "   ", 1);
	mprintf(dbglevel, "  ", 1);
	for (j = 0; j < first_dimn; j++) {
		sprintf(msg, "     %c", sequences[0][j]);
		mprintf (dbglevel, msg, 1);
	}
	mprintf (dbglevel, "\n", 1);
	for (j = 0; j < ident; j++) 
		mprintf(dbglevel, "   ", 1);
	mprintf(dbglevel, "  ", 1);
	for (j = 0; j < first_dimn; j++) {
		mprintf (dbglevel, "______", 1);
	}
	mprintf (dbglevel, "\n", 1);
	dim2cnt = 0;
	for (j = 0; j < ident; j++) 
		mprintf(dbglevel, "   ", 1);
	sprintf(msg, "%c|", sequences[1][dim2cnt]);
	mprintf (dbglevel, msg, 1);
	dim2cnt ++;	
	for (i = 0; i < total_elements; i++) {
		if (i>0 && (i%first_dimn) == 0) {
			mprintf (dbglevel, "\n", 1);
			for (j = 0; j < ident; j++) 
				mprintf(dbglevel, "   ", 1);
			sprintf(msg, "%c|", sequences[1][dim2cnt]);
			mprintf (dbglevel, msg, 1);
			dim2cnt ++;	
		}
		if (elm_type == 0) 
        		sprintf(msg, " %5lld", elm_value);
		else {
                    sprintf (msg, "{%5lld", MOA_indexes[i][0]);
                    for (k = 1; k < dimn; k++) 
                        sprintf(msg, "%s, %5lld", msg, MOA_indexes[i][k]);
                    sprintf (msg, "%s} ", msg);
                }
		mprintf (dbglevel, msg, 1);
	}
	mprintf (dbglevel, "\n\n", 1);
	for (j = 0; j < ident; j++) 
		mprintf(dbglevel, "___", 1);
	mprintf(dbglevel, "__", 1);
	for (j = 0; j < first_dimn; j++) {
		mprintf (dbglevel, "______", 1);
	}
	mprintf (dbglevel, "\n", 1);
}
#endif

#ifndef NDEBUG
int printMOA_Sequences(int dbglevel, int ident, MOATypeDimn dimn, MOATypeShape *ind, char * * sequences, MOA_rec * MOA, int elm_type) {
	MOATypeDimn i;
	MOATypeShape d, k;
	char msg[MID_MESSAGE_SIZE];
	
	if (dimn > 1) {
		for (i = 0; i < ident; i++) 
			mprintf(dbglevel, "   ", 1);
		sprintf(msg, "[%c]", sequences[dimn][ind[dimn]]);
		mprintf (dbglevel, msg, 1);
		mprintf (dbglevel, "\n", 1);
		ident++;
		ident = printMOA_Sequences(dbglevel, ident, dimn-1, ind, sequences, MOA, elm_type);
	}
	if (dimn == 1) {
		d = 1;
		k = 0;
		for (i = 0; i < MOA->dimn; i++) {
			k += d * ind[i];
			d *= MOA->shape[i];
		}
		if (elm_type == 0)
			printMOA_Matrix(dbglevel, MOA->dimn, ident, &(MOA->elements[k]), elm_type, MOA->shape[0] * MOA->shape[1], MOA->shape[0], sequences);
		else  
			printMOA_Matrix(dbglevel, MOA->dimn, ident, &(MOA->indexes[k]), elm_type, MOA->shape[0] * MOA->shape[1], MOA->shape[0], sequences);
                
	}
	return 0;
}
#endif

#ifndef NDEBUG
void printMOA_dimn(int dbglevel, MOATypeDimn dimn, MOATypeShape * ind, char * * sequences, MOA_rec * MOA, int elm_type) {
	MOATypeDimn i, j;
	int ident;
	
	for (i = 0; i < MOA->shape[2] ; i++) {
		ident = printMOA_Sequences(dbglevel, 0, dimn, ind, sequences, MOA, elm_type);
		ind[2]++;
	}
	ind[2] = 0;
	for (i = 3; i < MOA->dimn; i++) {
		while (ind[i]< MOA->shape[i] - 1) {
			ind[i]++;
			for (j = i; j > 3; j--) ind[j-1] = 0;
			printMOA_dimn (dbglevel, dimn, ind, sequences, MOA, elm_type);
		}
	}
}
#endif
#ifndef NDEBUG
void printMOA (int dbglevel, MOA_rec * MOA, char * * sequences, int elm_type) { 
	long i, dimn;
	MOATypeShape * ind;
	char msg[MID_MESSAGE_SIZE];

	/* elm_type = 0 to print MOA elements else (=1) to print MOA Indexes.*/
	
	/* allocate memory for indexes of tensor.*/
	ind = NULL;
	ind =  mmalloc (((MOATypeInd) MOA->dimn) * ((MOATypeInd) sizeof *ind));
	if (ind == NULL) 
		return;
	
	/* Print sequences incase of printing MOA elements only. */
	if ((elm_type == 0) && (sequences != NULL))
		PrintSequencies (dbglevel, MOA->dimn, sequences, MOA->shape);

	if (elm_type == 0)
		sprintf(msg, "\n The Tensor of Dimn %lld & elements %lld & Shape: ", MOA->dimn, MOA->elements_ub);
	else 
		sprintf(msg, "\n The Relative Indices of Tensor of Dimn %lld & elements %lld & Shape: ", MOA->dimn, MOA->elements_ub);
	mprintf (dbglevel, msg, 1);
	for (i = 0; i < MOA->dimn; i++) {
		ind[i] = 0;
		sprintf(msg, " %lld ", MOA->shape[i]);
		mprintf (dbglevel, msg, 1);
	}
	mprintf(dbglevel, ":\n\n", 1);
	dimn = MOA->dimn - 1;
	
	printMOA_dimn(dbglevel, dimn, ind, sequences, MOA, elm_type);

	/* Free memory of indexes */
	if (ind != NULL) 
		free (ind);
}
#endif

void printMOA1 (MOA_rec * MOA, int elm_type) { 
    MOATypeInd i, j, lcnt, dim1cnt, dim2cnt, dim3cnt, * ind;
    MOATypeDimn k;
    char msg[MID_MESSAGE_SIZE];

  
    ind = NULL;
    ind =  mmalloc (((MOATypeInd) MOA->dimn) * ((MOATypeInd) sizeof *ind));
    if (ind == NULL) 
        return;

  
    if (elm_type == 0)
        sprintf(msg, "\n The Tensor of Dimn %lld & elements %lld & Shape: ", MOA->dimn, MOA->elements_ub);
    else 
        sprintf(msg, "\n The Relative Indices of Tensor of Dimn %lld & elements %lld & Shape: ", MOA->dimn, MOA->elements_ub);
    mprintf (1, msg, threadnum);
  for (i = 0; i < MOA->dimn; i++) {
    if (i == 0)
      ind[i] = MOA->shape[i];
    else
      ind[i] = MOA->shape[i] * ind[i-1];
    sprintf(msg, " %lld ", MOA->shape[i]);
		mprintf (1, msg, threadnum);
  }

  dim1cnt = dim2cnt = dim3cnt = 0;
  /* print indices for dimn 1 */
  mprintf(1, ": \n", threadnum);
  /* leave 7 spaces for the 2nd dimension indices */
  mprintf(1, "       ", threadnum);	
  /* leave 2 spaces for each < for ewach dimension in the tensor*/
  for (j = 0; j < MOA->dimn; j++) 
    mprintf(1, "  ", threadnum);
   
  for (j = 0; j < MOA->shape[0]; j++) {
    sprintf(msg, "%5ld ", dim1cnt);
		mprintf (1, msg, threadnum);
    dim1cnt ++;
  }
  /* print 3rd dimension indices */
  sprintf(msg, " %9ld", dim3cnt);
	mprintf (1, msg, threadnum);
  dim3cnt ++;	
  if (dim3cnt > MOA->shape[2] - 1)
    dim3cnt = 0;
 
  /* print separator line */
  dim1cnt = 0;
  mprintf(1, "\n", threadnum);
  /* leave 7 spaces for the 2nd dimension indices & more for the end */
  mprintf(1, "___________________", threadnum);	
  /* leave 2 spaces for each < for ewach dimension in the tensor */
  for (j = 0; j < MOA->dimn; j++) 
    mprintf(1, "__", threadnum);
   
  for (j = 0; j < MOA->shape[0]; j++) {
    mprintf(1, "______", threadnum);	
    dim1cnt ++;
  }
  mprintf(1, "\n", threadnum);
  /* print 2nd dimension indices */
  sprintf(msg, "%5ld| ", dim2cnt);
	mprintf (1, msg, threadnum);
  dim2cnt ++;	
  for (j = 0; j < MOA->dimn; j++) 
    mprintf(1, "< ", threadnum);			
  
  for (i = 0; i < MOA->elements_ub; i++) {
	if (elm_type == 0)
    		sprintf(msg, "%5lld ", MOA->elements[i].val);
	else {
            sprintf (msg, "{%5lld", MOA->indexes[i][0]);
            for (k = 1; k < MOA->dimn; k++) 
    		sprintf(msg, "%s, %5lld", msg, MOA->indexes[i][k]);
            sprintf (msg, "%s} ", msg);
        }
        mprintf (1, msg, threadnum);
    for (j = MOA->dimn - 1; j >= 0; j--) 
      if (((i+1) %  ind[j]) == 0) 
	mprintf(1, " >", threadnum);
    lcnt = 0;
    for (j = MOA->dimn - 1; j >= 0; j--) 
      if ((i+1) < MOA->elements_ub - 1) 
	if (((i+1) %  ind[j]) == 0) 
	  lcnt ++;
	
    if (lcnt > 0) {
      if (lcnt == 2) {
	/* print 3rd dimension indices */
	sprintf(msg, " %5ld", dim3cnt);
	mprintf (1, msg, threadnum);
	dim3cnt ++;	
	if (dim3cnt > MOA->shape[2] - 1)
	  dim3cnt = 0;
      }
      if (dim2cnt == 0) {
	/* print separator line */
	dim1cnt = 0;
	mprintf(1, "\n", threadnum);
	/* leave 7 spaces for the 2nd dimension indices & more for the end */
	mprintf(1, "___________________", threadnum);	
	/* leave 2 spaces for each < for ewach dimension in the tensor */
	for (j = 0; j < MOA->dimn; j++) 
	  mprintf(1, "__", threadnum);
   
	for (j = 0; j < MOA->shape[0]; j++) {
	  mprintf(1, "______", threadnum);	
	  dim1cnt ++;
	}
      }
      mprintf(1, "\n", threadnum);
      /* print 2nd dimension indices */
      sprintf(msg, "%5ld| ", dim2cnt);
			mprintf (1, msg, threadnum);
      dim2cnt ++;	
      if (dim2cnt > MOA->shape[1] - 1)
				dim2cnt = 0;
      for (j = 0; j <  (MOA->dimn - lcnt); j++)
				mprintf(1, "  ", threadnum); 		
      for (j = 0; j <  lcnt; j++)
				mprintf(1, "< ", threadnum);
    }    
  }
  mprintf(1, "\n", threadnum);
  if (ind != NULL) 
		free (ind);

}

void printMOA_scr (MOA_rec * MOA, int elm_type) { 
  MOATypeInd i, j, lcnt, dim1cnt, dim2cnt, dim3cnt, * ind;

  MOATypeDimn k;
  
  ind = NULL;
  ind =  mmalloc (((MOATypeInd) MOA->dimn) * ((MOATypeInd) sizeof *ind));
  if (ind == NULL) 
    return;

  
if (elm_type == 0)
		printf("\n The Tensor of Dimn %lld & elements %lld & Shape: ", MOA->dimn, MOA->elements_ub);
	else 
		printf("\n The Relative Indices of Tensor of Dimn %lld & elements %lld & Shape: ", MOA->dimn, MOA->elements_ub);
  for (i = 0; i < MOA->dimn; i++) {
    if (i == 0)
      ind[i] = MOA->shape[i];
    else
      ind[i] = MOA->shape[i] * ind[i-1];
    printf(" %lld ", MOA->shape[i]);
  }

  dim1cnt = dim2cnt = dim3cnt = 0;
  /* print indices for dimn 1 */
  printf(": \n");
  /* leave 7 spaces for the 2nd dimension indices */
  printf("       ");	
  /* leave 2 spaces for each < for ewach dimension in the tensor*/
  for (j = 0; j < MOA->dimn; j++) 
    printf("  ");
   
  for (j = 0; j < MOA->shape[0]; j++) {
    printf("%5ld ", dim1cnt);
		
    dim1cnt ++;
  }
  /* print 3rd dimension indices */
  printf(" %9ld", dim3cnt);
  dim3cnt ++;	
  if (dim3cnt > MOA->shape[2] - 1)
    dim3cnt = 0;
 
  /* print separator line */
  dim1cnt = 0;
  printf("\n");
  /* leave 7 spaces for the 2nd dimension indices & more for the end */
  printf("___________________");	
  /* leave 2 spaces for each < for ewach dimension in the tensor */
  for (j = 0; j < MOA->dimn; j++) 
    printf("__");
   
  for (j = 0; j < MOA->shape[0]; j++) {
    printf("______");	
    dim1cnt ++;
  }
  printf("\n");
  /* print 2nd dimension indices */
  printf("%5ld| ", dim2cnt);
	
  dim2cnt ++;	
  for (j = 0; j < MOA->dimn; j++) 
    printf("< ");			
  
  for (i = 0; i < MOA->elements_ub; i++) {
	if (elm_type == 0)
    		printf("%5lld ", MOA->elements[i].val);
        else {
            printf ("{%5lld", MOA->indexes[i][0]);
            for (k = 1; k < MOA->dimn; k++) 
    		printf(", %5lld", MOA->indexes[i][k]);
            printf ("} ");
        }
		
    for (j = MOA->dimn - 1; j >= 0; j--) 
      if (((i+1) %  ind[j]) == 0) 
	printf(" >");
    lcnt = 0;
    for (j = MOA->dimn - 1; j >= 0; j--) 
      if ((i+1) < MOA->elements_ub - 1) 
	if (((i+1) %  ind[j]) == 0) 
	  lcnt ++;
	
    if (lcnt > 0) {
      if (lcnt == 2) {
	/* print 3rd dimension indices */
	printf(" %5ld", dim3cnt);
	
	dim3cnt ++;	
	if (dim3cnt > MOA->shape[2] - 1)
	  dim3cnt = 0;
      }
      if (dim2cnt == 0) {
	/* print separator line */
	dim1cnt = 0;
	printf("\n");
	/* leave 7 spaces for the 2nd dimension indices & more for the end */
	printf("___________________");	
	/* leave 2 spaces for each < for ewach dimension in the tensor */
	for (j = 0; j < MOA->dimn; j++) 
	  printf("__");
   
	for (j = 0; j < MOA->shape[0]; j++) {
	  printf("______");	
	  dim1cnt ++;
	}
      }
      printf("\n");
      /* print 2nd dimension indices */
      printf("%5ld| ", dim2cnt);
			
      dim2cnt ++;	
      if (dim2cnt > MOA->shape[1] - 1)
				dim2cnt = 0;
      for (j = 0; j <  (MOA->dimn - lcnt); j++)
				printf("  "); 		
      for (j = 0; j <  lcnt; j++)
				printf("< ");
    }    
  }
  printf("\n");
  if (ind != NULL) 
		free (ind);

}

#ifndef NDEBUG
void printMOAIndices (MOA_rec * MOA) { 
  MOATypeInd i, j, lcnt, dim1cnt, dim2cnt, dim3cnt, * ind;
  MOATypeDimn k;
  char msg[MID_MESSAGE_SIZE];
  
  ind = NULL;
  ind =  mmalloc (((MOATypeInd) MOA->dimn) * ((MOATypeInd) sizeof *ind));
  if (ind == NULL) 
    return;

  sprintf(msg, "\n The Relative Indices of Tensor of Dimn %lld & elements %lld & Shape: ", MOA->dimn, MOA->elements_ub);
	mprintf (1, msg, threadnum);
  for (i = 0; i < MOA->dimn; i++) {
    if (i == 0)
      ind[i] = MOA->shape[i];
    else
      ind[i] = MOA->shape[i] * ind[i-1];
    sprintf(msg, " %lld ", MOA->shape[i]);
		mprintf (1, msg, threadnum);
  }

  dim1cnt = dim2cnt = dim3cnt = 0;
  /* print indices for dimn 1 */
  mprintf(1, ": \n", threadnum);
  /* leave 7 spaces for the 2nd dimension indices */
  mprintf(1, "       ", threadnum);	
  /* leave 2 spaces for each < for ewach dimension in the tensor*/
  for (j = 0; j < MOA->dimn; j++) 
    mprintf(1, "  ", threadnum);
   
  for (j = 0; j < MOA->shape[0]; j++) {
    sprintf(msg, "%5ld ", dim1cnt);
		mprintf (1, msg, threadnum);
    dim1cnt ++;
  }
  /* print 3rd dimension indices */
  sprintf(msg, " %9ld", dim3cnt);
	mprintf (1, msg, threadnum);
  dim3cnt ++;	
  if (dim3cnt > MOA->shape[2] - 1)
    dim3cnt = 0;
 
  /* print separator line */
  dim1cnt = 0;
  mprintf(1, "\n", threadnum);
  /* leave 7 spaces for the 2nd dimension indices & more for the end */
  mprintf(1, "___________________", threadnum);	
  /* leave 2 spaces for each < for ewach dimension in the tensor */
  for (j = 0; j < MOA->dimn; j++) 
    mprintf(1, "__", threadnum);
   
  for (j = 0; j < MOA->shape[0]; j++) {
    mprintf(1, "______", threadnum);	
    dim1cnt ++;
  }
  mprintf(1, "\n", threadnum);
  /* print 2nd dimension indices */
  sprintf(msg, "%5ld| ", dim2cnt);
	mprintf (1, msg, threadnum);
  dim2cnt ++;	
  for (j = 0; j < MOA->dimn; j++) 
    mprintf(1, "< ", threadnum);			
  
  for (i = 0; i < MOA->elements_ub; i++) {
    sprintf(msg, "{%5lld", MOA->indexes[i][0]);
    for (k = 1; k < MOA->dimn; k++) 
        sprintf(msg, "%s, %5lld", msg, MOA->indexes[i][k]);
    sprintf(msg, "%s} ", msg);

    mprintf (1, msg, threadnum);
    for (j = MOA->dimn - 1; j >= 0; j--) 
      if (((i+1) %  ind[j]) == 0) 
	mprintf(1, " >", threadnum);
    lcnt = 0;
    for (j = MOA->dimn - 1; j >= 0; j--) 
      if ((i+1) < MOA->elements_ub - 1) 
	if (((i+1) %  ind[j]) == 0) 
	  lcnt ++;
	
    if (lcnt > 0) {
      if (lcnt == 2) {
	/* print 3rd dimension indices */
	sprintf(msg, " %5ld", dim3cnt);
	mprintf (1, msg, threadnum);
	dim3cnt ++;	
	if (dim3cnt > MOA->shape[2] - 1)
	  dim3cnt = 0;
      }
      if (dim2cnt == 0) {
	/* print separator line */
	dim1cnt = 0;
	mprintf(1, "\n", threadnum);
	/* leave 7 spaces for the 2nd dimension indices & more for the end */
	mprintf(1, "___________________", threadnum);	
	/* leave 2 spaces for each < for ewach dimension in the tensor */
	for (j = 0; j < MOA->dimn; j++) 
	  mprintf(1, "__", threadnum);
   
	for (j = 0; j < MOA->shape[0]; j++) {
	  mprintf(1, "______", threadnum);	
	  dim1cnt ++;
	}
      }
      mprintf(1, "\n", threadnum);
      /* print 2nd dimension indices */
      sprintf(msg, "%5ld| ", dim2cnt);
			mprintf (1, msg, threadnum);
      dim2cnt ++;	
      if (dim2cnt > MOA->shape[1] - 1)
				dim2cnt = 0;
      for (j = 0; j <  (MOA->dimn - lcnt); j++)
				mprintf(1, "  ", threadnum); 		
      for (j = 0; j <  lcnt; j++)
				mprintf(1, "< ", threadnum);
    }    
  }
  mprintf(1, "\n", threadnum);
  if (ind != NULL) 
		free (ind);
}
#endif

/*****************************************************************
	Function: isLowerBorderCell
*****************************************************************/
int isLowerBorderCell (MOATypeShape * index, MOATypeDimn dimn) {
  MOATypeDimn j;
  int borderCell = 0;
  for (j = 0; j< dimn; j++) {
    if (index[j] == 0)
      borderCell = 1;
  }
  return borderCell;
}
/***************************************************************
	Function: isHigherBorderCell
***************************************************************/
int isHigherBorderCell (MOATypeShape * index, MOATypeDimn dimn, MOATypeShape * shape) {
  MOATypeDimn j;
  int borderCell = 0;

  for (j = 0; j< dimn; j++) {
    if (index[j] == (shape[j] -1))
      borderCell = 1;
  }
  return borderCell;
}
/**************************************************************************
	Function: isHigherBorderCellandNotLower
**************************************************************************/
int isHigherBorderCellandNotLower (MOATypeShape * index, MOATypeDimn dimn, MOATypeShape * shape) {
  MOATypeDimn j;
  int HigherBorderCell = 0;
  int LowerBorderCell = 0;
  for (j = 0; j< dimn; j++) {
    if (index[j] == (shape[j] -1))
      HigherBorderCell = 1;
    if (index[j] == 0)
      LowerBorderCell = 1;
  }
  if ((HigherBorderCell == 1) && (LowerBorderCell == 0))
    return 1;
  else
    return 0;
}

/* ***************************************************************
 **********************   VecIsEqual   *************************
 ***************************************************************/
int VecIsEqual (MOATypeShape * Array_1, MOATypeDimn array1_ub, MOATypeShape * Array_2, MOATypeDimn array2_ub) {
	int rslt = 1;
	MOATypeDimn i;

	if (array1_ub != array2_ub)  {
		rslt = 0; /*Not Equal*/
		return (rslt);
	}
  
	for (i = 0 ; (i < array1_ub) && (rslt); i++ ) {
		if (Array_1[i] != Array_2[i]) {
			rslt = 0;/*Not Equal*/
		}
	}
	return (rslt);
}

/* ***************************************************************
 **************************** NextIndex ************************
 ***************************************************************/

void NextIndex (MOATypeShape * shape, MOATypeDimn shape_ub, MOATypeShape * * Prev_Index) {
	MOATypeInd i, j;

	for (i = 0; i < shape_ub; i++){
		if ((*Prev_Index)[i] < shape[i] - 1) {
			(*Prev_Index)[i] += 1;
			if (i > 0) {
				for (j = i - 1; j >= 0; j--) {
					(*Prev_Index)[j] = 0;					
				}
			}
			break;
		}
	}

} /* end NextIndex */

/* ***************************************************************
 **************************** copyIndices **********************
 ***************************************************************/
void analyzeCellPosition (MOATypeShape * cellIndex, MOATypeInd flatCellIndex, MOA_rec * * rslt, MOATypeShape * global_shape) {
    MOATypeDimn k;
    
    (*rslt)->elements[flatCellIndex].cp.llOC = 0;
    (*rslt)->elements[flatCellIndex].cp.lhOC = 0;
    (*rslt)->elements[flatCellIndex].cp.glOC = 0;
    (*rslt)->elements[flatCellIndex].cp.ghOC = 0;
    for (k = 0; k< (*rslt)->dimn; k++) {
        /*Checking local position*/
        if (cellIndex[k] == 0) 
            (*rslt)->elements[flatCellIndex].cp.llOC = 1;
        if (cellIndex[k] == (*rslt)->shape[k] - 1)
            (*rslt)->elements[flatCellIndex].cp.lhOC = 1;
        /*Checking global position*/
        if ((*rslt)->indexes[flatCellIndex][k] == 0)
            (*rslt)->elements[flatCellIndex].cp.glOC = 1;
        if ((*rslt)->indexes[flatCellIndex][k] == global_shape[k] - 1)
            (*rslt)->elements[flatCellIndex].cp.ghOC = 1;
    }
}

void copyIndicesElm (MOATypeShape * offset, MOA_rec * MOA_in, MOA_rec * * rslt, MOATypeShape * ref_dimn, MOATypeShape * Last_ref_dimn, MOATypeShape * valid_index, int copyElemFlag, int analyzeLBHB) {
    MOATypeInd counter;
    MOATypeDimn k;
    int l1_finished;
#ifndef NDEBUG
    char msg[MID_MESSAGE_SIZE];
    int dbglevel = 15;
#endif

    counter = 0;
    l1_finished = 0;
    while (l1_finished == 0) {
        //(*rslt)->indexes[counter] = offset + Gamma(ref_dimn, MOA_in->dimn, MOA_in->shape, MOA_in->dimn, 1);
        (*rslt)->indexes[counter] = mmalloc (((MOATypeInd) (*rslt)->dimn) * ((MOATypeInd) sizeof *(*rslt)->indexes[counter]));
#ifndef NDEBUG
        mprintf (dbglevel, "{", 1);
#endif
        for (k=0;k<MOA_in->dimn;k++) {
            if (offset == NULL)
                (*rslt)->indexes[counter][k] = ref_dimn[k];
            else
                (*rslt)->indexes[counter][k] = offset[k] + ref_dimn[k];
#ifndef NDEBUG
            if (k > 0)
                mprintf (dbglevel, ", ", 1);
            sprintf(msg, "%lld", (*rslt)->indexes[counter][k]);
            mprintf (dbglevel, msg, 1);
#endif
        }
#ifndef NDEBUG
        mprintf (dbglevel, "}", 1);
#endif
        if (copyElemFlag == 1)
            (*rslt)->elements[counter].val = MOA_in->elements[Gamma((*rslt)->indexes[counter], MOA_in->dimn, MOA_in->shape, MOA_in->dimn, 1)].val;
        else
            (*rslt)->elements[counter].val = 0;
        (*rslt)->elements[counter].prev = NULL;
        (*rslt)->elements[counter].prev_ub = 0;
        /*Begin AnalyzeCells*/
        if (analyzeLBHB == 1) {
            analyzeCellPosition (ref_dimn, counter, rslt, MOA_in->shape);
        }
        /*Finish AnalyzeCells*/
        counter ++;
        l1_finished = VecIsEqual(ref_dimn, MOA_in->dimn, Last_ref_dimn, MOA_in->dimn);
        NextIndex(/*shape*/valid_index,/*dimn*/MOA_in->dimn,/*current index*/&ref_dimn);
    }
}

/**************************************************************************
	Function: copyIndices
**************************************************************************/

void copyIndices (MOATypeShape * offset, MOATypeShape * global_shape, MOA_rec * * rslt, MOATypeShape * ref_dimn, MOATypeShape * Last_ref_dimn, MOATypeShape * valid_index, int analyzeLBHB) {
    MOATypeInd counter;
    MOATypeDimn k;
    int l1_finished;
#ifndef NDEBUG
    char msg[MID_MESSAGE_SIZE];
    int dbglevel = 4;
#endif

    counter = 0;
    l1_finished = 0;
#ifndef NDEBUG
    mprintf (dbglevel, ">>>>Copy Indices {\n", 1);
#endif
    while (l1_finished == 0) {
        //(*rslt)->indexes[counter] = offset + Gamma(ref_dimn, (*rslt)->dimn, global_shape, (*rslt)->dimn, 1);
        (*rslt)->indexes[counter] = mmalloc (((MOATypeInd) (*rslt)->dimn) * ((MOATypeInd) sizeof *(*rslt)->indexes[counter]));
#ifndef NDEBUG
        mprintf (dbglevel, "{", 1);
#endif
        for (k=0;k<(*rslt)->dimn;k++) {
            if (offset == NULL) {
                (*rslt)->indexes[counter][k] = ref_dimn[k];
            }
            else {
                (*rslt)->indexes[counter][k] = offset[k] + ref_dimn[k];
            }
#ifndef NDEBUG
            if (k > 0)
                mprintf (dbglevel, ", ", 1);
            sprintf(msg, "%lld", (*rslt)->indexes[counter][k]);
            mprintf (dbglevel, msg, 1);
#endif
        }
#ifndef NDEBUG
        mprintf (dbglevel, "}", 1);
#endif
        (*rslt)->elements[counter].val = 0;
        (*rslt)->elements[counter].prev = NULL;
        (*rslt)->elements[counter].prev_ub = 0;
        /* Begin AnalyzeCells =======================================*/
        if (analyzeLBHB == 1) {
            analyzeCellPosition (ref_dimn, counter, rslt, global_shape);
        }
        /* Finish AnalyzeCells =======================================*/
        counter ++;
        l1_finished = VecIsEqual(ref_dimn, (*rslt)->dimn, Last_ref_dimn, (*rslt)->dimn);
        NextIndex(valid_index, (*rslt)->dimn, &ref_dimn);
    }
#ifndef NDEBUG
    mprintf (dbglevel, "\n}\n", 1);
#endif
}

/* ***************************************************************
 **********************      Take      *************************
 ***************************************************************/

int Take (MOATypeShape * ind, MOATypeDimn ind_ub, MOA_rec * MOA_in, MOA_rec * * rslt) {
    MOATypeShape * valid_index = NULL; /*how much to take exactly after fixing the direction*/
    MOATypeShape * ref_dimn = NULL; /*starting index and iterator*/
    MOATypeShape * Last_ref_dimn = NULL; /* stopping index*/
    int positive; /*direction of take decided based on index passed values*/
    MOATypeDimn i;
    MOATypeInd offset;
#ifndef NDEBUG
    char msg[MID_MESSAGE_SIZE];
    int dbglevel = 4;
#endif

    positive = 0;
    /* decide from which direction to take */
    if (ind[0] >= 0) {
            positive = 1;
    }

    (*rslt)->dimn = MOA_in->dimn;

    /* test the index parameter to form a valid index to the Psi function */
    if (MOA_in->dimn == 0) {
            valid_index = mmalloc ((MOATypeInd) sizeof *valid_index);
            ref_dimn = mmalloc ((MOATypeInd) sizeof *ref_dimn);
            Last_ref_dimn = mmalloc ((MOATypeInd) sizeof *Last_ref_dimn);
            (*rslt)->shape = mmalloc ((MOATypeInd) sizeof *((*rslt)->shape));
            if ((ind[0] > 0) && (ind[0] <= MOA_in->shape[0])) {
                    if (positive == 0) {
                            mprintf (0,"Take Error: Invalid index!\n", 1);
                            return -1;
                    }
                    valid_index[0] = ind[0];
            }
            else if ((ind[0] < 0) && (abs(ind[0]) <= MOA_in->shape[0] + 1)) {
                    if (positive == 1) {
                            mprintf (0, "Take Error: Invalid index!\n", 1);
                            return -1;
                    }
                    valid_index[0] = MOA_in->shape[0] - abs(ind[0]);
            }
            else {
                    if (  !((ind[0] >= 0) && (ind[0] <= MOA_in->shape[0])) 
                            && !((ind[0] < 0) && (abs(ind[0]) <= MOA_in->shape[0] + 1 ))) {
                            mprintf (0, "Take Error: Invalid index!\n", 1);
                            return -1;
                    }
            }
    }
    else {
        valid_index = mmalloc (((MOATypeInd) MOA_in->dimn) * ((MOATypeInd) sizeof *valid_index));
        ref_dimn = mmalloc (((MOATypeInd) MOA_in->dimn) * ((MOATypeInd) sizeof *ref_dimn)); 
        Last_ref_dimn =mmalloc (((MOATypeInd) MOA_in->dimn) * ((MOATypeInd) sizeof *Last_ref_dimn));
        (*rslt)->shape = mmalloc (((MOATypeInd) MOA_in->dimn) * ((MOATypeInd) sizeof *((*rslt)->shape)));
        for (i = 0; i < MOA_in->dimn; i++) {
            if (ind_ub > i) {
                if ((ind[i] > 0) && (ind[i] <= MOA_in->shape[i])) {
                    if (positive == 0) {
                        mprintf (0, "Take Error: Invalid index!\n", 1);
                    return -1;
                    }
                    valid_index[i] = ind[i];
                }
                else if ((ind[i] < 0) && (abs(ind[i]) <= MOA_in->shape[i] + 1)) {
                    if (positive == 1) {
                        mprintf (0, "Take Error: Invalid index, can be all negative, or all positive!\n", 1);
                        return -1;
                    }
                    valid_index[i] = MOA_in->shape[i] - abs(ind[i]);
                }
                else {
                    if (  !((ind[i] >= 0) && (ind[i] <= MOA_in->shape[i])) 
                            && !((ind[i] < 0) && (abs(ind[i]) <= MOA_in->shape[i] + 1 ) ) ) {
                        if (ind[i] >= 0)
                            valid_index[i] = MOA_in->shape[i];
                        else
                            valid_index[i] = 0;
                    }
                }
            }
            else {
                valid_index[i] = MOA_in->shape[i];
            }
        }
    }
#ifndef NDEBUG
    sprintf(msg, "\npositive %d and valid_indexes: {\n", positive);
    mprintf (dbglevel, msg, 1);
    for (i = 0; i < MOA_in->dimn; i++) {
        sprintf(msg, " %lld ", valid_index[i]);
        mprintf (dbglevel, msg, 1);
    }
    mprintf (dbglevel, "\n}\n and input ind: {\n", 1);
    for (i = 0; i < MOA_in->dimn; i++) {
        sprintf(msg, " %lld ", ind[i]);
        mprintf (dbglevel, msg, 1);
    }
    mprintf (dbglevel, "\n}\n and shape: {\n", 1);
    for (i = 0; i < MOA_in->dimn; i++) {
        sprintf(msg, " %lld ", MOA_in->shape[i]);
        mprintf (dbglevel, msg, 1);
    }
    mprintf (dbglevel, "\n}\n", 1);
#endif
    /*call the Psi function with the valid index parameter to return the required parition of the array */
    if (positive == 1) {
        if (MOA_in->dimn == 0) {
            ref_dimn[0] = 0;
            (*rslt)->shape[0] = valid_index[0];
            if (valid_index[0] > 0)
                Last_ref_dimn[0] = valid_index[0] - 1;
            else
                Last_ref_dimn[0] = 0;
        }
        else {
            for (i=0; i< MOA_in->dimn; i++) {
                ref_dimn[i] = 0; 
                (*rslt)->shape[i] = valid_index[i]; 
                if (valid_index[i] > 0)
                    Last_ref_dimn[i] = valid_index[i] - 1; 
                else
                    Last_ref_dimn[i] = 0;
            }
        }
        (*rslt)->elements_ub = Tau((*rslt)->shape, (*rslt)->dimn);
        (*rslt)->indexes = mcalloc ((*rslt)->elements_ub, (MOATypeInd) sizeof *((*rslt)->indexes));
        (*rslt)->elements = mmalloc ((*rslt)->elements_ub * (MOATypeInd) sizeof *((*rslt)->elements));
        if ((*rslt)->elements_ub == 1) 
            (*rslt)->dimn = 0;
#ifndef NDEBUG
        mprintf (dbglevel, "and  rslt->shape: {\n", 1);
        for (i = 0; i < MOA_in->dimn; i++) {
            sprintf(msg, " %lld ", (*rslt)->shape[i]);
            mprintf (dbglevel, msg, 1);
        }
        mprintf (dbglevel, "\n}\n", 1);
#endif
        copyIndicesElm(NULL, MOA_in, rslt, ref_dimn, Last_ref_dimn, valid_index, 1, 0);
    }
    else {
        if (MOA_in->dimn == 0) {
            if (abs(ind[0]) < MOA_in->shape[0]) {
                ind[0] = abs(ind[0]);
            }
            else{
                valid_index[0] = 0;
                ind[0] = MOA_in->shape[0];
            }
            (*rslt)->shape[0] = ind[0];
            ref_dimn[0] = 0;
            if (ind[0] > 0)
                Last_ref_dimn[0] = ind[0] - 1;
            else
                Last_ref_dimn[0] = 0;
        }
        else {
            for (i = 0; i < MOA_in->dimn; i++) {
                if (abs(ind[i]) < MOA_in->shape[i]) {
                    ind[i] = abs(ind[i]);
                }
                else {
                    valid_index[i] = 0;
                    ind[i] = MOA_in->shape[i];
                }
                (*rslt)->shape[i] = ind[i];
                ref_dimn[i] = 0;
                if (ind[i] > 0)
                    Last_ref_dimn[i] = ind[i] - 1;
                else
                    Last_ref_dimn[i] = 0;
            }
        }
        (*rslt)->elements_ub = Tau((*rslt)->shape, (*rslt)->dimn);
        (*rslt)->indexes = mcalloc ((*rslt)->elements_ub, (MOATypeInd) sizeof *((*rslt)->indexes));
        (*rslt)->elements = mmalloc ((*rslt)->elements_ub * (MOATypeInd) sizeof *((*rslt)->elements));
        if ((*rslt)->elements_ub == 1) 
            (*rslt)->dimn = 0;	
#ifndef NDEBUG
        mprintf (dbglevel, "and  rslt->shape: {\n", 1);
        for (i = 0; i < MOA_in->dimn; i++) {
            sprintf(msg, " %lld ", (*rslt)->shape[i]);
            mprintf (dbglevel, msg, 1);
        }
        mprintf (dbglevel, "\n}\n", 1);
#endif
        //offset = Gamma(valid_index, MOA_in->dimn, MOA_in->shape, MOA_in->dimn, 1);
        copyIndicesElm(valid_index, MOA_in, rslt, ref_dimn, Last_ref_dimn, ind, 1, 0);
    }
  
    if (valid_index != NULL) 	
        free(valid_index);
    if (ref_dimn  != NULL) 
        free(ref_dimn);
    if (Last_ref_dimn != NULL)
        free(Last_ref_dimn);

    return 0;	
}

/***************************************************************
 *******  Take without creating the whole Tensor     ***********
 ***************************************************************/
int TakeInd (MOATypeShape * ind, MOATypeDimn dimn, MOATypeShape * global_shape, MOA_rec * * rslt, int analyzeLBHB) {
	MOATypeShape * valid_index = NULL;
	MOATypeShape * ref_dimn = NULL;
	MOATypeShape * Last_ref_dimn = NULL;
	int positive;
	MOATypeInd offset;
	MOATypeDimn i;
#ifndef NDEBUG
	char msg[MID_MESSAGE_SIZE];
	int dbglevel = 4;
#endif

	
	positive = 0;
	/* decide from which direction to take*/
	if (ind[0] >= 0) {
		positive = 1;
	}
	
  
	(*rslt)->dimn = dimn;

	/* test the index parameter to form a valid index to the Psi function*/
	if (dimn == 0) {
		valid_index = mmalloc ((MOATypeInd) sizeof *valid_index);
		ref_dimn = mmalloc ((MOATypeInd) sizeof *ref_dimn);
		Last_ref_dimn = mmalloc ((MOATypeInd) sizeof *Last_ref_dimn);
		(*rslt)->shape = mmalloc ((MOATypeInd) sizeof *(*rslt)->shape );
		if ((ind[0] > 0) && (ind[0] <= global_shape[0])) {
			if (positive == 0) {
				mprintf(0,"Take Error: Invalid index!\n", 1);
				return -1;
			}
			valid_index[0] = ind[0];
		}
		else if ((ind[0] < 0) && (abs(ind[0]) <= global_shape[0] + 1)) {
			if (positive == 1) {
				mprintf(0, "Take Error: Invalid index!\n", 1);
				return -1;
			}
			valid_index[0] = global_shape[0] - abs(ind[0]);
		}
		else {
			if (  !((ind[0] >= 0) && (ind[0] <= global_shape[0])) 
				&& !((ind[0] < 0) && (abs(ind[0]) <= global_shape[0] + 1 ) ) ) {
				mprintf(0, "Take Error: Invalid index!\n", 1);
				return -1;
			}
		}
	}
	else {
		valid_index = mmalloc (((MOATypeInd) dimn) * ((MOATypeInd) sizeof *valid_index));
		ref_dimn = mmalloc (((MOATypeInd) dimn) * ((MOATypeInd) sizeof *ref_dimn)); 
		Last_ref_dimn =mmalloc ((MOATypeInd) dimn * ((MOATypeInd) sizeof *Last_ref_dimn));
		(*rslt)->shape = mmalloc (((MOATypeInd) dimn) * ((MOATypeInd) sizeof *(*rslt)->shape) );
		for (i = 0; i < dimn; i++) {
			if (dimn > i) {
				if ((ind[i] >= 0) && (ind[i] <= global_shape[i])) {
					if (positive == 0) {
						mprintf(0, "Take Error: Invalid index!\n", 1);
						return -1;
					}
					valid_index[i] = ind[i];
				}
				else if ((ind[i] < 0) && (abs(ind[i]) <= global_shape[i] + 1)) {
					if (positive == 1) {
	   				mprintf(0, "Take Error: Invalid index!\n", 1);
						return -1;
					}
					valid_index[i] = global_shape[i] - abs(ind[i]);
				}
				else {
					if (  !((ind[i] >= 0) && (ind[i] <= global_shape[i])) 
						&& !((ind[i] < 0) && (abs(ind[i]) <= global_shape[i] + 1 ) ) ) {
						if (ind[i] >= 0)
							valid_index[i] = global_shape[i];
						else
							valid_index[i] = 0;
					}
				}
			} 
			else {
				valid_index[i] = global_shape[i];
			}
		}
	}
	/* call the Psi function with the valid index parameter to return the required parition of the array*/
	if (positive == 1) {
		if (dimn == 0) {
			ref_dimn[0] = 0;
			(*rslt)->shape[0] = valid_index[0];
			if (valid_index[0] > 0)
				Last_ref_dimn[0] = valid_index[0] - 1;
			else
				Last_ref_dimn[0] = 0;
		}
		else {
			for (i=0; i< dimn; i++) {
				ref_dimn[i] = 0;
				(*rslt)->shape[i] = valid_index[i];
				if (valid_index[i] > 0)
					Last_ref_dimn[i] = valid_index[i] - 1;
				else
					Last_ref_dimn[i] = 0;
			}
		}
		(*rslt)->elements_ub = Tau((*rslt)->shape, (*rslt)->dimn);
                (*rslt)->indexes = mcalloc ((*rslt)->elements_ub, (MOATypeInd) sizeof *((*rslt)->indexes));
		(*rslt)->elements = mmalloc ((*rslt)->elements_ub * (MOATypeInd) sizeof *((*rslt)->elements));
		if ((*rslt)->elements_ub == 1) 
			(*rslt)->dimn = 0;
    
   
#ifndef NDEBUG
		mprintf (dbglevel, "\nand  rslt->shape:{\n", 1);
		for (i = 0; i < dimn; i++) {
			sprintf (msg, " %lld ", (*rslt)->shape[i]);
			mprintf (dbglevel, msg, 1);
		}
		mprintf (dbglevel, "\n}\n", 1);
#endif
		copyIndices(NULL, global_shape, rslt, ref_dimn, Last_ref_dimn, valid_index, analyzeLBHB);
	}
	else {
		if (dimn == 0) {
			if (abs(ind[0]) < global_shape[0]) {
				ind[0] = abs(ind[0]);
			}
			else{
				valid_index[0] = 0;
				ind[0] = global_shape[0];
			}
			(*rslt)->shape[0] = ind[0];
			ref_dimn[0] = 0;
			if (ind[0] > 0)
				Last_ref_dimn[0] = ind[0] - 1;
			else
				Last_ref_dimn[0] = 0;
		}
		else {
			for (i = 0; i < dimn; i++) {
				if (abs(ind[i]) < global_shape[i]) {
					ind[i] = abs(ind[i]);
				}
				else {
					valid_index[i] = 0;
					ind[i] = global_shape[i];
				}
				(*rslt)->shape[i] = ind[i];
				ref_dimn[i] = 0;
				if (ind[i] > 0)
					Last_ref_dimn[i] = ind[i] - 1;
				else
					Last_ref_dimn[i] = 0;
			}
		}
		(*rslt)->elements_ub = Tau((*rslt)->shape, (*rslt)->dimn);
		(*rslt)->indexes = mcalloc ((*rslt)->elements_ub, (MOATypeInd) sizeof *((*rslt)->indexes));
		(*rslt)->elements = mmalloc ((*rslt)->elements_ub * (MOATypeInd) sizeof *((*rslt)->elements));
		if ((*rslt)->elements_ub == 1) 
			(*rslt)->dimn = 0;

#ifndef NDEBUG
		mprintf (dbglevel, "\nand  rslt->shape: {\n", 1);
		for (i = 0; i < dimn; i++) {
			sprintf (msg, " %lld ", (*rslt)->shape[i]);
			mprintf (dbglevel, msg, 1);
		}
		mprintf (dbglevel, "\n}\n", 1);
#endif
		//offset = Gamma(valid_index, dimn, global_shape, dimn, 1);
		copyIndices(valid_index, global_shape, rslt, ref_dimn, Last_ref_dimn, ind, analyzeLBHB);
	}
  
	if (valid_index != NULL) 	
		free(valid_index);
	if (ref_dimn  != NULL) 
		free(ref_dimn);
	if (Last_ref_dimn != NULL)
		free(Last_ref_dimn);
	return 0;	
}

/***************************************************************
 **********************      Drop      *************************
 ***************************************************************/
int Drop (MOATypeShape * ind, MOATypeDimn ind_ub, MOA_rec * MOA_in, MOA_rec * * rslt) {
	MOATypeShape * valid_index = NULL;
	int positive, ret;
	MOATypeDimn i;
	char msg[MID_MESSAGE_SIZE];
#ifndef NDEBUG
	int dbglevel = 4;
#endif
	positive = 0;
  
	if (ind[0] >= 0) {
		positive = 1;
	} 
	
  
	if (MOA_in->dimn == 0) {
		valid_index = mmalloc ((MOATypeInd)  sizeof *valid_index);
		if ((ind[0] >= 0) && ( ind[0] < MOA_in->shape[0])) {
			if (positive == 0) {
				mprintf(0, "Drop Error: Invalid index!\n", 1);
				return -1;
			}
			valid_index[0] = MOA_in->shape[0] - ind[0];
		}
		else if ((ind[0] < 0) && (abs(ind[0]) <= MOA_in->shape[0])) {
			if (positive == 1) {
				mprintf(0, "Drop Error: Invalid index!\n", 1);
				return -1;
			}
			valid_index[0] = MOA_in->shape[0] - abs(ind[0]);
		}
		else {
			mprintf(0, "Drop Error: Invalid index!\n", 1);
      	return -1;
		}
	}
	else {
		valid_index = mmalloc(((MOATypeInd)  MOA_in->dimn) * ((MOATypeInd) sizeof *valid_index));
		for (i = 0; i < MOA_in->dimn; i++) {
			if (ind_ub > i) {
				if ((ind[i] >= 0) && ( ind[i] < MOA_in->shape[i])) {
					if (positive == 0) {
						sprintf(msg, "Drop Error: Invalid index at dimension %lld!\n", i);
						mprintf (0, msg, 1);
						return -1;
					}
					valid_index[i] = MOA_in->shape[i] - ind[i];
				}
				else if ((ind[i] < 0) && (abs(ind[i]) <= MOA_in->shape[i])) {
					if (positive == 1) {
						sprintf(msg, "Drop Error: Invalid index at dimension %lld!\n", i);
						mprintf (0, msg, 1);
						return -1;
					}
					valid_index[i] = MOA_in->shape[i] - abs(ind[i]);
				}
				else {
					sprintf(msg, "Drop Error: Invalid index at dimension %lld!\n", i);
					mprintf (0, msg, 1);
					return -1;
				}
			}
			else {
				valid_index[i] = MOA_in->shape[i];
			}
		}
	}
  
	if (ind_ub == 0) {
		valid_index[0] = MOA_in->shape[0] - ind[0];
	}
  
	if (positive) {
		if (MOA_in->dimn  == 0) {
			valid_index[0] = valid_index[0] * (-1);
		}
		else {
			for (i = 0; i < MOA_in->dimn; i ++) {
				valid_index[i] = valid_index[i] * (-1);
			}
		}
	}    
#ifndef NDEBUG
	mprintf (dbglevel, "\nin drop before take with valid index =  {\n", 1);
	for (i = 0; i < MOA_in->dimn; i ++) {
		sprintf(msg, " %lld ", valid_index[i]);
		mprintf (dbglevel, msg, 1);
	}
	mprintf (dbglevel, "\n}\n and  MOA_in->shape =  {\n", 1);
	for (i = 0; i < MOA_in->dimn; i ++) {
		sprintf(msg, " %lld ", MOA_in->shape[i]);
		mprintf (dbglevel, msg, 1);
	}
	mprintf (dbglevel, "\n}\n and  ind  =  {\n", 1);
	for (i = 0; i < MOA_in->dimn; i ++) {
		sprintf(msg, " %lld ", ind[i]);
		mprintf (dbglevel, msg, 1);
	}
	mprintf (dbglevel, "\n}\n", 1);
#endif
	ret = Take(valid_index, MOA_in->dimn, MOA_in, rslt);
  
	if (valid_index != NULL) 	
		free (valid_index);
  
	return ret;
}

/* ***************************************************************
 *******  Drop without creating the whole Tensor  **************
 ***************************************************************/
int DropInd (MOATypeShape * ind, MOATypeDimn dimn, MOATypeShape * shape, MOA_rec * * rslt) {
	MOATypeShape * valid_index = NULL;
	int positive, ret;
	MOATypeDimn i;
	char msg[MID_MESSAGE_SIZE];
#ifndef NDEBUG
	int dbglevel = 4;
#endif
	
	positive = 0;
	if (ind[0] >= 0) {
		positive = 1;
	} 
	
	/*Decide the valid index to take (shape - index)*/
	if (dimn == 0) {
		valid_index = mmalloc ((MOATypeInd)  sizeof *valid_index);
		if ((ind[0] >= 0) && ( ind[0] < shape[0])) {
			if (positive == 0) {
				mprintf(0, "Drop Error: Invalid index!\n", 1);
				return -1;
      			}
	      		valid_index[0] = shape[0] - ind[0];
		}
		else if ((ind[0] < 0) && (abs(ind[0]) <= shape[0])) {
			if (positive == 1) {
				mprintf(0, "Drop Error: Invalid index!\n", 1);
				return -1;
			}
      			valid_index[0] = shape[0] - abs(ind[0]);
		}
		else {
			mprintf(0, "Drop Error: Invalid index!\n", 1);
      			return -1;
		}
  	}
	else {
		valid_index = mmalloc(((MOATypeInd)  dimn) * ((MOATypeInd)  sizeof *valid_index ));
		for (i = 0; i < dimn; i++) {
			if (dimn > i) {
				if ((ind[i] >= 0) && ( ind[i] < shape[i])) {
					if (positive == 0) {
						sprintf(msg, "Drop Error: Invalid index at dimension %lld!\n", i);
						mprintf (0, msg, 1);
						return -1;
					}
					valid_index[i] = shape[i] - ind[i];
				}
				else if ((ind[i] < 0) && (abs(ind[i]) <= shape[i])) {
					if (positive == 1) {
						sprintf(msg, "Drop Error: Invalid index at dimension %lld!\n", i);
						mprintf (0, msg, 1);
						return -1;
					}
					valid_index[i] = shape[i] - abs(ind[i]);
				}
				else {
					sprintf(msg, "Drop Error: Invalid index at dimension %lld!\n", i);
					mprintf (0, msg, 1);
					return -1;
				}
			}
			else {
				valid_index[i] = shape[i];
			}
		}
	}
  	/*make the valid index negative to send it to TakeInd*/
	if (dimn == 0) {
		valid_index[0] = shape[0] - ind[0];
	}
  
	if (positive) {
		if (dimn  == 0) {
			valid_index[0] = valid_index[0] * (-1);
		}
		else {
			for (i = 0; i < dimn; i ++) {
				/**/
				valid_index[i] = valid_index[i] * (-1);
			}
		}
	}    
  
#ifndef NDEBUG
	mprintf (dbglevel, "\nin drop before take with valid index =  {\n", 1);  
	for (i = 0; i < dimn; i ++) {
		sprintf(msg, " %lld ", valid_index[i]);
		mprintf (dbglevel, msg, 1);
	}
	mprintf (dbglevel, "}\nand  shape =  {\n", 1);
	for (i = 0; i < dimn; i ++) {
		sprintf(msg, " %lld ", shape[i]);
		mprintf (dbglevel, msg, 1);
	}
	mprintf (dbglevel, "\n}\nand  ind  =  {\n", 1);
	for (i = 0; i < dimn; i ++) {
		sprintf(msg, " %lld ", ind[i]);
		mprintf (dbglevel, msg, 1);
	}
	mprintf (dbglevel, "\n}\n", 1);
#endif
	ret = TakeInd(valid_index, dimn, shape, rslt, 0);
  
	if (valid_index != NULL) 	
		free (valid_index);
  
	return ret;
}

/* ***************************************************************
 **********************    MOAGetLocalLowerNeighbors ****************
 ***************************************************************/
MOATypeInd MOAGetLocalLowerNeighbors_old (MOATypeShape * startingIndex, MOA_rec * MOA_in, MOA_rec * * rslt) {
    MOATypeDimn i, k, idx;
    MOATypeShape * ind, * ind2, * ind3;
    int l1_finished, ret, noDrop;
#ifndef NDEBUG
    char msg[MID_MESSAGE_SIZE];
    int dbglevel = 5;
#endif
    ind = NULL;
    ind =  mcalloc ((MOATypeInd)  sizeof *ind, (MOATypeInd)  MOA_in->dimn);
    if (ind == NULL) 
        return 0;

    ind2 = NULL;
    ind2 =  mcalloc ((MOATypeInd) sizeof *ind2, (MOATypeInd) MOA_in->dimn);
    if (ind2 == NULL) 
        return 0;
    ind3 = NULL;
    ind3 =  mcalloc ((MOATypeInd) sizeof *ind3, (MOATypeInd) MOA_in->dimn);
    if (ind3 == NULL) 
        return 0;
#ifndef NDEBUG
    sprintf(msg, "in MOAGetLowerNeighbors dimn %lld\n", MOA_in->dimn);  
    mprintf(dbglevel, msg, 1);
#endif
    MOA_rec * rslt2;
    noDrop = 1;
#ifndef NDEBUG
    //printf("(");
#endif
    for (i = 0; i < MOA_in->dimn; i++) {
        /*ind is the vector passed to take to decide how much to take from each dimension from the starting position, because we are getting only 1 distance away neighbors, it should be = 2, to get the current cell and it's direct neighbor*/
        ind[i] = 2;
        /*ind2 here is the vector passed to drop function to decide how much to drop from each dimension, it might be nothing and drop operation ommitted, if one of the neighbors is at the all zeros coordinates, hence the noDrop flag.*/
        if (startingIndex[i] > 0)
            ind2[i] = startingIndex[i] - 1;
        else
            return 0;
        if (ind2[i] > 0)
            noDrop = 0;
#ifndef NDEBUG
	//	printf (" %ld ", ind2[i]);
#endif
	}
#ifndef NDEBUG
	//printf (") Drop Vector from shape: (");
#endif
    if (noDrop == 0) {
        createMOAStruct (&rslt2);
        ret = Drop(ind2, MOA_in->dimn, MOA_in, &rslt2);  
        //printMOA_scr(rslt2, 1);
        if (ret >= 0) {
            createMOAStruct (rslt);
            ret = Take(ind, rslt2->dimn, rslt2, rslt);
        }  
        if (ret < 0) 
            return 0;
  
        for (i = 0; i < MOA_in->dimn; i++) {
#ifndef NDEBUG
            //printf (" %ld ", MOA_in->shape[i]);
#endif
            ind[i] = 0; /*this is the starting index of the neighbors MoA structure, and used as the iterator index to */
            ind2[i] = 1; /*this is the ending index of the neighbors MoA structure, always 1 because of getting the direct neighbors only shape -1*/
            ind3[i] = 2;/*this is the shape vector of the neighbors MoA structure, always 2 because of getting the direct neighbors only*/
        }
        l1_finished = 0;
        i = 0;
		
        if ((*rslt)->elements_ub > 0) {
            while (l1_finished == 0) {
                idx = Gamma(ind, rslt2->dimn, rslt2->shape, rslt2->dimn, 1);
                (*rslt)->indexes[i] = mcalloc ((MOATypeInd) (*rslt)->dimn, (MOATypeInd) sizeof *(*rslt)->indexes[i]);
                for (k=0;k<MOA_in->dimn;k++)
                    (*rslt)->indexes[i][k] =  rslt2->indexes[idx][k];
                i++;
                l1_finished = VecIsEqual(ind, MOA_in->dimn, ind2, MOA_in->dimn);
                NextIndex(ind3, MOA_in->dimn, &ind);
            }
        }
        deleteMOA (rslt2);  
    }
    else {
        createMOAStruct (rslt);
        ret = Take(ind, MOA_in->dimn, MOA_in, rslt);
    } 
    if (ind != NULL)
        free(ind);
    if (ind2 != NULL)
        free(ind2);
    if (ind3 != NULL)
        free(ind3);
    return (*rslt)->elements_ub-1;
}

/* ***************************************************************
 **********************    MOAGetGlobalLowerNeighbors ****************
 ***************************************************************/
MOATypeInd MOAGetGlobalLowerNeighbors_old (MOATypeShape * startingIndex, MOATypeDimn dimn, MOATypeShape * shape, MOA_rec * * rslt) {
    MOATypeDimn k, i, idx;
    MOATypeShape * ind, * ind2, * ind3;
    int l1_finished, ret, doDrop;
#ifndef NDEBUG
    char msg[MID_MESSAGE_SIZE];
    int dbglevel = 5;
#endif
    ind = NULL;
    ind =  mcalloc ((MOATypeInd) sizeof *ind, (MOATypeInd) dimn);
    if (ind == NULL) 
            return 0;

    ind2 = NULL;
    ind2 =  mcalloc ((MOATypeInd) sizeof *ind2, (MOATypeInd) dimn);
    if (ind2 == NULL) 
            return 0;

    ind3 = NULL;
    ind3 =  mcalloc ((MOATypeInd) sizeof *ind3, (MOATypeInd) dimn);
    if (ind3 == NULL) 
            return 0;

#ifndef NDEBUG
    sprintf(msg, "in MOAGetLowerNeighbors dimn %lld\n", dimn);  
    mprintf(dbglevel, msg, 1);
#endif
    MOA_rec * rslt2;
    doDrop = 1;
#ifndef NDEBUG
    //printf("(");
#endif
    for (i = 0; i < dimn; i++) {
        /*ind is the vector passed to take to decide how much to take from each dimension from the starting position, because we are getting only 1 distance away neighbors, it should be = 2, to get the current cell and it's direct neighbor*/
        if (startingIndex[i] == 0)
            ind[i] = 0;
        else 
            ind[i] = 2;
        /*ind2 here is the vector passed to drop function to decide how much to drop from each dimension, it might be nothing and drop operation ommitted, if one of the neighbors is at the all zeros coordinates, hence the noDrop flag.*/
        if (startingIndex[i] > 0)
            ind2[i] = startingIndex[i] - 1;
        else
            ind2[i] = startingIndex[i]; //return 0;
        if (ind2[i] > 0)
            doDrop = 0;
#ifndef NDEBUG
        //printf (" %ld ", ind2[i]);
#endif
    }
#ifndef NDEBUG
    //printf (") Drop Vector from shape: (");
#endif

    if (doDrop == 0) {
        createMOAStruct (&rslt2);
        ret = DropInd(ind2, dimn, shape, &rslt2);  
        //printMOA_scr(rslt2, 1);
        if (ret >= 0) {
            createMOAStruct (rslt);
            ret = TakeInd(ind, rslt2->dimn, rslt2->shape, rslt, 0); 
        }
        if (ret < 0) 
            return 0;

        for (i = 0; i < dimn; i++) {
            //printf (" %ld ", MOA_in->shape[i]);
            ind[i] = 0; /*this is the starting index of the neighbors MoA structure, and used as the iterator index to */
            if ((*rslt)->shape[i] > 0)  
                ind2[i] = (*rslt)->shape[i] - 1; /*this is the ending index of the neighbors MoA structure, always 1 because of getting the direct neighbors only shape -1*/
            else
                ind2[i] = 0;
        }
        l1_finished = 0;
        i = 0;
        if ((*rslt)->elements_ub > 0) {
            while (l1_finished == 0) {
                idx = Gamma(ind, rslt2->dimn, rslt2->shape, rslt2->dimn, 1);
                (*rslt)->indexes[i] = mcalloc ((MOATypeInd) (*rslt)->dimn, (MOATypeInd) sizeof *(*rslt)->indexes[i]);
                for (k=0;k<dimn;k++)
                    (*rslt)->indexes[i][k] =  rslt2->indexes[idx][k];
                i++;
                l1_finished = VecIsEqual(ind, dimn, ind2, dimn);
                NextIndex((*rslt)->shape, dimn, &ind);
            }
        }
        deleteMOA (rslt2);  
    }
    else {
        createMOAStruct (rslt);
        ret = TakeInd(ind, dimn, shape, rslt, 0); 
    } 
    if (ind != NULL)
        free(ind);
    if (ind2 != NULL)
        free(ind2);
    if (ind3 != NULL)
        free(ind3);
  return (*rslt)->elements_ub-1;
}
/* ***************************************************************
 **********************    MOAGetHigherNeighbors ****************
 ***************************************************************/
MOATypeInd MOAGetHigherNeighbors_old (long stride, MOATypeShape * startingIndex, MOATypeDimn dimn, MOATypeShape * shape, MOA_rec * * rslt) {
    MOATypeDimn i, idx, k;
    MOATypeShape * ind, * ind2, * ind3, * ind4;
    int l1_finished, ret, noDrop, noTake;

    ind = NULL;
    ind =  mmalloc (((MOATypeInd) dimn) * ((MOATypeInd) sizeof *ind));
    if (ind == NULL) 
        return 0;

    ind2 = NULL;
    ind2 =  mmalloc (((MOATypeInd) dimn) * ((MOATypeInd) sizeof *ind2));
    if (ind2 == NULL) 
        return 0;

    ind3 = NULL;
    ind3 =  mmalloc (((MOATypeInd) dimn) * ((MOATypeInd) sizeof *ind3));
    if (ind3 == NULL) 
        return 0;

    ind4 = NULL;
    ind4 =  mmalloc (((MOATypeInd) dimn) * ((MOATypeInd) sizeof *ind4));
    if (ind4 == NULL) 
        return 0;

    MOA_rec * rslt2;
    noDrop = 1;
    noTake = 1;
    for (i = 0; i < dimn; i++) {
        if (startingIndex[i] > shape[i] - 1)
            return 0;
        else {
            ind[i] = stride; /*how much to take (stride)*/
            ind2[i] = startingIndex[i]; /*how much to drop(whats b4 current index)*/
            if (ind2[i] > 0)
                noDrop = 0; /*There will be a dropping*/
            if (startingIndex[i] < shape[i]-1)
                noTake = 0; /*There will be taking*/
            
        }
    }
    if (noTake == 1)
        return 0;
    if (noDrop == 0) {
        createMOAStruct (&rslt2);
        ret = DropInd(ind2, dimn, shape, &rslt2);  
        if (ret >= 0) {
            createMOAStruct (rslt);
            ret = TakeInd(ind, rslt2->dimn, rslt2->shape, rslt, 0); 
            if (ret < 0) 
                return 0;
        }
        for (i = 0; i < dimn; i++) {
            ind[i] = 0;
            ind3[i] = (*rslt)->shape[i];
            ind2[i] = ind3[i] - 1;
        }
        l1_finished = 0;
        i = 0;
        if ((*rslt)->elements_ub > 0) {
            while (l1_finished == 0) {
                idx = Gamma(ind, rslt2->dimn, rslt2->shape, rslt2->dimn, 1);
                //(*rslt)->indexes[i] = mcalloc ((MOATypeInd) (*rslt)->dimn, (MOATypeInd) sizeof *(*rslt)->indexes[i]);
                for (k=0;k<dimn;k++) 
                    (*rslt)->indexes[i][k] =  rslt2->indexes[idx][k];
                analyzeCellPosition (ind, i, rslt, shape);
                i++;
                l1_finished = VecIsEqual(ind, dimn, ind2, dimn);
                NextIndex(ind3, dimn, &ind);
            }
        } 
        deleteMOA (rslt2);  
    }
    else {
        createMOAStruct (rslt);
        ret = TakeInd(ind, dimn, shape, rslt, 1);
        if (ret < 0) 
            return 0;
    }

    if (ind != NULL)
        free(ind);
    if (ind2 != NULL)
        free(ind2);
    if (ind3 != NULL)
        free(ind3);
    if (ind4 != NULL)
        free(ind4);

    return (*rslt)->elements_ub-1;
}
/*********************************************************************************
	Function: getLowerNeighbors
		Construct Partition MOA record for higher neighbors of a starting index, up to the specified stride
	Input:
		startingIndex: First cell Index of the partition
		dimn: dimension of MOA array.
		shape: array of lengths of each dimention.
		partSize: Partition Size.
	Output:
		part: MOA Partition record.
*********************************************************************************/
int getLowerNeighbors (long stride, MOATypeShape * cellIndex, MOATypeDimn dimn, MOATypeShape * shape, MOA_rec * * part) {
    MOATypeInd i;
    MOATypeDimn j;
    MOATypeShape * startingIndex = NULL;
    MOATypeShape * endingIndex = NULL;
    /* contruct MOA record ====================*/
    deleteMOA((*part));
    createMOAStruct(part);
    (*part)->dimn = dimn;
    (*part)->shape = mcalloc ((MOATypeInd) dimn, (MOATypeInd) sizeof *(*part)->shape);
    startingIndex = mcalloc ((MOATypeInd) dimn, (MOATypeInd) sizeof *startingIndex);
    if (startingIndex == NULL) {
        printf ("Error creating memory for startingIndex in getLowerNeighbors. Exiting\n");
        return -1;
    }
    endingIndex = mcalloc ((MOATypeInd) dimn, (MOATypeInd) sizeof *endingIndex);
    if (endingIndex == NULL) {
        printf ("Error creating memory for endingIndex in getLowerNeighbors. Exiting\n");
        return -1;
    }
    for (i=0;i<dimn;i++) {
        endingIndex[i] = cellIndex[i];
        startingIndex[i] = endingIndex[i] - stride +1;
        if (startingIndex[i] < 0) 
            startingIndex[i] = 0;		
        (*part)->shape[i] = endingIndex[i] - startingIndex[i] + 1;
    }
    (*part)->elements_ub = Tau((*part)->shape, (*part)->dimn);   
    (*part)->elements = mmalloc ((*part)->elements_ub * ((MOATypeInd) sizeof *(*part)->elements));
    if ((*part)->elements == NULL) {
        printf ("Error creating memory for MoA elements. Exiting\n");
        return -1;
    }
    (*part)->indexes = mmalloc ((*part)->elements_ub * ((MOATypeInd) sizeof *(*part)->indexes));
    if ((*part)->indexes == NULL) {
        printf ("Error creating memory for MoA indexes. Exiting\n");
        return -1;
    }
    for (i=0;i<(*part)->elements_ub;i++) {
        (*part)->indexes[i] = mcalloc ((MOATypeInd) dimn, (MOATypeInd) sizeof *(*part)->indexes[i]);
        if ((*part)->indexes[i] == NULL) {
            printf ("Error creating memory for MoA indexes elements. Exiting\n");
            return -1;
        }
        (*part)->elements[i].prev_ub = 0;
        (*part)->elements[i].prev = NULL;
        Gamma_Inverse(i, (*part)->shape, (*part)->dimn, &endingIndex, 1);
        for (j=0;j<dimn;j++) {
            endingIndex[j] = endingIndex[j] + startingIndex[j]; /* add the starting offset to get the Global Index*/
            (*part)->indexes[i][j] = endingIndex[j];
        }
    }	
    if (startingIndex != NULL) 
        free (startingIndex);
    startingIndex = NULL;
    if (endingIndex != NULL) 
        free (endingIndex);
    endingIndex = NULL;
    
    return 0;    
}
/*********************************************************************************
	Function: getHigherNeighbors
		Construct Partition MOA record for higher neighbors of a starting index, up to the specified stride
	Input:
		startingIndex: First cell Index of the partition
		dimn: dimension of MOA array.
		shape: array of lengths of each dimention.
		partSize: Partition Size.
	Output:
		part: MOA Partition record.
*********************************************************************************/
int getHigherNeighbors (long stride, MOATypeShape * startingIndex, MOATypeDimn dimn, MOATypeShape * shape, MOA_rec * * part) {
    MOATypeInd i;
    MOATypeDimn j;
    MOATypeShape * endingIndex = NULL;

    /* contruct MOA record ====================*/
    deleteMOA((*part));
    createMOAStruct(part);
    (*part)->dimn = dimn;
    (*part)->shape = mcalloc ((MOATypeInd) dimn, (MOATypeInd) sizeof *(*part)->shape);
    endingIndex = mcalloc ((MOATypeInd) dimn, (MOATypeInd) sizeof *endingIndex);
    if (endingIndex == NULL) {
        printf ("Error creating memory for endingIndex in getLowerNeighbors. Exiting\n");
        return -1;
    }
    for (i=0;i<dimn;i++) {
        endingIndex[i] = startingIndex[i]+ stride -1;
        if (endingIndex[i] >= shape[i]) 
            endingIndex[i] = shape[i] - 1;		
        (*part)->shape[i] = endingIndex[i] - startingIndex[i] + 1;
    }
    (*part)->elements_ub = Tau((*part)->shape, (*part)->dimn);   
    (*part)->elements = mmalloc ((*part)->elements_ub * ((MOATypeInd) sizeof *(*part)->elements));
    if ((*part)->elements == NULL) {
        printf ("Error creating memory for MoA elements. Exiting\n");
        return -1;
    }
    (*part)->indexes = mmalloc ((*part)->elements_ub * ((MOATypeInd) sizeof *(*part)->indexes));
    if ((*part)->indexes == NULL) {
        printf ("Error creating memory for MoA indexes. Exiting\n");
        return -1;
    }
    for (i=0;i<(*part)->elements_ub;i++) {
        (*part)->indexes[i] = mcalloc ((MOATypeInd) dimn, (MOATypeInd) sizeof *(*part)->indexes[i]);
        if ((*part)->indexes[i] == NULL) {
            printf ("Error creating memory for MoA indexes elements. Exiting\n");
            return -1;
        }
        (*part)->elements[i].val = 0;
        (*part)->elements[i].prev_ub = 0;
        (*part)->elements[i].prev = NULL;
        Gamma_Inverse(i, (*part)->shape, (*part)->dimn, &endingIndex, 1);
        for (j=0;j<dimn;j++) {
            endingIndex[j] = endingIndex[j] + startingIndex[j]; /* add the starting offset to get the Global Index*/
            (*part)->indexes[i][j] = endingIndex[j];
        }
    }	
    if (endingIndex != NULL) 
        free (endingIndex);
    endingIndex = NULL;
    
    return 0;    
}
/*********************************************************************************
	Function: getNeighbors
		Construct Partition MOA record for all neighbors of a starting index, up to the specified stride, by decrementing 
                the index, and getting the higher neighbors, so that it encapsulates all surrouding neighbors.
	Input:
		startingIndex: First cell Index of the partition
		dimn: dimension of MOA array.
		shape: array of lengths of each dimention.
		partSize: Partition Size.
	Output:
		part: MOA Partition record.
*********************************************************************************/

int getNeighbors (long stride, MOATypeShape * startingIndex, MOATypeDimn dimn, MOATypeShape * shape, MOA_rec * * part) {
    MOATypeDimn j;
    MOATypeShape * newIndex = NULL;
    newIndex = mcalloc ((MOATypeInd) dimn, (MOATypeInd) sizeof *newIndex);
   for (j=0;j<dimn;j++)
        newIndex[j] = startingIndex[j] - 1;;
    getHigherNeighbors (stride+1, newIndex, dimn, shape, part);
    if (newIndex != NULL) 
        free (newIndex);
    newIndex = NULL;
    
    return 0;    
}

/* ***************************************************************
 **********************    MOAGetHigherPartitions ****************
 ***************************************************************/
MOATypeInd MOAGetHigherPartitions (long pSize, long stride, MOATypeShape * startingIndex, MOATypeDimn dimn, MOATypeShape * shape, MOA_rec * * rslt, int * * validPartitions) {
    MOATypeInd i; /*Number of Neighbors*/
    MOATypeDimn k;
    /* get all higher neighbors of the  Cell Index*/
    if (getHigherNeighbors (stride, startingIndex, dimn, shape, rslt) == 0) {
        for (i=1;i<(*rslt)->elements_ub;i++) {
            (*validPartitions)[i] = 1;
            for (k=0;k<dimn;k++) {
                /*add partition Size to get the starting index*/
                if (startingIndex[k] != (*rslt)->indexes[i][k] && (*rslt)->indexes[i][k] + (((MOATypeShape) pSize-1) - ((*rslt)->indexes[i][k] % ((MOATypeShape) pSize-1))) < shape[k]-1)
                    (*rslt)->indexes[i][k] += ((MOATypeShape) pSize-1) - ((*rslt)->indexes[i][k] % ((MOATypeShape) pSize-1));
                else if (startingIndex[k] == (*rslt)->indexes[i][k]) 
                    (*rslt)->indexes[i][k] -= (*rslt)->indexes[i][k] % ((MOATypeShape) pSize-1);         
                else
                    (*validPartitions)[i] = 0;            
            }
        }
        return (*rslt)->elements_ub;
    }
    return 0;
}
