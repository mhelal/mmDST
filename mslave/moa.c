#include <stdio.h>
#include <stdlib.h>
#include "globals.h"
#include "moa.h"
#include "utils.h"

void createMOAStruct(MOA_rec * * MOA_val) {
  (*MOA_val) = (MOA_rec *) mmalloc(sizeof(MOA_rec));
  (*MOA_val)->dimn = 1;
  (*MOA_val)->shape = NULL;
  (*MOA_val)->indexes = NULL;
  (*MOA_val)->elements = NULL;
}

void createMOA(long * shape, long dimn, MOA_rec * MOA_val, int callFlag, int cellValue)
{
  long i;
  
  MOA_val->elements_ub = Tau(shape, dimn);
  if (MOA_val->elements_ub == 0) {
    mprintf(1, "The smallest type of MOA array is a scalar that requires at least one element! NO MOA type supports zero element list!\n", threadnum);
    return;
  }
  MOA_val->dimn = dimn;
  if (dimn == 0) {
    MOA_val->shape = (long *) mcalloc (1, sizeof(long));
    MOA_val->shape[0] = 1;
  }
  else {
    MOA_val->shape = (long *) mcalloc(MOA_val->dimn, sizeof(long));
    for (i = 0; i<= dimn; i++) {
      MOA_val->shape[i] = shape[i];
    }
  }
  /* if callFlag is -1 means, don't allocate memory for the tensor, this is used with sparse arrays implementation. */

  if (callFlag != -1) {
    
    /*MOA_val->elements_ub =  Tau(shape, dimn); */
    MOA_val->elements = (MOA_elm * ) mcalloc (MOA_val->elements_ub, sizeof(MOA_elm));
    MOA_val->indexes = (unsigned long * ) mcalloc (MOA_val->elements_ub, sizeof(unsigned long));    /* if callFlag is 0, means initialize cell values with 0, which is done with mcalloc */
    if (callFlag == 1) {
      for (i = 0; i< MOA_val->elements_ub; i++) {
			/* if callFlag is 1 meaning, initialize the array values with the cellValue argument value */
				MOA_val->elements[i].prev_ub = 0;	
				MOA_val->elements[i].prev = NULL;
				MOA_val->indexes[i] = i;
				if (callFlag == 1) {
				  MOA_val->elements[i].val = cellValue;
				}
				else { /* else initialize the cell values with their flat index value */
				  MOA_val->elements[i].val = i;
				}
      }
    }
  }
}

void printMOA_old (MOA_rec * MOA) { 
  long i;
	char msg[MID_MESSAGE_SIZE];

  mprintf(1, "The Alignment Matrix:\n", threadnum);
  mprintf(1,"< ", threadnum);
  for (i=0; i<MOA->elements_ub;i++) {
    sprintf(msg, "%4ld  ", MOA->elements[i].val);
		mprintf (1, msg, threadnum);
    if (((i+1) %  MOA->shape[ MOA->dimn - 1]) == 0) {
      mprintf(1," >\n", threadnum);
      if ((i+1)<MOA->elements_ub)
	mprintf(1,"< ", threadnum);
    }
  }

}

int isLowerBorderCell (long * index, long dimn) {
  long j;
  int borderCell = 0;
  for (j = 0; j< dimn; j++) {
    if (index[j] == 0)
      borderCell = 1;
  }
  return borderCell;
}

int isHigherBorderCell (long * index, long dimn, long * shape) {
  long j;
	char msg[MID_MESSAGE_SIZE];
  int borderCell = 0;

  for (j = 0; j< dimn; j++) {
    if (index[j] == (shape[j] -1))
      borderCell = 1;
  }
  return borderCell;
}

int isHigherBorderCellandNotLower (long * index, long dimn, long * shape) {
  long j;
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


void printMOA_old2 (MOA_rec * MOA) { 
  long i, j, lcnt;
  long * stride = NULL; 
	char msg[MID_MESSAGE_SIZE];
  
  stride = (long *) mmalloc (MOA->dimn * sizeof(long));
  sprintf(msg, "\n The Tensor of Dimn %ld & elements %ld & Shape: ", MOA->dimn, MOA->elements_ub);
	mprintf (1, msg, threadnum);
  for (i = 0; i < MOA->dimn; i++) {
    if (i == 0)
      stride[i] = MOA->shape[i];
    else
      stride[i] = MOA->shape[i] * stride[i-1];
    sprintf(msg, " %ld ",MOA->shape[i]);
		mprintf (1, msg, threadnum);
  }
  mprintf(1,"\n", threadnum);
  for (j = 0; j < MOA->dimn; j++) 
    mprintf(1,"< ", threadnum);      
 
  for (i = 0; i < MOA->elements_ub; i++) {
    sprintf(msg, "%5ld  ", MOA->elements[i].val);
		mprintf (1, msg, threadnum);
    for (j = MOA->dimn - 1; j >= 0; j--)      
      if (((i+1) %  stride[j]) == 0) 
	mprintf(1," >", threadnum);          		   
    lcnt=0;
    for (j = MOA->dimn - 1; j >= 0; j--) 
      if ((i+1) < MOA->elements_ub - 1) 
	if (((i+1) %  stride[j]) == 0) 
	  lcnt ++;      
    if (lcnt > 0) {
      mprintf(1,"\n", threadnum);
      for (j = 0; j <  (MOA->dimn - lcnt); j++)
	mprintf(1,"  ", threadnum); 		
      for (j = 0; j <  lcnt; j++)
	mprintf(1,"< ", threadnum);
    }    
  }
  
  mprintf(1, "\n", threadnum);
  if (stride != NULL)	
    free(stride);
}
void printMOA (MOA_rec * MOA) { 
  long i, j, lcnt, dim1cnt, dim2cnt, dim3cnt;
  long * stride = NULL;
  char msg[MID_MESSAGE_SIZE];
  stride = (long *) mmalloc (MOA->dimn * sizeof(long));
  
  sprintf(msg, "\n The Tensor of Dimn %ld & elements %ld & Shape: ", MOA->dimn, MOA->elements_ub);
	mprintf (1, msg, threadnum);
  for (i = 0; i < MOA->dimn; i++) {
    if (i == 0)
      stride[i] = MOA->shape[i];
    else
      stride[i] = MOA->shape[i] * stride[i-1];
    sprintf(msg, " %ld ", MOA->shape[i]);
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
    sprintf(msg, "%5ld ", MOA->elements[i].val);
		mprintf (1, msg, threadnum);
    for (j = MOA->dimn - 1; j >= 0; j--) 
      if (((i+1) %  stride[j]) == 0) 
	mprintf(1, " >", threadnum);
    lcnt = 0;
    for (j = MOA->dimn - 1; j >= 0; j--) 
      if ((i+1) < MOA->elements_ub - 1) 
	if (((i+1) %  stride[j]) == 0) 
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
  //if (stride != NULL)	
    //free(stride);
}

void printMOAIndices (MOA_rec * MOA) { 
  long i, j, lcnt, dim1cnt, dim2cnt, dim3cnt;
  long * stride = NULL;
  char msg[MID_MESSAGE_SIZE];
  stride = (long *) mmalloc (MOA->dimn * sizeof(long));
  
  sprintf(msg, "\n The Relative Indices of Tensor of Dimn %ld & elements %ld & Shape: ", MOA->dimn, MOA->elements_ub);
	mprintf (1, msg, threadnum);
  for (i = 0; i < MOA->dimn; i++) {
    if (i == 0)
      stride[i] = MOA->shape[i];
    else
      stride[i] = MOA->shape[i] * stride[i-1];
    sprintf(msg, " %ld ", MOA->shape[i]);
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
    sprintf(msg, "%5ld ", MOA->indexes[i]);
		mprintf (1, msg, threadnum);
    for (j = MOA->dimn - 1; j >= 0; j--) 
      if (((i+1) %  stride[j]) == 0) 
	mprintf(1, " >", threadnum);
    lcnt = 0;
    for (j = MOA->dimn - 1; j >= 0; j--) 
      if ((i+1) < MOA->elements_ub - 1) 
	if (((i+1) %  stride[j]) == 0) 
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
  //if (stride != NULL)	
    //free(stride);
}

void deleteMOA(MOA_rec * MOA)
{
  long i;
  char msg[MID_MESSAGE_SIZE];

  if (MOA != NULL) {	
	  mprintf(10, "deleteMOA 1\n", threadnum);
    if (MOA->elements != NULL) {
		  mprintf(10, "deleteMOA 2\n", threadnum);
      for (i = 0; i< MOA->elements_ub; i++) {
			  mprintf(10, "deleteMOA 3\n", threadnum);
				if (MOA->elements[i].prev != NULL) {
				  free(MOA->elements[i].prev);
        }
      }
      free(MOA->elements);
		  mprintf(10, "deleteMOA 4\n", threadnum);
    }
    if (MOA->shape != NULL)
      free(MOA->shape);
	  mprintf(10, "deleteMOA 5\n", threadnum);
    if (MOA->indexes != NULL) {
      free(MOA->indexes);
    }
	  mprintf(10, "deleteMOA 6\n", threadnum);
    
    //free(MOA);
	  mprintf(10, "deleteMOA 7\n", threadnum);
  }
}

/* ***************************************************************
// **************************  Tau  ******************************
// ***************************************************************/
/*Number of elements in an MOA structure based on dimn and shape - regular MOA*/
long Tau (long * array_in, long array_ub)
{
  long rslt, i;
  rslt = array_in[0];
  for ( i = 1; i< array_ub; i++) {
    rslt = rslt * array_in[i];	
  }
  return (rslt);
} /* end Tau */

/* ***************************************************************
// **********************  Gamma        **************************
// ***************************************************************/
long Gamma (long * ind, long ind_ub, long * arr_shape, long shape_ub, int Front) {
  long rslt, i;
  rslt = ind[shape_ub - 1];
  for (i = shape_ub - 2; i>=0;i--) 
    /*rslt = rslt + (ind[i] * arr_shape[i]);*/
    rslt = (rslt * arr_shape[i]) + ind[i];
  return rslt;
}

long Gamma_old (long * ind, long ind_ub, long * arr_shape, long shape_ub, int Front) {
  long * new_ind = NULL;
  long rslt, i;
  if (shape_ub == 0) {
    new_ind = (long *) mmalloc(sizeof(long));
  }
  else {
    new_ind =  (long *) mmalloc(sizeof(long) * shape_ub);
  }
  
  if (Front == 1) {
    if (ind_ub == 0) {
      new_ind[0] = ind[0];
    }
    else {
      for (i = shape_ub - ind_ub; i< shape_ub; i++) {
	new_ind[i] = ind[i];
      }
    }
    /* fill the empty dimensions with zeros */
    if (ind_ub < shape_ub) {
      if (ind_ub + 1 == shape_ub) {
	new_ind[ind_ub] = 0;
      }
      else {
	for (i = 0; i < shape_ub - ind_ub; i++) {
	  new_ind[i] = 0;
	}
      }
    }
  }
  else {
    if (ind_ub == 0) {
      new_ind[0] = ind[0];
    }
    else {
      for (i = 0; i< ind_ub; i++) {
	new_ind[i] = ind[i];
      }
    }
    /* fill the empty dimensions with zeros */
    if (ind_ub < shape_ub) {
      if (ind_ub + 1 == shape_ub) {
	new_ind[ind_ub] = 0;
      }
      else {
	for (i = ind_ub; i< shape_ub; i++) {
	  new_ind[i] = 0;
	}
      }
    }
  }
  
  if (IsValidIndex(arr_shape, shape_ub, new_ind) == -1) {
    if (new_ind != NULL)	
      free(new_ind);
    return (-1);
  }
  else {
    rslt = new_ind[0];
    for (i = 1; i< shape_ub; i++) {
      /*for (i = shape_ub - 1; i>= 1; i--) { */
      rslt = (rslt * arr_shape[i]) + new_ind[i];
    }
    if (new_ind != NULL)	
      free(new_ind);
    return (rslt);
  }
} /* end Gamma*/
/***************************************************************
// **********************  Gamma Inverse************************
// *************************************************************/

int Gamma_Inverse (long ind, long * arr_shape, long shape_ub, long * rslt)
{
  long i, ind_temp;
  div_t div_result;
  if (ind >= Tau(arr_shape, shape_ub)) {
    mprintf(1, "The Index is beyond the upper bound of the array of the defined shape!", threadnum);
    return -1;
  }
  /*for (i = 0; i< shape_ub; i++) {
  //  rslt[i] = 0;
  }*/
  ind_temp = ind;
  /*for (i = shape_ub - 1; i>= 0; i --) { */
  for (i = 0; i< shape_ub; i ++) {
    div_result = div(ind_temp, arr_shape[i]);
    rslt[i] = div_result.rem; /* ind_temp Mod arr_shape->GetAt(i) */
    ind_temp = div_result.quot; /* ind_temp \ arr_shape->GetAt(i); */
  }
  return 0;
} /* end Gamma_Inverse */

/* ***************************************************************
// ****************************** PSI ****************************
// ***************************************************************/

void Psi (long * ind, long ind_ub, MOA_rec * MOA_in, MOA_rec * rslt) {
  long offset, i;
  long j;
	char msg[MID_MESSAGE_SIZE];
  
  if (ind_ub > MOA_in->dimn) {
    mprintf(1, "Psi Error: Invalid index dimension. Psi Can't Take index of dimension greater than the Array dimension!\n", threadnum);
    return;
  }
  else if (ind_ub == 0) {
    if (ind[0] > MOA_in->shape[0] - 1) {
      mprintf  (1, "Psi Error: Invalid index at dimension: 0\n", threadnum);
      return;
    }
  }
  else {
    for (j=0;j<ind_ub;j++) {
      if (ind[j] > MOA_in->shape[j] - 1) {
				sprintf(msg, "Psi Error: Invalid index at dimension: %ld\n", j);
				mprintf (1, msg, threadnum);
				return;
      }
    }
  }
  
  if (ind_ub < MOA_in->dimn) {
    rslt->shape = (long *) mmalloc (sizeof(long) * (MOA_in->dimn - ind_ub));
    rslt->dimn = MOA_in->dimn - ind_ub;
    if (rslt->dimn == 0){
      rslt->shape[0] = MOA_in->shape[MOA_in->dimn - 1];
    }
    else {
      for (j = 0; j < rslt->dimn; j++) {
	rslt->shape[j] = MOA_in->shape[(MOA_in->dimn - 1) - j];
      }
    }
  }
  else {
    rslt->shape = (long *) mmalloc(sizeof(long));
    rslt->shape[0] = 1;
    rslt->dimn = 0;
  }
  
  rslt->elements_ub = Tau (rslt->shape, rslt->dimn);
  rslt->elements = (MOA_elm *) mmalloc(sizeof(MOA_elm) *rslt->elements_ub);
  offset = Gamma (ind, ind_ub, MOA_in->shape, MOA_in->dimn, 0);
  for (i = 0; i < rslt->elements_ub;i++) {
    rslt->elements[i].val = MOA_in->elements[offset + i].val;
  }
  
} /* end Psi */

/* ***************************************************************
// **********************   VecIsEqual   *************************
// ***************************************************************/
int VecIsEqual (long * Array_1, long array1_ub, long * Array_2, long array2_ub) 
{
  long rslt, i;
  if (array1_ub != array2_ub) 
    rslt = 0;
  else
    rslt = 1;
  
  for (i = 0 ; (i < array1_ub) && (rslt); i++ ) {
    if (Array_1[i] != Array_2[i]) {
      rslt = 0;
    }
  }
  
  return (rslt);
}

/* ***************************************************************
// **********************  GetLowerNeighbors  ********************
// ***************************************************************/

long GetLowerNeighbors_test(long seqNum, long * seqLen,  long * m_index, long * * lNeighbors) {
  long score, i, j, k;
  long lnCount; /*lower neighbor processed Counter*/
  long combnum; /* number of combinations for decrementing i number of indices*/
  long * * combin = NULL; /* combinations of lower neighbor indices*/
  long flatIndex;
  long validIndex;
  char msg[MID_MESSAGE_SIZE];
 
  score = 0;
  lnCount = 0;
  (*lNeighbors) = (long *) mcalloc (1, sizeof(long));
  /*process all decremented index neighbor */
  validIndex = 0;
  for (k=0; k<seqNum; k++) {
    m_index[k]--;
    if ( m_index[k] < 0) {
      validIndex = 1;
    }
  }
  if (validIndex == 0) { /* add it to lower neighbors list */
    
  	mprintf(3, "The lower neighbor index: ", threadnum);
 
    for (k=0; k<seqNum; k++) {
			sprintf(msg, "%ld ", m_index[k]);
			mprintf (3, msg, threadnum);
    }
    flatIndex = Gamma(m_index, seqNum, seqLen, seqNum, 1);
    (*lNeighbors)[lnCount] = flatIndex;
    /*(*dd)[lnCount] = seqNum;
    (*ddpos)[lnCount] = (int *) mcalloc (seqNum, sizeof(int));
    */
    
    sprintf(msg," ln %ld ind %ld ", lnCount, flatIndex);
		mprintf (3, msg, threadnum);
    lnCount ++;
  }
  for (k=0; k<seqNum; k++) {
    /*(*ddpos)[lnCount][k] = k; */
    m_index[k]++;
  }
  
  /* process remaining neighbors by starting by decrementing only one index entry at a time and try all different locations, then 2 locations at a time, and so forth till decrement all*/
  for (i = 1; i< seqNum; i++)  {
    /* number of the different combinations of which indeces in the multidimensional index to decrement so that all lower neighbors can be processed */
    combnum = (Factorial(seqNum)/(Factorial(seqNum-i) * Factorial(i)));
    /* create memory for combinations matrix */
    combin = (long * *) mcalloc (combnum, sizeof(long *));
    for (k=0; k<combnum; k++)
      combin[k] = (long *) mcalloc (i, sizeof(long));
    if (combin == NULL) {
      mprintf(3,"Can not allocate memory to combinations matrix\n", threadnum);
      return -1;
    }
    /* get matrix of all possible combinations */
    Combinations(seqNum, i, &combin);
    /* loop to decrement the selected indices */
    for (j=0; j<combnum; j++) {
      validIndex = 0;
      for (k=0; k<i; k++) {
	m_index[combin[j][k]-1]--;
	if ( m_index[combin[j][k]-1] < 0) {
	  validIndex = 1; /* false */
	}
      }
      if (validIndex == 0) { /*  add it to lower neighbors list*/
			  mprintf(3, "The lower neighbor index: ", threadnum);
	  	for (k=0; k<seqNum; k++) {
	    	sprintf(msg, "%ld ", m_index[k]);
				mprintf (3, msg, threadnum);
	  	}
			flatIndex = Gamma(m_index, seqNum, seqLen, seqNum, 1);
			(*lNeighbors) = (long *) realloc ((*lNeighbors), (lnCount+1) * sizeof(long));
     	if ((*lNeighbors) == NULL) {
				sprintf(msg, "Could not reallocate memory for Lower Neighbors %ld!\n", (lnCount+1));
				mprintf (3, msg, threadnum);
				return -1;
      }
	(*lNeighbors)[lnCount] = flatIndex;
	lnCount++; /*increment lower neighbor counter */
	  sprintf(msg," ln %ld ind %ld ", lnCount, flatIndex);
		mprintf (3, msg, threadnum);
     }
      for (k=0; k<i; k++) {
	m_index[combin[j][k]-1]++;
	}
    }
    /*free combinations matrix*/
    if (combin != NULL) 
		for (k=0; k<combnum; k++) {
	        if (combin[k] != NULL)	
			  free(combin[k]);
		}
    free(combin);
  }
  return lnCount;
}

/* ***************************************************************
// **********************  IsValidIndex **************************
// ***************************************************************/

int IsValidIndex (long * shape, long shape_ub, long * index) {
  long i;
	if (shape_ub == 0) {
		if ((index[0] < 0) || (index[0] >= shape[0])) 
			return (-1);
	}
	else {
		for (i = 0; i < shape_ub; i++) {
			if ((index[i] < 0) || (index[i] >= shape[i])) 
				return (-1);
		}
	}
	return (1);
}
/* ***************************************************************
// **************************** NextIndex ************************
// ***************************************************************/

void NextIndex (long * shape, long shape_ub, long * Prev_Index) {
  long i, j;
  long * rslt = NULL;
  rslt = (long *) mmalloc( sizeof(long) * shape_ub);

  /* copy the previous index */
  for (i = 0; i < shape_ub; i++){
    rslt[i] = Prev_Index[i];
  }

  /*for (i = shape_ub - 1; i>= 0; i--){ */
  for (i = 0; i < shape_ub; i++){
    if (Prev_Index[i] < shape[i] - 1) {
      rslt[i] = Prev_Index[i] + 1;
      
      if (i > 0) {
	for (j = i - 1; j >= 0; j--) {
	  rslt[j] = 0;					
	}
      }
      
      break;
    }
  }

  for (i = 0; i< shape_ub; i++){
    Prev_Index[i] = rslt[i];
  }
  if (rslt != NULL)	
    free(rslt);
}/* end NextIndex */
/* ***************************************************************
// **********************      Take      *************************
// ***************************************************************/

int Take (long * ind, long ind_ub, MOA_rec * MOA_in, MOA_rec  * rslt)
{
  long * valid_index = NULL;
  long * ref_dimn = NULL;
  long * Last_ref_dimn = NULL;
  long negative, positive, l1_finished;
  long counter;
  long offset, i;
  char msg[MID_MESSAGE_SIZE];

  negative = 0;
  positive = 0;
  /* decide from which direction to take */
  if (ind[0] >= 0) {
    positive = 1;
  }
  else {
    negative = 1;
  }
  
  rslt->dimn = MOA_in->dimn;
  
  /* test the index parameter to form a valid index to the Psi function */
  if (MOA_in->dimn == 0) {
    valid_index = (long *) mmalloc (sizeof(long));
    ref_dimn = (long *) mmalloc (sizeof(long));
    Last_ref_dimn = (long * ) mmalloc (sizeof(long));
    rslt->shape = (long *) mmalloc (sizeof(long));
    if ((ind[0] > 0) && (ind[0] <= MOA_in->shape[0])) {
      if (negative == 1) {
	mprintf(3,"Take Error: Invalid index!\n", threadnum);
	return -1;
      }
      valid_index[0] = ind[0];
    }
    else if ((ind[0] < 0) && (abs(ind[0]) <= MOA_in->shape[0] + 1)) {
      if (positive == 1) {
	mprintf(3, "Take Error: Invalid index!\n", threadnum);
	return -1;
      }
      valid_index[0] = MOA_in->shape[0] - abs(ind[0]);
    }
    else {
      if (  !((ind[0] >= 0) && (ind[0] <= MOA_in->shape[0])) 
	    && !((ind[0] < 0) && (abs(ind[0]) <= MOA_in->shape[0] + 1 ) ) ) {
	mprintf(3, "Take Error: Invalid index!\n", threadnum);
	return -1;
      }
    }
  }
  else {
    valid_index = (long *) mmalloc (MOA_in->dimn * sizeof(long));
    ref_dimn = (long *) mmalloc (MOA_in->dimn * sizeof(long)); 
    Last_ref_dimn =(long *) mmalloc (MOA_in->dimn * sizeof(long));
    rslt->shape = (long *) mmalloc (MOA_in->dimn * sizeof(long));
    for (i = 0; i < MOA_in->dimn; i++) {
      if (ind_ub > i) {
	if ((ind[i] > 0) && (ind[i] <= MOA_in->shape[i])) {
	  if (negative == 1) {
	    mprintf(3, "Take Error: Invalid index!\n", threadnum);
	    return -1;
	  }
	  valid_index[i] = ind[i];
	}
	else if ((ind[i] < 0) && (abs(ind[i]) <= MOA_in->shape[i] + 1)) {
	  if (positive == 1) {
	    mprintf(3, "Take Error: Invalid index!\n", threadnum);
	    return -1;
	  }
	  
	  valid_index[i] = MOA_in->shape[i] - abs(ind[i]);
	}
	else {
	  if (  !((ind[i] >= 0) && (ind[i] <= MOA_in->shape[i])) 
		&& !((ind[i] < 0) && (abs(ind[i]) <= MOA_in->shape[i] + 1 ) ) ) {
	    /*mprintf(1, "Take Error: Invalid index!\n", threadnum);
	    //return -1; */
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
  
  l1_finished = 0;
  counter = 0;
  
  sprintf(msg, "\npositive %ld and valid_index:", positive);
	mprintf (7, msg, threadnum);
  for (i = 0; i < MOA_in->dimn; i++) {
    sprintf(msg, " %ld ", valid_index[i]);
		mprintf (7, msg, threadnum);
  }
  mprintf(7, " and input ind:", threadnum);
  for (i = 0; i < MOA_in->dimn; i++) {
    sprintf(msg, " %ld ", ind[i]);
		mprintf (7, msg, threadnum);
  }
  mprintf(7, " and shape:", threadnum);
  for (i = 0; i < MOA_in->dimn; i++) {
    sprintf(msg, " %ld ", MOA_in->shape[i]);
		mprintf (7, msg, threadnum);
  }
  /*call the Psi function with the valid index parameter to return the required parition of the array */
  if (positive == 1) {
    if (MOA_in->dimn == 0) {
      ref_dimn[0] = 0;
      Last_ref_dimn[0] = valid_index[0] - 1;
      rslt->shape[0] = valid_index[0];
    }
    else {
      for (i=0; i< MOA_in->dimn; i++) {
				ref_dimn[i] = 0;
				Last_ref_dimn[i] = valid_index[i] - 1;
				rslt->shape[i] = valid_index[i];
      }
    }
    rslt->elements_ub = Tau(rslt->shape, rslt->dimn);
    rslt->indexes = (unsigned long * ) mcalloc (rslt->elements_ub, sizeof(unsigned long));
    rslt->elements = (MOA_elm * ) mcalloc (rslt->elements_ub, sizeof(MOA_elm));
    if (rslt->elements_ub == 1) 
      rslt->dimn = 0;
    counter = 0;
    l1_finished = 0;
    
   
    mprintf(7, " and  rslt->shape:", threadnum);
    for (i = 0; i < MOA_in->dimn; i++) {
      sprintf(msg, " %ld ", rslt->shape[i]);
			mprintf (7, msg, threadnum);
    }
    offset = Gamma(valid_index, MOA_in->dimn, MOA_in->shape, MOA_in->dimn, 1);
    sprintf(msg, "\n  offset = %ld  & indices are : ", offset);
		mprintf (7, msg, threadnum);
    
    while (l1_finished == 0) {
      rslt->indexes[counter] = Gamma(ref_dimn, MOA_in->dimn, MOA_in->shape, MOA_in->dimn, 1);
      rslt->elements[counter] = MOA_in->elements[rslt->indexes[counter]];
      rslt->elements[counter].prev = NULL;
      rslt->elements[counter].prev_ub = 0;
      
      sprintf(msg, "  %ld ", rslt->indexes[counter]);
			mprintf (7, msg, threadnum);
      
      counter ++;
      l1_finished = VecIsEqual(ref_dimn, MOA_in->dimn, Last_ref_dimn, MOA_in->dimn);
      NextIndex(valid_index, MOA_in->dimn, ref_dimn);
    }
    
  }
  
  else {
    if (MOA_in->dimn == 0) {
      if (abs(ind[0]) < MOA_in->shape[0]) {
	ind[0] = abs(ind[0]);
	rslt->shape[0] = ind[0];
	ref_dimn[0] = 0;
	Last_ref_dimn[0] = ind[0] - 1;
      }
      else{
	valid_index[0] = 0;
	ind[0] = MOA_in->shape[0];
	rslt->shape[0] = ind[0];
	ref_dimn[0] = 0;
	Last_ref_dimn[0] = ind[0] - 1;
      }
    }
    else {
      for (i = 0; i < MOA_in->dimn; i++) {
	if (abs(ind[i]) < MOA_in->shape[i]) {
	  ind[i] = abs(ind[i]);
	  rslt->shape[i] = ind[i];
	  ref_dimn[i] = 0;
	  Last_ref_dimn[i] = ind[i] - 1;
	}
	else {
	  valid_index[i] = 0;
	  ind[i] = MOA_in->shape[i];
	  rslt->shape[i] = ind[i];
	  ref_dimn[i] = 0;
	  Last_ref_dimn[i] = ind[i] - 1;
	}
      }
    }
    rslt->elements_ub = Tau(rslt->shape, rslt->dimn);
    rslt->indexes = (unsigned long * ) mcalloc (rslt->elements_ub, sizeof(unsigned long));
    rslt->elements = (MOA_elm * ) mcalloc (rslt->elements_ub, sizeof(MOA_elm));
    if (rslt->elements_ub == 1) 
      rslt->dimn = 0;
  
    mprintf(7, " and  rslt->shape:", threadnum);
    for (i = 0; i < MOA_in->dimn; i++) {
      sprintf(msg, " %ld ", rslt->shape[i]);
			mprintf (7, msg, threadnum);
    }
    
    
    counter = 0;
    l1_finished = 0;
    offset = Gamma(valid_index, MOA_in->dimn, MOA_in->shape, MOA_in->dimn, 1);
    
    sprintf(msg, "\n  offset %ld  = & indices are : ", offset);
		mprintf (7, msg, threadnum);
      
    
    while (l1_finished == 0) {
      rslt->indexes[counter] = offset + Gamma(ref_dimn, MOA_in->dimn, MOA_in->shape, MOA_in->dimn, 1);
      rslt->elements[counter] = MOA_in->elements[rslt->indexes[counter]];	    	    
      rslt->elements[counter].prev = NULL;
      rslt->elements[counter].prev_ub = 0;
    
      sprintf(msg, "  %ld ", rslt->indexes[counter]);
			mprintf (7, msg, threadnum);
      
      counter ++;
      l1_finished = VecIsEqual(ref_dimn, MOA_in->dimn, Last_ref_dimn, MOA_in->dimn);
      NextIndex(ind, MOA_in->dimn, ref_dimn);
    }
  }
  
  if (valid_index != NULL) 	
    free(valid_index);
  if (ref_dimn  != NULL) 
    free(ref_dimn);
  if (Last_ref_dimn != NULL)
    free(Last_ref_dimn);
  return 0;	
}
/* ***************************************************************
// *******  Take without creating the whole Tensor     ***********
// ***************************************************************/
int TakeInd (long * ind, long ind_ub, MOA_rec * MOA_in, MOA_rec  * rslt) {
  long * valid_index = NULL;
  long * ref_dimn = NULL;
  long * Last_ref_dimn = NULL;
  long negative, positive, l1_finished;
  long counter;
  long offset, i;
  char msg[MID_MESSAGE_SIZE];

  negative = 0;
  positive = 0;
  /* decide from which direction to take */
  if (ind[0] >= 0) {
    positive = 1;
  }
  else {
    negative = 1;
  }
  
  rslt->dimn = MOA_in->dimn;
  
  /* test the index parameter to form a valid index to the Psi function */
  if (MOA_in->dimn == 0) {
    valid_index = (long *) mmalloc (sizeof(long));
    ref_dimn = (long *) mmalloc (sizeof(long));
    Last_ref_dimn = (long * ) mmalloc (sizeof(long));
    rslt->shape = (long *) mmalloc (sizeof(long));
    if ((ind[0] > 0) && (ind[0] <= MOA_in->shape[0])) {
      if (negative == 1) {
	mprintf(1,"Take Error: Invalid index!\n", threadnum);
	return -1;
      }
      valid_index[0] = ind[0];
    }
    else if ((ind[0] < 0) && (abs(ind[0]) <= MOA_in->shape[0] + 1)) {
      if (positive == 1) {
	mprintf(1, "Take Error: Invalid index!\n", threadnum);
	return -1;
      }
      valid_index[0] = MOA_in->shape[0] - abs(ind[0]);
    }
    else {
      if (  !((ind[0] >= 0) && (ind[0] <= MOA_in->shape[0])) 
	    && !((ind[0] < 0) && (abs(ind[0]) <= MOA_in->shape[0] + 1 ) ) ) {
	mprintf(1, "Take Error: Invalid index!\n", threadnum);
	return -1;
      }
    }
  }
  else {
    valid_index = (long *) mmalloc (MOA_in->dimn * sizeof(long));
    ref_dimn = (long *) mmalloc (MOA_in->dimn * sizeof(long)); 
    Last_ref_dimn =(long *) mmalloc (MOA_in->dimn * sizeof(long));
    rslt->shape = (long *) mmalloc (MOA_in->dimn * sizeof(long));
    for (i = 0; i < MOA_in->dimn; i++) {
      if (ind_ub > i) {
	if ((ind[i] > 0) && (ind[i] <= MOA_in->shape[i])) {
	  if (negative == 1) {
	    mprintf(1, "Take Error: Invalid index!\n", threadnum);
	    return -1;
	  }
	  valid_index[i] = ind[i];
	}
	else if ((ind[i] < 0) && (abs(ind[i]) <= MOA_in->shape[i] + 1)) {
	  if (positive == 1) {
	    mprintf(1, "Take Error: Invalid index!\n", threadnum);
	    return -1;
	  }
	  
	  valid_index[i] = MOA_in->shape[i] - abs(ind[i]);
	}
	else {
	  if (  !((ind[i] >= 0) && (ind[i] <= MOA_in->shape[i])) 
		&& !((ind[i] < 0) && (abs(ind[i]) <= MOA_in->shape[i] + 1 ) ) ) {
	    /*mprintf(1, "Take Error: Invalid index!\n", threadnum);
	    //return -1; */
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
  
  l1_finished = 0;
  counter = 0;
  
  sprintf(msg, "\npositive %ld and valid_index:", positive);
	mprintf (7, msg, threadnum);
  for (i = 0; i < MOA_in->dimn; i++) {
    sprintf(msg, " %ld ", valid_index[i]);
		mprintf (7, msg, threadnum);
  }
  mprintf(7, " and input ind:", threadnum);
  for (i = 0; i < MOA_in->dimn; i++) {
    sprintf(msg, " %ld ", ind[i]);
		mprintf (7, msg, threadnum);
  }
  mprintf(7, " and shape:", threadnum);
  for (i = 0; i < MOA_in->dimn; i++) {
    sprintf(msg, " %ld ", MOA_in->shape[i]);
		mprintf (7, msg, threadnum);
  }
  
  
  
  /*call the Psi function with the valid index parameter to return the required parition of the array */
  if (positive == 1) {
    if (MOA_in->dimn == 0) {
      ref_dimn[0] = 0;
      Last_ref_dimn[0] = valid_index[0] - 1;
      rslt->shape[0] = valid_index[0];
    }
    else {
      for (i=0; i< MOA_in->dimn; i++) {
				ref_dimn[i] = 0;
				Last_ref_dimn[i] = valid_index[i] - 1;
				rslt->shape[i] = valid_index[i];
      }
    }
    rslt->elements_ub = Tau(rslt->shape, rslt->dimn);
    rslt->indexes = (unsigned long * ) mcalloc (rslt->elements_ub, sizeof(unsigned long));
    rslt->elements = (MOA_elm * ) mcalloc (rslt->elements_ub, sizeof(MOA_elm));
    if (rslt->elements_ub == 1) 
      rslt->dimn = 0;
    counter = 0;
    l1_finished = 0;
    
   
    mprintf(7, " and  rslt->shape:", threadnum);
    for (i = 0; i < MOA_in->dimn; i++) {
      sprintf (msg, " %ld ", rslt->shape[i]);
			mprintf (7, msg, threadnum);
    }
    offset = Gamma(valid_index, MOA_in->dimn, MOA_in->shape, MOA_in->dimn, 1);
    sprintf(msg, "\n  offset %ld  = & indices are : ", offset);
		mprintf (7, msg, threadnum);
    
    while (l1_finished == 0) {
      rslt->indexes[counter] = Gamma(ref_dimn, MOA_in->dimn, MOA_in->shape, MOA_in->dimn, 1);
      rslt->elements[counter].val = 0;
      rslt->elements[counter].prev = NULL;
      rslt->elements[counter].prev_ub = 0;
      
      sprintf (msg, "  %ld ", rslt->indexes[counter]);
			mprintf (7, msg, threadnum);
      
      counter ++;
      l1_finished = VecIsEqual(ref_dimn, MOA_in->dimn, Last_ref_dimn, MOA_in->dimn);
      NextIndex(valid_index, MOA_in->dimn, ref_dimn);
    }
  }
  
  else {
    if (MOA_in->dimn == 0) {
      if (abs(ind[0]) < MOA_in->shape[0]) {
				ind[0] = abs(ind[0]);
				rslt->shape[0] = ind[0];
				ref_dimn[0] = 0;
				Last_ref_dimn[0] = ind[0] - 1;
      }
      else{
				valid_index[0] = 0;
				ind[0] = MOA_in->shape[0];
				rslt->shape[0] = ind[0];
				ref_dimn[0] = 0;
				Last_ref_dimn[0] = ind[0] - 1;
      }
    }
    else {
      for (i = 0; i < MOA_in->dimn; i++) {
	if (abs(ind[i]) < MOA_in->shape[i]) {
	  ind[i] = abs(ind[i]);
	  rslt->shape[i] = ind[i];
	  ref_dimn[i] = 0;
	  Last_ref_dimn[i] = ind[i] - 1;
	}
	else {
	  valid_index[i] = 0;
	  ind[i] = MOA_in->shape[i];
	  rslt->shape[i] = ind[i];
	  ref_dimn[i] = 0;
	  Last_ref_dimn[i] = ind[i] - 1;
	}
      }
    }
    rslt->elements_ub = Tau(rslt->shape, rslt->dimn);
    rslt->indexes = (unsigned long * ) mcalloc (rslt->elements_ub, sizeof(unsigned long));
    rslt->elements = (MOA_elm * ) mcalloc (rslt->elements_ub, sizeof(MOA_elm));
    if (rslt->elements_ub == 1) 
      rslt->dimn = 0;
   
    mprintf(7, " and  rslt->shape:", threadnum);
    for (i = 0; i < MOA_in->dimn; i++) {
      sprintf (msg, " %ld ", rslt->shape[i]);
			mprintf (7, msg, threadnum);
    }
    
    
    counter = 0;
    l1_finished = 0;
    offset = Gamma(valid_index, MOA_in->dimn, MOA_in->shape, MOA_in->dimn, 1);
    sprintf(msg, "\n  offset %ld  = & indices are : ", offset);
		mprintf (7, msg, threadnum);
    while (l1_finished == 0) {
      rslt->indexes[counter] = offset + Gamma(ref_dimn, MOA_in->dimn, MOA_in->shape, MOA_in->dimn, 1);
      
      sprintf(msg, "  %ld ", rslt->indexes[counter]);
			mprintf (7, msg, threadnum);
      
      counter ++;
      l1_finished = VecIsEqual(ref_dimn, MOA_in->dimn, Last_ref_dimn, MOA_in->dimn);
      NextIndex(ind, MOA_in->dimn, ref_dimn);
    }
  }
  
  if (valid_index != NULL) 	
    free(valid_index);
  if (ref_dimn  != NULL) 
    free(ref_dimn);
  if (Last_ref_dimn != NULL)
    free(Last_ref_dimn);
  return 0;	
}
/* ***************************************************************
// **********************      Drop      *************************
// ***************************************************************/
int Drop (long * ind, long ind_ub, MOA_rec * MOA_in, MOA_rec *  rslt) {
  long * valid_index = NULL;
  int positive, negative, ret;
  long i;
 	char msg[MID_MESSAGE_SIZE];
 
  negative = 0;
  positive = 0;
  
  if (ind[0] >= 0) {
    positive = 1;
  } 
  else {
    negative = 1;
  }
  
  if (MOA_in->dimn == 0) {
    valid_index = (long *) mmalloc (sizeof(long));
    if ((ind[0] >= 0) && ( ind[0] < MOA_in->shape[0])) {
      if (negative == 1) {
	mprintf(1, "Drop Error: Invalid index!\n", threadnum);
	return -1;
      }
      valid_index[0] = MOA_in->shape[0] - ind[0];
    }
    else if ((ind[0] < 0) && (abs(ind[0]) <= MOA_in->shape[0])) {
      if (positive == 1) {
	mprintf(1, "Drop Error: Invalid index!\n", threadnum);
	return -1;
      }
      valid_index[0] = MOA_in->shape[0] - abs(ind[0]);
    }
    else {
      mprintf(1, "Drop Error: Invalid index!\n", threadnum);
      return -1;
    }
  }
  else {
    valid_index = (long *) mmalloc(MOA_in->dimn * sizeof(long));
    for (i = 0; i < MOA_in->dimn; i++) {
      if (ind_ub > i) {
	if ((ind[i] >= 0) && ( ind[i] < MOA_in->shape[i])) {
	  if (negative == 1) {
	    sprintf(msg, "Drop Error: Invalid index at dimension %ld!\n", i);
			mprintf (1, msg, threadnum);
	    return -1;
	  }
	  valid_index[i] = MOA_in->shape[i] - ind[i];
	}
	else if ((ind[i] < 0) && (abs(ind[i]) <= MOA_in->shape[i])) {
	  if (positive == 1) {
	    sprintf(msg, "Drop Error: Invalid index at dimension %ld!\n", i);
			mprintf (1, msg, threadnum);
	    return -1;
	  }
	  valid_index[i] = MOA_in->shape[i] - abs(ind[i]);
	}
	else {
	  sprintf(msg, "Drop Error: Invalid index at dimension %ld!\n", i);
		mprintf (1, msg, threadnum);
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
 
  mprintf(7, "\nin drop before take with valid index =  ", threadnum);
  
  for (i = 0; i < MOA_in->dimn; i ++) {
    sprintf(msg, " %ld ", valid_index[i]);
		mprintf (7, msg, threadnum);
  }
  mprintf(7, " and  MOA_in->shape =  ", threadnum);
  for (i = 0; i < MOA_in->dimn; i ++) {
    sprintf(msg, " %ld ", MOA_in->shape[i]);
		mprintf (7, msg, threadnum);
  }
  mprintf(7, " and  ind  =  ", threadnum);
  for (i = 0; i < MOA_in->dimn; i ++) {
    sprintf(msg, " %ld ", ind[i]);
		mprintf (7, msg, threadnum);
  }
  
  ret = Take(valid_index, MOA_in->dimn, MOA_in, rslt);
  
  if (valid_index != NULL) 	
    free (valid_index);
  
  return ret;
}

/* ***************************************************************
// *******  Drop without creating the whole Tensor  **************
// ***************************************************************/
int DropInd (long * ind, long ind_ub, MOA_rec * MOA_in, MOA_rec *  rslt) {
  long * valid_index = NULL;
  int positive, negative, ret;
  long i;
	char msg[MID_MESSAGE_SIZE];
  
  negative = 0;
  positive = 0;
  
  if (ind[0] >= 0) {
    positive = 1;
  } 
  else {
    negative = 1;
  }
  
  if (MOA_in->dimn == 0) {
    valid_index = (long *) mmalloc (sizeof(long));
    if ((ind[0] >= 0) && ( ind[0] < MOA_in->shape[0])) {
      if (negative == 1) {
	mprintf(1, "Drop Error: Invalid index!\n", threadnum);
	return -1;
      }
      valid_index[0] = MOA_in->shape[0] - ind[0];
    }
    else if ((ind[0] < 0) && (abs(ind[0]) <= MOA_in->shape[0])) {
      if (positive == 1) {
	mprintf(1, "Drop Error: Invalid index!\n", threadnum);
	return -1;
      }
      valid_index[0] = MOA_in->shape[0] - abs(ind[0]);
    }
    else {
      mprintf(1, "Drop Error: Invalid index!\n", threadnum);
      return -1;
    }
  }
  else {
    valid_index = (long *) mmalloc(MOA_in->dimn * sizeof(long));
    for (i = 0; i < MOA_in->dimn; i++) {
      if (ind_ub > i) {
	if ((ind[i] >= 0) && ( ind[i] < MOA_in->shape[i])) {
	  if (negative == 1) {
	    sprintf(msg, "Drop Error: Invalid index at dimension %ld!\n", i);
			mprintf (1, msg, threadnum);
	    return -1;
	  }
	  valid_index[i] = MOA_in->shape[i] - ind[i];
	}
	else if ((ind[i] < 0) && (abs(ind[i]) <= MOA_in->shape[i])) {
	  if (positive == 1) {
	    sprintf(msg, "Drop Error: Invalid index at dimension %ld!\n", i);
			mprintf (1, msg, threadnum);
	    return -1;
	  }
	  valid_index[i] = MOA_in->shape[i] - abs(ind[i]);
	}
	else {
	  sprintf(msg, "Drop Error: Invalid index at dimension %ld!\n", i);
		mprintf (1, msg, threadnum);
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
 
  mprintf(7, "\nin drop before take with valid index =  ", threadnum);
    
  for (i = 0; i < MOA_in->dimn; i ++) {
    sprintf(msg, " %ld ", valid_index[i]);
		mprintf (7, msg, threadnum);
  }
  mprintf(7, " and  MOA_in->shape =  ", threadnum);
  for (i = 0; i < MOA_in->dimn; i ++) {
    sprintf(msg, " %ld ", MOA_in->shape[i]);
		mprintf (7, msg, threadnum);
  }
  mprintf(7, " and  ind  =  ", threadnum);
  for (i = 0; i < MOA_in->dimn; i ++) {
    sprintf(msg, " %ld ", ind[i]);
		mprintf (7, msg, threadnum);
  }
  
  ret = TakeInd(valid_index, MOA_in->dimn, MOA_in, rslt);
  
  if (valid_index != NULL) 	
    free (valid_index);
  
  return ret;
}


/* ***************************************************************
// **********************     Reshape    *************************
// ***************************************************************/

void Reshape (long * new_shape, long  n_shape_ub, MOA_rec * MOA_in, MOA_rec * rslt)
{

  long num_elements, i;
  div_t div_result;
  int k;
  
  num_elements = Tau(new_shape, n_shape_ub);
  rslt->dimn = n_shape_ub;
  
	
  if (rslt->dimn ==0) {
    rslt->shape = (long * ) mmalloc (sizeof(long));
    rslt->shape[0] = new_shape[0];
  }
  else {
    rslt->shape = (long *) mmalloc (sizeof(long) * rslt->dimn);
    
    for (k = 0; k < rslt->dimn; k++) {
      rslt->shape[k] = new_shape[k];
    }
  }
  /* rslt->lelm_size = MOA_in->lelm_size; */
  rslt->elements_ub = num_elements;
  rslt->elements = (MOA_elm *) mmalloc (sizeof(MOA_elm) * rslt->elements_ub);
  for (i = 0; i < rslt->elements_ub ; i++) {
    /* i Mod orig_num_elements */
    div_result = div(i, MOA_in->elements_ub);
    rslt->elements[i].val = MOA_in->elements[div_result.rem].val; /* i Mod orig_num_elements */
  }
      
}
/* ***************************************************************
// **********************    Catenate    *************************
// ***************************************************************/



void Catenate (MOA_rec * MOA_1, MOA_rec * MOA_2, long  Cat_DIM, MOA_rec * rslt)
{
  
  long * ref_dim = NULL;
  long * Last_ref_dimn = NULL;
  long current_idx, ref_idx, i;
  int l1_finished;
  char msg[MID_MESSAGE_SIZE];

  if (MOA_1->dimn != MOA_2->dimn) {
    mprintf(1, "Catenate Error: The Dimension of the Input Arrays must be equal!\n", threadnum);
    return;
  }
  if ((Cat_DIM < 0) || (Cat_DIM > MOA_1->dimn)) {
    sprintf(msg, "Catenation Dimension Can't be less than nor greater than the Dimension of the Input Arrays which is: %ld \n", MOA_1->dimn);
		mprintf (1, msg, threadnum);
    return;
  }
  
  
  rslt->dimn = MOA_1->dimn;
  if (MOA_1->dimn == 0) {
    ref_dim = (long *) mmalloc (sizeof(long));
    Last_ref_dimn = (long *) mmalloc (sizeof(long));
    rslt->shape = (long *) mmalloc (sizeof(long));
    rslt->shape[0] = MOA_1->shape[0];
  }
  else {
    ref_dim =(long *) mmalloc (sizeof(long) * MOA_1->dimn);
    Last_ref_dimn = (long *) mmalloc (sizeof(long) * MOA_1->dimn);
    rslt->shape = (long *) mmalloc (sizeof(long) * MOA_1->dimn);
    for (i = 0; i< MOA_1->dimn; i++) {
      rslt->shape[i] = MOA_1->shape[i];
    }
  }
  
  if (Cat_DIM == 0) 
    Cat_DIM = 1;
  
  rslt->shape[Cat_DIM - 1] = MOA_2->shape[Cat_DIM - 1] + MOA_1->shape[Cat_DIM - 1];
  
  if (MOA_1->dimn == 0) {
    ref_dim[0] = 0;
    Last_ref_dimn[0] = rslt->shape[0] - 1;
  }
  else {
    for (i = 0; i < MOA_1->dimn; i++ ) {
      ref_dim[i] = 0;
      Last_ref_dimn[i] = rslt->shape[i] - 1;
    }
  }
  
  rslt->elements_ub = Tau(rslt->shape, rslt->dimn);
  rslt->elements = (MOA_elm *) mmalloc (sizeof(MOA_elm) * rslt->elements_ub);
  l1_finished = 0;
  while (!l1_finished == 0) {
    current_idx = Gamma(ref_dim, MOA_1->dimn, rslt->shape, rslt->dimn, 1);
    if (ref_dim[Cat_DIM - 1] >= MOA_1->shape[Cat_DIM - 1]){
      ref_dim[Cat_DIM - 1] = ref_dim[Cat_DIM - 1] - MOA_1->shape[Cat_DIM - 1];
      ref_idx = Gamma(ref_dim, MOA_1->dimn, MOA_2->shape, MOA_2->dimn, 1);
      rslt->elements[current_idx] = MOA_2->elements[ref_idx];
      ref_dim[Cat_DIM - 1] = ref_dim[Cat_DIM - 1] + MOA_1->shape[Cat_DIM - 1];
    } 
    else {
      ref_idx = Gamma(ref_dim, MOA_1->dimn, MOA_1->shape, MOA_1->dimn, 1);
      rslt->elements[current_idx] = MOA_1->elements[ref_idx];
    }
    l1_finished = VecIsEqual(ref_dim, rslt->dimn, Last_ref_dimn, rslt->dimn);
    NextIndex(rslt->shape, rslt->dimn, ref_dim);
  }
  
  if (ref_dim != NULL)	
    free(ref_dim);
  if (Last_ref_dimn != NULL)	
    free(Last_ref_dimn);
  
}

/* ***************************************************************
// **********************    scalar_op    ************************
// ***************************************************************/

void scalar_op(char Op, MOA_rec * MOA_in, long scalar, MOA_rec * rslt)
{
  long i;

  rslt->dimn = MOA_in->dimn;
  if (rslt->dimn == 0) {
    rslt->shape = (long *) mmalloc (sizeof(long));
    rslt->shape[0] = MOA_in->shape[0];
  }
  else {
    rslt->shape = (long *) mmalloc (sizeof(long) * rslt->dimn);
    for (i = 0; i< rslt->dimn; i++) {
      rslt->shape[i] = MOA_in->shape[i];
    }
  }

  rslt->elements_ub = Tau(rslt->shape, rslt->dimn);
  rslt->elements = (MOA_elm *) mmalloc (sizeof(MOA_elm) * rslt->elements_ub);
  
  for (i = 0; i < rslt->elements_ub; i++) {
    switch (Op) {
    case '+':
      rslt->elements[i].val = MOA_in->elements[i].val + scalar;
      break;
    case '*':
      rslt->elements[i].val = MOA_in->elements[i].val * scalar;
      break;
    case '-':
      rslt->elements[i].val = MOA_in->elements[i].val - scalar;
      break;
    case '/':
      if (scalar == 0) 
	rslt->elements[i].val = MOA_in->elements[i].val;
      else
	rslt->elements[i].val = MOA_in->elements[i].val / scalar;
      
      break;
    case '<':
      if (MOA_in->elements[i].val < scalar) {
	rslt->elements[i].val = MOA_in->elements[i].val;
      }
      else {
	rslt->elements[i].val = scalar;
      }
      break;
    case '>':
      if (MOA_in->elements[i].val > scalar) {
	rslt->elements[i].val = MOA_in->elements[i].val;
      }
      else {
	rslt->elements[i].val = scalar;
      }
      break;
    default:
      rslt->elements[i].val = MOA_in->elements[i].val;
    }
  }
        
}

long getMaxOnLastBorder(MOA_rec * MOA_in, long * flatindex) {
  long i, j, maxValue, maxIndex;
  long * m_index = NULL; /* multidimensional index */
  
 
  maxIndex = -1;
  m_index = (long *)  mcalloc ((MOA_in->dimn), sizeof(long));
  
  maxValue = MOA_in->elements[0].val;
  for (i = 1; i< MOA_in->elements_ub; i++)  {
    maxValue = MOA_in->elements[i].val;
    Gamma_Inverse(i, MOA_in->shape, MOA_in->dimn, m_index);
    for (j = 0; j< MOA_in->dimn; j++) {
      if (m_index[j] == MOA_in->shape[j] - 1) {
	if (MOA_in->elements[i].val >= maxValue) {
	  maxValue = MOA_in->elements[i].val;
	  maxIndex = i;
	}
      }
    }
  }

  if (m_index != NULL)	
    free(m_index);
  (*flatindex) = maxIndex;
  return maxValue;
}

long MOA_max(MOA_rec *  MOA_in, int callFlag, long ubMaxVal, long * flatindex) {
  long i, maxValue, maxIndex;

  
  maxValue = MOA_in->elements[0].val;
  maxIndex = 0;
    
  for (i = 0; i< MOA_in->elements_ub; i++)  {
    if ((callFlag == 1) && (MOA_in->elements[i].val >= maxValue) && (MOA_in->elements[i].val < ubMaxVal)) {
      maxValue = MOA_in->elements[i].val;
      maxIndex = i;
    }
    else if ((callFlag == 0) && (MOA_in->elements[i].val >= maxValue)) {
      maxValue = MOA_in->elements[i].val;
      maxIndex = i;
    }
  }

  (*flatindex) = maxIndex;
  return maxValue;
}


/* ***************************************************************
// **********************    MOAGetLowerNeighbors ****************
// ***************************************************************/
int MOAGetLowerNeighbors (long * ind, MOA_rec * MOA_in, MOA_rec * rslt) {
  long * tmpInd1 = NULL;
  long * tmpInd2 = NULL;
  long * tmpInd3 = NULL;
  long i, idx;
  int l1_finished, ret;
  char msg[MID_MESSAGE_SIZE];

	sprintf(msg, "in MOAGetLowerNeighbors dimn %ld\n", MOA_in->dimn);  
	mprintf(10, msg, threadnum);
  MOA_rec * rslt2;
  tmpInd1 =  (long *) mmalloc (sizeof(long) * MOA_in->dimn);
  tmpInd2 =  (long *) mmalloc (sizeof(long) * MOA_in->dimn);
  tmpInd3 =  (long *) mmalloc (sizeof(long) * MOA_in->dimn);
  /*rslt2 = (MOA_rec *) mmalloc(sizeof(MOA_rec));*/
  createMOAStruct (&rslt2);
	sprintf(msg, "ind0 %ld 1 %ld 2 %ld\n", ind[0], ind[1], ind[2]);
	mprintf(10, msg, threadnum);
  for (i = 0; i < MOA_in->dimn; i++) {
    tmpInd1[i] = 2;
    if (ind[i] > 0)
      tmpInd2[i] = ind[i] - 1;
    else
      return -1;
  }
 
  ret = Drop(tmpInd2, MOA_in->dimn, MOA_in, rslt2);  
  if (ret >= 0)
    ret = Take(tmpInd1, rslt2->dimn, rslt2, rslt);

  if (ret < 0) 
      return -1;
  
  for (i = 0; i < MOA_in->dimn; i++) {
    tmpInd1[i] = 0;
    tmpInd2[i] = 1;
    tmpInd3[i] = 2;
  }
  
    l1_finished = 0;
    i = 0;
    if (rslt->elements_ub > 0) {
      while (l1_finished == 0) {
	idx = Gamma(tmpInd1, rslt2->dimn, rslt2->shape, rslt2->dimn, 1);
	rslt->indexes[i] =  rslt2->indexes[idx];
	i++;
	l1_finished = VecIsEqual(tmpInd1, MOA_in->dimn, tmpInd2, MOA_in->dimn);
	NextIndex(tmpInd3, MOA_in->dimn, tmpInd1);
      }
    }

  if (tmpInd1 != NULL)
    free(tmpInd1);
  if (tmpInd2 != NULL)
    free(tmpInd2);
  if (tmpInd3 != NULL)
    free(tmpInd3);
  deleteMOA (rslt2);  
  return i-1;
}


/* ***************************************************************
// **********************    MOAGetHigherNeighbors ****************
// ***************************************************************/
int MOAGetHigherNeighbors (long stride, long * ind, MOA_rec * MOA_in, MOA_rec * rslt) {
  long * tmpInd1 = NULL;
  long * tmpInd2 = NULL;
  long * tmpInd3 = NULL;
  long i, idx, ret;
  int l1_finished;

  MOA_rec * rslt2;
  tmpInd1 =  (long *) mmalloc (sizeof(long) * MOA_in->dimn);
  tmpInd2 =  (long *) mmalloc (sizeof(long) * MOA_in->dimn);
  tmpInd3 =  (long *) mmalloc (sizeof(long) * MOA_in->dimn);
  /*rslt2 = (MOA_rec *) mmalloc (sizeof(MOA_rec));*/
  createMOAStruct (&rslt2);
  for (i = 0; i < MOA_in->dimn; i++) {

    if (ind[i] > MOA_in->shape[i] - 1)
      return -1;
    else {
      tmpInd1[i] = stride;
      tmpInd2[i] = ind[i];
    }
  }
 

  ret = DropInd(tmpInd2, MOA_in->dimn, MOA_in, rslt2);  
  if (ret >= 0)
    ret = TakeInd(tmpInd1, rslt2->dimn, rslt2, rslt);

  if (ret < 0) 
      return -1;
  
  for (i = 0; i < MOA_in->dimn; i++) {
    tmpInd1[i] = 0;
    tmpInd3[i] = rslt->shape[i];
    tmpInd2[i] = tmpInd3[i] - 1;
  }
  
  l1_finished = 0;
  i = 0;
  if (rslt->elements_ub > 0) {
    while (l1_finished == 0) {
      idx = Gamma(tmpInd1, rslt2->dimn, rslt2->shape, rslt2->dimn, 1);
      rslt->indexes[i] =  rslt2->indexes[idx];
      i++;
      l1_finished = VecIsEqual(tmpInd1, MOA_in->dimn, tmpInd2, MOA_in->dimn);
      NextIndex(tmpInd3, MOA_in->dimn, tmpInd1);
    }
  }
  
  if (tmpInd1 != NULL)	
    free(tmpInd1);
  if (tmpInd2 != NULL)	
    free(tmpInd2);
  if (tmpInd3 != NULL)	
    free(tmpInd3);
  deleteMOA (rslt2);  
  return i-1;
}


/* ***************************************************************
// **********************    getDiagonals ************************
// ***************************************************************/
void getDiagonals (long waveNo, long waveSize, long * ind, MOA_rec * MOA_in) {
  MOA_rec * MOA_nghb;
  long i, j, k, l, border, flatIndex; 
  int moreneighb;
  long * nghb_ind;
  long combnum; /* number of combinations for decrementing i number of indices*/
  long * * combin; /* combinations of lower neighbor indices*/
	char msg[MID_MESSAGE_SIZE];

  createMOAStruct (&MOA_nghb);

  moreneighb = MOAGetHigherNeighbors (waveSize, ind, MOA_in, MOA_nghb);
  if (moreneighb <= 0)
    return;
  else {
    sprintf(msg, "\nWave test 5  %ld moreneighb %d MOA_nghb->elements_ub %ld from point : ", waveNo , moreneighb, MOA_nghb->elements_ub);
		mprintf (6, msg, threadnum);
    for (j = 0; j <  MOA_in->dimn; j++) {
      sprintf(msg, "%ld ", ind[j]);
			mprintf (6, msg, threadnum);
     
    }
    sprintf(msg, " has %d points with Indices: ", moreneighb);
		mprintf (6, msg, threadnum);
    for (i=0;i<MOA_nghb->elements_ub; i++) {
      sprintf(msg, "%ld ", MOA_nghb->indexes[i]);    
			mprintf (6, msg, threadnum);
      /*assign point to the next available processor*/
    }
    
    
    nghb_ind =  (long *) mmalloc (sizeof(long) * MOA_nghb->dimn);
      
    for (j = 1; j< MOA_nghb->dimn; j++)  {
      
      /* number of the different combinations of which indeces in the multidimensional index to give the shape bound so that all edge points can be processed */
      combnum = (Factorial(MOA_nghb->dimn)/(Factorial(MOA_nghb->dimn-j) * Factorial(j)));
      /* create memory for combinations matrix */
      combin = (long * *) mcalloc (combnum, sizeof(long *));
      for (k=0; k<combnum; k++)
	combin[k] = (long *) mcalloc (j, sizeof(long));
      if (combin == NULL) {
				mprintf (1, "getDiagonals Error: Can not allocate memory to combinations matrix\n", threadnum);
				return;
      }
      /* get matrix of all possible combinations */
      Combinations(MOA_nghb->dimn, j, &combin);
      /* loop to decrement the selected indices */
      for (k=0; k<combnum; k++) {
	for (l = 0; l <MOA_nghb->dimn; l++) {
	  nghb_ind[l] = 0;
	}
	for (l=0; l<j; l++) {
	  nghb_ind[combin[k][l]-1] = MOA_nghb->shape[combin[k][l]-1] - 1;
	}
	flatIndex = Gamma( nghb_ind, MOA_nghb->dimn,  MOA_nghb->shape, MOA_nghb->dimn, 1);
	
	Gamma_Inverse(MOA_nghb->indexes[flatIndex], MOA_in->shape, MOA_in->dimn, ind);
	
	border = 1;
	for (l=0; l<MOA_in->dimn; l++) {
	  if (ind[l] >= MOA_in->shape[l] - 1) {
	    border = 0;
	  }
	}
	if (border == 1) {
	  getDiagonals (waveNo + 1, waveSize, ind, MOA_in);
	}
      }
    }
    /* do the final edge where all at the shape limit*/
    for (l = 0; l <MOA_nghb->dimn; l++) {
      nghb_ind[l] = MOA_nghb->shape[l] - 1;
    }
    flatIndex = Gamma( nghb_ind, MOA_nghb->dimn,  MOA_nghb->shape, MOA_nghb->dimn, 1);
    
    Gamma_Inverse(MOA_nghb->indexes[flatIndex], MOA_in->shape, MOA_in->dimn, ind);
    
    border = 1;
    for (l=0; l<MOA_in->dimn; l++) {
      if (ind[l] >= MOA_in->shape[l] - 1) {
	border = 0;
      }
    }
    if (border == 1) {
      getDiagonals (waveNo + 1, waveSize, ind, MOA_in);
    }
    free(nghb_ind);
  }
  deleteMOA (MOA_nghb);  
}


/* ***************************************************************
// **********************    Navigate ****************************
// ***************************************************************/
void GetNextBorderCell (long * * ind,  MOA_rec * MOA_in) {
  
}

void Navigate (long * NavDir, long stride, long * stIndex, MOA_rec * MOA_in, long * * rslt, long * * rsltInd, long * rslt_ub) {
  long * ind = NULL;
  long i, count, rounds, error;
  int border, done;
  border = done = rounds = 0;

  ind =  (long *) mmalloc (sizeof(long) * MOA_in->dimn);
  count = 0;

  (*rslt) =  (long *) mmalloc (sizeof(long));
  (*rsltInd) =  (long *) mmalloc (sizeof(long));
  (*rsltInd)[count] = Gamma(stIndex, MOA_in->dimn, MOA_in->shape, MOA_in->dimn, 1);
  (*rslt)[count] = MOA_in->elements[(*rsltInd)[count]].val;
  for (i=0; i<MOA_in->dimn; i++) {
    ind[i] = stIndex[i] + NavDir[i];
    if ((ind[i] < 0) || (ind[i] >= MOA_in->shape[i]))
      border = 1;
  }
  while (done == 0) {
    while (border == 0) {
      count ++;    
      (*rslt) =  (long *) realloc ((*rslt), sizeof(long) * (count+1));
     if ((*rslt) == NULL) {
       mprintf(1, "Navigate Error: Could not reallocate memory!\n", threadnum);
       return;
     }
     (*rsltInd) =  (long *) realloc ((*rsltInd), sizeof(long) * (count+1));
     if ((*rsltInd) == NULL) {
       mprintf(1, "Navigate Error:Could not reallocate memory!\n", threadnum);
       return;
     }
     
     (*rsltInd)[count] = Gamma(ind, MOA_in->dimn, MOA_in->shape, MOA_in->dimn, 1);
     (*rslt)[count] = MOA_in->elements[(*rsltInd)[count]].val;
     
     for (i=0; i<MOA_in->dimn; i++) {
       ind[i] = ind[i] + NavDir[i];
       if ((ind[i] < 0) || (ind[i] >= MOA_in->shape[i]))
	 border = 1;
     }
    }
    border = 0;
    rounds ++;
    GetNextBorderCell(&ind, MOA_in);
    error = 0;
    for (i=0; i<MOA_in->dimn; i++) {
      if (ind[i] >= MOA_in->shape[i]) {
	done = 1;
	error = 1;
      }
      if (ind[i] >= MOA_in->shape[i] - 1)
	done = 1;
      NavDir[i] = NavDir[i] * (-1);
    }
    if (rounds == stride)  
	done = 1;
    if (error == 0) {
      count ++;    
      (*rslt) =  (long *) realloc ((*rslt), sizeof(long) * (count+1));
     if ((*rslt) == NULL) {
			mprintf(1, "Navigate Error: Could not reallocate memory!\n", threadnum);
			return;
      }
      (*rsltInd) =  (long *) realloc ((*rsltInd), sizeof(long) * (count+1));
     if ((*rsltInd) == NULL) {
			mprintf(1, "Navigate Error: Could not reallocate memory!\n", threadnum);
			return;
      }
      (*rsltInd)[count] = Gamma(ind, MOA_in->dimn, MOA_in->shape, MOA_in->dimn, 1);
      (*rslt)[count] = MOA_in->elements[(*rsltInd)[count]].val;
    }
  }
  (*rslt_ub) = count + 1;
  if (ind != NULL)	
    free(ind);
}
