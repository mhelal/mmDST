/**********************************************************
* Author: Manal Helal																			*
* Last Modification Date: Fri 12 Jan 2007 03:39:51 AM EST *
* Project : MMSA - Multiple Sequence Alignment Based on 	*
* 	Mathematics of Arrays - PhD Experimentation		*
* File: moa.h, header for the basic MOA functions used in *
* this project																						*
***********************************************************/
#ifndef MOAH_
#define MOAH_

typedef long long MOATypeDimn;
typedef long long MOATypeShape;
typedef long long MOATypeInd;
typedef long long MOATypeElmVal;

typedef struct cellProp_elm { /* Cell Properties */
    int llOC; /* Local Lower Overlapping Cell*/
    int glOC; /* Global Lower Overlapping Cell*/
    int lhOC; /* Local Higher Overlapping Cell*/
    int ghOC; /* Global Higher Overlapping Cell*/
} cellProp;

typedef struct tagMOA_elm {
    cellProp cp;
    MOATypeElmVal val;
    MOATypeDimn prev_ub;
    MOATypeShape * * prev;
} MOA_elm ;

typedef struct tagMOA_rec {
    MOATypeDimn dimn;
    MOATypeShape * shape;
    MOATypeInd elements_ub;
    MOATypeShape * * indexes; /*Multidimensional Global Index*/
    MOA_elm * elements;	
} MOA_rec;


void printMOA_Matrix(int dbglevel, MOATypeDimn dimn, int ident, void * elements, int elm_type, MOATypeInd total_elements, MOATypeDimn first_dimn, char * * sequences);
int printMOA_Sequences(int dbglevel, int ident, MOATypeDimn dimn, MOATypeShape *ind, char * * sequences, MOA_rec * MOA, int elm_type);
void printMOA_dimn(int dbglevel, MOATypeDimn dimn, MOATypeShape * ind, char * * sequences, MOA_rec * MOA, int elm_type);
void printMOA (int dbglevel, MOA_rec * MOA, char * * sequences, int elm_type);
void printMOA1 (MOA_rec * MOA, int elm_type);
void printMOA_scr (MOA_rec * MOA, int elm_type);
void printMOAIndices (MOA_rec * MOA);

void createMOAStruct(MOA_rec * * MOA_val);
int createMOA(MOATypeShape * shape, MOATypeDimn dimn, MOA_rec * MOA_val, int callFlag, int cellValue);
void deleteMOA(MOA_rec * MOA);
MOATypeInd Tau (MOATypeShape * array_in, MOATypeDimn array_ub);
MOATypeInd Sum (MOATypeShape * array_in, MOATypeDimn array_ub);
MOATypeInd Gamma (MOATypeShape * ind, MOATypeDimn ind_ub, MOATypeShape * arr_shape,  MOATypeDimn shape_ub, int Front);
int Gamma_Inverse (MOATypeInd ind, MOATypeShape * arr_shape, MOATypeDimn shape_ub, MOATypeShape * * rslt, int thrd);
void Psi (MOATypeShape * ind, MOATypeInd ind_ub, MOA_rec * MOA_in, MOA_rec * rslt);
int VecIsEqual (MOATypeShape * Array_1, MOATypeDimn array1_ub, MOATypeShape * Array_2, MOATypeDimn array2_ub);
int IsValidIndex (MOATypeShape * shape, MOATypeInd shape_ub, MOATypeShape * index);
int isLowerBorderCell (MOATypeShape * index, MOATypeDimn dimn);
int isHigherBorderCell (MOATypeShape * index, MOATypeDimn dimn, MOATypeShape * shape);
int isHigherBorderCellandNotLower (MOATypeShape * index, MOATypeDimn dimn, MOATypeShape * shape);
void NextIndex (MOATypeShape * shape, MOATypeDimn shape_ub, MOATypeShape * * Prev_Index);

void analyzeCellPosition (MOATypeShape * cellIndex, MOATypeInd flatCellIndex, MOA_rec * * rslt, MOATypeShape * global_shape);
void copyIndicesElm (MOATypeShape * offset, MOA_rec * MOA_in, MOA_rec * * rslt, MOATypeShape * ref_dimn, MOATypeShape * Last_ref_dimn, MOATypeShape * valid_index, int copyElemFlag, int analyzeLBHB);
void copyIndices (MOATypeShape * offset, MOATypeShape * global_shape, MOA_rec * * rslt, MOATypeShape * ref_dimn, MOATypeShape * Last_ref_dimn, MOATypeShape * valid_index, int analyzeLBHB);

int Take (MOATypeShape * ind, MOATypeDimn ind_ub, MOA_rec * MOA_in, MOA_rec * * rslt);
int Drop (MOATypeShape * ind, MOATypeDimn ind_ub, MOA_rec * MOA_in, MOA_rec * * rslt);
int TakeInd (MOATypeShape * ind, MOATypeDimn dimn, MOATypeShape * global_shape, MOA_rec * * rslt, int analyzeLBHB);
int DropInd (MOATypeShape * ind, MOATypeDimn dimn, MOATypeShape * Sshape, MOA_rec * * rslt); 

void Reshape (MOATypeShape * new_shape, MOATypeInd  n_shape_ub, MOA_rec * MOA_in, MOA_rec * rslt);
void Catenate (MOA_rec * MOA_1, MOA_rec * MOA_2, MOATypeDimn  Cat_DIM, MOA_rec * rslt);
void scalar_op(char Op, MOA_rec * MOA_in, long scalar, MOA_rec * rslt);
long getMaxOnLastBorder(MOA_rec * MOA_in, MOATypeInd * flatindex);
long MOA_max(MOA_rec *  MOA_in, int callFlag, MOATypeInd ubMaxVal, MOATypeInd * flatindex);
MOATypeInd MOAGetLocalLowerNeighbors_old (MOATypeShape * ind, MOA_rec * MOA_in, MOA_rec * * rslt);
MOATypeInd MOAGetGlobalLowerNeighbors_old (MOATypeShape * startingIndex, MOATypeDimn dimn, MOATypeShape * shape, MOA_rec * * rslt);
MOATypeInd MOAGetHigherNeighbors_old (long stride, MOATypeShape * startingIndex, MOATypeDimn dimn, MOATypeShape * shape, MOA_rec * * rslt);
MOATypeInd MOAGetHigherPartitions (long pSize, long stride, MOATypeShape * startingIndex, MOATypeDimn dimn, MOATypeShape * shape, MOA_rec * * rslt, int * * validPartitions);
int getLowerNeighbors (long stride, MOATypeShape * startingIndex, MOATypeDimn dimn, MOATypeShape * shape, MOA_rec * * part);
int getHigherNeighbors (long stride, MOATypeShape * startingIndex, MOATypeDimn dimn, MOATypeShape * shape, MOA_rec * * part);
int getNeighbors (long stride, MOATypeShape * startingIndex, MOATypeDimn dimn, MOATypeShape * shape, MOA_rec * * part);
#endif
