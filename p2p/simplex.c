/*
 * Copyright (c) 2000 Matteo Frigo
 * Copyright (c) 2000 Biomedin s. r. l.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* $Id: simplex.c,v 1.2 2001/05/08 18:20:15 matley Exp $ */
/* interface between C code and the LP solver */

#include "config.h"
#include "simplex.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "f2c.h"

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

static void panic(const char *format,...)
{
     va_list ap;

     fprintf(stderr, "fatal error: ");
     va_start(ap, format);	/* Variable argument begin */
     vfprintf(stderr, format, ap);
     va_end(ap);		/* Variable argument end */

     exit(1);
}

static void *xmalloc(size_t size)
{
     void *x = malloc(size);
     if (!x)
	  panic("Out of memory\n");

     return x;
}

/*
 * d_mod from libf2c
 */
double d_mod(double *x, double *y);
double d_mod(double *x, double *y)
{
#ifdef HAVE_DREM
        double xa, ya, z;
        if ((ya = *y) < 0.)
                ya = -ya;
        z = drem(xa = *x, ya);
        if (xa > 0) {
                if (z < 0)
                        z += ya;
                }
        else if (z > 0)
                z -= ya;
        return z;
#else
        double quotient;
	double result;
        if( (quotient = (*x) / (*y)) >= 0)
                quotient = floor(quotient);
        else
                quotient = -floor(-quotient);
	result = (*x) - ((*y) * quotient);
        return result;
#endif
}

/*************************************************************
 * SIMPLEX manipulation routines
 *************************************************************/
/* allocate a SIMPLEX */
SIMPLEX *simplex_make(int maxnrow, int ncol)
{
     SIMPLEX *self = (SIMPLEX *)xmalloc(sizeof(SIMPLEX));
     self->nrow = 0;
     self->maxnrow = maxnrow;
     self->iteration_limit = 10000;

     self->A = (double *)xmalloc(maxnrow * ncol * sizeof(double));
     self->ncol = ncol;
     self->rhs = (double *)xmalloc(maxnrow * sizeof(double));
     self->obj = (double *)xmalloc(ncol * sizeof(double));
     self->duals = (double *)xmalloc(maxnrow * sizeof(double));
     self->primals = (double *)xmalloc(ncol * sizeof(double));
     return self;
}

/* destroy a SIMPLEX */
void simplex_free(SIMPLEX *self)
{
     free(self->A);
     free(self->rhs);
     free(self->obj);
     free(self->duals);
     free(self->primals);
     free(self);
}

/* add an equation to a SIMPLEX. */
void simplex_add_equation(SIMPLEX *self, double *row, double rhs)
{
     int i, j;
     int N = self->ncol;
     int M = self->maxnrow;
     
     if (self->nrow == self->maxnrow)
	  panic("Too many equations\n");

     i = self->nrow;

     for (j = 0; j < N; j++) {
	  self->A[i + M * j] = row[j];
     }

     self->rhs[i] = rhs;
     self->nrow++;
}

/* set the SIMPLEX objective function */
void simplex_set_objective(SIMPLEX *self, double *row)
{
     int j;
     for (j = 0; j < self->ncol; ++j)
	  self->obj[j] = row[j];
}

void simplex_set_iteration_limit(SIMPLEX *self, int l)
{
     self->iteration_limit = l;
}

enum simplex_result simplex_solve(SIMPLEX *self)
{
     int M = self->nrow;
     int N = self->ncol;
     integer INPUT[19];
     double TOL[4];
     double OBJ;
     double *X = (doublereal *)xmalloc((M + 2) * sizeof(double));
     integer *JX = (integer *)xmalloc(M * sizeof(integer));
     /*
      * The size of DJ must be at least min(N, M).  The LPPRIM
      * manual is wrong.  DJ is used in LPPVER, aliased as
      * IS, which is a vector of size N. 
      */
     double *DJ = (double *)xmalloc((N + M) * sizeof(double));
     double ERR[4];
     integer IOUT[4];
     double *WS = (double *)xmalloc((M * M + 5 * M + 6) * sizeof(double));
     int i;

     INPUT[0] = self->maxnrow;       /* row dimension of A */
     INPUT[1] = N;                   /* column dimension of A */
     INPUT[2] = M;                   /* number of constraints */
     INPUT[3] = N;                   /* number of structural variables */
     INPUT[4] = 0;                   /* minimize objective */
     INPUT[5] = (M * M + 5 * M + 6); /* size of workspace */
     INPUT[6] = 0;                   /* solve new problem */
     INPUT[7] = 0;                   /* all-equality flag */
     INPUT[8] = 0;                   /* no costs assigned to slack variables */
     INPUT[9] = self->iteration_limit;  /* iteration limit */
     INPUT[10] = -1;                /* equal consecutive objective limit */
     INPUT[11] = -1;                 /* inversion frequency (30) */
     INPUT[12] = 1;                  /* re-invert at optimum */
     INPUT[13] = 0;                  /* no output */
     INPUT[14] = 0;
     INPUT[15] = 0;
     INPUT[16] = 0;
     INPUT[17] = 0;
     INPUT[18] = 72;
     TOL[0] = -1.0;                  /* pivot tolerance */
     TOL[1] = -1.0;                  /* reduced costs tolerance */
     TOL[2] = -1.0;                  /* feasibility tolerance */
     TOL[3] = -1.0;                  /* round-off tolerance */

     lpprim_(self->A, (int *)0, self->rhs, self->obj, INPUT, TOL,
	    &OBJ, X, JX, self->duals, DJ, ERR, IOUT, WS);

     switch (IOUT[0]) {
	 case 1:
	      self->result = G_OPTIMAL;
	      for (i = 0; i < M; ++i)
		   self->duals[i] = -self->duals[i];
	      for (i = 0; i < N; ++i)
		   self->primals[i] = 0.0;
	      for (i = 0; i < M; ++i) {
		   if (JX[i] >= 1 && JX[i] <= N)
			self->primals[JX[i]-1] = X[i];
	      }
	      self->objective_value = OBJ;
	      break;
	 case 2:
	      self->result = G_INFEASIBLE;
	      break;
	 case 3:
	      self->result = G_UNBOUNDED;
	      break;
	 case 4:
	      self->result = G_ITERATION_LIMIT;
	      break;
	 default:
	      panic("Unrecognized return value %d from LPPRIM\n", IOUT[0]);
	      break;
     }

     free(X);
     free(JX);
     free(DJ);
     free(WS);

     return self->result;
}

