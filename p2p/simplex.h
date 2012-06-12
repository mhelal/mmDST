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

/* $Id: simplex.h,v 1.2 2001/05/08 18:20:15 matley Exp $ */

enum simplex_result {
     G_OPTIMAL, G_INFEASIBLE, G_UNBOUNDED, G_ITERATION_LIMIT
};

typedef struct {
     /* public fields */
     enum simplex_result result;
     double *primals;                  /* primal variables */
     double *duals;                    /* dual variables */
     double objective_value;
     int iteration_limit;
     double *obj;                      /* objective function */

     /* private fields */
     int maxnrow;
     int ncol;
     double *A;
     double *rhs;
     int nrow;
} SIMPLEX;

SIMPLEX *simplex_make(int nrow, int ncol);
void simplex_free(SIMPLEX *self);
void simplex_add_equation(SIMPLEX *self, double *row, double rhs);
void simplex_set_objective(SIMPLEX *self, double *row);
enum simplex_result simplex_solve(SIMPLEX *self);
void simplex_set_iteration_limit(SIMPLEX *self, int l);
