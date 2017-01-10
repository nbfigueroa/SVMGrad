/*
 * Copyright (C) 2012 Learning Algorithms and Systems Laboratory, EPFL, Switzerland
 *
 *  util.h
 *
 *  Created on : Aug 14, 2012
 *  Author     : Ashwini Shukla and Saurav Aryan
 *  Email      : ashwini.shukla@epfl.ch
 *  Website    : lasa.epfl.ch
 *
 * Permission is granted to copy, distribute, and/or modify this program
 * under the terms of the GNU General Public License, version 2 or any
 * later version published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details
 */

#ifndef _util_H_
#define _util_H_

#include <math.h>
#include <string.h>
#include <iostream>

#ifndef KILL
#define KILL(a) {if(a!=0) {delete [] a; a = 0;}}
#endif


double getkernel(double *x1, double *x2, double lambda, const char* type, int n);
bool getfirstkernelderivative(double *x1, double *x2, double lambda, const char* type, int der_wrt, double* der_val, int n);
bool getsecondkernelderivative(double *x1, double *x2, int n, double lambda, const char *type, double **hesval);

double getkernel(double *x1, double *x2, double* lambda, const char* type, int n);
bool getfirstkernelderivative(double *x1, double *x2, double* lambda, const char* type, int der_wrt, double* der_val, int n);
bool getsecondkernelderivative(double *x1, double *x2, int n, double* lambda, const char *type, double **hesval);

void VectorMatrixMultipy(double *VectorA, double **MatrixB, double *Result, int n, int p);
void MatrixVectorMultipy(double **MatrixB, double *VectorA, double *Result, int cols, int rows);
double arraydot(double *x, double *y, int m);
double norm(double *x, int m);
double norm2(double *x, int m);
#endif
