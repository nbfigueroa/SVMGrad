/*
 * Copyright (C) 2012 Learning Algorithms and Systems Laboratory, EPFL, Switzerland
 *
 *  util.cpp
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

#include "util.h"
using namespace std;

void VectorMatrixMultipy(double *VectorA, double **MatrixB, double *Result, int cols, int rows)	//m=1
{
	for(int i=0; i<cols; i++)
	{
		Result[i] = 0;
		for(int j=0; j<rows; j++)
		{
			Result[i] += VectorA[j]*MatrixB[j][i];
		}
	}
}

void MatrixVectorMultipy(double **MatrixB, double *VectorA, double *Result, int cols, int rows)	//m=1
{
	for(int i=0; i<rows; i++)
	{
		Result[i] = 0;
		for(int j=0; j<cols; j++)
		{
			Result[i] += VectorA[j]*MatrixB[i][j];
		}
	}
}

double arraydot(double *x, double *y, int m)
{
	double dot = 0.0;
	for(int i =0; i<m; i++)
		dot += x[i]*y[i];
	return dot;
}

double norm(double *x, int m)
{
	double sum = 0.0;
	for(int i=0; i<m; i++)
		sum += x[i]*x[i];

	return sqrt(sum);
}

double norm2(double *x, int m)
{
	double sum = 0.0;
	for(int i=0; i<m; i++)
		sum += x[i]*x[i];

	return sum;
}

double getkernel(double *x1, double *x2, double lambda, const char* type, int n)
{
	double ker_val = 0.0;
	double *diff = new double[n];

	for(int i =0; i<n; i++)
	{
		diff[i] = x1[i] - x2[i];
	}
	if(strcmp(type, "poly") == 0)
		ker_val = pow((arraydot(x1, x2, n)+1),lambda);
	else if(strcmp(type, "rbf") == 0)
	{
		ker_val = exp(-lambda*norm2(diff, n));
	}
	else
		cout<<"\nInvalid kernel type specified in getkernel function!";
	delete diff;
	return ker_val;
}
double getkernel(double *x1, double *x2, double* lambda, const char* type, int n)
{
	double ker_val = 0.0;
	double diff=0;

	for(int i =0; i<n; i++)
	{
		diff -= lambda[i]*pow(x1[i] - x2[i],2);
	}
	if(strcmp(type, "rbf") == 0)
	{
		ker_val = exp(diff);
	}
	else
		cout<<"\nInvalid kernel type specified in getkernel function!";

	return ker_val;
}

bool getfirstkernelderivative(double *x1, double *x2, double* lambda, const char* type, int der_wrt, double* der_val, int n)
{
	double temp = 0.0;

	if(strcmp(type, "rbf") == 0)
	{
		double diff=0;
		for(int i =0; i<n; i++)
		{
			der_val[i] = 2*lambda[i]*(x1[i] - x2[i]);
			diff -= lambda[i]*pow(x1[i] - x2[i],2);
		}

		if(der_wrt == 1)
		{
			for(int i =0; i<n; i++)
				der_val[i] *= -exp(diff);
		}
		else
		{
			for(int i =0; i<n; i++)
				der_val[i] *= exp(diff);
		}


		return true;
	}
	else
	{
		cout<<"\nInvalid kernel type specified in getkernel function!";
		return false;
	}
}

bool getfirstkernelderivative(double *x1, double *x2, double lambda, const char* type, int der_wrt, double* der_val, int n)
{
	double temp = 0.0;

	if(strcmp(type, "poly") == 0)
	{
		temp = lambda*pow((arraydot(x1, x2, n)+1),(lambda-1));
		if(der_wrt == 1)
		{
			for(int i=0; i<n; i++)
				der_val[i] = temp*x2[i];
		}
		else
		{
			for(int i=0; i<n; i++)
				der_val[i] = temp*x1[i];
		}

		return true;
	}
	else if(strcmp(type, "rbf") == 0)
	{
		double *diff;
		diff = new double[n];
		for(int i =0; i<n; i++)
			diff[i] = x1[i] - x2[i];

		if(der_wrt == 1)
		{
			temp = -2*lambda*exp(-lambda*norm2(diff, n));
		}
		else
		{
			temp = 2*lambda*exp(-lambda*norm2(diff, n));
		}

		for(int i =0; i<n; i++)
			der_val[i] = temp*diff[i];
		delete diff;
		return true;
	}
	else
	{
		cout<<"\nInvalid kernel type specified in getkernel function!";
		return false;
	}
}


bool getsecondkernelderivative(double *x1, double *x2, int n, double* lambda, const char *type, double **hesval)
{
	if(strcmp(type, "rbf") == 0)
	{
		double *tmp1 = new double[n];
		double diff=0;
		for(int i =0; i<n; i++)
		{
			tmp1[i] = lambda[i]*(x1[i] - x2[i]);
			diff -=lambda[i]*pow(x1[i]-x2[i],2);
		}

		double diff2=2*exp(diff);

		for(int i = 0; i<n; i++)
			for(int j =0; j<n; j++)
				if(i == j)
					hesval[i][j] = diff2*(lambda[i] - 2*tmp1[i]*tmp1[j]);
				else
					hesval[i][j] = diff2*(- 2*tmp1[i]*tmp1[j]);
		delete tmp1;
		return true;
	}
	else
	{
		cout<<"\nInvalid type specified in the getsecondkernelderivative";
		return false;
	}
}


bool getsecondkernelderivative(double *x1, double *x2, int n, double lambda, const char *type, double **hesval)
{
	if(strcmp(type, "poly") == 0)
	{
		double tmp = arraydot(x1, x2, n) + 1;
		for(int i =0; i< n; i++)
			for(int j=0; j< n; j++)
			{
				if(i == j)
					hesval[i][j] = lambda* pow(tmp,(lambda-2))*( tmp + (lambda-1)*x2[i]*x1[j]);
				else
					hesval[i][j] = lambda* pow(tmp,(lambda-2))*((lambda-1)*x2[i]*x1[j]);
			}
		return true;
	}
	else if(strcmp(type, "rbf") == 0)
	{
		double *tmp1 = new double[n];
		for(int i =0; i<n; i++)
			tmp1[i] = x1[i] - x2[i];

		for(int i = 0; i<n; i++)
			for(int j =0; j<n; j++)
				if(i == j)
					hesval[i][j] = 2*lambda*exp(-lambda*norm2(tmp1,n)) * ( 1 - 2*lambda*(tmp1[i]*tmp1[j]));
				else
					hesval[i][j] = 2*lambda*exp(-lambda*norm2(tmp1, n)) * (- 2*lambda*(tmp1[i]*tmp1[j]));
		delete tmp1;
		return true;
	}
	else
	{
		cout<<"\nInvalid type specified in the getsecondkernelderivative";
		return false;
	}
}


