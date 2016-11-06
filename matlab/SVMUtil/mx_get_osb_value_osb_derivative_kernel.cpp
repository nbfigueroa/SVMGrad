/*
 * Copyright (C) 2012 Learning Algorithms and Systems Laboratory, EPFL, Switzerland
 * 
 *  mx_get_bsb_value_bsb_derivative_kernel.cpp
 *
 *  Created on : Jan 24, 2013
 *  Author     : Ashwini Shukla
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

#include "../../include/asvm.h"
#include "../SVMWidget/mxUtil.h"
#ifdef MX_API_VER
#if MX_API_VER < 0x07030000
typedef int mwIndex;
#endif
#endif

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{

	if(nrhs == 0)
	{
		printf("Arguments: input_data, target, labels, lambda, ktype.");
		return;
	}
	int M = mxGetN(prhs[0]);
	int N = mxGetM(prhs[0])/2;

	double* lambda;
	if(mxGetNumberOfElements(prhs[3]) == 1)
	{
		lambda = new double[N];
		double tmp = mxGetScalar(prhs[3]);
		for(int i=0;i<N;i++)
			lambda[i] = tmp;
	}
	else if(mxGetNumberOfElements(prhs[3]) == N)
	{
		lambda = (double*)mxGetPr(prhs[3]);
	}
	else
	{
		printf("ERROR: lambda can have length = 1 or equal to the number of dimensions \n");
		return;
	}

	double* labels = (double*)mxGetPr(prhs[2]);
	int P=0;
	for(int i=0;i<M;i++)
		if(labels[i] == 1)
			P++;
	int *pl=new int[P];
	int count=0;
	for(int i=0;i<M;i++) if(labels[i] == 1) pl[count++] = i;

	char type[1025];
	mxGetString(prhs[4],type, 1025 );
	double* input_data = (double*)mxGetPr(prhs[0]);
	double* target = (double*)mxGetPr(prhs[1]);

	double** Q = new double*[M+N+P];
	for(int i=0;i<M+N+P;i++)
		Q[i] = new double[M+N+P];



	double* tmp = new double[N];
	int tmpind, tmpind2;
	double sum;

	//K
	for(int i=0;i<M;i++)
	{
		tmpind=2*i*N;
		for(int j=i;j<M;j++)
		{
			Q[i][j] = labels[i]*labels[j]*getkernel(input_data + tmpind, input_data + 2*j*N, lambda, type,N);
			Q[j][i]=Q[i][j];
		}
	}



	//G
	for(int i=0;i<M;i++)
	{
		tmpind=2*i*N;
		for(int j=0;j<P;j++)
		{
			tmpind2 = 2*pl[j]*N;
			getfirstkernelderivative(input_data+tmpind, input_data + tmpind2, lambda, type, 2, tmp, N);
			Q[i][j+M]=0;
			for(int k=0;k<N;k++)
				Q[i][j+M] += tmp[k]*input_data[tmpind2+N+k];

			Q[i][j+M] *= labels[i];
			Q[j+M][i] = Q[i][j+M];
		}
	}

	//G_s
	for(int i=0;i<M;i++)
	{
		getfirstkernelderivative(input_data+2*i*N, target, lambda, type, 2, tmp, N);
		for(int j=0;j<N;j++)
		{
			Q[i][ j+M+P ] = -labels[i]*tmp[j];
			Q[j+M+P][i] = Q[i][ j+M+P ];
		}
	}

	double** tmpmat = new double*[N];
	for(int i=0;i<N;i++)
		tmpmat[i] = new double[N];

	//H
	for(int i=0;i<P;i++)
	{
		tmpind = 2*pl[i]*N;
		for(int j=i;j<P;j++)
		{
			tmpind2 = 2*pl[j]*N;
			getsecondkernelderivative(input_data+tmpind, input_data + tmpind2, N, lambda, type, tmpmat);
			sum=0;
			for(int k=0;k<N;k++)
				for(int l=0;l<N;l++)
					sum += tmpmat[k][l]*input_data[tmpind+N+k]*input_data[tmpind2+N+l];

			Q[i+M][j+M] = sum;
			Q[j+M][i+M] = Q[i+M][j+M];

		}
	}

	//H_s
	for(int i=0;i<P;i++)
	{
		getsecondkernelderivative(input_data+2*pl[i]*N, target, N, lambda, type, tmpmat);
		MatrixVectorMultipy(tmpmat, input_data+2*pl[i]*N+N, tmp, N,N);
		for(int j=0;j<N;j++)
		{
			Q[i+M][ j+M+P ] = -tmp[j];
			Q[j+M+P][i+M] = Q[i+M][ j+M+P ];
		}
	}

	getsecondkernelderivative(target, target, N, lambda, type, tmpmat);
	for(int i=0;i<N;i++)
		for(int j=i;j<N;j++)
		{
			Q[i+M+P][j+M+P] = tmpmat[i][j];
			Q[j+M+P][i+M+P] = Q[i+M+P][j+M+P];
		}

	plhs[0] = mxCreateDoubleMatrix(M+P+N, M+P+N, mxREAL);
		double* dptr = mxGetPr(plhs[0]);

		count=0;
		for(int i=0;i<M+P+N;i++)
			for(int j=0;j<M+P+N;j++)
				dptr[count++] = Q[j][i];

		delete tmp;

		for(int i=0;i<N;i++)
			delete [] tmpmat[i];
		delete [] tmpmat;

		for(int i=0;i<M+P+N;i++)
			delete [] Q[i];

		delete [] Q;

}

