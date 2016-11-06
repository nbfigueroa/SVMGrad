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

	int M = mxGetN(prhs[0]);
	int N = mxGetM(prhs[0])/2;

//	double lambda = mxGetScalar(prhs[1]);
//	printf("sfasdas");
	double* lambda;
	if(mxGetNumberOfElements(prhs[1]) == 1)
	{
		lambda = new double[N];
		double tmp = mxGetScalar(prhs[1]);
		for(int i=0;i<N;i++)
			lambda[i] = tmp;
	}
	else if(mxGetNumberOfElements(prhs[1]) == N)
	{
		lambda = (double*)mxGetPr(prhs[1]);
	}
	else
	{
		printf("ERROR: lambda can have length = 1 or equal to the number of dimensions \n");
		return;
	}

	char type[1025];
	mxGetString(prhs[2],type, 1025 );
	double* input_data = (double*)mxGetPr(prhs[0]);
	double** K = new double*[M];
	for(int i=0;i<M;i++)
		K[i] = new double[M];

	double** G = new double*[M];
	for(int i=0;i<M;i++)
		G[i] = new double[N*M];

	double** H = new double*[N*M];
	for(int i=0;i<N*M;i++)
		H[i] = new double[N*M];


	for(int i=0;i<M;i++)
	{
		for(int j=0;j<M;j++)
		{
			K[i][j] = getkernel(input_data+2*i*N, input_data + 2*j*N, lambda, type,N);
		}
	}

	double tmp[N];
	for(int b1=0;b1<N;b1++)
		for(int i=0;i<M;i++)
		{
			for(int j=0;j<M;j++)
			{
				getfirstkernelderivative(input_data+2*i*N, input_data + 2*j*N, lambda, type, 2, tmp, N);
				G[i][j+M*b1] = tmp[b1];
			}
		}

	double** tmpmat = new double*[N];
	for(int i=0;i<N;i++)
		tmpmat[i] = new double[N];
	for(int b1=0;b1<N;b1++)
		for(int b2=0;b2<N;b2++)
			for(int i=0;i<M;i++)
			{
				for(int j=0;j<M;j++)
				{

					getsecondkernelderivative(input_data+2*i*N, input_data + 2*j*N, N, lambda, type,tmpmat);
					H[i+b1*M][ j+M*b2 ] = tmpmat[b1][b2];
				}
			}

	double** Qt = new double*[M+N*M];
	for(int i=0;i<M+N*M;i++)
		Qt[i] = new double[M+N*M];

	for(int i=0;i<M;i++)
	{
		for(int j=0;j<M;j++)
			Qt[i][j] = K[i][j];
		for(int j=M;j<M+N*M;j++)
			Qt[i][j] = -G[i][j-M];
	}
	for(int i=M;i<M+N*M;i++)
	{
		for(int j=0;j<M;j++)
			Qt[i][j] = -G[j][i-M];
		for(int j=M;j<M+N*M;j++)
			Qt[i][j] = H[i-M][j-M];
	}

	for(int i=0;i<M;i++)
	{
		delete [] K[i];
		delete [] G[i];
	}
	delete [] K;
	delete [] G;

	for(int i=0;i<M*N;i++)
	{
		delete [] H[i];
	}
	delete [] H;


	double** Qt2 = new double*[(M+N*M)];
	for(int i=0;i<(M+N*M);i++)
		Qt2[i] = new double[2*(M+N*M)];


	for (int i=0;i<N+1;i++)
	{
		for(int j=0;j<M;j++)
		{
			for(int k=0;k<M+N*M;k++)
			{
					Qt2[k][2*i*M+j] = Qt[k][i*M+j];
					Qt2[k][2*i*M+j+M] = -Qt[k][i*M+j];
			}

		}
	}

	for(int i=0;i<M+N*M;i++)
		delete [] Qt[i];
	delete [] Qt;


	double** Q = new double*[2*(M+N*M)];
	for(int i=0;i<2*(M+N*M);i++)
		Q[i] = new double[2*(M+N*M)];

	for (int i=0;i<N+1;i++)
	{
		for(int j=0;j<M;j++)
		{
			for(int k=0;k<2*(M+N*M);k++)
				{
					Q[2*i*M+j][k] = Qt2[i*M+j][k];
					Q[2*i*M+j+M][k] = -Qt2[i*M+j][k];
				}
		}
	}

	for(int i=0;i<M+N*M;i++)
		delete [] Qt2[i];
	delete [] Qt2;

	plhs[0] = mxCreateDoubleMatrix(2*(M+N*M), 2*(M+N*M), mxREAL);
	double* dptr = mxGetPr(plhs[0]);

	int count=0;
	for(int i=0;i<2*(M+N*M);i++)
		for(int j=0;j<2*(M+N*M);j++)
			dptr[count++] = Q[j][i];

	for(int i=0;i<2*(M+N*M);i++)
		delete [] Q[i];
	delete [] Q;

//	plhs[1] = mxCreateDoubleMatrix((M+N*M), (M+N*M), mxREAL);
//	double* dptr2 = mxGetPr(plhs[1]);
//
//	count=0;
//	for(int i=0;i<(M+N*M);i++)
//		for(int j=0;j<(M+N*M);j++)
//			dptr2[count++] = Qt[j][i];



}

