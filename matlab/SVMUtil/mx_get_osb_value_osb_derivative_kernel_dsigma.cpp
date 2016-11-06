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

	double** K = new double*[M];
	for(int i=0;i<M;i++)
		K[i] = new double[M];

	double* Kt = new double[M];



	double* tmp = new double[N];
	double* tmp2 = new double[N];
	int tmpind, tmpind2;
	double sum;
	plhs[0] = mxCreateCellMatrix(N,1);
	mxArray* matptr;
	double kij;double tmpval;
	double** tmpmat = new double*[N];
	for(int i=0;i<N;i++)
		tmpmat[i] = new double[N];

	double* dptr;

	for(int i=0;i<M;i++)
		{
			tmpind=2*i*N;
			Kt[i] = getkernel(input_data + tmpind, target, lambda, type,N);
			for(int j=i;j<M;j++)
			{
				K[i][j] = getkernel(input_data + tmpind, input_data + 2*j*N, lambda, type,N);
				K[j][i]=K[i][j];
			}
		}

	for(int l=0;l<N;l++)
	{
	//K
	for(int i=0;i<M;i++)
	{
		tmpind=2*i*N;
		for(int j=i;j<M;j++)
		{
			Q[i][j] = -labels[i]*labels[j]*K[i][j]*pow(input_data[tmpind+l]-input_data[2*j*N+l],2);
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
			kij = K[i][pl[j]];
			for(int k=0;k<N;k++)
				tmp[k] = input_data[tmpind+k]-input_data[tmpind2+k];

			for(int k=0;k<N;k++)
				tmp2[k] = -2*tmp[l]*tmp[l]*kij*tmp[k]*lambda[k];

			tmp2[l] += 2*kij*tmp[l];

			Q[i][j+M]=0;
			for(int k=0;k<N;k++)
				Q[i][j+M] += tmp2[k]*input_data[tmpind2+N+k];
			Q[i][j+M] *= labels[i];

			Q[j+M][i] = Q[i][j+M];
		}
	}

	//G_s

	for(int i=0;i<M;i++)
	{
		kij=Kt[i];
		tmpind = 2*i*N;
		tmpval = input_data[tmpind+l]-target[l];
		for(int k=0;k<N;k++)
		{
			tmp2[k] = -2*tmpval*tmpval*kij*( input_data[tmpind+k]-target[k])*lambda[k];
		}
		tmp2[l] += 2*kij*tmpval;
		for(int j=0;j<N;j++)
		{
			Q[i][ j+M+P ] = -labels[i]*tmp2[j];
			Q[j+M+P][i] = Q[i][ j+M+P ];
		}
	}



	//H
	for(int i=0;i<P;i++)
	{
		tmpind = 2*pl[i]*N;
		for(int j=i;j<P;j++)
		{
			tmpind2 = 2*pl[j]*N;
			kij = K[pl[i]][pl[j]];
			tmpval = (input_data[tmpind + l] - input_data[tmpind2 + l])*(input_data[tmpind + l] - input_data[tmpind2 + l]);
			for(int k=0;k<N;k++)
				tmp[k] = (input_data[tmpind + k] - input_data[tmpind2 + k])*lambda[k];
			for(int k=0;k<N;k++)
				for(int p=0;p<N;p++)
				{
					if(k == p)
						tmpmat[k][p] = 2*kij*tmpval*(2*tmp[k]*tmp[p] - lambda[k]);
					else
						tmpmat[k][p] = 4*kij*tmpval*tmp[k]*tmp[p];
				}
			for(int k=0;k<N;k++)
			{
				tmpmat[l][k] = -4*(input_data[tmpind + l] - input_data[tmpind2 + l])*kij*(1-lambda[l]*tmpval)*tmp[k];
				tmpmat[k][l] = tmpmat[l][k];
			}

			tmpmat[l][l] = -4*tmpval*kij*lambda[l]*(2-lambda[l]*tmpval) + 2*kij*(1-lambda[l]*tmpval);

			sum=0;
			for(int k=0;k<N;k++)
				for(int p=0;p<N;p++)
					sum += tmpmat[k][p]*input_data[tmpind+N+k]*input_data[tmpind2+N+p];

			Q[i+M][j+M] = sum;
			Q[j+M][i+M] = Q[i+M][j+M];

		}
	}

	//H_s
	for(int i=0;i<P;i++)
	{
		tmpind = 2*pl[i]*N;
		kij = Kt[pl[i]];
		tmpval = (input_data[tmpind+l]-target[l])*(input_data[tmpind+l]-target[l]);
		for(int k=0;k<N;k++)
			tmp[k] = (input_data[tmpind+k]-target[k])*lambda[k];

		for(int k=0;k<N;k++)
			for(int p=0;p<N;p++)
			{
				if(k == p)
					tmpmat[k][p] = 2*kij*tmpval*(2*tmp[k]*tmp[p] - lambda[k]);
				else
					tmpmat[k][p] = 4*kij*tmpval*tmp[k]*tmp[p];
			}

		for(int k=0;k<N;k++)
		{
			tmpmat[l][k] = -4*(input_data[tmpind + l] - target[l])*kij*(1-lambda[l]*tmpval)*tmp[k];
			tmpmat[k][l] = tmpmat[l][k];
		}

		tmpmat[l][l] = -4*tmpval*kij*lambda[l]*(2-lambda[l]*tmpval) + 2*kij*(1-lambda[l]*tmpval);

		for(int k=0;k<N;k++)
		{
			tmp[k] = 0;
			for(int p=0;p<N;p++)
				tmp[k] += tmpmat[k][p]*input_data[tmpind + N + p];
		}

		for(int j=0;j<N;j++)
		{
			Q[i+M][ j+M+P ] = -tmp[j];
			Q[j+M+P][i+M] = Q[i+M][ j+M+P ];
		}
	}

	for(int i=0;i<N;i++)
		for(int j=i;j<N;j++)
		{
			Q[i+M+P][j+M+P] = 0;
			Q[j+M+P][i+M+P] = Q[i+M+P][j+M+P];
		}

	Q[l+M+P][l+M+P] = 2;


	matptr = mxCreateDoubleMatrix(M+P+N, M+P+N, mxREAL);
		dptr = mxGetPr(matptr);

		count=0;
		for(int i=0;i<M+P+N;i++)
			for(int j=0;j<M+P+N;j++)
				dptr[count++] = Q[j][i];

		mxSetCell(plhs[0], l, matptr);
	}
		delete tmp;

		for(int i=0;i<N;i++)
			delete [] tmpmat[i];
		delete [] tmpmat;

		for(int i=0;i<M+P+N;i++)
			delete [] Q[i];

		delete [] Q;

}

