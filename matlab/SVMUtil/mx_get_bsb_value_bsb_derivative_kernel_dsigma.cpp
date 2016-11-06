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

#include "mex.h"
#include "../../include/util.h"
#ifdef MX_API_VER
#if MX_API_VER < 0x07030000
typedef int mwIndex;
#endif
#endif

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{

	int M = mxGetN(prhs[0]);
	int N = mxGetM(prhs[0])/2;

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

	int size = M+N*M;
	double** Qp = new double*[size];
	for(int i=0;i<size;i++)
		Qp[i] = new double[size];
	double** tmpmat = new double*[N];
	for(int i=0;i<N;i++)
		tmpmat[i] = new double[N];

	double* tmp = new double[N];
	double* tmp2 = new double[N];

	double** Q = new double*[2*size];
	for(int i=0;i<2*size;i++)
		Q[i] = new double[2*size];

	double** Qt2 = new double*[size];
	for(int i=0;i<size;i++)
		Qt2[i] = new double[2*size];

	plhs[0] = mxCreateCellMatrix(N,1);
	mxArray* matptr; double* dptr;

	double kij, tmpval;
	int tmpind, tmpind2;

	for(int i=0;i<M;i++)
	{
		for(int j=0;j<M;j++)
		{
			K[i][j] = getkernel(input_data+2*i*N, input_data + 2*j*N, lambda, type,N);
		}
	}

	//	int l=0;
	for(int l=0;l<N;l++)
	{
		//K
		for(int i=0;i<M;i++)
		{
			tmpind=2*i*N;
			for(int j=i;j<M;j++)
			{
				Qp[i][j] = -K[i][j]*(input_data[tmpind+l]-input_data[2*j*N+l])*(input_data[tmpind+l]-input_data[2*j*N+l]);
				Qp[j][i]=Qp[i][j];
			}
		}

		//G
		for(int i=0;i<M;i++)
		{
			tmpind=2*i*N;
			for(int j=0;j<M;j++)
			{
				tmpind2 = 2*j*N;
				kij = K[i][j];
				for(int k=0;k<N;k++)
					tmp[k] = input_data[tmpind+k]-input_data[tmpind2+k];

				for(int k=0;k<N;k++)
					tmp2[k] = -2*tmp[l]*tmp[l]*kij*tmp[k]*lambda[k];

				tmp2[l] += 2*kij*tmp[l];

				for(int k=0;k<N;k++)
				{
					Qp[i][j+(k+1)*M] = -tmp2[k];
					Qp[j+(k+1)*M][i] = Qp[i][j+(k+1)*M];
				}
			}
		}




		//H
		for(int i=0;i<M;i++)
		{
			tmpind = 2*i*N;
			for(int j=i;j<M;j++)
			{
				tmpind2 = 2*j*N;
				kij = K[i][j];
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

				for(int k=0;k<N;k++)
					for(int p=0;p<N;p++)
					{
						Qp[i+(k+1)*M][j+(p+1)*M] = tmpmat[k][p];
						Qp[j+(p+1)*M][i+(k+1)*M] = tmpmat[k][p];
					}

			}
		}


		for (int i=0;i<N+1;i++)
		{
			for(int j=0;j<M;j++)
			{
				for(int k=0;k<M+N*M;k++)
				{
					Qt2[k][2*i*M+j] = Qp[k][i*M+j];
					Qt2[k][2*i*M+j+M] = -Qp[k][i*M+j];
				}

			}
		}

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

		matptr = mxCreateDoubleMatrix(2*(M+N*M), 2*(M+N*M), mxREAL);
		dptr = mxGetPr(matptr);

		int count=0;
		for(int i=0;i<2*(M+N*M);i++)
			for(int j=0;j<2*(M+N*M);j++)
				dptr[count++] = Q[j][i];

		mxSetCell(plhs[0], l, matptr);
	}

	for(int i=0;i<2*(M+N*M);i++)
		delete [] Q[i];
	delete [] Q;


	for(int i=0;i<M+N*M;i++)
		delete [] Qt2[i];
	delete [] Qt2;

	for(int i=0;i<M;i++)
		delete [] K[i];
	delete [] K;

	delete tmp;
	delete tmp2;
}

