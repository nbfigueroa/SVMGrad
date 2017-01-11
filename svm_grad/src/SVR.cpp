
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>

#include "SVR.h"

SVR::SVR(char *f_svrmodel)
{
	unsigned int i,s;
	ifstream fin(f_svrmodel);
	
	// Load SVR Model
	fin >> SVRModel.kernelType
	    >> SVRModel.nbDim     
	    >> SVRModel.totalSV   
	    >> SVRModel.gamma     
	    >> SVRModel.b;

	SVRModel.SVs   = zeros<mat>(SVRModel.nbDim, SVRModel.totalSV);
	SVRModel.alpha = zeros<colvec>(SVRModel.totalSV);

	for( s=0; s<SVRModel.totalSV; s++ ){
		fin >> SVRModel.alpha(s);
	}

	for( i=0; i<SVRModel.nbDim; i++ ){
		for( s=0; s<SVRModel.totalSV; s++ ){		
			fin >> SVRModel.SVs(i,s);
		}
	}
	fin >> SVRModel.mux;
	fin.close();

	diffx = zeros<vec>(SVRModel.nbDim);
}

double SVR::regression(vec xi)
{
	double out;
	unsigned int d,s;
	double k;

	out=0.0;
	
	// Kernel type
	//"	0 -- linear: u'*v\n"
	//"	1 -- polynomial: (gamma*u'*v + coef0)^degree\n"
	//"	2 -- radial basis function: exp(-gamma*|u-v|^2)\n"
	//"	3 -- sigmoid: tanh(gamma*u'*v + coef0)\n"
	//"	4 -- precomputed kernel (kernel values in training_instance_matrix)\n"

	// RBF kernel
	if(SVRModel.kernelType == 2){

		for( s=0; s<SVRModel.totalSV; s++ ){
			for(d=0; d<SVRModel.nbDim; d++) diffx(d) = xi(d)-SVRModel.SVs(d,s);

			k = exp(-SVRModel.gamma*dot(diffx,diffx));
			out += SVRModel.alpha(s)*k;
		}

	}
	// Polynomial kernel
	else if(SVRModel.kernelType = 1){
		//not implemented yet
	}
	
	return (out - SVRModel.b)*SVRModel.mux;
}



//double SVR::regression(vec xi)
//{
//	double out = 0.0;
//	unsigned int s;
//	double k, outr;
//	
//	outr=0;
//	//for( s=0; s<SVRModel.totalSV; s++ ){
//	for( s=0; s<50; s++ ){
//		k = dot(xi,SVRModel.SVs.col(s)) +1.0;
//		outr += SVRModel.alpha(s)*k;
//	}
//	
//	return (outr - SVRModel.b)*SVRModel.mux;
//}
