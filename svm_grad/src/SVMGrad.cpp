#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>

#include "SVMGrad.h"

SVMGrad::SVMGrad(char *f_SVMGradmodel)
{
//	unsigned int i,s;
//    ifstream fin(f_SVMGradmodel);
	
//    // Load SVMGrad Model
//    fin >> SVMGradModel.kernelType
//        >> SVMGradModel.D
//        >> SVMGradModel.nSV
//        >> SVMGradModel.gamma
//        >> SVMGradModel.b;

//    SVMGradModel.SVs   = zeros<mat>(SVMGradModel.D, SVMGradModel.totalSV);
//    SVMGradModel.yalphas = zeros<colvec>(SVMGradModel.totalSV);

//    for( s=0; s<SVMGradModel.totalSV; s++ ){
//        fin >> SVMGradModel.yalphas(s);
//	}

//    for( i=0; i<SVMGradModel.D; i++ ){
//        for( s=0; s<SVMGradModel.totalSV; s++ ){
//            fin >> SVMGradModel.SVs(i,s);
//		}
//	}
//    fin >> SVMGradModel.mux;
//	fin.close();

//    diffx = zeros<vec>(SVMGradModel.D);
}

double SVMGrad::calculateGamma(vec xi)
{
//	double out;
//	unsigned int d,s;
//	double k;

//	out=0.0;

//	// RBF kernel
//    if(SVMGradModel.kernelType == 2){

//        for( s=0; s<SVMGradModel.totalSV; s++ ){
//            for(d=0; d<SVMGradModel.D; d++) diffx(d) = xi(d)-SVMGradModel.SVs(d,s);

//            k = exp(-SVMGradModel.gamma*dot(diffx,diffx));
//            out += SVMGradModel.yalphas(s)*k;
//		}

//	}
	
//    return (out - SVMGradModel.b)*SVMGradModel.mux;
}

double SVMGrad::calculateGammaDerivative(vec xi)
{
//    double out;
//    unsigned int d,s;
//    double k;

//    out=0.0;

//    // RBF kernel
//    if(SVMGradModel.kernelType == 2){

//        for( s=0; s<SVMGradModel.totalSV; s++ ){
//            for(d=0; d<SVMGradModel.D; d++) diffx(d) = xi(d)-SVMGradModel.SVs(d,s);

//            k = exp(-SVMGradModel.gamma*dot(diffx,diffx));
//            out += SVMGradModel.yalphas(s)*k;
//        }

//    }
//    return (out - SVMGradModel.b)*SVMGradModel.mux;
}
