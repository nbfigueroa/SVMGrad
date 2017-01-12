#include "SVMGrad.h"

SVMGrad::SVMGrad(string& f_SVMGradmodel)
{

    cout << "SVM Model: " << f_SVMGradmodel << endl;
    ifstream fin(f_SVMGradmodel.c_str());
    unsigned int    d, s;

    // Load SVMGrad Model
    fin >> SVMGradModel.D
        >> SVMGradModel.nSV
        >> SVMGradModel.b
        >> SVMGradModel.sigma;

    cout << "model.D: " << SVMGradModel.D << endl
         << "model.nSV: " << SVMGradModel.nSV << endl
         << "model.b: " << SVMGradModel.b << endl
         << "model.sigma: " << SVMGradModel.sigma << endl;

    SVMGradModel.yalphas = zeros<colvec>(SVMGradModel.nSV);
    SVMGradModel.SVs     = zeros<mat>(SVMGradModel.D, SVMGradModel.nSV);

    for( s=0; s<SVMGradModel.nSV; s++ ){
        fin >> SVMGradModel.yalphas(s);
    }

    for( d=0; d<SVMGradModel.D; d++ ){
        for( s=0; s<SVMGradModel.nSV; s++ ){
            fin >> SVMGradModel.SVs(d,s);
        }
    }
    fin.close();

    lambda = 1/(2*SVMGradModel.sigma*SVMGradModel.sigma);
}


double SVMGrad::calculateClass(vec x)
{
    y = calculateGamma(x);
    return y<0? -1 : y>0;
}

double SVMGrad::calculateGamma(vec x)
{

    double gamma_val = 0.0000;

    for(unsigned int s=0; s<SVMGradModel.nSV; s++ ){
        // \alpha_i*y_i*exp(-lamdba*||x-x_i||^2))
        gamma_val += SVMGradModel.yalphas(s)*getKernel(x, s);
    }

    return (gamma_val + SVMGradModel.b);
}

double SVMGrad::getKernel(vec x, unsigned int s){

    double kernel = 0.0;

    // Compute vector diff
    diffx = zeros<vec>(SVMGradModel.D);
    for(unsigned int d=0; d<SVMGradModel.D; d++)
        diffx(d) = x(d) - SVMGradModel.SVs(d,s);

    // RBF Kernel Exponential
    kernel = exp(-lambda*pow(norm(diffx),2));

    return kernel;
}


