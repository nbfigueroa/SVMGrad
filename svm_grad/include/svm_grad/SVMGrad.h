#ifndef __SVMGrad_H__
#define __SVMGrad_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include "armadillo"

using namespace arma;
using namespace std;

/* SVMGrad: This structure holds the parameters of an RBF-SVM
 * D:       Datapoint Dimension
 * nSV:     # of SVs
 * b:       Offset
 * sigma:   kernel width
 * yalphas: alpha_i*yi     (1 X nSV)
 * SVs:     Support Vector (D X nSV)
*/

struct SVMGradModels{
    unsigned int D;
    unsigned int nSV;
    double sigma;
    double b;
    vec yalphas;
    mat SVs;
};

/* SVMGrad: This class computes the following function from an SMVGradModel
 *       y      = sign(Gamma(x))
 *       Gamma  = \sum_{i=1}^{N_sv}\alpha_iy_ik(x,x_i) + b
 *       DGamma = \sum_{i=1}^{N_sv}-1/2\sigma^2\alpha_iy_ik(x,x_i)(x-x_i)
*/

class SVMGrad
{
private:
    SVMGradModels SVMGradModel;
    vec    diffx;
    double lambda; // lambda = 1/(2*sigma^2)
    double y;      // y      = sign(Gamma(x))

public:       
    SVMGrad(string& f_SVMGradmodel);

    // Armadillo input
    double calculateClass(vec x);
    double calculateGamma(vec x);
    double getKernel(vec x, unsigned int s);

    vec  getKernelDerivative(vec x);
    vec calculateGammaDerivative(vec x);


    // Eigen input
//    void eigen2arma(Eigen::VectorXf x, vec& x);
//    void arma2eigen(vec x, Eigen::VectorXf& x);
//    double calculateGamma(Eigen::VectorXf xi);
//    Eigen::VectorXf double calculateGammaDerivative(Eigen::VectorXf xi);

};


#endif //__SVMGrad_H__
