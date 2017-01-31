#include <stdio.h>
#include <fstream>
#include <time.h>
#include "armadillo"
#include "eigen3/Eigen/Dense"
#include "svm_grad.h"

using namespace std;
using namespace arma;

struct SVMTestData{
    int D;         // dimensions
    int M;         // samples
    mat x;         // input
    vec y;         // output
    vec value;     // expected Gamma(X) value
    mat gradient;  // expected DGamma(X) value
};

void loadTestData(string& f_SVMTestData, SVMTestData& data)
{

    cout << "\nSVM Testing File: " << f_SVMTestData << endl;
    ifstream fin(f_SVMTestData.c_str());
    unsigned int d , s;

    // Load SVMTest Data
    fin >> data.D
        >> data.M;

    cout << "data.D: " << data.D << endl
         << "data.M: " << data.M << endl;

    data.x = zeros<mat>(data.D, data.M);
    data.y        = zeros<colvec>(data.M);
    data.value    = zeros<colvec>(data.M);
    data.gradient = zeros<mat>(data.D, data.M);


    for( d=0; d<data.D; d++ ){
        for( s=0; s<data.M; s++ ){
            fin >> data.x(d,s);
        }
    }

    for( s=0; s<data.M; s++ ){
        fin >> data.y(s);
    }

    for( s=0; s<data.M; s++ ){
        fin >> data.value(s);
    }

    for( d=0; d<data.D; d++ ){
        for( s=0; s<data.M; s++ ){
            fin >> data.gradient(d,s);
        }
    }

    fin.close();
}


int
main (int argc, char **argv)
{
    string package_path  = "/home/nbfigueroa/dev/catkin_ws/src/SVMGrad/";
    string modelFilename = package_path + "matlab/models/Fender/36D-054k-Optimal-Model-Fender.txt";
    string dataFilename  = package_path + "matlab/models/Fender/36D-054k-Data-Fender.txt";

    SVMGrad svm_(modelFilename);
    SVMTestData data_;
    loadTestData(dataFilename,data_);

    /* Testing Gamma and DGamma function with Armadillo Inputs */
    double gamma;
    vec x = zeros<vec>(data_.D);
    vec gamma_grad      = zeros<vec>(data_.D);
    vec gamma_grad_mat  = zeros<vec>(data_.D);
    vec gamma_diff      = zeros<vec>(data_.M);
    vec gamma_grad_diff = zeros<vec>(data_.M);
    vec gamma_time      = zeros<vec>(data_.M);
    int samples         = data_.M;

    for (unsigned int  i=0;i<samples;i++){

        x = data_.x.col(i);
        gamma_grad_mat = data_.gradient.col(i);

        svm_.preComputeKernel(true);

        clock_t t;
        t = clock();
        {
            gamma         = svm_.calculateGamma(x);
            gamma_grad    = svm_.calculateGammaDerivative(x);
        }
        t = clock() - t;

        gamma_time(i)      = ((float)t)/CLOCKS_PER_SEC;
        gamma_diff(i)      = abs(gamma - data_.value(i));
        gamma_grad_diff(i) = norm(gamma_grad_mat-gamma_grad);

    }

    cout << "\n---Testing Gamma and DGamma function with Armadillo Inputs---\n";
    cout << "Average Gamma Numerical Error: "
         << sum(gamma_diff)/(double)samples << endl;
    cout << "Average Norm DGamma Numerical Error: "
         << sum(gamma_grad_diff)/(double)samples << endl;
    cout << "Average Time -Gamma+DGamma Calculation-: "
         << sum(gamma_time)/(double)samples << endl;


    /* Testing Gamma and DGamma function with Eigen Inputs */
    Eigen::VectorXd x_eig;
    Eigen::VectorXd grad_eig;
    Eigen::VectorXd grad_eig_mat;
    Eigen::VectorXd grad_eig_diff;
    x_eig.resize(x.size());
    grad_eig.resize(gamma_grad.size());
    grad_eig_mat.resize(gamma_grad.size());
    grad_eig_diff.resize(gamma_grad.size());

    for (unsigned int  i=0;i<samples;i++){

        for(unsigned int d=0; d<data_.D; d++){
            x_eig(d)        = data_.x(d,i);
            grad_eig_mat(d) = data_.gradient(d,i);
        }

        svm_.preComputeKernel(true);

        clock_t t;
        t = clock();
        {
           gamma       = svm_.calculateGamma(x_eig);
           grad_eig    = svm_.calculateGammaDerivative(x_eig);

        }
        t = clock() - t;

        gamma_time(i)      = ((float)t)/CLOCKS_PER_SEC;
        gamma_diff(i)      = abs(gamma - data_.value(i));
        grad_eig_diff      = grad_eig_mat-grad_eig;
        gamma_grad_diff(i) = grad_eig_diff.norm();

    }

    cout << "\n---Testing Gamma and DGamma function with Eigen Inputs---\n";
    cout << "Average Gamma Numerical Error: "
         << sum(gamma_diff)/(double)samples << endl;
    cout << "Average Norm DGamma Numerical Error: "
         << sum(gamma_grad_diff)/(double)samples << endl;
    cout << "Average Time -Gamma+DGamma Calculation-: "
         << sum(gamma_time)/(double)samples << endl;


    return 0;
}
