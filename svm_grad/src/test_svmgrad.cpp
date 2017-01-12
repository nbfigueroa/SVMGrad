#include <stdio.h>
#include <fstream>
#include <time.h>
#include "armadillo"
#include "SVMGrad.h"

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

    cout << "SVM Testing File: " << f_SVMTestData << endl;
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
    std::string package_path  = "/home/nbfigueroa/dev/catkin_ws/src/SVMGrad/";
    std::string modelFilename = package_path + "matlab/models/36d-robotcollision-svm.txt";
    std::string dataFilename  = package_path + "matlab/models/36d-robotcollision-data.txt";

    SVMGrad svm_(modelFilename);
    SVMTestData data_;
    loadTestData(dataFilename,data_);

    // Test Gamma and DGamma
    double gamma;
    vec x = zeros<vec>(data_.D);
    vec gamma_diff = zeros<vec>(data_.M);
    vec gamma_time = zeros<vec>(data_.M);
    int samples = data_.M;

    for (unsigned int  i=0;i<samples;i++){

        for(unsigned int d=0; d<data_.D; d++)
            x(d) = data_.x(d,i);

        clock_t t;
        t = clock();
        gamma         = svm_.calculateGamma(x);
        t = clock() - t;
        gamma_time(i) = ((float)t)/CLOCKS_PER_SEC;
        gamma_diff(i) = abs(gamma - data_.value(i));

    }

    cout << "Average Gamma Numerical Error: "    << sum(gamma_diff)/(double)samples << endl;
    cout << "Average Time -Gamma Calculation-: " << sum(gamma_time)/(double)samples << endl;


    return 0;
}
