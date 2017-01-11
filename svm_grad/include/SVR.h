#ifndef __SVR_H__
#define __SVR_H__

#include "armadillo"

using namespace arma;
using namespace std;


struct SVRModels{
	unsigned int kernelType;
	unsigned int nbDim;
	unsigned int totalSV;
	double gamma;
	double b;
	mat SVs;  // nbDim X totalSV
	vec alpha; // 1 X totalSV
	double mux;
};

// SVR
class SVR
{
private:
	SVRModels SVRModel;

	vec diffx;
public:
	SVR(char *f_svrmodel);	
	double regression(vec xi);	
};


#endif //__SVR_H__
