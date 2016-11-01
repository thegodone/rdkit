#include <cstdlib>
#include <iostream>
#include "Data3Ddescriptors.h"

using namespace std;


double Data3Ddescriptors::mw[53]={0.084,0,0,0,0.900,1.000,1.166,1.332,1.582,0,0,0, 2.246, 2.339, 2.579,2.670, 2.952,0,0,0,0,0,0,0,0, 4.650, 4.907, 4.887, 5.291, 5.445,0,0,0,0, 6.653,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 9.884,0,0, 10.56};
double Data3Ddescriptors::vdW[53]={0.299,0,0,0,0.796,1.000,0.695,0.512,0.410,0,0,0, 1.626, 1.424, 1.181,1.088, 1.035,0,0,0,0,0,0,0,0, 1.829, 1.561, 0.764, 0.512, 1.708,0,0,0,0, 1.384,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 2.042,0,0, 1.728};
double Data3Ddescriptors::neg[53]={0.944,0,0,0,0.828,1.000,1.163,1.331,1.457,0,0,0, 0.624, 0.779, 0.916,1.077, 1.265,0,0,0,0,0,0,0,0, 0.728, 0.728, 0.728, 0.740, 0.810,0,0,0,0, 1.172,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0.837,0,0, 1.012};
double Data3Ddescriptors::pol[53]={0.379,0,0,0,1722,1.000,0.625,0.456,0.316,0,0,0, 3.864, 3.057, 2.063,1.648, 1.239,0,0,0,0,0,0,0,0, 4.773, 4.261, 3.864, 3.466, 4.034,0,0,0,0, 1.733,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 4.375,0,0, 3.040};
double Data3Ddescriptors::ionpol[83]={1.208,0,0.479,0.828,0.737,1,1.291,1.209,1.547,0,0.456,0.679,0.532,0.724,0.931,0.92,1.152,0,0.386,0.543,0,0,0,0.601,0.66,0.702,0.7,0.679,0.686,0.834,0.533,0.702,0.872,0.866,1.049,0,0.371,0.506,0,0,0,0.63,0,0,0,0,0.673,0.799,0.514,0.652,0.767,0.8,0.928,0,0,0,0,0,0,0,0,0,0,0.546,0,0,0,0,0,0,0,0,0,0,0,0,0,0.799,0.819,0.927,0.542,0.659,0.647};
double Data3Ddescriptors::rcov[83]={0.37,0,1.34,0.90,0.82,0.77,0.73,0.71,0,1.54,1.30,1.18,1.11,1.06,1.02,0.99,0,1.96,1.74,0,0,0,1.27,1.39,1.25,1.26,1.21,1.38,1.31,1.26,1.22,1.19,1.16,1.14,0,2.11,1.92,0,0,0,1.45,0,0,0,0,1.53,1.48,1.44,1.41,1.38,1.35,1.33,0,0,0,0,0,0,0,0,0,0,1.79,0,0,0,0,0,0,0,0,0,0,0,0,0,1.28,1.44,1.49,1.48,1.47, 1.46};


Data3Ddescriptors::Data3Ddescriptors(){}

double* Data3Ddescriptors::getMW( ){
  	return mw;
}

double* Data3Ddescriptors::getVDW( ){
	return vdW;
}

double* Data3Ddescriptors::getNEG( ){
	return neg;
}

double* Data3Ddescriptors::getPOL( ){
	return pol;
}

double* Data3Ddescriptors::getIonPOL( ){
	return ionpol;
}

double* Data3Ddescriptors::getRCOV( ){
	return rcov;
}

