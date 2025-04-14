#include <adaptation.h>

double computePSNR_modif(Mesh *MshInit, Mesh *MshIntern, double *solInit, double *solInterp, double *solQC);
void ifOriPixel(double *solInit, double x, double y);
