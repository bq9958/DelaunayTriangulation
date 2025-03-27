#include <triangulation.h>

void createPixelSignif(Mesh *MshInit, Mesh *MshDel, int NbrPixelSignif, double *solInit);
double *readVerSolution(const char *filename, int num_values);
void TriangulationDelaunay(Mesh *MshDel, HashTable *hshDel, 
    double box_xmin, double box_xmax, double box_ymin, double box_ymax,
    int NbrVerDel);
    void interpolateSolution(Mesh *MshInit, Mesh *MshDel, double *solInit, double *solInterp) ;
double computePSNR(Mesh *MshInit, Mesh *MshDel, double *solInit, double *solInterp, double *solPSNR);