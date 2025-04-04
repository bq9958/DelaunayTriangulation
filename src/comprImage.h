#include <triangulation.h>

void createPixelSignif(Mesh *MshInit, Mesh *MshDel, int NbrPixelSignif, 
                       double *solInit, int NbPtFr, double *solInterp, double *solQc);
double *readVerSolution(const char *filename, int num_values);
void TriangulationDelaunay(Mesh *MshDel, HashTable *hshDel, int NbrPtInsMin, int NbrPtInsMax);
void interpolateSolution(Mesh *MshInit, Mesh *MshDel, double *solInit, double *solInterp, int NbrPtFr) ;
void solFr(double *solInit, double *solInterp, double *solQC) ;      // only for h = 1
double computePSNR(Mesh *MshInit, Mesh *MshIntern, double *solInit, double *solInterp, double *solPSNR);
void FrPoints(Mesh *MshDel, double h);
int *find_connex_components_modif(Mesh *Msh);
void selectInternalTri(Mesh *MshIntern, Mesh *MshDel, int *color);
HashTable * msh_neighbors_modif(Mesh *Msh, const char* keyMode);