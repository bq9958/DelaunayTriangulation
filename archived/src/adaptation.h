#include <comprImage.h>
#include <eigen.h>

void metric(double3d *Hessien, Mesh *Msh, double3d *Metric, double *TriArea, int normalization);
void rebuildHessien(Mesh *Msh, double3d *Hessien);
void hessien(Mesh *Msh, double *uh, dynamArr *ArrBase, 
    double3d *Hessien, double2d *NablaRUh, double *TriArea);
double L2Projector(double *TriArea, dynamArr *ArrBoucle, double *ListQOI, int NbrBoucle);
void nabla2_uh(int iTri, Mesh *Msh, double2d *NablaRUh, 
    double *nablaxx_uh, double *nablaxy_uh, double *nablayy_uh);
void nabla_uh(int iTri, Mesh *Msh, double *uh, double *nabla_uhx, double *nabla_uhy);
dynamArr *listBase(Mesh *Msh);
dynamArr *boucleDetection(Mesh *Msh, int iPt, int iTri, int *NbrBoucle);
void export_NablaRUh_to_file(const char* filename, Mesh *Msh, double2d *NablaRUh);
void export_Hessien_to_file(const char* filename, Mesh *Msh, double3d *Hessien);