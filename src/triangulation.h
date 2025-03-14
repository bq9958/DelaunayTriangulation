#include <mesh.h>
#include <stdbool.h>

typedef struct {
    double x; 
    double y;
} Point;



void kernelDelaunay(double2d *points, int NbrPts, int *NbrTri, double3d *Tri, double3d *Voi);


///////////////////////////// Tool functions //////////////////////////////////////
bool inCircumcircle(double x1, double y1, double x2, double y2, double x3, double y3, double x, double y);
double rCircumcircle(double x1, double y1, double x2, double y2, double x3, double y3);
Point centerCircumcircle(double x1, double y1, double x2, double y2, double x3, double y3);
void barycenter(double x0, double y0, 
    double x1, double y1, 
    double x2, double y2, 
    double x, double y,
    double *beta0, double *beta1, double *beta2);
double aireSignee(double x0, double y0, double x1, double y1, double x2, double y2);
double2d *TestPointSet(int NbrPts);
void save_triangulation_to_file(const char *filename, Mesh *Msh, HashTable *hsh);
void plot_with_gnuplot(const char *filename);
void addElement(int NbrTri, int SizVector, double3d *Tri, double3d *Voi);