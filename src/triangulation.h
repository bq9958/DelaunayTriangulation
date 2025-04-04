#include <dynamicArray.h>
#include <time.h>
#include <stdbool.h>

typedef struct {
    double x; 
    double y;
} Point;



void Cavity(Mesh *Msh, HashTable *hsh, int iPtIns, int *iTriLocLast);
int location(Mesh *Msh, int iTri, double x, double y, int *move);
void boucleDetection(Mesh *Msh, int iPt, int iTri, int *mark, int *step);
int inTriangle(double x0, double y0, double x1, double y1, double x2, double y2, double x, double y);
int inCircumcircle(double x1, double y1, double x2, double y2, double x3, double y3, double x, double y);
double rCircumcircle(double x1, double y1, double x2, double y2, double x3, double y3);
Point centerCircumcircle(double x1, double y1, double x2, double y2, double x3, double y3);
void barycenter(double x0, double y0, 
    double x1, double y1, 
    double x2, double y2, 
    double x, double y,
    double *beta0, double *beta1, double *beta2);
double aireSignee(double x0, double y0, double x1, double y1, double x2, double y2);
double2d *TestPointSet(int NbrPts, double xmin, double xmax, double ymin, double ymax);
void save_triangulation_to_file(const char *filename, Mesh *Msh, int NbrPtFr);
int inList(int element, int *list, int SizList);
void print_pile_to_file(int *pile);
void addElement(int NbrTri, int SizVector, double3d *Tri, double3d *Voi);
void supprAreteCom(int iVerCom1, int iVerCom2, dynamArr *PilEdg);
void addTwoOtherSides(Mesh *Msh, int iTri, int iVerComLocal1, int iVerComLocal2, dynamArr *PilEdg);
void areteCommune(Mesh *Msh, int iTri, int jTri, int *iVerCom1, int *iVerCom2,
                  int *iVerComLocal1, int *iVerComLocal2);
int hash_find_modified(HashTable *hsh, int iVer1, int iVer2, const char* keyMode);
void save_test_points_to_file(const char *filename, Mesh *Msh);