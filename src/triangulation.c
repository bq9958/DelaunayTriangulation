#include <triangulation.h>

void kernelDelaunay(double2d *points, int NbrPts, int *NbrTri, double3d *Tri, double3d *Voi)
{
    

}



///////////////////////////////// Tool functions //////////////////////////////////////


double2d *TestPointSet(int NbrPts)
{
    srand(NbrPts);

    // Pre-allocate memory for the points of boiteGlobante
    double2d *Pts = (double2d *)malloc((NbrPts + 5) * sizeof(double2d));

    if (!Pts) {
        printf("Error: Memory allocation failed!\n");
        return NULL;
    }

    for (int i = 1; i <= NbrPts; i++) {
        Pts[i][0] = (double)rand() / RAND_MAX;
        Pts[i][1] = (double)rand() / RAND_MAX;
    }

    return Pts;
}

bool inCircumcircle(double x1, double y1, double x2, double y2, double x3, double y3, double x, double y)
{
    // determine if (x,y) is inside the circumcircle of the triangle (x1,y1), (x2,y2), (x3,y3)
    double radius = rCircumcircle(x1, y1, x2, y2, x3, y3);
    Point center = centerCircumcircle(x1, y1, x2, y2, x3, y3);
    double dist = distance(center.x, center.y, x, y);
    if (dist < radius)
        return true;
    return false;
}

double rCircumcircle(double x1, double y1, double x2, double y2, double x3, double y3)
{
    double l1 = distance(x1, y1, x2, y2);
    double l2 = distance(x2, y2, x3, y3);
    double l3 = distance(x3, y3, x1, y1);
    double area = triArea(x1, y1, x2, y2, x3, y3);

    return l1 * l2 * l3 / (4.0 * area); 
}

Point centerCircumcircle(double x1, double y1, double x2, double y2, double x3, double y3)
{
    Point center;

    center.x = ((y2-y1)*(y3*y3-y1*y1+x3*x3-x1*x1)-(y3-y1)* \ 
               (y2*y2-y1*y1+x2*x2-x1*x1))/(2.0*((x3-x1)*(y2-y1)-(x2-x1)*(y3-y1)));
    center.y = ((x2-x1)*(x3*x3-x1*x1+y3*y3-y1*y1)-(x3-x1)* \ 
               (x2*x2-x1*x1+y2*y2-y1*y1))/(2.0*((y3-y1)*(x2-x1)-(y2-y1)*(x3-x1)));

    return center;
}

void barycenter(double x0, double y0, 
                  double x1, double y1, 
                  double x2, double y2, 
                  double x, double y,
                  double *beta0, double *beta1, double *beta2)
{
    double aireNonSignee = triArea(x0, y0, x1, y1, x2, y2);
    double PP1P2 = aireSignee(x, y, x1, y1, x2, y2);
    double P0PP2 = aireSignee(x0, y0, x, y, x2, y2);
    double P0P1P = aireSignee(x0, y0, x1, y1, x, y);

    *beta0 = PP1P2 / aireNonSignee;
    *beta1 = P0PP2 / aireNonSignee;
    *beta2 = P0P1P / aireNonSignee;
}

double aireSignee(double x0, double y0, double x1, double y1, double x2, double y2)
{
    double aireSignee = 0.5 * ((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));
    return aireSignee;
}

void save_triangulation_to_file(const char *filename, Mesh *Msh, HashTable *hsh) 
{
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        printf("Error: Cannot open file %s for writing!\n", filename);
        return;
    }

    fprintf(fp, "# Vertex\n");
    for (int i = 1; i <= Msh->NbrVer + 4; i++) {
        fprintf(fp, "%lf %lf\n", Msh->Crd[i][0], Msh->Crd[i][1]);
    }

    fprintf(fp, "\n\n# Edges\n");
    for (int i = 1; i <= hsh->NbrObj; i++) {
        fprintf(fp, "%lf %lf\n", Msh->Crd[hsh->LstObj[i][0]][0], Msh->Crd[hsh->LstObj[i][0]][1]);
        fprintf(fp, "%lf %lf\n", Msh->Crd[hsh->LstObj[i][1]][0], Msh->Crd[hsh->LstObj[i][1]][1]);
    }

    fclose(fp);
    printf("[Ouput File] triangulation written to %s\n", filename);
}

void plot_with_gnuplot(const char *filename) {
    FILE *gnuplot = popen("gnuplot -persist", "w");
    if (!gnuplot) {
        printf("Error: Cannot open Gnuplot!\n");
        return;
    }

    fprintf(gnuplot, "set title 'Point Set Visualization'\n");
    fprintf(gnuplot, "set xlabel 'X'\n");
    fprintf(gnuplot, "set ylabel 'Y'\n");
    fprintf(gnuplot, "set grid\n");
    fprintf(gnuplot, "unset key\n");

    fprintf(gnuplot, "plot '%s' index 0 with points pt 7 lc rgb 'blue', \\\n", filename);
    fprintf(gnuplot, "     '%s' index 1 with lines lw 2 lc rgb 'red' \n", filename);
    fflush(gnuplot);
    pclose(gnuplot);
}

void addElement(int NbrTri, int SizVector, double3d *Tri, double3d *Voi)
{
    if (NbrTri >= SizVector) {
        SizVector *= 2;
        double3d *tmp1 =(double3d *)realloc(Tri, SizVector * sizeof(double3d));
        double3d *tmp2 = (double3d *)realloc(Voi, SizVector * sizeof(double3d));

        if (!tmp1 || !tmp2) {
            printf("Error: Memory reallocation failed!\n");
            return;
        }
        else {
            Tri = tmp1;
            Voi = tmp2;
        }
    }
}