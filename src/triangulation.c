#include <triangulation.h>

void kernelDelaunay(double2d *points, int NbrPts, int *NbrTri, double3d *Tri, double3d *Voi)
{
    

}

void boucleDetection(Mesh *Msh, int iPt, int iTri, int *mark, int *step)
{
    int search = 3;
    mark[iTri] = 1;

    printf("Step %d\n", *step); 
    (*step)++;

    for (int iEdglocal=0; iEdglocal < 3; iEdglocal++){
        int neighbor = Msh->TriVoi[iTri][iEdglocal];
        if (mark[neighbor] == 0 && (Msh->Tri[neighbor][0] == iPt || Msh->Tri[neighbor][1] == iPt 
            || Msh->Tri[neighbor][2] == iPt)){
            mark[neighbor] = 1;     // belongs to boucle
            boucleDetection(Msh, iPt, neighbor, mark, step);
        }
        else{
            if (mark[neighbor] != 1)
                mark[neighbor] = -1;    // not belong to boucle
            search--;
        }
        if (search == 0) {return;}
    }
}

int location(Mesh *Msh, int iTri, double x, double y, int *move)
{
    double x0, y0, x1, y1, x2, y2;
    int position;

    printf("Move %d\n", *move); 
    (*move)++;

    x0 = Msh->Crd[Msh->Tri[iTri][0]][0]; y0 = Msh->Crd[Msh->Tri[iTri][0]][1];
    x1 = Msh->Crd[Msh->Tri[iTri][1]][0]; y1 = Msh->Crd[Msh->Tri[iTri][1]][1];
    x2 = Msh->Crd[Msh->Tri[iTri][2]][0]; y2 = Msh->Crd[Msh->Tri[iTri][2]][1];

    position = inTriangle(x0, y0, x1, y1, x2, y2, x, y);
    if (position == -1)
    {
        return iTri;
    }
    else {
        iTri = Msh->TriVoi[iTri][position];
        location(Msh, iTri, x, y, move);
    }
}

int inTriangle(double x0, double y0, double x1, double y1, double x2, double y2, double x, double y)
{
    double beta[3];
    int choice[3] = {0, 0, 0};
    int count = 0;

    barycenter(x0, y0, x1, y1, x2, y2, x, y, &beta[0], &beta[1], &beta[2]);
    for (int i=0; i<3; i++){
        if (beta[i] < 0)
        {
            count++;
            choice[count] = i;
        }
    }
    if (count == 1)    // one beta is negative
    {
        return choice[count];
    }
    else if (count == 2)  // two betas are negative
    {
        srand(time(NULL));
        int random_number = (rand() % 2) ? 2 : 1;
        return choice[random_number];
    }
    else if (count == 0)  // all beta are positive
    {
        return -1;
    }
    else
    {
        printf("Error: Invalid location\n");
        return -2;
    }
}


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
    
    fprintf(gnuplot, "plot '%s' index 0 using 1:2 with points pt 7 lc rgb 'blue', \\\n", filename);
    fprintf(gnuplot, "     '%s' index 1 using 1:2 with lines lw 2 lc rgb 'red', \\\n", filename);
    fprintf(gnuplot, "     '%s' index 0 using 1:2:(sprintf('%%d', $0+1)) with labels offset 0.3,0.3 font ',10' textcolor rgb 'black'\n", filename);
    
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