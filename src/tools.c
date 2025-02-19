#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

double distance(double x1, double y1, double x2, double y2)
{
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

double triArea(double x0, double y0, double x1, double y1, double x2, double y2)
{
    return 0.5 * ((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));
}


// int main()
// {
//     double area = triArea(0, 0, 1, 0, 0, 1);
//     double dist = distance(0, 0, 1, 1);
//     printf("rhoK: %f\n", rho);
//     printf("Distance between (0, 0) and (1, 1): %f\n", dist);
//     printf("Area of triangle: %f\n", area);
//     return 0;
// }