#include <comprImage.h>


/////////// Global variables //////////
const char *pixelSignifMode = "aleatoire";
int SizPil = 100;       // cavity Tri Number max  //TODO Size Control
const char *keyMode = "sum";
const int m = 58;      // pixel width
const int n = 88;      // pixel height



int main(int argc, char *argv[])
{
    /////////// Read mesh of initial image //////////
    // if ( argc < 2 ) {
    //     printf(" usage : mesh file \n");
    //     return 0;
    // }

    int readEfr = 1;
    const char *meshFile = "../data/joconde.lowres.mesh";
    Mesh * MshInit = msh_read(meshFile, readEfr);

    if ( ! MshInit ) return 0;
  
    printf("  Vertices   %10d \n", MshInit->NbrVer);
    printf("  Triangles  %10d \n", MshInit->NbrTri);

    /////////// Define Mesh Delaunay ////////////
    int NbrPixelSignif = 3000;

    Mesh *MshDel = msh_init();
    MshDel->NbrVer = NbrPixelSignif;
    int SizTri = 4*MshDel->NbrVer;   
    MshDel->Crd = (double2d *)calloc(MshDel->NbrVer + 5, sizeof(double2d));
    MshDel->Tri = (int3d *)calloc(SizTri + 1, sizeof(int3d));;  
    MshDel->TriVoi = (int3d *)calloc(SizTri + 1, sizeof(int3d));;

    /////////// Define Hash Table //////////
    HashTable *hshInit = msh_neighbors(MshInit, keyMode);
    HashTable *hshDel = hash_init(0.8*SizTri, 5*SizTri+1);     //TODO MaxObj Siz control


    /////////// Define solution //////////
    const char *solFile = "../data/joconde.lowres.sol";
    double *solInit = sol_read(solFile, 2, MshInit->NbrVer);
    double *solInterp = (double *)calloc(MshDel->NbrVer + 5, sizeof(double));
    double *solPSNR = (double *)calloc(MshInit->NbrVer, sizeof(double));

    /////////// Create pixels significatifs ///////////
    printf(" NbrVer of MshDel = %d\n", MshDel->NbrVer);
    createPixelSignif(MshInit, MshDel, NbrPixelSignif, solInit);
    for (int i=MshDel->NbrVer; i>=1; i--)
    {
        MshDel->Crd[i+4][0] = MshDel->Crd[i][0];
        MshDel->Crd[i+4][1] = MshDel->Crd[i][1];
        MshDel->Crd[i][0] = 0;
        MshDel->Crd[i][1] = 0;
    }
    save_triangulation_to_file("../output/triangulation.mesh", MshDel, hshDel);
    printf("NbrVer of MshDel = %d\n", MshDel->NbrVer);

    // /////////// Triangulation Delaunay ///////////
    save_triangulation_to_file("../output/triangulation.mesh", MshDel, hshDel);
    TriangulationDelaunay(MshDel, hshDel,  -5, m+5, -5, n+5, NbrPixelSignif);
    printf("NbrVer of MshDel = %d\n", MshDel->NbrVer);

    /////////// Interpolation and PSNR ///////////
    interpolateSolution(MshInit, MshDel, solInit, solInterp);
    printf("NbrVer of MshDel = %d\n", MshDel->NbrVer);
    save_triangulation_to_file("../output/triangulation.mesh", MshDel, hshDel);
    double PSNR = computePSNR(MshInit, MshDel, solInit, solInterp, solPSNR);
    printf(" NbrVer of MshDel = %d\n", MshDel->NbrVer);
    printf("PSNR = %f (decibel)\n", PSNR);

    /////////// Free memory ///////////
    msh_free(MshInit);
    msh_free(MshDel);
    hash_free(hshDel);
    hash_free(hshInit);
    return 0;
}