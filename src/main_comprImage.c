#include <comprImage.h>


/////////// Global variables //////////
const char *pixelSignifMode = "aleatoire";   // "aleatoire" ou "bloc"
int SizPil = 100;       // cavity Tri Number max  //TODO Size Control
const char *keyMode = "sum";
const int m = 58;      // pixel width (58/289)    //! Update domain
const int n = 88;      // pixel height (88/440)  //! Update domain



int main(int argc, char *argv[])
{
    /////////// Read mesh of initial image //////////
    const char *meshFile = "../data/joconde.lowres.mesh";       //! Update here

    int readEfr = 1;
    Mesh * MshInit = msh_read(meshFile, readEfr);

    if ( ! MshInit ) return 0;
  
    printf("  Vertices   %10d \n", MshInit->NbrVer);
    printf("  Triangles  %10d \n", MshInit->NbrTri);

    /////////// Define Mesh Delaunay ////////////
    const int NbrPixelSignif = 3500;      //! Update here
    double h_fr = 1; int nx = m / h_fr; int ny = n /h_fr;
    printf("  nx, ny = %d, %d\n", nx, ny);
    int NbPtFr = 2*(nx + ny) - 4;
    printf("  NbPtFr = %d\n", NbPtFr);

    Mesh *MshDel = msh_init();
    MshDel->NbrVer = NbrPixelSignif + NbPtFr;
    int SizTri = 4*MshDel->NbrVer;  //TODO Size control 
    MshDel->Crd = (double2d *)calloc(MshDel->NbrVer + 5, sizeof(double2d));
    MshDel->Tri = (int3d *)calloc(SizTri + 1, sizeof(int3d));;  
    MshDel->TriVoi = (int3d *)calloc(SizTri + 1, sizeof(int3d));;

    /////////// Define Hash Table //////////
    HashTable *hshInit = msh_neighbors(MshInit, keyMode);
    HashTable *hshDel = hash_init(0.8*SizTri, 3*SizTri+1);     //TODO MaxObj Siz control


    /////////// Define solution //////////
    const char *solFile = "../data/joconde.lowres.sol";          //! Update here
    double *solInit = sol_read(solFile, 2, MshInit->NbrVer);
    double *solInterp = (double *)calloc(MshDel->NbrVer + 5, sizeof(double));
    double *solQC = (double *)calloc(MshInit->NbrVer, sizeof(double));

    /////////// Add frontier points ///////////
    FrPoints(MshDel, h_fr);
    save_triangulation_to_file("../output/triangulation.mesh", MshDel, hshDel); 

    /////////// Create pixels significatifs ///////////
    debug_printf(" NbrVer of MshDel = %d\n", MshDel->NbrVer);
    createPixelSignif(MshInit, MshDel, NbrPixelSignif, solInit, NbPtFr);
    save_triangulation_to_file("../output/triangulation.mesh", MshDel, hshDel);
    
    // // /////////// Triangulation Delaunay ///////////
    printf("NbrVer of MshDel : %d\n", MshDel->NbrVer);
    TriangulationDelaunay(MshDel, hshDel,  -0.1*m, 1.1*m, -0.1*n, 1.1*n, NbrPixelSignif);
    save_test_points_to_file("../output/test_points.txt", MshDel);
    save_triangulation_to_file("../output/triangulation.mesh", MshDel, hshDel);  
    debug_printf("NbrVer of MshDel = %d\n", MshDel->NbrVer);
    
    // /////////// Interpolation and PSNR ///////////
    // interpolateSolution(MshInit, MshDel, solInit, solInterp);
    // double PSNR = computePSNR(MshInit, MshDel, solInit, solInterp, solQC);
    // printf("NbrVer of MshDel = %d\n", MshDel->NbrVer);
    // printf("PSNR = %f (decibel)\n", PSNR);

    /////////// Free memory ///////////
    free(solInit);
    free(solInterp);
    free(solQC);
    msh_free(MshInit);
    msh_free(MshDel);
    hash_free(hshDel);
    hash_free(hshInit);
    return 0;
}
