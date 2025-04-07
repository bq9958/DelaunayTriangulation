#include <comprImage.h>

/////////// Global variables //////////
const char *pixelSignifMode = "bloc";   // "regulier" ou "aleatoire" ou "bloc"
int SizPil = 100;       // cavity Tri Number max  //TODO Size Control
const char *keyMode = "sum";
const int m = 289;      // pixel width (58/289)    //! Update domain
const int n = 440;      // pixel height (88/440)  //! Update domain
const int incrementRegulier = 9;   // increment for regular pixel significatif


int main(int argc, char *argv[])
{
    /////////// Read mesh of initial image //////////
    const char *meshFile = "../data/joconde.mesh";       //! Update here

    printf("[ [Read initial mesh] ]\n");
    int readEfr = 1;
    Mesh * MshInit = msh_read(meshFile, readEfr);

    if ( ! MshInit ) return 0;
  
    printf("  Vertices   %10d \n", MshInit->NbrVer);
    printf("  Triangles  %10d \n", MshInit->NbrTri);

    /////////// Define Mesh Delaunay ////////////
    int NbrPixelSignif = 10000;      //! Update here, NbrPixelSignif < MshInit->NbrVer
    double h_fr = 1; int nx = m / h_fr; int ny = n /h_fr;
    printf("  nx, ny = %d, %d\n", nx, ny);
    int NbrPtFr = 2*(nx + ny) - 4;
    printf("  NbPtFr = %d\n", NbrPtFr);
    if (strcmp(pixelSignifMode, "regulier") == 0)
    {
        NbrPixelSignif = (MshInit->NbrVer-2*m) / incrementRegulier;
    }

    Mesh *MshDel = msh_init();
    MshDel->NbrVer = NbrPixelSignif + NbrPtFr;
    int SizTri = 4*MshDel->NbrVer;  //TODO Size control 
    MshDel->Crd = (double2d *)calloc(MshDel->NbrVer + 5, sizeof(double2d));
    MshDel->Tri = (int3d *)calloc(SizTri + 1, sizeof(int3d));
    MshDel->TriVoi = (int3d *)calloc(SizTri + 1, sizeof(int3d));

    /////////// Define Hash Table //////////
    HashTable *hshInit = msh_neighbors(MshInit, keyMode);
    HashTable *hshDel = hash_init(0.8*SizTri, 3*SizTri+1);     //TODO MaxObj Siz control


    /////////// Define solution //////////
    const char *solFile = "../data/joconde.sol";          //! Update here
    double *solInit = sol_read(solFile, 2, MshInit->NbrVer);
    double *solInterp = (double *)calloc(MshDel->NbrVer + 5, sizeof(double));
    double *solQC = (double *)calloc(MshInit->NbrVer, sizeof(double));

    /////////// Add frontier points ///////////
    
    FrPoints(MshDel, h_fr);
    solFr(solInit, solInterp, solQC);

    /////////// Create pixels significatifs ///////////
    double PSNR;
    debug_printf(" NbrVer of MshDel = %d\n", MshDel->NbrVer);
    createPixelSignif(MshInit, MshDel, NbrPixelSignif, solInit, NbrPtFr, solInterp, solQC);
    
    
    ///////////// Triangulation Delaunay ///////////
    debug_printf("NbrVer of MshDel : %d\n", MshDel->NbrVer);
    msh_boundingbox(MshDel);
    printf("[ Finish creating bounding box ]\n");
    
    TriangulationDelaunay(MshDel, hshDel, 1, MshDel->NbrVer);
    //save_test_points_to_file("../output/test_points.txt", MshDel);
    save_triangulation_to_file("../output/triangulation.mesh", MshDel, NbrPtFr);  
    debug_printf("NbrVer of MshDel = %d\n", MshDel->NbrVer);
    

    Mesh *MshTmp = msh_read("../output/triangulation.mesh", readEfr);
    HashTable *hshTmp;
    hshTmp = msh_neighbors(MshTmp, keyMode);
    
    int *color = find_connex_components_modif(MshTmp);
    
    Mesh *MshIntern = msh_init();
    selectInternalTri(MshIntern, MshDel, color);
    HashTable *hshIntern = msh_neighbors(MshIntern, keyMode);
    //write_TriVoi_to_file("../output/TriVoi_MshIntern.txt", MshIntern);
    save_triangulation_to_file("../output/triangulation.mesh", MshIntern, 0); 
    
    // /////////// Interpolation and PSNR ///////////
    
    msh_write2dfield_Vertices("../output/solCompr.sol", MshIntern->NbrVer+4, solInterp);
    interpolateSolution(MshInit, MshDel, solInit, solInterp, NbrPtFr);
    msh_write2dfield_Vertices("../output/solCompr.sol", MshIntern->NbrVer+4, solInterp);
    PSNR = computePSNR(MshInit, MshIntern, solInit, solInterp, solQC);
    printf("MshIntern->NbrVer = %d\n", MshIntern->NbrVer);
    printf("PSNR = %f (decibel)\n", PSNR);

    /////////// Free memory ///////////
    free(solInit);
    free(solInterp);
    free(solQC);
    msh_free(MshInit);
    msh_free(MshDel);
    msh_free(MshIntern);
    msh_free(MshTmp);
    hash_free(hshDel);
    hash_free(hshInit);
    hash_free(hshIntern);
    hash_free(hshTmp);
    return 0;
}
