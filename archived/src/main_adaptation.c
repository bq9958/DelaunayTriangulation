#include <adaptation.h>

const char *keyMode = "sum";
int SizPil = 100; 
double Ntarget = 8177.505066;

int main(int argc, char *argv[])
{
    ///////// Read mesh //////////
    const char *meshFile = "../data/joconde.lowres.mesh";
    const char *solFile = "../data/joconde.lowres.sol";
    if (argc >= 2) meshFile = argv[1];
    if (argc >= 3) solFile  = argv[2];

    printf("[ Read mesh ]\n");
    int readEfr = 1;
    Mesh * Msh = msh_read(meshFile, readEfr);
    HashTable *hsh = msh_neighbors(Msh, keyMode);
    //write_TriVoi_to_file("../output/TriVoi.txt", Msh);
    double *uh = sol_read(solFile, 2, Msh->NbrVer);


    ///////// Create base list ////////
    dynamArr *ArrBase = listBase(Msh);
    //dyArr_export_to_file(ArrBase, "../output/Base.txt");

    ///////// Calculate Hessien //////////
    double3d *Hessien = (double3d *)malloc(sizeof(double3d) * (Msh->NbrVer+1));
    double2d *NablaRUh = (double2d *)malloc(sizeof(double2d) * (Msh->NbrVer+1));
    double *TriArea = (double *)malloc(sizeof(double) * (Msh->NbrTri+1));
    hessien(Msh, uh, ArrBase, Hessien, NablaRUh, TriArea);
    //export_NablaRUh_to_file("../output/NablaRUh.txt", Msh, NablaRUh);
    //export_Hessien_to_file("../output/Hessien.txt", Msh, Hessien);

    //////// Rebuild Hessien //////////
    rebuildHessien(Msh, Hessien);
    //export_Hessien_to_file("../output/HessienRebuild.txt", Msh, Hessien);

    //////// Calculate metric /////////
    int normalization = 1;
    if (strcmp(argv[1], "maillage.mesh") == 0) { normalization = 0; printf("Global normalization unenabled\n"); }
    else { printf("Global normalization enabled\n"); }
    double3d *Metric = (double3d *)malloc(sizeof(double3d) * (Msh->NbrVer+1));
    metric(Hessien, Msh, Metric, TriArea, normalization);
    //export_Hessien_to_file("../output/Metric.txt", Msh, Metric);
    printf("Msh NbrVer = %d\n", Msh->NbrVer);
    msh_write2dmetric("../output/maillage.met.sol", Msh->NbrVer, Metric);
    
    
    ///////// Test session /////////
    // int iPt = 2634; int iTri = ArrBase->data1d[iPt]; int NbrBoucle = 0;
    // dynamArr *ArrBoucle = boucleDetection(Msh, iPt, iTri, &NbrBoucle);
    // dyArr_print(ArrBoucle);
    // dyArr_free(ArrBoucle);

    ///////// Free memory /////////
    free(uh);
    free(TriArea);
    free(Hessien);
    free(Metric);
    free(NablaRUh);
    dyArr_free(ArrBase);
    msh_free(Msh);
    hash_free(hsh);
}