#include <adaptation.h>

const char *pixelSignifMode = "regulier";   // "regulier" ou "aleatoire" ou "bloc"
int SizPil = 100;       // cavity Tri Number max  //TODO Size Control
const char *keyMode = "sum";
const int m = 58;      // pixel width (58/289)    //! Update domain
const int n = 88;      // pixel height (88/440)  //! Update domain
const int incrementRegulier = 4;   // increment for regular pixel significatif
double Ntarget = 8177.505066;

int main(int argc, char *argv[])
{
    ////////// Read mesh //////////
    const char *FileMeshInit = "../data/joconde.lowres.mesh";
    const char *FileSolInit = "../data/joconde.lowres.sol";
    const char *FileMeshAdapte = "../dir/maillage.adapte.mesh";
    const char *FileSolAdapte = "../dir/maillage.niveaugris.itp.sol";
    if (argc >= 2) FileMeshInit = argv[1];
    if (argc >= 3) FileSolInit = argv[2];
    if (argc >= 4) FileMeshAdapte = argv[3];
    if (argc >= 5) FileSolAdapte = argv[4];

    Mesh *MshInit = msh_read(FileMeshInit, 1);
    Mesh *MshAdapte = msh_read(FileMeshAdapte, 1);
    double *solInit = sol_read(FileSolInit, 2, MshInit->NbrVer);
    double *solAdapte = sol_read(FileSolAdapte, 2, MshAdapte->NbrVer);
    HashTable *hshAdapte = msh_neighbors(MshAdapte, keyMode);
    
    ////////// Calculate error //////////
    double *solQC = (double *)calloc(MshInit->NbrVer, sizeof(double));
    double PSNR = computePSNR(MshInit, MshAdapte, solInit, solAdapte, solQC);
    printf("PSNR = %.6f\n", PSNR);

    ////////// Free memory //////////
    free(solInit);
    free(solAdapte);
    msh_free(MshInit);
    msh_free(MshAdapte);
    free(solQC);
    return 0;
}


