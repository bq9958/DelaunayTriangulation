#include <error.h>

void ifOriPixel(double *solInit, double x, double y)
{
    
}

double computePSNR_modif(Mesh *MshInit, Mesh *MshIntern, double *solInit, double *solInterp, double *solQC)
{
    printf("[[ Begin computePSNR ]]\n");
    double x, y; 
    int move = 0;
    int iTri, iTriLoc; 
    double u, uc;
    double qc = 0;
    double beta0, beta1, beta2;
    double PSNR;
    srand(time(NULL));
    for (int i=1; i<=MshInit->NbrVer; i++)
    {
        if (solQC[i] == -1)
        {
            debug_printf(" i = %d\n", i);
            solQC[i] = 0;
            continue;
        }
        
        move = 0;
        u = solInit[i];
        x = MshInit->Crd[i][0]; y = MshInit->Crd[i][1];
        iTri = rand() % MshIntern->NbrTri;
        iTriLoc = location(MshIntern, iTri, x, y, &move);
        
        debug_printf("iTri = %d, iTriLoc = %d\n", iTri, iTriLoc);
        barycenter(MshIntern->Crd[MshIntern->Tri[iTriLoc][0]][0], MshIntern->Crd[MshIntern->Tri[iTriLoc][0]][1], 
                MshIntern->Crd[MshIntern->Tri[iTriLoc][1]][0], MshIntern->Crd[MshIntern->Tri[iTriLoc][1]][1], 
                MshIntern->Crd[MshIntern->Tri[iTriLoc][2]][0], MshIntern->Crd[MshIntern->Tri[iTriLoc][2]][1], 
                x, y, &beta0, &beta1, &beta2);
        debug_printf("beta0 = %f, beta1 = %f, beta2 = %f\n", beta0, beta1, beta2);
        uc = beta0 * solInterp[MshIntern->Tri[iTriLoc][0]] + beta1 * solInterp[MshIntern->Tri[iTriLoc][1]]\ 
                                                        + beta2 * solInterp[MshIntern->Tri[iTriLoc][2]];
        
        solQC[i] = (u - uc) * (u - uc) / MshInit->NbrVer;
        qc += (u - uc) * (u - uc) / MshInit->NbrVer;
    }
    PSNR = 10 * log10(255*255 / qc);
    msh_write2dfield_Vertices("../output/solQC.sol", MshInit->NbrVer, solQC);
    printf("[Output File] PSNR solution written in solQC.sol \n");

    return PSNR;
}