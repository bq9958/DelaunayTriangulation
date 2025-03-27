#include <comprImage.h>

extern const char *pixelSignifMode;
extern const int m, n;

double computePSNR(Mesh *MshInit, Mesh *MshDel, double *solInit, double *solInterp, double *solPSNR)
{
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
        move = 0;
        u = solInit[i];
        x = MshInit->Crd[i][0]; y = MshInit->Crd[i][1];
        iTri = rand() % MshDel->NbrTri + 4;     // avoid the bounding box
        iTriLoc = location(MshDel, iTri, x, y, &move);
        printf("iTri = %d, iTriLoc = %d\n", iTri, iTriLoc);
        barycenter(MshDel->Crd[MshDel->Tri[iTriLoc][0]][0], MshDel->Crd[MshDel->Tri[iTriLoc][0]][1], 
                   MshDel->Crd[MshDel->Tri[iTriLoc][1]][0], MshDel->Crd[MshDel->Tri[iTriLoc][1]][1], 
                   MshDel->Crd[MshDel->Tri[iTriLoc][2]][0], MshDel->Crd[MshDel->Tri[iTriLoc][2]][1], 
                   x, y, &beta0, &beta1, &beta2);
        printf("beta0 = %f, beta1 = %f, beta2 = %f\n", beta0, beta1, beta2);
        uc = beta0 * solInterp[MshDel->Tri[iTriLoc][0]] + beta1 * solInterp[MshDel->Tri[iTriLoc][1]]\ 
                                                        + beta2 * solInterp[MshDel->Tri[iTriLoc][2]];
        qc += (u - uc) * (u - uc) / MshInit->NbrVer;
        solPSNR[i] = qc;
    }
    PSNR = 10 * log10(255*255 / qc);
    msh_write2dfield_Vertices("../output/solPSNR.sol", MshInit->NbrVer, solPSNR);
    printf("[Output File] PSNR solution written in solPSNR.sol \n");

    return PSNR;
}


void interpolateSolution(Mesh *MshInit, Mesh *MshDel, double *solInit, double *solInterp)       
{
    
    int iTriLoc = 0;
    double beta0, beta1, beta2;

    if (strcmp(pixelSignifMode, "regulier") == 0)
    {
        
    }
    else if (strcmp(pixelSignifMode, "aleatoire") == 0 || strcmp(pixelSignifMode, "gradient") == 0)
    {
        printf(" -------- Begin Interpolation -------- \n");
        srand(time(NULL));
        for (int i=1; i<=MshDel->NbrVer;i++)
        {
            int move = 0;
            int iTri = rand() % MshInit->NbrTri ;
            printf("iTri = %d\n", iTri);
            iTriLoc = location(MshInit, iTri, MshDel->Crd[4+i][0], MshDel->Crd[4+i][1], &move);
            printf("iTriLoc = %d\n", iTriLoc);
            barycenter(MshInit->Crd[MshInit->Tri[iTriLoc][0]][0], MshInit->Crd[MshInit->Tri[iTriLoc][0]][1], 
                       MshInit->Crd[MshInit->Tri[iTriLoc][1]][0], MshInit->Crd[MshInit->Tri[iTriLoc][1]][1], 
                       MshInit->Crd[MshInit->Tri[iTriLoc][2]][0], MshInit->Crd[MshInit->Tri[iTriLoc][2]][1], 
                       MshDel->Crd[4+i][0], MshDel->Crd[4+i][1], &beta0, &beta1, &beta2);
            printf("beta0 = %f, beta1 = %f, beta2 = %f\n", beta0, beta1, beta2);
            solInterp[4+i] = beta0 * solInit[MshInit->Tri[iTriLoc][0]] + beta1 * solInit[MshInit->Tri[iTriLoc][1]]\ 
                                                                     + beta2 * solInit[MshInit->Tri[iTriLoc][2]];
        }
        printf("NbrVer of MshDel = %d\n", MshDel->NbrVer);
        msh_write2dfield_Vertices("../output/solCompr.sol", MshDel->NbrVer+4, solInterp);
        printf("[Output File] Interpolated solution written in solCompr.sol \n");
    }
}

// TriangulationDelaunay is a packaged version for image compression
void TriangulationDelaunay(Mesh *MshDel, HashTable *hshDel, 
    double box_xmin, double box_xmax, double box_ymin, double box_ymax,
    int NbrVerDel)
{
    //////// Constants ////////
    int iTriLocLast = 0;

    //////// Create the bounding box ////////
    MshDel->Crd[1][0] = box_xmin; MshDel->Crd[1][1] = box_ymax; // (xmin, ymax)
    MshDel->Crd[2][0] = box_xmin; MshDel->Crd[2][1] = box_ymin; // (xmin, ymin)
    MshDel->Crd[3][0] = box_xmax; MshDel->Crd[3][1] = box_ymin; // (xmax, ymin)
    MshDel->Crd[4][0] = box_xmax; MshDel->Crd[4][1] = box_ymax; // (xmax, ymax)

    MshDel->Tri[1][0] = 1; MshDel->Tri[1][1] = 2; MshDel->Tri[1][2] = 4;
    MshDel->Tri[2][0] = 2; MshDel->Tri[2][1] = 3; MshDel->Tri[2][2] = 4;
    MshDel->NbrTri += 2;
    MshDel->TriVoi[1][0] = 2; MshDel->TriVoi[1][1] = 0; MshDel->TriVoi[1][2] = 0;
    MshDel->TriVoi[2][0] = 1; MshDel->TriVoi[2][1] = 0; MshDel->TriVoi[2][2] = 0;
    printf("[ Finish creating the bounding box ]\n");

    for (int iPtIns = 1; iPtIns <= MshDel->NbrVer; iPtIns++)
    {
        printf("[ iPtIns = %d ]\n", iPtIns);
        Cavity(MshDel, hshDel, iPtIns, &iTriLocLast);
    }
}

void createPixelSignif(Mesh *MshInit, Mesh *MshDel, int NbrPixelSignif, double *solInit)
{    
    dynamArr *pixelSignif = NULL;

    if (NbrPixelSignif > MshInit->NbrVer)
    {
        printf("Error: NbrPixelSignif is larger than the number of vertices in the mesh!\n");
        return;
    }

    else if (strcmp(pixelSignifMode, "regulier") == 0)  //! Not robust
    {
        int increment = MshInit->NbrVer / NbrPixelSignif;
        for (int i=1; i<=NbrPixelSignif; i++)
        {
            MshDel->Crd[i][0] = MshInit->Crd[i*increment][0];
            MshDel->Crd[i][1] = MshInit->Crd[i*increment][1];
        }
    }
    else if (strcmp(pixelSignifMode, "aleatoire") == 0)
    {
        MshDel->Crd = TestPointSet(NbrPixelSignif, 1, m-1, 1, n-1);
    }
    else if (strcmp(pixelSignifMode, "gradient") == 0)  // 现在用作 block-wise random
    {
        srand(time(NULL));
        MshDel->NbrVer = 0;
        int countTot = 0;

        int Nbr_bloc_m = 3;   // number of vertical blocks (rows)
        int Nbr_bloc_n = 3;   // number of horizontal blocks (cols)
        int NbPerBlock = NbrPixelSignif / (Nbr_bloc_m * Nbr_bloc_n) + 1; 
        printf("--------- Begin Block-wise Random Selection ---------\n");

        for (int bloc_m = 1; bloc_m <= Nbr_bloc_m; bloc_m++)
        {
            for (int bloc_n = 1; bloc_n <= Nbr_bloc_n; bloc_n++)
            {
                int Mprev = m / Nbr_bloc_m * (bloc_m - 1);
                int Nprev = n / Nbr_bloc_n * (bloc_n - 1);
                int M = m / Nbr_bloc_m * bloc_m;
                int N = n / Nbr_bloc_n * bloc_n;

                int H = M - Mprev;  
                int W = N - Nprev; 

                for (int k = 0; k < NbPerBlock && countTot < NbrPixelSignif; k++)
                {
                    double x = (double)rand() / RAND_MAX ;
                    double y = (double)rand() / RAND_MAX ;
                    MshDel->NbrVer++;
                    countTot++;
                    MshDel->Crd[MshDel->NbrVer][0] = Mprev + 1 + x * (M - Mprev);
                    MshDel->Crd[MshDel->NbrVer][1] = Nprev + 1 + y * (N - Nprev);
                }     
                if (countTot >= NbrPixelSignif)
                break;           
            }
            if (countTot >= NbrPixelSignif)
            break;
        }
        printf("Total randomly selected points: %d\n", countTot);
    }
    else{ printf("Error: Unable to create pixelSignif"); }

}