#include <comprImage.h>

extern const char *pixelSignifMode;
extern const int m, n;

double computePSNR(Mesh *MshInit, Mesh *MshDel, double *solInit, double *solInterp, double *solPSNR)
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
        move = 0;
        u = solInit[i];
        x = MshInit->Crd[i][0]; y = MshInit->Crd[i][1];
        iTri = rand() % MshDel->NbrTri + 4;     // avoid the bounding box
        iTriLoc = location(MshDel, iTri, x, y, &move);
        debug_printf("iTri = %d, iTriLoc = %d\n", iTri, iTriLoc);
        barycenter(MshDel->Crd[MshDel->Tri[iTriLoc][0]][0], MshDel->Crd[MshDel->Tri[iTriLoc][0]][1], 
                   MshDel->Crd[MshDel->Tri[iTriLoc][1]][0], MshDel->Crd[MshDel->Tri[iTriLoc][1]][1], 
                   MshDel->Crd[MshDel->Tri[iTriLoc][2]][0], MshDel->Crd[MshDel->Tri[iTriLoc][2]][1], 
                   x, y, &beta0, &beta1, &beta2);
        debug_printf("beta0 = %f, beta1 = %f, beta2 = %f\n", beta0, beta1, beta2);
        uc = beta0 * solInterp[MshDel->Tri[iTriLoc][0]] + beta1 * solInterp[MshDel->Tri[iTriLoc][1]]\ 
                                                        + beta2 * solInterp[MshDel->Tri[iTriLoc][2]];
        
        solPSNR[i] = (u - uc) * (u - uc) / MshInit->NbrVer;
        qc += (u - uc) * (u - uc) / MshInit->NbrVer;
    }
    PSNR = 10 * log10(255*255 / qc);
    msh_write2dfield_Vertices("../output/solPSNR.sol", MshInit->NbrVer, solPSNR);
    printf("[Output File] PSNR solution written in solPSNR.sol \n");

    return PSNR;
}


void interpolateSolution(Mesh *MshInit, Mesh *MshDel, double *solInit, double *solInterp)       
{
    printf("[[ Begin interpolateSolution ]]\n");
    int iTriLoc = 0;
    double beta0, beta1, beta2;

    if (strcmp(pixelSignifMode, "regulier") == 0)
    {
        
    }
    else if (strcmp(pixelSignifMode, "aleatoire") == 0 || strcmp(pixelSignifMode, "bloc") == 0)
    {
        debug_printf(" -------- Begin Interpolation -------- \n");
        srand(time(NULL));
        for (int i=1; i<=MshDel->NbrVer;i++)
        {
            debug_printf("i = %d\n", i);
            debug_printf("MshDel->NbrVer = %d\n", MshDel->NbrVer);
            int move = 0;
            int iTri = rand() % MshInit->NbrTri;
            debug_printf("iTri = %d\n", iTri);
            iTriLoc = location(MshInit, iTri, MshDel->Crd[4+i][0], MshDel->Crd[4+i][1], &move);
            debug_printf("iTriLoc = %d\n", iTriLoc);
            barycenter(MshInit->Crd[MshInit->Tri[iTriLoc][0]][0], MshInit->Crd[MshInit->Tri[iTriLoc][0]][1], 
                       MshInit->Crd[MshInit->Tri[iTriLoc][1]][0], MshInit->Crd[MshInit->Tri[iTriLoc][1]][1], 
                       MshInit->Crd[MshInit->Tri[iTriLoc][2]][0], MshInit->Crd[MshInit->Tri[iTriLoc][2]][1], 
                       MshDel->Crd[4+i][0], MshDel->Crd[4+i][1], &beta0, &beta1, &beta2);
            debug_printf("beta0 = %f, beta1 = %f, beta2 = %f\n", beta0, beta1, beta2);
            solInterp[4+i] = beta0 * solInit[MshInit->Tri[iTriLoc][0]] + beta1 * solInit[MshInit->Tri[iTriLoc][1]]\ 
                                                                     + beta2 * solInit[MshInit->Tri[iTriLoc][2]];
        }
        debug_printf("NbrVer of MshDel = %d\n", MshDel->NbrVer);
        msh_write2dfield_Vertices("../output/solCompr.sol", MshDel->NbrVer+4, solInterp);
        printf("[Output File] Interpolated solution written in solCompr.sol \n");
    }
}

// TriangulationDelaunay is a packaged version for image compression
void TriangulationDelaunay(Mesh *MshDel, HashTable *hshDel, 
    double box_xmin, double box_xmax, double box_ymin, double box_ymax,
    const int NbrVerDel)
{
    printf("[[ Begin Triangulation ]]\n");
    //////// Constants ////////
    int iTriLocLast = 0;

    //////// Create the bounding box ////////
    for (int i=MshDel->NbrVer; i>=1; i--)
    {
        MshDel->Crd[i+4][0] = MshDel->Crd[i][0];
        MshDel->Crd[i+4][1] = MshDel->Crd[i][1];
        MshDel->Crd[i][0] = 0;
        MshDel->Crd[i][1] = 0;
    }
    save_triangulation_to_file("../output/triangulation.mesh", MshDel, hshDel);
    MshDel->Crd[1][0] = box_xmin; MshDel->Crd[1][1] = box_ymax; // (xmin, ymax)
    MshDel->Crd[2][0] = box_xmin; MshDel->Crd[2][1] = box_ymin; // (xmin, ymin)
    MshDel->Crd[3][0] = box_xmax; MshDel->Crd[3][1] = box_ymin; // (xmax, ymin)
    MshDel->Crd[4][0] = box_xmax; MshDel->Crd[4][1] = box_ymax; // (xmax, ymax)

    MshDel->Tri[1][0] = 1; MshDel->Tri[1][1] = 2; MshDel->Tri[1][2] = 4;
    MshDel->Tri[2][0] = 2; MshDel->Tri[2][1] = 3; MshDel->Tri[2][2] = 4;
    MshDel->NbrTri += 2;
    MshDel->TriVoi[1][0] = 2; MshDel->TriVoi[1][1] = 0; MshDel->TriVoi[1][2] = 0;
    MshDel->TriVoi[2][0] = 1; MshDel->TriVoi[2][1] = 0; MshDel->TriVoi[2][2] = 0;
    debug_printf("[ Finish creating the bounding box ]\n");
    for (int iPtIns = 1; iPtIns <= MshDel->NbrVer; iPtIns++)
    {
        debug_printf("[ iPtIns = %d ]\n", iPtIns);
        Cavity(MshDel, hshDel, iPtIns, &iTriLocLast);
    }
}

void createPixelSignif(Mesh *MshInit, Mesh *MshDel, int NbrPixelSignif, double *solInit)
{    
    printf("[[ Begin createPixelSignif ]]\n");
    if (NbrPixelSignif > MshInit->NbrVer)
    {
        LOG_ERROR("NbrPixelSignif is larger than the number of vertices in the mesh!\n");
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
    else if (strcmp(pixelSignifMode, "bloc") == 0) 
    {
        srand(time(NULL));
        MshDel->NbrVer = 0;
        int countTot = 0;

        int Nbr_bloc_m = 3;   // number of vertical blocks (rows)
        int Nbr_bloc_n = 3;   // number of horizontal blocks (cols)
        int NbPerBlockMoy = NbrPixelSignif / (Nbr_bloc_m * Nbr_bloc_n) + 1; 
        int NbPerBlock[3][3] = {{(int)(NbPerBlockMoy*0.2 + 1), (int)(NbPerBlockMoy*0.1 + 1), (int)(NbPerBlockMoy*0.1 + 1)}, 
                                {(int)(NbPerBlockMoy*3.1 + 1), (int)(NbPerBlockMoy*2.6 + 1), (int)(NbPerBlockMoy*2.5 + 1)}, 
                                {(int)(NbPerBlockMoy*0.2 + 1), (int)(NbPerBlockMoy*0.1 + 1), (int)(NbPerBlockMoy*0.1 + 1)}};
        debug_printf("--------- Begin Block-wise Random Selection ---------\n");

        for (int bloc_m = 1; bloc_m <= Nbr_bloc_m; bloc_m++)
        {
            for (int bloc_n = 1; bloc_n <= Nbr_bloc_n; bloc_n++)
            {
                int Mprev = m / Nbr_bloc_m * (bloc_m - 1);
                int Nprev = n / Nbr_bloc_n * (bloc_n - 1);
                int M = m / Nbr_bloc_m * bloc_m;
                int N = n / Nbr_bloc_n * bloc_n;

                for (int k = 0; k < NbPerBlock[bloc_m-1][bloc_n-1] && countTot < NbrPixelSignif; k++)
                {
                    double x = (double)rand() / RAND_MAX ;
                    double y = (double)rand() / RAND_MAX ;
                    MshDel->NbrVer++;
                    countTot++;
                    MshDel->Crd[MshDel->NbrVer][0] = Mprev + x * (M - Mprev - 1 - 1e-5) + 1 + 1e-10;   // avoid exceeding boundary
                    MshDel->Crd[MshDel->NbrVer][1] = Nprev + y * (N - Nprev - 1 - 1e-5) + 1 + 1e-10;   // avoid exceeding boundary
                }     
                if (countTot >= NbrPixelSignif)
                break;           
            }
            if (countTot >= NbrPixelSignif)
            break;
        }
        printf("Total randomly selected points: %d\n", MshDel->NbrVer);
    }
    else{ LOG_ERROR("Unable to create pixelSignif"); }
}