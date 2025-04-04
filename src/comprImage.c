#include <comprImage.h>

extern const char *pixelSignifMode;
extern const int m, n;
extern int tri2edg[3][2];
extern int incrementRegulier;
#define EPSILON 1e-14


double computePSNR(Mesh *MshInit, Mesh *MshIntern, double *solInit, double *solInterp, double *solQC)
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
        iTri = rand() % MshIntern->NbrTri + 4;     // avoid the bounding box
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


void interpolateSolution(Mesh *MshInit, Mesh *Msh, double *solInit, double *solInterp, int NbrPtFr)       
{
    printf("[[ Begin interpolateSolution ]]\n");
    int iTriLoc = 0;
    double beta0, beta1, beta2;

    if (strcmp(pixelSignifMode, "regulier") == 0)
    {
        // implement in createPixelSignif
    }
    else if (strcmp(pixelSignifMode, "aleatoire") == 0 || strcmp(pixelSignifMode, "bloc") == 0)
    {
        debug_printf(" -------- Begin Interpolation -------- \n");
        srand(time(NULL));
        for (int i=NbrPtFr+1; i<=Msh->NbrVer;i++)
        {
            debug_printf("i = %d\n", i);
            debug_printf("MshDel->NbrVer = %d\n", Msh->NbrVer);
            int move = 0;
            int iTri = rand() % MshInit->NbrTri;
            debug_printf("iTri = %d\n", iTri);
            iTriLoc = location(MshInit, iTri, Msh->Crd[4+i][0], Msh->Crd[4+i][1], &move);
            debug_printf("iTriLoc = %d\n", iTriLoc);
            barycenter(MshInit->Crd[MshInit->Tri[iTriLoc][0]][0], MshInit->Crd[MshInit->Tri[iTriLoc][0]][1], 
                       MshInit->Crd[MshInit->Tri[iTriLoc][1]][0], MshInit->Crd[MshInit->Tri[iTriLoc][1]][1], 
                       MshInit->Crd[MshInit->Tri[iTriLoc][2]][0], MshInit->Crd[MshInit->Tri[iTriLoc][2]][1], 
                       Msh->Crd[4+i][0], Msh->Crd[4+i][1], &beta0, &beta1, &beta2);
            debug_printf("beta0 = %f, beta1 = %f, beta2 = %f\n", beta0, beta1, beta2);
            solInterp[4+i] = beta0 * solInit[MshInit->Tri[iTriLoc][0]] + beta1 * solInit[MshInit->Tri[iTriLoc][1]]\ 
                                                                     + beta2 * solInit[MshInit->Tri[iTriLoc][2]];
        }
        debug_printf("NbrVer of MshDel = %d\n", Msh->NbrVer);
        msh_write2dfield_Vertices("../output/solCompr.sol", Msh->NbrVer+4, solInterp);
        printf("[Output File] Interpolated solution written in solCompr.sol \n");
    }
}

void solFr(double *solInit, double *solInterp, double *solQC)      // only for h = 1
{
    solInterp[5] = solInit[1];     // upper-left corner
    solQC[1] = -1;

    for (int i=1; i<=m-2; i++)     // up
    {
        solInterp[5 + i] = solInit[1+i];
        solQC[1+i] = -1;
    }
    solInterp[5 + m - 1] = solInit[m];       // upper-right corner
    solQC[m] = -1;

    for (int i=1; i<=n-2; i++)      //right
    {
        solInterp[5 + m -1 + i] = solInit[(1+i)*m];
        solQC[(1+i)*m] = -1;
    }
    solInterp[5 + m - 1 + n - 1] = solInit[n*m];   // bottom-right corner
    solQC[n*m] = -1;

    for (int i=1; i<=m-2; i++)     // bottom
    {
        solInterp[5 + m - 1 + n - 1 + i] = solInit[n*m - i];
        solQC[n*m - i] = -1;
    }
    solInterp[5 + m-1 + n-1 + m-1] = solInit[n*m - (m-1)];   // bottom-left corner
    solQC[n*m - (m-1)] = -1;

    for (int i=1; i<=n-2; i++)      // left
    {
        solInterp[5 + m-1 + n-1 + m-1 + i] = solInit[n*m - (m-1) - i*m];
        solQC[n*m - (m-1) - i*m] = -1;
    }
}

// TriangulationDelaunay is a packaged version for image compression
void TriangulationDelaunay(Mesh *MshDel, HashTable *hshDel, int NbrPtInsMin, int NbrPtInsMax)
{
    printf("[[ Begin Triangulation ]]\n");
    //////// Constants ////////
    int iTriLocLast = 0;
    
    //////////////////////////
    if (NbrPtInsMin == 1)    // start from beginning
    {
        MshDel->Tri[1][0] = 1; MshDel->Tri[1][1] = 2; MshDel->Tri[1][2] = 4;
        MshDel->Tri[2][0] = 2; MshDel->Tri[2][1] = 3; MshDel->Tri[2][2] = 4;
        MshDel->NbrTri += 2;
        MshDel->TriVoi[1][0] = 2; MshDel->TriVoi[1][1] = 0; MshDel->TriVoi[1][2] = 0;
        MshDel->TriVoi[2][0] = 1; MshDel->TriVoi[2][1] = 0; MshDel->TriVoi[2][2] = 0;
    }
    // debug_printf("[ Finish creating the bounding box ]\n");
    for (int iPtIns = NbrPtInsMin; iPtIns <= NbrPtInsMax; iPtIns++)
    {
        debug_printf("[ iPtIns = %d ]\n", iPtIns);
        Cavity(MshDel, hshDel, iPtIns, &iTriLocLast);
    }
}

void createPixelSignif(Mesh *MshInit, Mesh *MshDel, int NbrPixelSignif, double *solInit, 
                       int NbPtFr, double *solInterp, double *solQC)
{    
    printf("[[ Begin createPixelSignif ]]\n");
    if (NbrPixelSignif > MshInit->NbrVer)
    {
        LOG_ERROR("NbrPixelSignif is larger than the number of vertices in the mesh!\n");
        return;
    }

    else if (strcmp(pixelSignifMode, "regulier") == 0)  
    {
        int Nbr = (MshInit->NbrVer-2*m) / incrementRegulier;      // exclude upper and lower side
        int count = 0;

        for (int i=1; i<=Nbr; i++)
        {
            debug_printf("i = %d\n", i);
            if ((i*incrementRegulier-1)%m > EPSILON && (i*incrementRegulier)%m > EPSILON)    // exclude left and right side
            {
                count ++;
                if (MshInit->Crd[m+i*incrementRegulier] == NULL)
                {
                    LOG_ERROR("MshInit->Crd[][0] is NULL!\n");
                    return;
                }
                MshDel->Crd[count+4+NbPtFr][0] = MshInit->Crd[m+i*incrementRegulier][0];
                MshDel->Crd[count+4+NbPtFr][1] = MshInit->Crd[m+i*incrementRegulier][1];
                solInterp[count+4+NbPtFr] = solInit[m+i*incrementRegulier];
                solQC[m+i*incrementRegulier] = -1;
            }
        }
        MshDel->NbrVer = count + NbPtFr;
        debug_printf("count = %d\n", count);

        double2d* tmp = realloc(MshDel->Crd, (count+NbPtFr+5) * sizeof(double2d));
        MshDel->Crd = tmp;
        if (MshDel->Crd == NULL)
        {
            LOG_ERROR("Memory allocation failed for MshDel->Crd!\n");
            return;
        }
        else{ printf("Msh->Crd resized to %d\n", count+NbPtFr+5); }
    }
    else if (strcmp(pixelSignifMode, "aleatoire") == 0)
    {
        double2d *Points = TestPointSet(NbrPixelSignif, 1, m-1, 1, n-1);
        for (int i=1; i<=NbrPixelSignif; i++)
        {
            debug_printf("i + NbPtFr = %d\n", i + NbPtFr);
            MshDel->Crd[i+4+NbPtFr][0] = Points[i][0];
            MshDel->Crd[i+4+NbPtFr][1] = Points[i][1];
        }
        free(Points);
    }
    else if (strcmp(pixelSignifMode, "bloc") == 0) 
    {
        srand(time(NULL));
        int position = 4;   // leave space for bounding box
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
                    position++;
                    countTot++;
                    MshDel->Crd[position + NbPtFr][0] = Mprev + x * (M - Mprev - 1 - 1e-5) + 1 + 1e-10;   // avoid exceeding boundary
                    MshDel->Crd[position + NbPtFr][1] = Nprev + y * (N - Nprev - 1 - 1e-5) + 1 + 1e-10;   // avoid exceeding boundary
                }     
                if (countTot >= NbrPixelSignif)
                break;           
            }
            if (countTot >= NbrPixelSignif)
            break;
        }
        printf("Total randomly selected points: %d\n", position);
    }
    else{ LOG_ERROR("Unable to create pixelSignif"); }
}

void FrPoints(Mesh *MshDel, double h)
{
    printf("[[ Add frontier points ]]\n");
    int nx = (m-1) / h - 1;
    int ny = (n-1) / h - 1;
    int count = 5;     // leave space for bounding box 

    MshDel->Crd[count][0] = 1; MshDel->Crd[count][1] = n;    // upper left
    count ++;
    for (int i=1; i<=nx; i++)   // up
    {
        MshDel->Crd[count][0] = 1 + i*h; MshDel->Crd[count][1] = n; 
        count ++;
    }
    
    MshDel->Crd[count][0] = m; MshDel->Crd[count][1] = n;    // upper right
    count ++;
    for (int i=1; i<=ny; i++)   // right
    {
        MshDel->Crd[count][0] = m; MshDel->Crd[count][1] = n - i*h; 
        count ++;
    }
    
    MshDel->Crd[count][0] = m; MshDel->Crd[count][1] = 1;    // bottom right
    count ++;
    for (int i=1; i<=nx; i++)   // bottom
    {
        MshDel->Crd[count][0] = m - i*h; MshDel->Crd[count][1] = 1; 
        count ++;
    }
    
    MshDel->Crd[count][0] = 1; MshDel->Crd[count][1] = 1;    // bottom left
    count ++;
    for (int i=1; i<=ny; i++)   // left
    {
        MshDel->Crd[count][0] = 1; MshDel->Crd[count][1] = 1 + i*h; 
        count ++;
    }
    
}

void selectInternalTri(Mesh *MshIntern, Mesh *MshDel, int *color)
{
    MshIntern->NbrVer = MshDel->NbrVer;
    int SizTri = 4*MshIntern->NbrVer;
    MshIntern->Crd = (double2d *)calloc(MshIntern->NbrVer + 5, sizeof(double2d));
    MshIntern->Tri = (int3d *)calloc(SizTri + 1, sizeof(int3d));
    MshIntern->TriVoi = (int3d *)calloc(SizTri + 1, sizeof(int3d));    

    for (int i=1; i<=4+MshIntern->NbrVer; i++)    // copy Crd
    {
        MshIntern->Crd[i][0] = MshDel->Crd[i][0];
        MshIntern->Crd[i][1] = MshDel->Crd[i][1];
    }

    int indicator;
    if ((MshDel->Crd[MshDel->Tri[1][0]][0] < EPSILON || MshDel->Crd[MshDel->Tri[1][0]][0] > m + EPSILON) || 
        (MshDel->Crd[MshDel->Tri[1][1]][0] < EPSILON || MshDel->Crd[MshDel->Tri[1][1]][0] > m + EPSILON) || 
        (MshDel->Crd[MshDel->Tri[1][2]][0] < EPSILON || MshDel->Crd[MshDel->Tri[1][2]][0] > m + EPSILON))     // bounding box triangular
    {
        indicator = color[1];    // the color of triangular to delete
    }

    for (int i=1; i<=MshDel->NbrTri; i++)
    {
        if (color[i] != indicator)
        {
            MshIntern->NbrTri ++;
            MshIntern->Tri[MshIntern->NbrTri][0] = MshDel->Tri[i][0];
            MshIntern->Tri[MshIntern->NbrTri][1] = MshDel->Tri[i][1];
            MshIntern->Tri[MshIntern->NbrTri][2] = MshDel->Tri[i][2];
        }
    }
    if (MshIntern->NbrTri == 0)
    {
        LOG_ERROR("No internal triangles found!\n");
        return;
    }
    else
    {
        printf("NbrTri of MshIntern = %d\n", MshIntern->NbrTri);
    }
}

int *find_connex_components_modif(Mesh *Msh)
{
  int *color = (int *)calloc(Msh->NbrTri+1, sizeof(int));
  int *stack = (int *)malloc((Msh->NbrTri+1) * sizeof(int));
  int ndomn = 0;
  int top;

  for (int iTri=1; iTri <= Msh->NbrTri; iTri++){
    if (color[iTri] == 0){
      ndomn++;
      color[iTri] = ndomn;

      top = 0;
      stack[top++] = iTri; 

      while (top > 0){
        int topTri = stack[--top];

        for (int iEdglocal=0; iEdglocal < 3; iEdglocal++){
          int neighbor = Msh->TriVoi[topTri][iEdglocal];
          if (neighbor != 0 && color[neighbor] == 0){
            color[neighbor] = ndomn;
            stack[top++] = neighbor;
          }
        }
      }
    }
  }

  int *NbrConnex = (int *)calloc(ndomn+1, sizeof(int));
  for (int iTri=1; iTri <= Msh->NbrTri; iTri++){
    NbrConnex[color[iTri]]++;
  }

  double *colorDouble = convertIntToDouble(color, Msh->NbrTri);
  //write_color_to_txt("../output/color.txt", color, Msh->NbrTri);

  for (int i=1; i<=ndomn; i++){
    printf("  Connex component %d has %d triangles\n", i, NbrConnex[i]);
  }

  //msh_write2dfield_Triangles("../output/connex.solb", Msh->NbrTri, colorDouble);
  printf("[Output File] connex components written in connex.solb \n");

  free(stack);
  free(NbrConnex);

  return color;
}

HashTable * msh_neighbors_modif(Mesh *Msh, const char* keyMode)
{
  int iTri, iEdg, iVer1, iVer2;
  
  
  if ( ! Msh ) return 0;
  
  if ( Msh->TriVoi == NULL )
    Msh->TriVoi = calloc( (Msh->NbrTri+1), sizeof(int3d) );  
  
  //--- initialize HashTable and set the hash table 
  int NbrMaxObj = 3 * Msh->NbrTri + 1; int SizHead = 0.8 * NbrMaxObj;
  HashTable *hsh = hash_init(SizHead, NbrMaxObj);
  
  //--- Compute the neighbors using the hash table
  for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
    for (iEdg=0; iEdg<3; iEdg++) {
      int boolEfr = 0;
      iVer1 = Msh->Tri[iTri][ tri2edg[iEdg][0] ];
      iVer2 = Msh->Tri[iTri][ tri2edg[iEdg][1] ];

      if (Msh->Efr != NULL){
        for (int iEfr=1; iEfr<=Msh->NbrEfr; iEfr++){
          if ((Msh->Efr[iEfr][0] == iVer1 && Msh->Efr[iEfr][1] == iVer2) \
           || (Msh->Efr[iEfr][0] == iVer2 && Msh->Efr[iEfr][1] == iVer1)){
            Msh->TriVoi[iTri][iEdg] = 0;
            boolEfr = 1;
          }
        }
      }

      if (boolEfr == 1) continue;

      // compute the key   
      // do we have objects as that key   hash_find () */
      //  if yes ===> look among objects and potentially update TriVoi */
      //  if no  ===> add to hash table   hash_add()   */
      int id = hash_find(hsh, iVer1, iVer2, keyMode);
      if (id != 0) {
        Msh->TriVoi[hsh->LstObj[id][2]][hsh->LstObj[id][3]] = iTri;
        hsh->LstObj[id][3] = iTri;
        Msh->TriVoi[iTri][iEdg] = hsh->LstObj[id][2];
      }
      else {
        hash_add(hsh, iVer1, iVer2, iTri, iEdg, keyMode);
        hash_add(hsh, iVer1, iVer2, iTri, iEdg, keyMode);
      }
    }
  }

  printf("  Nbr edges hash %10d \n", hsh->NbrObj);
  
  //write_Head_to_file("../output/Head.txt", hsh);
  //write_LstObj_to_file("../output/LstObj.txt", hsh);

  return hsh;
}  