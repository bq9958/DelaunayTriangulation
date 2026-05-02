#include <adaptation.h>

extern int SizPil;
extern double Ntarget;

void metric(double3d *Hessien, Mesh *Msh, double3d *Metric, double *TriArea, int normalization)
{
    ////////// Calculate metric ///////////
    double *sqrtDetM = (double *)malloc(sizeof(double) * (Msh->NbrVer+1));
    for (int iPt=1; iPt<=Msh->NbrVer; iPt++)
    {
        double det = Hessien[iPt][0] * Hessien[iPt][2] - Hessien[iPt][1] * Hessien[iPt][1];
        Metric[iPt][0] = pow(det, -1.0/6.0) * Hessien[iPt][0];
        Metric[iPt][1] = pow(det, -1.0/6.0) * Hessien[iPt][1];
        Metric[iPt][2] = pow(det, -1.0/6.0) * Hessien[iPt][2];
        sqrtDetM[iPt] = sqrt(Metric[iPt][0] * Metric[iPt][2] - Metric[iPt][1] * Metric[iPt][1]);
    }

    ////////// Calculate complexite du metric ////////////
    double C = 0;
    for (int iTri=1; iTri<=Msh->NbrTri; iTri++)
    {
        double tmp = 0;
        tmp += sqrtDetM[Msh->Tri[iTri][0]];
        tmp += sqrtDetM[Msh->Tri[iTri][1]];
        tmp += sqrtDetM[Msh->Tri[iTri][2]];
        C += TriArea[iTri]/3 * tmp;
    }
    printf("C = %lf\n", C);
    printf("Ntarget = %lf\n", Ntarget);
    /////////// Normalisation du metric //////////
    if (normalization)
    {
        for (int i=1; i<=Msh->NbrVer; i++)
        {
            Metric[i][0] = Metric[i][0] * Ntarget / C;
            Metric[i][1] = Metric[i][1] * Ntarget / C;
            Metric[i][2] = Metric[i][2] * Ntarget / C;
        }
    }
    
    /////////// Free memory //////////
    free(sqrtDetM); 
}


void rebuildHessien(Mesh *Msh, double3d *Hessien)
{
    double met[3];
    double eigVec[4]; double eigVal[2];
    for (int i=1; i<=Msh->NbrVer; i++)
    {
        met[0] = Hessien[i][0];
        met[1] = Hessien[i][1];
        met[2] = Hessien[i][2];
        Sol_Eigen2d(met, eigVal, eigVec);
        double v1 = fmax(fabs(eigVal[0]), 1e-10);
        double v2 = fmax(fabs(eigVal[1]), 1e-10);
        Hessien[i][0] = v1 * eigVec[0] * eigVec[0] + v2 * eigVec[2] * eigVec[2];
        Hessien[i][1] = v1 * eigVec[0] * eigVec[1] + v2 * eigVec[2] * eigVec[3];
        Hessien[i][2] = v1 * eigVec[1] * eigVec[1] + v2 * eigVec[3] * eigVec[3];
    }
}


void hessien(Mesh *Msh, double *uh, dynamArr *ArrBase, 
             double3d *Hessien, double2d *NablaRUh, double *TriArea)
{   
    ///////// Calculate TriArea //////////
    printf("Begin calculating TriArea\n");
    for (int iTri=1; iTri<=Msh->NbrTri; iTri++)
    {
        double x0 = Msh->Crd[Msh->Tri[iTri][0]][0];
        double y0 = Msh->Crd[Msh->Tri[iTri][0]][1];
        double x1 = Msh->Crd[Msh->Tri[iTri][1]][0];
        double y1 = Msh->Crd[Msh->Tri[iTri][1]][1];
        double x2 = Msh->Crd[Msh->Tri[iTri][2]][0];
        double y2 = Msh->Crd[Msh->Tri[iTri][2]][1];
        TriArea[iTri] = triArea(x0, y0, x1, y1, x2, y2);
    }

    ///////// Calculate TriArea and nablaR_uh //////////
    printf("Begin calculating nablaR_uh\n");
    for (int iPt=1; iPt<=Msh->NbrVer; iPt++)
    {
        int iTri = ArrBase->data1d[iPt]; 
        int NbrBoucle = 0;
        dynamArr *ArrBoucle = boucleDetection(Msh, iPt, iTri, &NbrBoucle);
        double *ListNablaRUhX = (double *)malloc(sizeof(double) * NbrBoucle);
        double *ListNablaRUhY = (double *)malloc(sizeof(double) * NbrBoucle);
        double nabla_uhx, nabla_uhy;
        double nablaR_uhx, nablaR_uhy;

        for (int iBoucle=0; iBoucle<NbrBoucle; iBoucle++)
        {
            int iTriBoucle = ArrBoucle->data1d[iBoucle];
            nabla_uh(iTriBoucle, Msh, uh, &nabla_uhx, &nabla_uhy);
            ListNablaRUhX[iBoucle] = nabla_uhx;
            ListNablaRUhY[iBoucle] = nabla_uhy;
        }
        nablaR_uhx = L2Projector(TriArea, ArrBoucle, ListNablaRUhX, NbrBoucle);
        nablaR_uhy = L2Projector(TriArea, ArrBoucle, ListNablaRUhY, NbrBoucle);
        
        // update NablaRUh
        NablaRUh[iPt][0] = nablaR_uhx; 
        NablaRUh[iPt][1] = nablaR_uhy;

        free(ListNablaRUhX);
        free(ListNablaRUhY);
        dyArr_free(ArrBoucle);
    }

    ///////// Calculate nabla2R_uh //////////
    printf("Begin calculating nabla2R_uh\n");
    for (int iPt=1; iPt<=Msh->NbrVer; iPt++)
    {
        int iTri = ArrBase->data1d[iPt]; 
        int NbrBoucle = 0;
        dynamArr *ArrBoucle = boucleDetection(Msh, iPt, iTri, &NbrBoucle);
        double *ListNabla2RUhXX = (double *)malloc(sizeof(double) * NbrBoucle);
        double *ListNabla2RUhXY = (double *)malloc(sizeof(double) * NbrBoucle);
        double *ListNabla2RUhYY = (double *)malloc(sizeof(double) * NbrBoucle);

        double nablaxx_uh, nablaxy_uh, nablayy_uh;
        double nablaRxx_uh, nablaRxy_uh, nablaRyy_uh;

        for (int iBoucle=0; iBoucle<NbrBoucle; iBoucle++)
        {
            int iTriBoucle = ArrBoucle->data1d[iBoucle];
            nabla2_uh(iTriBoucle, Msh, NablaRUh, 
                      &nablaxx_uh, &nablaxy_uh, &nablayy_uh);
            ListNabla2RUhXX[iBoucle] = nablaxx_uh;
            ListNabla2RUhXY[iBoucle] = nablaxy_uh;
            ListNabla2RUhYY[iBoucle] = nablayy_uh;
        }
        nablaRxx_uh = L2Projector(TriArea, ArrBoucle, ListNabla2RUhXX, NbrBoucle);
        nablaRxy_uh = L2Projector(TriArea, ArrBoucle, ListNabla2RUhXY, NbrBoucle);
        nablaRyy_uh = L2Projector(TriArea, ArrBoucle, ListNabla2RUhYY, NbrBoucle);

        // update Hessien 
        Hessien[iPt][0] = nablaRxx_uh; 
        Hessien[iPt][1] = nablaRxy_uh; 
        Hessien[iPt][2] = nablaRyy_uh;

        free(ListNabla2RUhXX);
        free(ListNabla2RUhXY);
        free(ListNabla2RUhYY);
        dyArr_free(ArrBoucle);
    }
}


double L2Projector(double *TriArea, dynamArr *ArrBoucle, double *ListQOI, int NbrBoucle)
{
    double L2 = 0; double sumTriArea = 0;
    for (int i=0; i<NbrBoucle; i++)
    {
        L2 += TriArea[ArrBoucle->data1d[i]] * ListQOI[i];
        sumTriArea += TriArea[ArrBoucle->data1d[i]];
    }
    L2 = L2 / sumTriArea;

    return L2;
}

void nabla2_uh(int iTri, Mesh *Msh, double2d *NablaRUh, 
               double *nablaxx_uh, double *nablaxy_uh, double *nablayy_uh)
{
    // output : matrix 2*2 : nabla2_uh = (nablaxx_uh, nablayx_uh; nablayx_uh, nablayy_uh)
    double x0 = Msh->Crd[Msh->Tri[iTri][0]][0];
    double y0 = Msh->Crd[Msh->Tri[iTri][0]][1];
    double x1 = Msh->Crd[Msh->Tri[iTri][1]][0];
    double y1 = Msh->Crd[Msh->Tri[iTri][1]][1];
    double x2 = Msh->Crd[Msh->Tri[iTri][2]][0];
    double y2 = Msh->Crd[Msh->Tri[iTri][2]][1];
    double area = triArea(x0, y0, x1, y1, x2, y2);

    double nabla_phi0x = (y1 - y2) / (2 * area);
    double nabla_phi0y = (x2 - x1) / (2 * area);
    double nabla_phi1x = (y2 - y0) / (2 * area);
    double nabla_phi1y = (x0 - x2) / (2 * area);
    double nabla_phi2x = (y0 - y1) / (2 * area);
    double nabla_phi2y = (x1 - x0) / (2 * area);
    double nablaR_uhp0x = NablaRUh[Msh->Tri[iTri][0]][0];
    double nablaR_uhp0y = NablaRUh[Msh->Tri[iTri][0]][1];
    double nablaR_uhp1x = NablaRUh[Msh->Tri[iTri][1]][0];
    double nablaR_uhp1y = NablaRUh[Msh->Tri[iTri][1]][1];
    double nablaR_uhp2x = NablaRUh[Msh->Tri[iTri][2]][0];
    double nablaR_uhp2y = NablaRUh[Msh->Tri[iTri][2]][1];

    *nablaxx_uh = nabla_phi0x * nablaR_uhp0x + nabla_phi1x * nablaR_uhp1x + nabla_phi2x * nablaR_uhp2x;
    *nablaxy_uh = nabla_phi0y * nablaR_uhp0x + nabla_phi1y * nablaR_uhp1x + nabla_phi2y * nablaR_uhp2x;
    *nablayy_uh = nabla_phi0y * nablaR_uhp0y + nabla_phi1y * nablaR_uhp1y + nabla_phi2y * nablaR_uhp2y;
}

void nabla_uh(int iTri, Mesh *Msh, double *uh, double *nabla_uhx, double *nabla_uhy)
{
    // output : vector nabla_uh = (nabla_uhx, nabla_uhy)
    double x0 = Msh->Crd[Msh->Tri[iTri][0]][0];
    double y0 = Msh->Crd[Msh->Tri[iTri][0]][1];
    double x1 = Msh->Crd[Msh->Tri[iTri][1]][0];
    double y1 = Msh->Crd[Msh->Tri[iTri][1]][1];
    double x2 = Msh->Crd[Msh->Tri[iTri][2]][0];
    double y2 = Msh->Crd[Msh->Tri[iTri][2]][1];

    double area = triArea(x0, y0, x1, y1, x2, y2);
    double nabla_phi0x = (y1 - y2) / (2 * area);
    double nabla_phi0y = (x2 - x1) / (2 * area);
    double nabla_phi1x = (y2 - y0) / (2 * area);
    double nabla_phi1y = (x0 - x2) / (2 * area);
    double nabla_phi2x = (y0 - y1) / (2 * area);
    double nabla_phi2y = (x1 - x0) / (2 * area);

    double uhp0 = uh[Msh->Tri[iTri][0]];
    double uhp1 = uh[Msh->Tri[iTri][1]];
    double uhp2 = uh[Msh->Tri[iTri][2]];

    *nabla_uhx = nabla_phi0x * uhp0 + nabla_phi1x * uhp1 + nabla_phi2x * uhp2;
    *nabla_uhy = nabla_phi0y * uhp0 + nabla_phi1y * uhp1 + nabla_phi2y * uhp2;
}



dynamArr *listBase(Mesh *Msh)
{
    printf("Begin listBase\n");
    dynamArr *ArrBase = dyArr_init(Msh->NbrVer + 1, 1, INT); 
    ArrBase->SizCur = Msh->NbrVer + 1;
    int idPt;

    for (int iTri=1; iTri<=Msh->NbrTri; iTri++)
    {
        debug_printf("iTri = %d\n", iTri);
        for (int iEdgeLocal = 0; iEdgeLocal <=2; iEdgeLocal++)
        {
            idPt = Msh->Tri[iTri][iEdgeLocal];
            if (ArrBase->data1d[idPt] == 0)
            {
                ArrBase->data1d[idPt] = iTri;
            }
        } 
    }

    return ArrBase;
}

dynamArr *boucleDetection(Mesh *Msh, int iPt, int iTri, int *NbrBoucle)
{
    debug_printf("----------- [ Begin boucleDetection ] ------------\n");
    ////////// Initialisation de pile //////////
    dynamArr *ArrBoucle = dyArr_init(SizPil, 1, INT);
    dyArr_addElement_INT(ArrBoucle, iTri, 0);
    int ptCur = 0; int ptSiz = 1;
    int depile; int jTri;
    int step = 0;
    *NbrBoucle = 1;

    while (ptCur < ptSiz)
    {
        debug_printf("step = %d\n", step);
        depile = ArrBoucle->data1d[ptCur]; ptCur ++;
        for (int iEdgeLocal = 0; iEdgeLocal <= 2; iEdgeLocal ++)
        {
            jTri = Msh->TriVoi[depile][iEdgeLocal];
            if (!dyInList(jTri, ArrBoucle))
            {
                if (Msh->Tri[jTri][0] == iPt || Msh->Tri[jTri][1] == iPt
                    || Msh->Tri[jTri][2] == iPt)
                {
                    dyArr_auto_resize(ArrBoucle);
                    ArrBoucle->data1d[ptSiz] = jTri; ptSiz ++; ArrBoucle->SizCur ++;
                    *NbrBoucle = *NbrBoucle + 1;
                }
            }
        }
        step ++;
    }
    debug_printf("Boucle : \n");
    //dyArr_print(ArrBoucle);
    debug_printf("NbrBoucle = %d\n", *NbrBoucle);
    return ArrBoucle;
}

void export_NablaRUh_to_file(const char* filename, Mesh *Msh, double2d *NablaRUh)
{
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "Error opening file %s\n", filename);
        return;
    }
    for (int i=1; i<=Msh->NbrVer; i++)
    {
        fprintf(fp, "%d %lf %lf\n", i, NablaRUh[i][0], NablaRUh[i][1]);
    }
    fclose(fp);
    printf("NablaRUh data exported to %s\n", filename);
}


void export_Hessien_to_file(const char* filename, Mesh *Msh, double3d *Hessien)
{
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "Error opening file %s\n", filename);
        return;
    }
    for (int i=1; i<=Msh->NbrVer; i++)
    {
        fprintf(fp, "%d %lf %lf %lf\n", i, Hessien[i][0], Hessien[i][1], Hessien[i][2]);
    }
    fclose(fp);
    printf("Hessien data exported to %s\n", filename);
}