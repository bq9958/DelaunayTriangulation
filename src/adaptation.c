#include <adaptation.h>

extern int SizPil;

dynamArr *listBase(Mesh *Msh)
{
    printf("--------- [ Begin listBase ] --------\n");
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

dynamArr *boucleDetection(Mesh *Msh, int iPt, int iTri)
{
    printf("----------- [Begin boucleDetection] ------------\n");
    ////////// Initialisation de pile //////////
    dynamArr *ArrBoucle = dyArr_init(SizPil, 1, INT);
    dyArr_addElement_INT(ArrBoucle, iTri, 0);
    int ptCur = 0; int ptSiz = 1;
    int depile; int jTri;
    int step = 0;

    while (ptCur < ptSiz)
    {
        printf("step = %d\n", step);
        depile = ArrBoucle->data1d[ptCur]; ptCur ++;
        for (int iEdgeLocal = 0; iEdgeLocal <= 2; iEdgeLocal ++)
        {
            jTri = Msh->TriVoi[depile][iEdgeLocal];
            if (!dyInList(jTri, ArrBoucle))
            {
                if (Msh->Tri[jTri][0] == iPt || Msh->Tri[jTri][1] == iPt
                    || Msh->Tri[jTri][2] == iPt)
                {
                    ArrBoucle->data1d[ptSiz] = jTri; ptSiz ++; ArrBoucle->SizCur ++;
                }
            }
        }
        step ++;
    }

    return ArrBoucle;
}