#include <triangulation.h>
#define EPSILON  1e-14

extern int SizPil;
extern const char *keyMode;    

void Cavity(Mesh *Msh, HashTable *hsh, int iPtIns, int *iTriLocLast)
{
    debug_printf("--------------Begin Cavity-------------\n");
    ////////// Initialization //////////
    dynamArr *PilElt = dyArr_init(SizPil, 1, INT);
    dynamArr *PilEdg = dyArr_init(SizPil, 2, INT);
    int ptrCur = 0; int ptrNxt = 1;
    int move = 0; int iTri = 0;
    
    ////////// Localisation ///////////
    debug_printf("--------------Begin location-------------\n");
    if (*iTriLocLast == 0)
        iTri = (rand()%(Msh->NbrTri-1) + 1);          
    else 
        iTri = *iTriLocLast;
    debug_printf("iTri = %d\n", iTri);
    int TriPtIns = location(Msh, iTri, Msh->Crd[4+iPtIns][0], Msh->Crd[4+iPtIns][1], &move);
    *iTriLocLast = TriPtIns;
    debug_printf("iTriLocLast = %d\n", *iTriLocLast);
    debug_printf("TriPtIns = %d\n", TriPtIns);
    PilElt->data1d[0] = TriPtIns; PilElt->SizCur = 1;
    PilEdg->data2d[0][0] = Msh->Tri[TriPtIns][1]; 
    PilEdg->data2d[0][1] = Msh->Tri[TriPtIns][2]; PilEdg->SizCur ++;
    PilEdg->data2d[1][0] = Msh->Tri[TriPtIns][2]; 
    PilEdg->data2d[1][1] = Msh->Tri[TriPtIns][0]; PilEdg->SizCur ++;
    PilEdg->data2d[2][0] = Msh->Tri[TriPtIns][0]; 
    PilEdg->data2d[2][1] = Msh->Tri[TriPtIns][1]; PilEdg->SizCur ++;

    ////////// Stack loop ////////////
    debug_printf("--------------Begin Stack loop-------------\n");
    int depile, jTri;
    int iVerCom1, iVerCom2, iVerComLocal1, iVerComLocal2;
    double x0, y0, x1, y1, x2, y2;
    double x = Msh->Crd[4+iPtIns][0]; double y = Msh->Crd[4+iPtIns][1]; 
    while (ptrCur < ptrNxt)
    {
        depile = PilElt->data1d[ptrCur]; ptrCur ++;
        debug_printf("ptrCur = %d\n", ptrCur);
        debug_printf("depile = %d\n", depile);
        for (int iEdglocal=0; iEdglocal < 3; iEdglocal++)
        {
            jTri = Msh->TriVoi[depile][iEdglocal];
            debug_printf("jTri = %d\n", jTri);
            //debug_printf("Press Enter to continue\n");
            debug_printf("PilElt[0] = %d, PilElt[1] = %d\n", PilElt->data1d[0], PilElt->data1d[1]);
            //getchar();
            if (! dyArr_ifInArr(PilElt, jTri, 0)){
                x0 = Msh->Crd[Msh->Tri[jTri][0]][0]; y0 = Msh->Crd[Msh->Tri[jTri][0]][1];
                x1 = Msh->Crd[Msh->Tri[jTri][1]][0]; y1 = Msh->Crd[Msh->Tri[jTri][1]][1];
                x2 = Msh->Crd[Msh->Tri[jTri][2]][0]; y2 = Msh->Crd[Msh->Tri[jTri][2]][1];
                if (inCircumcircle(x0, y0, x1, y1, x2, y2, x, y)) 
                {
                    debug_printf("[inCircumcircle]\n");
                    dyArr_auto_resize(PilElt);
                    PilElt->data1d[ptrNxt] = jTri; PilElt->SizCur ++; ptrNxt ++;
                    debug_printf("ptrNxt = %d\n", ptrNxt);
                    areteCommune(Msh, depile, jTri, &iVerCom1, &iVerCom2, &iVerComLocal1, &iVerComLocal2); // index for jTri
                    debug_printf("iVerCom1 = %d, iVerCom2 = %d\n", iVerCom1, iVerCom2);
                    debug_printf("iVerComLocal1 = %d, iVerComLocal2 = %d\n", iVerComLocal1, iVerComLocal2);
                    debug_printf("------- PilEdg before suppr ---------\n");
                    #if LOG_LEVEL == 4
                        dyArr_print(PilEdg);
                    #endif
                    supprAreteCom(iVerCom1, iVerCom2, PilEdg);
                    debug_printf("------- PilEdg after suppr ---------\n");
                    #if LOG_LEVEL == 4
                        dyArr_print(PilEdg);
                    #endif
                    debug_printf("-------- PilEdg after add ---------\n");
                    addTwoOtherSides(Msh, jTri, iVerComLocal1, iVerComLocal2, PilEdg);
                    #if LOG_LEVEL == 4
                        dyArr_print(PilEdg);
                    #endif
                    debug_printf("\n");
                }
            }
        }
    }

    debug_printf("PilEdg->SizCur = %d\n", PilEdg->SizCur);
    dynamArr *tmp = dyArr_sort_circle(PilEdg);    // sort the edges in a circle
    dyArr_free(PilEdg);
    PilEdg = tmp;
    debug_printf("-------- PilEdg after sort ---------\n");
    #if LOG_LEVEL == 4
        dyArr_print(PilEdg);
    #endif

    debug_printf("-----------------End Stack loop-----------------\n");

    ////////// Delete the triangles in the cavity //////////
    debug_printf("--------------Begin Delete the triangles in the cavity-------------\n");
    debug_printf("PilElt->SizCur = %d\n", PilElt->SizCur);
    for (int i=0; i < PilElt->SizCur; i++)
    {
        //Msh->Tri[PilElt->data1d[i]][0] = -Msh->Tri[PilElt->data1d[i]][0];
        Msh->NbrTri --;
        debug_printf("Triangle %d is deleted\n", PilElt->data1d[i]);
    }
    debug_printf("--------------End Delete the triangles in the cavity-------------\n");

    ////////////// Add new vertices and triangles //////////////
    debug_printf("--------------Begin Add new vertices and triangles-------------\n");
    dynamArr *CavTri = dyArr_init(SizPil*3, 1, INT);    // stocker iTri dans la cavité //TODO: Size control
    for (int i=0; i < PilEdg->SizCur; i++)
    {
        Msh->NbrTri ++;
        debug_printf("---- i = %d -----\n", i);
        debug_printf("PilElt->data1d[i] = %d\n", PilElt->data1d[i]);
        
        if (i < PilElt->SizCur)
        {
            // create new triangle
            Msh->Tri[PilElt->data1d[i]][0] = 4 + iPtIns;
            Msh->Tri[PilElt->data1d[i]][1] = PilEdg->data2d[i][0];
            Msh->Tri[PilElt->data1d[i]][2] = PilEdg->data2d[i][1];
            hash_add(hsh, Msh->Tri[PilElt->data1d[i]][0], Msh->Tri[PilElt->data1d[i]][1], PilElt->data1d[i], 2, keyMode);
            hash_add(hsh, Msh->Tri[PilElt->data1d[i]][2], Msh->Tri[PilElt->data1d[i]][0], PilElt->data1d[i], 1, keyMode);
            // update Edge 0
            int id = hash_find_modified(hsh, Msh->Tri[PilElt->data1d[i]][1], Msh->Tri[PilElt->data1d[i]][2], keyMode);
            if (id == 0) { LOG_INFO("Cannot find Edge 0 in hsh\n");}
            else { hsh->LstObj[id][2] = PilElt->data1d[i];}

            dyArr_addElement_INT(CavTri, PilElt->data1d[i], 0);
            debug_printf("-- CavTri -- : \n");
            #if LOG_LEVEL == 4
                dyArr_print(CavTri);
            #endif
            debug_printf("-- end CavTri --\n");

            // find another triangle that shares edge 0
            int Edgid = hash_find_modified(hsh, Msh->Tri[PilElt->data1d[i]][2], Msh->Tri[PilElt->data1d[i]][1], keyMode);
            debug_printf("Edgid = %d\n", Edgid);
            int neighbor;
            if (Edgid == 0)
            { LOG_INFO("Edgid is 0\n"); neighbor = 0;}
            else { neighbor = hsh->LstObj[Edgid][2]; }
            debug_printf("neighbor = %d\n", neighbor);
            if (neighbor != 0)
            {
                areteCommune(Msh, PilElt->data1d[i], neighbor, &iVerCom1, &iVerCom2, &iVerComLocal1, &iVerComLocal2);
                debug_printf("iVerCom1, iVerCom2 = %d, %d\n", iVerCom1, iVerCom2);
                Msh->TriVoi[neighbor][3-iVerComLocal1-iVerComLocal2] = PilElt->data1d[i];
                Msh->TriVoi[PilElt->data1d[i]][0] = neighbor;     
            } 
            else { LOG_INFO("Note: Neighbor not found\n"); Msh->TriVoi[PilElt->data1d[i]][0] = neighbor; }
        }
        else
        {
            Msh->Tri[Msh->NbrTri][0] = 4 + iPtIns;
            Msh->Tri[Msh->NbrTri][1] = PilEdg->data2d[i][0];
            Msh->Tri[Msh->NbrTri][2] = PilEdg->data2d[i][1];  
            hash_add(hsh, Msh->Tri[Msh->NbrTri][0], Msh->Tri[Msh->NbrTri][1], Msh->NbrTri, 2, keyMode);
            hash_add(hsh, Msh->Tri[Msh->NbrTri][2], Msh->Tri[Msh->NbrTri][0], Msh->NbrTri, 1, keyMode);
            int id = hash_find_modified(hsh, Msh->Tri[Msh->NbrTri][1], Msh->Tri[Msh->NbrTri][2], keyMode);
            if (id == 0) { LOG_INFO("Cannot find Edge 0 in hsh\n");}
            else { hsh->LstObj[id][2] = Msh->NbrTri;}
            dyArr_addElement_INT(CavTri, Msh->NbrTri, 0); 
            debug_printf("-- CavTri -- : \n");
            #if LOG_LEVEL == 4
                dyArr_print(CavTri);
            #endif
            debug_printf("-- end CavTri --\n");
            int Edgid = hash_find_modified(hsh, Msh->Tri[Msh->NbrTri][2], Msh->Tri[Msh->NbrTri][1], keyMode);
            debug_printf("Edgid = %d\n", Edgid);
            int neighbor;
            if (Edgid == 0)
            { LOG_INFO("Edgid is 0\n"); neighbor = 0;}
            else { neighbor = hsh->LstObj[Edgid][2]; }
            debug_printf("neighbor = %d\n", neighbor);
            if (neighbor != 0)
            {
                areteCommune(Msh, Msh->NbrTri, neighbor, &iVerCom1, &iVerCom2, &iVerComLocal1, &iVerComLocal2);
                debug_printf("iVerCom1, iVerCom2 = %d, %d\n", iVerCom1, iVerCom2);
                Msh->TriVoi[neighbor][3-iVerComLocal1-iVerComLocal2] = Msh->NbrTri; // update the neighbor of neighbor
                Msh->TriVoi[Msh->NbrTri][0] = neighbor;  // update one neighbor
            } 
            else { LOG_INFO("Note: Neighbor not found\n"); Msh->TriVoi[Msh->NbrTri][0] = neighbor; }
        }
        //getchar();
    }

    for (int i=0; i < CavTri->SizCur; i++)
    {
        if (i == 0)
        {
            Msh->TriVoi[CavTri->data1d[i]][2] = CavTri->data1d[CavTri->SizCur-1];
            Msh->TriVoi[CavTri->data1d[i]][1] = CavTri->data1d[i+1];
        }
        else if (i == CavTri->SizCur-1)
        {
            Msh->TriVoi[CavTri->data1d[i]][2] = CavTri->data1d[i-1];
            Msh->TriVoi[CavTri->data1d[i]][1] = CavTri->data1d[0];
        }
        else
        {
            Msh->TriVoi[CavTri->data1d[i]][2] = CavTri->data1d[i-1];
            Msh->TriVoi[CavTri->data1d[i]][1] = CavTri->data1d[i+1];
        }
    }
    debug_printf("--------------End Add new vertices and triangles-------------\n");
    debug_printf("-------------- End iPtIns = %d -------------\n", iPtIns);

    #if LOG_LEVEL == 4
        dyArr_export_to_file(PilElt, "../output/PilElt.txt");
        dyArr_export_to_file(PilEdg, "../output/PilEdg.txt");
    #endif

    dyArr_free(PilElt);
    dyArr_free(PilEdg);
    dyArr_free(CavTri);
}

void boucleDetection(Mesh *Msh, int iPt, int iTri, int *mark, int *step)
{
    int search = 3;
    mark[iTri] = 1;

    debug_printf("Step %d\n", *step); 
    (*step)++;

    for (int iEdglocal=0; iEdglocal < 3; iEdglocal++){
        int neighbor = Msh->TriVoi[iTri][iEdglocal];
        if (mark[neighbor] == 0 && (Msh->Tri[neighbor][0] == iPt || Msh->Tri[neighbor][1] == iPt 
            || Msh->Tri[neighbor][2] == iPt)){
            mark[neighbor] = 1;     // belongs to boucle
            boucleDetection(Msh, iPt, neighbor, mark, step);
        }
        else{
            if (mark[neighbor] != 1)
                mark[neighbor] = -1;    // not belong to boucle
            search--;
        }
        if (search == 0) {return;}
    }
}

int location(Mesh *Msh, int iTri, double x, double y, int *move)
{
    debug_printf("--------------Begin location Step [%d]-------------\n", *move);

    while (1) {
        (*move)++;
        if (*move > 2*Msh->NbrTri) {    //TODO : Optimize
            LOG_ERROR("Too many moves\n");
            return -1;
        }

        double x0 = Msh->Crd[Msh->Tri[iTri][0]][0];
        double y0 = Msh->Crd[Msh->Tri[iTri][0]][1];
        double x1 = Msh->Crd[Msh->Tri[iTri][1]][0];
        double y1 = Msh->Crd[Msh->Tri[iTri][1]][1];
        double x2 = Msh->Crd[Msh->Tri[iTri][2]][0];
        double y2 = Msh->Crd[Msh->Tri[iTri][2]][1];

        int position = inTriangle(x0, y0, x1, y1, x2, y2, x, y);

        if (position == -1) {
            return iTri;
        }
        else if (Msh->TriVoi[iTri][position] == 0) {
            LOG_ERROR("Reached boundary.\n");
            return -1;
        }

        iTri = Msh->TriVoi[iTri][position];
        debug_printf("Moved to iTri = %d\n", iTri);
    }

    return -1;
}

int inTriangle(double x0, double y0, double x1, double y1, double x2, double y2, double x, double y)
{
    double beta[3] = {0, 0, 0};
    int choice[3] = {0, 0, 0};
    int count = 0;

    barycenter(x0, y0, x1, y1, x2, y2, x, y, &beta[0], &beta[1], &beta[2]);
    for (int i=0; i<3; i++){
        if (beta[i] < -EPSILON)
        {
            count++;
            choice[count] = i;
        }
    }
    if (count == 1)    // one beta is negative
    {
        return choice[count];
    }
    else if (count == 2)  // two betas are negative
    {
        srand(time(NULL));
        int random_number = (rand() % 2) ? 2 : 1;
        return choice[random_number];
    }
    else if (count == 0)  // all beta are positive
    {
        return -1;
    }
    else
    {
        LOG_ERROR("Error: Invalid location\n");
        return -2;
    }
}


double2d *TestPointSet(int NbrPts, double xmin, double xmax, double ymin, double ymax)
{
    srand(NbrPts);

    // Pre-allocate memory for the points of boiteGlobante
    double2d *Pts = (double2d *)malloc((NbrPts + 5) * sizeof(double2d));

    if (!Pts) {
        LOG_ERROR("Memory allocation failed!\n");
        return NULL;
    }

    for (int i = 1; i <= NbrPts; i++) {
        double x = (double)rand() / RAND_MAX ;
        double y = (double)rand() / RAND_MAX ;
        Pts[i][0] = xmin + x * (xmax - xmin);
        Pts[i][1] = ymin + y * (ymax - ymin);
    }

    return Pts;
}

int inCircumcircle(double x1, double y1, double x2, double y2, double x3, double y3, double x, double y)
{
    // determine if (x,y) is inside the circumcircle of the triangle (x1,y1), (x2,y2), (x3,y3)
    double radius = rCircumcircle(x1, y1, x2, y2, x3, y3);
    Point center = centerCircumcircle(x1, y1, x2, y2, x3, y3);
    double dist = distance(center.x, center.y, x, y);
    if (dist < radius)
        return 1;
    return 0;
}

double rCircumcircle(double x1, double y1, double x2, double y2, double x3, double y3)
{
    double l1 = distance(x1, y1, x2, y2);
    double l2 = distance(x2, y2, x3, y3);
    double l3 = distance(x3, y3, x1, y1);
    double area = triArea(x1, y1, x2, y2, x3, y3);
    if (area < -EPSILON) { LOG_ERROR("area Negative\n");}

    return l1 * l2 * l3 / (4.0 * area); 
}

Point centerCircumcircle(double x1, double y1, double x2, double y2, double x3, double y3)
{
    Point center;

    center.x = ((y2-y1)*(y3*y3-y1*y1+x3*x3-x1*x1)-(y3-y1)*\
                (y2*y2-y1*y1+x2*x2-x1*x1))/(2.0*((x3-x1)*(y2-y1)-(x2-x1)*(y3-y1)));
    center.y = ((x2-x1)*(x3*x3-x1*x1+y3*y3-y1*y1)-(x3-x1)*\
                (x2*x2-x1*x1+y2*y2-y1*y1))/(2.0*((y3-y1)*(x2-x1)-(y2-y1)*(x3-x1)));

    return center;
}

void barycenter(double x0, double y0, 
                  double x1, double y1, 
                  double x2, double y2, 
                  double x, double y,
                  double *beta0, double *beta1, double *beta2)
{
    double aireNonSignee = triArea(x0, y0, x1, y1, x2, y2);
    if (aireNonSignee < -EPSILON) { LOG_ERROR("aireNonSignee Negative\n");}
    double PP1P2 = aireSignee(x, y, x1, y1, x2, y2);
    double P0PP2 = aireSignee(x0, y0, x, y, x2, y2);
    double P0P1P = aireSignee(x0, y0, x1, y1, x, y);

    *beta0 = PP1P2 / aireNonSignee;
    *beta1 = P0PP2 / aireNonSignee;
    *beta2 = P0P1P / aireNonSignee;
}

double aireSignee(double x0, double y0, double x1, double y1, double x2, double y2)
{
    double aireSignee = 0.5 * ((x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0));
    return aireSignee;
}

void areteCommune(Mesh *Msh, int iTri, int jTri, int *iVerCom1, int *iVerCom2, int *iVerComLocal1, int *iVerComLocal2)
{
    int list[3] = {Msh->Tri[iTri][0], Msh->Tri[iTri][1], Msh->Tri[iTri][2]};
    int count = 1;
    for (int i=0; i < 3; i++)
    {
        if (inList(Msh->Tri[jTri][i], list, 3))
        {
            if (count == 1) {*iVerCom1 = Msh->Tri[jTri][i]; *iVerComLocal1 = i; count ++;}
            else {*iVerCom2 = Msh->Tri[jTri][i]; *iVerComLocal2 = i;}   // index for jTri
        }
    }
}

void supprAreteCom(int iVerCom1, int iVerCom2, dynamArr *PilEdg)
{
    for (int i = 0; i < PilEdg->SizCur; i++){
        if ((PilEdg->data2d[i][0] == iVerCom1 && PilEdg->data2d[i][1] == iVerCom2) || 
            (PilEdg->data2d[i][0] == iVerCom2 && PilEdg->data2d[i][1] == iVerCom1))
        {
            dyArr_deleteElement(PilEdg, i);
        }
    }
}

void addTwoOtherSides(Mesh *Msh, int iTri, int iVerComLocal1, int iVerComLocal2, dynamArr *PilEdg)
{
    int iVerReste = 3 - iVerComLocal1 - iVerComLocal2;   // local index
    int iVerMin = iVerComLocal1 < iVerComLocal2 ? iVerComLocal1 : iVerComLocal2;  // local index
    int iVerMax = iVerComLocal1 > iVerComLocal2 ? iVerComLocal1 : iVerComLocal2;  // local index
    int iVerMinGl = Msh->Tri[iTri][iVerMin];  // global index
    int iVerMaxGl = Msh->Tri[iTri][iVerMax];  // global index
    int iVerResteGl = Msh->Tri[iTri][iVerReste];  // global index
    if (iVerReste == 1){
        if (!dyArr_ifInArr(PilEdg, iVerMinGl, iVerResteGl) && !dyArr_ifInArr(PilEdg, iVerResteGl, iVerMaxGl)){
            dyArr_addElement_INT(PilEdg, iVerMinGl, iVerResteGl);   // 0 -> 1
            dyArr_addElement_INT(PilEdg, iVerResteGl, iVerMaxGl);   // 1 -> 2
        }
    }
    else if (iVerReste == 2){
        if (!dyArr_ifInArr(PilEdg, iVerMaxGl, iVerResteGl) && !dyArr_ifInArr(PilEdg, iVerResteGl, iVerMinGl)){
            dyArr_addElement_INT(PilEdg, iVerMaxGl, iVerResteGl);   // 1 -> 2
            dyArr_addElement_INT(PilEdg, iVerResteGl, iVerMinGl);   // 2 -> 0
        }
    }
    else if (iVerReste == 0){
        if (!dyArr_ifInArr(PilEdg, iVerResteGl, iVerMaxGl) && !dyArr_ifInArr(PilEdg, iVerMinGl, iVerMaxGl)){
            dyArr_addElement_INT(PilEdg, iVerMaxGl, iVerResteGl);     // 2 -> 0
            dyArr_addElement_INT(PilEdg, iVerResteGl, iVerMinGl);   // 0 -> 1
        }
    }
}

int inList(int element, int *List, int Siz)
{
    for (int i=0; i < Siz; i++)
    {
        if (element == List[i])
            return 1;
    }
    return 0;
}

void save_test_points_to_file(const char *filename, Mesh *Msh)
{
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        LOG_ERROR("Cannot open file %s for writing!\n", filename);
        return;
    }
    for (int i = 1; i <= Msh->NbrVer + 4; i++) {
        if (i != Msh->NbrVer + 4)
            fprintf(fp, "%lf %lf\n", Msh->Crd[i][0], Msh->Crd[i][1]);
        else
        {
            fprintf(fp, "%lf %lf", Msh->Crd[i][0], Msh->Crd[i][1]);
        }
    }
    fclose(fp);
    printf("[Ouput File] test points written to %s\n", filename);
}

void save_triangulation_to_file(const char *filename, Mesh *Msh, HashTable *hsh) 
{
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        LOG_ERROR("Cannot open file %s for writing!\n", filename);
        return;
    }
    fprintf(fp, "MeshVersionFormatted 2\n\n");
    fprintf(fp, "Dimension 2\n\n");
    fprintf(fp, "Vertices\n");
    fprintf(fp, "%d\n", Msh->NbrVer+4);
    for (int i = 1; i <= Msh->NbrVer + 4; i++) {
        fprintf(fp, "%lf %lf %d\n", Msh->Crd[i][0], Msh->Crd[i][1], 0);
    }

    fprintf(fp, "\n\nTriangles\n");
    fprintf(fp, "%d\n", Msh->NbrTri);
    for (int i = 1; i <= Msh->NbrTri; i++) {
        fprintf(fp, "%d %d %d %d\n", Msh->Tri[i][0], Msh->Tri[i][1], Msh->Tri[i][2], 1);
    }
    fprintf(fp, "End\n");
    fclose(fp);
    printf("[Ouput File] triangulation written to %s\n", filename);
}


void addElement(int NbrTri, int SizVector, double3d *Tri, double3d *Voi)
{
    if (NbrTri >= SizVector) {
        SizVector *= 2;
        double3d *tmp1 =(double3d *)realloc(Tri, SizVector * sizeof(double3d));
        double3d *tmp2 = (double3d *)realloc(Voi, SizVector * sizeof(double3d));

        if (!tmp1 || !tmp2) {
            LOG_ERROR("Memory reallocation failed!\n");
            return;
        }
        else {
            Tri = tmp1;
            Voi = tmp2;
        }
    }
}

int hash_find_modified(HashTable *hsh, int iVer1, int iVer2, const char* keyMode)
{
  
	// to be implemented
	
	// return the id found (in LstObj ), if 0 the object is not in the list

  // ------------------------- Two hash keys -------------------------
  int key;
  if (strcmp(keyMode, "min") == 0){
    key = iVer1 > iVer2 ? iVer2 : iVer1;
  }
  else if (strcmp(keyMode, "sum") == 0){
    key = iVer1 + iVer2;
  }
  else if (strcmp(keyMode, "divide") == 0){
    key = (int)((1 + \ 
            (iVer1 > iVer2 ? iVer2 : iVer1) / (double)(iVer1 > iVer2 ? iVer1 : iVer2) ) * 1e9);
  }
  else {
    LOG_ERROR("Invalid keyMode\n");
    return 0;
  }
  //----------------------------------------------------------------

  int id = hsh->Head[key%(hsh->SizHead)];
  while (id != 0) {
    if ((hsh->LstObj[id][0] == iVer1 && hsh->LstObj[id][1] == iVer2)) 
    {
      return id;
    }
    id = hsh->LstObj[id][4];
  }

	return 0;
}