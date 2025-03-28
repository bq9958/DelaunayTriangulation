#include <triangulation.h>

//////////// Global varialbes //////////
int SizPil = 100;       // cavity Tri Number max
const char *keyMode = "sum";

int main(int argc, char *argv[])
{
    //////// Mesh Definition //////
    Mesh *Msh = msh_init();
    Msh->NbrVer = 2395;
    int SizTri = 4*Msh->NbrVer;   
    Msh->Crd = TestPointSet(Msh->NbrVer, 0, 1, 0, 1);
    Msh->Tri = (int3d *)calloc(SizTri + 1, sizeof(int3d));      // le tableau actif
    Msh->TriVoi = (int3d *)calloc(SizTri + 1, sizeof(int3d));
    
    //////// Hash Table Definition ////////
    HashTable *hsh = hash_init(0.8*SizTri, 3*SizTri+1);

    
    //////// Other constants /////////
    int iTriLocLast = 0;


    //////// Initialize the bounding box /////////
    // leave speace to insert bounding box points
    for (int i=Msh->NbrVer; i>=1; i--)
    {
        Msh->Crd[i+4][0] = Msh->Crd[i][0];
        Msh->Crd[i+4][1] = Msh->Crd[i][1];
        Msh->Crd[i][0] = 0;
        Msh->Crd[i][1] = 0;
    }
    msh_boundingbox(Msh);
    Msh->Tri[1][0] = 1; Msh->Tri[1][1] = 2; Msh->Tri[1][2] = 4;
    Msh->Tri[2][0] = 2; Msh->Tri[2][1] = 3; Msh->Tri[2][2] = 4;
    Msh->NbrTri += 2;

    Msh->TriVoi[1][0] = 2; Msh->TriVoi[1][1] = 0; Msh->TriVoi[1][2] = 0;
    Msh->TriVoi[2][0] = 1; Msh->TriVoi[2][1] = 0; Msh->TriVoi[2][2] = 0;
    save_test_points_to_file("../output/test_points.txt", Msh);

    #if LOG_LEVEL == 4
        write_Tri_to_file("../output/Tri.txt", Msh);
        save_triangulation_to_file("../output/triangulation.mesh", Msh, hsh);
    #endif
    
    // ///////// Test session ////////
    double to, ti;
    to = GetWallClock();
    for (int iPtIns = 1; iPtIns<=Msh->NbrVer; iPtIns++)
    {
        printf("[[--------iPtIns-------]] = %d\n", iPtIns);
        Cavity(Msh, hsh, iPtIns, &iTriLocLast);
    }
    ti = GetWallClock();

    ////////// Output /////////
    #if LOG_LEVEL == 4
        write_LstObj_to_file("../output/LstObj.txt", hsh);
        write_Tri_to_file("../output/Tri.txt", Msh);
        write_TriVoi_to_file("../output/TriVoi.txt", Msh);
    #endif
    save_triangulation_to_file("../output/triangulation.mesh", Msh, hsh);
    printf("Time to insert all points: %f (s)\n", ti-to);

    // ///////// Free memory ////////
    msh_free(Msh);
    hash_free(hsh);
}