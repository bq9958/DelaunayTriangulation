#include <triangulation.h>

//////////// Global varialbes //////////
int SizPil = 20;
const char *keyMode = "sum";

int main(int argc, char *argv[])
{
    //////// Mesh Definition //////
    Mesh *Msh = msh_init();
    Msh->NbrVer = 10;
    int SizTri = 4*Msh->NbrVer;   
    Msh->Crd = TestPointSet(Msh->NbrVer);
    Msh->Tri = (int3d *)malloc((SizTri + 1) * sizeof(int3d));      // le tableau actif
    Msh->TriVoi = (int3d *)malloc((SizTri + 1) * sizeof(int3d));
    int *pile = (int *)malloc(SizPil * sizeof(int));         // la pile
    
    //////// Hash Table Definition ////////
    HashTable *hsh;
    
    
    //////// Other constants /////////
    int move = 1; int step = 0;


    //////// Initialize the bounding box /////////
    msh_boundingbox(Msh);
    Msh->Tri[1][0] = Msh->NbrVer+1; Msh->Tri[1][1] = Msh->NbrVer+2; Msh->Tri[1][2] = Msh->NbrVer+4;
    Msh->Tri[2][0] = Msh->NbrVer+2; Msh->Tri[2][1] = Msh->NbrVer+3; Msh->Tri[2][2] = Msh->NbrVer+4;
    Msh->NbrTri += 2;

    hsh = msh_neighbors(Msh, keyMode);

    
    ///////// Test session ////////


    ////////// Output /////////
    write_LstObj_to_file("../output/LstObj.txt", hsh);
    write_Tri_to_file("../output/Tri.txt", Msh);
    write_TriVoi_to_file("../output/TriVoi.txt", Msh);
    save_triangulation_to_file("../output/triangulation.txt", Msh, hsh);
}