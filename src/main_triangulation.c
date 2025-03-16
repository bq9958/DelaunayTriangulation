#include <triangulation.h>

int main(int argc, char *argv[])
{
    //////// Mesh Definition //////
    Mesh *Msh = msh_init();
    Msh->NbrVer = 10;
    int SizTri = 3*Msh->NbrVer;
    Msh->Crd = TestPointSet(Msh->NbrVer);
    Msh->Tri = (int3d *)malloc((SizTri + 1) * sizeof(int3d));
    Msh->TriVoi = (int3d *)malloc((SizTri + 1) * sizeof(int3d));
    
    //////// Hash Table Definition ////////
    HashTable *hsh;
    HashTable *hsh_CavEdg = hash_init(10, Msh->NbrVer);  // Msh->NbrTri is a conservative choice
    const char *keyMode = "sum";
    
    //////// Other constants /////////
    int move = 1; int step = 0;


    //////// Initialize the bounding box /////////
    msh_boundingbox(Msh);
    Msh->Tri[1][0] = Msh->NbrVer+1; Msh->Tri[1][1] = Msh->NbrVer+2; Msh->Tri[1][2] = Msh->NbrVer+4;
    Msh->Tri[2][0] = Msh->NbrVer+2; Msh->Tri[2][1] = Msh->NbrVer+3; Msh->Tri[2][2] = Msh->NbrVer+4;
    Msh->NbrTri += 2;

    hsh = msh_neighbors(Msh, keyMode);

    
    ///////// Test session ////////
    int *mark = (int *)calloc(Msh->NbrTri+1, sizeof(int));
    printf("%f\n", Msh->Crd[1][0]);
    
    EdgBuffer(Msh, 1, 1, hsh_CavEdg, mark, &step, keyMode);   
    for (int i=0; i<Msh->NbrTri+1; i++){
        printf("mark[%d] %d\n", i, mark[i]);
    }
    removeEdgBuffer(Msh, hsh, hsh_CavEdg, keyMode);


    ////////// Output /////////
    write_LstObj_to_file("../output/LstObj_CavEdg.txt", hsh_CavEdg);
    write_LstObj_to_file("../output/LstObj.txt", hsh);
    write_Tri_to_file("../output/Tri.txt", Msh);
    write_TriVoi_to_file("../output/TriVoi.txt", Msh);
    save_triangulation_to_file("../output/triangulation.txt", Msh, hsh);
}