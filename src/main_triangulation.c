#include <triangulation.h>

int main(int argc, char *argv[])
{
    Mesh *Msh = msh_init();
    HashTable *hsh;

    const char *keymode = "sum";
    Msh->NbrVer = 10;
    int SizTri = 3*Msh->NbrVer;
    int move = 0; int step = 0;

    // Allocate memory for Msh
    Msh->Crd = TestPointSet(Msh->NbrVer);
    Msh->Tri = (int3d *)malloc((SizTri + 1) * sizeof(int3d));
    Msh->TriVoi = (int3d *)malloc((SizTri + 1) * sizeof(int3d));

    // Initialize the bounding box
    msh_boundingbox(Msh);
    Msh->Tri[1][0] = Msh->NbrVer+1; Msh->Tri[1][1] = Msh->NbrVer+2; Msh->Tri[1][2] = Msh->NbrVer+4;
    Msh->Tri[2][0] = Msh->NbrVer+2; Msh->Tri[2][1] = Msh->NbrVer+3; Msh->Tri[2][2] = Msh->NbrVer+4;
    Msh->NbrTri += 2;

    hsh = msh_neighbors(Msh, keymode);

    write_TriVoi_to_file("../output/TriVoi.txt", Msh);

    //////////////////////////////////// Test session ///////////////////////////////////////
    int *mark = (int *)calloc(Msh->NbrTri+1, sizeof(int));
    boucleDetection(Msh, 14, 1, mark, &step);    
    for (int i=0; i<Msh->NbrTri+1; i++){
        printf("mark[%d] %d\n", i, mark[i]);
    }

    save_triangulation_to_file("../output/triangulation.txt", Msh, hsh);
    plot_with_gnuplot("../output/triangulation.txt");

}