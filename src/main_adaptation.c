#include <adaptation.h>

const char *keyMode = "sum";
int SizPil = 100; 

int main()
{
    const char *meshFile = "../data/squarecircle.mesh";

    printf("[ [Read mesh] ]\n");
    int readEfr = 1;
    Mesh * Msh = msh_read(meshFile, readEfr);
    HashTable *hsh = msh_neighbors(Msh, keyMode);
    write_TriVoi_to_file("../output/TriVoi.txt", Msh);

    dynamArr *ArrBase = listBase(Msh);
    printf("ArrBase[1] = %d\n", ArrBase->SizCur);
    dyArr_export_to_file(ArrBase, "../output/Base.txt");

    int iPt = 10; int iTri = ArrBase->data1d[iPt]; 
    dynamArr *ArrBoucle = boucleDetection(Msh, iPt, iTri);
    dyArr_print(ArrBoucle);
}