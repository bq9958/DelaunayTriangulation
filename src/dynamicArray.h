#include <mesh.h>

typedef struct {
    int *data1d;
    int2d *data2d;
    int SizCur;     
    int MaxSiz;  
    int dim;        
} dynamArr;

dynamArr *dyArr_init(int MaxSiz, int dim);
void dyArr_deleteElement(dynamArr *dyArr, int id);
void dyArr_addElement(dynamArr *dyArr, int element1, int element2);
int dyArr_ifInArr(dynamArr *dyArr, int element1, int element2);


