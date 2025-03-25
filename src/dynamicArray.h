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
void dyArr_export_to_file(dynamArr *dyArr, const char *filename);
void dyArr_print(dynamArr *dyArr);
dynamArr *dyArr_sort_circle(dynamArr *list_to_sort);
void dyArr_free(dynamArr *dyArr);
