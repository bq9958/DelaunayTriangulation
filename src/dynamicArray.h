#include <mesh.h>
typedef enum {
    INT,
    DOUBLE
} DataType;

typedef struct {
    int *data1d;
    int2d *data2d;
    double *data1dd;
    double2d *data2dd;
    int SizCur;     
    int MaxSiz;  
    DataType type;
    int dim;        
} dynamArr;

dynamArr *dyArr_init(int MaxSiz, int dim, DataType type);
void dyArr_deleteElement(dynamArr *dyArr, int id);
void dyArr_addElement_INT(dynamArr *dyArr, int element1, int element2);
void dyArr_addElement_DOUBLE(dynamArr *dyArr, double element1, double element2);
int dyArr_ifInArr(dynamArr *dyArr, int element1, int element2);
void dyArr_export_to_file(dynamArr *dyArr, const char *filename);
void dyArr_print(dynamArr *dyArr);
dynamArr *dyArr_sort_circle(dynamArr *list_to_sort);
void dyArr_free(dynamArr *dyArr);
int dyArr_auto_resize(dynamArr* dyArr);
