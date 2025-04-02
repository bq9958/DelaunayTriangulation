#include <dynamicArray.h>

extern int SizPil;

dynamArr *dyArr_init(int MaxSiz, int dim, DataType type)
{
    dynamArr *dyArr = (dynamArr *)malloc(sizeof(dynamArr));
    if ( ! dyArr ) return NULL;

    dyArr->dim = dim;
    dyArr->MaxSiz = MaxSiz;
    dyArr->SizCur = 0;
    dyArr->type = type;
    if (dyArr->dim == 1)
    {
        if (type == INT)
        {
            dyArr->data1d = (int *)calloc(dyArr->MaxSiz, sizeof(int));
            dyArr->data2d = NULL;
            dyArr->data1dd = NULL;
            dyArr->data2dd = NULL;
        }
        else if (type == DOUBLE)
        {
            dyArr->data1dd = (double *)calloc(dyArr->MaxSiz, sizeof(double));
            dyArr->data2dd = NULL;
            dyArr->data1d = NULL;
            dyArr->data2d = NULL;
        }
    }
    else if (dyArr->dim == 2)
    {
        if (type == INT)
        {
            dyArr->data2d = (int2d *)calloc(dyArr->MaxSiz, sizeof(int2d));
            dyArr->data1d = NULL;
            dyArr->data1dd = NULL;
            dyArr->data2dd = NULL;
        }
        else if (type == DOUBLE)
        {
            dyArr->data2dd = (double2d *)calloc(dyArr->MaxSiz, sizeof(double2d));
            dyArr->data1dd = NULL;
            dyArr->data1d = NULL;
            dyArr->data2d = NULL;
        }
    }

    return dyArr;
}

void dyArr_deleteElement(dynamArr *dyArr, int id)     // delete element at position id
{                                                     // for INT or DOUBLE
    if (dyArr->dim == 1)
    {
        if (id < 0 || id >= dyArr->MaxSiz) {
            LOG_ERROR("Out of Bound when deleting element\n");
            return;
        }
        
        for (int i = id; i < dyArr->MaxSiz - 1; i++) {
            dyArr->data1d[i] = dyArr->data1d[i + 1];
        }
        
        dyArr->data1d[dyArr->MaxSiz - 1] = 0;

        dyArr->SizCur --;
    }
    else if (dyArr->dim == 2)
    {
        if (id < 0 || id >= dyArr->MaxSiz) {
            LOG_ERROR("Out of Bound when deleting element\n");
            return;
        }
        
        for (int i = id; i < dyArr->MaxSiz - 1; i++) {
            dyArr->data2d[i][0] = dyArr->data2d[i + 1][0];
            dyArr->data2d[i][1] = dyArr->data2d[i + 1][1];
        }
        
        dyArr->data2d[dyArr->MaxSiz - 1][0] = 0; dyArr->data2d[dyArr->MaxSiz - 1][1] = 0;

        dyArr->SizCur --;
    }
}

void dyArr_addElement_DOUBLE(dynamArr *dyArr, double element1, double element2) // for DOUBLE
{
    int resize = dyArr_auto_resize(dyArr);
    if (resize < 0) {
        fprintf(stderr, "Resize failed!\n");
        return;
    }

    if (dyArr->SizCur >= dyArr->MaxSiz) {
        LOG_ERROR("Array is full, cannot add more elements!\n");
        return;
    }

    if (dyArr->dim == 1) {
        dyArr->data1dd[dyArr->SizCur] = element1;
    } else if (dyArr->dim == 2) {
        dyArr->data2dd[dyArr->SizCur][0] = element1;
        dyArr->data2dd[dyArr->SizCur][1] = element2;
    }

    dyArr->SizCur++;
}

void dyArr_addElement_INT(dynamArr *dyArr, int element1, int element2) // for INT
{
    if (dyArr_auto_resize(dyArr) < 0) {
        fprintf(stderr, "Resize failed!\n");
        return;
    }

    if (dyArr->SizCur >= dyArr->MaxSiz) {
        LOG_ERROR("Array is full, cannot add more elements!\n");
        return;
    }

    if (dyArr->dim == 1) {
        dyArr->data1d[dyArr->SizCur] = element1;
    } else if (dyArr->dim == 2) {
        dyArr->data2d[dyArr->SizCur][0] = element1;
        dyArr->data2d[dyArr->SizCur][1] = element2;
    }

    dyArr->SizCur++;
}

int dyArr_ifInArr(dynamArr *dyArr, int element1, int element2)   //! only for INT
{
    if (dyArr->dim == 1)
    {
        for (int i=0; i < dyArr->SizCur; i++)
        {
            if (element1 == dyArr->data1d[i])
                return 1;
        }
        return 0;
    }
    else if (dyArr->dim == 2)
    {
        for (int i=0; i < dyArr->SizCur; i++)
        {
            if (element1 == dyArr->data2d[i][0] && element2 == dyArr->data2d[i][1])
                return 1;
        }
        return 0;
    }

    return 0;
}

void dyArr_export_to_file(dynamArr *dyArr, const char *filename)
{
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        LOG_ERROR("Cannot open file %s for writing!\n", filename);
        return;
    }

    if (dyArr->dim == 1)
    {
        if (dyArr->type == INT)
        {
            for (int i=0; i < dyArr->SizCur; i++)
            {
                fprintf(fp, "%d\n", dyArr->data1d[i]);
            }
        }
        else if (dyArr->type == DOUBLE)
        {
            for (int i=0; i < dyArr->SizCur; i++)
            {
                fprintf(fp, "%lf\n", dyArr->data1dd[i]);
            }
        }
    }
    else if (dyArr->dim == 2)
    {
        if (dyArr->type == INT)
        {
            for (int i=0; i < dyArr->SizCur; i++)
            {
                fprintf(fp, "%d %d\n", dyArr->data2d[i][0], dyArr->data2d[i][1]);
            }
        }
        else if (dyArr->type == DOUBLE)
        {
            for (int i=0; i < dyArr->SizCur; i++)
            {
                fprintf(fp, "%lf %lf\n", dyArr->data2dd[i][0], dyArr->data2dd[i][1]);
            }
        }
    }
    fclose(fp);
}

void dyArr_print(dynamArr *dyArr)
{
    if (dyArr->dim == 1)
    {
        if (dyArr->type == INT)
        {
            for (int i=0; i < dyArr->SizCur; i++)
            {
                printf("%d\n", dyArr->data1d[i]);
            }
        }
        else if (dyArr->type == DOUBLE)
        {
            for (int i=0; i < dyArr->SizCur; i++)
            {
                printf("%lf\n", dyArr->data1dd[i]);
            }
        }
    }
    else if (dyArr->dim == 2)
    {
        if (dyArr->type == INT)
        {
            for (int i=0; i < dyArr->SizCur; i++)
            {
                printf("%d %d\n", dyArr->data2d[i][0], dyArr->data2d[i][1]);
            }
        }
        else if (dyArr->type == DOUBLE)
        {
            for (int i=0; i < dyArr->SizCur; i++)
            {
                printf("%f %f\n", dyArr->data2dd[i][0], dyArr->data2dd[i][1]);
            }
        }
    }
}

dynamArr *dyArr_sort_circle(dynamArr *list_to_sort)    //! only for INT
{
    if (list_to_sort->dim == 1)
        LOG_ERROR("The input list is 1D !\n");
    else if (list_to_sort->dim == 2)
    {
        int *used = calloc(list_to_sort->MaxSiz, sizeof(int));
        dynamArr *list_sorted = dyArr_init(list_to_sort->MaxSiz, 2, INT);

        int ObjTmp[2];
        ObjTmp[0] = list_to_sort->data2d[0][0]; ObjTmp[1] = list_to_sort->data2d[0][1];
        dyArr_addElement_INT(list_sorted, ObjTmp[0], ObjTmp[1]); used[0] = 1;
        for (int i = 1; i < list_to_sort->SizCur; i++) {
            for (int j = 1; j < list_to_sort->SizCur; j++) {
                if (!used[j] && list_to_sort->data2d[j][0] == ObjTmp[1]) {
                    dyArr_addElement_INT(list_sorted, list_to_sort->data2d[j][0], list_to_sort->data2d[j][1]);
                    ObjTmp[0] = list_to_sort->data2d[j][0]; ObjTmp[1] = list_to_sort->data2d[j][1];
                    used[j] = 1;
                    break;
                }
            }
        }
        free(used);
        return list_sorted;
    }
}

void dyArr_free(dynamArr *dyArr) {
    if (!dyArr) return;

    if (dyArr->dim == 1 && dyArr->data1d) {
        free(dyArr->data1d);
        dyArr->data1d = NULL;
    }
    else if (dyArr->dim == 2 && dyArr->data2d) {
        free(dyArr->data2d);
        dyArr->data2d = NULL;
    }

    free(dyArr);
}

int dyArr_auto_resize(dynamArr* dyArr) 
{
    if (dyArr->SizCur < dyArr->MaxSiz)
        return 0; 

    int newSize = dyArr->MaxSiz * 2;
    if (newSize == 0) newSize = 4;

    if (dyArr->dim == 1) {
        if (dyArr->type == INT) {
            int* newData = realloc(dyArr->data1d, newSize * sizeof(int));
            if (!newData) return -1;
            dyArr->data1d = newData;
        } else if (dyArr->type == DOUBLE) {
            double* newData = realloc(dyArr->data1dd, newSize * sizeof(double));
            if (!newData) return -1;
            dyArr->data1dd = newData;
        } else {
            return -2;
        }
        printf("dyArr resized to %d\n", newSize);
    }

    else if (dyArr->dim == 2) {
        if (dyArr->type == INT) {
            int2d *oldData = dyArr->data2d;
            int2d *newData = realloc(oldData, newSize * sizeof(int2d));
            if (!newData) return -1;
            dyArr->data2d = newData;
        } else if (dyArr->type == DOUBLE) {
            double2d *oldData = dyArr->data2dd;
            double2d *newData = realloc(oldData, newSize * sizeof(double2d));
            if (!newData) return -1;
            dyArr->data2dd = newData;
        } else {
            return -2; 
        }
        printf("dyArr resized to %d\n", newSize);
    }

    dyArr->MaxSiz = newSize;
    return 1;  
}