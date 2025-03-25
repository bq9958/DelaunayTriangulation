#include <dynamicArray.h>

extern int SizPil;

dynamArr *dyArr_init(int MaxSiz, int dim)
{
    dynamArr *dyArr = (dynamArr *)malloc(sizeof(dynamArr));
    if ( ! dyArr ) return NULL;

    dyArr->dim = dim;
    dyArr->MaxSiz = SizPil;
    dyArr->SizCur = 0;
    if (dyArr->dim == 1)
    {
        dyArr->data1d = (int *)malloc(dyArr->MaxSiz * sizeof(int));
        dyArr->data2d = NULL;
    }
    else if (dyArr->dim == 2)
    {
        dyArr->data2d = (int2d *)malloc(dyArr->MaxSiz * sizeof(int2d));
        dyArr->data1d = NULL;
    }

    return dyArr;
}

void dyArr_deleteElement(dynamArr *dyArr, int id)     // delete element at position id
{
    if (dyArr->dim == 1)
    {
        if (id < 0 || id >= dyArr->MaxSiz) {
            printf("Error : Out of Bound when deleting element\n");
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
            printf("Error : Out of Bound when deleting element\n");
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

void dyArr_addElement(dynamArr *dyArr, int element1, int element2) // add element to the end
{
    if (dyArr->SizCur >= dyArr->MaxSiz) {
        printf("Error: Array is full, cannot add more elements!\n");
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

int dyArr_ifInArr(dynamArr *dyArr, int element1, int element2)
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
        printf("Error: Cannot open file %s for writing!\n", filename);
        return;
    }

    if (dyArr->dim == 1)
    {
        for (int i=0; i < dyArr->SizCur; i++)
        {
            fprintf(fp, "%d\n", dyArr->data1d[i]);
        }
    }
    else if (dyArr->dim == 2)
    {
        for (int i=0; i < dyArr->SizCur; i++)
        {
            fprintf(fp, "%d %d\n", dyArr->data2d[i][0], dyArr->data2d[i][1]);
        }
    }

    fclose(fp);
}

void dyArr_print(dynamArr *dyArr)
{
    if (dyArr->dim == 1)
    {
        for (int i=0; i < dyArr->SizCur; i++)
        {
            printf("%d\n", dyArr->data1d[i]);
        }
    }
    else if (dyArr->dim == 2)
    {
        for (int i=0; i < dyArr->SizCur; i++)
        {
            printf("%d %d\n", dyArr->data2d[i][0], dyArr->data2d[i][1]);
        }
    }
}

dynamArr *dyArr_sort_circle(dynamArr *list_to_sort)
{
    if (list_to_sort->dim == 1)
        printf(" Warning: The input list is 1D !\n");
    else if (list_to_sort->dim == 2)
    {
        int *used = calloc(list_to_sort->SizCur, sizeof(int));
        dynamArr *list_sorted = dyArr_init(list_to_sort->MaxSiz, 2);

        int ObjTmp[2];
        ObjTmp[0] = list_to_sort->data2d[0][0]; ObjTmp[1] = list_to_sort->data2d[0][1];
        dyArr_addElement(list_sorted, ObjTmp[0], ObjTmp[1]); used[0] = 1;
        for (int i = 1; i < list_to_sort->SizCur; i++) {
            for (int j = 1; j < list_to_sort->SizCur; j++) {
                if (!used[j] && list_to_sort->data2d[j][0] == ObjTmp[1]) {
                    dyArr_addElement(list_sorted, list_to_sort->data2d[j][0], list_to_sort->data2d[j][1]);
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