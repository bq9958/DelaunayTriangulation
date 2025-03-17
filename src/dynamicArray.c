#include <dynamicArray.h>

extern SizPil;

dynamArr *dyArr_init(int MaxSiz, int dim)
{
    dynamArr *dyArr;
    dyArr->dim = dim;
    dyArr->MaxSiz = SizPil;
    dyArr->SizCur = 0;
    if (dyArr->dim == 1)
    {
        dyArr->data1d = (int *)malloc(malloc(dyArr->MaxSiz * sizeof(int)));
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
}