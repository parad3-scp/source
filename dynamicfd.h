//
// Created by Francisco Cruz on 29/12/16.
//
#include "fgstruct.h"


#ifndef AD3NEW_DYNAMICFD_H
#define AD3NEW_DYNAMICFD_H

int __dynamicfd_currentSize;

typedef struct active_set_def_ {
    uint32_t allocated;
    uint32_t lenActiveSet;
    uint32_t *activeSet;
    uint32_t *tmpR;
    double   *varB;
    double   *tau;
    double   *z;
    double   *distribution;
    double   *matrixA;
    double   *inverseA;
    int      updatedInverseA;


} dynamicFactorData;

dynamicFactorData * createDynamicFactorData();
dynamicFactorData * expandDynamicFactorData(dynamicFactorData* dynamicdata);
uint32_t countAllocatedMemoryDynamicFD();


#endif //AD3NEW_DYNAMICFD_H
