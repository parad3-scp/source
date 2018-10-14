//
// Created by Francisco Cruz on 29/12/16.
//

#include <stdlib.h>
#include <memory.h>
#include <printf.h>
#include <assert.h>
#include "dynamicfd.h"


#define initial_allocation 50
#define growing_allocation_step 10

dynamicFactorData * createDynamicFactorData() {

    dynamicFactorData* dynamicData = (dynamicFactorData*) malloc(sizeof(dynamicFactorData));
    dynamicData->activeSet    = (uint32_t *) malloc(sizeof(uint32_t)*initial_allocation);
    dynamicData->tmpR         = (uint32_t *) malloc(sizeof(uint32_t)*(initial_allocation+1));
    dynamicData->distribution = (double *)   malloc(sizeof(double)*initial_allocation);
    dynamicData->tau          = (double *)   malloc(sizeof(double)*(initial_allocation+1));
    dynamicData->z            = &dynamicData->tau[1];
    dynamicData->varB         = (double *)   malloc(sizeof(double)*(initial_allocation+1));
    dynamicData->matrixA      = (double *)   malloc(sizeof(double)*(initial_allocation+1)*(initial_allocation+1));
    dynamicData->inverseA     = (double *)   malloc(sizeof(double)*(initial_allocation+1)*(initial_allocation+1));
    dynamicData->allocated=initial_allocation;
    dynamicData->lenActiveSet=0;
    __dynamicfd_currentSize = initial_allocation;
    return dynamicData;
}

dynamicFactorData * expandDynamicFactorData(dynamicFactorData* dynamicdata) {

    if (dynamicdata->lenActiveSet == dynamicdata->allocated) {
        uint32_t  newsize = dynamicdata->allocated + growing_allocation_step;
        uint32_t* newDataActiveSet  = (uint32_t *) malloc(sizeof(uint32_t)* newsize);
        uint32_t* newTmpR           = (uint32_t *) malloc(sizeof(uint32_t)* (newsize+1));
        double * newDistribution    = (double *)   malloc(sizeof(double)*   newsize);
        double * newTau             = (double *)   malloc(sizeof(double)*   (newsize+1));
        double * newVarB            = (double *)   malloc(sizeof(double)*   (newsize+1));
        double * newMatrixA         = (double *)   malloc(sizeof(double)*   (newsize+1)*(newsize+1));
        double * newInverseA        = (double *)   malloc(sizeof(double)*   (newsize+1)*(newsize+1));

        memcpy(newDataActiveSet ,  dynamicdata->activeSet,   sizeof(uint32_t)* dynamicdata->allocated);
        // No need to memcpy tmpR, new distribution
        // memcpy(newDistribution,    dynamicdata->distribution,sizeof(double)*   dynamicdata->allocated);
        memcpy(newTau,             dynamicdata->tau,         sizeof(double)*   (dynamicdata->allocated+1));
        memcpy(newVarB,            dynamicdata->varB,        sizeof(double)*   (dynamicdata->allocated+1));
        memcpy(newMatrixA,         dynamicdata->matrixA,     sizeof(double)*   (dynamicdata->allocated+1)*(dynamicdata->allocated+1));
        memcpy(newInverseA,        dynamicdata->inverseA,    sizeof(double)*   (dynamicdata->allocated+1)*(dynamicdata->allocated+1));


        free(dynamicdata->activeSet);
        free(dynamicdata->tmpR);
        free(dynamicdata->distribution);
        free(dynamicdata->tau);
        free(dynamicdata->varB);
        free(dynamicdata->matrixA);
        free(dynamicdata->inverseA);

        dynamicdata->activeSet    = newDataActiveSet;
        dynamicdata->tmpR         = newTmpR;
        dynamicdata->distribution = newDistribution;
        dynamicdata->tau          = newTau;
        dynamicdata->z            = &dynamicdata->tau[1];
        dynamicdata->varB         = newVarB;
        dynamicdata->matrixA      = newMatrixA;
        dynamicdata->inverseA     = newInverseA;
        dynamicdata->allocated=newsize;
        __dynamicfd_currentSize = newsize;
        //printf("Extra allocation of memory for activeSet: new size %i\n",newsize);
    }

    return dynamicdata;
}

uint32_t countAllocatedMemoryDynamicFD() {
    uint32_t count=0;
    count += sizeof(uint32_t) * __dynamicfd_currentSize;     // activeSet
    count += sizeof(uint32_t) * (__dynamicfd_currentSize+1); // tmpR
    count += sizeof(double) * __dynamicfd_currentSize;       // distribution
    count += sizeof(double) * (__dynamicfd_currentSize+1)*2; // tau && varB
    count += sizeof(double) * (__dynamicfd_currentSize+1)*(__dynamicfd_currentSize+1)*2; // matrices

    return count;
}