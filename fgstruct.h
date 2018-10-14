//
// Created by Francisco Cruz on 27/12/16.
//

#ifndef AD3NEW_FGSTRUCT_H_H
#define AD3NEW_FGSTRUCT_H_H
#include <ctype.h>
#include <stdint.h>
#include "dynamicfd.h"

#define int_size 4

typedef struct buffer_read_def_ {
    uint32_t    id;
    uint32_t    n_configurations;
    double*     configurations;
} FileReadBuffer;

typedef struct variable_info_def_ {
    uint32_t    cMultiVariables;
    uint32_t    cBinaryVariables;
    double*     logPotentials;
    uint32_t*   variablesFrom;
    uint32_t*   variablesTo;
    uint32_t*   degree;
    double*     diff;
    double*     p;
    double*     sumQiAlpha;
    double*     outProbability;
    uint32_t*   binaryVariablesActive;
    uint32_t*   binaryBelongsTo;
    uint32_t**  connectedFactors;

} VarInfo;

typedef struct factor_info_def_ {
    double*     initialFactorPotentials;
    double*     currentFactorPotentials;
    double*     outputFactorPotentials;
    uint32_t*   potentialsStart;
    uint32_t*   potentialsEnd;
    uint32_t*   edgesStart;
    uint32_t*   edgesEnd;
    uint32_t*   factorsActive;
    uint32_t    cFactors;
    uint32_t    cFactorPotentials;
} FactorInfo;

typedef struct factor_def_ {
    uint32_t    variableLeft;
    uint32_t    variableRight;
    uint32_t    sizeI;
    uint32_t    sizeJ;
    dynamicFactorData* dynamicFactorData;
    double*     outVarLogPotentials;  // Could be moved to the edge

} Factor;

typedef struct edge_info_def_ {
    uint32_t    cEdges;
    double*     lamdas;
    double*     qs;
    double*     initialLogPotentials;
    double*     currentLogPotentials;
    double*     tmpDataComputeDual;
    uint32_t*   factors;
    uint32_t*   multiVariables;
    uint32_t*   binaryVariables;
    int*        edgeFactorActive;
    double **   logPotentialInFactor;
} EdgeInfo;

typedef struct factor_graph_def_ {

    double      optimal_energy;
    uint32_t *  optimal_configuration;
    VarInfo *   varInfo;
    FactorInfo* factorInfo;
    Factor*     factors;
    EdgeInfo*   edgeInfo;

} FactorGraph;





#endif //AD3NEW_FGSTRUCT_H_H
