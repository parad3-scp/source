//
// Created by Francisco Cruz on 29/12/16.
//

#ifndef AD3NEW_MATRIXINV_H
#define AD3NEW_MATRIXINV_H


#include <stdint.h>
#include "definitions.h"

void invertMatrix(double * src, double*dest, uint32_t dimension, double* scratch);
void invertMatrixMar(uint32_t *r, double *inverse_A,int size_A);
void InvertAfterRemovalMar(double * inverse_A_, int size_A, int removed_index);
void eigenDecomposition(double *matrix, double *dest, int dimension, double*scratch);



#endif //AD3NEW_MATRIXINV_H
