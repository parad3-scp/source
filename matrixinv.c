//
// Created by Francisco Cruz on 29/12/16.
//
#include "definitions.h"
#include <memory.h>
#include <stdint.h>
#include <stdlib.h>
#include <lapacke.h>
#include "matrixinv.h"




void invertMatrixOld(double * src, double*dest, uint32_t dimension) {
    int i, j, k, temp;
    double temporary, r ;
    double *src_copy = (double*) malloc(sizeof(double)*dimension*dimension);
    memcpy(src_copy,src,sizeof(double)*dimension*dimension);

    // Init output
    for(i=0;i<dimension; i++)
        for(j=0; j<dimension; j++)
            if(i==j%dimension)  dest[i+j*dimension] = 1;
            else dest[i+j*dimension] = 0;


    for(j=0; j<dimension; j++)
    {
        temp=j;
        /* finding maximum jth column element in last (dimension-j) rows */
        for(i=j+1; i<dimension; i++)
            if(src_copy[i+j*dimension]>src_copy[temp+j*dimension])
                temp=i;

        /* swapping row which has maximum jth column element */
        if(temp!=j)
            for(k=0; k<dimension; k++)
            {
                temporary=src_copy[j+k*dimension] ;
                src_copy[j+k*dimension]=src_copy[temp+k*dimension] ;
                src_copy[temp+k*dimension]=temporary ;
                temporary = dest[j+k*dimension];
                dest[j+k*dimension] = dest[temp+k*dimension];
                dest[temp+k*dimension] = temporary;
            }
        /* performing row operations to form required identity matrix out of the input matrix */

        for(i=0; i<dimension; i++)
            if(i!=j) {
                r=src_copy[i+j*dimension];
                for(k=0; k<dimension; k++) {
                    src_copy[i+k*dimension] -= (src_copy[j+k*dimension] / src_copy[j+j*dimension]) * r;
                    dest[i+k*dimension]-=(dest[j+k*dimension]/src_copy[j+j*dimension])*r ;
                }
            } else {
                r=src_copy[i+j*dimension];
                for(k=0; k<dimension; k++) {
                    src_copy[i+k*dimension] /= r;
                    dest[i+k*dimension]/=r ;
                }
            }
    }
}

void invertMatrix(double * src, double*dest, uint32_t dimension,double *scratch) {
    int i, j, k, temp;
    double temporary, r ;
    int scratch_size=SCRATCH_SIZE;
    int * ipiv = (int * ) malloc(sizeof(int)*dimension);
    memcpy(dest,src,sizeof(double)*dimension*dimension);
    int n=dimension;
    int lda=dimension;
    int info=dimension;
    dgetrf_(&n,&n,dest,&lda,ipiv,&info);
    dgetri_(&n,dest,&lda,ipiv,scratch,&scratch_size,&info);
    free(ipiv);    
}


void eigenDecomposition(double *matrix, double *dest, int dimension, double *scratch) {
    int info;
    int scratch_size=SCRATCH_SIZE;
    dsyev_("V", "U", &dimension, matrix, &dimension, dest,scratch,&scratch_size,&info);
}

void invertMatrixMar(uint32_t *r, double *inverse_A,int size_A) {


    double s=2.0;
    for (int i = 0; i < size_A; ++i) {
        if (r[i] == 0.0) continue;
        s -= r[i] * r[i] * inverse_A[i * size_A + i];
        for (int j = i+1; j < size_A; ++j) {
            if (r[j] == 0.0) continue;
            s -= 2 * r[i] * r[j] * inverse_A[i * size_A + j];
        }
    }

    double invs = 1.0 / s;

    double * d = (double*) malloc(sizeof(double)*size_A);
    for (int i=0;i<size_A;i++) {
        d[i]=0.0;
    }

    for (int i = 0; i < size_A; ++i) {
        if (r[i] == 0.0) continue;
        for (int j = 0; j < size_A; ++j) {
            d[j] += inverse_A[i * size_A + j] * r[i];
        }
    }

    int size_A_after = size_A + 1;
    double *dest = (double*) malloc(sizeof(double)*size_A_after*size_A_after);

    for (int i = 0; i < size_A; ++i) {
        for (int j = 0; j < size_A; ++j) {
            dest[i * size_A_after + j] = inverse_A[i * size_A + j] +
                                               invs * d[i] * d[j];
        }
        dest[i * size_A_after + size_A] = -invs * d[i];
        dest[size_A * size_A_after + i] = -invs * d[i];
    }
    dest[size_A * size_A_after + size_A] = invs;

    memcpy(inverse_A,dest,sizeof(double)*size_A_after*size_A_after);
    free(d);
    free(dest);
}


void InvertAfterRemovalMar(double * inverse_A_, int size_A, int removed_index) {



    double * inverse_A = (double *) malloc(sizeof(double)*size_A*size_A);
    memcpy (inverse_A,inverse_A_,sizeof(double)*size_A*size_A);
    double *r  = (double *) malloc(sizeof(double)*size_A);
    double *d  = (double *) malloc(sizeof(double)*size_A-1);
    for (int i=0;i<size_A-1;i++) d[i]=0.0;


    ++removed_index; // Index in A has an offset of 1.
    double invs = inverse_A[removed_index * size_A + removed_index];

    double s = 1.0 / invs;

    int k = 0;
    for (int i = 0; i < size_A; ++i) {
        if (i == removed_index) continue;
        d[k] = -s * inverse_A[removed_index * size_A + i];
        ++k;
    }

    int size_A_after = size_A - 1;

    k = 0;
    for (int i = 0; i < size_A; ++i) {
        if (i == removed_index) continue;
        int l = 0;
        for (int j = 0; j < size_A; ++j) {
            if (j == removed_index) continue;
            inverse_A_[k * size_A_after + l] = inverse_A[i * size_A + j] - invs * d[k] * d[l];
            ++l;
        }
        ++k;
    }

}