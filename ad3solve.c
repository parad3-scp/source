//
// Created by Francisco Cruz on 28/12/16.
//

#define VERSION2

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "ad3solve.h"
#include "util.h"
#include "matrixinv.h"
#include "definitions.h"

#ifndef SEQ
#include "omp.h"
#endif

#define true 1
#define false 0

double  ad3_eta = 0.1;
int     ad3_adapt_eta = true;
int     ad3_num_iterations_adapt_eta = 10;
double  ad3_residual_threshold = 1e-6;
int     ad3_num_iterations_compute_dual = 50;
double  ad3_lower_bound = -1e100;
int     ad3_num_iterations_reset = 50;
double  ad3_cache_tolerance = 1e-12;
// Stepsize adjustment parameters
double  ad3_max_eta = 100.0;
double  ad3_min_eta = 1e-3;
double  ad3_gamma_primal = 100.0;
double  ad3_gamma_dual = 10.0;
double  ad3_tau = 1.0;
double  ad3_factor_step = 2.0;


double computeExtraScore(FactorGraph * fg) {
    double extra_score=0.0;
    for (int mv=0;mv<fg->varInfo->cMultiVariables;mv++) {
        int vardegree=fg->varInfo->degree[mv];

        double max=-1e100;
        if (vardegree==0) {
            for (int bv=fg->varInfo->variablesFrom[mv]; bv< fg->varInfo->variablesTo[mv];bv++) {
                if (fg->varInfo->logPotentials[bv] > max ) max = fg->varInfo->logPotentials[bv];
            }
            extra_score+=max;
        }
    }

    return extra_score;
}


void updateEdgeAndVar(FactorGraph *fg, int recompute_everything, int eta_changed, int iteration) {
    
    int threads=omp_get_num_threads();
    int thread_id=omp_get_thread_num(); 

    #pragma omp for
    for (int edge=0; edge< fg->edgeInfo->cEdges;edge++) {    
        if (recompute_everything || fg->edgeInfo->edgeFactorActive[edge]) {
            int bv = fg->edgeInfo->binaryVariables[edge];
            double vp = *(fg->edgeInfo->logPotentialInFactor[edge]);
            fg->varInfo->sumQiAlpha[bv*threads+thread_id] += vp - fg->edgeInfo->qs[edge]; 
            double p_i = fg->varInfo->p[bv];
             if ((eta_changed) ||
                (!NEARLY_BINARY(vp,ad3_cache_tolerance)) ||
                (!NEARLY_EQ_TOL(vp,fg->edgeInfo->qs[edge],ad3_cache_tolerance)) ||
                (!NEARLY_EQ_TOL(vp,p_i,ad3_cache_tolerance)))
                fg->varInfo->binaryVariablesActive[bv] = true; // true
            fg->edgeInfo->qs[edge]= vp;
        } 
    }
}


void propagateFactorActivationsToEdges(FactorGraph *fg) {

/*
    for (int factor=0; factor<fg->factorInfo->cFactors;factor++) {
        int activated=false;

        for (int edge=fg->factorInfo->edgesStart[factor]; edge< fg->factorInfo->edgesEnd[factor];edge++ ) {
            activated = activated || fg->edgeInfo->edgeFactorActive[edge];
        }
        fg->factorInfo->factorsActive[factor]=activated;
    } */
    #pragma omp parallel for    
    for (int factor=0; factor<fg->factorInfo->cFactors;factor++) {
        int activated=fg->factorInfo->factorsActive[factor];
        #pragma simd 
        for (int edge=fg->factorInfo->edgesStart[factor]; edge< fg->factorInfo->edgesEnd[factor];edge++ ) {
            fg->edgeInfo->edgeFactorActive[edge] =activated;
        }
    }
}

void findMax(double * col, uint32_t sizeCol, double * row, uint32_t sizeRow, double* matrix, uint32_t * ret_pos, double * ret_val) {
    *ret_pos = 0;
    *ret_val = matrix[0]+ col[0] + row[0];
    #pragma ivdep
    for (uint32_t idx=0;idx < sizeCol*sizeRow; idx ++) {
        uint32_t pos_row = idx % sizeRow;
        uint32_t pos_col = idx / sizeRow;
        double cell_value=matrix[idx]+ col[pos_col] + row[pos_row];
        if (cell_value > *ret_val) {
            *ret_val=cell_value;
            *ret_pos=idx;
        }
    }
}

void findMaxMinusExtra(double * col, uint32_t sizeCol, double * row, uint32_t sizeRow, double* matrix,  double *extra_col, double* extra_row,uint32_t * ret_pos, double * ret_val) {
    *ret_pos = 0;
    *ret_val = matrix[0]+ col[0] + row[0] - extra_col[0] - extra_row[0];
    #pragma ivdep
    for (uint32_t idx=0;idx < sizeCol*sizeRow; idx ++) {
        uint32_t pos_row = idx % sizeRow;
        uint32_t pos_col = idx / sizeRow;
        double tmp1;
        double tmp2;
        tmp1=row[pos_row]- extra_row[pos_row];
        tmp2=col[pos_col] -extra_col[pos_col];
        double cell_value=matrix[idx]+ tmp1 + tmp2 ;

        if (cell_value > *ret_val) {
            *ret_val=cell_value;
            *ret_pos=idx;
        }
    }
}

void findMaxMinusExtraPrint(double * col, uint32_t sizeCol, double * row, uint32_t sizeRow, double* matrix,  double *extra_col, double* extra_row,uint32_t * ret_pos, double * ret_val) {
    *ret_pos = 0;
    *ret_val = matrix[0]+ col[0] + row[0] - extra_col[0] - extra_row[0];
#pragma ivdep
    for (uint32_t idx=0;idx < sizeCol*sizeRow; idx ++) {
        uint32_t pos_row = idx % sizeRow;
        uint32_t pos_col = idx / sizeRow;
        double cell_value=matrix[idx]+ col[pos_col] + row[pos_row] - extra_col[pos_col] - extra_row[pos_row];
        printf("pos %i (%i,%i) total %5.20f: %5.20f %5.20f %5.20f %5.20f %5.20f\n",idx,pos_col,pos_row,cell_value, matrix[idx], col[pos_col], row[pos_row] , extra_col[pos_col], extra_row[pos_row]);
        if (cell_value > *ret_val) {
            *ret_val=cell_value;
            *ret_pos=idx;
        }
    }
}

double computeSingle (double * col, uint32_t sizeCol, double * row, uint32_t sizeRow, double* matrix,  double *extra_col, double* extra_row,uint32_t pos) {

    uint32_t pos_row = pos % sizeRow;
    uint32_t pos_col = pos / sizeRow;
    return matrix[pos]+ col[pos_col] + row[pos_row] - extra_col[pos_col] - extra_row[pos_row];
}


void computeInverseA(dynamicFactorData* data, double *scratch) {

    int i,j,n,k;
    n = data->lenActiveSet +1;
    double *inverseA = data->inverseA;
    double *matrixA = data->matrixA;
    invertMatrix(matrixA,inverseA,n, scratch);
    data->updatedInverseA = true;

}

void printMatrixRect(double *data, int size_x, int size_y) {

    for (int x = 0; x < size_x; x++) {
        for (int y = 0; y < size_y; y++) {
            printf("%g ", data[x+y*size_x]);
        }
        printf("\n");
    }
}

void printMatrix(double *data, int sizesqrt) {

    for (int x = 0; x < sizesqrt; x++) {
        for (int y = 0; y < sizesqrt; y++) {
            printf("%g ", data[x+y*sizesqrt]);
        }
        printf("\n");
    }
}
void printMatrixInt(int *data, int sizesqrt) {

    for (int x = 0; x < sizesqrt; x++) {
        for (int y = 0; y < sizesqrt; y++) {
            printf("%i ", data[x+y*sizesqrt]);
        }
        printf("\n");
    }
}

void printVectorInt(int *data, int sizes) {

    for (int x = 0; x < sizes; x++) {
        if (x>0) printf(", ");
        else printf("[");
        printf("%i", data[x]);

    }
    printf("]\n");
}

void printVector(double *data, int sizes) {

    for (int x = 0; x < sizes; x++) {
        printf("%5.16f ", data[x]);

    }
    printf("\n");
}


void    matrixMultVector(dynamicFactorData * data) {
    // Mutplies inverseA matrix per varB vector
    // Output is stored as follows
    // [0] : tau
    // [1..n] : z
    // in order to save computation tau and z are allocated toghether in the memory space.
    double *output = &data->tau[0];
    int width=data->lenActiveSet+1;
    for (int x=0;x<width;x++) {
        output[x]=0;
        for (int y =0; y < width; y++ ){
            output[x]+=data->inverseA[x*width + y] * data->varB[y];
        }
    }
}

void computeMarginalsFromSparseDistribution(Factor * factor, FactorGraph *fg, int factor_id){
    for (int i=0;i<factor->sizeI+factor->sizeJ;i++) {
        factor->outVarLogPotentials[i]=0.0;
    }

    for (int f=fg->factorInfo->potentialsStart[factor_id];f<fg->factorInfo->potentialsEnd[factor_id];f++) {
        fg->factorInfo->outputFactorPotentials[f]=0.0;
    }
    double *outputFactorPotentials = &fg->factorInfo->outputFactorPotentials[fg->factorInfo->potentialsStart[factor_id]];

    for (int i=0;i<factor->dynamicFactorData->lenActiveSet;i++) {
        int pos_x = factor->dynamicFactorData->activeSet[i] / factor->sizeJ;
        int pos_y = factor->dynamicFactorData->activeSet[i] % factor->sizeJ;
        factor->outVarLogPotentials[pos_x] += factor->dynamicFactorData->distribution[i];
        factor->outVarLogPotentials[pos_y+factor->sizeI] += factor->dynamicFactorData->distribution[i];
        outputFactorPotentials[pos_y+pos_x*factor->sizeJ] = factor->dynamicFactorData->distribution[i];
    }
}

void insertRowColumnIntoMatrixA(dynamicFactorData * data, uint32_t *r) {
    // Since InverseA matrix is going to be invalidated we use that allocated memory to skip the allocation of a new matrix
    // InvA will we overwritten anyways

    int dim = data->lenActiveSet+1;
    for (int row=0; row<dim;row++) {
        for (int col=0; col<dim;col++) {
            data->inverseA[row*(dim+1)+col] = data->matrixA[row*dim + col];
        }
        data->inverseA[row*(dim+1)+dim] = r[row];
    }
    for (int col=0;col<(dim);col++) {
        data->inverseA[dim*(dim+1)+col] = r[col];
    }
    data->inverseA[(dim+1)*(dim+1)-1]=2;
    double * tmp = data->matrixA;
    data->matrixA= data->inverseA;
    data->inverseA = tmp;
    data->updatedInverseA=false;
}
void removeVarFromEverywhere(dynamicFactorData * data, int rmv, int keepdistribution, int keepz) {

    int matrix_dest_pos=0;
    int vector_dest_pos=0;

    double * distributionplusone=&data->distribution[0];
    double * zplusone= &data->z[0];
    vector_dest_pos=0;
    for (int i=0; i < (data->lenActiveSet);i++) {
        if (i!=rmv) {
            data->activeSet[vector_dest_pos] = data->activeSet[i];
            if (!keepz) zplusone[vector_dest_pos] = zplusone[i];
            if (!keepdistribution) distributionplusone[vector_dest_pos] = distributionplusone[i];
            vector_dest_pos++;
        }
    }
    if (!keepdistribution) distributionplusone[data->lenActiveSet-1] =0.0;
    for (int i=0; i < (data->lenActiveSet+1)*(data->lenActiveSet+1);i++) {
        int row= i / (data->lenActiveSet +1);
        int column = i % (data->lenActiveSet +1);
        if ((row!=rmv+1) && (column!=rmv+1)) {
            data->matrixA[matrix_dest_pos] = data->matrixA[i];
            matrix_dest_pos++;
        }

    }
    data->updatedInverseA=false;
    data->lenActiveSet -=1;

}

void removeBlockingFromEverywhere(dynamicFactorData * data, int blocking) {

    int matrix_dest_pos=0;
    int vector_dest_pos=0;

    for (int i=0; i < (data->lenActiveSet+1)*(data->lenActiveSet+1);i++) {
        int row= i / (data->lenActiveSet +1);
        int column = i % (data->lenActiveSet +1);
        if ((row==0) && (column!=(blocking))) {
            data->activeSet[vector_dest_pos] =  data->activeSet[i];
            data->z[vector_dest_pos] =  data->z[i];
            data->distribution[vector_dest_pos] =  data->distribution[i];
            vector_dest_pos++;
        }

        if ((row!=blocking+1) && (column!=blocking+1)) {
            data->matrixA[matrix_dest_pos] = data->matrixA[i];
            matrix_dest_pos++;
        }

    }
    data->updatedInverseA=false;
    data->lenActiveSet -=1;

}

void printActiveSetFile(FILE * fout, dynamicFactorData* data) {

    fprintf(fout, "ActiveSet ");
    for (int i=0;i<data->lenActiveSet;i++) {
        fprintf(fout,"%i ",data->activeSet[i]);
    }
    fprintf(fout,"\n");
    fprintf(fout, "Distribution ");
    double acc=0.0;
    for (int i=0;i<data->lenActiveSet;i++) {
        fprintf(fout,"%.20f ",data->distribution[i]);
        acc+=data->distribution[i];
    }
    fprintf(fout,"Total %f \n",acc);

}

void printMapsFile(const char * s,FactorGraph *fg, FILE *fout) {
    for (int i=0;i<fg->edgeInfo->cEdges; i++) {
        fprintf(fout,"MAPS %s,%i->%f\n",s, i, fg->edgeInfo->qs[i]);
    }
};

void printMapsAVFile(const char * s,FactorGraph *fg, FILE *fout) {
    for (int i=0;i<fg->varInfo->cBinaryVariables; i++) {
        fprintf(fout,"MAPS_AV %s,%i->%.9f\n" , s, i, fg->varInfo->p[i]);
    }
};

void printMapsSumFile(const char * s,FactorGraph *fg, FILE *fout) {
    for (int i=0;i<fg->varInfo->cBinaryVariables; i++) {
        fprintf(fout,"MAPS_SUM %s,%i->%.09f\n" , s, i, fg->varInfo->sumQiAlpha[i]);
    }
};

void printLambdasFile(FactorGraph *fg, FILE *fout) {
    for (int i=0;i<fg->edgeInfo->cEdges; i++) {
        double value=fg->edgeInfo->lamdas[i];
        if (NEARLY_ZERO_TOL(value,1e-9)) value=0.0;
        fprintf(fout,"LAMBDAS %i->%.09f\n",i, value);
    }
};

void printFacOutLogPotential(FactorGraph *fg, int factor_id, FILE *fout) {
    uint32_t varsI = fg->factors[factor_id].sizeI;
    uint32_t varsJ = fg->factors[factor_id].sizeJ;
    for (int prt=0;prt<varsI+varsJ;prt++) {
        fprintf(fout,"factor var output f%i %i %f\n",factor_id,prt,fg->factors[factor_id].outVarLogPotentials[prt]);
    }
};

void printVectorDblFile(const char * s, double * vector, int len , FILE *fout){
    for (int i=0;i<len;i++){
        fprintf(fout,"%s %i %f \n",s,i,vector[i]);
    }
}

void printVectorIntFile(const char * s, int * vector, int len , FILE *fout){
    for (int i=0;i<len;i++){
        fprintf(fout,"%s %i %i \n",s,i,vector[i]);
    }
}

int willProduceSingularMatrix(double* data,uint32_t *r, int len) {
    double s=2.0;
    for (int i = 0; i < len; ++i) {
        if (r[i] == 0.0) continue;
        s -= r[i] * r[i] * data[i * len + i];
        for (int j = i+1; j < len; ++j) {
            if (r[j] == 0.0) continue;
            s -= 2 * r[i] * r[j] * data[i * len + j];
        }
    }
    return (NEARLY_ZERO_TOL(s,1e-6));
}

int countCommon(int a, int b, int sizeJ) {
    if (a==b) return 2;
    int a_col = a % sizeJ;
    int a_row = a / sizeJ;
    int b_col = b % sizeJ;
    int b_row = b / sizeJ;
    if (a_col==b_col) return 1;
    if (a_row==b_row) return 1;
    return 0;

}

void printActiveSetPairs(dynamicFactorData * data, int width) {
    for(int i=0; i<data->lenActiveSet;i++) {
        printf("%i(%i,%i) ",data->activeSet[i],data->activeSet[i]/width, data->activeSet[i] % width);
    }
    printf("\n");
}

void printActiveSetPairs2(dynamicFactorData * data, int width) {
    for(int i=0; i<data->lenActiveSet;i++) {
        printf("(%i,%i) ", data->activeSet[i]/width, data->activeSet[i] % width);
    }
    printf("\n");
}


void solveQP (FactorGraph *fg, int factor_id, double *scratch) {

    Factor *factor =&fg->factors[factor_id];

    double *varLogPotentials = &(fg->edgeInfo->currentLogPotentials[fg->factorInfo->edgesStart[factor_id]]);
    double *currentFactorPotentials = &(fg->factorInfo->currentFactorPotentials[fg->factorInfo->potentialsStart[factor_id]]);
    dynamicFactorData* data = fg->factors[factor_id].dynamicFactorData;
    uint32_t varsI = fg->factors[factor_id].sizeI;
    uint32_t varsJ = fg->factors[factor_id].sizeJ;

    if (data->lenActiveSet == 0) {
        uint32_t max_pos;
        double   max_val;
        findMax(&varLogPotentials[0]     ,varsI,
                &varLogPotentials[varsI] ,varsJ,
                currentFactorPotentials,  &max_pos,&max_val);
        data->inverseA[0]=-2; data->inverseA[1]= 1;
        data->inverseA[2]= 1; data->inverseA[3]= 0;
        data->matrixA[0]= 0;  data->matrixA[1]= 1;
        data->matrixA[2]= 1;  data->matrixA[3]= 2;
        data->activeSet[0]=max_pos;
        data->distribution[0]=1.0;
        data->lenActiveSet=1;
        data->updatedInverseA = true;
    }
    int changed_active_set=true;
    int num_max_iterationsQP = 10;
    for (int qp_iteration=0; qp_iteration<num_max_iterationsQP;qp_iteration++) {
        int same_as_before = true;
        if (changed_active_set) {
            data->varB[0] = 1.0;
            for (int activeSetIterator = 0; activeSetIterator < data->lenActiveSet; activeSetIterator++) {
                double value;
                uint32_t as = data->activeSet[activeSetIterator];
                value = varLogPotentials[as / varsJ];
                value += varLogPotentials[as % varsJ + varsI];
                value += currentFactorPotentials[as];
                data->varB[activeSetIterator + 1] = value;
            }
            if (!data->updatedInverseA)
                computeInverseA(data,scratch);
            matrixMultVector(data);  // updates data->z && data->tau
            same_as_before = false;
        }
        if (same_as_before) {
            computeMarginalsFromSparseDistribution(&fg->factors[factor_id],fg,factor_id);
            uint32_t max_pos;
            double   max_val;
            findMaxMinusExtra(
                    &varLogPotentials[0]     ,varsI,
                    &varLogPotentials[varsI] ,varsJ,
                    currentFactorPotentials,
                    &factor->outVarLogPotentials[0],
                    &factor->outVarLogPotentials[varsI],
                    &max_pos,
                    &max_val);


            int pos_y = max_pos % varsJ;
            int pos_x = max_pos / varsJ;
            if (max_val <= data->tau[0] + 1e-9 /**Small threshold */) {
                return;
            }

            for (int i= 0; i< data->lenActiveSet; i++) {
                if (max_pos==data->activeSet[i]) {
                    data->lenActiveSet=0;
                    return;
                }
            }
            // Adding a new variable, first we check if it fits in memory and we allocate more memory just in case
            uint32_t *r = (uint32_t *) malloc(sizeof(uint32_t)*(data->lenActiveSet+1));
            r[0]=1;
            int max_pos_x = max_pos / varsJ;
            int max_pos_y = max_pos % varsJ;
            for (int i=0; i<data->lenActiveSet;i++) {
                int cx= data->activeSet[i] / varsJ;
                int cy= data->activeSet[i] % varsJ;
                r[i+1] = ((cx==max_pos_x) || (cy==max_pos_y))?1:0;
            }
            if (!data->updatedInverseA) computeInverseA(data, scratch) ;
            int singular=willProduceSingularMatrix(data->inverseA,r,data->lenActiveSet+1) ;
            expandDynamicFactorData(data);
            data->z[data->lenActiveSet] =0.0;
            for (int i=0; i < data->lenActiveSet+1;i++) {
                data->distribution[i] = data->z[i];
            }
            if (singular) {
                int asp2=data->lenActiveSet+2;
                int as=data->lenActiveSet;
                double * padded_similarities=(double *) malloc(sizeof(double)*(asp2)*(asp2));
                int countpadded=asp2*asp2;
                for (int j=0;j<asp2;j++) {
                    for (int i=0;i<asp2;i++) {
                        int pos=j*asp2+i;
                        // (0,0) <- 0
                        if (pos==0) {padded_similarities[pos]=0; continue; }
                        // diagonal <- 2
                        if (i==j)   {padded_similarities[pos] =2; continue; }
                        // first row/column <-1
                        if ((i==0) || (j==0)) {padded_similarities[pos] =1; continue; }
                        // last row/column <- similarity vector (active_set vs max_pos)
                        if (i==asp2-1)  {
                            padded_similarities[pos]=countCommon(data->activeSet[j-1],max_pos,varsJ);
                            continue;
                        }
                        if (j==asp2-1)  {
                            padded_similarities[pos]=countCommon(data->activeSet[i-1],max_pos,varsJ);
                            continue;
                        }
                        // middle positions <- similarity between active_set and active_set
                        padded_similarities[pos]=countCommon(data->activeSet[i-1], data->activeSet[j-1],varsJ);
                    }
                }

                double * eigenvalues = (double*) malloc (sizeof(double)*asp2);
                eigenDecomposition(padded_similarities, eigenvalues,asp2, scratch);

                
                int zero_eigenvalue = -1;
                for (int i = 0; i < asp2; ++i) {
                    if (NEARLY_EQ_TOL(eigenvalues[i], 0.0, 1e-9)) {
                        if (zero_eigenvalue >= 0) {
                            data->lenActiveSet=0;
                            return;
                        }
                        zero_eigenvalue = i;
                    }
                } 
                
                int configuration_to_remove;
                for (int j = 1; j < data->lenActiveSet+1; ++j) {
                    if (zero_eigenvalue*(data->lenActiveSet+2) + j>=countpadded) {
                        printf("Factor %i countpaded %i searching %i\n",factor_id,countpadded,zero_eigenvalue*(data->lenActiveSet+2) + j);
                        printf("zero %i data->len+2 %i j %i \n",zero_eigenvalue,(data->lenActiveSet+2) ,j);

                    }

                    double value = padded_similarities[zero_eigenvalue*(data->lenActiveSet+2) + j];
                    if (!NEARLY_EQ_TOL(value, 0.0, 1e-9)) {
                        configuration_to_remove=(j-1);
                        break;
                    }
                }

#ifdef LAPACK_INVERSION
                removeVarFromEverywhere(data, configuration_to_remove,1,1);
                changed_active_set=true;
#else
                // Version con inversion de martins
                InvertAfterRemovalMar(data->inverseA,data->lenActiveSet+1,configuration_to_remove);
                removeVarFromEverywhere(data, configuration_to_remove,1,1);
                changed_active_set=true;
                data->updatedInverseA=true;
#endif
                free(padded_similarities);
                free(eigenvalues);
                r[0]=1;
                for (int i=0; i<data->lenActiveSet;i++) {
                    r[i+1] = (countCommon(data->activeSet[i],max_pos,varsJ)>0)?1:0;
                }

            }

#ifdef LAPACK_INVERSION
            insertRowColumnIntoMatrixA(data,r);
            data->activeSet[data->lenActiveSet] = max_pos;
            data->lenActiveSet += 1;
            data->updatedInverseA = false;
#else
            data->activeSet[data->lenActiveSet] = max_pos;
            data->lenActiveSet += 1;
            invertMatrixMar(r, data->inverseA, data->lenActiveSet);
            data->updatedInverseA = true;
#endif
            changed_active_set = true;
            free(r);


        } else {

            int exist_blocking = false;
            int blocking  = -1;
            double alpha = 1.0;
            changed_active_set = false;
            for (int i=0; i< data->lenActiveSet;i++) {
                if (data->z[i]< data->distribution[i]) {
                    if (data->z[i] < 0) exist_blocking = true;
                    double tmp = data->distribution[i] / (data->distribution[i] - data->z[i]);
                    if ((blocking < 0) || (tmp < alpha)) {
                        alpha = tmp;
                        blocking = i;
                    }
                }
            }

            if (exist_blocking) {

                if (alpha > 1.0) alpha =1.0;
                for (int i=0; i<data->lenActiveSet;i++) {
                    data->distribution[i] = (1.0-alpha) * data->distribution[i] + alpha*data->z[i];
                    data->z[i] = data->distribution[i];
                }
#ifdef LAPACK_INVERSION
                removeVarFromEverywhere(data, blocking,0,0);
                changed_active_set=true;
#else


                InvertAfterRemovalMar(data->inverseA,data->lenActiveSet+1,blocking);
                removeVarFromEverywhere(data, blocking,0,0);
                changed_active_set=true;
                data->updatedInverseA=true;


#endif
            } else { // !exist_blocking
                for (int i=0; i<data->lenActiveSet;i++) {
                    data->distribution[i] = data->z[i];
                }
            }
        }
    }
    computeMarginalsFromSparseDistribution(&fg->factors[factor_id],fg,factor_id);


}

double getAssignment(FactorGraph *fg) {
    int *choosen = (int *) malloc(sizeof(int)*fg->varInfo->cMultiVariables);

    double sum=0.0;
    for (int mv = 0; mv < fg->varInfo->cMultiVariables; mv++) {
        double best_probability = -1;
        int position=0;

        for (int bv = fg->varInfo->variablesFrom[mv]; bv < fg->varInfo->variablesTo[mv]; bv++) {
            if (fg->varInfo->p[bv] > best_probability) {
                best_probability = fg->varInfo->p[bv];
                position=bv;
            }
        }
        if (best_probability==0.0) {
            double min=-fg->varInfo->logPotentials[fg->varInfo->variablesFrom[mv]];
            for (int bv=fg->varInfo->variablesFrom[mv];bv < fg->varInfo->variablesTo[mv];bv++) {
                if (-fg->varInfo->logPotentials[bv] < min) {
                    position = bv;
                    min = -fg->varInfo->logPotentials[bv];
                }

            }
        }
        choosen[mv] = position -fg->varInfo->variablesFrom[mv];
        sum+=-fg->varInfo->logPotentials[position];
    }

    for (int factor=0;factor< fg->factorInfo->cFactors; factor++){
        int varleft=fg->factors[factor].variableLeft;
        int varright=fg->factors[factor].variableRight;
        int sizeI=fg->factors[factor].sizeI;
        int sizeJ=fg->factors[factor].sizeJ;
        int pos= choosen[varright] + choosen[varleft]*sizeJ;
        sum+= -fg->factorInfo->initialFactorPotentials[fg->factorInfo->potentialsStart[factor] + pos];
    }

    return sum;
}

double getPrimal(FactorGraph *fg) {
    int *choosen = (int *) malloc(sizeof(int)*fg->varInfo->cMultiVariables);

    double sum=0.0;
    for (int mv = 0; mv < fg->varInfo->cMultiVariables; mv++) {
        for (int bv = fg->varInfo->variablesFrom[mv]; bv < fg->varInfo->variablesTo[mv]; bv++) {
            sum+=fg->varInfo->p[bv]*fg->varInfo->logPotentials[bv];
        }
    }

    for (int factor=0;factor< fg->factorInfo->cFactors; factor++){
        int varleft=fg->factors[factor].variableLeft;
        int varright=fg->factors[factor].variableRight;
        int sizeI=fg->factors[factor].sizeI;
        int sizeJ=fg->factors[factor].sizeJ;
        int pos= choosen[varright] + choosen[varleft]*sizeJ;
        sum+= -fg->factorInfo->initialFactorPotentials[fg->factorInfo->potentialsStart[factor] + pos];
    }

    return sum;
}



void solveAD3(FactorGraph * fg, ad3_params * params) {
    
    int    ad3_max_iterations     = params->ad3_max_iterations;
    double ad3_eta                = params->ad3_eta;
    double ad3_residual_threshold = params->ad3_residual_threshold;

#ifndef SEQ
    printf("threads %i\n", omp_get_max_threads());
#endif
    printf("Running AD3\n");
    int eta_changed = true;
    double dual_obj_best = 1e100;
    double primal_rel_obj = -1e100;
    double primal_rel_obj_best = -1e100;
    double primal_obj = 0.0;
    int reached_lower_bound = false;
    int optimal = false;
    double dual_obj=1e100;
    double num_inactive_factors = 0;
    double extra_score=computeExtraScore(fg);

    //extra_score = 0.0;
    printf("extra score %f\n",extra_score);
    timerStart();
    double dual_residual=0.0;
    double primal_residual =0.0;
    int do_break = false;
    int recompute_everything = true;
    int total_iterations =0;
    
    float min_value =1e100;
    int min_pos=0;

    #pragma omp parallel  
    {
    int threads=omp_get_num_threads();    
    double scratch[SCRATCH_SIZE];
    #pragma omp for
    for (int variable=0; variable<fg->varInfo->cBinaryVariables ;variable++) {
        fg->varInfo->binaryVariablesActive[variable] = false;
    }

    int iteration=0;
    for (; iteration<ad3_max_iterations;iteration++ ) {

        #pragma omp barrier
        #pragma omp master 
        {
            recompute_everything =((eta_changed) || ((iteration % ad3_num_iterations_reset) ==0));
            num_inactive_factors = 0;
        }   
        if (eta_changed) {
            #pragma omp for
            for (int edge = 0; edge < fg->edgeInfo->cEdges;edge++) {
                int bv = fg->edgeInfo->binaryVariables[edge];
                int factor = fg->edgeInfo->factors[edge];
                if (!fg->factorInfo->factorsActive[factor]) continue;
                double Xi = fg->edgeInfo->initialLogPotentials[edge]+ 2.0*fg->edgeInfo->lamdas[edge];
                double p_i = fg->varInfo->p[bv];
                fg->edgeInfo->currentLogPotentials[edge] = p_i + Xi / (2*ad3_eta);
            }
            #pragma omp for
            for (int factor=0;factor<fg->factorInfo->cFactors;factor++) {
                if (!fg->factorInfo->factorsActive[factor]) continue;
                for (int flp=fg->factorInfo->potentialsStart[factor]; flp < fg->factorInfo->potentialsEnd[factor]; flp++) {
                    fg->factorInfo->currentFactorPotentials[flp] = fg->factorInfo->initialFactorPotentials[flp] / (2*ad3_eta);
                }
            }
        }


        #pragma omp  for
        for (int factor=0;factor < fg->factorInfo->cFactors;factor++) {
            if ((!recompute_everything) && (!fg->factorInfo->factorsActive[factor])) {
                num_inactive_factors++;
                continue;
            }
            solveQP(fg,factor,scratch);
        }
        #pragma omp barrier

  
        updateEdgeAndVar(fg,recompute_everything,eta_changed,iteration);

        // Deactivate Factors
        #pragma omp  for   
        for (int factor=0; factor<fg->factorInfo->cFactors;factor++) {
            fg->factorInfo->factorsActive[factor] = false;
        }
        #pragma omp  for 
        for (int edge=0; edge<fg->edgeInfo->cEdges;edge++) {
            fg->edgeInfo->edgeFactorActive[edge] = false;
        }
        #pragma omp barrier
        
         // computePi, iterating over binaryVariables

        //// ----
        #pragma omp master
        {
        dual_residual=0.0;
        primal_residual = 0.0;
        }
        #pragma omp barrier
        //#pragma omp master
        #pragma omp for reduction(+:dual_residual)    
        for (int bv=0;bv<fg->varInfo->cBinaryVariables;bv++) {
            int mv= fg->varInfo->binaryBelongsTo[bv];
            double p_prev = fg->varInfo->p[bv];
            if ((fg->varInfo->binaryVariablesActive[bv]) && (fg->varInfo->degree[mv] !=0)) {
                double sumQiAlpha=0.0;
                for (int i=0;i<threads;i++) {
                    sumQiAlpha+=fg->varInfo->sumQiAlpha[bv*threads+i];
                }
                fg->varInfo->p[bv] = (sumQiAlpha)/ ((1.0)* fg->varInfo->degree[mv]);

            }
            double diff =fg->varInfo->p[bv]- p_prev;
            fg->varInfo->diff[bv] = diff;
            dual_residual += fg->varInfo->degree[mv] * diff * diff;
        }

        // computeEdgeLogPotential
        
        //#pragma omp master
        #pragma omp for reduction(+:primal_residual)  
        for (int edge=0; edge < fg->edgeInfo->cEdges; edge++) {
            int bv = fg->edgeInfo->binaryVariables[edge];
            if (!fg->varInfo->binaryVariablesActive[bv]) {
                continue;
            }
            double diff_penalty = fg->edgeInfo->qs[edge] - fg->varInfo->p[bv];
            double diff=fg->varInfo->diff[bv];
            fg->edgeInfo->lamdas[edge] -= ad3_tau * ad3_eta * diff_penalty;
            fg->factorInfo->factorsActive[fg->edgeInfo->factors[edge]] = true;
            fg->edgeInfo->currentLogPotentials[edge] += diff - ad3_tau * diff_penalty;
            if (NEARLY_ZERO_TOL(fg->edgeInfo->currentLogPotentials[edge],1e-12)) {
                fg->edgeInfo->currentLogPotentials[edge] = 0.0;
            }
            primal_residual += diff_penalty*diff_penalty;
        } 
        

        #pragma omp barrier
        #pragma omp master 
        {
           dual_residual= sqrt(dual_residual / fg->edgeInfo->cEdges); 
           primal_residual = sqrt(primal_residual / fg->edgeInfo->cEdges);      
        }        
        
        
        #pragma omp for    
        for (int factor=0; factor<fg->factorInfo->cFactors;factor++) {
            int activated=fg->factorInfo->factorsActive[factor];
            #pragma simd 
            for (int edge=fg->factorInfo->edgesStart[factor]; edge< fg->factorInfo->edgesEnd[factor];edge++ ) {
                fg->edgeInfo->edgeFactorActive[edge] =activated;
            }
        }
        #pragma omp for 
        for (int variable=0; variable<fg->varInfo->cBinaryVariables ;variable++) {
            fg->varInfo->binaryVariablesActive[variable] = false;
        }
       

        #pragma omp barrier
        double compute_dual = ((primal_residual< ad3_residual_threshold) || ((iteration > 0) && (iteration % ad3_num_iterations_compute_dual == 0)));
        #pragma omp master
        dual_obj = dual_obj_best;
        if (compute_dual) {
            #pragma omp master
            dual_obj = 0.0;
            double maxVarLogPotential;
            uint32_t tmpMaxPos;
            #pragma omp barrier
            #pragma omp for reduction (+:dual_obj)
            for (int factor=0;factor<fg->factorInfo->cFactors;factor++) {
                double delta=0;
                for (int edge=fg->factorInfo->edgesStart[factor]; edge< fg->factorInfo->edgesEnd[factor]; edge++) {
                    int bv = fg->edgeInfo->binaryVariables[edge];
                    double varLogPotential = fg->varInfo->logPotentials[bv];
                    int varDegree = fg->varInfo->degree[fg->varInfo->binaryBelongsTo[bv]];
                    double tmpVarLogPotential = varLogPotential / varDegree + (2*fg->edgeInfo->lamdas[edge]);
                    fg->edgeInfo->tmpDataComputeDual[edge]=tmpVarLogPotential;
                    if (tmpVarLogPotential > maxVarLogPotential) {
                        maxVarLogPotential = tmpVarLogPotential;
                    }
                    delta -= fg->edgeInfo->lamdas[edge];
                }

                int sizeI=fg->factors[factor].sizeI;
                int sizeJ=fg->factors[factor].sizeJ;
                int edgeStart = fg->factorInfo->edgesStart[factor];
                findMax(&fg->edgeInfo->tmpDataComputeDual[edgeStart + 0   ],sizeI,
                        &fg->edgeInfo->tmpDataComputeDual[edgeStart +sizeI],sizeJ,
                        &fg->factorInfo->initialFactorPotentials[fg->factorInfo->potentialsStart[factor]],
                        &tmpMaxPos,&maxVarLogPotential);
                dual_obj += maxVarLogPotential + delta;
            }
            #pragma omp master
            {
                dual_obj += extra_score;
                if (dual_obj_best > dual_obj) {
                    dual_obj_best = dual_obj;
                }
                primal_obj = getAssignment(fg);
                if (min_value > primal_obj) {
                    min_value = primal_obj;
                    min_pos=iteration;
                } 
                // primal_rel = getPrimal(fg);

                printf("Iteration=%i Dual obj=%f Dual residual=%f Primal obj=%f Primal residual=%f "
                       "Best dual obj=%f Cached factors=%f eta =%f Changed eta=%s time=%f\n",
                       iteration, dual_obj, dual_residual,primal_obj,primal_residual, dual_obj_best,
                       num_inactive_factors*1.0/fg->factorInfo->cFactors, ad3_eta,
                       (eta_changed)?"True":"False",getRunningTimer() );
            }
        }
        #pragma omp barrier
        
        if (dual_obj_best > dual_obj) {
            #pragma omp master
            dual_obj_best = dual_obj;
            #pragma omp for
            for (int bv=0; bv<fg->varInfo->cBinaryVariables;bv++) {
                fg->varInfo->outProbability[bv] = fg->varInfo->p[bv];
            }
            #pragma omp master 
            if (dual_obj_best < ad3_lower_bound) {
                reached_lower_bound = true;
                do_break=true;
                
            }
        }
        
        #pragma omp barrier
        if ((dual_residual < ad3_residual_threshold) && (primal_residual <ad3_residual_threshold)) {
            #pragma omp for
            for (int bv=0; bv<fg->varInfo->cBinaryVariables;bv++) {
                fg->varInfo->outProbability[bv] = fg->varInfo->p[bv];
            }
            #pragma omp master
            {
                optimal = true;
                do_break=true;
                printf("Exiting residual %5.20f %5.20f \n",dual_residual,primal_residual);
            }
        }

        #pragma omp master 
        {
            eta_changed=false;
            if ((ad3_adapt_eta) && ((iteration % ad3_num_iterations_adapt_eta)==0)) {
                if (primal_residual > ad3_gamma_primal * dual_residual) {
                    if (ad3_eta < ad3_max_eta) {
                        ad3_eta *=ad3_factor_step;
                        eta_changed= true;
                        // TODO : Ver cómo cambia el flujo en base a este cambio de factor. 
                        printf("Iteracion : %i * ad3_factor_step %5.20f %5.20f\n",iteration,primal_residual,dual_residual);
                    }
                } else
                    if (dual_residual> ad3_gamma_dual * primal_residual) {
                        if (ad3_eta > ad3_min_eta) {
                            ad3_eta /= ad3_factor_step;
                            eta_changed = true;
                            printf("Iteracion : %i / ad3_factor_step %5.20f %5.20f\n",iteration,primal_residual,dual_residual);
                        }
                        // TODO : Ver cómo cambia el flujo en base a este cambio de factor. 
                    }
            }
            total_iterations = iteration+1;
        }
        #pragma omp barrier
        if (do_break) {
            iteration= ad3_max_iterations;
        }

        // #pragma omp master
        // printf("%i  "      "%s "                           "%s "                          "%5.16f "          "%5.16f\n",
        //        iteration, (compute_dual)?"True": "False", (eta_changed)?"True": "False", primal_residual,dual_residual); 


    } // Main for

    } // omp parallel
    printf("Exiting residual %5.20f %5.20f \n",dual_residual,primal_residual);
    int fractional=false;
    for (int bv=0; bv<fg->varInfo->cBinaryVariables;bv++) {
        fg->varInfo->outProbability[bv] = fg->varInfo->p[bv];
        if (!NEARLY_BINARY(fg->varInfo->p[bv], 1e-12)) fractional=true;
    }
    if (optimal) {
        if (fractional) printf("The solution is fractional\n");
        else printf("The solution is optimal\n");
    } else {
        printf("The solution is approximate\n");
    }

    printf("Primal objective = %f\nFinished AD3 in %i iterations in %f seconds\n\n", getAssignment(fg),total_iterations, getRunningTimer() );
    printf("Best primal %f at iteration %i\n",min_value,min_pos);
};