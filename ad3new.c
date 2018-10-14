#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "fgstruct.h"
#include "util.h"
#include "ad3solve.h"
#include "dynamicfd.h"
#include "omp.h"



FactorGraph * createFactorGraph() {
    return (FactorGraph *) malloc(sizeof(FactorGraph));

}

void fileImport(FILE * file, FactorGraph *fg, int threads) {

    int printinfo=0;


    fg->varInfo = (VarInfo*) malloc (sizeof(VarInfo));
    uint32_t *cMultiVariables = &fg->varInfo->cMultiVariables;
    uint32_t *cBinaryVariables = &fg->varInfo->cBinaryVariables;

    fread(cMultiVariables, int_size,1,file);

    *cBinaryVariables =0;


    fg->varInfo->variablesFrom      = (uint32_t*  ) malloc((*cMultiVariables)*sizeof(uint32_t));
    fg->varInfo->variablesTo        = (uint32_t*  ) malloc((*cMultiVariables)*sizeof(uint32_t));
    fg->varInfo->degree             = (uint32_t*  ) malloc((*cMultiVariables)*sizeof(uint32_t));
    fg->varInfo->connectedFactors   = (uint32_t** ) malloc((*cMultiVariables)*sizeof(uint32_t*));


    FileReadBuffer * bufferVars = (FileReadBuffer * ) malloc((*cMultiVariables)*sizeof(FileReadBuffer));
    for (int idx=0;idx<*cMultiVariables;idx++) {
        fread(&(bufferVars[idx].id), int_size,1,file);
        fread(&(bufferVars[idx].n_configurations), int_size,1,file);
        bufferVars[idx].configurations=(double * ) malloc(bufferVars[idx].n_configurations*sizeof(double));
        fread(bufferVars[idx].configurations,sizeof(double),bufferVars[idx].n_configurations,file);
        fg->varInfo->variablesFrom[idx]=*cBinaryVariables;
        *cBinaryVariables+=bufferVars[idx].n_configurations;
        fg->varInfo->variablesTo[idx]=*cBinaryVariables;
        fg->varInfo->degree[idx] = 0;
    }
    printf("Threads %i\n ---",threads);

    fg->varInfo->logPotentials  = (double *) malloc((*cBinaryVariables)*sizeof(double));
    fg->varInfo->diff           = (double *) malloc((*cBinaryVariables)*sizeof(double));
    fg->varInfo->p              = (double *) malloc((*cBinaryVariables)*sizeof(double));
    fg->varInfo->sumQiAlpha     = (double *) malloc((*cBinaryVariables)*sizeof(double)*threads);
    fg->varInfo->outProbability = (double *) malloc((*cBinaryVariables)*sizeof(double));
    fg->varInfo->binaryVariablesActive = (uint32_t *) malloc((*cBinaryVariables)*sizeof(uint32_t));
    fg->varInfo->binaryBelongsTo = (uint32_t *) malloc((*cBinaryVariables)*sizeof(uint32_t));

    for (int i=0;i<*cBinaryVariables;i++) {
        fg->varInfo->diff[i] = 0.0;
        fg->varInfo->p[i] = 0.5;
        
        fg->varInfo->outProbability[i] = 0.0;
        fg->varInfo->binaryVariablesActive[i] = 1; // True
        for (int j=0;j<threads;j++) {
            fg->varInfo->sumQiAlpha[i*threads+j] = 0.0;
        }
    }

    int i=0;
    for (unsigned int mv=0;mv<*cMultiVariables;mv++) {
        for (int bv=0;bv<bufferVars[mv].n_configurations;bv++) {
            fg->varInfo->logPotentials[i] = -bufferVars[mv].configurations[bv];
            fg->varInfo->binaryBelongsTo[i] = mv;
            i++;
        }
        free(bufferVars[mv].configurations);
    }
    free(bufferVars);


    fg->factorInfo = (FactorInfo*) malloc(sizeof(FactorInfo));
    uint32_t *cFactors = &fg->factorInfo->cFactors;
    fread(cFactors, int_size,1,file);

    fg->factors                     =  (Factor * )  malloc ((*cFactors)*sizeof(Factor));
    fg->factorInfo->potentialsStart =  (uint32_t *) malloc ((*cFactors)*sizeof(uint32_t));
    fg->factorInfo->potentialsEnd   =  (uint32_t *) malloc ((*cFactors)*sizeof(uint32_t));
    fg->factorInfo->edgesStart      =  (uint32_t *) malloc ((*cFactors)*sizeof(uint32_t));
    fg->factorInfo->edgesEnd        =  (uint32_t *) malloc ((*cFactors)*sizeof(uint32_t));
    fg->factorInfo->factorsActive   =  (uint32_t *) malloc ((*cFactors)*sizeof(uint32_t));
    fg->factorInfo->cFactorPotentials =0;
    FileReadBuffer * bufferFactors = (FileReadBuffer*) malloc ((*cFactors)*sizeof(FileReadBuffer));
    for (int idx=0;idx<*cFactors;idx++) {
        int varleft,varright;
        fread(&varleft, int_size,1,file);
        fread(&varright, int_size,1,file);
        //varleft-=1;
        //varright-=1;

        fg->factors[idx].variableLeft = (uint32_t) varleft;
        fg->factors[idx].variableRight = (uint32_t) varright;
        fg->varInfo->degree[varleft]++;
        fg->varInfo->degree[varright]++;
        fread(&fg->factors[idx].sizeI, int_size,1,file);
        fread(&fg->factors[idx].sizeJ, int_size,1,file);
        fg->factors[idx].outVarLogPotentials = (double*) malloc(sizeof(double)*(fg->factors[idx].sizeI+fg->factors[idx].sizeJ));
        unsigned int nconf = fg->factors[idx].sizeI*fg->factors[idx].sizeJ;
        bufferFactors[idx].n_configurations = nconf;
        bufferFactors[idx].configurations = (double*) malloc(sizeof(double)*nconf);
        fg->factorInfo->factorsActive[idx] = 1;  // true
        fg->factorInfo->potentialsStart[idx] = fg->factorInfo->cFactorPotentials;
        for (int lp=0; lp<nconf;lp++) {
            double rdouble;
            fread(&rdouble,sizeof(double),1,file);
            bufferFactors[idx].configurations[lp] = rdouble;
            fg->factorInfo->cFactorPotentials++;
        }
        fg->factorInfo->potentialsEnd[idx] = fg->factorInfo->cFactorPotentials;
        fg->factors[idx].dynamicFactorData=createDynamicFactorData();


    }

    fg->factorInfo->initialFactorPotentials = (double *) malloc((fg->factorInfo->cFactorPotentials)*sizeof(double));
    fg->factorInfo->currentFactorPotentials = (double *) malloc((fg->factorInfo->cFactorPotentials)*sizeof(double));
    fg->factorInfo->outputFactorPotentials  = (double *) malloc((fg->factorInfo->cFactorPotentials)*sizeof(double));
    i=0;
    for (unsigned int f=0;f<*cFactors;f++) {
        for (int conf=0;conf<bufferFactors[f].n_configurations;conf++) {
            fg->factorInfo->initialFactorPotentials[i] = -bufferFactors[f].configurations[conf];
            i++;
        }
        free(bufferFactors[f].configurations);
    }
    free(bufferFactors);

    /**** create factor-variable->connection*/
    uint32_t countEdges=0;
    for (i=0;i<fg->varInfo->cMultiVariables;i++) {
        countEdges += fg->varInfo->degree[i] * (fg->varInfo->variablesTo[i] - fg->varInfo->variablesFrom[i]);
        fg->varInfo->connectedFactors[i]= (uint32_t *) malloc(sizeof(uint32_t)*fg->varInfo->degree[i]);
        fg->varInfo->degree[i]=0;
    }


    for (i=0;i<fg->factorInfo->cFactors;i++) {

        uint32_t varleft = fg->factors[i].variableLeft;
        uint32_t varright = fg->factors[i].variableRight;
        fg->varInfo->connectedFactors[varleft][fg->varInfo->degree[varleft]]= (uint32_t) i;
        fg->varInfo->connectedFactors[varright][fg->varInfo->degree[varright]]= (uint32_t) i;
        fg->varInfo->degree[varleft] ++;
        fg->varInfo->degree[varright] ++;
    }

    /*** Precompute probabilites for isolated vars **/
    for (int bv=0;bv<fg->varInfo->cBinaryVariables;bv++) {
        int vardegree= fg->varInfo->degree[fg->varInfo->binaryBelongsTo[bv]];
        if (vardegree==0) {
            if (fg->varInfo->logPotentials[bv] > 0) fg->varInfo->p[bv] =1.0;
            else                                    fg->varInfo->p[bv] =0.0;
        }
    }
    /*** Compute Edge information */
    fg->edgeInfo = (EdgeInfo*) malloc (sizeof(EdgeInfo));
    fg->edgeInfo->cEdges = countEdges;
    fg->edgeInfo->factors= (uint32_t *) malloc(sizeof(uint32_t)*countEdges);
    fg->edgeInfo->multiVariables= (uint32_t *) malloc(sizeof(uint32_t)*countEdges);
    fg->edgeInfo->binaryVariables= (uint32_t *) malloc(sizeof(uint32_t)*countEdges);
    fg->edgeInfo->lamdas = (double *) malloc(sizeof(double)*countEdges);
    fg->edgeInfo->qs= (double *) malloc(sizeof(double)*countEdges);
    fg->edgeInfo->initialLogPotentials= (double *) malloc(sizeof(double)*countEdges);
    fg->edgeInfo->currentLogPotentials= (double *) malloc(sizeof(double)*countEdges);
    fg->edgeInfo->tmpDataComputeDual = (double *) malloc(sizeof(double)*countEdges);
    fg->edgeInfo->edgeFactorActive = (int  *) malloc(sizeof(int)*countEdges);
    fg->edgeInfo->logPotentialInFactor = (double **) malloc(sizeof (double*)*countEdges);
    uint32_t counted=0;
    for (i=0;i<fg->factorInfo->cFactors;i++) {
        fg->factorInfo->edgesStart[i] = counted;
        uint32_t mv;
        for (int leftright=0; leftright<2; leftright++) {

            if (leftright==0)   mv=fg->factors[i].variableLeft;
            else                mv=fg->factors[i].variableRight;
            uint32_t binaryVarStart = fg->varInfo->variablesFrom[mv];
            uint32_t binaryVarEnd =   fg->varInfo->variablesTo[mv];
            uint32_t varDegree=       fg->varInfo->degree[mv];
            for (uint32_t bv=binaryVarStart;bv<binaryVarEnd;bv++) {
                fg->edgeInfo->factors[counted]= (uint32_t) i;
                fg->edgeInfo->multiVariables[counted]=mv;
                fg->edgeInfo->binaryVariables[counted]=bv;
                fg->edgeInfo->lamdas[counted]=0.0;
                fg->edgeInfo->qs[counted]=0.0;
                fg->edgeInfo->currentLogPotentials[counted]=0.0;
                fg->edgeInfo->initialLogPotentials[counted]=(fg->varInfo->logPotentials[bv] / varDegree);


                int offset = (leftright==0)?0:fg->factors[i].sizeI;
                int posInFactor = bv- fg->varInfo->variablesFrom[mv] +offset;

                fg->edgeInfo->logPotentialInFactor[counted] = &fg->factors[i].outVarLogPotentials[posInFactor];


                counted += 1;
            }
        }
        fg->factorInfo->edgesEnd[i] = counted;
    }

    fread(&fg->optimal_energy,sizeof(double),1,file);
    fg->optimal_configuration = (uint32_t *) malloc (sizeof(uint32_t) * fg->varInfo->cMultiVariables);
    fread(fg->optimal_configuration, int_size, fg->varInfo->cMultiVariables, file);
    for(int idx=0; idx<fg->varInfo->cMultiVariables; idx++) {
            fg->optimal_configuration[idx]-=1;
    }
}

void showCounts(FactorGraph * fg) {

    int accDegree=0;
    double avgDegree;
    for(int mv=0;mv<fg->varInfo->cMultiVariables;mv++)
        accDegree += fg->varInfo->degree[mv];

    printf("Factors: %i\n", fg->factorInfo->cFactors);
    printf("Edges %i\n", fg->edgeInfo->cEdges);
    printf("Sum (Factor Values) %i\n", fg->factorInfo->cFactorPotentials);
    printf("Multi Variables %i\n", fg->varInfo->cMultiVariables);
    printf("Binary Variables %i\n", fg->varInfo->cBinaryVariables);
    printf("Variable Average Degree %d\n", accDegree/fg->varInfo->cMultiVariables);

    int memoryOnVars =  sizeof(uint32_t) * fg->varInfo->cMultiVariables * 4 +
                        sizeof(double) * fg->varInfo->cBinaryVariables * 5 +
                        sizeof(uint32_t) * fg->varInfo->cBinaryVariables * 2+
                        sizeof(uint32_t) * accDegree;

    int memoryOnFactors = sizeof(uint32_t) * fg->factorInfo->cFactors * 2 +
                          sizeof(double) * fg->factorInfo->cFactorPotentials * 2 +
                          fg->factorInfo->cFactors * countAllocatedMemoryDynamicFD();

    int memoryOnEdges =     sizeof(double) * fg->edgeInfo->cEdges * 5 +
                            sizeof(uint32_t) * fg->edgeInfo->cEdges * 3;

    int memoryTotal= memoryOnVars+memoryOnFactors+memoryOnEdges;

    char printingTotal[20],printingVars[20],printingFactors[20],printingEdges[20];
    printf("Memory consumption : %s\n\tVars:    %s\n\tFactors: %s\n\tEdges:   %s\n",
           readable_fs(memoryTotal,printingTotal),
           readable_fs(memoryOnVars,printingVars),
           readable_fs(memoryOnFactors,printingFactors),
           readable_fs(memoryOnEdges,printingEdges));

}


int main (int argc, char **argv) {

    ad3_params ad3_parameters;

    ArgsInfo *args=parse(argc,argv);
    FactorGraph *fg = createFactorGraph();
    // Parse file
    struct timeval  tv1, tv2;
    int threads=omp_get_max_threads();
    timerStart();
    fileImport(args->fdin, fg, threads);
    timerStop();
    showCounts(fg);
    fclose(args->fdin);
    printf("File loaded %f secs\n\n",getTimer());

    ad3_parameters.ad3_max_iterations     = args->ad3_max_iterations;
    ad3_parameters.ad3_eta                = args->ad3_eta;
    ad3_parameters.ad3_residual_threshold = args->ad3_residual_threshold;
    solveAD3(fg,&ad3_parameters);
    printf("Minimal energy found %f\n",getAssignment(fg));


    return 0;
}