//
// Created by Francisco Cruz on 28/12/16.
//

#ifndef AD3NEW_AD3SOLVE_H
#define AD3NEW_AD3SOLVE_H
#include "fgstruct.h"

typedef struct ad3_parameter_def_ {
    uint32_t ad3_max_iterations;
    double   ad3_eta;
    double   ad3_residual_threshold;

} ad3_params;

double getAssignment(FactorGraph *fg);
void solveAD3(FactorGraph * fg, ad3_params *ad3_parameters);




#endif //AD3NEW_AD3SOLVE_H
