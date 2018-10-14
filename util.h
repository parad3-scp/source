//
// Created by Francisco Cruz on 27/12/16.
//
#include <stdio.h>
#include <sys/time.h>


#ifndef AD3NEW_UTIL_H
#define AD3NEW_UTIL_H

#define NEARLY_EQ_TOL(a,b,tol) (((a)-(b))*((a)-(b))<=(tol))
#define NEARLY_BINARY(a,tol) (NEARLY_EQ_TOL((a),1.0,(tol)) || NEARLY_EQ_TOL((a),0.0,(tol)))
#define NEARLY_ZERO_TOL(a,tol) (((a)<=(tol)) && ((a)>=(-(tol))))

/*** Argument parsing **/
typedef struct _arguments_information_ {
    char *filein;
    FILE *fdin;
    char *fileout;
    FILE *fdout;
    int ad3_max_iterations;
    double ad3_eta;
    double ad3_residual_threshold;
} ArgsInfo;

ArgsInfo * parse(char argc, char**argv);

/**** Measuring time **/
struct timeval  tv1, tv2;
void timerStart();
void timerStop();
double getTimer();
double getRunningTimer();

/*** print helpers ***/
char* readable_fs(double size/*in bytes*/, char *buf);



#endif //AD3NEW_UTIL_H
