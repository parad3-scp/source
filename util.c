//
// Created by Francisco Cruz on 28/12/16.
//

#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <ctype.h>
#include <memory.h>
#include <sys/time.h>


ArgsInfo* parse(char argc, char**argv) {

    ArgsInfo *args = (ArgsInfo * ) malloc(sizeof(ArgsInfo));
    args->filein = (char *) malloc(sizeof(char)*255);
    char *fvalue = NULL;
    char *nvalue = NULL;

    int index;
    int c;

    opterr = 0;
    
    args->ad3_max_iterations=1000;
    args->ad3_eta=0.1;
    args->ad3_residual_threshold=1e-6;

    while ((c = getopt(argc, argv, "f:n:e:r:h")) != -1)
        switch (c) {
            case 'f':
                fvalue = optarg;
                break;
            case 'n':
                args->ad3_max_iterations = atoi(optarg);
                if (args->ad3_max_iterations==0) args->ad3_max_iterations=1000;
                break;

            case 'e':
                args->ad3_eta = atof(optarg);
                if (args->ad3_eta==0.0) args->ad3_eta=0.1;
                break;

            case 'r':
                args->ad3_residual_threshold = atof(optarg);
                if (args->ad3_residual_threshold==0.0) args->ad3_residual_threshold=1e-6;
                break;

            case 'h':
                printf("Usage:\n %s [OPTIONS] -f file.in  \n", argv[0]);
                printf("\nAvailable OPTIONS:\n");
                printf("\t-h        : Print this help and exit\n");
                printf("\t-n number : Number of AD3 iterations (Default: 1000)\n");
                printf("\t-e number : Initial ETA (Default: 0.1)\n");
                printf("\t-r number : Initial Residual threshold (Default: 1e-6)\n");
                printf("\n");    
                
                exit(0);    
            case '?':
                if (optopt == 'c')
                    fprintf(stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint(optopt))
                    fprintf(stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf(stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                exit(1);
            default:
                abort();
        }

    if (fvalue == NULL) {
        printf("Please specify a file name %s -f <filename>\n ", argv[0]);
        exit(0);
    }

    // Open file
    args->fdin  = fopen(fvalue, "r");
    strcpy(args->filein,fvalue);

    if (args->fdin == NULL) {
        printf("error opening file  %s\n", fvalue);
        exit(1);
    }

    return args;
}

void timerStart() {
    gettimeofday(&tv1, NULL);
};
void timerStop() {
    gettimeofday(&tv2, NULL);
}
double getRunningTimer() {
    struct timeval  tv2;
    gettimeofday(&tv2, NULL);
    return (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
}

double getTimer() {
    return (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
}

char* readable_fs(double size/*in bytes*/, char *buf) {
    int i = 0;
    const char* units[] = {"B", "kB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"};
    while (size > 1024) {
        size /= 1024;
        i++;
    }
    sprintf(buf, "%.*f %s", i, size, units[i]);
    return buf;
}