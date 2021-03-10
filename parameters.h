//
// Created by vojta on 10.11.20.
//

#ifndef BCPROJECT_PARAMETERS_H
#define BCPROJECT_PARAMETERS_H

#define Matrix double*

#include <stdbool.h>

#define UNIFORM 0
#define LT_FOS 1
#define MP_FOS 2
#define LT_FOS_f 3
#define MP_FOS_f 4

#define MIC 1
#define LINC 2

#define LT_FOS_MIC 11
#define LT_FOS_LINC 12
#define MP_FOS_MIC 21
#define MP_FOS_LINC 22

#define LT_FOS_f_MIC 31
#define LT_FOS_f_LINC 32
#define MP_FOS_f_MIC 41
#define MP_FOS_f_LINC 42


extern int DIMENSION;
extern double CR;
extern double F;
extern double (*func)(void*);
extern double TERMINAL_CONDITION;

int size_of_fraction;


double** matrix_1D;
double** matrix_2D;
double** matrix_3D;
double** matrix_4D;
double** matrix_5D;

typedef struct Chromosome_t {
    double* vector;
    double fitness;
}Chromosome;

typedef struct Population_t {
    Chromosome* individuals;
    int size;
} Population;

typedef struct Change_t {
    int index;
    Chromosome chromosome;
    int size;
} Change;

typedef struct Minim_t {
    double value;
    int index;
} Minim;

typedef struct Statistics_t {
    bool success;
    int num_of_evaluations;
    Minim minimum;
} Statistics;

typedef struct Dependency_t {
    int* sorted_indexes;
    int num_of_classes;
    int* size_of_classes;
} Dependency;

typedef struct subset_t {
    int number_of_variables;
    int* indexes;
} subset;

typedef struct FOS_t {
    int number_of_subsets;
    subset* parts;
} FOS;

typedef struct Mack_t {
    int* variables;
    int lenght;
} Mask;

typedef struct Node_t {
    Mask information;
    struct Node_t* left;
    struct Node_t* right;
} Node;


#endif //BCPROJECT_PARAMETERS_H
