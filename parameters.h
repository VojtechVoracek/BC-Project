//
// Created by vojta on 10.11.20.
//

#ifndef BCPROJECT_PARAMETERS_H
#define BCPROJECT_PARAMETERS_H

#define Matrix double*

#include <stdbool.h>

extern int DIMENSION;
extern double CR;
extern double F;
extern double (*func)(void*);
extern double TERMINAL_CONDITION;

int size_of_fraction;

//radians
double angle;

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
