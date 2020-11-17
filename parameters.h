//
// Created by vojta on 10.11.20.
//

#ifndef BCPROJECT_PARAMETERS_H
#define BCPROJECT_PARAMETERS_H


#include <stdbool.h>

extern int DIMENSION;
extern double CR;
extern double F;
extern double (*func)(void*);
extern double TERMINAL_CONDITION;

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


#endif //BCPROJECT_PARAMETERS_H
