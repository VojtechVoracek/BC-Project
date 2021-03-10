//
// Created by vojta on 11.11.20.
//

#ifndef BCPROJECT_SEARCH_H
#define BCPROJECT_SEARCH_H


#include "parameters.h"

Chromosome uniformMutate(Chromosome* x, Chromosome* a, Chromosome* b, Chromosome* c);

Chromosome fosMutate(Chromosome* x, Chromosome* a, Chromosome* b, Population* p, FOS fos);

void ltMutate(int index, Chromosome* x, Population* p, Mask* traversal_array);

Minim step(Population* population, int flag, FOS* fos);

FOS initFos(int flag);

FOS learnFos(double** D, int flag);

Statistics DEAlgorithm(Population* population, int fos_flag, int dist_flag);

#endif //BCPROJECT_SEARCH_H
