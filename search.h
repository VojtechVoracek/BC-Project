//
// Created by vojta on 11.11.20.
//

#ifndef BCPROJECT_SEARCH_H
#define BCPROJECT_SEARCH_H


#include "parameters.h"

Chromosome uniformMutate(Chromosome* x, Chromosome* a, Chromosome* b, Chromosome* c);

Chromosome mpMutate(Chromosome* x, Chromosome* a, Chromosome* b, Chromosome* c, Dependency* dependency);

void ltMutate(int index, Chromosome* x, Population* p, Mask* traversal_array);

Minim step(Population* population, Dependency* dependency, int flag, Mask* traversal_array);

Statistics DEAlgorithm(Population* population, Dependency* dependency, int flag, Mask* traversal_array);


#endif //BCPROJECT_SEARCH_H
