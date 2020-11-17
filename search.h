//
// Created by vojta on 11.11.20.
//

#ifndef BCPROJECT_SEARCH_H
#define BCPROJECT_SEARCH_H


#include "parameters.h"

Chromosome uniformMutate(Chromosome* x, Chromosome* a, Chromosome* b, Chromosome* c);

Chromosome informMutate(Chromosome* x, Chromosome* a, Chromosome* b, Chromosome* c, Dependency* dependency);

Minim step(Population* population, Dependency* dependency, int flag);

Statistics DEAlgorithm(Population* population, Dependency* dependency, int flag);


#endif //BCPROJECT_SEARCH_H
