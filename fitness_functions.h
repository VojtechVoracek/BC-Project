//
// Created by vojta on 10.11.20.
//

#ifndef BCPROJECT_FITNESS_FUNCTIONS_H
#define BCPROJECT_FITNESS_FUNCTIONS_H

#include "parameters.h"

double oneMax(Chromosome* ch);

int size_of_fraction;

double ellipticFunction(int start, Chromosome* ch);

double myBenchMarkFunction(Chromosome* ch);

#endif //BCPROJECT_FITNESS_FUNCTIONS_H
