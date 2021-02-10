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

double* makeRotatedVector(Chromosome* ch, int start, int size);

double ellipsoid(double* vector, int size);

double SoREB(Chromosome* ch);

#endif //BCPROJECT_FITNESS_FUNCTIONS_H
