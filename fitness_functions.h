//
// Created by vojta on 10.11.20.
//

#ifndef BCPROJECT_FITNESS_FUNCTIONS_H
#define BCPROJECT_FITNESS_FUNCTIONS_H

#include "parameters.h"

double oneMax(Chromosome* ch);

double sphere(Chromosome* ch);

int size_of_fraction;

double* makeRotatedVector(Chromosome* ch, int start, int size);

double ellipsoid(double* vector, int size);

double SoREB(Chromosome* ch);

double Michalewicz(Chromosome* ch);

double Rastrigin(Chromosome* ch);

double Rosenbrock(Chromosome* ch);

double OSoREB(Chromosome* ch);

#endif //BCPROJECT_FITNESS_FUNCTIONS_H
