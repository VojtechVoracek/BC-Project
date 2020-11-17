//
// Created by vojta on 10.11.20.
//

#include "fitness_functions.h"
#include <math.h>
#include <stdio.h>

double oneMax(Chromosome* ch) {
    double ret = 0;
    for (int i = 0; i < DIMENSION; ++i) {
        ret += ch->vector[i];
    }

    return fabs(ret);
}

int size_of_fraction = 4;

double ellipticFunction(int start, Chromosome* ch) {
    if (start >= DIMENSION) return 0;

    double sum = 0;
    double product = 1;
    for (int i = start; i < start + size_of_fraction; ++i) {
        if (i >= DIMENSION) break;
        //sum += ch->vector[i] / 10;

        product *= ch->vector[i];
    }
    //return fabs(product);
    return fabs(sum + product);
}

double myBenchMarkFunction(Chromosome* ch) {

    int num_of_fractions = DIMENSION / size_of_fraction;

    double result = 0;

    for (int i = 0; i < num_of_fractions + 1; ++i) {
        result += ellipticFunction(i * size_of_fraction, ch);
    }
    return fabs(result);
}