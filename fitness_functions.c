//
// Created by vojta on 10.11.20.
//

#include "fitness_functions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double oneMax(Chromosome* ch) {
    double ret = 0;
    for (int i = 0; i < DIMENSION; ++i) {
        ret += ch->vector[i];
    }

    return fabs(ret);
}

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

double* makeRotatedVector(Chromosome* ch, int start, int size) {

    double* result = (double*) malloc(size * sizeof(double));

    double** matrix;

    switch (size) {
        case 1:
            matrix = matrix_1D;
            break;
        case 2:
            matrix = matrix_2D;
            break;
        case 3:
            matrix = matrix_3D;
            break;
        case 4:
            matrix = matrix_4D;
            break;
        case 5:
            matrix = matrix_5D;
            break;
        default:
            return -1;
    }

    for (int i = 0; i < size; ++i) {
        double tmp_res = 0;
        for (int j = start; j < start+size; ++j) {
            tmp_res += matrix[i][j-start] * ch->vector[j];
        }
        result[i] = tmp_res;
    }

    return result;
}

double ellipsoid(double* vector, int size) {
    double result = 0;
    for (int i = 0; i < size; ++i) {
        double tmp = pow(vector[i],2);
        tmp *= pow(10, ((double )(6*i)) / ((double) (size-1)));
        result += tmp;
    }

    //result *= pow(10, ((double )(6*size)) / ((double) (size - 1)));
    return result;
}

double SoREB(Chromosome* ch) {

    double result = 0;

    for (int i = 0; i < DIMENSION; i += size_of_fraction) {
        int real_size = size_of_fraction;
        if (i + real_size >= DIMENSION) real_size = DIMENSION - i;
        double* vector = makeRotatedVector(ch, i, real_size);
        //for (int j = 0; j < real_size; ++j) {
        //    printf("%f, ", vector[j]);
        //}
        //printf("--- rotated vector\n");
        //printf("%d\n", real_size);
        double tmp = ellipsoid(vector, real_size);
        free(vector);
        //printf("Vysledek od sud: %f\n", tmp);
        result += tmp;
    }
    //printf("OVerall: %5f\n", result);
    return result;

}







