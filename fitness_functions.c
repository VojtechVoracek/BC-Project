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

double sphere(Chromosome* ch) {
    double result = 0;
    for (int i = 0; i < DIMENSION; ++i) {
        result += pow(ch->vector[i], 2);
    }

    return result;
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

double Michalewicz(Chromosome* ch) {
    double result = 0;
    for (int i = 0; i < DIMENSION; ++i) {
        result += -sin(ch->vector[i]) * pow(sin((i + 1) * pow(ch->vector[i], 2) / /*M_PI*/ acos(-1)), 20);
    }
    return fabs(result);
}

double Rastrigin(Chromosome* ch) {
    double result = 10 * DIMENSION;
    for (int i = 0; i < DIMENSION; ++i) {
        result += pow(ch->vector[i], 2) - 10 * cos(2 * /*M_PI*/acos(-1) * ch->vector[i]);
    }
    return fabs(result);
}

double Rosenbrock(Chromosome* ch) {
    double result = 0;
    for (int i = 0; i < DIMENSION - 1; ++i) {
        result += 100 * pow(ch->vector[i + 1] - pow(ch->vector[i], 2), 2) + pow(1 - ch->vector[i], 2);
    }
    return result;
}

double OSoREB(Chromosome* ch) {
    double result = SoREB(ch);

    for (int i = size_of_fraction; i < DIMENSION; i += size_of_fraction) {
        if (i >= DIMENSION) break;
        double* vector = makeRotatedVector(ch, i-1, 2);
        double tmp = ellipsoid(vector, 2);
        free(vector);
        result += tmp;
    }

    return result;
}

