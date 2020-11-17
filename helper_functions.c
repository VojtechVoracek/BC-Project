//
// Created by vojta on 10.11.20.
//

#include "helper_functions.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double RandomGenerator() {
    return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
}

double normalRandom() {
    double y1=RandomGenerator();
    double y2=RandomGenerator();
    return cos(2*3.14*y2)*sqrt(-2.*log(y1));
}

double* randomSamples(int size) {
    double* samples = (double *) malloc(size * sizeof(double ));
    for (int i = 0; i < size; ++i) {
        samples[i] = normalRandom();
    }

    return samples;
}

Chromosome initChromosome() {
    Chromosome ch;
    ch.vector = (double *) malloc(DIMENSION * sizeof(double));
    for (int i = 0; i < DIMENSION; ++i) {
        ch.vector[i] = normalRandom() * 100;
    }
    ch.fitness = func(&ch);

    return ch;
}

Population initPopulation(int size) {
    Population p;
    p.size = size;
    p.individuals = (Chromosome*) malloc(size * sizeof(Chromosome));
    for (int i = 0; i < size; ++i) {
        p.individuals[i] = initChromosome();
    }

    return p;
}

int* pickThree(int size, int my_idx) {

    int* triplet = (int*) malloc(3 * sizeof(int));
    int a;

    while ((a = rand() % size) == my_idx);
    triplet[0] = a;
    while ((a = rand() % size) == my_idx || a == triplet[0]);
    triplet[1] = a;
    while ((a = rand() % size) == my_idx || (a == triplet[0]) || (a == triplet[1]));
    triplet[2] = a;

    if (triplet[0] == triplet[1] || triplet[1] == triplet[2] || triplet[2] == triplet[0]) {
        printf("Chyba trojic!!!!!\n");
    }

    return triplet;

}

Statistics initStatistics(Population* population) {
    Statistics ret;
    Minim minimum;
    minimum.value = 999999999;
    ret.num_of_evaluations = population->size;
    ret.success = false;

    for (int i = 0; i < population->size; ++i) {
        if (population->individuals[i].fitness < minimum.value) {
            minimum.value = population->individuals[i].fitness;
            minimum.index = i;
        }
    }
    ret.minimum = minimum;
    if (ret.minimum.value <= TERMINAL_CONDITION) ret.success = true;

    return ret;
}

void printChromosome(Chromosome ch) {
    for (int i = 0; i < DIMENSION; ++i) {
        printf("%5f ", ch.vector[i]);
    }
    printf("%5f\n", ch.fitness);
}

void printPopulation(Population* p) {
    printf("---------------------------------------------\n");
    for (int i = 0; i < p->size; ++i) {
        printChromosome(p->individuals[i]);
    }
    printf("---------------------------------------------\n");
}

void printStatistics(Statistics* stats, Population* population) {
    if (stats->success) printf("Successfully finished, ");
    else printf("Unsuccessfully finished, ");
    printf("number of evaluations: %d, ", stats->num_of_evaluations);
    printf("minimal value: %5f, vector: ", stats->minimum.value);
    printChromosome(population->individuals[stats->minimum.index]);
    printf("\n");

}

bool equalPopulation(Population* population) {
    Chromosome first = population->individuals[0];

    double epsilon = 0.00000001;

    for (int i = 1; i < population->size; ++i) {
        for (int j = 0; j < DIMENSION; ++j) {
            if (fabs(first.vector[j] - population->individuals[i].vector[j]) > epsilon) return false;
        }
    }

    return true;
}

void makeChanges(Change* changes, Population* population) {

    for (int i = 0; i < changes->size; ++i) {
        free(population->individuals[changes[i].index].vector);
        population->individuals[changes[i].index] = changes[i].chromosome;
    }
}

int findMax(int* vector) {

    int maximum = -1;
    for (int i = 0; i < DIMENSION; ++i) {
        if (vector[i] > maximum) maximum = vector[i];
    }

    return maximum;
}

int* findSizeOfClasses(int* vector, int num_of_classes) {
    int* result = (int*) calloc(num_of_classes, sizeof(int));
    for (int i = 0; i < DIMENSION; ++i) {
        result[vector[i]]++;
    }
    return result;
}

int* makeCumulative(int* sizes, int num_of_classes) {
    int* cumulative = (int*) malloc(num_of_classes * sizeof(int));
    cumulative[0] = 0;
    for (int i = 1; i < num_of_classes; ++i) {
        cumulative[i] = cumulative[i-1] + sizes[i-1];
    }
    return cumulative;
}

Dependency makeDependencies(int* vector) {

    Dependency result;
    result.num_of_classes = findMax(vector) + 1;
    result.size_of_classes = findSizeOfClasses(vector, result.num_of_classes);
    result.sorted_indexes = (int*) malloc(DIMENSION * sizeof(int));

    int* offsets = (int*) calloc(result.num_of_classes, sizeof(int));
    int* cumulative = makeCumulative(result.size_of_classes, result.num_of_classes);


    for (int i = 0; i < DIMENSION; ++i) {
        int index = cumulative[vector[i]] + offsets[vector[i]]++;
        result.sorted_indexes[index] = i;
    }

    free(cumulative);
    free(vector);
    free(offsets);

    return result;
}

void freePopulation(Population population) {
    for (int i = 0; i < population.size; ++i) {
        free(population.individuals[i].vector);
    }
    free(population.individuals);
}

