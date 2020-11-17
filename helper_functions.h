//
// Created by vojta on 10.11.20.
//

#ifndef BCPROJECT_HELPER_FUNCTIONS_H
#define BCPROJECT_HELPER_FUNCTIONS_H

#include <stdio.h>
#include "parameters.h"

double RandomGenerator();

double normalRandom();

double* randomSamples(int size);

Chromosome initChromosome();

Population initPopulation(int size);

int* pickThree(int size, int my_idx);

Statistics initStatistics(Population* population);

void printChromosome(Chromosome ch);

void printPopulation(Population* p);

void printStatistics(Statistics* stats, Population* population);

bool equalPopulation(Population* population);

void makeChanges(Change* changes, Population* population);

int findMax(int* vector);

int* findSizeOfClasses(int* vector, int num_of_classes);

int* makeCumulative(int* sizes, int num_of_classes);

Dependency makeDependencies(int* vector);

void freePopulation(Population population);


#endif //BCPROJECT_HELPER_FUNCTIONS_H
