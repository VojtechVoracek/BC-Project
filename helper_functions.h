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

double** allocMatrix(int size);

double** makeOriginRotationMatrix(int dimension, double angle);

double** makeTranslationMatrix(int dimension, double* point);

double** makeInverseTranslationMatrix(int dimension, double* point);

double** makeRotationMatrix(int dimension, double angle, double* point);

double** make1D (double angle);

double** make2D (double angle);

double** make3D (double angle);

double** make4D (double angle);

double** make5D(double angle);

void makeMatrices(double angle);

void printSquareMatrix(double** M, int dimension, char* name);

void freeMatrices();

void freeSquareMatrix(double** M, int dimension);

Node* makeSubTree(int from, int to);

void printTree(Node* root);

Mask* makeArrayFromTree();

void freeArrayFromTree(Node* root);

void masksFromTree(Node* root, Mask** mask_array);

void freeMaskArray(Mask* m);

int* individuals2Learning (Population * p);

void freeFOS(FOS fos, double** D);

double** createDependencyMatrix(Population* population);

double** createMIMatrix(Population* population);

void printFOS(FOS fos);

int* findNBest(Population* population, int n);

double** createMICMatrix(Population* population);

double** getOptimalDMatrix();

#endif //BCPROJECT_HELPER_FUNCTIONS_H
