//
// Created by voracvoj on 10.03.21.
//

#ifndef BCPROJECT_COMPARISON_H
#define BCPROJECT_COMPARISON_H

#include <stdbool.h>
#include <stdio.h>

bool performN(int n, int size, FILE* f, int fos_flag, int dist_flag);


int findOptimalSize(FILE* f, int num_of_iterations, int fos_flag, int dist_flag);


int*** findOptimalPopulationSize(void** problems, int number_of_problems, int* dimensions, int number_of_dimension, int number_of_iterations, char* name_of_file);


double findMean(int iterations, int size_of_population, int fos_flag, int dist_flag);


double*** compare(void** problems, int number_of_problems, int* dimensions, int number_of_dimension, int*** optimal_population_size, int number_of_iterations, char* name_of_file);

#endif //BCPROJECT_COMPARISON_H
