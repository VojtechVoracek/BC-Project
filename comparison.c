//
// Created by voracvoj on 10.03.21.
//

#include <stdio.h>
#include <stdlib.h>

#include "comparison.h"
#include "search.h"
#include "helper_functions.h"



bool performN(int n, int size, FILE* f, int fos_flag, int dist_flag) {
    printf("\n");
    printf("Trying %d: ", size);
    fprintf(f,"\nTrying %d: ", size);

    int failed_count = 0;
    double mean = 0;

    for (int i = 0; i < n; ++i) {
        Population population = initPopulation(size);
        Statistics stats = DEAlgorithm(&population, fos_flag, dist_flag);
        if (stats.success) {
            fprintf(f, "*");
            printf("*");
            freePopulation(population);
            mean += stats.num_of_evaluations;

        } else {
            //printStatistics(&stats, &population);
            if (equalPopulation(&population)) {
                fprintf(f, "@(%5f)", stats.minimum.value);
                printf("@(%5f)", stats.minimum.value);
            }
            else {
                fprintf(f, "#(%5f)", stats.minimum.value);
                printf("#(%5f)", stats.minimum.value);
            }
            freePopulation(population);
            failed_count++;
            if (failed_count > n*0.04 ) {
                return false;
            }
        }

    }
    mean /= n-failed_count;

    fprintf(f, "\nmean: %5f\n", mean);
    printf("\nmean: %5f\n", mean);

    return true;
}


int findOptimalSize(FILE* f, int num_of_iterations, int fos_flag, int dist_flag) {
    fprintf(f, "Finding optimal size of population for dimension: %d\n", DIMENSION);
    printf("Finding optimal size of population for dimension: %d\n", DIMENSION);
    int current = 5;
    if (dist_flag == MIC && fos_flag != UNIFORM) current = 30;

    while (performN(num_of_iterations, current, f, fos_flag, dist_flag) == false) {

        if (dist_flag == MIC && fos_flag != UNIFORM) {
            current += 10;
            if (current > 1000) return -1;
        }

        if (dist_flag == LINC || fos_flag == UNIFORM) {
            current += 2;
            if (current > 30) return -1;
        }
    }

    int lower = current - 2;
    int high = current;
    current = lower + (high-lower) / 2;

    while (current != lower) {
        bool res = performN(num_of_iterations, current, f, fos_flag, dist_flag);
        if (res) {
            high = current;
            current = lower + (high-lower) / 2;
        } else {
            lower = current;
            current = lower + (high-lower) / 2;
        }
    }
    fprintf(f, "Optimal size for dimension %d is %d\n", DIMENSION, high);
    printf("Optimal size for dimension %d is %d\n", DIMENSION, high);
    FILE* a = fopen("minimum.txt", "a");
    fprintf(a, "MP: Optimal size for dimension %d and size of fraction %d is %d\n", DIMENSION,size_of_fraction, high);
    fclose(a);
    return high;
}


int*** findOptimalPopulationSize(void** problems, int number_of_problems, int* dimensions, int number_of_dimension, int number_of_iterations, char* name_of_file) {
    FILE* f = fopen(name_of_file, "w");
    int algorithms[9] = {UNIFORM, LT_FOS_MIC, LT_FOS_LINC, MP_FOS_MIC, MP_FOS_LINC, LT_FOS_f_MIC, LT_FOS_f_LINC, MP_FOS_f_MIC, MP_FOS_f_LINC};
    int*** results = malloc(9 * sizeof(int **));

    for (int a = 0; a < 9; ++a) {
        results[a] = malloc(number_of_problems * sizeof(int *));
        for (int p = 0; p < number_of_problems; ++p) {
            results[a][p] = malloc(number_of_problems * sizeof(int));
            for (int d = 0; d < number_of_dimension; ++d) {
                func = problems[p];
                DIMENSION = dimensions[d];
                int dist_flag = (a+1) % 2 + 1;
                if (a == 0) dist_flag = 0;
                int fos_flag = (a + 1) / 2;
                printf("Zkousim s %d %d, problem cislo %d\n", fos_flag, dist_flag, p);
                results[a][p][d] = findOptimalSize(f, number_of_iterations, fos_flag, dist_flag);
            }
        }
    }

    return results;
}


double findMean(int iterations, int size_of_population, int fos_flag, int dist_flag) {
    printf("Trying dim: %d, size of fraction: %d, size of population: %d\n", DIMENSION, size_of_fraction, size_of_population);
    double mean = 0;

    for (int i = 0; i < iterations; ++i) {
        Population p = initPopulation(size_of_population);
        Statistics stats = DEAlgorithm(&p, fos_flag, dist_flag);
        freePopulation(p);
        if (!stats.success) {
            printf("#");
            i -= 1;
        } else {
            printf("*");
            mean += stats.num_of_evaluations;
        }
    }
    printf("\n");
    mean /= iterations;

    return mean;
}


double*** compare(void** problems, int number_of_problems, int* dimensions, int number_of_dimension, int*** optimal_population_size, int number_of_iterations, char* name_of_file) {

    FILE* f = fopen(name_of_file, "w");
    int algorithms[9] = {UNIFORM, LT_FOS_MIC, LT_FOS_LINC, MP_FOS_MIC, MP_FOS_LINC, LT_FOS_f_MIC, LT_FOS_f_LINC, MP_FOS_f_MIC, MP_FOS_f_LINC};
    double*** results = malloc(9 * sizeof(double **));

    for (int a = 0; a < 9; ++a) {
        results[a] = malloc(number_of_problems * sizeof(double *));
        for (int p = 0; p < number_of_problems; ++p) {
            results[a] = malloc(number_of_problems * sizeof(double));
            for (int d = 0; d < number_of_dimension; ++d) {
                func = problems[p];
                DIMENSION = dimensions[d];
                int dist_flag = a % 2;
                int fos_flag = (a + 1) / 2;
                results[a][p][d] = findMean(number_of_iterations, optimal_population_size[a][p][d], fos_flag, dist_flag);
            }
        }
    }

    return results;
}
