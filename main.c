#include <stdio.h>
#include <time.h>
#include <stdlib.h>
//#include "parameters.h"
#include "helper_functions.h"
#include "fitness_functions.h"
//#include "uniform_search.h"
//#include "inform_search.h"
#include "search.h"

FILE* f;
bool LOG = true;

bool perform25Inform(int size) {
    printf("\n");
    printf("Trying %d: ", size);
    if (LOG) fprintf(f,"\nTrying %d: ", size);

    int* dep_vector = (int*) malloc(DIMENSION * sizeof(int));
    for (int i = 0; i < DIMENSION; ++i) {
        dep_vector[i] = i / size_of_fraction;
    }

    Dependency dependency = makeDependencies(dep_vector);

    int failed_count = 0;

    for (int i = 0; i < 25; ++i) {
        Population population = initPopulation(size);
        Statistics stats = DEAlgorithm(&population, &dependency, 1);
        if (stats.success) {
            if (LOG) fprintf(f, "*");
            printf("*");
            freePopulation(population);

        } else {
            //printStatistics(&stats, &population);
            if (equalPopulation(&population)) {
                if (LOG) fprintf(f, "@(%5f)", stats.minimum.value);
                printf("@(%5f)", stats.minimum.value);
            }
            else {
                if (LOG) fprintf(f, "#(%5f)", stats.minimum.value);
                printf("#(%5f)", stats.minimum.value);
            }
            freePopulation(population);
            failed_count++;
            if (failed_count >= 3) {
                free(dependency.size_of_classes);
                free(dependency.sorted_indexes);
                return false;
            }
        }

    }

    free(dependency.size_of_classes);
    free(dependency.sorted_indexes);

    return true;
}

int findOptimalSizeInformed() {
    if (LOG) fprintf(f, "Finding optimal size of population for dimension: %d - inform\n", DIMENSION);
    printf("Finding optimal size of population for dimension: %d - inform\n", DIMENSION);
    int current = 4;

    while (perform25Inform(current) == false) {
        current += 10;
    }

    int lower = current - 10;
    int high = current;
    current = lower + (high-lower) / 2;

    while (current != lower) {
        bool res = perform25Inform(current);
        if (res) {
            high = current;
            current = lower + (high-lower) / 2;
        } else {
            lower = current;
            current = lower + (high-lower) / 2;
        }
    }
    if (LOG) fprintf(f, "Optimal size for dimension %d is %d\n", DIMENSION, high);
    printf("Optimal size for dimension %d is %d\n", DIMENSION, high);
    return high;
}

bool perform25Uniform(int size) {
    printf("\n");
    printf("Trying %d: ", size);
    if (LOG) fprintf(f,"\nTrying %d: ", size);

    int failed_count = 0;

    for (int i = 0; i < 25; ++i) {
        Population population = initPopulation(size);
        Statistics stats = DEAlgorithm(&population, NULL, 0);
        if (stats.success) {
            printf("*");
            freePopulation(population);

        } else {
            if (equalPopulation(&population)) {
                if (LOG) fprintf(f,"@(%5f)", stats.minimum.value);
                printf("@(%5f)", stats.minimum.value);
            }
            else {
                if (LOG) fprintf(f, "#(%5f)", stats.minimum.value);
                printf("#(%5f)", stats.minimum.value);
            }
            freePopulation(population);
            failed_count++;
            if (failed_count >= 3) return false;
        }

    }

    return true;
}

int findOptimalSizeUniformed() {
    if (LOG) fprintf(f, "Finding optimal size of population for dimension: %d - uniform\n", DIMENSION);
    printf("Finding optimal size of population for dimension: %d - uniform\n", DIMENSION);
    int current = 4;

    while (perform25Uniform(current) == false) {
        current += 10;
    }

    int lower = current - 10;
    int high = current;
    current = lower + (high-lower) / 2;

    while (current != lower) {
        bool res = perform25Uniform(current);
        if (res) {
            high = current;
            current = lower + (high-lower) / 2;
        } else {
            lower = current;
            current = lower + (high-lower) / 2;
        }
    }
    if (LOG) fprintf(f, "Optimal size for dimension %d is %d\n", DIMENSION, high);
    printf("Optimal size for dimension %d is %d\n", DIMENSION, high);
    return high;
}


double findMeanInform(int iterations, int size_of_population) {

    double mean = 0;

    int* dep_vector = (int*) malloc(DIMENSION * sizeof(int));
    for (int i = 0; i < DIMENSION; ++i) {
        dep_vector[i] = i / size_of_fraction;
    }

    Dependency dependency = makeDependencies(dep_vector);

    for (int i = 0; i < iterations; ++i) {
        Population p = initPopulation(size_of_population);
        Statistics stats = DEAlgorithm(&p, &dependency, 1);
        freePopulation(p);
        if (!stats.success) {
            printf("#");
            i -= 1;
        } else {
            printf("*");
            mean += stats.num_of_evaluations;
        }
    }

    mean /= iterations;

    return mean;
}

double findMeanUniform(int iterations, int size_of_population) {

    double mean = 0;
    printf("\n");
    for (int i = 0; i < iterations; ++i) {
        Population p = initPopulation(size_of_population);
        Statistics stats = DEAlgorithm(&p, NULL, 0);
        freePopulation(p);
        if (!stats.success) {
            printf("#");
            i -= 1;
        } else {
            printf("*");
            mean += stats.num_of_evaluations;
        }
    }

    mean /= iterations;

    return mean;
}

int main() {

    srand(time(NULL));
    DIMENSION = 20;
    func = myBenchMarkFunction;
    size_of_fraction = 2;

    //LOG = false;

    int dimensions[10] = {10, 13, 17, 22, 29, 36, 46, 60, 77, 100};
    int inform_dimensions[10] = {24, 30, 33, 47, 34, 41, 41, 56};
    int uniform_dimensions[10] = {7, 10, 8, 6, 10, 5, 9, 5};
    f = fopen("Results.txt", "w");


    for (int i = 0; i < 1; ++i) {
        DIMENSION = dimensions[i];
        //findOptimalSizeInformed();
        findOptimalSizeUniformed();
    }


    /*
    for (int i = 6; i < 8; ++i) {
        if (LOG) fprintf(f, "Dimension: %d\n", dimensions[i]);
        printf("Dimension: %d\n", dimensions[i]);
        DIMENSION = dimensions[i];
        size_of_fraction = DIMENSION / 4;
        //size_of_fraction = 1;
        double I_mean = findMeanInform(25, inform_dimensions[i]);
        printf("Mean inform: %5f\n", I_mean);
        double U_mean = findMeanUniform(25, uniform_dimensions[i]);
        printf("Mean uniform: %5f\n", U_mean);
    }
    */

    fclose(f);




    /*
    int* dep_vector = (int*) malloc(DIMENSION * sizeof(int));
    for (int i = 0; i < DIMENSION; ++i) {
        dep_vector[i] = i / size_of_fraction;
    }

    Dependency dependency = makeDependencies(dep_vector);

    informAlgorithm(&p, &dependency);
    */

    return 0;
}
