#include <stdio.h>
#include <time.h>
#include <stdlib.h>
//#include "parameters.h"
#include "helper_functions.h"
#include "fitness_functions.h"
//#include "uniform_search.h"
//#include "inform_search.h"
//#include "search.h"
#include <math.h>
#include "comparison.h"



/*
FILE* f;
bool LOG = true;

bool perform25LT(int size) {
    printf("\n");
    printf("Trying %d: ", size);
    if (LOG) fprintf(f,"\nTrying %d: ", size);

    Mask* traversal_array = makeArrayFromTree();

    int failed_count = 0;

    double mean = 0;

    for (int i = 0; i < 25; ++i) {
        Population population = initPopulation(size);
        Statistics stats = DEAlgorithm(&population, 2);
        if (stats.success) {
            if (LOG) fprintf(f, "*");
            printf("*");
            freePopulation(population);
            mean += stats.num_of_evaluations;

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
            if (failed_count >= 2) {
                freeMaskArray(traversal_array);
                free(traversal_array);
                return false;
            }
        }

    }
    mean /= 25-failed_count;
    freeMaskArray(traversal_array);
    free(traversal_array);

    if (LOG) fprintf(f, "\nmean: %5f\n", mean);
    printf("\nmean: %5f\n", mean);

    return true;
}

int findOptimalSizeLT() {
    if (LOG) fprintf(f, "Finding optimal size of population for dimension: %d - LT\n", DIMENSION);
    printf("Finding optimal size of population for dimension: %d - LT\n", DIMENSION);
    int current = 10;

    while (perform25LT(current) == false) {
        current += 2;
        if (current > 100) return -1;
    }

    int lower = current - 2;
    int high = current;
    current = lower + (high-lower) / 2;

    while (current != lower) {
        bool res = perform25LT(current);
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
    FILE* a = fopen("minimum.txt", "a");
    fprintf(a, "LT: Optimal size for dimension %d and size of fraction %d is %d\n", DIMENSION,size_of_fraction, high);
    fclose(a);
    return high;
}

bool perform25MP(int size) {
    printf("\n");
    printf("Trying %d: ", size);
    if (LOG) fprintf(f,"\nTrying %d: ", size);

    int* dep_vector = (int*) malloc(DIMENSION * sizeof(int));
    for (int i = 0; i < DIMENSION; ++i) {
        dep_vector[i] = i / size_of_fraction;
    }

    Dependency dependency = makeDependencies(dep_vector);

    int failed_count = 0;

    double mean = 0;

    for (int i = 0; i < 5; ++i) {
        Population population = initPopulation(size);
        Statistics stats = DEAlgorithm(&population, 1);
        if (stats.success) {
            if (LOG) fprintf(f, "*");
            printf("*");
            freePopulation(population);
            mean += stats.num_of_evaluations;

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
            if (failed_count >= 1) {
                free(dependency.size_of_classes);
                free(dependency.sorted_indexes);
                return false;
            }
        }

    }
    mean /= 5-failed_count;
    free(dependency.size_of_classes);
    free(dependency.sorted_indexes);
    if (LOG) fprintf(f, "\nmean: %5f\n", mean);
    printf("\nmean: %5f\n", mean);

    return true;
}

int findOptimalSizeMP() {
    if (LOG) fprintf(f, "Finding optimal size of population for dimension: %d - MP\n", DIMENSION);
    printf("Finding optimal size of population for dimension: %d - MP\n", DIMENSION);
    int current = 8;

    while (perform25MP(current) == false) {
        current += 2;
        if (current > 30) return -1;
    }

    int lower = current - 2;
    int high = current;
    current = lower + (high-lower) / 2;

    while (current != lower) {
        bool res = perform25MP(current);
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
    FILE* a = fopen("minimum.txt", "a");
    fprintf(a, "MP: Optimal size for dimension %d and size of fraction %d is %d\n", DIMENSION,size_of_fraction, high);
    fclose(a);
    return high;
}

bool perform25Uniform(int size) {
    printf("\n");
    printf("Trying %d: ", size);
    if (LOG) fprintf(f,"\nTrying %d: ", size);

    int failed_count = 0;

    double mean = 0;

    for (int i = 0; i < 25; ++i) {
        Population population = initPopulation(size);
        Statistics stats = DEAlgorithm(&population, 0);
        if (stats.success) {
            printf("*");
            freePopulation(population);
            mean += stats.num_of_evaluations;

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
            if (failed_count >= 2) return false;
        }

    }
    mean /= 25-failed_count;
    if (LOG) fprintf(f, "\nmean: %5f\n", mean);
    printf("\nmean: %5f\n", mean);

    return true;
}

int findOptimalSizeUniformed() {
    if (LOG) fprintf(f, "Finding optimal size of population for dimension: %d - uniform\n", DIMENSION);
    printf("Finding optimal size of population for dimension: %d - uniform\n", DIMENSION);
    int current = 8;

    while (perform25Uniform(current) == false) {
        current += 2;
        if (current > 40) return -1;
    }

    int lower = current - 2;
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
    FILE* a = fopen("minimum.txt", "a");
    fprintf(a, "Uniformaed: Optimal size for dimension %d and size of fraction %d is %d\n", DIMENSION,size_of_fraction, high);
    fclose(a);
    printf("Optimal size for dimension %d is %d\n", DIMENSION, high);
    return high;
}

double findMeanLT(int iterations, int size_of_population) {
    printf("Trying dim: %d, size of fraction: %d, size of population: %d\n", DIMENSION, size_of_fraction, size_of_population);
    double mean = 0;

    Mask* traversal_array = makeArrayFromTree();

    for (int i = 0; i < iterations; ++i) {
        Population p = initPopulation(size_of_population);
        Statistics stats = DEAlgorithm(&p, 2);
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

double findMeanMP(int iterations, int size_of_population) {
    printf("Trying dim: %d, size of fraction: %d, size of population: %d\n", DIMENSION, size_of_fraction, size_of_population);
    double mean = 0;

    int* dep_vector = (int*) malloc(DIMENSION * sizeof(int));
    for (int i = 0; i < DIMENSION; ++i) {
        dep_vector[i] = i / size_of_fraction;
    }

    Dependency dependency = makeDependencies(dep_vector);

    for (int i = 0; i < iterations; ++i) {
        Population p = initPopulation(size_of_population);
        Statistics stats = DEAlgorithm(&p, 1);
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

double findMeanUniform(int iterations, int size_of_population) {

    double mean = 0;
    printf("\n");
    for (int i = 0; i < iterations; ++i) {
        Population p = initPopulation(size_of_population);
        Statistics stats = DEAlgorithm(&p, 0);
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
*/


int main() {

    srand(time(NULL));
    double angle = acos(-1)/2.5;
    makeMatrices(angle);
    F = 0.7;
    size_of_fraction = 5;


    void** problems = malloc(5 * sizeof(void*));
    problems[0] = (void*)sphere;
    problems[2] = (void*)Michalewicz;
    problems[1] = (void*)oneMax;
    problems[3] = (void*)Rastrigin;
    problems[4] = (void*)Rosenbrock;

    int* dimension = malloc(5 * sizeof(int));
    dimension[0] = 5;
    dimension[1] = 10;
    dimension[2] = 15;
    dimension[3] = 20;
    dimension[4] = 25;


    int*** sizes = findOptimalPopulationSize(problems, 5, dimension, 5, 5, "testik.txt");

    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < 5; ++j) {
            for (int k = 0; k < 5; ++k) {
                printf("%d ", sizes[i][j][k]);
            }
            printf("\n");
        }
        printf("----------------\n");
    }

    freeMatrices();
    return 0;
}
