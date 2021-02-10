//
// Created by vojta on 11.11.20.
//

#include <stdlib.h>
#include "search.h"
#include "helper_functions.h"
#include <math.h>

Chromosome uniformMutate(Chromosome* x, Chromosome* a, Chromosome* b, Chromosome* c) {

    Chromosome new;
    new.vector = (double *) malloc(DIMENSION * sizeof(double));

    double* r_is = randomSamples(DIMENSION);
    int R = rand() % DIMENSION;

    for (int i = 0; i < DIMENSION; ++i) {

        if (r_is[i] < CR || i == R) {
            new.vector[i] = x->vector[i] + (b->vector[i] - c->vector[i]) * F;
        } else {
            new.vector[i] = a->vector[i];
        }
    }

    free(r_is);

    //if (fabs(new.vector[0]) > 9999999) printChromosome(new);

    return new;
}

Chromosome mpMutate(Chromosome* x, Chromosome* a, Chromosome* b, Chromosome* c, Dependency* dependency) {

    Chromosome new;
    new.vector = (double *) malloc(DIMENSION * sizeof(double));

    double* r_is = randomSamples(dependency->num_of_classes);
    int R = rand() % DIMENSION;//dependency->num_of_classes;

    int current_index = 0;

    for (int i = 0; i < dependency->num_of_classes; ++i) {

        if (r_is[i] < CR || i == R) {

            for (int j = 0; j < dependency->size_of_classes[i]; ++j) {
                int index = dependency->sorted_indexes[current_index++];
                new.vector[index] = x->vector[index] + (b->vector[index] - c->vector[index]) * F;
            }

        } else {

            for (int j = 0; j < dependency->size_of_classes[i]; ++j) {
                int index = dependency->sorted_indexes[current_index++];
                new.vector[index] = a->vector[index];
            }

        }
    }


    free(r_is);

    //if (fabs(new.vector[0]) > 9999999) printChromosome(new);

    return new;
}

void ltMutate(int index, Chromosome* x, Population* population, Mask* traversal_array) {

    Chromosome new;
    new.vector = (double *) malloc(DIMENSION * sizeof(double));

    for (int i = 0; i < DIMENSION * 2 -1; ++i) {

        Mask m = traversal_array[i];

        int* triplet = pickThree(population->size, i);
        Chromosome* a = &population->individuals[triplet[0]];
        Chromosome* b = &population->individuals[triplet[1]];
        Chromosome* c = &population->individuals[triplet[2]];

        int from = m.variables[0];
        int to = m.variables[m.lenght-1];
        for (int j = 0; j < from; ++j) {
            new.vector[j] = x->vector[j] + (b->vector[j] - c->vector[j]) * F;
        }
        for (int j = from; j <= to; ++j) {
            new.vector[j] = a->vector[j];
        }
        for (int j = to+1; j < DIMENSION; ++j) {
            new.vector[j] = x->vector[j] + (b->vector[j] - c->vector[j]) * F;
        }

        new.fitness = func(&new);

        if (x->fitness > new.fitness) {
            for (int j = 0; j < DIMENSION; ++j) {
                x->vector[j] = new.vector[j];

            }
            x->fitness = new.fitness;
        }
        free(triplet);
    }
    free(new.vector);

}

Minim step(Population* population, Dependency* dependency, int flag, Mask* traversal_array) {

    Minim minimum;
    minimum.value = 9999999999999999;
    Change* changes = (Change*) malloc(population->size * sizeof(Change));
    changes->size = 0;

    for (int i = 0; i < population->size; ++i) {
        Chromosome individual = population->individuals[i];
        int* triplet = pickThree(population->size, i);
        Chromosome new_individual;
        if (flag == 0) {
            new_individual = uniformMutate(&individual, &population->individuals[triplet[0]],
                                                     &population->individuals[triplet[1]],
                                                     &population->individuals[triplet[2]]);
        } else if (flag == 1) {
            new_individual = mpMutate(&individual, &population->individuals[triplet[0]],
                                                     &population->individuals[triplet[1]],
                                                     &population->individuals[triplet[2]], dependency);
        } else if (flag == 2) {
            ltMutate(i, &individual, population, traversal_array);
            if (individual.fitness <= minimum.value) {
                minimum.value = individual.fitness;
                minimum.index = i;
            }
            free(triplet);
            continue;

        } else {
            printf("Wrong flag\n");
        }

        free(triplet);
        double new_value = func(&new_individual);
        new_individual.fitness = new_value;

        if (new_value <= individual.fitness) {
            changes[changes->size].index = i;
            changes[changes->size++].chromosome = new_individual;

            if (new_value <= minimum.value) {
                minimum.value = new_value;
                minimum.index = i;
            }
        } else {
            free(new_individual.vector);
        }


    }

    makeChanges(changes, population);
    free(changes);

    return minimum;
}

Statistics DEAlgorithm(Population* population, Dependency* dependency, int flag, Mask* traversal_array) {
    //printPopulation(population);
    Statistics stats = initStatistics(population);
    int max_num_of_evals = 100000 * DIMENSION;

    while (fabs(stats.minimum.value) >= TERMINAL_CONDITION) {
        //printPopulation(population);
        //sleep(1);
        Minim m = step(population, dependency, flag, traversal_array);

        if (m.value < stats.minimum.value) {
            //printf("New minimum: %5f\n", m.value);
            stats.minimum.value = m.value;
            stats.minimum.index = m.index;
        }
        if (flag == 2) stats.num_of_evaluations += population->size * (2*DIMENSION - 1);
        else stats.num_of_evaluations += population->size;

        if (stats.num_of_evaluations >= max_num_of_evals) {
            //printf("Too much\n");
            if (stats.minimum.value <= TERMINAL_CONDITION) {
                stats.success = true;
            } else stats.success = false;

            return stats;
        }

        if (equalPopulation(population)) {
            stats.success = false;
            //printPopulation(population);
            return stats;
        }

        if (stats.num_of_evaluations % (population->size * 1000) == 0) {
            //printPopulation(population);
        }

        //printChromosome(population->individuals[stats.minimum.index]);
        //printPopulation(population);

    }
    stats.success = true;
    return stats;
}