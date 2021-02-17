//
// Created by vojta on 11.11.20.
//

#include <stdlib.h>
#include "search.h"
#include "helper_functions.h"
#include <math.h>

#define LT_FOS 1
#define MP_FOS 2

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

Chromosome fosMutate(Chromosome* x, Chromosome* a, Chromosome* b, Population* p, FOS fos) {

    Chromosome new;
    new.vector = (double *) malloc(DIMENSION * sizeof(double));

    double* r_is = randomSamples(fos.number_of_subsets);
    int R = rand() % DIMENSION;

    for (int i = 0; i < fos.number_of_subsets; ++i) {

        if (r_is[i] < CR || i == R) {

            for (int j = 0; j < fos.parts[i].number_of_variables; ++j) {
                int index = fos.parts[i].indexes[j];
                new.vector[index] = x->vector[index] + (a->vector[index] - b->vector[index]) * F;
            }

        } else {

            int donor_idx = rand() % p->size;

            for (int j = 0; j < fos.parts[i].number_of_variables; ++j) {
                int index = fos.parts[i].indexes[j];
                new.vector[index] = p->individuals[donor_idx].vector[index];
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


Minim step(Population* population, int flag, FOS* fos) {

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
        } else {
            new_individual = fosMutate(&individual, &population->individuals[triplet[0]],
                                                     &population->individuals[triplet[1]],
                                                     population, *fos);
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


FOS learnFos(Population* population, double** D, int flag) {

    FOS fos;
    if (flag == LT_FOS) {
        fos.number_of_subsets = 2 * DIMENSION - 1;
        fos.parts = malloc(fos.number_of_subsets * sizeof(subset));
        bool* isNodeActive = calloc(2*DIMENSION-1, sizeof(bool));
        for (int i = 0; i < DIMENSION; ++i) {
            isNodeActive[i] = true;
            subset variable;
            variable.number_of_variables = 1;
            variable.indexes = malloc(sizeof(int));
            variable.indexes[0] = i;
            fos.parts[i] = variable;
        }

        for (int i = 0; i < DIMENSION - 1; ++i) {
            double maxMI = -1;
            int nodes[2];

            for (int j = 0; j < DIMENSION+i; ++j) {
                if (!isNodeActive[j]) continue;
                for (int k = j+1; k < DIMENSION+i; ++k) {
                    if (!isNodeActive[k]) continue;

                    double dependency = 0;
                    for (int l = 0; l < fos.parts[j].number_of_variables; ++l) {
                        for (int m = 0; m < fos.parts[k].number_of_variables; ++m) {
                            dependency += D[fos.parts[j].indexes[l]][fos.parts[k].indexes[m]];
                        }
                    }
                    dependency /= fos.parts[j].number_of_variables * fos.parts[k].number_of_variables;

                    if (dependency > maxMI) {
                        maxMI = dependency;
                        nodes[0] = j;
                        nodes[1] = k;
                    }
                }
            }

            fos.parts[DIMENSION + i].number_of_variables = fos.parts[nodes[0]].number_of_variables + fos.parts[nodes[1]].number_of_variables;
            fos.parts[DIMENSION + i].indexes = malloc(fos.parts[DIMENSION + 1].number_of_variables * sizeof(int));

            for (int j = 0; j < fos.parts[DIMENSION + i].number_of_variables; ++j) {

                if (j < fos.parts[nodes[0]].number_of_variables) fos.parts[DIMENSION + i].indexes[j] = fos.parts[nodes[0]].indexes[j];
                else fos.parts[DIMENSION + i].indexes[j] = fos.parts[nodes[1]].indexes[j - fos.parts[nodes[0]].number_of_variables];

            }

            isNodeActive[nodes[0]] = false;
            isNodeActive[nodes[1]] = false;
            isNodeActive[DIMENSION+i] = true;
        }

    } else if (flag == MP_FOS) {

        fos.number_of_subsets = 2 * DIMENSION - 1;
        fos.parts = malloc(fos.number_of_subsets * sizeof(subset));
        bool* isNodeActive = calloc(2*DIMENSION-1, sizeof(bool));
        for (int i = 0; i < DIMENSION; ++i) {
            isNodeActive[i] = true;
            subset variable;
            variable.number_of_variables = 1;
            variable.indexes = malloc(sizeof(int));
            variable.indexes[0] = i;
            fos.parts[i] = variable;
        }

        for (int i = 0; i < DIMENSION - 1; ++i) {
            double maxMI = -1;
            int nodes[2];

            for (int j = 0; j < DIMENSION+i; ++j) {
                if (!isNodeActive[j]) continue;
                for (int k = j+1; k < DIMENSION+i; ++k) {
                    if (!isNodeActive[k]) continue;

                    double dependency = 0;
                    for (int l = 0; l < fos.parts[j].number_of_variables; ++l) {
                        for (int m = 0; m < fos.parts[k].number_of_variables; ++m) {
                            dependency += D[fos.parts[j].indexes[l]][fos.parts[k].indexes[m]];
                        }
                    }
                    dependency /= fos.parts[j].number_of_variables * fos.parts[k].number_of_variables;

                    if (dependency > maxMI) {
                        maxMI = dependency;
                        nodes[0] = j;
                        nodes[1] = k;
                    }
                }
            }

            if(fabs(maxMI) < 0.00001) break;

            fos.parts[DIMENSION + i].number_of_variables = fos.parts[nodes[0]].number_of_variables + fos.parts[nodes[1]].number_of_variables;
            fos.parts[DIMENSION + i].indexes = malloc(fos.parts[DIMENSION + 1].number_of_variables * sizeof(int));

            for (int j = 0; j < fos.parts[DIMENSION + i].number_of_variables; ++j) {

                if (j < fos.parts[nodes[0]].number_of_variables) fos.parts[DIMENSION + i].indexes[j] = fos.parts[nodes[0]].indexes[j];
                else fos.parts[DIMENSION + i].indexes[j] = fos.parts[nodes[1]].indexes[j - fos.parts[nodes[0]].number_of_variables];

            }

            isNodeActive[nodes[0]] = false;
            isNodeActive[nodes[1]] = false;
            isNodeActive[DIMENSION+i] = true;
        }

        int active_count = 0;
        for (int i = 0; i < 2*DIMENSION-1; ++i) {
            if (isNodeActive[i]) active_count++;
        }
        FOS reduced_fos;
        reduced_fos.number_of_subsets = active_count;
        reduced_fos.parts = malloc(reduced_fos.number_of_subsets * sizeof(subset));
        int idx = 0;
        for (int i = 0; i < 2*DIMENSION - 1; ++i) {
            if (isNodeActive[i]) {
                reduced_fos.parts[idx].number_of_variables = fos.parts[i].number_of_variables;
                reduced_fos.parts[idx].indexes = malloc(reduced_fos.parts[idx].number_of_variables * sizeof(int));
                for (int j = 0; j < reduced_fos.parts[idx].number_of_variables; ++j) {
                    reduced_fos.parts[idx].indexes[j] = fos.parts[i].indexes[j];
                }
                idx++;
            }
        }
        printFOS(reduced_fos);
        return reduced_fos;
    }

    return fos;
}

Statistics DEAlgorithm(Population* population, int flag) {
    //printPopulation(population);
    Statistics stats = initStatistics(population);
    int max_num_of_evals = 100000 * DIMENSION;



    while (fabs(stats.minimum.value) >= TERMINAL_CONDITION) {
        //printPopulation(population);
        //sleep(1);

        //double** D = createDependencyMatrix(population->individuals[rand() % population->size]);

        double** D = createDependencyMatrix(population);
        FOS fos = learnFos(population, D, LT_FOS);
        printFOS(fos);
        Minim m = step(population, flag, &fos);
        freeFOS(fos, D);

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