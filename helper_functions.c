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
    minimum.value = 999999999999999;
    ret.num_of_evaluations = population->size;
    ret.success = false;

    for (int i = 0; i < population->size; ++i) {
        if (population->individuals[i].fitness < minimum.value) {
            minimum.value = population->individuals[i].fitness;
            minimum.index = i;
        }
    }
    ret.minimum = minimum;
    if (fabs(ret.minimum.value) <= TERMINAL_CONDITION) ret.success = true;
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

    double epsilon = 0.0000000000000001;

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

double** allocMatrix(int size) {
    double** m = (double **) malloc(size * sizeof(double*));
    for (int i = 0; i < size; ++i) {
        m[i] = (double *) malloc(size * sizeof(double));
    }
    return m;
}

double** make1D (double angle) {
    double** matrix = allocMatrix(1);
    matrix[0][0] = cos(angle);
    return matrix;
}

double** make2D (double angle) {
    double** matrix = allocMatrix(2);
    matrix[0][0] = cos(angle);  matrix[0][1] = -sin(angle);
    matrix[1][0] = sin(angle);  matrix[1][1] = cos(angle);
    return matrix;
}

double** make3D (double angle) {
    double** matrix = allocMatrix(3);
    matrix[0][0] = cos(angle);  matrix[0][1] = -sin(angle); matrix[0][2] = 0;
    matrix[1][0] = sin(angle);  matrix[1][1] = cos(angle);  matrix[1][2] = 0;
    matrix[2][0] = 0;           matrix[2][1] = 0;           matrix[2][2] = 1;
    return matrix;
}

double** make4D (double angle) {
    double** matrix = allocMatrix(4);
    matrix[0][0] = cos(angle);  matrix[0][1] = -sin(angle); matrix[0][2] = 0;   matrix[0][3] = 0;
    matrix[1][0] = sin(angle);  matrix[1][1] = cos(angle);  matrix[1][2] = 0;   matrix[1][3] = 0;
    matrix[2][0] = 0;           matrix[2][1] = 0;           matrix[2][2] = 1;   matrix[2][3] = 0;
    matrix[3][0] = 0;           matrix[3][1] = 0;           matrix[3][2] = 0;   matrix[3][3] = 1;
    return matrix;
}

double** make5D () {
    double** matrix = allocMatrix(5);
    matrix[0][0] = 0.954328;    matrix[0][1] = 0.0236714;   matrix[0][2] = 0.0930151;   matrix[0][3] = 0.162359;    matrix[0][4] = 0.231703;
    matrix[1][0] = 0.0845678;   matrix[1][1] = 0.977164;    matrix[1][2] = 0.0388955;   matrix[1][3] = 0.100627;    matrix[1][4] = 0.162359;
    matrix[2][0] = -0.123463;   matrix[2][1] = -0.0693437;  matrix[2][2] = 0.984776;    matrix[2][3] = 0.0388955;   matrix[2][4] = 0.0930151;
    matrix[3][0] = -0.162359;   matrix[3][1] = -0.115851;   matrix[3][2] = -0.0693437;  matrix[3][3] = 0.977164;    matrix[3][4] = 0.0236714;
    matrix[4][0] = -0.201254;   matrix[4][1] = -0.162359;   matrix[4][2] = -0.123463;   matrix[4][3] = -0.0845678;  matrix[4][4] = 0.954328;
    return matrix;
}

void makeMatrices() {
    matrix_1D = make1D(angle);
    matrix_2D = make2D(angle);
    matrix_3D = make3D(angle);
    matrix_4D = make4D(angle);
    matrix_5D = make5D();
}

void freeMatrices() {
    double** a[5] = {matrix_1D, matrix_2D, matrix_3D, matrix_4D, matrix_5D};
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < i+1; ++j) {
            free(a[i][j]);
        }
        free(a[i]);
    }
}


Node* makeSubTree(int from, int to) {
    if (to >= DIMENSION) to = DIMENSION - 1;
    if (from == to) {
        Node* leaf = malloc(sizeof(Node));
        leaf->left = NULL;
        leaf->right = NULL;
        Mask m;
        m.variables = malloc(sizeof(int));
        m.variables[0] = from;
        m.lenght = 1;
        leaf->information = m;
        return leaf;
    }

    Node* inner = malloc(sizeof(Node));

    inner->right = makeSubTree(to, to);
    inner->left = makeSubTree(from, to-1);

    Mask m;
    m.lenght = to - from + 1;
    m.variables = malloc(sizeof (int) * m.lenght);
    int idx = 0;
    for (int i = from; i <= to; ++i) {
        m.variables[idx] = i;
        idx++;
    }
    inner->information = m;
    return inner;

}

void printTree(Node* root) {

    if (root->right != NULL) {
        printTree(root->right);
    }
    if (root->left != NULL) {
        printTree(root->left);
    }

    for (int i = 0; i < root->information.lenght; ++i) {
        printf("%d ", root->information.variables[i]);
    }
    printf("\n");
}

int index_for_mask = 0;

void masksFromTree(Node* root, Mask** mask_array) {

    if (root->right != NULL) {
        mask_array[0][index_for_mask] = root->information;
        index_for_mask++;
        masksFromTree(root->left, mask_array);
        masksFromTree(root->right, mask_array);
    } else {
        mask_array[0][root->information.variables[0] + DIMENSION - 1] = root->information;
    }
}

Mask* makeArrayFromTree() {
    int num_of_leafs = (int) ceil(((double)DIMENSION) / size_of_fraction);

    Node ** nodes = malloc((2 * num_of_leafs - 1) * sizeof(Node*));
    int from = 0;
    int idx = 0;
    for (int i = 0; i < num_of_leafs; ++i) {
        nodes[idx] = makeSubTree(from, from + size_of_fraction - 1);
        idx++;
        from += size_of_fraction;
        if (i != 0) {
            Node *inner = malloc(sizeof(Node));
            Mask m;
            m.lenght = nodes[idx - 1]->information.lenght + nodes[idx - 2]->information.lenght;
            m.variables = malloc(m.lenght * sizeof(int));
            for (int j = 0; j < m.lenght; ++j) {
                m.variables[j] = j;
            }
            inner->right = nodes[idx - 1];
            inner->left = nodes[idx - 2];
            inner->information = m;
            nodes[idx] = inner;
            idx++;
        }
    }
    index_for_mask = 0;
    Mask* mask_array = malloc((2*DIMENSION - 1) * sizeof(Mask));
    masksFromTree(nodes[2*num_of_leafs-2], &mask_array);
    free(nodes);


    Mask* rotated = malloc((2*DIMENSION - 1) * sizeof(Mask));
    for (int i = 0; i < 2*DIMENSION - 1; ++i) {
        rotated[i] = mask_array[2*DIMENSION-2-i];
    }

    return mask_array;


    return mask_array;
}

void freeArrayFromTree(Node* root) {
    if (root->left != NULL) freeArrayFromTree(root->left);
    if (root->right != NULL) freeArrayFromTree(root->right);
    free(root->information.variables);
    free(root);
}

void freeMaskArray(Mask* m) {
    for (int i = 0; i < 2*DIMENSION-1; ++i) {
        free(m[i].variables);
    }
}

