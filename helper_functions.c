//
// Created by vojta on 10.11.20.
//

#include "helper_functions.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "minepy-1.2.5/libmine/mine.h"
#include "fitness_functions.h"

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
        m[i] = (double *) calloc(size, sizeof(double));
    }
    return m;
}

double** makeOriginRotationMatrix(int dimension, double angle) {
    double** rot_matrix = allocMatrix(dimension);
    for (int i = 0; i < dimension / 2; ++i) {
        rot_matrix[2*i][2*i] = cos(angle);  rot_matrix[2*i][2*i+1] = -sin(angle);
        rot_matrix[2*i+1][2*i] = sin(angle);  rot_matrix[2*i+1][2*i+1] = cos(angle);
    }
    if (dimension % 2 == 1) rot_matrix[dimension-1][dimension-1] = cos(angle);
    return rot_matrix;
}

double** makeTranslationMatrix(int dimension, double* point) {
    double** transition_matrix = allocMatrix(dimension);
    for (int i = 0; i < dimension; ++i) {
        transition_matrix[i][i] = 1;
    }
    for (int i = 0; i < dimension - 1; ++i) {
        transition_matrix[i][dimension - 1] = -point[i];
    }


    return transition_matrix;
}

double** makeInverseTranslationMatrix(int dimension, double* point) {
    double** inv_transition_matrix = allocMatrix(dimension);
    for (int i = 0; i < dimension; ++i) {
        inv_transition_matrix[i][i] = 1;
    }
    for (int i = 0; i < dimension - 1; ++i) {
        inv_transition_matrix[i][dimension - 1] = point[i];
    }

    return inv_transition_matrix;
}

double** squareMatrixMultiplication(double** M1, double** M2, int dimension) {
    double** result = allocMatrix(dimension);
    for (int r = 0; r < dimension; ++r) {
        for (int c = 0; c < dimension; ++c) {
            double tmp = 0;
            for (int i = 0; i < dimension; ++i) {
                tmp += M1[r][i] * M2[i][c];
            }
            result[r][c] = tmp;
        }
    }

    return result;

}

double** makeRotationMatrix(int dimension, double angle, double* point) {
    double** rotation = makeOriginRotationMatrix(dimension, angle);
    double** translation = makeTranslationMatrix(dimension, point);
    double** inv_translation = makeInverseTranslationMatrix(dimension, point);

    double** tmp = squareMatrixMultiplication(rotation, translation, dimension);
    freeSquareMatrix(rotation, dimension);
    freeSquareMatrix(translation, dimension);

    double** result = squareMatrixMultiplication(inv_translation, tmp, dimension);
    freeSquareMatrix(inv_translation, dimension);
    freeSquareMatrix(tmp, dimension);

    return result;
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
    /*
    double** matrix = allocMatrix(3);
    matrix[0][0] = pow(cos(angle), 3) - pow(sin(angle), 2);                    matrix[0][1] = -sin(angle) * pow(cos(angle), 2) - sin(angle)*cos(angle);    matrix[0][2] = sin(angle) * cos(angle);
    matrix[1][0] = sin(angle) * cos(angle) + sin(angle) * pow(cos(angle), 2);     matrix[1][1] = pow(cos(angle), 2) - pow(sin(angle), 2)*cos(angle);       matrix[1][2] = pow(sin(angle), 2);
    matrix[2][0] = -sin(angle) * cos(angle);                                         matrix[2][1] = pow(sin(angle), 2);                                          matrix[2][2] = cos(angle);

    double** rot_matrix = allocMatrix(3);
    rot_matrix[0][0] = cos(angle);  rot_matrix[0][1] = -sin(angle);  rot_matrix[0][2] = 0;
    rot_matrix[1][0] = sin(angle);  rot_matrix[1][1] = cos(angle);   rot_matrix[1][2] = 0;
    rot_matrix[2][0] = 0;           rot_matrix[2][1] = 0;            rot_matrix[2][2] = 1;

    double** new_base_matrix = allocMatrix(3);
    new_base_matrix[0][0] = 1;  new_base_matrix[0][1] = 2;   new_base_matrix[0][2] = 3;
    new_base_matrix[1][0] = 0;  new_base_matrix[1][1] = 1;   new_base_matrix[1][2] = 0;
    new_base_matrix[2][0] = 1;  new_base_matrix[2][1] = 3;   new_base_matrix[2][2] = 1;

    double** inv_new_base_matrix = allocMatrix(3);
    inv_new_base_matrix[0][0] = 1;  inv_new_base_matrix[0][1] = 0;   inv_new_base_matrix[0][2] = 0;
    inv_new_base_matrix[1][0] = 0;  inv_new_base_matrix[1][1] = 1;   inv_new_base_matrix[1][2] = 0;
    inv_new_base_matrix[2][0] = -1;  inv_new_base_matrix[2][1] = 0;   inv_new_base_matrix[2][2] = 1;
    */

    double** rot_matrix_1 = allocMatrix(3);
    rot_matrix_1[0][0] = cos(angle);  rot_matrix_1[0][1] = -sin(angle);   rot_matrix_1[0][2] = 0;
    rot_matrix_1[1][0] = sin(angle);  rot_matrix_1[1][1] = cos(angle);    rot_matrix_1[1][2] = 0;
    rot_matrix_1[2][0] = 0;           rot_matrix_1[2][1] = 0;             rot_matrix_1[2][2] = 1;

    double** rot_matrix_2 = allocMatrix(3);
    rot_matrix_2[0][0] = cos(angle);    rot_matrix_2[0][1] = 0;   rot_matrix_2[0][2] = sin(angle);
    rot_matrix_2[1][0] = 0;             rot_matrix_2[1][1] = 1;   rot_matrix_2[1][2] = 0;
    rot_matrix_2[2][0] = -sin(angle);   rot_matrix_2[2][1] = 0;   rot_matrix_2[2][2] = cos(angle);

    double** matrix = allocMatrix(3);
    //matrix[0][0] = 1;    matrix[0][1] = 0;            matrix[0][2] = 0;
    //matrix[1][0] = 0;    matrix[1][1] = cos(angle);   matrix[1][2] = -sin(angle);
    //matrix[2][0] = 0;    matrix[2][1] = sin(angle);   matrix[2][2] = cos(angle);

    matrix[0][0] = cos(angle) * sin(angle); matrix[0][1] = -pow(sin(angle), 2); matrix[0][2] = pow(cos(angle), 2);
    matrix[1][0] = pow(sin(angle), 2); matrix[1][1] = cos(angle) * sin(angle); matrix[1][2] = cos(angle) * sin(angle);
    matrix[2][0] = cos(angle); matrix[2][1] = 0; matrix[2][2] = -sin(angle);

    return matrix;
    //return squareMatrixMultiplication(squareMatrixMultiplication(rot_matrix_2, rot_matrix_1, 3), matrix, 3);
    //return squareMatrixMultiplication(rot_matrix_1, rot_matrix_2, 3);
}

double** make4D (double angle) {
    /*
    double** matrix1 = allocMatrix(4);
    matrix1[0][0] = cos(angle); matrix1[0][1] = -sin(angle); matrix1[0][2] = 0; matrix1[0][3] = 0;
    matrix1[1][0] = sin(angle); matrix1[1][1] = cos(angle); matrix1[1][2] = 0; matrix1[1][3] = 0;
    matrix1[2][0] = 0; matrix1[2][1] = 0; matrix1[2][2] = cos(angle); matrix1[2][3] = -sin(angle);
    matrix1[3][0] = 0; matrix1[3][1] = 0; matrix1[3][2] = sin(angle); matrix1[3][3] = cos(angle);

    double** matrix2 = allocMatrix(4);
    matrix2[0][0] = cos(angle);   matrix2[0][1] = 0;            matrix2[0][2] = sin(angle);   matrix2[0][3] = 0;
    matrix2[1][0] = 0;            matrix2[1][1] = cos(angle);   matrix2[1][2] = 0;            matrix2[1][3] = sin(angle);
    matrix2[2][0] = -sin(angle);  matrix2[2][1] = 0;            matrix2[2][2] = cos(angle);   matrix2[2][3] = 0;
    matrix2[3][0] = 0;            matrix2[3][1] = -sin(angle);  matrix2[3][2] = 0;            matrix2[3][3] = cos(angle);

    double** matrix3 = allocMatrix(4);
    matrix3[0][0] = cos(angle);  matrix3[0][1] = 0;           matrix3[0][2] = 0;             matrix3[0][3] = -sin(angle);
    matrix3[1][0] = 0;           matrix3[1][1] = cos(angle);  matrix3[1][2] = -sin(angle);   matrix3[1][3] = 0;
    matrix3[2][0] = 0;           matrix3[2][1] = sin(angle);  matrix3[2][2] = cos(angle);    matrix3[2][3] = 0;
    matrix3[3][0] = sin(angle);  matrix3[3][1] = 0;           matrix3[3][2] = 0;             matrix3[3][3] = cos(angle);

    return squareMatrixMultiplication(squareMatrixMultiplication(matrix3, matrix2, 4), matrix1, 4);

    */
    double** matrix1 = allocMatrix(4);
    matrix1[0][0] = pow(cos(angle), 3);                matrix1[0][1] = -pow(cos(angle), 2) * sin(angle);    matrix1[0][2] = -pow(cos(angle), 2) * sin(angle);    matrix1[0][3] = -pow(cos(angle), 2) * sin(angle);
    matrix1[1][0] = pow(cos(angle), 2) * sin(angle);   matrix1[1][1] = pow(cos(angle), 3);    matrix1[1][2] = -pow(sin(angle), 2) * cos(angle);    matrix1[1][3] = -pow(sin(angle), 2) * cos(angle);
    matrix1[2][0] = cos(angle) * sin(angle);              matrix1[2][1] = 0;                 matrix1[2][2] = pow(cos(angle),2);                   matrix1[2][3] = -pow(sin(angle), 2);
    matrix1[3][0] = sin(angle);                           matrix1[3][1] = 0;                                      matrix1[3][2] = 0;                                      matrix1[3][3] = cos(angle);

    return matrix1;

}

double** make5D (double angle) {
    double** matrix = allocMatrix(5);
    matrix[0][0] = 0.954328;    matrix[0][1] = 0.0236714;   matrix[0][2] = 0.0930151;   matrix[0][3] = 0.162359;    matrix[0][4] = 0.231703;
    matrix[1][0] = 0.0845678;   matrix[1][1] = 0.977164;    matrix[1][2] = 0.0388955;   matrix[1][3] = 0.100627;    matrix[1][4] = 0.162359;
    matrix[2][0] = -0.123463;   matrix[2][1] = -0.0693437;  matrix[2][2] = 0.984776;    matrix[2][3] = 0.0388955;   matrix[2][4] = 0.0930151;
    matrix[3][0] = -0.162359;   matrix[3][1] = -0.115851;   matrix[3][2] = -0.0693437;  matrix[3][3] = 0.977164;    matrix[3][4] = 0.0236714;
    matrix[4][0] = -0.201254;   matrix[4][1] = -0.162359;   matrix[4][2] = -0.123463;   matrix[4][3] = -0.0845678;  matrix[4][4] = 0.954328;

    matrix[0][0] = cos(angle);  matrix[0][1] = -sin(angle);  matrix[0][2] = 0;                         matrix[0][3] = 0;                        matrix[0][4] = 0;
    matrix[1][0] = sin(angle);  matrix[1][1] = cos(angle);   matrix[1][2] = 0;                         matrix[1][3] = 0;                        matrix[1][4] = 0;
    matrix[2][0] = -0;          matrix[2][1] = -0;           matrix[2][2] = cos(angle) * sin(angle);   matrix[2][3] = -pow(sin(angle), 2);   matrix[2][4] = pow(cos(angle), 2);
    matrix[3][0] = -0;          matrix[3][1] = -0;           matrix[3][2] = pow(sin(angle), 2);     matrix[3][3] = cos(angle) * sin(angle);  matrix[3][4] = cos(angle) * sin(angle);
    matrix[4][0] = -0;          matrix[4][1] = -0;           matrix[4][2] = cos(angle);                matrix[4][3] = 0;                       matrix[4][4] = -sin(angle);

    double** matrix2 = allocMatrix(5);
    matrix2[0][0] = cos(angle) * sin(angle);   matrix2[0][1] = 0;             matrix2[0][2] = -pow(sin(angle), 2);    matrix2[0][3] = 0;             matrix2[0][4] = pow(cos(angle), 2);
    matrix2[1][0] = 0;                         matrix2[1][1] = cos(angle);    matrix2[1][2] = 0;                         matrix2[1][3] = -sin(angle);   matrix2[1][4] = 0;
    matrix2[2][0] = pow(sin(angle), 2);     matrix2[2][1] = 0;             matrix2[2][2] = cos(angle) * sin(angle);   matrix2[2][3] = 0;             matrix2[2][4] = cos(angle) * sin(angle);
    matrix2[3][0] = 0;                         matrix2[3][1] = sin(angle);    matrix2[3][2] = 0;                         matrix2[3][3] = cos(angle);    matrix2[3][4] = cos(angle) * sin(angle);;
    matrix2[4][0] = cos(angle);                matrix2[4][1] = 0;             matrix2[4][2] = 0;                         matrix2[4][3] = 0;             matrix2[4][4] = -sin(angle);

    double** matrix3 = allocMatrix(5);
    matrix3[0][0] = cos(angle) * sin(angle);   matrix3[0][1] = 0;             matrix3[0][2] = 0;                         matrix3[0][3] = -pow(sin(angle), 2);   matrix3[0][4] = pow(cos(angle), 2);
    matrix3[1][0] = 0;                         matrix3[1][1] = cos(angle);    matrix3[1][2] = sin(angle);                matrix3[1][3] = 0;                        matrix3[1][4] = 0;
    matrix3[2][0] = 0;                         matrix3[2][1] = -sin(angle);   matrix3[2][2] = cos(angle);                matrix3[2][3] = 0;                        matrix3[2][4] = 0;
    matrix3[3][0] = pow(sin(angle), 2);     matrix3[3][1] = 0;             matrix3[3][2] = 0;                         matrix3[3][3] = cos(angle) * sin(angle);  matrix3[3][4] = cos(angle) * sin(angle);
    matrix3[4][0] = cos(angle);                matrix3[4][1] = 0;             matrix3[4][2] = 0;                         matrix3[4][3] = 0;                        matrix3[4][4] = -sin(angle);



    matrix[0][0] = pow(cos(angle), 4);                matrix[0][1] = -sin(angle) * pow(cos(angle), 3);   matrix[0][2] = -sin(angle) * pow(cos(angle), 3);              matrix[0][3] = -sin(angle) * pow(cos(angle), 3);              matrix[0][4] = -sin(angle) * pow(cos(angle), 3);
    matrix[1][0] = sin(angle) * pow(cos(angle), 3);   matrix[1][1] = pow(cos(angle), 4);                 matrix[1][2] = -pow(sin(angle), 2) * pow(cos(angle), 2);   matrix[1][3] = -pow(sin(angle), 2) * pow(cos(angle), 2);   matrix[1][4] = -pow(sin(angle), 2) * pow(cos(angle), 2);
    matrix[2][0] = sin(angle) * pow(cos(angle), 2);   matrix[2][1] = 0;                                     matrix[2][2] = pow(cos(angle), 3);                            matrix[2][3] = -pow(sin(angle), 2) * cos(angle);              matrix[2][4] = -pow(sin(angle), 2) * cos(angle);
    matrix[3][0] = sin(angle) * pow(cos(angle), 1);   matrix[3][1] = 0;                                     matrix[3][2] = 0;                                                matrix[3][3] = pow(cos(angle), 2);                            matrix[3][4] = -pow(sin(angle), 2);
    matrix[4][0] = sin(angle);                           matrix[4][1] = 0;                                     matrix[4][2] = 0;                                                matrix[4][3] = 0;                                                matrix[4][4] = cos(angle);

    return matrix;
    //return squareMatrixMultiplication(squareMatrixMultiplication(matrix, matrix2, 5), matrix3, 5);
}

void makeMatrices(double angle) {

    double* point = malloc(5 * sizeof(double));

    for (int i = 0; i < 5; ++i) {
        point[i] = 1;
    }

    matrix_1D = makeRotationMatrix(1, angle, point);
    matrix_2D = makeRotationMatrix(2, angle, point);
    matrix_3D = makeRotationMatrix(3, angle, point);
    matrix_4D = makeRotationMatrix(4, angle, point);
    matrix_5D = makeRotationMatrix(5, angle, point);
    matrix_5D = make5D(angle);
    matrix_3D = make3D(angle);
    matrix_4D = make4D(angle);
}

void printSquareMatrix(double** M, int dimension, char* name) {
    printf("%s:\n", name);
    for (int i = 0; i < dimension; ++i) {
        for (int j = 0; j < dimension; ++j) {
            printf("%5f ", M[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void freeSquareMatrix(double** M, int dimension) {
    for (int i = 0; i < dimension; ++i) {
        free(M[i]);
    }
    free(M);
}

void freeMatrices() {
    double** a[5] = {matrix_1D, matrix_2D, matrix_3D, matrix_4D, matrix_5D};
    for (int i = 0; i < 5; ++i) {
        freeSquareMatrix(a[i], i + 1);
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

int* individuals2Learning (Population* population) {
    int* indexes = malloc(population->size * sizeof(int));
    int idx = 0;

    for (int i = 0; i < population->size && idx < population->size; ++i) {
        int rn = population->size - i;
        int rm = population->size - idx;
        if (rand() % rn < rm) indexes[idx++] = i;
    }

    int* selected_indexes = malloc((population->size/2) * sizeof(int));

    for (int i = 0; i < population->size/2; ++i) {
        if (population->individuals[indexes[i]].fitness < population->individuals[indexes[i+1]].fitness) {
            selected_indexes[i] = indexes[i];
        } else {
            selected_indexes[i] = indexes[i+1];
        }
    }

    free(indexes);

    return selected_indexes;
}

void freeFOS(FOS fos, double** D) {
    for (int i = 0; i < fos.number_of_subsets; ++i) {
        free(fos.parts[i].indexes);
    }
    free(fos.parts);
    for (int i = 0; i < DIMENSION; ++i) {
        free(D[i]);
    }
    free(D);
}

double** createDependencyMatrix(Population* population) {

    //int* selected_individuals = individuals2Learning(population);
    int n = floor(population->size * 0.5);
    int* selected_individuals = findNBest(population, n);
    double* a = calloc(DIMENSION, sizeof(double));
    double* b = calloc(DIMENSION, sizeof(double));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < DIMENSION; ++j) {
            double value = population->individuals[selected_individuals[i]].vector[j];
            if (value > a[j]) a[j] = value;
            if (value < b[j]) b[j] = value;
        }
    }

    for (int i = 0; i < DIMENSION; ++i) {
        double tmp = b[i];
        b[i] = 0.35 * (a[i] - b[i]);
        a[i] = tmp + b[i];
    }

    double** D = malloc(DIMENSION * sizeof(double*));

    for (int i = 0; i < DIMENSION; ++i) {
        D[i] = malloc(DIMENSION * sizeof(double));

    }
    int z = rand() % n;
    //printf("A vyslo jako: %d\n", z);
    int idx = selected_individuals[z];
    //printf("Index je: %d\n", idx);
    Chromosome test_individual;
    test_individual.vector = malloc(DIMENSION * sizeof(double));
    for (int i = 0; i < DIMENSION; ++i) {
        test_individual.vector[i] = population->individuals[idx].vector[i];
    }


    for (int i = 0; i < DIMENSION; ++i) {
        for (int j = i+1; j < DIMENSION; ++j) {

            test_individual.vector[i] = a[i];
            test_individual.vector[j] = a[j];
            double diff_i = func(&test_individual);

            test_individual.vector[i] += b[i];
            diff_i -= func(&test_individual);

            test_individual.vector[i] = a[i];
            test_individual.vector[j] = a[j] + b[j];
            double diff_ij = func(&test_individual);

            test_individual.vector[i] += b[i];
            diff_ij -= func(&test_individual);

            diff_ij = fabs(diff_ij);
            diff_i = fabs(diff_i);
            if (diff_i - diff_ij >= 0) {
                D[i][j] = 1 - fabs(diff_ij/diff_i);
                D[j][i] = D[i][j];
            } else {
                D[i][j] = 1 - fabs(diff_i/diff_ij);
                D[j][i] = D[i][j];

            }
        }
    }

    return D;
}

double** createMIMatrix(Population* population) {
    //int* selected_individuals = individuals2Learning(population);
    int n = population->size * 0.2;
    int* selected_individuals = findNBest(population, n);
    double* means = calloc(DIMENSION, sizeof(double));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < DIMENSION; ++j) {
            means[j] += population->individuals[selected_individuals[i]].vector[j];
        }
    }
    for (int i = 0; i < DIMENSION; ++i) {
        means[i] /= n;
    }
    double** covariance_matrix = allocMatrix(DIMENSION);
    for (int r = 0; r < DIMENSION; ++r) {
        for (int c = r; c < DIMENSION; ++c) {
            for (int i = 0; i < n; ++i) {
                covariance_matrix[r][c] += (population->individuals[selected_individuals[i]].vector[r] - means[r]) *
                        (population->individuals[selected_individuals[i]].vector[c] - means[c]);
            }
            covariance_matrix[r][c] /= n;
            covariance_matrix[c][r] = covariance_matrix[r][c];
        }
    }
    free(means);

    double** mutual_information_matrix = allocMatrix(DIMENSION);
    for (int i = 0; i < DIMENSION; ++i) {
        for (int j = i; j < DIMENSION; ++j) {
            double R = covariance_matrix[i][j] / sqrt(covariance_matrix[i][i] * covariance_matrix[j][j]);
            mutual_information_matrix[i][j] = log2(sqrt(1 / (1 - pow(R, 2))));
            mutual_information_matrix[j][i] = mutual_information_matrix[i][j];
        }
    }
    freeSquareMatrix(covariance_matrix, DIMENSION);

    return mutual_information_matrix;
}

double** createMICMatrix(Population* population) {

    int n = population->size * 0.2;
    //printf("n: %d\n", n);
    int* selected_individuals = findNBest(population, n);
    
    mine_parameter param;
    param.alpha = 9;
    param.c = 5;
    param.est = EST_MIC_E;
    
    mine_matrix X;
    X.n = DIMENSION;
    X.m = n;
    X.data = (double *) malloc ((X.n * X.m) * sizeof(double));

    for (int j = 0; j < X.m; ++j) {
        for (int i = 0; i < X.n; ++i) {
            X.data[X.m * i + j] = population->individuals[selected_individuals[j]].vector[i];
        }
    }

    mine_pstats *pstats = mine_compute_pstats(&X, &param);


    int i, j, k, z;

    if (pstats == NULL)
    {
        printf("ERROR: pstats()\n");
        return 1;
    }
    /*
    printf("     ");
    for (j=1; j<X.n; j++) printf("X[%d]  ", j);
    printf("\n");
    k = 0;
    for (i=0; i<X.n-1; i++)
    {
        printf("X[%d] ", i);
        for (z=0; z<i; z++)
            printf("      ");

        for (j=i+1; j<X.n; j++)
        {
            printf("%.3lf ", pstats->mic[k]);
            k++;
        }
        printf("\n");
    }
    */

    double** mutual_information_matrix = allocMatrix(DIMENSION);
    for (int i = 0; i < DIMENSION; ++i) {
        mutual_information_matrix[i][i] = 0;
        for (int j = i + 1; j < DIMENSION; ++j) {
            mutual_information_matrix[i][j] = pstats->mic[DIMENSION * i - ((i*i+i)/2) - i - 1 + j];
            mutual_information_matrix[j][i] = mutual_information_matrix[i][j];
        }
    }
    free(pstats->mic);
    free(pstats->tic);
    free(pstats);

    free(X.data);

    return mutual_information_matrix;
}

void printFOS(FOS fos) {
    for (int i = 0; i < fos.number_of_subsets; ++i) {
        for (int j = 0; j < fos.parts[i].number_of_variables; ++j) {
            printf("%2d ", fos.parts[i].indexes[j]);
        }
        printf(" | ");
    }
    printf("\n");
}

int* argSort(double* values, int size) {
    int* sorted_indexes = malloc(size * sizeof(int));
    for (int i = 0; i < size; ++i) {
        sorted_indexes[i] = i;
    }

    for (int i = 0; i < size - 1; ++i) {
        int min_idx = i;
        for (int j = i+1; j < size; ++j) {
            if(values[j] < values[min_idx]) min_idx = j;
        }

        double min = values[min_idx];
        values[min_idx] = values[i];
        values[i] = min;
        int tmp = sorted_indexes[min_idx];
        sorted_indexes[min_idx] = sorted_indexes[i];
        sorted_indexes[i] = tmp;
    }

    return sorted_indexes;
}

int* findNBest(Population* population, int n) {
    double* fitness = malloc(population->size * sizeof(double));
    for (int i = 0; i < population->size; ++i) {
        fitness[i] = population->individuals[i].fitness;
    }

    int* sorted_idxs = argSort(fitness, population->size);
    free(fitness);
    int* ret_array = malloc(n * sizeof(int));
    for (int i = 0; i < n; ++i) {
        ret_array[i] = sorted_idxs[i];
    }
    free(sorted_idxs);

    return ret_array;
}

double** getOptimalDMatrix() {
    double** D = allocMatrix(DIMENSION);
    if (func == oneMax || func == sphere || func == Michalewicz || func == Rastrigin) return D;
    if (func == Rosenbrock) {
        for (int i = 0; i < DIMENSION - 1; ++i) {
            D[i][i+1] = 1;
            D[i+1][i] = 1;
        }
    }
    if (func == SoREB || func == OSoREB) {
        for (int i = 0; i < DIMENSION; i += size_of_fraction) {
            int real_size = size_of_fraction;
            if (i + real_size >= DIMENSION) real_size = DIMENSION - i;

            for (int j = 0; j < real_size; ++j) {
                for (int k = 0; k < real_size; ++k) {
                    if (j == k) continue;
                    D[i+j][i+k] = 1;
                }
            }
            if (func == OSoREB && i != 0) {
                D[i][i-1] = 1;
                D[i-1][i] = 1;
            }


        }
    }

    return D;
}



