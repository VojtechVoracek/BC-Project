//
// Created by vojta on 10.11.20.
//

#include "parameters.h"
#include "fitness_functions.h"
#include <math.h>


#include <stdbool.h>

int DIMENSION = 20;
double CR = 0.9;
double F = 0.7;
double (*func)(void*) = &oneMax;
double TERMINAL_CONDITION = 0.00000001;



