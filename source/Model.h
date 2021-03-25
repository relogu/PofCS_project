#include "Initializers.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>
#include <fstream>
#include <string>
#include <chrono>

using namespace std;

void compute3DAdjMatrix(float**, float*, float*, float*, int, float, float, float);
float compute3DEclideanD(float, float, float, float, float, float);
void computeAdjDegree(float*, float**,  int);
void computeAdjLaplacian(float**, float**, float*, int);
void spatial1DDiffEq(float*, float*, float*, float, int);
void velocity1DDiffEq(float*, float**, float*, int);
float* computeSpatialConvergence(float*, float*, float*, float*, float*, float*, float*, int);
float computeV3DConvergence(float*, float*, float*, int);
bool checkConvergence(int, float*, float, int);
void DumpResult(string, float, int, float, int, float*, float);
void run3DFlockingModel(int, int, float, float);
void draw3DFlockingModel(int, int, float, float);