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

void compute3DAdjMatrix(double**, double*, double*, double*, int, double, double, double);
double compute3DEclideanD(double, double, double, double, double, double);
void computeAdjDegree(double*, double**,  int);
void computeAdjLaplacian(double**, double**, double*, int);
void computeDisplacements(double*, double*, int);
void spatial1DDiffEq(double*, double*, double*, double, int);
void velocity1DDiffEq(double*, double**, double*, int);
double* computeSpatialConvergence1(double*, double*, double*, double*, double*, double*, double*, int);
double* computeSpatialConvergence(double*, double*, double*, double*, double*, double*, double*, int);
double computeV3DConvergence1(double*, double*, double*, int);
double* computeV3DConvergence(double*, double*, double*, double*, int);
bool checkConvergence1(int, double*, double, int);
bool checkConvergence(int, double*, double*, int, double);
void DumpResult1(string, double, int, double, int, double*, double);
void DumpResult(string, double, int, double, int, double*, double*, double*, double*, double*, double*, double*, double*);
void run3DFlockingModel(string, int, double, double, double, int, double);