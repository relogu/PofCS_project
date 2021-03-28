#include "Initializers.h"
#include <random>
#include <cmath>

using namespace std;

void initializeRandomMatrix(double*** matrix, int N, mt19937* generator){

    (*matrix) = new double*[N];
    for(int i=0; i<N; i++){
        (*matrix)[i] = new double[N];
    }

    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
        uniform_real_distribution<double> unif(0.0,10.0);
            (*matrix)[i][j] = unif(*generator);
        }
    }
    
}

void initializeNullMatrix(double*** matrix, int N){

    (*matrix) = new double*[N];
    for(int i=0; i<N; i++){
        (*matrix)[i] = new double[N];
    }

    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            (*matrix)[i][j] = 0.0;
        }
    }
    
}

void initializeRandomVector(double** vector, int N, mt19937* generator){

    (*vector) = new double[N];
    for(int i=0; i<N; i++){
        uniform_real_distribution<double> unif(-5.0,5.0);
        (*vector)[i] = unif(*generator);
    }
    
}

void initializeNullVector(double** vector, int N){

    (*vector) = new double[N];
    for(int i=0; i<N; i++){
        (*vector)[i] = 0.0;
    }
    
}