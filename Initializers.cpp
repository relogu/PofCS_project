#include "Initializers.h"
#include <random>
#include <cmath>

using namespace std;

void initializeRandomMatrix(float*** matrix, int N, mt19937* generator){

    (*matrix) = new float*[N];
    for(int i=0; i<N; i++){
        (*matrix)[i] = new float[N];
    }

    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
        uniform_real_distribution<float> unif(0.0,10.0);
            (*matrix)[i][j] = unif(*generator);
        }
    }
    
}

void initializeNullMatrix(float*** matrix, int N){

    (*matrix) = new float*[N];
    for(int i=0; i<N; i++){
        (*matrix)[i] = new float[N];
    }

    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            (*matrix)[i][j] = 0.0f;
        }
    }
    
}

void initializeRandomVector(float** vector, int N, mt19937* generator){

    (*vector) = new float[N];
    for(int i=0; i<N; i++){
        uniform_real_distribution<float> unif(-1.0,1.0);
        (*vector)[i] = unif(*generator);
    }
    
}

void initializeNullVector(float** vector, int N){

    (*vector) = new float[N];
    for(int i=0; i<N; i++){
        (*vector)[i] = 0.0f;
    }
    
}