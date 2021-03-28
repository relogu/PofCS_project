#include "Model.h"
#include <cmath>

using namespace std;

int main(){
    double dt = 1.0;
    double sigma = 1.0;
    double beta = 1.0;
    int N = 10;
    // data for studying the dependence on beta
    #pragma omp parallel for
    for (int i=1; i < 200; i++){
        double b = -0.005*i + beta;
        run3DFlockingModel(N, dt, sigma, b);
    }
    beta = 0.4;
    sigma = 1.0;
    //TODO: data for studying the dependence on N
    N = 2;
    #pragma omp parallel for
    for (int i=0; i < 100; i+=5) run3DFlockingModel(N*(1+i), dt, sigma, beta);
    // data for studying the dependence on sigma
    #pragma omp parallel for
    for (int i=1; i < 200; i++){
        double s = exp(0.1*i) - 1.0 + sigma;
        run3DFlockingModel(N, dt, s, beta);
        run3DFlockingModel(N, dt, 1/s, beta);
    }
    /*
    for (int i=N; i<=100; i+=10) run3DFlockingModel(i, rep, dt, sigma, beta);
    for (int i=0; i<10; i++) run3DFlockingModel(100, rep, dt, 0.1f+(i*0.1f), beta);
    */   
}