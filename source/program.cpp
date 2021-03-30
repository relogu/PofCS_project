#include "Model.h"
#include <cmath>
#include <string>

using namespace std;

int main(){
    double dt = 1.0;
    double sigma = 1.0;
    double beta = 1.0;
    int N = 10;
    int max_iter = 50*N;//old 100*N
    double epsilon = 1.0e-5f;
    string path = "../results1/";

    // data for studying the dependence on beta
    #pragma omp parallel for
    for (int i=0; i < 200; i++){
        double b = -0.005*i + beta;
        run3DFlockingModel(path, N, dt, sigma, b, max_iter, epsilon);
    }
    path = "../results4/";
    beta = 0.5;
    #pragma omp parallel for
    for (int i=0; i < 200; i++){
        int m = 50*N - i;
        double e = epsilon*m/(50*N);
        run3DFlockingModel(path, N, dt, sigma, beta, m, e);
    }
    path = "../results5/";
    beta = 0.5;
    N = 2;
    #pragma omp parallel for
    for (int i=0; i < 200; i++) run3DFlockingModel(path, N+i, dt, sigma, beta, 50*(N), epsilon);


    // data for studying the acceptance of the model

    beta = 0.4;
    sigma = 1.0;
    
    //data for studying the dependence on N
    N = 2;
    path = "../results2/";
    #pragma omp parallel for
    for (int i=0; i < 200; i++) run3DFlockingModel(path, N+i, dt, sigma, beta, 50*(N+i), epsilon);
    
    // data for studying the dependence on sigma
    N = 10;
    path = "../results3/";
    #pragma omp parallel for
    for (int i=1; i < 100; i++){
        double s = exp(0.1*i) - 1.0 + sigma;
        run3DFlockingModel(path, N, dt, s, beta, 100*(N), epsilon);
        run3DFlockingModel(path, N, dt, 1/s, beta, 100*(N), epsilon);
    }

    /*
    for (int i=N; i<=100; i+=10) run3DFlockingModel(i, rep, dt, sigma, beta);
    for (int i=0; i<10; i++) run3DFlockingModel(100, rep, dt, 0.1f+(i*0.1f), beta);
    */   
}