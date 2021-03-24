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

#ifndef EPSILON
#define EPSILON 1.2e-7f
#endif

#ifndef MAX_ITER
#define MAX_ITER 100
#endif

#ifndef BETA_MAX
#define BETA_MAX 1.1009f
#endif

#ifndef BETA_MIN
#define BETA_MIN -0.01f
#endif

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

void computeAdjMatrix(float** adj, float* x, int N, float beta, float sigma_sqrd, float K){
      
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            float d = sigma_sqrd + pow(abs(x[i]-x[j]), 2.0);
            float val = K / pow(d, beta);
            adj[i][j] = val;
        }
    }

}

float compute3DEclideanD(float x_1, float x_2, float y_1, float y_2, float z_1, float z_2){
    float dx = x_1-x_2;
    float dy = y_1-y_2;
    float dz = z_1-z_2;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

void compute3DAdjMatrix(float** adj, float* x, float* y, float* z, int N, float beta, float sigma_sqrd, float K){
      
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            float d = sigma_sqrd + compute3DEclideanD(x[i], x[j], y[i], y[j], z[i], z[j]);
            float val = K / pow(d, beta);
            adj[i][j] = val;
        }
    }

}

void computeAdjDegree(float* degree, float** adj,  int N){

    for(int i=0; i<N; i++){
        float degree_i = 0.0f;
        for(int j=0; j<N; j++){
            degree_i += adj[i][j];
        }
        degree[i] = degree_i;
    }

}

void computeAdjLaplacian(float** lpc, float** adj, float* degree, int N){
    
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            if(i==j) lpc[i][j] = degree[i] - adj[i][j];
            else lpc[i][j] = 0.0f - adj[i][j];
        }
    }

}


void spatialDiffEq(float* x_t1, float* x_t, float* v_t, float dt, int N){

    for(int i=0; i<N; i++){
        x_t1[i] = x_t[i] + dt * v_t[i];
    }

}

void velocityDiffEq(float* v_t1, float** lplc, float* v_t, int N){

    for(int i=0; i<N; i++){
        float val = 0.0f;
        for(int j=0; j<N; j++){
            val += lplc[i][j] * v_t[j];
        }
        v_t1[i] = v_t[i] - val;
    }

}

float* computeSpatialConvergence(float* conv, float*x_1, float*y_1, float*z_1, float*x, float*y, float*z, int N){
    float conv_i = 0.0f;
    float conv_j = 0.0f;
    float conv_k = 0.0f;
    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            conv_i = abs(x_1[i] - x_1[j] - x[i] + x[j]);
            conv_j = abs(y_1[i] - y_1[j] - y[i] + y[j]);
            conv_k = abs(z_1[i] - z_1[j] - z[i] + z[j]);
            conv[0] = max(0.0f, conv_i);
            conv[1] = max(0.0f, conv_j);
            conv[2] = max(0.0f, conv_k);
        }
    }
    return conv;
}

float computeV3DConvergence(float*vx, float*vy, float*vz, int N){

    float max_convergence = 0.0f;
    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            float val = compute3DEclideanD(vx[i], vx[j], vy[i], vy[j], vz[i], vz[j]);
            if (val > max_convergence) max_convergence = val;
        }
    }
    return max_convergence;
}

bool checkConvergence(int iter, float* Xprecs, float Vprec, int max_iter){
    bool ret = true;
    if ((Vprec < EPSILON && Xprecs[0] < EPSILON && Xprecs[1] < EPSILON && Xprecs[2] < EPSILON) ||
    iter > max_iter) ret = false;
    return ret;
}

void DumpPrecision(string filename, float *beta, int *counter, float *precision, int N){
    // substitute this encoded path to yours
    filename = "/home/relogu/OneDrive/UNIBO/Magistrale/Physics of Complex Systems/project/convergence/" + filename;
    ofstream filepp;
    filepp.open(filename);
    //cout << "Dumping the file " << filename << endl;
    filepp << "beta,iter,prec" << endl;
    for(int i=0;i<N;i++){//for-loop on the i-th row
        filepp << beta[i] << setprecision(4);
        filepp << ",";
        filepp << counter[i];
        filepp << ",";
        filepp << precision[i] << setprecision(4);
        if (i != N-1) filepp << endl;
    }
    filepp.close();
}

void DumpPerformance(string filename, float time, int N_sim){
    // substitute this encoded path to yours
    filename = "/home/relogu/OneDrive/UNIBO/Magistrale/Physics of Complex Systems/project/performance/" + filename;
    ofstream filepp;
    filepp.open(filename);
    //cout << "Dumping the file " << filename << endl;
    filepp << "time_elapsed (ms),number_of_iterations;" << endl;
    filepp << time << setprecision(8);
    filepp << ",";
    filepp << N_sim;
    filepp.close();
}

void DumpResult(string filename, float beta, int iter, float* Xprec, float Vprec){
    // substitute this encoded path to yours
    ofstream filepp;
    filepp.open(filename, ofstream::app);
    ostringstream tmp;
    
    tmp << beta << setprecision(4);
    tmp << ", ";
    tmp << iter;
    tmp << ", ";
    for (int i=0; i<3; i++) {
        tmp << Xprec[i] << setprecision(4);
        if (i<2) tmp << ",";
    }
    tmp << ", ";
    tmp << Vprec << setprecision(4);
    tmp << endl;
    cout << tmp.str().c_str();
    filepp << tmp.str().c_str();
        
    filepp.close();
}

void run3DFlockingModel(int birds, int jumps, float dt, float sigma_sqrd, float K){
    // Seed with a real random value, if available
    //random_device r;
    // Retrieve a seed sequence fror next RNGs
    //seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
        
    //definition of chronometers
    chrono::steady_clock::time_point begin, end;
    begin = chrono::steady_clock::now();

    int N;
    int reps;
    float m_dt;
    float m_sigma_sqrd;
    float m_K;
    float time;
    int max_iter;
    ostringstream tmp;
    reps = jumps;

    // first step computation
    m_dt = dt;
    float beta;
    beta = BETA_MAX;
    float d = (BETA_MAX - BETA_MIN) / reps;
    m_sigma_sqrd = sigma_sqrd;
    m_K = K;
    N = birds;
    max_iter = 10*N;
    tmp << "convergence/3D_N" << N << ".out";

    ofstream filepp;
    filepp.open(tmp.str().c_str(), ofstream::trunc);
    filepp << "beta,iter,Xprec_i,Xprec_j,Xprec_k,Vprec" << endl;
    filepp.close();
    cout << "Preliminary operations done" << endl;

    #pragma omp parallel for
    for (int i=0; i<reps; i++){
        beta -= d;
        //cout << "Computing simulation for beta " << beta << endl;
        // instantiate initial (random) conditions
        float *x, *y, *z;
        float *vx, *vy, *vz;
        float *x_1, *y_1, *z_1;
        float *vx_1, *vy_1, *vz_1;
        float **adj;
        float *degree;
        float **lapl;

        int seed = 51550;
        // Instantiate a the Marsenne-Twister PRNG with a seed of the sequence
        mt19937 generator(seed);

        initializeRandomVector(&x, N, &generator);
        initializeRandomVector(&y, N, &generator);
        initializeRandomVector(&z, N, &generator);
        initializeRandomVector(&vx, N, &generator);
        initializeRandomVector(&vy, N, &generator);
        initializeRandomVector(&vz, N, &generator);
        initializeNullVector(&x_1, N);
        initializeNullVector(&y_1, N);
        initializeNullVector(&z_1, N);
        initializeNullVector(&vx_1, N);
        initializeNullVector(&vy_1, N);
        initializeNullVector(&vz_1, N);
        initializeNullVector(&degree, N);
        initializeNullMatrix(&adj, N);
        initializeNullMatrix(&lapl, N);

        float* Xmax_convergence = new float[3];
        float Vmax_convergence = 0.0f;
        int iter = 0;
        do {
            compute3DAdjMatrix(adj, x, y, z, N, beta, m_sigma_sqrd, m_K);
            computeAdjDegree(degree, adj, N);
            computeAdjLaplacian(lapl, adj, degree, N);
            spatialDiffEq(x_1, x, vx, m_dt, N);
            spatialDiffEq(y_1, y, vy, m_dt, N);
            spatialDiffEq(z_1, z, vz, m_dt, N);
            velocityDiffEq(vx_1, lapl, vx, N);
            velocityDiffEq(vy_1, lapl, vy, N);
            velocityDiffEq(vz_1, lapl, vz, N);
            Vmax_convergence = computeV3DConvergence(vx_1, vy_1, vz_1, N);
            Xmax_convergence = computeSpatialConvergence(Xmax_convergence, x_1, y_1, z_1, x, y, z, N);
            
            swap(x, x_1);
            swap(y, y_1);
            swap(z, z_1);
            swap(vx, vx_1);
            swap(vy, vy_1);
            swap(vz, vz_1);

            iter++;
        } while (checkConvergence(iter, Xmax_convergence, Vmax_convergence, max_iter));
        
        cout << "Computed simulation for beta " << beta << endl;

        if (iter >= MAX_ITER) {
            cout << "The convergence NOT reached after " << MAX_ITER << " iterations" << endl;
        } else {
            cout << "The FINAL number of iterations reached is " << iter << endl;
        }
        
        cout << "The FINAL max convergence reached is V: " << Vmax_convergence 
        << " -- X: (" << Xmax_convergence[0] << ", " << Xmax_convergence[1] << ", " << Xmax_convergence[2] << ")" << endl;
        cout << "Dumping results" << endl;
        DumpResult(tmp.str().c_str(), beta, iter, Xmax_convergence, Vmax_convergence);
    }
}

int main()
{
    int N = 10;
    float K = -0.00001f + 1 / (sqrt(N)*(N-1));
    run3DFlockingModel(N, 200, 1.0f, 1.0f, K);

    N = 30;
    K = -0.00001f + 1 / (sqrt(N)*(N-1));
    run3DFlockingModel(N, 200, 1.0f, 1.0f, K);

    N = 100;
    K = -0.00001f + 1 / (sqrt(N)*(N-1));
    run3DFlockingModel(N, 200, 1.0f, 1.0f, K);
}