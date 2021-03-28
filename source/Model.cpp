#include "Model.h"
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

#ifndef EPSILON
#define EPSILON 1.2e-7f
#endif

#ifndef BETA_MAX
#define BETA_MAX 1.1009f
#endif

#ifndef BETA_MIN
#define BETA_MIN -0.01f
#endif

void computeAdjMatrix(double** adj, double* x, int N, double beta, double sigma_sqrd, double K){
      
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            double d = sigma_sqrd + pow(abs(x[i]-x[j]), 2.0);
            double val = K / pow(d, beta);
            adj[i][j] = val;
        }
    }

}

double compute3DEclideanD(double x_1, double x_2, double y_1, double y_2, double z_1, double z_2){
    double dx = x_1-x_2;
    double dy = y_1-y_2;
    double dz = z_1-z_2;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

void compute3DAdjMatrix(double** adj, double* x, double* y, double* z, int N, double beta, double sigma_sqrd, double K){
      
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            double d = sigma_sqrd + compute3DEclideanD(x[i], x[j], y[i], y[j], z[i], z[j]);
            double val = K / pow(d, beta);
            adj[i][j] = val;
        }
    }

}

void computeAdjDegree(double* degree, double** adj,  int N){

    for(int i=0; i<N; i++){
        double degree_i = 0.0;
        for(int j=0; j<N; j++){
            degree_i += adj[i][j];
        }
        degree[i] = degree_i;
    }

}

void computeAdjLaplacian(double** lpc, double** adj, double* degree, int N){
    
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            if(i==j) lpc[i][j] = degree[i] - adj[i][j];
            else lpc[i][j] = - adj[i][j];
        }
    }

}

void computeDisplacements(double* displacements, double* x,  int N){

    for(int i=0; i<N; i++){
        double degree_i = 0.0;
        for(int j=i+1; j<N; j++){
            displacements[i+j] = x[i] - x[j];
        }
    }

}


void spatialDiffEq(double* x_t1, double* x_t, double* v_t, double dt, int N){

    for(int i=0; i<N; i++){
        x_t1[i] = x_t[i] + dt * v_t[i];
    }

}

void velocityDiffEq(double* v_t1, double** lplc, double* v_t, int N){

    for(int i=0; i<N; i++){
        double val = 0.0;
        for(int j=0; j<N; j++){
            val += lplc[i][j] * v_t[j];
        }
        v_t1[i] = v_t[i] - val;
    }

}

double* computeSpatialConvergence1(double* conv, double*x_1, double*y_1, double*z_1, double*x, double*y, double*z, int N){
    double conv_i = 0.0;
    double conv_j = 0.0;
    double conv_k = 0.0;
    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            conv_i = x_1[i] - x_1[j] - x[i] + x[j];
            conv_j = y_1[i] - y_1[j] - y[i] + y[j];
            conv_k = z_1[i] - z_1[j] - z[i] + z[j];
            if(conv_i < 0.0) conv_i = - conv_i;
            if(conv_j < 0.0) conv_j = - conv_j;
            if(conv_k < 0.0) conv_k = - conv_k;
            conv[0] = max(0.0, conv_i);
            conv[1] = max(0.0, conv_j);
            conv[2] = max(0.0, conv_k);
        }
    }
    return conv;
}

double* computeSpatialConvergence(double* conv, double* disp_x_1, double* disp_y_1, double* disp_z_1, double* disp_x, double* disp_y, double* disp_z, int N){
    double conv_i = 0.0;
    double conv_j = 0.0;
    double conv_k = 0.0;
    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            conv_i = disp_x_1[i+j] - disp_x[i+j];
            conv_j = disp_y_1[i+j] - disp_y[i+j];
            conv_k = disp_z_1[i+j] - disp_z[i+j];
            if(conv_i < 0.0) conv_i = - conv_i;
            if(conv_j < 0.0) conv_j = - conv_j;
            if(conv_k < 0.0) conv_k = - conv_k;
            conv[0] = max(0.0, conv_i);
            conv[1] = max(0.0, conv_j);
            conv[2] = max(0.0, conv_k);
        }
    }
    return conv;
}

double* computeV3DConvergence(double* conv, double*vx, double*vy, double*vz, int N){
    double conv_i = 0.0;
    double conv_j = 0.0;
    double conv_k = 0.0;
    double max_convergence = 0.0;
    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            conv_i = vx[i] - vx[j];
            conv_j = vy[i] - vy[j];
            conv_k = vz[i] - vz[j];
            if(conv_i < 0.0) conv_i = - conv_i;
            if(conv_j < 0.0) conv_j = - conv_j;
            if(conv_k < 0.0) conv_k = - conv_k;
            conv[0] = max(0.0, conv_i);
            conv[1] = max(0.0, conv_j);
            conv[2] = max(0.0, conv_k);
        }
    }
    return conv;
}

bool checkConvergence1(int iter, double* Xprecs, double Vprec, int max_iter){
    bool ret = true;
    if ((Vprec < EPSILON && Xprecs[0] < EPSILON && Xprecs[1] < EPSILON && Xprecs[2] < EPSILON) ||
    iter > max_iter) ret = false;
    return ret;
}

bool checkConvergence(int iter, double* Xprecs, double* Vprec, int max_iter){
    bool ret = true;
    if ((Vprec[0] < EPSILON && Vprec[1] < EPSILON && Vprec[2] < EPSILON
            && Xprecs[0] < EPSILON && Xprecs[1] < EPSILON && Xprecs[2] < EPSILON) ||
        iter > max_iter) ret = false;
    return ret;
}

void DumpResult1(string filename, double sigma, int N, double beta, int iter, double* Xprec, double Vprec){
    // substitute this encoded path to yours
    ofstream filepp;
    filepp.open(filename, ofstream::app);
    ostringstream tmp;
    
    tmp << sigma << setprecision(4);
    tmp << ", ";
    tmp << N;
    tmp << ", ";
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

void DumpResult(string filename, double sigma, int N, double beta, int iter, double* Xprec, double* Vprec,
    double* x, double* y, double* z, double* vx, double* vy, double* vz){
    // substitute this encoded path to yours
    ofstream filepp;
    filepp.open(filename, ofstream::app);
    ostringstream tmp;
    
    tmp << sigma << setprecision(4);
    tmp << ",";
    tmp << N;
    tmp << ",";
    tmp << beta << setprecision(4);
    tmp << ",";
    tmp << iter;
    tmp << ",";
    for (int i=0; i<3; i++) {
        tmp << Xprec[i] << setprecision(4);
        if (i != 2) tmp << " ";
    }
    tmp << ",";
    for (int i=0; i<3; i++) {
        tmp << Vprec[i] << setprecision(4);
        if (i != 2) tmp << " ";
    }
    tmp << ",";
    for (int j=0; j<N; j++){
        tmp << x[j] << setprecision(4);
        tmp << " ";
        tmp << y[j] << setprecision(4);
        tmp << " ";
        tmp << z[j] << setprecision(4);
        if (j != N*(N-1)-1) tmp << " ";
    }
    tmp << ",";
    for (int j=0; j<N; j++){
        tmp << vx[j] << setprecision(4);
        tmp << " ";
        tmp << vy[j] << setprecision(4);
        tmp << " ";
        tmp << vz[j] << setprecision(4);
        tmp << "";
        if (j != N-1) tmp << " ";
    }
    tmp << "";
    tmp << endl;
    //cout << tmp.str().c_str();
    filepp << tmp.str().c_str();
        
    filepp.close();
}

void run3DFlockingModel(int birds, double dt, double sigma, double beta){
    // Seed with a real random value, if available
    //random_device r;
    // Retrieve a seed sequence fror next RNGs
    //seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
        
    //definition of chronometers
    chrono::steady_clock::time_point begin, end;
    begin = chrono::steady_clock::now();

    int N;
    int reps;
    double m_dt;
    double m_sigma_sqrd;
    double m_beta;
    double time;
    int max_iter;
    double K_den;
    double m_K;
    ostringstream tmp1;

    // first step computation
    m_dt = dt;
    m_beta = beta;
    m_sigma_sqrd = sigma*sigma;
    N = birds;
    K_den = (N-1)*sqrt(N);
    m_K = pow(m_sigma_sqrd, m_beta)/K_den; // - EPSILON;
    max_iter = 100*N;
    tmp1 << "../results/" << N << "_" << sigma << "_" << m_beta << ".dat";

    ofstream filepp;
    filepp.open(tmp1.str().c_str(), ofstream::trunc);
    filepp << "sigma,N,beta,iter,[space_conv],[vel_conv],[[positions]],[[velocities]]" << endl;
    filepp.close();
    cout << "Preliminary operations done" << endl;

    //cout << "Computing simulation for beta " << beta << endl;
    // instantiate initial (random) conditions
    double *x, *y, *z;
    double *vx, *vy, *vz;
    double *x_1, *y_1, *z_1;
    double *vx_1, *vy_1, *vz_1;
    double **adj;
    double *degree;
    double **lapl;
    double *displacements_x, *displacements_y, *displacements_z;
    double *displacements_x_1, *displacements_y_1, *displacements_z_1;

    // Instantiate a the Marsenne-Twister PRNG with a seed of the sequence
    int seed = 51550;
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
    initializeNullVector(&displacements_x, (int)N*(N-1));
    initializeNullVector(&displacements_y, (int)N*(N-1));
    initializeNullVector(&displacements_z, (int)N*(N-1));
    initializeNullVector(&displacements_x_1, (int)N*(N-1));
    initializeNullVector(&displacements_y_1, (int)N*(N-1));
    initializeNullVector(&displacements_z_1, (int)N*(N-1));
    initializeNullMatrix(&adj, N);
    initializeNullMatrix(&lapl, N);

    double* X_convergence = new double[3];
    double* V_convergence = new double[3];
    int iter = 0;
    do {
        iter++;
        
        compute3DAdjMatrix(adj, x, y, z, N, m_beta, m_sigma_sqrd, m_K);
        computeAdjDegree(degree, adj, N);
        computeAdjLaplacian(lapl, adj, degree, N);
        spatialDiffEq(x_1, x, vx, m_dt, N);
        spatialDiffEq(y_1, y, vy, m_dt, N);
        spatialDiffEq(z_1, z, vz, m_dt, N);
        velocityDiffEq(vx_1, lapl, vx, N);
        velocityDiffEq(vy_1, lapl, vy, N);
        velocityDiffEq(vz_1, lapl, vz, N);
        computeDisplacements(displacements_x_1, x_1, N);
        computeDisplacements(displacements_y_1, y_1, N);
        computeDisplacements(displacements_z_1, z_1, N);

        X_convergence = computeSpatialConvergence(X_convergence, displacements_x_1,displacements_y_1,
            displacements_z_1,displacements_x,displacements_y,displacements_z,N);
        V_convergence = computeV3DConvergence(V_convergence, vx_1, vy_1, vz_1, N);
        DumpResult(tmp1.str().c_str(), sigma, N, m_beta, iter, X_convergence, V_convergence,
            x_1,y_1,z_1,vx_1,vy_1,vz_1);

        swap(x, x_1);
        swap(y, y_1);
        swap(z, z_1);
        swap(vx, vx_1);
        swap(vy, vy_1);
        swap(vz, vz_1);
        swap(displacements_x, displacements_x_1);
        swap(displacements_y, displacements_y_1);
        swap(displacements_z, displacements_z_1);

    } while (checkConvergence(iter, X_convergence, V_convergence, max_iter));


    ostringstream tmp;
    
    tmp << sigma << setprecision(4);
    tmp << ",";
    tmp << N;
    tmp << ",";
    tmp << beta << setprecision(4);
    tmp << ",";
    tmp << iter;
    tmp << ",[";
    for (int i=0; i<3; i++) {
        tmp << X_convergence[i] << setprecision(4);
        if (i<2) tmp << " ";
    }
    tmp << "],[";
    for (int i=0; i<3; i++) {
        tmp << V_convergence[i] << setprecision(4);
        if (i<2) tmp << " ";
    }
    tmp << "],[";
    for (int j=0; j<N;j++){
        tmp << "[";
        tmp << x[j] << setprecision(4);
        tmp << " ";
        tmp << y[j] << setprecision(4);
        tmp << " ";
        tmp << z[j] << setprecision(4);
        tmp << "]";
        if (j != N-1) tmp << " ";
    }
    tmp << "],[";
    for (int j=0; j<N; j++){
        tmp << "[";
        tmp << vx[j] << setprecision(4);
        tmp << " ";
        tmp << vy[j] << setprecision(4);
        tmp << " ";
        tmp << vz[j] << setprecision(4);
        tmp << "]";
        if (j != N-1) tmp << " ";
    }
    tmp << "]";
    tmp << endl;
    //cout << tmp.str().c_str();
    

}