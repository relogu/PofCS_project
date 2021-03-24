#include "Model.h"
#include <cmath>

using namespace std;

int main(){
    int rep = 200;
    float dt = 1.0f;
    float sigma = 1.0f;
    int N = 10;
    for (int i=N; i<=100; i+=10) run3DFlockingModel(i, rep, dt, sigma);
    for (int i=0; i<10; i++) run3DFlockingModel(100, rep, dt, 0.1f+(i*0.1f));
}