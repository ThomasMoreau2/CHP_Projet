
#include "fonction.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

// Cas stationnaire 1
double f1(double x, double y, double t, double Lx, double Ly) {
    return 2 * (y - y*y + x - x*x);
}

double g1(double, double) {
    return 0.0;
}

double h1(double, double) {
    return 0.0;
}

// Cas stationnaire 2
double f2(double x, double y, double t, double Lx, double Ly) {
    return sin(x) + cos(y);
}

double g2(double x, double y) {
    return sin(x) + cos(y);
}

double h2(double x, double y) {
    return sin(x) + cos(y);
}

// Cas instationnaire périodique
double f3(double x, double y, double t, double Lx, double Ly) {
    double dx = x - (Lx / 2.0);
    double dy = y - (Ly / 2.0);
    return exp(-dx*dx) * exp(-dy*dy) * cos(M_PI * t / 2.0);
}

double g3(double, double) {
    return 0.0;
}

double h3(double, double) {
    return 1.0;
}

//Solution exacte cas 1
double f1_ex(double x,double y){
    return x*(1-x)*y*(1-y);
}

//Solution exacte cas 2

double f2_ex(double x,double y){
    return cos(y)+sin(x);
}


double norme2(vector<double>& u) {
    double somme = 0.0;
    for (double val : u) {
        somme += val * val;
    }
    return sqrt(somme);
}