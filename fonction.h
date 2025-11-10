#ifndef FONCTION_H
#define FONCTION_H

#include <vector>

// Cas stationnaire 1
double f1(double x, double y, double t, double Lx, double Ly);
double g1(double x, double y);
double h1(double x, double y);

// Cas stationnaire 2
double f2(double x, double y, double t, double Lx, double Ly);
double g2(double x, double y);
double h2(double x, double y);

// Cas instationnaire périodique
double f3(double x, double y, double t, double Lx, double Ly);
double g3(double x, double y);
double h3(double x, double y);

//Solution exacte 
double f1_ex(double x,double y);
double f2_ex(double x,double y);

double norme2(std::vector<double>& u);

#endif
