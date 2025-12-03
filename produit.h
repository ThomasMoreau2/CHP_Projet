#ifndef PRODUIT_H
#define PRODUIT_H

#include <vector>
#include "fonction.h"

std::vector<double> prod_mat_vec(int Nx, int Ny, double D, double dx, double dy, double dt, std::vector<double> x);

std::vector<double> make_b(std::vector<double> u, int Nx, int Ny, double dx, double dy, double Lx, double Ly, double D, double dt, double t, int cas,
                    int rank, int nproc, int iBeg, std::vector<double> u_gauche, std::vector<double> u_droit, int r);


#endif