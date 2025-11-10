#ifndef GRAD_CONJ_H
#define GRAD_CONJ_H

#include <vector>

std::vector<double> gradient_conjugue(
    int Nx, int Ny, double D, double dx, double dy, double dt,
    const std::vector<double>& b, double tol = 1e-8, int max_iter = 1000
);

#endif
