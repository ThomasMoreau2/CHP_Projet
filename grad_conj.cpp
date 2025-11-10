#include "grad_conj.h"
#include "produit.h"
#include <cmath>
#include <iostream>

using namespace std;

// Produit scalaire entre deux vecteurs
static double dot(const vector<double>& a, const vector<double>& b) {
    double res = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        res += a[i] * b[i];
    }
    return res;
}

// Gradient conjugué
vector<double> gradient_conjugue(
    int Nx, int Ny, double D, double dx, double dy, double dt,
    const vector<double>& b, double tol, int max_iter
) {
    int n = Nx * Ny;
    vector<double> x(n, 0.0); // x₀
    vector<double> r = b;
    vector<double> p = r;
    vector<double> Ap(n);

    double rsold = dot(r, r);

    for (int k = 0; k < max_iter; ++k) {
        Ap = prod_mat_vec(Nx, Ny, D, dx, dy, dt, p);

        double alpha = rsold / dot(p, Ap);

        for (int i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double rsnew = dot(r, r);
        if (sqrt(rsnew) < tol) {
            // std::cout << "Convergence en " << k + 1 << " itérations." << std::endl;
            break;
        }

        double beta = rsnew / rsold;
        for (int i = 0; i < n; ++i) {
            p[i] = r[i] + beta * p[i];
        }

        rsold = rsnew;
    }

    return x;
}
