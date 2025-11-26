#include "produit.h"
#include <vector>
#include <fstream>

using namespace std;


vector<double> prod_mat_vec(int Nx, int Ny, double D, double dx, double dy, double dt, vector<double> x) 
{   
    int n;
    double a, lambda_x, lambda_y;

    n = Nx*Ny;
    lambda_x = D*dt/(dx*dx);
    lambda_y = D*dt/(dy*dy);
    a = 1 +2*lambda_x+2*lambda_y;

    vector<double> y(n);

    y[0] = a*x[0] - lambda_y*x[1] - lambda_x*x[Ny];
    for (int i=1; i<Ny-1; i++)
    {
        y[i] = -lambda_y*x[i-1] + a*x[i] - lambda_y*x[i+1]-lambda_x*x[Ny+i];
    }
    y[Ny-1] = -lambda_y*x[Ny-2] + a*x[Ny-1] - lambda_x*x[Ny-1+Ny];

    for (int i=Ny; i<=n-Ny-1; i++)
    {   
        if ((i+1)%Ny==0)
        {
            y[i] = -lambda_x*x[i-Ny] - lambda_y*x[i-1] + a*x[i] -lambda_x*x[Ny+i];
        } 
        else if (i%Ny==0)
        {
            y[i] = -lambda_x*x[i-Ny] + a*x[i] - lambda_y*x[i+1]-lambda_x*x[Ny+i];
        }
        else 
        {
            y[i] = -lambda_x*x[i-Ny]-lambda_y*x[i-1] + a*x[i] - lambda_y*x[i+1]-lambda_x*x[Ny+i];
        }
    }

    y[n-Ny] = -lambda_x*x[n-Ny-Ny] + a*x[n-Ny] - lambda_y*x[n-Ny+1];
    for (int i=n-Ny+1; i<n-1; i++)
    {
        y[i] = -lambda_x*x[i-Ny]-lambda_y*x[i-1] + a*x[i] - lambda_y*x[i+1];
    }
    y[n-1] = -lambda_x*x[n-1-Ny]-lambda_y*x[n-2] + a*x[n-1];

    return y;
}

vector<double> make_b(vector<double> u, int Nx, int Ny, double dx, double dy, double Lx, double Ly, double D, double dt, double t, int cas, 
                int rank, int nproc, vector<double> u_gauche, vector<double> u_droit)
{   
    double (*f)(double, double, double, double, double);
    double (*g)(double, double);
    double (*h_g)(double, double);
    double (*h_d)(double, double);

   if (cas == 1) 
    {
        f = f1;
        g = g1;
        h_g = h1;
        h_d = h1;
    } 
    else if (cas == 2) 
    {
        f = f2;
        g = g2;
        h_g = h2;
        h_d = h2;
    } 
    else if (cas == 3) 
    {
        f = f3;
        g = g3;
        h_g = h3;
        h_d = h3;
    } 

    int n;
    n = Nx*Ny;

    double lambda_x, lambda_y;
    lambda_x = D*dt/(dx*dx);
    lambda_y = D*dt/(dy*dy);

    vector<double> b(n);

    if (rank == 0)
    {
        for(int i=1; i<=Ny; i++){u[n - i] += u_droit[n - i];}
        h_d = null_f;

    }

    else if ((rank >0) && (rank <nproc-1))
    {
        h_d = null_f;
        h_g = null_f;
        for(int i=0; i<Ny; i++){u[i] += u_gauche[i];}
        for(int i=1; i<=Ny; i++){u[n - i] += u_droit[n - i];}
    }

    else if (rank == nproc-1)
    {
        h_g = null_f;
        for(int i=0; i<Ny; i++){u[i] += u_gauche[i];}
    }
    


    b[0] = u[0] + dt*f(dx, dy, t, Lx, Ly) + lambda_x*h_g(0, dy) + lambda_y*g(dx, 0);
    for (int i=1; i<Ny-1; i++)
    {
        b[i] = u[i] + dt*f((i/Ny+1)*dx, (i%Ny+1)*dy, t, Lx, Ly) + lambda_x*h_g(0, (i%Ny+1)*dy);
    }
    b[Ny-1] = u[Ny-1] + dt*f(dx, Ny*dy, t, Lx, Ly) + lambda_x*h_g(0, Ny*dy) + lambda_y*g(dx, (Ny+1)*dy);

    for (int i=Ny; i<=n-Ny-1; i++)
    {   
        if ((i+1)%Ny==0)
        {
            b[i] = u[i] + dt*f((i/Ny+1)*dx, (i%Ny+1)*dy, t, Lx, Ly) + lambda_y*g((i/Ny+1)*dx, (Ny+1)*dy);
        } 
        else if (i%Ny==0)
        {
            b[i] = u[i] + dt*f((i/Ny+1)*dx, (i%Ny+1)*dy, t, Lx, Ly) + lambda_y*g((i/Ny+1)*dx, 0);
        }
        else 
        {
            b[i] = u[i] + dt*f((i/Ny+1)*dx, (i%Ny+1)*dy, t, Lx, Ly);
        }
    }

    b[n-Ny] = u[n-Ny] + dt*f(((n-Ny)/Ny+1)*dx, ((n-Ny)%Ny+1)*dy, t, Lx, Ly) + lambda_x*h_d((Nx+1)*dx, ((n-Ny)%Ny+1)*dy) + lambda_y*g(((n-Ny)/Ny+1)*dx, 0);
    for (int i=n-Ny+1; i<n-1; i++)
    {
        b[i] = u[i] + dt*f((i/Ny+1)*dx, (i%Ny+1)*dy, t, Lx, Ly) + lambda_x*h_d((Nx+1)*dx, (i%Ny+1)*dy);
    }
    b[n-1] = u[n-1] + dt*f(Nx*dx, Ny*dy, t, Lx, Ly) + lambda_x*h_d((Nx+1)*dx, Ny*dy) + lambda_y*g(Nx*dx, (Ny+1)*dy);

    return b;
}
