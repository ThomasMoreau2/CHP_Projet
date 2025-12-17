#include "fonction.h"
#include "produit.h"
#include "grad_conj.h"
#include "charge.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[])
{   
    int Nx, Ny, n_iter, N, cas, r, nproc, iBeg, iEnd, rank, N_loc, iter_max, Nx_loc;
    double Lx, Ly, D, dt, dx, dy, t, eps, ecart_loc, ecart;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    ifstream file("valeurs.txt");
    file >> Nx >> Ny >> Lx >> Ly >> D >> dt >> r >> cas ;
    file.close();

    double (*f_ex)(double, double);
    if (cas == 1){f_ex = f1_ex;}
    else if (cas == 2) {f_ex = f2_ex;}

    charge(rank, Nx, nproc, &iBeg, &iEnd);


    vector<double> U0, U, b;
    // if ((rank == 0) || (rank == nproc-1))
    // {
    //     N_loc = (iEnd-iBeg+1+r)*Ny;
    //     U0.resize(N_loc); U.resize(N_loc); b.resize(N_loc);
    //     Nx_loc = N_loc/Ny;
    //     printf("N_loc = %d\n", N_loc);
    //     printf("Nx_loc = %d\n", Nx_loc);
    // }
    // else 
    // {   
    //     N_loc = (iEnd-iBeg+1+2*r)*Ny;
    //     U0.resize(N_loc); U.resize(N_loc); b.resize(N_loc);
    //     Nx_loc = N_loc/Ny;
    //     printf("N_loc = %d\n", N_loc);
    //     printf("Nx_loc = %d\n", Nx_loc);
    // }

    if ((rank == 0)) 
    {
        N_loc = (iEnd-iBeg+1)*Ny;
        U0.resize(N_loc); U.resize(N_loc); b.resize(N_loc);
        Nx_loc = N_loc/Ny;
        printf("N_loc = %d\n", N_loc);
        printf("Nx_loc = %d\n", Nx_loc);
    }
    else 
    {
        N_loc = (iEnd-iBeg+1+r)*Ny;
        U0.resize(N_loc); U.resize(N_loc); b.resize(N_loc);
        Nx_loc = N_loc/Ny;
        printf("N_loc = %d\n", N_loc);
        printf("Nx_loc = %d\n", Nx_loc);
    }

    printf("Proc %d: iBeg = %d, iEnd = %d, N_loc = %d\n", rank, iBeg, iEnd, N_loc);

    // vector<double> sol(Nx*Ny), u_Global(Nx*Ny);

    N = Nx*Ny;
    dx = Lx/(Nx+1);
    dy = Ly/(Ny+1);
    n_iter = 10;    
    t = 0;  
  
    iter_max = 100;
    eps = 0.000001;

    for (int i=0; i<N_loc; i++){U0[i]=0;}

    vector<double> U_g(Ny), U_d(Ny), U_g_0(Ny), U_d_0(Ny);
    for (int i=0; i<Ny; i++){U_g_0[i]=1; U_d_0[i]=1;}

    for (int i=0; i<n_iter; i++)
    {   
        for (int j=0; j<iter_max; j++)
        {    
            int left = rank-1;
            int right = rank+1;

            if (rank==0){left = MPI_PROC_NULL;}
            if (rank==nproc-1){right = MPI_PROC_NULL;}

            // Envoie à droite, reçoit de gauche
            MPI_Sendrecv(&U0[N_loc-2*r*Ny], Ny, MPI_DOUBLE, right, 0,
            &U_g[0], Ny, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // printf("%d\n", N_loc-2*r*Ny);
            
            // Envoie à gauche, reçoit de droite
            MPI_Sendrecv(&U0[2*(r-1)*Ny], Ny, MPI_DOUBLE, left,  1, 
            &U_d[0], Ny, MPI_DOUBLE, right, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // printf("%d\n", (2*r-1)*Ny);
            
            for (int k=0; k<Ny; k++){U_g_0[k] -= U_g[k]; U_d_0[k] -= U_d[k];}
            ecart_loc = (norme2(U_g_0) + norme2(U_d_0))/2;

            MPI_Allreduce(&ecart_loc, &ecart, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            // if (ecart < eps){printf("n_iter = %d\n", j); break;}

            b = make_b(U0, Nx_loc, Ny, dx, dy, Lx, Ly, D, dt, t, cas, rank, nproc, iBeg, U_g, U_d, r);
            U = gradient_conjugue(Nx_loc, Ny, D, dx, dy, dt, b, 0.0001, 1000);

            U0=U;
            U_g_0 = U_g;
            U_d_0 = U_d;

        }

        t += dt;
    }
 
    ofstream fichier("affichage.dat", ios::app);
    int i_off;
    if (rank==0){i_off = 0;}
    else {i_off = (iBeg - r);}
    for (int i=0; i<N_loc; i++)
    {
        fichier << (i_off + i/Ny+1)*dx << " " << (i%Ny+1)*dy << " " << U[i] << endl;
    }
    fichier.close();

    MPI_Finalize();

    return 0;
}
