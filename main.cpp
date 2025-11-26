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
    int Nx, Ny, n_iter, N, cas, r, nproc, iBeg, iEnd, rank, N_loc;
    double Lx, Ly, D, dt, dx, dy, t;

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
    if ((iBeg == 0) && (iEnd == Nx-1))
    {
        N_loc = (iEnd-iBeg+1+(r-1))*Ny;
        U0.resize(N_loc); U.resize(N_loc); b.resize(N_loc);
    }
    else 
    {   
        N_loc = (iEnd-iBeg+1+2*(r-1))*Ny;
        U0.resize(N_loc); U.resize(N_loc); b.resize(N_loc);
    }

    printf("Proc %d: iBeg = %d, iEnd = %d, N_loc = %d\n", rank, iBeg, iEnd, N_loc);

    vector<double> sol(Nx*Ny), u_Global(Nx*Ny);

    N = Nx*Ny;
    dx = Lx/(Nx+1);
    dy = Ly/(Ny+1);
    n_iter = 10;    
    t = 0;  
  

    n_iter = 1;
    t = 0;


    for (int i=0; i<N_loc; i++){U0[i]=0;
    }

    vector<double> U_g(Ny), U_d(Ny);

    for (int i=0; i<n_iter; i++)
    {      
        int left = rank-1;
        int right = rank+1;
        printf("Proc %d: itération %d\n", rank, i);

        if (rank==0){left = MPI_PROC_NULL;}
        if (rank==nproc-1){right = MPI_PROC_NULL;}

        // Envoie à droite, reçoit de gauche
        MPI_Sendrecv(&U0[N_loc-r*Ny], Ny, MPI_DOUBLE, right, 0,
        &U_g[0], Ny, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        
        
        // Envoie à gauche, reçoit de droite
        MPI_Sendrecv(&U0[(r-1)*Ny], Ny, MPI_DOUBLE, left,  1, 
        &U_d[0], Ny, MPI_DOUBLE, right, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        t += dt;
        b = make_b(U0, N_loc/Ny, Ny, dx, dy, Lx, Ly, D, dt, t, cas, rank, nproc, U_g, U_d);
        for (int i = 0; i < N_loc; i++)
        {
            printf("Proc %d: b[%d] = %f\n", rank, i, b[i]);
        }
        
        U = gradient_conjugue(N_loc/Ny, Ny, D, dx, dy, dt, b, 0.0001, 1000);
        for (int i = 0; i < N_loc; i++)
        {
            printf("Proc %d: U[%d] = %f\n", rank, i, U[i]);
        }
        U0=U;
    }
    vector<double> U_final(Nx*Ny);
    // Rassembler les résultats de tous les processus
    MPI_Allgather(U0.data() + ((rank == 0) ? 0 : (r-1)*Ny), 
                  (iEnd - iBeg + 1)*Ny, MPI_DOUBLE, 
                  U_final.data(), 
                  (iEnd - iBeg + 1)*Ny, MPI_DOUBLE, 
                  MPI_COMM_WORLD);
                
    
    
    ofstream fichier("affichage.dat");
    for (int i=0; i<N; i++)
    {
        fichier << (i/Ny+1)*dx << " " << (i%Ny+1)*dy << " " << U_final[i] << endl;
    }
    fichier.close();

    MPI_Finalize();

    return 0;
}
