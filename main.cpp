#include "fonction.h"
#include "produit.h"
#include "grad_conj.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include "charge.h"
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[])
{   

    int Nx, Ny, n_iter, N, cas, r, nproc, iBeg, iEnd, rank;
    double Lx, Ly, D, dt, dx, dy, t, err;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        printf("Choissisez le cas: 1, 2 ou 3\n");
        cin >> cas;

        printf("Choisissez le recouvrement\n");
        cin >> r;
    }

    MPI_Bcast(&cas, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&r, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double (*f_ex)(double, double);
    if (cas == 1){f_ex = f1_ex;}
    else if (cas == 2) {f_ex = f2_ex;}
  
    ifstream file("valeurs.txt");
    file >> Nx >> Ny >> Lx >> Ly >> D >> dt ;
    file.close();

    for(int i = 0; i<nproc; i++){charge(rank, Nx*Ny, nproc, &iBeg, &iEnd);}

    vector<double> U0, U, b;
    if ((iBeg == 0) && (iEnd = Nx*Ny-1)){U0.resize(iEnd-iBeg + 1 + (r-1)); U.resize(iEnd-iBeg + 1 + (r-1)); b.resize(iEnd-iBeg + 1 + (r-1));}
    else {U0.resize(iEnd-iBeg + 1 + 2*(r-1)); U.resize(iEnd-iBeg + 1 + 2*(r-1)); b.resize(iEnd-iBeg + 1 + 2*(r-1));}

    vector<double> sol(Nx*Ny), u_Global(Nx*Ny);

    N = Nx*Ny;
    dx = Lx/(Nx+1);
    dy = Ly/(Ny+1);
    n_iter = 10;    
    t = 0;  
    // --------------------------------------
    // Fichier pour écire l'erreur
    // ofstream fichier("convergence.dat");
    // --------------------------------------


    for (int j=1; j<=1; j++)
    {   
        // --------------------------------------
        // Paramètre pour le calcul d'erreur   

        // Nx = pow(2, j);
        // Ny = Nx;
        // N = Nx*Ny;
        // dx = Lx/(Nx+1);
        // dy = Ly/(Ny+1);
        // --------------------------------------

        n_iter = 50;
        t = 0;


        if ((iBeg == 0) && (iEnd = Nx*Ny-1)){for (int i=0; i<iEnd-iBeg + 1 + (r-1); i++){U0[i]=0;}}
        else {for (int i=0; i<iEnd-iBeg + 1 + 2*(r-1); i++){U0[i]=0;}}

        vector<double> U_g(Ny), U_d(Ny);

        for (int i=0; i<n_iter; i++)
        {      
            int left = rank-1;
            int right = rank+1;

            if (rank==0){left = MPI_PROC_NULL;}
            if (rank==nproc-1){right = MPI_PROC_NULL;}

            // Envoie à droite, reçoit de gauche
            MPI_Sendrecv(&U0[N-r*Ny], Ny, MPI_DOUBLE, right, 0,
            &U_g[0], Ny, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            // Envoie à gauche, reçoit de droite
            MPI_Sendrecv(&U0[(r-1)*Ny], Ny, MPI_DOUBLE, left,  1, 
            &U_d[0], Ny, MPI_DOUBLE, right, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            t += dt;
            b = make_b(U0, Nx, Ny, dx, dy, Lx, Ly, D, dt, t, cas, rank, nproc, U_g, U_d);
            U = gradient_conjugue(Nx, Ny, D, dx, dy, dt, b, 0.0001, 1000);

            U0=U;
        }

        ofstream fichier("affichage.dat");
        for (int i=0; i<N; i++)
        {
            fichier << (i/Ny+1)*dx << " " << (i%Ny+1)*dy << " " << U[i] << endl;
        }
        fichier.close();

        // --------------------------------------
        // Affichage de l'erreur dans le cas 2
        // if ((cas == 2))
        // {
        // err = 0;
        // for (int i=0; i<N; i++)
        // {
        //     sol[i] = f_ex((i/Ny+1)*dx, (i%Ny+1)*dy);
        //     err += (sol[i]-U[i])*(sol[i]-U[i])*dx*dy;
        // }   
        // err = sqrt(err);
        // }
        // fichier << dx << " " << err << endl;
        // printf("dx = %f, err = %f\n", dx, err);
        // --------------------------------------


    }
    // --------------------------------------
    // fermer fichier erreur
    // fichier.close();
    // --------------------------------------


    MPI_Finalize();

    return 0;
}
