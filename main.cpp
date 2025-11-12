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

int main(int argc, char *argv[]){   
    int Nx, Ny, n_iter, N, cas,r, nproc, iBeg, iEnd, rank;
    double Lx, Ly, D, dt, dx, dy, t, err;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank==0)
    {
            printf("Choissisez le cas: 1, 2 ou 3\n");
            cin >> cas;

    }
    if (rank==0)
    {
            printf("Choississez le nombre de recouvrement\n");
            cin >> r;

    }
    MPI_Bcast(&cas, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double (*f_ex)(double, double);
    if (cas == 1){f_ex = f1_ex;}
    else if (cas == 2) {f_ex = f2_ex;}
  

    ifstream file("valeurs.txt");
    file >> Nx >> Ny >> Lx >> Ly >> D >> dt ;
    file.close();

    for(int i=0; i<nproc; i++){charge(rank, Nx*Ny, nproc, &iBeg, &iEnd);}
    printf("Proc %d: iBeg = %d, iEnd = %d\n", rank, iBeg, iEnd);

    vector<double>  sol(Nx*Ny), U_global(Nx*Ny);
    if (iBeg==0)
    {
        vector<double> U0(iEnd-iBeg+1+r), U(iEnd-iBeg+1+r), b(iEnd-iBeg+1+r);
    }
    else if (iEnd==Nx*Ny-1)
    {
        vector<double> U0(iEnd-iBeg+1+r), U(iEnd-iBeg+1+r), b(iEnd-iBeg+1+r) ;  

    }
    else
    {
        vector<double> U0(iEnd-iBeg+1+2*r), U(iEnd-iBeg+1+2*r), b(iEnd-iBeg+1+2*r) ; 
    }
    
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

        vector<double> U0(N), U(N), b(N), sol(N);

        for (int i=0; i<N; i++){U0[i]=0;}

        for (int i=0; i<n_iter; i++)
        {   
            t += dt;
            b = make_b(U0, Nx, Ny, dx, dy, Lx, Ly, D, dt, t, cas);
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
