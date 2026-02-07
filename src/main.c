#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "mesh_p1.h"
#include "p2.h"
#include "fem_stokes.h"
#include "gmres.h"
#include "csr.h"

int main(int argc, char** argv){
    MPI_Init(&argc,&argv);
    MPI_Comm comm=MPI_COMM_WORLD;

    int rank=0,size=1;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&size);

    int N=16;
    if(argc>=2) N=atoi(argv[1]);
    if(N<4) N=4;

    double nu=1.0;

    if(rank==0){
        printf("Stokes Taylor–Hood P2/P1 + MPI\n");
        printf("N=%d | ranks=%d | nu=%.3g\n", N, size, nu);
    }

    MeshP1 m = make_mesh_p1(N);
    P2Data p2 = build_p2_from_p1(&m);

    int nP2=p2.nP2;
    int nV=m.nV;
    int nTot=2*nP2 + nV;

    if(rank==0){
        printf("Vertices (P1 pressure): %d\n", nV);
        printf("Edges (P2 mid nodes): %d\n", p2.nE);
        printf("Scalar P2 nodes: %d\n", nP2);
        printf("Total unknowns (Ux,Uy,P): %d\n", nTot);
        printf("Triangles: %d\n", m.nElem);
    }

    CSR A = build_stokes_pattern(&m,&p2);

    double* rhs_local=(double*)calloc((size_t)nTot,sizeof(double));
    double* val_local=(double*)calloc((size_t)A.nnz,sizeof(double));
    if(!rhs_local||!val_local) die("calloc local failed");

    CSR Aloc=A; Aloc.val=val_local;
    assemble_stokes_local(&m,&p2,&Aloc,rhs_local,rank,size,nu);

    double* rhs=(double*)malloc((size_t)nTot*sizeof(double));
    if(!rhs) die("malloc rhs failed");

    MPI_Reduce(val_local, A.val, A.nnz, MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(rhs_local, rhs, nTot, MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Bcast(A.val, A.nnz, MPI_DOUBLE, 0, comm);
    MPI_Bcast(rhs, nTot, MPI_DOUBLE, 0, comm);

    apply_dirichlet_velocity(&A,rhs,&m,&p2);
    fix_pressure_dof(&A,rhs,0,nP2);

    double* x=(double*)calloc((size_t)nTot,sizeof(double));
    if(!x) die("calloc x failed");

    double t0=MPI_Wtime();
    int it = gmres_restart_mpi(&A,rhs,x,40,2000,1e-8,comm);
    double t1=MPI_Wtime();

    if(rank==0){
        printf("GMRES iterations: %d\n", it);
        printf("Solve time: %.6f s\n", t1-t0);

        int ic=N/2, jc=N/2;
        int vcenter=vid(ic,jc,N);
        printf("Sample @ center vertex: ux=%.6e uy=%.6e p=%.6e\n",
               x[vcenter], x[nP2+vcenter], x[2*nP2+vcenter]);
    }

    free(x);
    free(rhs);
    free(rhs_local);
    free(val_local);
    free_csr(&A);
    free_p2(&p2);
    free_mesh_p1(&m);

    MPI_Finalize();
    return 0;
}
