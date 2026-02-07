#ifndef GMRES_H
#define GMRES_H

#include "csr.h"
#include <mpi.h>

int gmres_restart_mpi(
    const CSR* A, const double* b, double* x,
    int restart_m, int maxit, double tol,
    MPI_Comm comm
);

#endif
