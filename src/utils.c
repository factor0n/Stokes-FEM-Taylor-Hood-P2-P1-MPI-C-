#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void die(const char* msg) {
    fprintf(stderr, "%s\n", msg);
    exit(EXIT_FAILURE);
}

double dot_mpi(const double* a, const double* b, int n, MPI_Comm comm) {
    double local = 0.0;
    for (int i=0;i<n;++i) local += a[i]*b[i];
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, comm);
    return global;
}

double norm2_mpi(const double* a, int n, MPI_Comm comm) {
    return sqrt(dot_mpi(a,a,n,comm));
}

void axpy(double* y, const double* x, double alpha, int n) {
    for (int i=0;i<n;++i) y[i] += alpha*x[i];
}

void scal(double* x, double alpha, int n) {
    for (int i=0;i<n;++i) x[i] *= alpha;
}

void copy_vec(double* dst, const double* src, int n) {
    memcpy(dst, src, (size_t)n*sizeof(double));
}

void zero_vec(double* x, int n) {
    memset(x, 0, (size_t)n*sizeof(double));
}
