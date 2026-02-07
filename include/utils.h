#ifndef UTILS_H
#define UTILS_H

#include <mpi.h>

void die(const char* msg);

double dot_mpi(const double* a, const double* b, int n, MPI_Comm comm);
double norm2_mpi(const double* a, int n, MPI_Comm comm);

void axpy(double* y, const double* x, double alpha, int n);
void scal(double* x, double alpha, int n);
void copy_vec(double* dst, const double* src, int n);
void zero_vec(double* x, int n);

#endif
