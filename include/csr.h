#ifndef CSR_H
#define CSR_H

typedef struct {
    int n;
    int nnz;
    int* rowptr;
    int* colind;
    double* val;
} CSR;

void csr_matvec(const CSR* A, const double* x, double* y);
int  csr_find(const CSR* A, int row, int col);
void free_csr(CSR* A);

#endif
