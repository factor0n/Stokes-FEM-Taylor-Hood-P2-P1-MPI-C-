#include "csr.h"
#include <stdlib.h>

int csr_find(const CSR* A, int row, int col){
    int lo=A->rowptr[row], hi=A->rowptr[row+1]-1;
    while(lo<=hi){
        int mid=(lo+hi)/2;
        int c=A->colind[mid];
        if(c==col) return mid;
        if(c<col) lo=mid+1; else hi=mid-1;
    }
    return -1;
}

void csr_matvec(const CSR* A, const double* x, double* y){
    for(int i=0;i<A->n;++i){
        double s=0.0;
        for(int k=A->rowptr[i]; k<A->rowptr[i+1]; ++k)
            s += A->val[k]*x[A->colind[k]];
        y[i]=s;
    }
}

void free_csr(CSR* A){
    free(A->rowptr); free(A->colind); free(A->val);
    A->rowptr=NULL; A->colind=NULL; A->val=NULL;
    A->n=0; A->nnz=0;
}
