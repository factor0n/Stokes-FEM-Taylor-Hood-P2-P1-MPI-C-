#include "gmres.h"
#include "utils.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

int gmres_restart_mpi(
    const CSR* A, const double* b, double* x,
    int restart_m, int maxit, double tol,
    MPI_Comm comm
){
    int n=A->n;

    double* r=(double*)malloc((size_t)n*sizeof(double));
    double* Ax=(double*)malloc((size_t)n*sizeof(double));
    double* w=(double*)malloc((size_t)n*sizeof(double));
    if(!r||!Ax||!w) die("malloc gmres");

    double* V=(double*)malloc((size_t)(restart_m+1)*(size_t)n*sizeof(double));
    double* H=(double*)malloc((size_t)(restart_m+1)*(size_t)restart_m*sizeof(double));
    double* cs=(double*)malloc((size_t)restart_m*sizeof(double));
    double* sn=(double*)malloc((size_t)restart_m*sizeof(double));
    double* g =(double*)malloc((size_t)(restart_m+1)*sizeof(double));
    double* y =(double*)malloc((size_t)restart_m*sizeof(double));
    if(!V||!H||!cs||!sn||!g||!y) die("malloc gmres work");

    int it_total=0;

    while(it_total<maxit){
        csr_matvec(A,x,Ax);
        for(int i=0;i<n;++i) r[i]=b[i]-Ax[i];

        double beta = norm2_mpi(r,n,comm);
        if(beta<tol) break;

        for(int i=0;i<n;++i) V[i]=r[i]/beta;

        memset(H,0,(size_t)(restart_m+1)*(size_t)restart_m*sizeof(double));
        for(int i=0;i<restart_m+1;++i) g[i]=0.0;
        g[0]=beta;

        int m_used=0;
        for(int j=0;j<restart_m && it_total<maxit; ++j, ++it_total){
            m_used=j+1;

            csr_matvec(A, &V[(size_t)j*(size_t)n], w);

            for(int i=0;i<=j;++i){
                double hij = dot_mpi(&V[(size_t)i*(size_t)n], w, n, comm);
                H[(size_t)i*(size_t)restart_m + (size_t)j] = hij;
                axpy(w, &V[(size_t)i*(size_t)n], -hij, n);
            }
            double hnext = norm2_mpi(w,n,comm);
            H[(size_t)(j+1)*(size_t)restart_m + (size_t)j] = hnext;

            if(hnext>1e-14){
                for(int i=0;i<n;++i) V[(size_t)(j+1)*(size_t)n + (size_t)i]=w[i]/hnext;
            } else {
                for(int i=0;i<n;++i) V[(size_t)(j+1)*(size_t)n + (size_t)i]=0.0;
            }

            for(int i=0;i<j;++i){
                double h1=H[(size_t)i*(size_t)restart_m + (size_t)j];
                double h2=H[(size_t)(i+1)*(size_t)restart_m + (size_t)j];
                H[(size_t)i*(size_t)restart_m + (size_t)j]     =  cs[i]*h1 + sn[i]*h2;
                H[(size_t)(i+1)*(size_t)restart_m + (size_t)j] = -sn[i]*h1 + cs[i]*h2;
            }

            double h_jj=H[(size_t)j*(size_t)restart_m + (size_t)j];
            double h_j1j=H[(size_t)(j+1)*(size_t)restart_m + (size_t)j];
            double denom=sqrt(h_jj*h_jj + h_j1j*h_j1j);
            if(denom<1e-30){ cs[j]=1.0; sn[j]=0.0; }
            else { cs[j]=h_jj/denom; sn[j]=h_j1j/denom; }

            H[(size_t)j*(size_t)restart_m + (size_t)j] = cs[j]*h_jj + sn[j]*h_j1j;
            H[(size_t)(j+1)*(size_t)restart_m + (size_t)j] = 0.0;

            double g_j=g[j], g_j1=g[j+1];
            g[j]   = cs[j]*g_j + sn[j]*g_j1;
            g[j+1] = -sn[j]*g_j + cs[j]*g_j1;

            if(fabs(g[j+1]) < tol) break;
        }

        for(int i=0;i<m_used;++i) y[i]=0.0;
        for(int i=m_used-1;i>=0;--i){
            double s=g[i];
            for(int k=i+1;k<m_used;++k) s -= H[(size_t)i*(size_t)restart_m + (size_t)k]*y[k];
            double hii=H[(size_t)i*(size_t)restart_m + (size_t)i];
            if(fabs(hii)<1e-30) hii=1e-30;
            y[i]=s/hii;
        }

        for(int j=0;j<m_used;++j){
            axpy(x, &V[(size_t)j*(size_t)n], y[j], n);
        }
    }

    free(r); free(Ax); free(w);
    free(V); free(H); free(cs); free(sn); free(g); free(y);
    return it_total;
}
