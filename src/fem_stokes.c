#include "fem_stokes.h"
#include "utils.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ---- small intlist for pattern build ---- */
typedef struct { int* data; int size; int cap; } IntList;
static void il_init(IntList* L){ L->data=NULL; L->size=0; L->cap=0; }
static void il_free(IntList* L){ free(L->data); L->data=NULL; L->size=L->cap=0; }
static void il_push_unique(IntList* L, int v){
    for(int i=0;i<L->size;++i) if(L->data[i]==v) return;
    if(L->size==L->cap){
        int nc=(L->cap==0)?8:(2*L->cap);
        int* nd=(int*)realloc(L->data,(size_t)nc*sizeof(int));
        if(!nd) die("realloc ilist failed");
        L->data=nd; L->cap=nc;
    }
    L->data[L->size++]=v;
}
static int int_cmp(const void* a, const void* b){
    int ia=*(const int*)a, ib=*(const int*)b;
    return (ia>ib)-(ia<ib);
}

/* ---- quadrature ---- */
static void tri_quad_3pt(double xi[3], double eta[3], double w[3]){
    xi[0]=1.0/6.0; eta[0]=1.0/6.0;
    xi[1]=2.0/3.0; eta[1]=1.0/6.0;
    xi[2]=1.0/6.0; eta[2]=2.0/3.0;
    w[0]=w[1]=w[2]=1.0/6.0;
}

/* ---- P1 basis for pressure ---- */
static void p1_psi(double xi, double eta, double psi[3]){
    double l1=1.0-xi-eta, l2=xi, l3=eta;
    psi[0]=l1; psi[1]=l2; psi[2]=l3;
}

/* ---- P2 basis for velocity (scalar) ---- */
static void p2_phi_and_grad_ref(double xi, double eta, double phi[6], double gref[6][2]){
    double l1=1.0-xi-eta, l2=xi, l3=eta;

    phi[0]=l1*(2*l1-1);
    phi[1]=l2*(2*l2-1);
    phi[2]=l3*(2*l3-1);
    phi[3]=4*l1*l2;
    phi[4]=4*l2*l3;
    phi[5]=4*l3*l1;

    double gl1[2]={-1.0,-1.0};
    double gl2[2]={ 1.0, 0.0};
    double gl3[2]={ 0.0, 1.0};

    double c1=4*l1-1, c2=4*l2-1, c3=4*l3-1;

    gref[0][0]=c1*gl1[0]; gref[0][1]=c1*gl1[1];
    gref[1][0]=c2*gl2[0]; gref[1][1]=c2*gl2[1];
    gref[2][0]=c3*gl3[0]; gref[2][1]=c3*gl3[1];

    gref[3][0]=4*(l2*gl1[0] + l1*gl2[0]);
    gref[3][1]=4*(l2*gl1[1] + l1*gl2[1]);

    gref[4][0]=4*(l3*gl2[0] + l2*gl3[0]);
    gref[4][1]=4*(l3*gl2[1] + l2*gl3[1]);

    gref[5][0]=4*(l1*gl3[0] + l3*gl1[0]);
    gref[5][1]=4*(l1*gl3[1] + l3*gl1[1]);
}

/* ---- MPI ownership ---- */
void elem_owner_centroid_x(const MeshP1* m, int e, int size, int* owner){
    int v0=m->triV[3*e+0], v1=m->triV[3*e+1], v2=m->triV[3*e+2];
    double cx=(m->xV[v0]+m->xV[v1]+m->xV[v2])/3.0;
    int r=(int)floor(cx*(double)size);
    if(r<0) r=0;
    if(r>=size) r=size-1;
    *owner=r;
}

/* ---- build CSR pattern for global saddle-point system ---- */
CSR build_stokes_pattern(const MeshP1* m, const P2Data* p2){
    int nP2=p2->nP2;
    int nV=m->nV;
    int nTot=2*nP2 + nV;

    #define IDX_UX(i) ( (i) )
    #define IDX_UY(i) ( (nP2) + (i) )
    #define IDX_P(j)  ( (2*nP2) + (j) )

    IntList* rows=(IntList*)malloc((size_t)nTot*sizeof(IntList));
    if(!rows) die("malloc rows failed");
    for(int i=0;i<nTot;++i) il_init(&rows[i]);

    for(int e=0;e<m->nElem;++e){
        int v0=m->triV[3*e+0], v1=m->triV[3*e+1], v2=m->triV[3*e+2];

        int e01=p2->triE[3*e+0];
        int e12=p2->triE[3*e+1];
        int e20=p2->triE[3*e+2];

        int s[6]={ v0, v1, v2, m->nV+e01, m->nV+e12, m->nV+e20 };
        int pj[3]={ v0, v1, v2 };

        for(int a=0;a<6;++a){
            for(int b=0;b<6;++b){
                il_push_unique(&rows[IDX_UX(s[a])], IDX_UX(s[b]));
                il_push_unique(&rows[IDX_UY(s[a])], IDX_UY(s[b]));
            }
        }
        for(int a=0;a<6;++a){
            for(int q=0;q<3;++q){
                il_push_unique(&rows[IDX_UX(s[a])], IDX_P(pj[q]));
                il_push_unique(&rows[IDX_UY(s[a])], IDX_P(pj[q]));
            }
        }
        for(int q=0;q<3;++q){
            for(int a=0;a<6;++a){
                il_push_unique(&rows[IDX_P(pj[q])], IDX_UX(s[a]));
                il_push_unique(&rows[IDX_P(pj[q])], IDX_UY(s[a]));
            }
        }
    }

    int nnz=0;
    for(int i=0;i<nTot;++i){
        qsort(rows[i].data,(size_t)rows[i].size,sizeof(int),int_cmp);
        nnz += rows[i].size;
    }

    CSR A;
    A.n=nTot; A.nnz=nnz;
    A.rowptr=(int*)malloc((size_t)(nTot+1)*sizeof(int));
    A.colind=(int*)malloc((size_t)nnz*sizeof(int));
    A.val=(double*)malloc((size_t)nnz*sizeof(double));
    if(!A.rowptr||!A.colind||!A.val) die("malloc CSR failed");

    int p=0;
    A.rowptr[0]=0;
    for(int i=0;i<nTot;++i){
        for(int k=0;k<rows[i].size;++k) A.colind[p++]=rows[i].data[k];
        A.rowptr[i+1]=p;
    }
    for(int k=0;k<nnz;++k) A.val[k]=0.0;

    for(int i=0;i<nTot;++i) il_free(&rows[i]);
    free(rows);

    #undef IDX_UX
    #undef IDX_UY
    #undef IDX_P
    return A;
}

/* ---- assembly ---- */
void assemble_stokes_local(
    const MeshP1* m, const P2Data* p2,
    CSR* A, double* rhs,
    int rank, int size, double nu
){
    int nP2=p2->nP2;

    #define IDX_UX(i) ( (i) )
    #define IDX_UY(i) ( (nP2) + (i) )
    #define IDX_P(j)  ( (2*nP2) + (j) )

    double qxi[3], qeta[3], qw[3];
    tri_quad_3pt(qxi,qeta,qw);

    // body force f=(1,1) for demo
    const double fx=1.0, fy=1.0;

    for(int e=0;e<m->nElem;++e){
        int owner=0; elem_owner_centroid_x(m,e,size,&owner);
        if(owner!=rank) continue;

        int v0=m->triV[3*e+0], v1=m->triV[3*e+1], v2=m->triV[3*e+2];

        int e01=p2->triE[3*e+0];
        int e12=p2->triE[3*e+1];
        int e20=p2->triE[3*e+2];

        int s[6]={ v0, v1, v2, m->nV+e01, m->nV+e12, m->nV+e20 };
        int pj[3]={ v0, v1, v2 };

        double x1=m->xV[v0], y1=m->yV[v0];
        double x2=m->xV[v1], y2=m->yV[v1];
        double x3=m->xV[v2], y3=m->yV[v2];

        double J00=x2-x1, J01=x3-x1;
        double J10=y2-y1, J11=y3-y1;
        double detJ = J00*J11 - J01*J10;
        double absDetJ=fabs(detJ);
        if(absDetJ<1e-30) die("degenerate element");

        double invJT00 =  J11/detJ;
        double invJT01 = -J01/detJ;
        double invJT10 = -J10/detJ;
        double invJT11 =  J00/detJ;

        double Kloc[6][6]; memset(Kloc,0,sizeof(Kloc));
        double Bxloc[3][6]; memset(Bxloc,0,sizeof(Bxloc));
        double Byloc[3][6]; memset(Byloc,0,sizeof(Byloc));
        double Flocx[6]; memset(Flocx,0,sizeof(Flocx));
        double Flocy[6]; memset(Flocy,0,sizeof(Flocy));

        for(int qp=0; qp<3; ++qp){
            double xi=qxi[qp], eta=qeta[qp], w=qw[qp];

            double phi[6], gref[6][2], psi[3];
            p2_phi_and_grad_ref(xi,eta,phi,gref);
            p1_psi(xi,eta,psi);

            double gphys[6][2];
            for(int a=0;a<6;++a){
                double gx = invJT00*gref[a][0] + invJT01*gref[a][1];
                double gy = invJT10*gref[a][0] + invJT11*gref[a][1];
                gphys[a][0]=gx; gphys[a][1]=gy;
            }

            double dV = w * absDetJ;

            for(int a=0;a<6;++a){
                Flocx[a] += dV * fx * phi[a];
                Flocy[a] += dV * fy * phi[a];
                for(int b=0;b<6;++b){
                    Kloc[a][b] += dV * (gphys[a][0]*gphys[b][0] + gphys[a][1]*gphys[b][1]);
                }
            }

            for(int j=0;j<3;++j){
                for(int a=0;a<6;++a){
                    Bxloc[j][a] += dV * psi[j] * gphys[a][0];
                    Byloc[j][a] += dV * psi[j] * gphys[a][1];
                }
            }
        }

        // assemble
        for(int a=0;a<6;++a){
            int rowUx=IDX_UX(s[a]);
            int rowUy=IDX_UY(s[a]);

            rhs[rowUx] += Flocx[a];
            rhs[rowUy] += Flocy[a];

            for(int b=0;b<6;++b){
                int colUx=IDX_UX(s[b]);
                int colUy=IDX_UY(s[b]);

                int posUx = csr_find(A,rowUx,colUx);
                int posUy = csr_find(A,rowUy,colUy);
                if(posUx<0||posUy<0) die("CSR pos K not found");

                A->val[posUx] += nu*Kloc[a][b];
                A->val[posUy] += nu*Kloc[a][b];
            }

            for(int j=0;j<3;++j){
                int colP=IDX_P(pj[j]);
                int posUxP=csr_find(A,rowUx,colP);
                int posUyP=csr_find(A,rowUy,colP);
                if(posUxP<0||posUyP<0) die("CSR pos BT not found");
                A->val[posUxP] += -Bxloc[j][a];
                A->val[posUyP] += -Byloc[j][a];
            }
        }

        for(int j=0;j<3;++j){
            int rowP=IDX_P(pj[j]);
            for(int a=0;a<6;++a){
                int colUx=IDX_UX(s[a]);
                int colUy=IDX_UY(s[a]);
                int posPUx=csr_find(A,rowP,colUx);
                int posPUy=csr_find(A,rowP,colUy);
                if(posPUx<0||posPUy<0) die("CSR pos B not found");
                A->val[posPUx] += Bxloc[j][a];
                A->val[posPUy] += Byloc[j][a];
            }
        }
    }

    #undef IDX_UX
    #undef IDX_UY
    #undef IDX_P
}

/* ---- BCs ---- */
void apply_dirichlet_velocity(CSR* A, double* rhs, const MeshP1* m, const P2Data* p2){
    int nP2=p2->nP2;

    for(int i=0;i<nP2;++i){
        if(is_boundary_xy(p2->xP2[i], p2->yP2[i])){
            int rowUx=i;
            int rowUy=nP2+i;

            for(int k=A->rowptr[rowUx]; k<A->rowptr[rowUx+1]; ++k) A->val[k]=0.0;
            int pos=csr_find(A,rowUx,rowUx);
            if(pos<0) die("diag Ux not found");
            A->val[pos]=1.0; rhs[rowUx]=0.0;

            for(int k=A->rowptr[rowUy]; k<A->rowptr[rowUy+1]; ++k) A->val[k]=0.0;
            pos=csr_find(A,rowUy,rowUy);
            if(pos<0) die("diag Uy not found");
            A->val[pos]=1.0; rhs[rowUy]=0.0;
        }
    }
    (void)m;
}

void fix_pressure_dof(CSR* A, double* rhs, int pressure_vertex_id, int nP2){
    int rowP=2*nP2 + pressure_vertex_id;
    for(int k=A->rowptr[rowP]; k<A->rowptr[rowP+1]; ++k) A->val[k]=0.0;
    int pos=csr_find(A,rowP,rowP);
    if(pos<0) die("diag P not found");
    A->val[pos]=1.0;
    rhs[rowP]=0.0;
}
