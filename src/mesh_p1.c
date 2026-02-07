#include "mesh_p1.h"
#include "utils.h"
#include <stdlib.h>

int vid(int i, int j, int N) { return j*(N+1)+i; }

MeshP1 make_mesh_p1(int N) {
    MeshP1 m;
    m.N = N;
    m.nV = (N+1)*(N+1);
    m.nElem = 2*N*N;

    m.xV = (double*)malloc((size_t)m.nV*sizeof(double));
    m.yV = (double*)malloc((size_t)m.nV*sizeof(double));
    m.triV = (int*)malloc((size_t)3*(size_t)m.nElem*sizeof(int));
    if (!m.xV || !m.yV || !m.triV) die("malloc make_mesh_p1 failed");

    for (int j=0;j<=N;++j) {
        for (int i=0;i<=N;++i) {
            int id = vid(i,j,N);
            m.xV[id] = (double)i/(double)N;
            m.yV[id] = (double)j/(double)N;
        }
    }

    int e=0;
    for (int j=0;j<N;++j) {
        for (int i=0;i<N;++i) {
            int v0 = vid(i,  j,  N);
            int v1 = vid(i+1,j,  N);
            int v2 = vid(i,  j+1,N);
            int v3 = vid(i+1,j+1,N);

            m.triV[3*e+0]=v0; m.triV[3*e+1]=v1; m.triV[3*e+2]=v2; e++;
            m.triV[3*e+0]=v1; m.triV[3*e+1]=v3; m.triV[3*e+2]=v2; e++;
        }
    }
    return m;
}

void free_mesh_p1(MeshP1* m) {
    free(m->xV); free(m->yV); free(m->triV);
    m->xV=NULL; m->yV=NULL; m->triV=NULL;
    m->nV=0; m->nElem=0; m->N=0;
}
