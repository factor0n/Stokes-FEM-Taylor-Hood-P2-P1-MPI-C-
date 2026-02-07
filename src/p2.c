#include "p2.h"
#include "utils.h"
#include <stdlib.h>
#include <math.h>

static int edge_lower(int a,int b){ return a<b?a:b; }
static int edge_upper(int a,int b){ return a<b?b:a; }

static int edge_cmp(const void* A, const void* B) {
    const Edge* e1=(const Edge*)A;
    const Edge* e2=(const Edge*)B;
    if (e1->a != e2->a) return (e1->a>e2->a)-(e1->a<e2->a);
    return (e1->b>e2->b)-(e1->b<e2->b);
}

int is_boundary_xy(double x, double y){
    const double eps=1e-12;
    return (x<eps || x>1.0-eps || y<eps || y>1.0-eps);
}

P2Data build_p2_from_p1(const MeshP1* m) {
    int nAll = 3*m->nElem;
    Edge* all = (Edge*)malloc((size_t)nAll*sizeof(Edge));
    int* triE = (int*)malloc((size_t)3*(size_t)m->nElem*sizeof(int));
    if (!all || !triE) die("malloc build_p2_from_p1 failed");

    for (int e=0;e<m->nElem;++e) {
        int v0=m->triV[3*e+0], v1=m->triV[3*e+1], v2=m->triV[3*e+2];
        all[3*e+0]=(Edge){edge_lower(v0,v1), edge_upper(v0,v1)};
        all[3*e+1]=(Edge){edge_lower(v1,v2), edge_upper(v1,v2)};
        all[3*e+2]=(Edge){edge_lower(v2,v0), edge_upper(v2,v0)};
        triE[3*e+0]=triE[3*e+1]=triE[3*e+2]=-1;
    }

    qsort(all, (size_t)nAll, sizeof(Edge), edge_cmp);

    Edge* uniq = (Edge*)malloc((size_t)nAll*sizeof(Edge));
    if (!uniq) die("malloc uniq failed");
    int nE=0;
    for (int i=0;i<nAll;++i) {
        if (i==0 || all[i].a!=all[i-1].a || all[i].b!=all[i-1].b) uniq[nE++]=all[i];
    }

    for (int e=0;e<m->nElem;++e) {
        int v0=m->triV[3*e+0], v1=m->triV[3*e+1], v2=m->triV[3*e+2];
        Edge E[3] = {
            {edge_lower(v0,v1), edge_upper(v0,v1)},
            {edge_lower(v1,v2), edge_upper(v1,v2)},
            {edge_lower(v2,v0), edge_upper(v2,v0)}
        };
        for (int k=0;k<3;++k) {
            int lo=0, hi=nE-1, pos=-1;
            while (lo<=hi) {
                int mid=(lo+hi)/2;
                if (uniq[mid].a==E[k].a && uniq[mid].b==E[k].b){ pos=mid; break; }
                if (uniq[mid].a<E[k].a || (uniq[mid].a==E[k].a && uniq[mid].b<E[k].b)) lo=mid+1;
                else hi=mid-1;
            }
            if (pos<0) die("edge not found");
            triE[3*e+k]=pos;
        }
    }

    int nP2 = m->nV + nE;
    double* xP2 = (double*)malloc((size_t)nP2*sizeof(double));
    double* yP2 = (double*)malloc((size_t)nP2*sizeof(double));
    if (!xP2 || !yP2) die("malloc xP2/yP2 failed");

    for (int i=0;i<m->nV;++i){ xP2[i]=m->xV[i]; yP2[i]=m->yV[i]; }
    for (int ei=0;ei<nE;++ei){
        int a=uniq[ei].a, b=uniq[ei].b;
        xP2[m->nV+ei] = 0.5*(m->xV[a]+m->xV[b]);
        yP2[m->nV+ei] = 0.5*(m->yV[a]+m->yV[b]);
    }

    free(all);

    P2Data p2;
    p2.nE=nE;
    p2.edges=(Edge*)realloc(uniq, (size_t)nE*sizeof(Edge));
    p2.triE=triE;
    p2.nP2=nP2;
    p2.xP2=xP2;
    p2.yP2=yP2;
    return p2;
}

void free_p2(P2Data* p2) {
    free(p2->edges); free(p2->triE); free(p2->xP2); free(p2->yP2);
    p2->edges=NULL; p2->triE=NULL; p2->xP2=NULL; p2->yP2=NULL;
    p2->nE=0; p2->nP2=0;
}
