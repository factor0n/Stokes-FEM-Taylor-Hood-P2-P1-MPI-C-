#ifndef P2_H
#define P2_H

#include "mesh_p1.h"

typedef struct { int a,b; } Edge; // a<b

typedef struct {
    int nE;
    Edge* edges;   // unique edges
    int* triE;     // 3*nElem edge indices for each triangle: (v0,v1),(v1,v2),(v2,v0)
    int nP2;       // scalar P2 nodes = nV + nE
    double* xP2;
    double* yP2;
} P2Data;

P2Data build_p2_from_p1(const MeshP1* m);
void free_p2(P2Data* p2);

int is_boundary_xy(double x, double y);

#endif
