#ifndef MESH_P1_H
#define MESH_P1_H

typedef struct {
    int N;
    int nV;
    int nElem;
    double* xV;
    double* yV;
    int* triV; // 3*nElem : vertex ids per triangle
} MeshP1;

int vid(int i, int j, int N);

MeshP1 make_mesh_p1(int N);
void free_mesh_p1(MeshP1* m);

#endif
