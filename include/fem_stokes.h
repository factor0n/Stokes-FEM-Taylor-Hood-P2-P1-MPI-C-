#ifndef FEM_STOKES_H
#define FEM_STOKES_H

#include "mesh_p1.h"
#include "p2.h"
#include "csr.h"
#include <mpi.h>

CSR build_stokes_pattern(const MeshP1* m, const P2Data* p2);

void assemble_stokes_local(
    const MeshP1* m, const P2Data* p2,
    CSR* A_local_same_pattern, double* rhs_local,
    int rank, int size, double nu
);

void apply_dirichlet_velocity(CSR* A, double* rhs, const MeshP1* m, const P2Data* p2);
void fix_pressure_dof(CSR* A, double* rhs, int pressure_vertex_id, int nP2);

// MPI partition helper
void elem_owner_centroid_x(const MeshP1* m, int e, int size, int* owner);

#endif
