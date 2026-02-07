# Stokes FEM Taylor–Hood (P2/P1) + MPI (C)

## Overview
This project solves the 2D stationary Stokes problem on the unit square using:
- Velocity: Taylor–Hood P2 (quadratic) on triangles
- Pressure: P1 (linear) on triangle vertices
- MPI for parallel element assembly (stripe partitioning)
- Restarted GMRES for the saddle-point linear system

## Model
\[
-\nu \Delta \mathbf{u} + \nabla p = \mathbf{f},\quad \nabla\cdot \mathbf{u}=0,\quad \mathbf{u}=0 \text{ on } \partial\Omega
\]

## Build
```bash
make
