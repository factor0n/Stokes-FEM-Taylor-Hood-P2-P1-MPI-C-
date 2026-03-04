#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Minimal 3D incompressible fluid simulation (stable fluids style):
 * - Semi-Lagrangian advection for velocity and density
 * - Jacobi pressure projection for incompressibility
 * - Simple smoke source and buoyancy force
 *
 * Build:
 *   gcc -O2 -std=c11 examples/fluid3d_semi_lagrangian.c -lm -o fluid3d
 * Run:
 *   ./fluid3d
 */

typedef struct {
    int nx, ny, nz;
    float dt;
    float viscosity;
    int pressure_iters;

    float *u, *v, *w;
    float *u0, *v0, *w0;
    float *density, *density0;
    float *pressure, *div;
} Fluid3D;

static int idx3(const Fluid3D *f, int x, int y, int z) {
    return x + f->nx * (y + f->ny * z);
}

static int clampi(int a, int lo, int hi) {
    return a < lo ? lo : (a > hi ? hi : a);
}

static float sample_trilinear(const Fluid3D *f, const float *field, float x, float y, float z) {
    x = fmaxf(0.0f, fminf((float)(f->nx - 1), x));
    y = fmaxf(0.0f, fminf((float)(f->ny - 1), y));
    z = fmaxf(0.0f, fminf((float)(f->nz - 1), z));

    int x0 = (int)floorf(x), y0 = (int)floorf(y), z0 = (int)floorf(z);
    int x1 = clampi(x0 + 1, 0, f->nx - 1);
    int y1 = clampi(y0 + 1, 0, f->ny - 1);
    int z1 = clampi(z0 + 1, 0, f->nz - 1);

    float tx = x - x0, ty = y - y0, tz = z - z0;

    float c000 = field[idx3(f, x0, y0, z0)];
    float c100 = field[idx3(f, x1, y0, z0)];
    float c010 = field[idx3(f, x0, y1, z0)];
    float c110 = field[idx3(f, x1, y1, z0)];
    float c001 = field[idx3(f, x0, y0, z1)];
    float c101 = field[idx3(f, x1, y0, z1)];
    float c011 = field[idx3(f, x0, y1, z1)];
    float c111 = field[idx3(f, x1, y1, z1)];

    float c00 = c000 * (1.0f - tx) + c100 * tx;
    float c10 = c010 * (1.0f - tx) + c110 * tx;
    float c01 = c001 * (1.0f - tx) + c101 * tx;
    float c11 = c011 * (1.0f - tx) + c111 * tx;

    float c0 = c00 * (1.0f - ty) + c10 * ty;
    float c1 = c01 * (1.0f - ty) + c11 * ty;

    return c0 * (1.0f - tz) + c1 * tz;
}

static void set_boundary_velocity(Fluid3D *f) {
    int nx = f->nx, ny = f->ny, nz = f->nz;
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                if (x == 0 || x == nx - 1 || y == 0 || y == ny - 1 || z == 0 || z == nz - 1) {
                    int id = idx3(f, x, y, z);
                    f->u[id] = f->v[id] = f->w[id] = 0.0f;
                }
            }
        }
    }
}

static void add_source_and_forces(Fluid3D *f) {
    int cx = f->nx / 2;
    int cy = f->ny / 6;
    int cz = f->nz / 2;

    for (int z = cz - 2; z <= cz + 2; ++z) {
        for (int y = cy - 1; y <= cy + 1; ++y) {
            for (int x = cx - 2; x <= cx + 2; ++x) {
                int xx = clampi(x, 1, f->nx - 2);
                int yy = clampi(y, 1, f->ny - 2);
                int zz = clampi(z, 1, f->nz - 2);
                int id = idx3(f, xx, yy, zz);
                f->density[id] += 5.0f * f->dt;
                f->w[id] += 1.5f * f->dt;
            }
        }
    }

    for (int i = 0; i < f->nx * f->ny * f->nz; ++i) {
        f->w[i] += 0.05f * f->density[i] * f->dt;
    }
}

static void advect_scalar(Fluid3D *f, float *dst, const float *src) {
    for (int z = 1; z < f->nz - 1; ++z) {
        for (int y = 1; y < f->ny - 1; ++y) {
            for (int x = 1; x < f->nx - 1; ++x) {
                int id = idx3(f, x, y, z);
                float xb = x - f->dt * f->u[id];
                float yb = y - f->dt * f->v[id];
                float zb = z - f->dt * f->w[id];
                dst[id] = sample_trilinear(f, src, xb, yb, zb);
            }
        }
    }
}

static void advect_velocity(Fluid3D *f) {
    memcpy(f->u0, f->u, (size_t)f->nx * f->ny * f->nz * sizeof(float));
    memcpy(f->v0, f->v, (size_t)f->nx * f->ny * f->nz * sizeof(float));
    memcpy(f->w0, f->w, (size_t)f->nx * f->ny * f->nz * sizeof(float));

    for (int z = 1; z < f->nz - 1; ++z) {
        for (int y = 1; y < f->ny - 1; ++y) {
            for (int x = 1; x < f->nx - 1; ++x) {
                int id = idx3(f, x, y, z);
                float xb = x - f->dt * f->u0[id];
                float yb = y - f->dt * f->v0[id];
                float zb = z - f->dt * f->w0[id];
                f->u[id] = sample_trilinear(f, f->u0, xb, yb, zb);
                f->v[id] = sample_trilinear(f, f->v0, xb, yb, zb);
                f->w[id] = sample_trilinear(f, f->w0, xb, yb, zb);
            }
        }
    }
}

static void diffuse_velocity(Fluid3D *f) {
    if (f->viscosity <= 0.0f) return;

    const float a = f->viscosity * f->dt;
    const int n = f->nx * f->ny * f->nz;
    memcpy(f->u0, f->u, (size_t)n * sizeof(float));
    memcpy(f->v0, f->v, (size_t)n * sizeof(float));
    memcpy(f->w0, f->w, (size_t)n * sizeof(float));

    for (int iter = 0; iter < 8; ++iter) {
        for (int z = 1; z < f->nz - 1; ++z) {
            for (int y = 1; y < f->ny - 1; ++y) {
                for (int x = 1; x < f->nx - 1; ++x) {
                    int id = idx3(f, x, y, z);
                    int xp = idx3(f, x + 1, y, z), xm = idx3(f, x - 1, y, z);
                    int yp = idx3(f, x, y + 1, z), ym = idx3(f, x, y - 1, z);
                    int zp = idx3(f, x, y, z + 1), zm = idx3(f, x, y, z - 1);

                    f->u[id] = (f->u0[id] + a * (f->u[xp] + f->u[xm] + f->u[yp] + f->u[ym] + f->u[zp] + f->u[zm])) / (1.0f + 6.0f * a);
                    f->v[id] = (f->v0[id] + a * (f->v[xp] + f->v[xm] + f->v[yp] + f->v[ym] + f->v[zp] + f->v[zm])) / (1.0f + 6.0f * a);
                    f->w[id] = (f->w0[id] + a * (f->w[xp] + f->w[xm] + f->w[yp] + f->w[ym] + f->w[zp] + f->w[zm])) / (1.0f + 6.0f * a);
                }
            }
        }
        set_boundary_velocity(f);
    }
}

static void project(Fluid3D *f) {
    for (int z = 1; z < f->nz - 1; ++z) {
        for (int y = 1; y < f->ny - 1; ++y) {
            for (int x = 1; x < f->nx - 1; ++x) {
                int id = idx3(f, x, y, z);
                int xp = idx3(f, x + 1, y, z), xm = idx3(f, x - 1, y, z);
                int yp = idx3(f, x, y + 1, z), ym = idx3(f, x, y - 1, z);
                int zp = idx3(f, x, y, z + 1), zm = idx3(f, x, y, z - 1);

                f->div[id] = -0.5f * (f->u[xp] - f->u[xm] + f->v[yp] - f->v[ym] + f->w[zp] - f->w[zm]);
                f->pressure[id] = 0.0f;
            }
        }
    }

    for (int it = 0; it < f->pressure_iters; ++it) {
        for (int z = 1; z < f->nz - 1; ++z) {
            for (int y = 1; y < f->ny - 1; ++y) {
                for (int x = 1; x < f->nx - 1; ++x) {
                    int id = idx3(f, x, y, z);
                    int xp = idx3(f, x + 1, y, z), xm = idx3(f, x - 1, y, z);
                    int yp = idx3(f, x, y + 1, z), ym = idx3(f, x, y - 1, z);
                    int zp = idx3(f, x, y, z + 1), zm = idx3(f, x, y, z - 1);

                    f->pressure[id] = (f->div[id] + f->pressure[xp] + f->pressure[xm] + f->pressure[yp] + f->pressure[ym] + f->pressure[zp] + f->pressure[zm]) / 6.0f;
                }
            }
        }
    }

    for (int z = 1; z < f->nz - 1; ++z) {
        for (int y = 1; y < f->ny - 1; ++y) {
            for (int x = 1; x < f->nx - 1; ++x) {
                int id = idx3(f, x, y, z);
                int xp = idx3(f, x + 1, y, z), xm = idx3(f, x - 1, y, z);
                int yp = idx3(f, x, y + 1, z), ym = idx3(f, x, y - 1, z);
                int zp = idx3(f, x, y, z + 1), zm = idx3(f, x, y, z - 1);

                f->u[id] -= 0.5f * (f->pressure[xp] - f->pressure[xm]);
                f->v[id] -= 0.5f * (f->pressure[yp] - f->pressure[ym]);
                f->w[id] -= 0.5f * (f->pressure[zp] - f->pressure[zm]);
            }
        }
    }

    set_boundary_velocity(f);
}

static void step(Fluid3D *f) {
    add_source_and_forces(f);
    diffuse_velocity(f);
    advect_velocity(f);
    project(f);

    memcpy(f->density0, f->density, (size_t)f->nx * f->ny * f->nz * sizeof(float));
    advect_scalar(f, f->density, f->density0);
}


static void write_vtk(const Fluid3D *f, const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Cannot open %s for writing\n", filename);
        return;
    }

    int n = f->nx * f->ny * f->nz;
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "Fluid3D output\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS\n");
    fprintf(fp, "DIMENSIONS %d %d %d\n", f->nx, f->ny, f->nz);
    fprintf(fp, "ORIGIN 0 0 0\n");
    fprintf(fp, "SPACING 1 1 1\n");
    fprintf(fp, "POINT_DATA %d\n", n);

    fprintf(fp, "SCALARS density float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n; ++i) {
        fprintf(fp, "%g\n", f->density[i]);
    }

    fprintf(fp, "VECTORS velocity float\n");
    for (int i = 0; i < n; ++i) {
        fprintf(fp, "%g %g %g\n", f->u[i], f->v[i], f->w[i]);
    }

    fclose(fp);
}

static float density_sum(const Fluid3D *f) {
    float s = 0.0f;
    int n = f->nx * f->ny * f->nz;
    for (int i = 0; i < n; ++i) s += f->density[i];
    return s;
}

static void init_fluid(Fluid3D *f, int nx, int ny, int nz) {
    memset(f, 0, sizeof(*f));
    f->nx = nx;
    f->ny = ny;
    f->nz = nz;
    f->dt = 0.1f;
    f->viscosity = 0.001f;
    f->pressure_iters = 40;

    int n = nx * ny * nz;
    f->u = (float *)calloc((size_t)n, sizeof(float));
    f->v = (float *)calloc((size_t)n, sizeof(float));
    f->w = (float *)calloc((size_t)n, sizeof(float));
    f->u0 = (float *)calloc((size_t)n, sizeof(float));
    f->v0 = (float *)calloc((size_t)n, sizeof(float));
    f->w0 = (float *)calloc((size_t)n, sizeof(float));
    f->density = (float *)calloc((size_t)n, sizeof(float));
    f->density0 = (float *)calloc((size_t)n, sizeof(float));
    f->pressure = (float *)calloc((size_t)n, sizeof(float));
    f->div = (float *)calloc((size_t)n, sizeof(float));

    if (!f->u || !f->v || !f->w || !f->u0 || !f->v0 || !f->w0 ||
        !f->density || !f->density0 || !f->pressure || !f->div) {
        fprintf(stderr, "Allocation failed\n");
        exit(EXIT_FAILURE);
    }
}

static void free_fluid(Fluid3D *f) {
    free(f->u); free(f->v); free(f->w);
    free(f->u0); free(f->v0); free(f->w0);
    free(f->density); free(f->density0);
    free(f->pressure); free(f->div);
}

int main(void) {
    Fluid3D fluid;
    init_fluid(&fluid, 32, 32, 32);

    for (int frame = 0; frame < 120; ++frame) {
        step(&fluid);
        if (frame % 10 == 0) {
            printf("step=%d total_density=%.3f\n", frame, density_sum(&fluid));
        }
        if (frame % 20 == 0) {
            char filename[64];
            snprintf(filename, sizeof(filename), "fluid3d_%04d.vtk", frame);
            write_vtk(&fluid, filename);
            printf("wrote %s\n", filename);
        }
    }

    free_fluid(&fluid);
    return 0;
}

