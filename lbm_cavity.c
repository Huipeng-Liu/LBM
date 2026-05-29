#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NX 101
#define NY 101
#define Q 9

#define RE 100.0
#define U_LID 0.1
#define MAX_ITERS 20000
#define OUTPUT_INTERVAL 1000
#define TOL 1.0e-8

#define IDX(x, y, k) ((((y) * NX + (x)) * Q) + (k))
#define NODE(x, y) ((y) * NX + (x))

static const int cx[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
static const int cy[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
static const int opp[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
static const double w[Q] = {
    4.0 / 9.0,
    1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
    1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
};

static double equilibrium(int k, double rho, double ux, double uy)
{
    double cu = 3.0 * (cx[k] * ux + cy[k] * uy);
    double u2 = ux * ux + uy * uy;
    return w[k] * rho * (1.0 + cu + 0.5 * cu * cu - 1.5 * u2);
}

static int is_solid_wall(int x, int y)
{
    return x == 0 || x == NX - 1 || y == 0;
}

static void initialize(double *f, double *rho, double *ux, double *uy)
{
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int n = NODE(x, y);
            rho[n] = 1.0;
            ux[n] = (y == NY - 1 && x > 0 && x < NX - 1) ? U_LID : 0.0;
            uy[n] = 0.0;

            for (int k = 0; k < Q; ++k) {
                f[IDX(x, y, k)] = equilibrium(k, rho[n], ux[n], uy[n]);
            }
        }
    }
}

static void collide_stream(const double *f, double *f_next, double omega)
{
    for (int i = 0; i < NX * NY * Q; ++i) {
        f_next[i] = 0.0;
    }

    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            if (is_solid_wall(x, y)) {
                continue;
            }

            double rho = 0.0;
            double ux = 0.0;
            double uy = 0.0;

            for (int k = 0; k < Q; ++k) {
                double fk = f[IDX(x, y, k)];
                rho += fk;
                ux += fk * cx[k];
                uy += fk * cy[k];
            }

            ux /= rho;
            uy /= rho;

            if (y == NY - 1) {
                ux = U_LID;
                uy = 0.0;
            }

            for (int k = 0; k < Q; ++k) {
                double feq = equilibrium(k, rho, ux, uy);
                double f_post = f[IDX(x, y, k)] - omega * (f[IDX(x, y, k)] - feq);
                int xn = x + cx[k];
                int yn = y + cy[k];

                if (xn < 0 || xn >= NX || yn < 0 || yn >= NY || is_solid_wall(xn, yn)) {
                    f_next[IDX(x, y, opp[k])] += f_post;
                } else {
                    f_next[IDX(xn, yn, k)] += f_post;
                }
            }
        }
    }
}

static void apply_lid_velocity(double *f)
{
    int y = NY - 1;

    for (int x = 1; x < NX - 1; ++x) {
        double rho = f[IDX(x, y, 0)] + f[IDX(x, y, 1)] + f[IDX(x, y, 3)]
                   + 2.0 * (f[IDX(x, y, 2)] + f[IDX(x, y, 5)] + f[IDX(x, y, 6)]);

        f[IDX(x, y, 4)] = f[IDX(x, y, 2)];
        f[IDX(x, y, 7)] = f[IDX(x, y, 5)] - 0.5 * (f[IDX(x, y, 1)] - f[IDX(x, y, 3)])
                        - (1.0 / 6.0) * rho * U_LID;
        f[IDX(x, y, 8)] = f[IDX(x, y, 6)] + 0.5 * (f[IDX(x, y, 1)] - f[IDX(x, y, 3)])
                        + (1.0 / 6.0) * rho * U_LID;
    }

    for (int k = 0; k < Q; ++k) {
        f[IDX(0, y, k)] = equilibrium(k, 1.0, 0.0, 0.0);
        f[IDX(NX - 1, y, k)] = equilibrium(k, 1.0, 0.0, 0.0);
    }
}

static double update_macroscopic(const double *f, double *rho, double *ux, double *uy,
                                 double *ux_old, double *uy_old)
{
    double diff2 = 0.0;
    double norm2 = 0.0;

    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int n = NODE(x, y);

            if (is_solid_wall(x, y)) {
                rho[n] = 1.0;
                ux[n] = 0.0;
                uy[n] = 0.0;
            } else if (y == NY - 1) {
                rho[n] = 0.0;
                for (int k = 0; k < Q; ++k) {
                    rho[n] += f[IDX(x, y, k)];
                }
                ux[n] = U_LID;
                uy[n] = 0.0;
            } else {
                rho[n] = 0.0;
                ux[n] = 0.0;
                uy[n] = 0.0;

                for (int k = 0; k < Q; ++k) {
                    double fk = f[IDX(x, y, k)];
                    rho[n] += fk;
                    ux[n] += fk * cx[k];
                    uy[n] += fk * cy[k];
                }

                ux[n] /= rho[n];
                uy[n] /= rho[n];
            }

            double du = ux[n] - ux_old[n];
            double dv = uy[n] - uy_old[n];
            diff2 += du * du + dv * dv;
            norm2 += ux[n] * ux[n] + uy[n] * uy[n];

            ux_old[n] = ux[n];
            uy_old[n] = uy[n];
        }
    }

    if (norm2 <= 0.0) {
        return sqrt(diff2);
    }

    return sqrt(diff2 / norm2);
}

static int has_invalid_values(const double *rho, const double *ux, const double *uy)
{
    for (int i = 0; i < NX * NY; ++i) {
        if (!isfinite(rho[i]) || !isfinite(ux[i]) || !isfinite(uy[i])) {
            return 1;
        }
    }
    return 0;
}

static int write_csv(const char *path, const double *rho, const double *ux, const double *uy)
{
    FILE *fp = fopen(path, "w");
    if (!fp) {
        return 0;
    }

    fprintf(fp, "x,y,rho,u,v,velocity_magnitude\n");
    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int n = NODE(x, y);
            double xn = (double)x / (double)(NX - 1);
            double yn = (double)y / (double)(NY - 1);
            double speed = sqrt(ux[n] * ux[n] + uy[n] * uy[n]);
            fprintf(fp, "%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
                    xn, yn, rho[n], ux[n], uy[n], speed);
        }
    }

    fclose(fp);
    return 1;
}

int main(void)
{
    const double length = (double)(NX - 1);
    const double nu = U_LID * length / RE;
    const double tau = 3.0 * nu + 0.5;
    const double omega = 1.0 / tau;
    const char *output_path = "lid_driven_cavity.csv";

    if (tau <= 0.5) {
        fprintf(stderr, "Invalid relaxation time: tau = %.12g must be greater than 0.5\n", tau);
        return EXIT_FAILURE;
    }

    double *f = (double *)malloc((size_t)NX * NY * Q * sizeof(double));
    double *f_next = (double *)malloc((size_t)NX * NY * Q * sizeof(double));
    double *rho = (double *)malloc((size_t)NX * NY * sizeof(double));
    double *ux = (double *)malloc((size_t)NX * NY * sizeof(double));
    double *uy = (double *)malloc((size_t)NX * NY * sizeof(double));
    double *ux_old = (double *)calloc((size_t)NX * NY, sizeof(double));
    double *uy_old = (double *)calloc((size_t)NX * NY, sizeof(double));

    if (!f || !f_next || !rho || !ux || !uy || !ux_old || !uy_old) {
        fprintf(stderr, "Memory allocation failed.\n");
        free(f);
        free(f_next);
        free(rho);
        free(ux);
        free(uy);
        free(ux_old);
        free(uy_old);
        return EXIT_FAILURE;
    }

    printf("2D lid-driven cavity LBM (D2Q9 BGK)\n");
    printf("NX=%d NY=%d Re=%.6g U_lid=%.6g nu=%.12g tau=%.12g omega=%.12g\n",
           NX, NY, RE, U_LID, nu, tau, omega);

    initialize(f, rho, ux, uy);
    for (int i = 0; i < NX * NY; ++i) {
        ux_old[i] = ux[i];
        uy_old[i] = uy[i];
    }

    double residual = 0.0;
    int iter = 0;
    for (iter = 1; iter <= MAX_ITERS; ++iter) {
        collide_stream(f, f_next, omega);
        apply_lid_velocity(f_next);

        double *tmp = f;
        f = f_next;
        f_next = tmp;

        residual = update_macroscopic(f, rho, ux, uy, ux_old, uy_old);

        if (has_invalid_values(rho, ux, uy)) {
            fprintf(stderr, "Invalid numeric value detected at iteration %d.\n", iter);
            free(f);
            free(f_next);
            free(rho);
            free(ux);
            free(uy);
            free(ux_old);
            free(uy_old);
            return EXIT_FAILURE;
        }

        if (iter % OUTPUT_INTERVAL == 0 || residual < TOL) {
            printf("iter=%d residual=%.12e center_u=%.12g center_v=%.12g\n",
                   iter, residual, ux[NODE(NX / 2, NY / 2)], uy[NODE(NX / 2, NY / 2)]);
        }

        if (residual < TOL) {
            break;
        }
    }

    if (!write_csv(output_path, rho, ux, uy)) {
        fprintf(stderr, "Failed to write %s\n", output_path);
        free(f);
        free(f_next);
        free(rho);
        free(ux);
        free(uy);
        free(ux_old);
        free(uy_old);
        return EXIT_FAILURE;
    }

    printf("Finished at iteration %d with residual %.12e\n", iter <= MAX_ITERS ? iter : MAX_ITERS, residual);
    printf("Wrote %s (%d data rows)\n", output_path, NX * NY);

    free(f);
    free(f_next);
    free(rho);
    free(ux);
    free(uy);
    free(ux_old);
    free(uy_old);

    return EXIT_SUCCESS;
}
