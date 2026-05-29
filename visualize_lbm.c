#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NX 101
#define NY 101
#define WIDTH 900
#define HEIGHT 900
#define MARGIN 70
#define ARROW_STEP 5

#define NODE(x, y) ((y) * NX + (x))

static void turbo_color(double t, int *r, int *g, int *b)
{
    double x = t;
    if (x < 0.0) {
        x = 0.0;
    } else if (x > 1.0) {
        x = 1.0;
    }

    double rr = 34.61 + x * (1172.33 + x * (-10793.56 + x * (33300.12 + x * (-38394.49 + x * 14825.05))));
    double gg = 23.31 + x * (557.33 + x * (1225.33 + x * (-3574.96 + x * (1073.77 + x * 707.56))));
    double bb = 27.20 + x * (3211.10 + x * (-15327.97 + x * (27814.00 + x * (-22569.18 + x * 6838.66))));

    if (rr < 0.0) rr = 0.0;
    if (rr > 255.0) rr = 255.0;
    if (gg < 0.0) gg = 0.0;
    if (gg > 255.0) gg = 255.0;
    if (bb < 0.0) bb = 0.0;
    if (bb > 255.0) bb = 255.0;

    *r = (int)(rr + 0.5);
    *g = (int)(gg + 0.5);
    *b = (int)(bb + 0.5);
}

static int read_csv(const char *path, double *rho, double *ux, double *uy, double *speed)
{
    FILE *fp = fopen(path, "r");
    if (!fp) {
        return 0;
    }

    char line[512];
    if (!fgets(line, sizeof(line), fp)) {
        fclose(fp);
        return 0;
    }

    int n = 0;
    while (fgets(line, sizeof(line), fp) && n < NX * NY) {
        double x_norm, y_norm;
        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf",
                   &x_norm, &y_norm, &rho[n], &ux[n], &uy[n], &speed[n]) != 6) {
            fclose(fp);
            return 0;
        }
        ++n;
    }

    fclose(fp);
    return n == NX * NY;
}

static void write_arrow(FILE *fp, double x0, double y0, double dx, double dy, double scale)
{
    double len = sqrt(dx * dx + dy * dy);
    if (len <= 1.0e-14) {
        return;
    }

    double x1 = x0 + dx * scale;
    double y1 = y0 - dy * scale;
    double ux = (x1 - x0) / sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
    double uy = (y1 - y0) / sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
    double head = 6.0;
    double wing = 3.5;
    double px = -uy;
    double py = ux;
    double ax = x1 - ux * head + px * wing;
    double ay = y1 - uy * head + py * wing;
    double bx = x1 - ux * head - px * wing;
    double by = y1 - uy * head - py * wing;

    fprintf(fp, "<path d=\"M %.3f %.3f L %.3f %.3f M %.3f %.3f L %.3f %.3f L %.3f %.3f\" "
                "stroke=\"#111827\" stroke-width=\"1.15\" fill=\"none\" stroke-linecap=\"round\" "
                "stroke-linejoin=\"round\" opacity=\"0.72\"/>\n",
            x0, y0, x1, y1, ax, ay, x1, y1, bx, by);
}

static int write_svg(const char *path, const double *rho, const double *ux,
                     const double *uy, const double *speed)
{
    (void)rho;

    double max_speed = 0.0;
    for (int i = 0; i < NX * NY; ++i) {
        if (speed[i] > max_speed) {
            max_speed = speed[i];
        }
    }

    FILE *fp = fopen(path, "w");
    if (!fp) {
        return 0;
    }

    const double plot = WIDTH - 2.0 * MARGIN;
    const double cell = plot / (double)(NX - 1);
    const double arrow_scale = 0.42 * cell * ARROW_STEP / (max_speed > 0.0 ? max_speed : 1.0);

    fprintf(fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\">\n",
            WIDTH, HEIGHT, WIDTH, HEIGHT);
    fprintf(fp, "<rect width=\"100%%\" height=\"100%%\" fill=\"#ffffff\"/>\n");
    fprintf(fp, "<text x=\"%d\" y=\"32\" font-family=\"Arial, sans-serif\" font-size=\"22\" "
                "font-weight=\"700\" fill=\"#111827\">LBM Lid-Driven Cavity Velocity Field</text>\n",
            MARGIN);
    fprintf(fp, "<text x=\"%d\" y=\"55\" font-family=\"Arial, sans-serif\" font-size=\"13\" "
                "fill=\"#374151\">D2Q9 BGK, Re=100, lid velocity=0.1. Color shows speed; arrows show velocity direction.</text>\n",
            MARGIN);

    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            int n = NODE(x, y);
            double t = speed[n] / (max_speed > 0.0 ? max_speed : 1.0);
            int r, g, b;
            turbo_color(t, &r, &g, &b);

            double sx = MARGIN + x * cell - 0.5 * cell;
            double sy = MARGIN + (NY - 1 - y) * cell - 0.5 * cell;
            fprintf(fp, "<rect x=\"%.3f\" y=\"%.3f\" width=\"%.3f\" height=\"%.3f\" "
                        "fill=\"rgb(%d,%d,%d)\"/>\n",
                    sx, sy, cell + 0.8, cell + 0.8, r, g, b);
        }
    }

    fprintf(fp, "<rect x=\"%d\" y=\"%d\" width=\"%.3f\" height=\"%.3f\" fill=\"none\" "
                "stroke=\"#111827\" stroke-width=\"1.5\"/>\n",
            MARGIN, MARGIN, plot, plot);

    for (int y = ARROW_STEP; y < NY - 1; y += ARROW_STEP) {
        for (int x = ARROW_STEP; x < NX - 1; x += ARROW_STEP) {
            int n = NODE(x, y);
            double sx = MARGIN + x * cell;
            double sy = MARGIN + (NY - 1 - y) * cell;
            write_arrow(fp, sx, sy, ux[n], uy[n], arrow_scale);
        }
    }

    const int bar_x = WIDTH - MARGIN + 24;
    const int bar_y = MARGIN;
    const int bar_w = 22;
    const int bar_h = WIDTH - 2 * MARGIN;
    for (int i = 0; i < bar_h; ++i) {
        double t = 1.0 - (double)i / (double)(bar_h - 1);
        int r, g, b;
        turbo_color(t, &r, &g, &b);
        fprintf(fp, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"1\" fill=\"rgb(%d,%d,%d)\"/>\n",
                bar_x, bar_y + i, bar_w, r, g, b);
    }

    fprintf(fp, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"none\" "
                "stroke=\"#111827\" stroke-width=\"1\"/>\n", bar_x, bar_y, bar_w, bar_h);
    fprintf(fp, "<text x=\"%d\" y=\"%d\" font-family=\"Arial, sans-serif\" font-size=\"12\" "
                "fill=\"#111827\">%.4f</text>\n", bar_x + 30, bar_y + 5, max_speed);
    fprintf(fp, "<text x=\"%d\" y=\"%d\" font-family=\"Arial, sans-serif\" font-size=\"12\" "
                "fill=\"#111827\">0</text>\n", bar_x + 30, bar_y + bar_h);
    fprintf(fp, "<text x=\"%d\" y=\"%d\" font-family=\"Arial, sans-serif\" font-size=\"12\" "
                "fill=\"#111827\" transform=\"rotate(90 %d %d)\">speed magnitude</text>\n",
            bar_x + 54, bar_y + bar_h / 2 - 45, bar_x + 54, bar_y + bar_h / 2 - 45);

    fprintf(fp, "<text x=\"%d\" y=\"%d\" font-family=\"Arial, sans-serif\" font-size=\"13\" "
                "fill=\"#374151\">x</text>\n", MARGIN + (int)(plot / 2.0), HEIGHT - 22);
    fprintf(fp, "<text x=\"22\" y=\"%d\" font-family=\"Arial, sans-serif\" font-size=\"13\" "
                "fill=\"#374151\" transform=\"rotate(-90 22 %d)\">y</text>\n",
            MARGIN + (int)(plot / 2.0), MARGIN + (int)(plot / 2.0));

    fprintf(fp, "</svg>\n");
    fclose(fp);
    return 1;
}

int main(void)
{
    double *rho = (double *)malloc((size_t)NX * NY * sizeof(double));
    double *ux = (double *)malloc((size_t)NX * NY * sizeof(double));
    double *uy = (double *)malloc((size_t)NX * NY * sizeof(double));
    double *speed = (double *)malloc((size_t)NX * NY * sizeof(double));

    if (!rho || !ux || !uy || !speed) {
        fprintf(stderr, "Memory allocation failed.\n");
        free(rho);
        free(ux);
        free(uy);
        free(speed);
        return EXIT_FAILURE;
    }

    if (!read_csv("lid_driven_cavity.csv", rho, ux, uy, speed)) {
        fprintf(stderr, "Failed to read lid_driven_cavity.csv. Run lbm_cavity first.\n");
        free(rho);
        free(ux);
        free(uy);
        free(speed);
        return EXIT_FAILURE;
    }

    if (!write_svg("velocity_field.svg", rho, ux, uy, speed)) {
        fprintf(stderr, "Failed to write velocity_field.svg.\n");
        free(rho);
        free(ux);
        free(uy);
        free(speed);
        return EXIT_FAILURE;
    }

    printf("Wrote velocity_field.svg\n");

    free(rho);
    free(ux);
    free(uy);
    free(speed);
    return EXIT_SUCCESS;
}
