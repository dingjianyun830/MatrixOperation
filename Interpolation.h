#pragma once
#include <math.h>
double Lagrange(double *x, double *y, int n, double t);
double Pqrat(double *x, double *y, int n, double t);
double Hermite(double *x, double *y, double *dy, int n, double t);

double eps;
double Aitken(double *x, double*y, int n, double t);

// 3-nd ployline
double Akima(double *x, double *y, int n, double t);

// 2d lagrange interpolation.
double Lagrange2(double *x, double *y, double **z, int n, int m, double u, double v);

// Least Squares
bool FitLS(double *x, double *y, int n, double *a, int m, double *dt);

// 2d Least Squares
bool FitLS2(double *x, double *y, double **z, int n, int m, double **a, int p, int q, double *dt);