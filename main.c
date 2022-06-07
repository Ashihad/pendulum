#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define g 9.81
#define R 1.0
#define m 1.0

void pochodne(double t, double *s, double *k) {
    k[0]=s[1];
    k[1]=-g/R*sin(s[0]);
}

void rk4_vec(double t, double dt, int n, double *s, void (*f)(double , double *, double  *)) {
    #define M 1000
    static double k1[M], k2[M], k3[M], k4[M], w[M];
    int i;
    // Kopia tablicy
    for(i = 0; i < n; i++) {
        w[i] = s[i];
    }
    // get rk1
    f(t, w, k1);
    for(i = 0; i < n; i++) {
        w[i] = s[i] + dt/2*k1[i];
    }
    // get rk2
    f(t+dt/2, w, k2);
    for(i = 0; i < n; i++) {
        w[i] = s[i] + dt/2*k2[i];
    }
    // get rk3
    f(t+dt/2, w, k3);
    for(i = 0; i < n; i++) {
        w[i] = s[i] + dt*k3[i];
    }
    // get rk4
    f(t+dt, w, k4);
    // do s[]
    for(i = 0; i < n; i++) {
        s[i]= s[i] + dt/6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    }
}

int main() {
    // inicjalizacja  parametrów:
    int n = 2;              // ilość zmiennych w układzie  RRZ1
    double dt = 0.01;
    double tmax = 10;
    int N = (int)tmax/dt;    // ilość kroków  czasowych
    double t = 0;

    void (*f)(double, double *, double *);
    f = pochodne;

    double *s;
    s = (double *)malloc(n*sizeof(double));    // tablica  rozwiązań

    // warunki początkowe
    s[0] = 4 * 3.14/180; // degrees to radians
    s[1] = 0;

    // plik
    FILE *fp;
    fp = fopen("wyniki.txt", "w");
    if(fp == NULL) {
        printf("Error: No file found\n");
        return EXIT_FAILURE;
    }

    for(size_t i = 1; i <= N; i++) {
        rk4_vec(t, dt, n, s, f);
        t=t+dt;
        fprintf(fp, "%f %f\n", t, s[0]);
    }
    fclose(fp);
    return EXIT_SUCCESS;
}