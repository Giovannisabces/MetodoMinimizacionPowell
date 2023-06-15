#include <stdio.h>
#include <math.h>

double min_unid(double (*fun)(double,double*,double*), double* x, double a, double* pk);
void triplete_inicial(double (*fun)(double,double*,double*), double* x0, double* alpha, double* pk, double* a, double* b, double* c);
double rosenb(double* x);
void vectorUpdate(double* xini, double* x0, int m);

double f(double alfa, double *x0, double *pk);
#define M 2

void mini_powell_2(double (*fun)(double,double*,double*), double* x0, int dim, double* result, int* iter)
{
    int m = dim;
    printf("Entrada x0[0]=%lf\tx0[1]=%lf\n",x0[0],x0[1]);
    double direc[m][m];
    for(int i=0; i<m ;i++)                  //  Equivalente eye
        for (int j=0;j<m;j++){
            direc[i][j] = (i==j)?1.0 : 0.0;
        }

    int fin = 0;                            //  Variable para indicar el fin del while, cuando la diferencia es inferior a  1e-5
    *iter = 0;                              //  Variable utilizada para contar las iteraciones.
    double x1[m];
    double xini[m];
    double pk[m];
    while (fin == 0) {
        vectorUpdate(xini,x0,m);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                pk[j] = direc[j][i];
            }

            double alfa_b = min_unid(fun, x0, 0, pk);
            for (int j = 0; j < m; j++) {
                x1[j] = x0[j] + alfa_b * pk[j];
                x0[j] = x1[j];
            }
        }

        for (int i = 0; i < m - 1; i++) {
            for (int j = 0; j < m; j++) {
                direc[j][i] = direc[j][i + 1];
            }
        }

        for (int i = 0; i < m; i++) {
            direc[i][m - 1] = x1[i] - xini[i];
        }


        for (int j = 0; j < m; j++) {
            pk[j] = direc[j][m - 1];
        }

        double alfa_b = min_unid(fun, x0, 0, pk);
        for (int j = 0; j < m; j++) {
            x1[j] = x0[j] + alfa_b * pk[j];
        }

        *iter = *iter + 1;

        double norm = 0.0;
        for (int j = 0; j < m; j++) {
            norm += pow(x1[j] - xini[j], 2);
        }

        norm = sqrt(norm);
        if (norm < 1e-5&& *iter>8) {
            fin = 1;
        }

        for (int j = 0; j < m; j++) {
            x0[j] = x1[j];
        }
    }
    for (int i = 0; i < m; i++) {
        result[i] = x0[i];
    }
}

double min_unid(double (*fun)(double, double*,double*), double* x, double alpha, double* pk)
{
    double a, b, c, d;
    triplete_inicial(fun, x, &alpha, pk, &a, &b, &c);

    double tol = 1e-7;
    double resphi = 2 - (1 + sqrt(5)) / 2.0;
    c = a + ( b - a) * resphi;
    d = b - ( b - a) * resphi;
    double fc = fun(c,x,pk);
    double fd = fun(d,x,pk);
    while (fabs(b - a) > tol * (fabs(c) + fabs(d))) {
        if (fc < fd) {
            b = d;
            d = c;
            c = a + (b - a) * resphi;
            fd = fc;
            fc = fun(c,x,pk);
        } else {
            a = c;
            c = d;
            d = b - (b - a) * resphi;
            fc = fd;
            fd = fun(d,x,pk);
        }
    }

    return (a + b) / 2.0;
}
/***********************************
 * Llamado de la función tipo:
 * [b, c, a]=triplete_inicial(fun,a);
 *
 * Recepcion de la función:
 * function [a,b,c] = triplete_inicial(fun,x0)
 * Equivalencias dentro de la fun trip a salida:
 * CodEjemplo = Valor Real
 * x0=a;     a=b;    b=c;     c=a(dist input);
 **********************************/
void triplete_inicial(double (*fun)(double, double*, double*), double* x0, double* alpha, double* pk,double* a, double* b, double* c)

    *b= *alpha;                 // *b = *x0;
    *a = *alpha - 0.01;
    *c = *alpha + 0.01;
    double fb = fun(*b,x0,pk);  //alfa, *x0, *pk
    double fa = fun(*a,x0,pk);
    double fc = fun(*c,x0,pk);

    if ((fb < fa) && (fb < fc)) {
        return;
    } else {                    //  Busca hacia abajo desde a hacia b
        double dir = *b - *a;
        if (fa < fb) {          //  Intercambia roles para ir hacia abajo
            double temp = *a;
            *a = *b;
            *b = temp;
            double fx = fa;
            fa = fb;
            fb = fx;
            dir = -dir;
        }

        *b = *a + dir;
        fb = fun(*b,x0,pk);
        *c = *a + dir * 2;
        fc = fun(*c,x0,pk);

        while (fc < fb) {
            *a = *b;
            *b = *c;
            fb = fc;
            *c = *c + dir;
            fc = fun(*c,x0,pk);
        }
    }
}

double rosenb(double* x)
{
    double y = 100 * pow(x[0] - x[1] * x[1], 2) + pow(1 - x[1], 2);
    return y;
}

double f(double alfa, double *x0, double *pk) {
    double temp[2];
    for (int i = 0; i < M; i++) {
        temp[i] = x0[i] + alfa * pk[i];
    }
    return rosenb(temp);
}
void vectorUpdate(double* xini, double* x0, int m){
    for (int i = 0; i < m; i++) {
        xini[i] = x0[i];
    }
}

int main()
{
    double x0[M] = {8.0, 2.0};
    double result[2];
    int iter;

    mini_powell_2(f, x0, M, result, &iter);

    printf("Resultado:\n");
    for (int i = 0; i < 2; i++) {
        printf("%f ", result[i]);
    }
    printf("\n");
    printf("Iteraciones: %d\n", iter);

    return 0;
}