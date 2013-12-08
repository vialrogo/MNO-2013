/*
 * Victor Alberto Romero
 * Métodos numéricos em otimização
 * 12 de Setembro de 2013
 * IME - USP
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*
 * Global variables
*/
double M = 2.0;
double gama = 0.0001;
double epsilon = 0.0000000001; /* 10^(-10) its practicaly 0 */

/*
 * Objective function: Rosenbrock function
 * f(x1,x2) = 100(x2-x1^2)^2 + (1-x1)^2
 * x* = (1,1)^T is the only local minimizator
*/
void Rosenbrock(double* xk, double* output)
{
    *output = 100*(xk[1]-xk[0]*xk[0])*(xk[1]-xk[0]*xk[0]) + (1-xk[0])*(1-xk[0]);
}

/*
 * First derivate of the Rosenbrock function in R2
 * f'(x1,x2) = [400*x1^3 -400*x2*x1 +2*x1 - 2 ; 200(x2 - x1^2) ]
*/
void RosenbrockD(double* xk, double* output)
{
    output[0] = 400*xk[0]*xk[0]*xk[0] - 400*xk[1]*xk[0] + 2*xk[0] -2;
    output[1] = 200*(xk[1]-xk[0]*xk[0]);
}

/*
 * Function that return the minimum valor with linear search
*/
void returnMin(double* xk)
{
    double* d     = (double*)malloc(2*sizeof(double));
    double* fD    = (double*)malloc(2*sizeof(double));
    double* xkNew = (double*)malloc(2*sizeof(double));
    double f, fNew, alfa;
    int iterations = 0;

    while(iterations < 1000000000) /* Stop by iterations */
    {
        /* Basic parameters */
        alfa = 1;
        Rosenbrock(xk, &f);
        RosenbrockD(xk, fD);

        /* Calculate the gradient and multiply by -1 */
        d[0]=(-1)*fD[0];
        d[1]=(-1)*fD[1];

        /* Check stop conditions: The two norm */
        if( sqrt(d[0]*d[0] + d[1]*d[1]) < epsilon )
            break;

        /* Calculate a possible new xk */
        xkNew[0]= xk[0]+alfa*d[0];
        xkNew[1]= xk[1]+alfa*d[1];
        Rosenbrock(xkNew, &fNew);

        /* Calculate a new alfa */
        while( fNew > ( f + gama*alfa*(fD[0]*d[0] + fD[1]*d[1])  ))
        {            
            alfa = alfa / M;
            xkNew[0]= xk[0]+alfa*d[0];
            xkNew[1]= xk[1]+alfa*d[1];
            Rosenbrock(xkNew, &fNew);
        }

        /* Update xk */
        xk[0]=xkNew[0];
        xk[1]=xkNew[1];

        iterations++;
    }

    /* Print the result */
    printf("(%f,%f) %d\n",xk[0],xk[1],iterations);

    /* Memory free */
    free(fD);
    free(d);
    free(xkNew);
}

int main(void)
{
    double i,j;
    double limit = 2;
    double step = 1;
    double* x0 = (double*)malloc(2*sizeof(double));

    x0[0]=1.2;
    x0[1]=1.2;
    printf("x0=(%f,%f) ",x0[0],x0[1]);
    returnMin(x0);

    x0[0]=-1.2;
    x0[1]=1;
    printf("x0=(%f,%f) ",x0[0],x0[1]);
    returnMin(x0);


    for (i = -limit; i <= limit; i+=step)
    {
        for (j = -limit; j <= limit; j+=step)
        {
            x0[0]=i;
            x0[1]=j;
            printf("(%f,%f) ",x0[0],x0[1]);
            returnMin(x0);
        }
    }

    x0[0]=-47.5;
    x0[1]=-39;
    printf("x0=(%f,%f)  ",x0[0],x0[1]);
    returnMin(x0);

    x0[0]=-48;
    x0[1]=-40;
    printf("x0=(%f,%f)  ",x0[0],x0[1]);
    returnMin(x0);

    x0[0]=-46;
    x0[1]=-37;
    printf("x0=(%f,%f)  ",x0[0],x0[1]);
    returnMin(x0);

    return 0;
}
