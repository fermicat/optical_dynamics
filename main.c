///@title Main program
///@Auther QCat
///@dev Nt may subject to change

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "physics_constants.h"
#include "laser_dynamics.c"
#include "initialization.c"


/******************* MAIN PROGRAM START HERE ********************/

int main(void)
{
    //cavity parameters:
    double L, rl, rr;                                   //size and reflection coef
    L = 1e-3;   rl = 1; rr = 0.60;                      // L = 3mm, Rl > 95%, Rr^2 = 36%

    // numerical parameters:
    int Nz, Nt;                                         //mesh size TBD
    double dz, dt;                                      //step
    dz = 1e-6; dt = 1e-14; Nt = 700000;                 // dt >= 1e-13: integration divergence
    Nz = (int) (L/dz);

    //laser state
    double *n0;
    double _Complex *Ep, *En, *n2, *etap, *etan;
    double _Complex *Eoutput;
    double _Complex *n2out;
    double *n0out;

    double _Complex *EpC, *EnC, *n2C;
    double *n0C;
    
    //allocate the laser state:
    n0 = (double*) malloc(sizeof(double) * (Nz+1));
    n2 = (double _Complex*) malloc(sizeof(double _Complex) * (Nz+1));
    etap = (double _Complex*) malloc(sizeof(double _Complex) * (Nz+1));
    etan = (double _Complex*) malloc(sizeof(double _Complex) * (Nz+1));
    Ep = (double _Complex*) malloc(sizeof(double _Complex) * (Nz+1));
    En = (double _Complex*) malloc(sizeof(double _Complex) * (Nz+1));
    
    //allocate the output:
    Eoutput = (double _Complex*) malloc(sizeof(double _Complex) * (Nt+1));
    n2out = (double _Complex*) malloc(sizeof(double _Complex) * (Nt+1));
    n0out = (double*) malloc(sizeof(double) * (Nt+1));

    n0C = (double*) malloc(sizeof(double) * (Nz+1));
    n2C = (double _Complex*) malloc(sizeof(double _Complex) * (Nz+1));
    EpC = (double _Complex*) malloc(sizeof(double _Complex) * (Nz+1));
    EnC = (double _Complex*) malloc(sizeof(double _Complex) * (Nz+1));



    //initialization of boundary conditions: laser_sate at t = 0
    field_initialization(Ep, En, rl, rr, Nz);
    population_initialization(n0, n2, Nz);
    eta_initialization(etap, etan, Nz);

    //main caculation
    evolution(Eoutput, Ep, En, n0, n2, etap, etan, Nz, Nt, dz, dt, L, rl, rr, n2out, n0out, n0C, n2C, EpC, EnC);
    
    //print to file
    char userdescribe[40] = "initial n0 = 0.5, E~5e2, n2 = 0";
    //printf("Input the description:");
    //scanf("%s", userdescribe);
    FILE *fp = NULL;
    fp = fopen("./data.dat", "w+");
    fprintf(fp, "#%s\n", userdescribe);
    fprintf(fp, "#Output (z = L) field:\n");
    fprintf(fp, "#time(ps)\t Re[E]\t Im[E] \t |E| \t |E|^2 \t Re[n2] \t n0 \n");
    for(int j = 0; j <= Nt; j++)
    {
        fprintf(fp, "%f\t%f\t%f\t", j*dt/1e-12, creal(Eoutput[j]), cimag(Eoutput[j]));
        double intensity = (double) (conj(Eoutput[j]) * Eoutput[j]);
        fprintf(fp, "%f\t%f\t%f\t%f\n", sqrt(intensity), intensity, creal(n2out[j]), n0out[j] );
    }

    fclose(fp);

    fp = fopen("./cavity.dat", "w+");
    fprintf(fp, "z(10^-6 m) \t n0 \t Re[n2] \t Im[n2] \t Re[Ep] \t Im[Ep] \t Re[En] \t Im[En] \n");
    for (int i = 0; i <= Nz; i++) {
        fprintf(fp, "%d\t%f\t%1.9f\t%f\t", i, n0C[i], creal(n2C[i]), cimag(n2C[i]) );
        fprintf(fp, "%f\t%f\t%f\t%f\n", creal(EpC[i]), cimag(EpC[i]), creal(EnC[i]), cimag(EnC[i]));
    }

    //delocate
    free(n0);
    free(n2);
    free(etap);
    free(etan);
    free(Ep);
    free(En);
    free(Eoutput);

    return 0;
}