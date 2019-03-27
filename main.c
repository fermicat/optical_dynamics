#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "physics_constants.h"


/************************************** parameters ******************************************/
    //cavity parameters:
    double L  = 3e-3, rl = 1, rr = 0.60;                      // L = 3mm, Rl > 95%, Rr^2 = 36%
    double dz = 1e-6;                   
    int    Nz;
    double h  = 3e-4, J  = 0.8;

    // laser parameter
    double T2     = 50e-15,    Diff = 0.0046;          // T1 = 1ns, T2 = 50fs, Diff = 46cm2s^-1
    double overlap = 0.1,       dipole = 0.8e-9 * echarge;
    double nr      = 3.5,       loss   = 500.0;
    double rho_e   = 4.34e22;                                   // for testing 
    double lambda0 = 3.2e-6;
    double kwave, G;

    // time parameter
    int    Nt;    
    double dt = 1e-14;              // dt > 1e-13 integral divergence

/************************************ laser state ********************************************/
    double *n0;
    double _Complex *Ep, *En, *n2, *etap, *etan, *dzEp, *dzEn;
    // slope
    double **k0;
    double _Complex **k2, **kap, **kan, **kp, **kn;

    //copy
    double *n0_old;
    double _Complex *Ep_old, *En_old, *n2_old, *etap_old, *etan_old, *dzEp_old, *dzEn_old;
    
    double _Complex *Eout;

/****************************************** PDE **********************************************/
 //Neq is equi-population-inversion
    //z = iz*dz, 0<z<h is gain region
    double Neq(double z) {
        //return 2 * J - 1;
        if (z < L - h) return 0.6;
        else return 0.2; 
    }

    double T1(double z) {
        if (z < L - h) return 1e-12;
        else return 1e-13;
    }

    /**
     *  six dynamic equation 
     *  Schordinger Maxwell equation
     */
    void dt_n0(int round) {
        for (int j = 0; j < Nz; j++) {
            double z = j * dz;
            double temp = cimag( (Ep[j] * conj(etap[j]) + En[j] * conj(etan[j])) );
            k0[round][j] = - (n0[j] - Neq(z)) / T1(z) - 2.0 * dipole / hbar * temp;
        } 
    }

    void dt_n2(int round) {
        for (int j = 0; j < Nz; j++) {
            complex temp = Ep[j] * conj(etan[j]) - conj(En[j]) * etap[j];
            k2[round][j] = - (1 / T1(j * dz) + 4 * kwave * kwave * Diff) * n2[j] + I * dipole / hbar * temp;
        }
    }

    void dt_etap(int round) {
        for (int j = 0; j < Nz; j++) {
            complex temp = n0[j] * Ep[j] + n2[j] * En[j];
            kap[round][j] = - I * dipole / (2*hbar) * temp - etap[j] / T2;
        }
    }

    void dt_etan(int round) {
        for (int j = 0; j < Nz; j++) {
            complex temp = n0[j] * En[j] + conj(n2[j]) * Ep[j];
            kan[round][j] = - I * dipole / (2*hbar) * temp - etan[j] / T2;
        }
    }
  
    void dt_Ep(int round) {
        for (int j = 0; j < Nz; j++) {
            kp[round][j] = (c0/nr) * (I * G * etap[j] - loss * Ep[j] - dzEp[j]);
        }
    }

    void dt_En(int round) {
        for (int j = 0; j < Nz; j++) {
            kn[round][j] = (c0/nr) * (I * G * etan[j] - loss * En[j] + dzEn[j]);
        }
    }
    
    // dE/dz 
    void gradE(complex *dzE, complex *E) {
        dzE[0] = (E[1] - E[0]) / dz;
        dzE[Nz - 1] = (E[Nz - 1] - E[Nz - 2]) / dz;
        for (int j = 1; j < Nz - 1; j++) {
            dzE[j] = (E[j + 1] - E[j - 1]) / (2 * dz);
        }
    }


/********************************* Runge - Kutta *****************************************/

    void copyArray() {
        for (int j = 0; j < Nz; j++) {
            n0_old[j] = n0[j];      
            n2_old[j] = n2[j];
            Ep_old[j] = Ep[j];     
            En_old[j] = En[j];
            etap_old[j] = etap[j];  
            etan_old[j] = etan[j];
        }
    }

    void updateTempState(int round, double dt) {
        for (int i = 0; i < Nz; i++) {
            n0[i]   = n0_old[i]   + k0[round][i]  * dt;
            n2[i]   = n2_old[i]   + k2[round][i]  * dt;
            //n2[i] = 0;
            etap[i] = etap_old[i] + kap[round][i] * dt;
            etan[i] = etan_old[i] + kan[round][i] * dt;
            Ep[i]   = Ep_old[i]   + kp[round][i]  * dt;
            En[i]   = En_old[i]   + kn[round][i]  * dt;
        }
    }

    void updateState() {
        for (int i = 0; i < Nz; i++) {
            k0[0][i]  = k0[1][i]  + 2.0 * k0[2][i]  + 2.0 * k0[3][i]  + k0[4][i];
            k2[0][i]  = k2[1][i]  + 2.0 * k2[2][i]  + 2.0 * k2[3][i]  + k2[4][i];
            //k2[0][i] = 0;
            kap[0][i] = kap[1][i] + 2.0 * kap[2][i] + 2.0 * kap[3][i] + kap[4][i];
            kan[0][i] = kan[1][i] + 2.0 * kan[2][i] + 2.0 * kan[3][i] + kan[4][i];
            kn[0][i]  = kn[1][i]  + 2.0 * kn[2][i]  + 2.0 * kn[3][i]  + kn[4][i];
            kp[0][i]  = kp[1][i]  + 2.0 * kp[2][i]  + 2.0 * kp[3][i]  + kp[4][i];
        }

        for (int i = 0; i < Nz; i++) {
            n0[i]   = n0_old[i]   + k0[0][i]  * dt/6.0;     // [(j+1)dt] = [jdt] + kdt
            n2[i]   = n2_old[i]   + k2[0][i]  * dt/6.0;
            etap[i] = etap_old[i] + kap[0][i] * dt/6.0;
            etan[i] = etan_old[i] + kan[0][i] * dt/6.0;
            Ep[i]   = Ep_old[i]   + kp[0][i]  * dt/6.0;
            En[i]   = En_old[i]   + kn[0][i]  * dt/6.0;
        }
    }

    // update value
    void rungeKutta(int round) {
        if (round == 1 || round == 2) {
            updateTempState(round, dt / 2.0);
        }
        else if (round == 3) {
            updateTempState(round, dt);
        }
        else if (round == 4) {
            updateState();
        }
        else {
            printf("error!\n");
        }
    }


/********************************** State Initialization *********************************/

    // initialize the E(z, t=0)
    void field_initialization(complex Epl, complex Epr) {
        complex Enl, Enr;
        Enr = Epr * rr;                                     //reflection
        Enl = Epl / rl;
        complex kEp = (Epr - Epl) / (Nz - 1.0);
        complex kEn = (Enr - Enl) / (Nz - 1.0);
        for (int i = 0; i < Nz; i++) {
            Ep[i] = Epl + (double)(i) * kEp;
            En[i] = Enl + (double)(i) * kEn;
        }
    }

    //initialize the n(z, t=0)
    void population_initialization(double _n0, complex _n2) {
        for (int i = 0; i < Nz; i++) {
            n0[i] = _n0;       
            n2[i] = _n2;
        }
    }

    void eta_initialization(complex _etap, complex _etan) {
        for (int i = 0; i < Nz; i++) {
            etap[i] = _etap;
            etan[i] = _etan;
        }
    }


/*********************************** Dynamics **************************************/
    void evolution(int round) {
        dt_n0(round);   
        dt_n2(round);
        dt_etap(round); 
        dt_etan(round);
        dt_Ep(round);   
        dt_En(round);
    }

    void model() {
        void showSteps(int t);

        for (int i = 0; i < Nt; i++) {
            copyArray();
            Eout[i] = (1 - rr) * Ep[Nz - 1];                      //output field

            // Runge-Kutta
            for (int round = 1; round <= 4; round++) {
                gradE(dzEp, Ep);
                gradE(dzEn, En);
                evolution(round);
                rungeKutta(round);
            }
            
            Ep[0] = rl * En[0];                                 
            En[Nz - 1] = rr * Ep[Nz - 1];
        
            showSteps(i);
        }
    }

/**************************** system & malloc & free ****************************/
    void showSteps(int t) {
        if (t % 500 == 0) {
            printf("t = %d\n", t);
            printf("Eout = %f, \tn0 = %f, \tn2 = %f\n", creal(Eout[t]), n0[Nz - 1], creal(n2[Nz - 1]));
        }
    }

    void creatArray() {
        // new for laser state
        n0   = (double*) malloc(sizeof(double) * Nz);
        n2   = (double _Complex*) malloc(sizeof(double _Complex) * Nz);
        etap = (double _Complex*) malloc(sizeof(double _Complex) * Nz);
        etan = (double _Complex*) malloc(sizeof(double _Complex) * Nz);
        Ep   = (double _Complex*) malloc(sizeof(double _Complex) * Nz);
        En   = (double _Complex*) malloc(sizeof(double _Complex) * Nz);
        dzEp = (double _Complex*) malloc(sizeof(double _Complex) * Nz);
        dzEn = (double _Complex*) malloc(sizeof(double _Complex) * Nz);

        // new for copy state
        n0_old   = (double*) malloc(sizeof(double) * Nz);
        n2_old   = (double _Complex*) malloc(sizeof(double _Complex) * Nz);
        etap_old = (double _Complex*) malloc(sizeof(double _Complex) * Nz);
        etan_old = (double _Complex*) malloc(sizeof(double _Complex) * Nz);
        Ep_old   = (double _Complex*) malloc(sizeof(double _Complex) * Nz);
        En_old   = (double _Complex*) malloc(sizeof(double _Complex) * Nz);
        dzEp_old = (double _Complex*) malloc(sizeof(double _Complex) * Nz);
        dzEn_old = (double _Complex*) malloc(sizeof(double _Complex) * Nz);

        Eout     = (complex*) malloc(sizeof(complex) * Nt);

        // new for slope
        k0  = (double**)  malloc(sizeof(double*) * 5);
        k2  = (complex**) malloc(sizeof(complex*) * 5);
        kap = (complex**) malloc(sizeof(complex*) * 5);
        kan = (complex**) malloc(sizeof(complex*) * 5);
        kp  = (complex**) malloc(sizeof(complex*) * 5);
        kn  = (complex**) malloc(sizeof(complex*) * 5);

        for (int i = 0; i < 5; i++) {
            k0[i]  = (double*)  malloc(sizeof(double) * Nz);
            k2[i]  = (complex*) malloc(sizeof(complex) * Nz);
            kap[i] = (complex*) malloc(sizeof(complex) * Nz);
            kan[i] = (complex*) malloc(sizeof(complex) * Nz);
            kp[i]  = (complex*) malloc(sizeof(complex) * Nz);
            kn[i]  = (complex*) malloc(sizeof(complex) * Nz);
        }

    }

    void deleteArray() {
        free(n0);
        free(n2);
        free(etap);
        free(etan);
        free(Ep);
        free(En);
        //free(dzEp);
        free(dzEn);
        
        free(n0_old);
        free(n2_old);
        free(etap_old);
        free(etan_old);
        
        free(Ep_old);
        free(En_old);
        //free(dzEp_old);
        free(dzEn_old);
        
        free(Eout);
        
        for (int i = 0; i < 5; i++) {
            free(k0[i]);
            free(k2[i]);
            free(kap[i]);
            free(kan[i]);
            free(kp[i]);
            free(kn[i]);
        }
        free(k0);
        free(k2);
        free(kap);
        free(kan);
        free(kp);
        free(kn);
    }


/**************************************** Main ********************************/

int main() {
    // quantity
    double Etransition = hbar * c0 * 2.0 * PI / lambda0; 
    kwave = nr * Etransition / (c0 * hbar);
    G     = overlap * dipole * (Etransition/hbar) * rho_e / (nr * epsilon0 * c0);
    Nz    = (int) (L / dz);
    Nt    = 500000;

    // malloc resource
    creatArray();
    
    // initialization
    complex Epl = 5e4 + I, Epr = 2e4 + 2.0 * I;
    field_initialization(Epl, Epr);
    
    double  _n0 = 0.56;  
    complex _n2 = 0.0;
    population_initialization(_n0, _n2);
    
    complex _etap = 0.001, _etan = 0.001;
    eta_initialization(_etap, _etan);

    // dynamics
    model();

    // file
    FILE *fp = NULL;
    fp = fopen("./data.dat", "w+");
    
    fprintf(fp, "#time(ps)\t Re[E]\t Im[E] \t |E| \t |E|^2 \n");
    for(int j = 0; j <= Nt; j++)
    {
        fprintf(fp, "%f\t%f\t%f\t", j*dt/1e-12, creal(Eout[j]), cimag(Eout[j]));
        double intensity = (double) (conj(Eout[j]) * Eout[j]);
        fprintf(fp, "%f\t%f\n", sqrt(intensity), intensity );
    }
    
    fclose(fp);

    //free
    deleteArray();
    return 0;
}