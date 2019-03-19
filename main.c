#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "physics_constants.h"

/// @dev change ii -> I
/// @dev complex<double> -> double_Complex
/// @dev slope
/// @dev debug print

/************************************** parameters ******************************************/
    //cavity parameters:
    double L  = 1e-3, rl = 1, rr = 0.60;                      // L = 3mm, Rl > 95%, Rr^2 = 36%
    double dz = 1e-6;                   // dt >= 1e-13: integration divergence
    int    Nz = (int) (L/dz);
    double h  = 3e-4, J  = 0.8;

    // laser parameter
    double T1      = 1e-12,     T2     = 50e-15,    Diff = 1;          // T1 = 1ns, T2 = 50fs, Diff = 46cm2s^-1
    double overlap = 0.1,       dipole = 0.8e-9 * echarge;
    double nr      = 3.5,       loss    = 500.0;
    double rho_e   = 5.62e22;                                   // for testing 
    double lambda0 = 3.2e-6;
    double Etransition = hbar * c0 * 2.0 * PI / lambda0; 
    double kwave       = nr * Etransition / (c0 * hbar);
    double G = overlap * dipole * (Etransition/hbar) * rho_e / (nr * epsilon0 * c0);

    // time parameter
    int    Nt = 500;    
    double dt = 1e-14;

/************************************ laser state ********************************************/
    double *n0;
    double _Complex *Ep, *En, *n2, *etap, *etan, *dzEp, *dzEn;
    // slope
    double *n0_old;
    complex<double> *Ep_old, *En_old, *n2_old, *etap_old, *etan_old;

    // copy


/****************************************** PDE **********************************************/
 //Neq is equi-population-inversion
    //z = iz*dz, 0<z<h is gain region
    double Neq(double z) {
        return 2 * J - 1;
        //if (z < h) return 0.6;
        //else return 0.1; 
    }

    /**
     *  six dynamic equation 
     *  Schordinger Maxwell equation
     */
    void dt_n0(int round) {
        for (int j = 0; j < Nz; j++) {
            double z = j * dz;
            double temp = (Ep[j] * conj(etap[j]) + En[j] * conj(etan[j])).imag();
            k0[round][j] = - (n0[j] - Neq(z)) / T1 - 2.0 * dipole / hbar * temp;
        } 
    }

    void dt_n2(int round) {
        for (int j = 0; j < Nz; j++) {
            complex<double> temp = Ep[j] * conj(etan[j]) - conj(En[j]) * etap[j];
            k2[round][j] = - (1 / T1 + 4 * kwave * kwave * Diff) * n2[j] + ii * dipole / hbar * temp;
        }
    }

    void dt_etap(int round) {
        for (int j = 0; j < Nz; j++) {
            complex<double> temp = n0[j] * Ep[j] + n2[j] * En[j];
            kap[round][j] = - ii * dipole / (2*hbar) * temp - etap[j] / T2;
        }
    }

    void dt_etan(int round) {
        for (int j = 0; j < Nz; j++) {
            complex<double> temp = n0[j] * En[j] + conj(n2[j]) * Ep[j];
            kan[round][j] - ii * dipole / (2*hbar) * temp - etan[j] / T2;
        }
    }
  
    void dt_Ep(int round) {
        for (int j = 0; j < Nz; j++) {
            kp[round][j] = (c0/nr) * (ii * G * etap[j] - loss * Ep[j] - dzEp[j]);
        }
    }

    void dt_En(int round) {
        for (int j = 0; j < Nz; j++) {
            kn[round][j] = (c0/nr) * (ii * G * etan[j] - loss * En[j] + dzEn[j]);
        }
    }
    
    // dE/dz 
    void gradE(complex<double> dzE[], complex<double> E[]) {
        dzE[0] = (E[1] - E[0]) / dz;
        dzE[Nz] = (E[Nz] - E[Nz-1]) / dz;
        for (int j = 1; j < Nz - 1; j++) {
            dzE[j] = (E[j + 1] - E[j - 1]) / (2*dz);
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
    void field_initialization(complex<double> Epl, complex<double> Epr) {
        complex<double> Enl, Enr;
        Enr = Epr * rr;                                     //reflection
        Enl = Epl / rl;
        complex<double> kEp = (Epr - Epl) / double(Nz - 1);
        complex<double> kEn = (Enr - Enl) / double(Nz - 1);
        for (int i = 0; i < Nz; i++) {
            Ep[i] = Epl + double(i) * kEp;
            En[i] = Enl + double(i) * kEn;
        }
    }

    //initialize the n(z, t=0)
    void population_initialization(double _n0, complex<double> _n2) {
        for (int i = 0; i < Nz; i++) {
            n0[i] = _n0;       
            n2[i] = _n2;
        }
    }

    void eta_initialization(complex<double> _etap, complex<double> _etan) {
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
        for (int i = 0; i < Nt; i++) {
            copyArray();
            Eout[i] = (1 - rr) * Ep[Nz];                      //output field

            // Runge-Kutta
            for (int round = 1; round <= 4; round++) {
                gradE(dzEp, Ep);
                gradE(dzEn, En);
                evolution(round);
                rungeKutta(round);
            }
            
            Ep[0] = rl * En[0];                                 
            En[Nz] = rr * Ep[Nz];
        
            //debugPrint(i);
        }
    }



/**************************************** Main ********************************/

int main() {
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

    Eout = new complex<double>[Nt];

    // new for slope


    // initialization
    complex<double> Epl = 5e2 + ii;
    complex<double> Epr = 2e2 + 2.0 * ii;
    field_initialization(Epl, Epr);
    
    double _n0 = 0.56;  
    complex<double> _n2 = 0.0;
    population_initialization(_n0, _n2);
    
    complex<double> _etap = 0.001, _etan = 0.001;
    eta_initialization(_etap, _etan);

    // dynamics
    model();

    FILE *fp = NULL;
    fp = fopen("./data.dat", "w+");
    
    fclose(fp);

    //free

}