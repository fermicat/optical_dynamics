/// @title Subprogram - the laser_dynamics evolution, mainly Maxwell-Bloch evolution
/// @author QCat
/*
 * @dev the functions before evolution() are to caculate the first order time derivative of each laser parameters.
 *Neq(z, h, J) need to change with absorber.
 *In the further program, we need to include "gain_threshold.c" to determine J and Neq
 *****************************
 *The evolution() function is the main subprogram. We use 4th order R-K algorithm to caculate the evolution.
 *The method is f( x(t + dt), t + dt) = f(x(t), t) + K * dt + o(dt^4)
 *where K = (k1 + 2 * k2 + 2 * k3 + k4)/6

 *k1 is the slope predicted at t, i.e. k1 = f'(x(t), t)
 *k2 is the slope predicted at t + dt/2 with k1, i.e. k2 = f'(x(t+dt/2), t+dt/2)
 *k3 is the slope predicted at t + dt/2 with k2, i.e. k3 = f'(x(t+dt/2), t+dt/2)
 *k4 is the slope predicted at t + dt with k3, i.e. k2 = f'(x(t+dt), t+dt)
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "physics_constants.h"
#include "initialization.c"

// laser parameter
// temporory set for testing:
double T1 = 1e-9; double T2 = 50e-15; double Diff = 0.0046;          // T1 = 1ns, T2 = 50fs, Diff = 46cm-2
double overlap = 0.1;  double dipole = 0.8e-9 * echarge;  double nr = 3.5;   double loss = 500.0;
double rho_e = 5.62e22;                                   // for testing 

//Neq is equi-population-inversion
//z = iz*dz, 0<z<h is gain region
double Neq(double z, double h, double J) {
    //no absorber testing:
    return 2 * J - 1;
    /*
    if (z < h)
        return 0.6;
        else
        return 0.1;*/
}

// time derivative dt_n0(n0, etap, etan, Ep, En, i*Nz, h, J)
double dt_n0(double n0, double _Complex etap, double _Complex etan, double _Complex Ep, double _Complex En, double z, double h, double J) {
    double temp = cimag( Ep * conj(etap) + En * conj(etan));
    return -(n0-Neq(z, h, J))/T1 - 2 * dipole / hbar * temp;
}

// time derivative dt_n2(n2, etap, etan, Ep, En, kwave)
double _Complex dt_n2(double _Complex n2, double _Complex etap, double _Complex etan, double _Complex Ep, double _Complex En, double k) {
    double _Complex temp = Ep * conj(etan) - conj(En) * etap;
    return - (1/T1 + 4*k*k*Diff) * n2 + I * dipole / hbar * temp;
}

// time derivative dt_etap(n0, n2, etap, Ep, En)
double _Complex dt_etap(double n0, double _Complex n2, double _Complex etap, double _Complex Ep, double _Complex En) {
    double _Complex temp = n0 * Ep + n2 * En;
    return - I * dipole / (2*hbar) * temp - etap/T2;
}

// time derivative dt_etap(n0, n2, etap, Ep, En)
double _Complex dt_etan(double n0, double _Complex n2, double _Complex etan, double _Complex Ep, double _Complex En) {
    double _Complex temp = n0 * En + conj(n2) * Ep;
    return - I * dipole / (2*hbar) * temp - etan/T2;
}

// time derivative dt_Ep(etap, Ep, dzEp, G)
double _Complex dt_Ep(double _Complex etap, double _Complex Ep, double _Complex dzE, double _G) {
    return (c0/nr) * (I * _G * etap - loss * Ep - dzE);
}

// time derivative dt_En(etan, En, dzEn, G)
double _Complex dt_En(double _Complex etan, double _Complex En, double _Complex dzE, double _G) {
    return (c0/nr) * (I * _G * etan - loss * En + dzE);
}


/************************ MAJOR FUNCTION START HERE *********************/

// Dynamic evolution
void evolution(double _Complex* Eout, double _Complex* Ep, double _Complex* En, double* n0, double _Complex* n2, double _Complex* etap, double _Complex* etan, int Nz, int Nt, double dz, double dt, double L, double _rl, double _rr, double _Complex *n2out, double *n0out, double *n0C, double _Complex *n2C, double _Complex *EpC, double _Complex *EnC)
{
    // 1D temporary sequence in this procedure
    double n0temp[Nz];
    double _Complex n2temp[Nz], Eptemp[Nz], Entemp[Nz]; 
    double _Complex etaptemp[Nz], etantemp[Nz];
    double _Complex dzEp[Nz], dzEn[Nz];
    
    for (int i = 0; i <= Nz; i++)                       //initial for temp
    {
        n0temp[i] = n0[i];
        n2temp[i] = n2[i];
        etaptemp[i] = etap[i];
        etantemp[i] = etan[i];
        Eptemp[i] = Ep[i];
        Entemp[i] = En[i];
    }

    //store the 1st derivative for R-K algorithm as Euler ajust
    double k0[5][Nz+1];
    double _Complex k2[5][Nz+1], kap[5][Nz+1], kan[5][Nz+1], kp[5][Nz+1], kn[5][Nz+1];

    //cavity parameters:
    double h = 3e-4;
    double J = 0.8;
    double gain_th = loss - 1/L * log(_rl*_rr);
    
    //physical quantities in lasing process: (to save caculation resource)
    double lambda0 = 3.2e-6;
    double Etransition = hbar * c0 * 2.0 * PI / lambda0; 
    double kwave = nr * Etransition / (c0 * hbar);
    double G = overlap *dipole * (Etransition/hbar) * rho_e / (nr * epsilon0 * c0);
    double gain_coef = overlap * dipole * dipole * (Etransition/hbar) * T2 * rho_e / (2.0*hbar*epsilon0*nr*c0);

    // first output field
    Eout[0] = (1 - _rr) * Ep[Nz];
    
    //Using 4th-order Runge-Kutta Algorithm to caculate the laser state for each time point
    //laser state ~ o(dt^4)
    for (int j = 0; j < Nt; j++) {
        
        printf("t:j=%d\t", j); //debug

        //caculating k1 for each z = i*dz
        grad_caculate(dzEp, Eptemp, Nz, dz);
        grad_caculate(dzEn, Entemp, Nz, dz);
        for (int i = 0; i <= Nz; i++) {
            k0[1][i] = dt_n0(n0temp[i], etaptemp[i], etantemp[i], Eptemp[i], Entemp[i], i*dz, h, J);
            k2[1][i] = dt_n2(n2temp[i], etaptemp[i], etantemp[i], Eptemp[i], Entemp[i], kwave);
            kap[1][i] = dt_etap(n0temp[i], n2temp[i], etaptemp[i], Eptemp[i], Entemp[i]);
            kan[1][i] = dt_etan(n0temp[i], n2temp[i], etantemp[i], Eptemp[i], Entemp[i]);;
            kp[1][i] = dt_Ep(etaptemp[i], Eptemp[i], dzEp[i], G);
            kn[1][i] = dt_En(etantemp[i], Entemp[i], dzEn[i], G);
            
            // temp here means z = i*dz, t = j*dt + dt/2 (pre) for computing k2
            n0temp[i] = n0[i] + k0[1][i]*dt/2;
            n2temp[i] = n2[i] + k2[1][i]*dt/2;
            //n2temp[i] = 0;
            etaptemp[i] = etap[i] +kap[1][i]*dt/2;
            etantemp[i] = etan[i] +kan[1][i]*dt/2;
            Eptemp[i] = Ep[i] + kp[1][i]*dt/2;
            Entemp[i] = En[i] + kn[1][i]*dt/2;
        }
    
        //debug
        //int tk0 = (int)(k0[1][Nz] / 1e7); 
        //int tk2 = (int)(creal(k2[1][Nz]) / 1e5);
        //int tkp = (int)(creal(kp[1][Nz])/ 1e10);
        //int tkn = (int)(creal(kn[1][Nz])/ 1e10);
        //int tdzE = (int) (creal(dzEp[Nz]) / 1e7);
        //printf("k0=%de7,\tk2=%de5,\tkp=%de10,\tkn=%de10\tdE/dz=%de7\n", tk0, tk2, tkp, tkn,tdzE);

        //caculating k2 for each z = i*dz
        grad_caculate(dzEp, Eptemp, Nz, dz);
        grad_caculate(dzEn, Entemp, Nz, dz);
        for (int i = 0; i <= Nz; i++) {
            k0[2][i] = dt_n0(n0temp[i], etaptemp[i], etantemp[i], Eptemp[i], Entemp[i], i*dz, h, J);
            k2[2][i] = dt_n2(n2temp[i], etaptemp[i], etantemp[i], Eptemp[i], Entemp[i], kwave);
            kap[2][i] = dt_etap(n0temp[i], n2temp[i], etaptemp[i], Eptemp[i], Entemp[i]);
            kan[2][i] = dt_etan(n0temp[i], n2temp[i], etantemp[i], Eptemp[i], Entemp[i]);;
            kp[2][i] = dt_Ep(etaptemp[i], Eptemp[i], dzEp[i], G);
            kn[2][i] = dt_En(etantemp[i], Entemp[i], dzEn[i], G);
            
            // temp here means z = i*dz, t = j*dt + dt/2 (post) for computing k3
            n0temp[i] = n0[i] + k0[2][i]*dt/2;
            n2temp[i] = n2[i] + k2[2][i]*dt/2;
            //n2temp[i] = 0;
            etaptemp[i] = etap[i] +kap[2][i]*dt/2;
            etantemp[i] = etan[i] +kan[2][i]*dt/2;
            Eptemp[i] = Ep[i] + kp[2][i]*dt/2;
            Entemp[i] = En[i] + kn[2][i]*dt/2;
        }

        //caculating k3 for each z = i*dz
        grad_caculate(dzEp, Eptemp, Nz, dz);
        grad_caculate(dzEn, Entemp, Nz, dz);
        for (int i = 0; i <= Nz; i++) {
            k0[3][i] = dt_n0(n0temp[i], etaptemp[i], etantemp[i], Eptemp[i], Entemp[i], i*dz, h, J);
            k2[3][i] = dt_n2(n2temp[i], etaptemp[i], etantemp[i], Eptemp[i], Entemp[i], kwave);
            kap[3][i] = dt_etap(n0temp[i], n2temp[i], etaptemp[i], Eptemp[i], Entemp[i]);
            kan[3][i] = dt_etan(n0temp[i], n2temp[i], etantemp[i], Eptemp[i], Entemp[i]);;
            kp[3][i] = dt_Ep(etaptemp[i], Eptemp[i], dzEp[i], G);
            kn[3][i] = dt_En(etantemp[i], Entemp[i], dzEn[i], G);
            
            // temp here means z = i*dz, t = j*dt + dt for computing k4 
            n0temp[i] = n0[i] + k0[3][i]*dt;
            n2temp[i] = n2[i] + k2[3][i]*dt;
            //n2temp[i] = 0;
            etaptemp[i] = etap[i] +kap[3][i]*dt;
            etantemp[i] = etan[i] +kan[3][i]*dt;
            Eptemp[i] = Ep[i] + kp[3][i]*dt;
            Entemp[i] = En[i] + kn[3][i]*dt;
        }
        
        //caculating k4 for each z = i*dz
        grad_caculate(dzEp, Eptemp, Nz, dz);
        grad_caculate(dzEn, Entemp, Nz, dz);
        for (int i = 0; i <= Nz; i++) {
            k0[4][i] = dt_n0(n0temp[i], etaptemp[i], etantemp[i], Eptemp[i], Entemp[i], i*dz, h, J);
            k2[4][i] = dt_n2(n2temp[i], etaptemp[i], etantemp[i], Eptemp[i], Entemp[i], kwave);
            kap[4][i] = dt_etap(n0temp[i], n2temp[i], etaptemp[i], Eptemp[i], Entemp[i]);
            kan[4][i] = dt_etan(n0temp[i], n2temp[i], etantemp[i], Eptemp[i], Entemp[i]);;
            kp[4][i] = dt_Ep(etaptemp[i], Eptemp[i], dzEp[i], G);
            kn[4][i] = dt_En(etantemp[i], Entemp[i], dzEn[i], G);
        }
        
        //caculating the final slope for each z = i*dz
        for (int i = 0; i <= Nz; i++)
        {
            k0[0][i] = k0[1][i] + 2 * k0[2][i] + 2 * k0[3][i] + k0[4][i];
            k2[0][i] = k2[1][i] + 2 * k2[2][i] + 2 * k2[3][i] + k2[4][i];
            kap[0][i] = kap[1][i] + 2 * kap[2][i] + 2 * kap[3][i] + kap[4][i];
            kan[0][i] = kan[1][i] + 2 * kan[2][i] + 2 * kan[3][i] + kan[4][i];
            kn[0][i] = kn[1][i] + 2 * kn[2][i] + 2 * kn[3][i] + kn[4][i];
            kp[0][i] = kp[1][i] + 2 * kp[2][i] + 2 * kp[3][i] + kp[4][i];
        }

        // caculate the laser state of time (j+1)dt for each z =i*dz
        for (int i = 0; i <= Nz; i++) {
        
            n0[i] += k0[0][i]*dt/6;                         // [(j+1)dt] = [jdt] + kdt
            n2[i] += k2[0][i]*dt/6;
            etap[i] += kap[0][i]*dt/6;
            etan[i] += kan[0][i]*dt/6;
            Ep[i] += kp[0][i]*dt/6;
            En[i] += kn[0][i]*dt/6;
            
            n0temp[i] = n0[i];
            n2temp[i] = n2[i];
            //n2temp[i] = 0;
            etaptemp[i] = etap[i];
            etantemp[i] = etan[i];
            Eptemp[i] = Ep[i];
            Entemp[i] = En[i];
        }

        //reflection correction
        Ep[0] = _rl * En[0];                                 
        En[Nz] = _rr * Ep[Nz];
        Eptemp[0] = Ep[0];
        Entemp[Nz] = En[Nz];

        Eout[j+1] = (1 - _rr) * Ep[Nz];                      //output field
        n0out[j+1] = n0[Nz];
        n2out[j+1] = n2[Nz];
            printf("\tEout[t] = %f\tn0 = %f\t n2 = %f\n", creal(Eout[j+1]),n0[Nz], creal(n2[Nz]));
    }

    for (int i = 0; i <= Nz; i++) {
        n0C[i] = n0[i];
        n2C[i] = n2[i];
        EpC[i] = Ep[i];
        EnC[i] = En[i];
    }

}