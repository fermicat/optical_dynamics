///@title functions - initialization, delocation and dE/dz caculation
///@author QCat
///@dev 2D delocation is no use in current program

#include <stdlib.h>
#include <complex.h>

#ifndef INITIALIZATION_C
#define INITIALIZATION_C

// initialize the E(z, t=0)
void field_initialization(double _Complex Ep[], double _Complex En[], double _rl, double _rr, int _Nz) {
    double _Complex Epl, Epr, Enl, Enr;
    double kp, kn;                                      //slope
    Epl = 5e2 + 1*I;                                  //initial E
    Epr = 2e2 + 2*I;
    Enr = Epr * _rr;                                     //reflection
    Enl = Epl / _rl;
    kp = (Epr - Epl) / _Nz;
    kn = (Enr - Enl) / _Nz;
    for (int i = 0; i <= _Nz; i++) {
        Ep[i] = Epl + i * kp;
        En[i] = Enl + i * kn;
    }
}

//initialize the n(z, t=0)
void population_initialization(double n0[], double _Complex n2[], int _Nz) {
    for (int i = 0; i <= _Nz; i++) {
        n0[i] = +0.45;                                    //TBD
        n2[i] = 0;
        //n2[i] = +0.2 + 0.2*I;
    }
}

void eta_initialization(double _Complex etap[], double _Complex etan[], int _Nz) {
    for (int i = 0; i <= _Nz; i++) {
        etap[i] = 0.01;
        etan[i] = 0.01;
    }
}

// caculate dE/dz
void grad_caculate(double _Complex dzE[], double _Complex E[], int _Nz, double _dz) {
    dzE[0] = (E[1] - E[0]) / _dz;
    dzE[_Nz] = (E[_Nz] - E[_Nz-1]) / _dz;
    for (int i = 1; i < _Nz; i++)
        dzE[i] = (E[i+1] - E[i-1]) / (2*_dz);
}

//delocate 2D Array
void delocate(double** Array, int n) {
    for(int i = 0; i <= n; i++)
        free(Array[i]);
    free(Array);
}

void delocateC(double _Complex** Array, int n) {
    for(int i = 0; i <= n; i++)
        free(Array[i]);
    free(Array);
}

#endif