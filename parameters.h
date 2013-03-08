#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>

const double hbarC = 0.197327053;  //GeV*fm

// physical parameters for prefactors of emission rates 
const double e_sq = 1./137.*4.*M_PI;
const double q_sq = pow(2./3., 2) + pow(-1./3., 2) + pow(-1./3., 2); //u, d, s
const double d_F = 3; //number of flavor
const double N_c = 3; //number of color
const double C_F = (N_c*N_c - 1)/2/N_c;
const double g_s = 2.0; // strong coupling constant

// parameters for non-thermal equilibrium delta f correction
const double deltaf_alpha = 2;   //the exponent of the p dependence of the delta f correction deltaf_alpha = 2 for democratic choice

// parameters for numerical integrations
const double q_cutoff = 0.0;
const int n_s = 100;
const double s_max = 50.0;
const int n_t = 20;

const int n_E1 = 20;
const int n_E2 = 20;

// parameters for output photon emission rate table
const int n_Eq = 80;
const double Eq_i = 0.05;
const double dEq = 0.05;
const int n_Temp = 1;
const double T_i = 0.3;
const double dT = 0.2;

#endif
