#include <iostream>
#include <cmath>

#include "Matrix_elements_sq.h"
using namespace std;

void Matrix_elements_sq(int channel, double s, double t, double E1, double E2, double Eq, double T, double* result_ptr)
{
   switch(channel)
   {
      case 1:  //Compton scattering
         Matrix_elements_sq_Compton(s, t, E1, E2, Eq, T, result_ptr);
         break;
      case 2:  //pair annihilation
         Matrix_elements_sq_Annihilation(s, t, E1, E2, Eq, T, result_ptr);
         break;
      default:
         cout << "Matrix elements squarted ERROR: can not find the corresponding channel! (channel =  " << channel << ")" << endl;
         exit(1);
   }
   return;
}

/*************************************************************************************/
//Compton scattering: q + g -> q + gamma and qbar + g -> qbar + gamma
/*************************************************************************************/
void Matrix_elements_sq_Compton(double s, double t, double E1, double E2, double Eq, double T, double* result_ptr)
{
    double prefactor = 2.*e_sq*q_sq*d_F*C_F*g_s*g_s; //2 counts for fermions and anti-fermions

    double pprime = E1;
    double p = E2;
    double k = Eq;
    double kprime = p + pprime - k;

    double p_minus_k_0 = p - k;
    double p_minus_k_i = sqrt(p_minus_k_0*p_minus_k_0 - t);

    //Calculate self-energy coefficients
    Selfenergy_coefficients Sigma02;
    Selfenergy_coefficients* Sigma02ptr = &Sigma02;
    get_quark_selfenergy_coefficients(p_minus_k_0, p_minus_k_i, T, Sigma02ptr);
    double ReA02 = Sigma02.Re_A0;
    double ImA02 = Sigma02.Im_A0;
    double ReB02 = Sigma02.Re_B0;
    double ImB02 = Sigma02.Im_B0;
//    ReA02 = 0.0;  // for test
//    ImA02 = 0.0;
//    ReB02 = 0.0;
//    ImB02 = 0.0;

    //equlibrium contribution
    double trace1_eq, trace2_eq, trace4_eq;
    
    //trace of the 1st term
    trace1_eq = -8*t/(s);
 
    //trace of the 4th term
    double Re_kprime_dot_Q02prime = (-2*kprime*ReB02 + t - ReA02*t)/2.;
    double Im_kprime_dot_Q02prime = (-2*ImB02*kprime - ImA02*t)/2.;
    double Re_p_dot_Q02primestar = (-2*p*ReB02 + t - ReA02*t)/2.;
    double Im_p_dot_Q02primestar = ImB02*p + (ImA02*t)/2.;
    double Q02prime_dot_Q02primestar = Power(ImB02,2) + 2*ImA02*ImB02*(-k + p) 
                        - 2*p*ReB02 - 2*k*(-1 + ReA02)*ReB02 + 2*p*ReA02*ReB02 
                        + Power(ReB02,2) + t + Power(ImA02,2)*t - 2*ReA02*t 
                        + Power(ReA02,2)*t;
    double Re_Q02prime_dot_Q02prime = -Power(ImB02,2) + 2*ImA02*ImB02*(k - p) 
                        - 2*p*ReB02 - 2*k*(-1 + ReA02)*ReB02 + 2*p*ReA02*ReB02 
                        + Power(ReB02,2) + t - Power(ImA02,2)*t - 2*ReA02*t 
                        + Power(ReA02,2)*t;
    double Im_Q02prime_dot_Q02prime = 2*(ImB02*(k + p*(-1 + ReA02) - k*ReA02 + ReB02) 
                        + ImA02*(-(k*ReB02) + p*ReB02 + (-1 + ReA02)*t));
    double Mag_Q02prime_dot_Q02prime_sq = Re_Q02prime_dot_Q02prime*Re_Q02prime_dot_Q02prime
                        + Im_Q02prime_dot_Q02prime*Im_Q02prime_dot_Q02prime;
    double trace4_numerator1 = 2.*Repart_ComplexMultiply(Re_kprime_dot_Q02prime, Im_kprime_dot_Q02prime, Re_p_dot_Q02primestar, Im_p_dot_Q02primestar);
    double trace4_numerator2 = (s + t)/2.*Q02prime_dot_Q02primestar;
    trace4_eq = 16./Mag_Q02prime_dot_Q02prime_sq*(trace4_numerator1 - trace4_numerator2);

    //trace of the 2nd + the 3rd terms
    double Re_Q0_dot_Q02primestar = -((p + pprime)*ReB02);
    double Im_Q0_dot_Q02primestar = ImB02*(p + pprime);
    double Re_Q02primestar_dot_Q02primestar = Re_Q02prime_dot_Q02prime;
    double Im_Q02primestar_dot_Q02primestar = - Im_Q02prime_dot_Q02prime;
    
    double Re_trace2_temp = 1./s*Repart_ComplexDivide(Re_Q0_dot_Q02primestar, Im_Q0_dot_Q02primestar, Re_Q02primestar_dot_Q02primestar, Im_Q02primestar_dot_Q02primestar);
    double Im_trace2_temp = 1./s*Impart_ComplexDivide(Re_Q0_dot_Q02primestar, Im_Q0_dot_Q02primestar, Re_Q02primestar_dot_Q02primestar, Im_Q02primestar_dot_Q02primestar);

    trace2_eq = -64.*(s + t)/2.*Re_trace2_temp;

    result_ptr[0] = prefactor*(trace1_eq + trace2_eq + trace4_eq);
    //if(result<0.0)
    //{
    //  cout << result << endl;
    //  exit(0);
    //}

    //Viscous corrections
    double trace1_vis, trace4_vis, trace2_vis;

    //Self-energy Viscous corrections coefficients 
    double ReA12 = Sigma02.Re_A1;
    double ImA12 = Sigma02.Im_A1;
    double ReB12 = Sigma02.Re_B1;
    double ImB12 = Sigma02.Im_B1;
    double ReC12 = Sigma02.Re_C1;
    double ImC12 = Sigma02.Im_C1;
//    ReA12 = 0.0e0;
//    ImA12 = 0.0e0;
//    ReB12 = 0.0e0;
//    ImB12 = 0.0e0;
//    ReC12 = 0.0e0;
//    ImC12 = 0.0e0;

    // trace of the 1st terms
    trace1_vis = 0.0;

    // trace of the 4th terms
    double Re_Q02primestar_dot_Sigma1prime = -((8*Power(k,4) - 16*Power(k,3)*p 
                 + 8*Power(k,2)*(Power(p,2) - t) + 12*k*p*t + 3*Power(t,2))*(ImB02*(ImB12 
                 + ImA12*(-k + p)) - k*ReA12*ReB02 + p*ReA12*ReB02 + k*ReB12 - p*ReB12 
                 - k*ReA02*ReB12 + p*ReA02*ReB12 + ReB02*ReB12 - 2*ReC12 + 2*ReA02*ReC12 
                 - ReA12*t + ReA02*ReA12*t + ImA02*(2*ImC12 - ImB12*k + ImB12*p + ImA12*t)))
                 /(8.*Power(k,2));
    double Im_Q02primestar_dot_Sigma1prime = ((-2*ImC12*(-1 + ReA02) - ImB02*k*ReA12 
                 + ImB02*p*ReA12 + ImA12*k*ReB02 - ImA12*p*ReB02 - ImB12*(k + p*(-1 + ReA02) 
                 - k*ReA02 + ReB02) + ImB02*ReB12 - ImA02*k*ReB12 + ImA02*p*ReB12 
                 + 2*ImA02*ReC12 + ImA12*t - ImA12*ReA02*t + ImA02*ReA12*t)*(8*Power(k,4) 
                 - 16*Power(k,3)*p + 8*Power(k,2)*(Power(p,2) - t) + 12*k*p*t 
                 + 3*Power(t,2)))/(8.*Power(k,2));
    double Re_p_dot_Sigma1primestar = (8*Power(k,4)*(2*p*ReB12 + ReA12*t) 
                 - 16*Power(k,3)*p*(2*p*ReB12 + 2*ReC12 + ReA12*t) + 12*k*p*t*(2*p*ReB12 
                 + 4*ReC12 + ReA12*t) + 3*Power(t,2)*(2*p*ReB12 + 4*ReC12 + ReA12*t) 
                 + 8*Power(k,2)*(2*Power(p,3)*ReB12 - 2*p*ReB12*t - t*(2*ReC12 + ReA12*t) 
                 + Power(p,2)*(4*ReC12 + ReA12*t)))/(16.*Power(k,2));
    double Im_p_dot_Sigma1primestar = (-((2*ImB12*p + ImA12*t)*(8*Power(k,4) 
                 - 16*Power(k,3)*p + 8*Power(k,2)*(Power(p,2) - t) + 12*k*p*t 
                 + 3*Power(t,2))) + 4*ImC12*(8*Power(k,3)*p - 12*k*p*t - 3*Power(t,2) 
                 + Power(k,2)*(-8*Power(p,2) + 4*t)))/(16.*Power(k,2));
    double Re_kprime_dot_Sigma1prime = ((2*kprime*ReB12 + ReA12*t)*(8*Power(k,4) 
                 - 16*Power(k,3)*p + 8*Power(k,2)*(Power(p,2) - t) + 12*k*p*t 
                 + 3*Power(t,2)) - 4*ReC12*(8*Power(k,3)*kprime + 3*s*t 
                 - 2*Power(k,2)*(4*kprime*p + 3*s + t) + 6*k*(p*s - kprime*t)))
                 /(16.*Power(k,2));
    double Im_kprime_dot_Sigma1prime = ((2*ImB12*kprime + ImA12*t)*(8*Power(k,4) 
                 - 16*Power(k,3)*p + 8*Power(k,2)*(Power(p,2) - t) + 12*k*p*t 
                 + 3*Power(t,2)) - 4*ImC12*(8*Power(k,3)*kprime + 3*s*t 
                 - 2*Power(k,2)*(4*kprime*p + 3*s + t) + 6*k*(p*s - kprime*t)))
                 /(16.*Power(k,2));
    double Re_Q02prime_dot_Sigma1primestar = Re_Q02primestar_dot_Sigma1prime;

    double trace4_vis_temp1 = Repart_ComplexMultiply(Re_Q02prime_dot_Q02prime, Im_Q02prime_dot_Q02prime, Re_Q02primestar_dot_Sigma1prime, Im_Q02primestar_dot_Sigma1prime);
    double trace4_vis_temp2 = 2.*(Repart_ComplexMultiply(Re_kprime_dot_Q02prime, Im_kprime_dot_Q02prime, Re_p_dot_Sigma1primestar, Im_p_dot_Sigma1primestar)
                                + Repart_ComplexMultiply(Re_p_dot_Q02primestar, Im_p_dot_Q02primestar, Re_kprime_dot_Sigma1prime, Im_kprime_dot_Sigma1prime));
    double trace4_vis_temp3 = 2.*(s + t)/2.*Re_Q02prime_dot_Sigma1primestar;

    double trace4_vis_1 = trace4_eq*4./Mag_Q02prime_dot_Q02prime_sq*trace4_vis_temp1;
    double trace4_vis_2 = 16./Mag_Q02prime_dot_Q02prime_sq*trace4_vis_temp2;
    double trace4_vis_3 = 16./Mag_Q02prime_dot_Q02prime_sq*trace4_vis_temp3;
    trace4_vis = trace4_vis_1 - trace4_vis_2 + trace4_vis_3;

    //trace of the 2nd terms
    double Re_Q02primestar_dot_Sigma1primestar = ((8*Power(k,4) - 16*Power(k,3)*p 
                 + 8*Power(k,2)*(Power(p,2) - t) + 12*k*p*t + 3*Power(t,2))*(ImB02*(ImB12 
                 + ImA12*(-k + p)) + k*ReA12*ReB02 - p*ReA12*ReB02 - k*ReB12 + p*ReB12 
                 + k*ReA02*ReB12 - p*ReA02*ReB12 - ReB02*ReB12 + 2*ReC12 - 2*ReA02*ReC12 
                 + ReA12*t - ReA02*ReA12*t + ImA02*(2*ImC12 - ImB12*k + ImB12*p + ImA12*t)))
                 /(8.*Power(k,2));
    double Im_Q02primestar_dot_Sigma1primestar = ((2*ImC12*(-1 + ReA02) - ImB02*k*ReA12 
                 + ImB02*p*ReA12 - ImA12*k*ReB02 + ImA12*p*ReB02 + ImB12*(k + p*(-1 + ReA02) 
                 - k*ReA02 + ReB02) + ImB02*ReB12 - ImA02*k*ReB12 + ImA02*p*ReB12 
                 + 2*ImA02*ReC12 - ImA12*t + ImA12*ReA02*t + ImA02*ReA12*t)*(8*Power(k,4) 
                 - 16*Power(k,3)*p + 8*Power(k,2)*(Power(p,2) - t) + 12*k*p*t 
                 + 3*Power(t,2)))/(8.*Power(k,2));
    double Re_Q0_dot_Sigma1primestar = (8*Power(k,4)*(p + pprime)*ReB12 - 16*Power(k,3)
                 *(p + pprime)*(p*ReB12 + ReC12) + 3*t*(-2*ReC12*s + (p + pprime)*ReB12*t) 
                 + 4*Power(k,2)*(2*Power(p,3)*ReB12 + 2*Power(p,2)*(pprime*ReB12 + 2*ReC12) 
                 + 3*ReC12*s - 2*pprime*ReB12*t + p*(4*pprime*ReC12 - 2*ReB12*t)) 
                 + 12*k*(Power(p,2)*ReB12*t + pprime*ReC12*t + p*(-(ReC12*s) 
                 + pprime*ReB12*t + ReC12*t)))/(8.*Power(k,2));
    double Im_Q0_dot_Sigma1primestar = -(ImB12*(p + pprime)*(8*Power(k,4) - 16*Power(k,3)*p 
                 + 8*Power(k,2)*(Power(p,2) - t) + 12*k*p*t + 3*Power(t,2)) 
                 - 2*ImC12*(8*Power(k,3)*(p + pprime) - 2*Power(k,2)*(4*Power(p,2) 
                 + 4*p*pprime + 3*s) + 3*s*t + 6*k*(p*s - p*t - pprime*t)))/(8.*Power(k,2));
    
    double trace2_vis_Re_temp1 = 2.*Repart_ComplexDivide(Re_Q02primestar_dot_Sigma1primestar, Im_Q02primestar_dot_Sigma1primestar, Re_Q02primestar_dot_Q02primestar, Im_Q02primestar_dot_Q02primestar);
    double trace2_vis_Im_temp1 = 2.*Impart_ComplexDivide(Re_Q02primestar_dot_Sigma1primestar, Im_Q02primestar_dot_Sigma1primestar, Re_Q02primestar_dot_Q02primestar, Im_Q02primestar_dot_Q02primestar);
    double trace2_vis_temp1 = 1./s*Repart_ComplexMultiply(Re_trace2_temp, Im_trace2_temp, trace2_vis_Re_temp1, trace2_vis_Im_temp1);
    double trace2_vis_temp2 = 1./s*Repart_ComplexDivide(Re_Q0_dot_Sigma1primestar, Im_Q0_dot_Sigma1primestar, Re_Q02primestar_dot_Q02primestar, Im_Q02primestar_dot_Q02primestar);
    trace2_vis = -64.*(s + t)/2.*(trace2_vis_temp1 - trace2_vis_temp2);

    result_ptr[1] = prefactor*(trace1_vis + trace2_vis + trace4_vis);

    return;
}

/*************************************************************************************/
//Pair annihilation: q + qbar -> g + gamma
/*************************************************************************************/
void Matrix_elements_sq_Annihilation(double s, double t, double E1, double E2, double Eq, double T, double* result_ptr)
{
    double u = - s - t;
    double prefactor = e_sq*q_sq*d_F*C_F*g_s*g_s;
    
    double pprime = E1;
    double p = E2;
    double k = Eq;
    double kprime = p + pprime - k;

    double p_minus_k_0 = p - k;
    double p_minus_k_i = sqrt(p_minus_k_0*p_minus_k_0 - t);
    double p_minus_kprime_0 = p - kprime;
    double p_minus_kprime_i = sqrt(p_minus_kprime_0*p_minus_kprime_0 - u);

    //Calculate self-energy coefficients
    Selfenergy_coefficients Sigma02;
    Selfenergy_coefficients* Sigma02ptr = &Sigma02;
    get_quark_selfenergy_coefficients(p_minus_k_0, p_minus_k_i, T, Sigma02ptr);
    double ReA02 = Sigma02.Re_A0;
    double ImA02 = Sigma02.Im_A0;
    double ReB02 = Sigma02.Re_B0;
    double ImB02 = Sigma02.Im_B0;
//    ReA02 = 0.0e0;
//    ImA02 = 0.0e0;
//    ReB02 = 0.0e0;
//    ImB02 = 0.0e0;

    Selfenergy_coefficients Sigma03;
    Selfenergy_coefficients* Sigma03ptr = &Sigma03;
    get_quark_selfenergy_coefficients(p_minus_kprime_0, p_minus_kprime_i, T, Sigma03ptr);
    double ReA03 = Sigma03.Re_A0;
    double ImA03 = Sigma03.Im_A0;
    double ReB03 = Sigma03.Re_B0;
    double ImB03 = Sigma03.Im_B0;
//    ReA03 = 0.0e0;
//    ImA03 = 0.0e0;
//    ReB03 = 0.0e0;
//    ImB03 = 0.0e0;
    
    //equlibrium contribution
    double trace1_eq, trace2_eq, trace4_eq;
    
    //trace of the 1st term
    double Re_pprime_dot_Q02prime = (-2*pprime*ReB02 + (-1 + ReA02)*t)/2.;
    double Im_pprime_dot_Q02prime = -(ImB02*pprime) + (ImA02*t)/2.;
    double Re_p_dot_Q02primestar = (-2*p*ReB02 + t - ReA02*t)/2.;
    double Im_p_dot_Q02primestar = ImB02*p + (ImA02*t)/2.;
    double Q02prime_dot_Q02primestar = Power(ImB02,2) + 2*ImA02*ImB02*(-k + p) 
                        - 2*p*ReB02 - 2*k*(-1 + ReA02)*ReB02 + 2*p*ReA02*ReB02 
                        + Power(ReB02,2) + t + Power(ImA02,2)*t - 2*ReA02*t 
                        + Power(ReA02,2)*t;
    double Re_Q02prime_dot_Q02prime = -Power(ImB02,2) + 2*ImA02*ImB02*(k - p) 
                        - 2*p*ReB02 - 2*k*(-1 + ReA02)*ReB02 + 2*p*ReA02*ReB02 
                        + Power(ReB02,2) + t - Power(ImA02,2)*t - 2*ReA02*t 
                        + Power(ReA02,2)*t;
    double Im_Q02prime_dot_Q02prime = 2*(ImB02*(k + p*(-1 + ReA02) - k*ReA02 + ReB02) 
                        + ImA02*(-(k*ReB02) + p*ReB02 + (-1 + ReA02)*t));
    double Mag_Q02prime_dot_Q02prime_sq = Re_Q02prime_dot_Q02prime*Re_Q02prime_dot_Q02prime
                        + Im_Q02prime_dot_Q02prime*Im_Q02prime_dot_Q02prime;

    double trace1_numerator1 = 2.*Repart_ComplexMultiply(Re_pprime_dot_Q02prime, Im_pprime_dot_Q02prime, Re_p_dot_Q02primestar, Im_p_dot_Q02primestar);
    double trace1_numerator2 = s/2.*Q02prime_dot_Q02primestar;
    trace1_eq = 16./Mag_Q02prime_dot_Q02prime_sq*(trace1_numerator1 - trace1_numerator2);

    //trace of the 4th term
    double Re_pprime_dot_Q03Tilde = (-2*pprime*ReB03 - (-1 + ReA03)*(s + t))/2.;
    double Im_pprime_dot_Q03Tilde = (-2*ImB03*pprime - ImA03*s - ImA03*t)/2.;
    double Re_p_dot_Q03Tildestar = (-2*p*ReB03 + (-1 + ReA03)*(s + t))/2.;
    double Im_p_dot_Q03Tildestar = ImB03*p - (ImA03*(s + t))/2.;
    double Q03Tilde_dot_Q03TIildestar = Power(ImB03,2) - 2*ImA03*ImB03*kprime + 2*ImA03*ImB03*p 
                        + 2*kprime*ReB03 - 2*p*ReB03 - 2*kprime*ReA03*ReB03 + 2*p*ReA03*ReB03 
                        + Power(ReB03,2) - s - t - Power(ImA03,2)*(s + t) + 2*ReA03*(s + t) 
                        - Power(ReA03,2)*(s + t);
    double Re_Q03Tilde_dot_Q03Tilde = -Power(ImB03,2) + 2*ImA03*ImB03*kprime - 2*ImA03*ImB03*p 
                        + 2*kprime*ReB03 - 2*p*ReB03 - 2*kprime*ReA03*ReB03 + 2*p*ReA03*ReB03 
                        + Power(ReB03,2) - s - t + Power(ImA03,2)*(s + t) + 2*ReA03*(s + t) 
                        - Power(ReA03,2)*(s + t);
    double Im_Q03Tilde_dot_Q03Tilde = 2*(ImB03*(kprime + p*(-1 + ReA03) - kprime*ReA03 + ReB03) 
                        + ImA03*(-(kprime*ReB03) + p*ReB03 - (-1 + ReA03)*(s + t)));
    double Mag_Q03Tilde_dot_Q03Tilde_sq = Re_Q03Tilde_dot_Q03Tilde*Re_Q03Tilde_dot_Q03Tilde
                        + Im_Q03Tilde_dot_Q03Tilde*Im_Q03Tilde_dot_Q03Tilde;

    double trace4_numerator1 = 2.*Repart_ComplexMultiply(Re_pprime_dot_Q03Tilde, Im_pprime_dot_Q03Tilde, Re_p_dot_Q03Tildestar, Im_p_dot_Q03Tildestar);
    double trace4_numerator2 = s/2.*Q03Tilde_dot_Q03TIildestar;
    trace4_eq = 16./Mag_Q03Tilde_dot_Q03Tilde_sq*(trace4_numerator1 - trace4_numerator2);

    //trace of the 2nd + 3rd terms
    double Re_Q03Tilde_dot_Q02primestar = ImA02*ImB03*(-k + p) + ImB02*(ImB03 + ImA03*(-kprime + p))
                        + kprime*ReB02 - p*ReB02 - kprime*ReA03*ReB02 + p*ReA03*ReB02 
                        + k*ReB03 - p*ReB03 - k*ReA02*ReB03 + p*ReA02*ReB03 + ReB02*ReB03;
    double Im_Q03Tilde_dot_Q02primestar = -(ImA03*kprime*ReB02) + ImA03*p*ReB02 
                        + ImB03*(k + p*(-1 + ReA02) - k*ReA02 + ReB02) + ImA02*k*ReB03 
                        - ImA02*p*ReB03 - ImB02*(kprime + p*(-1 + ReA03) - kprime*ReA03 
                        + ReB03);
    double Re_Q02primestar_dot_Q02primestar = Re_Q02prime_dot_Q02prime;
    double Im_Q02primestar_dot_Q02primestar = - Im_Q02prime_dot_Q02prime;

    double Re_trace2_eq_denominator = Repart_ComplexMultiply(Re_Q02primestar_dot_Q02primestar, Im_Q02primestar_dot_Q02primestar, Re_Q03Tilde_dot_Q03Tilde, Im_Q03Tilde_dot_Q03Tilde);
    double Im_trace2_eq_denominator = Impart_ComplexMultiply(Re_Q02primestar_dot_Q02primestar, Im_Q02primestar_dot_Q02primestar, Re_Q03Tilde_dot_Q03Tilde, Im_Q03Tilde_dot_Q03Tilde);
    double Re_trace2_eq = Repart_ComplexDivide(Re_Q03Tilde_dot_Q02primestar, Im_Q03Tilde_dot_Q02primestar, Re_trace2_eq_denominator, Im_trace2_eq_denominator);
    double Im_trace2_eq = Impart_ComplexDivide(Re_Q03Tilde_dot_Q02primestar, Im_Q03Tilde_dot_Q02primestar, Re_trace2_eq_denominator, Im_trace2_eq_denominator);
    trace2_eq = -64.*s/2.*Re_trace2_eq;
    
    result_ptr[0] = prefactor*(trace1_eq + trace2_eq + trace4_eq);
    //result = prefactor*16.*(u/t);

    //Viscous corrections
    double trace1_vis, trace2_vis, trace4_vis;
    //Self-energy Viscous corrections coefficients 
    double ReA12 = Sigma02.Re_A1;
    double ImA12 = Sigma02.Im_A1;
    double ReB12 = Sigma02.Re_B1;
    double ImB12 = Sigma02.Im_B1;
    double ReC12 = Sigma02.Re_C1;
    double ImC12 = Sigma02.Im_C1;

    double ReA13 = Sigma03.Re_A1;
    double ImA13 = Sigma03.Im_A1;
    double ReB13 = Sigma03.Re_B1;
    double ImB13 = Sigma03.Im_B1;
    double ReC13 = Sigma03.Re_C1;
    double ImC13 = Sigma03.Im_C1;
//    ReA12 = 0.0e0;
//    ImA12 = 0.0e0;
//    ReB12 = 0.0e0;
//    ImB12 = 0.0e0;
//    ReC12 = 0.0e0;
//    ImC12 = 0.0e0;
//    ReA13 = 0.0e0;
//    ImA13 = 0.0e0;
//    ReB13 = 0.0e0;
//    ImB13 = 0.0e0;
//    ReC13 = 0.0e0;
//    ImC13 = 0.0e0;

    //trace of the 1st terms
    double Re_Q02primestar_dot_Sigma1primestar = ((8*Power(k,4) - 16*Power(k,3)*p 
                          + 8*Power(k,2)*(Power(p,2) - t) + 12*k*p*t + 3*Power(t,2))
                          *(ImB02*(ImB12 + ImA12*(-k + p)) + k*ReA12*ReB02 
                          - p*ReA12*ReB02 - k*ReB12 + p*ReB12 + k*ReA02*ReB12 
                          - p*ReA02*ReB12 - ReB02*ReB12 + 2*ReC12 - 2*ReA02*ReC12 
                          + ReA12*t - ReA02*ReA12*t + ImA02*(2*ImC12 - ImB12*k + ImB12*p 
                          + ImA12*t)))/(8.*Power(k,2));
    double Im_Q02primestar_dot_Sigma1primestar = ((2*ImC12*(-1 + ReA02) - ImB02*k*ReA12 
                          + ImB02*p*ReA12 - ImA12*k*ReB02 + ImA12*p*ReB02 
                          + ImB12*(k + p*(-1 + ReA02) - k*ReA02 + ReB02) + ImB02*ReB12 
                          - ImA02*k*ReB12 + ImA02*p*ReB12 + 2*ImA02*ReC12 - ImA12*t 
                          + ImA12*ReA02*t + ImA02*ReA12*t)*(8*Power(k,4) - 16*Power(k,3)*p 
                          + 8*Power(k,2)*(Power(p,2) - t) + 12*k*p*t + 3*Power(t,2)))
                          /(8.*Power(k,2));
    double Re_p_dot_Sigma1primestar = (8*Power(k,4)*(2*p*ReB12 + ReA12*t) 
                          - 16*Power(k,3)*p*(2*p*ReB12 + 2*ReC12 + ReA12*t) 
                          + 12*k*p*t*(2*p*ReB12 + 4*ReC12 + ReA12*t) 
                          + 3*Power(t,2)*(2*p*ReB12 + 4*ReC12 + ReA12*t) 
                          + 8*Power(k,2)*(2*Power(p,3)*ReB12 - 2*p*ReB12*t 
                          - t*(2*ReC12 + ReA12*t) + Power(p,2)*(4*ReC12 + ReA12*t)))
                          /(16.*Power(k,2));
    double Im_p_dot_Sigma1primestar = (-((2*ImB12*p + ImA12*t)*(8*Power(k,4) 
                          - 16*Power(k,3)*p + 8*Power(k,2)*(Power(p,2) - t) + 12*k*p*t 
                          + 3*Power(t,2))) + 4*ImC12*(8*Power(k,3)*p - 12*k*p*t 
                          - 3*Power(t,2) + Power(k,2)*(-8*Power(p,2) + 4*t)))
                          /(16.*Power(k,2));
    double Re_pprime_dot_Sigma1prime = (8*Power(k,4)*(2*pprime*ReB12 - ReA12*t) 
                          - 16*Power(k,3)*(2*p*pprime*ReB12 + 2*pprime*ReC12 - p*ReA12*t) 
                          + 8*Power(k,2)*(4*p*pprime*ReC12 + 3*ReC12*s - 2*pprime*ReB12*t 
                          + 2*ReC12*t + ReA12*Power(t,2) + Power(p,2)*(2*pprime*ReB12 
                          - ReA12*t)) - 3*t*(4*ReC12*(s + t) + t*(-2*pprime*ReB12 
                          + ReA12*t)) - 12*k*(-2*pprime*ReC12*t + p*(2*ReC12*(s + t) 
                          + t*(-2*pprime*ReB12 + ReA12*t))))/(16.*Power(k,2));
    double Im_pprime_dot_Sigma1prime = ((2*ImB12*pprime - ImA12*t)*(8*Power(k,4) 
                          - 16*Power(k,3)*p + 8*Power(k,2)*(Power(p,2) - t) + 12*k*p*t 
                          + 3*Power(t,2)) - 4*ImC12*(8*Power(k,3)*pprime + 3*t*(s + t) 
                          - 2*Power(k,2)*(4*p*pprime + 3*s + 2*t) + 6*k*(-(pprime*t) 
                          + p*(s + t))))/(16.*Power(k,2));
    double Re_Q02prime_dot_Sigma1primestar = -((8*Power(k,4) - 16*Power(k,3)*p 
                          + 8*Power(k,2)*(Power(p,2) - t) + 12*k*p*t 
                          + 3*Power(t,2))*(ImB02*(ImB12 + ImA12*(-k + p)) - k*ReA12*ReB02 
                          + p*ReA12*ReB02 + k*ReB12 - p*ReB12 - k*ReA02*ReB12 
                          + p*ReA02*ReB12 + ReB02*ReB12 - 2*ReC12 + 2*ReA02*ReC12 - ReA12*t 
                          + ReA02*ReA12*t + ImA02*(2*ImC12 - ImB12*k + ImB12*p + ImA12*t)))
                          /(8.*Power(k,2));
    
    double trace1_vis_temp1 = Repart_ComplexMultiply(Re_Q02prime_dot_Q02prime, Im_Q02prime_dot_Q02prime, Re_Q02primestar_dot_Sigma1primestar, Im_Q02primestar_dot_Sigma1primestar);
    double trace1_vis_temp2 = 2*(Repart_ComplexMultiply(Re_pprime_dot_Q02prime, Im_pprime_dot_Q02prime, Re_p_dot_Sigma1primestar, Im_p_dot_Sigma1primestar) 
                              + Repart_ComplexMultiply(Re_p_dot_Q02primestar, Im_p_dot_Q02primestar, Re_pprime_dot_Sigma1prime, Im_pprime_dot_Sigma1prime));
    double trace1_vis_temp3 = 2.*s/2.*Re_Q02prime_dot_Sigma1primestar;

    double trace1_vis_1 = trace1_eq*4./Mag_Q02prime_dot_Q02prime_sq*trace1_vis_temp1;
    double trace1_vis_2 = 16./Mag_Q02prime_dot_Q02prime_sq*trace1_vis_temp2;
    double trace1_vis_3 = 16./Mag_Q02prime_dot_Q02prime_sq*trace1_vis_temp3;
    trace1_vis = trace1_vis_1 - trace1_vis_2 + trace1_vis_3;

    //trace of the 4th terms
    double Re_Q03Tildestar_dot_Sigma1Tildestar = ((4*Power(k,2)*(2*Power(kprime,2) - 4*kprime*p 
                          + 2*Power(p,2) - s - t) - 12*k*(kprime - p)*(s + t) 
                          + 3*Power(s + t,2))*(ImB03*(ImB13 + ImA13*(-kprime + p)) 
                          + kprime*ReA13*ReB03 - p*ReA13*ReB03 - kprime*ReB13 + p*ReB13 
                          + kprime*ReA03*ReB13 - p*ReA03*ReB13 - ReB03*ReB13 + 2*ReC13 
                          - 2*ReA03*ReC13 - ReA13*s + ReA03*ReA13*s - ReA13*t 
                          + ReA03*ReA13*t + ImA03*(2*ImC13 + ImB13*(-kprime + p) 
                          - ImA13*(s + t))))/(8.*Power(k,2));
    double Im_Q03Tildestar_dot_Sigma1Tildestar = ((2*ImC13*(-1 + ReA03) - ImB03*kprime*ReA13 
                          + ImB03*p*ReA13 - ImA13*kprime*ReB03 + ImA13*p*ReB03 
                          + ImB13*(kprime + p*(-1 + ReA03) - kprime*ReA03 + ReB03) 
                          + ImB03*ReB13 - ImA03*kprime*ReB13 + ImA03*p*ReB13 + 2*ImA03*ReC13 
                          + ImA13*s - ImA13*ReA03*s - ImA03*ReA13*s + ImA13*t - ImA13*ReA03*t 
                          - ImA03*ReA13*t)*(4*Power(k,2)*(2*Power(kprime,2) - 4*kprime*p + 2*Power(p,2) 
                          - s - t) - 12*k*(kprime - p)*(s + t) + 3*Power(s + t,2)))
                          /(8.*Power(k,2));
    double Re_p_dot_Sigma1Tildestar = (-3*(s + t)*(-4*ReC13*t - 2*p*ReB13*(s + t) 
                          + ReA13*Power(s + t,2)) + 4*Power(k,2)*(4*Power(p,3)*ReB13 
                          - 2*p*ReB13*(s + t) + Power(kprime,2)*(4*p*ReB13 - 2*ReA13*(s + t)) 
                          - 2*Power(p,2)*(-4*ReC13 + ReA13*(s + t)) + (s + t)*(-2*ReC13 
                          + ReA13*(s + t)) + 4*kprime*p*(-2*p*ReB13 - 2*ReC13 
                          + ReA13*(s + t))) - 12*k*(kprime*(2*ReC13*t + 2*p*ReB13*(s + t) 
                          - ReA13*Power(s + t,2)) + p*(-2*p*ReB13*(s + t) 
                          + ReA13*Power(s + t,2) - 2*ReC13*(s + 2*t))))/(16.*Power(k,2));
    double Im_p_dot_Sigma1Tildestar = (-((2*ImB13*p - ImA13*(s + t))*(4*Power(k,2)
                          *(2*Power(kprime,2) - 4*kprime*p + 2*Power(p,2) - s - t) 
                          - 12*k*(kprime - p)*(s + t) + 3*Power(s + t,2))) 
                          + 4*ImC13*(-3*t*(s + t) + 2*Power(k,2)*(4*kprime*p - 4*Power(p,2) 
                          + s + t) + 6*k*(kprime*t - p*(s + 2*t))))/(16.*Power(k,2));
    double Re_pprime_dot_Sigma1Tilde = (3*Power(s + t,2)*(2*pprime*ReB13 - 4*ReC13 
                          + ReA13*(s + t)) + 12*k*(s + t)*(2*pprime*ReC13 
                          - kprime*(2*pprime*ReB13 - 2*ReC13 + ReA13*(s + t)) 
                          + p*(2*pprime*ReB13 - 2*ReC13 + ReA13*(s + t))) 
                          + 4*Power(k,2)*(8*p*pprime*ReC13 + 2*Power(kprime,2)
                          *(2*pprime*ReB13 + ReA13*(s + t)) + 2*Power(p,2)*(2*pprime*ReB13 
                          + ReA13*(s + t)) - (s + t)*(2*pprime*ReB13 - 2*ReC13 
                          + ReA13*(s + t)) - 4*kprime*(2*pprime*ReC13 + p*(2*pprime*ReB13 
                          + ReA13*(s + t)))))/(16.*Power(k,2));
    double Im_pprime_dot_Sigma1Tilde =((2*ImB13*pprime + ImA13*(s + t))*(4*Power(k,2)
                          *(2*Power(kprime,2) - 4*kprime*p + 2*Power(p,2) - s - t) 
                          - 12*k*(kprime - p)*(s + t) + 3*Power(s + t,2)) 
                          - 4*ImC13*(-6*k*(kprime - p + pprime)*(s + t) + 3*Power(s + t,2) 
                          + Power(k,2)*(8*kprime*pprime - 2*(4*p*pprime + s + t))))
                          /(16.*Power(k,2));
    double Re_Q03Tilde_dot_Sigma1Tildestar = -((4*Power(k,2)*(2*Power(kprime,2) - 4*kprime*p 
                          + 2*Power(p,2) - s - t) - 12*k*(kprime - p)*(s + t) 
                          + 3*Power(s + t,2))*(ImB03*(ImB13 + ImA13*(-kprime + p)) 
                          - kprime*ReA13*ReB03 + p*ReA13*ReB03 + kprime*ReB13 - p*ReB13 
                          - kprime*ReA03*ReB13 + p*ReA03*ReB13 + ReB03*ReB13 - 2*ReC13 
                          + 2*ReA03*ReC13 + ReA13*s - ReA03*ReA13*s + ReA13*t - ReA03*ReA13*t 
                          + ImA03*(2*ImC13 + ImB13*(-kprime + p) - ImA13*(s + t))))
                          /(8.*Power(k,2));
    
    double trace4_vis_temp1 = Repart_ComplexMultiply(Re_Q03Tilde_dot_Q03Tilde, Im_Q03Tilde_dot_Q03Tilde, Re_Q03Tildestar_dot_Sigma1Tildestar, Im_Q03Tildestar_dot_Sigma1Tildestar);
    double trace4_vis_temp2 = 2*(Repart_ComplexMultiply(Re_pprime_dot_Q03Tilde, Im_pprime_dot_Q03Tilde, Re_p_dot_Sigma1Tildestar, Im_p_dot_Sigma1Tildestar)
                              + Repart_ComplexMultiply(Re_p_dot_Q03Tildestar, Im_p_dot_Q03Tildestar, Re_pprime_dot_Sigma1Tilde, Im_pprime_dot_Sigma1Tilde));
    double trace4_vis_temp3 = 2.*s/2.*Re_Q03Tilde_dot_Sigma1Tildestar;

    double trace4_vis_1 = trace4_eq*4./Mag_Q03Tilde_dot_Q03Tilde_sq*trace4_vis_temp1;
    double trace4_vis_2 = 16./Mag_Q03Tilde_dot_Q03Tilde_sq*trace4_vis_temp2;
    double trace4_vis_3 = 16./Mag_Q03Tilde_dot_Q03Tilde_sq*trace4_vis_temp3;
    trace4_vis = trace4_vis_1 - trace4_vis_2 + trace4_vis_3;

    //trace of the 2nd + the 3rd terms
    double Re_Q03Tilde_dot_Sigma1Tilde = ((4*Power(k,2)*(2*Power(kprime,2) - 4*kprime*p 
                          + 2*Power(p,2) - s - t) - 12*k*(kprime - p)*(s + t) 
                          + 3*Power(s + t,2))*(ImB03*(ImB13 + ImA13*(-kprime + p)) 
                          + kprime*ReA13*ReB03 - p*ReA13*ReB03 - kprime*ReB13 + p*ReB13 
                          + kprime*ReA03*ReB13 - p*ReA03*ReB13 - ReB03*ReB13 + 2*ReC13 
                          - 2*ReA03*ReC13 - ReA13*s + ReA03*ReA13*s - ReA13*t 
                          + ReA03*ReA13*t + ImA03*(2*ImC13 + ImB13*(-kprime + p) 
                          - ImA13*(s + t))))/(8.*Power(k,2));
    double Im_Q03Tilde_dot_Sigma1Tilde = -((2*ImC13*(-1 + ReA03) - ImB03*kprime*ReA13 
                          + ImB03*p*ReA13 - ImA13*kprime*ReB03 + ImA13*p*ReB03 
                          + ImB13*(kprime + p*(-1 + ReA03) - kprime*ReA03 + ReB03) 
                          + ImB03*ReB13 - ImA03*kprime*ReB13 + ImA03*p*ReB13 + 2*ImA03*ReC13 
                          + ImA13*s - ImA13*ReA03*s - ImA03*ReA13*s + ImA13*t - ImA13*ReA03*t 
                          - ImA03*ReA13*t)*(4*Power(k,2)*(2*Power(kprime,2) - 4*kprime*p 
                          + 2*Power(p,2) - s - t) - 12*k*(kprime - p)*(s + t) 
                          + 3*Power(s + t,2)))/(8.*Power(k,2));
    double Re_Q03Tilde_dot_Sigma1primestar = -(-8*Power(k,5)*ReA12*ReB03 
                          + 24*Power(k,4)*p*ReA12*ReB03 - 24*Power(k,3)*Power(p,2)*ReA12*ReB03 
                          + 8*Power(k,2)*Power(p,3)*ReA12*ReB03 + 8*Power(k,4)*kprime*ReB12 
                          - 8*Power(k,4)*p*ReB12 - 16*Power(k,3)*kprime*p*ReB12 + 
      16*Power(k,3)*Power(p,2)*ReB12 + 
      8*Power(k,2)*kprime*Power(p,2)*ReB12 - 
      8*Power(k,2)*Power(p,3)*ReB12 - 
      8*Power(k,4)*kprime*ReA03*ReB12 + 
      8*Power(k,4)*p*ReA03*ReB12 + 
      16*Power(k,3)*kprime*p*ReA03*ReB12 - 
      16*Power(k,3)*Power(p,2)*ReA03*ReB12 - 
      8*Power(k,2)*kprime*Power(p,2)*ReA03*ReB12 + 
      8*Power(k,2)*Power(p,3)*ReA03*ReB12 + 
      8*Power(k,4)*ReB03*ReB12 - 16*Power(k,3)*p*ReB03*ReB12 + 
      8*Power(k,2)*Power(p,2)*ReB03*ReB12 - 
      16*Power(k,3)*kprime*ReC12 + 16*Power(k,3)*p*ReC12 + 
      16*Power(k,2)*kprime*p*ReC12 - 
      16*Power(k,2)*Power(p,2)*ReC12 + 
      16*Power(k,3)*kprime*ReA03*ReC12 - 
      16*Power(k,3)*p*ReA03*ReC12 - 
      16*Power(k,2)*kprime*p*ReA03*ReC12 + 
      16*Power(k,2)*Power(p,2)*ReA03*ReC12 + 
      12*Power(k,2)*ReC12*s - 12*k*p*ReC12*s - 
      12*Power(k,2)*ReA03*ReC12*s + 12*k*p*ReA03*ReC12*s + 
      8*Power(k,3)*ReA12*ReB03*t - 20*Power(k,2)*p*ReA12*ReB03*t + 
      12*k*Power(p,2)*ReA12*ReB03*t - 8*Power(k,2)*kprime*ReB12*t + 
      8*Power(k,2)*p*ReB12*t + 12*k*kprime*p*ReB12*t - 
      12*k*Power(p,2)*ReB12*t + 8*Power(k,2)*kprime*ReA03*ReB12*t - 
      8*Power(k,2)*p*ReA03*ReB12*t - 12*k*kprime*p*ReA03*ReB12*t + 
      12*k*Power(p,2)*ReA03*ReB12*t - 8*Power(k,2)*ReB03*ReB12*t + 
      12*k*p*ReB03*ReB12*t + 12*Power(k,2)*ReC12*t + 
      12*k*kprime*ReC12*t - 24*k*p*ReC12*t - 
      12*Power(k,2)*ReA03*ReC12*t - 12*k*kprime*ReA03*ReC12*t + 
      24*k*p*ReA03*ReC12*t - 6*ReC12*s*t + 6*ReA03*ReC12*s*t - 
      3*k*ReA12*ReB03*Power(t,2) + 3*p*ReA12*ReB03*Power(t,2) + 
      3*kprime*ReB12*Power(t,2) - 3*p*ReB12*Power(t,2) - 
      3*kprime*ReA03*ReB12*Power(t,2) + 
      3*p*ReA03*ReB12*Power(t,2) + 3*ReB03*ReB12*Power(t,2) - 
      6*ReC12*Power(t,2) + 6*ReA03*ReC12*Power(t,2) + 
      ImB03*(ImB12 + ImA12*(-k + p))*
       (8*Power(k,4) - 16*Power(k,3)*p + 
         8*Power(k,2)*(Power(p,2) - t) + 12*k*p*t + 3*Power(t,2)) + 
      ImA03*(-(ImB12*(kprime - p)*
            (8*Power(k,4) - 16*Power(k,3)*p + 
              8*Power(k,2)*(Power(p,2) - t) + 12*k*p*t + 
              3*Power(t,2))) + 
         2*ImC12*(8*Power(k,3)*(kprime - p) + 3*t*(s + t) - 
            2*Power(k,2)*(4*kprime*p - 4*Power(p,2) + 3*(s + t)) + 
            6*k*(-(kprime*t) + p*(s + 2*t)))))/(8.*Power(k,2));
    double Im_Q03Tilde_dot_Sigma1primestar =(8*ImB03*Power(k,5)*ReA12 - 24*ImB03*Power(k,4)*p*ReA12 + 
     24*ImB03*Power(k,3)*Power(p,2)*ReA12 - 
     8*ImB03*Power(k,2)*Power(p,3)*ReA12 - 
     8*ImA12*Power(k,5)*ReB03 + 24*ImA12*Power(k,4)*p*ReB03 - 
     24*ImA12*Power(k,3)*Power(p,2)*ReB03 + 
     8*ImA12*Power(k,2)*Power(p,3)*ReB03 - 
     8*ImB03*Power(k,4)*ReB12 + 8*ImA03*Power(k,4)*kprime*ReB12 + 
     16*ImB03*Power(k,3)*p*ReB12 - 8*ImA03*Power(k,4)*p*ReB12 - 
     16*ImA03*Power(k,3)*kprime*p*ReB12 - 
     8*ImB03*Power(k,2)*Power(p,2)*ReB12 + 
     16*ImA03*Power(k,3)*Power(p,2)*ReB12 + 
     8*ImA03*Power(k,2)*kprime*Power(p,2)*ReB12 - 
     8*ImA03*Power(k,2)*Power(p,3)*ReB12 - 
     16*ImA03*Power(k,3)*kprime*ReC12 + 
     16*ImA03*Power(k,3)*p*ReC12 + 
     16*ImA03*Power(k,2)*kprime*p*ReC12 - 
     16*ImA03*Power(k,2)*Power(p,2)*ReC12 + 
     12*ImA03*Power(k,2)*ReC12*s - 12*ImA03*k*p*ReC12*s - 
     8*ImB03*Power(k,3)*ReA12*t + 20*ImB03*Power(k,2)*p*ReA12*t - 
     12*ImB03*k*Power(p,2)*ReA12*t + 8*ImA12*Power(k,3)*ReB03*t - 
     20*ImA12*Power(k,2)*p*ReB03*t + 
     12*ImA12*k*Power(p,2)*ReB03*t + 8*ImB03*Power(k,2)*ReB12*t - 
     8*ImA03*Power(k,2)*kprime*ReB12*t - 12*ImB03*k*p*ReB12*t + 
     8*ImA03*Power(k,2)*p*ReB12*t + 12*ImA03*k*kprime*p*ReB12*t - 
     12*ImA03*k*Power(p,2)*ReB12*t + 12*ImA03*Power(k,2)*ReC12*t + 
     12*ImA03*k*kprime*ReC12*t - 24*ImA03*k*p*ReC12*t - 
     6*ImA03*ReC12*s*t + 3*ImB03*k*ReA12*Power(t,2) - 
     3*ImB03*p*ReA12*Power(t,2) - 3*ImA12*k*ReB03*Power(t,2) + 
     3*ImA12*p*ReB03*Power(t,2) - 3*ImB03*ReB12*Power(t,2) + 
     3*ImA03*kprime*ReB12*Power(t,2) - 3*ImA03*p*ReB12*Power(t,2) - 
     6*ImA03*ReC12*Power(t,2) + 
     ImB12*(kprime + p*(-1 + ReA03) - kprime*ReA03 + ReB03)*
      (8*Power(k,4) - 16*Power(k,3)*p + 
        8*Power(k,2)*(Power(p,2) - t) + 12*k*p*t + 3*Power(t,2)) + 
     2*ImC12*(-1 + ReA03)*(8*Power(k,3)*(kprime - p) + 
        3*t*(s + t) - 2*Power(k,2)*
         (4*kprime*p - 4*Power(p,2) + 3*(s + t)) + 
        6*k*(-(kprime*t) + p*(s + 2*t))))/(8.*Power(k,2));
    double Re_Q02primestar_dot_Sigma1Tilde = -(-8*Power(k,2)*Power(kprime,3)*ReA13*ReB02 + 
      24*Power(k,2)*Power(kprime,2)*p*ReA13*ReB02 - 
      24*Power(k,2)*kprime*Power(p,2)*ReA13*ReB02 + 
      8*Power(k,2)*Power(p,3)*ReA13*ReB02 + 
      8*Power(k,3)*Power(kprime,2)*ReB13 - 
      16*Power(k,3)*kprime*p*ReB13 - 
      8*Power(k,2)*Power(kprime,2)*p*ReB13 + 
      8*Power(k,3)*Power(p,2)*ReB13 + 
      16*Power(k,2)*kprime*Power(p,2)*ReB13 - 
      8*Power(k,2)*Power(p,3)*ReB13 - 
      8*Power(k,3)*Power(kprime,2)*ReA02*ReB13 + 
      16*Power(k,3)*kprime*p*ReA02*ReB13 + 
      8*Power(k,2)*Power(kprime,2)*p*ReA02*ReB13 - 
      8*Power(k,3)*Power(p,2)*ReA02*ReB13 - 
      16*Power(k,2)*kprime*Power(p,2)*ReA02*ReB13 + 
      8*Power(k,2)*Power(p,3)*ReA02*ReB13 + 
      8*Power(k,2)*Power(kprime,2)*ReB02*ReB13 - 
      16*Power(k,2)*kprime*p*ReB02*ReB13 + 
      8*Power(k,2)*Power(p,2)*ReB02*ReB13 - 
      16*Power(k,3)*kprime*ReC13 + 16*Power(k,3)*p*ReC13 + 
      16*Power(k,2)*kprime*p*ReC13 - 
      16*Power(k,2)*Power(p,2)*ReC13 + 
      16*Power(k,3)*kprime*ReA02*ReC13 - 
      16*Power(k,3)*p*ReA02*ReC13 - 
      16*Power(k,2)*kprime*p*ReA02*ReC13 + 
      16*Power(k,2)*Power(p,2)*ReA02*ReC13 + 
      4*Power(k,2)*kprime*ReA13*ReB02*s + 
      12*k*Power(kprime,2)*ReA13*ReB02*s - 
      4*Power(k,2)*p*ReA13*ReB02*s - 24*k*kprime*p*ReA13*ReB02*s + 
      12*k*Power(p,2)*ReA13*ReB02*s - 4*Power(k,3)*ReB13*s - 
      12*Power(k,2)*kprime*ReB13*s + 16*Power(k,2)*p*ReB13*s + 
      12*k*kprime*p*ReB13*s - 12*k*Power(p,2)*ReB13*s + 
      4*Power(k,3)*ReA02*ReB13*s + 
      12*Power(k,2)*kprime*ReA02*ReB13*s - 
      16*Power(k,2)*p*ReA02*ReB13*s - 12*k*kprime*p*ReA02*ReB13*s + 
      12*k*Power(p,2)*ReA02*ReB13*s - 4*Power(k,2)*ReB02*ReB13*s - 
      12*k*kprime*ReB02*ReB13*s + 12*k*p*ReB02*ReB13*s + 
      12*Power(k,2)*ReC13*s - 12*k*p*ReC13*s - 
      12*Power(k,2)*ReA02*ReC13*s + 12*k*p*ReA02*ReC13*s - 
      3*kprime*ReA13*ReB02*Power(s,2) + 
      3*p*ReA13*ReB02*Power(s,2) + 3*k*ReB13*Power(s,2) - 
      3*p*ReB13*Power(s,2) - 3*k*ReA02*ReB13*Power(s,2) + 
      3*p*ReA02*ReB13*Power(s,2) + 3*ReB02*ReB13*Power(s,2) + 
      4*Power(k,2)*kprime*ReA13*ReB02*t + 
      12*k*Power(kprime,2)*ReA13*ReB02*t - 
      4*Power(k,2)*p*ReA13*ReB02*t - 24*k*kprime*p*ReA13*ReB02*t + 
      12*k*Power(p,2)*ReA13*ReB02*t - 4*Power(k,3)*ReB13*t - 
      12*Power(k,2)*kprime*ReB13*t + 16*Power(k,2)*p*ReB13*t + 
      12*k*kprime*p*ReB13*t - 12*k*Power(p,2)*ReB13*t + 
      4*Power(k,3)*ReA02*ReB13*t + 
      12*Power(k,2)*kprime*ReA02*ReB13*t - 
      16*Power(k,2)*p*ReA02*ReB13*t - 12*k*kprime*p*ReA02*ReB13*t + 
      12*k*Power(p,2)*ReA02*ReB13*t - 4*Power(k,2)*ReB02*ReB13*t - 
      12*k*kprime*ReB02*ReB13*t + 12*k*p*ReB02*ReB13*t + 
      12*Power(k,2)*ReC13*t + 12*k*kprime*ReC13*t - 
      24*k*p*ReC13*t - 12*Power(k,2)*ReA02*ReC13*t - 
      12*k*kprime*ReA02*ReC13*t + 24*k*p*ReA02*ReC13*t - 
      6*kprime*ReA13*ReB02*s*t + 6*p*ReA13*ReB02*s*t + 
      6*k*ReB13*s*t - 6*p*ReB13*s*t - 6*k*ReA02*ReB13*s*t + 
      6*p*ReA02*ReB13*s*t + 6*ReB02*ReB13*s*t - 6*ReC13*s*t + 
      6*ReA02*ReC13*s*t - 3*kprime*ReA13*ReB02*Power(t,2) + 
      3*p*ReA13*ReB02*Power(t,2) + 3*k*ReB13*Power(t,2) - 
      3*p*ReB13*Power(t,2) - 3*k*ReA02*ReB13*Power(t,2) + 
      3*p*ReA02*ReB13*Power(t,2) + 3*ReB02*ReB13*Power(t,2) - 
      6*ReC13*Power(t,2) + 6*ReA02*ReC13*Power(t,2) + 
      ImB02*(ImB13 + ImA13*(-kprime + p))*
       (4*Power(k,2)*(2*Power(kprime,2) - 4*kprime*p + 
            2*Power(p,2) - s - t) - 12*k*(kprime - p)*(s + t) + 
         3*Power(s + t,2)) + 
      ImA02*(-(ImB13*(k - p)*
            (4*Power(k,2)*(2*Power(kprime,2) - 4*kprime*p + 
                 2*Power(p,2) - s - t) - 
              12*k*(kprime - p)*(s + t) + 3*Power(s + t,2))) + 
         2*ImC13*(8*Power(k,3)*(kprime - p) + 3*t*(s + t) - 
            2*Power(k,2)*(4*kprime*p - 4*Power(p,2) + 3*(s + t)) + 
            6*k*(-(kprime*t) + p*(s + 2*t)))))/(8.*Power(k,2));
    double Im_Q02primestar_dot_Sigma1Tilde = (-8*ImB02*Power(k,2)*Power(kprime,3)*ReA13 + 
     24*ImB02*Power(k,2)*Power(kprime,2)*p*ReA13 - 
     24*ImB02*Power(k,2)*kprime*Power(p,2)*ReA13 + 
     8*ImB02*Power(k,2)*Power(p,3)*ReA13 + 
     8*ImA13*Power(k,2)*Power(kprime,3)*ReB02 - 
     24*ImA13*Power(k,2)*Power(kprime,2)*p*ReB02 + 
     24*ImA13*Power(k,2)*kprime*Power(p,2)*ReB02 - 
     8*ImA13*Power(k,2)*Power(p,3)*ReB02 + 
     8*ImB02*Power(k,2)*Power(kprime,2)*ReB13 - 
     8*ImA02*Power(k,3)*Power(kprime,2)*ReB13 - 
     16*ImB02*Power(k,2)*kprime*p*ReB13 + 
     16*ImA02*Power(k,3)*kprime*p*ReB13 + 
     8*ImA02*Power(k,2)*Power(kprime,2)*p*ReB13 + 
     8*ImB02*Power(k,2)*Power(p,2)*ReB13 - 
     8*ImA02*Power(k,3)*Power(p,2)*ReB13 - 
     16*ImA02*Power(k,2)*kprime*Power(p,2)*ReB13 + 
     8*ImA02*Power(k,2)*Power(p,3)*ReB13 + 
     16*ImA02*Power(k,3)*kprime*ReC13 - 
     16*ImA02*Power(k,3)*p*ReC13 - 
     16*ImA02*Power(k,2)*kprime*p*ReC13 + 
     16*ImA02*Power(k,2)*Power(p,2)*ReC13 + 
     4*ImB02*Power(k,2)*kprime*ReA13*s + 
     12*ImB02*k*Power(kprime,2)*ReA13*s - 
     4*ImB02*Power(k,2)*p*ReA13*s - 24*ImB02*k*kprime*p*ReA13*s + 
     12*ImB02*k*Power(p,2)*ReA13*s - 
     4*ImA13*Power(k,2)*kprime*ReB02*s - 
     12*ImA13*k*Power(kprime,2)*ReB02*s + 
     4*ImA13*Power(k,2)*p*ReB02*s + 24*ImA13*k*kprime*p*ReB02*s - 
     12*ImA13*k*Power(p,2)*ReB02*s - 4*ImB02*Power(k,2)*ReB13*s + 
     4*ImA02*Power(k,3)*ReB13*s - 12*ImB02*k*kprime*ReB13*s + 
     12*ImA02*Power(k,2)*kprime*ReB13*s + 12*ImB02*k*p*ReB13*s - 
     16*ImA02*Power(k,2)*p*ReB13*s - 12*ImA02*k*kprime*p*ReB13*s + 
     12*ImA02*k*Power(p,2)*ReB13*s - 12*ImA02*Power(k,2)*ReC13*s + 
     12*ImA02*k*p*ReC13*s - 3*ImB02*kprime*ReA13*Power(s,2) + 
     3*ImB02*p*ReA13*Power(s,2) + 3*ImA13*kprime*ReB02*Power(s,2) - 
     3*ImA13*p*ReB02*Power(s,2) + 3*ImB02*ReB13*Power(s,2) - 
     3*ImA02*k*ReB13*Power(s,2) + 3*ImA02*p*ReB13*Power(s,2) + 
     4*ImB02*Power(k,2)*kprime*ReA13*t + 
     12*ImB02*k*Power(kprime,2)*ReA13*t - 
     4*ImB02*Power(k,2)*p*ReA13*t - 24*ImB02*k*kprime*p*ReA13*t + 
     12*ImB02*k*Power(p,2)*ReA13*t - 
     4*ImA13*Power(k,2)*kprime*ReB02*t - 
     12*ImA13*k*Power(kprime,2)*ReB02*t + 
     4*ImA13*Power(k,2)*p*ReB02*t + 24*ImA13*k*kprime*p*ReB02*t - 
     12*ImA13*k*Power(p,2)*ReB02*t - 4*ImB02*Power(k,2)*ReB13*t + 
     4*ImA02*Power(k,3)*ReB13*t - 12*ImB02*k*kprime*ReB13*t + 
     12*ImA02*Power(k,2)*kprime*ReB13*t + 12*ImB02*k*p*ReB13*t - 
     16*ImA02*Power(k,2)*p*ReB13*t - 12*ImA02*k*kprime*p*ReB13*t + 
     12*ImA02*k*Power(p,2)*ReB13*t - 12*ImA02*Power(k,2)*ReC13*t - 
     12*ImA02*k*kprime*ReC13*t + 24*ImA02*k*p*ReC13*t - 
     6*ImB02*kprime*ReA13*s*t + 6*ImB02*p*ReA13*s*t + 
     6*ImA13*kprime*ReB02*s*t - 6*ImA13*p*ReB02*s*t + 
     6*ImB02*ReB13*s*t - 6*ImA02*k*ReB13*s*t + 
     6*ImA02*p*ReB13*s*t + 6*ImA02*ReC13*s*t - 
     3*ImB02*kprime*ReA13*Power(t,2) + 3*ImB02*p*ReA13*Power(t,2) + 
     3*ImA13*kprime*ReB02*Power(t,2) - 3*ImA13*p*ReB02*Power(t,2) + 
     3*ImB02*ReB13*Power(t,2) - 3*ImA02*k*ReB13*Power(t,2) + 
     3*ImA02*p*ReB13*Power(t,2) + 6*ImA02*ReC13*Power(t,2) + 
     ImB13*(p + k*(-1 + ReA02) - p*ReA02 - ReB02)*
      (4*Power(k,2)*(2*Power(kprime,2) - 4*kprime*p + 
           2*Power(p,2) - s - t) - 12*k*(kprime - p)*(s + t) + 
        3*Power(s + t,2)) - 
     2*ImC13*(-1 + ReA02)*(8*Power(k,3)*(kprime - p) + 
        3*t*(s + t) - 2*Power(k,2)*
         (4*kprime*p - 4*Power(p,2) + 3*(s + t)) + 
        6*k*(-(kprime*t) + p*(s + 2*t))))/(8.*Power(k,2));

    double trace2_vis_Re_temp1 = 2.*(Repart_ComplexDivide(Re_Q02primestar_dot_Sigma1primestar, Im_Q02primestar_dot_Sigma1primestar, Re_Q02primestar_dot_Q02primestar, Im_Q02primestar_dot_Q02primestar) 
                                   + Repart_ComplexDivide(Re_Q03Tilde_dot_Sigma1Tilde, Im_Q03Tilde_dot_Sigma1Tilde, Re_Q03Tilde_dot_Q03Tilde, Im_Q03Tilde_dot_Q03Tilde));
    double trace2_vis_Im_temp1 = 2.*(Impart_ComplexDivide(Re_Q02primestar_dot_Sigma1primestar, Im_Q02primestar_dot_Sigma1primestar, Re_Q02primestar_dot_Q02primestar, Im_Q02primestar_dot_Q02primestar) 
                                   + Impart_ComplexDivide(Re_Q03Tilde_dot_Sigma1Tilde, Im_Q03Tilde_dot_Sigma1Tilde, Re_Q03Tilde_dot_Q03Tilde, Im_Q03Tilde_dot_Q03Tilde));
    double trace2_vis_Re_temp2 = Repart_ComplexMultiply(Re_Q03Tilde_dot_Q03Tilde, Im_Q03Tilde_dot_Q03Tilde, Re_Q02primestar_dot_Q02primestar, Im_Q02primestar_dot_Q02primestar);
    double trace2_vis_Im_temp2 = Impart_ComplexMultiply(Re_Q03Tilde_dot_Q03Tilde, Im_Q03Tilde_dot_Q03Tilde, Re_Q02primestar_dot_Q02primestar, Im_Q02primestar_dot_Q02primestar);

    double trace2_vis_1 = Repart_ComplexMultiply(Re_trace2_eq, Im_trace2_eq, trace2_vis_Re_temp1, trace2_vis_Im_temp1);
    double trace2_vis_2 = Repart_ComplexDivide(Re_Q03Tilde_dot_Sigma1primestar, Im_Q03Tilde_dot_Sigma1primestar, trace2_vis_Re_temp2, trace2_vis_Im_temp2) 
                        + Repart_ComplexDivide(Re_Q02primestar_dot_Sigma1Tilde, Im_Q02primestar_dot_Sigma1Tilde, trace2_vis_Re_temp2, trace2_vis_Im_temp2);
    trace2_vis = -64.*s/2.*(trace2_vis_1 - trace2_vis_2);

    result_ptr[1] = prefactor*(trace1_vis + trace2_vis + trace4_vis);
    
    return;
}

void get_quark_selfenergy_coefficients(double p_0, double p_i, double T, Selfenergy_coefficients* Sigma_ptr)
{
    double omega_0_sq = g_s*g_s*C_F*T*T/8.0;

    double p0_sq = p_0*p_0;
    double p_sq = p_i*p_i;
    double p0_cubic = p0_sq*p_0;
    double p_cubic = p_sq*p_i;
    double Re_Qx = quark_selfenergy_Q(p_0/p_i);
    double Im_Qx = - 0.5*(p_0/p_i)*M_PI;
    
    //Equilibrium quark self energy from Hard Thermal Loop Approximation
    Sigma_ptr->Re_A0 = omega_0_sq/p_sq*(Re_Qx - 1.); 
    Sigma_ptr->Re_B0 = omega_0_sq*((- p_0/p_sq + 1./p_0)*Re_Qx + p_0/p_sq);

    double neq_coeff = (4./(M_PI*M_PI))*12.62159748;  // for delta_alpha = 2
    //Off-equilibrium corrections coefficients from Hard Loop Approximation
    Sigma_ptr->Re_A1 = omega_0_sq*neq_coeff*1./(2.*p_cubic*p_cubic)
                  *((5.*p0_sq - 3.*p_sq)*Re_Qx - 5.*p0_sq + 4./3.*p_sq);
    Sigma_ptr->Re_B1 = omega_0_sq*neq_coeff*1./(2.*p_cubic*p_cubic)
                  *((- 5.*p0_cubic + 6.*p_0*p_sq - p_cubic*p_i/p_0)*Re_Qx
                     + 5.*p0_cubic - 13./3.*p_0*p_sq);
    Sigma_ptr->Re_C1 = omega_0_sq*neq_coeff*1./(2.*p_sq*p_sq)
                  *((p0_sq - p_sq)*Re_Qx - p0_sq + 2./3.*p_sq);
    
    //imaginary part
    if(p_0 > p_i || p_0 < - p_i)
    {
       Sigma_ptr->Im_A0 = 0.0e0;
       Sigma_ptr->Im_B0 = 0.0e0;
       Sigma_ptr->Im_A1 = 0.0e0;
       Sigma_ptr->Im_B1 = 0.0e0;
       Sigma_ptr->Im_C1 = 0.0e0;
    }
    else
    {
       Sigma_ptr->Im_A0 = omega_0_sq/p_sq*Im_Qx; 
       Sigma_ptr->Im_B0 = omega_0_sq*(- p_0/p_sq + 1./p_0)*Im_Qx;
       Sigma_ptr->Im_A1 = omega_0_sq*neq_coeff*1./(2.*p_cubic*p_cubic)
                          *(5.*p0_sq - 3.*p_sq)*Im_Qx;
       Sigma_ptr->Im_B1 = omega_0_sq*neq_coeff*1./(2.*p_cubic*p_cubic)
                          *(- 5.*p0_cubic + 6.*p_0*p_sq - p_cubic*p_i/p_0)*Im_Qx;
       Sigma_ptr->Im_C1 = omega_0_sq*neq_coeff*1./(2.*p_sq*p_sq)
                          *(p0_sq - p_sq)*Im_Qx;
    }

    return;
}

double quark_selfenergy_Q(double x)
{
    double result; 
    double eps = 1e-100;
    if(fabs(x) > 1.)
      result = 0.5*x*log((x + 1. - eps)/(x - 1. + eps));
    else
      result = 0.5*x*log((x + 1. + eps)/(1. - x + eps));

    return(result);
}

inline double Repart_ComplexMultiply(double Re1, double Im1, double Re2, double Im2)
{
    return(Re1*Re2 - Im1*Im2);
}

inline double Impart_ComplexMultiply(double Re1, double Im1, double Re2, double Im2)
{
    return(Re1*Im2 + Im1*Re2);
}

inline double Repart_ComplexDivide(double Re1, double Im1, double Re2, double Im2)
{
    return((Re1*Re2 + Im1*Im2)/(Re2*Re2 + Im2*Im2));
}

inline double Impart_ComplexDivide(double Re1, double Im1, double Re2, double Im2)
{
    return((-Re1*Im2 + Im1*Re2)/(Re2*Re2 + Im2*Im2));
}

inline double Power(double x, int a)
{
    double result = 1.0;
    for(int i=0; i<a; i++)
       result *= x;
    return(result);
}
