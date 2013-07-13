#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "QGP_2to2_Scattering_Kinetic.h"
#include "Physicalconstants.h"
#include "ParameterReader.h"
#include "gauss_quadrature.h"
#include "Arsenal.h"

using namespace std;

QGP_2to2_Scattering_Kinetic::QGP_2to2_Scattering_Kinetic(ParameterReader* paraRdr_in)
{
   paraRdr = paraRdr_in;
   
   n_Eq = paraRdr->getVal("n_Eq");
   double Eq_i = paraRdr->getVal("Eq_min");
   double dEq = paraRdr->getVal("dEq");
   Eq_tb = new double [n_Eq];
   for(int i=0; i<n_Eq; i++)
      Eq_tb[i] = Eq_i + i*dEq;

   n_Temp = paraRdr->getVal("n_Temp");
   double T_i = paraRdr->getVal("T_min");
   double dT = paraRdr->getVal("dT");
   T_tb = new double [n_Temp];
   for(int i=0; i<n_Temp; i++)
      T_tb[i] = T_i + i*dT;

   equilibrium_results = new double* [n_Eq];
   viscous_results = new double*[n_Eq];
   for(int i=0; i<n_Eq; i++)
   {
      equilibrium_results[i] = new double [n_Temp];
      viscous_results[i] = new double [n_Temp];
   }

   n_qtilde = paraRdr->getVal("n_qtilde");
   double qtilde_i = paraRdr->getVal("qtilde_i");
   double qtilde_f = paraRdr->getVal("qtilde_max");
   double dqtilde = (qtilde_f - qtilde_i)/(n_qtilde - 1 + 1e-10);
   qtilde_pt = new double [n_qtilde];
   for(int i=0; i<n_qtilde; i++)
      qtilde_pt[i] = qtilde_i + i*dqtilde;
   equilibriumTilde_results = new double [n_qtilde];
   viscousTilde_results = new double [n_qtilde];
   
   //initialize the Gaussian quadrature lattices
   n_s = paraRdr->getVal("n_s");
   s_pt = new double [n_s];
   s_weight = new double [n_s];
   s_pt_standard = new double [n_s];
   s_weight_standard = new double [n_s];

   t_pt = new double* [n_s];
   t_weight = new double* [n_s];
   t_pt_standard = new double* [n_s];
   t_weight_standard = new double* [n_s];
   Matrix_elements_sq_ptr = new double* [n_s];
   n_t = paraRdr->getVal("n_t");
   for(int i=0; i<n_s; i++)
   {
       t_pt[i] = new double [n_t];
       t_weight[i] = new double [n_t];
       t_pt_standard[i] = new double [n_t];
       t_weight_standard[i] = new double [n_t];
       Matrix_elements_sq_ptr[i] = new double [n_t];
   }

   n_E1 = paraRdr->getVal("n_E1");
   E1_pt_standard = new double [n_E1];
   E1_weight_standard = new double [n_E1];
   n_E2 = paraRdr->getVal("n_E2");
   E2_pt_standard = new double* [3];
   E2_weight_standard = new double* [3];
   for(int i=0; i<3; i++)
   {
      E2_pt_standard[i] = new double [n_E2];
      E2_weight_standard[i] = new double [n_E2];
   }

   deltaf_alpha = paraRdr->getVal("deltaf_alpha");

   set_gausspoints();
   
}

QGP_2to2_Scattering_Kinetic::~QGP_2to2_Scattering_Kinetic()
{
   delete[] Eq_tb;
   delete[] T_tb;
   for(int i=0; i<n_Eq; i++)
   {
      delete[] equilibrium_results[i];
      delete[] viscous_results[i];
   }
   delete[] equilibrium_results;
   delete[] viscous_results;
   delete[] equilibriumTilde_results;
   delete[] viscousTilde_results;

   delete[] qtilde_pt;

   delete[] s_pt;
   delete[] s_weight;
   delete[] s_pt_standard;
   delete[] s_weight_standard;
   for(int i=0; i<n_s; i++)
   {
      delete[] t_pt[i];
      delete[] t_weight[i];
      delete[] t_pt_standard[i];
      delete[] t_weight_standard[i];
      delete[] Matrix_elements_sq_ptr[i];
   }
   delete[] t_pt;
   delete[] t_weight;
   delete[] t_pt_standard;
   delete[] t_weight_standard;
   delete[] Matrix_elements_sq_ptr;
   
   delete[] E1_pt_standard;
   delete[] E1_weight_standard;
   for(int i=0; i<3; i++)
   {
      delete[] E2_pt_standard[i];
      delete[] E2_weight_standard[i];
   }
   delete[] E2_pt_standard;
   delete[] E2_weight_standard;

}

void QGP_2to2_Scattering_Kinetic::set_gausspoints()
{
   gauss_quadrature_standard(n_s, 1, 0.0, 0.0, 0.0, 1.0, s_pt_standard, s_weight_standard);
  
   for(int i=0; i<n_s; i++)
      gauss_quadrature_standard(n_t, 1, 0.0, 0.0, 0.0, 1.0, t_pt_standard[i], t_weight_standard[i]);
    
   gauss_quadrature_standard(n_E1, 5, 0.0, 0.0, 0.0, 1.0, E1_pt_standard, E1_weight_standard);

   // use Chebyshev–Gauss quadrature for channels: pi + rho, pi + Kstar, rho + K, and K + Kstar
   gauss_quadrature_standard(n_E2, 2, 0.0, 0.0, 0.0, 1.0, E2_pt_standard[0], E2_weight_standard[0]);

   // use Jacobi-Gauss quadrature for channels: pi + pi, pi + K
   gauss_quadrature_standard(n_E2, 4, 0.0, -0.5, 0.0, 1.0, E2_pt_standard[1], E2_weight_standard[1]);
   gauss_quadrature_standard(n_E2, 4, -0.5, 0.0, 0.0, 1.0, E2_pt_standard[2], E2_weight_standard[2]);

}

void QGP_2to2_Scattering_Kinetic::calculateEmissionrates(int channel_in, string filename_in)
{
   double hbarC = Phycons.get_hbarC();
   channel = channel_in;
   filename = filename_in;

   double *results = new double [2];

   for(int i=0; i<n_qtilde; i++)
   {
       double qtilde = qtilde_pt[i];

       double prefactor = 1./(16.*pow(2.0*M_PI, 7)*qtilde);

       scale_gausspoints_st(qtilde);

       double equilibrium_result_s = 0.0;
       double viscous_result_s = 0.0;
       for(int k=0; k<n_s; k++)
       {
          double equilibrium_result_t = 0.0;
          double viscous_result_t = 0.0;
          for(int l=0; l<n_t; l++)
          {
             Integrate_E1(qtilde, s_pt[k], t_pt[k][l], results);
             equilibrium_result_t += results[0]*t_weight[k][l];
             viscous_result_t += results[1]*t_weight[k][l];
          }
          equilibrium_result_s += equilibrium_result_t*s_weight[k];
          viscous_result_s += viscous_result_t*s_weight[k];
       }
       equilibriumTilde_results[i] = equilibrium_result_s*prefactor/pow(hbarC, 4); // convert units to 1/(GeV^2 fm^4) for the emission rates

       viscousTilde_results[i] = viscous_result_s*prefactor/(qtilde*qtilde)/pow(hbarC, 4); // convert units to 1/(GeV^2 fm^4) for the emission rates
   }
   
   buildupEmissionrate2DTable();
   output_emissionrateTable();
   delete [] results;
}

void QGP_2to2_Scattering_Kinetic::buildupEmissionrate2DTable()
{
   double *qtildeT = new double [n_Eq];
   double *rawResult_eq = new double [n_Eq];
   double *rawResult_vis = new double [n_Eq];
   double *log_eq = new double [n_qtilde];
   double *scaled_vis = new double [n_qtilde];
   for(int i = 0; i < n_qtilde; i++)
   {
      log_eq[i] = log(equilibriumTilde_results[i]);
      scaled_vis[i] = viscousTilde_results[i]/equilibriumTilde_results[i];
   }
   for(int j = 0; j < n_Temp; j++)
   {
      double T = T_tb[j];
      for(int i = 0; i < n_Eq; i++)
         qtildeT[i] = Eq_tb[i]/T;
      interpolation1D_linear(qtilde_pt, log_eq, qtildeT, rawResult_eq, n_qtilde, n_Eq);
      interpolation1D_linear(qtilde_pt, scaled_vis, qtildeT, rawResult_vis, n_qtilde, n_Eq);
      for(int i = 0; i < n_Eq; i++)
      {
         double temp = exp(rawResult_eq[i]);
         equilibrium_results[i][j] = T*T*temp;
         viscous_results[i][j] = rawResult_vis[i]*temp;
      }
   }
   delete [] qtildeT;
   delete [] rawResult_eq;
   delete [] rawResult_vis;
   delete [] log_eq;
   delete [] scaled_vis;
}

void QGP_2to2_Scattering_Kinetic::output_emissionrateTable()
{
   // output emission rate tables
   ostringstream output_file_eqrate;
   ostringstream output_file_viscous;
   output_file_eqrate << "rate_" << filename << "_eqrate.dat";
   output_file_viscous << "rate_" << filename << "_viscous.dat";
   ofstream of_eqrate(output_file_eqrate.str().c_str());
   ofstream of_viscous(output_file_viscous.str().c_str());
   for(int j=0; j<n_Temp; j++)
   {
      for(int i=0; i<n_Eq; i++)
      {
         of_eqrate << scientific << setw(20) << setprecision(8)
                   << equilibrium_results[i][j] << "   ";
         of_viscous << scientific << setw(20) << setprecision(8)
                   << viscous_results[i][j] << "   ";
      }
      of_eqrate << endl;
      of_viscous << endl;
   }

   of_eqrate.close();
   of_viscous.close();
}

void QGP_2to2_Scattering_Kinetic::scale_gausspoints_st(double qtilde)
{
   p_cutoff = 0.0e0;
   double s_min = 2*p_cutoff*p_cutoff;
   double s_max = paraRdr->getVal("s_max");
   if(s_max < qtilde*200) s_max = qtilde*200;
  
   for(int i=0; i < n_s; i++)
   {
      s_pt[i] = s_pt_standard[i];
      s_weight[i] = s_weight_standard[i];
      for(int j=0; j<n_t; j++)
      {
         t_pt[i][j] = t_pt_standard[i][j];
         t_weight[i][j] = t_weight_standard[i][j];
      }
   }
   scale_gausspoints(n_s, 1, 0.0, 0.0, s_min, s_max, s_pt, s_weight);
  
   for(int i=0; i<n_s; i++)
   {
      double s = s_pt[i];
      double t_min;
      double t_max;
      t_min = - s + p_cutoff*p_cutoff;
      t_max = - p_cutoff*p_cutoff;

      scale_gausspoints(n_t, 1, 0.0, 0.0, t_min, t_max, t_pt[i], t_weight[i]);
    }
}

void QGP_2to2_Scattering_Kinetic::Integrate_E1(double qtilde, double s_tilde, double t_tilde, double* results)
{
   double equilibrium_result = 0.0e0;
   double viscous_result = 0.0e0;
   double E1_min;
   double u_tilde = - s_tilde - t_tilde;
   E1_min = (- u_tilde)/(4.*qtilde);

   double* E1_pt = new double [n_E1];
   double* E1_weight = new double [n_E1];
   for(int i=0; i<n_E1; i++)
   {
      E1_pt[i] = E1_pt_standard[i];
      E1_weight[i] = E1_weight_standard[i];
   }
   
   //double slope = qtilde;
   double slope = 1.0;
   scale_gausspoints(n_E1, 5, 0.0, 0.0, E1_min, slope, E1_pt, E1_weight);

   for(int i=0; i<n_E1; i++)
   {
      Integrate_E2(qtilde, s_tilde, t_tilde, E1_pt[i], results);
      equilibrium_result += results[0]*E1_weight[i];
      viscous_result += results[1]*E1_weight[i];
   }
   
   results[0] = equilibrium_result;
   results[1] = viscous_result;

   delete[] E1_pt;
   delete[] E1_weight;
}

void QGP_2to2_Scattering_Kinetic::Integrate_E2(double qtilde, double s_tilde, double t_tilde, double E1tilde, double* results)
{
   double eps = 1e-100;
   double equilibrium_result = 0.0;
   double viscous_result = 0.0;
   double E2_min;
   double E2_max;
   double min_1 = (- t_tilde)/(4.*qtilde);

   double a = - (s_tilde + t_tilde)*(s_tilde + t_tilde);
   double b = qtilde*((s_tilde + t_tilde)*(s_tilde)) + E1tilde*(-t_tilde)*(s_tilde + t_tilde);
   double c = - (qtilde*s_tilde + E1tilde*t_tilde)*(qtilde*s_tilde + E1tilde*t_tilde) + s_tilde*t_tilde*(s_tilde+t_tilde);

   if((b*b - a*c) >= 0) 
   {
      double min_2 = (-b + sqrt(b*b - a*c))/(a + eps);
      if(min_1 < min_2)
         E2_min = min_2;
      else
         E2_min = min_1;
      E2_max = (-b - sqrt(b*b - a*c))/(a + eps);

      if(E2_max < E2_min)
      {
         results[0] = 0.0e0;
         results[1] = 0.0e0;
         return;
      }
   
      double common_factor;

      double* E2_pt = new double [n_E2];
      double* E2_weight = new double [n_E2];
      
      for(int i=0; i<n_E2; i++)
      {
         E2_pt[i] = E2_pt_standard[0][i];
         E2_weight[i] = E2_weight_standard[0][i];
      }
      scale_gausspoints(n_E2, 2, 0.0, 0.0, E2_min, E2_max, E2_pt, E2_weight);
      for(int i=0; i<n_E2; i++)
      {
         double* Matrix_sq_ptr = new double [2];
         Matrix_elements_sq(s_tilde, t_tilde, E1tilde, E2_pt[i], qtilde, Matrix_sq_ptr);
         double Matrix_sq_eq = Matrix_sq_ptr[0];
         double Matrix_sq_noneq = Matrix_sq_ptr[1];

         double f0_E1, f0_E2, f0_E3;
         if(channel == 1) //Compton Scattering
         {
            f0_E1 = Bose_distribution(E1tilde);
            f0_E2 = Fermi_distribution(E2_pt[i]);
            f0_E3 = Fermi_distribution(E1tilde + E2_pt[i] - qtilde);
            common_factor = f0_E1*f0_E2*(1 - f0_E3)/(sqrt(a*E2_pt[i]*E2_pt[i] + 2*b*E2_pt[i] + c) + eps);
         }
         else if(channel ==2) //pair annihilation
         {
            f0_E1 = Fermi_distribution(E1tilde);
            f0_E2 = Fermi_distribution(E2_pt[i]);
            f0_E3 = Bose_distribution(E1tilde + E2_pt[i] - qtilde);
            common_factor = f0_E1*f0_E2*(1 + f0_E3)/(sqrt(a*E2_pt[i]*E2_pt[i] + 2*b*E2_pt[i] + c) + eps);
         }
         equilibrium_result += common_factor*1.*Matrix_sq_eq*E2_weight[i];
         viscous_result += common_factor*(Matrix_sq_eq*viscous_integrand(s_tilde, t_tilde, E1tilde, E2_pt[i], qtilde, f0_E1, f0_E2, f0_E3) + 1.*Matrix_sq_noneq)*E2_weight[i];
         delete [] Matrix_sq_ptr;
      }

      delete[] E2_pt;
      delete[] E2_weight;
   }
   else  // no kinematic phase space
   {
      equilibrium_result = 0.0e0;
      viscous_result = 0.0e0;
   }
   results[0] = equilibrium_result;
   results[1] = viscous_result;
}

double QGP_2to2_Scattering_Kinetic::viscous_integrand(double s_tilde, double t_tilde, double E1tilde, double E2tilde, double qtilde, double f0_E1, double f0_E2, double f0_E3)
{
   double eps = 1e-100;
   double E3tilde = E1tilde + E2tilde - qtilde;
   double p1tilde = E1tilde;
   double p2tilde = E2tilde;
   double p3tilde = E3tilde;
   double costheta1 = (- s_tilde - t_tilde + 2*E1tilde*qtilde)/(2*p1tilde*qtilde + eps);
   double costheta2 = (t_tilde + 2*E2tilde*qtilde)/(2*p2tilde*qtilde + eps);
   double p3_z_tilde = p1tilde*costheta1 + p2tilde*costheta2 - qtilde; 
   double integrand;
   if(channel == 1)
   {
      integrand = (1. + f0_E1)*deltaf_chi(p1tilde)*0.5*(-1. + 3.*costheta1*costheta1) + (1. - f0_E2)*deltaf_chi(p2tilde)*0.5*(-1. + 3.*costheta2*costheta2) - f0_E3*deltaf_chi(p3tilde)/p3tilde/p3tilde*(-0.5*p3tilde*p3tilde + 1.5*p3_z_tilde*p3_z_tilde);
   }
   else if(channel == 2)
   {
      integrand = (1. - f0_E1)*deltaf_chi(p1tilde)*0.5*(-1. + 3.*costheta1*costheta1) + (1. - f0_E2)*deltaf_chi(p2tilde)*0.5*(-1. + 3.*costheta2*costheta2) + f0_E3*deltaf_chi(p3tilde)/p3tilde/p3tilde*(-0.5*p3tilde*p3tilde + 1.5*p3_z_tilde*p3_z_tilde);
   }
   else
   {
      cout << "viscous_integrand:: channel can not be found! (channel = " << channel << ")" << endl;
      exit(1);
   }

   return(integrand);
}

void QGP_2to2_Scattering_Kinetic::get_quark_selfenergy_coefficients(double p_0_tilde, double p_i_tilde, Selfenergy_coefficients* Sigma_ptr)
{
    double g_s = Phycons.get_g_s_const();
    double C_F = Phycons.get_C_F();

    double omega_0_sq = g_s*g_s*C_F/8.0;

    double p0_tilde_sq = p_0_tilde*p_0_tilde;
    double p_tilde_sq = p_i_tilde*p_i_tilde;
    double p0_tilde_cubic = p0_tilde_sq*p_0_tilde;
    double p_tilde_cubic = p_tilde_sq*p_i_tilde;
    double Re_Qx = quark_selfenergy_Q(p_0_tilde/p_i_tilde);
    double Im_Qx = - 0.5*(p_0_tilde/p_i_tilde)*M_PI;
    
    //Equilibrium quark self energy from Hard Thermal Loop Approximation
    Sigma_ptr->Re_A0 = omega_0_sq/p_tilde_sq*(Re_Qx - 1.); 
    Sigma_ptr->Re_B0 = omega_0_sq*((- p_0_tilde/p_tilde_sq + 1./p_0_tilde)*Re_Qx + p_0_tilde/p_tilde_sq);

    double neq_coeff = (4./(M_PI*M_PI))*12.62159748;  // for delta_alpha = 2
    //Off-equilibrium corrections coefficients from Hard Loop Approximation
    Sigma_ptr->Re_A1 = omega_0_sq*neq_coeff*1./(2.*p_tilde_cubic*p_tilde_cubic)
                  *((5.*p0_tilde_sq - 3.*p_tilde_sq)*Re_Qx - 5.*p0_tilde_sq + 4./3.*p_tilde_sq);
    Sigma_ptr->Re_B1 = omega_0_sq*neq_coeff*1./(2.*p_tilde_cubic*p_tilde_cubic)
                  *((- 5.*p0_tilde_cubic + 6.*p_0_tilde*p_tilde_sq - p_tilde_cubic*p_i_tilde/p_0_tilde)*Re_Qx
                     + 5.*p0_tilde_cubic - 13./3.*p_0_tilde*p_tilde_sq);
    Sigma_ptr->Re_C1 = omega_0_sq*neq_coeff*1./(2.*p_tilde_sq*p_tilde_sq)
                  *((p0_tilde_sq - p_tilde_sq)*Re_Qx - p0_tilde_sq + 2./3.*p_tilde_sq);
    
    //imaginary part
    if(p_0_tilde > p_i_tilde || p_0_tilde < - p_i_tilde)
    {
       Sigma_ptr->Im_A0 = 0.0e0;
       Sigma_ptr->Im_B0 = 0.0e0;
       Sigma_ptr->Im_A1 = 0.0e0;
       Sigma_ptr->Im_B1 = 0.0e0;
       Sigma_ptr->Im_C1 = 0.0e0;
    }
    else
    {
       Sigma_ptr->Im_A0 = omega_0_sq/p_tilde_sq*Im_Qx; 
       Sigma_ptr->Im_B0 = omega_0_sq*(- p_0_tilde/p_tilde_sq + 1./p_0_tilde)*Im_Qx;
       Sigma_ptr->Im_A1 = omega_0_sq*neq_coeff*1./(2.*p_tilde_cubic*p_tilde_cubic)
                          *(5.*p0_tilde_sq - 3.*p_tilde_sq)*Im_Qx;
       Sigma_ptr->Im_B1 = omega_0_sq*neq_coeff*1./(2.*p_tilde_cubic*p_tilde_cubic)
                          *(- 5.*p0_tilde_cubic + 6.*p_0_tilde*p_tilde_sq - p_tilde_cubic*p_i_tilde/p_0_tilde)*Im_Qx;
       Sigma_ptr->Im_C1 = omega_0_sq*neq_coeff*1./(2.*p_tilde_sq*p_tilde_sq)
                          *(p0_tilde_sq - p_tilde_sq)*Im_Qx;
    }
    return;
}

double QGP_2to2_Scattering_Kinetic::quark_selfenergy_Q(double x)
{
    double result; 
    double eps = 1e-100;
    if(fabs(x) > 1.)
      result = 0.5*x*log((x + 1. - eps)/(x - 1. + eps));
    else
      result = 0.5*x*log((x + 1. + eps)/(1. - x + eps));

    return(result);
}

double QGP_2to2_Scattering_Kinetic::Bose_distribution(double Etilde)
{
   return(1.0/(exp(Etilde)-1.0));
}

double QGP_2to2_Scattering_Kinetic::Fermi_distribution(double Etilde)
{
    return(1.0/(exp(Etilde) + 1.0));
}

double QGP_2to2_Scattering_Kinetic::deltaf_chi(double ptilde)
{ 
    return(pow(ptilde, deltaf_alpha));
}

void QGP_2to2_Scattering_Kinetic::Matrix_elements_sq(double s, double t, double E1, double E2, double qtilde, double* result_ptr)
{
   switch(channel)
   {
      case 1:  //Compton scattering
         Matrix_elements_sq_Compton(s, t, E1, E2, qtilde, result_ptr);
         break;
      case 2:  //pair annihilation
         Matrix_elements_sq_Annihilation(s, t, E1, E2, qtilde, result_ptr);
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
void QGP_2to2_Scattering_Kinetic::Matrix_elements_sq_Compton(double s, double t, double E1, double E2, double qtilde, double* result_ptr)
{
    double e_sq = Phycons.get_e_sq();
    double q_sq = Phycons.get_q_sq();
    double d_F = Phycons.get_d_F();
    double C_F = Phycons.get_C_F();
    double g_s = Phycons.get_g_s_const();

    double prefactor = 2.*e_sq*q_sq*d_F*C_F*g_s*g_s; //2 counts for fermions and anti-fermions

    double pprime = E1;
    double p = E2;
    double k = qtilde;
    double kprime = p + pprime - k;

    double p_minus_k_0 = p - k;
    double p_minus_k_i = sqrt(p_minus_k_0*p_minus_k_0 - t);

    //Calculate self-energy coefficients
    Selfenergy_coefficients Sigma02;
    Selfenergy_coefficients* Sigma02ptr = &Sigma02;
    get_quark_selfenergy_coefficients(p_minus_k_0, p_minus_k_i, Sigma02ptr);
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
void QGP_2to2_Scattering_Kinetic::Matrix_elements_sq_Annihilation(double s, double t, double E1, double E2, double qtilde, double* result_ptr)
{
    double u = - s - t;

    double e_sq = Phycons.get_e_sq();
    double q_sq = Phycons.get_q_sq();
    double d_F = Phycons.get_d_F();
    double C_F = Phycons.get_C_F();
    double g_s = Phycons.get_g_s_const();
    
    double prefactor = e_sq*q_sq*d_F*C_F*g_s*g_s;
    
    double pprime = E1;
    double p = E2;
    double k = qtilde;
    double kprime = p + pprime - k;

    double p_minus_k_0 = p - k;
    double p_minus_k_i = sqrt(p_minus_k_0*p_minus_k_0 - t);
    double p_minus_kprime_0 = p - kprime;
    double p_minus_kprime_i = sqrt(p_minus_kprime_0*p_minus_kprime_0 - u);

    //Calculate self-energy coefficients
    Selfenergy_coefficients Sigma02;
    Selfenergy_coefficients* Sigma02ptr = &Sigma02;
    get_quark_selfenergy_coefficients(p_minus_k_0, p_minus_k_i, Sigma02ptr);
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
    get_quark_selfenergy_coefficients(p_minus_kprime_0, p_minus_kprime_i, Sigma03ptr);
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
    
    result_ptr[0] = prefactor*(2.0*trace1_eq + trace2_eq + 0.0*trace4_eq);
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

/*
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
    */
    trace2_vis = 0.0;

    result_ptr[1] = prefactor*(2.0*trace1_vis + trace2_vis + 0.0*trace4_vis);
    
    return;
}


double QGP_2to2_Scattering_Kinetic::Repart_ComplexMultiply(double Re1, double Im1, double Re2, double Im2)
{
    return(Re1*Re2 - Im1*Im2);
}

double QGP_2to2_Scattering_Kinetic::Impart_ComplexMultiply(double Re1, double Im1, double Re2, double Im2)
{
    return(Re1*Im2 + Im1*Re2);
}

double QGP_2to2_Scattering_Kinetic::Repart_ComplexDivide(double Re1, double Im1, double Re2, double Im2)
{
    return((Re1*Re2 + Im1*Im2)/(Re2*Re2 + Im2*Im2));
}

double QGP_2to2_Scattering_Kinetic::Impart_ComplexDivide(double Re1, double Im1, double Re2, double Im2)
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
