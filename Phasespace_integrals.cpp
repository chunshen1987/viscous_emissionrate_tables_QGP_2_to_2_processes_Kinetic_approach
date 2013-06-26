#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Phasespace_integrals.h"
#include "gauss_quadrature.h"

using namespace std;

int Calculate_emissionrates(int channel, string filename)
{
   double* results = new double [3];
   double* s_pt = new double [n_s];
   double* s_weight = new double [n_s];

   double** t_pt = new double* [n_s];
   double** t_weight = new double* [n_s];
   for(int i=0; i<n_s; i++)
   {
       t_pt[i] = new double [n_t];
       t_weight[i] = new double [n_t];
   }
   
   double* E1_pt_standard = new double [n_E1];
   double* E1_weight_standard = new double [n_E1];
   double** E2_pt_standard = new double* [3];
   double** E2_weight_standard = new double* [3];
   for(int i=0; i<3; i++)
   {
      E2_pt_standard[i] = new double [n_E2];
      E2_weight_standard[i] = new double [n_E2];
   }

   set_gausspoints(s_pt, s_weight, t_pt, t_weight, E1_pt_standard, E1_weight_standard, E2_pt_standard, E2_weight_standard);
   
   double* Eq_tb = new double [n_Eq];
   double* T_tb = new double [n_Temp];
   for(int i=0; i<n_Eq; i++)
      Eq_tb[i] = Eq_i + i*dEq;
   for(int i=0; i<n_Temp; i++)
      T_tb[i] = T_i + i*dT;

   double** equilibrium_results = new double* [n_Eq];
   double** viscous_results1 = new double*[n_Eq];
   double** viscous_results2 = new double*[n_Eq];
   for(int i=0; i<n_Eq; i++)
   {
      equilibrium_results[i] = new double [n_Temp];
      viscous_results1[i] = new double [n_Temp];
      viscous_results2[i] = new double [n_Temp];
   }

   double Eq;
   double T;
   for(int j=0; j<n_Temp; j++)
   {
      T = T_tb[j];
      for(int i=0; i<n_Eq; i++)
      {
          Eq = Eq_tb[i];
          double prefactor = 1./(16.*pow(2.0*M_PI, 7)*Eq);

          double equilibrium_result_s = 0.0;
          double viscous_result1_s = 0.0;
          double viscous_result2_s = 0.0;
          for(int k=0; k<n_s; k++)
          {
             double equilibrium_result_t = 0.0;
             double viscous_result1_t = 0.0;
             double viscous_result2_t = 0.0;
             for(int l=0; l<n_t; l++)
             {
                Integrate_E1(Eq, T, channel, s_pt[k], t_pt[k][l], E1_pt_standard, E1_weight_standard, E2_pt_standard, E2_weight_standard, results);
                equilibrium_result_t += results[0]*t_weight[k][l];
                viscous_result1_t += results[1]*t_weight[k][l];
                viscous_result2_t += results[2]*t_weight[k][l];
             }
             equilibrium_result_s += equilibrium_result_t*s_weight[k];
             viscous_result1_s += viscous_result1_t*s_weight[k];
             viscous_result2_s += viscous_result2_t*s_weight[k];
          }
          equilibrium_results[i][j] = equilibrium_result_s*prefactor/pow(hbarC, 4); // convert units to 1/(GeV^2 fm^4) for the emission rates
          viscous_results1[i][j] = viscous_result1_s*prefactor/(Eq*Eq)/pow(hbarC, 4); // convert units to 1/(GeV^4 fm^4) for the emission rates
          viscous_results2[i][j] = viscous_result2_s*prefactor/(Eq*Eq)/pow(hbarC, 4); // convert units to 1/(GeV^4 fm^4) for the emission rates
          cout << Eq << "  " << T << "  " << equilibrium_results[i][j] << "   " << viscous_results1[i][j]*1. + 1.*viscous_results2[i][j]<< endl;
          exit(0);
      }
   }

   // output emission rate tables
   ostringstream output_file_eqrate;
   ostringstream output_file_viscous1;
   ostringstream output_file_viscous2;
   output_file_eqrate << "rate_" << filename << "_eqrate.dat";
   output_file_viscous1 << "rate_" << filename << "_viscous.dat";
   output_file_viscous2 << "rate_" << filename << "_viscous2.dat";
   ofstream of_eqrate(output_file_eqrate.str().c_str());
   ofstream of_viscous1(output_file_viscous1.str().c_str());
   ofstream of_viscous2(output_file_viscous2.str().c_str());
   for(int j=0; j<n_Temp; j++)
   {
      for(int i=0; i<n_Eq; i++)
      {
         of_eqrate << scientific << setw(20) << setprecision(8)
                   << equilibrium_results[i][j] << "   ";
         of_viscous1 << scientific << setw(20) << setprecision(8)
                   << viscous_results1[i][j] << "   ";
         of_viscous2 << scientific << setw(20) << setprecision(8)
                   << viscous_results2[i][j] << "   ";
      }
      of_eqrate << endl;
      of_viscous1 << endl;
      of_viscous2 << endl;
   }

   of_eqrate.close();
   of_viscous1.close();
   of_viscous2.close();

   ofstream check("check.dat");
   for(int i=0; i<n_Eq; i++)
      check << Eq_tb[i] << "   " << viscous_results2[i][0] << endl;
   check.close();

   delete[] results;
   delete[] s_pt;
   delete[] s_weight;
   for(int i=0; i<n_s; i++)
   {
      delete[] t_pt[i];
      delete[] t_weight[i];
   }
   delete[] t_pt;
   delete[] t_weight;
   
   delete[] E1_pt_standard;
   delete[] E1_weight_standard;
   for(int i=0; i<3; i++)
   {
      delete[] E2_pt_standard[i];
      delete[] E2_weight_standard[i];
   }
   delete[] E2_pt_standard;
   delete[] E2_weight_standard;

   delete[] Eq_tb;
   delete[] T_tb;
   for(int i=0; i<n_Eq; i++)
   {
      delete[] equilibrium_results[i];
      delete[] viscous_results1[i];
      delete[] viscous_results2[i];
   }
   delete[] equilibrium_results;
   delete[] viscous_results1;
   delete[] viscous_results2;

   return(0);
}

double set_gausspoints(double* s_pt, double* s_weight, double** t_pt, double** t_weight, double* E1_pt_standard, double* E1_weight_standard, double** E2_pt_standard, double** E2_weight_standard)
{
   double s_min = 2*q_cutoff*q_cutoff;
  
   gauss_quadrature(n_s, 1, 0.0, 0.0, s_min, s_max, s_pt, s_weight);
  
   for(int i=0; i<n_s; i++)
   {
      double s = s_pt[i];
      double t_min;
      double t_max;
      t_min = - s + q_cutoff*q_cutoff;
      t_max = - q_cutoff*q_cutoff;

      gauss_quadrature(n_t, 1, 0.0, 0.0, t_min, t_max, t_pt[i], t_weight[i]);
    }
    
    gauss_quadrature_standard(n_E1, 5, 0.0, 0.0, 0.0, 1.0, E1_pt_standard, E1_weight_standard);

    // use Chebyshev–Gauss quadrature for channels: pi + rho, pi + Kstar, rho + K, and K + Kstar
    gauss_quadrature_standard(n_E2, 2, 0.0, 0.0, 0.0, 1.0, E2_pt_standard[0], E2_weight_standard[0]);

    // use Jacobi-Gauss quadrature for channels: pi + pi, pi + K
    gauss_quadrature_standard(n_E2, 4, 0.0, -0.5, 0.0, 1.0, E2_pt_standard[1], E2_weight_standard[1]);
    gauss_quadrature_standard(n_E2, 4, -0.5, 0.0, 0.0, 1.0, E2_pt_standard[2], E2_weight_standard[2]);

    return(0);

}

double Integrate_E1(double Eq, double T, int channel, double s, double t, double* E1_pt_standard, double* E1_weight_standard, double** E2_pt_standard, double** E2_weight_standard, double* results)
{
   double equilibrium_result = 0.0e0;
   double viscous_result1 = 0.0e0;
   double viscous_result2 = 0.0e0;
   double u = - s - t;
   double E1_min;
   E1_min = (- u)/(4.*Eq);

   double* E1_pt = new double [n_E1];
   double* E1_weight = new double [n_E1];
   for(int i=0; i<n_E1; i++)
   {
      E1_pt[i] = E1_pt_standard[i];
      E1_weight[i] = E1_weight_standard[i];
   }
   
   double slope = 1./T;
   slope = 1.0;
   scale_gausspoints(n_E1, 5, 0.0, 0.0, E1_min, slope, E1_pt, E1_weight);
   
   for(int i=0; i<n_E1; i++)
   {
      Integrate_E2(Eq, T, channel, s, t, E1_pt[i], E2_pt_standard, E2_weight_standard, results);
      equilibrium_result += results[0]*E1_weight[i];
      viscous_result1 += results[1]*E1_weight[i];
      viscous_result2 += results[2]*E1_weight[i];
   }

   results[0] = equilibrium_result;
   results[1] = viscous_result1;
   results[2] = viscous_result2;

   delete[] E1_pt;
   delete[] E1_weight;

   return(0);
}


double Integrate_E2(double Eq, double T, int channel, double s, double t, double E1, double** E2_pt_standard, double** E2_weight_standard, double* results)
{
   double eps = 1e-100;
   double equilibrium_result = 0.0;
   double viscous_result1 = 0.0;
   double viscous_result2 = 0.0;
   double E2_min;
   double E2_max;
   double min_1 = (- t)/(4.*Eq);

   double a = - (s + t)*(s + t);
   double b = Eq*((s + t)*(s)) + E1*(-t)*(s + t);
   double c = - (Eq*s + E1*t)*(Eq*s + E1*t) + s*t*(s+t);

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
         results[2] = 0.0e0;
         return (0.0);
      }
   
      double common_factor;
      double mu1 = 0.0e0;
      double mu2 = 0.0e0;
      double mu3 = 0.0e0;

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
         Matrix_elements_sq(channel, s, t, E1, E2_pt[i], Eq, T, Matrix_sq_ptr);
         double Matrix_sq_eq = Matrix_sq_ptr[0];
         double Matrix_sq_noneq = Matrix_sq_ptr[1];

         double f0_E1, f0_E2, f0_E3;
         if(channel == 1) //Compton Scattering
         {
            f0_E1 = Bose_distribution(E1, T, mu1);
            f0_E2 = Fermi_distribution(E2_pt[i], T, mu2);
            f0_E3 = Fermi_distribution(E1 + E2_pt[i] - Eq, T, mu3);
            common_factor = f0_E1*f0_E2*(1 - f0_E3)/(sqrt(a*E2_pt[i]*E2_pt[i] + 2*b*E2_pt[i] + c) + eps);
         }
         else if (channel == 2) // pair annihilation
         {
            f0_E1 = Fermi_distribution(E1, T, mu1);
            f0_E2 = Fermi_distribution(E2_pt[i], T, mu2);
            f0_E3 = Bose_distribution(E1 + E2_pt[i] - Eq, T, mu3);
            common_factor = f0_E1*f0_E2*(1 + f0_E3)/(sqrt(a*E2_pt[i]*E2_pt[i] + 2*b*E2_pt[i] + c) + eps);
         }
         equilibrium_result += common_factor*1.*Matrix_sq_eq*E2_weight[i];
         viscous_result1 += common_factor*(Matrix_sq_eq*viscous_integrand(channel, s, t, E1, E2_pt[i], Eq, T, f0_E1, f0_E2, f0_E3))*E2_weight[i];
         viscous_result2 += common_factor*(Matrix_sq_noneq*1.0)*E2_weight[i];

         delete[] Matrix_sq_ptr;
      }

      delete[] E2_pt;
      delete[] E2_weight;
   }
   else  // no kinematic phase space
   {
      equilibrium_result = 0.0e0;
      viscous_result1 = 0.0e0;
      viscous_result2 = 0.0e0;
   }
   results[0] = equilibrium_result;
   results[1] = viscous_result1;
   results[2] = viscous_result2;

   return(0);
}

inline double viscous_integrand(double channel, double s, double t, double E1, double E2, double Eq, double T, double f0_E1, double f0_E2, double f0_E3)
{
   double eps = 1e-100;
   double E3 = E1 + E2 - Eq;
   double p1 = E1;
   double p2 = E2;
   double p3 = E3;
   double costheta1 = (- s - t + 2*E1*Eq)/(2*p1*Eq + eps);
   double costheta2 = (t + 2*E2*Eq)/(2*p2*Eq + eps);
   double p3_z = p1*costheta1 + p2*costheta2 - Eq; 
   double integrand;
   if(channel == 1)
   {
      integrand = (1. + f0_E1)*deltaf_chi(p1/T)*0.5*(-1. + 3.*costheta1*costheta1) + (1. - f0_E2)*deltaf_chi(p2/T)*0.5*(-1. + 3.*costheta2*costheta2) - f0_E3*deltaf_chi(p3/T)/p3/p3*(-0.5*p3*p3 + 1.5*p3_z*p3_z);
   }
   else if(channel == 2)
   {
      integrand = (1. - f0_E1)*deltaf_chi(p1/T)*0.5*(-1. + 3.*costheta1*costheta1) + (1. - f0_E2)*deltaf_chi(p2/T)*0.5*(-1. + 3.*costheta2*costheta2) + f0_E3*deltaf_chi(p3/T)/p3/p3*(-0.5*p3*p3 + 1.5*p3_z*p3_z);
   }
   else
   {
      cout << "viscous_integrand:: channel can not be found! (channel = " << channel << ")" << endl;
      exit(1);
   }

   return(integrand);
}

inline double Bose_distribution(double E, double T, double mu)
{
    return(1.0/(exp((E-mu)/T) - 1.0));
}

inline double Fermi_distribution(double E, double T, double mu)
{
    return(1.0/(exp((E-mu)/T) + 1.0));
}

inline double deltaf_chi(double p)
{ 
    return(pow(p, deltaf_alpha));
}
