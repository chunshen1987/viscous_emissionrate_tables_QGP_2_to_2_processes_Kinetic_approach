#ifndef QGP_2to2_Scattering_Kinetic_H
#define QGP_2to2_Scattering_Kinetic_H

#include <string>
#include "ParameterReader.h"
#include "Physicalconstants.h"

typedef struct
{
   double Re_A0;
   double Im_A0;
   double Re_B0;
   double Im_B0;
   double Re_A1;
   double Im_A1;
   double Re_B1;
   double Im_B1;
   double Re_C1;
   double Im_C1;
}Selfenergy_coefficients;

class QGP_2to2_Scattering_Kinetic;

struct CCallbackHolder
{
   QGP_2to2_Scattering_Kinetic* clsPtr;
   void *params;
};

inline double Power(double x, int a);

class QGP_2to2_Scattering_Kinetic
{
   private:
      ParameterReader *paraRdr;

      Physicalconstants Phycons;

      int n_Eq, n_Temp;
      double *Eq_tb, *T_tb;
      double **equilibrium_results, **viscous_results;

      int channel;
      string filename;
      double p_cutoff;

      int n_qtilde;
      double *qtilde_pt;
      double *equilibriumTilde_results, *viscousTilde_results;

      // Gaussian quadrature points for phase space integrations 
      int n_s;
      double *s_pt, *s_weight, *s_pt_standard, *s_weight_standard;
      int n_t;
      double **t_pt, **t_weight, **t_pt_standard, **t_weight_standard;
      double **Matrix_elements_sq_ptr;
      int n_E1;
      double *E1_pt_standard, *E1_weight_standard;
      int n_E2;
      double **E2_pt_standard, **E2_weight_standard;

      double deltaf_alpha;

   public:
      QGP_2to2_Scattering_Kinetic(ParameterReader* paraRdr_in);
      ~QGP_2to2_Scattering_Kinetic();

      void buildupEmissionrate2DTable();
      void output_emissionrateTable();

      void set_gausspoints();
      void calculateEmissionrates(int channel, string filename);
      void scale_gausspoints_st(double qtilde);
      void Integrate_E1(double qtilde, double s_tilde, double t_tilde, double* results);
      void Integrate_E2(double qtilde, double s_tilde, double t_tilde, double E1tilde, double* results);
      double testFunc(double x, void* params);
      double viscous_integrand(double s_tilde, double t_tilde, double E1tilde, double E2tilde, double qtilde, double f0_E1, double f0_E2, double f0_E3);

      void get_quark_selfenergy_coefficients(double p_0_tilde, double p_i_tilde, Selfenergy_coefficients* Sigma_ptr);
      double quark_selfenergy_Q(double x);

      void Matrix_elements_sq(double s, double t, double E1, double E2, double qtilde, double* result_ptr);
      void Matrix_elements_sq_Compton(double s, double t, double E1, double E2, double qtilde, double* result_ptr);
      void Matrix_elements_sq_Annihilation(double s, double t, double E1, double E2, double qtilde, double* result_ptr);
      
      double Bose_distribution(double Etilde);
      double Fermi_distribution(double Etilde);
      double deltaf_chi(double ptilde);

      double Repart_ComplexMultiply(double Re1, double Im1, double Re2, double Im2);
      double Impart_ComplexMultiply(double Re1, double Im1, double Re2, double Im2);
      double Repart_ComplexDivide(double Re1, double Im1, double Re2, double Im2);
      double Impart_ComplexDivide(double Re1, double Im1, double Re2, double Im2);
      
      double eqRateintegrands(double E1tilde, void *params);
      static double CCallback_eqRateintegrands(double x, void* params)
      {
         CCallbackHolder* h = static_cast<CCallbackHolder*>(params);
         return h->clsPtr->eqRateintegrands(x, h->params);
      }
      double eqRateintegrandt(double E1tilde, void *params);
      static double CCallback_eqRateintegrandt(double x, void* params)
      {
         CCallbackHolder* h = static_cast<CCallbackHolder*>(params);
         return h->clsPtr->eqRateintegrandt(x, h->params);
      }
      double eqRateintegrandE1(double E1tilde, void *params);
      static double CCallback_eqRateintegrandE1(double x, void* params)
      {
         CCallbackHolder* h = static_cast<CCallbackHolder*>(params);
         return h->clsPtr->eqRateintegrandE1(x, h->params);
      }
      double eqRateintegrandE2(double E2, void *params);
      static double CCallback_eqRateintegrandE2(double x, void* params)
      {
         CCallbackHolder* h = static_cast<CCallbackHolder*>(params);
         return h->clsPtr->eqRateintegrandE2(x, h->params);
      }
      static double CCallback_test(double x, void* params)
      {
         CCallbackHolder* h = static_cast<CCallbackHolder*>(params);
         return h->clsPtr->testFunc(x, h->params);
      }
};


#endif
