#ifndef MATRIX_ELEMENTS_SQ_H
#define MATRIX_ELEMENTS_SQ_H

#include "parameters.h"
#include "Arsenal.h"
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

void Matrix_elements_sq(int channel, double s, double t, double E1, double E2, double Eq, double T, double* result_ptr);
void Matrix_elements_sq_Compton(double s, double t, double E1, double E2, double Eq, double T, double* result_ptr);
void Matrix_elements_sq_Annihilation(double s, double t, double E1, double E2, double Eq, double T, double* result_ptr);

void get_quark_selfenergy_coefficients(double p_0, double p_i, double T, Selfenergy_coefficients* Sigma_ptr);
double quark_selfenergy_Q(double x);

inline double Repart_ComplexMultiply(double Re1, double Im1, double Re2, double Im2);
inline double Impart_ComplexMultiply(double Re1, double Im1, double Re2, double Im2);
inline double Repart_ComplexDivide(double Re2, double Im1, double Re2, double Im2);
inline double Impart_ComplexDivide(double Re1, double Im1, double Re2, double Im2);
inline double Power(double x, int a);
#endif
