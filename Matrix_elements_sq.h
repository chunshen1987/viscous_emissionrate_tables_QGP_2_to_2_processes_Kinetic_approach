#ifndef MATRIX_ELEMENTS_SQ_H
#define MATRIX_ELEMENTS_SQ_H

#include "parameters.h"
#include "Arsenal.h"
#include "QGP_2to2_Scattering_Kinetic.h"


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
