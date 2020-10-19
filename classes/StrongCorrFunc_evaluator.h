#ifndef StrongCorrFunc_evaluator_H
#define StrongCorrFunc_evaluator_H

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <cmath>

class StrongCorrFunc_evaluator
{
 public:
  StrongCorrFunc_evaluator(const char* filename);
  virtual ~StrongCorrFunc_evaluator();
  float get_value_CC_wout(const double alpha, const double k, const double Rcc) const;
  float get_value_CC_with(const double alpha, const double k, const double Rcc) const;
  float read_table_wout(const double alpha, const double k, const double Rcc) const;
  float read_table_with(const double alpha, const double k, const double Rcc) const;
  float interpolate_in_table_wout(const double alpha, const double k, const double Rcc) const;
  float interpolate_in_table_with(const double alpha, const double k, const double Rcc) const;
  float getValue_wout(const double alpha, const double k, const double Rcc) const; // the user should use this!!!
  float getValue_with(const double alpha, const double k, const double Rcc) const; // the user should use this!!!

  float Levy_func(const double alpha, const double k, const double Rcc) const;
 
 private:
  float*** CC_array_wout;
  float*** CC_array_with;
  double alpha_min;
  double alpha_max;
  double k_min;
  double k_max;
  double Rcc_min;
  double Rcc_max;
  int Nalpha;
  int Nk;
  int NRcc;
  double d_alpha;
  double d_k;
  double d_Rcc;

  void Init(const char* filename);
};

#endif // StrongCorrFunc_evaluator_H
