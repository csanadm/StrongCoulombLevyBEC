#ifndef CorrFunc_reader_H
#define CorrFunc_reader_H

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <cmath>
using namespace std;

class CorrFunc_reader
{
 public:
  CorrFunc_reader(const char* filename);
  virtual ~CorrFunc_reader();
  double get_value_CC1(const double alpha, const double k, const double Rcc) const;
  double get_value_CC2(const double alpha, const double k, const double Rcc) const;
  double read_table1(const double alpha, const double k, const double Rcc) const;
  double read_table2(const double alpha, const double k, const double Rcc) const;
  double getValue(const double alpha, const double k, const double Rcc) const; // the user should use this!!!
 
 private:
  double*** CC1_array;
  double*** CC2_array;
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

  void InitTable(const char* filename);
};

#endif // CorrFunc_reader_H
