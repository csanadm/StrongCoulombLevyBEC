#ifndef HypSInt_reader_H
#define HypSInt_reader_H

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <cmath>
using namespace std;

class HypSInt_reader
{
 public:
  HypSInt_reader(const char* filename, bool coulomb);
  virtual ~HypSInt_reader();
  double read_table1(const double eta, const double kr) const;
  double read_table2(const double eta, const double kr) const;
  double get_F2int_1(const double eta, const double kr) const; // these are the proper one, with linear interpolation
  double get_F2int_2(const double eta, const double kr) const; // these are the proper one, with linear interpolation
  double get_F2int_3_real(const double eta, const double kr) const;
  double get_F2int_3_imag(const double eta, const double kr) const;
 
 private:
  float** F2int_1_array; // int_-1^1 dy abs(F^2)) 
  float** F2int_2_array; // int_-1^1 dy e^2ikr F_1 F*_2 
  float** F2int_3_real_array; //  int_-1^1 dy F
  float** F2int_3_imag_array;
  
  double eta_min;
  double eta_max;
  double kr_min;
  double kr_max;
  int Neta;
  int Nkr;
  double d_eta;
  double d_kr;

  void InitTable(const char* filename);
  void InitTableCoulomb(const char* filename);
  void InitTableStrong(const char* filename);
};

#endif // HypSInt_reader_H
