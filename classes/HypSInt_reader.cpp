#include "HypSInt_reader.h"

HypSInt_reader::HypSInt_reader(const char* filename, bool coulomb)
{
  F2int_1_array = NULL;
  F2int_2_array = NULL;
  F2int_3_real_array = NULL;
  F2int_3_imag_array = NULL;
//  InitTable(filename);
  if(coulomb) InitTableCoulomb(filename);
  else        InitTableStrong(filename);
}

HypSInt_reader::~HypSInt_reader()
{
  if(F2int_1_array) delete [] F2int_1_array;
  if(F2int_2_array) delete [] F2int_2_array;
  if(F2int_3_real_array) delete [] F2int_3_real_array;
  if(F2int_3_imag_array) delete [] F2int_3_imag_array;
}

double HypSInt_reader::read_table1(const double eta, const double kr) const
{
  if(!F2int_1_array)
    return 0;
  int ieta = floor((eta - eta_min) / d_eta); 
  int ikr = floor((kr - kr_min) / d_kr); 
  if(ieta >= 0 && ieta < (Neta + 1))
    if(ikr >= 0 && ikr < (Nkr + 1))
      return F2int_1_array[ieta][ikr];
  
  return 0;
}

double HypSInt_reader::read_table2(const double eta, const double kr) const
{
  if(!F2int_2_array)
    return 0;
  int ieta = floor((eta - eta_min) / d_eta); 
  int ikr = floor((kr - kr_min) / d_kr); 
  if(ieta >= 0 && ieta < (Neta + 1))
    if(ikr >= 0 && ikr < (Nkr + 1))
      return F2int_2_array[ieta][ikr];
  
  return 0;
}

double HypSInt_reader::get_F2int_1(const double eta, const double kr) const
{
  if(!F2int_1_array)
    return 0;
  int ikr = floor((kr - kr_min) / d_kr);
  double krlo = kr_min + ikr * d_kr;
  int ieta = floor((eta - eta_min) / d_eta);
  double etalo = eta_min + ieta * d_eta;
  if(ikr < 0 || ikr > Nkr || ieta < 0 || ieta > Neta)
    return 0;
  if(ikr < Nkr && ieta < Neta)
  { 
    double krlo_lint = 0;
    double krhi_lint = 0;
    krlo_lint = F2int_1_array[ieta][ikr]     * (krlo - kr + d_kr) / d_kr + F2int_1_array[ieta][ikr + 1]     * (kr - krlo) / d_kr;
    krhi_lint = F2int_1_array[ieta + 1][ikr] * (krlo - kr + d_kr) / d_kr + F2int_1_array[ieta + 1][ikr + 1] * (kr - krlo) / d_kr;
    return krlo_lint * (etalo - eta + d_eta) / d_eta + krhi_lint * (eta - etalo) / d_eta;
  }
  else if(ikr == Nkr && ieta < Neta) // Iterpolation only in eta
    return F2int_1_array[ieta][ikr] * (etalo - eta + d_eta) / d_eta + F2int_1_array[ieta + 1][ikr] * (eta - etalo) / d_eta;
  else if(ieta == Neta && ikr < Nkr) // Iterpolation only in kr
    return F2int_1_array[ieta][ikr] * (krlo - kr + d_kr) / d_kr + F2int_1_array[ieta][ikr + 1] * (kr - krlo) / d_kr;
  else if(ikr == Nkr && ieta == Neta)
    return F2int_1_array[ieta][ikr];
  return 0;
}

double HypSInt_reader::get_F2int_2(const double eta, const double kr) const
{
  if(!F2int_2_array)
    return 0;
  int ikr = floor((kr - kr_min) / d_kr);
  double krlo = kr_min + ikr * d_kr;
  int ieta = floor((eta - eta_min) / d_eta);
  double etalo = eta_min + ieta * d_eta;
  if(ikr < 0 || ikr > Nkr || ieta < 0 || ieta > Neta)
    return 0;
  if(ikr < Nkr && ieta < Neta)
  { 
    double krlo_lint = 0;
    double krhi_lint = 0;
    krlo_lint = F2int_2_array[ieta][ikr]     * (krlo - kr + d_kr) / d_kr + F2int_2_array[ieta][ikr + 1]     * (kr - krlo) / d_kr;
    krhi_lint = F2int_2_array[ieta + 1][ikr] * (krlo - kr + d_kr) / d_kr + F2int_2_array[ieta + 1][ikr + 1] * (kr - krlo) / d_kr;
    return krlo_lint * (etalo - eta + d_eta) / d_eta + krhi_lint * (eta - etalo) / d_eta;
  }
  else if(ikr == Nkr && ieta < Neta) // Iterpolation only in eta
    return F2int_2_array[ieta][ikr] * (etalo - eta + d_eta) / d_eta + F2int_2_array[ieta + 1][ikr] * (eta - etalo) / d_eta;
  else if(ieta == Neta && ikr < Nkr) // Iterpolation only in kr
    return F2int_2_array[ieta][ikr] * (krlo - kr + d_kr) / d_kr + F2int_2_array[ieta][ikr + 1] * (kr - krlo) / d_kr;
  else if(ikr == Nkr && ieta == Neta)
    return F2int_2_array[ieta][ikr];
  return 0;
}

double HypSInt_reader::get_F2int_3_real(const double eta, const double kr) const
{
  if(!F2int_3_real_array)
    return 0;
  int ikr = floor((kr - kr_min) / d_kr);
  double krlo = kr_min + ikr * d_kr;
  int ieta = floor((eta - eta_min) / d_eta);
  double etalo = eta_min + ieta * d_eta;
  if(ikr < 0 || ikr > Nkr || ieta < 0 || ieta > Neta)
    return 0;
  if(ikr < Nkr && ieta < Neta)
  {
    double krlo_lint = 0;
    double krhi_lint = 0;
    krlo_lint = F2int_3_real_array[ieta][ikr]     * (krlo - kr + d_kr) / d_kr + F2int_3_real_array[ieta][ikr + 1]     * (kr - krlo) / d_kr;
    krhi_lint = F2int_3_real_array[ieta + 1][ikr] * (krlo - kr + d_kr) / d_kr + F2int_3_real_array[ieta + 1][ikr + 1] * (kr - krlo) / d_kr;
    return krlo_lint * (etalo - eta + d_eta) / d_eta + krhi_lint * (eta - etalo) / d_eta;
  }
  else if(ikr == Nkr && ieta < Neta) // Iterpolation only in eta
    return F2int_3_real_array[ieta][ikr] * (etalo - eta + d_eta) / d_eta + F2int_3_real_array[ieta + 1][ikr] * (eta - etalo) / d_eta;
  else if(ieta == Neta && ikr < Nkr) // Iterpolation only in kr
    return F2int_3_real_array[ieta][ikr] * (krlo - kr + d_kr) / d_kr + F2int_3_real_array[ieta][ikr + 1] * (kr - krlo) / d_kr;
  else if(ikr == Nkr && ieta == Neta)
    return F2int_3_real_array[ieta][ikr];
  return 0;
}

double HypSInt_reader::get_F2int_3_imag(const double eta, const double kr) const
{
  if(!F2int_3_imag_array)
    return 0;
  int ikr = floor((kr - kr_min) / d_kr);
  double krlo = kr_min + ikr * d_kr;
  int ieta = floor((eta - eta_min) / d_eta);
  double etalo = eta_min + ieta * d_eta;
  if(ikr < 0 || ikr > Nkr || ieta < 0 || ieta > Neta)
    return 0;
  if(ikr < Nkr && ieta < Neta)
  {
    double krlo_lint = 0;
    double krhi_lint = 0;
    krlo_lint = F2int_3_imag_array[ieta][ikr]     * (krlo - kr + d_kr) / d_kr + F2int_3_imag_array[ieta][ikr + 1]     * (kr - krlo) / d_kr;
    krhi_lint = F2int_3_imag_array[ieta + 1][ikr] * (krlo - kr + d_kr) / d_kr + F2int_3_imag_array[ieta + 1][ikr + 1] * (kr - krlo) / d_kr;
    return krlo_lint * (etalo - eta + d_eta) / d_eta + krhi_lint * (eta - etalo) / d_eta;
  }
  else if(ikr == Nkr && ieta < Neta) // Iterpolation only in eta
    return F2int_3_imag_array[ieta][ikr] * (etalo - eta + d_eta) / d_eta + F2int_3_imag_array[ieta + 1][ikr] * (eta - etalo) / d_eta;
  else if(ieta == Neta && ikr < Nkr) // Iterpolation only in kr
    return F2int_3_imag_array[ieta][ikr] * (krlo - kr + d_kr) / d_kr + F2int_3_imag_array[ieta][ikr + 1] * (kr - krlo) / d_kr;
  else if(ikr == Nkr && ieta == Neta)
    return F2int_3_imag_array[ieta][ikr];
  return 0;
}

void HypSInt_reader::InitTable(const char* filename)
{ 
  ifstream infile;
  infile.open(filename, ios::in | ios::binary);
  if(infile.is_open()) cout << "File OK" << endl;
  else cout << "Error!" << endl;
  cout << "Initializing HypSInt_reader table in memory using file: " << filename << endl;
  infile.read((char*)&eta_min, sizeof(eta_min));  cout << "  eta_min = " << eta_min << endl;
  infile.read((char*)&eta_max, sizeof(eta_max));  cout << "  eta_max = " << eta_max << endl;
  infile.read((char*)&Neta,    sizeof(Neta));     cout << "  Neta = " << Neta << endl;
  infile.read((char*)&kr_min,  sizeof(kr_min));   cout << "  kr_min = " << kr_min << endl;
  infile.read((char*)&kr_max,  sizeof(kr_max));   cout << "  kr_max = " << kr_max << endl;
  infile.read((char*)&Nkr,     sizeof(Nkr));      cout << "  Nkr = " << Nkr << endl;
  d_eta = (eta_max - eta_min) * 1.0 / Neta;     
  d_kr = (kr_max - kr_min) * 1.0 / Nkr;
  cout << "    d_eta: " << d_eta << " d_kr: " << d_kr << endl;

  float value1 = 0, value2 = 0, value3_real = 0, value3_imag = 0;
  F2int_1_array = new float*[Neta + 1];
  F2int_2_array = new float*[Neta + 1];
  F2int_3_real_array = new float*[Neta + 1];
  F2int_3_imag_array = new float*[Neta + 1];
  for(int ieta = 0; ieta < (Neta + 1); ieta++)
  {
    F2int_1_array[ieta] = new float[Nkr + 1];
    F2int_2_array[ieta] = new float[Nkr + 1];
    F2int_3_real_array[ieta] = new float[Nkr + 1];
    F2int_3_imag_array[ieta] = new float[Nkr + 1];
  }

  int _ikr = 0;

  for(int ikr = 0; ikr < (Nkr + 1); ikr++)
  {
    infile.read((char*)&_ikr, sizeof(_ikr)); // cout << " _ikr: " << _ikr << "\n";

    for(int ieta = 0; ieta < (Neta + 1); ieta++)
    {
      infile.read((char*)&value1, sizeof(value1));
      F2int_1_array[ieta][_ikr] = value1;
//      cout << "reader, eta: " << 0.1*ieta/(Neta+1) << " kr: " << 900.*ikr/(Nkr+1) << " value1: " << value1 << "\t";
//      cout << " value1: " << value1 << "\t";
      infile.read((char*)&value2, sizeof(value2));
      F2int_2_array[ieta][_ikr] = value2;
//      cout << " value2: " << value2 << endl;
      infile.read((char*)&value3_real, sizeof(value3_real));
      F2int_3_real_array[ieta][_ikr] = value3_real;
      infile.read((char*)&value3_imag, sizeof(value3_imag));
      F2int_3_imag_array[ieta][_ikr] = value3_imag;
    }
  }
/*  for(int ikr = 0; ikr < (Nkr + 1); ikr++)
  {
    cout << "good: " << _ikr << "\t" << F2int_1_array[1][ikr] << "\n";
    cout << "bad:  " << _ikr << "\t" << F2int_2_array[1][ikr] << "\n";
  }
*/
  infile.close();
}

void HypSInt_reader::InitTableCoulomb(const char* filename)
{ 
  ifstream infile;
  infile.open(filename, ios::in | ios::binary);
  if(infile.is_open()) cout << "File OK" << endl;
  else cout << "Error!" << endl;
  cout << "Initializing HypSInt_reader table in memory using file: " << filename << endl;
  infile.read((char*)&eta_min, sizeof(eta_min));  cout << "  eta_min = " << eta_min << endl;
  infile.read((char*)&eta_max, sizeof(eta_max));  cout << "  eta_max = " << eta_max << endl;
  infile.read((char*)&Neta,    sizeof(Neta));     cout << "  Neta = " << Neta << endl;
  infile.read((char*)&kr_min,  sizeof(kr_min));   cout << "  kr_min = " << kr_min << endl;
  infile.read((char*)&kr_max,  sizeof(kr_max));   cout << "  kr_max = " << kr_max << endl;
  infile.read((char*)&Nkr,     sizeof(Nkr));      cout << "  Nkr = " << Nkr << endl;
  d_eta = (eta_max - eta_min) * 1.0 / Neta;     
  d_kr = (kr_max - kr_min) * 1.0 / Nkr;
  cout << "    d_eta: " << d_eta << " d_kr: " << d_kr << endl;

  double value1 = 0., value2 = 0., value3_real = 0., value3_imag = 0.;
  F2int_1_array = new float*[Neta + 1];
  F2int_2_array = new float*[Neta + 1];
//  F2int_3_real_array = new float*[Neta + 1];
//  F2int_3_imag_array = new float*[Neta + 1];
  for(int ieta = 0; ieta < (Neta + 1); ieta++)
  {
    F2int_1_array[ieta] = new float[Nkr + 1];
    F2int_2_array[ieta] = new float[Nkr + 1];
//    F2int_3_real_array[ieta] = new float[Nkr + 1];
//    F2int_3_imag_array[ieta] = new float[Nkr + 1];
  }

  int _ikr = 0;
  for(int ikr = 0; ikr < (Nkr + 1); ikr++)
  {
    infile.read((char*)&_ikr, sizeof(_ikr)); // cout << " _ikr: " << _ikr << "\n";
//    if(ikr % 10000 == 0) cout << "ikr: " << ikr << endl;
    for(int ieta = 0; ieta < (Neta + 1); ieta++)
    {
      infile.read((char*)&value1, sizeof(value1));
      F2int_1_array[ieta][_ikr] = (float)value1;
//      cout << "reader, eta: " << 0.1*ieta/(Neta+1) << " kr: " << 900.*ikr/(Nkr+1) << " value1: " << value1 << "\t";
//      cout << " value1: " << value1 << "\t";
      infile.read((char*)&value2, sizeof(value2));
      F2int_2_array[ieta][_ikr] = (float)value2;
//      cout << " value2: " << value2 << endl;
      infile.read((char*)&value3_real, sizeof(value3_real));
//      F2int_3_real_array[ieta][_ikr] = value3_real;
      infile.read((char*)&value3_imag, sizeof(value3_imag));
//      F2int_3_imag_array[ieta][_ikr] = value3_imag;
    }
  }
/*  for(int ikr = 0; ikr < (Nkr + 1); ikr++)
  {
    cout << "good: " << _ikr << "\t" << F2int_1_array[1][ikr] << "\n";
    cout << "bad:  " << _ikr << "\t" << F2int_2_array[1][ikr] << "\n";
  }
*/
  cout << "Initialization done." << endl;
  infile.close();
}

void HypSInt_reader::InitTableStrong(const char* filename)
{ 
  ifstream infile;
  infile.open(filename, ios::in | ios::binary);
  if(infile.is_open()) cout << "File OK" << endl;
  else cout << "Error!" << endl;
  cout << "Initializing HypSInt_reader table in memory using file: " << filename << endl;
  infile.read((char*)&eta_min, sizeof(eta_min));  cout << "  eta_min = " << eta_min << endl;
  infile.read((char*)&eta_max, sizeof(eta_max));  cout << "  eta_max = " << eta_max << endl;
  infile.read((char*)&Neta,    sizeof(Neta));     cout << "  Neta = " << Neta << endl;
  infile.read((char*)&kr_min,  sizeof(kr_min));   cout << "  kr_min = " << kr_min << endl;
  infile.read((char*)&kr_max,  sizeof(kr_max));   cout << "  kr_max = " << kr_max << endl;
  infile.read((char*)&Nkr,     sizeof(Nkr));      cout << "  Nkr = " << Nkr << endl;
  d_eta = (eta_max - eta_min) * 1.0 / Neta;     
  d_kr = (kr_max - kr_min) * 1.0 / Nkr;
  cout << "    d_eta: " << d_eta << " d_kr: " << d_kr << endl;

  double value1 = 0, value2 = 0, value3_real = 0, value3_imag = 0;
//  F2int_1_array = new float*[Neta + 1];
//  F2int_2_array = new float*[Neta + 1];
  F2int_3_real_array = new float*[Neta + 1];
  F2int_3_imag_array = new float*[Neta + 1];
  for(int ieta = 0; ieta < (Neta + 1); ieta++)
  {
//    F2int_1_array[ieta] = new float[Nkr + 1];
//    F2int_2_array[ieta] = new float[Nkr + 1];
    F2int_3_real_array[ieta] = new float[Nkr + 1];
    F2int_3_imag_array[ieta] = new float[Nkr + 1];
  }

  int _ikr = 0;

  for(int ikr = 0; ikr < (Nkr + 1); ikr++)
  {
    infile.read((char*)&_ikr, sizeof(_ikr)); // cout << " _ikr: " << _ikr << "\n";
    for(int ieta = 0; ieta < (Neta + 1); ieta++)
    {
      infile.read((char*)&value1, sizeof(value1));
//      F2int_1_array[ieta][_ikr] = value1;
//      cout << "reader, eta: " << 0.1*ieta/(Neta+1) << " kr: " << 900.*ikr/(Nkr+1) << " value1: " << value1 << "\t";
//      cout << " value1: " << value1 << "\t";
      infile.read((char*)&value2, sizeof(value2));
//      F2int_2_array[ieta][_ikr] = value2;
//      cout << " value2: " << value2 << endl;
      infile.read((char*)&value3_real, sizeof(value3_real));
      F2int_3_real_array[ieta][_ikr] = (float)value3_real;
      infile.read((char*)&value3_imag, sizeof(value3_imag));
      F2int_3_imag_array[ieta][_ikr] = (float)value3_imag;
    }
  }
/*  for(int ikr = 0; ikr < (Nkr + 1); ikr++)
  {
    cout << "good: " << _ikr << "\t" << F2int_1_array[1][ikr] << "\n";
    cout << "bad:  " << _ikr << "\t" << F2int_2_array[1][ikr] << "\n";
  }
*/
  infile.close();
  cout << "Initialization done." << endl;
}
