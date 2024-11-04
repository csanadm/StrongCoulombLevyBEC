#include <functions.h>
using namespace std;
//#define K_CONST (139.57018 / 2. / 137.036)
#define K_CONST (493.677 / 2. / 137.036)
#define HBARC 197.3269788
#define MASS_PI 139.57018
#define MASS2_PI 19479.8351452
#define MASS_KA 493.677
#define MASS2_KA 243716.98
const char* infile_name = "tables/CorrFunc_raw_table_kaons.dat";
const char* outfile_name = "tables/CorrFunc_table_kaons.dat";

const complex<double> I(0., 1.);
const double a00 = 0.220;
const double a02 = -0.0444;
const double b02 = -0.0857;
const double c02 = -0.00221;
const double d02 = -0.000129;
const double s02 = -21.62*MASS2_PI;
const double f0 = -0.063;//(2.*a00 + a02)/3.;// = 0.132 in femtometers (Lednicky value = 0.2)
const double d0 = -10.; // in femtometers
//const double B0 = -80.4;
//const double B1 = -73.6;
const double B0 = -79.4;
const double B1 = -63.0;
const double z2 = 143.5;//MASS_PI;
double K_calc_kaons(const double k);
double K_calc_CGL(const double k);
double K_calc_GM(const double k);
complex<double> chi(const double eta);

int main(int argc, char** argv)
{
  if(argc != 1)
  {
    cout << "Usage: " << argv[0] << " and nothing else." << endl;
    return 1;
  }

  ifstream infile;
  ofstream outfile;
  infile.open(infile_name, ios::in | ios::binary);
  outfile.open(outfile_name, ios::out | ios::binary); 

  cout << "Reading aux info..." << endl;
  double alpha_min = 0.; infile.read((char*)&alpha_min, sizeof(alpha_min));  cout << "alpha_min = " << alpha_min << endl;
  double alpha_max = 0.; infile.read((char*)&alpha_max, sizeof(alpha_max));  cout << "alpha_max = " << alpha_max << endl;
  int    Nalpha    = 0;  infile.read((char*)&Nalpha,    sizeof(Nalpha));     cout << "Nalpha = "    << Nalpha    << endl;
  double k_min     = 0.; infile.read((char*)&k_min,     sizeof(k_min));      cout << "k_min = "     << k_min     << endl;
  double k_max     = 0.; infile.read((char*)&k_max,     sizeof(k_max));      cout << "k_max = "     << k_max     << endl;
  int    Nk        = 0;  infile.read((char*)&Nk,        sizeof(Nk));         cout << "Nk = "        << Nk        << endl;
  double Rcc_min   = 0.; infile.read((char*)&Rcc_min,   sizeof(Rcc_min));    cout << "Rcc_min = "   << Rcc_min   << endl;
  double Rcc_max   = 0.; infile.read((char*)&Rcc_max,   sizeof(Rcc_max));    cout << "Rcc_max = "   << Rcc_max   << endl;
  int    NRcc      = 0;  infile.read((char*)&NRcc,      sizeof(NRcc));       cout << "NRcc = "      << NRcc      << endl;
  cout << "Read aux info." << endl;

  outfile.write((char*)&alpha_min, sizeof(alpha_min));
  outfile.write((char*)&alpha_max, sizeof(alpha_max));
  outfile.write((char*)&Nalpha,    sizeof(Nalpha));
  outfile.write((char*)&k_min,     sizeof(k_min));
  outfile.write((char*)&k_max,     sizeof(k_max));
  outfile.write((char*)&Nk,        sizeof(Nk));
  outfile.write((char*)&Rcc_min,   sizeof(Rcc_min));
  outfile.write((char*)&Rcc_max,   sizeof(Rcc_max));
  outfile.write((char*)&NRcc,      sizeof(NRcc));
  
  double d_k = (k_max - k_min) * 1.0 / Nk;
  double d_Rcc = (Rcc_max - Rcc_min) * 1.0 / NRcc;

  for(int ialpha = 0; ialpha < (Nalpha + 1); ialpha++)  
  {
    double k = k_min;
    double** result1 = new double*[Nk+1]; for(int _ik = 0; _ik < Nk+1; _ik++) result1[_ik] = new double[NRcc+1];
    double** result2 = new double*[Nk+1]; for(int _ik = 0; _ik < Nk+1; _ik++) result2[_ik] = new double[NRcc+1];
    double** result3 = new double*[Nk+1]; for(int _ik = 0; _ik < Nk+1; _ik++) result3[_ik] = new double[NRcc+1];
    double** result4 = new double*[Nk+1]; for(int _ik = 0; _ik < Nk+1; _ik++) result4[_ik] = new double[NRcc+1];
    cout << "result arrays created." << endl;
    for(int ik = 0; ik < (Nk + 1); ik++)  
    {
      double eta = K_CONST / k;
      double Gamow_factor = (2. * M_PI * eta) / (exp(2. * M_PI * eta) - 1.); // == |\c N|^2
      double delta_coul = arg(Gamma(1. + I * eta));
      complex<double> f_c = 1. / (1. / K_calc_kaons(k) - 2. * k / HBARC * eta * chi(eta)); 
      // calculated in units of fm, notation and formula by Lednicky, Physics of Particles and Nuclei, 2009, Vol. 40, No. 3, pp. 307-352.
      // sin(Delta)*exp(I*Delta) = k/HBARC*Gamow_factor*f_c .
      double result3_factor = -2. * M_PI * 16. * exp(M_PI * eta) * SQR(k / HBARC * Gamow_factor) * norm(f_c);
      double Rcc = Rcc_min;
      for(int iRcc = 0; iRcc < (NRcc + 1); iRcc++)
      {
        infile.read((char*)&result1[ik][iRcc], sizeof(result1[ik][iRcc]));
        infile.read((char*)&result2[ik][iRcc], sizeof(result2[ik][iRcc]));
        infile.read((char*)&result3[ik][iRcc], sizeof(result3[ik][iRcc]));
        result1[ik][iRcc] = result1[ik][iRcc] * 2. * M_PI * Gamow_factor;
        result2[ik][iRcc] = result2[ik][iRcc] * 2. * M_PI * Gamow_factor;
        result3[ik][iRcc] = result3[ik][iRcc] * result3_factor;
        Rcc += d_Rcc;
      }
      k += d_k;
    }
    k = k_min;
    for(int ik = 0; ik < (Nk + 1); ik++)  
    {
      double eta = K_CONST / k;
      double Gamow_factor = (2. * M_PI * eta) / (exp(2. * M_PI * eta) - 1.); // == |\c N|^2
      double delta_coul = arg(Gamma(1. + I * eta));
      complex<double> f_c = 1. / (1. / K_calc_kaons(k) - 2. * k / HBARC * eta * chi(eta)); 
      // calculated in units of fm, notation and formula by Lednicky, Physics of Particles and Nuclei, 2009, Vol. 40, No. 3, pp. 307-352.
      // sin(Delta)*exp(I*Delta) = k/HBARC*Gamow_factor*f_c .
      complex<double> result4_factor = 2. * M_PI * 8. * I * Gamma(1. + I * eta) * exp(-2. * I * delta_coul) * k / HBARC * Gamow_factor * conj(f_c);
      double Rcc = Rcc_min;
      for(int iRcc = 0; iRcc < (NRcc + 1); iRcc++)
      {
        double result4_real; infile.read((char*)&result4_real, sizeof(result4_real));
        double result4_imag; infile.read((char*)&result4_imag, sizeof(result4_imag));
        complex<double> result4complex(result4_real, result4_imag);
        result4complex = result4complex * result4_factor;
        result4[ik][iRcc] = real(result4complex);
        Rcc += d_Rcc;
      }
      k += d_k;
    }
    for(int ik = 0; ik < (Nk + 1); ik++)
      for(int iRcc = 0; iRcc < (NRcc + 1); iRcc++)
      {
        float CorrFunc_wout = (float)(result1[ik][iRcc] + result2[ik][iRcc]);
        float CorrFunc_with = (float)(result1[ik][iRcc] + result2[ik][iRcc] + result3[ik][iRcc] + result4[ik][iRcc]);
        outfile.write((char*)&CorrFunc_wout, sizeof(CorrFunc_wout));
        outfile.write((char*)&CorrFunc_with, sizeof(CorrFunc_with));
      }

    cout << "writing CorrFunc w/o ialpha: " << ialpha << endl;
    for(int _ik = 0; _ik < Nk+1; _ik++) delete result1[_ik]; delete result1;
    for(int _ik = 0; _ik < Nk+1; _ik++) delete result2[_ik]; delete result2;
    for(int _ik = 0; _ik < Nk+1; _ik++) delete result3[_ik]; delete result3;
    for(int _ik = 0; _ik < Nk+1; _ik++) delete result4[_ik]; delete result4;
  }
  
  infile.close();
  outfile.close();
  return 0;
}

double K_calc_kaons(const double k)
{
  double f0_K = -0.0646;
  double d0_K = -10.;
  return f0_K - 0.5 * d0_K * SQR(f0_K * k / HBARC);
}

double K_calc_CGL(const double k)
{
//  return f0 - 0.5 * d0 * SQR(f0 * k / HBARC);
  double s = 4.*(MASS2_PI + k*k);
  double x = k/MASS_PI;
  return HBARC * (2./sqrt(s)) * ((4.*MASS2_PI - s02)/(s-s02)) * (a02 + b02*x*x + c02*x*x*x*x + d02*x*x*x*x*x*x);
}

double K_calc_GM(const double k)
{
  double s = 4.*(MASS2_PI + k*k);
  double s0 = 1050.*1050.;
  double w = ((sqrt(s) - sqrt(s0-s))/(sqrt(s) + sqrt(s0-s)));
  return HBARC * (2./sqrt(s)) * (s/MASS2_PI - 2.*z2*z2/MASS2_PI) * (1./(B0 + B1 * w));
}

complex<double> chi(const double eta)
{
//  complex<double> result(real(Digamma(I*eta)-log(eta)), 0.);
  complex<double> result(real(Digamma(I*eta)) - log(eta), 0.);
  return result + I * 0.5 / eta * exp(- M_PI * eta / 2.) * Gamma(1. + I*eta);
}
