#include "functions.h"
#define NMIN_SIN 60 // at least this much points in a quarter period of the sine function
#define NMIN_EXP 50 // at least this much points in one 1/e change of the exponential
#define LIMIT 30    // if the exponent in the exponential in the iterand is less than minus LIMIT, stop iteration
#define EPSILON 1e-20 // this is for the terms in the series ofdthe hypergeometric function.
#define EPSILON_ASYM 1e-6 // this is for the terms in the asymptotic series ofdthe hypergeometric function.

using namespace std;
using std::complex;

// for the Lanczos approximation of the gamma function
const double lanczos_coeff[9] = {
0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
const double Euler_gamma = 0.57721566490153286060;
const double g_lanczos = 7.;
const complex<double> I(0., 1.);

double levy_calc(double x, double R, double alpha)
{
  double result = 0.;  
  double t = 0.;
  if(x < 0. || alpha < 0.5) return 0;
  if(x < 0.001)
    return  pow(2., 3 / alpha - 1.) / alpha / SQR(M_PI) / CUBE(R) * Gamma(3. / alpha);
  
  double dt_sin = M_PI / 2. / NMIN_SIN;
  double dt_exp = pow(2., 1.0 / alpha) * x / NMIN_EXP;
  double dt = (dt_sin < dt_exp) ? dt_sin : dt_exp;
  double exponent_dt = 0., exponent_mean = 0.;
  double exponential = 1., exponential_dt = 1., exponential_mean = 1.;
  double sinust = 0., sinust_dt = 0., sinust_mean = 0.;
  for(;;)
  {
    exponent_dt = -pow(t + dt, alpha ) / ( 2 * pow(x, alpha) );
    exponent_mean =  -pow((2 * t + dt) / 2, alpha ) / ( 2 * pow(x, alpha) );
    exponential = exponential_dt;
    exponential_mean = exp(exponent_mean);
    exponential_dt = exp(exponent_dt);
    sinust = sinust_dt;
    sinust_mean = sin((2. * t + dt) / 2.);
    sinust_dt = sin(t + dt);
    if(exponent_dt + log(1. + t) < -LIMIT) // change was here
      break;
    double parabolic_sum = (dt / 6.) * ( (exponential * sinust * t) + 4. * (exponential_mean * sinust_mean * ((2. * t + dt) / 2.)) + (exponential_dt * sinust_dt * (t + dt)));
    result += parabolic_sum;
    t += dt; 
  }
  return result * 1.0/(2. * SQR(M_PI) * CUBE(R * x));
}

double lorentz_calc(double x, double R)
{
  return 8. / (SQR(M_PI) * CUBE(R) * SQR(4. * SQR(x) + 1.));
}

double gauss_calc(double x, double R)
{
  return 1. / sqrt( CUBE(2. * M_PI)) / CUBE(R) * exp (-SQR(x) / 2.);
}

complex<double> Gamma(complex<double> z)  //Utilizing the "Lanczos approximation"
{
  complex<double> A_g = lanczos_coeff[0]; 
  for(int i=1; i <= (int)(g_lanczos + 1.); i++)
    A_g += lanczos_coeff[i] / (z - 1. + double(i));
  return sqrt(2. * M_PI) * pow( z + g_lanczos - 0.5, z - 0.5 ) * exp(- (z + g_lanczos - 0.5)) * A_g;
}

double Gamma(double x)
{
  complex<double> z(x, 0.);
  return real(Gamma(z));
}

complex<double> Digamma(const complex<double> z) // We use Lanczos approximation; series in trial_Digamma(), ditched.
{
  complex<double> A_g_denom = lanczos_coeff[0];
  for(int i=1; i <= (int)(g_lanczos + 1.); i++)
    A_g_denom += lanczos_coeff[i] / (z - 1. + double(i));
  complex<double> A_g_numer(0., 0.);
  for(int i=1; i <= (int)(g_lanczos + 1.); i++)
    A_g_numer += lanczos_coeff[i] / SQR(z - 1. + double(i));
  return log(z + g_lanczos - 0.5) - g_lanczos / (z + g_lanczos - 0.5) - A_g_numer / A_g_denom;
}

complex<double> U1_series(complex<double> a, complex<double> z) // b=1
{
  complex<double> beta_n(-2.*Euler_gamma, 0.);
  beta_n = beta_n - Digamma(a) - log(z); // calculating the first term's coefficient.
  complex<double> alpha_n(1. , 0.);
  complex<double> term = beta_n;
  complex<double> result(0., 0.);
  double n = 0;
  while(abs(term) > EPSILON)
  {
    result += term;
    term = term * (a + n) / (n + 1.) / (n + 1.) * z / beta_n;
    beta_n = beta_n + 2./(n + 1.) - 1./(n + a); 
    term = term * beta_n;
    n = n + 1.;
  }
  return result / Gamma(a);
}

complex<double> U2_series(complex<double> a, complex<double> z) // b=2
{
  complex<double> beta_n(-2.*Euler_gamma + 1. , 0.);
  beta_n = beta_n - Digamma(a) - log(z); // calculating the first term's coefficient.
  complex<double> alpha_n(1. , 0.);
  complex<double> term = beta_n;
  complex<double> result(0., 0.);
  double n = 0;
  while(abs(term) > EPSILON)
  {
    result += term;
    term = term * (a + n) / (n + 1.) / (n + 2.) * z / beta_n;
    beta_n = beta_n + 1./(n+1.) + 1./(n+2.) - 1./(n+a); 
    term = term * beta_n;
    n = n + 1.;
  }
  return -1. * (result - 1. / z / (a - 1.)) / Gamma(a - 1.);
}

complex<double> Hyp_conf_series(complex<double> a, complex<double> b, complex<double> z)
{
  complex<double> term(1. , 0.);
  complex<double> result(0., 0.);
  double n = 0;
  while(abs(term) > EPSILON)
  {
    result += term; 
    term = term * (a + n) / (b + n) * z / (n + 1.);
    n = n + 1.;
  }
  return result;
//  complex<double> term(1. , 0.);
//  complex<float> result(0., 0.);
//  double n = 0;
//  while(abs(term) > EPSILON)
//  {
//    result += (complex<double>)term; 
//    term = term * (complex<double>)( (a + n) / (b + n) * z / (n + 1.) );
//    n = n + 1.;
//  }
//  return result;
}

complex<double> U_asymp(complex<double> a, complex<double> b, complex<double> z) // this is NOT the function denoted by "G" by Landau, F needs to be calculated in a bit more complicated way.
{ 
  complex<double> result(0, 0);
  if(abs(z) < 20)
    return result;

  complex<double> term(1. , 0.);
  double n = 0;
  while((abs(term) > EPSILON_ASYM))  // BEWARE!!! Asymptotic series, be careful.
//  while(n < 5.)  // BEWARE!!! Asymptotic series, be careful.
  {
    result += term; 
    term = -1. * term * ((a + n) * (a - b + 1. + n)) / (z * (n + 1.));
    n = n + 1.;
  }
  return result * pow(z, -a);
}

complex<double> Hyp_conf_asymp(complex<double> a, complex<double> b, complex<double> z) 
{
  double Pi_phase = (arg(z) > 0.) ? M_PI : -M_PI;
  return Gamma(b) * exp(I * Pi_phase * a) / Gamma(b - a) * U_asymp(a, b, -z) + exp(I * Pi_phase * (a - b)) * Gamma(b) / Gamma(a) * exp(z) * U_asymp(b - a, b, -z);
}

complex<double> U2(complex<double> a, complex<double> z)
{
  return (abs(z) < 30.) ? U2_series(a, z) : U_asymp(a, 2., z);
}

complex<double> Hyp_conf(complex<double> a, complex<double> b, complex<double> z)
{
  return (abs(z) < 30.) ? Hyp_conf_series(a, b, z) : Hyp_conf_asymp(a, b, z); // beware!!! just for "normal" a and b!!!
}

complex<double> trial_Digamma(const complex<double> z, const int Nterms)
{
  complex<double> result(0., 0.);
  if(abs(z - 1.) < 1e-20) return result;
  if(Nterms < 0)
  {
    complex<double> A_g_denom = lanczos_coeff[0];
    for(int i=1; i <= (int)(g_lanczos + 1.); i++)
      A_g_denom += lanczos_coeff[i] / (z - 1. + double(i));
    complex<double> A_g_numer(0., 0.);
    for(int i=1; i <= (int)(g_lanczos + 1.); i++)
      A_g_numer += lanczos_coeff[i] / SQR(z - 1. + double(i));

    return log(z + g_lanczos - 0.5) - g_lanczos / (z + g_lanczos - 0.5) - A_g_numer / A_g_denom;
  }
  double n = 1.;
  while((int)n < Nterms)
  {
    result = result + 1. / n / (n + z - 1.);
    n = n + 1.;
  }
  return result * (z - 1.) - Euler_gamma;
}

