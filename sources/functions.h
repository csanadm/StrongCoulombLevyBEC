#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <complex>
#include <vector>
#include <list>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <sstream>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TNamed.h>
#include <TLegend.h>

#define SQR(x)  ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

//using namespace std;
using std::complex;
double levy_calc(double x, double R, double alpha);
double lorentz_calc(double x, double R);
double gauss_calc(double x, double R);
complex<double> Gamma(complex<double> z);
double Gamma(double x);
complex<double> Digamma(const complex<double> z);
complex<double> Hyp_conf_series(complex<double> a, complex<double> b, complex<double> z);
complex<double> Hyp_conf_asymp(complex<double> a, complex<double> b, complex<double> z);
complex<double> Hyp_conf(complex<double> a, complex<double> b, complex<double> z); // just for ,,normal'' a and b!!!
complex<double> U1_series(complex<double> a, complex<double> z); // b=1
complex<double> U2_series(complex<double> a, complex<double> z); // b=2
complex<double> U_asymp(complex<double> a, complex<double> b, complex<double> z);
complex<double> U2(complex<double> a, complex<double> z);

complex<double> trial_Digamma(const complex<double> z, const int Nterms);
