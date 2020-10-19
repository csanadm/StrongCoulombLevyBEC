#include <functions.h>
//#define NJOBS 90
#define NJOBS 18000

using namespace std;

complex<double> I(0, 1.);

int main(int argc, char** argv)
{
  if(argc != 2)
  {
    cout << "Usage: " << argv[0] << " ijob (out of " << NJOBS << ")." << endl;
    return 1;
  }

//  double eta_min = atof(argv[1]);
//  double eta_max = atof(argv[2]);
//  int Neta = atoi(argv[3]);
//  double kr_min = atof(argv[4]);
//  double kr_max = atof(argv[5]);
//  int Nkr = atoi(argv[6]);

  double eta_min = 0.;
  double eta_max = 1.00;// Qmin[MeV] = 1.02/etamax
  int    Neta    = 1000;
  double kr_min  = 0.;
//  double kr_max  = 900;  // kr = 15*rt, if hbarc=200, k=100, Rcc=30
  double kr_max  = 6000;  // trying out a more precise calc.
//  int    Nkr     = 18000;
  int    Nkr     = 180000;// trying out a more precise calc.
//krmax*2, Nkr*2*1.5
  int ijob  = atoi(argv[1]);
  if(ijob < 0 || ijob > NJOBS)
  {
    cout << "Wrong job number, exiting." << endl; 
    return 1;
  }

  cout << "Filling HypSInt_table: " << endl;
  cout << eta_min << "<=eta<" << eta_max << ", " << Neta << "+1 steps; for each eta: " << kr_min << "<=kr<" << kr_max << ", " << Nkr << "+1 steps." << endl;

  ofstream outfile; // the integral of abs(F)^2 wrt y is the first value, the integral of e^2ikr * F1  F2 is the second.
  stringstream outfile_name(""); outfile_name << "HypSInt_table_" << ijob << ".dat";
  outfile.open(outfile_name.str().c_str(), ios::out | ios::binary); 
  if(ijob == 0)
  {
    outfile.write((char*)&eta_min , sizeof(eta_min));
    outfile.write((char*)&eta_max , sizeof(eta_max)); 
    outfile.write((char*)&Neta   , sizeof(Neta)); 
    outfile.write((char*)&kr_min , sizeof(kr_min));
    outfile.write((char*)&kr_max , sizeof(kr_max)); 
    outfile.write((char*)&Nkr   , sizeof(Nkr)); 
    cout << eta_min << endl;
    cout << eta_max << endl;
    cout << Neta << endl;
    cout << kr_min << endl;
    cout << kr_max << endl;
    cout << Nkr << endl;
  }
  
  double d_kr = (kr_max - kr_min) * 1.0 / Nkr; 
  for(int ikr = 0; ikr < (Nkr + 1); ikr++)  
  {
    if(ikr % NJOBS != ijob) continue;
    double eta = eta_min;
    double d_eta = (eta_max - eta_min) * 1.0 / Neta; 
    double kr = kr_min + (kr_max - kr_min) * ikr / Nkr;
    outfile.write((char*)&ikr, sizeof(ikr));
//    int Ny = floor(20 * kr); // to catch the oscillations of both F and e^2ikr : approx. 16 points in e^2ikry's quarter period. 
    int Ny = floor(25 * kr); //trying out a more precise calculation 
    if(Ny < 100) Ny = 100;
    for(int ieta = 0; ieta < (Neta + 1); ieta++)
    {
      double hyp_integral1 = 0.;
      double hyp_integral2 = 0.;
      complex<double> hyp_integral3 = 0.;
      double y = -1.;
      double d_y = 2. / Ny;
      for(int iy = 0; iy < Ny; iy++)
      {
        double term1_lo  = norm(Hyp_conf(-I*eta, 1., I*kr*(1. - y)));
        double term1_mid = norm(Hyp_conf(-I*eta, 1., I*kr*(1. - (y + .5*d_y))));
        double term1_hi  = norm(Hyp_conf(-I*eta, 1., I*kr*(1. - (y +    d_y))));
        
        double term2_lo  = real(Hyp_conf(1.+I*eta, 1., -I*kr*(1. - y)) * Hyp_conf(1.-I*eta, 1., I*kr*(1.+y)));
        double term2_mid = real(Hyp_conf(1.+I*eta, 1., -I*kr*(1. - (y+.5*d_y))) * Hyp_conf(1.-I*eta, 1., I*kr*(1.+(y+.5*d_y))));
        double term2_hi  = real(Hyp_conf(1.+I*eta, 1., -I*kr*(1. - (y+   d_y))) * Hyp_conf(1.-I*eta, 1., I*kr*(1.+(y+   d_y))));

        complex<double> term3_lo  = Hyp_conf(1.+I*eta, 1., -I*kr*(1. + y));
        complex<double> term3_mid = Hyp_conf(1.+I*eta, 1., -I*kr*(1. + (y + .5*d_y)));
        complex<double> term3_hi  = Hyp_conf(1.+I*eta, 1., -I*kr*(1. + (y +    d_y)));

        hyp_integral1 += d_y / 6. * (term1_lo + 4. * term1_mid + term1_hi);
        hyp_integral2 += d_y / 6. * (term2_lo + 4. * term2_mid + term2_hi);
        hyp_integral3 += d_y / 6. * (term3_lo + 4. * term3_mid + term3_hi);
        y += d_y;
      }
//      cout << "kr: " << 900.*ikr/(Nkr+1) << "\t eta: " << 0.1*ieta/(Neta+1) << "\t int1: " << hyp_integral1 << "\t int2: " << hyp_integral2 << endl;
      outfile.write((char*)&hyp_integral1, sizeof(hyp_integral1));
      outfile.write((char*)&hyp_integral2, sizeof(hyp_integral2));
      outfile.write((char*)&real(hyp_integral3), sizeof(real(hyp_integral3)));
      outfile.write((char*)&imag(hyp_integral3), sizeof(imag(hyp_integral3)));
      eta += d_eta;
    }
    cerr << "  ikr = " << ikr << " finished." << endl;
  }

  outfile.close();
  cout << "Output file closed." << endl;
  return 0;
}
