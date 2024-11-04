#include <functions.h>
#include <HypSInt_reader.h>
#include <Levy_reader.h>

#define K_CONST (493.677 / 2. / 137.036)
#define NJOBS 251
#define HBARC 197.3269788

const double R_cutoff = 20.;

const complex<double> I(0., 1.);
const complex<double> ZERO(0., 0.);

int main(int argc, char** argv)
{
  if(argc != 2)
  {
    cout << "Usage: " << argv[0] << " ijob (out of " << NJOBS << " jobs.)" << endl;
    return 1;
  }
  int ijob = atoi(argv[1]);

  double alpha_min = 0.5;
  double alpha_max = 2.0;
  int Nalpha = 250;//240;
  
//kr = 0 - 6000, 
//rt = r/Rcc, r = rt*Rcc, kr = k*rt*Rcc
//6000 = k*rt*Rcc = 200/HBARC*rt*30 = 30*rt --> rt max = 200
//kr = 0 - 1500, 
//rt = r/Rcc, r = rt*Rcc, kr = k*rt*Rcc
//1500 = k*rt*Rcc = 100/HBARC*rt*30 = 15*rt --> rt max = 100

  double k_min = 1.81;
  double k_max = 200.;
//  int Nk = 8;
  int Nk = 800;//600

  double Rcc_min = 2.;
  double Rcc_max = 29.;
  int NRcc = 540;  // I changed the original ratio to: 0.015 . It gives finer steps then the original (S.L.)
//  int NRcc = 27;  // I changed the original ratio to: 0.015 . It gives finer steps then the original (S.L.)

  double rt_min = 0.;
//  double rt_max = 60.;
//  int Nrt = 12000;
  double rt_max = 200.;//trying out a more precise calc.
  int Nrt = 60000;//26000
  
  if(ijob < 0 || ijob >= NJOBS || ((Nalpha + 1) % NJOBS != 0))//supposed to be Nalpha + 1 % NJOBS. If it is Nalpha % NJOBS, and Nalpha = NJOBS, the last job contains two ialphas.
  {
    cout << "Check again the job distribution!" << endl;
    return 1;
  }

  ofstream outfile; // int dr (Levy *  (int dy abs(F)^2))    is the 1st value,    int dr (Levy * (int dy e^2ik * F1  F2))   is the 2nd
  ofstream logfile; // int dr (Levy *  (int dy abs(F)^2))    is the 1st value,    int dr (Levy * (int dy e^2ik * F1  F2))   is the 2nd
  stringstream outfile_name(""); outfile_name << "CorrFunc_kaons_raw_table_" << ijob << ".dat";
  stringstream logfile_name(""); logfile_name << "log_CorrFunc_kaons_raw_table_" << ijob << ".log";
  outfile.open(outfile_name.str().c_str(), ios::out | ios::binary); 
  logfile.open(logfile_name.str().c_str());
  if(ijob == 0)
  {
    cout << "Job 0, writing aux info." << endl;
    cerr << outfile.tellp() << endl;
    outfile.write((char*)&alpha_min , sizeof(alpha_min));
    outfile.write((char*)&alpha_max , sizeof(alpha_max)); 
    outfile.write((char*)&Nalpha   , sizeof(Nalpha)); 
    outfile.write((char*)&k_min , sizeof(k_min));
    outfile.write((char*)&k_max , sizeof(k_max)); 
    outfile.write((char*)&Nk   , sizeof(Nk)); 
    outfile.write((char*)&Rcc_min , sizeof(Rcc_min));
    outfile.write((char*)&Rcc_max , sizeof(Rcc_max)); 
    outfile.write((char*)&NRcc   , sizeof(NRcc)); 
    cerr << outfile.tellp() << endl;
    cout << "aux writing done." << endl;
  }
  
//  HypSInt_reader* myHypSInt_reader = new HypSInt_reader("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/kincsesd/Strong_Coulomb_Levy/tables/HypSInt_table_mostprecise.dat");
  Levy_reader* myLevy_reader         = new Levy_reader("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/kincsesd/Strong_Coulomb_Levy/tables/Levy_table_mostprecise.dat");
  HypSInt_reader* myHypSInt_reader_C = new HypSInt_reader("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/kincsesd/Strong_Coulomb_Levy/tables/HypSInt_table_final.dat",true);
  cout << "Not reached" << endl;
  
  double alpha = alpha_min;
  double d_alpha = (alpha_max - alpha_min) * 1.0 / Nalpha;
  for(int ialpha = 0; ialpha < (Nalpha + 1); ialpha++)  
  {
    if( (ialpha / ((Nalpha + 1)/ NJOBS) != ijob) && !( (ijob == (NJOBS - 1)) && ialpha == Nalpha ) )
    {
      alpha += d_alpha;
      continue; // performing job distribution here...
    }
    cerr << "ialpha = " << ialpha << endl;
    double k = k_min;
    double d_k = (k_max - k_min) * 1.0 / Nk; 
    for(int ik = 0; ik < (Nk + 1); ik++)  
    {
      cerr << "ik = " << ik << endl;
      double eta = K_CONST / k;
      double Rcc = Rcc_min;
      double d_Rcc = (Rcc_max - Rcc_min) * 1.0 / NRcc; 
      for(int iRcc = 0; iRcc < (NRcc + 1); iRcc++)
      {
        cerr << "iRcc = " << iRcc << endl;
        double rt = rt_min;
        double d_rt = (rt_max - rt_min) * 1.0 / Nrt; 

        double result1 = 0.;
        double result2 = 0.;
        double result3 = 0.;
//        complex<double> result4 = 0.;
      
        for(int irt = 0; irt < (Nrt + 1); irt++)
        {
//          double levy_val = (myLevy_reader->getValue(alpha, rt))*(exp(-rt/R_cutoff));
          double levy_val = (myLevy_reader->getValue(alpha, rt));
          result1 += levy_val * myHypSInt_reader_C->get_F2int_1(eta, (k * rt * Rcc) / HBARC) * SQR(rt) * d_rt;
          result2 += levy_val * myHypSInt_reader_C->get_F2int_2(eta, (k * rt * Rcc) / HBARC) * SQR(rt) * d_rt;
          complex<double> U2_val = ((rt == 0.) ? ZERO : U2(1. - I * eta, 2. * I * (k * rt * Rcc) / HBARC));
          result3 += levy_val * norm(U2_val) * SQR(rt) * d_rt;
//          complex<double> F_int4(myHypSInt_reader->get_F2int_3_real(eta, (k * rt * Rcc) / HBARC), myHypSInt_reader->get_F2int_3_imag(eta, (k * rt * Rcc) / HBARC));
//          result4 += levy_val * U2_val * F_int4 *  SQR(rt) * d_rt;
          rt += d_rt;
        }
//        double result4_real = real(result4);
//        double result4_imag = imag(result4);
        outfile.write((char*)&result1, sizeof(result1));
        outfile.write((char*)&result2, sizeof(result2));
        outfile.write((char*)&result3, sizeof(result3));
//        outfile.write((char*)&result4_real, sizeof(result4_real));
//        outfile.write((char*)&result4_imag, sizeof(result4_imag));
        Rcc += d_Rcc;
      }
      logfile << "ik = " << ik << " finished." << endl;
      k += d_k;
    }
    alpha += d_alpha;
    logfile << "  ialpha = " << ialpha << " finished." << endl;
  }
  delete myHypSInt_reader_C;
  HypSInt_reader* myHypSInt_reader_S = new HypSInt_reader("/gpfs/mnt/gpfs02/phenix/plhf/plhf1/kincsesd/Strong_Coulomb_Levy/tables/HypSInt_table_final.dat",false);
  alpha = alpha_min;
  for(int ialpha = 0; ialpha < (Nalpha + 1); ialpha++)  
  {
    if( (ialpha / ((Nalpha + 1) / NJOBS) != ijob) && !( (ijob == (NJOBS - 1)) && ialpha == Nalpha ) )
    {
      alpha += d_alpha;
      continue; // performing job distribution here...
    }
    cerr << "ialpha = " << ialpha << endl;
    double k = k_min;
    double d_k = (k_max - k_min) * 1.0 / Nk; 
    for(int ik = 0; ik < (Nk + 1); ik++)  
    {
      cerr << "ik = " << ik << endl;
      double eta = K_CONST / k;
      double Rcc = Rcc_min;
      double d_Rcc = (Rcc_max - Rcc_min) * 1.0 / NRcc; 
      for(int iRcc = 0; iRcc < (NRcc + 1); iRcc++)
      {
        cerr << "iRcc = " << iRcc << endl;
        double rt = rt_min;
        double d_rt = (rt_max - rt_min) * 1.0 / Nrt; 

//        double result1 = 0.;
//        double result2 = 0.;
//        double result3 = 0.;
        complex<double> result4 = 0.;
      
        for(int irt = 0; irt < (Nrt + 1); irt++)
        {
          double levy_val = (myLevy_reader->getValue(alpha, rt));
//          double levy_val = (myLevy_reader->getValue(alpha, rt))*(exp(-rt/R_cutoff));
//          result1 += levy_val * myHypSInt_reader->get_F2int_1(eta, (k * rt * Rcc) / HBARC) * SQR(rt) * d_rt;
//          result2 += levy_val * myHypSInt_reader->get_F2int_2(eta, (k * rt * Rcc) / HBARC) * SQR(rt) * d_rt;
          complex<double> U2_val = ((rt == 0.) ? ZERO : U2(1. - I * eta, 2. * I * (k * rt * Rcc) / HBARC));
//          result3 += levy_val * norm(U2_val) * SQR(rt) * d_rt;
          
          complex<double> F_int4(myHypSInt_reader_S->get_F2int_3_real(eta, (k * rt * Rcc) / HBARC), myHypSInt_reader_S->get_F2int_3_imag(eta, (k * rt * Rcc) / HBARC));
          
          result4 += levy_val * U2_val * F_int4 *  SQR(rt) * d_rt;
          rt += d_rt;
        }
        double result4_real = real(result4);
        double result4_imag = imag(result4);
//        outfile.write((char*)&result1, sizeof(result1));
//        outfile.write((char*)&result2, sizeof(result2));
//        outfile.write((char*)&result3, sizeof(result3));
        outfile.write((char*)&result4_real, sizeof(result4_real));
        outfile.write((char*)&result4_imag, sizeof(result4_imag));
        Rcc += d_Rcc;
      }
      logfile << "ik = " << ik << " finished." << endl;
      k += d_k;
    }
    alpha += d_alpha;
    logfile << "  ialpha = " << ialpha << " finished." << endl;
  }
  
  delete myHypSInt_reader_S;
  delete myLevy_reader;
  outfile.close();
  logfile << "Output file closed." << endl;
  logfile.close();
  return 0;
}
