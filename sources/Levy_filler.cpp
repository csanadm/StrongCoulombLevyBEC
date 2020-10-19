#include <functions.h>
#include <ctime>
//#define NJOBS 400
#define NJOBS 3000
using namespace std;

int main(int argc, char** argv)
{
  if(argc != 2)
  {
    cout << "Usage: " << argv[0] << " ijob (out of " << NJOBS << ")" << endl;
    return 1;
  }

//  double xmin = atof(argv[1]);
//  double xmax = atof(argv[2]);
//  int Nx   = atoi(argv[3]);
//  double alphamin = atof(argv[4]);
//  double alphamax = atof(argv[5]);
//  int Nalpha = atoi(argv[6]);

  int ijob = atoi(argv[1]);
  if(ijob < 0 || ijob >= NJOBS)
  {
    cout << "Wrong job number!" << endl;
    return 1;
  }

  double xmin = 0.;
  double xmax = 300.;  // if kr=900, xmax=60 is enough
  int Nx   = 60000.;
  double alphamin = 0.5;
  double alphamax = 2.;
  int Nalpha = 300;

  cout << "Filling levy table: " << endl;
  cout << alphamin << "<=alpha<" << alphamax << ", " << Nalpha << "+1 steps, for each alpha for " << xmin <<  "<=x<" << xmax << ", " << Nx << "+1 steps." << endl;

  ofstream outfile;
  stringstream outfile_name(""); outfile_name << "Levy_table_job" << ijob << ".dat";
  outfile.open(outfile_name.str().c_str(), ios::out | ios::binary);
  if(ijob == 0)
  {
    outfile.write((char*)&alphamin , sizeof(alphamin));
    outfile.write((char*)&alphamax , sizeof(alphamax)); 
    outfile.write((char*)&Nalpha   , sizeof(Nalpha)); 
    outfile.write((char*)&xmin , sizeof(xmin));
    outfile.write((char*)&xmax , sizeof(xmax)); 
    outfile.write((char*)&Nx   , sizeof(Nx)); 
  }
  
  int full_comp_time;
  for(int ix = 0; ix < (Nx + 1); ix++)  
  {
    clock_t begin = clock();
    if(ix % NJOBS != ijob) continue;
    double x = xmin + (xmax - xmin) * ix / Nx;
    outfile.write((char*)&ix, sizeof(ix));
    for(int ialpha = 0; ialpha < (Nalpha + 1); ialpha++)  
    {
      double alpha = alphamin + (alphamax - alphamin) * ialpha / Nalpha; 
      double levy = levy_calc(x, 1., alpha);
      outfile.write((char*)&levy, sizeof(levy));
    }
    clock_t end = clock();
    int seconds = (end-begin)/CLOCKS_PER_SEC;
    cerr << "  ix = " << ix << " finished in " << ( seconds / 60 ) << " minutes and " << ( seconds % 60 ) << " seconds" << endl;
    full_comp_time += seconds;
  }
  cout << "Full computation time: " << ( full_comp_time / 3600 ) << " hours " << ( full_comp_time % 3600 ) << " minutes " << endl;

  outfile.close();
  cout << "Output file closed." << endl;
  return 0;
}

