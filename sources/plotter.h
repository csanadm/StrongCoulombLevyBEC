#ifndef Plotter_hh
#define Plotter_hh
#include "functions.h"

//const int colors[18]  = {1, 2, 3, 4, 5, 6, 7, 46, 12, 28, 8, 9, 29, 30, 32, 38, 41, 42};
//const int markers[10] = {20, 21, 22, 23, 34, 24, 25, 26, 32, 28};

class Plotter
{
 public:
  Plotter();  
  Plotter(const int wx, const int wy);
  ~Plotter();  

  void Update();
  void plot(const char* path, TNamed* named);
  void plot(const char* path, const char* name);
  void SetHistTitleSize(const float TitleH, const float TitleW = .9, const float TitleX=.5);
  void SetHistLabelsX(const float TitleH, const float LabelH);
  void SetHistLabelsY(const float TitleH, const float LabelH, const float LeftMargin = -1.);
  void SetMargins(const float MarginTop, const float MarginRight, const float MarginBottom, const float MarginLeft);
  void SetLog(bool isLogx, bool isLogy, bool isLogz);
 private:
  TCanvas* c1;
  TStyle* gStyle;
};

#endif // Plotter_hh
