#include "plotter.h"
//const int colors[18]  = {1, 2, 3, 4, 5, 6, 7, 46, 12, 28, 8, 9, 29, 30, 32, 38, 41, 42};
//const int markers[10] = {20, 21, 22, 23, 34, 24, 25, 26, 32, 28};
using namespace std;

Plotter::Plotter()
{
  gStyle = new TStyle();
  gStyle->SetHistLineWidth(2.);
  gStyle->SetFuncWidth(2.);
  c1 = NULL;
  c1 = new TCanvas("c1", "Canvas", 10, 10, 700, 500);
}

Plotter::Plotter(const int wx, const int wy)
{
  gStyle = new TStyle();
  gStyle->SetHistLineWidth(2.);
  gStyle->SetFuncWidth(2.);
  c1 = NULL;
  c1 = new TCanvas("c1", "Canvas", 10, 10, wx, wy);
  c1->SetRightMargin(0.02);
}

Plotter::~Plotter()
{
  delete gStyle;
  if(c1) delete c1;
}

void Plotter::Update()
{
  c1->Update();
  c1->Modified();
}

void Plotter::plot(const char* path, TNamed* named)
{
//  stringstream gifname("");
//  gifname << path << "/" << named->GetName() << ".gif";
//  c1->Print(gifname.str().c_str());
  stringstream pngname("");
  pngname << path << "/" << named->GetName() << ".pdf";
  c1->Print(pngname.str().c_str());
  c1->Clear();
}

void Plotter::plot(const char* path, const char* name)
{
//  stringstream gifname("");
//  gifname << path << "/" << name << ".gif";
//  c1->Print(gifname.str().c_str());
  stringstream pngname("");
  pngname << path << "/" << name << ".pdf";
  c1->Print(pngname.str().c_str());
  c1->Clear();
}

void Plotter::SetLog(bool isLogx, bool isLogy, bool isLogz)
{
  c1->SetLogx(isLogx ? 1 : 0);
  c1->SetLogy(isLogy ? 1 : 0);
  c1->SetLogz(isLogz ? 1 : 0);
}

void Plotter::SetHistTitleSize(const float TitleH, const float TitleW, const float TitleX)
{
  gStyle->SetTitleX(TitleX);
  gStyle->SetTitleY(.99);
  gStyle->SetTitleW(TitleW);
  gStyle->SetTitleH(TitleH);
  c1->SetTopMargin(TitleH + .02);
}

void Plotter::SetHistLabelsX(const float TitleH, const float LabelH)
{
  gStyle->SetLabelSize(1.25 * LabelH, "X");
  gStyle->SetLabelOffset(LabelH * 0.05, "X");
  gStyle->SetTitleSize(TitleH, "X");
  gStyle->SetTitleOffset( (1.3 * LabelH + 0.01) / 1.6 / TitleH + 0.32, "X");
  c1->SetBottomMargin(LabelH + TitleH + .02);
}

void Plotter::SetHistLabelsY(const float TitleH, const float LabelH, const float LeftMargin)
{
  gStyle->SetLabelSize(1.25 * LabelH, "Y");
  gStyle->SetLabelOffset(0.05 * LabelH, "Y");
  gStyle->SetTitleSize(TitleH, "Y");
  gStyle->SetTitleOffset(c1->GetWh() * 1. / c1->GetWw() * (1.5 * LabelH / 1.6 / TitleH + 0.35), "Y");
  c1->SetLeftMargin((c1->GetWh() * 1. / c1->GetWw() * (LabelH + TitleH + .02)));
//  c1->SetLeftMargin(LeftMargin < 0. ? (c1->GetWh() / c1->GetWw() * (LabelH + TitleH + .03)) : LeftMargin);
}

void Plotter::SetMargins(const float MarginTop, const float MarginRight, const float MarginBottom, const float MarginLeft)
{
  c1->SetTopMargin(MarginTop);
  c1->SetRightMargin(MarginRight);
  c1->SetBottomMargin(MarginBottom);
  c1->SetLeftMargin(MarginLeft);
}

/*
Labels, axes...

LabelX:
  absolute box height = 0.8 * LabelSize() * CanvasH, 
  upper box edge from axis (absolute): Offset * CanvasH
TitleX:
  absolute box height = CanvasH * TitleSize()
  middle of box from axis: TitleOffset() * 1.6 * TitleSize() * CanvasH
LabelY:
  absolute box height = 0.8 * LabelSize() * CanvasH, 
  right box edge from axis (absolute): Offset * CanvasW
TitleY:
  absolute box height = CanvasH * TitleSize()
  middle of box from axis: TitleOffset() * 1.6 * CanvasW * TitleSize()
  
Inside boxes, texts resize themselves if Canvas too asymmetric...

*/
