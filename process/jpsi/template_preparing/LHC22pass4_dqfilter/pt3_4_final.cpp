#include "MALICE.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "RooAddPdf.h"
#include "RooCrystalBall.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "TApplication.h"
#include "TColor.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TStyle.h"
#include "math.h"

void pt3_4_final() {
  gROOT->SetBatch(true);
  TFile *file_input = new TFile("/home/szhu/work/alice/analysis/QA/input/jpsi/"
                                "JpsiQA_LHC22pass4_dqfilter.root");
  auto fPosZ_fMass_fPT_NumContribCalib_Binned = MRootIO::GetObjectDiectly<THnD>(
      file_input, "fPosZ_fMass_fPT_NumContribCalib_Binned");

  MHnTool hn_tool(fPosZ_fMass_fPT_NumContribCalib_Binned);
  hn_tool.PrintAllAxis();
  hn_tool.Rebin(0, 200);
  hn_tool.Rebin(3, 2);
  hn_tool.PrintAllAxis();

  gPublisherCanvas = new MPublisherCanvas("./pt3_4_final.pdf", 1, 1, 600, 600);

  TF1 *signal =
      new TF1("signal",
              "ROOT::Math::crystalball_function(x,[Alpha],[N],[Sigma],[Mean])",
              1.55, 5.);
  cout << signal->GetExpFormula() << endl;
  signal->SetParameters(2.8186e-01, 3.0632e+00, 1.9479e+00, 4.4484e-02);
  // signal->SetParameters(5.12185, 3.06213, 0.0491487, 0.234623);
  gPublisherCanvas->Draw(signal);

  TF1 *bkg = new TF1(
      "bkg", "pow(x - [bkg0], [bkg1]) / (exp([bkg2] * x) + [bkg3])", 1.5, 5.);

  // 3  p2           1.08717e+00   8.70323e-03   8.34182e-06  -3.47160e-02
  //    4  p3           1.47485e+00   2.43288e-01   2.78783e-04  -2.33890e-03
  //    5  p4          -1.16754e-05   5.48934e-02   4.08062e-05  -3.23054e+01
  //    6  p5           1.26424e+00   2.90258e-02   2.24681e-05   1.16494e-02

#define parfast -1.16754e-05, 1.26424e+00, 1.08717e+00, 1.47485e+00
  // #define parfast 0.8, 10.5, -6.3, 0., -0.448589
  bkg->SetParameters(parfast);
  gPublisherCanvas->Draw(bkg);

  TF1 *bkg_qa =
      new TF1("bkg_qa",
              "exp([bkg0]*x) +exp([bkg1]+[bkg2]/pow(x-[bkg3],[bkg4]))", 1., 5.);
  bkg_qa->SetParameters(parfast);
  gPublisherCanvas->Draw(bkg_qa);
  gPad->SetLogy();

  auto lFormulaTotal = [signal, bkg](double *x, double *par) {
    return par[0] * signal->Eval(x[0]) + par[1] * bkg->Eval(x[0]);
  };

  auto lFormulaTotal2 = [signal, bkg](double *x, double *par) {
    double x1 = x[0];
    double bkg0 = par[2];
    double bkg1 = par[3];

    return par[0] * signal->Eval(x1) +
           par[1] * (1 + bkg0 * x1 + bkg1 * x1 * x1);
  };

  auto lFormulaTotal3 = [signal, bkg](double *x, double *par) {
    double x1 = x[0];
    double bkg0 = par[2];
    double bkg1 = par[3];
    double bkg2 = par[4];
    double bkg3 = par[5];
    // double bkg4 = par[6];
    return par[0] * signal->Eval(x1) +
           par[1] * pow(x1 - bkg2, bkg3) / (exp(bkg0 * x1) + bkg1);
  };
  gPublisherCanvas->SetCanvasNwNh(3, 1);
  for (int i_mult = 1; i_mult <= 5; i_mult++) {
    auto h1_lowMult = hn_tool.Project(1, {0, 4, i_mult});
    TF1 *f1_total =
        new TF1(Form("f1_total_%d", i_mult), lFormulaTotal, 1., 5., 2);
    f1_total->SetParNames("a", "b");
    h1_lowMult->Fit(f1_total, "0M", "same", 1., 5.);
    gPublisherCanvas->Draw(h1_lowMult)->DrawSame(f1_total);

    double n_sig = f1_total->GetParameter(0);
    double n_bkg = f1_total->GetParameter(1);
    TF1 *f1_total_2 =
        new TF1(Form("f1_total_2_%d", i_mult), lFormulaTotal2, 1., 5., 8);
    f1_total_2->SetParameters(n_sig, n_bkg, parfast, 1.);
    f1_total_2->FixParameter(0, n_sig);
    h1_lowMult->Fit(f1_total_2, "0M", "same", 1.5, 5.);
    gPublisherCanvas->Draw(h1_lowMult)->DrawSame(f1_total_2);

    TF1 *f1_total_3 =
        new TF1(Form("f1_total_3_%d", i_mult), lFormulaTotal3, 1.5, 5., 6);
    f1_total_3->SetParameters(n_sig, n_bkg);
    f1_total_3->FixParameter(0, n_sig);
    h1_lowMult->Fit(f1_total_3, "0M", "same", 1.5, 5.);
    gPublisherCanvas->Draw(h1_lowMult)->DrawSame(f1_total_3);

    h1_lowMult->SetTitle(Form("pT: 1, Mult: %d", i_mult));
  }
  gPublisherCanvas->finalize();

  TFile *file = new TFile("/home/szhu/work/alice/analysis/QA/output/jpsi/"
                          "fit_template/LHC22pass4_dqfilter/pt3_4.root",
                          "RECREATE");
  file->cd();
  signal->Write();
  bkg->Write();
  file->Close();
}