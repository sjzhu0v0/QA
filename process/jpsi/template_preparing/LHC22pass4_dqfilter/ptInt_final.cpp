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

void ptInt_final() {
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

  gPublisherCanvas =
      new MPublisherCanvas("./mass_ptInt_final.pdf", 1, 1, 600, 600);

  TF1 *signal =
      new TF1("signal",
              "ROOT::Math::crystalball_function(x,[Alpha],[N],[Sigma],[Mean])",
              1.55, 5.);
  cout << signal->GetExpFormula() << endl;
  signal->SetParameters(2.5858e-01, 3.0586e+00, 3.1221e+00, 5.3778e-02);
  // signal->SetParameters(5.12185, 3.06213, 0.0491487, 0.234623);
  gPublisherCanvas->Draw(signal);

  TF1 *bkg = new TF1(
      "bkg", "1./(exp([bkg0]*x) + exp([bkg1]+[bkg2]/pow(x-[bkg3],[bkg4])))",
      1.2, 5.);
// 3  p2           1.12754e+00   3.97184e-03  -2.98829e-03   1.52103e+03
//    4  p3           4.32104e+01   2.00620e-01   1.14192e-01   3.69477e+02
//    5  p4          -1.39842e+01   1.69123e-01   9.76459e-02   1.06373e+03
//    6  p5          -1.60591e+01   2.80002e-01  -1.56060e-01   2.89888e+02
//    7  p6          -3.60336e-01   4.13338e-03  -2.46396e-03   4.37906e+04
#define parfast 1.0972e+00, 5.2326e+01, -1.0819e+01, -2.4148e+01, -4.6344e-01

  bkg->SetParameters(parfast);
  gPublisherCanvas->Draw(bkg);

  TF1 *bkg_qa = new TF1(
      "bkg_qa", "exp([bkg0]*x) +exp([bkg1]+[bkg2]/pow(x-[bkg3],[bkg4]))", 1.2,
      5.);
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
    double bkg2 = par[4];
    double bkg3 = par[5];
    double bkg4 = par[6];
    return par[0] * signal->Eval(x1) +
           par[1] * (1. / (exp(bkg0 * x1) +
                           exp(bkg1 + bkg2 / pow(x1 - bkg3, bkg4))));
  };

  for (int i_mult = 1; i_mult <= 5; i_mult++) {
    gPublisherCanvas->NewPad()->cd();
    auto h1_lowMult = hn_tool.Project(1, {0, 0, i_mult});
    TF1 *f1_total =
        new TF1(Form("f1_total_%d", i_mult), lFormulaTotal, 1.5, 5., 2);
    f1_total->SetParNames("a", "b");
    h1_lowMult->Fit(f1_total, "M", "same", 1.5, 5.);
    double n_sig = f1_total->GetParameter(0);
    double n_bkg = f1_total->GetParameter(1);
    TF1 *f1_total_2 =
        new TF1(Form("f1_total_2_%d", i_mult), lFormulaTotal2, 1., 5., 7);
    f1_total_2->SetParameters(n_sig, n_bkg, parfast);
    h1_lowMult->Fit(f1_total_2, "M", "same", 1., 5.);

    h1_lowMult->SetTitle(Form("pT: 1, Mult: %d", i_mult));
  }
  gPublisherCanvas->finalize();

  TFile *file = new TFile("/home/szhu/work/alice/analysis/QA/output/jpsi/"
                          "fit_template/LHC22pass4_dqfilter/ptInt.root",
                          "RECREATE");
  file->cd();
  signal->Write();
  bkg->Write();
  file->Close();
}