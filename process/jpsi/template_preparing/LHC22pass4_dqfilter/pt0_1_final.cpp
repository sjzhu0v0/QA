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

void pt0_1_final() {
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

  gPublisherCanvas = new MPublisherCanvas("./mass_pt0_1.pdf", 1, 1, 600, 600);

  TF1 *signal =
      new TF1("signal",
              "ROOT::Math::crystalball_function(x,[Alpha],[N],[Sigma],[Mean])",
              1.55, 5.);
  cout << signal->GetExpFormula() << endl;
  signal->SetParameters(0.234623, 3.06213, 5.12185, 0.0491487);
  // signal->SetParameters(5.12185, 3.06213, 0.0491487, 0.234623);
  gPublisherCanvas->Draw(signal);

  TF1 *bkg = new TF1(
      "bkg", "1./(exp([bkg0]*x) + [bkg1]*exp([bkg2]/pow(x-[bkg3],[bkg4])))",
      1.55, 5.);
  cout << bkg->GetExpFormula() << endl;
  bkg->SetParameters(1.39635, exp(9.7430275), -6.34629, 1.54991, -0.448589);

  auto lFormulaTotal = [signal, bkg](double *x, double *par) {
    return par[0] * signal->Eval(x[0]) + par[1] * bkg->Eval(x[0]);
  };

  for (int i_mult = 1; i_mult <= 5; i_mult++) {
    gPublisherCanvas->NewPad()->cd();
    auto h1_lowMult = hn_tool.Project(1, {0, 1, i_mult});
    TF1 *f1_total =
        new TF1(Form("f1_total_%d", i_mult), lFormulaTotal, 1.55, 5., 2);
    f1_total->SetParNames("a", "b");
    // f1_total->GetParName(0);
    h1_lowMult->Fit(f1_total, "M", "same", 1.55, 5.);
    h1_lowMult->SetTitle(Form("pT: 1, Mult: %d", i_mult));
  }
  gPublisherCanvas->finalize();

  TFile *file = new TFile("/home/szhu/work/alice/analysis/QA/output/jpsi/"
                          "fit_template/LHC22pass4_dqfilter/pt0_1.root",
                          "RECREATE");
  file->cd();
  signal->Write();
  bkg->Write();
  file->Close();
}