#include "MFit.h"
#include "MHist.h"
#include "MRootIO.h"

void JpsiFit_Example() {
  TFile *file_input =
      new TFile("/home/szhu/work/alice/analysis/QA/test/jpsi_fit_results.root");
  TFile *file_input2 = new TFile("/home/szhu/work/alice/analysis/QA/input/jpsi/"
                                 "JpsiQA_LHC22pass4_dqfilter.root");
  auto fPosZ_fMass_fPT_NumContribCalib_Binned = MRootIO::GetObjectDiectly<THnD>(
      file_input2, "fPosZ_fMass_fPT_NumContribCalib_Binned");

  MHnTool hn_tool(fPosZ_fMass_fPT_NumContribCalib_Binned);
  hn_tool.Rebin(3, 2);
  auto hist = hn_tool.Project(1, {0, 1, 1});
  TF1 *signal = MRootIO::GetObjectDiectly<TF1>(file_input, "signal");
  TF1 *bkg = MRootIO::GetObjectDiectly<TF1>(file_input, "bkg");
  signal->SetNormalized(true);
  bkg->SetNormalized(true);

  gPublisherCanvas =
      new MPublisherCanvas("./jpsi_fit_public.pdf", 1, 1, 600, 600);
  (*gPublisherCanvas) << hist << signal <<= bkg;
  MRootGraphic::StyleCommon();

  MSignalFit signal_fit("JpsiFit", signal, bkg);
  signal_fit << hist;
  signal_fit.chi2Fit();
  signal_fit.chi2Fit(false);
  signal_fit >> (gPublisherCanvas->NewPad());

  gPublisherCanvas->finalize();
}