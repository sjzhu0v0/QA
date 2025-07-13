#include "MFit.h"
#include "MHist.h"
#include "MRootIO.h"

vector<array<TF1 *, 2>> *gTemplate;

vector<array<TF1 *, 2>> GetTemplate(TString path_file) {
  TF1 *signal0 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/ptInt.root" + ":signal");
  TF1 *bkg0 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/ptInt.root" + ":bkg");
  signal0->SetNormalized(true);
  bkg0->SetNormalized(true);
  array<TF1 *, 2> arr0 = {signal0, bkg0};

  TF1 *signal1 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/pt0_1.root" + ":signal");
  TF1 *bkg1 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/pt0_1.root" + ":bkg");
  signal1->SetNormalized(true);
  bkg1->SetNormalized(true);
  array<TF1 *, 2> arr1 = {signal1, bkg1};

  TF1 *signal2 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/pt1_2.root" + ":signal");
  TF1 *bkg2 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/pt1_2.root" + ":bkg");
  signal2->SetNormalized(true);
  bkg2->SetNormalized(true);
  array<TF1 *, 2> arr2 = {signal2, bkg2};

  TF1 *signal3 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/pt2_3.root" + ":signal");
  TF1 *bkg3 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/pt2_3.root" + ":bkg");
  signal3->SetNormalized(true);
  bkg3->SetNormalized(true);
  array<TF1 *, 2> arr3 = {signal3, bkg3};

  TF1 *signal4 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/pt3_4.root" + ":signal");
  TF1 *bkg4 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/pt3_4.root" + ":bkg");
  signal4->SetNormalized(true);
  bkg4->SetNormalized(true);
  array<TF1 *, 2> arr4 = {signal4, bkg4};

  return {arr0, arr1, arr2, arr3, arr4};
}

TF1 *ProxyTemplate(int i, int j) {
  int j_temp = j >= 2 ? 1 : j;
  int i_temp = i >= gTemplate->size() ? gTemplate->size() - 1 : i;
  return (*gTemplate)[i_temp][j_temp];
}

void JpsiFitPtMult() {
  TFile *file_input2 = new TFile("/home/szhu/work/alice/analysis/QA/input/jpsi/"
                                 "JpsiQA_LHC22pass4_dqfilter.root");
  auto vec_arr_template =
      GetTemplate("/home/szhu/work/alice/analysis/QA/output/jpsi/fit_template/"
                  "LHC22pass4_dqfilter");
  gTemplate = &vec_arr_template;
  auto fPosZ_fMass_fPT_NumContribCalib_Binned = MRootIO::GetObjectDiectly<THnD>(
      file_input2, "fPosZ_fMass_fPT_NumContribCalib_Binned");

  MHnTool hn_tool(fPosZ_fMass_fPT_NumContribCalib_Binned);
  hn_tool.Rebin(3, 5);
  hn_tool.PrintAllAxis();

  auto hist = hn_tool.Project(1, {0, 0, 0});

  gPublisherCanvas = new MPublisherCanvas(
      "/home/szhu/work/alice/analysis/QA/plot/jpsi/jpsi_fit_pt_mult.pdf", 1, 1,
      600, 600);
  MRootGraphic::StyleCommon();
  // RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  MSignalFit signal_fit_total("JpsiFit1", ProxyTemplate(0, 0),
                              ProxyTemplate(0, 1), 1.6, 5.);

  auto chi2Fit = [](MSignalFit &signal_fit_total) {
    signal_fit_total.FixBkg();
    signal_fit_total.FixSignal();
    signal_fit_total.chi2Fit();
    signal_fit_total.FixBkg(false);
    signal_fit_total.FixSignal();
    signal_fit_total.chi2Fit();
    signal_fit_total.FixBkg(false);
    signal_fit_total.FixSignal(false);
    signal_fit_total.chi2Fit();
    signal_fit_total.RemoveLimit();
    signal_fit_total.chi2Fit();
  };

  auto fit = [](MSignalFit &signal_fit_total) {
    signal_fit_total.FixBkg();
    signal_fit_total.FixSignal();
    signal_fit_total.Fit();
    signal_fit_total.FixBkg(false);
    signal_fit_total.FixSignal();
    signal_fit_total.Fit();
    signal_fit_total.FixBkg(false);
    signal_fit_total.FixSignal(false);
    signal_fit_total.Fit();
    signal_fit_total.RemoveLimit();
    signal_fit_total.Fit();
  };

  signal_fit_total.InputData(hist);
  chi2Fit(signal_fit_total);
  signal_fit_total >> (gPublisherCanvas->NewPad());

  for (int i_pt = 1; i_pt <= hn_tool.GetNbins(2); i_pt++) {
    MSignalFit signal_fit_pt(Form("JpsiFit_pt_%d", i_pt),
                             ProxyTemplate(i_pt, 0), ProxyTemplate(i_pt, 1),
                             1.7, 5.);
    auto hist_pt = hn_tool.Project(1, {0, i_pt, 1});
    signal_fit_pt.InputData(hist_pt);
    chi2Fit(signal_fit_pt);
    signal_fit_pt >> (gPublisherCanvas->NewPad());
    gPublisherCanvas->AddText(Form("pT: %d, Mult: 1", i_pt), 0.55,
                              0.86 - 0.045 * 6, kBlack, 0.04);
    // for (int i_mult = 1.; i_mult <= hn_tool.GetNbins(3); i_mult++) {
    //   MSignalFit signal_fit(Form("JpsiFit_%d_%d", i_pt, i_mult),
    //                         ProxyTemplate(i_pt, 0), ProxyTemplate(i_pt,
    //                         1), 1., 5.);
    //   auto hist_pt = hn_tool.Project(1, {0, i_pt, i_mult});
    //   signal_fit << hist_pt;
    //   signal_fit.chi2Fit();
    //   signal_fit.chi2Fit(false);
    //   signal_fit.RemoveLimit();
    //   signal_fit.chi2Fit(false);
    //   signal_fit >> (gPublisherCanvas->NewPad());
    // }
    gPublisherCanvas->SetCanvasNwNh(2, 1);
    for (int i_mult = 1; i_mult <= hn_tool.GetNbins(3); i_mult++) {
      auto hist_pt_mult = hn_tool.Project(1, {0, i_pt, i_mult});
      gPublisherCanvas->Draw(hist_pt_mult);
    }
    gPublisherCanvas->SetCanvasNwNh(1, 1);
  }
  gPublisherCanvas->finalize();
}
