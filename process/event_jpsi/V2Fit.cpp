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

void V2Fit() {
  gPublisherCanvas = new MPublisherCanvas(
      "/home/szhu/work/alice/analysis/QA/plot/event_jpsi/V2FitTest.pdf", 1, 1,
      600, 600);
  MRootGraphic::StyleCommon();
  TFile *file_v2 = new TFile("/home/szhu/work/alice/analysis/QA/output/"
                             "event_jpsi/AssoYieldGroupQAEtaGap.root");

  auto MassUS_EtaGap_v22 =
      MRootIO::GetObjectDiectly<TH2D>(file_v2, "MassUS_EtaGap_v22");
  gPublisherCanvas->Draw(MassUS_EtaGap_v22);

  gPublisherCanvas->SetCanvasNwNh(3, 2);

  StrVar4Hist var_EtaGap("EtaGap", "#Delta#eta_{gap}", "", 6, {-0.4, 2.});
  StrVar4Hist var_MassJpsiCandidate("MassUS", "M_{ee}", "GeV^{2}/c^{4}", 90,
                                    {1.8, 5.4});
  MIndexHist indexHistEtaGap(var_EtaGap, 1, 1);
  MIndexHist indexHistMass(var_MassJpsiCandidate, 1, 2);
  MVec<MHist1D> vec_v22(indexHistEtaGap);

  for (int i_etaGap = 1; i_etaGap <= MassUS_EtaGap_v22->GetNbinsY();
       i_etaGap++) {
    auto h1_v22 = MassUS_EtaGap_v22->ProjectionX(Form("h1_v22_%d", i_etaGap),
                                                 i_etaGap, i_etaGap);
    h1_v22->GetYaxis()->SetTitle("V_{2}");
    h1_v22->SetTitle(
        Form("V_{2} vs M_{ee} for #Delta#eta_{gap} = %.2f",
             MassUS_EtaGap_v22->GetYaxis()->GetBinUpEdge(i_etaGap)));
    MRootGraphic::StyleHistCommonHist(h1_v22);
    vec_v22.fVec.emplace_back(indexHistMass, h1_v22);
  }

  ////////////////////////// fitting part /////////////////////////////
  TFile *file_input2 = new TFile("/home/szhu/work/alice/analysis/QA/input/jpsi/"
                                 "JpsiQA_LHC22pass4_dqfilter.root");
  auto vec_arr_template =
      GetTemplate("/home/szhu/work/alice/analysis/QA/output/jpsi/fit_template/"
                  "LHC22pass4_dqfilter");
  gTemplate = &vec_arr_template;
  auto fPosZ_fMass_fPT_NumContribCalib_Binned = MRootIO::GetObjectDiectly<THnD>(
      file_input2, "fPosZ_MassUS_PtUS_NumContribCalib");

  MHnTool hn_tool(fPosZ_fMass_fPT_NumContribCalib_Binned);
  hn_tool.Rebin(3, 5);
  hn_tool.PrintAllAxis();

  auto hist = hn_tool.Project(1, {0, 0, 0});

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

  gPublisherCanvas->SetCanvasNwNh(1, 1);
  signal_fit_total.InputData(hist);
  chi2Fit(signal_fit_total);
  signal_fit_total >> (gPublisherCanvas->NewPad());

  auto gr_s2b = signal_fit_total.GetSignalToBackgroundCurve(100);
  gPublisherCanvas->Draw(gr_s2b);
  auto gr_s2t = signal_fit_total.GetSignalFractionCurve(100);
  gPublisherCanvas->Draw(gr_s2t);
  auto gr_b2t = signal_fit_total.GetBkgFractionCurve(100);
  gPublisherCanvas->Draw(gr_b2t);

  gPublisherCanvas->SetCanvasNwNh(3, 2);
  for (auto i_etaGap : indexHistEtaGap) {
    auto h1_v22 = vec_v22.currentObject();
    h1_v22.fHisto->GetXaxis()->SetRangeUser(1.4, 5.);
    h1_v22.fHisto->SetTitle(Form("V_{2} vs M_{ee} for #Delta#eta_{gap} = %.2f",
                                 indexHistEtaGap.GetBinUpperEdge()));
    gPublisherCanvas->DrawClone(h1_v22);
  }

  double jpsiV2 = 2e-3;
  for (auto i_etaGap : indexHistEtaGap) {
    auto h1_v22 = vec_v22.currentObject();
    for (auto i_etaMass : indexHistMass) {
      auto value = h1_v22.GetBinInfo();
      double center_bin = indexHistMass.GetBinCenter();
      double s2t = gr_s2t->Eval(center_bin);
      MDouble value_new(value.fValue - jpsiV2 * s2t, value.fError);
      h1_v22.SetBinInfo(value_new);
    }
    gPublisherCanvas->DrawClone(h1_v22);
  }

  gPublisherCanvas->finalize();
}
