#include "MHist.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "Math/ChebyshevPol.h"

TF1 *f1_signal;

double mass_total(double *x, double *par) {
  return par[0] * f1_signal->Eval(x[0]) + par[1] + par[2] * x[0] +
         par[3] * x[0] * x[0] + par[4] * x[0] * x[0] * x[0];
}

double ratio_s2t(double *x, double *par) {
  return par[0] * f1_signal->Eval(x[0]) /
         (par[0] * f1_signal->Eval(x[0]) + par[1] + par[2] * x[0] +
          par[3] * x[0] * x[0] + par[4] * x[0] * x[0] * x[0]);
}

void CorrelationQA(
    TString path_input = "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
                         "AssoYeild_24pass1/mass_rebin/"
                         "AssoYeildGroupEtagapPt_bInt_f_BS_bs.root",
    TString path_input_mass = "/home/szhu/work/alice/analysis/QA/input/jpsi/"
                              "JpsiMass_LHC24_apass1_DiElectron.root",
    TString path_pdf = "/home/szhu/work/alice/analysis/QA/plot/event_jpsi/"
                       "LHC24_pass1_DiElectron/"
                       "AssoYeildGroupEtagapPt_bInt_f_BS_bs_massRebin.pdf") {
  TFile *file_input = new TFile(path_input);
  TFile *file_input_mass = new TFile(path_input_mass);

  f1_signal = new TF1(Form("signal"),
                      "ROOT::Math::crystalball_function(x,[Alpha],[N],["
                      "Sigma],[Mean])",
                      2., 4.);
  f1_signal->SetParameters(0.234623, 3.06213, 5.12185, 0.0491487);
  // f1_signal->SetParameters(2.8186e-01, 3.0632e+00, 1.9479e+00, 4.4484e-02);

  auto h_mass =
      (THnD *)file_input_mass->Get("fPosZ_MassUS_PtUS_NumContribCalib");
  MHnTool hnTool_mass(h_mass);
  hnTool_mass.Rebin(0, 200);
  hnTool_mass.Rebin(3, 5);
  hnTool_mass.PrintAllAxis();

  gPublisherCanvas = new MPublisherCanvas(path_pdf, 1, 1, 600, 600);
  TH2D *mass_pt = hnTool_mass.Project(1, 2, {0, 2});
  MRootGraphic::StyleHistCommonHist(mass_pt);
  MRootGraphic::StyleCommon();
  gPublisherCanvas->Draw(mass_pt);
  gPublisherCanvas->Draw(f1_signal);

  gPublisherCanvas->SetCanvasNwNh(5, 4);

  for (int i_pt = 1; i_pt <= 10; i_pt++) {
    auto h = mass_pt->ProjectionY(Form("_pt%d", i_pt), i_pt, i_pt);
    gPublisherCanvas->NewPad()->cd();
    MRootGraphic::StyleHistCommonHist(h);
    double max = h->GetMaximum();
    h->GetYaxis()->SetRangeUser(0, 1.2 * max);
    h->GetXaxis()->SetRangeUser(2., 4.);
    h->SetTitle(Form("p_{T} bin: %d", i_pt));

    TF1 f_total("f_total", mass_total, 2., 4., 5);
    f_total.SetParameters(5000., 0.1, -0.1, 0.01, -0.001);
    h->Fit(&f_total, "", "", 2., 4.);
    h->Draw();

    TF1 *f_ratio = new TF1("f_ratio", ratio_s2t, 2., 4., 5);
    f_ratio->SetParameters(f_total.GetParameters());
    gPublisherCanvas->Draw(f_ratio);
  }

  struct StrAny_ptV2 {
    const vector<vector<int>> bins = {{1},
                                      {2},
                                      {3},
                                      {4},
                                      {5},
                                      {2, 3},
                                      {4, 5},
                                      {1, 2, 3},
                                      {1, 2, 3, 4, 5},
                                      {6, 7, 8, 9, 10},
                                      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}};
    const int fNbins = bins.size();
    const TString fName = "ptV2";

    vector<int> operator[](int index) { return bins[index]; }
  } strAny_ptV2;

  StrVar4Hist var_PtV2Jpsi("PtV2Jpsi", "p_{T}", "GeV/c", strAny_ptV2.fNbins,
                           {0., 1.});
  StrVar4Hist var_EtaGap("EtaGap", "#Delta#eta_{gap}", "", 6, {-0.4, 2.});
  MIndexAny indexAnyPtV2Jpsi(strAny_ptV2, 1);
  MIndexHist indexHistEtaGap(var_EtaGap, 1, 1);

  MHGroupTool2D *hgroupTool2d_a0PlusB = new MHGroupTool2D(
      file_input, "MassUS_EtaGap_a0PlusB_ptV2_%d", {var_PtV2Jpsi}, {1});
  MHGroupTool2D *hgroupTool2d_a0 = new MHGroupTool2D(
      file_input, "MassUS_EtaGap_a0_ptV2_%d", {var_PtV2Jpsi}, {1});
  MHGroupTool2D *hgroupTool2d_a1 = new MHGroupTool2D(
      file_input, "MassUS_EtaGap_a1_ptV2_%d", {var_PtV2Jpsi}, {1});
  MHGroupTool2D *hgroupTool2d_a2 = new MHGroupTool2D(
      file_input, "MassUS_EtaGap_a2_ptV2_%d", {var_PtV2Jpsi}, {1});
  MHGroupTool2D *hgroupTool2d_a3 = new MHGroupTool2D(
      file_input, "MassUS_EtaGap_a3_ptV2_%d", {var_PtV2Jpsi}, {1});
  MHGroupTool2D *hgroupTool2d_b = new MHGroupTool2D(
      file_input, "MassUS_EtaGap_b_ptV2_%d", {var_PtV2Jpsi}, {1});
  MHGroupTool2D *hgroupTool2d_v22 = new MHGroupTool2D(
      file_input, "MassUS_EtaGap_v22_ptV2_%d", {var_PtV2Jpsi}, {1});

  for (auto iPtV2 : indexAnyPtV2Jpsi) {
    auto vec_iPtV2 = indexAnyPtV2Jpsi[iPtV2 - 1];

    gPublisherCanvas->SetCanvasNwNh(2, 1);
    auto h2_v22 = hgroupTool2d_v22->GetHist(vector<int>{iPtV2});
    // gPublisherCanvas->Draw(h2_v22);
    auto bins_pt = indexAnyPtV2Jpsi[iPtV2 - 1];
    auto mass_ptSpecific =
        mass_pt->ProjectionY(Form("_pT_%d_%d", bins_pt.front(), bins_pt.back()),
                             bins_pt.front(), bins_pt.back());
    MRootGraphic::StyleHistCommonHist(mass_ptSpecific);
    mass_ptSpecific->GetXaxis()->SetRangeUser(2., 4.);
    double max = mass_ptSpecific->GetMaximum();
    mass_ptSpecific->GetYaxis()->SetRangeUser(0, 1.2 * max);
    gPublisherCanvas->NewPad()->cd();
    TF1 f_total("f_total", mass_total, 2., 4., 5);
    f_total.SetParameters(5000., 0.1, -0.1, 0.01, -0.001);
    mass_ptSpecific->Fit(&f_total, "", "", 2., 4.);
    mass_ptSpecific->Draw();
    TF1 *f_ratio = new TF1("f_ratio", ratio_s2t, 2., 4., 5);
    f_ratio->SetParameters(f_total.GetParameters());
    gPublisherCanvas->Draw(f_ratio);

    // TF1 *signal = new TF1(
    //     Form("signal_ptV2_%d", iPtV2),
    //     "ROOT::Math::crystalball_function(x,[Alpha],[N],[Sigma],[Mean])", 2.,
    //     4.);
    // signal->SetParameters(2.8186e-01, 3.0632e+00, 1.9479e+00, 4.4484e-02);

    // TF1 *bkg = new TF1(Form("bkg_ptV2_%d", iPtV2), "pol4", 2., 4.);

    // TF1 *total = new TF1(Form("total_ptV2_%d", iPtV2), lTotalFit, 2., 4., 6);

    TString text_bins = "p_{T} bins: ";
    for (auto bin : vec_iPtV2) {
      text_bins += Form("%d ", bin);
    }
    gPublisherCanvas->AddText(text_bins, 0., 0.);
    gPublisherCanvas->SetCanvasNwNh(3, 2);
    for (auto iEtaGap : indexHistEtaGap) {
      auto h1_v22 = h2_v22->ProjectionX(Form("h1_v22_%d_%d", iPtV2, iEtaGap),
                                        iEtaGap, iEtaGap);
      // h1_v22->GetYaxis()->SetRangeUser(-0.05, 0.05);
      h1_v22->GetXaxis()->SetRangeUser(2., 4.);
      gPublisherCanvas->DrawClone(h1_v22);
    }
  }

  gPublisherCanvas->finalize();
}