#include "MHist.h"
#include "MRootGraphic.h"

void CorrelationQA(
    TString path_input =
        "/home/szhu/work/alice/analysis/QA/input/event_jpsi/AssoYeild_24pass1/"
        "AssoYeildGroupEtagapPt_bInt_f_BS_bs.root",
    // TString path_output = "/home/szhu/work/alice/analysis/QA/test/"
    //                       "AssoYeildGroupEtagapPt_bLow.root",
    TString path_pdf = "/home/szhu/work/alice/analysis/QA/plot/event_jpsi/"
                       "LHC24_pass1_DiElectron/"
                       "AssoYeildGroupEtagapPt_bInt_f_BS_bs.pdf") {
  TFile *file_input = new TFile(path_input);

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

  gPublisherCanvas = new MPublisherCanvas(path_pdf, 1, 1, 600, 600);

  for (auto iPtV2 : indexAnyPtV2Jpsi) {
    auto vec_iPtV2 = indexAnyPtV2Jpsi[iPtV2 - 1];

    if (iPtV2 != 1)
      gPublisherCanvas->SetCanvasNwNh(1, 1);
    auto h2_v22 = hgroupTool2d_v22->GetHist(vector<int>{iPtV2});
    gPublisherCanvas->Draw(h2_v22);

    gPublisherCanvas->SetCanvasNwNh(3, 2);
    TString text_bins = "p_{T} bins: ";
    for (auto bin : vec_iPtV2) {
      text_bins += Form("%d ", bin);
    }
    gPublisherCanvas->AddText(text_bins, 0., 0.);

    for (auto iEtaGap : indexHistEtaGap) {
      auto h1_v22 = h2_v22->ProjectionX(Form("h1_v22_%d_%d", iPtV2, iEtaGap),
                                        iEtaGap, iEtaGap);
      h1_v22->GetYaxis()->SetRangeUser(-0.05, 0.05);
      gPublisherCanvas->DrawClone(h1_v22);
    }
  }

  gPublisherCanvas->finalize();
}