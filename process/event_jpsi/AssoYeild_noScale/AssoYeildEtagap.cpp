#include "MFit.h"
#include "MHelper.h"
#include "MMath.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"

funcWithJson(void, AssoYeildEtagap)(
    TString path_input = "/home/szhu/work/alice/analysis/QA/test/"
                         "AssoYeildGroupEtagap_NoScale.root",
    TString path_output = "/home/szhu/work/alice/analysis/QA/test/"
                          "AssoYeildFit_noScale.root") {
  gErrorIgnoreLevel = kWarning;
  SetUpJson();
  Configurable<int> config_n_rebin_mass("n_rebin_mass", 3);
  int n_rebin_mass = config_n_rebin_mass.data;

  TFile *file_input = new TFile(path_input);
  TFile *file_output = new TFile(path_output, "RECREATE");

  struct StrAny_ptV2 {
    const vector<vector<int>> bins = {{1},
                                      {2},
                                      {3},
                                      {4},
                                      {5},
                                      {6},
                                      {7},
                                      {8},
                                      {9},
                                      {10},
                                      {1, 2},
                                      {2, 3},
                                      {3, 4},
                                      {4, 5},
                                      {5, 6},
                                      {6, 7},
                                      {7, 8},
                                      {8, 9},
                                      {9, 10},
                                      {1, 2, 3},
                                      {2, 3, 4},
                                      {3, 4, 5},
                                      {4, 5, 6},
                                      {5, 6, 7},
                                      {6, 7, 8},
                                      {7, 8, 9},
                                      {8, 9, 10},
                                      {1, 2, 3, 4},
                                      {2, 3, 4, 5},
                                      {3, 4, 5, 6},
                                      {4, 5, 6, 7},
                                      {5, 6, 7, 8},
                                      {6, 7, 8, 9},
                                      {7, 8, 9, 10},
                                      {1, 2, 3, 4, 5},
                                      {2, 3, 4, 5, 6},
                                      {3, 4, 5, 6, 7},
                                      {4, 5, 6, 7, 8},
                                      {5, 6, 7, 8, 9},
                                      {6, 7, 8, 9, 10},
                                      {1, 2, 3, 4, 5, 6},
                                      {2, 3, 4, 5, 6, 7},
                                      {3, 4, 5, 6, 7, 8},
                                      {4, 5, 6, 7, 8, 9},
                                      {5, 6, 7, 8, 9, 10},
                                      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}};
    const int fNbins = bins.size();
    const TString fName = "ptV2";

    vector<int> operator[](int index) { return bins[index]; }
  } strAny_ptV2;

  StrVar4Hist var_fPosZ("PosZUS", "#it{V}_{Z}", "cm", 8, {-10, 10});
  StrVar4Hist var_NumContribCalibBinned(
      "NumContribCalibUS", "N_{vtx contrib} Calibrated", "", 10,
      {0, 7, 12, 17, 22, 29, 37, 46, 57, 73, 300});
  StrVar4Hist var_MassJpsiCandidate("MassUS", "M_{ee}", "GeV^{2}/c^{4}", 90,
                                    {1.8, 5.4});
  StrVar4Hist var_PtJpsiCandidate("PtUS", "p_{T}", "GeV/c", 10, {0., 10.});
  StrVar4Hist var_DeltaEtaUS("DeltaEtaUS", "#Delta#eta_{J/#psi, track}", "", 80,
                             {-4., 4.});
  StrVar4Hist var_DeltaPhiUS("DeltaPhiUS", "#Delta#phi_{J/#psi, track}", "", 10,
                             {-M_PI_2, M_PI + M_PI_2});
  StrVar4Hist var_EtaGap("EtaGap", "#Delta#eta_{gap}", "", 6, {-0.1, 1.});
  StrVar4Hist var_PtV2Jpsi("PtV2Jpsi", "p_{T}", "GeV/c", strAny_ptV2.fNbins,
                           {0., 1.});

  MIndexHist indexHistMass(var_MassJpsiCandidate, 1, n_rebin_mass);
  MIndexHist indexHistPtJpsiCandidate(var_PtJpsiCandidate, 1, 1);
  MIndexHist indexHistDeltaPhiUS(var_DeltaPhiUS, 1, 1);
  MIndexHist indexHistDeltaEtaUS(var_DeltaEtaUS, 1, 2);
  MIndexHist indexHistEtaGap(var_EtaGap, 1, 1);
  MIndexHist indexHistPtV2Jpsi(var_PtV2Jpsi, 1, 1);
  MIndexAny indexAnyPtV2Jpsi(strAny_ptV2, 1);
  gDirectory = nullptr;
  MHGroupTool1D assoYeild_sub(file_input,
                              "DeltaPhiUS_AssoYeild_sub_DeltaEtaUS_%d_ptV2_%d",
                              {var_DeltaEtaUS, var_PtV2Jpsi}, {2, 1});
  MHGroupTool1D assoYeild_low(file_input,
                              "DeltaPhiUS_AssoYeild_low_DeltaEtaUS_%d_ptV2_%d",
                              {var_DeltaEtaUS, var_PtV2Jpsi}, {2, 1});
  MHGroupTool1D assoYeild_high(
      file_input, "DeltaPhiUS_AssoYeild_high_DeltaEtaUS_%d_ptV2_%d",
      {var_DeltaEtaUS, var_PtV2Jpsi}, {2, 1});

  MHist1D assoYeild_sub_int(indexHistDeltaPhiUS, "AssoYeild_sub_int");
  MHist1D assoYeild_low_int(indexHistDeltaPhiUS, "AssoYeild_low_int");
  MHist1D assoYeild_high_int(indexHistDeltaPhiUS, "AssoYeild_high_int");

  MVec<MHist1D> assoYeild_sub_DeltaEta(indexHistEtaGap, assoYeild_sub_int);
  MVec<MHist1D> assoYeild_low_DeltaEta(indexHistEtaGap, assoYeild_low_int);
  MVec<MHist1D> assoYeild_high_DeltaEta(indexHistEtaGap, assoYeild_high_int);
  MVec<MVec<MHist1D>, MIndexAny<StrAny_ptV2>> assoYeild_sub_EtaGap(
      indexAnyPtV2Jpsi, assoYeild_sub_DeltaEta);
  MVec<MVec<MHist1D>, MIndexAny<StrAny_ptV2>> assoYeild_low_EtaGap(
      indexAnyPtV2Jpsi, assoYeild_low_DeltaEta);
  MVec<MVec<MHist1D>, MIndexAny<StrAny_ptV2>> assoYeild_high_EtaGap(
      indexAnyPtV2Jpsi, assoYeild_high_DeltaEta);

#define MH1DGetBin(...) GetHist(vector<int>{__VA_ARGS__})

  for (auto iPtV2 : indexAnyPtV2Jpsi) {
    for (auto iEtaGap : indexHistEtaGap) {
      for (int i_deltaEta = 11; i_deltaEta <= 20 - iEtaGap + 1; i_deltaEta++)
        assoYeild_sub_EtaGap.currentObject().fHisto->Add(
            assoYeild_sub.MH1DGetBin(i_deltaEta, iPtV2));
      for (int i_deltaEta = 20 + iEtaGap - 1; i_deltaEta <= 30; i_deltaEta++)
        assoYeild_sub_EtaGap.currentObject().fHisto->Add(
            assoYeild_sub.MH1DGetBin(i_deltaEta, iPtV2));
    }
  }

  file_output->cd();
  for (auto iPtV2 : indexAnyPtV2Jpsi) {
    for (auto iEtaGap : indexHistEtaGap) {
      assoYeild_sub_EtaGap.currentObject().fHisto->Write();
      assoYeild_low_EtaGap.currentObject().fHisto->Write();
      assoYeild_high_EtaGap.currentObject().fHisto->Write();
    }
  }
  file_output->Close();
}

int main(int argc, char **argv) {
  gROOT->SetBatch(kTRUE);
  TString path_input = argc > 1 ? argv[1]
                                : "/home/szhu/work/alice/analysis/QA/test/"
                                  "AssoYeildGroupEtagap_NoScale.root";
  TString path_output = argc > 2 ? argv[2]
                                 : "/home/szhu/work/alice/analysis/QA/test/"
                                   "AssoYeildFit_noScale.root";

  AssoYeildEtagap(path_input, path_output);

  return 0;
}