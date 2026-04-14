#include "MFit.h"
#include "MHelper.h"
#include "MMath.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"
#include "yaml-cpp/yaml.h"
#include <algorithm>

// Helper function to parse StrVar4Hist from YAML config
void AssoYieldEtagap(
    TString path_input = "/home/szhu/work/alice/analysis/QA/test/"
                         "AssoYieldGroupEtagap_NoScale.root",
    TString path_output = "/home/szhu/work/alice/analysis/QA/test/"
                          "AssoYieldFit_noScale.root",
    TString path_config = "config_new.yaml") {
  gErrorIgnoreLevel = kWarning;
  YAML::Node config = YAML::LoadFile(path_config.Data());
  
  TFile *file_input = new TFile(path_input);
  TFile *file_output = new TFile(path_output, "RECREATE");

  YAML::Node hist_config = config["hist_binning"];
  
  StrVar4Hist var_fPosZ = ParseStrVar4Hist(hist_config["fPosZ"]);
  StrVar4Hist var_NumContribCalibBinned = ParseStrVar4Hist(hist_config["NumContribCalibBinned"]);
  StrVar4Hist var_MassJpsiCandidate = ParseStrVar4Hist(hist_config["MassJpsiCandidate"]);
  StrVar4Hist var_PtJpsiCandidate = ParseStrVar4Hist(hist_config["PtJpsiCandidate"]);
  StrVar4Hist var_DeltaEtaUS = ParseStrVar4Hist(hist_config["DeltaEtaUS"]);
  StrVar4Hist var_DeltaPhiUS = ParseStrVar4Hist(hist_config["AbsDeltaPhiUS"]);
  
  // Calculate EtaGap binning based on DeltaEtaUS configuration
  int n_bins_deltaEta_assoYield = var_DeltaEtaUS.fNbins;
  double min_deltaEta = var_DeltaEtaUS.fBins[0];
  double max_deltaEta = var_DeltaEtaUS.fBins.back();
  double bin_width_etaGap = (max_deltaEta - min_deltaEta) / n_bins_deltaEta_assoYield;
  int n_bins_etaGap = n_bins_deltaEta_assoYield / 2 - 2;
  std::vector<double> etaGap_bins = {-1.0 * bin_width_etaGap, 
                                      (n_bins_deltaEta_assoYield / 2 - 3) * bin_width_etaGap};
  StrVar4Hist var_EtaGap("EtaGap", "#Delta#eta_{gap}", "", n_bins_etaGap, etaGap_bins);
  
  // strAny_ptV2 structure for custom pt V2 binning (kept hardcoded due to complex structure)
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
  
  StrVar4Hist var_PtV2Jpsi("PtV2Jpsi", "p_{T}", "GeV/c", strAny_ptV2.fNbins, {0., 1.});
  int n_rebin_mass_assoYield =
      config["hist_binning"]["n_rebin_mass_assoYield"].as<int>();
  MIndexHist indexHistMass(var_MassJpsiCandidate, 1, n_rebin_mass_assoYield);
  MIndexHist indexHistPtJpsiCandidate(var_PtJpsiCandidate, 1, 1);
  MIndexHist indexHistDeltaPhiUS(var_DeltaPhiUS, 1, 1);
  int n_rebin_deltaEta_assoYield =
      config["hist_binning"]["n_rebin_deltaEta_assoYield"].as<int>();
  MIndexHist indexHistDeltaEtaUS(var_DeltaEtaUS, 1, n_rebin_deltaEta_assoYield);
  MIndexHist indexHistEtaGap(var_EtaGap, 1, 1);
  MIndexHist indexHistPtV2Jpsi(var_PtV2Jpsi, 1, 1);
  MIndexAny indexAnyPtV2Jpsi(strAny_ptV2, 1);
  gDirectory = nullptr;
  MHGroupTool1D assoYield_sub(
      file_input, "DeltaPhi_AssoYield_sub_DeltaEta_%d_ptV2_%d",
      {var_DeltaEtaUS, var_PtV2Jpsi}, {indexHistDeltaEtaUS, 1});
  MHGroupTool1D assoYield_low(
      file_input, "DeltaPhi_AssoYield_low_DeltaEta_%d_ptV2_%d",
      {var_DeltaEtaUS, var_PtV2Jpsi}, {indexHistDeltaEtaUS, 1});
  MHGroupTool1D assoYield_high(
      file_input, "DeltaPhi_AssoYield_high_DeltaEta_%d_ptV2_%d",
      {var_DeltaEtaUS, var_PtV2Jpsi}, {indexHistDeltaEtaUS, 1});

  MHist1D assoYield_sub_int(indexHistDeltaPhiUS, "AssoYield_sub_int");
  MHist1D assoYield_low_int(indexHistDeltaPhiUS, "AssoYield_low_int");
  MHist1D assoYield_high_int(indexHistDeltaPhiUS, "AssoYield_high_int");

  MHist2D assoYeild_sub_int2(indexHistDeltaEtaUS, indexHistDeltaPhiUS,
                             "AssoYield_sub_int2");
  MHist2D assoYeild_low_int2(indexHistDeltaEtaUS, indexHistDeltaPhiUS,
                             "AssoYield_low_int2");
  MHist2D assoYeild_high_int2(indexHistDeltaEtaUS, indexHistDeltaPhiUS,
                              "AssoYield_high_int2");

  auto assoYield_sub_EtaGap =
      MakeMVec(assoYield_sub_int, indexHistEtaGap, indexAnyPtV2Jpsi);
  auto assoYield_low_EtaGap =
      MakeMVec(assoYield_low_int, indexHistEtaGap, indexAnyPtV2Jpsi);
  auto assoYield_high_EtaGap =
      MakeMVec(assoYield_high_int, indexHistEtaGap, indexAnyPtV2Jpsi);

  auto assoYield_sub_PtV2 = MakeMVec(assoYeild_sub_int2, indexAnyPtV2Jpsi);
  auto assoYield_low_PtV2 = MakeMVec(assoYeild_low_int2, indexAnyPtV2Jpsi);
  auto assoYield_high_PtV2 = MakeMVec(assoYeild_high_int2, indexAnyPtV2Jpsi);

#define MH1DGetBin(...) GetHist(vector<int>{__VA_ARGS__})

  for (auto iPtV2 : indexAnyPtV2Jpsi) {
    for (auto iEtaGap : indexHistEtaGap) {
      for (int i_deltaEta = 1;
           i_deltaEta <= n_bins_deltaEta_assoYield / 2 - iEtaGap + 1;
           i_deltaEta++) {
        assoYield_sub_EtaGap.currentObject().fHisto->Add(
            assoYield_sub.MH1DGetBin(i_deltaEta, iPtV2));
        assoYield_low_EtaGap.currentObject().fHisto->Add(
            assoYield_low.MH1DGetBin(i_deltaEta, iPtV2));
        assoYield_high_EtaGap.currentObject().fHisto->Add(
            assoYield_high.MH1DGetBin(i_deltaEta, iPtV2));
      }
      for (int i_deltaEta = n_bins_deltaEta_assoYield / 2 + iEtaGap - 1;
           i_deltaEta <= n_bins_deltaEta_assoYield; i_deltaEta++) {
        assoYield_sub_EtaGap.currentObject().fHisto->Add(
            assoYield_sub.MH1DGetBin(i_deltaEta, iPtV2));
        assoYield_low_EtaGap.currentObject().fHisto->Add(
            assoYield_low.MH1DGetBin(i_deltaEta, iPtV2));
        assoYield_high_EtaGap.currentObject().fHisto->Add(
            assoYield_high.MH1DGetBin(i_deltaEta, iPtV2));
      }

      double low_edge_deltaPhi = var_DeltaPhiUS.fBins[0];
      double up_edge_deltaPhi = var_DeltaPhiUS.fBins.back();
      
      TF1 func_modulation(
          "f1_modulation", "[a0]+2*([a1]*cos(x)+[a2]*cos(2*x)+[a3]*cos(3*x))",
          low_edge_deltaPhi, up_edge_deltaPhi);
      assoYield_sub_EtaGap.currentObject().fHisto->Fit(
          &func_modulation, "QR", "", low_edge_deltaPhi, up_edge_deltaPhi);
      assoYield_low_EtaGap.currentObject().fHisto->Fit(
          &func_modulation, "QR", "", low_edge_deltaPhi, up_edge_deltaPhi);
      assoYield_high_EtaGap.currentObject().fHisto->Fit(
          &func_modulation, "QR", "", low_edge_deltaPhi, up_edge_deltaPhi);
    }
  }

  for (auto iPtV2 : indexAnyPtV2Jpsi)
    for (auto iDeltaEta : indexHistDeltaEtaUS)
      for (auto iDeltaPhi : indexHistDeltaPhiUS) {
        auto assoYield_sub_temp = assoYield_sub.MH1DGetBin(iDeltaEta, iPtV2);
        auto assoYield_low_temp = assoYield_low.MH1DGetBin(iDeltaEta, iPtV2);
        auto assoYield_high_temp = assoYield_high.MH1DGetBin(iDeltaEta, iPtV2);
        assoYield_sub_PtV2.currentObject().SetBinInfo(
            assoYield_sub_temp->GetBinContent(iDeltaPhi),
            assoYield_sub_temp->GetBinError(iDeltaPhi));
        assoYield_low_PtV2.currentObject().SetBinInfo(
            assoYield_low_temp->GetBinContent(iDeltaPhi),
            assoYield_low_temp->GetBinError(iDeltaPhi));
        assoYield_high_PtV2.currentObject().SetBinInfo(
            assoYield_high_temp->GetBinContent(iDeltaPhi),
            assoYield_high_temp->GetBinError(iDeltaPhi));
      }

#define ModulationMultClass(mult_class)                                        \
  MHist2D a0_##mult_class(indexHistEtaGap, indexHistPtV2Jpsi,                  \
                          "a0_" #mult_class);                                  \
  MHist2D a1_##mult_class(indexHistEtaGap, indexHistPtV2Jpsi,                  \
                          "a1_" #mult_class);                                  \
  MHist2D a2_##mult_class(indexHistEtaGap, indexHistPtV2Jpsi,                  \
                          "a2_" #mult_class);                                  \
  MHist2D a3_##mult_class(indexHistEtaGap, indexHistPtV2Jpsi,                  \
                          "a3_" #mult_class);                                  \
                                                                               \
  for (auto iPtV2 : indexHistPtV2Jpsi)                                         \
    for (auto iEtaGap : indexHistEtaGap) {                                     \
      TF1 *fitFunc = assoYield_##mult_class##_EtaGap[iPtV2 - 1][iEtaGap - 1]   \
                         .fHisto->GetFunction("f1_modulation");                \
      a0_##mult_class.SetBinInfo(fitFunc->GetParameter(0),                     \
                                 fitFunc->GetParError(0));                     \
      a1_##mult_class.SetBinInfo(fitFunc->GetParameter(1),                     \
                                 fitFunc->GetParError(1));                     \
      a2_##mult_class.SetBinInfo(fitFunc->GetParameter(2),                     \
                                 fitFunc->GetParError(2));                     \
      a3_##mult_class.SetBinInfo(fitFunc->GetParameter(3),                     \
                                 fitFunc->GetParError(3));                     \
    }                                                                          \
  file_output->cd();                                                           \
  a0_##mult_class.Write();                                                     \
  a1_##mult_class.Write();                                                     \
  a2_##mult_class.Write();                                                     \
  a3_##mult_class.Write();

  ModulationMultClass(high);
  ModulationMultClass(low);
  ModulationMultClass(sub);

  file_output->cd();
  for (auto iPtV2 : indexAnyPtV2Jpsi) {
    for (auto iEtaGap : indexHistEtaGap) {
      assoYield_sub_EtaGap.currentObject().fHisto->Write();
      assoYield_low_EtaGap.currentObject().fHisto->Write();
      assoYield_high_EtaGap.currentObject().fHisto->Write();
    }
  }

  for (auto _ : indexAnyPtV2Jpsi) {
    assoYield_sub_PtV2.currentObject().fHisto->Write();
    assoYield_low_PtV2.currentObject().fHisto->Write();
    assoYield_high_PtV2.currentObject().fHisto->Write();
  }

  file_output->Close();
}

int main(int argc, char **argv) {
  gROOT->SetBatch(kTRUE);
  TString path_input = argc > 1 ? argv[1]
                                : "/home/szhu/work/alice/analysis/QA/test/"
                                  "AssoYieldGroupEtagap_NoScale.root";
  TString path_output = argc > 2 ? argv[2]
                                 : "/home/szhu/work/alice/analysis/QA/test/"
                                   "AssoYieldFit_noScale.root";
  TString path_config = argc > 3 ? argv[3] : "config_new.yaml";

  AssoYieldEtagap(path_input, path_output, path_config);

  return 0;
}
