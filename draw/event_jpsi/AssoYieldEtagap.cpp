#include "MHist.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "yaml-cpp/yaml.h"

void AssoYieldEtagap(
    TString path_input = "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
                         "JpsiAssocYield_24apass1.root",
    TString path_output = "/home/szhu/work/alice/analysis/QA/output/event_jpsi/"
                          "AssoYieldEtagap.root") {
  TFile *file_input = new TFile(path_input);
  TFile *file_output = new TFile(path_output, "RECREATE");
  YAML::Node config = YAML::LoadFile("config.yaml");
  // int n_rebin_mass_assoYield =
  //     config["hist_binning"]["n_rebin_mass_assoYield"].as<int>();
  const double low_edge_deltaPhiToPi =
      config["hist_binning"]["low_edge_deltaPhiToPi"].as<double>();
  const double up_edge_deltaPhiToPi =
      config["hist_binning"]["up_edge_deltaPhiToPi"].as<double>();

  // TFile *file_input = new TFile(path_input);

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
  int n_bins_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["n_bins"].as<int>();
  double min_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["min"].as<double>();
  double max_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["max"].as<double>();
  double bin_width_etaGap = (max_deltaEta_assoYield - min_deltaEta_assoYield) /
                            (double)n_bins_deltaEta_assoYield / 2.;
  StrVar4Hist var_DeltaEtaUS("DeltaEtaUS", "#Delta#eta_{J/#psi, track}", "",
                             n_bins_deltaEta_assoYield,
                             {min_deltaEta_assoYield, max_deltaEta_assoYield});
  int n_bins_deltaPhi_assoYield =
      config["hist_binning"]["n_bins_deltaPhi_assoYield"].as<int>();
  StrVar4Hist var_DeltaPhiUS(
      "DeltaPhiUS", "#Delta#phi_{J/#psi, track}", "", n_bins_deltaPhi_assoYield,
      {low_edge_deltaPhiToPi * M_PI, up_edge_deltaPhiToPi * M_PI});
  StrVar4Hist var_EtaGap(
      "EtaGap", "#Delta#eta_{gap}", "", n_bins_deltaEta_assoYield / 2 - 2,
      {-1. * bin_width_etaGap,
       (n_bins_deltaEta_assoYield / 2 - 3) * bin_width_etaGap});
  StrVar4Hist var_PtV2Jpsi("PtV2Jpsi", "p_{T}", "GeV/c", strAny_ptV2.fNbins,
                           {0., 1.});

  MHGroupTool1D DeltaPhi_sub(file_input,
                             "DeltaPhiUS_AssoYield_sub_int_EtaGap_%d_ptV2_%d",
                             {var_EtaGap, var_PtV2Jpsi}, {1, 1});
  MHGroupTool1D DeltaPhi_low(file_input,
                             "DeltaPhiUS_AssoYield_low_int_EtaGap_%d_ptV2_%d",
                             {var_EtaGap, var_PtV2Jpsi}, {1, 1});
  MHGroupTool1D DeltaPhi_high(file_input,
                              "DeltaPhiUS_AssoYield_high_int_EtaGap_%d_ptV2_%d",
                              {var_EtaGap, var_PtV2Jpsi}, {1, 1});

  uint n_bins_varPtV2Jpsi = var_PtV2Jpsi.fNbins;
  uint n_bins_varEtaGap = var_EtaGap.fNbins;

#define MH1DGetBin(...) GetHist(vector<int>{__VA_ARGS__})
  for (int iPtV2 = 1; iPtV2 <= n_bins_varPtV2Jpsi; iPtV2++) {
    for (int iEtaGap = 1; iEtaGap <= n_bins_varEtaGap; iEtaGap++) {
      auto h_deltaPhi_low = DeltaPhi_low.MH1DGetBin(iEtaGap, iPtV2);
    }
  }
}