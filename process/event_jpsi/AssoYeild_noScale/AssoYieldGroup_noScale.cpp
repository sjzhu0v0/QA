#include "MFit.h"
#include "MHelper.h"
#include "MMath.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"
#include "TSystem.h"
#include "yaml-cpp/yaml.h"

void AssoYieldGroup_noScale(
    TString path_input = "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
                         "AssoYield_24pass1/AssoYield_noScale.root",
    TString path_output = "/home/szhu/work/alice/analysis/QA/test/"
                          "AssoYieldGroupEtagap_NoScale.root",
    TString path_pdf = "/home/szhu/work/alice/analysis/QA/test/"
                       "AssoYieldGroup_noScale.pdf") {
  gROOT->SetBatch(kTRUE);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  YAML::Node config = YAML::LoadFile("config.yaml");
  const double low_edge_deltaPhiToPi =
      config["hist_binning"]["low_edge_deltaPhiToPi"].as<double>();
  const double up_edge_deltaPhiToPi =
      config["hist_binning"]["up_edge_deltaPhiToPi"].as<double>();
  int n_rebin_mass_assoYield =
      config["hist_binning"]["n_rebin_mass_assoYield"].as<int>();

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
  int n_bins_mass_assoYield =
      config["hist_binning"]["n_bins_mass_assoYield"].as<int>();
  StrVar4Hist var_MassJpsiCandidate("MassUS", "M_{ee}", "GeV^{2}/c^{4}",
                                    n_bins_mass_assoYield, {1.8, 5.4});
  StrVar4Hist var_PtJpsiCandidate("PtUS", "p_{T}", "GeV/c", 10, {0., 5.});
  int n_bins_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["n_bins"].as<int>();
  double min_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["min"].as<double>();
  double max_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["max"].as<double>();
  StrVar4Hist var_DeltaEtaUS("DeltaEtaUS", "#Delta#eta_{J/#psi, track}", "",
                             n_bins_deltaEta_assoYield,
                             {min_deltaEta_assoYield, max_deltaEta_assoYield});
  int n_bins_deltaPhi_assoYield =
      config["hist_binning"]["n_bins_deltaPhi_assoYield"].as<int>();
  StrVar4Hist var_DeltaPhiUS(
      "DeltaPhiUS", "#Delta#phi_{J/#psi, track}", "", n_bins_deltaPhi_assoYield,
      {low_edge_deltaPhiToPi * M_PI, up_edge_deltaPhiToPi * M_PI});
  StrVar4Hist var_EtaGap("EtaGap", "#Delta#eta_{gap}", "", 6, {-0.4, 2.});
  StrVar4Hist var_PtV2Jpsi("PtV2Jpsi", "p_{T}", "GeV/c", strAny_ptV2.fNbins,
                           {0., 1.});

  MIndexHist indexHistMass(var_MassJpsiCandidate, 1, n_rebin_mass_assoYield);
  MIndexHist indexHistPtJpsiCandidate(var_PtJpsiCandidate, 1, 1);
  MIndexHist indexHistDeltaPhiUS(var_DeltaPhiUS, 1, 1);
  int n_rebin_deltaEta_assoYield =
      config["hist_binning"]["n_rebin_deltaEta_assoYield"].as<int>();
  MIndexHist indexHistDeltaEtaUS(var_DeltaEtaUS, 1, n_rebin_deltaEta_assoYield);
  MIndexAny indexAnyPtV2Jpsi(strAny_ptV2, 1);

  gDirectory = nullptr;
  MHGroupTool2D *hgroupTool2d_total_mass = new MHGroupTool2D(
      file_input, "h2_total_mass_pt_%d_%d",
      {var_MassJpsiCandidate, var_PtV2Jpsi}, {n_rebin_mass_assoYield, 1});
  MHGroupTool2D *hgroupTool2d_lowMult_mass = new MHGroupTool2D(
      file_input, "h2_lowMult_mass_pt_%d_%d",
      {var_MassJpsiCandidate, var_PtV2Jpsi}, {n_rebin_mass_assoYield, 1});
  MHGroupTool2D *hgroupTool2d_highMult_mass = new MHGroupTool2D(
      file_input, "h2_highMult_mass_pt_%d_%d",
      {var_MassJpsiCandidate, var_PtV2Jpsi}, {n_rebin_mass_assoYield, 1});

  MHist1D h1_AssoYield_lowMult(indexHistMass, "AssoYield_lowMult");
  MHist1D h1_AssoYield_highMult(indexHistMass, "AssoYield_highMult");

  gDirectory = nullptr;

  auto h1Vec_AssoYield_lowMult =
      MakeMVec(h1_AssoYield_lowMult, indexHistDeltaEtaUS, indexHistDeltaPhiUS,
               indexAnyPtV2Jpsi);
  auto h1Vec_AssoYield_highMult =
      MakeMVec(h1_AssoYield_highMult, indexHistDeltaEtaUS, indexHistDeltaPhiUS,
               indexAnyPtV2Jpsi);

  gPublisherCanvas = new MPublisherCanvas(path_pdf, 2, 1, 600, 600);

  for (auto i_ptV2 : indexAnyPtV2Jpsi)
    for (auto i_mass : indexHistMass) {
      auto hist_high = hgroupTool2d_highMult_mass->GetHist({i_mass, i_ptV2});
      auto hist_low = hgroupTool2d_lowMult_mass->GetHist({i_mass, i_ptV2});
      gPublisherCanvas->Draw(hist_high)->Draw(hist_low);
      for (auto i_deltaEta : indexHistDeltaEtaUS)
        for (auto i_deltaPhi : indexHistDeltaPhiUS) {
          MDouble valHigh(hist_high->GetBinContent(i_deltaEta, i_deltaPhi),
                          hist_high->GetBinError(i_deltaEta, i_deltaPhi));
          MDouble valLow(hist_low->GetBinContent(i_deltaEta, i_deltaPhi),
                         hist_low->GetBinError(i_deltaEta, i_deltaPhi));
          h1Vec_AssoYield_lowMult.currentObject().SetBinInfo(valLow);
          h1Vec_AssoYield_highMult.currentObject().SetBinInfo(valHigh);
        }
    }

  file_output->cd();
  for (auto i_deltaEta : indexHistDeltaEtaUS)
    for (auto i_deltaPhi : indexHistDeltaPhiUS)
      for (auto i_ptV2 : indexAnyPtV2Jpsi) {
        h1Vec_AssoYield_lowMult.currentObject().fHisto->Write();
        h1Vec_AssoYield_highMult.currentObject().fHisto->Write();
      }

  gPublisherCanvas->finalize();
  file_output->Write();
  file_output->Close();
}

int main(int argc, char **argv) {
  if (argc < 1) {
    cerr << "Usage: ./AssoYieldGroup_noScale [path_input] [path_output] "
            "[path_pdf]"
         << endl;
    return 1;
  }
  if (argc == 1) {
    AssoYieldGroup_noScale();
  } else if (argc == 4) {
    AssoYieldGroup_noScale(TString(argv[1]), TString(argv[2]),
                           TString(argv[3]));
  } else {
    cerr << "Usage: ./AssoYieldGroup_noScale [path_input] [path_output] "
            "[path_pdf]"
         << endl;
    return 1;
  }
  return 0;
}