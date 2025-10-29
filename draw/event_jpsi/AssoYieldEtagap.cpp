#include "MHist.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "yaml-cpp/yaml.h"

void AssoYieldEtagap(
    TString path_input = "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
                         "JpsiAssocYield_24apass1.root",
    TString path_output = "/home/szhu/work/alice/analysis/QA/output/event_jpsi/"
                          "AssoYieldEtagap.root",
    TString path_pdf = "/home/szhu/work/alice/analysis/QA/plot/event_jpsi/"
                       "AssoYieldEtagap.pdf") {
  gErrorIgnoreLevel = kWarning;
  TFile *file_input = new TFile(path_input);
  TFile *file_output = new TFile(path_output, "RECREATE");
  YAML::Node config = YAML::LoadFile("config.yaml");
  const double low_edge_deltaPhiToPi =
      config["hist_binning"]["low_edge_deltaPhiToPi"].as<double>();
  const double up_edge_deltaPhiToPi =
      config["hist_binning"]["up_edge_deltaPhiToPi"].as<double>();

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

  MIndexHist indexEtagap(var_EtaGap);
  MIndexHist indexPtV2Jpsi(var_PtV2Jpsi);
  file_output->cd();
  MHist2D a0_high(indexEtagap, indexPtV2Jpsi, "a0_high");
  MHist2D a0_low(indexEtagap, indexPtV2Jpsi, "a0_low");
  MHist2D a0_sub(indexEtagap, indexPtV2Jpsi, "a0_sub");
  MHist2D a1_high(indexEtagap, indexPtV2Jpsi, "a1_high");
  MHist2D a1_low(indexEtagap, indexPtV2Jpsi, "a1_low");
  MHist2D a1_sub(indexEtagap, indexPtV2Jpsi, "a1_sub");
  MHist2D a2_high(indexEtagap, indexPtV2Jpsi, "a2_high");
  MHist2D a2_low(indexEtagap, indexPtV2Jpsi, "a2_low");
  MHist2D a2_sub(indexEtagap, indexPtV2Jpsi, "a2_sub");
  MHist2D a3_high(indexEtagap, indexPtV2Jpsi, "a3_high");
  MHist2D a3_low(indexEtagap, indexPtV2Jpsi, "a3_low");
  MHist2D a3_sub(indexEtagap, indexPtV2Jpsi, "a3_sub");
  MHist2D V2(indexEtagap, indexPtV2Jpsi, "V2");

  gDirectory = file_output;
  gPublisherCanvas = new MPublisherCanvas(path_pdf, 3, 1);

#define MH1DGetBin(...) GetHist(vector<int>{__VA_ARGS__})
  for (auto iPt : indexPtV2Jpsi)
    for (auto iEtagpa : indexEtagap) {
      auto f_sub = (TF1 *)(DeltaPhi_sub.MH1DGetBin(iEtagpa, iPt)
                               ->GetFunction("f1_modulation"));
      auto f_high = (TF1 *)(DeltaPhi_high.MH1DGetBin(iEtagpa, iPt)
                                ->GetFunction("f1_modulation"));
      auto f_low = (TF1 *)(DeltaPhi_low.MH1DGetBin(iEtagpa, iPt)
                               ->GetFunction("f1_modulation"));
      gPublisherCanvas->Draw(f_sub)->Draw(f_low)->Draw(f_high);

      MDouble val_a0_sub(f_sub->GetParameter(0), f_sub->GetParError(0));
      a0_sub.SetBinInfo(val_a0_sub);
      MDouble val_a1_sub(f_sub->GetParameter(1), f_sub->GetParError(1));
      a1_sub.SetBinInfo(val_a1_sub);
      MDouble val_a2_sub(f_sub->GetParameter(2), f_sub->GetParError(2));
      a2_sub.SetBinInfo(val_a2_sub);
      MDouble val_a3_sub(f_sub->GetParameter(3), f_sub->GetParError(3));
      a3_sub.SetBinInfo(val_a3_sub);
      MDouble val_a0_low(f_low->GetParameter(0), f_low->GetParError(0));
      a0_low.SetBinInfo(val_a0_low);
      MDouble val_a1_low(f_low->GetParameter(1), f_low->GetParError(1));
      a1_low.SetBinInfo(val_a1_low);
      MDouble val_a2_low(f_low->GetParameter(2), f_low->GetParError(2));
      a2_low.SetBinInfo(val_a2_low);
      MDouble val_a3_low(f_low->GetParameter(3), f_low->GetParError(3));
      a3_low.SetBinInfo(val_a3_low);
      MDouble val_a0_high(f_high->GetParameter(0), f_high->GetParError(0));
      a0_high.SetBinInfo(val_a0_high);
      MDouble val_a1_high(f_high->GetParameter(1), f_high->GetParError(1));
      a1_high.SetBinInfo(val_a1_high);
      MDouble val_a2_high(f_high->GetParameter(2), f_high->GetParError(2));
      a2_high.SetBinInfo(val_a2_high);
      MDouble val_a3_high(f_high->GetParameter(3), f_high->GetParError(3));
      a3_high.SetBinInfo(val_a3_high);
      MDouble val_v2 = (val_a2_high - val_a2_low) / val_a0_high;
      V2.SetBinInfo(val_v2);
    }

  gPublisherCanvas->finalize();
  file_output->Write();
  file_output->Close();
}