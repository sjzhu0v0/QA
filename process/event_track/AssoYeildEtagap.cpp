#include "MFit.h"
#include "MHelper.h"
#include "MMath.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"
#include "yaml-cpp/yaml.h"

void AssoYieldEtagap(
    TString path_input = "/home/szhu/work/alice/analysis/QA/test/"
                         "AssoYieldGroupEtagap_NoScale.root",
    TString path_output = "/home/szhu/work/alice/analysis/QA/test/"
                          "AssoYieldFit_noScale.root") {
  gErrorIgnoreLevel = kWarning;
  TFile *file_input = new TFile(path_input);
  TFile *file_output = new TFile(path_output, "RECREATE");
  YAML::Node config = YAML::LoadFile("config.yaml");
  // int n_rebin_mass_assoYield =
  //     config["hist_binning"]["n_rebin_mass_assoYield"].as<int>();
  const double low_edge_deltaPhiToPi =
      config["hist_binning"]["low_edge_deltaPhiToPi"].as<double>();
  const double up_edge_deltaPhiToPi =
      config["hist_binning"]["up_edge_deltaPhiToPi"].as<double>();

  int n_bins_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["n_bins"].as<int>();
  double min_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["min"].as<double>();
  double max_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["max"].as<double>();
  double bin_width_etaGap = (max_deltaEta_assoYield - min_deltaEta_assoYield) /
                            (double)n_bins_deltaEta_assoYield;
  StrVar4Hist var_EtaGap(
      "EtaGap", "#Delta#eta_{gap}", "", n_bins_deltaEta_assoYield / 2 - 2,
      {-1. * bin_width_etaGap,
       (n_bins_deltaEta_assoYield / 2 - 3) * bin_width_etaGap});
  MIndexHist indexHistEtaGap(var_EtaGap, 1, 1);

  auto h2_total = (TH2D *)file_input->Get("h2_total");
  auto h2_lowMult = (TH2D *)file_input->Get("h2_lowMult");
  auto h2_highMult = (TH2D *)file_input->Get("h2_highMult");
  auto h2_sub = (TH2D *)file_input->Get("h2_sub");

  auto hist_etaGap = [](double eta_gap, TH2D *h2) {
    h2->GetXaxis()->SetRangeUser(-1.8, -abs(eta_gap));
    auto h1_1 = h2->ProjectionY(Form("%s_negEtaGap%.2f", h2->GetName()));
    h2->GetXaxis()->SetRangeUser(abs(eta_gap), 1.8);
    auto h1_2 = h2->ProjectionY(Form("%s_posEtaGap%.2f", h2->GetName()));
    h1_1->Add(h1_2);
    h1_1->SetTitle(
        Form("%s with #Delta#eta_{gap}=%.2f", h2->GetTitle(), eta_gap));
    h1_2->Delete();
    return h1_1;
  };
  file_output->cd();

  gPublisherCanvas = new MPublisherCanvas("test.pdf", 2, 2);

#define ModulationMultClass(mult_class)                                        \
  MHist1D a0_##mult_class(indexHistEtaGap, "a0_" #mult_class);                 \
  MHist1D a1_##mult_class(indexHistEtaGap, "a1_" #mult_class);                 \
  MHist1D a2_##mult_class(indexHistEtaGap, "a2_" #mult_class);                 \
  MHist1D a3_##mult_class(indexHistEtaGap, "a3_" #mult_class);                 \
  for (auto iEtaGap : indexHistEtaGap) {                                       \
    double val_etaGap = var_EtaGap.GetBinUpperEdge(iEtaGap - 1);               \
    cout << "Fitting eta gap: " << val_etaGap << endl;                         \
    h2_##mult_class->GetXaxis()->SetRangeUser(-1.8, -abs(val_etaGap));         \
    auto h1_##mult_class = h2_##mult_class->ProjectionY(                       \
        Form("%s_EtaGap%d", h2_##mult_class->GetName(), iEtaGap));             \
    gPublisherCanvas->DrawClone(h1_##mult_class);                              \
    h2_##mult_class->GetXaxis()->SetRangeUser(abs(val_etaGap), 1.8);           \
    auto h1_temp = h2_##mult_class->ProjectionY(                               \
        Form("%s_posGap%d", h2_##mult_class->GetName(), iEtaGap));             \
    gPublisherCanvas->DrawClone(h1_temp);                                      \
    h1_##mult_class->Add(h1_temp);                                             \
    gPublisherCanvas->DrawClone(h1_##mult_class);                              \
    h1_temp->Delete();                                                         \
                                                                               \
    h2_##mult_class->GetXaxis()->SetRangeUser(-1.8, 1.8);                      \
    h1_##mult_class->SetName(                                                  \
        Form("%s_EtaGap%d", h2_##mult_class->GetName(), iEtaGap));             \
    h1_##mult_class->Write();                                                  \
    TF1 func_modulation(                                                       \
        "f1_modulation", "[a0]+2*([a1]*cos(x)+[a2]*cos(2*x)+[a3]*cos(3*x))",   \
        low_edge_deltaPhiToPi * M_PI, up_edge_deltaPhiToPi * M_PI);            \
    gPublisherCanvas->NewPad()->cd();                                          \
    h1_##mult_class->Fit(&func_modulation, "QR", "",                           \
                         low_edge_deltaPhiToPi * M_PI,                         \
                         up_edge_deltaPhiToPi * M_PI);                         \
    a0_##mult_class.SetBinInfo(func_modulation.GetParameter(0),                \
                               func_modulation.GetParError(0));                \
    a1_##mult_class.SetBinInfo(func_modulation.GetParameter(1),                \
                               func_modulation.GetParError(1));                \
    a2_##mult_class.SetBinInfo(func_modulation.GetParameter(2),                \
                               func_modulation.GetParError(2));                \
    a3_##mult_class.SetBinInfo(func_modulation.GetParameter(3),                \
                               func_modulation.GetParError(3));                \
  }                                                                            \
  a0_##mult_class.Write();                                                     \
  a1_##mult_class.Write();                                                     \
  a2_##mult_class.Write();                                                     \
  a3_##mult_class.Write();

  ModulationMultClass(highMult);
  ModulationMultClass(lowMult);
  ModulationMultClass(sub);
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

  AssoYieldEtagap(path_input, path_output);

  return 0;
}