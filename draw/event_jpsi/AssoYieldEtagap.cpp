#include "MHist.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TFitResult.h"
#include "yaml-cpp/yaml.h"

#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMatrixDSym.h"
#include <cmath>
#include <stdexcept>
#include <utility>

std::pair<double, double>
compute_V2_and_error(const TFitResultPtr &result_sub,
                     const TFitResultPtr &result_low,
                     int idx_a0_sub = 0, // result_sub 中 a0 的参数索引
                     int idx_a2_sub = 2, // result_sub 中 a2 的参数索引
                     int idx_a0_low = 0) // result_low 中 a0 的参数索引
{
  if (!result_sub.Get() || !result_low.Get()) {
    throw std::runtime_error("compute_v2_and_error: empty TFitResultPtr.");
  }

  const double a0s = result_sub->Parameter(idx_a0_sub);
  const double sa0s = result_sub->ParError(idx_a0_sub);
  const double a2s = result_sub->Parameter(idx_a2_sub);
  const double sa2s = result_sub->ParError(idx_a2_sub);

  const double a0l = result_low->Parameter(idx_a0_low);
  const double sa0l = result_low->ParError(idx_a0_low);

  const TMatrixDSym cov_sub = result_sub->GetCovarianceMatrix();
  const double cov_a2s_a0s =
      cov_sub(idx_a2_sub, idx_a0_sub); // = (idx_a0_sub, idx_a2_sub)

  const double D = a0s + a0l;
  if (D == 0.0 || !std::isfinite(D)) {
    throw std::runtime_error(
        "compute_v2_and_error: D=a0_sub+a0_low is zero or not finite.");
  }

  const double v2 = a2s / D;

  // 线性误差传播（独立近似）：
  // Var(v2) = (sa2s^2)/D^2 + (a2s^2)(sa0s^2 + sa0l^2)/D^4 -
  // 2*a2s*Cov(a2s,a0s)/D^3
  double var_v2 =
      (sa2s * sa2s) / (D * D) +
      (a2s * a2s) * ((sa0s * sa0s) + (sa0l * sa0l)) / (D * D * D * D) -
      2.0 * a2s * cov_a2s_a0s / (D * D * D);

  // 数值稳健性：负的极小值视作 0（浮点误差可能造成）
  if (var_v2 < 0 && var_v2 > -1e-30)
    var_v2 = 0.0;

  const double sv2 = std::sqrt(var_v2);
  return {v2, sv2};
}

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
                            (double)n_bins_deltaEta_assoYield;
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

  // selected bins for ptV2
  // 10, 12, 14, 16, 18
  // 10, 12, 14, 33

  MHGroupTool1D DeltaPhi_sub(file_input,
                             "DeltaPhiUS_AssoYield_sub_int_EtaGap_%d_ptV2_%d",
                             {var_EtaGap, var_PtV2Jpsi}, {1, 1});
  MHGroupTool1D DeltaPhi_low(file_input,
                             "DeltaPhiUS_AssoYield_low_int_EtaGap_%d_ptV2_%d",
                             {var_EtaGap, var_PtV2Jpsi}, {1, 1});
  MHGroupTool1D DeltaPhi_high(file_input,
                              "DeltaPhiUS_AssoYield_high_int_EtaGap_%d_ptV2_%d",
                              {var_EtaGap, var_PtV2Jpsi}, {1, 1});

  gStyle->SetOptStat(0);
  TCanvas *c_assocYield =
      new TCanvas("c_assocYield", "c_assocYield", 1800, 600);
  c_assocYield->Divide(3, 1);
  auto assocYield_highMult_ptInt =
      (TH1D *)DeltaPhi_high.GetHist({1, 46})->Clone(
          "assocYield_highMult_ptInt");
  auto assocYield_lowMult_ptInt =
      (TH1D *)DeltaPhi_low.GetHist({1, 46})->Clone("assocYield_lowMult_ptInt");
  auto assocYield_sub_ptInt =
      (TH1D *)DeltaPhi_sub.GetHist({1, 46})->Clone("assocYield_sub_ptInt");
  MRootGraphic::StyleHistCommon(assocYield_highMult_ptInt);
  MRootGraphic::StyleHistCommon(assocYield_lowMult_ptInt);
  MRootGraphic::StyleHistCommon(assocYield_sub_ptInt);
  // MRootGraphic::StyleCommon();
  c_assocYield->cd(1);
  assocYield_highMult_ptInt->SetTitle("Associated yield in high mult.");
  assocYield_highMult_ptInt->GetYaxis()->SetTitle(
      "Associated yield per J/#psi");
  assocYield_highMult_ptInt->DrawClone();
  c_assocYield->cd(2);
  assocYield_lowMult_ptInt->SetTitle("Associated yield in low mult.");
  assocYield_lowMult_ptInt->GetYaxis()->SetTitle("Associated yield per J/#psi");
  assocYield_lowMult_ptInt->DrawClone();
  c_assocYield->cd(3);
  assocYield_sub_ptInt->SetTitle("Substracted associated yield");
  assocYield_sub_ptInt->GetYaxis()->SetTitle("Associated yield per J/#psi");
  assocYield_sub_ptInt->DrawClone();

  MIndexHist indexEtagap(var_EtaGap);
  MIndexHist indexPtV2Jpsi(var_PtV2Jpsi);

  gPublisherCanvas = new MPublisherCanvas(path_pdf, 3, 1);

#define MH1DGetBin(...) GetHist(vector<int>{__VA_ARGS__})

  auto v2REF_etaGap = MRootIO::GetObjectDiectly<TH1D>(
      "/home/szhu/work/alice/analysis/QA/output/event_track/"
      "v2_etagap_24pass1.root:EtaGap_v2");

  file_output->cd();

  auto V2_pT_etaGap =
      new TH2D("V2Jpsi_pT_etaGap",
               "V2Jpsi_pT_etaGap;p_{T} "
               "(GeV/c);#Delta#eta;V_{2}^{J/#psi}",
               4, vector<double>{0., 1., 2., 3., 5.}.data(),
               n_bins_deltaEta_assoYield / 2 - 2, -1. * bin_width_etaGap,
               (n_bins_deltaEta_assoYield / 2 - 3) * bin_width_etaGap);
  auto v2_pT_etaGap =
      new TH2D("v2Jpsi_pT_etaGap",
               "v2Jpsi_pT_etaGap;p_{T} "
               "(GeV/c);#Delta#eta;v_{2}^{J/#psi}",
               4, vector<double>{0., 1., 2., 3., 5.}.data(),
               n_bins_deltaEta_assoYield / 2 - 2, -1. * bin_width_etaGap,
               (n_bins_deltaEta_assoYield / 2 - 3) * bin_width_etaGap);
  gDirectory = nullptr;

  int iPt_v2_pT_etaGap = 0;
  for (auto iPt : {11, 13, 15, 34}) {
    iPt_v2_pT_etaGap++;
    for (auto iEtagpa : indexEtagap) {
      auto h_sub = DeltaPhi_sub.MH1DGetBin(iEtagpa, iPt);
      auto h_high = DeltaPhi_high.MH1DGetBin(iEtagpa, iPt);
      auto h_low = DeltaPhi_low.MH1DGetBin(iEtagpa, iPt);
      auto f_sub = (TF1 *)(h_sub->GetFunction("f1_modulation"));
      auto f_high = (TF1 *)(h_high->GetFunction("f1_modulation"));
      auto f_low = (TF1 *)(h_low->GetFunction("f1_modulation"));
      // gPublisherCanvas->Draw(f_sub)->Draw(f_low)->Draw(f_high);
      auto result_sub = h_sub->Fit(f_sub, "S Q N R");
      auto result_low = h_low->Fit(f_low, "S Q N R");
      auto result_high = h_high->Fit(f_high, "S Q N R");
      auto [val_V2, err_V2] = compute_V2_and_error(result_sub, result_low);
      V2_pT_etaGap->SetBinContent(iPt_v2_pT_etaGap, iEtagpa, val_V2);
      V2_pT_etaGap->SetBinError(iPt_v2_pT_etaGap, iEtagpa, err_V2);
      double val_v2REF = v2REF_etaGap->GetBinContent(iEtagpa);
      double err_v2REF = v2REF_etaGap->GetBinError(iEtagpa);

      double v2Jpsi = val_V2 / val_v2REF;
      double err_v2Jpsi =
          std::sqrt((err_V2 * err_V2) / (val_v2REF * val_v2REF) +
                    (val_V2 * val_V2) * (err_v2REF * err_v2REF) /
                        (val_v2REF * val_v2REF * val_v2REF * val_v2REF));
      v2_pT_etaGap->SetBinContent(iPt_v2_pT_etaGap, iEtagpa, v2Jpsi);
      v2_pT_etaGap->SetBinError(iPt_v2_pT_etaGap, iEtagpa, err_v2Jpsi);
    }
  }

  file_output->cd();
  for (int i = 1; i <= v2_pT_etaGap->GetNbinsX(); ++i) {
    auto h_proj = v2_pT_etaGap->ProjectionY(
        TString::Format("v2Jpsi_etaGap_pTbin%d", i), i, i);
    MRootGraphic::StyleHistCommon(h_proj);
    h_proj->GetYaxis()->SetTitle("v_{2}^{J/#psi}");
    auto H_proj = V2_pT_etaGap->ProjectionY(
        TString::Format("V2Jpsi_etaGap_pTbin%d", i), i, i);
    MRootGraphic::StyleHistCommon(H_proj);
    H_proj->GetYaxis()->SetTitle("V_{2}^{J/#psi}");
  }
  for (int i = 1; i <= v2_pT_etaGap->GetNbinsY(); ++i) {
    auto h_proj = v2_pT_etaGap->ProjectionX(
        TString::Format("v2Jpsi_pT_etaGapbin%d", i), i, i);
    MRootGraphic::StyleHistCommon(h_proj);
    h_proj->GetYaxis()->SetTitle("v_{2}^{J/#psi}");
    auto H_proj = V2_pT_etaGap->ProjectionX(
        TString::Format("V2Jpsi_pT_etaGapbin%d", i), i, i);
    MRootGraphic::StyleHistCommon(H_proj);
    H_proj->GetYaxis()->SetTitle("V_{2}^{J/#psi}");
  }

  TCanvas *c_v2_pT_etaGap =
      new TCanvas("c_v2_pT_etaGap", "c_v2_pT_etaGap", 1200, 1200);
  c_v2_pT_etaGap->Divide(2, 2);

  for (int i = 1; i <= 4; ++i) {
    c_v2_pT_etaGap->cd(i);
    auto h = (TH1D *)file_output->Get(
        TString::Format("v2Jpsi_pT_etaGapbin%d", i * 2));
    double etaGap = (-1. + (i * 2 - 1)) * bin_width_etaGap;
    h->SetTitle(TString::Format("v_{2}^{J/#psi} at #Delta#eta = %.2f", etaGap));
    h->Draw();
  }

  TCanvas *c_V2_pT_etaGap =
      new TCanvas("c_V2_pT_etaGap", "c_V2_pT_etaGap", 1200, 1200);
  c_V2_pT_etaGap->Divide(2, 2);
  for (int i = 1; i <= 4; ++i) {
    c_V2_pT_etaGap->cd(i);
    auto h = (TH1D *)file_output->Get(
        TString::Format("V2Jpsi_pT_etaGapbin%d", i * 2));
    double etaGap = (-1. + (i * 2 - 1)) * bin_width_etaGap;
    h->SetTitle(TString::Format("V_{2}^{J/#psi} at #Delta#eta = %.2f", etaGap));
    h->Draw();
  }

  // file_output->ls();

  gPublisherCanvas->finalize();
  file_output->Write();
  // file_output->Close();
}