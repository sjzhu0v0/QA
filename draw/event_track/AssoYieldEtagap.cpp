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

// 计算 v2 及误差（独立情形：Cov(a0_low, a2_sub)=Cov(a0_low, a0_sub)=0）
std::pair<double, double>
compute_v2_and_error(const TFitResultPtr &result_sub,
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
  double v2_final = sqrt(v2);
  double err_final = 0.5 * sv2 / sqrt(v2);
  return {v2_final, err_final};
}

void AssoYieldEtagap(
    TString path_input = "/home/szhu/work/alice/analysis/QA/input/event_track/"
                         "TrackAssocYeild_24apass1.root",
    TString path_output = "/home/szhu/work/alice/analysis/QA/output/"
                          "event_track/v2_etagap_24pass1.root",
    TString path_pdf = "/home/szhu/work/alice/analysis/QA/plot/event_track/"
                       "v2_etagap_24pass1.pdf") {
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;
  TFile *file_input = new TFile(path_input);
  TFile *file_output = new TFile(path_output, "RECREATE");
  YAML::Node config = YAML::LoadFile("config.yaml");
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

  MIndexHist index_etaGap(var_EtaGap);

  MHGroupTool1D Hs_highMult(file_input, "h2_highMult_EtaGap%d", {var_EtaGap},
                            {1});
  MHGroupTool1D Hs_lowMult(file_input, "h2_lowMult_EtaGap%d", {var_EtaGap},
                           {1});
  MHGroupTool1D Hs_sub(file_input, "h2_sub_EtaGap%d", {var_EtaGap}, {1});
  file_output->cd();
  MHist1D v2_etaGap(index_etaGap, "v2");
  v2_etaGap.fHisto->GetYaxis()->SetTitle("v_{2}^{REF}");
  gDirectory = nullptr;
  MRootGraphic::StyleCommon();
  gPublisherCanvas = new MPublisherCanvas(path_pdf, 2, 2);

  for (auto iEtagap : index_etaGap) {
    auto h_sub = Hs_sub.GetHist(vector<int>{iEtagap});
    auto h_high = Hs_highMult.GetHist(vector<int>{iEtagap});
    auto h_low = Hs_lowMult.GetHist(vector<int>{iEtagap});
    auto f_sub =
        (TF1 *)h_sub->GetListOfFunctions()->FindObject("f1_modulation");
    auto f_high =
        (TF1 *)h_high->GetListOfFunctions()->FindObject("f1_modulation");
    auto f_low =
        (TF1 *)h_low->GetListOfFunctions()->FindObject("f1_modulation");

    auto result_sub = h_sub->Fit(f_sub, "S Q N R");
    // The line 'auto result_sub = h_sub->Fit(f_sub, "S Q N R");' performs a fit
    // of the histogram 'h_sub' using the function 'f_sub'. The fit options
    // specified are: "S": This option indicates that the fit should return a
    // 'TFitResultPtr' object, which contains detailed information about the fit
    // results, including parameter values, errors, and fit quality. "Q": This
    // option suppresses the printing of fit results to the console. It is
    // useful when you want to perform the fit silently without displaying
    // output. "N": This option prevents the histogram from being drawn after
    // the fit. It is useful when you want to fit the data without visualizing
    // it immediately. "R": This option restricts the fit to the range defined
    // by the function 'f_sub'. It ensures that the fit is performed only within
    // the specified limits of the function.
    auto result_high = h_high->Fit(f_high, "S Q N R");
    auto result_low = h_low->Fit(f_low, "S Q N R");

    auto pair_v2 = compute_v2_and_error(result_sub, result_low);
    v2_etaGap.SetBinInfo(pair_v2.first, pair_v2.second);
    cout << "Etagap bin " << iEtagap << ": v2 = " << pair_v2.first << " +/- "
         << pair_v2.second << endl;
  }
  MRootGraphic::StyleHistCommon(v2_etaGap.fHisto.get());
  gPublisherCanvas->DrawClone(v2_etaGap.fHisto.get(), "Text");
  gPublisherCanvas->finalize();
  file_output->Write();
  file_output->Close();
}