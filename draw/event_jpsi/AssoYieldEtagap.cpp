#include "MHist.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TFitResult.h"
#include "TVectorD.h"
#include "yaml-cpp/yaml.h"

#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMatrixDSym.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

std::pair<double, double> EvalError(TF1* f, const TFitResultPtr& r, double x) {
  if (!f || !r.Get()) {
    return {0.0, 0.0};
  }

  // 1. 函数值
  double value = f->Eval(x);

  // 2. 协方差矩阵
  TMatrixDSym cov = r->GetCovarianceMatrix();
  int npar = f->GetNpar();

  // 3. 计算对参数的偏导数 ∂f/∂p_i
  std::vector<double> grad(npar);
  f->GradientPar(x, grad.data()); // ✔ TF1*（非 const）才能调用

  // 4. 误差传播
  double err2 = 0;
  for (int i = 0; i < npar; i++) {
    for (int j = 0; j < npar; j++) {
      err2 += grad[i] * cov(i, j) * grad[j];
    }
  }

  return {value, std::sqrt(std::max(err2, 0.0))};
}

static double PropagateWithCov(const TMatrixDSym& Cov, const TVectorD& grad) {
  // var = grad^T * Cov * grad
  double var = 0.0;
  const int n = grad.GetNrows();
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      var += grad(i) * Cov(i, j) * grad(j);
    }
  }
  return (var >= 0.0) ? std::sqrt(var) : 0.0; // 数值误差下 var 可能出现极小负数
}

std::vector<int> parseBinningString(const TString& str_binning) {
  std::vector<int> result;
  std::string s = str_binning.Data();
  std::stringstream ss(s);
  std::string item;

  while (std::getline(ss, item, ',')) {
    try {
      int value = std::stoi(item);
      result.push_back(value);
    } catch (const std::invalid_argument& e) {
      std::cerr << "Invalid integer in binning string: " << item << std::endl;
    } catch (const std::out_of_range& e) {
      std::cerr << "Integer out of range in binning string: " << item << std::endl;
    }
  }

  return result;
}

void ComputePrimeErrorsWithFullCov(const TFitResultPtr& result_high,
                                   const TFitResultPtr& result_low, double& a0_prime,
                                   double& ea0_prime, double& a2_prime, double& ea2_prime,
                                   double& ratio, double& eratio) {
  if (!result_high.Get() || !result_low.Get()) {
    cout << "Fit result is null! Cannot compute errors." << endl;
    a0_prime = ea0_prime = a2_prime = ea2_prime = ratio = eratio = NAN;
    return;
  }

  // ---- 取参数（假设: par0=a0, par1=a1, par2=a2）----
  const double a0_h = result_high->Parameter(0);
  const double a1_h = result_high->Parameter(1);
  const double a2_h = result_high->Parameter(2);

  const double a0_l = result_low->Parameter(0);
  const double a1_l = result_low->Parameter(1);
  const double a2_l = result_low->Parameter(2);

  // ---- 取协方差矩阵（每个 fit 内部 3x3，含相关性）----
  // ROOT: GetCovarianceMatrix() 返回 TMatrixDSym
  const TMatrixDSym CovH = result_high->GetCovarianceMatrix();
  const TMatrixDSym CovL = result_low->GetCovarianceMatrix();

  // ---- 拼成总协方差矩阵 6x6： (a0_h,a1_h,a2_h,a0_l,a1_l,a2_l) ----
  TMatrixDSym Cov(6);
  Cov.Zero();

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      Cov(i, j) = CovH(i, j);
      Cov(i + 3, j + 3) = CovL(i, j);
    }
  }

  // ---- 定义 a0' 与 a2' ----
  a0_prime = a0_h - a1_h * a0_l / a1_l;
  a2_prime = a2_h - a1_h * a2_l / a1_l;

  // ---- 梯度（对 6 个参数）----
  // a0' = a0_h - a1_h*a0_l/a1_l
  // grad_a0' = (∂/∂a0_h, ∂/∂a1_h, ∂/∂a2_h, ∂/∂a0_l, ∂/∂a1_l, ∂/∂a2_l)
  TVectorD g_a0p(6);
  g_a0p = 0.0;
  g_a0p(0) = 1.0;
  g_a0p(1) = -a0_l / a1_l;
  g_a0p(2) = 0.0;
  g_a0p(3) = -a1_h / a1_l;
  g_a0p(4) = a1_h * a0_l / (a1_l * a1_l);
  g_a0p(5) = 0.0;

  // a2' = a2_h - a1_h*a2_l/a1_l
  TVectorD g_a2p(6);
  g_a2p = 0.0;
  g_a2p(0) = 0.0;
  g_a2p(1) = -a2_l / a1_l;
  g_a2p(2) = 1.0;
  g_a2p(3) = 0.0;
  g_a2p(4) = a1_h * a2_l / (a1_l * a1_l);
  g_a2p(5) = -a1_h / a1_l;

  // ---- 误差（含协方差）----
  ea0_prime = PropagateWithCov(Cov, g_a0p);
  ea2_prime = PropagateWithCov(Cov, g_a2p);

  // ---- ratio = a2'/a0' ----
  ratio = a2_prime / a0_prime;

  // 用链式法则： r = a2'/a0'
  // ∂r/∂p = (1/a0') * ∂a2'/∂p  - (a2'/a0'^2) * ∂a0'/∂p
  TVectorD g_r(6);
  const double inv_a0p = 1.0 / a0_prime;
  const double coeff = a2_prime / (a0_prime * a0_prime);
  for (int k = 0; k < 6; ++k) {
    g_r(k) = inv_a0p * g_a2p(k) - coeff * g_a0p(k);
  }

  eratio = PropagateWithCov(Cov, g_r);

  // ---- 输出 ----
  cout << "a0' = " << a0_prime << " +/- " << ea0_prime << endl;
  cout << "a2' = " << a2_prime << " +/- " << ea2_prime << endl;
  cout << "a2'/a0' = " << ratio << " +/- " << eratio << endl;
}

std::pair<double, double> compute_V2_and_error(const TFitResultPtr& result_sub,
                                               const TFitResultPtr& result_low, TF1* f_low,
                                               int idx_a0_sub = 0, // result_sub 中 a0 的参数索引
                                               int idx_a2_sub = 2 // result_sub 中 a2 的参数索引
) {
  if (!result_sub.Get() || !result_low.Get()) {
    throw std::runtime_error("compute_v2_and_error: empty TFitResultPtr.");
  }

  const double a0s = result_sub->Parameter(idx_a0_sub);
  const double sa0s = result_sub->ParError(idx_a0_sub);
  const double a2s = result_sub->Parameter(idx_a2_sub);
  const double sa2s = result_sub->ParError(idx_a2_sub);

  auto [a0l, sa0l] = EvalError(f_low, result_low, M_PI / 2.0);

  const TMatrixDSym cov_sub = result_sub->GetCovarianceMatrix();
  const double cov_a2s_a0s = cov_sub(idx_a2_sub, idx_a0_sub); // = (idx_a0_sub, idx_a2_sub)

  const double D = a0s + a0l;
  if (D == 0.0 || !std::isfinite(D)) {
    throw std::runtime_error("compute_v2_and_error: D=a0_sub+a0_low is zero or not finite.");
  }

  const double v2 = a2s / D;

  // 线性误差传播（独立近似）：
  // Var(v2) = (sa2s^2)/D^2 + (a2s^2)(sa0s^2 + sa0l^2)/D^4 -
  // 2*a2s*Cov(a2s,a0s)/D^3
  double var_v2 = (sa2s * sa2s) / (D * D) +
                  (a2s * a2s) * ((sa0s * sa0s) + (sa0l * sa0l)) / (D * D * D * D) -
                  2.0 * a2s * cov_a2s_a0s / (D * D * D);

  // 数值稳健性：负的极小值视作 0（浮点误差可能造成）
  if (var_v2 < 0 && var_v2 > -1e-30)
    var_v2 = 0.0;

  const double sv2 = std::sqrt(var_v2);
  return {v2, sv2};
}

void AssoYieldEtagap(TString path_input = "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
                                          "JpsiAssocYield_24apass1.root",
                     TString path_hist_ref = "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
                                             "JpsiAssocYield_24apass1.root:hist",
                     TString path_hist_ref_new =
                         "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
                         "JpsiAssocYield_24apass1.root:hist_new",
                     TString str_binning = "11,13,15,34",
                     TString path_output = "/home/szhu/work/alice/analysis/QA/output/event_jpsi/"
                                           "AssoYieldEtagap.root",
                     TString path_pdf = "/home/szhu/work/alice/analysis/QA/plot/event_jpsi/"
                                        "AssoYieldEtagap") {
  gErrorIgnoreLevel = kWarning;
  TFile* file_input = new TFile(path_input);
  TFile* file_output = new TFile(path_output, "RECREATE");
  YAML::Node config = YAML::LoadFile("config.yaml");
  const double low_edge_deltaPhiToPi = config["hist_binning"]["low_edge_deltaPhiToPi"].as<double>();
  const double up_edge_deltaPhiToPi = config["hist_binning"]["up_edge_deltaPhiToPi"].as<double>();

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
  StrVar4Hist var_NumContribCalibBinned("NumContribCalibUS", "N_{vtx contrib} Calibrated", "", 10,
                                        {0, 7, 12, 17, 22, 29, 37, 46, 57, 73, 300});
  StrVar4Hist var_MassJpsiCandidate("MassUS", "M_{ee}", "GeV^{2}/c^{4}", 90, {1.8, 5.4});
  StrVar4Hist var_PtJpsiCandidate("PtUS", "p_{T}", "GeV/c", 10, {0., 10.});
  int n_bins_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["n_bins"].as<int>();
  double min_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["min"].as<double>();
  double max_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["max"].as<double>();
  double bin_width_etaGap =
      (max_deltaEta_assoYield - min_deltaEta_assoYield) / (double)n_bins_deltaEta_assoYield;
  StrVar4Hist var_DeltaEtaUS("DeltaEtaUS", "#Delta#eta_{J/#psi, track}", "",
                             n_bins_deltaEta_assoYield,
                             {min_deltaEta_assoYield, max_deltaEta_assoYield});
  int n_bins_deltaPhi_assoYield = config["hist_binning"]["n_bins_deltaPhi_assoYield"].as<int>();
  StrVar4Hist var_DeltaPhiUS("DeltaPhiUS", "#Delta#phi_{J/#psi, track}", "",
                             n_bins_deltaPhi_assoYield,
                             {low_edge_deltaPhiToPi * M_PI, up_edge_deltaPhiToPi * M_PI});
  StrVar4Hist var_EtaGap(
      "EtaGap", "#Delta#eta_{gap}", "", n_bins_deltaEta_assoYield / 2 - 2,
      {-1. * bin_width_etaGap, (n_bins_deltaEta_assoYield / 2 - 3) * bin_width_etaGap});
  StrVar4Hist var_PtV2Jpsi("PtV2Jpsi", "p_{T}", "GeV/c", strAny_ptV2.fNbins, {0., 1.});

  // selected bins for ptV2
  // 10, 12, 14, 16, 18
  // 10, 12, 14, 33

  MHGroupTool1D DeltaPhi_sub(file_input, "DeltaPhiUS_AssoYield_sub_int_EtaGap_%d_ptV2_%d",
                             {var_EtaGap, var_PtV2Jpsi}, {1, 1});
  MHGroupTool1D DeltaPhi_low(file_input, "DeltaPhiUS_AssoYield_low_int_EtaGap_%d_ptV2_%d",
                             {var_EtaGap, var_PtV2Jpsi}, {1, 1});
  MHGroupTool1D DeltaPhi_high(file_input, "DeltaPhiUS_AssoYield_high_int_EtaGap_%d_ptV2_%d",
                              {var_EtaGap, var_PtV2Jpsi}, {1, 1});

  gStyle->SetOptStat(0);
  TCanvas* c_assocYield = new TCanvas("c_assocYield", "c_assocYield", 1800, 600);
  c_assocYield->Divide(3, 1);
  auto assocYield_highMult_ptInt =
      (TH1D*)DeltaPhi_high.GetHist({1, 46})->Clone("assocYield_highMult_ptInt");
  auto assocYield_lowMult_ptInt =
      (TH1D*)DeltaPhi_low.GetHist({1, 46})->Clone("assocYield_lowMult_ptInt");
  auto assocYield_sub_ptInt = (TH1D*)DeltaPhi_sub.GetHist({1, 46})->Clone("assocYield_sub_ptInt");
  MRootGraphic::StyleHistCommon(assocYield_highMult_ptInt);
  MRootGraphic::StyleHistCommon(assocYield_lowMult_ptInt);
  MRootGraphic::StyleHistCommon(assocYield_sub_ptInt);
  // MRootGraphic::StyleCommon();
  c_assocYield->cd(1);
  assocYield_highMult_ptInt->SetTitle("Associated yield in high mult.");
  assocYield_highMult_ptInt->GetYaxis()->SetTitle("Associated yield per J/#psi");
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

  gPublisherCanvas = new MPublisherCanvas(path_pdf + ".pdf", 3, 1);

#define MH1DGetBin(...) GetHist(vector<int>{__VA_ARGS__})

  auto v2REF_etaGap = MRootIO::GetObjectDiectly<TH1D>(path_hist_ref);
  auto v2REF_etaGap_new = MRootIO::GetObjectDiectly<TH1D>(path_hist_ref_new);

  file_output->cd();

  vector<int> vec_binning_sel = parseBinningString(str_binning);
  vector<double> vec_binning_pt;

  cout << "vec_binning_sel: ";
  for (auto i : vec_binning_sel) {
    cout << i << ": ";
    for (auto j : strAny_ptV2.bins[i - 1])
      cout << j << ", ";
    cout << ";";
  }
  cout << endl;

  int ibin_last = 0;
  for (auto bin : vec_binning_sel) {
    vector<int> bin_sel = strAny_ptV2.bins[bin - 1];
    if (ibin_last != 0 && ibin_last != bin_sel[0] - 1) {
      cerr << "str_binning is wrong. exit" << endl;
      exit(1);
    }
    vec_binning_pt.push_back(bin_sel[0]);
    ibin_last = bin_sel[bin_sel.size() - 1];
  }
  cout << "raw vec_binning_pt: ";
  for (auto i : vec_binning_pt)
    cout << i << ", ";
  cout << endl;

  for (int i = 0; i < vec_binning_pt.size(); i++) {
    vec_binning_pt[i] = vec_binning_pt[i] * 0.5 - 0.5;
  }
  vec_binning_pt.push_back(ibin_last * 0.5);
  cout << "vec_binning_pt: ";
  for (auto i : vec_binning_pt)
    cout << i << ", ";
  cout << endl;

  auto V2_pT_etaGap =
      new TH2D("V2Jpsi_pT_etaGap",
               "V2Jpsi_pT_etaGap;p_{T} "
               "(GeV/c);#Delta#eta;V_{2}^{J/#psi}",
               vec_binning_pt.size() - 1, vec_binning_pt.data(), n_bins_deltaEta_assoYield / 2 - 2,
               -1. * bin_width_etaGap, (n_bins_deltaEta_assoYield / 2 - 3) * bin_width_etaGap);
  auto v2_pT_etaGap =
      new TH2D("v2Jpsi_pT_etaGap",
               "v2Jpsi_pT_etaGap;p_{T} "
               "(GeV/c);#Delta#eta;v_{2}^{J/#psi}",
               vec_binning_pt.size() - 1, vec_binning_pt.data(), n_bins_deltaEta_assoYield / 2 - 2,
               -1. * bin_width_etaGap, (n_bins_deltaEta_assoYield / 2 - 3) * bin_width_etaGap);

  auto V2_pT_etaGap_new =
      new TH2D("V2Jpsi_pT_etaGap_new",
               "V2Jpsi_pT_etaGap;p_{T} "
               "(GeV/c);#Delta#eta;V_{2}^{J/#psi}",
               vec_binning_pt.size() - 1, vec_binning_pt.data(), n_bins_deltaEta_assoYield / 2 - 2,
               -1. * bin_width_etaGap, (n_bins_deltaEta_assoYield / 2 - 3) * bin_width_etaGap);
  auto v2_pT_etaGap_new =
      new TH2D("v2Jpsi_pT_etaGap_new",
               "v2Jpsi_pT_etaGap;p_{T} "
               "(GeV/c);#Delta#eta;v_{2}^{J/#psi}",
               vec_binning_pt.size() - 1, vec_binning_pt.data(), n_bins_deltaEta_assoYield / 2 - 2,
               -1. * bin_width_etaGap, (n_bins_deltaEta_assoYield / 2 - 3) * bin_width_etaGap);
  gDirectory = nullptr;

  int iPt_v2_pT_etaGap = 0;
  for (auto iPt : vec_binning_sel) {
    iPt_v2_pT_etaGap++;
    for (auto iEtagap : indexEtagap) {
      auto h_sub = DeltaPhi_sub.MH1DGetBin(iEtagap, iPt);
      auto h_high = DeltaPhi_high.MH1DGetBin(iEtagap, iPt);
      auto h_low = DeltaPhi_low.MH1DGetBin(iEtagap, iPt);
      auto f_sub = (TF1*)(h_sub->GetFunction("f1_modulation"));
      auto f_high = (TF1*)(h_high->GetFunction("f1_modulation"));
      auto f_low = (TF1*)(h_low->GetFunction("f1_modulation"));
      // gPublisherCanvas->Draw(f_sub)->Draw(f_low)->Draw(f_high);
      auto result_sub = h_sub->Fit(f_sub, "S Q N R");
      auto result_low = h_low->Fit(f_low, "S Q N R");
      auto result_high = h_high->Fit(f_high, "S Q N R");
      auto [val_V2, err_V2] = compute_V2_and_error(result_sub, result_low, f_low);
      V2_pT_etaGap->SetBinContent(iPt_v2_pT_etaGap, iEtagap, val_V2);
      V2_pT_etaGap->SetBinError(iPt_v2_pT_etaGap, iEtagap, err_V2);
      double val_v2REF = v2REF_etaGap->GetBinContent(iEtagap);
      double err_v2REF = v2REF_etaGap->GetBinError(iEtagap);

      double v2Jpsi = val_V2 / val_v2REF;
      double err_v2Jpsi = std::sqrt((err_V2 * err_V2) / (val_v2REF * val_v2REF) +
                                    (val_V2 * val_V2) * (err_v2REF * err_v2REF) /
                                        (val_v2REF * val_v2REF * val_v2REF * val_v2REF));
      v2_pT_etaGap->SetBinContent(iPt_v2_pT_etaGap, iEtagap, v2Jpsi);
      v2_pT_etaGap->SetBinError(iPt_v2_pT_etaGap, iEtagap, err_v2Jpsi);

      double a0prime, a2prime, ratio;
      double err_a0prime, err_a2prime, err_ratio;
      ComputePrimeErrorsWithFullCov(result_high, result_low, a0prime, err_a0prime, a2prime,
                                    err_a2prime, ratio, err_ratio);

      V2_pT_etaGap_new->SetBinContent(iPt_v2_pT_etaGap, iEtagap, ratio);
      V2_pT_etaGap_new->SetBinError(iPt_v2_pT_etaGap, iEtagap, err_ratio);

      double v2_ref_new, err_v2_ref_new;
      v2_ref_new = v2REF_etaGap_new->GetBinContent(iEtagap);
      err_v2_ref_new = v2REF_etaGap_new->GetBinError(iEtagap);

      double v2_poi_new;
      double err_v2_poi_new;

      if (v2_ref_new == 0) {
        v2_poi_new = 0;
        err_v2_poi_new = 0;
      } else {
        v2_poi_new = ratio / v2_ref_new;
        err_v2_poi_new =
            err_ratio * err_ratio / v2_ref_new / v2_ref_new +
            (err_v2_ref_new * ratio / v2_ref_new) * (err_v2_ref_new * ratio / v2_ref_new);
        err_v2_poi_new = sqrt(err_v2_poi_new);
      }

      v2_pT_etaGap_new->SetBinContent(iPt_v2_pT_etaGap, iEtagap, v2_poi_new);
      v2_pT_etaGap_new->SetBinError(iPt_v2_pT_etaGap, iEtagap, err_v2_poi_new);
    }
  }

  file_output->cd();
  for (int i = 1; i <= v2_pT_etaGap->GetNbinsX(); ++i) {
    auto h_proj = v2_pT_etaGap->ProjectionY(TString::Format("v2Jpsi_etaGap_pTbin%d", i), i, i);
    MRootGraphic::StyleHistCommon(h_proj);
    h_proj->GetYaxis()->SetTitle("v_{2}^{J/#psi}");
    auto H_proj = V2_pT_etaGap->ProjectionY(TString::Format("V2Jpsi_etaGap_pTbin%d", i), i, i);
    MRootGraphic::StyleHistCommon(H_proj);
    H_proj->GetYaxis()->SetTitle("V_{2}^{J/#psi}");
  }
  for (int i = 1; i <= v2_pT_etaGap->GetNbinsY(); ++i) {
    auto h_proj = v2_pT_etaGap->ProjectionX(TString::Format("v2Jpsi_pT_etaGapbin%d", i), i, i);
    MRootGraphic::StyleHistCommon(h_proj);
    h_proj->GetYaxis()->SetTitle("v_{2}^{J/#psi}");
    auto H_proj = V2_pT_etaGap->ProjectionX(TString::Format("V2Jpsi_pT_etaGapbin%d", i), i, i);
    MRootGraphic::StyleHistCommon(H_proj);
    H_proj->GetYaxis()->SetTitle("V_{2}^{J/#psi}");
  }

  TCanvas* c_v2_pT_etaGap = new TCanvas("c_v2_pT_etaGap", "c_v2_pT_etaGap", 1200, 1200);
  c_v2_pT_etaGap->Divide(2, 2);

  for (int i = 1; i <= 4; ++i) {
    c_v2_pT_etaGap->cd(i);
    auto h = (TH1D*)file_output->Get(TString::Format("v2Jpsi_pT_etaGapbin%d", i * 2));
    double etaGap = (-1. + (i * 2 - 1)) * bin_width_etaGap;
    h->SetTitle(TString::Format("v_{2}^{J/#psi} at #Delta#eta = %.2f", etaGap));
    h->Draw();
  }
  c_v2_pT_etaGap->SaveAs(path_pdf + "_v2Jpsi_etaGap.pdf");

  TCanvas* c_v2_pT_etaGap1p2 = new TCanvas("c_v2_pT_etaGap1p2", "c_v2_pT_etaGap1p2", 600, 600);
  c_v2_pT_etaGap1p2->cd();
  auto h_v2_etaGap1p2 = (TH1D*)file_output->Get("v2Jpsi_pT_etaGapbin8");
  // h_v2_etaGap1p2->GetYaxis()->SetRangeUser();
  h_v2_etaGap1p2->SetTitle("v_{2}^{J/#psi} at #Delta#eta = 1.2");
  h_v2_etaGap1p2->Draw();
  c_v2_pT_etaGap1p2->SaveAs(path_pdf + "_v2_pT_etaGap1p2.pdf");

  TCanvas* c_V2_pT_etaGap = new TCanvas("c_V2_pT_etaGap", "c_V2_pT_etaGap", 1200, 1200);
  c_V2_pT_etaGap->Divide(2, 2);
  for (int i = 1; i <= 4; ++i) {
    c_V2_pT_etaGap->cd(i);
    auto h = (TH1D*)file_output->Get(TString::Format("V2Jpsi_pT_etaGapbin%d", i * 2));
    double etaGap = (-1. + (i * 2 - 1)) * bin_width_etaGap;
    h->SetTitle(TString::Format("V_{2}^{J/#psi} at #Delta#eta = %.2f", etaGap));
    h->Draw();
  }

  // file_output->ls();

  gPublisherCanvas->finalize();
  file_output->Write();
  // file_output->Close();
}

int main(int argc, char** argv) {
  AssoYieldEtagap(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);
  return 0;
}
