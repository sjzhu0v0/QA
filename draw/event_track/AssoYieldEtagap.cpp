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

// 计算 v2 及误差（独立情形：Cov(a0_low, a2_sub)=Cov(a0_low, a0_sub)=0）
std::pair<double, double> compute_v2_and_error(const TFitResultPtr& result_sub,
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
  double v2_final = sqrt(v2);
  double err_final = 0.5 * sv2 / sqrt(v2);
  return {v2_final, err_final};
}

void AssoYieldEtagap(TString path_input = "/home/szhu/work/alice/analysis/QA/input/event_track/"
                                          "TrackAssocYeild_24apass1.root",
                     TString path_output = "/home/szhu/work/alice/analysis/QA/output/"
                                           "event_track/v2_etagap_24pass1.root",
                     TString path_pdf = "/home/szhu/work/alice/analysis/QA/plot/event_track/"
                                        "v2_etagap_24pass1.pdf") {
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;
  TFile* file_input = new TFile(path_input);
  TFile* file_output = new TFile(path_output, "RECREATE");
  YAML::Node config = YAML::LoadFile("config.yaml");
  const double low_edge_deltaPhiToPi = config["hist_binning"]["low_edge_deltaPhiToPi"].as<double>();
  const double up_edge_deltaPhiToPi = config["hist_binning"]["up_edge_deltaPhiToPi"].as<double>();

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

  MIndexHist index_etaGap(var_EtaGap);

  MHGroupTool1D Hs_highMult(file_input, "h2_highMult_EtaGap%d", {var_EtaGap}, {1});
  MHGroupTool1D Hs_lowMult(file_input, "h2_lowMult_EtaGap%d", {var_EtaGap}, {1});
  MHGroupTool1D Hs_sub(file_input, "h2_sub_EtaGap%d", {var_EtaGap}, {1});
  file_output->cd();
  MHist1D v2_etaGap(index_etaGap, "v2");
  v2_etaGap.fHisto->GetYaxis()->SetTitle("v_{2}^{REF}");

  // Create histograms to store a0_prime, a2_prime, and their ratio
  MHist1D a0prime_etaGap(index_etaGap, "a0_prime");
  a0prime_etaGap.fHisto->GetYaxis()->SetTitle("a0'");

  MHist1D a2prime_etaGap(index_etaGap, "a2_prime");
  a2prime_etaGap.fHisto->GetYaxis()->SetTitle("a2'");

  MHist1D ratio_etaGap(index_etaGap, "ratio_a2prime_a0prime");
  ratio_etaGap.fHisto->GetYaxis()->SetTitle("a2'/a0'");

  gDirectory = nullptr;
  MRootGraphic::StyleCommon();
  gPublisherCanvas = new MPublisherCanvas(path_pdf, 2, 2);

  for (auto iEtagap : index_etaGap) {
    auto h_sub = Hs_sub.GetHist(vector<int>{iEtagap});
    auto h_high = Hs_highMult.GetHist(vector<int>{iEtagap});
    auto h_low = Hs_lowMult.GetHist(vector<int>{iEtagap});
    auto f_sub = (TF1*)h_sub->GetListOfFunctions()->FindObject("f1_modulation");
    auto f_high = (TF1*)h_high->GetListOfFunctions()->FindObject("f1_modulation");
    auto f_low = (TF1*)h_low->GetListOfFunctions()->FindObject("f1_modulation");

    auto result_sub = h_sub->Fit(f_sub, "S Q N R");
    auto result_high = h_high->Fit(f_high, "S Q N R");
    auto result_low = h_low->Fit(f_low, "S Q N R");

    double a0_prime, a2_prime, ratio;
    double err_a0_prime, err_a2_prime, err_ratio;
    ComputePrimeErrorsWithFullCov(result_high, result_low, a0_prime, err_a0_prime, a2_prime,
                                  err_a2_prime, ratio, err_ratio);

    a0prime_etaGap.SetBinInfo(a0_prime, err_a0_prime);
    a2prime_etaGap.SetBinInfo(a2_prime, err_a2_prime);
    ratio_etaGap.SetBinInfo(a2_prime / a0_prime, err_ratio);

    auto pair_v2 = compute_v2_and_error(result_sub, result_low, f_low);
    v2_etaGap.SetBinInfo(pair_v2.first, pair_v2.second);
    cout << "  v2 = " << pair_v2.first << " +/- " << pair_v2.second << endl;
  }
  MRootGraphic::StyleHistCommon(v2_etaGap.fHisto.get());
  gPublisherCanvas->DrawClone(v2_etaGap.fHisto.get(), "Text");
  gPublisherCanvas->finalize();
  file_output->Write();
  file_output->Close();
}

int main(int argc, char** argv) {
  AssoYieldEtagap(argv[1], argv[2], argv[3]);
  return 0;
}
