#include "MHist.h"
#include "MRootGraphic.h"
#include "yaml-cpp/yaml.h"

#include <TMinuit.h>
#include <cmath>

static TH1D* hHigh = nullptr;
static TH1D* hLow = nullptr;

struct StrResult {
  double F, F_err;
  double A0, A0_err;
  double A2, A2_err;
  double A3, A3_err;
  double V2, V2_err;
  double V3, V3_err;
  double chi2;
  int ndf;
  void Print() const {
    std::cout << std::fixed << std::setprecision(4);

    std::cout << "===== Fit Result =====" << std::endl;
    std::cout << "F   = " << F << " ± " << F_err << std::endl;
    std::cout << "A0  = " << A0 << " ± " << A0_err << std::endl;
    std::cout << "A2  = " << A2 << " ± " << A2_err << std::endl;
    std::cout << "A3  = " << A3 << " ± " << A3_err << std::endl;

    std::cout << std::endl;
    std::cout << "V2  = " << V2 << " ± " << V2_err << std::endl;
    std::cout << "V3  = " << V3 << " ± " << V3_err << std::endl;

    std::cout << std::endl;
    std::cout << "chi2 / ndf = " << chi2 << " / " << ndf << " = " << chi2 / ndf << std::endl;

    std::cout << "======================" << std::endl;
  }
};

double chi2(const double* p) {
  double F = p[0];
  double A0 = p[1];
  double A2 = p[2];
  double A3 = p[3];

  double chi2 = 0.0;

  for (int i = 1; i <= hHigh->GetNbinsX(); ++i) {
    double phi = hHigh->GetBinCenter(i);
    double dphi = hHigh->GetBinWidth(i);

    double H = hHigh->GetBinContent(i);
    double L = hLow->GetBinContent(i);

    double sH = hHigh->GetBinError(i);
    double sL = hLow->GetBinError(i);

    double model = F * L + A0 + A2 * cos(2 * phi) + A3 * cos(3 * phi);

    double sigma2 = sH * sH + F * F * sL * sL;
    if (sigma2 <= 0)
      continue;

    chi2 += (H - model) * (H - model) / sigma2;
  }
  return chi2;
}

void fcn(int& npar, double*, double& f, double* par, int) { f = chi2(par); }

StrResult fit(TH1D* hist_high, TH1D* hist_low) {
  hHigh = hist_high;
  hLow = hist_low;

  // ---------- 初值估计 ----------
  double initF = hist_high->Integral() / hist_low->Integral();
  double initA0 = hist_high->Integral() / (2.0 * TMath::Pi());
  double initA2 = 0.0;
  double initA3 = 0.0;

  TMinuit minuit(4);
  minuit.SetFCN(fcn);
  minuit.SetPrintLevel(-1);

  int ierflg = 0;
  double arglist[10];

  // ---------- 参数定义 + 约束 ----------
  minuit.DefineParameter(0, "F", initF, 0.01 * initF, 0.2 * initF, 5.0 * initF);

  minuit.DefineParameter(1, "A0", initA0, 0.01 * initA0, 1e-6, 10.0 * initA0); // A0 > 0 强制

  minuit.DefineParameter(2, "A2", initA2, 0.001 * initA0, -5.0 * initA0, 5.0 * initA0);

  minuit.DefineParameter(3, "A3", initA3, 0.001 * initA0, -5.0 * initA0, 5.0 * initA0);

  // ---------- MIGRAD ----------
  arglist[0] = 5000;
  arglist[1] = 1e-6;
  minuit.mnexcm("MIGRAD", arglist, 2, ierflg);

  // ---------- 取结果 ----------
  StrResult res;
  minuit.GetParameter(0, res.F, res.F_err);
  minuit.GetParameter(1, res.A0, res.A0_err);
  minuit.GetParameter(2, res.A2, res.A2_err);
  minuit.GetParameter(3, res.A3, res.A3_err);

  double amin, edm, errdef;
  int nvpar, nparx, icstat;
  minuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);

  res.chi2 = amin;
  res.ndf = hist_high->GetNbinsX() - nvpar;

  // ---------- Vn & 误差 ----------
  res.V2 = res.A2 / (2.0 * res.A0);
  res.V3 = res.A3 / (2.0 * res.A0);

  res.V2_err = std::abs(res.V2) *
               std::sqrt(std::pow(res.A2_err / res.A2, 2) + std::pow(res.A0_err / res.A0, 2));

  res.V3_err = std::abs(res.V3) *
               std::sqrt(std::pow(res.A3_err / res.A3, 2) + std::pow(res.A0_err / res.A0, 2));

  return res;
}
void AssoYieldEtagap(
    TString path_input = "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
                         "JpsiAssocYield_24apass1.root",
    TString path_hist_ref = "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
                         "JpsiAssocYield_24apass1.root:hist",
    TString path_output = "/home/szhu/work/alice/analysis/QA/output/event_jpsi/"
                          "AssoYieldEtagap.root",
    TString path_pdf = "/home/szhu/work/alice/analysis/QA/plot/event_jpsi/"
                       "AssoYieldEtagap") {
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
                                      {4, 5,6, 7, 8},
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

  gPublisherCanvas = new MPublisherCanvas(path_pdf+".pdf", 3, 1);

#define MH1DGetBin(...) GetHist(vector<int>{__VA_ARGS__})

  auto v2REF_etaGap = MRootIO::GetObjectDiectly<TH1D>(path_hist_ref);

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
      auto h_high = DeltaPhi_high.MH1DGetBin(iEtagpa, iPt);
      auto h_low = DeltaPhi_low.MH1DGetBin(iEtagpa, iPt);
      auto result = fit(h_high, h_low);
      double val_V2 = result.V2;
      double err_V2 = result.V2_err;

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
  c_v2_pT_etaGap->SaveAs(path_pdf + "_v2Jpsi_etaGap.pdf");

  TCanvas *c_v2_pT_etaGap1p2 =
      new TCanvas("c_v2_pT_etaGap1p2", "c_v2_pT_etaGap1p2", 600, 600);
  c_v2_pT_etaGap1p2->cd();
  auto h_v2_etaGap1p2 = (TH1D *)file_output->Get("v2Jpsi_pT_etaGapbin8");
  // h_v2_etaGap1p2->GetYaxis()->SetRangeUser();
  h_v2_etaGap1p2->SetTitle("v_{2}^{J/#psi} at #Delta#eta = 1.2");
  h_v2_etaGap1p2->Draw();
  c_v2_pT_etaGap1p2->SaveAs(
      path_pdf + "_v2_pT_etaGap1p2.pdf");

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

int main(int argc, char **argv) {
  AssoYieldEtagap(argv[1], argv[2], argv[3], argv[4]);
  return 0;
}
