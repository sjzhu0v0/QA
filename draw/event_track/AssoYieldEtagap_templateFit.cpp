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
  gDirectory = nullptr;
  MRootGraphic::StyleCommon();
  gPublisherCanvas = new MPublisherCanvas(path_pdf, 2, 2);

  for (auto iEtagap : index_etaGap) {
    auto h_high = Hs_highMult.GetHist(vector<int>{iEtagap});
    auto h_low = Hs_lowMult.GetHist(vector<int>{iEtagap});
    auto result = fit(h_high, h_low);
    double val_V2 = result.V2;
    // if (val_V2 < 0) {
    //   std::cerr << "Warning: V2 < 0 for eta gap bin " << iEtagap << ", setting V2 to 0."
    //             << std::endl;
    //   exit(1);
    // }
    double err_V2 = result.V2_err;
    double val_v2 = sqrt(val_V2);
    double err_v2 = 0.5 * err_V2 / val_v2;

    v2_etaGap.SetBinInfo(val_v2, err_v2);
    cout << "Etagap bin " << iEtagap << ": v2 = " << val_v2 << " +/- " << err_v2 << endl;
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
