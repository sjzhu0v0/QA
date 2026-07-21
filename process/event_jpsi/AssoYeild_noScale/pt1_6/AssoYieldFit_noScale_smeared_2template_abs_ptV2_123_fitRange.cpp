#include "MFit.h"
#include "MHelper.h"
#include "MMath.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"
#include "yaml-cpp/yaml.h"

void AssoYieldEtagap(TString path_input = "/home/szhu/work/alice/analysis/QA/test/"
                                          "AssoYieldGroupEtagap_NoScale.root",
                     TString path_input_mass =
                         "/home/szhu/work/alice/analysis/QA/input/jpsi/"
                         "JpsiMass_LHC24_apass1_DiElectron.root:PosZUSSingle_MassUSSingle_"
                         "PtUSSingle_NumContribCalibUSSingle",
                     TString path_input_tf1 = "", TString path_input_tf1_assoYield = "",
                     TString path_output = "/home/szhu/work/alice/analysis/QA/test/"
                                           "AssoYieldFit_noScale.root",
                     TString path_pdf = "/home/szhu/work/alice/analysis/QA/test/"
                                        "AssoYieldFit_noScale.pdf",
                     TString path_config = "config_new.yaml") {
  gErrorIgnoreLevel = kWarning;
  YAML::Node config = YAML::LoadFile(path_config.Data());
  YAML::Node hist_config = config["hist_binning"];
  const int poly_order = config["poly_order"].as<int>();

  TFile* file_input = new TFile(path_input);
  TFile* file_input_tf1 = new TFile(path_input_tf1);
  TFile* file_input_tf1_assoYield = new TFile(path_input_tf1_assoYield);
  TFile* file_output = new TFile(path_output, "RECREATE");

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
                                      {11},
                                      {12},
                                      {1, 2},
                                      {2, 3},
                                      {3, 4},
                                      {4, 5},
                                      {5, 6},
                                      {6, 7},
                                      {7, 8},
                                      {8, 9},
                                      {9, 10},
                                      {10, 11},
                                      {11, 12},
                                      {1, 2, 3},
                                      {2, 3, 4},
                                      {3, 4, 5},
                                      {4, 5, 6},
                                      {5, 6, 7},
                                      {6, 7, 8},
                                      {7, 8, 9},
                                      {8, 9, 10},
                                      {9, 10, 11},
                                      {10, 11, 12},
                                      {1, 2, 3, 4},
                                      {2, 3, 4, 5},
                                      {3, 4, 5, 6},
                                      {4, 5, 6, 7},
                                      {5, 6, 7, 8},
                                      {6, 7, 8, 9},
                                      {7, 8, 9, 10},
                                      {8, 9, 10, 11},
                                      {9, 10, 11, 12},
                                      {1, 2, 3, 4, 5},
                                      {2, 3, 4, 5, 6},
                                      {3, 4, 5, 6, 7},
                                      {4, 5, 6, 7, 8},
                                      {5, 6, 7, 8, 9},
                                      {6, 7, 8, 9, 10},
                                      {7, 8, 9, 10, 11},
                                      {8, 9, 10, 11, 12},
                                      {1, 2, 3, 4, 5, 6},
                                      {2, 3, 4, 5, 6, 7},
                                      {3, 4, 5, 6, 7, 8},
                                      {4, 5, 6, 7, 8, 9},
                                      {5, 6, 7, 8, 9, 10},
                                      {6, 7, 8, 9, 10, 11},
                                      {7, 8, 9, 10, 11, 12},
                                      {1, 2, 3, 4, 5, 6, 7},
                                      {2, 3, 4, 5, 6, 7, 8},
                                      {3, 4, 5, 6, 7, 8, 9},
                                      {4, 5, 6, 7, 8, 9, 10},
                                      {5, 6, 7, 8, 9, 10, 11},
                                      {6, 7, 8, 9, 10, 11, 12},
                                      {1, 2, 3, 4, 5, 6, 7, 8},
                                      {2, 3, 4, 5, 6, 7, 8, 9},
                                      {3, 4, 5, 6, 7, 8, 9, 10},
                                      {4, 5, 6, 7, 8, 9, 10, 11},
                                      {5, 6, 7, 8, 9, 10, 11, 12},
                                      {1, 2, 3, 4, 5, 6, 7, 8, 9},
                                      {2, 3, 4, 5, 6, 7, 8, 9, 10},
                                      {3, 4, 5, 6, 7, 8, 9, 10, 11},
                                      {4, 5, 6, 7, 8, 9, 10, 11, 12},
                                      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
                                      {2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
                                      {3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                                      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
                                      {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
                                      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}};
    const int fNbins = bins.size();
    const TString fName = "ptV2";

    vector<int> operator[](int index) { return bins[index]; }
  } strAny_ptV2;

  StrVar4Hist var_fPosZ = ParseStrVar4Hist(hist_config["fPosZ"]);
  StrVar4Hist var_NumContribCalibBinned = ParseStrVar4Hist(hist_config["NumContribCalibBinned"]);
  StrVar4Hist var_MassJpsiCandidate = ParseStrVar4Hist(hist_config["MassJpsiCandidate"]);
  StrVar4Hist var_PtJpsiCandidate = ParseStrVar4Hist(hist_config["PtJpsiCandidate"]);
  StrVar4Hist var_DeltaEtaUS = ParseStrVar4Hist(hist_config["DeltaEtaUS"]);
  StrVar4Hist var_DeltaPhiUS = ParseStrVar4Hist(hist_config["AbsDeltaPhiUS"]);

  // Calculate EtaGap binning based on DeltaEtaUS configuration
  int n_bins_deltaEta_assoYield = var_DeltaEtaUS.fNbins;
  double min_deltaEta = var_DeltaEtaUS.fBins[0];
  double max_deltaEta = var_DeltaEtaUS.fBins.back();
  double bin_width_etaGap = (max_deltaEta - min_deltaEta) / n_bins_deltaEta_assoYield;
  int n_bins_etaGap = n_bins_deltaEta_assoYield / 2 - 2;
  std::vector<double> etaGap_bins = {-1.0 * bin_width_etaGap,
                                     (n_bins_deltaEta_assoYield / 2 - 3) * bin_width_etaGap};
  StrVar4Hist var_EtaGap("EtaGap", "#Delta#eta_{gap}", "", n_bins_etaGap, etaGap_bins);
  StrVar4Hist var_PtV2Jpsi("PtV2Jpsi", "p_{T}", "GeV/c", strAny_ptV2.fNbins, {0., 1.});

  int n_rebin_mass_assoYield = config["hist_binning"]["n_rebin_mass_assoYield"].as<int>();
  MIndexHist indexHistMass(var_MassJpsiCandidate, 1, n_rebin_mass_assoYield);
  MIndexHist indexHistPtJpsiCandidate(var_PtJpsiCandidate, 1, 1);
  MIndexHist indexHistDeltaPhiUS(var_DeltaPhiUS, 1, 1);
  int n_rebin_deltaEta_assoYield = config["hist_binning"]["n_rebin_deltaEta_assoYield"].as<int>();
  MIndexHist indexHistDeltaEtaUS(var_DeltaEtaUS, 1, n_rebin_deltaEta_assoYield);
  MIndexHist indexHistEtaGap(var_EtaGap, 1, 1);
  MIndexHist indexHistPtV2Jpsi(var_PtV2Jpsi, 1, 1);
  MIndexAny indexAnyPtV2Jpsi(strAny_ptV2, 1);

  MHGroupTool<TH1D> g_tf1_input(file_input_tf1, "mass_template_pt_%d", {var_PtV2Jpsi}, {1});
  MHGroupTool<TH1D> g_tf1_assoYield_input(file_input_tf1_assoYield, "mass_template_pt_%d",
                                          {var_PtV2Jpsi}, {1});

  auto h_mass = MRootIO::GetObjectDiectly<THnD>(path_input_mass);
  MHnTool hnTool_mass(h_mass);
  // hnTool_mass.Rebin(0, 200);
  hnTool_mass.Rebin(3, 2);
  hnTool_mass.PrintAllAxis();

  TH2D* mass_pt_lowMult = hnTool_mass.Project(1, 2, {0, 1});
  mass_pt_lowMult->Add(hnTool_mass.Project(1, 2, {0, 2}));
  mass_pt_lowMult->Add(hnTool_mass.Project(1, 2, {0, 3}));
  mass_pt_lowMult->Add(hnTool_mass.Project(1, 2, {0, 4}));
  TH2D* mass_pt_highMult = hnTool_mass.Project(1, 2, {0, 5});

  MHGroupTool1D assoYield_lowMult(
      file_input, "jpsi_mass_AssoYield_lowMult_DeltaEta_%d_" + var_DeltaPhiUS.fName + "_%d_ptV2_%d",
      {var_DeltaEtaUS, var_DeltaPhiUS, var_PtV2Jpsi}, {n_rebin_deltaEta_assoYield, 1, 1});
  MHGroupTool1D assoYield_highMult(
      file_input,
      "jpsi_mass_AssoYield_highMult_DeltaEta_%d_" + var_DeltaPhiUS.fName + "_%d_ptV2_%d",
      {var_DeltaEtaUS, var_DeltaPhiUS, var_PtV2Jpsi}, {n_rebin_deltaEta_assoYield, 1, 1});

  gDirectory = nullptr;
  MHist1D h1_assoYield_sub(indexHistDeltaPhiUS, "AssoYield_sub");
  MHist1D h1_assoYield_low(indexHistDeltaPhiUS, "AssoYield_low");
  MHist1D h1_assoYield_high(indexHistDeltaPhiUS, "AssoYield_high");

  gDirectory = file_output;
  TH1D* nsignal_highMult_ptV2 = new TH1D("nsignal_highMult_ptV2", "nsignal_highMult_ptV2",
                                         strAny_ptV2.fNbins, 0, strAny_ptV2.fNbins);
  TH1D* nsignal_lowMult_ptV2 = new TH1D("nsignal_lowMult_ptV2", "nsignal_lowMult_ptV2",
                                        strAny_ptV2.fNbins, 0, strAny_ptV2.fNbins);

  auto vec_assoYield_sub = MakeMVec(h1_assoYield_sub, indexHistDeltaEtaUS, indexAnyPtV2Jpsi);
  auto vec_assoYield_low = MakeMVec(h1_assoYield_low, indexHistDeltaEtaUS, indexAnyPtV2Jpsi);
  auto vec_assoYield_high = MakeMVec(h1_assoYield_high, indexHistDeltaEtaUS, indexAnyPtV2Jpsi);

  for (auto _ : indexAnyPtV2Jpsi)
    for (auto __ : indexHistDeltaEtaUS) {
      vec_assoYield_sub.currentObject().fHisto->SetDirectory(file_output);
      vec_assoYield_low.currentObject().fHisto->SetDirectory(file_output);
      vec_assoYield_high.currentObject().fHisto->SetDirectory(file_output);
    }

  gDirectory = nullptr;

  for (auto iPtV2 : indexAnyPtV2Jpsi) {
    auto bins_pt = indexAnyPtV2Jpsi[iPtV2 - 1];

    TString str_bins_pt = Form("%d p_{T} (GeV/c): ", iPtV2);
    for (auto bin : bins_pt)
      str_bins_pt += Form("%d_", bin);

    auto mass_highMult = mass_pt_highMult->ProjectionY(Form("mass_highMult_ptV2_%d", iPtV2),
                                                       bins_pt.front(), bins_pt.back());
    auto mass_lowMult = mass_pt_lowMult->ProjectionY(Form("mass_lowMult_ptV2_%d", iPtV2),
                                                     bins_pt.front(), bins_pt.back());

    const bool use_narrow_fit_range =
        bins_pt.size() == 1 && bins_pt.front() >= 1 && bins_pt.front() <= 3;
    const double fit_min = use_narrow_fit_range ? 2.6 : 1.88;
    const double fit_max = 4.32;

    auto signal_tf1 = g_tf1_input.GetHist(vector<int>{iPtV2});
    MFitterPoly fitterPoly_mass(mass_highMult, fit_min, fit_max);
    fitterPoly_mass.initializeBasis(poly_order);
    fitterPoly_mass.inputSignal(signal_tf1);

    auto signal_tf1_assoYield = g_tf1_assoYield_input.GetHist(vector<int>{iPtV2});
    MFitterPoly fitterPoly_asso(assoYield_highMult.GetHist(vector<int>{1, 1, iPtV2}), fit_min,
                                fit_max);
    fitterPoly_asso.initializeBasis(poly_order);
    fitterPoly_asso.inputSignal(signal_tf1_assoYield);

    fitterPoly_mass.setHisto(mass_highMult);
    fitterPoly_mass.fitWithSignal();
    double nsignal_highMult = fitterPoly_mass.fNSignal;
    fitterPoly_mass.setHisto(mass_lowMult);
    fitterPoly_mass.fitWithSignal();
    double nsignal_lowMult = fitterPoly_mass.fNSignal;

    // An underconstrained mass fit reports a zero signal yield.  Keep the
    // corresponding associated-yield point at zero rather than propagating a
    // division by zero (or a non-finite value) to the output histogram.
    const auto normalizeYield = [](double yield, double nsignal) {
      if (!(nsignal > 0.) || !std::isfinite(yield)) return 0.;
      return yield / nsignal;
    };

    for (auto i_deltaEta : indexHistDeltaEtaUS)
      for (auto i_deltaPhi : indexHistDeltaPhiUS) {
        auto assoYield_lowMult_diff =
            assoYield_lowMult.GetHist(vector<int>{i_deltaEta, i_deltaPhi, iPtV2});
        auto assoYield_highMult_diff =
            assoYield_highMult.GetHist(vector<int>{i_deltaEta, i_deltaPhi, iPtV2});

        fitterPoly_asso.setHisto(assoYield_highMult_diff);
        fitterPoly_asso.fitWithSignal();
        double nyield_highmult = fitterPoly_asso.fNSignal;

        fitterPoly_asso.setHisto(assoYield_lowMult_diff);
        fitterPoly_asso.fitWithSignal();
        double nyield_lowmult = fitterPoly_asso.fNSignal;

        double normalized_nyield_highmult = normalizeYield(nyield_highmult, nsignal_highMult);
        double normalized_nyield_lowmult = normalizeYield(nyield_lowmult, nsignal_lowMult);

        double normalized_nyield_sub = normalized_nyield_highmult - normalized_nyield_lowmult;
        vec_assoYield_sub.currentObject().SetBinInfo(normalized_nyield_sub, 0);
        vec_assoYield_low.currentObject().SetBinInfo(normalized_nyield_lowmult, 0);
        vec_assoYield_high.currentObject().SetBinInfo(normalized_nyield_highmult, 0);
      }
  }

  file_output->Write();
  file_output->Close();
}

int main(int argc, char** argv) {
  TString path_input = argc > 1 ? argv[1]
                                : "/home/szhu/work/alice/analysis/QA/test/"
                                  "AssoYieldGroupEtagap_NoScale.root";
  TString path_input_mass = argc > 2 ? argv[2]
                                     : "/home/szhu/work/alice/analysis/QA/input/jpsi/"
                                       "JpsiMass_LHC24_apass1_DiElectron.root";
  TString path_input_tf1 = argc > 3 ? argv[3] : "";
  TString path_input_assoYield_tf1 = argc > 4 ? argv[4] : "";
  TString path_output = argc > 5 ? argv[5]
                                 : "/home/szhu/work/alice/analysis/QA/test/"
                                   "AssoYeilFit_noScale.root";
  TString path_pdf = argc > 6 ? argv[6]
                              : "/home/szhu/work/alice/analysis/QA/test/"
                                "AssoYieldFit_noScale.pdf";
  TString path_config = argc > 7 ? argv[7] : "config_new.yaml";

  AssoYieldEtagap(path_input, path_input_mass, path_input_tf1, path_input_assoYield_tf1,
                  path_output, path_pdf, path_config);

  return 0;
}
