#include "MFit.h"

// This executable uses the inverse-variance implementation supplied by MFit.h.
#define MFitterPoly MFitterPolyInvSigma2
#include "MHelper.h"
#include "MMath.h"
#include "MRootIO.h"
#include "yaml-cpp/yaml.h"

#include <algorithm>
#include <map>
#include <vector>

void AssoYieldEtagap(TString path_input = "/home/szhu/work/alice/analysis/QA/test/"
                                          "AssoYieldGroupEtagap_NoScale.root",
                     TString path_input_mass =
                         "/home/szhu/work/alice/analysis/QA/input/jpsi/"
                         "JpsiMass_LHC24_apass1_DiElectron.root:PosZUSSingle_MassUSSingle_"
                         "PtUSSingle_NumContribCalibUSSingle",
                     TString path_input_tf1 = "", TString path_input_tf1_assoYield = "",
                     TString path_output = "/home/szhu/work/alice/analysis/QA/test/"
                                           "AssoYieldFit_noScale.root",
                     TString path_config = "config_new.yaml") {
  gErrorIgnoreLevel = kWarning;
  YAML::Node config = YAML::LoadFile(path_config.Data());
  YAML::Node hist_config = config["hist_binning"];

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

  struct FitOrderSelection {
    std::vector<int> ptBins;
    std::vector<int> lowMultOrders;
    std::vector<int> highMultOrders;
  };
  const std::vector<FitOrderSelection> orderSelections = {
      {{1, 2, 3, 4}, {3, 4, 5, 6}, {5, 6, 7, 8}},   // 0 < pT < 2 GeV/c
      {{3, 4, 5, 6}, {2, 3, 4, 5}, {5, 6, 7, 8}},   // 1 < pT < 3 GeV/c
      {{5, 6, 7, 8, 9}, {2, 3, 4, 5}, {5, 6, 7, 8}},      // 2 < pT < 4.5 GeV/c
      {{7, 8, 9, 10, 11, 12}, {3, 4, 5, 6}, {5, 6, 7, 8}} // 3 < pT < 6 GeV/c
  };
  const std::vector<int> massSignalOrders = {6, 7, 8};

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
  hnTool_mass.Rebin(3, 2);

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

  for (auto iPtV2 : indexAnyPtV2Jpsi) {
    auto bins_pt = indexAnyPtV2Jpsi[iPtV2 - 1];
    const auto selection = std::find_if(
        orderSelections.begin(), orderSelections.end(),
        [&bins_pt](const FitOrderSelection& candidate) { return candidate.ptBins == bins_pt; });
    if (selection == orderSelections.end()) continue;

    auto mass_highMult = mass_pt_highMult->ProjectionY(Form("mass_highMult_ptV2_%d", iPtV2),
                                                       bins_pt.front(), bins_pt.back());
    auto mass_lowMult = mass_pt_lowMult->ProjectionY(Form("mass_lowMult_ptV2_%d", iPtV2),
                                                     bins_pt.front(), bins_pt.back());

    // The current study removes the highest associated-yield mass bin and
    // keeps the first physical bin.  The mass normalization removes its four
    // highest bins.  These endpoints select bins [1.80, 4.60] and [1.80, 4.36]
    // respectively, without including the associated histogram underflow bin.
    const double associatedFitMin = 1.80;
    const double associatedFitMax = 4.56;
    const double massFitMin = 1.80;
    const double massFitMax = 4.32;
    auto signal_tf1 = g_tf1_input.GetHist(vector<int>{iPtV2});
    auto signal_tf1_assoYield = g_tf1_assoYield_input.GetHist(vector<int>{iPtV2});

    std::map<int, double> massSignalLow;
    std::map<int, double> massSignalHigh;
    for (const int massOrder : massSignalOrders) {
      MFitterPoly massFitter(mass_highMult, massFitMin, massFitMax);
      massFitter.initializeBasis(massOrder);
      massFitter.inputSignal(signal_tf1);
      massFitter.setHisto(mass_highMult);
      massFitter.fitWithSignal();
      massSignalHigh[massOrder] = massFitter.fNSignal;
      massFitter.setHisto(mass_lowMult);
      massFitter.fitWithSignal();
      massSignalLow[massOrder] = massFitter.fNSignal;
    }

    // An underconstrained mass fit reports a zero signal yield.  Keep the
    // corresponding associated-yield point at zero rather than propagating a
    // division by zero (or a non-finite value) to the output histogram.
    const auto normalizeYield = [](double yield, double nsignal) {
      if (!(nsignal > 0.) || !std::isfinite(yield)) return 0.;
      return yield / nsignal;
    };

    for (auto iEtaGap : indexHistEtaGap) {
      std::map<int, std::map<int, TH1D*>> lowOutputs;
      std::map<int, std::map<int, TH1D*>> highOutputs;
      for (const int massOrder : massSignalOrders)
        for (const int associatedOrder : selection->lowMultOrders) {
          auto* output = new TH1D(
              Form("AssoYield_low_ptV2_%d_EtaGap_%d_assocPoly%d_massPoly%d", iPtV2, iEtaGap,
                   associatedOrder, massOrder),
              Form("low multiplicity, p_{T} bin %d, #Delta#eta gap %d;|#Delta#varphi|;"
                   "normalized associated yield",
                   iPtV2, iEtaGap),
              var_DeltaPhiUS.fNbins, var_DeltaPhiUS.fBins.data());
          output->SetDirectory(file_output);
          lowOutputs[massOrder][associatedOrder] = output;
        }
      for (const int massOrder : massSignalOrders)
        for (const int associatedOrder : selection->highMultOrders) {
          auto* output = new TH1D(
              Form("AssoYield_high_ptV2_%d_EtaGap_%d_assocPoly%d_massPoly%d", iPtV2, iEtaGap,
                   associatedOrder, massOrder),
              Form("high multiplicity, p_{T} bin %d, #Delta#eta gap %d;|#Delta#varphi|;"
                   "normalized associated yield",
                   iPtV2, iEtaGap),
              var_DeltaPhiUS.fNbins, var_DeltaPhiUS.fBins.data());
          output->SetDirectory(file_output);
          highOutputs[massOrder][associatedOrder] = output;
        }

      for (auto iDeltaPhi : indexHistDeltaPhiUS) {
        // Build the gap-selected mass spectra before fitting.  The selected
        // bins are the two signed Δη regions outside the central gap.
        auto* mass_low_gap = static_cast<TH1D*>(
            assoYield_lowMult.GetHist(vector<int>{1, iDeltaPhi, iPtV2})
                ->Clone(Form("mass_low_etaGap_%d_deltaPhi_%d_ptV2_%d", iEtaGap, iDeltaPhi,
                             iPtV2)));
        auto* mass_high_gap = static_cast<TH1D*>(
            assoYield_highMult.GetHist(vector<int>{1, iDeltaPhi, iPtV2})
                ->Clone(Form("mass_high_etaGap_%d_deltaPhi_%d_ptV2_%d", iEtaGap, iDeltaPhi,
                             iPtV2)));
        mass_low_gap->Reset("ICES");
        mass_high_gap->Reset("ICES");
        mass_low_gap->SetDirectory(nullptr);
        mass_high_gap->SetDirectory(nullptr);

        const int last_negative_bin = n_bins_deltaEta_assoYield / 2 - iEtaGap + 1;
        const int first_positive_bin = n_bins_deltaEta_assoYield / 2 + iEtaGap - 1;
        for (int iDeltaEta = 1; iDeltaEta <= last_negative_bin; ++iDeltaEta) {
          mass_low_gap->Add(
              assoYield_lowMult.GetHist(vector<int>{iDeltaEta, iDeltaPhi, iPtV2}));
          mass_high_gap->Add(
              assoYield_highMult.GetHist(vector<int>{iDeltaEta, iDeltaPhi, iPtV2}));
        }
        for (int iDeltaEta = first_positive_bin; iDeltaEta <= n_bins_deltaEta_assoYield;
             ++iDeltaEta) {
          mass_low_gap->Add(
              assoYield_lowMult.GetHist(vector<int>{iDeltaEta, iDeltaPhi, iPtV2}));
          mass_high_gap->Add(
              assoYield_highMult.GetHist(vector<int>{iDeltaEta, iDeltaPhi, iPtV2}));
        }

        for (const int associatedOrder : selection->lowMultOrders) {
          MFitterPoly associatedFitter(mass_low_gap, associatedFitMin, associatedFitMax);
          associatedFitter.initializeBasis(associatedOrder);
          associatedFitter.inputSignal(signal_tf1_assoYield);
          associatedFitter.fitWithSignal();
          for (const int massOrder : massSignalOrders)
            lowOutputs[massOrder][associatedOrder]->SetBinContent(
                iDeltaPhi, normalizeYield(associatedFitter.fNSignal, massSignalLow[massOrder]));
        }
        for (const int associatedOrder : selection->highMultOrders) {
          MFitterPoly associatedFitter(mass_high_gap, associatedFitMin, associatedFitMax);
          associatedFitter.initializeBasis(associatedOrder);
          associatedFitter.inputSignal(signal_tf1_assoYield);
          associatedFitter.fitWithSignal();
          for (const int massOrder : massSignalOrders)
            highOutputs[massOrder][associatedOrder]->SetBinContent(
                iDeltaPhi, normalizeYield(associatedFitter.fNSignal, massSignalHigh[massOrder]));
        }

        delete mass_low_gap;
        delete mass_high_gap;
      }
    }
  }

  file_output->Write();
  file_output->Close();
}

#undef MFitterPoly

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
  TString path_config = argc > 6 ? argv[6] : "config_new.yaml";

  AssoYieldEtagap(path_input, path_input_mass, path_input_tf1, path_input_assoYield_tf1,
                  path_output, path_config);

  return 0;
}
