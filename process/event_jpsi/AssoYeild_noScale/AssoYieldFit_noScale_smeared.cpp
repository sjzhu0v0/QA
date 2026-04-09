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
    TString path_input_mass =
        "/home/szhu/work/alice/analysis/QA/input/jpsi/"
        "JpsiMass_LHC24_apass1_DiElectron.root:PosZUSSingle_MassUSSingle_"
        "PtUSSingle_NumContribCalibUSSingle",
    TString path_input_tf1 = "",
    TString path_output = "/home/szhu/work/alice/analysis/QA/test/"
                          "AssoYieldFit_noScale.root",
    TString path_pdf = "/home/szhu/work/alice/analysis/QA/test/"
                       "AssoYieldFit_noScale.pdf",
    TString path_config = "config_new.yaml") {
  gErrorIgnoreLevel = kWarning;
  YAML::Node config = YAML::LoadFile(path_config.Data());
  YAML::Node hist_config = config["hist_binning"];

  TFile *file_input = new TFile(path_input);
  TFile *file_input_tf1 = new TFile(path_input_tf1);
  TFile *file_output = new TFile(path_output, "RECREATE");

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

  StrVar4Hist var_fPosZ = ParseStrVar4Hist(hist_config["fPosZ"]);
  StrVar4Hist var_NumContribCalibBinned = ParseStrVar4Hist(hist_config["NumContribCalibBinned"]);
  StrVar4Hist var_MassJpsiCandidate = ParseStrVar4Hist(hist_config["MassJpsiCandidate"]);
  StrVar4Hist var_PtJpsiCandidate = ParseStrVar4Hist(hist_config["PtJpsiCandidate"]);
  StrVar4Hist var_DeltaEtaUS = ParseStrVar4Hist(hist_config["DeltaEtaUS"]);
  StrVar4Hist var_DeltaPhiUS = ParseStrVar4Hist(hist_config["DeltaPhiUS"]);
  
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

  int n_rebin_mass_assoYield =
      config["hist_binning"]["n_rebin_mass_assoYield"].as<int>();
  MIndexHist indexHistMass(var_MassJpsiCandidate, 1, n_rebin_mass_assoYield);
  MIndexHist indexHistPtJpsiCandidate(var_PtJpsiCandidate, 1, 1);
  MIndexHist indexHistDeltaPhiUS(var_DeltaPhiUS, 1, 1);
  int n_rebin_deltaEta_assoYield =
      config["hist_binning"]["n_rebin_deltaEta_assoYield"].as<int>();
  MIndexHist indexHistDeltaEtaUS(var_DeltaEtaUS, 1, n_rebin_deltaEta_assoYield);
  MIndexHist indexHistEtaGap(var_EtaGap, 1, 1);
  MIndexHist indexHistPtV2Jpsi(var_PtV2Jpsi, 1, 1);
  MIndexAny indexAnyPtV2Jpsi(strAny_ptV2, 1);

  MHGroupTool<TH1D> g_tf1_input(file_input_tf1, "mass_template_pt_%d",
                               {var_PtV2Jpsi}, {1});

  auto h_mass = MRootIO::GetObjectDiectly<THnD>(path_input_mass);
  MHnTool hnTool_mass(h_mass);
  // hnTool_mass.Rebin(0, 200);
  hnTool_mass.Rebin(3, 2);
  hnTool_mass.PrintAllAxis();

  TH2D *mass_pt_lowMult = hnTool_mass.Project(1, 2, {0, 1});
  mass_pt_lowMult->Add(hnTool_mass.Project(1, 2, {0, 2}));
  mass_pt_lowMult->Add(hnTool_mass.Project(1, 2, {0, 3}));
  mass_pt_lowMult->Add(hnTool_mass.Project(1, 2, {0, 4}));
  TH2D *mass_pt_highMult = hnTool_mass.Project(1, 2, {0, 5});
  MFitterPoly fitterPoly_mass(mass_pt_highMult->ProjectionY(), 1.88, 4.32);
  fitterPoly_mass.initializeBasis(6);

  MHGroupTool1D assoYield_lowMult(
      file_input,
      "MassUS_AssoYield_lowMult_DeltaEtaUS_%d_DeltaPhiUS_%d_ptV2_%d",
      {var_DeltaEtaUS, var_DeltaPhiUS, var_PtV2Jpsi},
      {n_rebin_deltaEta_assoYield, 1, 1});
  MHGroupTool1D assoYield_highMult(
      file_input,
      "MassUS_AssoYield_highMult_DeltaEtaUS_%d_DeltaPhiUS_%d_ptV2_%d",
      {var_DeltaEtaUS, var_DeltaPhiUS, var_PtV2Jpsi},
      {n_rebin_deltaEta_assoYield, 1, 1});
  MFitterPoly fitterPoly_asso(assoYield_highMult.GetHist(vector<int>{1, 1, 1}),
                              1.88, 4.32);
  fitterPoly_asso.initializeBasis(6);

  gDirectory = nullptr;
  MHist1D h1_assoYield_sub(indexHistDeltaPhiUS, "AssoYield_sub");
  MHist1D h1_assoYield_low(indexHistDeltaPhiUS, "AssoYield_low");
  MHist1D h1_assoYield_high(indexHistDeltaPhiUS, "AssoYield_high");

  gDirectory = file_output;
  TH1D *nsignal_highMult_ptV2 =
      new TH1D("nsignal_highMult_ptV2", "nsignal_highMult_ptV2",
               strAny_ptV2.fNbins, 0, strAny_ptV2.fNbins);
  TH1D *nsignal_lowMult_ptV2 =
      new TH1D("nsignal_lowMult_ptV2", "nsignal_lowMult_ptV2",
               strAny_ptV2.fNbins, 0, strAny_ptV2.fNbins);

  auto vec_assoYield_sub =
      MakeMVec(h1_assoYield_sub, indexHistDeltaEtaUS, indexAnyPtV2Jpsi);
  auto vec_assoYield_low =
      MakeMVec(h1_assoYield_low, indexHistDeltaEtaUS, indexAnyPtV2Jpsi);
  auto vec_assoYield_high =
      MakeMVec(h1_assoYield_high, indexHistDeltaEtaUS, indexAnyPtV2Jpsi);

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

    auto mass_highMult = mass_pt_highMult->ProjectionY(
        Form("mass_highMult_ptV2_%d", iPtV2), bins_pt.front(), bins_pt.back());
    auto mass_lowMult = mass_pt_lowMult->ProjectionY(
        Form("mass_lowMult_ptV2_%d", iPtV2), bins_pt.front(), bins_pt.back());

    auto signal_tf1 = g_tf1_input.GetHist(vector<int>{iPtV2});
    fitterPoly_mass.inputSignal(signal_tf1);
    fitterPoly_asso.inputSignal(signal_tf1);

    fitterPoly_mass.setHisto(mass_highMult);
    fitterPoly_mass.fitWithSignal();
    double nsignal_highMult = fitterPoly_mass.fNSignal;
    fitterPoly_mass.setHisto(mass_lowMult);
    fitterPoly_mass.fitWithSignal();
    double nsignal_lowMult = fitterPoly_mass.fNSignal;

    for (auto i_deltaEta : indexHistDeltaEtaUS)
      for (auto i_deltaPhi : indexHistDeltaPhiUS) {
        auto assoYield_lowMult_diff = assoYield_lowMult.GetHist(
            vector<int>{i_deltaEta, i_deltaPhi, iPtV2});
        auto assoYield_highMult_diff = assoYield_highMult.GetHist(
            vector<int>{i_deltaEta, i_deltaPhi, iPtV2});

        fitterPoly_asso.setHisto(assoYield_highMult_diff);
        fitterPoly_asso.fitWithSignal();
        double nyield_highmult = fitterPoly_asso.fNSignal;

        fitterPoly_asso.setHisto(assoYield_lowMult_diff);
        fitterPoly_asso.fitWithSignal();
        double nyield_lowmult = fitterPoly_asso.fNSignal;

        double normalized_nyield_highmult = nyield_highmult / nsignal_highMult;
        double normalized_nyield_lowmult = nyield_lowmult / nsignal_lowMult;

        double normalized_nyield_sub =
            normalized_nyield_highmult - normalized_nyield_lowmult;
        vec_assoYield_sub.currentObject().SetBinInfo(normalized_nyield_sub, 0);
        vec_assoYield_low.currentObject().SetBinInfo(normalized_nyield_lowmult,
                                                     0);
        vec_assoYield_high.currentObject().SetBinInfo(
            normalized_nyield_highmult, 0);
      }
  }

  file_output->Write();
  file_output->Close();
}

int main(int argc, char **argv) {
  TString path_input = argc > 1 ? argv[1]
                                : "/home/szhu/work/alice/analysis/QA/test/"
                                  "AssoYieldGroupEtagap_NoScale.root";
  TString path_input_mass =
      argc > 2 ? argv[2]
               : "/home/szhu/work/alice/analysis/QA/input/jpsi/"
                 "JpsiMass_LHC24_apass1_DiElectron.root";
  TString path_input_tf1 = argc > 3 ? argv[3] : "";
  TString path_output = argc > 4 ? argv[4]
                                 : "/home/szhu/work/alice/analysis/QA/test/"
                                   "AssoYeilFit_noScale.root";
  TString path_pdf = argc > 5 ? argv[5]
                              : "/home/szhu/work/alice/analysis/QA/test/"
                                "AssoYieldFit_noScale.pdf";
  TString path_config = argc > 6 ? argv[6] : "config_new.yaml";

  AssoYieldEtagap(path_input, path_input_mass, path_input_tf1, path_output,
                  path_pdf, path_config);

  return 0;
}
