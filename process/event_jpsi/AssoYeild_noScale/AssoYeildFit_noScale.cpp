#include "MFit.h"
#include "MHelper.h"
#include "MMath.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"

funcWithJson(void, AssoYeildFit_noScale)(
    TString path_input = "/home/szhu/work/alice/analysis/QA/test/"
                         "AssoYeildGroupEtagap_NoScale.root",
    TString path_input_mass = "/home/szhu/work/alice/analysis/QA/input/jpsi/"
                              "JpsiMass_LHC24_apass1_DiElectron.root",
    TString path_input_tf1 = "",
    TString path_output = "/home/szhu/work/alice/analysis/QA/test/"
                          "AssoYeildFit_noScale.root",
    TString path_pdf = "/home/szhu/work/alice/analysis/QA/test/"
                       "AssoYeildFit_noScale.pdf") {
  gErrorIgnoreLevel = kWarning;
  SetUpJson();
  Configurable<int> config_n_rebin_mass("n_rebin_mass", 3);
  int n_rebin_mass = config_n_rebin_mass.data;

  TFile *file_input = new TFile(path_input);
  TFile *file_input_tf1 = new TFile(path_input_tf1);
  TFile *file_input_mass = new TFile(path_input_mass);
  TFile *file_output = new TFile(path_output, "RECREATE");

  gDirectory = nullptr;

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
  StrVar4Hist var_DeltaEtaUS("DeltaEtaUS", "#Delta#eta_{J/#psi, track}", "", 80,
                             {-4., 4.});
  StrVar4Hist var_DeltaPhiUS("DeltaPhiUS", "#Delta#phi_{J/#psi, track}", "", 10,
                             {-M_PI_2, M_PI + M_PI_2});
  StrVar4Hist var_EtaGap("EtaGap", "#Delta#eta_{gap}", "", 6, {-0.4, 2.});
  StrVar4Hist var_PtV2Jpsi("PtV2Jpsi", "p_{T}", "GeV/c", strAny_ptV2.fNbins,
                           {0., 1.});
  cout << var_PtV2Jpsi.fNbins << endl;

  MIndexHist indexHistMass(var_MassJpsiCandidate, 1, n_rebin_mass);
  MIndexHist indexHistPtJpsiCandidate(var_PtJpsiCandidate, 1, 1);
  MIndexHist indexHistDeltaPhiUS(var_DeltaPhiUS, 1, 1);
  MIndexHist indexHistDeltaEtaUS(var_DeltaEtaUS, 1, 2);
  MIndexHist indexHistEtaGap(var_EtaGap, 1, 1);
  MIndexHist indexHistPtV2Jpsi(var_PtV2Jpsi, 1, 1);
  MIndexAny indexAnyPtV2Jpsi(strAny_ptV2, 1);

  MHGroupTool<TF1> g_tf1_input(file_input_tf1, "fitted_signal_poly6_pt%d",
                               {var_PtV2Jpsi}, {1});

  auto h_mass =
      (THnD *)file_input_mass->Get("fPosZ_MassUS_PtUS_NumContribCalib");
  MHnTool hnTool_mass(h_mass);
  hnTool_mass.Rebin(0, 200);
  hnTool_mass.Rebin(3, 5);
  hnTool_mass.PrintAllAxis();

  TH2D *mass_pt_lowMult = hnTool_mass.Project(1, 2, {0, 1});
  TH2D *mass_pt_highMult = hnTool_mass.Project(1, 2, {0, 2});
  MFitterPoly fitterPoly_mass(mass_pt_highMult->ProjectionY(), 1.88, 4.32);
  fitterPoly_mass.initializeBasis(6);

  MHGroupTool1D assoYeild_lowMult(
      file_input,
      "MassUS_AssoYeild_lowMult_DeltaEtaUS_%d_DeltaPhiUS_%d_ptV2_%d",
      {var_DeltaEtaUS, var_DeltaPhiUS, var_PtV2Jpsi}, {2, 1, 1});
  MHGroupTool1D assoYeild_highMult(
      file_input,
      "MassUS_AssoYeild_highMult_DeltaEtaUS_%d_DeltaPhiUS_%d_ptV2_%d",
      {var_DeltaEtaUS, var_DeltaPhiUS, var_PtV2Jpsi}, {2, 1, 1});
  MFitterPoly fitterPoly_asso(assoYeild_highMult.GetHist(vector<int>{1, 1, 1}),
                              2., 4.32);
  fitterPoly_asso.initializeBasis(6);

  // MHist1D h1_assoYeild_lowMult(indexHistDeltaPhiUS, "AssoYeild_lowMult");
  // MHist1D h1_assoYeild_highMult(indexHistDeltaPhiUS, "AssoYeild_highMult");
  // MHist1D h1_assoYeild_temp(indexHistDeltaPhiUS, "AssoYeild_sub");
  // MVec<MHist1D> assoYeild_lowMult_DeltaEta(indexHistDeltaEtaUS,
  //                                          h1_assoYeild_lowMult);
  // MVec<MHist1D> assoYeild_highMult_DeltaEta(indexHistDeltaEtaUS,
  //                                           h1_assoYeild_highMult);
  MHist1D h1_assoYeild_sub(indexHistDeltaPhiUS, "AssoYeild_sub");
  MVec<MHist1D> h1_assoYeild_sub_DeltaEta(indexHistDeltaEtaUS,
                                          h1_assoYeild_sub);

  gDirectory = file_output;
  TH1D *nsignal_highMult_ptV2 =
      new TH1D("nsignal_highMult_ptV2", "nsignal_highMult_ptV2",
               strAny_ptV2.fNbins, 0, strAny_ptV2.fNbins);
  TH1D *nsignal_lowMult_ptV2 =
      new TH1D("nsignal_lowMult_ptV2", "nsignal_lowMult_ptV2",
               strAny_ptV2.fNbins, 0, strAny_ptV2.fNbins);

  MVec<MVec<MHist1D>, MIndexAny<StrAny_ptV2>> vec_assoYeild_sub(
      indexAnyPtV2Jpsi, h1_assoYeild_sub_DeltaEta);
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
        auto assoYeild_lowMult_diff = assoYeild_lowMult.GetHist(
            vector<int>{i_deltaEta, i_deltaPhi, iPtV2});
        auto assoYeild_highMult_diff = assoYeild_highMult.GetHist(
            vector<int>{i_deltaEta, i_deltaPhi, iPtV2});

        fitterPoly_asso.setHisto(assoYeild_highMult_diff);
        fitterPoly_asso.fitWithSignal();
        double nyeild_highmult = fitterPoly_asso.fNSignal;

        fitterPoly_asso.setHisto(assoYeild_lowMult_diff);
        fitterPoly_asso.fitWithSignal();
        double nyeild_lowmult = fitterPoly_asso.fNSignal;

        double normalized_nyeild_highmult = nyeild_highmult / nsignal_highMult;
        double normalized_nyeild_lowmult = nyeild_lowmult / nsignal_lowMult;

        double normalized_nyeild_sub =
            normalized_nyeild_highmult - normalized_nyeild_lowmult;
        vec_assoYeild_sub.currentObject().SetBinInfo(normalized_nyeild_sub, 0);
      }
  }

  file_output->Write();
  file_output->Close();
}

int main(int argc, char **argv) {
  TString path_input = argc > 1 ? argv[1]
                                : "/home/szhu/work/alice/analysis/QA/test/"
                                  "AssoYeildGroupEtagap_NoScale.root";
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
                                "AssoYeildFit_noScale.pdf";

  AssoYeildFit_noScale(path_input, path_input_mass, path_input_tf1, path_output,
                       path_pdf);

  return 0;
}