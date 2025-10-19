#include "MFit.h"
#include "MHelper.h"
#include "MMath.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"

funcWithJson(void, AssoYeildProj_noScale)(
    TString path_input = "/home/szhu/work/alice/analysis/QA/test/"
                         "AssoYeildGroupEtagap_NoScale.root",
    TString path_input_mass = "/home/szhu/work/alice/analysis/QA/input/jpsi/"
                              "JpsiMass_LHC24_apass1_DiElectron.root",
    TString path_output = "/home/szhu/work/alice/analysis/QA/test/"
                          "AssoYeildFit_noScale.root",
    TString path_pdf = "/home/szhu/work/alice/analysis/QA/test/"
                       "AssoYeildFit_noScale.pdf") {
  gErrorIgnoreLevel = kWarning;
  SetUpJson();
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  Configurable<int> config_n_rebin_mass("n_rebin_mass", 3);
  int n_rebin_mass = config_n_rebin_mass.data;

  TFile *file_input = new TFile(path_input);
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

  MIndexHist indexHistMass(var_MassJpsiCandidate, 1, n_rebin_mass);
  MIndexHist indexHistPtJpsiCandidate(var_PtJpsiCandidate, 1, 1);
  MIndexHist indexHistDeltaPhiUS(var_DeltaPhiUS, 1, 1);
  MIndexHist indexHistDeltaEtaUS(var_DeltaEtaUS, 1, 2);
  MIndexHist indexHistEtaGap(var_EtaGap, 1, 1);
  MIndexHist indexHistPtV2Jpsi(var_PtV2Jpsi, 1, 1);
  MIndexAny indexAnyPtV2Jpsi(strAny_ptV2, 1);

  MHGroupTool3D *hg3_assoYeild_highMult = new MHGroupTool3D(
      file_input, "MassUS_DeltaEtaUS_DeltaPhiUS_AssoYeild_highMult_ptV2_%d",
      {var_PtV2Jpsi}, {1});
  MHGroupTool3D *hg3_assoYeild_lowMult = new MHGroupTool3D(
      file_input, "MassUS_DeltaEtaUS_DeltaPhiUS_AssoYeild_lowMult_ptV2_%d",
      {var_PtV2Jpsi}, {1});

  auto h_mass =
      (THnD *)file_input_mass->Get("fPosZ_MassUS_PtUS_NumContribCalib");
  MHnTool hnTool_mass(h_mass);
  hnTool_mass.Rebin(0, 200);
  hnTool_mass.Rebin(3, 5);
  hnTool_mass.PrintAllAxis();

  TH2D *mass_pt_lowMult = hnTool_mass.Project(1, 2, {0, 1});
  TH2D *mass_pt_highMult = hnTool_mass.Project(1, 2, {0, 2});

  for (auto iPtV2 : indexAnyPtV2Jpsi) {
    auto bins_pt = indexAnyPtV2Jpsi[iPtV2 - 1];

    auto mass_highMult = mass_pt_highMult->ProjectionY(
        Form("mass_highMult_ptV2_%d", iPtV2), bins_pt.front(), bins_pt.back());
    auto mass_lowMult = mass_pt_lowMult->ProjectionY(
        Form("mass_lowMult_ptV2_%d", iPtV2), bins_pt.front(), bins_pt.back());

    file_output->WriteObject(mass_highMult,
                             Form("mass_highMult_ptV2_%d", iPtV2));
    file_output->WriteObject(mass_lowMult, Form("mass_lowMult_ptV2_%d", iPtV2));

    auto assoYeild_highMult =
        hg3_assoYeild_highMult->GetHist(vector<int>{iPtV2});
    auto assoYeild_lowMult = hg3_assoYeild_lowMult->GetHist(vector<int>{iPtV2});

    file_output->WriteObject(assoYeild_highMult,
                             Form("assoYeild_highMult_ptV2_%d", iPtV2));
    file_output->WriteObject(assoYeild_lowMult,
                             Form("assoYeild_lowMult_ptV2_%d", iPtV2));

    auto assoYeild_highMult_proj = assoYeild_highMult->ProjectionX(
        Form("assoYeild_highMult_proj_ptV2_%d", iPtV2));
    auto assoYeild_lowMult_proj = assoYeild_lowMult->ProjectionX(
        Form("assoYeild_lowMult_proj_ptV2_%d", iPtV2));

    file_output->WriteObject(assoYeild_highMult_proj,
                             Form("assoYeild_highMult_proj_ptV2_%d", iPtV2));
    file_output->WriteObject(assoYeild_lowMult_proj,
                             Form("assoYeild_lowMult_proj_ptV2_%d", iPtV2));

    for (auto i_deltaEta : indexHistDeltaEtaUS) {
      for (auto i_deltaPhi : indexHistDeltaPhiUS) {
        auto assoYeild_diff_highMult = assoYeild_highMult->ProjectionX(
            Form("assoYeild_diff_highMult_ptV2_%d_dEta_%d_dPhi_%d", iPtV2,
                 i_deltaEta, i_deltaPhi),
            i_deltaEta, i_deltaEta, i_deltaPhi, i_deltaPhi);
        auto assoYeild_diff_lowMult = assoYeild_lowMult->ProjectionX(
            Form("assoYeild_diff_lowMult_ptV2_%d_dEta_%d_dPhi_%d", iPtV2,
                 i_deltaEta, i_deltaPhi),
            i_deltaEta, i_deltaEta, i_deltaPhi, i_deltaPhi);
        file_output->WriteObject(
            assoYeild_diff_highMult,
            Form("assoYeild_diff_highMult_ptV2_%d_dEta_%d_dPhi_%d", iPtV2,
                 i_deltaEta, i_deltaPhi));
        file_output->WriteObject(
            assoYeild_diff_lowMult,
            Form("assoYeild_diff_lowMult_ptV2_%d_dEta_%d_dPhi_%d", iPtV2,
                 i_deltaEta, i_deltaPhi));
      }
    }
  }

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
  TString path_output = argc > 3 ? argv[3]
                                 : "/home/szhu/work/alice/analysis/QA/test/"
                                   "AssoYeilFit_noScale.root";
  TString path_pdf = argc > 4 ? argv[4]
                              : "/home/szhu/work/alice/analysis/QA/test/"
                                "AssoYeildFit_noScale.pdf";

  AssoYeildProj_noScale(path_input, path_input_mass, path_output, path_pdf);

  return 0;
}