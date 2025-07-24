#include "MHelper.h"
#include "MMath.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"

void AssoYeildGroupQA(
    TString path_input = "/home/szhu/work/alice/analysis/QA/output/event_jpsi/"
                         "AssoYeildQA_LHC22pass4.root",
    TString path_output = "/home/szhu/work/alice/analysis/QA/output/event_jpsi/"
                          "AssoYeildGroupQA.root") {
  TFile *file_input = new TFile(path_input);
  TFile *file_output = new TFile(path_output, "RECREATE");

  StrVar4Hist var_fPosZ("PosZUS", "#it{V}_{Z}", "cm", 8, {-10, 10});
  StrVar4Hist var_NumContribCalibBinned(
      "NumContribCalibUS", "N_{vtx contrib} Calibrated", "", 10,
      {0, 23, 31, 37, 43, 48, 54, 61, 69, 81, 297});
  StrVar4Hist var_MassJpsiCandidate("MassUS", "M_{ee}", "GeV^{2}/c^{4}", 90,
                                    {1.4, 5.});
  StrVar4Hist var_PtJpsiCandidate("PtUS", "p_{T}", "GeV/c", 10, {0., 10.});
  StrVar4Hist var_DeltaEtaUS("DeltaEtaUS", "#Delta#eta_{J/#psi, track}", "", 80,
                             {-4., 4.});
  StrVar4Hist var_DeltaPhiUS("DeltaPhiUS", "#Delta#phi_{J/#psi, track}", "", 10,
                             {-M_PI_2, M_PI + M_PI_2});

  MHGroupTool2D *hgroupTool2d_total_mass = new MHGroupTool2D(
      file_input, "h2_total_mass_%d", {var_MassJpsiCandidate}, {1});
  MHGroupTool2D *hgroupTool2d_lowMult_mass = new MHGroupTool2D(
      file_input, "h2_lowMult_mass_%d", {var_MassJpsiCandidate}, {1});
  MHGroupTool2D *hgroupTool2d_highMult_mass = new MHGroupTool2D(
      file_input, "h2_highMult_mass_%d", {var_MassJpsiCandidate}, {1});
  MHGroupTool2D *hgroupTool2d_highSubLow_mass = new MHGroupTool2D(
      file_input, "h2_highSubLow_mass_%d", {var_MassJpsiCandidate}, {1});

  gPublisherCanvas = new MPublisherCanvas(
      "/home/szhu/work/alice/analysis/QA/plot/event_jpsi/AssoYeildQA.pdf", 2, 3,
      600, 600);
  MRootGraphic::StyleCommon();

#define LGetBinContent(hist, vec_index, ...)                                   \
  hist->GetHist(vec_index)->GetBinContent(__VA_ARGS__)

#define LGetBinError(hist, vec_index, ...)                                     \
  hist->GetHist(vec_index)->GetBinError(__VA_ARGS__)

  MIndexHist indexHistMass(var_MassJpsiCandidate);
  MHist1D h_b(indexHistMass, "b");
  MHist1D h_a0(indexHistMass, "a0");
  MHist1D h_a1(indexHistMass, "a1");
  MHist1D h_a2(indexHistMass, "a2");
  MHist1D h_a3(indexHistMass, "a3");
  MHist1D h_v22(indexHistMass, "v22");

  for (auto i_mass : indexHistMass) {
    auto h2_highSubLow_mass =
        hgroupTool2d_highSubLow_mass->GetHist(vector<int>{i_mass});

    h2_highSubLow_mass->GetXaxis()->SetRangeUser(-2, 2);
    StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(),
                                 h2_highSubLow_mass);

    auto h1_highSubLow_mass = (TH1D *)h2_highSubLow_mass->ProjectionY(
        Form("h1_highSubLow_mass_%d", i_mass));
    gPublisherCanvas->NewPad()->cd();
    TF1 f1_modu("f1_modu", "[0]+[1]*cos(x)+[2]*cos(2*x)+[3]*cos(3*x)", -M_PI_2,
                M_PI + M_PI_2);
    h1_highSubLow_mass->Fit(&f1_modu, "Q", "", -M_PI_2, M_PI + M_PI_2);
    double results_modu[4];
    f1_modu.GetParameters(results_modu);

    auto h2_lowMult_mass =
        hgroupTool2d_lowMult_mass->GetHist(vector<int>{i_mass});
    auto h1_lowMult_mass = (TH1D *)h2_lowMult_mass->ProjectionY(
        Form("h1_lowMult_mass_%d", i_mass));
    double b_lowMult_mass = h1_lowMult_mass->GetBinContent(1);
    h1_highSubLow_mass->Draw();
    TF1 f1_modu_0("f1_modu_0", "[0]", -M_PI_2, M_PI + M_PI_2);
    TF1 f1_modu_1("f1_modu_1", "[0]+[1]*cos(x)", -M_PI_2, M_PI + M_PI_2);
    TF1 f1_modu_2("f1_modu_2", "[0]+[1]*cos(2*x)", -M_PI_2, M_PI + M_PI_2);
    TF1 f1_modu_3("f1_modu_3", "[0]+[1]*cos(3*x)", -M_PI_2, M_PI + M_PI_2);
    f1_modu_0.SetParameter(0, results_modu[0]);
    f1_modu_1.SetParameter(0, results_modu[0]);
    f1_modu_1.SetParameter(1, results_modu[1]);
    f1_modu_2.SetParameter(0, results_modu[0]);
    f1_modu_2.SetParameter(1, results_modu[2]);
    f1_modu_3.SetParameter(0, results_modu[0]);
    f1_modu_3.SetParameter(1, results_modu[3]);
    f1_modu_0.SetLineColor(kBlue);
    f1_modu_1.SetLineColor(kGreen + 2);
    f1_modu_2.SetLineColor(kYellow + 2);
    f1_modu_3.SetLineColor(kOrange + 2);
    f1_modu_0.DrawClone("same");
    f1_modu_1.DrawClone("same");
    f1_modu_2.DrawClone("same");

    MDouble bValue(h1_lowMult_mass->GetBinContent(1),
                   h1_lowMult_mass->GetBinError(1));
    MDouble a0Value(results_modu[0], f1_modu_0.GetParError(0));
    MDouble a1Value(results_modu[1], f1_modu_1.GetParError(1));
    MDouble a2Value(results_modu[2], f1_modu_2.GetParError(1));
    MDouble a3Value(results_modu[3], f1_modu_3.GetParError(1));
    MDouble v22Value = a2Value / (a0Value + bValue);

    h_b.SetBinInfo(bValue);
    h_a0.SetBinInfo(a0Value);
    h_a1.SetBinInfo(a1Value);
    h_a2.SetBinInfo(a2Value);
    h_a3.SetBinInfo(a3Value);
    h_v22.SetBinInfo(v22Value);
  }

  gPublisherCanvas->finalize();
  file_output->cd();
  h_b.Write();
  h_a0.Write();
  h_a1.Write();
  h_a2.Write();
  h_a3.Write();
  h_v22.Write();
  file_output->Close();
}