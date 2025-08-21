#include "MHelper.h"
#include "MMath.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"

void AssoYeildGroupEtaGap(
    TString path_input =
        "/home/szhu/work/alice/analysis/QA/test/AssoYeildQA.root",
    TString path_output =
        "/home/szhu/work/alice/analysis/QA/test/AssoYeildGroupEtaGap_bInt.root",
    TString path_pdf = "/home/szhu/work/alice/analysis/QA/test/"
                       "AssoYeildGroupEtaGap_bInt.pdf",
    bool is_bInt = true) {
  TFile *file_input = new TFile(path_input);
  TFile *file_output = new TFile(path_output, "RECREATE");

  StrVar4Hist var_fPosZ("PosZUS", "#it{V}_{Z}", "cm", 8, {-10, 10});
  StrVar4Hist var_NumContribCalibBinned(
      "NumContribCalibUS", "N_{vtx contrib} Calibrated", "", 10,
      {0, 23, 31, 37, 43, 48, 54, 61, 69, 81, 297});
  StrVar4Hist var_MassJpsiCandidate("MassUS", "M_{ee}", "GeV^{2}/c^{4}", 90,
                                    {1.8, 5.4});
  StrVar4Hist var_PtJpsiCandidate("PtUS", "p_{T}", "GeV/c", 10, {0., 10.});
  StrVar4Hist var_DeltaEtaUS("DeltaEtaUS", "#Delta#eta_{J/#psi, track}", "", 80,
                             {-4., 4.});
  StrVar4Hist var_DeltaPhiUS("DeltaPhiUS", "#Delta#phi_{J/#psi, track}", "", 10,
                             {-M_PI_2, M_PI + M_PI_2});
  StrVar4Hist var_EtaGap("EtaGap", "#Delta#eta_{gap}", "", 6, {-0.4, 2.});

  MHGroupTool2D *hgroupTool2d_total_mass = new MHGroupTool2D(
      file_input, "h2_total_mass_%d", {var_MassJpsiCandidate}, {2});
  MHGroupTool2D *hgroupTool2d_lowMult_mass = new MHGroupTool2D(
      file_input, "h2_lowMult_mass_%d", {var_MassJpsiCandidate}, {2});
  MHGroupTool2D *hgroupTool2d_highMult_mass = new MHGroupTool2D(
      file_input, "h2_highMult_mass_%d", {var_MassJpsiCandidate}, {2});
  MHGroupTool2D *hgroupTool2d_highSubLow_mass = new MHGroupTool2D(
      file_input, "h2_highSubLow_mass_%d", {var_MassJpsiCandidate}, {2});

  gPublisherCanvas = new MPublisherCanvas(path_pdf, 3, 6, 600, 600);
  MRootGraphic::StyleCommon();

#define LGetBinContent(hist, vec_index, ...)                                   \
  hist->GetHist(vec_index)->GetBinContent(__VA_ARGS__)

#define LGetBinError(hist, vec_index, ...)                                     \
  hist->GetHist(vec_index)->GetBinError(__VA_ARGS__)

  MIndexHist indexHistMass(var_MassJpsiCandidate, 1, 2);
  MIndexHist indexHistEtaGap(var_EtaGap, 1, 1);
  file_output->cd();
  MHist2D h_b(indexHistMass, indexHistEtaGap, "b");
  MHist2D h_a0(indexHistMass, indexHistEtaGap, "a0");
  MHist2D h_a1(indexHistMass, indexHistEtaGap, "a1");
  MHist2D h_a2(indexHistMass, indexHistEtaGap, "a2");
  MHist2D h_a3(indexHistMass, indexHistEtaGap, "a3");
  MHist2D h_v22(indexHistMass, indexHistEtaGap, "v22");
  MHist2D h_v22part(indexHistMass, indexHistEtaGap, "v22part");
  MHist2D h_a0PlusB(indexHistMass, indexHistEtaGap, "a0PlusB");
  //   MHist1D h1d_b(indexHistEtaGap, "b");
  //   MVec<MHist1D> hVec_b(indexHistEtaGap, h1d_b);
  // quit file_output directory
  gDirectory = nullptr;

  auto ApplyEtaGap = [](TH2D *h2, double gap_eta) {
    h2->GetXaxis()->SetRangeUser(-1.8, -gap_eta / 2.);
    auto h1_highSubLow_mass1 =
        (TH1D *)h2->ProjectionY(Form("h1_highSubLow_mass1_%d", GenerateUID()));
    h2->GetXaxis()->SetRangeUser(gap_eta / 2., 1.8);
    auto h1_highSubLow_mass =
        (TH1D *)h2->ProjectionY(Form("h1_highSubLow_mass_%d", GenerateUID()));
    h1_highSubLow_mass->Add(h1_highSubLow_mass1);
    return h1_highSubLow_mass;
  };

  for (auto i_mass : indexHistMass) {
    for (auto i_etaGap : indexHistEtaGap) {
      double deltaEta = indexHistEtaGap.GetBinUpperEdge();
      auto h2_highSubLow_mass =
          hgroupTool2d_highSubLow_mass->GetHist(vector<int>{i_mass});
      auto h2_lowMult_mass =
          hgroupTool2d_lowMult_mass->GetHist(vector<int>{i_mass});
      auto h2_highMult_mass =
          hgroupTool2d_highMult_mass->GetHist(vector<int>{i_mass});

      auto h1_highSubLow_mass = ApplyEtaGap(h2_highSubLow_mass, deltaEta);
      auto h1_lowMult_mass = ApplyEtaGap(h2_lowMult_mass, deltaEta);
      auto h1_highMult_mass = ApplyEtaGap(h2_highMult_mass, deltaEta);
      /////////////////////////////////////////////////
      gPublisherCanvas->NewPad()->cd();
      TF1 f1_modu("f1_modu", "[0]+2*([1]*cos(x)+[2]*cos(2*x)+[3]*cos(3*x))",
                  -M_PI_2, M_PI + M_PI_2);
      h1_highSubLow_mass->Fit(&f1_modu, "Q", "", -M_PI_2, M_PI + M_PI_2);
      double results_modu[4];
      f1_modu.GetParameters(results_modu);

      double b_lowMult_mass = h1_lowMult_mass->GetBinContent(1);
      h1_highSubLow_mass->SetTitle(Form("i_mass = %d, #Delta#eta_{gap} = %.2f",
                                        i_mass,
                                        indexHistEtaGap.GetBinUpperEdge()));
      h1_highSubLow_mass->Draw();
      TF1 f1_modu_0("f1_modu_0", "[0]", -M_PI_2, M_PI + M_PI_2);
      TF1 f1_modu_1("f1_modu_1", "[0]+2*[1]*cos(x)", -M_PI_2, M_PI + M_PI_2);
      TF1 f1_modu_2("f1_modu_2", "[0]+2*[1]*cos(2*x)", -M_PI_2, M_PI + M_PI_2);
      TF1 f1_modu_3("f1_modu_3", "[0]+2*[1]*cos(3*x)", -M_PI_2, M_PI + M_PI_2);
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

#define FillHist(name, ...)                                                    \
  MDouble name##Value __VA_ARGS__;                                             \
  h_##name.SetBinInfo(name##Value);

      double b_value, b_error;
      if (is_bInt) {
        b_value = h1_highSubLow_mass->IntegralAndError(
            1, h1_highSubLow_mass->GetNbinsX(), b_error);
        b_value /= (double)h1_highSubLow_mass->GetNbinsX();
        b_error /= (double)h1_highSubLow_mass->GetNbinsX();
      } else {
        b_value = h1_highSubLow_mass->GetBinContent(1);
        b_error = h1_highSubLow_mass->GetBinError(1);
      }
      FillHist(b, (b_value, b_error));

      FillHist(a0, (results_modu[0], f1_modu.GetParError(0)));
      FillHist(a1, (results_modu[1], f1_modu.GetParError(1)));
      FillHist(a2, (results_modu[2], f1_modu.GetParError(2)));
      FillHist(a3, (results_modu[3], f1_modu.GetParError(3)));
      FillHist(v22, = a2Value / (a0Value + bValue));
      FillHist(v22part, = a2Value / a0Value);
      FillHist(a0PlusB, = a0Value + bValue);

      gPublisherCanvas->Draw(h1_highMult_mass)->Draw(h1_lowMult_mass);
    }
  }

  gPublisherCanvas->SetCanvasNwNh(1, 1);
  gPublisherCanvas->Draw(h_v22, "surf2")
      /* ->Draw(h_a2)
      ->Draw(h_b)
      ->Draw(h_a0)
      ->Draw(h_a0PlusB)
      ->Draw(h_v22part)
      ->Draw(h_v22) */
      ;

  gPublisherCanvas->SetCanvasNwNh(3, 2);
  for (int i_etaGap = 1; i_etaGap <= h_v22.fHisto->GetNbinsY(); i_etaGap++) {
    auto h1_v22 = h_v22.fHisto->ProjectionX(Form("h1_v22_%d", i_etaGap),
                                            i_etaGap, i_etaGap);
    h1_v22->GetYaxis()->SetTitle("V_{2}");
    h1_v22->SetTitle(Form("V_{2} vs M_{ee} for #Delta#eta_{gap} = %.2f",
                          h_v22.fHisto->GetYaxis()->GetBinUpEdge(i_etaGap)));
    MRootGraphic::StyleHistCommonHist(h1_v22);
    gPublisherCanvas->Draw(h1_v22);
  }

  gPublisherCanvas->finalize();

  file_output->Write();
  file_output->Close();
}

int main(int argc, char **argv) {
  gROOT->SetBatch(kTRUE);
  TString path_input = "/home/szhu/work/alice/analysis/QA/output/event_jpsi/"
                       "AssoYeildQA_LHC22pass4.root";
  TString path_output = "/home/szhu/work/alice/analysis/QA/output/event_jpsi/"
                        "AssoYeildGroupQAEtaGap.root";
  TString path_pdf = "/home/szhu/work/alice/analysis/QA/plot/event_jpsi/"
                     "AssoYeildGroupQAEtaGap.pdf";
  bool is_bInt = false;

  if (argc > 1) {
    path_input = argv[1];
  }
  if (argc > 2) {
    path_output = argv[2];
  }
  if (argc > 3) {
    path_pdf = argv[3];
  }
  if (argc > 3) {
    is_bInt = atoi(argv[4]);
  }

  AssoYeildGroupEtaGap(path_input, path_output, path_pdf, is_bInt);
  return 0;
}