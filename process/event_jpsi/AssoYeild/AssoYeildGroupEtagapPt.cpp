#include "MHelper.h"
#include "MMath.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"

void AssoYeildGroupEtagapPt(
    TString path_input = "/home/szhu/work/alice/analysis/QA/test/"
                         "AssoYeildPt.root",
    TString path_output = "/home/szhu/work/alice/analysis/QA/test/"
                          "AssoYeildGroupEtagapPt_bLow.root",
    TString path_pdf = "/home/szhu/work/alice/analysis/QA/test/"
                       "AssoYeildGroupEtagapPt_bInt.pdf",
    bool is_bInt = false) {
  TFile *file_input = new TFile(path_input);
  TFile *file_output = new TFile(path_output, "RECREATE");

  struct StrAny_ptV2 {
    const vector<vector<int>> bins = {{1},
                                      {2},
                                      {3},
                                      {4},
                                      {5},
                                      {1, 2},
                                      {2, 3},
                                      {3, 4},
                                      {4, 5},
                                      {5, 6},
                                      {1, 2, 3},
                                      {2, 3, 4},
                                      {3, 4, 5},
                                      {4, 5, 6},
                                      {1, 2, 3, 4},
                                      {2, 3, 4, 5},
                                      {3, 4, 5, 6},
                                      {1, 2, 3, 4, 5},
                                      {6, 7, 8, 9, 10},
                                      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}};
    const int fNbins = bins.size();
    const TString fName = "ptV2";

    vector<int> operator[](int index) { return bins[index]; }
  } strAny_ptV2;

  StrVar4Hist var_fPosZ("PosZUS", "#it{V}_{Z}", "cm", 8, {-10, 10});
  StrVar4Hist var_NumContribCalibBinned(
      "NumContribCalibUS", "N_{vtx contrib} Calibrated", "", 10,
      {0,7,12,17,22,29,37,46,57,73,300});
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
  // StrVar4Hist var_PtV2Jpsi("PtV2Jpsi", "p_{T}", "GeV/c", 3,
  // {1., 3., 5., 10.});

  MHGroupTool2D *hgroupTool2d_total_mass =
      new MHGroupTool2D(file_input, "h2_total_mass_pt_%d_%d",
                        {var_MassJpsiCandidate, var_PtV2Jpsi}, {2, 1});
  MHGroupTool2D *hgroupTool2d_lowMult_mass =
      new MHGroupTool2D(file_input, "h2_lowMult_mass_pt_%d_%d",
                        {var_MassJpsiCandidate, var_PtV2Jpsi}, {2, 1});
  MHGroupTool2D *hgroupTool2d_highMult_mass =
      new MHGroupTool2D(file_input, "h2_highMult_mass_pt_%d_%d",
                        {var_MassJpsiCandidate, var_PtV2Jpsi}, {2, 1});
  MHGroupTool2D *hgroupTool2d_highSubLow_mass =
      new MHGroupTool2D(file_input, "h2_highSubLow_mass_pt_%d_%d",
                        {var_MassJpsiCandidate, var_PtV2Jpsi}, {2, 1});

  gPublisherCanvas = new MPublisherCanvas(path_pdf, 3, 3, 600, 600);
  MRootGraphic::StyleCommon();

#define LGetBinContent(hist, vec_index, ...)                                   \
  hist->GetHist(vec_index)->GetBinContent(__VA_ARGS__)

#define LGetBinError(hist, vec_index, ...)                                     \
  hist->GetHist(vec_index)->GetBinError(__VA_ARGS__)

  MIndexHist indexHistMass(var_MassJpsiCandidate, 1, 2);
  MIndexHist indexHistEtaGap(var_EtaGap, 1, 1);
  MIndexHist indexHistPtJpsiCandidate(var_PtJpsiCandidate, 1, 1);
  // MIndexHist indexHistPtV2Jpsi(var_PtV2Jpsi, 1, 1);
  MIndexAny indexAnyPtV2Jpsi(strAny_ptV2, 1);
  // using MIndexAny1 = MIndexAny<StrAny_ptV2>;
  gDirectory = nullptr;
  MHist2D h_b(indexHistMass, indexHistEtaGap, "b");
  MHist2D h_a0(indexHistMass, indexHistEtaGap, "a0");
  MHist2D h_a1(indexHistMass, indexHistEtaGap, "a1");
  MHist2D h_a2(indexHistMass, indexHistEtaGap, "a2");
  MHist2D h_a3(indexHistMass, indexHistEtaGap, "a3");
  MHist2D h_v22(indexHistMass, indexHistEtaGap, "v22");
  MHist2D h_v22part(indexHistMass, indexHistEtaGap, "v22part");
  MHist2D h_a0PlusB(indexHistMass, indexHistEtaGap, "a0PlusB");

  // MHist3D h_b(indexHistPtV2Jpsi, indexHistMass, indexHistEtaGap, "b");
  // MHist3D h_a0(indexHistPtV2Jpsi, indexHistMass, indexHistEtaGap, "a0");
  // MHist3D h_a1(indexHistPtV2Jpsi, indexHistMass, indexHistEtaGap, "a1");
  // MHist3D h_a2(indexHistPtV2Jpsi, indexHistMass, indexHistEtaGap, "a2");
  // MHist3D h_a3(indexHistPtV2Jpsi, indexHistMass, indexHistEtaGap, "a3");
  // MHist3D h_v22(indexHistPtV2Jpsi, indexHistMass, indexHistEtaGap, "v22");
  // MHist3D h_v22part(indexHistPtV2Jpsi, indexHistMass, indexHistEtaGap,
  //                   "v22part");
  // MHist3D h_a0PlusB(indexHistPtV2Jpsi, indexHistMass, indexHistEtaGap,
  //                   "a0PlusB");

  file_output->cd();
  using MVec1 = MVec<MHist2D, MIndexAny<StrAny_ptV2>>;
  MVec1 hVec_b(indexAnyPtV2Jpsi, h_b);
  MVec1 hVec_a0(indexAnyPtV2Jpsi, h_a0);
  MVec1 hVec_a1(indexAnyPtV2Jpsi, h_a1);
  MVec1 hVec_a2(indexAnyPtV2Jpsi, h_a2);
  MVec1 hVec_a3(indexAnyPtV2Jpsi, h_a3);
  MVec1 hVec_v22(indexAnyPtV2Jpsi, h_v22);
  MVec1 hVec_v22part(indexAnyPtV2Jpsi, h_v22part);
  MVec1 hVec_a0PlusB(indexAnyPtV2Jpsi, h_a0PlusB);

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
    for (auto i_ptV2 : indexAnyPtV2Jpsi) {
      for (auto i_etaGap : indexHistEtaGap) {
        double deltaEta = indexHistEtaGap.GetBinUpperEdge();
        auto h2_highSubLow_mass =
            hgroupTool2d_highSubLow_mass->GetHist(vector<int>{i_mass, i_ptV2});
        auto h2_lowMult_mass =
            hgroupTool2d_lowMult_mass->GetHist(vector<int>{i_mass, i_ptV2});
        auto h2_highMult_mass =
            hgroupTool2d_highMult_mass->GetHist(vector<int>{i_mass, i_ptV2});

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
        h1_highSubLow_mass->SetTitle(
            Form("i_mass = %d, #Delta#eta_{gap} = %.2f", i_mass,
                 indexHistEtaGap.GetBinUpperEdge()));
        h1_highSubLow_mass->Draw();
        TF1 f1_modu_0("f1_modu_0", "[0]", -M_PI_2, M_PI + M_PI_2);
        TF1 f1_modu_1("f1_modu_1", "[0]+2*[1]*cos(x)", -M_PI_2, M_PI + M_PI_2);
        TF1 f1_modu_2("f1_modu_2", "[0]+2*[1]*cos(2*x)", -M_PI_2,
                      M_PI + M_PI_2);
        TF1 f1_modu_3("f1_modu_3", "[0]+2*[1]*cos(3*x)", -M_PI_2,
                      M_PI + M_PI_2);
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
  hVec_##name.current().SetBinInfo(name##Value);

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
  }

  gPublisherCanvas->SetCanvasNwNh(4, 3);
  for (auto i_ptv2 : indexAnyPtV2Jpsi) {
    auto h2_v22 = (TH2D *)hVec_v22.current();
    auto h2_b = (TH2D *)hVec_b.current();
    auto h2_a0 = (TH2D *)hVec_a0.current();
    auto h2_a2 = (TH2D *)hVec_a2.current();
    for (int i_etaGap = 1; i_etaGap <= h2_v22->GetNbinsY(); i_etaGap++) {
      auto h1_v22 = h2_v22->ProjectionX(Form("h1_v22_%d", GenerateUID()),
                                        i_etaGap, i_etaGap);
      h1_v22->GetYaxis()->SetTitle("V_{2}");
      h1_v22->SetTitle(
          Form("V_{2} vs M_{ee} for #Delta#eta_{gap} = %.2f, p_{T} bin: %d",
               h_v22.fHisto->GetYaxis()->GetBinUpEdge(i_etaGap), i_ptv2));
      MRootGraphic::StyleHistCommonHist(h1_v22);
      gPublisherCanvas->DrawClone(h1_v22);
      auto h1_b =
          h2_b->ProjectionX(Form("h1_b_%d", GenerateUID()), i_etaGap, i_etaGap);
      h1_b->GetYaxis()->SetTitle("b");
      h1_b->SetTitle(
          Form("b vs M_{ee} for #Delta#eta_{gap} = %.2f, p_{T} bin: %d",
               h_v22.fHisto->GetYaxis()->GetBinUpEdge(i_etaGap), i_ptv2));
      MRootGraphic::StyleHistCommonHist(h1_b);
      gPublisherCanvas->DrawClone(h1_b);
      auto h1_a0 = h2_a0->ProjectionX(Form("h1_a0_%d", GenerateUID()), i_etaGap,
                                      i_etaGap);
      h1_a0->GetYaxis()->SetTitle("a0");
      h1_a0->SetTitle(
          Form("a0 vs M_{ee} for #Delta#eta_{gap} = %.2f, p_{T} bin: %d",
               h_v22.fHisto->GetYaxis()->GetBinUpEdge(i_etaGap), i_ptv2));
      MRootGraphic::StyleHistCommonHist(h1_a0);
      gPublisherCanvas->DrawClone(h1_a0);
      auto h1_a2 = h2_a2->ProjectionX(Form("h1_a2_%d", GenerateUID()), i_etaGap,
                                      i_etaGap);
      h1_a2->GetYaxis()->SetTitle("a2");
      h1_a2->SetTitle(
          Form("a2 vs M_{ee} for #Delta#eta_{gap} = %.2f, p_{T} bin: %d",
               h_v22.fHisto->GetYaxis()->GetBinUpEdge(i_etaGap), i_ptv2));
      MRootGraphic::StyleHistCommonHist(h1_a2);
      gPublisherCanvas->DrawClone(h1_a2);
    }
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
  if (argc > 4) {
    is_bInt = atoi(argv[4]);
  }

  AssoYeildGroupEtagapPt(path_input, path_output, path_pdf, is_bInt);
  return 0;
}