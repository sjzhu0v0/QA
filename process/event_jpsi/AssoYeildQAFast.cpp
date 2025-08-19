#include "MHelper.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"
#include "TLegend.h"

void fitPad(TPad *pad, TH1D *h1) {
  pad->cd();
  TF1 f1_modu("f1_modu", "[0]+2*[1]*cos(x)+2*[2]*cos(2*x)+2*[3]*cos(3*x)", -M_PI_2,
              M_PI + M_PI_2);
  h1->Fit(&f1_modu, "Q", "", -M_PI_2, M_PI + M_PI_2);
  double results_modu[4];
  f1_modu.GetParameters(results_modu);
  h1->Draw("E1");
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
  f1_modu_1.DrawClone("same");
  f1_modu_2.DrawClone("same");
  f1_modu_3.DrawClone("same");

  TLegend *legend = new TLegend(0.6, 0.62, 0.93, 0.88);
	legend->SetLineColor(0);
  legend->AddEntry(f1_modu.Clone(), "a_{0}+2#Sigma^{3}_{i=1}a_{i}cos(i*#Delta#phi)",
                   "l");
  legend->AddEntry(f1_modu_1.Clone(), "cos(1*#Delta#phi) term", "l");
  legend->AddEntry(f1_modu_2.Clone(), "cos(2*#Delta#phi) term", "l");
  legend->AddEntry(f1_modu_3.Clone(), "cos(3*#Delta#phi) term", "l");

  legend->Draw("same");
};

void AssoYeildQAFast(
    TString input_se_pr =
        "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
        "JpsiAsso_cluster1_LHC22pass4_dqfilter.root:DeltaEtaUS_"
        "DeltaPhiUS_PosZUS_MassUS_PtUS_NumContribCalibUS",
    TString input_se_raw = "/home/szhu/work/alice/analysis/QA/input/jpsi/"
                           "JpsiQA_LHC22pass4_dqfilter.root:fPosZ_"
                           "MassUS_PtUS_NumContribCalib",
    TString input_me_pr = "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
                          "MixEventReading_LHC22pass4_dqfilter.root:DeltaEtaUS_"
                          "DeltaPhiUS_PosZUS_MassUS_PtUS_NumContribCalibUS",
    TString path_output = "/home/szhu/work/alice/analysis/QA/output/event_jpsi/"
                          "AssoYeildQAFast") {
  auto hist_se_pr = MRootIO::GetObjectDiectly<THnD>(input_se_pr);
  auto hist_se_raw = MRootIO::GetObjectDiectly<THnD>(input_se_raw);
  auto hist_me_pr = MRootIO::GetObjectDiectly<THnD>(input_me_pr);

  //   TFile *file_output =
  //       new TFile(Form("%s_LHC22pass4.root", path_output.Data()),
  //       "RECREATE");

  MHnTool hnTool_se_pr(hist_se_pr);
  MHnTool hnTool_se_raw(hist_se_raw);
  hnTool_se_raw.Rebin(0, 25); // Rebin DeltaPhiUS
  MHnTool hnTool_me_pr(hist_me_pr);

  AssocYeildHelper_v2 assoYeild(&hnTool_se_pr, &hnTool_me_pr, &hnTool_se_raw);
  assoYeild.Rebin(gtype_vars::kNumContrib, 5);
  assoYeild.Rebin(gtype_vars::kDeltaEta, 2);
  assoYeild.Rebin(gtype_vars::kMass, 2);

  hnTool_se_pr.PrintAllAxis();
  // Axis 0: axis0, title: #Delta#eta_{J/#psi, track}  nbins:80
  // Axis 1: axis1, title: #Delta#phi_{J/#psi, track}  nbins:10
  // Axis 2: axis2, title: #it{V}_{Z} (cm)  nbins:8
  // Axis 3: axis3, title: M_{ee} (GeV^{2}/c^{4})  nbins:90
  // Axis 4: axis4, title: p_{T} (GeV/c)  nbins:10
  // Axis 5: axis5, title: N_{vtx contrib} Calibrated  nbins:10

  hnTool_se_raw.PrintAllAxis();
  // Axis 0: axis0, title: #it{V}_{Z} (cm)  nbins:8
  // Axis 1: axis1, title: M_{ee} (GeV^{2}/c^{4})  nbins:100
  // Axis 2: axis2, title: p_{T} (GeV/c)  nbins:10
  // Axis 3: axis3, title: N_{vtx contrib} Calibrated  nbins:10

  hnTool_me_pr.PrintAllAxis();
  // Axis 0: axis0, title: #Delta#eta_{J/#psi, track}  nbins:80
  // Axis 1: axis1, title: #Delta#phi_{J/#psi, track}  nbins:10
  // Axis 2: axis2, title: #it{V}_{Z} (cm)  nbins:8
  // Axis 3: axis3, title: M_{ee} (GeV^{2}/c^{4})  nbins:90
  // Axis 4: axis4, title: p_{T} (GeV/c)  nbins:10
  // Axis 5: axis5, title: N_{vtx contrib} Calibrated  nbins:10

  gPublisherCanvas = new MPublisherCanvas(
      "/home/szhu/work/alice/analysis/QA/plot/event_jpsi/AssoYeildQAFast.pdf",
      3, 1, 600, 600);
  MRootGraphic::StyleCommon();
  gStyle->SetPalette(kRainBow);

  auto h2_se_pr = hnTool_se_pr.Project(1, 0, {0, 0, 0, 0});
  StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(), h2_se_pr);
  auto h2_me_pr = hnTool_me_pr.Project(1, 0, {0, 0, 0, 0});
  StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(), h2_me_pr);

  assoYeild.SetMixMultInt(false);

  auto h2_total = assoYeild.AssociatedYeild(0, 0, 0);

  StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(), h2_total);

  gPublisherCanvas->SetCanvasNwNh(2, 1);
  StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(), h2_total);
  auto h_total_proj = (TH1D *)h2_total->ProjectionY("h_total_proj");
  h_total_proj->GetYaxis()->SetTitle("Y(#Delta#phi)");
  MRootGraphic::StyleHistCommonHist(h_total_proj);

  gPublisherCanvas->NewPad()->cd();
  gPad->SetTopMargin(0.15);
  gPad->SetBottomMargin(0.1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  h_total_proj->Draw("E1");

  auto h2_lowMult = assoYeild.AssociatedYeild(0, 0, 1);
  auto h2_highMult = assoYeild.AssociatedYeild(0, 0, 2);
  auto h2_highSubLow = (TH2D *)h2_highMult->Clone("h2_highSubLow");
  HistSubstraction2D(h2_highSubLow, h2_highMult, h2_lowMult);

  h2_lowMult->SetTitle("Low Multiplicity");
  h2_highMult->SetTitle("High Multiplicity");
  h2_highSubLow->SetTitle("High Mult - Low Mult");

  gPublisherCanvas->SetCanvasNwNh(3, 1);
  //   auto pad_total = gPublisherCanvas->NewPad();
  //   StyleFlow::DeltaPhi_DeltaEta(pad_total, h2_total);
  auto pad_lowMult = gPublisherCanvas->NewPad();
  auto pad_highMult = gPublisherCanvas->NewPad();
  auto pad_highSubLow = gPublisherCanvas->NewPad();

  //   h2_total->GetZaxis()->SetRangeUser(3.e-3, 11.e-3);
  //   h2_lowMult->GetZaxis()->SetRangeUser(3.e-3, 11.e-3);
  //   h2_highMult->GetZaxis()->SetRangeUser(3.e-3, 11.e-3);
  //   h2_highSubLow->GetZaxis()->SetRangeUser(3.e-3, 11.e-3);

  StyleFlow::DeltaPhi_DeltaEta(pad_lowMult, h2_lowMult);
  StyleFlow::DeltaPhi_DeltaEta(pad_highMult, h2_highMult);
  StyleFlow::DeltaPhi_DeltaEta(pad_highSubLow, h2_highSubLow);

  gPublisherCanvas->SetCanvasNwNh(2, 1);
  h2_highSubLow->SetTitle("High Mult - Low Mult");
  StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(), h2_highSubLow);
  h2_highSubLow->GetXaxis()->SetRangeUser(-2, 2);
  h2_highSubLow->GetZaxis()->SetRangeUser(4.5e-3, 6e-3);
  auto h_highSubLow_proj =
      (TH1D *)h2_highSubLow->ProjectionY("h_highSubLow_proj");
  MRootGraphic::StyleHistCommonHist(h_highSubLow_proj);
  h_highSubLow_proj->GetYaxis()->SetTitle("Y(#Delta#phi)");
  auto pad_highSubLow_proj = gPublisherCanvas->NewPad();
  pad_highSubLow_proj->SetLeftMargin(0.15);
  fitPad(pad_highSubLow_proj, h_highSubLow_proj);
  h_highSubLow_proj->SetTitle("#eta_{gap}=0, |#Delta#eta_{gap}|<2");
  h_highSubLow_proj->GetYaxis()->SetTitleOffset(2.);

  auto h2_highSubLow_subEta1 =
      (TH2D *)h2_highSubLow->Clone("h2_highSubLow_subEta1");
  h2_highSubLow_subEta1->GetXaxis()->SetRangeUser(-2, -0.4);
  auto h2_highSubLow_subEta2 =
      (TH2D *)h2_highSubLow->Clone("h2_highSubLow_subEta2");
  h2_highSubLow_subEta2->GetXaxis()->SetRangeUser(0.4, 2);
  auto h2_highSubLow_subEta1_proj =
      (TH1D *)h2_highSubLow_subEta1->ProjectionY("h2_highSubLow_subEta1_proj");
  auto h2_highSubLow_subEta2_proj =
      (TH1D *)h2_highSubLow_subEta2->ProjectionY("h2_highSubLow_subEta2_proj");
  h2_highSubLow_subEta1_proj->Add(h2_highSubLow_subEta2_proj);

  // clear h2_highSubLow_subEta1 from -0.4 to 2
  for (int i = 1; i <= h2_highSubLow_subEta1->GetNbinsX(); i++) {
    if (h2_highSubLow_subEta1->GetXaxis()->GetBinCenter(i) > -0.4 &&
        h2_highSubLow_subEta1->GetXaxis()->GetBinCenter(i) < 0.4) {
      for (int j = 1; j <= h2_highSubLow_subEta1->GetNbinsY(); j++) {
        h2_highSubLow_subEta1->SetBinContent(i, j, 0);
        h2_highSubLow_subEta1->SetBinError(i, j, 0);
      }
    }
  }

  h2_highSubLow_subEta1->GetXaxis()->SetRangeUser(-2, 2);
  h2_highSubLow_subEta2->GetXaxis()->SetRangeUser(-2, 2);
  h2_highSubLow_subEta1->GetZaxis()->SetRangeUser(4.5e-3, 6e-3);
  h2_highSubLow_subEta2->GetZaxis()->SetRangeUser(4.5e-3, 6e-3);

  auto pad_subEta = gPublisherCanvas->NewPad();
  h2_highSubLow_subEta1->SetTitle("High Mult - Low Mult, #eta_{gap}=0.8");
  StyleFlow::DeltaPhi_DeltaEta(pad_subEta, h2_highSubLow_subEta1);

  auto pad_subEta_proj = gPublisherCanvas->NewPad();
  pad_subEta_proj->SetLeftMargin(0.15);
  fitPad(pad_subEta_proj, h2_highSubLow_subEta1_proj);
  h2_highSubLow_subEta1_proj->SetTitle("#eta_{gap}=0.8, |#Delta#eta_{gap}|<2");
  h2_highSubLow_subEta1_proj->GetYaxis()->SetTitleOffset(2.);

  // // gPublisherCanvas->SetCanvasNwNh(2, 2);
  //   file_output->cd();
  //   for (int i_mass = 1; i_mass <= hnTool_se_pr.GetNbins(3); i_mass++) {
  //     auto h2_total_mass = assoYeild.AssociatedYeild(i_mass, 0, 0);
  //     auto h2_lowMult_mass = assoYeild.AssociatedYeild(i_mass, 0, 1);
  //     auto h2_highMult_mass = assoYeild.AssociatedYeild(i_mass, 0, 2);
  //     auto h2_highSubLow_mass =
  //         (TH2D *)h2_highMult_mass->Clone(Form("h2_highSubLow_mass_%d",
  //         i_mass));

  //     h2_total_mass->SetName(Form("h2_total_mass_%d", i_mass));
  //     h2_lowMult_mass->SetName(Form("h2_lowMult_mass_%d", i_mass));
  //     h2_highMult_mass->SetName(Form("h2_highMult_mass_%d", i_mass));
  //     h2_highSubLow_mass->SetName(Form("h2_highSubLow_mass_%d", i_mass));
  //     h2_total_mass->Write();
  //     h2_lowMult_mass->Write();
  //     h2_highMult_mass->Write();
  //     h2_highSubLow_mass->Write();

  //     HistSubstraction2D(h2_highSubLow_mass, h2_highMult_mass,
  //     h2_lowMult_mass);
  //     StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(),
  //     h2_total_mass);
  //     StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(),
  //     h2_lowMult_mass);
  //     StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(),
  //     h2_highMult_mass);
  //     StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(),
  //                                  h2_highSubLow_mass);
  //     for (int i_pt = 1; i_pt <= hnTool_se_pr.GetNbins(4); i_pt++) {
  //       auto h2_total_mass_pt = assoYeild.AssociatedYeild(i_mass, i_pt, 0);
  //       auto h2_lowMult_mass_pt = assoYeild.AssociatedYeild(i_mass, i_pt, 1);
  //       auto h2_highMult_mass_pt = assoYeild.AssociatedYeild(i_mass, i_pt,
  //       2); auto h2_highSubLow_mass_pt = (TH2D *)h2_highMult_mass_pt->Clone(
  //           Form("h2_highSubLow_mass_pt_%d_%d", i_mass, i_pt));
  //       HistSubstraction2D(h2_highSubLow_mass_pt, h2_highMult_mass_pt,
  //                          h2_lowMult_mass_pt);
  //       StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(),
  //                                    h2_total_mass_pt);
  //       StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(),
  //                                    h2_lowMult_mass_pt);
  //       StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(),
  //                                    h2_highMult_mass_pt);
  //       StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(),
  //                                    h2_highSubLow_mass_pt);
  //       h2_total_mass_pt->SetName(Form("h2_total_mass_pt_%d_%d", i_mass,
  //       i_pt)); h2_lowMult_mass_pt->SetName(
  //           Form("h2_lowMult_mass_pt_%d_%d", i_mass, i_pt));
  //       h2_highMult_mass_pt->SetName(
  //           Form("h2_highMult_mass_pt_%d_%d", i_mass, i_pt));
  //       h2_highSubLow_mass_pt->SetName(
  //           Form("h2_highSubLow_mass_pt_%d_%d", i_mass, i_pt));
  //       h2_total_mass_pt->Write();
  //       h2_lowMult_mass_pt->Write();
  //       h2_highMult_mass_pt->Write();
  //       h2_highSubLow_mass_pt->Write();
  //     }
  //   }

  gPublisherCanvas->finalize();
  //   file_output->Close();
}