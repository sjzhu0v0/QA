#include "MHelper.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"

void AssoYeildQA(
    TString input_se_pr = "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
                          "JpsiAsso_LHC22pass4_dqfilter.root:DeltaEtaUS_"
                          "DeltaPhiUS_PosZUS_MassUS_PtUS_NumContribCalibUS",
    TString input_se_raw = "/home/szhu/work/alice/analysis/QA/input/jpsi/"
                           "JpsiQA_LHC22pass4_dqfilter.root:fPosZ_"
                           "MassUS_PtUS_NumContribCalib",
    TString input_me_pr = "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
                          "MixEventReading_LHC22pass4_dqfilter.root:DeltaEtaUS_"
                          "DeltaPhiUS_PosZUS_MassUS_PtUS_NumContribCalibUS",
    TString path_output = "/home/szhu/work/alice/analysis/QA/output/event_jpsi/"
                          "AssoYeildQA") {
  auto hist_se_pr = MRootIO::GetObjectDiectly<THnD>(input_se_pr);
  auto hist_se_raw = MRootIO::GetObjectDiectly<THnD>(input_se_raw);
  auto hist_me_pr = MRootIO::GetObjectDiectly<THnD>(input_me_pr);

  TFile *file_output =
      new TFile(Form("%s.root", path_output.Data()), "RECREATE");

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

  gPublisherCanvas = new MPublisherCanvas(path_output + ".pdf", 3, 1, 600, 600);
  MRootGraphic::StyleCommon();
  gStyle->SetPalette(kRainBow);

  auto h2_se_pr = hnTool_se_pr.Project(1, 0, {0, 0, 0, 0});
  StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(), h2_se_pr);
  auto h2_me_pr = hnTool_me_pr.Project(1, 0, {0, 0, 0, 0});
  StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(), h2_me_pr);

  assoYeild.SetMixMultInt(false);

  auto h2_total = assoYeild.AssociatedYeild(0, 0, 0);

  StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(), h2_total);

  auto h2_lowMult = assoYeild.AssociatedYeild(0, 0, 1);
  auto h2_highMult = assoYeild.AssociatedYeild(0, 0, 2);
  auto h2_highSubLow = (TH2D *)h2_highMult->Clone("h2_highSubLow");
  HistSubstraction2D(h2_highSubLow, h2_highMult, h2_lowMult);

  // h2_total->GetZaxis()->SetRangeUser(1., 4.5);
  // h2_lowMult->GetZaxis()->SetRangeUser(1., 4.5);
  // h2_highMult->GetZaxis()->SetRangeUser(1., 4.5);

  gPublisherCanvas->SetCanvasNwNh(2, 2);
  auto pad_total = gPublisherCanvas->NewPad();
  StyleFlow::DeltaPhi_DeltaEta(pad_total, h2_total);
  auto pad_lowMult = gPublisherCanvas->NewPad();
  StyleFlow::DeltaPhi_DeltaEta(pad_lowMult, h2_lowMult);
  auto pad_highMult = gPublisherCanvas->NewPad();
  StyleFlow::DeltaPhi_DeltaEta(pad_highMult, h2_highMult);
  auto pad_highSubLow = gPublisherCanvas->NewPad();
  StyleFlow::DeltaPhi_DeltaEta(pad_highSubLow, h2_highSubLow);

  // // gPublisherCanvas->SetCanvasNwNh(2, 2);
  file_output->cd();
  for (int i_mass = 1; i_mass <= hnTool_se_pr.GetNbins(3); i_mass++) {
    auto h2_total_mass = assoYeild.AssociatedYeild(i_mass, 0, 0);
    auto h2_lowMult_mass = assoYeild.AssociatedYeild(i_mass, 0, 1);
    auto h2_highMult_mass = assoYeild.AssociatedYeild(i_mass, 0, 2);
    auto h2_highSubLow_mass =
        (TH2D *)h2_highMult_mass->Clone(Form("h2_highSubLow_mass_%d", i_mass));

    h2_total_mass->SetName(Form("h2_total_mass_%d", i_mass));
    h2_lowMult_mass->SetName(Form("h2_lowMult_mass_%d", i_mass));
    h2_highMult_mass->SetName(Form("h2_highMult_mass_%d", i_mass));
    h2_highSubLow_mass->SetName(Form("h2_highSubLow_mass_%d", i_mass));
    h2_total_mass->Write();
    h2_lowMult_mass->Write();
    h2_highMult_mass->Write();
    h2_highSubLow_mass->Write();

    HistSubstraction2D(h2_highSubLow_mass, h2_highMult_mass, h2_lowMult_mass);
    
    cout << "Mass: " << i_mass << " integral " << h2_total_mass->Integral()
         << " lowMult: " << h2_lowMult_mass->Integral()
         << " highMult: " << h2_highMult_mass->Integral()
         << " highSubLow: " << h2_highSubLow_mass->Integral() << endl;

    StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(), h2_total_mass);
    StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(), h2_lowMult_mass);
    StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(), h2_highMult_mass);
    StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(),
                                 h2_highSubLow_mass);
    // for (int i_pt = 1; i_pt <= hnTool_se_pr.GetNbins(4); i_pt++) {
    //   auto h2_total_mass_pt = assoYeild.AssociatedYeild(i_mass, i_pt, 0);
    //   auto h2_lowMult_mass_pt = assoYeild.AssociatedYeild(i_mass, i_pt, 1);
    //   auto h2_highMult_mass_pt = assoYeild.AssociatedYeild(i_mass, i_pt, 2);
    //   auto h2_highSubLow_mass_pt = (TH2D *)h2_highMult_mass_pt->Clone(
    //       Form("h2_highSubLow_mass_pt_%d_%d", i_mass, i_pt));
    //   HistSubstraction2D(h2_highSubLow_mass_pt, h2_highMult_mass_pt,
    //                      h2_lowMult_mass_pt);
    //   StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(),
    //                                h2_total_mass_pt);
    //   StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(),
    //                                h2_lowMult_mass_pt);
    //   StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(),
    //                                h2_highMult_mass_pt);
    //   StyleFlow::DeltaPhi_DeltaEta(gPublisherCanvas->NewPad(),
    //                                h2_highSubLow_mass_pt);
    //   h2_total_mass_pt->SetName(Form("h2_total_mass_pt_%d_%d", i_mass, i_pt));
    //   h2_lowMult_mass_pt->SetName(
    //       Form("h2_lowMult_mass_pt_%d_%d", i_mass, i_pt));
    //   h2_highMult_mass_pt->SetName(
    //       Form("h2_highMult_mass_pt_%d_%d", i_mass, i_pt));
    //   h2_highSubLow_mass_pt->SetName(
    //       Form("h2_highSubLow_mass_pt_%d_%d", i_mass, i_pt));
    //   h2_total_mass_pt->Write();
    //   h2_lowMult_mass_pt->Write();
    //   h2_highMult_mass_pt->Write();
    //   h2_highSubLow_mass_pt->Write();
    // }
  }

  gPublisherCanvas->finalize();
  file_output->Close();
}

int main(int argc, char **argv) {
  gROOT->SetBatch(kTRUE);
  TString path_input_se_pr =
      "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
      "JpsiAsso_LHC22pass4_dqfilter.root:DeltaEtaUS_DeltaPhiUS_PosZUS_"
      "MassUS_"
      "PtUS_NumContribCalibUS";
  TString path_input_se_raw =
      "/home/szhu/work/alice/analysis/QA/input/jpsi/"
      "JpsiQA_LHC22pass4_dqfilter.root:fPosZ_MassUS_PtUS_NumContribCalib";
  TString path_input_me_pr =
      "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
      "MixEventReading_LHC22pass4_dqfilter.root:DeltaEtaUS_DeltaPhiUS_"
      "PosZUS_"
      "MassUS_PtUS_NumContribCalibUS";
  TString path_output =
      "/home/szhu/work/alice/analysis/QA/output/event_jpsi/AssoYeildQA";

  if (argc > 1) {
    path_input_se_pr = argv[1];
  }
  if (argc > 2) {
    path_input_se_raw = argv[2];
  }
  if (argc > 3) {
    path_input_me_pr = argv[3];
  }
  if (argc > 4) {
    path_output = argv[4];
  }
  AssoYeildQA(path_input_se_pr, path_input_se_raw, path_input_me_pr,
              path_output);

  return 0;
}