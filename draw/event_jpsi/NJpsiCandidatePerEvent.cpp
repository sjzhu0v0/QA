#include "MALICE.h"
#include "MHead.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"
#include "TColor.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TStyle.h"

void NJpsiCandidatePerEvent(
    TString path_input =
        "/home/szhu/data/DoubleJpsi/22pass4_highIR/doubleJpsi.root") {
  MRootGraphic::StyleCommon();

  auto NJpsiCandidata = MRootIO::GetTH1D(path_input + ":NJpsiCandidata");
  auto pair_mass = MRootIO::GetTH2D(path_input + ":pair_mass");
  pair_mass->Rebin2D(5, 5);

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  gPad->SetLogy();
  MRootGraphic::StyleHistCommonHist(NJpsiCandidata);
  NJpsiCandidata->GetXaxis()->SetTitle("N_{J/#psi candidates}");
  NJpsiCandidata->SetTitle(
      "N_{J/#psi candidates}: 2.5 GeV/c^{2} < M_{e^{+}e^{-}} < 3.2 GeV/c^{2}");
  NJpsiCandidata->GetYaxis()->SetTitle("Events");
  NJpsiCandidata->GetXaxis()->SetNdivisions(10);
  NJpsiCandidata->Draw();
  c1->SaveAs("/home/szhu/work/alice/analysis/QA/plot/event_jpsi/"
             "NJpsiCandidatePerEvent_NJpsiCandidata.pdf");
  c1->SaveAs("/home/szhu/work/alice/analysis/QA/plot/event_jpsi/"
             "NJpsiCandidatePerEvent_NJpsiCandidata.json");
  cout << NJpsiCandidata->Integral() << endl;
  cout << NJpsiCandidata->GetBinContent(1) << endl;
  cout << NJpsiCandidata->GetBinContent(2) << endl;

  TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
  StyleFlow::DeltaPhi_DeltaEta(c2, pair_mass);
  pair_mass->GetXaxis()->SetTitle("M_{e^{+}e^{-}} (GeV/c^{2})");
  pair_mass->GetYaxis()->SetTitle("M_{e^{+}e^{-}} (GeV/c^{2})");
  pair_mass->GetZaxis()->SetTitle("Counts");
  pair_mass->GetZaxis()->SetTitleOffset(1.2);
  c2->SaveAs("/home/szhu/work/alice/analysis/QA/plot/event_jpsi/"
             "NJpsiCandidatePerEvent_pair_mass.pdf");
  c2->SaveAs("/home/szhu/work/alice/analysis/QA/plot/event_jpsi/"
             "NJpsiCandidatePerEvent_pair_mass.json");
}