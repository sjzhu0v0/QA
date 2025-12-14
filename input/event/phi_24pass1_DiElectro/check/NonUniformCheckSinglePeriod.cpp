#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TLegend.h"

void NonUniformCheckSinglePeriod(
    TString path_input = "NonUniform_551013_24DiElectron.root") {
  TFile *fInput = TFile::Open(path_input, "READ");

  // KEY: TH2D     fPTREF_fPhiREF;1
  // KEY: TH2D     fEtaREF_fPhiREF;1
  // KEY: TH2D     fPosZ_fPhiREF;1
  // KEY: TH2D     NumContribCalib_fPhiREF;1
  // KEY: TH2D     fEtaREF_fPosZ;1
  auto fPTREF_fPhiREF = (TH2D *)fInput->Get("fPTREF_fPhiREF");
  auto fEtaREF_fPhiREF = (TH2D *)fInput->Get("fEtaREF_fPhiREF");
  auto fPosZ_fPhiREF = (TH2D *)fInput->Get("fPosZ_fPhiREF");
  auto NumContribCalib_fPhiREF = (TH2D *)fInput->Get("NumContribCalib_fPhiREF");
  auto fEtaREF_fPosZ = (TH2D *)fInput->Get("fEtaREF_fPosZ");

  MRootGraphic::StyleHistCommon(fPTREF_fPhiREF);
  MRootGraphic::StyleHistCommon(fEtaREF_fPhiREF);
  MRootGraphic::StyleHistCommon(fPosZ_fPhiREF);
  MRootGraphic::StyleHistCommon(NumContribCalib_fPhiREF);
  MRootGraphic::StyleHistCommon(fEtaREF_fPosZ);

  TCanvas *c_phi_integral =
      new TCanvas("c_phi_integral", "c_phi_integral", 800, 800);
  TH1D *phi_integral = fPTREF_fPhiREF->ProjectionY("phi_integral");
  phi_integral->GetYaxis()->SetRangeUser(0, phi_integral->GetMaximum() * 1.2);
  phi_integral->GetYaxis()->SetMaxDigits(1);
  phi_integral->GetYaxis()->SetTitle("Normalized counts");
  phi_integral->DrawNormalized();
  c_phi_integral->SaveAs("phi_integral_af.pdf");

  TCanvas *c_phi_eta = new TCanvas("c_phi_eta", "c_phi_eta", 1680, 800);
  c_phi_eta->Divide(2, 1);
  c_phi_eta->cd(1);
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.15);
  fEtaREF_fPhiREF->GetZaxis()->SetMaxDigits(1);
  fEtaREF_fPhiREF->Draw("colz");

  c_phi_eta->cd(2);
  gStyle->SetOptStat(0);
  // draw it differentially and normalized, draw an legend
  for (int i = 1; i <= fEtaREF_fPhiREF->GetNbinsX() / 4; i++) {
    TH1D *h_proj = fEtaREF_fPhiREF->ProjectionY(Form("h_proj_%d", i),
                                                (i - 1) * 4 + 1, i * 4);
    h_proj->GetYaxis()->SetRangeUser(0, h_proj->GetMaximum() * 2.);
    h_proj->GetYaxis()->SetMaxDigits(1);
    h_proj->GetYaxis()->SetTitle("Normalized counts");
    h_proj->SetLineColor(i);
    if (i == 1) {
      h_proj->DrawNormalized();
    } else {
      h_proj->DrawNormalized("same");
    }
  }
  TLegend *leg_phi_eta = new TLegend(0.4, 0.6, 0.9, 0.9);
  for (int i = 1; i <= fEtaREF_fPhiREF->GetNbinsX() / 4; i++) {
    leg_phi_eta->AddEntry(
        Form("h_proj_%d", i),
        Form("Eta: %.2f - %.2f",
             fEtaREF_fPhiREF->GetXaxis()->GetBinLowEdge((i - 1) * 4 + 1),
             fEtaREF_fPhiREF->GetXaxis()->GetBinUpEdge(i * 4)),
        "l");
  }
  leg_phi_eta->Draw();
  c_phi_eta->SaveAs("phi_eta_diff_af.pdf");

  TCanvas *c_phi_z = new TCanvas("c_phi_z", "c_phi_z", 1680, 800);
  c_phi_z->Divide(2, 1);
  c_phi_z->cd(1);
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.15);
  fPosZ_fPhiREF->GetZaxis()->SetMaxDigits(1);
  fPosZ_fPhiREF->Draw("colz");
  c_phi_z->cd(2);
  gStyle->SetOptStat(0);
  // draw it differentially and normalized, draw an legend
  for (int i = 1; i <= fPosZ_fPhiREF->GetNbinsX() / 5; i++) {
    TH1D *h_proj = fPosZ_fPhiREF->ProjectionY(Form("h_proj_z_%d", i),
                                              (i - 1) * 5 + 1, i * 5);
    h_proj->GetYaxis()->SetRangeUser(0, h_proj->GetMaximum() * 2.);
    h_proj->GetYaxis()->SetMaxDigits(1);
    h_proj->GetYaxis()->SetTitle("Normalized counts");
    h_proj->SetLineColor(i);
    if (i == 1) {
      h_proj->DrawNormalized();
    } else {
      h_proj->DrawNormalized("same");
    }
  }
  TLegend *leg_phi_z = new TLegend(0.4, 0.6, 0.9, 0.9);
  for (int i = 1; i <= fPosZ_fPhiREF->GetNbinsX() / 5; i++) {
    leg_phi_z->AddEntry(
        Form("h_proj_z_%d", i),
        Form("PosZ: %.2f - %.2f cm",
             fPosZ_fPhiREF->GetXaxis()->GetBinLowEdge((i - 1) * 5 + 1),
             fPosZ_fPhiREF->GetXaxis()->GetBinUpEdge(i * 5)),
        "l");
  }
  leg_phi_z->Draw();
  c_phi_z->SaveAs("phi_z_diff_af.pdf");
}
