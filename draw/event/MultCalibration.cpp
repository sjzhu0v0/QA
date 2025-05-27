#include "MRootGraphic.h"
#include "MRootIO.h"

void MultCalibration(
    TString path = "/home/szhu/work/alice/analysis/QA/input/event/"
                   "MultRaw_LHC22pass4_dqfilter.root",
    TString path_output = "/home/szhu/work/alice/analysis/QA/output/event/"
                          "MultCalibration_LHC22pass4_dqfilter.root") {
  // fNumContribfPosZ
  MRootGraphic::StyleCommon();
  TProfile *fNumContribfPosZ = MRootIO::GetTProfile(path + ":fNumContribfPosZ");
  TProfile *fNumContribRun = MRootIO::GetTProfile(path + ":fNumContribRun");
  TFile *file_output = new TFile(path_output, "RECREATE");

  MRootGraphic::StyleHistCommon1D(fNumContribfPosZ);

  TCanvas *c_fNumContribfPosZ =
      new TCanvas("c_fNumContribfPosZ", "c_fNumContribfPosZ", 600, 600);
  c_fNumContribfPosZ->cd();
  c_fNumContribfPosZ->SetGrid();
  fNumContribfPosZ->SetTitle("fNumContribfPosZ;fPosZ [cm];<fNumContrib>");
  fNumContribfPosZ->GetYaxis()->SetRangeUser(0, 60);
  fNumContribfPosZ->Draw();
  c_fNumContribfPosZ->Modified();
  c_fNumContribfPosZ->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/"
      "MultCalibration_fNumContribfPosZ_LHC22pass4_dqfilter.pdf");
  c_fNumContribfPosZ->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/"
      "MultCalibration_fNumContribfPosZ_LHC22pass4_dqfilter.json");

  TCanvas *c_fNumContribRun =
      new TCanvas("c_fNumContribRun", "c_fNumContribRun", 1200, 600);
  c_fNumContribRun->cd();
  c_fNumContribRun->SetGrid();
  fNumContribRun->GetXaxis()->SetTitle("");
  fNumContribRun->GetYaxis()->SetTitle("<fNumContrib>");
  fNumContribRun->GetYaxis()->SetLabelSize(0.02);
  MRootGraphic::StyleHistCommon1D(fNumContribRun);
  fNumContribRun->GetYaxis()->SetRangeUser(30, 60);
  fNumContribRun->Draw();
  c_fNumContribRun->Modified();
  c_fNumContribRun->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/"
      "MultCalibration_fNumContribRun_LHC22pass4_dqfilter.pdf");
  c_fNumContribRun->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/"
      "MultCalibration_fNumContribRun_LHC22pass4_dqfilter.json");

  TCanvas *c_h_fNumContribfPosZ =
      new TCanvas("c_h_fNumContribfPosZ", "c_h_fNumContribfPosZ", 600, 600);
  c_h_fNumContribfPosZ->cd();
  c_h_fNumContribfPosZ->SetGrid();
  TH1D *h_fNumContribfPosZ =
      fNumContribfPosZ->ProjectionX("h_NumContribPosZ_calibration");
  h_fNumContribfPosZ->GetYaxis()->SetTitle("1/(Normlized <fNumContrib>)");
  h_fNumContribfPosZ->Sumw2(false);
  h_fNumContribfPosZ->Scale(h_fNumContribfPosZ->GetNbinsX() /
                            h_fNumContribfPosZ->Integral());
  for (int i = 1; i <= h_fNumContribfPosZ->GetNbinsX(); i++) {
    double content = h_fNumContribfPosZ->GetBinContent(i);
    double error = h_fNumContribfPosZ->GetBinError(i);
    double content1 = 1.0 / content;
    double error1 = content1 * content1 * error;
    h_fNumContribfPosZ->SetBinContent(i, content1);
    h_fNumContribfPosZ->SetBinError(i, error1);
  }
  h_fNumContribfPosZ->Fit("pol4", "Q", "", -10, 10);
  h_fNumContribfPosZ->GetYaxis()->SetRangeUser(0, 1.3);
  h_fNumContribfPosZ->Draw();
  c_h_fNumContribfPosZ->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/"
      "MultCalibration_fNumContribfPosZ_Calibration_LHC22pass4_dqfilter.pdf");
  c_h_fNumContribfPosZ->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/"
      "MultCalibration_fNumContribfPosZ_Calibration_LHC22pass4_dqfilter.json");

  file_output->cd();
  h_fNumContribfPosZ->Write();
  fNumContribRun->Write();
  file_output->Close();
}