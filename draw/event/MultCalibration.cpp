#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TProfile2D.h"

void MultCalibration(
    TString path = "/home/szhu/work/alice/analysis/QA/input/event/"
                   "MultRaw_LHC22pass4_dqfilter.root",
    TString path_output = "/home/szhu/work/alice/analysis/QA/output/event/"
                          "MultCalibration_LHC22pass4_dqfilter.root") {
  // fNumContribfPosZ
  gROOT->SetBatch(true);
  MRootGraphic::StyleCommon();
  TProfile *fNumContribfPosZ = MRootIO::GetTProfile(path + ":fNumContribfPosZ");
  TProfile *fNumContribRun = MRootIO::GetTProfile(path + ":fNumContribRun");
  auto fNumContribfPosZRun =
      MRootIO::GetObjectDiectly<TProfile2D>(path + ":fNumContribfPosZRun");
  TFile *file_output = new TFile(path_output, "RECREATE");

  MRootGraphic::StyleHistCommonHist(fNumContribfPosZ);

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
  MRootGraphic::StyleHistCommonHist(fNumContribRun);
  fNumContribRun->GetYaxis()->SetRangeUser(30, 60);
  fNumContribRun->Draw();
  c_fNumContribRun->Modified();
  c_fNumContribRun->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/"
      "MultCalibration_fNumContribRun_LHC22pass4_dqfilter.pdf");
  c_fNumContribRun->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/"
      "MultCalibration_fNumContribRun_LHC22pass4_dqfilter.json");

  TCanvas *c_fNumContribfPosZRun =
      new TCanvas("c_fNumContribfPosZRun", "c_fNumContribfPosZRun", 800, 600);
  c_fNumContribfPosZRun->cd();
  c_fNumContribfPosZRun->SetGrid();
  fNumContribfPosZRun->SetTitle(
      "fNumContribfPosZRun;fPosZ [cm];Run;<fNumContrib>");
  fNumContribfPosZRun->GetXaxis()->SetTitle("fPosZ [cm]");
  fNumContribfPosZRun->GetYaxis()->SetTitle("Run");
  fNumContribfPosZRun->Draw("colz");
  c_fNumContribfPosZRun->Modified();
  c_fNumContribfPosZRun->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/"
      "MultCalibration_fNumContribfPosZRun_LHC22pass4_dqfilter.pdf");
  c_fNumContribfPosZRun->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/"
      "MultCalibration_fNumContribfPosZRun_LHC22pass4_dqfilter.json");

  TCanvas *c_fNumContribfPosZRun_calib = new TCanvas(
      "c_fNumContribfPosZRun_calib", "c_fNumContribfPosZRun_calib", 600, 600);
  c_fNumContribfPosZRun_calib->cd();
  for (int i_run = 1; i_run <= fNumContribfPosZRun->GetNbinsY(); i_run++) {
    TString label = fNumContribfPosZRun->GetYaxis()->GetBinLabel(i_run);
    TH1D *h_NumContribPosZ_calibration = fNumContribfPosZRun->ProjectionX(
        Form("h_NumContribPosZ_calibration_%s", label.Data()), i_run, i_run);
    for (int i_vtxz = 1; i_vtxz <= h_NumContribPosZ_calibration->GetNbinsX();
         i_vtxz++) {
      double content = h_NumContribPosZ_calibration->GetBinContent(i_vtxz);
      double error = h_NumContribPosZ_calibration->GetBinError(i_vtxz);
      double content1 = 50.0 / content;
      double error1 = content1 * content1 * error;
      h_NumContribPosZ_calibration->SetBinContent(i_vtxz, content1);
      h_NumContribPosZ_calibration->SetBinError(i_vtxz, error1);
    }
    TF1 *fit_func = new TF1(Form("pol4_%s", label.Data()), "pol4", -10, 10);
    h_NumContribPosZ_calibration->Fit(fit_func, "Q", "", -10, 10);
    // fNumContribfPosZRun->GetListOfFunctions()->Add((TF1 *)fit_func->Clone());
    file_output->WriteObject<TF1>(
        fit_func, Form("fNumContribfPosZRun_calib_%s", label.Data()));
    MRootGraphic::StyleHistCommonHist(h_NumContribPosZ_calibration);
    gPad->SetTopMargin(0.1);
    h_NumContribPosZ_calibration->SetTitle(Form("Run %s", label.Data()));
    h_NumContribPosZ_calibration->Draw();
    cout << "Run: " << label.Data() << endl;

    if (i_run == 1) {
      c_fNumContribfPosZRun_calib->SaveAs(
          "/home/szhu/work/alice/analysis/QA/plot/event/"
          "MultCalibration_fNumContribfPosZRun_calib_LHC22pass4_dqfilter.pdf(",
          Form("Title: Run%s", label.Data()));
    } else if (i_run == fNumContribfPosZRun->GetNbinsY()) { // last run
      c_fNumContribfPosZRun_calib->SaveAs(
          "/home/szhu/work/alice/analysis/QA/plot/event/"
          "MultCalibration_fNumContribfPosZRun_calib_LHC22pass4_dqfilter."
          "pdf)",
          Form("Title: Run%s", label.Data()));
    } else {
      c_fNumContribfPosZRun_calib->SaveAs(
          "/home/szhu/work/alice/analysis/QA/plot/event/"
          "MultCalibration_fNumContribfPosZRun_calib_LHC22pass4_dqfilter."
          "pdf",
          Form("Title: Run%s", label.Data()));
    }
  }

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
  fNumContribfPosZRun->Write();
  file_output->Close();
}