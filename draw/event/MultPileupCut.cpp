#include "MALICE.h"
#include "MHead.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"
#include "TColor.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TStyle.h"
// KEY: TH1D     fNumContrib;1   fNumContrib
// KEY: TH2D     fNumContrib_fMultTPC;1  fNumContrib_fMultTPC
// KEY: TH2D     NumContribCalib_fMultTPC;1      NumContribCalib_fMultTPC
// KEY: TH2D     fMultNTracksPV_fMultTPC;1       fMultNTracksPV_fMultTPC
// KEY: TProfile fNumContribRun;1        fNumContribRun
// KEY: TProfile fNumContribfPosZ;1      fNumContribVtxZ
// KEY: TProfile NumContribCalibPosZ;1   NumContribCalib
// KEY: TProfile NumContribCalibRun;1    fNumContribCalibRun

void MultPileupCut(TString path_input =
                       "/home/szhu/work/alice/analysis/QA/input/event/"
                       "MultPileupCut_LHC22pass4_dqfilter.root",
                   TString tag_period = "LHC22pass4_dqfilter") {
  gROOT->SetBatch(true);
  MRootGraphic::StyleCommon();
  TCanvas *c_NumContrib = new TCanvas("c_NumContrib", "c_NumContrib", 800, 400);
  c_NumContrib->Divide(2, 1);
  c_NumContrib->cd(1);
  auto fNumContribRun =
      MRootIO::GetObjectDiectly<TProfile>(path_input + ":fNumContribRun");
  auto NumContribCalibRun =
      MRootIO::GetObjectDiectly<TProfile>(path_input + ":NumContribCalibRun");
  MRootGraphic::StyleHistCommonHist(fNumContribRun);
  MRootGraphic::StyleHistCommonHist(NumContribCalibRun);
  fNumContribRun->SetLineColor(kRed);
  fNumContribRun->SetMarkerColor(kRed);
  NumContribCalibRun->SetLineColor(kBlue);
  NumContribCalibRun->SetMarkerColor(kBlue);

  fNumContribRun->SetTitle("fNumContribRun;Run;<fNumContrib>");
  fNumContribRun->GetYaxis()->SetRangeUser(0, 65);

  fNumContribRun->Draw();
  NumContribCalibRun->Draw("same");

  TLegend *legend = new TLegend(0.4, 0.3, 0.7, 0.5);
  legend->SetHeader("<Number of Vtx Contributors>");
  legend->AddEntry(fNumContribRun, "Raw", "l");
  legend->AddEntry(NumContribCalibRun, "Calibrated", "l");
  legend->SetLineColor(0);
  legend->SetTextSize(0.04);
  legend->Draw();

  c_NumContrib->cd(2);
  auto fNumContribfPosZ =
      MRootIO::GetObjectDiectly<TProfile>(path_input + ":fNumContribfPosZ");
  auto NumContribCalibPosZ =
      MRootIO::GetObjectDiectly<TProfile>(path_input + ":NumContribCalibPosZ");
  MRootGraphic::StyleHistCommonHist(fNumContribfPosZ);
  MRootGraphic::StyleHistCommonHist(NumContribCalibPosZ);
  fNumContribfPosZ->SetLineColor(kRed);
  fNumContribfPosZ->SetMarkerColor(kRed);
  NumContribCalibPosZ->SetLineColor(kBlue);
  NumContribCalibPosZ->SetMarkerColor(kBlue);

  fNumContribfPosZ->SetTitle("fNumContribfPosZ;fPosZ [cm];<fNumContrib>");
  fNumContribfPosZ->GetYaxis()->SetRangeUser(0, 65);

  fNumContribfPosZ->Draw();
  NumContribCalibPosZ->Draw("same");
  c_NumContrib->SaveAs("/home/szhu/work/alice/analysis/QA/plot/event/"
                       "MultPileupCut_NumContrib_RawCalibComparison_" +
                       tag_period + ".pdf");
  c_NumContrib->SaveAs("/home/szhu/work/alice/analysis/QA/plot/event/"
                       "MultPileupCut_NumContrib_RawCalibComparison_" +
                       tag_period + ".json");
}