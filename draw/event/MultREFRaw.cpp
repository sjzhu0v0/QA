#include "MALICE.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"
#include "TColor.h"
#include "TLegend.h"
#include "TPaveStats.h"
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

void MultREFRaw(
    TString path_input = "/home/szhu/work/alice/analysis/QA/input/event/"
                         "MultREFRaw_LHC22pass4_dqfilter.root",
    TString path_output = "/home/szhu/work/alice/analysis/QA/output/event/"
                          "MultREFRaw_LHC22pass4_dqfilter.root",
    TString tag_period = "LHC22pass4_dqfilter") {
  //   gROOT->SetBatch(true);

  StrVar4Hist var_fNumContrib("fNumContrib", "fNumContrib", "", 300, {0, 300});
  StrVar4Hist var_NumContribCalib("NumContribCalib", "NumContrib Calibrated",
                                  "", 300, {0, 300});
  StrVar4Hist var_fMultTPC("fMultTPC", "fMultTPC", "", 600, {0, 600});
  StrVar4Hist var_fMultREF("fMultREF", "fMultREF", "", 100, {0, 100});
  StrVar4Hist var_fMultFT0C("fMultFT0C", "fMultFT0C", "", 130,
                            {-1000., 12000.});
  StrVar4Hist var_fMultNTracksPV("fMultNTracksPV", "fMultNTracksPV", "a.u.",
                                 150, {0, 150});

  vector<array<string, 2>> conditions_samePileup = {
      {"isntSameBunchPileup || !isntSameBunchPileup", "NoSameBunchCut"},
      {"isntSameBunchPileup", "NoSameBunchPileup"},
      {"!isntSameBunchPileup", "AllSameBunchPileup"}};

  auto GetHistName1D = [](const StrVar4Hist &var, TString tag) {
    return var.fName + "_" + (TString)tag[1] + "_" + tag;
  };

  auto GetHistName2D = [](const StrVar4Hist &var1, const StrVar4Hist &var2,
                          TString tag) {
    return var1.fName + "_" + var2.fName + "_" + tag;
  };

  TFile *file_input = TFile::Open(path_input, "READ");
  TFile *file_output = new TFile(path_output, "RECREATE");
  MRootGraphic::StyleCommon();

  TCanvas *c_fNumContrib_fMultREF = new TCanvas(
      "c_fNumContrib_fMultREF", "c_fNumContrib_fMultREF", 1200, 800);
  c_fNumContrib_fMultREF->Divide(3, 3);
  gStyle->SetOptStat("e");

  vector<TH1D *> vec_proj_lowNumContrib;
  for (int i_CondSameBunchPileup = 0; i_CondSameBunchPileup < 3;
       i_CondSameBunchPileup++) {
    auto fNumContrib_fMultREF = MRootIO::GetObjectDiectly<TH2D>(
        file_input,
        GetHistName2D(
            var_NumContribCalib, var_fMultREF,
            (TString)conditions_samePileup[i_CondSameBunchPileup][1]));
    MRootGraphic::StyleHistCommonHist(fNumContrib_fMultREF);
    c_fNumContrib_fMultREF->cd(i_CondSameBunchPileup + 1);
    gPad->SetLogz();
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.11);
    gPad->SetTopMargin(0.1);
    // get the pavestats from the canvas
    fNumContrib_fMultREF->Draw("colz");
    gPad->Update();
    TPaveStats *st = (TPaveStats *)fNumContrib_fMultREF->FindObject("stats");
    if (st) {
      st->SetX1NDC(0.15);
      st->SetX2NDC(0.45);
      st->SetY1NDC(0.8);
      st->SetY2NDC(0.9);
      st->SetTextSize(0.03);
    }
    gPad->Modified();
    auto proj_lowNumContrib_y = (TH1D *)fNumContrib_fMultREF->ProjectionY(
        GetHistName1D(var_NumContribCalib,
                      conditions_samePileup[i_CondSameBunchPileup][1]),
        1, 4);
    MRootGraphic::StyleHistCommonHist(proj_lowNumContrib_y);
    c_fNumContrib_fMultREF->cd(i_CondSameBunchPileup + 4);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.11);
    gPad->SetTopMargin(0.1);
    proj_lowNumContrib_y->Draw();
    gPad->Update();
    st = (TPaveStats *)proj_lowNumContrib_y->FindObject("stats");
    if (st) {
      st->SetTextSize(0.03);
      st->SetX1NDC(0.2);
      st->SetX2NDC(0.5);
      st->SetY1NDC(0.8);
      st->SetY2NDC(0.9);
    }
    gPad->Modified();
    vec_proj_lowNumContrib.push_back(proj_lowNumContrib_y);
  }
}