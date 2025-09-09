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

  StrVar4Hist var_fPosX("fPosX", "#it{V}_{x}", "cm", 200, {-10, 10});
  StrVar4Hist var_fPosY("fPosY", "#it{V}_{Y}", "cm", 200, {-10, 10});
  StrVar4Hist var_fPosZ("fPosZ", "#it{V}_{Z}", "cm", 200, {-10, 10});
  StrVar4Hist var_fNumContrib("fNumContrib", "#it{N}_{vtx contrib} ", "", 300,
                              {0, 300});
  StrVar4Hist var_NumContribCalib(
      "NumContribCalib", "N_{vtx contrib} Calibrated", "", 300, {0, 300});
  StrVar4Hist var_fMultTPC("fMultTPC", "Mult_{TPC}", "", 600, {0, 600});
  StrVar4Hist var_fMultREF("fMultREF", "Mult_{REF}", "", 100, {0, 100});
  StrVar4Hist var_fMultFT0C("fMultFT0C", "Mult_{FT0C}", "", 130,
                            {-1000., 12000.});
  StrVar4Hist var_MultNTracksPV("fMultNTracksPV", "#it{N}_{Tracks PV}", "", 150,
                                {0, 150});
  StrVar4Hist var_MassJpsiCandidate("fMass", "M_{ee}", "GeV^{2}/c^{4}", 100,
                                    {1.8, 5.4});
  StrVar4Hist var_PtJpsiCandidate("fPT", "M_{ee}", "GeV/c", 10, {0., 10.});

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

  TCanvas *c_fNumContrib_fMultREF =
      new TCanvas("c_fNumContrib_fMultREF", "c_fNumContrib_fMultREF", 900, 900);
  c_fNumContrib_fMultREF->Divide(3, 3);
  gStyle->SetOptStat("e");

  vector<TH1D *> vec_proj_lowNumContrib;

  auto setTitleSize = [](TH1 *hist, double size) {
    hist->GetXaxis()->SetTitleSize(size);
    hist->GetYaxis()->SetTitleSize(size);
    hist->GetYaxis()->SetTitleOffset(0.9);
    hist->GetZaxis()->SetTitleSize(size);
  };
  vector<TH2D *> vec_fNumContrib_fMultREF;
  vector<TH1D *> vec_proj_lowNumContrib_x;
  for (int i_CondSameBunchPileup = 0; i_CondSameBunchPileup < 3;
       i_CondSameBunchPileup++) {
    auto fNumContrib_fMultREF = MRootIO::GetObjectDiectly<TH2D>(
        file_input,
        GetHistName2D(
            var_NumContribCalib, var_fMultREF,
            (TString)conditions_samePileup[i_CondSameBunchPileup][1]));
    vec_fNumContrib_fMultREF.push_back(fNumContrib_fMultREF);
    TString title_x = var_NumContribCalib.fTitle;
    if (var_NumContribCalib.fUnit != "")
      title_x += " (" + var_NumContribCalib.fUnit + ")";
    TString title_y = var_fMultREF.fTitle;
    if (var_fMultREF.fUnit != "")
      title_y += " (" + var_fMultREF.fUnit + ")";
    fNumContrib_fMultREF->GetXaxis()->SetTitle(title_x);
    fNumContrib_fMultREF->GetYaxis()->SetTitle(title_y);
    MRootGraphic::StyleHistCommonHist(fNumContrib_fMultREF);
    c_fNumContrib_fMultREF->cd(i_CondSameBunchPileup + 1);
    gPad->SetLogz();
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.11);
    gPad->SetTopMargin(0.1);
    // get the pavestats from the canvas
    setTitleSize(fNumContrib_fMultREF, 0.05);
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
    setTitleSize(proj_lowNumContrib_y, 0.05);
    proj_lowNumContrib_y->Draw();
    TLatex *latex = new TLatex();
    latex->SetTextSize(0.05);
    latex->SetTextAlign(12);
    latex->SetNDC();
    latex->DrawLatex(0.55, 0.85, "#it{N}_{contrib}<4");
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
    c_fNumContrib_fMultREF->cd(i_CondSameBunchPileup + 7);
    auto proj_lowNumContrib_x = (TH1D *)fNumContrib_fMultREF->ProjectionX(
        GetHistName1D(var_fMultREF,
                      conditions_samePileup[i_CondSameBunchPileup][1]),
        1, 4);
    vec_proj_lowNumContrib_x.push_back(proj_lowNumContrib_x);
    MRootGraphic::StyleHistCommonHist(proj_lowNumContrib_x);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.11);
    gPad->SetTopMargin(0.1);
    setTitleSize(proj_lowNumContrib_x, 0.05);
    proj_lowNumContrib_x->Draw();
    latex->DrawLatex(0.55, 0.85, "MultREF<4");
    gPad->Update();
    st = (TPaveStats *)proj_lowNumContrib_x->FindObject("stats");
    if (st) {
      st->SetTextSize(0.03);
      st->SetX1NDC(0.2);
      st->SetX2NDC(0.5);
      st->SetY1NDC(0.8);
      st->SetY2NDC(0.9);
    }
    gPad->Modified();
  }
  c_fNumContrib_fMultREF->SaveAs("/home/szhu/work/alice/analysis/QA/plot/event/"
                                 "MultREFRaw_NumContrib_vs_MultREF_" +
                                 tag_period + ".pdf");
  c_fNumContrib_fMultREF->SaveAs("/home/szhu/work/alice/analysis/QA/plot/event/"
                                 "MultREFRaw_NumContrib_vs_MultREF_" +
                                 tag_period + ".json");

  TCanvas *c_fNumContrib_fMultREF_atLowMultREF =
      new TCanvas("c_fNumContrib_fMultREF_atLowMultREF",
                  "c_fNumContrib_fMultREF_atLowMultREF", 400, 400);
  c_fNumContrib_fMultREF_atLowMultREF->cd();
  gPad->SetLogy();
  for (auto hist : vec_proj_lowNumContrib_x) {
    hist->Scale(1. / hist->Integral("width"));
    TPaveStats *st = (TPaveStats *)hist->FindObject("stats");
    if (st) {
      hist->GetListOfFunctions()->Remove(st);
    }
  }
  gStyle->SetOptStat(0);

  auto fNumContrib_noCut_lowMult = vec_proj_lowNumContrib_x[0];
  auto fNumContrib_onlySameBunch_lowMult = vec_proj_lowNumContrib_x[2];
  auto fNumContrib_noSameBunch_lowMult = vec_proj_lowNumContrib_x[1];

  fNumContrib_noCut_lowMult->SetLineColor(kBlack);
  fNumContrib_noSameBunch_lowMult->SetLineColor(kRed);
  fNumContrib_onlySameBunch_lowMult->SetLineColor(kBlue);

  fNumContrib_onlySameBunch_lowMult->Draw();
  fNumContrib_noCut_lowMult->Draw("same");
  fNumContrib_noSameBunch_lowMult->Draw("same");

  TLegend *legend_proj_lowNumContrib_x = new TLegend(0.5, 0.6, 0.8, 0.8);
  legend_proj_lowNumContrib_x->SetHeader("Mult_{REF}<4");
  legend_proj_lowNumContrib_x->AddEntry(fNumContrib_noCut_lowMult, "No cut",
                                        "l");
  legend_proj_lowNumContrib_x->AddEntry(fNumContrib_noSameBunch_lowMult,
                                        "No same bunch pileup", "l");
  legend_proj_lowNumContrib_x->AddEntry(fNumContrib_onlySameBunch_lowMult,
                                        "Only same bunch pileup", "l");
  legend_proj_lowNumContrib_x->SetLineColor(0);
  legend_proj_lowNumContrib_x->SetTextSize(0.04);
  legend_proj_lowNumContrib_x->Draw();
  c_fNumContrib_fMultREF_atLowMultREF->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/"
      "MultREFRaw_NumContrib_atLowMultREF_" +
      tag_period + ".pdf");
  c_fNumContrib_fMultREF_atLowMultREF->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/"
      "MultREFRaw_NumContrib_atLowMultREF_" +
      tag_period + ".json");

  TCanvas *c_NumContrib = new TCanvas("c_NumContrib", "c_NumContrib", 800, 800);
  auto fNumContrib_isntSameBunchPileup =
      (TH1D *)vec_fNumContrib_fMultREF[1]->ProjectionX(
          GetHistName1D(var_NumContribCalib, "isntSameBunchPileup"));
  c_NumContrib->cd();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  MRootGraphic::StyleHistCommonHist(fNumContrib_isntSameBunchPileup);
  fNumContrib_isntSameBunchPileup->SetTitle("");
  fNumContrib_isntSameBunchPileup->Draw();

  double sum = fNumContrib_isntSameBunchPileup->Integral();
  double sum_per_bin = sum / 10;
  double mean = fNumContrib_isntSameBunchPileup->GetMean();

  vector<double> vec_bin;

  vec_bin.push_back(0.);
  double sum_bin = 0.;
  for (int i = 1; i <= fNumContrib_isntSameBunchPileup->GetNbinsX(); i++) {
    sum_bin += fNumContrib_isntSameBunchPileup->GetBinContent(i);
    if (sum_bin >= sum_per_bin * vec_bin.size()) {
      vec_bin.push_back(fNumContrib_isntSameBunchPileup->GetBinLowEdge(i + 1));
    }
  }

  for (auto value : vec_bin) {
    cout << value << ",";
  }

  for (auto &value : vec_bin) {
    double fraction = value / mean;
    cout << "Bin value: " << value << ", Fraction of mean: " << fraction
         << endl;
  }

  for (int i = 0; i < vec_bin.size() - 1; i++) {
    double x = vec_bin[i];
    double y = fNumContrib_isntSameBunchPileup->GetBinContent(
        fNumContrib_isntSameBunchPileup->FindBin(x));
    TLine *line = new TLine(x, 0, x, y);
    line->SetLineColor(kRed);
    line->Draw();
  }
  c_NumContrib->SaveAs("/home/szhu/work/alice/analysis/QA/plot/event/"
                       "MultREFRaw_NumContribCalibBinning_" +
                       tag_period + ".pdf");
  c_NumContrib->SaveAs("/home/szhu/work/alice/analysis/QA/plot/event/"
                       "MultREFRaw_NumContribCalibBinning_" +
                       tag_period + ".json");
}