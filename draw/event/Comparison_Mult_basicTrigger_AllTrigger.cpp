#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TLegend.h"

void Comparison_Mult_basicTrigger_AllTrigger(
    TString path_file_input = "/home/szhu/work/alice/analysis/QA/input/event/"
                              "MultQA_AllCut_LHC22pass4_dqfilter.root",
    TString tag_period = "LHC22pass4_dqfilter") {
  TFile *file_input = TFile::Open(path_file_input);
  MRootGraphic::StyleCommon();

  TCanvas *c_levelTrigger =
      new TCanvas("c_levelTrigger", "c_levelTrigger", 800, 600);
  TH1D *LevelTrigger =
      MRootIO::GetObjectDiectly<TH1D>(file_input, "LevelTrigger");
  MRootGraphic::StyleHistCommonHist(LevelTrigger);
  LevelTrigger->Scale(1. / LevelTrigger->GetBinContent(1));
  //  {"isTriggerTVX", "isntTimeFrameBorder", "isntITSROFrameBorder",
  //              "isntSameBunchPileup", "isntSelfDefinedPileup"};
  LevelTrigger->GetXaxis()->SetTitle("");
  LevelTrigger->GetYaxis()->SetTitle("Events normaliezed to all events (a.u.)");
  LevelTrigger->GetXaxis()->SetBinLabel(1, "All event");
  LevelTrigger->GetXaxis()->SetBinLabel(2, "TVX Trigger");
  LevelTrigger->GetXaxis()->SetBinLabel(3, "No time frame border");
  LevelTrigger->GetXaxis()->SetBinLabel(4, "No ITS RO frame border");
  LevelTrigger->GetXaxis()->SetBinLabel(5, "No same bunch pileup");
  LevelTrigger->GetXaxis()->SetBinLabel(6, "MultTPC tail removed");
  LevelTrigger->GetXaxis()->SetLabelSize(0.04);
  LevelTrigger->GetYaxis()->SetRangeUser(0, 1.2);
  gPad->SetGridy();
  gPad->SetRightMargin(0.150376);
  LevelTrigger->SetTitle("");
  LevelTrigger->GetYaxis()->SetTitle("Events normaliezed to all events");
  LevelTrigger->Draw("hist text");
  c_levelTrigger->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/"
      "Comparison_Mult_basicTrigger_AllTrigger_LevelTrigger.pdf");
  c_levelTrigger->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/"
      "Comparison_Mult_basicTrigger_AllTrigger_LevelTrigger.json");

  TCanvas *c_levelTrigger_0_5 =
      new TCanvas("c_levelTrigger_0_5", "c_levelTrigger", 800, 600);
  gPad->SetRightMargin(0.15);
  gPad->SetGridy();
  LevelTrigger->GetXaxis()->SetRangeUser(-0.5, 4.5);
  LevelTrigger->Draw("hist text");
  c_levelTrigger_0_5->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/"
      "Comparison_Mult_basicTrigger_AllTrigger_LevelTrigger_0to5.pdf");
  c_levelTrigger_0_5->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/"
      "Comparison_Mult_basicTrigger_AllTrigger_LevelTrigger_0to5.json");

  vector<string> vec_name_hist;
  vector<string> vec_name_hist_sort1d;
  vector<string> vec_name_hist_sort2d;

  TList *list_hist = file_input->GetListOfKeys();
  for (int i = 0; i < list_hist->GetEntries(); i++) {
    TKey *key = (TKey *)list_hist->At(i);
    TString name_hist = key->GetName();
    // if the end of name isn't fullTrigger or basicTrigger, skip it
    if (!name_hist.EndsWith("fullTrigger") &&
        !name_hist.EndsWith("basicTrigger")) {
      continue;
    }
    vec_name_hist.push_back(name_hist.Data());
    // remove the suffix "_fullTrigger" or "_basicTrigger"
    TH1 *hist = MRootIO::GetObjectDiectly<TH1>(file_input, name_hist);
    int dim = hist->GetDimension();
    if (name_hist.EndsWith("fullTrigger")) {
      name_hist.ReplaceAll("_fullTrigger", "");
    } else if (name_hist.EndsWith("basicTrigger")) {
      name_hist.ReplaceAll("_basicTrigger", "");
    }
    // check if the name_hist is 1D or 2D histogram
    if (dim == 1) {
      vec_name_hist_sort1d.push_back(name_hist.Data());
    } else if (dim == 2) {
      vec_name_hist_sort2d.push_back(name_hist.Data());
    }
  }
  // unique sort the vector
  sort(vec_name_hist_sort1d.begin(), vec_name_hist_sort1d.end());
  sort(vec_name_hist_sort2d.begin(), vec_name_hist_sort2d.end());
  // unique the vector
  vec_name_hist_sort1d.erase(
      unique(vec_name_hist_sort1d.begin(), vec_name_hist_sort1d.end()),
      vec_name_hist_sort1d.end());
  vec_name_hist_sort2d.erase(
      unique(vec_name_hist_sort2d.begin(), vec_name_hist_sort2d.end()),
      vec_name_hist_sort2d.end());

  TCanvas *c_comparison_2d =
      new TCanvas("c_comparison_2d", "c_comparison_2d", 800, 400);

  for (int i = 0; i < vec_name_hist_sort2d.size(); i++) {
    TString name_hist = vec_name_hist_sort2d[i];
    TH1 *hist_2d_fullTrigger =
        MRootIO::GetObjectDiectly<TH2D>(file_input, name_hist + "_fullTrigger");
    TH2D *hist_2d_basicTrigger = MRootIO::GetObjectDiectly<TH2D>(
        file_input, name_hist + "_basicTrigger");
    MRootGraphic::StyleHistCommonHist(hist_2d_fullTrigger);
    MRootGraphic::StyleHistCommonHist(hist_2d_basicTrigger);
    c_comparison_2d->Clear();
    c_comparison_2d->Divide(2, 1);
    c_comparison_2d->cd(1);
    gPad->SetLogz();
    gPad->SetRightMargin(0.12);
    hist_2d_fullTrigger->SetTitle(
        Form("%s, %s", name_hist.Data(), tag_period.Data()));
    hist_2d_fullTrigger->Draw("colz");
    TLatex *latex = new TLatex();
    latex->SetTextSize(0.04);
    latex->SetTextAlign(12);
    latex->SetNDC();
    latex->DrawLatex(0.6, 0.9, "Full Trigger");
    c_comparison_2d->cd(2);
    gPad->SetLogz();
    gPad->SetRightMargin(0.12);
    hist_2d_basicTrigger->SetTitle(
        Form("%s, %s", name_hist.Data(), tag_period.Data()));
    hist_2d_basicTrigger->Draw("colz");
    latex->DrawLatex(0.6, 0.9, "Basic Trigger");
    // show z axis
    // save the canvas
    if (i == 0) {
      c_comparison_2d->SaveAs(
          "/home/szhu/work/alice/analysis/QA/plot/event/"
          "Comparison_Mult_basicTrigger_AllTrigger_comparison2d.pdf(",
          Form("Title: %s", name_hist.Data()));
    } else if (i == vec_name_hist_sort2d.size() - 1) {
      c_comparison_2d->SaveAs(
          "/home/szhu/work/alice/analysis/QA/plot/event/"
          "Comparison_Mult_basicTrigger_AllTrigger_comparison2d.pdf)",
          Form("Title: %s", name_hist.Data()));
    } else {
      c_comparison_2d->SaveAs(
          "/home/szhu/work/alice/analysis/QA/plot/event/"
          "Comparison_Mult_basicTrigger_AllTrigger_comparison2d.pdf",
          Form("Title: %s", name_hist.Data()));
    }
  }

  TCanvas *c_comparison_1d =
      new TCanvas("c_comparison_1d", "c_comparison_1d", 800, 800);
  for (int i = 0; i < vec_name_hist_sort1d.size(); i++) {
    TString name_hist = vec_name_hist_sort1d[i];
    TH1D *hist_fullTrigger =
        MRootIO::GetObjectDiectly<TH1D>(file_input, name_hist + "_fullTrigger");
    TH1D *hist_basicTrigger = MRootIO::GetObjectDiectly<TH1D>(
        file_input, name_hist + "_basicTrigger");
    MRootGraphic::StyleHistCommonHist(hist_fullTrigger);
    MRootGraphic::StyleHistCommonHist(hist_basicTrigger);
    c_comparison_1d->Clear();
    double range_user_max, range_user_min;
    double range_max =
        max(hist_fullTrigger->GetMaximum(), hist_basicTrigger->GetMaximum());
    double range_min =
        min(hist_fullTrigger->GetMinimum(), hist_basicTrigger->GetMinimum());
    if (range_min < 0) {
      range_min = 0;
    }
    if (range_min == 0) {
      range_user_max = range_max * 10;
      range_user_min = 0.1;
    } else {
      range_user_max = range_max * 10.;
      range_user_min = range_min / 10.;
    }
    c_comparison_1d->cd();
    gPad->SetGridy();
    gPad->SetLogy();
    hist_fullTrigger->GetYaxis()->SetRangeUser(range_user_min, range_user_max);
    hist_fullTrigger->Draw("hist");
    hist_basicTrigger->SetLineColor(kRed);
    hist_basicTrigger->Draw("hist same");
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->AddEntry(hist_fullTrigger, "Full Trigger", "l");
    legend->AddEntry(hist_basicTrigger, "Basic Trigger", "l");
    legend->Draw();

    if (i == 0) {
      c_comparison_1d->SaveAs(
          "/home/szhu/work/alice/analysis/QA/plot/event/"
          "Comparison_Mult_basicTrigger_AllTrigger_comparison1d.pdf(",
          Form("Title: %s", name_hist.Data()));
    } else if (i == vec_name_hist_sort1d.size() - 1) {
      c_comparison_1d->SaveAs(
          "/home/szhu/work/alice/analysis/QA/plot/event/"
          "Comparison_Mult_basicTrigger_AllTrigger_comparison1d.pdf)",
          Form("Title: %s", name_hist.Data()));
    } else {
      c_comparison_1d->SaveAs(
          "/home/szhu/work/alice/analysis/QA/plot/event/"
          "Comparison_Mult_basicTrigger_AllTrigger_comparison1d.pdf",
          Form("Title: %s", name_hist.Data()));
    }
  }
  for (auto name : vec_name_hist_sort1d)
    cout << name << endl;
}