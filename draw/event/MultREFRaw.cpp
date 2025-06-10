#include "MALICE.h"
#include "MHead.h"
#include "MHist.h"
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

void MultREFRaw(
    TString path_input = "/home/szhu/work/alice/analysis/QA/input/event/"
                         "MultREFRaw_LHC22pass4_dqfilter.root",
    TString path_output = "/home/szhu/work/alice/analysis/QA/output/event/"
                          "MultREFRaw_LHC22pass4_dqfilter.root",
    TString tag_period = "LHC22pass4_dqfilter") {
  gROOT->SetBatch(true);

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

  auto GetHistName1D = [](const StrVar4Hist &var, const array<string, 2> &tag) {
    return var.fName + "_" + tag[1] + "_" + tag[1];
    };

  MRootGraphic::StyleCommon();
  TFile *file_output = new TFile(path_output, "RECREATE");

  TCanvas *c_fNumContrib_fMultREF =
      new TCanvas("c_fNumContrib_fMultREF", "c_fNumContrib_fMultREF", 800, 400);
  c_fNumContrib_fMultREF->Divide(2, 1);
  c_fNumContrib_fMultREF->cd(1);
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.12);
  auto fNumContrib_fMultREF =
      MRootIO::GetObjectDiectly<TH2D>(path_input + ":fNumContrib_fMultREF");
  auto NumContribCalib_fMultREF =
      MRootIO::GetObjectDiectly<TH2D>(path_input + ":NumContribCalib_fMultREF");
  MRootGraphic::StyleHistCommonHist(fNumContrib_fMultREF);
  MRootGraphic::StyleHistCommonHist(NumContribCalib_fMultREF);
  fNumContrib_fMultREF->SetTitle("fNumContrib_fMultREF;fNumContrib;fMultREF");
  NumContribCalib_fMultREF->SetTitle(
      "NumContribCalib_fMultREF;fNumContrib Calib;fMultREF");
  fNumContrib_fMultREF->GetZaxis()->SetTitle("");
  NumContribCalib_fMultREF->GetZaxis()->SetTitle("");
  gPad->SetLogz();
  fNumContrib_fMultREF->Draw("colz");
  c_fNumContrib_fMultREF->cd(2);
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.12);
  gPad->SetLogz();
  NumContribCalib_fMultREF->Draw("colz");
  c_fNumContrib_fMultREF->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/"
      "MultCalibrationResult_NumContrib_fMultREF_RawCalibComparison_" +
      tag_period + ".pdf");
  c_fNumContrib_fMultREF->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/"
      "MultCalibrationResult_NumContrib_fMultREF_RawCalibComparison_" +
      tag_period + ".json");

  //   TCanvas *c_fMultNTracksPV_fMultTPC_diff_fineBin =
  //       new TCanvas("c_fMultNTracksPV_fMultTPC_diff_fineBin",
  //                   "c_fMultNTracksPV_fMultTPC_diff_fineBin", 400, 400);
  //   c_fMultNTracksPV_fMultTPC_diff_fineBin->cd();
  //   gPad->SetLeftMargin(0.1);
  //   gPad->SetRightMargin(0.12);
  //   gPad->SetLogz();
  //   auto NumContribCalib_fMultTPC_clone2 =
  //       (TH2D *)NumContribCalib_fMultTPC->Clone(
  //           Form("NumContribCalib_fMultTPC_clone2_%d", GenerateUID()));
  //   // NumContribCalib_fMultTPC_clone->RebinX(2);
  //   // NumContribCalib_fMultTPC_clone->Draw("colz");
  //   for (int i = 1; i <= NumContribCalib_fMultTPC_clone2->GetNbinsX(); i++) {
  //     TH1D *h_proj = NumContribCalib_fMultTPC_clone2->ProjectionY(
  //         Form("h_proj_%d", i), i, i);
  //     h_proj->SetTitle(
  //         Form("NumContribCalibDiff_fMultTPC_clone;fMultTPC(%.0f<NumContribCalib<"
  //              "%.0f);Counts",
  //              NumContribCalib_fMultTPC_clone2->GetXaxis()->GetBinLowEdge(i),
  //              NumContribCalib_fMultTPC_clone2->GetXaxis()->GetBinUpEdge(i)));
  //     MRootGraphic::StyleHistCommonHist(h_proj);
  //     h_proj->Draw("hist");
  //     if (i == 1) {
  //       c_fMultNTracksPV_fMultTPC_diff_fineBin->SaveAs(
  //           "/home/szhu/work/alice/analysis/QA/plot/event/"
  //           "MultCalibrationResult_NumContribCalib_fMultTPC_Diff_fineBin_" +
  //               tag_period + ".pdf(",
  //           Form("Title: NumContribCalib%d", i));
  //     } else if (i == NumContribCalib_fMultTPC_clone2->GetNbinsX()) {
  //       c_fMultNTracksPV_fMultTPC_diff_fineBin->SaveAs(
  //           "/home/szhu/work/alice/analysis/QA/plot/event/"
  //           "MultCalibrationResult_NumContribCalib_fMultTPC_Diff_fineBin_" +
  //               tag_period + ".pdf)",
  //           Form("Title: NumContribCalib%d", i));
  //     } else {
  //       c_fMultNTracksPV_fMultTPC_diff_fineBin->SaveAs(
  //           "/home/szhu/work/alice/analysis/QA/plot/event/"
  //           "MultCalibrationResult_NumContribCalib_fMultTPC_Diff_fineBin_" +
  //               tag_period + ".pdf",
  //           Form("Title: NumContribCalib%d", i));
  //     }
  //   }

  //   TCanvas *c_fMultNTracksPV_fMultTPC_diff =
  //       new TCanvas("c_fMultNTracksPV_fMultTPC_diff",
  //                   "c_fMultNTracksPV_fMultTPC_diff", 400, 400);
  //   c_fMultNTracksPV_fMultTPC_diff->cd();
  //   gPad->SetLeftMargin(0.1);
  //   gPad->SetRightMargin(0.12);
  //   gPad->SetLogz();
  //   auto NumContribCalib_fMultTPC_clone =
  //       (TH2D *)NumContribCalib_fMultTPC->Clone();
  //   NumContribCalib_fMultTPC_clone->RebinX(3);
  //   // NumContribCalib_fMultTPC_clone->Draw("colz");
  //   vector<double> x_bins_upedge;
  //   vector<double> y_bins_upedge;
  //   int index_edge_up = 60;
  //   for (int i = 1; i <= NumContribCalib_fMultTPC_clone->GetNbinsX(); i++) {
  //     TH1D *h_proj =
  //         NumContribCalib_fMultTPC_clone->ProjectionY(Form("h_proj_%d", i),
  //         i, i);
  //     h_proj->SetTitle(
  //         Form("NumContribCalibDiff_fMultTPC_clone;fMultTPC(%.0f<NumContribCalib<"
  //              "%.0f);Counts",
  //              NumContribCalib_fMultTPC_clone->GetXaxis()->GetBinLowEdge(i),
  //              NumContribCalib_fMultTPC_clone->GetXaxis()->GetBinUpEdge(i)));
  //     MRootGraphic::StyleHistCommonHist(h_proj);
  //     double mean_fit = -100;
  //     double sigma_fit = 0;
  //     if (i <= index_edge_up) {
  //       // fit tow gaussian
  //       double mean1 =
  //           NumContribCalib_fMultTPC_clone->GetXaxis()->GetBinCenter(i);
  //       double sigma1 = 2.;
  //       TF1 *fit_func = new TF1(Form("fit_func_%d", i), "gaus", 0, 300);
  //       fit_func->SetParameters(h_proj->GetMaximum(), mean1, sigma1);
  //       h_proj->Fit(fit_func, "Q", "", 0, mean1 + 1.);
  //       double mean_new = fit_func->GetParameter(1);
  //       double sigma_new = fit_func->GetParameter(2);
  //       h_proj->Fit(fit_func, "Q", "", 0, mean_new + 1. * sigma_new);
  //       mean_fit = fit_func->GetParameter(1);
  //       sigma_fit = fit_func->GetParameter(2);
  //       x_bins_upedge.push_back(
  //           NumContribCalib_fMultTPC_clone->GetXaxis()->GetBinCenter(i));
  //       y_bins_upedge.push_back(mean_fit + 2.5 * sigma_fit);
  //     }
  //     h_proj->Draw("");
  //     double high = mean_fit + 2.5 * sigma_fit;
  //     TLine *line_high = new TLine(high, 0, high, h_proj->GetMaximum());
  //     line_high->SetLineColor(kGreen + 2);
  //     line_high->SetLineWidth(2);
  //     line_high->Draw("same");
  //     if (i == 1) {
  //       c_fMultNTracksPV_fMultTPC_diff->SaveAs(
  //           "/home/szhu/work/alice/analysis/QA/plot/event/"
  //           "MultCalibrationResult_NumContribCalib_fMultTPC_Diff_" +
  //               tag_period + ".pdf(",
  //           Form("Title: NumContribCalib%d", i));
  //     } else if (i == NumContribCalib_fMultTPC_clone->GetNbinsX()) {
  //       c_fMultNTracksPV_fMultTPC_diff->SaveAs(
  //           "/home/szhu/work/alice/analysis/QA/plot/event/"
  //           "MultCalibrationResult_NumContribCalib_fMultTPC_Diff_" +
  //               tag_period + ".pdf)",
  //           Form("Title: NumContribCalib%d", i));
  //     } else {
  //       c_fMultNTracksPV_fMultTPC_diff->SaveAs(
  //           "/home/szhu/work/alice/analysis/QA/plot/event/"
  //           "MultCalibrationResult_NumContribCalib_fMultTPC_Diff_" +
  //               tag_period + ".pdf",
  //           Form("Title: NumContribCalib%d", i));
  //     }
  //   }

  //   // remove the first two entries of x_bins_upedge and y_bins_upedge
  //   if (x_bins_upedge.size() > 2) {
  //     x_bins_upedge.erase(x_bins_upedge.begin(), x_bins_upedge.begin() + 2);
  //     y_bins_upedge.erase(y_bins_upedge.begin(), y_bins_upedge.begin() + 2);
  //   }

  //   TCanvas *c_fMultNTracksPV_fMultTPC_upedge =
  //       new TCanvas("c_fMultNTracksPV_fMultTPC_upedge",
  //                   "c_fMultNTracksPV_fMultTPC_upedge", 800, 400);
  //   c_fMultNTracksPV_fMultTPC_upedge->Divide(2, 1);
  //   c_fMultNTracksPV_fMultTPC_upedge->cd(1);
  //   gPad->SetLeftMargin(0.1);
  //   gPad->SetRightMargin(0.12);
  //   gPad->SetLogz();
  //   TGraph *graph_upedge =
  //       new TGraph(x_bins_upedge.size(), &x_bins_upedge[0],
  //       &y_bins_upedge[0]);
  //   NumContribCalib_fMultTPC->Draw("colz");
  //   graph_upedge->SetLineColor(kBlack);
  //   graph_upedge->SetLineWidth(2);
  //   graph_upedge->Draw("same");

  //   TLine *line_lowlimit =
  //       new TLine(6., 0, 6.,
  //       NumContribCalib_fMultTPC->GetYaxis()->GetXmax());
  //   line_lowlimit->SetLineColor(kRed);
  //   line_lowlimit->SetLineWidth(2);
  //   line_lowlimit->Draw("same");

  //   c_fMultNTracksPV_fMultTPC_upedge->cd(2);
  //   gPad->SetLeftMargin(0.1);
  //   gPad->SetRightMargin(0.12);
  //   gPad->SetLogz();
  //   TF1 *fit_func_upedge =
  //       new TF1("fit_func_upedge", "pol1",
  //       graph_upedge->GetXaxis()->GetXmin(),
  //               graph_upedge->GetXaxis()->GetXmax());
  //   ((TGraph *)graph_upedge->Clone())
  //       ->Fit(fit_func_upedge, "Q", "", graph_upedge->GetXaxis()->GetXmin(),
  //             graph_upedge->GetXaxis()->GetXmax());
  //   file_output->WriteObject<TF1>(fit_func_upedge, "fit_func_upedge");
  //   NumContribCalib_fMultTPC->Draw("colz");
  //   line_lowlimit->Draw("same");
  //   fit_func_upedge->SetLineColor(kBlue);
  //   fit_func_upedge->Draw("same");
  //   c_fMultNTracksPV_fMultTPC_upedge->cd(1);
  //   TLegend *legend_upedge = new TLegend(0.438649, 0.163482, 0.738878,
  //   0.315458); legend_upedge->SetMargin(0.3);
  //   // legend_upedge->SetHeader("MultTPC up edge");
  //   legend_upedge->AddEntry(graph_upedge, "#mu+2.5#sigma", "l");
  //   legend_upedge->AddEntry(fit_func_upedge, "pol1 fit", "l");
  //   legend_upedge->AddEntry(line_lowlimit, "NumContrib low limit", "l");
  //   legend_upedge->SetLineColor(0);
  //   legend_upedge->SetLineStyle(0);
  //   legend_upedge->SetFillStyle(0);
  //   legend_upedge->SetTextSize(0.04);
  //   legend_upedge->Draw();
  //   c_fMultNTracksPV_fMultTPC_upedge->SaveAs(
  //       "/home/szhu/work/alice/analysis/QA/plot/event/"
  //       "MultCalibrationResult_NumContribCalib_fMultTPC_UpEdge_" +
  //       tag_period + ".pdf");
  //   c_fMultNTracksPV_fMultTPC_upedge->SaveAs(
  //       "/home/szhu/work/alice/analysis/QA/plot/event/"
  //       "MultCalibrationResult_NumContribCalib_fMultTPC_UpEdge_" +
  //       tag_period + ".json");
  //   file_output->Close();
}