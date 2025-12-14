#include "MALICE.h"
#include "MHead.h"
#include "MRootGraphic.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "fstream"
#include "iostream"
TLatex *gLatex;
TLatex *gLatex2;

int gIndex = 0;

void GetLatex(TGraph *gr, double x) {
  double y = gr->Eval(x);
  TString content = TString::Format("%.2f%%", y * 100);
  gLatex->DrawLatex(x - 0.0002, y + 0.015, content);
}

void draw_proportion_main_vs_threshold() {
  TFile *file =
      new TFile("/home/szhu/work/alice/analysis/QA/plot/event/grouping/"
                "proportion_main_vs_threshold_distance_fine.root");
  gLatex = new TLatex();
  gLatex->SetTextSize(0.03);
  gLatex->SetTextAlign(12);
  gLatex->SetTextColor(kRed);
  gLatex2 = new TLatex();
  gLatex2->SetTextSize(0.03);
  gLatex2->SetTextAlign(12);

  auto gr = (TGraph *)file->Get("Graph");
  gr->SetTitle("Proportion of the largest group vs threshold "
               "dissimilarity;Dissimilarity threshold;"
               "Proportion of the largest group");
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
  c1->cd();
  MRootGraphic::StyleCommon();
  MRootGraphic::StyleHistCommon(gr);
  gr->Draw("ALP");
  gPad->SetGrid();
  gr->GetXaxis()->SetMaxDigits(1);

  // 0.0078 0.00625 0.00515 0.00485
  TLine *line1 = new TLine(0.00485, gr->Eval(0.00485) - 0.1, 0.00485,
                           gr->Eval(0.00485) + 0.1);
  gLatex2->DrawLatex(0.00485, gr->Eval(0.00485) - 0.1, "Thres 1");
  line1->SetLineColor(kRed);
  line1->SetLineStyle(2);
  line1->SetLineWidth(2);
  line1->Draw("same");

  TLine *line2 = new TLine(0.00515, gr->Eval(0.00515) - 0.1, 0.00515,
                           gr->Eval(0.00515) + 0.1);
  gLatex2->DrawLatex(0.00515, gr->Eval(0.00515) - 0.1, "Thres 2");
  line2->SetLineColor(kRed);
  line2->SetLineStyle(2);
  line2->SetLineWidth(2);
  line2->Draw("same");

  TLine *line3 = new TLine(0.00625, gr->Eval(0.00625) - 0.1, 0.00625,
                           gr->Eval(0.00625) + 0.1);
  gLatex2->DrawLatex(0.00625, gr->Eval(0.00625) - 0.1, "Thres 3");
  line3->SetLineColor(kRed);
  line3->SetLineStyle(2);
  line3->SetLineWidth(2);
  line3->Draw("same");

  TLine *line4 = new TLine(0.00785, gr->Eval(0.00785) - 0.1, 0.00785,
                           gr->Eval(0.00785) + 0.1);
  gLatex2->DrawLatex(0.00785, gr->Eval(0.00785) - 0.1, "Thres 4");
  line4->SetLineColor(kRed);
  line4->SetLineStyle(2);
  line4->SetLineWidth(2);
  line4->Draw("same");

  GetLatex(gr, 0.0049);
  GetLatex(gr, 0.0056);
  GetLatex(gr, 0.00666);
  GetLatex(gr, 0.0082);

  gLatex->DrawLatex(0.012, 0.25,
                    "Chosen threshold dissimilarity values for grouping");

  c1->SaveAs("/home/szhu/work/alice/analysis/QA/plot/event/grouping/"
             "proportion_main_vs_threshold_distance_fine.pdf");
}