#include "MALICE.h"
#include "MHead.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"
#include "TColor.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TStyle.h"

void StyleCommon1D(TH1 *hist_mb, Color_t color_line = kBlack,
                   Style_t style_marker = 1) {
  hist_mb->SetLineColor(color_line);
  hist_mb->SetMarkerStyle(style_marker);
  hist_mb->GetXaxis()->SetLabelSize(0.03);
  hist_mb->GetXaxis()->SetTitleSize(0.04);
  hist_mb->GetXaxis()->SetTitleOffset(1.2);
  hist_mb->GetYaxis()->SetLabelSize(0.03);
  hist_mb->GetYaxis()->SetTitleSize(0.04);
  hist_mb->GetYaxis()->SetTitleOffset(1.2);
  hist_mb->GetXaxis()->SetNdivisions(505);
  hist_mb->GetXaxis()->SetTickLength(0.02);
  hist_mb->GetYaxis()->SetNdivisions(10);
  hist_mb->GetYaxis()->SetTickLength(0.02);
  hist_mb->GetXaxis()->SetLabelOffset(0.01);
  hist_mb->GetYaxis()->SetLabelOffset(0.01);
  hist_mb->GetXaxis()->SetLabelFont(42);
  hist_mb->GetYaxis()->SetLabelFont(42);
  hist_mb->GetXaxis()->SetTitleFont(42);
  hist_mb->GetYaxis()->SetTitleFont(42);
  hist_mb->GetXaxis()->SetTickLength(0.02);
  hist_mb->GetYaxis()->SetTickLength(0.02);
  hist_mb->GetXaxis()->CenterTitle();
  hist_mb->GetYaxis()->CenterTitle();
}

TCanvas *EventSelectionDraw(
    TString path_hist = "/home/szhu/work/alice/analysis/MultStudy/output/"
                        "TriggerStudy.root:map_trigger") {
  TH1D *hist = MRootIO::GetTH1D(path_hist);
  gStyle->SetOptStat(0);
  // hist->Draw();
  // hist->Print();
  // hist->GetXaxis()->Print();
  // cout << "Number of bins: " << hist->GetXaxis()->GetNbins() << endl;
  // cout << "maximum bin: " << hist->GetXaxis()->GetXmin() << endl;
  // cout << "minimum bin: " << hist->GetXaxis()->GetXmax() << endl;

  int n_bins_not_empty = 0;
  double max_value = hist->GetBinContent(51);
  int n_bins_max = 0;

  vector<TString> vec_name;
  vector<int> bins;

  for (int i = 1; i < MALICE::kNsel; i++) {
    if (hist->GetBinContent(i) > 0) {
      n_bins_not_empty++;
      bins.push_back(i);
      vec_name.push_back(MALICE::name_eventselection[i - 1]);
    }
  }

  // for (auto i : bins) {
  //   cout << i << endl;
  // }

  // for (int i = 1; i <= MALICE::kNsel; i++) {
  //   if (hist->GetBinContent(i) == 0)
  //     continue;
  //   if (hist->GetBinContent(i) == max_value || hist->GetBinContent(i) < 5000)
  //   {
  //     n_bins_max++;
  //     // cout << MALICE::name_eventselection[i - 1] << endl;
  //   }
  // }
  // for (int i = 1; i <= hist->GetNbinsX(); i++) {
  //   double value = hist->GetBinContent(i);
  //   cout << value << endl;
  //   if (value != max_value && value > 5000) {
  //     vec_name.push_back(MALICE::name_eventselection[i - 1]);
  //     bins.push_back(i);
  //   }
  // }

  // cout << "Number of bins not empty: " << n_bins_not_empty << endl;
  // cout << "Max value: " << max_value << endl;
  // cout << "Number of bins with max value: " << n_bins_max << endl;

  // n_bins_not_empty -=
  //     n_bins_max; // remove the bins with max value from the count of not
  //     empty

  TH1D *hist2 =
      new TH1D(Form("hist%d", GenerateUID()), "hist2", n_bins_not_empty, 0, 1);
  for (int i = 0; i < n_bins_not_empty; i++) {
    hist2->GetXaxis()->SetBinLabel(i + 1, vec_name[i]);
    hist2->SetBinContent(i + 1, hist->GetBinContent(bins[i]) / max_value);
  }

  TCanvas *c1 = new TCanvas(Form("canvas%d", GenerateUID()),
                            Form("canvas%d", GenerateUID()), 1200, 600);
  c1->SetMargin(0.1, 0.1, 0.15, 0.05);
  hist2->GetYaxis()->SetRangeUser(0, 1.2);
  hist2->Draw("hist text");
  // TLatex *latex = new TLatex();
  // latex->SetTextSize(0.03);
  // latex->SetTextAlign(12);
  // latex->SetTextColor(kBlack);
  // latex->SetTextFont(42);
  // latex->DrawLatexNDC(0.1, 0.95, Form("Run %d", 526641));
  cout << n_bins_not_empty << endl;
  return c1;
}

template <typename T> T *GetObjectFromCanvas(TCanvas *c, TString name) {
  TList *list = c->GetListOfPrimitives();
  TIter next(list);
  TObject *obj;
  while ((obj = next())) {
    if (obj->InheritsFrom(name)) {
      return (T *)obj;
    }
  }
  cerr << "GetObjectFromCanvas: Object with name " << name
       << " does not inherit from the specified class." << endl;
  exit(1);
  return nullptr;
}

template <typename T> T *GetObjectFromCanvas(TCanvas *c, int index) {
  TList *list = c->GetListOfPrimitives();
  TObject *obj = list->At(index);
  if (obj->InheritsFrom(T::Class())) {
    return (T *)obj;
  }
  cerr << "GetObjectFromCanvas: Object at index " << index
       << " does not inherit from the specified class." << endl;
  exit(1);
  return nullptr;
}

void Comparison_mb_dqfilter(
    TString path_hist_mb = "/home/szhu/work/alice/analysis/QA/output/event/"
                           "triggerStudy_LHC22o_pass4_thin_526641_mb.root",
    TString path_hist_dqfilter =
        "/home/szhu/work/alice/analysis/QA/output/event/"
        "triggerStudy_LHC22o_pass4_thin_526641_dqfilter.root") {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(1);
  gStyle->SetNumberContours(255);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);

  /* #region :event selection comparison */
  TCanvas *c_mb = EventSelectionDraw(path_hist_mb + ":map_trigger");
  c_mb->SaveAs("/home/szhu/work/alice/analysis/QA/plot/event/"
               "Comparison_mb_dqfilter_effectiveTrigger_mb.pdf");
  c_mb->SaveAs("/home/szhu/work/alice/analysis/QA/plot/event/"
               "Comparison_mb_dqfilter_effectiveTrigger_mb.json");
  TCanvas *c_dqfilter = EventSelectionDraw(path_hist_dqfilter + ":map_trigger");
  c_dqfilter->SaveAs("/home/szhu/work/alice/analysis/QA/plot/event/"
                     "Comparison_mb_dqfilter_effectiveTrigger_dqfilter.pdf");
  c_dqfilter->SaveAs("/home/szhu/work/alice/analysis/QA/plot/event/"
                     "Comparison_mb_dqfilter_effectiveTrigger_dqfilter.json");
  // TH1D *hist_mb = GetObjectFromCanvas<TH1D>(c_mb, 1);
  // TH1D *hist_dqfilter = GetObjectFromCanvas<TH1D>(c_dqfilter, 0);
  // TCanvas *c_mb_dqfilter =
  //     new TCanvas(Form("canvas%d", GenerateUID()),
  //                 Form("canvas%d", GenerateUID()), 800, 600);
  // hist_mb->SetLineColor(kRed);
  // hist_mb->SetLineWidth(2);
  // hist_mb->SetTitle("Event Selection Comparison");
  // hist_dqfilter->SetLineColor(kBlue);
  // hist_dqfilter->SetLineWidth(2);

  // hist_mb->GetYaxis()->SetTitle("Normalized Counts");
  // hist_mb->GetXaxis()->SetTitle("Event Selection");
  // StyleCommon1D(hist_mb);
  // StyleCommon1D(hist_dqfilter);

  // hist_mb->Draw("hist text");
  // hist_dqfilter->Draw("hist text same");

  TCanvas *c_mb_tvx =
      EventSelectionDraw(path_hist_mb + ":map_trigger_isTriggerTVX");
  TCanvas *c_dqfilter_tvx =
      EventSelectionDraw(path_hist_dqfilter + ":map_trigger_isTriggerTVX");
  /* #endregion */

  TH1D *fNumContrib_mb = MRootIO::GetTH1D(path_hist_mb + ":fNumContrib");
  TH1D *fNumContrib_dqfilter =
      MRootIO::GetTH1D(path_hist_dqfilter + ":fNumContrib");

  TCanvas *c_mb_dqfilter =
      new TCanvas(Form("canvas%d", GenerateUID()),
                  Form("canvas%d", GenerateUID()), 800, 600);
  StyleCommon1D(fNumContrib_mb, kRed);
  StyleCommon1D(fNumContrib_dqfilter, kBlue);
  gPad->SetLogy();
  fNumContrib_mb->Draw();
  fNumContrib_dqfilter->Draw("same");
  fNumContrib_mb->GetXaxis()->SetTitle("Number of Vtx Contributors");
  fNumContrib_mb->GetYaxis()->SetTitle("Normalized Counts (a.u.)");
  fNumContrib_mb->Scale(fNumContrib_dqfilter->Integral() /
                        fNumContrib_mb->Integral());
  TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9, "NumContrib Comparison");
  leg->AddEntry(fNumContrib_mb, "minbias thin data", "l");
  leg->AddEntry(fNumContrib_dqfilter, "DQfiltered thin data", "l");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->Draw("same");
  c_mb_dqfilter->SaveAs("/home/szhu/work/alice/analysis/QA/plot/event/"
                        "Comparison_mb_dqfilter_NumContributor.pdf");
  c_mb_dqfilter->SaveAs("/home/szhu/work/alice/analysis/QA/plot/event/"
                        "Comparison_mb_dqfilter_NumContributor.json");
}
