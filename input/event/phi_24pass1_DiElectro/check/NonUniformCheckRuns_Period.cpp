#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TLegend.h"
#include "TLegendEntry.h"

double euclideanDistance(TH1 *h1, TH1 *h2) {
  double distance = 0.0;
  for (int i = 1; i <= h1->GetNbinsX(); i++) {
    distance += TMath::Power(h1->GetBinContent(i) - h2->GetBinContent(i), 2);
  }
  return sqrt(distance);
}

void NonUniformCheckRuns_af_sub(TString name_tag, TString path_list_base,
                                TString path_list_period,
                                TString path_list_cluster) {
  vector<TString> vec_base = ReadListFromFile(path_list_base);
  vector<TString> vec_period = ReadListFromFile(path_list_period);
  vector<TString> vec_output = ReadListFromFile(path_list_cluster);

  vector<TString> vec_intersection =
      GetIntersectionOfTwoLists({vec_base, vec_period, vec_output});

  if (vec_intersection.size() == 0) {
    cout << "No intersection found for " << name_tag.Data() << endl;
    return;
  }
  cout << "tag: " << name_tag.Data() << endl;
  int index = 0;
  for (const auto &run : vec_intersection)
    cout << ", " << run.Data();
  cout << " " << endl;

  TCanvas *c_phi_period =
      new TCanvas(Form("c_phi_period_%s", name_tag.Data()),
                  Form("c_phi_period_%s", name_tag.Data()), 1600, 800);
  c_phi_period->Divide(2, 1);
  c_phi_period->cd(1);
  vector<TH1D *> vec_phi;

  for (auto &name_period : vec_intersection) {
    TString path_input = Form("../%s/phi.root",
                              name_period.Data());
    TFile *fInput = TFile::Open(path_input, "READ");

    TH1D *phi_integral = (TH1D*)fInput->Get("Phi");
    phi_integral->SetDirectory(0);
    fInput->Close();
    phi_integral->GetXaxis()->SetRangeUser(0.,2.*TMath::Pi());
    phi_integral->SetName(Form("hist_%d", GenerateUID()));
    vec_phi.push_back(phi_integral);
    phi_integral->Scale(1.0 / phi_integral->Integral());
    phi_integral->GetYaxis()->SetRangeUser(0, phi_integral->GetMaximum() * 2.);
    phi_integral->GetYaxis()->SetMaxDigits(1);
    phi_integral->GetYaxis()->SetTitle("Normalized counts");

    phi_integral->SetLineColor(GetColorIndice(index));
    if (index == 1) {
      phi_integral->DrawNormalized();
    } else {
      phi_integral->DrawNormalized("same");
    }
    index++;
  }

  TLegend *leg_phi_period = new TLegend(0.1, 0.6, 0.9, 0.9);
  leg_phi_period->SetNColumns(5);
  index = 0;
  for (auto &name_period : vec_intersection) {
    auto entry = leg_phi_period->AddEntry(Form("phi_integral;1"),
                                          name_period.Data(), "l");
    //  set the line color
    entry->SetLineColor(GetColorIndice(index));
    index++;
  }
  leg_phi_period->Draw();

  c_phi_period->cd(2);
  // draw ratio to phi_af
  auto baseline = (TH1D *)vec_phi[0]->Clone(Form("hist_%d", GenerateUID()));
  baseline->Scale(1. / baseline->Integral());
  auto baseline2 = (TH1D *)vec_phi[1]->Clone(Form("hist_%d", GenerateUID()));
  baseline2->Scale(1. / baseline2->Integral());
  for (int index = 0; index < vec_phi.size(); index++) {
    auto hist = (TH1D *)vec_phi[index]->Clone(Form("hist_%d", GenerateUID()));
    hist->Scale(1. / hist->Integral());
    cout << euclideanDistance(hist, baseline) << ", "
         << euclideanDistance(hist, baseline2) << endl;
    hist->Divide(baseline);
    hist->GetYaxis()->SetRangeUser(0, 2.);
    hist->GetYaxis()->SetMaxDigits(1);
    hist->GetYaxis()->SetTitle("Ratio to first run in list");
    hist->SetLineColor(GetColorIndice(index));
    if (index == 0) {
      hist->Draw();
    } else {
      hist->Draw("same");
    }
  }
  c_phi_period->SaveAs(Form("phi_integral_runs_%s.pdf", name_tag.Data()));
}

void NonUniformCheckRuns_Period() {
  TString path_list =
      "./list/";
  gStyle->SetOptStat(0);
  for (auto period : {"af", "ag", "aj", "al", "am", "an", "ao"}) {
    NonUniformCheckRuns_af_sub(
        "cluster1_" + TString(period),
        path_list + "local_24pass1",
        path_list + "period/" + period, path_list + "cluster_1_985_inRow.list");
    NonUniformCheckRuns_af_sub(
        "cluster2_" + TString(period),
        path_list + "local_24pass1",
        path_list + "period/" + period, path_list + "cluster_2_951_inRow.list");
    NonUniformCheckRuns_af_sub(
        "cluster3_" + TString(period),
        path_list + "local_24pass1",
        path_list + "period/" + period, path_list + "cluster_3_922_inRow.list");
    NonUniformCheckRuns_af_sub(
        "cluster4_" + TString(period),
        path_list + "local_24pass1",
        path_list + "period/" + period, path_list + "cluster_4_471_inRow.list");
    NonUniformCheckRuns_af_sub(
        "cluster5_part1_" + TString(period),
        path_list + "local_24pass1",
        path_list + "period/" + period, path_list + "cluster5_part1");
    NonUniformCheckRuns_af_sub(
        "cluster5_part2_" + TString(period),
        path_list + "local_24pass1",
        path_list + "period/" + period, path_list + "cluster5_part2");
  }
  const char *period = "ao";
  NonUniformCheckRuns_af_sub(
      "cluster4_" + TString(period),
      path_list + "run_24pass1_minbias_derived_inCluster1.list",
      path_list + "period/" + period, path_list + "cluster_4_471_inRow.list");
  // NonUniformCheckRuns_af_sub(
  //     "cluster1", path_list + "run_24pass1_minbias_derived_inCluster1.list",
  //     path_list + "period/af", path_list + "cluster_1_985_inRow.list");
  // NonUniformCheckRuns_af_sub(
  //     "cluster2", path_list + "run_24pass1_minbias_derived_inCluster1.list",
  //     path_list + "period/af", path_list + "cluster_2_951_inRow.list");
  // NonUniformCheckRuns_af_sub(
  //     "cluster3", path_list + "run_24pass1_minbias_derived_inCluster1.list",
  //     path_list + "period/af", path_list + "cluster_3_922_inRow.list");
  // NonUniformCheckRuns_af_sub(
  //     "cluster4", path_list + "run_24pass1_minbias_derived_inCluster1.list",
  //     path_list + "period/af", path_list + "cluster_4_471.list");
}
