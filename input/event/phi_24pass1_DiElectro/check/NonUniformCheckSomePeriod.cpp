#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TLegend.h"
#include "TLegendEntry.h"

void NonUniformCheckSomePeriod_sub(TString name_cluster = "2") {
  // vector<TString> vec_name_period = {"LHC16d", "LHC16e", "LHC16f", "LHC16g",
  // "LHC16h"}; af  ag  aj  al  am  an  ao
  vector<TString> vec_name_period = {"af", "ag", "aj", "al", "am", "an", "ao"};
  TCanvas* c_phi_period = new TCanvas("c_phi_period", "c_phi_period", 1600, 800);
  c_phi_period->Divide(2, 1);
  c_phi_period->cd(1);
  vector<TH1D*> vec_phi;
  for (auto& name_period : vec_name_period) {
    TString path_input = Form("./merge/%s_%s.root", name_cluster.Data(), name_period.Data());
    TFile* fInput = TFile::Open(path_input, "READ");
    if (!fInput || fInput->IsZombie()) {
      cout << "Cannot open file: " << path_input << endl;
      vec_phi.push_back(nullptr);
      continue;
    }
    auto phi_integral = (TH1D*)fInput->Get("Phi");
    phi_integral->SetDirectory(0);
    fInput->Close();
    phi_integral->GetXaxis()->SetRangeUser(0., 2 * TMath::Pi());
    // MRootGraphic::StyleHistCommon(fPTREF_fPhiREF);
    MRootGraphic::StyleHistCommon(phi_integral);
    phi_integral->SetName(Form("hist_%d", GenerateUID()));
    vec_phi.push_back(phi_integral);
    phi_integral->Scale(1.0 / phi_integral->Integral());
    phi_integral->GetYaxis()->SetRangeUser(0, phi_integral->GetMaximum() * 2.);
    phi_integral->GetYaxis()->SetMaxDigits(1);
    phi_integral->GetYaxis()->SetTitle("Normalized counts");
    // phi_integral->DrawNormalized();
    static int index = 1;
    // phi_integral->SetLineWidth(9 - index);
    phi_integral->SetLineColor(index + 1);
    if (index == 1) {
      phi_integral->SetLineColor(index);
      phi_integral->DrawNormalized();
    } else {
      phi_integral->SetLineColor(index);
      phi_integral->DrawNormalized("same");
    }
    index++;
  }

  TLegend* leg_phi_period = new TLegend(0.4, 0.6, 0.9, 0.9);
  int index = 1;
  for (auto& name_period : vec_name_period) {
    if (vec_phi[index - 1] == nullptr) {
      index++;
      continue;
    }
    auto entry =
        leg_phi_period->AddEntry(Form("phi_integral;1"), Form("LHC24%s", name_period.Data()), "l");
    //  set the line color
    entry->SetLineColor(index);
    index++;
  }
  leg_phi_period->Draw();

  c_phi_period->cd(2);
  // draw ratio to phi_af
  TH1D* baseline;
  TString baseline_name;
  int index_baseline = 0;
  for (auto& hist : vec_phi) {
    if (hist != nullptr) {
      baseline = (TH1D*)hist->Clone(Form("baseline_%d", GenerateUID()));
      baseline_name = vec_name_period[index_baseline];
      break;
    }
    index_baseline++;
  }
  baseline->Scale(1. / baseline->Integral());
  // vec_phi[0]->Scale(1. / vec_phi[0]->Integral());
  for (int index = 0; index < vec_phi.size(); index++) {
    if (vec_phi[index] == nullptr) {
      continue;
    }
    auto hist = (TH1D*)vec_phi[index]->Clone(Form("hist_%d", GenerateUID()));
    hist->Scale(1. / hist->Integral());
    hist->Divide(baseline);
    // hist->Scale(1./hist->Integral());
    hist->GetYaxis()->SetRangeUser(0, 2.);
    hist->GetYaxis()->SetMaxDigits(1);
    hist->GetYaxis()->SetTitle(Form("Ratio to LHC24%s", baseline_name.Data()));
    hist->SetTitle(Form("Ratio to LHC24%s", baseline_name.Data()));
    // hist->SetLineWidth(9 - index);
    hist->SetLineColor(index + 1);
    if (index == 0) {
      hist->Draw();
    } else {
      hist->Draw("same");
      cout << hist->GetName() << endl;
    }
  }
  // c_phi_period->SaveAs("phi_integral_periods_cluster2.pdf");
  c_phi_period->SaveAs(Form("phi_integral_periods_cluster%s.pdf", name_cluster.Data()));
}

void NonUniformCheckSomePeriod() {
  gStyle->SetOptStat(0);
  for (auto cluster : {"1", "2", "3", "4"}) {
    NonUniformCheckSomePeriod_sub(cluster);
  }
}
