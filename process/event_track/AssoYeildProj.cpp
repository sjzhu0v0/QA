#include "MHelper.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"
#include "yaml-cpp/yaml.h"

class AssocYieldHelperREF {
private:
public:
  MHnTool* fHnSame;
  MHnTool* fHnMix;
  MHnTool* fHnTrigger;
  bool doMixMultInt = false;

  AssocYieldHelperREF(MHnTool* hnSame, MHnTool* hnMix, MHnTool* hnTrigger) {
    fHnSame = hnSame;
    fHnMix = hnMix;
    fHnTrigger = hnTrigger;
  }

  void SetMixMultInt(bool doMixMultInt_ = true) { doMixMultInt = doMixMultInt_; }

  void Rebin(int dimTarget, int n) {
    switch (dimTarget) {

    case kDeltaEta:
      fHnSame->Rebin(0, n);
      fHnMix->Rebin(0, n);
      break;
    case kDeltaPhi:
      fHnSame->Rebin(1, n);
      fHnMix->Rebin(1, n);
      break;
    case kVtxZ:
      fHnSame->Rebin(2, n);
      fHnMix->Rebin(2, n);
      fHnTrigger->Rebin(0, n);
      break;
    case kNumContrib:
      fHnSame->Rebin(3, n);
      fHnMix->Rebin(3, n);
      fHnTrigger->Rebin(1, n);
      break;
    default:
      cerr << "Error: AssocYieldHelper_v1::Rebin: dimTarget is out of range" << endl;
      exit(1);
    }
  }

  void SetRangeUser(int dim, double min, double max) {
    if (dim < 0 || dim >= fHnSame->hN->GetNdimensions()) {
      cerr << "Error: AssocYieldHelperREF::SetRangeUser: dim is out of range" << endl;
      exit(1);
    }
    switch (dim) {
    case kDeltaEta:
      fHnSame->SetRangeUser(0, min, max);
      fHnMix->SetRangeUser(0, min, max);
      break;
    case kDeltaPhi:
      fHnSame->SetRangeUser(1, min, max);
      fHnMix->SetRangeUser(1, min, max);
      break;
    case kVtxZ:
      fHnSame->SetRangeUser(2, min, max);
      fHnMix->SetRangeUser(2, min, max);
      fHnTrigger->SetRangeUser(0, min, max);
      break;
    case kNumContrib:
      fHnSame->SetRangeUser(3, min, max);
      fHnMix->SetRangeUser(3, min, max);
      fHnTrigger->SetRangeUser(1, min, max);
      break;
    default:
      cerr << "Error: AssocYieldHelperREF::SetRangeUser: dim is out of range" << endl;
      exit(1);
    }
  }

  TH2D* AssociatedYieldVtxZ(int iVtxZ, int iMult, bool doNTrigScale = true) {
    TH2D* h2D = fHnSame->Project(gtype_vars::kDeltaPhi, gtype_vars::kDeltaEta, {iVtxZ, iMult});
    TH2D* h2DMix = fHnMix->Project(gtype_vars::kDeltaPhi, gtype_vars::kDeltaEta,
                                   {iVtxZ, doMixMultInt ? iMult : 0});

    vector<int> vec_idTrigger_new = {iVtxZ, iMult};
    double number_triggered = fHnTrigger->GetBinContent(vec_idTrigger_new);
    TH2D* h_assoYield = (TH2D*)h2D->Clone(Form("h_assoYield_%d_%d", iVtxZ, iMult));
    DensityHisto2DNoWeight(h2DMix);

    h_assoYield->Divide(h2DMix);
    // h_assoYield->Scale(1.0 / number_triggered);
    if (doNTrigScale) {
      if (number_triggered != 0) {
        h_assoYield->Scale(1.0 / number_triggered);
      } else {
        cerr << "AssocYieldHelperREF:AssociatedYieldVtxZ: Error: "
                "number_triggered is zero, not scaling the histogram"
             << endl;
        exit(1);
      }
    }

    h2DMix->Delete();
    h2D->Delete();
    return h_assoYield;
  }

  TH2D* AssociatedYieldVtxZPtSum(int iVtxZ, int iMult, bool doNTrigScale = true) {
    TH2D* h2D = fHnSame->Project(gtype_vars::kDeltaPhi, gtype_vars::kDeltaEta, {iVtxZ, iMult});
    TH2D* h2DMix = fHnMix->Project(gtype_vars::kDeltaPhi, gtype_vars::kDeltaEta,
                                   {iVtxZ, doMixMultInt ? iMult : 0});

    vector<int> vec_idTrigger_new = {iVtxZ, iMult};
    double number_triggered = fHnTrigger->GetBinContent(vec_idTrigger_new);

    TH2D* h_assoYield = (TH2D*)h2D->Clone(Form("h_assoYield_%d_%d", iVtxZ, iMult));
    AccCorrHisto2DNoWeight(h2DMix);

    h_assoYield->Divide(h2DMix);
    // h_assoYield->Scale(1.0 / number_triggered);
    if (doNTrigScale) {
      if (number_triggered != 0) {
        h_assoYield->Scale(1.0 / number_triggered);
      } else {
        cerr << "AssocYieldHelperREF:AssociatedYieldVtxZ: Error: "
                "number_triggered is zero, not scaling the histogram"
             << endl;
        exit(1);
      }
    }

    h2DMix->Delete();
    h2D->Delete();
    return h_assoYield;
  }

  TH2D* AssociatedYieldVtxZPtSum(int iVtxZ, vector<int> vec_iMult, bool doNTrigScale = true) {
    TH2D* h2D =
        fHnSame->Project(gtype_vars::kDeltaPhi, gtype_vars::kDeltaEta, {iVtxZ, vec_iMult[0]});
    for (int i = 1; i < vec_iMult.size(); i++) {
      TH2D* h2D_temp =
          fHnSame->Project(gtype_vars::kDeltaPhi, gtype_vars::kDeltaEta, {iVtxZ, vec_iMult[i]});
      h2D->Add(h2D_temp);
      h2D_temp->Delete();
    }
    TH2D* h2DMix = fHnMix->Project(gtype_vars::kDeltaPhi, gtype_vars::kDeltaEta,
                                   {iVtxZ, doMixMultInt ? vec_iMult[0] : 0});
    for (int i = 1; i < vec_iMult.size(); i++) {
      TH2D* h2DMix_temp = fHnMix->Project(gtype_vars::kDeltaPhi, gtype_vars::kDeltaEta,
                                          {iVtxZ, doMixMultInt ? vec_iMult[i] : 0});
      h2DMix->Add(h2DMix_temp);
      h2DMix_temp->Delete();
    }

    vector<int> vec_idTrigger_new = {iVtxZ, vec_iMult[0]};
    double number_triggered = fHnTrigger->GetBinContent(vec_idTrigger_new);
    for (int i = 1; i < vec_iMult.size(); i++) {
      vector<int> vec_idTrigger_new_temp = {iVtxZ, vec_iMult[i]};
      number_triggered += fHnTrigger->GetBinContent(vec_idTrigger_new_temp);
    }

    TH2D* h_assoYield = (TH2D*)h2D->Clone(Form("h_assoYield_%d_%d", iVtxZ, vec_iMult[0]));
    AccCorrHisto2DNoWeight(h2DMix);

    h_assoYield->Divide(h2DMix);
    // h_assoYield->Scale(1.0 / number_triggered);
    if (doNTrigScale) {
      if (number_triggered != 0) {
        h_assoYield->Scale(1.0 / number_triggered);
      } else {
        cerr << "AssocYieldHelperREF:AssociatedYieldVtxZ: Error: "
                "number_triggered is zero, not scaling the histogram"
             << endl;
        exit(1);
      }
    }

    h2DMix->Delete();
    h2D->Delete();
    return h_assoYield;
  }

  TH2D* AssociatedYieldPtSum(int iMult, bool doNTrigScale = true) {
    int nVtxZ = fHnSame->GetNbins(2);
    TH2D* h2_first = AssociatedYieldVtxZPtSum(1, iMult, false);
    // int iMass_new;
    // if (iMass != 0)
    //   iMass_new = fHnTrigger->hN->GetAxis(1)->FindBin(
    //       fHnMix->hN->GetAxis(3)->GetBinCenter(iMass));
    // else
    //   iMass_new = 0; // default to 1 if iMass is 0

    TH1D* h1_trigger = fHnTrigger->Project(0, {iMult});
    h1_trigger->SetName(Form("h1_trigger_%d", GenerateUID()));
    double sum_number_triggered = 0;
    for (int i = 1; i <= nVtxZ; i++) {
      double number_triggered = h1_trigger->GetBinContent(i);
      sum_number_triggered += number_triggered;
    }
    for (int i = 2; i <= nVtxZ; i++) {
      TH2D* h2D = AssociatedYieldVtxZPtSum(i, iMult, false);
      h2_first->Add(h2D);
      h2D->Delete();
    }
    if (sum_number_triggered == 0) {
      return h2_first;
    }
    // ScaleHisto2D(h2_first, 1.0 / sum_number_triggered);
    if (doNTrigScale) {
      if (sum_number_triggered != 0) {
        h2_first->Scale(1.0 / sum_number_triggered);
      } else {
        cerr << "AssocYieldHelperREF:AssociatedYield: Error: "
                "sum_number_triggered is zero, not scaling the histogram"
             << endl;
        exit(1);
      }
    }
    h1_trigger->Delete();
    return h2_first;
  }

  TH2D* AssociatedYieldPtSum(vector<int> vec_iMult, bool doNTrigScale = true) {
    int nVtxZ = fHnSame->GetNbins(2);
    TH2D* h2_first = AssociatedYieldVtxZPtSum(1, vec_iMult[0], false);
    // int iMass_new;
    // if (iMass != 0)
    //   iMass_new = fHnTrigger->hN->GetAxis(1)->FindBin(
    //       fHnMix->hN->GetAxis(3)->GetBinCenter(iMass));
    // else
    //   iMass_new = 0; // default to 1 if iMass is 0
    TH1D* h1_trigger = fHnTrigger->Project(0, {vec_iMult[0]});
    h1_trigger->SetName(Form("h1_trigger_%d", GenerateUID()));
    for (int i = 1; i < vec_iMult.size(); i++) {
      auto h1_trigger_temp = fHnTrigger->Project(0, {vec_iMult[i]});
      h1_trigger->Add(h1_trigger_temp);
      h1_trigger_temp->Delete();
    }

    double sum_number_triggered = 0;
    for (int i = 1; i <= nVtxZ; i++) {
      double number_triggered = h1_trigger->GetBinContent(i);
      sum_number_triggered += number_triggered;
    }
    for (int i = 2; i <= nVtxZ; i++) {
      TH2D* h2D = AssociatedYieldVtxZPtSum(i, vec_iMult, false);
      h2_first->Add(h2D);
      h2D->Delete();
    }
    if (sum_number_triggered == 0) {
      return h2_first;
    }
    // ScaleHisto2D(h2_first, 1.0 / sum_number_triggered);
    if (doNTrigScale) {
      if (sum_number_triggered != 0) {
        h2_first->Scale(1.0 / sum_number_triggered);
      } else {
        cerr << "AssocYieldHelperREF:AssociatedYield: Error: "
                "sum_number_triggered is zero, not scaling the histogram"
             << endl;
        exit(1);
      }
    }
    h1_trigger->Delete();
    return h2_first;
  }

  TH2D* AssociatedYield(int iMult, bool doNTrigScale = true) {
    int nVtxZ = fHnSame->GetNbins(2);
    TH2D* h2_first = AssociatedYieldVtxZ(1, iMult, false);
    // int iMass_new;
    // if (iMass != 0)
    //   iMass_new = fHnTrigger->hN->GetAxis(1)->FindBin(
    //       fHnMix->hN->GetAxis(3)->GetBinCenter(iMass));
    // else
    //   iMass_new = 0; // default to 1 if iMass is 0

    TH1D* h1_trigger = fHnTrigger->Project(0, {iMult});
    h1_trigger->SetName(Form("h1_trigger_%d", GenerateUID()));
    double sum_number_triggered = 0;
    for (int i = 1; i <= nVtxZ; i++) {
      double number_triggered = h1_trigger->GetBinContent(i);
      sum_number_triggered += number_triggered;
    }
    for (int i = 2; i <= nVtxZ; i++) {
      TH2D* h2D = AssociatedYieldVtxZ(i, iMult, false);
      h2_first->Add(h2D);
      h2D->Delete();
    }
    if (sum_number_triggered == 0) {
      return h2_first;
    }
    // ScaleHisto2D(h2_first, 1.0 / sum_number_triggered);
    if (doNTrigScale) {
      if (sum_number_triggered != 0) {
        h2_first->Scale(1.0 / sum_number_triggered);
      } else {
        cerr << "AssocYieldHelperREF:AssociatedYield: Error: "
                "sum_number_triggered is zero, not scaling the histogram"
             << endl;
        exit(1);
      }
    }
    h1_trigger->Delete();
    return h2_first;
  }
};

void AssoYieldProj(
    TString input_se_pr = "/home/szhu/work/alice/analysis/QA/test/JpsiAsso.root:DeltaEtaUS_"
                          "DeltaPhiUS_PosZUS_MassUS_PtUS_NumContribCalibUS",
    TString input_se_raw = "/home/szhu/work/alice/analysis/QA/test/JpsiQA.root:fPosZ_"
                           "MassUS_PtUS_NumContribCalib",
    TString input_me_pr = "/home/szhu/work/alice/analysis/QA/test/"
                          "MixEventReading.root:DeltaEtaUS_"
                          "DeltaPhiUS_PosZUS_MassUS_PtUS_NumContribCalibUS",
    TString path_output = "/home/szhu/work/alice/analysis/QA/test/"
                          "AssoYieldPt") {
  YAML::Node config = YAML::LoadFile("config.yaml");

  auto hist_se_pr = MRootIO::GetObjectDiectly<THnD>(input_se_pr);
  auto hist_se_raw = MRootIO::GetObjectDiectly<THnD>(input_se_raw);
  auto hist_me_pr = MRootIO::GetObjectDiectly<THnD>(input_me_pr);

  TFile* file_output = new TFile(Form("%s.root", path_output.Data()), "RECREATE");

  MHnTool hnTool_se_pr(hist_se_pr);
  MHnTool hnTool_se_raw(hist_se_raw);

  int nbins_posz_raw = hnTool_se_raw.GetNbins(0);
  if (nbins_posz_raw == 200)
    hnTool_se_raw.Rebin(0, 25); // Rebin PosZ
  else if (nbins_posz_raw == 8) {
    // Do nothing, already rebinned
  } else {
    cerr << "Error: Unexpected number of bins for PosZ: " << nbins_posz_raw << endl;
    exit(1);
  }
  MHnTool hnTool_me_pr(hist_me_pr);

  AssocYieldHelperREF assoYield(&hnTool_se_pr, &hnTool_me_pr, &hnTool_se_raw);
  assoYield.Rebin(gtype_vars::kNumContrib, 2);
  int n_rebin_deltaEta_assoYield = config["hist_binning"]["n_rebin_deltaEta_assoYield"].as<int>();
  assoYield.Rebin(gtype_vars::kDeltaEta, n_rebin_deltaEta_assoYield);

  hnTool_se_pr.PrintAllAxis();
  hnTool_se_raw.PrintAllAxis();
  hnTool_me_pr.PrintAllAxis();

  MRootGraphic::StyleCommon();
  gStyle->SetPalette(kRainBow);
  file_output->cd();
  auto h2_total_mass_pt = assoYield.AssociatedYieldPtSum(0, true);
  auto h2_lowMult_mass_pt = assoYield.AssociatedYieldPtSum({1, 2, 3, 4}, true);
  auto h2_highMult_mass_pt = assoYield.AssociatedYieldPtSum({5}, true);
  h2_total_mass_pt->SetName("h2_total");
  h2_lowMult_mass_pt->SetName("h2_lowMult");
  h2_highMult_mass_pt->SetName("h2_highMult");

  auto h2_sub_highLow_mass_pt = (TH2D*)h2_highMult_mass_pt->Clone("h2_sub_highLow_mass_pt");
  HistSubstraction2D(h2_sub_highLow_mass_pt, h2_highMult_mass_pt, h2_lowMult_mass_pt);
  h2_sub_highLow_mass_pt->SetName("h2_sub");

  h2_total_mass_pt->Write();
  h2_lowMult_mass_pt->Write();
  h2_highMult_mass_pt->Write();
  h2_sub_highLow_mass_pt->Write();

  file_output->Close();
}

int main(int argc, char** argv) {
  gROOT->SetBatch(kTRUE);
  TString path_input_se_pr = "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
                             "JpsiAsso_LHC22pass4_dqfilter.root:DeltaEtaUS_DeltaPhiUS_PosZUS_"
                             "MassUS_"
                             "PtUS_NumContribCalibUS";
  TString path_input_se_raw = "/home/szhu/work/alice/analysis/QA/input/jpsi/"
                              "JpsiQA_LHC22pass4_dqfilter.root:fPosZ_MassUS_PtUS_NumContribCalib";
  TString path_input_me_pr = "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
                             "MixEventReading_LHC22pass4_dqfilter.root:DeltaEtaUS_DeltaPhiUS_"
                             "PosZUS_"
                             "MassUS_PtUS_NumContribCalibUS";
  TString path_output = "/home/szhu/work/alice/analysis/QA/output/event_jpsi/AssoYieldQA";

  if (argc > 1) {
    path_input_se_pr = argv[1];
  }
  if (argc > 2) {
    path_input_se_raw = argv[2];
  }
  if (argc > 3) {
    path_input_me_pr = argv[3];
  }
  if (argc > 4) {
    path_output = argv[4];
  }
  AssoYieldProj(path_input_se_pr, path_input_se_raw, path_input_me_pr, path_output);

  return 0;
}
