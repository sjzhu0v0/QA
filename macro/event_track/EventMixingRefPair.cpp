#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MEventMixing.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "opt/EventData.h"
#include <ROOT/RDataFrame.hxx>
#include "TTreeReaderArray.h"



template <typename T> std::vector<T> makeVec(const TTreeReaderArray<T> &arr) {
  return std::vector<T>(arr.begin(), arr.end());
}
void EventMixingRef(TString path_input_flowVecd = "../input1.root",
                             TString path_input_mult = "../input2.root",
                             TString path_input_index = "input3.root",
                             TString path_output_tree = "output_mix.root") {
  TChain *tree_flowVecd =
      MRootIO::OpenChain(path_input_flowVecd, "O2dqflowvecd");
  TChain *tree_mult = MRootIO::OpenChain(path_input_mult, "MultCalib");
  TChain *tree_index = MRootIO::OpenChain(path_input_index, "EventMixing");
  tree_flowVecd->AddFriend(tree_mult);
  TTreeReader rEvt(tree_flowVecd);
  TTreeReader rPairs(tree_index);
  TTreeReaderValue<std::vector<std::pair<ULong64_t, ULong64_t>>> abPair(
      rPairs, "MixedEvent");
  TTreeReaderValue<double> NumContribCalib(rEvt, "NumContribCalib");
  TTreeReaderValue<int> fMultTPC(rEvt, "fMultTPC");
  TTreeReaderValue<int> fMultTracklets(rEvt, "fMultTracklets");
  TTreeReaderValue<int> fMultNTracksPV(rEvt, "fMultNTracksPV");
  TTreeReaderValue<float> fMultFT0C(rEvt, "fMultFT0C");
  TTreeReaderValue<float> fPosX(rEvt, "fPosX");
  TTreeReaderValue<float> fPosY(rEvt, "fPosY");
  TTreeReaderValue<float> fPosZ(rEvt, "fPosZ");
  TTreeReaderValue<ULong64_t> fSelection(rEvt, "fSelection");
  TTreeReaderValue<float> fHadronicRate(rEvt, "fHadronicRate");
  TTreeReaderValue<int> fEta_size(rEvt, "fEta_size");
  TTreeReaderArray<float> fPT(rEvt, "fPT");
  TTreeReaderArray<float> fEta(rEvt, "fEta");
  TTreeReaderArray<float> fPhi(rEvt, "fPhi");
  TTreeReaderArray<float> fMass(rEvt, "fMass");
  TTreeReaderArray<float> fSign(rEvt, "fSign");
  TTreeReaderArray<float> fPt1(rEvt, "fPt1");
  TTreeReaderArray<float> fEta1(rEvt, "fEta1");
  TTreeReaderArray<float> fPhi1(rEvt, "fPhi1");
  TTreeReaderArray<int> fSign1(rEvt, "fSign1");
  TTreeReaderArray<float> fITSChi2NCl1(rEvt, "fITSChi2NCl1");
  TTreeReaderArray<float> fTPCNClsCR1(rEvt, "fTPCNClsCR1");
  TTreeReaderArray<float> fTPCNClsFound1(rEvt, "fTPCNClsFound1");
  TTreeReaderArray<float> fTPCChi2NCl1(rEvt, "fTPCChi2NCl1");
  TTreeReaderArray<float> fTPCSignal1(rEvt, "fTPCSignal1");
  TTreeReaderArray<float> fTPCNSigmaEl1(rEvt, "fTPCNSigmaEl1");
  TTreeReaderArray<float> fTPCNSigmaPi1(rEvt, "fTPCNSigmaPi1");
  TTreeReaderArray<float> fTPCNSigmaPr1(rEvt, "fTPCNSigmaPr1");
  TTreeReaderArray<float> fPt2(rEvt, "fPt2");
  TTreeReaderArray<float> fEta2(rEvt, "fEta2");
  TTreeReaderArray<float> fPhi2(rEvt, "fPhi2");
  TTreeReaderArray<int> fSign2(rEvt, "fSign2");
  TTreeReaderArray<float> fITSChi2NCl2(rEvt, "fITSChi2NCl2");
  TTreeReaderArray<float> fTPCNClsCR2(rEvt, "fTPCNClsCR2");
  TTreeReaderArray<float> fTPCNClsFound2(rEvt, "fTPCNClsFound2");
  TTreeReaderArray<float> fTPCChi2NCl2(rEvt, "fTPCChi2NCl2");
  TTreeReaderArray<float> fTPCSignal2(rEvt, "fTPCSignal2");
  TTreeReaderArray<float> fTPCNSigmaEl2(rEvt, "fTPCNSigmaEl2");
  TTreeReaderArray<float> fTPCNSigmaPi2(rEvt, "fTPCNSigmaPi2");
  TTreeReaderArray<float> fTPCNSigmaPr2(rEvt, "fTPCNSigmaPr2");
  TTreeReaderValue<int> fPTREF_size(rEvt, "fPTREF_size");
  TTreeReaderArray<float> fPTREF(rEvt, "fPTREF");
  TTreeReaderArray<float> fEtaREF(rEvt, "fEtaREF");
  TTreeReaderArray<float> fPhiREF(rEvt, "fPhiREF");
  TTreeReaderArray<float> fITSChi2NCl_ref(rEvt, "fITSChi2NCl");
  TTreeReaderArray<float> fTPCNClsCR_ref(rEvt, "fTPCNClsCR");
  TTreeReaderArray<float> fTPCNClsFound_ref(rEvt, "fTPCNClsFound");
  TTreeReaderArray<float> fTPCChi2NCl_ref(rEvt, "fTPCChi2NCl");
  TTreeReaderArray<float> fTPCSignal_ref(rEvt, "fTPCSignal");
  TTreeReaderArray<float> fTPCNSigmaEl_ref(rEvt, "fTPCNSigmaEl");
  TTreeReaderArray<float> fTPCNSigmaPi_ref(rEvt, "fTPCNSigmaPi");
  TTreeReaderArray<float> fTPCNSigmaPr_ref(rEvt, "fTPCNSigmaPr");

  TFile fout(path_output_tree, "RECREATE");
  TTree out("jpsi_ref_pairs", "mixed jpsi(A) x ref(B) pairs");

  double o_NumContribCalib;
  int o_fMultTPC, o_fMultTracklets, o_fMultNTracksPV;
  float o_fMultFT0C;
  float o_fPosX, o_fPosY, o_fPosZ;
  int o_fSelection;
  float o_fHadronicRate;
  out.Branch("fMultTPC", &o_fMultTPC);
  out.Branch("fMultTracklets", &o_fMultTracklets);
  out.Branch("fMultNTracksPV", &o_fMultNTracksPV);
  out.Branch("fMultFT0C", &o_fMultFT0C);
  out.Branch("fPosX", &o_fPosX);
  out.Branch("fPosY", &o_fPosY);
  out.Branch("fPosZ", &o_fPosZ);
  out.Branch("fSelection", &o_fSelection);
  out.Branch("fHadronicRate", &o_fHadronicRate);
  float o_ref1_pt, o_ref1_eta, o_ref1_phi;
  float o_ref1_ITSChi2NCl, o_ref1_TPCNClsCR, o_ref1_TPCNClsFound;
  float o_ref1_TPCChi2NCl, o_ref1_TPCSignal;
  float o_ref1_nsig_el, o_ref1_nsig_pi, o_ref1_nsig_pr;
  out.Branch("NumContribCalib", &o_NumContribCalib);
  out.Branch("ref1_pt", &o_ref1_pt);
  out.Branch("ref1_eta", &o_ref1_eta);
  out.Branch("ref1_phi", &o_ref1_phi);
  out.Branch("ref1_ITSChi2NCl", &o_ref1_ITSChi2NCl);
  out.Branch("ref1_TPCNClsCR", &o_ref1_TPCNClsCR);
  out.Branch("ref1_TPCNClsFound", &o_ref1_TPCNClsFound);
  out.Branch("ref1_TPCChi2NCl", &o_ref1_TPCChi2NCl);
  out.Branch("ref1_TPCSignal", &o_ref1_TPCSignal);
  out.Branch("ref1_nsig_el", &o_ref1_nsig_el);
  out.Branch("ref1_nsig_pi", &o_ref1_nsig_pi);
  out.Branch("ref1_nsig_pr", &o_ref1_nsig_pr);
  float o_ref2_pt, o_ref2_eta, o_ref2_phi;
  float o_ref2_ITSChi2NCl, o_ref2_TPCNClsCR, o_ref2_TPCNClsFound;
  float o_ref2_TPCChi2NCl, o_ref2_TPCSignal;
  float o_ref2_nsig_el, o_ref2_nsig_pi, o_ref2_nsig_pr;
  out.Branch("ref2_pt", &o_ref2_pt);
  out.Branch("ref2_eta", &o_ref2_eta);
  out.Branch("ref2_phi", &o_ref2_phi);
  out.Branch("ref2_ITSChi2NCl", &o_ref2_ITSChi2NCl);
  out.Branch("ref2_TPCNClsCR", &o_ref2_TPCNClsCR);
  out.Branch("ref2_TPCNClsFound", &o_ref2_TPCNClsFound);
  out.Branch("ref2_TPCChi2NCl", &o_ref2_TPCChi2NCl);
  out.Branch("ref2_TPCSignal", &o_ref2_TPCSignal);
  out.Branch("ref2_nsig_el", &o_ref2_nsig_el);
  out.Branch("ref2_nsig_pi", &o_ref2_nsig_pi);
  out.Branch("ref2_nsig_pr", &o_ref2_nsig_pr);
  long long nWritten = 0;

  bool isInteractive = is_interactive();
  long long nEntries = rPairs.GetEntries();

  long long iEntry = -1;
  while (rPairs.Next()) {
    rPairs.SetEntry(iEntry);

    for (const auto &abPair_single : *abPair) {
      ULong64_t entryA = abPair_single.first;
      ULong64_t entryB = abPair_single.second;
      int nPairs = fPT.GetSize() * fPTREF.GetSize();
      if (nPairs == 0) {
        continue;
      }
      rEvt.SetEntry(entryA);
      o_NumContribCalib = *NumContribCalib;
      o_fMultTPC = *fMultTPC;
      o_fMultTracklets = *fMultTracklets;
      o_fMultNTracksPV = *fMultNTracksPV;
      o_fMultFT0C = *fMultFT0C;
      o_fPosX = *fPosX;
      o_fPosY = *fPosY;
      o_fPosZ = *fPosZ;
      o_fSelection = *fSelection;
      o_fHadronicRate = *fHadronicRate;
      auto v_ref1_pt = makeVec(fPTREF);
      auto v_ref1_eta = makeVec(fEtaREF);
      auto v_ref1_phi = makeVec(fPhiREF);
      auto v_ref1_its = makeVec(fITSChi2NCl_ref);
      auto v_ref1_cr = makeVec(fTPCNClsCR_ref);
      auto v_ref1_found = makeVec(fTPCNClsFound_ref);
      auto v_ref1_chi2 = makeVec(fTPCChi2NCl_ref);
      auto v_ref1_sig = makeVec(fTPCSignal_ref);
      auto v_ref1_nel = makeVec(fTPCNSigmaEl_ref);
      auto v_ref1_npi = makeVec(fTPCNSigmaPi_ref);
      auto v_ref1_npr = makeVec(fTPCNSigmaPr_ref);


      rEvt.SetEntry(entryB);
      auto v_ref2_pt = makeVec(fPTREF);
      auto v_ref2_eta = makeVec(fEtaREF);
      auto v_ref2_phi = makeVec(fPhiREF);
      auto v_ref2_its = makeVec(fITSChi2NCl_ref);
      auto v_ref2_cr = makeVec(fTPCNClsCR_ref);
      auto v_ref2_found = makeVec(fTPCNClsFound_ref);
      auto v_ref2_chi2 = makeVec(fTPCChi2NCl_ref);
      auto v_ref2_sig = makeVec(fTPCSignal_ref);
      auto v_ref2_nel = makeVec(fTPCNSigmaEl_ref);
      auto v_ref2_npi = makeVec(fTPCNSigmaPi_ref);
      auto v_ref2_npr = makeVec(fTPCNSigmaPr_ref);

      // long long filled = 0;
      for (size_t iref1 = 0; iref1 < v_ref1_pt.size(); ++iref1) {
        o_ref1_pt = v_ref1_pt[iref1];
        o_ref1_eta = v_ref1_eta[iref1];
        o_ref1_phi = v_ref1_phi[iref1];
        o_ref1_ITSChi2NCl = v_ref1_its[iref1];
        o_ref1_TPCNClsCR = v_ref1_cr[iref1];
        o_ref1_TPCNClsFound = v_ref1_found[iref1];
        o_ref1_TPCChi2NCl = v_ref1_chi2[iref1];
        o_ref1_TPCSignal = v_ref1_sig[iref1];
        o_ref1_nsig_el = v_ref1_nel[iref1];
        o_ref1_nsig_pi = v_ref1_npi[iref1];
        o_ref1_nsig_pr = v_ref1_npr[iref1];
        for (size_t iref2 = 0; iref2 < v_ref2_pt.size(); ++iref2) {
          o_ref2_pt = v_ref2_pt[iref2];
          o_ref2_eta = v_ref2_eta[iref2];
          o_ref2_phi = v_ref2_phi[iref2];
          o_ref2_ITSChi2NCl = v_ref2_its[iref2];
          o_ref2_TPCNClsCR = v_ref2_cr[iref2];
          o_ref2_TPCNClsFound = v_ref2_found[iref2];
          o_ref2_TPCChi2NCl = v_ref2_chi2[iref2];
          o_ref2_TPCSignal = v_ref2_sig[iref2];
          o_ref2_nsig_el = v_ref2_nel[iref2];
          o_ref2_nsig_pi = v_ref2_npi[iref2];
          o_ref2_nsig_pr = v_ref2_npr[iref2];
          out.Fill();
        }
      }
    }
  }
  out.Write();
  fout.Close();
}



int main(int argc, char **argv) {
  TString path_input_flowVecd = "../input.root";
  TString path_input_mult = "../input2.root";
  TString path_output = "output.root";
  TString path_output_tree = "output_tree.root";

  if (argc > 1) {
    path_input_flowVecd = argv[1];
  }
  if (argc > 2) {
    path_input_mult = argv[2];
  }
  if (argc > 3) {
    path_output = argv[3];
  }
  if (argc > 4) {
    path_output_tree = argv[4];
  }

  gROOT->SetBatch(kTRUE); // Disable interactive graphics
  EventMixingRef(path_input_flowVecd, path_input_mult, path_output,
                 path_output_tree);

  return 0;
}
