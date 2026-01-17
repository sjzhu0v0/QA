#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MEventMixing.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "TTreeReaderArray.h"
#include "opt/EventData.h"
#include <ROOT/RDataFrame.hxx>

template <typename T> std::vector<T> makeVec(const TTreeReaderArray<T>& arr) {
  return std::vector<T>(arr.begin(), arr.end());
}
void EventMixingRef(TString path_input_flowVecd = "../input1.root",
                    TString path_input_mult = "../input2.root",
                    TString path_input_index = "input3.root",
                    TString path_output_tree = "output_mix.root") {
  TChain* tree_flowVecd = MRootIO::OpenChain(path_input_flowVecd, "O2dqflowvecd");
  TChain* tree_flowVecd2 = MRootIO::OpenChain(path_input_flowVecd, "O2dqflowvecd");
  // TChain* tree_mult = MRootIO::OpenChain(path_input_mult, "MultCalib");
  // TChain* tree_mult2 = MRootIO::OpenChain(path_input_mult, "MultCalib");
  TChain* tree_index = MRootIO::OpenChain(path_input_index, "EventMixing");
  // tree_flowVecd->AddFriend(tree_mult);
  // tree_flowVecd2->AddFriend(tree_mult2);

  tree_flowVecd->SetBranchStatus("*", 0);
  tree_flowVecd2->SetBranchStatus("*", 0);

  // tree_flowVecd->SetBranchStatus("fMultTPC", 1);
  // tree_flowVecd->SetBranchStatus("NumContribCalib", 1);
  // tree_flowVecd->SetBranchStatus("fMultTracklets", 1);
  // tree_flowVecd->SetBranchStatus("fMultNTracksPV", 1);
  // tree_flowVecd->SetBranchStatus("fMultFT0C", 1);
  // tree_flowVecd->SetBranchStatus("fPosX", 1);
  // tree_flowVecd->SetBranchStatus("fPosY", 1);
  // tree_flowVecd->SetBranchStatus("fPosZ", 1);
  // tree_flowVecd->SetBranchStatus("fSelection", 1);
  // tree_flowVecd->SetBranchStatus("fHadronicRate", 1);
  tree_flowVecd->SetBranchStatus("fPTREF", 1);
  tree_flowVecd->SetBranchStatus("fEtaREF", 1);
  tree_flowVecd->SetBranchStatus("fPhiREF", 1);
  tree_flowVecd->SetBranchStatus("fITSChi2NCl", 1);
  // tree_flowVecd->SetBranchStatus("fTPCNClsCR", 1);
  tree_flowVecd->SetBranchStatus("fTPCNClsFound", 1);
  tree_flowVecd->SetBranchStatus("fTPCChi2NCl", 1);
  tree_flowVecd->SetBranchStatus("fITSClusterMap", 1);
  // tree_flowVecd->SetBranchStatus("fTPCSignal", 1);
  // tree_flowVecd->SetBranchStatus("fTPCNSigmaEl", 1);
  // tree_flowVecd->SetBranchStatus("fTPCNSigmaPi", 1);
  // tree_flowVecd->SetBranchStatus("fTPCNSigmaPr", 1);
  tree_flowVecd->SetBranchStatus("fDcaXY", 1);
  tree_flowVecd->SetBranchStatus("fDcaZ", 1);

  tree_flowVecd2->SetBranchStatus("fPTREF", 1);
  tree_flowVecd2->SetBranchStatus("fEtaREF", 1);
  tree_flowVecd2->SetBranchStatus("fPhiREF", 1);
  tree_flowVecd2->SetBranchStatus("fITSChi2NCl", 1);
  // tree_flowVecd2->SetBranchStatus("fTPCNClsCR", 1);
  tree_flowVecd2->SetBranchStatus("fTPCNClsFound", 1);
  tree_flowVecd2->SetBranchStatus("fTPCChi2NCl", 1);
  tree_flowVecd2->SetBranchStatus("fITSClusterMap", 1);
  // tree_flowVecd2->SetBranchStatus("fTPCSignal", 1);
  // tree_flowVecd2->SetBranchStatus("fTPCNSigmaEl", 1);
  // tree_flowVecd2->SetBranchStatus("fTPCNSigmaPi", 1);
  // tree_flowVecd2->SetBranchStatus("fTPCNSigmaPr", 1);
  tree_flowVecd2->SetBranchStatus("fDcaXY", 1);
  tree_flowVecd2->SetBranchStatus("fDcaZ", 1);

  TTreeReader rEvt(tree_flowVecd);
  TTreeReader rEvt2(tree_flowVecd2);
  TTreeReader rPairs(tree_index);
  TTreeReaderValue<std::vector<std::pair<ULong64_t, ULong64_t>>> abPair(rPairs, "MixedEvent");
  TTreeReaderValue<int> iMultPair(rPairs, "IndexMixing_NumContribCalib");
  TTreeReaderValue<int> iVtxZPair(rPairs, "IndexMixing_PosZ");
  // TTreeReaderValue<double> NumContribCalib(rEvt, "NumContribCalib");
  // TTreeReaderValue<int> fMultTPC(rEvt, "fMultTPC");
  // TTreeReaderValue<int> fMultTracklets(rEvt, "fMultTracklets");
  // TTreeReaderValue<int> fMultNTracksPV(rEvt, "fMultNTracksPV");
  // TTreeReaderValue<float> fMultFT0C(rEvt, "fMultFT0C");
  // TTreeReaderValue<float> fPosX(rEvt, "fPosX");
  // TTreeReaderValue<float> fPosY(rEvt, "fPosY");
  // TTreeReaderValue<float> fPosZ(rEvt, "fPosZ");
  // TTreeReaderValue<ULong64_t> fSelection(rEvt, "fSelection");
  // TTreeReaderValue<float> fHadronicRate(rEvt, "fHadronicRate");
  TTreeReaderArray<float> fPTREF(rEvt, "fPTREF");
  TTreeReaderArray<float> fEtaREF(rEvt, "fEtaREF");
  TTreeReaderArray<float> fPhiREF(rEvt, "fPhiREF");
  TTreeReaderArray<float> fITSChi2NCl_ref(rEvt, "fITSChi2NCl");
  // TTreeReaderArray<float> fTPCNClsCR_ref(rEvt, "fTPCNClsCR");
  TTreeReaderArray<float> fTPCNClsFound_ref(rEvt, "fTPCNClsFound");
  TTreeReaderArray<float> fTPCChi2NCl_ref(rEvt, "fTPCChi2NCl");
  TTreeReaderArray<uint8_t> fITSClusterMap_ref(rEvt, "fITSClusterMap");
  // TTreeReaderArray<float> fTPCSignal_ref(rEvt, "fTPCSignal");
  // TTreeReaderArray<float> fTPCNSigmaEl_ref(rEvt, "fTPCNSigmaEl");
  // TTreeReaderArray<float> fTPCNSigmaPi_ref(rEvt, "fTPCNSigmaPi");
  // TTreeReaderArray<float> fTPCNSigmaPr_ref(rEvt, "fTPCNSigmaPr");
  TTreeReaderArray<float> fDcaXY_ref(rEvt, "fDcaXY");
  TTreeReaderArray<float> fDcaZ_ref(rEvt, "fDcaZ");

  TTreeReaderArray<float> fPTREF_2(rEvt2, "fPTREF");
  TTreeReaderArray<float> fEtaREF_2(rEvt2, "fEtaREF");
  TTreeReaderArray<float> fPhiREF_2(rEvt2, "fPhiREF");
  TTreeReaderArray<float> fITSChi2NCl_ref_2(rEvt2, "fITSChi2NCl");
  // TTreeReaderArray<float> fTPCNClsCR_ref_2(rEvt2, "fTPCNClsCR");
  TTreeReaderArray<float> fTPCNClsFound_ref_2(rEvt2, "fTPCNClsFound");
  TTreeReaderArray<float> fTPCChi2NCl_ref_2(rEvt2, "fTPCChi2NCl");
  TTreeReaderArray<uint8_t> fITSClusterMap_ref_2(rEvt2, "fITSClusterMap");
  // TTreeReaderArray<float> fTPCSignal_ref_2(rEvt2, "fTPCSignal");
  // TTreeReaderArray<float> fTPCNSigmaEl_ref_2(rEvt2, "fTPCNSigmaEl");
  // TTreeReaderArray<float> fTPCNSigmaPi_ref_2(rEvt2, "fTPCNSigmaPi");
  // TTreeReaderArray<float> fTPCNSigmaPr_ref_2(rEvt2, "fTPCNSigmaPr");
  TTreeReaderArray<float> fDcaXY_ref_2(rEvt2, "fDcaXY");
  TTreeReaderArray<float> fDcaZ_ref_2(rEvt2, "fDcaZ");

  TFile fout(path_output_tree, "RECREATE");
  fout.SetCompressionAlgorithm(ROOT::RCompressionSetting::EAlgorithm::kLZ4); // use LZ4 compression
  fout.SetCompressionLevel(1); // set compression level to 1 (fastest)
  TTree out("ref_pairs", "mixed ref(A) x ref(B) pairs");
  out.SetAutoSave(0);      // disable autosave
  out.SetAutoFlush(50000); // flush every 50000 bytes

  double o_NumContribCalib;
  int o_fMultTPC, o_fMultTracklets, o_fMultNTracksPV;
  float o_fMultFT0C;
  float o_fPosX, o_fPosY, o_fPosZ;
  int o_fSelection;
  float o_fHadronicRate;
  unsigned char o_iMult, o_iVtxZ;
  uint8_t o_ref1_ITSClusterMap, o_ref2_ITSClusterMap;
  out.Branch("iMult", &o_iMult);
  out.Branch("iVtxZ", &o_iVtxZ);
  // out.Branch("fMultTPC", &o_fMultTPC);
  // out.Branch("fMultTracklets", &o_fMultTracklets);
  // out.Branch("fMultNTracksPV", &o_fMultNTracksPV);
  // out.Branch("fMultFT0C", &o_fMultFT0C);
  // out.Branch("fPosX", &o_fPosX);
  // out.Branch("fPosY", &o_fPosY);
  // out.Branch("fPosZ", &o_fPosZ);
  // out.Branch("fSelection", &o_fSelection);
  // out.Branch("fHadronicRate", &o_fHadronicRate);
  float o_ref1_pt, o_ref1_eta, o_ref1_phi;
  float o_ref1_ITSChi2NCl, o_ref1_TPCNClsCR, o_ref1_TPCNClsFound;
  float o_ref1_TPCChi2NCl, o_ref1_TPCSignal;
  float o_ref1_nsig_el, o_ref1_nsig_pi, o_ref1_nsig_pr;
  float o_ref1_dcaxy, o_ref1_dcaz;
  // out.Branch("NumContribCalib", &o_NumContribCalib);
  out.Branch("ref1_pt", &o_ref1_pt);
  out.Branch("ref1_eta", &o_ref1_eta);
  out.Branch("ref1_phi", &o_ref1_phi);
  out.Branch("ref1_ITSChi2NCl", &o_ref1_ITSChi2NCl);
  // out.Branch("ref1_TPCNClsCR", &o_ref1_TPCNClsCR);
  out.Branch("ref1_TPCNClsFound", &o_ref1_TPCNClsFound);
  out.Branch("ref1_TPCChi2NCl", &o_ref1_TPCChi2NCl);
  out.Branch("ref1_ITSClusterMap", &o_ref1_ITSClusterMap);
  // out.Branch("ref1_TPCSignal", &o_ref1_TPCSignal);
  // out.Branch("ref1_nsig_el", &o_ref1_nsig_el);
  // out.Branch("ref1_nsig_pi", &o_ref1_nsig_pi);
  // out.Branch("ref1_nsig_pr", &o_ref1_nsig_pr);
  out.Branch("ref1_dcaxy", &o_ref1_dcaxy);
  out.Branch("ref1_dcaz", &o_ref1_dcaz);
  float o_ref2_pt, o_ref2_eta, o_ref2_phi;
  float o_ref2_ITSChi2NCl, o_ref2_TPCNClsCR, o_ref2_TPCNClsFound;
  float o_ref2_TPCChi2NCl, o_ref2_TPCSignal;
  float o_ref2_nsig_el, o_ref2_nsig_pi, o_ref2_nsig_pr;
  float o_ref2_dcaxy, o_ref2_dcaz;
  out.Branch("ref2_pt", &o_ref2_pt);
  out.Branch("ref2_eta", &o_ref2_eta);
  out.Branch("ref2_phi", &o_ref2_phi);
  out.Branch("ref2_ITSChi2NCl", &o_ref2_ITSChi2NCl);
  // out.Branch("ref2_TPCNClsCR", &o_ref2_TPCNClsCR);
  out.Branch("ref2_TPCNClsFound", &o_ref2_TPCNClsFound);
  out.Branch("ref2_TPCChi2NCl", &o_ref2_TPCChi2NCl);
  out.Branch("ref2_ITSClusterMap", &o_ref2_ITSClusterMap);
  // out.Branch("ref2_TPCSignal", &o_ref2_TPCSignal);
  // out.Branch("ref2_nsig_el", &o_ref2_nsig_el);
  // out.Branch("ref2_nsig_pi", &o_ref2_nsig_pi);
  // out.Branch("ref2_nsig_pr", &o_ref2_nsig_pr);
  out.Branch("ref2_dcaxy", &o_ref2_dcaxy);
  out.Branch("ref2_dcaz", &o_ref2_dcaz);
  out.SetBasketSize("*", 256 * 1024 * 1024); // set basket size to 256 M
  long long nWritten = 0;

  bool isInteractive = is_interactive();
  long long nEntries = rPairs.GetEntries();

  bool isntFirst = false;
  ULong64_t lastEventA = 0;

  int tagProcess = 0;
  while (rPairs.Next()) {
    tagProcess++;
    if (tagProcess != 9)
      continue;
    tagProcess = 0;
    o_iMult = *iMultPair;
    o_iVtxZ = *iVtxZPair;
    for (const auto& abPair_single : *abPair) {
      if (isInteractive) {
        if (isntFirst) {
          if (lastEventA != abPair_single.first) {
            rEvt.SetEntry(abPair_single.first);
            lastEventA = abPair_single.first;
          }
        } else {
          isntFirst = true;
          rEvt.SetEntry(abPair_single.first);
          lastEventA = abPair_single.first;
        }
      } else {
        rEvt.SetEntry(abPair_single.first);
      }
      rEvt2.SetEntry(abPair_single.second);
      // o_NumContribCalib = *NumContribCalib;
      // o_fMultTPC = *fMultTPC;
      // o_fMultTracklets = *fMultTracklets;
      // o_fMultNTracksPV = *fMultNTracksPV;
      // o_fMultFT0C = *fMultFT0C;
      // o_fPosX = *fPosX;
      // o_fPosY = *fPosY;
      // o_fPosZ = *fPosZ;
      // o_fSelection = *fSelection;
      // o_fHadronicRate = *fHadronicRate;

      for (int iRef1 = 0; iRef1 < fPTREF.GetSize(); iRef1++) {
        o_ref1_pt = fPTREF[iRef1];
        o_ref1_eta = fEtaREF[iRef1];
        o_ref1_phi = fPhiREF[iRef1];
        o_ref1_ITSChi2NCl = fITSChi2NCl_ref[iRef1];
        // o_ref1_TPCNClsCR = fTPCNClsCR_ref[iRef1];
        o_ref1_TPCNClsFound = fTPCNClsFound_ref[iRef1];
        o_ref1_TPCChi2NCl = fTPCChi2NCl_ref[iRef1];
        o_ref1_ITSClusterMap = fITSClusterMap_ref[iRef1];
        // o_ref1_TPCSignal = fTPCSignal_ref[iRef1];
        // o_ref1_nsig_el = fTPCNSigmaEl_ref[iRef1];
        // o_ref1_nsig_pi = fTPCNSigmaPi_ref[iRef1];
        // o_ref1_nsig_pr = fTPCNSigmaPr_ref[iRef1];
        o_ref1_dcaxy = fDcaXY_ref[iRef1];
        o_ref1_dcaz = fDcaZ_ref[iRef1];

        for (int iRef2 = 0; iRef2 < fPTREF_2.GetSize(); iRef2++) {
          o_ref2_pt = fPTREF_2[iRef2];
          o_ref2_eta = fEtaREF_2[iRef2];
          o_ref2_phi = fPhiREF_2[iRef2];
          o_ref2_ITSChi2NCl = fITSChi2NCl_ref_2[iRef2];
          // o_ref2_TPCNClsCR = fTPCNClsCR_ref_2[iRef2];
          o_ref2_TPCNClsFound = fTPCNClsFound_ref_2[iRef2];
          o_ref2_TPCChi2NCl = fTPCChi2NCl_ref_2[iRef2];
          o_ref2_ITSClusterMap = fITSClusterMap_ref_2[iRef2];
          // o_ref2_TPCSignal = fTPCSignal_ref_2[iRef2];
          // o_ref2_nsig_el = fTPCNSigmaEl_ref_2[iRef2];
          // o_ref2_nsig_pi = fTPCNSigmaPi_ref_2[iRef2];
          // o_ref2_nsig_pr = fTPCNSigmaPr_ref_2[iRef2];
          o_ref2_dcaxy = fDcaXY_ref_2[iRef2];
          o_ref2_dcaz = fDcaZ_ref_2[iRef2];

          out.Fill();
        }
      }
    }
  }
  out.Write();
  fout.Close();
}

int main(int argc, char** argv) {
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
  EventMixingRef(path_input_flowVecd, path_input_mult, path_output, path_output_tree);

  return 0;
}
