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
#include <TStopwatch.h>
#include <fstream>
#include <iomanip>
#include <iostream>

template <typename T> std::vector<T> makeVec(const TTreeReaderArray<T>& arr) {
  return std::vector<T>(arr.begin(), arr.end());
}

template <typename T> const T* makePtr(const TTreeReaderArray<T>& arr, size_t* size = nullptr) {
  if (size) {
    *size = arr.GetSize();
  }
  return arr.GetSize() > 0 ? &arr[0] : nullptr;
}

vector<pair<ULong64_t, ULong64_t>> MixEvent(unsigned int, const int id, const ULong64_t& event_id) {
  return MixVec<pair<ULong64_t, ULong64_t>, ULong64_t>(
      id, event_id, [](const ULong64_t& a, const ULong64_t& b) { return make_pair(a, b); }, 100);
}

void EventMixingIndexGen(TString path_input_flowVecd = "../input.root",
                         TString path_input_mult = "../input2.root",
                         TString path_output = "output.root",
                         TString path_output_tree = "output_tree.root") {
  // close multi-thread
  // ROOT::EnableImplicitMT(4);
  ROOT::DisableImplicitMT();

  StrVar4Hist var_fPosX("fPosX", "#it{V}_{x}", "cm", 200, {-10, 10});
  StrVar4Hist var_fPosY("fPosY", "#it{V}_{Y}", "cm", 200, {-10, 10});
  StrVar4Hist var_fPosZ("fPosZ", "#it{V}_{Z}", "cm", 200, {-10, 10});
  StrVar4Hist var_fPosZMix("fPosZ", "#it{V}_{Z}", "cm", 10, {-10, 10});
  StrVar4Hist var_fNumContrib("fNumContrib", "#it{N}_{vtx contrib} ", "", 300, {0, 300});
  StrVar4Hist var_NumContribCalib("NumContribCalib", "N_{vtx contrib} Calibrated", "", 300,
                                  {0, 300});
  StrVar4Hist var_NumContribCalibBinned("NumContribCalib", "N_{vtx contrib} Calibrated", "", 10,
                                        {0, 5, 8, 11, 14, 18, 23, 28, 36, 48, 300});
  StrVar4Hist var_fMultTPC("fMultTPC", "Mult_{TPC}", "", 600, {0, 600});
  StrVar4Hist var_fMultREF("fMultREF", "Mult_{REF}", "", 100, {0, 100});
  StrVar4Hist var_fMultFT0C("fMultFT0C", "Mult_{FT0C}", "", 130, {-1000., 12000.});
  StrVar4Hist var_MultNTracksPV("fMultNTracksPV", "#it{N}_{Tracks PV}", "", 150, {0, 150});
  StrVar4Hist var_MassJpsiCandidate("fMass", "M_{ee}", "GeV^{2}/c^{4}", 100, {1.8, 5.4});
  StrVar4Hist var_PtJpsiCandidate("fPT", "p_{T}", "GeV/c", 10, {0., 5.});

  const vector<double> bins_mix_numContrib = var_NumContribCalibBinned.fBins;
  const vector<double> bins_mix_posZ = var_fPosZMix.fBins;

  TChain* tree_flowVecd = MRootIO::OpenChain(path_input_flowVecd, "O2dqflowvecd");
  TChain* tree_mult = MRootIO::OpenChain(path_input_mult, "MultCalib");

  tree_flowVecd->AddFriend(tree_mult);

  cout << "Input file: " << path_input_flowVecd << endl;
  cout << "Input file: " << path_input_mult << endl;
  cout << "Output file: " << path_output << endl;
  cout << "Output tree file: " << path_output_tree << endl;

  ROOT::RDataFrame rdf(*tree_flowVecd);

  auto rdf_witTrigger =
      rdf.Define("map_trigger", MALICE::triggermapRVec, {"fSelection"})
          .Define("isntSPDPileup", MALICE::IsntSPDPileup, {"fSelection"})
          .Define("isntTPCPileup", MALICE::IsntTPCPileup, {"fSelection"})
          .Define("isntSameBunchPileup", MALICE::IsntSameBunchPileup_NoSlot, {"fSelection"})
          .Define("isntITSROFrameBorder", MALICE::IsntITSROFrameBorder, {"fSelection"})
          .Define("isntTimeFrameBorder", MALICE::IsntTimeFrameBorder, {"fSelection"})
          .Define("isTriggerTVX", MALICE::IsTriggerTVX, {"fSelection"})
          .Alias("fMultREF", "fPTREF_size")
          .Define("RunNumber", [] { return float(0.5); })
      /*   .DefineSlot("isntSelfDefinedPileup",
                    Cut_MultTPC_NumContrib::isInCutSlot,
                    {"NumContribCalib", "fMultTPC"}) */
      ;
  auto rdf_isntITSROFrameBorder =
      rdf_witTrigger.Filter("isntITSROFrameBorder", "no ITS RO Frame border");
  auto rdf_isntTimeFrameBorder =
      rdf_witTrigger.Filter("isntTimeFrameBorder", "no Time Frame border");
  auto rdf_isTriggerTVX = rdf_witTrigger.Filter("isTriggerTVX", "is Trigger TVX");
  auto rdf_PartTrigger = rdf_witTrigger.Filter("isTriggerTVX", "is Trigger TVX")
                             .Filter("isntITSROFrameBorder", "no ITS RO Frame border")
                             .Filter("isntTimeFrameBorder", "no Time Frame border")
                             .Filter("isntSameBunchPileup", "no Time Frame border")
      /*  .Filter("isntSelfDefinedPileup", "no self defined pileup") */;

  auto rdf_PartTriggerWithJpsi = rdf_PartTrigger.Filter("fEta_size>=1", "has Jpsi");

  auto rdf_PartTriggerWithJpsiWithEvent =
      rdf_PartTriggerWithJpsi
          .Define("EventData", CreateEventData,
                  {"fMultTPC", "fMultTracklets", "fMultNTracksPV", "fMultFT0C", "fNumContrib",
                   "NumContribCalib", "fPosX", "fPosY", "fPosZ", "fSelection", "fHadronicRate",
                   "fPT", "fEta", "fPhi", "fMass", "fSign", "fPTREF", "fEtaREF", "fPhiREF"})
          .Define("isEventGood", [](const EventData& event) { return event.isGood(); },
                  {"EventData"})
          .Filter("isEventGood", "Event is good");

  auto rdf_PartTriggerWithJpsiWithEventWithEventMixing =
      rdf_PartTriggerWithJpsiWithEvent
          .Define("IndexMixing_NumContribCalib",
                  [bins_mix_numContrib](double numContrib) {
                    int index = -1;
                    for (int i = 0; i < bins_mix_numContrib.size() - 1; i++) {
                      if (numContrib >= bins_mix_numContrib[i] &&
                          numContrib < bins_mix_numContrib[i + 1]) {
                        index = i;
                        break;
                      }
                    }
                    if (index == -1) {
                      return -1; // Invalid index
                    }
                    return int(index);
                  },
                  {"NumContribCalib"})
          .Filter("IndexMixing_NumContribCalib>=0", "valid NumContribCalib index")
          .Define("IndexMixing_PosZ",
                  [bins_mix_posZ](float posZ) {
                    int index = -1;
                    for (int i = 0; i < bins_mix_posZ.size() - 1; i++) {
                      if (posZ >= bins_mix_posZ[i] && posZ < bins_mix_posZ[i + 1]) {
                        index = i;
                        break;
                      }
                    }
                    if (index == -1) {
                      return -1; // Invalid index
                    }
                    return int(index);
                  },
                  {"fPosZ"})
          .Filter("IndexMixing_PosZ>=0", "valid PosZ index")
          .Define("IndexMixing",
                  [bins_mix_posZ](int indexNumContrib, int indexPosZ) {
                    if (indexNumContrib < 0 || indexPosZ < 0) {
                      return -1; // Invalid index
                    }
                    return int(indexPosZ + indexNumContrib * bins_mix_posZ.size());
                  },
                  {"IndexMixing_NumContribCalib", "IndexMixing_PosZ"})
          .DefineSlot("MixedEvent", MixEvent, {"IndexMixing", "rdfentry_"});

  rdf_PartTriggerWithJpsiWithEventWithEventMixing.Snapshot(
      "EventMixing", path_output_tree,
      {"MixedEvent", "IndexMixing_NumContribCalib", "IndexMixing_PosZ"});
  if (is_interactive())
    ROOT::RDF::Experimental::AddProgressBar(rdf_PartTriggerWithJpsiWithEventWithEventMixing);
}

void EventMixingJpsiAssoPair(TString path_input_flowVecd = "../input1.root",
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
  tree_flowVecd->SetBranchStatus("fPT", 1);
  tree_flowVecd->SetBranchStatus("fEta", 1);
  tree_flowVecd->SetBranchStatus("fPhi", 1);
  tree_flowVecd->SetBranchStatus("fMass", 1);
  // tree_flowVecd->SetBranchStatus("fSign", 1);
  tree_flowVecd->SetBranchStatus("fPt1", 1);
  tree_flowVecd->SetBranchStatus("fEta1", 1);
  tree_flowVecd->SetBranchStatus("fPhi1", 1);
  // tree_flowVecd->SetBranchStatus("fSign1", 1);
  // tree_flowVecd->SetBranchStatus("fITSChi2NCl1", 1);
  // tree_flowVecd->SetBranchStatus("fTPCNClsCR1", 1);
  // tree_flowVecd->SetBranchStatus("fTPCNClsFound1", 1);
  // tree_flowVecd->SetBranchStatus("fTPCChi2NCl1", 1);
  // tree_flowVecd->SetBranchStatus("fTPCSignal1", 1);
  tree_flowVecd->SetBranchStatus("fTPCNSigmaEl1", 1);
  tree_flowVecd->SetBranchStatus("fTPCNSigmaPi1", 1);
  tree_flowVecd->SetBranchStatus("fTPCNSigmaPr1", 1);
  tree_flowVecd->SetBranchStatus("fPt2", 1);
  tree_flowVecd->SetBranchStatus("fEta2", 1);
  tree_flowVecd->SetBranchStatus("fPhi2", 1);
  //  tree_flowVecd->SetBranchStatus("fSign2", 1);
  //  tree_flowVecd->SetBranchStatus("fITSChi2NCl2", 1);
  //  tree_flowVecd->SetBranchStatus("fTPCNClsCR2", 1);
  //  tree_flowVecd->SetBranchStatus("fTPCNClsFound2", 1);
  //  tree_flowVecd->SetBranchStatus("fTPCChi2NCl2", 1);
  //  tree_flowVecd->SetBranchStatus("fTPCSignal2", 1);
  tree_flowVecd->SetBranchStatus("fTPCNSigmaEl2", 1);
  tree_flowVecd->SetBranchStatus("fTPCNSigmaPi2", 1);
  tree_flowVecd->SetBranchStatus("fTPCNSigmaPr2", 1);
  tree_flowVecd2->SetBranchStatus("fPTREF", 1);
  tree_flowVecd2->SetBranchStatus("fEtaREF", 1);
  tree_flowVecd2->SetBranchStatus("fPhiREF", 1);
  tree_flowVecd2->SetBranchStatus("fITSChi2NCl", 1);
  // tree_flowVecd2->SetBranchStatus("fTPCNClsCR", 1);
  tree_flowVecd2->SetBranchStatus("fTPCNClsFound", 1);
  tree_flowVecd2->SetBranchStatus("fTPCChi2NCl", 1);
  tree_flowVecd2->SetBranchStatus("fITSClusterMap", 1);
  // tree_flowVecd2->SetBranchStatus("fTPCSignal", 1);
  tree_flowVecd2->SetBranchStatus("fTPCNSigmaEl", 1);
  tree_flowVecd2->SetBranchStatus("fTPCNSigmaPi", 1);
  tree_flowVecd2->SetBranchStatus("fTPCNSigmaPr", 1);
  tree_flowVecd2->SetBranchStatus("fDcaXY", 1);
  tree_flowVecd2->SetBranchStatus("fDcaZ", 1);
  TTreeReader rPairs(tree_index);
  TTreeReaderValue<std::vector<std::pair<ULong64_t, ULong64_t>>> abPair(rPairs, "MixedEvent");
  TTreeReaderValue<int> iMultPair(rPairs, "IndexMixing_NumContribCalib");
  TTreeReaderValue<int> iVtxZPair(rPairs, "IndexMixing_PosZ");

  TTreeReader rEvt(tree_flowVecd);
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
  TTreeReaderArray<float> fPT(rEvt, "fPT");
  TTreeReaderArray<float> fEta(rEvt, "fEta");
  TTreeReaderArray<float> fPhi(rEvt, "fPhi");
  TTreeReaderArray<float> fMass(rEvt, "fMass");
  // TTreeReaderArray<float> fSign(rEvt, "fSign");
  TTreeReaderArray<float> fPt1(rEvt, "fPt1");
  TTreeReaderArray<float> fEta1(rEvt, "fEta1");
  TTreeReaderArray<float> fPhi1(rEvt, "fPhi1");
  // TTreeReaderArray<int> fSign1(rEvt, "fSign1");
  // TTreeReaderArray<float> fITSChi2NCl1(rEvt, "fITSChi2NCl1");
  // TTreeReaderArray<float> fTPCNClsCR1(rEvt, "fTPCNClsCR1");
  // TTreeReaderArray<float> fTPCNClsFound1(rEvt, "fTPCNClsFound1");
  // TTreeReaderArray<float> fTPCChi2NCl1(rEvt, "fTPCChi2NCl1");
  // TTreeReaderArray<float> fTPCSignal1(rEvt, "fTPCSignal1");
  TTreeReaderArray<float> fTPCNSigmaEl1(rEvt, "fTPCNSigmaEl1");
  TTreeReaderArray<float> fTPCNSigmaPi1(rEvt, "fTPCNSigmaPi1");
  TTreeReaderArray<float> fTPCNSigmaPr1(rEvt, "fTPCNSigmaPr1");
  TTreeReaderArray<float> fPt2(rEvt, "fPt2");
  TTreeReaderArray<float> fEta2(rEvt, "fEta2");
  TTreeReaderArray<float> fPhi2(rEvt, "fPhi2");
  // TTreeReaderArray<int> fSign2(rEvt, "fSign2");
  // TTreeReaderArray<float> fITSChi2NCl2(rEvt, "fITSChi2NCl2");
  // TTreeReaderArray<float> fTPCNClsCR2(rEvt, "fTPCNClsCR2");
  // TTreeReaderArray<float> fTPCNClsFound2(rEvt, "fTPCNClsFound2");
  // TTreeReaderArray<float> fTPCChi2NCl2(rEvt, "fTPCChi2NCl2");
  // TTreeReaderArray<float> fTPCSignal2(rEvt, "fTPCSignal2");
  TTreeReaderArray<float> fTPCNSigmaEl2(rEvt, "fTPCNSigmaEl2");
  TTreeReaderArray<float> fTPCNSigmaPi2(rEvt, "fTPCNSigmaPi2");
  TTreeReaderArray<float> fTPCNSigmaPr2(rEvt, "fTPCNSigmaPr2");

  TTreeReader rEvt2(tree_flowVecd2);
  TTreeReaderArray<float> fPTREF(rEvt2, "fPTREF");
  TTreeReaderArray<float> fEtaREF(rEvt2, "fEtaREF");
  TTreeReaderArray<float> fPhiREF(rEvt2, "fPhiREF");
  TTreeReaderArray<float> fITSChi2NCl_ref(rEvt2, "fITSChi2NCl");
  // TTreeReaderArray<float> fTPCNClsCR_ref(rEvt2, "fTPCNClsCR");
  TTreeReaderArray<float> fTPCNClsFound_ref(rEvt2, "fTPCNClsFound");
  TTreeReaderArray<float> fTPCChi2NCl_ref(rEvt2, "fTPCChi2NCl");
  TTreeReaderArray<uint8_t> fITSClusterMap_ref(rEvt2, "fITSClusterMap");
  // TTreeReaderArray<float> fTPCSignal_ref(rEvt2, "fTPCSignal");
  // TTreeReaderArray<float> fTPCNSigmaEl_ref(rEvt2, "fTPCNSigmaEl");
  // TTreeReaderArray<float> fTPCNSigmaPi_ref(rEvt2, "fTPCNSigmaPi");
  // TTreeReaderArray<float> fTPCNSigmaPr_ref(rEvt2, "fTPCNSigmaPr");
  TTreeReaderArray<float> fDcaXY_ref(rEvt2, "fDcaXY");
  TTreeReaderArray<float> fDcaZ_ref(rEvt2, "fDcaZ");

  TFile fout(path_output_tree, "RECREATE");
  fout.SetCompressionAlgorithm(ROOT::RCompressionSetting::EAlgorithm::kLZ4); // use LZ4 compression
  fout.SetCompressionLevel(1); // set compression level to 1 (fastest)
  TTree out("jpsi_ref_pairs", "mixed jpsi(A) x ref(B) pairs");
  out.SetAutoSave(0);      // disable autosave
  out.SetAutoFlush(50000); // flush every 50000 bytes

  double o_NumContribCalib;
  int o_fMultTPC, o_fMultTracklets, o_fMultNTracksPV;
  float o_fMultFT0C;
  float o_fPosX, o_fPosY, o_fPosZ;
  int o_fSelection;
  float o_fHadronicRate;
  float o_jpsi_pt, o_jpsi_eta, o_jpsi_phi, o_jpsi_mass, o_jpsi_sign;
  float o_e1_pt, o_e1_eta, o_e1_phi;
  int o_e1_sign;
  float o_e1_ITSChi2NCl, o_e1_TPCNClsCR, o_e1_TPCNClsFound;
  float o_e1_TPCChi2NCl, o_e1_TPCSignal;
  float o_e1_nsig_el, o_e1_nsig_pi, o_e1_nsig_pr;
  float o_e2_pt, o_e2_eta, o_e2_phi;
  int o_e2_sign;
  uint8_t o_ref_ITSClusterMap;
  float o_e2_ITSChi2NCl, o_e2_TPCNClsCR, o_e2_TPCNClsFound;
  float o_e2_TPCChi2NCl, o_e2_TPCSignal;
  float o_e2_nsig_el, o_e2_nsig_pi, o_e2_nsig_pr;
  float o_ref_pt, o_ref_eta, o_ref_phi;
  float o_ref_ITSChi2NCl, o_ref_TPCNClsCR, o_ref_TPCNClsFound;
  float o_ref_TPCChi2NCl, o_ref_TPCSignal;
  float o_ref_nsig_el, o_ref_nsig_pi, o_ref_nsig_pr;
  float o_ref_dcaxy, o_ref_dcaz;
  unsigned char o_iMult, o_iVtxZ;
  // out.Branch("NumContribCalib", &o_NumContribCalib);
  // out.Branch("fMultTPC", &o_fMultTPC);
  // out.Branch("fMultTracklets", &o_fMultTracklets);
  // out.Branch("fMultNTracksPV", &o_fMultNTracksPV);
  // out.Branch("fMultFT0C", &o_fMultFT0C);
  // out.Branch("fPosX", &o_fPosX);
  // out.Branch("fPosY", &o_fPosY);
  // out.Branch("fPosZ", &o_fPosZ);
  // out.Branch("fSelection", &o_fSelection);
  // out.Branch("fHadronicRate", &o_fHadronicRate);

  out.Branch("iMult", &o_iMult);
  out.Branch("iVtxZ", &o_iVtxZ);
  out.Branch("jpsi_pt", &o_jpsi_pt);
  out.Branch("jpsi_eta", &o_jpsi_eta);
  out.Branch("jpsi_phi", &o_jpsi_phi);
  out.Branch("jpsi_mass", &o_jpsi_mass);
  // out.Branch("jpsi_sign", &o_jpsi_sign);
  out.Branch("e1_pt", &o_e1_pt);
  out.Branch("e1_eta", &o_e1_eta);
  out.Branch("e1_phi", &o_e1_phi);
  // out.Branch("e1_sign", &o_e1_sign);
  // out.Branch("e1_ITSChi2NCl", &o_e1_ITSChi2NCl);
  // out.Branch("e1_TPCNClsCR", &o_e1_TPCNClsCR);
  // out.Branch("e1_TPCNClsFound", &o_e1_TPCNClsFound);
  // out.Branch("e1_TPCChi2NCl", &o_e1_TPCChi2NCl);
  // out.Branch("e1_TPCSignal", &o_e1_TPCSignal);
  out.Branch("e1_nsig_el", &o_e1_nsig_el);
  out.Branch("e1_nsig_pi", &o_e1_nsig_pi);
  out.Branch("e1_nsig_pr", &o_e1_nsig_pr);

  out.Branch("e2_pt", &o_e2_pt);
  out.Branch("e2_eta", &o_e2_eta);
  out.Branch("e2_phi", &o_e2_phi);
  // out.Branch("e2_sign", &o_e2_sign);
  // out.Branch("e2_ITSChi2NCl", &o_e2_ITSChi2NCl);
  // out.Branch("e2_TPCNClsCR", &o_e2_TPCNClsCR);
  // out.Branch("e2_TPCNClsFound", &o_e2_TPCNClsFound);
  // out.Branch("e2_TPCChi2NCl", &o_e2_TPCChi2NCl);
  // out.Branch("e2_TPCSignal", &o_e2_TPCSignal);
  out.Branch("e2_nsig_el", &o_e2_nsig_el);
  out.Branch("e2_nsig_pi", &o_e2_nsig_pi);
  out.Branch("e2_nsig_pr", &o_e2_nsig_pr);

  out.Branch("ref_pt", &o_ref_pt);
  out.Branch("ref_eta", &o_ref_eta);
  out.Branch("ref_phi", &o_ref_phi);

  out.Branch("ref_ITSChi2NCl", &o_ref_ITSChi2NCl);
  out.Branch("ref_TPCNClsFound", &o_ref_TPCNClsFound);
  out.Branch("ref_TPCChi2NCl", &o_ref_TPCChi2NCl);
  out.Branch("ref_ITSClusterMap", &o_ref_ITSClusterMap);

  // out.Branch("ref_TPCSignal", &o_ref_TPCSignal);
  // out.Branch("ref_nsig_el", &o_ref_nsig_el);
  // out.Branch("ref_nsig_pi", &o_ref_nsig_pi);
  // out.Branch("ref_nsig_pr", &o_ref_nsig_pr);
  out.Branch("ref_dcaxy", &o_ref_dcaxy);
  out.Branch("ref_dcaz", &o_ref_dcaz);
  out.SetBasketSize("*", 16 * 1024 * 1024); // set basket size to 256 KB
  long long nWritten = 0;

  bool isInteractive = is_interactive();
  long long nEntries = rPairs.GetEntries();

  bool isntFirst = false;
  ULong64_t lastEventA = 0;
  while (rPairs.Next()) {
    o_iMult = iMultPair;
    o_iVtxZ = iVtxZPair;
    for (const auto& abPair_single : *abPair) {
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
      rEvt2.SetEntry(abPair_single.second);
      if (fPT.GetSize() == 0)
        continue;
      if (fPTREF.GetSize() == 0)
        continue;
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

      for (int iJpsi = 0; iJpsi < fEta.GetSize(); iJpsi++) {
        o_jpsi_pt = fPT[iJpsi];
        o_jpsi_eta = fEta[iJpsi];
        o_jpsi_phi = fPhi[iJpsi];
        o_jpsi_mass = fMass[iJpsi];
        // o_jpsi_sign = fSign[iJpsi];
        o_e1_pt = fPt1[iJpsi];
        o_e1_eta = fEta1[iJpsi];
        o_e1_phi = fPhi1[iJpsi];
        // o_e1_sign = fSign1[iJpsi];
        // o_e1_ITSChi2NCl = fITSChi2NCl1[iJpsi];
        // o_e1_TPCNClsCR = fTPCNClsCR1[iJpsi];
        // o_e1_TPCNClsFound = fTPCNClsFound1[iJpsi];
        // o_e1_TPCChi2NCl = fTPCChi2NCl1[iJpsi];
        // o_e1_TPCSignal = fTPCSignal1[iJpsi];
        o_e1_nsig_el = fTPCNSigmaEl1[iJpsi];
        o_e1_nsig_pi = fTPCNSigmaPi1[iJpsi];
        o_e1_nsig_pr = fTPCNSigmaPr1[iJpsi];
        o_e2_pt = fPt2[iJpsi];
        o_e2_eta = fEta2[iJpsi];
        o_e2_phi = fPhi2[iJpsi];
        // o_e2_sign = fSign2[iJpsi];
        // o_e2_ITSChi2NCl = fITSChi2NCl2[iJpsi];
        // o_e2_TPCNClsCR = fTPCNClsCR2[iJpsi];
        // o_e2_TPCNClsFound = fTPCNClsFound2[iJpsi];
        // o_e2_TPCChi2NCl = fTPCChi2NCl2[iJpsi];
        // o_e2_TPCSignal = fTPCSignal2[iJpsi];
        o_e2_nsig_el = fTPCNSigmaEl2[iJpsi];
        o_e2_nsig_pi = fTPCNSigmaPi2[iJpsi];
        o_e2_nsig_pr = fTPCNSigmaPr2[iJpsi];

        for (int iRef = 0; iRef < fPTREF.GetSize(); iRef++) {
          o_ref_pt = fPTREF[iRef];
          o_ref_eta = fEtaREF[iRef];
          o_ref_phi = fPhiREF[iRef];
          o_ref_ITSChi2NCl = fITSChi2NCl_ref[iRef];
          // o_ref_TPCNClsCR = fTPCNClsCR_ref[iRef];
          o_ref_TPCNClsFound = fTPCNClsFound_ref[iRef];
          o_ref_TPCChi2NCl = fTPCChi2NCl_ref[iRef];
          o_ref_ITSClusterMap = fITSClusterMap_ref[iRef];
          // o_ref_TPCSignal = fTPCSignal_ref[iRef];
          // o_ref_nsig_el = fTPCNSigmaEl_ref[iRef];
          // o_ref_nsig_pi = fTPCNSigmaPi_ref[iRef];
          // o_ref_nsig_pr = fTPCNSigmaPr_ref[iRef];
          o_ref_dcaxy = fDcaXY_ref[iRef];
          o_ref_dcaz = fDcaZ_ref[iRef];

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
  TString path_output_mix = "output_mix.root";
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
    path_output_mix = argv[4];
  }
  if (argc > 5) {
    path_output_tree = argv[5];
  }

  gROOT->SetBatch(kTRUE); // Disable interactive graphic
  cout << "Starting Event Mixing J/psi-Associated Pair Analysis..." << endl;
  EventMixingIndexGen(path_input_flowVecd, path_input_mult, path_output, path_output_mix);
  cout << "Event Mixing Index Generation Completed." << endl;
  EventMixingJpsiAssoPair(path_input_flowVecd, path_input_mult, path_output_mix, path_output_tree);
  cout << "Event Mixing J/psi-Associated Pair Analysis Completed." << endl;
  return 0;
}
