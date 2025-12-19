#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MEventMixing.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "opt/EventData.h"
#include <ROOT/RDataFrame.hxx>

vector<pair<ULong64_t, ULong64_t>> MixEvent(unsigned int, const int id,
                                            const ULong64_t &event_id) {
  return MixVec<pair<ULong64_t, ULong64_t>, ULong64_t>(
      id, event_id,
      [](const ULong64_t &a, const ULong64_t &b) { return make_pair(a, b); },
      100);
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
  StrVar4Hist var_fNumContrib("fNumContrib", "#it{N}_{vtx contrib} ", "", 300,
                              {0, 300});
  StrVar4Hist var_NumContribCalib(
      "NumContribCalib", "N_{vtx contrib} Calibrated", "", 300, {0, 300});
  StrVar4Hist var_NumContribCalibBinned(
      "NumContribCalib", "N_{vtx contrib} Calibrated", "", 10,
      {0, 7, 12, 17, 22, 29, 37, 46, 57, 73, 300});
  StrVar4Hist var_fMultTPC("fMultTPC", "Mult_{TPC}", "", 600, {0, 600});
  StrVar4Hist var_fMultREF("fMultREF", "Mult_{REF}", "", 100, {0, 100});
  StrVar4Hist var_fMultFT0C("fMultFT0C", "Mult_{FT0C}", "", 130,
                            {-1000., 12000.});
  StrVar4Hist var_MultNTracksPV("fMultNTracksPV", "#it{N}_{Tracks PV}", "", 150,
                                {0, 150});
  StrVar4Hist var_MassJpsiCandidate("fMass", "M_{ee}", "GeV^{2}/c^{4}", 100,
                                    {1.8, 5.4});
  StrVar4Hist var_PtJpsiCandidate("fPT", "p_{T}", "GeV/c", 10, {0., 5.});

  const vector<double> bins_mix_numContrib = var_NumContribCalibBinned.fBins;
  const vector<double> bins_mix_posZ = var_fPosZMix.fBins;

  TChain *tree_flowVecd =
      MRootIO::OpenChain(path_input_flowVecd, "O2dqflowvecd");
  TChain *tree_mult = MRootIO::OpenChain(path_input_mult, "MultCalib");

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
          .Define("isntSameBunchPileup", MALICE::IsntSameBunchPileup_NoSlot,
                  {"fSelection"})
          .Define("isntITSROFrameBorder", MALICE::IsntITSROFrameBorder,
                  {"fSelection"})
          .Define("isntTimeFrameBorder", MALICE::IsntTimeFrameBorder,
                  {"fSelection"})
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
  auto rdf_isTriggerTVX =
      rdf_witTrigger.Filter("isTriggerTVX", "is Trigger TVX");
  auto rdf_PartTrigger =
      rdf_witTrigger.Filter("isTriggerTVX", "is Trigger TVX")
          .Filter("isntITSROFrameBorder", "no ITS RO Frame border")
          .Filter("isntTimeFrameBorder", "no Time Frame border")
          .Filter("isntSameBunchPileup", "no Time Frame border")
      /*  .Filter("isntSelfDefinedPileup", "no self defined pileup") */;

  auto rdf_PartTriggerWithJpsi =
      rdf_PartTrigger.Filter("fEta_size>=1", "has Jpsi");

  auto rdf_PartTriggerWithJpsiWithEvent =
      rdf_PartTriggerWithJpsi
          .Define("EventData", CreateEventData,
                  {"fMultTPC", "fMultTracklets", "fMultNTracksPV", "fMultFT0C",
                   "fNumContrib", "NumContribCalib", "fPosX", "fPosY", "fPosZ",
                   "fSelection", "fHadronicRate", "fPT", "fEta", "fPhi",
                   "fMass", "fSign", "fPTREF", "fEtaREF", "fPhiREF"})
          .Define("isEventGood",
                  [](const EventData &event) { return event.isGood(); },
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
          .Define("IndexMixing_PosZ",
                  [bins_mix_posZ](float posZ) {
                    int index = -1;
                    for (int i = 0; i < bins_mix_posZ.size() - 1; i++) {
                      if (posZ >= bins_mix_posZ[i] &&
                          posZ < bins_mix_posZ[i + 1]) {
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
          .Define("IndexMixing",
                  [bins_mix_posZ](int indexNumContrib, int indexPosZ) {
                    if (indexNumContrib < 0 || indexPosZ < 0) {
                      return -1; // Invalid index
                    }
                    return int(indexPosZ +
                               indexNumContrib * bins_mix_posZ.size());
                  },
                  {"IndexMixing_NumContribCalib", "IndexMixing_PosZ"})
          .DefineSlot("MixedEvent", MixEvent, {"IndexMixing", "rdfentry_"});

  rdf_PartTriggerWithJpsiWithEventWithEventMixing.Snapshot(
      "EventMixing", path_output_tree,
      {"MixedEvent", "IndexMixing_NumContribCalib", "IndexMixing_PosZ"});
  if (is_interactive())
    ROOT::RDF::Experimental::AddProgressBar(
        rdf_PartTriggerWithJpsiWithEventWithEventMixing);
}

void EventMixingJpsiAssoPair(TString path_input_flowVecd = "../input1.root",
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
  // event-level（从事件A拷贝到输出；你 Snapshot 输出里有哪些就读哪些）
  TTreeReaderValue<int> NumContribCalib(rEvt, "NumContribCalib");
  TTreeReaderValue<float> fMultTPC(rEvt, "fMultTPC");
  TTreeReaderValue<float> fMultTracklets(rEvt, "fMultTracklets");
  TTreeReaderValue<float> fMultNTracksPV(rEvt, "fMultNTracksPV");
  TTreeReaderValue<float> fMultFT0C(rEvt, "fMultFT0C");
  TTreeReaderValue<float> fPosX(rEvt, "fPosX");
  TTreeReaderValue<float> fPosY(rEvt, "fPosY");
  TTreeReaderValue<float> fPosZ(rEvt, "fPosZ");
  TTreeReaderValue<int> fSelection(rEvt, "fSelection");
  TTreeReaderValue<float> fHadronicRate(rEvt, "fHadronicRate");

  // J/psi candidate 容器（事件A）
  TTreeReaderValue<std::vector<float>> fPT(rEvt, "fPT");
  TTreeReaderValue<std::vector<float>> fEta(rEvt, "fEta");
  TTreeReaderValue<std::vector<float>> fPhi(rEvt, "fPhi");
  TTreeReaderValue<std::vector<float>> fMass(rEvt, "fMass");
  TTreeReaderValue<std::vector<float>> fSign(rEvt, "fSign");

  // electron 1/2（事件A）
  TTreeReaderValue<std::vector<float>> fPt1(rEvt, "fPt1");
  TTreeReaderValue<std::vector<float>> fEta1(rEvt, "fEta1");
  TTreeReaderValue<std::vector<float>> fPhi1(rEvt, "fPhi1");
  TTreeReaderValue<std::vector<int>> fSign1(rEvt, "fSign1");
  TTreeReaderValue<std::vector<float>> fITSChi2NCl1(rEvt, "fITSChi2NCl1");
  TTreeReaderValue<std::vector<float>> fTPCNClsCR1(rEvt, "fTPCNClsCR1");
  TTreeReaderValue<std::vector<float>> fTPCNClsFound1(rEvt, "fTPCNClsFound1");
  TTreeReaderValue<std::vector<float>> fTPCChi2NCl1(rEvt, "fTPCChi2NCl1");
  TTreeReaderValue<std::vector<float>> fTPCSignal1(rEvt, "fTPCSignal1");
  TTreeReaderValue<std::vector<float>> fTPCNSigmaEl1(rEvt, "fTPCNSigmaEl1");
  TTreeReaderValue<std::vector<float>> fTPCNSigmaPi1(rEvt, "fTPCNSigmaPi1");
  TTreeReaderValue<std::vector<float>> fTPCNSigmaPr1(rEvt, "fTPCNSigmaPr1");

  TTreeReaderValue<std::vector<float>> fPt2(rEvt, "fPt2");
  TTreeReaderValue<std::vector<float>> fEta2(rEvt, "fEta2");
  TTreeReaderValue<std::vector<float>> fPhi2(rEvt, "fPhi2");
  TTreeReaderValue<std::vector<int>> fSign2(rEvt, "fSign2");
  TTreeReaderValue<std::vector<float>> fITSChi2NCl2(rEvt, "fITSChi2NCl2");
  TTreeReaderValue<std::vector<float>> fTPCNClsCR2(rEvt, "fTPCNClsCR2");
  TTreeReaderValue<std::vector<float>> fTPCNClsFound2(rEvt, "fTPCNClsFound2");
  TTreeReaderValue<std::vector<float>> fTPCChi2NCl2(rEvt, "fTPCChi2NCl2");
  TTreeReaderValue<std::vector<float>> fTPCSignal2(rEvt, "fTPCSignal2");
  TTreeReaderValue<std::vector<float>> fTPCNSigmaEl2(rEvt, "fTPCNSigmaEl2");
  TTreeReaderValue<std::vector<float>> fTPCNSigmaPi2(rEvt, "fTPCNSigmaPi2");
  TTreeReaderValue<std::vector<float>> fTPCNSigmaPr2(rEvt, "fTPCNSigmaPr2");

  // ref track 容器（事件B）
  TTreeReaderValue<std::vector<float>> fPTREF(rEvt, "fPTREF");
  TTreeReaderValue<std::vector<float>> fEtaREF(rEvt, "fEtaREF");
  TTreeReaderValue<std::vector<float>> fPhiREF(rEvt, "fPhiREF");
  TTreeReaderValue<std::vector<float>> fITSChi2NCl_ref(rEvt, "fITSChi2NCl");
  TTreeReaderValue<std::vector<float>> fTPCNClsCR_ref(rEvt, "fTPCNClsCR");
  TTreeReaderValue<std::vector<float>> fTPCNClsFound_ref(rEvt, "fTPCNClsFound");
  TTreeReaderValue<std::vector<float>> fTPCChi2NCl_ref(rEvt, "fTPCChi2NCl");
  TTreeReaderValue<std::vector<float>> fTPCSignal_ref(rEvt, "fTPCSignal");
  TTreeReaderValue<std::vector<float>> fTPCNSigmaEl_ref(rEvt, "fTPCNSigmaEl");
  TTreeReaderValue<std::vector<float>> fTPCNSigmaPi_ref(rEvt, "fTPCNSigmaPi");
  TTreeReaderValue<std::vector<float>> fTPCNSigmaPr_ref(rEvt, "fTPCNSigmaPr");
  // ----------- 输出文件与输出树（RDF Snapshot 的输出结构） -----------
  TFile fout(path_output_tree, "RECREATE");
  TTree out("jpsi_ref_pairs", "mixed jpsi(A) x ref(B) pairs (Snapshot-like)");

  // 输出的 event-level（这里用 A 的 event-level；如果你想同时存 B
  // 的，额外加分支）
  // event-level（来自 A）
  int o_NumContribCalib;
  float o_fMultTPC, o_fMultTracklets, o_fMultNTracksPV, o_fMultFT0C;
  float o_fPosX, o_fPosY, o_fPosZ;
  int o_fSelection;
  float o_fHadronicRate;

  // jpsi（来自 A）
  float o_jpsi_pt, o_jpsi_eta, o_jpsi_phi, o_jpsi_mass, o_jpsi_sign;

  // e1
  float o_e1_pt, o_e1_eta, o_e1_phi;
  int o_e1_sign;
  float o_e1_ITSChi2NCl, o_e1_TPCNClsCR, o_e1_TPCNClsFound;
  float o_e1_TPCChi2NCl, o_e1_TPCSignal;
  float o_e1_nsig_el, o_e1_nsig_pi, o_e1_nsig_pr;

  // e2
  float o_e2_pt, o_e2_eta, o_e2_phi;
  int o_e2_sign;
  float o_e2_ITSChi2NCl, o_e2_TPCNClsCR, o_e2_TPCNClsFound;
  float o_e2_TPCChi2NCl, o_e2_TPCSignal;
  float o_e2_nsig_el, o_e2_nsig_pi, o_e2_nsig_pr;

  // ref（来自 B）
  float o_ref_pt, o_ref_eta, o_ref_phi;
  float o_ref_ITSChi2NCl, o_ref_TPCNClsCR, o_ref_TPCNClsFound;
  float o_ref_TPCChi2NCl, o_ref_TPCSignal;
  float o_ref_nsig_el, o_ref_nsig_pi, o_ref_nsig_pr;

  out.Branch("NumContribCalib", &o_NumContribCalib);
  out.Branch("fMultTPC", &o_fMultTPC);
  out.Branch("fMultTracklets", &o_fMultTracklets);
  out.Branch("fMultNTracksPV", &o_fMultNTracksPV);
  out.Branch("fMultFT0C", &o_fMultFT0C);
  out.Branch("fPosX", &o_fPosX);
  out.Branch("fPosY", &o_fPosY);
  out.Branch("fPosZ", &o_fPosZ);
  out.Branch("fSelection", &o_fSelection);
  out.Branch("fHadronicRate", &o_fHadronicRate);

  out.Branch("jpsi_pt", &o_jpsi_pt);
  out.Branch("jpsi_eta", &o_jpsi_eta);
  out.Branch("jpsi_phi", &o_jpsi_phi);
  out.Branch("jpsi_mass", &o_jpsi_mass);
  out.Branch("jpsi_sign", &o_jpsi_sign);

  out.Branch("e1_pt", &o_e1_pt);
  out.Branch("e1_eta", &o_e1_eta);
  out.Branch("e1_phi", &o_e1_phi);
  out.Branch("e1_sign", &o_e1_sign);
  out.Branch("e1_ITSChi2NCl", &o_e1_ITSChi2NCl);
  out.Branch("e1_TPCNClsCR", &o_e1_TPCNClsCR);
  out.Branch("e1_TPCNClsFound", &o_e1_TPCNClsFound);
  out.Branch("e1_TPCChi2NCl", &o_e1_TPCChi2NCl);
  out.Branch("e1_TPCSignal", &o_e1_TPCSignal);
  out.Branch("e1_nsig_el", &o_e1_nsig_el);
  out.Branch("e1_nsig_pi", &o_e1_nsig_pi);
  out.Branch("e1_nsig_pr", &o_e1_nsig_pr);

  out.Branch("e2_pt", &o_e2_pt);
  out.Branch("e2_eta", &o_e2_eta);
  out.Branch("e2_phi", &o_e2_phi);
  out.Branch("e2_sign", &o_e2_sign);
  out.Branch("e2_ITSChi2NCl", &o_e2_ITSChi2NCl);
  out.Branch("e2_TPCNClsCR", &o_e2_TPCNClsCR);
  out.Branch("e2_TPCNClsFound", &o_e2_TPCNClsFound);
  out.Branch("e2_TPCChi2NCl", &o_e2_TPCChi2NCl);
  out.Branch("e2_TPCSignal", &o_e2_TPCSignal);
  out.Branch("e2_nsig_el", &o_e2_nsig_el);
  out.Branch("e2_nsig_pi", &o_e2_nsig_pi);
  out.Branch("e2_nsig_pr", &o_e2_nsig_pr);

  out.Branch("ref_pt", &o_ref_pt);
  out.Branch("ref_eta", &o_ref_eta);
  out.Branch("ref_phi", &o_ref_phi);
  out.Branch("ref_ITSChi2NCl", &o_ref_ITSChi2NCl);
  out.Branch("ref_TPCNClsCR", &o_ref_TPCNClsCR);
  out.Branch("ref_TPCNClsFound", &o_ref_TPCNClsFound);
  out.Branch("ref_TPCChi2NCl", &o_ref_TPCChi2NCl);
  out.Branch("ref_TPCSignal", &o_ref_TPCSignal);
  out.Branch("ref_nsig_el", &o_ref_nsig_el);
  out.Branch("ref_nsig_pi", &o_ref_nsig_pi);
  out.Branch("ref_nsig_pr", &o_ref_nsig_pr);

  // ----------- 主循环：对每个 (A,B) 做笛卡尔积配对并输出 -----------
  long long nWritten = 0;
  while (rPairs.Next())
    for (const auto &abPair_single : *abPair) {
      const ULong64_t entryA = abPair_single.first;
      const ULong64_t entryB = abPair_single.second;

      // 读 A（J/psi + e1/e2 + event-level）
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

      const auto &A_jpsi_pt = *fPT;
      const auto &A_jpsi_eta = *fEta;
      const auto &A_jpsi_phi = *fPhi;
      const auto &A_jpsi_mass = *fMass;
      const auto &A_jpsi_sign = *fSign;

      const auto &A_e1_pt = *fPt1;
      const auto &A_e1_eta = *fEta1;
      const auto &A_e1_phi = *fPhi1;
      const auto &A_e1_sign = *fSign1;
      const auto &A_e1_its = *fITSChi2NCl1;
      const auto &A_e1_cr = *fTPCNClsCR1;
      const auto &A_e1_found = *fTPCNClsFound1;
      const auto &A_e1_chi2 = *fTPCChi2NCl1;
      const auto &A_e1_sig = *fTPCSignal1;
      const auto &A_e1_nel = *fTPCNSigmaEl1;
      const auto &A_e1_npi = *fTPCNSigmaPi1;
      const auto &A_e1_npr = *fTPCNSigmaPr1;

      const auto &A_e2_pt = *fPt2;
      const auto &A_e2_eta = *fEta2;
      const auto &A_e2_phi = *fPhi2;
      const auto &A_e2_sign = *fSign2;
      const auto &A_e2_its = *fITSChi2NCl2;
      const auto &A_e2_cr = *fTPCNClsCR2;
      const auto &A_e2_found = *fTPCNClsFound2;
      const auto &A_e2_chi2 = *fTPCChi2NCl2;
      const auto &A_e2_sig = *fTPCSignal2;
      const auto &A_e2_nel = *fTPCNSigmaEl2;
      const auto &A_e2_npi = *fTPCNSigmaPi2;
      const auto &A_e2_npr = *fTPCNSigmaPr2;

      // 读 B（ref tracks）
      rEvt.SetEntry(entryB);

      const auto &B_ref_pt = *fPTREF;
      const auto &B_ref_eta = *fEtaREF;
      const auto &B_ref_phi = *fPhiREF;
      const auto &B_ref_its = *fITSChi2NCl_ref;
      const auto &B_ref_cr = *fTPCNClsCR_ref;
      const auto &B_ref_found = *fTPCNClsFound_ref;
      const auto &B_ref_chi2 = *fTPCChi2NCl_ref;
      const auto &B_ref_sig = *fTPCSignal_ref;
      const auto &B_ref_nel = *fTPCNSigmaEl_ref;
      const auto &B_ref_npi = *fTPCNSigmaPi_ref;
      const auto &B_ref_npr = *fTPCNSigmaPr_ref;

      const long long nA = (long long)A_jpsi_pt.size();
      const long long nB = (long long)B_ref_pt.size();
      long long nPairs = nA * nB;
      if (nPairs == 0) {
        continue;
      }

      // 笛卡尔积：把 A 的每个 jpsi 与 B 的每个 ref 组成一对
      long long filled = 0;
      for (size_t ia = 0; ia < A_jpsi_pt.size(); ++ia) {
        for (size_t ib = 0; ib < B_ref_pt.size(); ++ib) {

          // jpsi (A)
          o_jpsi_pt = A_jpsi_pt[ia];
          o_jpsi_eta = A_jpsi_eta[ia];
          o_jpsi_phi = A_jpsi_phi[ia];
          o_jpsi_mass = A_jpsi_mass[ia];
          o_jpsi_sign = A_jpsi_sign[ia];

          // e1
          o_e1_pt = A_e1_pt[ia];
          o_e1_eta = A_e1_eta[ia];
          o_e1_phi = A_e1_phi[ia];
          o_e1_sign = A_e1_sign[ia];
          o_e1_ITSChi2NCl = A_e1_its[ia];
          o_e1_TPCNClsCR = A_e1_cr[ia];
          o_e1_TPCNClsFound = A_e1_found[ia];
          o_e1_TPCChi2NCl = A_e1_chi2[ia];
          o_e1_TPCSignal = A_e1_sig[ia];
          o_e1_nsig_el = A_e1_nel[ia];
          o_e1_nsig_pi = A_e1_npi[ia];
          o_e1_nsig_pr = A_e1_npr[ia];

          // e2
          o_e2_pt = A_e2_pt[ia];
          o_e2_eta = A_e2_eta[ia];
          o_e2_phi = A_e2_phi[ia];
          o_e2_sign = A_e2_sign[ia];
          o_e2_ITSChi2NCl = A_e2_its[ia];
          o_e2_TPCNClsCR = A_e2_cr[ia];
          o_e2_TPCNClsFound = A_e2_found[ia];
          o_e2_TPCChi2NCl = A_e2_chi2[ia];
          o_e2_TPCSignal = A_e2_sig[ia];
          o_e2_nsig_el = A_e2_nel[ia];
          o_e2_nsig_pi = A_e2_npi[ia];
          o_e2_nsig_pr = A_e2_npr[ia];

          // ref (B)
          o_ref_pt = B_ref_pt[ib];
          o_ref_eta = B_ref_eta[ib];
          o_ref_phi = B_ref_phi[ib];
          o_ref_ITSChi2NCl = B_ref_its[ib];
          o_ref_TPCNClsCR = B_ref_cr[ib];
          o_ref_TPCNClsFound = B_ref_found[ib];
          o_ref_TPCChi2NCl = B_ref_chi2[ib];
          o_ref_TPCSignal = B_ref_sig[ib];
          o_ref_nsig_el = B_ref_nel[ib];
          o_ref_nsig_pi = B_ref_npi[ib];
          o_ref_nsig_pr = B_ref_npr[ib];

          out.Fill(); // 🔥 每一个 (jpsi_i, ref_j) 一条 entry
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

  gROOT->SetBatch(kTRUE); // Disable interactive graphics
  EventMixingIndexGen(path_input_flowVecd, path_input_mult, path_output,
                      path_output_mix);
  EventMixingJpsiAssoPair(path_input_flowVecd, path_input_mult, path_output_mix,
                          path_output_tree);
  return 0;
}