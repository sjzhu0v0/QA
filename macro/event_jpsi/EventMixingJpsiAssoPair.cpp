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
#include <iomanip>
#include <iostream>
#include <iostream>
#include <fstream>

long getRSS() {
    std::ifstream statm("/proc/self/statm");
    long size, resident;
    statm >> size >> resident;
    return resident * sysconf(_SC_PAGESIZE) / 1024 / 1024;
}


template <typename T> std::vector<T> makeVec(const TTreeReaderArray<T> &arr) {
  return std::vector<T>(arr.begin(), arr.end());
}

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
  out.SetAutoFlush(10000);   // 每 1e4 行写一次
  out.SetAutoSave(30000000); // 30MB 保存一次

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
  float o_e2_ITSChi2NCl, o_e2_TPCNClsCR, o_e2_TPCNClsFound;
  float o_e2_TPCChi2NCl, o_e2_TPCSignal;
  float o_e2_nsig_el, o_e2_nsig_pi, o_e2_nsig_pr;
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
  long long nWritten = 0;

  bool isInteractive = is_interactive();
  long long nEntries = rPairs.GetEntries();
  for (long long iEntry = 0; iEntry < nEntries; ++iEntry) {
    if (isInteractive)
      // print progress bar
      if (iEntry % (100) == 0) {
          printf(
              "\rProcessing entry %lld / %lld  (RSS = %ld MB) ", iEntry,
              nEntries, getRSS());
        /*float progress = (float)iEntry / nEntries;
        int barWidth = 70;
        std::cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
          if (i < pos)
            std::cout << "=";
          else if (i == pos)
            std::cout << ">";
          else
            std::cout << " ";
        }*/
        // std::cout << "] " << int(progress * 100.0) << " %\r";
        // std::cout.flush();
      }
   // rPairs.SetEntry(iEntry);
    for (const auto &abPair_single : *abPair) {
      ULong64_t entryA = abPair_single.first;
      ULong64_t entryB = abPair_single.second;
      int nPairs = fPT.GetSize()*fPTREF.GetSize();
      if (nPairs == 0) {
        continue;
      }
      // rEvt.SetEntry(entryA);
      /*o_NumContribCalib = *NumContribCalib;
      o_fMultTPC = *fMultTPC;
      o_fMultTracklets = *fMultTracklets;
      o_fMultNTracksPV = *fMultNTracksPV;
      o_fMultFT0C = *fMultFT0C;
      o_fPosX = *fPosX;
      o_fPosY = *fPosY;
      o_fPosZ = *fPosZ;
      o_fSelection = *fSelection;
      o_fHadronicRate = *fHadronicRate;
      auto A_jpsi_pt = makeVec(fPT); 
      auto A_jpsi_eta = makeVec(fEta);
      auto A_jpsi_phi = makeVec(fPhi);
      auto A_jpsi_mass = makeVec(fMass);
      auto A_jpsi_sign = makeVec(fSign);

      auto A_e1_pt = makeVec(fPt1);
      auto A_e1_eta = makeVec(fEta1);
      auto A_e1_phi = makeVec(fPhi1);
      auto A_e1_sign = makeVec(fSign1);
      auto A_e1_its = makeVec(fITSChi2NCl1);
      auto A_e1_cr = makeVec(fTPCNClsCR1);
      auto A_e1_found = makeVec(fTPCNClsFound1);
      auto A_e1_chi2 = makeVec(fTPCChi2NCl1);
      auto A_e1_sig = makeVec(fTPCSignal1);
      auto A_e1_nel = makeVec(fTPCNSigmaEl1);
      auto A_e1_npi = makeVec(fTPCNSigmaPi1);
      auto A_e1_npr = makeVec(fTPCNSigmaPr1);

      auto A_e2_pt = makeVec(fPt2);
      auto A_e2_eta = makeVec(fEta2);
      auto A_e2_phi = makeVec(fPhi2);
      auto A_e2_sign = makeVec(fSign2);
      auto A_e2_its = makeVec(fITSChi2NCl2);
      auto A_e2_cr = makeVec(fTPCNClsCR2);
      auto A_e2_found = makeVec(fTPCNClsFound2);
      auto A_e2_chi2 = makeVec(fTPCChi2NCl2);
      auto A_e2_sig = makeVec(fTPCSignal2);
      auto A_e2_nel = makeVec(fTPCNSigmaEl2);
      auto A_e2_npi = makeVec(fTPCNSigmaPi2);
      auto A_e2_npr = makeVec(fTPCNSigmaPr2);

      rEvt.SetEntry(entryB);
      auto B_ref_pt = makeVec(fPTREF);
      auto B_ref_eta = makeVec(fEtaREF);
      auto B_ref_phi = makeVec(fPhiREF);
      auto B_ref_its = makeVec(fITSChi2NCl_ref);
      auto B_ref_cr = makeVec(fTPCNClsCR_ref);
      auto B_ref_found = makeVec(fTPCNClsFound_ref);
      auto B_ref_chi2 = makeVec(fTPCChi2NCl_ref);
      auto B_ref_sig = makeVec(fTPCSignal_ref);
      auto B_ref_nel = makeVec(fTPCNSigmaEl_ref);
      auto B_ref_npi = makeVec(fTPCNSigmaPi_ref);
      auto B_ref_npr = makeVec(fTPCNSigmaPr_ref);*/

      // long long filled = 0;
     /* for (size_t ia = 0; ia < A_jpsi_pt.size(); ++ia) {
        for (size_t ib = 0; ib < B_ref_pt.size(); ++ib) {
          o_jpsi_pt = A_jpsi_pt[ia];
          o_jpsi_eta = A_jpsi_eta[ia];
          o_jpsi_phi = A_jpsi_phi[ia];
          o_jpsi_mass = A_jpsi_mass[ia];
          o_jpsi_sign = A_jpsi_sign[ia];
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
          out.Fill();
        }
      }*/
     /* A_jpsi_pt.clear();
        A_jpsi_eta.clear();
        A_jpsi_phi.clear();
        A_jpsi_mass.clear();
        A_jpsi_sign.clear();
        A_e1_pt.clear();
        A_e1_eta.clear();
        A_e1_phi.clear();
        A_e1_sign.clear();
        A_e1_its.clear();
        A_e1_cr.clear();
        A_e1_found.clear();
        A_e1_chi2.clear();
        A_e1_sig.clear();
        A_e1_nel.clear();
        A_e1_npi.clear();
        A_e1_npr.clear();
        A_e2_pt.clear();
        A_e2_eta.clear();
        A_e2_phi.clear();
        A_e2_sign.clear();
        A_e2_its.clear();
        A_e2_cr.clear();
        A_e2_found.clear();
        A_e2_chi2.clear();
        A_e2_sig.clear();
        A_e2_nel.clear();
        A_e2_npi.clear();
        A_e2_npr.clear();
        B_ref_pt.clear();
        B_ref_eta.clear();
        B_ref_phi.clear();
        B_ref_its.clear();
        B_ref_cr.clear();
        B_ref_found.clear();
        B_ref_chi2.clear();
        B_ref_sig.clear();
        B_ref_nel.clear();
        B_ref_npi.clear();
        B_ref_npr.clear();*/
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
