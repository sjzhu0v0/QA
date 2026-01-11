#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "opt/EventData.h"
#include "yaml-cpp/yaml.h"
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <utility>
#include <vector>

using ROOT::RDataFrame;
using ROOT::RVec;
using ROOT::VecOps::Take;

void JpsiAsso(TString path_input_flowVecd = "../input.root",
              TString path_input_mult = "../input2.root",
              TString path_output = "output.root", /* int runNumber = 0, */
              double threshold_bs = 1.) {
  TFile* fOutput = new TFile(path_output, "RECREATE");

  TChain* tree_flowVecd = MRootIO::OpenChain(path_input_flowVecd.Data(), "O2dqflowvecd");
  TChain* tree_mult = MRootIO::OpenChain(path_input_mult.Data(), "MultCalib");
  cout << "Input file: " << path_input_flowVecd << endl;
  cout << "Input file: " << path_input_mult << endl;
  cout << "Output file: " << path_output << endl;
  cout << "Threshold BS: " << threshold_bs << endl;
  YAML::Node config = YAML::LoadFile("config.yaml");
  const double low_edge_deltaPhiToPi = config["hist_binning"]["low_edge_deltaPhiToPi"].as<double>();
  const double up_edge_deltaPhiToPi = config["hist_binning"]["up_edge_deltaPhiToPi"].as<double>();

  tree_flowVecd->AddFriend(tree_mult);

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
      /*  .DefineSlot("NumContribCalib",
                   Calib_NumContrib_fPosZ_Run::NumContribCalibratedFloat,
                   {"fNumContrib", "fPosZ"}) */
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
                             .Define("randTag",
                                     []() {
                                       thread_local TRandom3 random_gen(0);
                                       return random_gen.Uniform(0, 1);
                                     })
      /*  .Filter("isntSelfDefinedPileup", "no self defined pileup") */;

  auto rdf_pair_PartTriggerWithJpsiWithEvent =
      rdf_PartTrigger.Define("pair_indices",
                             [](int njpsi, int nref) {
                               std::vector<std::pair<int, int>> out;
                               out.reserve(njpsi * nref);
                               for (int i = 0; i < njpsi; ++i)
                                 for (int j = 0; j < nref; ++j)
                                   out.emplace_back(i, j);
                               return out;
                             },
                             {"fPT_size", "fPTREF_size"});

  auto df_idx = rdf_pair_PartTriggerWithJpsiWithEvent
                    .Define("jpsi_idx",
                            [](const std::vector<std::pair<int, int>>& v) {
                              RVec<int> o;
                              for (auto& p : v)
                                o.push_back(p.first);
                              return o;
                            },
                            {"pair_indices"})
                    .Define("ref_idx",
                            [](const std::vector<std::pair<int, int>>& v) {
                              RVec<int> o;
                              for (auto& p : v)
                                o.push_back(p.second);
                              return o;
                            },
                            {"pair_indices"});

  auto df_take =
      df_idx
          // J/psi
          .Define("jpsi_pt", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPT", "jpsi_idx"})
          .Define("jpsi_eta", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fEta", "jpsi_idx"})
          .Define("jpsi_phi", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPhi", "jpsi_idx"})
          .Define("jpsi_mass", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fMass", "jpsi_idx"})
          .Define("jpsi_sign", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fSign", "jpsi_idx"})

          // electron 1
          .Define("e1_pt", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPt1", "jpsi_idx"})
          .Define("e1_eta", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fEta1", "jpsi_idx"})
          .Define("e1_phi", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPhi1", "jpsi_idx"})
          .Define("e1_sign", [](const RVec<int>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fSign1", "jpsi_idx"})
          .Define("e1_ITSChi2NCl",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fITSChi2NCl1", "jpsi_idx"})
          .Define("e1_TPCNClsCR",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNClsCR1", "jpsi_idx"})
          .Define("e1_TPCNClsFound",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNClsFound1", "jpsi_idx"})
          .Define("e1_TPCChi2NCl",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCChi2NCl1", "jpsi_idx"})
          .Define("e1_TPCSignal",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCSignal1", "jpsi_idx"})
          .Define("e1_nsig_el", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaEl1", "jpsi_idx"})
          .Define("e1_nsig_pi", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaPi1", "jpsi_idx"})
          .Define("e1_nsig_pr", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaPr1", "jpsi_idx"})

          // electron 2
          .Define("e2_pt", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPt2", "jpsi_idx"})
          .Define("e2_eta", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fEta2", "jpsi_idx"})
          .Define("e2_phi", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPhi2", "jpsi_idx"})
          .Define("e2_sign", [](const RVec<int>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fSign2", "jpsi_idx"})
          .Define("e2_ITSChi2NCl",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fITSChi2NCl2", "jpsi_idx"})
          .Define("e2_TPCNClsCR",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNClsCR2", "jpsi_idx"})
          .Define("e2_TPCNClsFound",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNClsFound2", "jpsi_idx"})
          .Define("e2_TPCChi2NCl",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCChi2NCl2", "jpsi_idx"})
          .Define("e2_TPCSignal",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCSignal2", "jpsi_idx"})
          .Define("e2_nsig_el", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaEl2", "jpsi_idx"})
          .Define("e2_nsig_pi", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaPi2", "jpsi_idx"})
          .Define("e2_nsig_pr", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaPr2", "jpsi_idx"});
  ;
  auto df_all =
      df_take
          .Define("ref_pt", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPTREF", "ref_idx"})
          .Define("ref_eta", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fEtaREF", "ref_idx"})
          .Define("ref_phi", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPhiREF", "ref_idx"})
          .Define("ref_ITSChi2NCl",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fITSChi2NCl", "ref_idx"})
          .Define("ref_TPCNClsCR",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNClsCR", "ref_idx"})
          .Define("ref_TPCNClsFound",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNClsFound", "ref_idx"})
          .Define("ref_TPCChi2NCl",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCChi2NCl", "ref_idx"})
          .Define("ref_TPCSignal",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCSignal", "ref_idx"})
          .Define("ref_nsig_el",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaEl", "ref_idx"})
          .Define("ref_nsig_pi",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaPi", "ref_idx"})
          .Define("ref_nsig_pr",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaPr", "ref_idx"})
          .Define("ref_dcaxy", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fDcaXY", "ref_idx"})
          .Define("ref_dcaz", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fDcaZ", "ref_idx"});

  if (is_interactive())
    ROOT::RDF::Experimental::AddProgressBar(df_all);
  df_all.Snapshot("jpsi_ref_pairs", path_output,
                  {"randTag",        "NumContribCalib",  "fMultTPC",        "fMultTracklets",
                   "fMultNTracksPV", "fMultFT0C",        "fPosX",           "fPosY",
                   "fPosZ",          "fSelection",       "fHadronicRate",   "jpsi_pt",
                   "jpsi_eta",       "jpsi_phi",         "jpsi_mass",       "jpsi_sign",
                   "e1_pt",          "e1_eta",           "e1_phi",          "e1_sign",
                   "e1_ITSChi2NCl",  "e1_TPCNClsCR",     "e1_TPCNClsFound", "e1_TPCChi2NCl",
                   "e1_TPCSignal",   "e1_nsig_el",       "e1_nsig_pi",      "e1_nsig_pr",
                   "e2_pt",          "e2_eta",           "e2_phi",          "e2_sign",
                   "e2_ITSChi2NCl",  "e2_TPCNClsCR",     "e2_TPCNClsFound", "e2_TPCChi2NCl",
                   "e2_TPCSignal",   "e2_nsig_el",       "e2_nsig_pi",      "e2_nsig_pr",
                   "ref_pt",         "ref_eta",          "ref_phi",         "ref_ITSChi2NCl",
                   "ref_TPCNClsCR",  "ref_TPCNClsFound", "ref_TPCChi2NCl",  "ref_TPCSignal",
                   "ref_nsig_el",    "ref_nsig_pi",      "ref_nsig_pr",     "ref_dcaxy",
                   "ref_dcaz"});
}

int main(int argc, char** argv) {
  TString path_input_flowVecd = "../input.root";
  TString path_input_mult = "../input2.root";
  TString path_output = "output.root";

  if (argc > 1) {
    path_input_flowVecd = argv[1];
  }
  if (argc > 2) {
    path_input_mult = argv[2];
  }
  if (argc > 3) {
    path_output = argv[3];
  }

  JpsiAsso(path_input_flowVecd, path_input_mult, path_output);

  return 0;
}
