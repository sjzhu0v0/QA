#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "opt/EventData.h"
#include "yaml-cpp/yaml.h"
#include <ROOT/RDataFrame.hxx>

void Ref_BS(TString path_input_flowVecd = "../input.root",
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

  auto CutTrackInfo = [](const TrackInfo& track_info, const int& index) {
    bool ptCut = track_info.fPTREF[index] > 0.4 && track_info.fPTREF[index] < 4.0;
    return ptCut;
  };

  auto rdf_PartTriggerWithJpsiWithEvent =
      rdf_PartTrigger.Define("pair_indices",
                             [](int ntrk) {
                               std::vector<std::pair<int, int>> out;
                               out.reserve((size_t)ntrk * (ntrk - 1) / 2);
                               for (int i = 0; i < ntrk; ++i)
                                 for (int j = i + 1; j < ntrk; ++j)
                                   out.emplace_back(i, j);
                               return out;
                             },
                             {"fPTREF_size"});

  auto df_idx = rdf_PartTriggerWithJpsiWithEvent
                    .Define("ref1_idx",
                            [](const std::vector<std::pair<int, int>>& v) {
                              RVec<int> o;
                              o.reserve(v.size());
                              for (auto& p : v)
                                o.push_back(p.first);
                              return o;
                            },
                            {"pair_indices"})
                    .Define("ref2_idx",
                            [](const std::vector<std::pair<int, int>>& v) {
                              RVec<int> o;
                              o.reserve(v.size());
                              for (auto& p : v)
                                o.push_back(p.second);
                              return o;
                            },
                            {"pair_indices"});
  auto df_all =
      df_idx
          .Define("ref1_pt", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPTREF", "ref1_idx"})
          .Define("ref1_eta", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fEtaREF", "ref1_idx"})
          .Define("ref1_phi", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPhiREF", "ref1_idx"})
          .Define("ref1_ITSChi2NCl",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fITSChi2NCl", "ref1_idx"})
          .Define("ref1_TPCNClsCR",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNClsCR", "ref1_idx"})
          .Define("ref1_TPCNClsFound",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNClsFound", "ref1_idx"})
          .Define("ref1_TPCChi2NCl",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCChi2NCl", "ref1_idx"})
          .Define("ref1_ITSClusterMap",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fITSClusterMap", "ref1_idx"})
          .Define("ref1_TPCSignal",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCSignal", "ref1_idx"})
          .Define("ref1_nsig_el",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaEl", "ref1_idx"})
          .Define("ref1_nsig_pi",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaPi", "ref1_idx"})
          .Define("ref1_nsig_pr",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaPr", "ref1_idx"})
          .Define("ref1_dcaxy", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fDcaXY", "ref1_idx"})
          .Define("ref1_dcaz", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fDcaZ", "ref1_idx"})
          .Define("ref2_pt", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPTREF", "ref2_idx"})
          .Define("ref2_eta", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fEtaREF", "ref2_idx"})
          .Define("ref2_phi", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPhiREF", "ref2_idx"})
          .Define("ref2_ITSChi2NCl",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fITSChi2NCl", "ref2_idx"})
          .Define("ref2_TPCNClsCR",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNClsCR", "ref2_idx"})
          .Define("ref2_TPCNClsFound",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNClsFound", "ref2_idx"})
          .Define("ref2_TPCChi2NCl",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCChi2NCl", "ref2_idx"})
          .Define("ref2_ITSClusterMap",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fITSClusterMap", "ref2_idx"})
          .Define("ref2_TPCSignal",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCSignal", "ref2_idx"})
          .Define("ref2_nsig_el",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaEl", "ref2_idx"})
          .Define("ref2_nsig_pi",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaPi", "ref2_idx"})
          .Define("ref2_nsig_pr",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaPr", "ref2_idx"})
          .Define("ref2_dcaxy", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fDcaXY", "ref2_idx"})
          .Define("ref2_dcaz", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fDcaZ", "ref2_idx"});

  if (is_interactive())
    ROOT::RDF::Experimental::AddProgressBar(rdf);
  df_all.Snapshot("jpsi_ref_pairs", path_output,
                  {"randTag",
                   "NumContribCalib",
                   "fMultTPC",
                   "fMultTracklets",
                   "fMultNTracksPV",
                   "fMultFT0C",
                   "fPosX",
                   "fPosY",
                   "fPosZ",
                   "fSelection",
                   "fHadronicRate",
                   "ref1_pt",
                   "ref1_eta",
                   "ref1_phi",
                   "ref1_ITSChi2NCl",
                   "ref1_TPCNClsCR",
                   "ref1_TPCNClsFound",
                   "ref1_TPCChi2NCl",
                   "ref1_ITSClusterMap",
                   "ref1_TPCSignal",
                   "ref1_nsig_el",
                   "ref1_nsig_pi",
                   "ref1_nsig_pr",
                   "ref1_dcaxy",
                   "ref1_dcaz",
                   "ref2_pt",
                   "ref2_eta",
                   "ref2_phi",
                   "ref2_ITSChi2NCl",
                   "ref2_TPCNClsCR",
                   "ref2_TPCNClsFound",
                   "ref2_TPCChi2NCl",
                   "ref2_ITSClusterMap",
                   "ref2_TPCSignal",
                   "ref2_nsig_el",
                   "ref2_nsig_pi",
                   "ref2_nsig_pr",
                   "ref2_dcaxy",
                   "ref2_dcaz"});
}

int main(int argc, char** argv) {
  TString path_input_flowVecd = "../input.root";
  TString path_input_mult = "../input2.root";
  TString path_output = "output.root";
  double threshold_bs = 1.;

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
    threshold_bs = atof(argv[4]);
  }

  Ref_BS(path_input_flowVecd, path_input_mult, path_output, threshold_bs);

  return 0;
}
