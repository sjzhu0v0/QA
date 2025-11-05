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
  TFile *fOutput = new TFile(path_output, "RECREATE");

  TChain *tree_flowVecd =
      MRootIO::OpenChain(path_input_flowVecd.Data(), "O2dqflowvecd");
  TChain *tree_mult = MRootIO::OpenChain(path_input_mult.Data(), "MultCalib");
  cout << "Input file: " << path_input_flowVecd << endl;
  cout << "Input file: " << path_input_mult << endl;
  cout << "Output file: " << path_output << endl;
  cout << "Threshold BS: " << threshold_bs << endl;
  YAML::Node config = YAML::LoadFile("config.yaml");
  const double low_edge_deltaPhiToPi =
      config["hist_binning"]["low_edge_deltaPhiToPi"].as<double>();
  const double up_edge_deltaPhiToPi =
      config["hist_binning"]["up_edge_deltaPhiToPi"].as<double>();

  tree_flowVecd->AddFriend(tree_mult);

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
  auto rdf_isTriggerTVX =
      rdf_witTrigger.Filter("isTriggerTVX", "is Trigger TVX");
  auto rdf_PartTrigger =
      rdf_witTrigger.Filter("isTriggerTVX", "is Trigger TVX")
          .Filter("isntITSROFrameBorder", "no ITS RO Frame border")
          .Filter("isntTimeFrameBorder", "no Time Frame border")
          .Filter("isntSameBunchPileup", "no Time Frame border")
      /*  .Filter("isntSelfDefinedPileup", "no self defined pileup") */;

  auto CutTrackInfo = [](const TrackInfo &track_info, const int &index) {
    bool ptCut =
        track_info.fPTREF[index] > 0.4 && track_info.fPTREF[index] < 4.0;
    return ptCut;
  };

  auto rdf_PartTriggerWithJpsiWithEvent =
      rdf_PartTrigger
          .Define("EventData", CreateEventData,
                  {"fMultTPC", "fMultTracklets", "fMultNTracksPV", "fMultFT0C",
                   "fNumContrib", "NumContribCalib", "fPosX", "fPosY", "fPosZ",
                   "fSelection", "fHadronicRate", "fPT", "fEta", "fPhi",
                   "fMass", "fSign", "fPTREF", "fEtaREF", "fPhiREF"})
          .DefineSlot("Cut_BS",
                      [threshold_bs](unsigned int) {
                        thread_local TRandom3 random_gen(0);
                        return random_gen.Uniform(0, 1) < threshold_bs;
                      })
          .Define(
              "RefInfo",
              [&CutTrackInfo, &low_edge_deltaPhiToPi,
               &up_edge_deltaPhiToPi](const EventData &eventData) {
                ROOT::VecOps::RVec<array<float, 4>> vec2return;
                for (size_t i = 0; i < eventData.track_info.fEtaREF.size(); ++i)
                  for (size_t j = i + 1;
                       j < eventData.track_info.fEtaREF.size(); ++j) {
                    if (!CutTrackInfo(eventData.track_info, j))
                      continue;
                    if (!CutTrackInfo(eventData.track_info, i))
                      continue;
                    float delta_eta = eventData.track_info.fEtaREF[i] -
                                      eventData.track_info.fEtaREF[j];
                    float delta_phi = eventData.track_info.fPhiREF[i] -
                                      eventData.track_info.fPhiREF[j];
                    int n = 0;
                    while (delta_phi > up_edge_deltaPhiToPi * M_PI && n < 10) {
                      n++;
                      delta_phi -= 2 * M_PI;
                    }
                    while (delta_phi < low_edge_deltaPhiToPi * M_PI && n < 10) {
                      n++;
                      delta_phi += 2 * M_PI;
                    }
                    if (n >= 10)
                      delta_phi = -999.;
                    vec2return.push_back({eventData.event_info.fPosZ,
                                          eventData.event_info.fNumContribCalib,
                                          delta_eta, delta_phi});
                  }
                return vec2return;
              },
              {"EventData"})
          .Define("RefInfoSingle",
                  [&CutTrackInfo](const EventData &eventData) {
                    ROOT::VecOps::RVec<array<float, 2>> vec2return;
                    for (size_t i = 0; i < eventData.track_info.fEtaREF.size();
                         ++i) {
                      if (!CutTrackInfo(eventData.track_info, i))
                        continue;
                      vec2return.push_back(
                          {eventData.event_info.fPosZ,
                           eventData.event_info.fNumContribCalib});
                    }
                    return vec2return;
                  },
                  {"EventData"})
          .Define("DeltaPhiRef",
                  [](const ROOT::VecOps::RVec<array<float, 4>> &refInfo) {
                    ROOT::VecOps::RVec<float> deltaPhiRef;
                    for (const auto &pair : refInfo) {
                      deltaPhiRef.push_back(pair[3]);
                    }
                    return deltaPhiRef;
                  },
                  {"RefInfo"})
          .Define("DeltaEtaRef",
                  [](const ROOT::VecOps::RVec<array<float, 4>> &refInfo) {
                    ROOT::VecOps::RVec<float> deltaEtaRef;
                    for (const auto &pair : refInfo) {
                      deltaEtaRef.push_back(pair[2]);
                    }
                    return deltaEtaRef;
                  },
                  {"RefInfo"})
          .Define("PosZRef",
                  [](const ROOT::VecOps::RVec<array<float, 4>> &refInfo) {
                    ROOT::VecOps::RVec<float> posZRef;
                    for (const auto &pair : refInfo) {
                      posZRef.push_back(pair[0]);
                    }
                    return posZRef;
                  },
                  {"RefInfo"})
          .Define("NumContribCalibRef",
                  [](const ROOT::VecOps::RVec<array<float, 4>> &refInfo) {
                    ROOT::VecOps::RVec<float> numContribCalibRef;
                    for (const auto &pair : refInfo) {
                      numContribCalibRef.push_back(pair[1]);
                    }
                    return numContribCalibRef;
                  },
                  {"RefInfo"})
          .Define("PosZRefSingle",
                  [](const ROOT::VecOps::RVec<array<float, 2>> &refInfo) {
                    ROOT::VecOps::RVec<float> posZRefSingle;
                    for (const auto &pair : refInfo) {
                      posZRefSingle.push_back(pair[0]);
                    }
                    return posZRefSingle;
                  },
                  {"RefInfoSingle"})
          .Define("NumContribCalibRefSingle",
                  [](const ROOT::VecOps::RVec<array<float, 2>> &refInfo) {
                    ROOT::VecOps::RVec<float> numContribCalibRefSingle;
                    for (const auto &pair : refInfo) {
                      numContribCalibRefSingle.push_back(pair[1]);
                    }
                    return numContribCalibRefSingle;
                  },
                  {"RefInfoSingle"})
          .Filter("Cut_BS", "Bootstrap cut for Jpsi association");
  if (is_interactive())
    ROOT::RDF::Experimental::AddProgressBar(rdf);
  StrVar4Hist var_fPosZ("PosZRef", "#it{V}_{Z}", "cm", 8, {-10, 10});
  StrVar4Hist var_NumContribCalibBinned(
      "NumContribCalibRef", "N_{vtx contrib} Calibrated", "", 10,
      {0, 7, 12, 17, 22, 29, 37, 46, 57, 73, 300});
  int n_bins_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["n_bins"].as<int>();
  double min_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["min"].as<double>();
  double max_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["max"].as<double>();
  StrVar4Hist DeltaEtaRef("DeltaEtaRef", "#Delta#eta_{track, track}", "",
                          n_bins_deltaEta_assoYield,
                          {min_deltaEta_assoYield, max_deltaEta_assoYield});
  int n_bins_deltaPhi_assoYield =
      config["hist_binning"]["n_bins_deltaPhi_assoYield"].as<int>();
  StrVar4Hist DeltaPhiRef(
      "DeltaPhiRef", "#Delta#phi_{track, track}", "", n_bins_deltaPhi_assoYield,
      {low_edge_deltaPhiToPi * M_PI, up_edge_deltaPhiToPi * M_PI});
  StrVar4Hist var_fPosZSingle("PosZRefSingle", "#it{V}_{Z}", "cm", 8,
                              {-10, 10});
  StrVar4Hist var_NumContribCalibBinnedSingle(
      "NumContribCalibRefSingle", "N_{vtx contrib} Calibrated", "", 10,
      {0, 7, 12, 17, 22, 29, 37, 46, 57, 73, 300});

#define obj2push_thnd(rdf2push, ...)                                           \
  do {                                                                         \
    TupleTHnDModel tuple_thnd = GetTHnDModelWithTitle(__VA_ARGS__);            \
    gRResultHandles.push_back(                                                 \
        rdf2push.HistoND(get<0>(tuple_thnd), get<1>(tuple_thnd)));             \
  } while (0)
  // cout << "Start pushing THnD objects..." << endl;
  obj2push_thnd(
      rdf_PartTriggerWithJpsiWithEvent,
      {DeltaEtaRef, DeltaPhiRef, var_fPosZ, var_NumContribCalibBinned});
  obj2push_thnd(rdf_PartTriggerWithJpsiWithEvent,
                {var_fPosZSingle, var_NumContribCalibBinnedSingle});

  RunGraphs(gRResultHandles);

  fOutput->cd();
  RResultWrite(gRResultHandles);
  fOutput->Close();
}

int main(int argc, char **argv) {
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
