#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MEventMixing.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "opt/EventData.h"
#include "yaml-cpp/yaml.h"
#include <ROOT/RDataFrame.hxx>

void EventMixingReading(TString path_input_flowVecd = "../input.root",
                        TString path_output = "output.root") {
  TFile *file_flowVecd = TFile::Open(path_input_flowVecd);
  TFile *fOutput = new TFile(path_output, "RECREATE");
  YAML::Node config = YAML::LoadFile("config.yaml");

  const double low_edge_deltaPhiToPi =
      config["hist_binning"]["low_edge_deltaPhiToPi"].as<double>();
  const double up_edge_deltaPhiToPi =
      config["hist_binning"]["up_edge_deltaPhiToPi"].as<double>();

  TTree *tree_input = (TTree *)file_flowVecd->Get("EventMixing");

  ROOT::RDataFrame rdf(*tree_input);

  auto CutTrackInfo = [](const TrackInfo &track_info, const int &index) {
    bool ptCut =
        track_info.fPTREF[index] > 0.4 && track_info.fPTREF[index] < 4.0;
    return ptCut;
  };

  auto rdf_AllVar =
      rdf.Define("JpsiInfoUS",
                 [&CutTrackInfo, &low_edge_deltaPhiToPi, &up_edge_deltaPhiToPi](
                     const vector<EventDataREF> &vec_eventData) {
                   ROOT::VecOps::RVec<array<float, 4>> vec2return;
                   for (const auto &eventData : vec_eventData) {
                     for (size_t i = 0; i < eventData.track_info.fEtaREF.size();
                          ++i)
                       for (size_t j = 0;
                            j < eventData.track_info2.fEtaREF.size(); ++j) {
                         if (!CutTrackInfo(eventData.track_info, j))
                           continue;
                         if (!CutTrackInfo(eventData.track_info2, i))
                           continue;
                         float delta_eta = eventData.track_info2.fEtaREF[i] -
                                           eventData.track_info.fEtaREF[j];
                         float delta_phi = eventData.track_info2.fPhiREF[i] -
                                           eventData.track_info.fPhiREF[j];
                         int n = 0;
                         while (delta_phi > up_edge_deltaPhiToPi * M_PI &&
                                n < 10) {
                           n++;
                           delta_phi -= 2 * M_PI;
                         }
                         while (delta_phi < low_edge_deltaPhiToPi * M_PI &&
                                n < 10) {
                           n++;
                           delta_phi += 2 * M_PI;
                         }
                         if (n >= 10)
                           delta_phi = -999.;
                         vec2return.push_back(
                             {eventData.event_info.fPosZ,
                              eventData.event_info.fNumContribCalib, delta_eta,
                              delta_phi});
                       }
                   }
                 },
                 {"MixedEvent"})
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
                  {"RefInfo"});
  // ROOT::RDF::Experimental::AddProgressBar(rdf);

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
      "NumContribCalibUSSingle", "N_{vtx contrib} Calibrated", "", 10,
      {0, 7, 12, 17, 22, 29, 37, 46, 57, 73, 300});

#define obj2push_thnd(rdf2push, ...)                                           \
  do {                                                                         \
    TupleTHnDModel tuple_thnd = GetTHnDModelWithTitle(__VA_ARGS__);            \
    gRResultHandles.push_back(                                                 \
        rdf2push.HistoND(get<0>(tuple_thnd), get<1>(tuple_thnd)));             \
  } while (0)
  obj2push_thnd(rdf_AllVar, {DeltaEtaRef, DeltaPhiRef, var_fPosZ,
                             var_NumContribCalibBinned});
  fOutput->cd();
  RResultWrite(gRResultHandles);
  fOutput->Close();
  cout << "Output written to " << path_output << endl;
}

int main(int argc, char **argv) {
  TString path_input_flowVecd = argv[1];
  TString path_output = argv[2];
  EventMixingReading(path_input_flowVecd, path_output);
  return 0;
}