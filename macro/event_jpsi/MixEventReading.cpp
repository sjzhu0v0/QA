#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MEventMixing.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "opt/EventData.h"
#include <ROOT/RDataFrame.hxx>

void EventMixingReading(TString path_input_flowVecd = "../input.root",
                        TString path_output = "output.root") {
  TFile *file_flowVecd = TFile::Open(path_input_flowVecd);
  TFile *fOutput = new TFile(path_output, "RECREATE");

  TTree *tree_input = (TTree *)file_flowVecd->Get("EventMixing");

  ROOT::RDataFrame rdf(*tree_input);

  auto CutTrackInfo = [](const TrackInfo &track_info, const int &index) {
    bool ptCut =
        track_info.fPTREF[index] > 0.4 && track_info.fPTREF[index] < 4.0;
    return ptCut;
  };

  auto rdf_AllVar =
      rdf.Define(
             "JpsiInfoUS",
             [&CutTrackInfo](const vector<EventData> &vec_eventData) {
               ROOT::VecOps::RVec<array<float, 6>> vec2return;
               for (auto eventData : vec_eventData)
                 for (size_t i = 0; i < eventData.jpsi_info.fPT.size(); ++i) {
                   if (eventData.jpsi_info.fSign[i] == 0) {
                     for (size_t j = 0; j < eventData.track_info.fEtaREF.size();
                          ++j) {
                       if (!CutTrackInfo(eventData.track_info, j))
                         continue;

                       float delta_eta = eventData.jpsi_info.fEta[i] -
                                         eventData.track_info.fEtaREF[j];
                       float delta_phi = eventData.jpsi_info.fPhi[i] -
                                         eventData.track_info.fPhiREF[j];
                       int n = 0;
                       while (delta_phi > 1.5 * M_PI && n < 10) {
                         n++;
                         delta_phi -= 2 * M_PI;
                       }
                       while (delta_phi < -0.5 * M_PI && n < 10) {
                         n++;
                         delta_phi += 2 * M_PI;
                       }
                       if (n >= 10)
                         delta_phi = -999.;
                       vec2return.push_back(
                           {eventData.event_info.fPosZ,
                            eventData.event_info.fNumContribCalib,
                            eventData.jpsi_info.fMass[i],
                            eventData.jpsi_info.fPT[i], delta_eta, delta_phi});
                     }
                   }
                 }
               return vec2return;
             },
             {"MixedEvent"})
          .Define("PosZUS",
                  [](const ROOT::VecOps::RVec<array<float, 6>> &jpsiInfoUS) {
                    ROOT::VecOps::RVec<float> posZUS;
                    for (const auto &pair : jpsiInfoUS) {
                      posZUS.push_back(pair[0]);
                    }
                    return posZUS;
                  },
                  {"JpsiInfoUS"})
          .Define("NumContribCalibUS",
                  [](const ROOT::VecOps::RVec<array<float, 6>> &jpsiInfoUS) {
                    ROOT::VecOps::RVec<float> numContribCalibUS;
                    for (const auto &pair : jpsiInfoUS) {
                      numContribCalibUS.push_back(pair[1]);
                    }
                    return numContribCalibUS;
                  },
                  {"JpsiInfoUS"})
          .Define("MassUS",
                  [](const ROOT::VecOps::RVec<array<float, 6>> &jpsiInfoUS) {
                    ROOT::VecOps::RVec<float> massUS;
                    for (const auto &pair : jpsiInfoUS) {
                      massUS.push_back(pair[2]);
                    }
                    return massUS;
                  },
                  {"JpsiInfoUS"})
          .Define("PtUS",
                  [](const ROOT::VecOps::RVec<array<float, 6>> &jpsiInfoUS) {
                    ROOT::VecOps::RVec<float> ptUS;
                    for (const auto &pair : jpsiInfoUS) {
                      ptUS.push_back(pair[3]);
                    }
                    return ptUS;
                  },
                  {"JpsiInfoUS"})
          .Define("DeltaEtaUS",
                  [](const ROOT::VecOps::RVec<array<float, 6>> &jpsiInfoUS) {
                    ROOT::VecOps::RVec<float> deltaEtaUS;
                    for (const auto &pair : jpsiInfoUS) {
                      deltaEtaUS.push_back(pair[4]);
                    }
                    return deltaEtaUS;
                  },
                  {"JpsiInfoUS"})
          .Define("DeltaPhiUS",
                  [](const ROOT::VecOps::RVec<array<float, 6>> &jpsiInfoUS) {
                    ROOT::VecOps::RVec<float> deltaPhiUS;
                    for (const auto &pair : jpsiInfoUS) {
                      deltaPhiUS.push_back(pair[5]);
                    }
                    return deltaPhiUS;
                  },
                  {"JpsiInfoUS"});

  StrVar4Hist var_fPosZ("PosZUS", "#it{V}_{Z}", "cm", 8, {-10, 10});
  StrVar4Hist var_NumContribCalibBinned(
      "NumContribCalibUS", "N_{vtx contrib} Calibrated", "", 10,
      {0, 23, 31, 37, 43, 48, 54, 61, 69, 81, 297});
  StrVar4Hist var_MassJpsiCandidate("MassUS", "M_{ee}", "GeV^{2}/c^{4}", 90,
                                    {1.8, 5.4});
  StrVar4Hist var_PtJpsiCandidate("PtUS", "p_{T}", "GeV/c", 10, {0., 10.});
  StrVar4Hist var_DeltaEtaUS("DeltaEtaUS", "#Delta#eta_{J/#psi, track}", "", 80,
                             {-4., 4.});
  StrVar4Hist var_DeltaPhiUS("DeltaPhiUS", "#Delta#phi_{J/#psi, track}", "", 10,
                             {-M_PI_2, M_PI + M_PI_2});

#define obj2push_thnd(rdf2push, ...)                                           \
  do {                                                                         \
    TupleTHnDModel tuple_thnd = GetTHnDModelWithTitle(__VA_ARGS__);            \
    gRResultHandles.push_back(                                                 \
        rdf2push.HistoND(get<0>(tuple_thnd), get<1>(tuple_thnd)));             \
  } while (0)

  obj2push_thnd(rdf_AllVar, {var_DeltaEtaUS, var_DeltaPhiUS, var_fPosZ,
                             var_MassJpsiCandidate, var_PtJpsiCandidate,
                             var_NumContribCalibBinned});
  fOutput->cd();
  RResultWrite(gRResultHandles);
  fOutput->Close();
}

int main(int argc, char **argv) {
  TString path_input_flowVecd = argv[1];
  TString path_output = argv[2];
  EventMixingReading(path_input_flowVecd, path_output);
  return 0;
}