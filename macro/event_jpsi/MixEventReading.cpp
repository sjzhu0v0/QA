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
      {0, 23, 31, 37, 43, 48, 54, 61, 69, 81, 297});
  StrVar4Hist var_fMultTPC("fMultTPC", "Mult_{TPC}", "", 600, {0, 600});
  StrVar4Hist var_fMultREF("fMultREF", "Mult_{REF}", "", 100, {0, 100});
  StrVar4Hist var_fMultFT0C("fMultFT0C", "Mult_{FT0C}", "", 130,
                            {-1000., 12000.});
  StrVar4Hist var_MultNTracksPV("fMultNTracksPV", "#it{N}_{Tracks PV}", "", 150,
                                {0, 150});
  StrVar4Hist var_MassJpsiCandidate("fMass", "M_{ee}", "GeV^{2}/c^{4}", 100,
                                    {1., 5.});
  StrVar4Hist var_PtJpsiCandidate("fPT", "p_{T}", "GeV/c", 10, {0., 10.});

  const vector<double> bins_mix_numContrib = var_NumContribCalibBinned.fBins;
  const vector<double> bins_mix_posZ = var_fPosZMix.fBins;

  TFile *file_flowVecd = TFile::Open(path_input_flowVecd);
  TFile *fOutput = new TFile(path_output, "RECREATE");

  TTree *tree_input = (TTree *)file_flowVecd->Get("EventMixing");

  ROOT::RDataFrame rdf(*tree_input);

  auto rdf_AllVar =
      rdf.Define("fPosZ1",
                 [](const vector<EventData> &events) {
                   vector<float> vec2return;
                   for (const auto &event : events) {
                     vec2return.push_back(event.event_info.fPosZ);
                   }
                   return ROOT::VecOps::RVec<float>(vec2return);
                 },
                 {"MixedEvent"})
          .Define("fPosZ2",
                  [](const vector<EventData> &events) {
                    vector<float> vec2return;
                    for (const auto &event : events) {
                      vec2return.push_back(event.event_info2.fPosZ);
                    }
                    return ROOT::VecOps::RVec<float>(vec2return);
                  },
                  {"MixedEvent"})
          .Define("NumContribCalib1",
                  [](const vector<EventData> &events) {
                    vector<float> vec2return;
                    for (const auto &event : events) {
                      vec2return.push_back(event.event_info.fNumContrib);
                    }
                    return ROOT::VecOps::RVec<float>(vec2return);
                  },
                  {"MixedEvent"})
          .Define("NumContribCalib2",
                  [](const vector<EventData> &events) {
                    vector<float> vec2return;
                    for (const auto &event : events) {
                      vec2return.push_back(event.event_info2.fNumContrib);
                    }
                    return ROOT::VecOps::RVec<float>(vec2return);
                  },
                  {"MixedEvent"})
          .Define("vector_Eta_Mass_Pt_Sign_DeltaPhi_DeltaEta",
                  [](const vector<EventData> &events) {
                    vector<array<float, 6>> vec2return;
                    for (const auto &event : events) {
                      for (size_t i = 0; i < event.jpsi_info.fPT.size(); ++i) {
                        for (size_t j = i + 1; j < event.jpsi_info.fPT.size();
                             ++j) {
                          if (event.jpsi_info.fSign[i] == 0 &&
                              event.jpsi_info.fSign[j] == 0) {
                            float deltaPhi = event.jpsi_info.fPhi[i] -
                                             event.jpsi_info.fPhi[j];
                            if (deltaPhi > M_PI)
                              deltaPhi -= 2 * M_PI;
                            else if (deltaPhi < -M_PI)
                              deltaPhi += 2 * M_PI;
                            float deltaEta = event.jpsi_info.fEta[i] -
                                             event.jpsi_info.fEta[j];
                            vec2return.push_back({event.jpsi_info.fEta[i],
                                                  event.jpsi_info.fMass[i],
                                                  event.jpsi_info.fPT[i],
                                                  event.jpsi_info.fSign[i],
                                                  deltaPhi, deltaEta});
                          }
                        }
                      }
                    }
                  },
                  {"MixedEvent"});
}

int main(int argc, char **argv) {
  TString path_input_flowVecd = argv[1];
  TString path_output = argv[2];
  EventMixingReading(path_input_flowVecd, path_output);
  return 0;
}