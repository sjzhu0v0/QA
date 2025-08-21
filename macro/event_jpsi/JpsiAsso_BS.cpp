#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "opt/EventData.h"
#include <ROOT/RDataFrame.hxx>

void JpsiAsso(
    TString path_input_flowVecd = "../input.root", TString path_input_mult = "../input2.root",
    TString path_output = "output.root", /* int runNumber = 0, */
    double threshold_bs = 1.
  /*   ,TString path_calib =
        "/lustre/alice/users/szhu/work/Analysis/InfoRun/MultCalib/"
        "MultCalibration_LHC22pass4_dqfilter.root:fNumContribfPosZRun_calib_",
    TString path_pileup =
        " /home/szhu/work/alice/analysis/QA/output/event/"
        "MultCalibrationResult_LHC22pass4_dqfilter.root:fit_func_upedge" */) {
  TFile *file_flowVecd = TFile::Open(path_input_flowVecd);
  TFile *file_mult = TFile::Open(path_input_mult);
  TFile *fOutput = new TFile(path_output, "RECREATE");

  cout << "Input file: " << path_input_flowVecd << endl;
  cout << "Input file: " << path_input_mult << endl;
  cout << "Output file: " << path_output << endl;
  cout << "Threshold BS: " << threshold_bs << endl;

  // Calib_NumContrib_fPosZ_Run::GetHistCali(path_calib, runNumber);
  // Cut_MultTPC_NumContrib::init(path_pileup);

  TChain *tree_flowVecd = MRootIO::OpenChain(file_flowVecd, "O2dqflowvecd");
  TChain *tree_mult = MRootIO::OpenChain(file_mult, "MultCalib");

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
          .Define(
              "JpsiInfoUS",
              [&CutTrackInfo](const EventData &eventData) {
                ROOT::VecOps::RVec<array<float, 6>> vec2return;
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
              {"EventData"})
          .Define(
              "JpsiInfoUSSingle",
              [](const EventData &eventData) {
                ROOT::VecOps::RVec<array<float, 4>> vec2return;
                for (size_t i = 0; i < eventData.jpsi_info.fPT.size(); ++i) {
                  if (eventData.jpsi_info.fSign[i] == 0) {
                    vec2return.push_back({eventData.event_info.fPosZ,
                                          eventData.event_info.fNumContribCalib,
                                          eventData.jpsi_info.fMass[i],
                                          eventData.jpsi_info.fPT[i]});
                  }
                }
                return vec2return;
              },
              {"EventData"})
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
                  {"JpsiInfoUS"})
          .Define("PosZUSSingle",
                  [](const ROOT::VecOps::RVec<array<float, 4>> &jpsiInfoUS) {
                    ROOT::VecOps::RVec<float> posZUSSingle;
                    for (const auto &pair : jpsiInfoUS) {
                      posZUSSingle.push_back(pair[0]);
                    }
                    return posZUSSingle;
                  },
                  {"JpsiInfoUSSingle"})
          .Define("NumContribCalibUSSingle",
                  [](const ROOT::VecOps::RVec<array<float, 4>> &jpsiInfoUS) {
                    ROOT::VecOps::RVec<float> numContribCalibUSSingle;
                    for (const auto &pair : jpsiInfoUS) {
                      numContribCalibUSSingle.push_back(pair[1]);
                    }
                    return numContribCalibUSSingle;
                  },
                  {"JpsiInfoUSSingle"})
          .Define("MassUSSingle",
                  [](const ROOT::VecOps::RVec<array<float, 4>> &jpsiInfoUS) {
                    ROOT::VecOps::RVec<float> massUSSingle;
                    for (const auto &pair : jpsiInfoUS) {
                      massUSSingle.push_back(pair[2]);
                    }
                    return massUSSingle;
                  },
                  {"JpsiInfoUSSingle"})
          .Define("PtUSSingle",
                  [](const ROOT::VecOps::RVec<array<float, 4>> &jpsiInfoUS) {
                    ROOT::VecOps::RVec<float> ptUSSingle;
                    for (const auto &pair : jpsiInfoUS) {
                      ptUSSingle.push_back(pair[3]);
                    }
                    return ptUSSingle;
                  },
                  {"JpsiInfoUSSingle"})
          .DefineSlot("Cut_BS",
                      [threshold_bs](unsigned int) {
                        thread_local TRandom3 random_gen(0);
                        return random_gen.Uniform(0, 1) < threshold_bs;
                      })
          .Filter("Cut_BS", "Bootstrap cut for Jpsi association");

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
  StrVar4Hist var_fPosZSingle("PosZUSSingle", "#it{V}_{Z}", "cm", 8, {-10, 10});
  StrVar4Hist var_NumContribCalibBinnedSingle(
      "NumContribCalibUSSingle", "N_{vtx contrib} Calibrated", "", 10,
      {0, 23, 31, 37, 43, 48, 54, 61, 69, 81, 297});
  StrVar4Hist var_MassJpsiCandidateSingle("MassUSSingle", "M_{ee}",
                                          "GeV^{2}/c^{4}", 90, {1.8, 5.4});
  StrVar4Hist var_PtJpsiCandidateSingle("PtUSSingle", "p_{T}", "GeV/c", 10,
                                        {0., 10.});

#define obj2push_thnd(rdf2push, ...)                                           \
  do {                                                                         \
    TupleTHnDModel tuple_thnd = GetTHnDModelWithTitle(__VA_ARGS__);            \
    gRResultHandles.push_back(                                                 \
        rdf2push.HistoND(get<0>(tuple_thnd), get<1>(tuple_thnd)));             \
  } while (0)
  // cout << "Start pushing THnD objects..." << endl;
  obj2push_thnd(rdf_PartTriggerWithJpsiWithEvent,
                {var_DeltaEtaUS, var_DeltaPhiUS, var_fPosZ,
                 var_MassJpsiCandidate, var_PtJpsiCandidate,
                 var_NumContribCalibBinned});
  obj2push_thnd(rdf_PartTriggerWithJpsiWithEvent,
                {var_fPosZSingle, var_MassJpsiCandidateSingle,
                 var_PtJpsiCandidateSingle, var_NumContribCalibBinnedSingle});

  RunGraphs(gRResultHandles);
  // Long64_t totalSize = gRResultHandles[0].GetPtr<THnD>()->GetNbins() *
  //                      sizeof(Double_t); // Approximate size
  // std::cout << "Approx. memory used: " << totalSize << " bytes ("
  //           << totalSize / (1024. * 1024.) << " MB)" << std::endl;
  // cout << "Start writing THnD objects..." << endl;

  // MHnTool hnTool(gRResultHandles[0].GetPtr<THnD>());
  // hnTool.PrintAllAxis();

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

  JpsiAsso(path_input_flowVecd, path_output, threshold_bs);

  return 0;
}
