#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "opt/EventData.h"
#include "yaml-cpp/yaml.h"
#include <ROOT/RDataFrame.hxx>

void JpsiAsso(TString path_input_flowVecd = "../input.root",
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
/* 
*Br    0 :fMultTPC  : fMultTPC/I                                             *
*Br    1 :fMultTracklets : fMultTracklets/I                                  *
*Br    2 :fMultNTracksPV : fMultNTracksPV/I                                  *
*Br    3 :fMultFT0C : fMultFT0C/F                                            *
*Br    4 :fNumContrib : fNumContrib/s                                        *
*Br    5 :fPosX     : fPosX/F                                                *
*Br    6 :fPosY     : fPosY/F                                                *
*Br    7 :fPosZ     : fPosZ/F                                                *
*Br    8 :fSelection : fSelection/l                                          *
*Br    9 :fHadronicRate : fHadronicRate/F                                    *
*Br   10 :fPT_size  : fPT_size/I                                             *
*Br   11 :fPT       : fPT[fPT_size]/F                                        *
*Br   12 :fEta_size : fEta_size/I                                            *
*Br   13 :fEta      : fEta[fEta_size]/F                                      *
*Br   14 :fPhi_size : fPhi_size/I                                            *
*Br   15 :fPhi      : fPhi[fPhi_size]/F                                      *
*Br   16 :fMass_size : fMass_size/I                                          *
*Br   17 :fMass     : fMass[fMass_size]/F                                    *
*Br   18 :fSign_size : fSign_size/I                                          *
*Br   19 :fSign     : fSign[fSign_size]/F                                    *
*Br   20 :fPTREF_size : fPTREF_size/I                                        *
*Br   21 :fPTREF    : fPTREF[fPTREF_size]/F                                  *
*Br   22 :fEtaREF_size : fEtaREF_size/I                                      *
*Br   23 :fEtaREF   : fEtaREF[fEtaREF_size]/F                                *
*Br   24 :fPhiREF_size : fPhiREF_size/I                                      *
*Br   25 :fPhiREF   : fPhiREF[fPhiREF_size]/F                                *
*Br   26 :fITSChi2NCl_size : fITSChi2NCl_size/I                              *
*Br   27 :fITSChi2NCl : fITSChi2NCl[fITSChi2NCl_size]/F                      *
*Br   28 :fTPCNClsCR_size : fTPCNClsCR_size/I                                *
*Br   29 :fTPCNClsCR : fTPCNClsCR[fTPCNClsCR_size]/F                         *
*Br   30 :fTPCNClsFound_size : fTPCNClsFound_size/I                          *
*Br   31 :fTPCNClsFound : fTPCNClsFound[fTPCNClsFound_size]/F                *
*Br   32 :fTPCChi2NCl_size : fTPCChi2NCl_size/I                              *
*Br   33 :fTPCChi2NCl : fTPCChi2NCl[fTPCChi2NCl_size]/F                      *
*Br   34 :fTPCSignal_size : fTPCSignal_size/I                                *
*Br   35 :fTPCSignal : fTPCSignal[fTPCSignal_size]/F                         *
*Br   36 :fTPCNSigmaEl_size : fTPCNSigmaEl_size/I                            *
*Br   37 :fTPCNSigmaEl : fTPCNSigmaEl[fTPCNSigmaEl_size]/F                   *
*Br   38 :fTPCNSigmaPi_size : fTPCNSigmaPi_size/I                            *
*Br   39 :fTPCNSigmaPi : fTPCNSigmaPi[fTPCNSigmaPi_size]/F                   *
*Br   40 :fTPCNSigmaPr_size : fTPCNSigmaPr_size/I                            *
*Br   41 :fTPCNSigmaPr : fTPCNSigmaPr[fTPCNSigmaPr_size]/F                   *
*Br   42 :fPt1_size : fPt1_size/I                                            *
*Br   43 :fPt1      : fPt1[fPt1_size]/F                                      *
*Br   44 :fEta1_size : fEta1_size/I                                          *
*Br   45 :fEta1     : fEta1[fEta1_size]/F                                    *
*Br   46 :fPhi1_size : fPhi1_size/I                                          *
*Br   47 :fPhi1     : fPhi1[fPhi1_size]/F                                    *
*Br   48 :fSign1_size : fSign1_size/I                                        *
*Br   49 :fSign1    : fSign1[fSign1_size]/I                                  *
*Br   50 :fITSChi2NCl1_size : fITSChi2NCl1_size/I                            *
*Br   51 :fITSChi2NCl1 : fITSChi2NCl1[fITSChi2NCl1_size]/F                   *
*Br   52 :fTPCNClsCR1_size : fTPCNClsCR1_size/I                              *
*Br   53 :fTPCNClsCR1 : fTPCNClsCR1[fTPCNClsCR1_size]/F                      *
*Br   54 :fTPCNClsFound1_size : fTPCNClsFound1_size/I                        *
*Br   55 :fTPCNClsFound1 : fTPCNClsFound1[fTPCNClsFound1_size]/F             *
*Br   56 :fTPCChi2NCl1_size : fTPCChi2NCl1_size/I                            *
*Br   57 :fTPCChi2NCl1 : fTPCChi2NCl1[fTPCChi2NCl1_size]/F                   *
*Br   58 :fTPCSignal1_size : fTPCSignal1_size/I                              *
*Br   59 :fTPCSignal1 : fTPCSignal1[fTPCSignal1_size]/F                      *
*Br   60 :fTPCNSigmaEl1_size : fTPCNSigmaEl1_size/I                          *
*Br   61 :fTPCNSigmaEl1 : fTPCNSigmaEl1[fTPCNSigmaEl1_size]/F                *
*Br   62 :fTPCNSigmaPi1_size : fTPCNSigmaPi1_size/I                          *
*Br   63 :fTPCNSigmaPi1 : fTPCNSigmaPi1[fTPCNSigmaPi1_size]/F                *
*Br   64 :fTPCNSigmaPr1_size : fTPCNSigmaPr1_size/I                          *
*Br   65 :fTPCNSigmaPr1 : fTPCNSigmaPr1[fTPCNSigmaPr1_size]/F                *
*Br   66 :fPt2_size : fPt2_size/I                                            *
*Br   67 :fPt2      : fPt2[fPt2_size]/F                                      *
*Br   68 :fEta2_size : fEta2_size/I                                          *
*Br   69 :fEta2     : fEta2[fEta2_size]/F                                    *
*Br   70 :fPhi2_size : fPhi2_size/I                                          *
*Br   71 :fPhi2     : fPhi2[fPhi2_size]/F                                    *
*Br   72 :fSign2_size : fSign2_size/I                                        *
*Br   73 :fSign2    : fSign2[fSign2_size]/I                                  *
*Br   74 :fITSChi2NCl2_size : fITSChi2NCl2_size/I                            *
*Br   75 :fITSChi2NCl2 : fITSChi2NCl2[fITSChi2NCl2_size]/F                   *
*Br   76 :fTPCNClsCR2_size : fTPCNClsCR2_size/I                              *
*Br   77 :fTPCNClsCR2 : fTPCNClsCR2[fTPCNClsCR2_size]/F                      *
*Br   78 :fTPCNClsFound2_size : fTPCNClsFound2_size/I                        *
*Br   79 :fTPCNClsFound2 : fTPCNClsFound2[fTPCNClsFound2_size]/F             *
*Br   80 :fTPCChi2NCl2_size : fTPCChi2NCl2_size/I                            *
*Br   81 :fTPCChi2NCl2 : fTPCChi2NCl2[fTPCChi2NCl2_size]/F                   *
*Br   82 :fTPCSignal2_size : fTPCSignal2_size/I                              *
*Br   83 :fTPCSignal2 : fTPCSignal2[fTPCSignal2_size]/F                      *
*Br   84 :fTPCNSigmaEl2_size : fTPCNSigmaEl2_size/I                          *
*Br   85 :fTPCNSigmaEl2 : fTPCNSigmaEl2[fTPCNSigmaEl2_size]/F                *
*Br   86 :fTPCNSigmaPi2_size : fTPCNSigmaPi2_size/I                          *
*Br   87 :fTPCNSigmaPi2 : fTPCNSigmaPi2[fTPCNSigmaPi2_size]/F                *
*Br   88 :fTPCNSigmaPr2_size : fTPCNSigmaPr2_size/I                          *
*Br   89 :fTPCNSigmaPr2 : fTPCNSigmaPr2[fTPCNSigmaPr2_size]/F                *
*/
  auto rdf_PartTriggerWithJpsiWithEvent =
      rdf_PartTrigger
          .Define("EventData", CreateEventData,
                  {"fMultTPC", "fMultTracklets", "fMultNTracksPV", "fMultFT0C",
                   "fNumContrib", "NumContribCalib", "fPosX", "fPosY", "fPosZ",
                   "fSelection", "fHadronicRate", "fPT", "fEta", "fPhi",
                   "fMass", "fSign", "fPTREF", "fEtaREF", "fPhiREF"})
          .Define(
              "JpsiInfoUS",
              [&CutTrackInfo, &low_edge_deltaPhiToPi,
               &up_edge_deltaPhiToPi](const EventData &eventData) {
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
      {0, 7, 12, 17, 22, 29, 37, 46, 57, 73, 300});
  int n_bins_mass_assoYield =
      config["hist_binning"]["n_bins_mass_assoYield"].as<int>();
  StrVar4Hist var_MassJpsiCandidate("MassUS", "M_{ee}", "GeV^{2}/c^{4}",
                                    n_bins_mass_assoYield, {1.8, 5.4});
  StrVar4Hist var_PtJpsiCandidate("PtUS", "p_{T}", "GeV/c", 10, {0., 5.});
  int n_bins_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["n_bins"].as<int>();
  double min_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["min"].as<double>();
  double max_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["max"].as<double>();
  StrVar4Hist var_DeltaEtaUS("DeltaEtaUS", "#Delta#eta_{J/#psi, track}", "",
                             n_bins_deltaEta_assoYield,
                             {min_deltaEta_assoYield, max_deltaEta_assoYield});
  int n_bins_deltaPhi_assoYield =
      config["hist_binning"]["n_bins_deltaPhi_assoYield"].as<int>();
  StrVar4Hist var_DeltaPhiUS(
      "DeltaPhiUS", "#Delta#phi_{J/#psi, track}", "", n_bins_deltaPhi_assoYield,
      {low_edge_deltaPhiToPi * M_PI, up_edge_deltaPhiToPi * M_PI});
  StrVar4Hist var_fPosZSingle("PosZUSSingle", "#it{V}_{Z}", "cm", 8, {-10, 10});
  StrVar4Hist var_NumContribCalibBinnedSingle(
      "NumContribCalibUSSingle", "N_{vtx contrib} Calibrated", "", 10,
      {0, 7, 12, 17, 22, 29, 37, 46, 57, 73, 300});
  StrVar4Hist var_MassJpsiCandidateSingle("MassUSSingle", "M_{ee}",
                                          "GeV^{2}/c^{4}", 90, {1.8, 5.4});
  StrVar4Hist var_PtJpsiCandidateSingle("PtUSSingle", "p_{T}", "GeV/c", 10,
                                        {0., 5.});

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

  JpsiAsso(path_input_flowVecd, path_input_mult, path_output, threshold_bs);

  return 0;
}
