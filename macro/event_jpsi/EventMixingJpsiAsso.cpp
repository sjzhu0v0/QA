#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "MixJpsiAsso.h"
#include <ROOT/RDataFrame.hxx>

void EventMixingJpsiAsso(
    TString path_input_flowVecd = "../input.root",
    TString path_output = "output.root", int runNumber = 0,
    TString path_calib =
        "/lustre/alice/users/szhu/work/Analysis/InfoRun/MultCalib/"
        "MultCalibration_LHC22pass4_dqfilter.root:fNumContribfPosZRun_calib_",
    TString path_pileup =
        " /home/szhu/work/alice/analysis/QA/output/event/"
        "MultCalibrationResult_LHC22pass4_dqfilter.root:fit_func_upedge") {
  TFile *file_flowVecd = TFile::Open(path_input_flowVecd);
  TFile *fOutput = new TFile(path_output, "RECREATE");

  Calib_NumContrib_fPosZ_Run::GetHistCali(path_calib, runNumber);
  Cut_MultTPC_NumContrib::init(path_pileup);

  TChain *tree_flowVecd = MRootIO::OpenChain(file_flowVecd, "O2dqflowvecd");

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
          .DefineSlot("NumContribCalib",
                      Calib_NumContrib_fPosZ_Run::NumContribCalibratedFloat,
                      {"fNumContrib", "fPosZ"})
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
      rdf_witTrigger.Filter("fEta_size>=1", "has Jpsi");

  auto rdf_PartTriggerWithJpsiWithEvent = rdf_PartTriggerWithJpsi.Define(
      "EventData", CreateEventData,
      {"fMultTPC", "fMultTracklets", "fMultNTracksPV", "fMultFT0C",
       "fNumContrib", "fPosX", "fPosY", "fPosZ", "fSelection", "fHadronicRate",
       "fPT", "fEta", "fPhi", "fMass", "fSign", "fPTREF", "fEtaREF",
       "fPhiREF"});

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

  const vector<double> bins_mix_numContrib = var_fPosZMix.fBins;
  const vector<double> bins_mix_posZ = var_fPosZMix.fBins;

  auto rdf_PartTriggerWithJpsiWithEventWithEventMixing =
      rdf_PartTriggerWithJpsiWithEvent
          .Define("IndexMixing",
                  [bins_mix_numContrib, bins_mix_posZ](double numContrib,
                                                       float posZ) {
                    int index1 = -1;
                    for (int i = 0; i < bins_mix_numContrib.size() - 1; i++) {
                      if (numContrib >= bins_mix_numContrib[i] &&
                          numContrib < bins_mix_numContrib[i + 1]) {
                        index1 = i;
                        break;
                      }
                    }
                    int index2 = -1;
                    for (int i = 0; i < bins_mix_posZ.size() - 1; i++) {
                      if (posZ >= bins_mix_posZ[i] &&
                          posZ < bins_mix_posZ[i + 1]) {
                        index2 = i;
                        break;
                      }
                    }
                    if (index1 == -1 || index2 == -1) {
                      return -1; // Invalid index
                    }
                    return int(index1 * bins_mix_posZ.size() + index2);
                  },
                  {"NumContribCalib", "fPosZ"})
          .Define("MixedEvent", MixEvent, {"IndexMixing", "EventData"})
          .Snapshot("EventMixing", "EventMixingJpsiAsso.root", {"MixedEvent"});

#define obj2push_thnd(rdf2push, ...)                                           \
  do {                                                                         \
    TupleTHnDModel tuple_thnd = GetTHnDModelWithTitle(__VA_ARGS__);            \
    gRResultHandles.push_back(                                                 \
        rdf2push.HistoND(get<0>(tuple_thnd), get<1>(tuple_thnd)));             \
  } while (0)

  obj2push_thnd(rdf_PartTrigger, {var_fPosZ, var_MassJpsiCandidate,
                                  var_PtJpsiCandidate, var_NumContribCalib});
  obj2push_thnd(rdf_PartTrigger, {var_fPosZ});
  obj2push_thnd(rdf_PartTrigger, {var_MassJpsiCandidate});
  obj2push_thnd(rdf_PartTrigger, {var_PtJpsiCandidate});
  obj2push_thnd(rdf_PartTrigger, {var_NumContribCalib}, "title_test");
  obj2push_thnd(rdf_PartTrigger,
                {var_fPosZ, var_MassJpsiCandidate, var_PtJpsiCandidate,
                 var_NumContribCalibBinned},
                "", "Binned");

  fOutput->cd();
  RResultWrite(gRResultHandles);
  fOutput->Close();
}

int main(int argc, char **argv) {
  TString path_input_flowVecd = "../input.root";
  TString path_output = "output.root";
  int runNumber = 0;
  TString path_calib =
      "/lustre/alice/users/szhu/work/Analysis/InfoRun/MultCalib/"
      "MultCalibration_LHC22pass4_dqfilter.root:fNumContribfPosZRun_calib_";
  TString path_pileup =
      "/lustre/alice/users/szhu/work/Analysis/InfoRun/MultCalib/"
      "MultPileup_LHC22pass4_dqfilter.root:fit_func_upedge";

  if (argc > 1) {
    path_input_flowVecd = argv[1];
  }
  if (argc > 2) {
    path_output = argv[2];
  }
  if (argc > 3) {
    runNumber = atoi(argv[3]);
  }
  if (argc > 4) {
    path_calib = argv[4];
  }
  if (argc > 5) {
    path_pileup = argv[5];
  }
  EventMixingJpsiAsso(path_input_flowVecd, path_output, runNumber, path_calib,
                      path_pileup);

  return 0;
}