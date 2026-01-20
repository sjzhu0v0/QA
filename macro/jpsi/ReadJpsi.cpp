#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "opt/EventData.h"
#include <ROOT/RDataFrame.hxx>

void JpsiQA(
    TString path_input_flowVecd = "../input.root",
    TString path_output = "output.root",
    TString path_output_tree = "output_tree.root", int runNumber = 0,
    TString path_calib =
        "/lustre/alice/users/szhu/work/Analysis/InfoRun/MultCalib/"
        "MultCalibration_LHC22pass4_dqfilter.root:fNumContribfPosZRun_calib_"/* ,
    TString path_pileup =
        " /home/szhu/work/alice/analysis/QA/output/event/"
        "MultCalibrationResult_LHC22pass4_dqfilter.root:fit_func_upedge" */) {
  ROOT::DisableImplicitMT();
  TFile* file_flowVecd = TFile::Open(path_input_flowVecd);
  TFile* fOutput = new TFile(path_output, "RECREATE");

  Calib_NumContrib_fPosZ_Run::GetHistCali(path_calib, runNumber);
  // Cut_MultTPC_NumContrib::init(path_pileup);

  TChain* tree_flowVecd = MRootIO::OpenChain(file_flowVecd, "O2dqflowvecd");

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
          .DefineSlot("NumContribCalib", Calib_NumContrib_fPosZ_Run::NumContribCalibratedFloat,
                      {"fNumContrib", "fPosZ"})
      /*   .DefineSlot("isntSelfDefinedPileup",
                    Cut_MultTPC_NumContrib::isInCutSlot,
                    {"NumContribCalib", "fMultTPC"}) */
      ;
  rdf_witTrigger.Snapshot("MultCalib", path_output_tree, {"NumContribCalib"});
  auto rdf_isntITSROFrameBorder =
      rdf_witTrigger.Filter("isntITSROFrameBorder", "no ITS RO Frame border");
  auto rdf_isntTimeFrameBorder =
      rdf_witTrigger.Filter("isntTimeFrameBorder", "no Time Frame border");
  auto rdf_isTriggerTVX = rdf_witTrigger.Filter("isTriggerTVX", "is Trigger TVX");
  auto rdf_PartTrigger = rdf_witTrigger.Filter("isTriggerTVX", "is Trigger TVX")
                             .Filter("isntITSROFrameBorder", "no ITS RO Frame border")
                             .Filter("isntTimeFrameBorder", "no Time Frame border")
                             .Filter("isntSameBunchPileup", "no Time Frame border")
      /*  .Filter("isntSelfDefinedPileup", "no self defined pileup") */;

  StrVar4Hist var_fPosZSel("fPosZ", "#it{V}_{Z}", "cm", 200, {-10, 10});
  StrVar4Hist var_NumContribCalibBinned("NumContribCalib", "N_{vtx contrib} Calibrated", "", 10,
                                        {0, 5, 8, 11, 14, 18, 23, 28, 36, 48, 300});
  StrVar4Hist var_MassJpsiCandidate("sel_Mass", "M_{ee}", "GeV^{2}/c^{4}", 90, {1.8, 5.4});
  StrVar4Hist var_PtJpsiCandidate("sel_Pt", "p_{T}", "GeV/c", 10, {0., 5.});

#define obj2push_thnd(rdf2push, ...)                                                               \
  do {                                                                                             \
    TupleTHnDModel tuple_thnd = GetTHnDModelWithTitle(__VA_ARGS__);                                \
    gRResultHandles.push_back(rdf2push.HistoND(get<0>(tuple_thnd), get<1>(tuple_thnd)));           \
  } while (0)

#define pushCut(name, ...)                                                                         \
  rdf_PartTrigger##name = rdf_PartTrigger.Define("sel" #name, __VA_ARGS__)                         \
                              .Define("sel_Mass", "fMass[sel" #name "]")                           \
                              .Define("sel_Pt", "fPT[sel" #name "]");                              \
  obj2push_thnd(                                                                                   \
      rdf_PartTrigger##name,                                                                       \
      {var_fPosZ, var_MassJpsiCandidate, var_PtJpsiCandidate, var_NumContribCalibBinned}, #name,   \
      #name);

  pushCut(default,
          "fTPCNSigmaPi1>3&&fTPCNsigmaPi2>3&&fTPCNSigmaPr1>3&&fTPCNSigmaPr2>3&&fTPCNSigmaEl1>"
          "-2&&fTPCNSigmaEl2>-2&&fTPCNSigmaEl1<3&&fTPCNSigmaEl2<3");
  pushCut(pid1, "fTPCNSigmaPi1>2.7&&fTPCNsigmaPi2>2.7&&fTPCNSigmaPr1>2.7&&fTPCNSigmaPr2>2.7&&"
                "fTPCNSigmaEl1>-3&&fTPCNSigmaEl2>-3&&fTPCNSigmaEl1<3&&fTPCNSigmaEl2<3");
  pushCut(pid2, "fTPCNSigmaPi1>2.5&&fTPCNsigmaPi2>2.5&&fTPCNSigmaPr1>2.5&&fTPCNSigmaPr2>2.5&&"
                "fTPCNSigmaEl1>-3&&fTPCNSigmaEl2>-3&&fTPCNSigmaEl1<3&&fTPCNSigmaEl2<3");
  pushCut(pid3, "fTPCNSigmaPi1>3.5&&fTPCNsigmaPi2>3.5&&fTPCNSigmaPr1>3.5&&fTPCNSigmaPr2>3.5&&"
                "fTPCNSigmaEl1>-2&&fTPCNSigmaEl2>-2&&fTPCNSigmaEl1<3&&fTPCNSigmaEl2<3");
  pushCut(vz1,
          "fTPCNSigmaPi1>3&&fTPCNsigmaPi2>3&&fTPCNSigmaPr1>3&&fTPCNSigmaPr2>3&&fTPCNSigmaEl1>-2&&"
          "fTPCNSigmaEl2>-2&&fTPCNSigmaEl1<3&&fTPCNSigmaEl2<3&&abs(fPosZ)<8");
  pushCut(vz2,
          "fTPCNSigmaPi1>3&&fTPCNsigmaPi2>3&&fTPCNSigmaPr1>3&&fTPCNSigmaPr2>3&&fTPCNSigmaEl1>-2&&"
          "fTPCNSigmaEl2>-2&&fTPCNSigmaEl1<3&&fTPCNSigmaEl2<3&&abs(fPosZ)<6");

  RunGraphs(gRResultHandles);

  fOutput->cd();
  RResultWrite(gRResultHandles);
  fOutput->Close();
}

int main(int argc, char** argv) {
  TString path_input_flowVecd = "../input.root";
  TString path_output = "output.root";
  TString path_output_tree = "output_tree.root";
  int runNumber = 0;
  TString path_calib = "/lustre/alice/users/szhu/work/Analysis/InfoRun/MultCalib/"
                       "MultCalibration_LHC22pass4_dqfilter.root:fNumContribfPosZRun_calib_";
  TString path_pileup = "/lustre/alice/users/szhu/work/Analysis/InfoRun/MultCalib/"
                        "MultPileup_LHC22pass4_dqfilter.root:fit_func_upedge";

  if (argc > 1) {
    path_input_flowVecd = argv[1];
  }
  if (argc > 2) {
    path_output = argv[2];
  }
  if (argc > 3) {
    path_output_tree = argv[3];
  }
  if (argc > 4) {
    runNumber = atoi(argv[4]);
  }
  if (argc > 5) {
    path_calib = argv[5];
  }
  if (argc > 6) {
    path_pileup = argv[6];
  }
  gROOT->SetBatch(kTRUE); // Disable interactive graphics
  JpsiQA(path_input_flowVecd, path_output, path_output_tree, runNumber,
         path_calib /* , path_pileup */);

  return 0;
}
