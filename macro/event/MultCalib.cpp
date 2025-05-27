#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include <ROOT/RDataFrame.hxx>

void MultCalib(
    TString path_input = "../input.root", TString path_output = "output.root",
    int runNumber = 0,
    TString path_hist_calib1 =
        "/lustre/alice/users/szhu/work/Analysis/InfoRun/MultCalib/"
        "MultCalibration_LHC22pass4_dqfilter.root:h_NumContribPosZ_calibration",
    TString path_hist_calib2 =
        "/lustre/alice/users/szhu/work/Analysis/InfoRun/MultCalib/"
        "MultCalibration_LHC22pass4_dqfilter.root:fNumContribRun") {
  TFile *file_event = TFile::Open(path_input);
  TFile *fOutput = new TFile(path_output, "RECREATE");

  Calib_NumContrib_fPosZ_Run::GetHistCali(path_hist_calib1, path_hist_calib2,
                                          runNumber);

  TTree *tree_event = (TTree *)file_event->Get("O2reducedevent");
  TTree *tree_event_ext = (TTree *)file_event->Get("O2reextended");

  tree_event->AddFriend(tree_event_ext);
  vector<RResultHandle> gRResultHandlesFast;
  ROOT::RDataFrame rdf(*tree_event);

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
          .Define("RunNumber", [] { return float(0.5); })
          .DefineSlot("NumContribCalib",
                      Calib_NumContrib_fPosZ_Run::NumContribCalibratedFloat,
                      {"fNumContrib", "fPosZ"});
  auto rdf_isntITSROFrameBorder =
      rdf_witTrigger.Filter("isntITSROFrameBorder", "no ITS RO Frame border");
  auto rdf_isntTimeFrameBorder =
      rdf_witTrigger.Filter("isntTimeFrameBorder", "no Time Frame border");
  auto rdf_isTriggerTVX =
      rdf_witTrigger.Filter("isTriggerTVX", "is Trigger TVX");
  auto rdf_fullTrigger =
      rdf_witTrigger.Filter("isTriggerTVX", "is Trigger TVX")
          .Filter("isntITSROFrameBorder", "no ITS RO Frame border")
          .Filter("isntTimeFrameBorder", "no Time Frame border")
          .Filter("isntSameBunchPileup", "no same bunch pileup");
  ROOT::RDF::Experimental::AddProgressBar(rdf_fullTrigger);

  /* #region mult */
  gRResultHandlesFast.push_back(rdf_fullTrigger.Histo1D(
      {"fNumContrib", "fNumContrib;Counts;NumContrib", 300, 0, 300},
      "fNumContrib"));
  gRResultHandlesFast.push_back(rdf_fullTrigger.Histo2D(
      {"fNumContrib_fMultTPC",
       "fNumContrib_fMultTPC;fNumContrib;fMultTPC;Counts", 300, 0, 300, 300, 0,
       300},
      "fNumContrib", "fMultTPC"));
  gRResultHandlesFast.push_back(rdf_fullTrigger.Histo2D(
      {"NumContribCalib_fMultTPC",
       "NumContribCalib_fMultTPC;NumContribCalib;fMultTPC;Counts", 300, 0, 300,
       300, 0, 300},
      "NumContribCalib", "fMultTPC"));
  gRResultHandlesFast.push_back(rdf_fullTrigger.Histo2D(
      {"fMultNTracksPV_fMultTPC",
       "fMultNTracksPV_fMultTPC;fMultNTracksPV;fMultTPC;Counts", 150, 0, 150,
       300, 0, 300},
      "fMultNTracksPV", "fMultTPC"));

  auto profile_fNumContribRun = rdf_fullTrigger.Profile1D(
      {"fNumContribRun", "fNumContribRun;run; fNumContrib", 1, 0., 1.},
      "RunNumber", "fNumContrib");
  gRResultHandlesFast.push_back(profile_fNumContribRun);
  auto profile_fNumContribVtxZ = rdf_fullTrigger.Profile1D(
      {"fNumContribfPosZ", "fNumContribVtxZ;fPosZ [cm]; fNumContrib", 10, -10,
       10.},
      "fPosZ", "fNumContrib");
  gRResultHandlesFast.push_back(profile_fNumContribVtxZ);
  auto profile_NumContribCalib = rdf_fullTrigger.Profile1D(
      {"NumContribCalib", "NumContribCalib;fPosZ [cm]; fNumContrib Calib", 10,
       -10, 10.},
      "fPosZ", "NumContribCalib");
  gRResultHandlesFast.push_back(profile_NumContribCalib);
  auto profile_fNumContribCalibRun = rdf_fullTrigger.Profile1D(
      {"fNumContribCalibRun", "fNumContribCalibRun;run; fNumContrib Calib", 1,
       0., 1.},
      "RunNumber", "NumContribCalib");
  gRResultHandlesFast.push_back(profile_fNumContribCalibRun);
  /* #endregion */
  RunGraphs(gRResultHandlesFast);
  profile_fNumContribRun->GetXaxis()->SetBinLabel(1, Form("%d", runNumber));

  fOutput->cd();
  RResultWrite(gRResultHandlesFast);
  fOutput->Close();
}

int main(int argc, char **argv) {
  TString path_input = "../input.root";
  TString path_output = "output.root";
  TString path_calib_hist1 =
      "/lustre/alice/users/szhu/work/Analysis/InfoRun/MultCalib/"
      "MultCalibration_LHC22pass4_dqfilter.root:h_NumContribPosZ_calibration";
  TString path_calib_hist2 =
      "/lustre/alice/users/szhu/work/Analysis/InfoRun/MultCalib/"
      "MultCalibration_LHC22pass4_dqfilter.root:fNumContribRun";
  int runNumber = 0;

  if (argc > 1) {
    path_input = argv[1];
  }
  if (argc > 2) {
    path_output = argv[2];
  }
  if (argc > 3) {
    runNumber = atoi(argv[3]);
  }
  if (argc > 4) {
    path_calib_hist1 = argv[4];
  }
  if (argc > 5) {
    path_calib_hist2 = argv[5];
  }

  MultCalib(path_input, path_output, runNumber, path_calib_hist1,
            path_calib_hist2);
  return 0;
}