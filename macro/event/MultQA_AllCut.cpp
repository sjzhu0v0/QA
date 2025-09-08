#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include <ROOT/RDataFrame.hxx>

void MultPileupCut(
    TString path_input = "../input.root", TString path_output = "output.root",
    int runNumber = 0,
    TString path_calib =
        "/lustre/alice/users/szhu/work/Analysis/InfoRun/MultCalib/"
        "MultCalibration_LHC22pass4_dqfilter.root:fNumContribfPosZRun_calib_",
    TString path_pileup =
        " /home/szhu/work/alice/analysis/QA/output/event/"
        "MultCalibrationResult_LHC22pass4_dqfilter.root:fit_func_upedge") {
  TFile *file_event = TFile::Open(path_input);
  TFile *fOutput = new TFile(path_output, "RECREATE");

  Calib_NumContrib_fPosZ_Run::GetHistCali(path_calib, runNumber);
  Cut_MultTPC_NumContrib::init(path_pileup);

  // TTree *tree_event = (TTree *)file_event->Get("O2reducedevent");
  // TTree *tree_event_ext = (TTree *)file_event->Get("O2reextended");
  TTree *tree_event = MRootIO::OpenChain(file_event, "O2reducedevent");
  TTree *tree_event_ext = MRootIO::OpenChain(file_event, "O2reextended");

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
                      {"fNumContrib", "fPosZ"})
          .DefineSlot("isntSelfDefinedPileup",
                      Cut_MultTPC_NumContrib::isInCutSlot,
                      {"NumContribCalib", "fMultTPC"})
          .Define(
              "LevelTrigger",
              [](bool bool1, bool bool2, bool bool3, bool bool4, bool bool5) {
                // if only bool1 but bool2 is false, return 1
                ROOT::RVec<int> vec_levelTrigger;
                vec_levelTrigger.push_back(0);
                if (bool1)
                  vec_levelTrigger.push_back(1);
                if (bool1 && bool2)
                  vec_levelTrigger.push_back(2);
                if (bool1 && bool2 && bool3)
                  vec_levelTrigger.push_back(3);
                if (bool1 && bool2 && bool3 && bool4)
                  vec_levelTrigger.push_back(4);
                if (bool1 && bool2 && bool3 && bool4 && bool5)
                  vec_levelTrigger.push_back(5);

                return vec_levelTrigger;
              },
              {"isTriggerTVX", "isntTimeFrameBorder", "isntITSROFrameBorder",
               "isntSameBunchPileup", "isntSelfDefinedPileup"});
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
          .Filter("isntSameBunchPileup", "no same bunch pileup")
          .Filter("isntSelfDefinedPileup", "no self defined pileup");
  auto rdf_basicTrigger =
      rdf_witTrigger.Filter("isTriggerTVX", "is Trigger TVX")
          .Filter("isntITSROFrameBorder", "no ITS RO Frame border")
          .Filter("isntTimeFrameBorder", "no Time Frame border")
          .Filter("isntSameBunchPileup", "no same bunch pileup");
  ROOT::RDF::Experimental::AddProgressBar(rdf_fullTrigger);

  gRResultHandlesFast.push_back(rdf_witTrigger.Histo1D(
      {"LevelTrigger", "Trigger Level; Trigger Level", 6, -0.5, 5.5},
      "LevelTrigger"));

  StrVar4Hist var_fPosX("fPosX", "fPosX", "cm", 200, {-10, 10});
  StrVar4Hist var_fPosY("fPosY", "fPosY", "cm", 200, {-10, 10});
  StrVar4Hist var_fPosZ("fPosZ", "fPosZ", "cm", 200, {-10, 10});
  StrVar4Hist var_fNumContrib("fNumContrib", "fNumContrib", "a.u.", 300,
                              {0, 300});
  StrVar4Hist var_fMultTPC("fMultTPC", "fMultTPC", "a.u.", 600, {0, 600});
  StrVar4Hist var_fMultFV0A("fMultFV0A", "fMultFV0A", "a.u.", 100, {0, 50.e3});
  StrVar4Hist var_fMultFT0A("fMultFT0A", "fMultFT0A", "a.u.", 310,
                            {-1000., 30000.});
  StrVar4Hist var_fMultFT0C("fMultFT0C", "fMultFT0C", "a.u.", 130,
                            {-1000., 12000.});
  StrVar4Hist var_fMultFDDA("fMultFDDA", "fMultFDDA", "a.u.", 360,
                            {-1000, 35000});
  StrVar4Hist var_fMultNTracksPV("fMultNTracksPV", "fMultNTracksPV", "a.u.",
                                 150, {0, 150});
  vector<StrVar4Hist> vec_var_pos = {var_fPosX, var_fPosY, var_fPosZ};
  vector<StrVar4Hist> vec_var_mult = {
      var_fNumContrib, var_fMultFV0A,      var_fMultFT0A, var_fMultFT0C,
      var_fMultFDDA,   var_fMultNTracksPV, var_fMultTPC};

#define localDefineHis1D(fVar, fTag)                                           \
  gRResultHandlesFast.push_back(                                               \
      rdf_##fTag.Histo1D(var_##fVar.GetTH1DModel(#fTag), var_##fVar.fName));

  localDefineHis1D(fPosX, basicTrigger);
  localDefineHis1D(fPosY, basicTrigger);
  localDefineHis1D(fPosZ, basicTrigger);
  localDefineHis1D(fNumContrib, basicTrigger);
  localDefineHis1D(fMultFV0A, basicTrigger);
  localDefineHis1D(fMultFT0A, basicTrigger);
  localDefineHis1D(fMultFT0C, basicTrigger);
  localDefineHis1D(fMultFDDA, basicTrigger);
  localDefineHis1D(fMultNTracksPV, basicTrigger);
  localDefineHis1D(fMultTPC, basicTrigger);

  localDefineHis1D(fPosX, fullTrigger);
  localDefineHis1D(fPosY, fullTrigger);
  localDefineHis1D(fPosZ, fullTrigger);
  localDefineHis1D(fNumContrib, fullTrigger);
  localDefineHis1D(fMultTPC, fullTrigger);
  localDefineHis1D(fMultFV0A, fullTrigger);
  localDefineHis1D(fMultFT0A, fullTrigger);
  localDefineHis1D(fMultFT0C, fullTrigger);
  localDefineHis1D(fMultFDDA, fullTrigger);
  localDefineHis1D(fMultNTracksPV, fullTrigger);

  for (auto &var_pos : vec_var_pos) {
    for (auto &var_mult : vec_var_mult) {
      gRResultHandlesFast.push_back(rdf_fullTrigger.Histo2D(
          GetTH2DModel(var_pos, var_mult, "fullTrigger"), var_pos.fName,
          var_mult.fName));
      gRResultHandlesFast.push_back(rdf_basicTrigger.Histo2D(
          GetTH2DModel(var_pos, var_mult, "basicTrigger"), var_pos.fName,
          var_mult.fName));
    }
  }

  for (int i_mult1 = 0; i_mult1 < vec_var_mult.size(); i_mult1++) {
    for (int i_mult2 = i_mult1 + 1; i_mult2 < vec_var_mult.size(); i_mult2++) {
      gRResultHandlesFast.push_back(rdf_fullTrigger.Histo2D(
          GetTH2DModel(vec_var_mult[i_mult1], vec_var_mult[i_mult2],
                       "fullTrigger"),
          vec_var_mult[i_mult1].fName, vec_var_mult[i_mult2].fName));
      gRResultHandlesFast.push_back(rdf_basicTrigger.Histo2D(
          GetTH2DModel(vec_var_mult[i_mult1], vec_var_mult[i_mult2],
                       "basicTrigger"),
          vec_var_mult[i_mult1].fName, vec_var_mult[i_mult2].fName));
    }
  }
  RunGraphs(gRResultHandlesFast);

  fOutput->cd();
  RResultWrite(gRResultHandlesFast);
  fOutput->Close();
}

int main(int argc, char **argv) {
  TString path_input = "../input.root";
  TString path_output = "output.root";
  int runNumber = 0;
  TString path_calib =
      "/lustre/alice/users/szhu/work/Analysis/InfoRun/MultCalib/"
      "MultCalibration_LHC22pass4_dqfilter.root:fNumContribfPosZRun_calib_";
  TString path_pileup =
      "/lustre/alice/users/szhu/work/Analysis/InfoRun/MultCalib/"
      "MultPileup_LHC22pass4_dqfilter.root:fit_func_upedge";

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
    path_calib = argv[4];
  }
  if (argc > 5) {
    path_pileup = argv[5];
  }

  MultPileupCut(path_input, path_output, runNumber, path_calib, path_pileup);
  return 0;
}