#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include <ROOT/RDataFrame.hxx>

void MultREF(
    TString path_input_flowVecd = "../input.root",
    TString path_output = "output.root", int runNumber = 0,
    TString path_calib =
        "/lustre/alice/users/szhu/work/Analysis/InfoRun/MultCalib/"
        "MultCalibration_LHC22pass4_dqfilter.root:fNumContribfPosZRun_calib_",
   /*  TString path_pileup =
        " /home/szhu/work/alice/analysis/QA/output/event/"
        "MultCalibrationResult_LHC22pass4_dqfilter.root:fit_func_upedge" */) {
  TFile *file_flowVecd = TFile::Open(path_input_flowVecd);
  TFile *fOutput = new TFile(path_output, "RECREATE");

  Calib_NumContrib_fPosZ_Run::GetHistCali(path_calib, runNumber);
  // Cut_MultTPC_NumContrib::init(path_pileup);

  // TTree *tree_flowVecd = (TTree *)file_flowVecd->Get("O2dqflowvecd");
  TTree *tree_flowVecd = MRootIO::OpenChain(file_flowVecd, "O2dqflowvecd");

  vector<RResultHandle> gRResultHandlesFast;
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
      /*  .Filter("isntSelfDefinedPileup", "no self defined pileup") */;

  StrVar4Hist var_fPosX("fPosX", "fPosX", "cm", 200, {-10, 10});
  StrVar4Hist var_fPosY("fPosY", "fPosY", "cm", 200, {-10, 10});
  StrVar4Hist var_fPosZ("fPosZ", "fPosZ", "cm", 200, {-10, 10});
  StrVar4Hist var_fNumContrib("fNumContrib", "fNumContrib", "", 300, {0, 300});
  StrVar4Hist var_NumContribCalib("NumContribCalib", "NumContrib Calibrated",
                                  "", 300, {0, 300});
  StrVar4Hist var_fMultTPC("fMultTPC", "fMultTPC", "", 600, {0, 600});
  StrVar4Hist var_fMultREF("fMultREF", "fMultREF", "", 100, {0, 100});
  StrVar4Hist var_fMultFT0C("fMultFT0C", "fMultFT0C", "", 130,
                            {-1000., 12000.});
  StrVar4Hist var_fMultNTracksPV("fMultNTracksPV", "fMultNTracksPV", "a.u.",
                                 150, {0, 150});

  vector<StrVar4Hist> vec_var_mult = {var_fNumContrib,    var_NumContribCalib,
                                      var_fMultTPC,       var_fMultFT0C,
                                      var_fMultNTracksPV, var_fMultREF};

  vector<array<string, 2>> conditions = {
      {"isntSameBunchPileup || !isntSameBunchPileup", "NoSameBunchCut"},
      {"isntSameBunchPileup", "NoSameBunchPileup"},
      {"!isntSameBunchPileup", "AllSameBunchPileup"}};

  for (const auto &condition : conditions) {
    auto rdf_PartTrigger_cond =
        rdf_PartTrigger.Filter(condition[0].c_str(), condition[1].c_str());
    TString tag_cond = condition[1];
    for (int i_mult = 0; i_mult < vec_var_mult.size(); i_mult++) {
      TString title = condition[0];
      TString title_x = vec_var_mult[i_mult].fTitle;
      if (vec_var_mult[i_mult].fUnit != "") {
        title_x += " (" + vec_var_mult[i_mult].fUnit + ")";
      }
      TString title1d = title + ";" + title_x + ";Counts";
      gRResultHandlesFast.push_back(rdf_PartTrigger_cond.Histo1D(
          GetTH1DModelWithTitle(vec_var_mult[i_mult], title1d, tag_cond),
          vec_var_mult[i_mult].fName));
      for (int j_mult = i_mult + 1; j_mult < vec_var_mult.size(); j_mult++) {
        TString title_y = vec_var_mult[j_mult].fTitle;
        if (vec_var_mult[j_mult].fUnit != "") {
          title_y += " (" + vec_var_mult[j_mult].fUnit + ")";
        }
        TString title2d = title + ";" + title_x + ";" + title_y;
        gRResultHandlesFast.push_back(rdf_PartTrigger_cond.Histo2D(
            GetTH2DModelWithTitle(vec_var_mult[i_mult], vec_var_mult[j_mult],
                                  title2d, tag_cond),
            vec_var_mult[i_mult].fName, vec_var_mult[j_mult].fName));
      }
    }
  }

  RunGraphs(gRResultHandlesFast);

  fOutput->cd();
  RResultWrite(gRResultHandlesFast);
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
  MultREF(path_input_flowVecd, path_output, runNumber, path_calib, path_pileup);

  return 0;
}