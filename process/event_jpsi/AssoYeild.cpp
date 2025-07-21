#include "MHelper.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"

void AssoYeild(
    TString input_se_pr = "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
                          "JpsiAsso_LHC22pass4_dqfilter.root:DeltaEtaUS_"
                          "DeltaPhiUS_PosZUS_MassUS_PtUS_NumContribCalibUS",
    TString input_se_raw =
        "/home/szhu/work/alice/analysis/QA/input/jpsi/"
        "JpsiQA_LHC22pass4_dqfilter.root:fPosZ_MassUS_PtUS_NumContribCalib",
    TString input_me_pr = "/home/szhu/work/alice/analysis/QA/input/event_jpsi/"
                          "MixEventReading_LHC22pass4_dqfilter.root:DeltaEtaUS_"
                          "DeltaPhiUS_PosZUS_MassUS_PtUS_NumContribCalibUS") {
  auto hist_se_pr = MRootIO::GetObjectDiectly<THnD>(input_se_pr);
  auto hist_se_raw = MRootIO::GetObjectDiectly<THnD>(input_se_raw);
  auto hist_me_pr = MRootIO::GetObjectDiectly<THnD>(input_me_pr);

  MHnTool hnTool_se_pr(hist_se_pr);
  MHnTool hnTool_se_raw(hist_se_raw);
  hnTool_se_raw.Rebin(0, 25); // Rebin DeltaPhiUS
  MHnTool hnTool_me_pr(hist_me_pr);

  hnTool_se_pr.PrintAllAxis();
  hnTool_se_raw.PrintAllAxis();
  hnTool_me_pr.PrintAllAxis();
}