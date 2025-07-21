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
  // Axis 0: axis0, title: #Delta#eta_{J/#psi, track}  nbins:80
  // Axis 1: axis1, title: #Delta#phi_{J/#psi, track}  nbins:10
  // Axis 2: axis2, title: #it{V}_{Z} (cm)  nbins:8
  // Axis 3: axis3, title: M_{ee} (GeV^{2}/c^{4})  nbins:90
  // Axis 4: axis4, title: p_{T} (GeV/c)  nbins:10
  // Axis 5: axis5, title: N_{vtx contrib} Calibrated  nbins:10
  hnTool_se_raw.PrintAllAxis();
  // Axis 0: axis0, title: #it{V}_{Z} (cm)  nbins:8
  // Axis 1: axis1, title: M_{ee} (GeV^{2}/c^{4})  nbins:100
  // Axis 2: axis2, title: p_{T} (GeV/c)  nbins:10
  // Axis 3: axis3, title: N_{vtx contrib} Calibrated  nbins:10
  hnTool_me_pr.PrintAllAxis();
  // Axis 0: axis0, title: #Delta#eta_{J/#psi, track}  nbins:80
  // Axis 1: axis1, title: #Delta#phi_{J/#psi, track}  nbins:10
  // Axis 2: axis2, title: #it{V}_{Z} (cm)  nbins:8
  // Axis 3: axis3, title: M_{ee} (GeV^{2}/c^{4})  nbins:90
  // Axis 4: axis4, title: p_{T} (GeV/c)  nbins:10
  // Axis 5: axis5, title: N_{vtx contrib} Calibrated  nbins:10

  
}