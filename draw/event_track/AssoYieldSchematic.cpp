#include "MHist.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TFitResult.h"
#include "yaml-cpp/yaml.h"

#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMatrixDSym.h"
#include <cmath>
#include <stdexcept>
#include <utility>

void AssoYieldSchematic(
    TString path_se = "/home/szhu/work/alice/analysis/QA/input/event_track/"
                      "SameEventReading_BS_1_24pass1.root:DeltaEtaRef_"
                      "DeltaPhiRef_PosZRef_NumContribCalibRef",
    TString path_me = "/home/szhu/work/alice/analysis/QA/input/event_track/"
                      "MixEventReading_24pass1.root:DeltaEtaRef_"
                      "DeltaPhiRef_PosZRef_NumContribCalibRef",
    TString path_trigger =
        "/home/szhu/work/alice/analysis/QA/input/event_track/"
        "SameEventReading_BS_1_24pass1.root:PosZRefSingle_"
        "NumContribCalibRefSingle") {
  auto hn_se = MRootIO::GetObjectDiectly<THnD>(path_se);
  auto hn_me = MRootIO::GetObjectDiectly<THnD>(path_me);
  auto hn_trigger = MRootIO::GetObjectDiectly<THnD>(path_trigger);

  MHnTool hnTool_se(hn_se);
  MHnTool hnTool_me(hn_me);
  hnTool_se.PrintAllAxis();
  hnTool_me.PrintAllAxis();

  auto deltaEta_deltaPhi_se = hnTool_se.Project(0, 1, {0, 0});
  auto deltaEta_deltaPhi_me = hnTool_me.Project(0, 1, {0, 0});

  auto localStyle = [](TH2D *h) {
    h->GetYaxis()->SetRangeUser(-1.799999, 1.799999);
  };

  localStyle(deltaEta_deltaPhi_se);
  localStyle(deltaEta_deltaPhi_me);

  MRootGraphic::StyleCommon();
  MRootGraphic::StyleHistCommon(deltaEta_deltaPhi_se);
  MRootGraphic::StyleHistCommon(deltaEta_deltaPhi_me);

  TCanvas *c1 = new TCanvas("c1", "c1", 1500, 500);
  c1->Divide(3, 1);
  c1->cd(1);
  StyleFlow::DeltaPhi_DeltaEta((TPad *)gPad, deltaEta_deltaPhi_se);

  c1->cd(2);
  StyleFlow::DeltaPhi_DeltaEta((TPad *)gPad, deltaEta_deltaPhi_me);

  auto deltaEta_deltaPhi_se_scaled = (TH2D *)deltaEta_deltaPhi_se->Clone();
  deltaEta_deltaPhi_se_scaled->Scale(1. / hn_trigger->GetEntries());

  auto deltaEta_deltaPhi_me_normalized = (TH2D *)deltaEta_deltaPhi_me->Clone();
  deltaEta_deltaPhi_me_normalized->Scale(deltaEta_deltaPhi_me->GetNbinsX() *
                                         deltaEta_deltaPhi_me->GetNbinsY() /
                                         deltaEta_deltaPhi_me->GetEntries());
  auto assoYield = (TH2D *)deltaEta_deltaPhi_se_scaled->Clone("assoYield");
  assoYield->Divide(deltaEta_deltaPhi_me_normalized);
  assoYield->GetZaxis()->SetRangeUser(12.e-3, 14.e-3);
  localStyle(assoYield);
  MRootGraphic::StyleHistCommon(assoYield);
  c1->cd(3);
  StyleFlow::DeltaPhi_DeltaEta((TPad *)gPad, assoYield);

  // TCanvas *c2 = new TCanvas("c2", "c2", 800, 400);
  // c2->Divide(2, 1);
  // c2->cd(1);
  // StyleFlow::DeltaPhi_DeltaEta((TPad *)gPad, deltaEta_deltaPhi_se_scaled);
  // c2->cd(2);
  // StyleFlow::DeltaPhi_DeltaEta((TPad *)gPad,
  // deltaEta_deltaPhi_me_normalized);
}