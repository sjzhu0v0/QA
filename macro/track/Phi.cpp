#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "opt/EventData.h"
#include "yaml-cpp/yaml.h"
#include <RCTSelectionFlags.h>
#include <ROOT/RDataFrame.hxx>

void Ref_BS(TString path_input_flowVecd = "../input.root", TString path_output = "output.root") {
  TFile* fOutput = new TFile(path_output, "RECREATE");

  TChain* tree_flowVecd = MRootIO::OpenChain(path_input_flowVecd.Data(), "O2dqflowvecd");
  TChain* tree_mult = MRootIO::OpenChain(path_input_flowVecd.Data(), "O2rctrawdq");
  cout << "Input file: " << path_input_flowVecd << endl;
  cout << "Output file: " << path_output << endl;

  tree_flowVecd->AddFriend(tree_mult);

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
          .Define("isCBT",
                  [](const uint32_t& fRct) {
                    thread_local o2::aod::rctsel::RCTFlagsChecker rctChecker{"CBT"};
                    thread_local MCollision mCollision;
                    mCollision.RctRaw = fRct;
                    return rctChecker(mCollision);
                  },
                  {"fRct"})
          .Define("RunNumber", [] { return float(0.5); });
  auto rdf_isntITSROFrameBorder =
      rdf_witTrigger.Filter("isntITSROFrameBorder", "no ITS RO Frame border");
  auto rdf_isntTimeFrameBorder =
      rdf_witTrigger.Filter("isntTimeFrameBorder", "no Time Frame border");
  auto rdf_isTriggerTVX = rdf_witTrigger.Filter("isTriggerTVX", "is Trigger TVX");
  auto rdf_PartTrigger = rdf_witTrigger.Filter("isTriggerTVX", "is Trigger TVX")
                             .Filter("isntITSROFrameBorder", "no ITS RO Frame border")
                             .Filter("isntTimeFrameBorder", "no Time Frame border")
                             .Filter("isntSameBunchPileup", "no Time Frame border")
                             .Filter("isCBT")
      /*  .Filter("isntSelfDefinedPileup", "no self defined pileup") */;

  //         .Define("EventData", CreateEventData,
  //                 {"fMultTPC", "fMultTracklets", "fMultNTracksPV",
  //                 "fMultFT0C",
  //                  "fNumContrib", "NumContribCalib", "fPosX", "fPosY",
  //                  "fPosZ", "fSelection", "fHadronicRate", "fPT", "fEta",
  //                  "fPhi", "fMass", "fSign", "fPTREF", "fEtaREF", "fPhiREF"})

  if (is_interactive())
    ROOT::RDF::Experimental::AddProgressBar(rdf);

  // StrVar4Hist var_fPosZ("fPosZ", "#it{V}_{Z}", "cm", 40, {-10, 10});
  // StrVar4Hist var_NumContribCalib("NumContribCalib", "N_{vtx contrib} Calibrated", "", 300,
  //                                 {0, 300});
  // StrVar4Hist var_pt("fPTREF", "p_{T}", "GeV/c", 29, {0.2, 6});
  // StrVar4Hist var_eta("fEtaREF", "#eta", "", 36, {-0.9, 0.9});
  StrVar4Hist var_phi("fPhiREF", "#phi", "rad", 200, {0, 2 * M_PI});
  // StrVar4Hist var_dcaxy("fDcaXY", "DCA_{XY}", "cm", 200, {-0.1, 0.1});
  // StrVar4Hist var_dcaz("fDcaZ", "DCA_{XY}", "cm", 200, {-0.1, 0.1});

#define obj2push_thnd(rdf2push, ...)                                                               \
  do {                                                                                             \
    TupleTHnDModel tuple_thnd = GetTHnDModelWithTitle(__VA_ARGS__);                                \
    gRResultHandles.push_back(rdf2push.HistoND(get<0>(tuple_thnd), get<1>(tuple_thnd)));           \
  } while (0)

#define obj2push_th2d(rdf2push, varX, varY, ...)                                                   \
  do {                                                                                             \
    gRResultHandles.push_back(                                                                     \
        rdf2push.Histo2D(GetTH2DM(varX, varY, __VA_ARGS__), varX.fName, varY.fName));              \
  } while (0)

#define obj2push_th3d(rdf2push, varX, varY, varZ, ...)                                             \
  do {                                                                                             \
    gRResultHandles.push_back(rdf2push.Histo3D(GetTH3DM(varX, varY, varZ, __VA_ARGS__),            \
                                               varX.fName, varY.fName, varZ.fName));               \
  } while (0)

  // obj2push_th2d(rdf_PartTrigger, var_pt, var_phi, "", "");
  // obj2push_th2d(rdf_PartTrigger, var_eta, var_phi, "", "");
  // obj2push_th2d(rdf_PartTrigger, var_fPosZ, var_phi, "", "");
  // obj2push_th2d(rdf_PartTrigger, var_dcaxy, var_phi, "", "");
  //
  // obj2push_th3d(rdf_PartTrigger, var_dcaz, var_phi, var_pt, "", "");
  // obj2push_th3d(rdf_PartTrigger, var_dcaz, var_phi, var_eta, "", "");
  // obj2push_th3d(rdf_PartTrigger, var_dcaxy, var_phi, var_pt, "", "");
  // obj2push_th3d(rdf_PartTrigger, var_dcaxy, var_phi, var_eta, "", "");
  //
  // obj2push_th2d(rdf_PartTrigger, var_NumContribCalib, var_phi, "", "");
  // obj2push_th2d(rdf_PartTrigger, var_eta, var_fPosZ, "", "");

  gRResultHandles.push_back(rdf_PartTrigger.Histo1D(GetTH1DM(var_phi, "", ""), var_phi.fName));

  RunGraphs(gRResultHandles);

  fOutput->cd();
  RResultWrite(gRResultHandles);
  fOutput->Close();
}

int main(int argc, char** argv) {
  TString path_input_flowVecd = "../input.root";
  TString path_output = "output.root";

  if (argc > 1) {
    path_input_flowVecd = argv[1];
  }
  if (argc > 2) {
    path_output = argv[3];
  }

  Ref_BS(path_input_flowVecd, path_output);

  return 0;
}
