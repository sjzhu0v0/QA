#define MRDF
#include "MALICE.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include <ROOT/RDataFrame.hxx>

void Trigger_Mult(TString path_input = "../input.root",
                  TString path_output = "output.root") {
  TFile *file_event = TFile::Open(path_input);
  TFile *fOutput = new TFile(path_output, "RECREATE");

  // find if O2reextended or O2reextended_001 is in the file
  TString treeName = "O2reextended";
  if (file_event->Get("O2reextended_001")) {
    treeName = "O2reextended_001";
  }

  TTree *tree_event = (TTree *)file_event->Get("O2reducedevent");
  TTree *tree_event_ext = (TTree *)file_event->Get(treeName);

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
          .Define("isTriggerTVX", MALICE::IsTriggerTVX, {"fSelection"});
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

  gRResultHandlesFast.push_back(rdf_witTrigger.Histo1D(
      {"map_trigger", "map_trigger", 51, 0, 51}, "map_trigger"));
  /* #region trigger */
  gRResultHandlesFast.push_back(rdf_witTrigger.Histo1D("isntSPDPileup"));
  gRResultHandlesFast.push_back(rdf_witTrigger.Histo1D("isntTPCPileup"));
  gRResultHandlesFast.push_back(rdf_witTrigger.Histo1D("isntSameBunchPileup"));
  gRResultHandlesFast.push_back(rdf_witTrigger.Histo1D(
      {"map_trigger", "map_trigger", 51, 0, 51}, "map_trigger"));
  gRResultHandlesFast.push_back(rdf_isntITSROFrameBorder.Histo1D(
      {"map_trigger_isntITSROFrameBorder", "map_trigger:isntITSROFrameBorder",
       51, 0, 51},
      "map_trigger"));
  gRResultHandlesFast.push_back(rdf_isntTimeFrameBorder.Histo1D(
      {"map_trigger_isntTimeFrameBorder", "map_trigger:isntTimeFrameBorder", 51,
       0, 51},
      "map_trigger"));
  gRResultHandlesFast.push_back(rdf_isTriggerTVX.Histo1D(
      {"map_trigger_isTriggerTVX", "map_trigger:isTriggerTVX", 51, 0, 51},
      "map_trigger"));
  /* #endregion */

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
      {"fMultNTracksPV_fMultTPC",
       "fMultNTracksPV_fMultTPC;fMultNTracksPV;fMultTPC;Counts", 150, 0, 150,
       300, 0, 300},
      "fMultNTracksPV", "fMultTPC"));
  /* #endregion */
  RunGraphs(gRResultHandlesFast);

  fOutput->cd();
  RResultWrite(gRResultHandlesFast);
  fOutput->Close();
}