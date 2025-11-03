#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MEventMixing.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "opt/EventData.h"
#include <ROOT/RDataFrame.hxx>

vector<EventDataREF> MixEvent(unsigned int, const int id,
                              const EventData &event_info) {
  return MixVec<EventDataREF, EventData>(
      id, event_info,
      [](const EventData &a, const EventData &b) {
        EventDataREF event;
        event.event_info.Copy(a.event_info);
        event.event_info2.Copy(b.event_info);
        event.track_info.Copy(a.track_info);
        event.track_info2.Copy(b.track_info);
        return event;
      },
      100);
}

void EventMixingRef(TString path_input_flowVecd = "../input.root",
                    TString path_input_mult = "../input2.root",
                    TString path_output = "output.root",
                    TString path_output_tree = "output_tree.root") {
  // close multi-thread
  ROOT::EnableImplicitMT(1);

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
      {0, 7, 12, 17, 22, 29, 37, 46, 57, 73, 300});
  StrVar4Hist var_fMultTPC("fMultTPC", "Mult_{TPC}", "", 600, {0, 600});
  StrVar4Hist var_fMultREF("fMultREF", "Mult_{REF}", "", 100, {0, 100});
  StrVar4Hist var_fMultFT0C("fMultFT0C", "Mult_{FT0C}", "", 130,
                            {-1000., 12000.});
  StrVar4Hist var_MultNTracksPV("fMultNTracksPV", "#it{N}_{Tracks PV}", "", 150,
                                {0, 150});
  StrVar4Hist var_MassJpsiCandidate("fMass", "M_{ee}", "GeV^{2}/c^{4}", 100,
                                    {1.8, 5.4});
  StrVar4Hist var_PtJpsiCandidate("fPT", "p_{T}", "GeV/c", 10, {0., 5.});

  const vector<double> bins_mix_numContrib = var_NumContribCalibBinned.fBins;
  const vector<double> bins_mix_posZ = var_fPosZMix.fBins;

  // TFile *file_flowVecd = TFile::Open(path_input_flowVecd);
  // TFile *file_mult = TFile::Open(path_input_mult);

  // Calib_NumContrib_fPosZ_Run::GetHistCali(path_calib, runNumber);
  // Cut_MultTPC_NumContrib::init(path_pileup);

  // TChain *tree_flowVecd = MRootIO::OpenChain(file_flowVecd, "O2dqflowvecd");
  // TChain *tree_mult = MRootIO::OpenChain(file_mult, "MultCalib");

  cout << "Input file: " << path_input_flowVecd << endl;
  cout << "Input file: " << path_input_mult << endl;
  cout << "Output file: " << path_output << endl;
  cout << "Output tree file: " << path_output_tree << endl;
  TChain *tree_flowVecd =
      MRootIO::OpenChain(path_input_flowVecd, "O2dqflowvecd");
  TChain *tree_mult = MRootIO::OpenChain(path_input_mult, "MultCalib");

  // tree_flowVecd->AddFriend(tree_mult);

  tree_flowVecd->Print();
  tree_mult->Print();

  ROOT::RDataFrame rdf(*tree_flowVecd);
  auto report = rdf.Report();
  report->Print();

  // auto rdf_witTrigger =
  //     rdf.Define("map_trigger", MALICE::triggermapRVec, {"fSelection"});

  // rdf_witTrigger.Histo1D("map_trigger")->Draw();
}

int main(int argc, char **argv) {
  TString path_input_flowVecd = "../input.root";
  TString path_input_mult = "../input2.root";
  TString path_output = "output.root";
  TString path_output_tree = "output_tree.root";

  if (argc > 1) {
    path_input_flowVecd = argv[1];
  }
  if (argc > 2) {
    path_input_mult = argv[2];
  }
  if (argc > 3) {
    path_output = argv[3];
  }
  if (argc > 4) {
    path_output_tree = argv[4];
  }

  gROOT->SetBatch(kTRUE); // Disable interactive graphics
  EventMixingRef(path_input_flowVecd, path_input_mult, path_output,
                 path_output_tree);

  return 0;
}