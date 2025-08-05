#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MEventMixing.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include <ROOT/RDataFrame.hxx>

void TrackInfoMC(
    TString path_input = "/home/szhu/work/alice/analysis/QA/input/track/"
                         "trackInfoMC_24fd4b_550367.root",
    TString path_output = "/home/szhu/work/alice/analysis/QA/output/track/"
                          "TrackInfoMC_24fd4b_550367.root") {

  TFile *file_input = new TFile(path_input);
  TFile *file_output = new TFile(path_output, "RECREATE");
  TTree *tree_cut = (TTree *)file_input->Get("O2barreltrackcuts");
  TTree *tree_info = (TTree *)file_input->Get("O2mctrackinfo");

  tree_info->AddFriend(tree_cut);
  ROOT::RDataFrame rdf(*tree_info);

#define push1d(name) gRResultHandles.push_back(rdf.Histo1D(#name));
  push1d(fPosX);
  push1d(fPosY);
  push1d(fPosZ);
  push1d(fNumContrib);
  push1d(fMCPosX);
  push1d(fMCPosY);
  push1d(fMCPosZ);
  push1d(fPt);
  push1d(fEta);
  push1d(fPhi);
  push1d(fSign);
  push1d(fDcaXY);
  push1d(fDcaZ);
  push1d(fITSClusterMap);
  push1d(fMcPt);
  push1d(fMcEta);
  push1d(fMcPhi);
  push1d(fPdgCode);
  push1d(fVx);
  push1d(fVy);
  push1d(fVz);
  push1d(fVt);

  RunGraphs(gRResultHandles);
  file_output->cd();
  RResultWrite(gRResultHandles);
  file_output->Close();
}