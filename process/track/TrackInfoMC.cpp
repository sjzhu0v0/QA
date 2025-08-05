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
  TTree *tree_info = MRootIO::OpenChain(file_input, "O2mctrackinfo");

  StrVar4Hist var_InversePt("InversePt", "1/p_{T}", "GeV/c", 90, {0.2, 7});
  StrVar4Hist var_DCAxy("fDcaXY", "DCA_{XY}", "cm", 100, {-0.04, 0.04});
  StrVar4Hist var_DCAz("fDcaZ", "DCA_{Z}", "cm", 100, {-0.04, 0.04});
  StrVar4Hist var_Pt("fPt", "p_{T}", "GeV/c", 100, {0.1, 6});
  StrVar4Hist var_DeltaPtOverPt("DeltaPtOverPt",
                                "(p^{rec}_{T}-p^{MC}_{T})/p^{rec}_{T}", "", 100,
                                {-0.3, 0.3});

  ROOT::RDataFrame rdf(*tree_info);

  auto rdf_full =
      rdf.Define("InversePt",
                 [](const float &pt) { return pt > 0 ? 1. / pt : 0.; }, {"fPt"})
          .Define("DeltaPtOverPt",
                  [](const float &pt, const float &mcpt) {
                    return pt > 0 ? (pt - mcpt) / pt : 0.;
                  },
                  {"fPt", "fMcPt"});

  gRResultHandles.push_back(
      rdf_full.Histo2D(GetTH2DModel(var_InversePt, var_DCAxy),
                       var_InversePt.fName, var_DCAxy.fName));
  gRResultHandles.push_back(
      rdf_full.Histo2D(GetTH2DModel(var_InversePt, var_DCAz),
                       var_InversePt.fName, var_DCAz.fName));
  gRResultHandles.push_back(
      rdf_full.Histo2D(GetTH2DModel(var_Pt, var_DeltaPtOverPt), var_Pt.fName,
                       var_DeltaPtOverPt.fName));

#define th1dPush(name)                                                         \
  gRResultHandles.push_back(                                                   \
      rdf_full.Histo1D(GetTH1DModelWithTitle2(var_##name), var_##name.fName));

  th1dPush(InversePt);
  th1dPush(DCAxy);
  th1dPush(DCAz);
  th1dPush(Pt);
  th1dPush(DeltaPtOverPt);

  //  OBJ: TLeafF    fPosX   fPosX : 0 at: 0x3ea36f0
  //  OBJ: TLeafF    fPosY   fPosY : 0 at: 0x3fbdfb0
  //  OBJ: TLeafF    fPosZ   fPosZ : 0 at: 0x3fbe2b0
  //  OBJ: TLeafS    fNumContrib     fNumContrib : 0 at: 0x3fbe5b0
  //  OBJ: TLeafF    fMCPosX fMCPosX : 0 at: 0x39bcfb0
  //  OBJ: TLeafF    fMCPosY fMCPosY : 0 at: 0x39bd2b0
  //  OBJ: TLeafF    fMCPosZ fMCPosZ : 0 at: 0x39bd5b0
  //  OBJ: TLeafF    fPt     fPt : 0 at: 0x39bdb10
  //  OBJ: TLeafF    fEta    fEta : 0 at: 0x39be070
  //  OBJ: TLeafF    fPhi    fPhi : 0 at: 0x39be5d0
  //  OBJ: TLeafI    fSign   fSign : 0 at: 0x39beb30
  //  OBJ: TLeafF    fDcaXY  fDcaXY : 0 at: 0x3fc1820
  //  OBJ: TLeafF    fDcaZ   fDcaZ : 0 at: 0x3fc1d80
  //  OBJ: TLeafB    fITSClusterMap  fITSClusterMap : 0 at: 0x3fc2310
  //  OBJ: TLeafF    fMcPt   fMcPt : 0 at: 0x3fc2870
  //  OBJ: TLeafF    fMcEta  fMcEta : 0 at: 0x3fc2dd0
  //  OBJ: TLeafF    fMcPhi  fMcPhi : 0 at: 0x3fc3330
  //  OBJ: TLeafI    fPdgCode        fPdgCode : 0 at: 0x3b370a0
  //  OBJ: TLeafF    fVx     fVx : 0 at: 0x3b37600
  //  OBJ: TLeafF    fVy     fVy : 0 at: 0x3b37b60
  //  OBJ: TLeafF    fVz     fVz : 0 at: 0x3b380c0
  //  OBJ: TLeafF    fVt     fVt : 0 at: 0x3b38620

  RunGraphs(gRResultHandles);
  file_output->cd();
  RResultWrite(gRResultHandles);
  file_output->Close();
}


int main(int argc, char **argv) {
  TString path_input = "/home/szhu/work/alice/analysis/QA/input/track/"
                       "trackInfoMC_24fd4b_550367.root";
  TString path_output = "/home/szhu/work/alice/analysis/QA/output/track/"
                        "TrackInfoMC_24fd4b_550367.root";

  if (argc > 1) {
    path_input = argv[1];
  }
  if (argc > 2) {
    path_output = argv[2];
  }

  TrackInfoMC(path_input, path_output);
  return 0;
}