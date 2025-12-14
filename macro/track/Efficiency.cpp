#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "opt/EventData.h"
#include "yaml-cpp/yaml.h"
#include <ROOT/RDataFrame.hxx>

/* root [4] O2mctrackinfo->Print()
******************************************************************************
*Tree    :O2mctrackinfo: O2mctrackinfo * *Entries :    11823 : Total = 1137413
bytes  File  Size =     435733 *
*        :          : Tree compression factor =   2.59                       *
******************************************************************************
*Br    0 :fPosX     : fPosX/F                                                *
*Entries :    11823 : Total  Size=      47858 bytes  File Size  =       3708 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=  12.78     *
*............................................................................*
*Br    1 :fPosY     : fPosY/F                                                *
*Entries :    11823 : Total  Size=      47858 bytes  File Size  =       3667 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=  12.92     *
*............................................................................*
*Br    2 :fPosZ     : fPosZ/F                                                *
*Entries :    11823 : Total  Size=      47858 bytes  File Size  =       4368 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=  10.85     *
*............................................................................*
*Br    3 :fSelection : fSelection/l                                          *
*Entries :    11823 : Total  Size=      95183 bytes  File Size  =       2038 *
*Baskets :        1 : Basket Size=      95608 bytes  Compression=  46.45     *
*............................................................................*
*Br    4 :fNumContrib : fNumContrib/s                                        *
*Entries :    11823 : Total  Size=      24238 bytes  File Size  =       1966 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=  12.07     *
*............................................................................*
*Br    5 :fMCPosX   : fMCPosX/F                                              *
*Entries :    11823 : Total  Size=      47868 bytes  File Size  =       3701 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=  12.80     *
*............................................................................*
*Br    6 :fMCPosY   : fMCPosY/F                                              *
*Entries :    11823 : Total  Size=      47868 bytes  File Size  =       3706 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=  12.78     *
*............................................................................*
*Br    7 :fMCPosZ   : fMCPosZ/F                                              *
*Entries :    11823 : Total  Size=      47868 bytes  File Size  =       4367 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=  10.85     *
*............................................................................*
*Br    8 :fPt       : fPt/F                                                  *
*Entries :    11823 : Total  Size=      47848 bytes  File Size  =      42047 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=   1.13     *
*............................................................................*
*Br    9 :fEta      : fEta/F                                                 *
*Entries :    11823 : Total  Size=      47853 bytes  File Size  =      43709 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=   1.08     *
*............................................................................*
*Br   10 :fPhi      : fPhi/F                                                 *
*Entries :    11823 : Total  Size=      47853 bytes  File Size  =      42036 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=   1.13     *
*............................................................................*
*Br   11 :fSign     : fSign/I                                                *
*Entries :    11823 : Total  Size=      47858 bytes  File Size  =       3061 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=  15.48     *
*............................................................................*
*Br   12 :fDcaXY    : fDcaXY/F                                               *
*Entries :    11823 : Total  Size=      47863 bytes  File Size  =      43492 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=   1.09     *
*............................................................................*
*Br   13 :fDcaZ     : fDcaZ/F                                                *
*Entries :    11823 : Total  Size=      47858 bytes  File Size  =      38900 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=   1.22     *
*............................................................................*
*Br   14 :fITSClusterMap : fITSClusterMap/b                                  *
*Entries :    11823 : Total  Size=      12428 bytes  File Size  =       3877 *
*Baskets :        1 : Basket Size=      32000 bytes  Compression=   3.07     *
*............................................................................*
*Br   15 :fMcPt     : fMcPt/F                                                *
*Entries :    11823 : Total  Size=      47858 bytes  File Size  =      42011 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=   1.13     *
*............................................................................*
*Br   16 :fMcEta    : fMcEta/F                                               *
*Entries :    11823 : Total  Size=      47863 bytes  File Size  =      43672 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=   1.08     *
*............................................................................*
*Br   17 :fMcPhi    : fMcPhi/F                                               *
*Entries :    11823 : Total  Size=      47863 bytes  File Size  =      41871 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=   1.13     *
*............................................................................*
*Br   18 :fPdgCode  : fPdgCode/I                                             *
*Entries :    11823 : Total  Size=      47873 bytes  File Size  =       6164 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=   7.69     *
*............................................................................*
*Br   19 :fVx       : fVx/F                                                  *
*Entries :    11823 : Total  Size=      47848 bytes  File Size  =      14984 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=   3.16     *
*............................................................................*
*Br   20 :fVy       : fVy/F                                                  *
*Entries :    11823 : Total  Size=      47848 bytes  File Size  =      15010 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=   3.16     *
*............................................................................*
*Br   21 :fVz       : fVz/F                                                  *
*Entries :    11823 : Total  Size=      47848 bytes  File Size  =      14756 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=   3.21     *
*............................................................................*
*Br   22 :fVt       : fVt/F                                                  *
*Entries :    11823 : Total  Size=      47848 bytes  File Size  =      11161 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression=   4.24     *
*............................................................................*
*Br   23 :fMcDecision : fMcDecision/i                                        *
*Entries :    11823 : Total  Size=      47888 bytes  File Size  =        118 *
*Baskets :        1 : Basket Size=      48316 bytes  Compression= 401.52     *
*............................................................................*

root [6] O2mctrktruth->Print()
******************************************************************************
*Tree    :O2mctrktruth: O2mctrktruth                                           *
*Entries :   230465 : Total =        16142699 bytes  File  Size =    3389593 *
*        :          : Tree compression factor =   4.76                       *
******************************************************************************
*Br    0 :fPosX     : fPosX/F                                                *
*Entries :   230465 : Total  Size=     922425 bytes  File Size  =       4697 *
*Baskets :        1 : Basket Size=     922884 bytes  Compression= 196.28     *
*............................................................................*
*Br    1 :fPosY     : fPosY/F                                                *
*Entries :   230465 : Total  Size=     922425 bytes  File Size  =       4688 *
*Baskets :        1 : Basket Size=     922884 bytes  Compression= 196.66     *
*............................................................................*
*Br    2 :fPosZ     : fPosZ/F                                                *
*Entries :   230465 : Total  Size=     922425 bytes  File Size  =       5534 *
*Baskets :        1 : Basket Size=     922884 bytes  Compression= 166.60     *
*............................................................................*
*Br    3 :fSelection : fSelection/l                                          *
*Entries :   230465 : Total  Size=    1844318 bytes  File Size  =       3227 *
*Baskets :        1 : Basket Size=    1844744 bytes  Compression= 571.37     *
*............................................................................*
*Br    4 :fNumContrib : fNumContrib/s                                        *
*Entries :   230465 : Total  Size=     461521 bytes  File Size  =       2276 *
*Baskets :        1 : Basket Size=     461954 bytes  Compression= 202.56     *
*............................................................................*
*Br    5 :fMCPosX   : fMCPosX/F                                              *
*Entries :   230465 : Total  Size=     922435 bytes  File Size  =       4689 *
*Baskets :        1 : Basket Size=     922884 bytes  Compression= 196.62     *
*............................................................................*
*Br    6 :fMCPosY   : fMCPosY/F                                              *
*Entries :   230465 : Total  Size=     922435 bytes  File Size  =       4695 *
*Baskets :        1 : Basket Size=     922884 bytes  Compression= 196.37     *
*............................................................................*
*Br    7 :fMCPosZ   : fMCPosZ/F                                              *
*Entries :   230465 : Total  Size=     922435 bytes  File Size  =       5522 *
*Baskets :        1 : Basket Size=     922884 bytes  Compression= 166.96     *
*............................................................................*
*Br    8 :fPt       : fPt/F                                                  *
*Entries :   230465 : Total  Size=     922415 bytes  File Size  =     830320 *
*Baskets :        1 : Basket Size=     922884 bytes  Compression=   1.11     *
*............................................................................*
*Br    9 :fEta      : fEta/F                                                 *
*Entries :   230465 : Total  Size=     922420 bytes  File Size  =     842871 *
*Baskets :        1 : Basket Size=     922884 bytes  Compression=   1.09     *
*............................................................................*
*Br   10 :fPhi      : fPhi/F                                                 *
*Entries :   230465 : Total  Size=     922420 bytes  File Size  =     812235 *
*Baskets :        1 : Basket Size=     922884 bytes  Compression=   1.14     *
*............................................................................*
*Br   11 :fPdgCode  : fPdgCode/I                                             *
*Entries :   230465 : Total  Size=     922440 bytes  File Size  =      88287 *
*Baskets :        1 : Basket Size=     922884 bytes  Compression=  10.44     *
*............................................................................*
*Br   12 :fVx       : fVx/F                                                  *
*Entries :   230465 : Total  Size=     922415 bytes  File Size  =     153274 *
*Baskets :        1 : Basket Size=     922884 bytes  Compression=   6.01     *
*............................................................................*
*Br   13 :fVy       : fVy/F                                                  *
*Entries :   230465 : Total  Size=     922415 bytes  File Size  =     155170 *
*Baskets :        1 : Basket Size=     922884 bytes  Compression=   5.94     *
*............................................................................*
*Br   14 :fVz       : fVz/F                                                  *
*Entries :   230465 : Total  Size=     922415 bytes  File Size  =     145555 *
*Baskets :        1 : Basket Size=     922884 bytes  Compression=   6.33     *
*............................................................................*
*Br   15 :fVt       : fVt/F                                                  *
*Entries :   230465 : Total  Size=     922415 bytes  File Size  =     325278 *
*Baskets :        1 : Basket Size=     922884 bytes  Compression=   2.83     *
*............................................................................*
*Br   16 :fMcDecision : fMcDecision/i                                        *
*Entries :   230465 : Total  Size=     922455 bytes  File Size  =        195 *
*Baskets :        1 : Basket Size=     922884 bytes  Compression= 4727.93     *
*............................................................................*
*/

void Efficiency(TString path_input, TString path_output) {
  TFile *file_input = new TFile(path_input, "READ");
  TFile *file_output = new TFile(path_output, "RECREATE");

  TChain *file_mc = MRootIO::OpenChain(path_input.Data(), "O2mctrktruth");
  TChain *file_reco = MRootIO::OpenChain(path_input.Data(), "O2mctrackinfo");

  ROOT::RDataFrame rdf_mc(*file_mc);
  ROOT::RDataFrame rdf_reco(*file_reco);

  auto rdf_mc_withTrigger =
      rdf_mc.Define("map_trigger", MALICE::triggermapRVec, {"fSelection"})
          .Define("isntSPDPileup", MALICE::IsntSPDPileup, {"fSelection"})
          .Define("isntTPCPileup", MALICE::IsntTPCPileup, {"fSelection"})
          .Define("isntSameBunchPileup", MALICE::IsntSameBunchPileup_NoSlot,
                  {"fSelection"})
          .Define("isntITSROFrameBorder", MALICE::IsntITSROFrameBorder,
                  {"fSelection"})
          .Define("isntTimeFrameBorder", MALICE::IsntTimeFrameBorder,
                  {"fSelection"})
          .Define("isTriggerTVX", MALICE::IsTriggerTVX, {"fSelection"});
  auto rdf_reco_withTrigger =
      rdf_reco.Define("map_trigger", MALICE::triggermapRVec, {"fSelection"})
          .Define("isntSPDPileup", MALICE::IsntSPDPileup, {"fSelection"})
          .Define("isntTPCPileup", MALICE::IsntTPCPileup, {"fSelection"})
          .Define("isntSameBunchPileup", MALICE::IsntSameBunchPileup_NoSlot,
                  {"fSelection"})
          .Define("isntITSROFrameBorder", MALICE::IsntITSROFrameBorder,
                  {"fSelection"})
          .Define("isntTimeFrameBorder", MALICE::IsntTimeFrameBorder,
                  {"fSelection"})
          .Define("isTriggerTVX", MALICE::IsTriggerTVX, {"fSelection"});

  auto rdf_mc_selected =
      rdf_mc_withTrigger.Filter("isTriggerTVX", "is Trigger TVX")
          .Filter("isntITSROFrameBorder", "no ITS RO Frame border")
          .Filter("isntTimeFrameBorder", "no Time Frame border")
          .Filter("isntSameBunchPileup", "no Time Frame border")
          .Filter(
              "abs(fPdgCode)==2212 || abs(fPdgCode) == 321 || abs(fPdgCode) == "
              "211 || abs(fPdgCode) == 13 || abs(fPdgCode) == 11",
              "is Final");
  auto rdf_reco_selected =
      rdf_reco_withTrigger.Filter("isTriggerTVX", "is Trigger TVX")
          .Filter("isntITSROFrameBorder", "no ITS RO Frame border")
          .Filter("isntTimeFrameBorder", "no Time Frame border")
          .Filter("isntSameBunchPileup", "no Time Frame border");

  auto rdf_mc_selected_kine =
      rdf_mc_selected.Filter("abs(fEta)<0.9 && fPt < 5", "abs(eta)<0.9 && fPt < 5");
  auto rdf_reco_selected_kine =
      rdf_reco_selected.Filter("abs(fEta)<0.9 && fPt < 5", "abs(eta)<0.9 && fPt < 5");

  StrVar4Hist var_mult("fNumContrib", "N_{vtx contrib}", "", 10,
                       {0, 5, 8, 11, 14, 18, 23, 28, 36, 48, 300});
  StrVar4Hist var_vz("fPosZ", "V_{z}", "cm", 20, {-10, 10});
  StrVar4Hist var_pt(
      "fPt", "P_{T}", "GeV/c", 14,
      {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.5, 2., 3., 4., 5.});
  StrVar4Hist var_eta("fEta", "#eta", "", 18, {-0.9, 0.9});

#define obj2push_thnd(rdf2push, ...)                                           \
  do {                                                                         \
    TupleTHnDModel tuple_thnd = GetTHnDModelWithTitle(__VA_ARGS__);            \
    gRResultHandles.push_back(                                                 \
        rdf2push.HistoND(get<0>(tuple_thnd), get<1>(tuple_thnd)));             \
  } while (0)

  gRResultHandles.push_back(rdf_mc_selected_kine.Histo3D(
      GetTH3DM(var_vz, var_pt, var_eta, "mc", "mc"), var_vz.fName.Data(),
      var_pt.fName.Data(), var_eta.fName.Data()));
  gRResultHandles.push_back(rdf_reco_selected_kine.Histo3D(
      GetTH3DM(var_vz, var_pt, var_eta, "reco", "reco"), var_vz.fName.Data(),
      var_pt.fName.Data(), var_eta.fName.Data()));

  gRResultHandles.push_back(rdf_reco_selected_kine.Histo2D(
      GetTH2DM(var_vz, var_mult, "reco", "reco"), var_vz.fName.Data(),
      var_mult.fName.Data());
  gRResultHandles.push_back(rdf_mc_selected_kine.Histo2D(
      GetTH2DM(var_vz, var_mult, "mc", "mc"), var_vz.fName.Data(),
      var_mult.fName.Data());

  gRResultHandles.push_back(rdf_reco_selected_kine.Histo2D(
      GetTH2DM(var_pt, var_mult, "reco", "reco"), var_pt.fName.Data(),
      var_mult.fName.Data()));
  gRResultHandles.push_back(
      rdf_mc_selected_kine.Histo2D(GetTH2DM(var_pt, var_mult, "mc", "mc"),
                                   var_pt.fName.Data(), var_mult.fName.Data()));
  gRResultHandles.push_back(rdf_mc_selected_kine.Histo2D(
      GetTH2DM(var_eta, var_mult, "mc", "mc"), var_eta.fName.Data(),
      var_mult.fName.Data()));
  gRResultHandles.push_back(rdf_reco_selected_kine.Histo2D(
      GetTH2DM(var_eta, var_mult, "reco", "reco"), var_eta.fName.Data(),
      var_mult.fName.Data()));

  gRResultHandles.push_back(rdf_mc_selected_kine.Histo1D(
      GetTH1DM(var_mult, "mc_mult", "mc_mult, |#eta|<0.9,1GeV/c<p_{T}<5GeV/c"),
      var_mult.fName.Data()));
  gRResultHandles.push_back(rdf_reco_selected_kine.Histo1D(
      GetTH1DM(var_mult, "reco_mult",
               "reco_mult, |#eta|<0.9,1GeV/c<p_{T}<5GeV/c"),
      var_mult.fName.Data()));

  RunGraphs(gRResultHandles);

  file_output->cd();
  RResultWrite(gRResultHandles);
  file_output->Close();
}

int main(int argc, char **argv) {
  TString path_input = "../input.root";
  TString path_output = "output.root";

  if (argc > 1) {
    path_input = argv[1];
  }
  if (argc > 2) {
    path_output = argv[2];
  }

  Efficiency(path_input, path_output);
  return 0;
}
