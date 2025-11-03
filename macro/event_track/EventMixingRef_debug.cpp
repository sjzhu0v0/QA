#define MRDF
#include <ROOT/RDataFrame.hxx>

#include "iostream"

using namespace std;

void EventMixingRef(TString path_input_flowVecd = "../input.root",
                    TString path_input_mult = "../input2.root",
                    TString path_output = "output.root",
                    TString path_output_tree = "output_tree.root") {
  // close multi-thread
  ROOT::EnableImplicitMT(1);

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
  // TChain *tree_mult = MRootIO::OpenChain(path_input_mult, "MultCalib");
  // TTree *tree_mult = GetObjectSingle<TTree>(
  //     "/lustre/alice/users/szhu/job/QA/LHC24pass1_DiElectron_Group/JpsiQA_mult/"
  //     "JpsiQA_mult_554207.root :MultCalib");
  TFile *file_mult =
      TFile::Open("/lustre/alice/users/szhu/job/QA/LHC24pass1_DiElectron_Group/"
                  "JpsiQA_mult/JpsiQA_mult_554207.root");
  TTree *tree_mult = (TTree *)file_mult->Get("MultCalib");

  // tree_flowVecd->AddFriend(tree_mult);

  tree_mult->Print();

  ROOT::RDataFrame rdf_flowVecd(*tree_mult);
  ROOT::RDataFrame rdf_mult(
      "MultCalib",
      "/lustre/alice/users/szhu/job/QA/LHC24pass1_DiElectron_Group/"
      "JpsiQA_mult/JpsiQA_mult_554207.root");
  // auto report = rdf_flowVecd.Report();
  // report->Print();

  // auto rdf_witTrigger =
  //     rdf_flowVecd.Define("map_trigger", MALICE::triggermapRVec,
  //     {"fSelection"});

  auto hist_mult = rdf_mult.Histo1D("NumContribCalib");
  hist_mult->Draw();
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