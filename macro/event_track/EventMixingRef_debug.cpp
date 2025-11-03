#define MRDF
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TH1.h>
#include <iostream>

using namespace std;

void EventMixingRef(TString path_input_flowVecd = "../input.root",
                    TString path_input_mult = "../input2.root",
                    TString path_output = "output.root",
                    TString path_output_tree = "output_tree.root") {
  // Disable multi-threading
  ROOT::EnableImplicitMT(1);

  cout << "Input file (flowVecd): " << path_input_flowVecd << endl;
  cout << "Input file (mult):     " << path_input_mult << endl;
  cout << "Output file:           " << path_output << endl;
  cout << "Output tree file:      " << path_output_tree << endl;

  // Use RDataFrame directly from file path (recommended)
  ROOT::RDataFrame rdf_mult("MultCalib", path_input_mult.Data());

  // Create histogram
  auto hist_mult = rdf_mult.Histo1D({"NumContribCalib_hist", "NumContribCalib;N_{contrib};Entries", 1000, 0, 10000}, "NumContribCalib");

  // Save to output file
  TFile outFile(path_output, "RECREATE");
  hist_mult->Write();
  outFile.Close();

  cout << "Histogram saved to " << path_output << endl;
}

int main(int argc, char **argv) {
  TString path_input_flowVecd = "../input.root";
  TString path_input_mult = "../input2.root";
  TString path_output = "output.root";
  TString path_output_tree = "output_tree.root";

  if (argc > 1) path_input_flowVecd = argv[1];
  if (argc > 2) path_input_mult     = argv[2];
  if (argc > 3) path_output         = argv[3];
  if (argc > 4) path_output_tree    = argv[4];

  gROOT->SetBatch(kTRUE); // No graphics

  EventMixingRef(path_input_flowVecd, path_input_mult, path_output, path_output_tree);

  return 0;
}