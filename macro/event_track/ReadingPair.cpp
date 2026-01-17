#define MRDF
#include "MALICE.h"
#include "MCalibration.h"
#include "MEventMixing.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "yaml-cpp/yaml.h"
#include <ROOT/RDataFrame.hxx>

void EventMixingReadingPair(TString path_input_flowVecd = "../input.root",
                            TString path_output = "output.root") {
  TFile* file_flowVecd = TFile::Open(path_input_flowVecd);
  TFile* fOutput = new TFile(path_output, "RECREATE");
  YAML::Node config = YAML::LoadFile("config.yaml");

  const double low_edge_deltaPhiToPi = config["hist_binning"]["low_edge_deltaPhiToPi"].as<double>();
  const double up_edge_deltaPhiToPi = config["hist_binning"]["up_edge_deltaPhiToPi"].as<double>();

  TTree* tree_input = nullptr;
  TList* list_keys = file_flowVecd->GetListOfKeys();
  for (int i = 0; i < list_keys->GetEntries(); i++) {
    TKey* key = (TKey*)list_keys->At(i);
    if (strcmp(key->GetClassName(), "TTree") == 0) {
      tree_input = (TTree*)file_flowVecd->Get(key->GetName());
      break;
    }
  }

  ROOT::RDataFrame rdf(*tree_input);

  auto rdf_AllVar =
      rdf.Define("DeltaPhi", "ref1_phi - ref2_phi").Define("DeltaEta", "ref1_eta - ref2_eta");
  if (is_interactive())
    ROOT::RDF::Experimental::AddProgressBar(rdf);

  StrVar4Hist var_fPosZ("fPosZ", "#it{V}_{Z}", "cm", 8, {-10, 10});
  StrVar4Hist var_NumContribCalibBinned("NumContribCalib", "N_{vtx contrib} Calibrated", "", 10,
                                        {0, 5, 8, 11, 14, 18, 23, 28, 36, 48, 300});
  int n_bins_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["n_bins"].as<int>();
  double min_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["min"].as<double>();
  double max_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["max"].as<double>();
  StrVar4Hist DeltaEtaRef("DeltaEta", "#Delta#eta_{track, track}", "", n_bins_deltaEta_assoYield,
                          {min_deltaEta_assoYield, max_deltaEta_assoYield});
  int n_bins_deltaPhi_assoYield = config["hist_binning"]["n_bins_deltaPhi_assoYield"].as<int>();
  StrVar4Hist DeltaPhiRef("DeltaPhi", "#Delta#phi_{track, track}", "", n_bins_deltaPhi_assoYield,
                          {low_edge_deltaPhiToPi * M_PI, up_edge_deltaPhiToPi * M_PI});

#define obj2push_thnd(rdf2push, ...)                                                               \
  do {                                                                                             \
    TupleTHnDModel tuple_thnd = GetTHnDModelWithTitle(__VA_ARGS__);                                \
    gRResultHandles.push_back(rdf2push.HistoND(get<0>(tuple_thnd), get<1>(tuple_thnd)));           \
  } while (0)
  obj2push_thnd(rdf_AllVar, {DeltaEtaRef, DeltaPhiRef, var_fPosZ, var_NumContribCalibBinned});
  fOutput->cd();
  RResultWrite(gRResultHandles);
  fOutput->Close();
  cout << "Output written to " << path_output << endl;
}

int main(int argc, char** argv) {
  TString path_input_flowVecd = argv[1];
  TString path_output = argv[2];
  EventMixingReadingPair(path_input_flowVecd, path_output);
  return 0;
}
