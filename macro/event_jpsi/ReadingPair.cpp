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

  TTree* tree_input = (TTree*)file_flowVecd->Get("EventMixing");

  ROOT::RDataFrame rdf(*tree_input);

  auto rdf_AllVar =
      rdf.Define("DeltaPhi", "jpsi_phi - ref_phi").Define("DeltaEta", "jpsi_eta - ref_eta");
  if (is_interactive())
    ROOT::RDF::Experimental::AddProgressBar(rdf);

  StrVar4Hist var_fPosZ("fPosZ", "#it{V}_{Z}", "cm", 8, {-10, 10});
  StrVar4Hist var_NumContribCalibBinned("NumContribCalib", "N_{vtx contrib} Calibrated", "", 10,
                                        {0, 5, 8, 11, 14, 18, 23, 28, 36, 48, 300});
  int n_bins_mass_assoYield = config["hist_binning"]["n_bins_mass_assoYield"].as<int>();
  StrVar4Hist var_MassJpsiCandidate("jpsi_mass", "M_{ee}", "GeV^{2}/c^{4}", n_bins_mass_assoYield,
                                    {1.8, 5.4});
  StrVar4Hist var_PtJpsiCandidate("jpsi_pt", "p_{T}", "GeV/c", 10, {0., 5.});
  StrVar4Hist var_PtJpsiCandidateFine("jpsi_pt", "p_{T}", "GeV/c", 10, {0., 5.});
  int n_bins_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["n_bins"].as<int>();
  double min_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["min"].as<double>();
  double max_deltaEta_assoYield =
      config["hist_binning"]["binning_deltaEta_assoYield"]["max"].as<double>();
  StrVar4Hist var_DeltaEtaUS("DeltaEta", "#Delta#eta_{J/#psi, track}", "",
                             n_bins_deltaEta_assoYield,
                             {min_deltaEta_assoYield, max_deltaEta_assoYield});
  int n_bins_deltaPhi_assoYield = config["hist_binning"]["n_bins_deltaPhi_assoYield"].as<int>();
  StrVar4Hist var_DeltaPhiUS("DeltaPhi", "#Delta#phi_{J/#psi, track}", "",
                             n_bins_deltaPhi_assoYield,
                             {low_edge_deltaPhiToPi * M_PI, up_edge_deltaPhiToPi * M_PI});

#define obj2push_thnd(rdf2push, ...)                                                               \
  do {                                                                                             \
    TupleTHnDModel tuple_thnd = GetTHnDModelWithTitle(__VA_ARGS__);                                \
    gRResultHandles.push_back(rdf2push.HistoND(get<0>(tuple_thnd), get<1>(tuple_thnd)));           \
  } while (0)

  obj2push_thnd(rdf_AllVar, {var_DeltaEtaUS, var_DeltaPhiUS, var_fPosZ, var_MassJpsiCandidate,
                             var_PtJpsiCandidate, var_NumContribCalibBinned});
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
