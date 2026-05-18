#define main JpsiAssoMEPoiEffMain
#include "JpsiAsso_me_poi_eff.cpp"
#undef main

#include <TStopwatch.h>
#include <cstdlib>

struct FillTimeStats {
  Long64_t index_entries_total = 0;
  Long64_t index_entries_processed = 0;
  Long64_t mixed_event_pairs = 0;
  Long64_t jpsi_candidates = 0;
  Long64_t ref_tracks = 0;
  Long64_t pair_fills = 0;
  Long64_t single_fills = 0;
};

FillTimeStats SampleMixedEventHistogramFill(TString path_input_flow, TString path_input_index,
                                            TString path_config, Long64_t sample_index_entries,
                                            TString only_cut = "") {
  if (sample_index_entries <= 0)
    throw std::invalid_argument("sample_index_entries must be positive");

  YAML::Node config = YAML::LoadFile(path_config.Data());
  LoadEfficiency(config);

  auto var_posz = ParseStrVar4Hist(config["hist_binning"]["fPosZ"]);
  auto var_mult = ParseStrVar4Hist(config["hist_binning"]["NumContribCalibBinned"]);
  auto var_mass = ParseStrVar4Hist(config["hist_binning"]["MassJpsiCandidate"]);
  auto var_pt = ParseStrVar4Hist(config["hist_binning"]["PtJpsiCandidate"]);
  auto var_deta = ParseStrVar4Hist(config["hist_binning"]["DeltaEtaUS"]);
  auto var_dphi = ParseStrVar4Hist(config["hist_binning"]["DeltaPhiUS"]);
  vector<StrVar4Hist> pair_vars = {var_deta, var_dphi, var_posz, var_mass, var_pt, var_mult};
  vector<StrVar4Hist> single_vars = {var_posz, var_mass, var_pt, var_mult};

  std::vector<HistSet> hist_sets;
  for (auto it = config["cuts"].begin(); it != config["cuts"].end(); ++it) {
    const std::string cut_name = it->first.as<std::string>();
    if (only_cut != "" && cut_name != only_cut.Data())
      continue;
    hist_sets.push_back({cut_name, FormulaExpr(it->second.as<std::string>()),
                         EfficiencySetupForCut(cut_name), MakeTHnD(pair_vars, "", cut_name.c_str()),
                         MakeTHnD(single_vars, "Single", cut_name.c_str())});
  }

  TChain* tree_a = MRootIO::OpenChain(path_input_flow.Data(), "O2dqflowvecd");
  TChain* tree_b = MRootIO::OpenChain(path_input_flow.Data(), "O2dqflowvecd");
  TChain* tree_index = MRootIO::OpenChain(path_input_index.Data(), "EventMixing");

  FillTimeStats stats;
  stats.index_entries_total = tree_index->GetEntries();

  TTreeReader pairs_reader(tree_index);
  TTreeReaderValue<std::vector<std::pair<ULong64_t, ULong64_t>>> mixed_events(pairs_reader,
                                                                              "MixedEvent");
  TTreeReaderValue<int> i_mult(pairs_reader, "IndexMixing_NumContribCalib");
  TTreeReaderValue<int> i_posz(pairs_reader, "IndexMixing_PosZ");

  TTreeReader reader_a(tree_a);
  TTreeReaderArray<float> fPT(reader_a, "fPT");
  TTreeReaderArray<float> fEta(reader_a, "fEta");
  TTreeReaderArray<float> fPhi(reader_a, "fPhi");
  TTreeReaderArray<float> fMass(reader_a, "fMass");
  TTreeReaderArray<float> fTPCNSigmaEl1(reader_a, "fTPCNSigmaEl1");
  TTreeReaderArray<float> fTPCNSigmaPi1(reader_a, "fTPCNSigmaPi1");
  TTreeReaderArray<float> fTPCNSigmaPr1(reader_a, "fTPCNSigmaPr1");
  TTreeReaderArray<float> fTPCNSigmaEl2(reader_a, "fTPCNSigmaEl2");
  TTreeReaderArray<float> fTPCNSigmaPi2(reader_a, "fTPCNSigmaPi2");
  TTreeReaderArray<float> fTPCNSigmaPr2(reader_a, "fTPCNSigmaPr2");

  TTreeReader reader_b(tree_b);
  TTreeReaderArray<float> fPTREF(reader_b, "fPTREF");
  TTreeReaderArray<float> fEtaREF(reader_b, "fEtaREF");
  TTreeReaderArray<float> fPhiREF(reader_b, "fPhiREF");
  TTreeReaderArray<float> fITSChi2NCl(reader_b, "fITSChi2NCl");
  TTreeReaderArray<float> fTPCNClsFound(reader_b, "fTPCNClsFound");
  TTreeReaderArray<uint8_t> fITSClusterMap(reader_b, "fITSClusterMap");
  TTreeReaderArray<float> fDcaXY(reader_b, "fDcaXY");
  TTreeReaderArray<float> fDcaZ(reader_b, "fDcaZ");

  const auto& mult_bins = var_mult.fBins;
  const auto& posz_bins = var_posz.fBins;
  while (pairs_reader.Next() && stats.index_entries_processed < sample_index_entries) {
    ++stats.index_entries_processed;
    const float mult_value = MixBinCenter(mult_bins, *i_mult);
    const float posz_value = MixBinCenter(posz_bins, *i_posz);
    for (const auto& [event_a, event_b] : *mixed_events) {
      ++stats.mixed_event_pairs;
      reader_a.SetEntry(event_a);
      reader_b.SetEntry(event_b);
      for (int i_jpsi = 0; i_jpsi < fPT.GetSize(); ++i_jpsi) {
        ++stats.jpsi_candidates;
        std::vector<bool> single_filled(hist_sets.size(), false);
        for (int i_ref = 0; i_ref < fPTREF.GetSize(); ++i_ref) {
          ++stats.ref_tracks;
          const int n_its = CountSetBits(fITSClusterMap[i_ref]);
          const float n_dcaz = NDca2Dev(fPTREF[i_ref], fDcaZ[i_ref]);
          const float n_dcaxy = NDca2Dev(fPTREF[i_ref], fDcaXY[i_ref]);
          const double pair_values[] = {fEta[i_jpsi] - fEtaREF[i_ref],
                                        AbsDeltaPhi(fPhi[i_jpsi], fPhiREF[i_ref]),
                                        posz_value,
                                        fMass[i_jpsi],
                                        fPT[i_jpsi],
                                        mult_value};
          const double single_values[] = {posz_value, fMass[i_jpsi], fPT[i_jpsi], mult_value};
          for (size_t i_hist = 0; i_hist < hist_sets.size(); ++i_hist) {
            auto& hist = hist_sets[i_hist];
            if (!PassCut(hist, fTPCNSigmaEl1[i_jpsi], fTPCNSigmaPi1[i_jpsi], fTPCNSigmaPr1[i_jpsi],
                         fTPCNSigmaEl2[i_jpsi], fTPCNSigmaPi2[i_jpsi], fTPCNSigmaPr2[i_jpsi],
                         fITSChi2NCl[i_ref], fTPCNClsFound[i_ref], n_its, n_dcaz, n_dcaxy,
                         fPTREF[i_ref], posz_value))
              continue;
            const float cut_weight =
                GetEfficiencyWeight(fPT[i_jpsi], fEta[i_jpsi], hist.efficiency_setup);
            hist.pair_hist->Fill(pair_values, cut_weight);
            ++stats.pair_fills;
            if (!single_filled[i_hist]) {
              hist.single_hist->Fill(single_values, cut_weight);
              single_filled[i_hist] = true;
              ++stats.single_fills;
            }
          }
        }
      }
    }
  }
  return stats;
}

void EstimateJpsiAssoMEPoiEffFillTime(TString path_input_flow, TString path_input_index,
                                      TString path_config, Long64_t sample_index_entries = 1000,
                                      TString only_cut = "") {
  TStopwatch timer;
  timer.Start();
  const FillTimeStats stats = SampleMixedEventHistogramFill(
      path_input_flow, path_input_index, path_config, sample_index_entries, only_cut);
  timer.Stop();

  if (stats.index_entries_processed == 0)
    throw std::runtime_error("no index entries were processed");

  const double elapsed_s = timer.RealTime();
  const double scale = static_cast<double>(stats.index_entries_total) / stats.index_entries_processed;
  const double estimated_total_s = elapsed_s * scale;
  cout << "sampled index entries: " << stats.index_entries_processed << " / "
       << stats.index_entries_total << endl;
  cout << "mixed event pairs: " << stats.mixed_event_pairs << endl;
  cout << "J/psi candidates: " << stats.jpsi_candidates << endl;
  cout << "reference tracks visited: " << stats.ref_tracks << endl;
  cout << "pair histogram fills: " << stats.pair_fills << endl;
  cout << "single histogram fills: " << stats.single_fills << endl;
  cout << "sample wall time [s]: " << elapsed_s << endl;
  cout << "estimated full histogram fill wall time [s]: " << estimated_total_s << endl;
  cout << "estimated full histogram fill wall time [min]: " << estimated_total_s / 60. << endl;
}

int main(int argc, char** argv) {
  TString path_input_flow = "../input.root";
  TString path_input_index = "output_mix_index.root";
  TString path_config = "config.yaml";
  Long64_t sample_index_entries = 1000;
  TString only_cut = "";
  if (argc > 1)
    path_input_flow = argv[1];
  if (argc > 2)
    path_input_index = argv[2];
  if (argc > 3)
    path_config = argv[3];
  if (argc > 4)
    sample_index_entries = std::atoll(argv[4]);
  if (argc > 5)
    only_cut = argv[5];
  EstimateJpsiAssoMEPoiEffFillTime(path_input_flow, path_input_index, path_config,
                                   sample_index_entries, only_cut);
  return 0;
}
