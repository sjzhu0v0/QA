#define MRDF
#include "MALICE.h"
#include "MEventMixing.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "opt/EventData.h"
#include "yaml-cpp/yaml.h"
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TFormula.h>
#include <TH2D.h>
#include <THn.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {
TH2D* g_eff_default = nullptr;
TH2D* g_eff_pid1 = nullptr;
TH2D* g_eff_pid2 = nullptr;
TH2D* g_eff_pid3 = nullptr;

int CountSetBits(uint8_t value) {
  int count = 0;
  while (value) {
    count += value & 1;
    value >>= 1;
  }
  return count;
}

float NDca2Dev(float pt, float dca) {
  const double dev_dca = 0.00179344 + 0.000924651 * pow(abs(pt), -1.4062);
  return abs(dca) / dev_dca;
}

float AbsDeltaPhi(float phi1, float phi2) {
  float delta = phi1 - phi2;
  while (delta > M_PI)
    delta -= 2.f * M_PI;
  while (delta < -M_PI)
    delta += 2.f * M_PI;
  return abs(delta);
}

float GetEfficiencyWeight(float pt, float eta, const std::string& setup) {
  TH2D* hist = g_eff_default;
  if (setup == "pid1" && g_eff_pid1)
    hist = g_eff_pid1;
  else if (setup == "pid2" && g_eff_pid2)
    hist = g_eff_pid2;
  else if (setup == "pid3" && g_eff_pid3)
    hist = g_eff_pid3;
  if (!hist)
    return 1.f;
  const int pt_bin = hist->GetXaxis()->FindFixBin(pt);
  const int eta_bin = hist->GetYaxis()->FindFixBin(eta);
  if (pt_bin < 1 || pt_bin > hist->GetNbinsX() || eta_bin < 1 || eta_bin > hist->GetNbinsY())
    return 0.f;
  const float eff = hist->GetBinContent(pt_bin, eta_bin);
  return eff > 0.f ? 1.f / eff : 0.f;
}

void LoadEfficiency(const YAML::Node& config) {
  const std::string path = config["efficiency_correction"]["file"].as<std::string>();
  auto file = TFile::Open(path.c_str());
  if (!file || file->IsZombie())
    throw std::runtime_error("cannot open efficiency file: " + path);
  const std::unordered_map<std::string, std::string> names = {
      {"default", "jpsi_reconstruction_efficiency_pt_eta_default_low_eff_removed"},
      {"pid1", "jpsi_reconstruction_efficiency_pt_eta_pid1_low_eff_removed"},
      {"pid2", "jpsi_reconstruction_efficiency_pt_eta_pid2_low_eff_removed"},
      {"pid3", "jpsi_reconstruction_efficiency_pt_eta_pid3_low_eff_removed"}};
  auto load = [&](const std::string& setup) {
    auto* hist = dynamic_cast<TH2D*>(file->Get(names.at(setup).c_str()));
    if (!hist)
      throw std::runtime_error("missing efficiency histogram for " + setup);
    auto* clone = dynamic_cast<TH2D*>(hist->Clone(("eff_" + setup).c_str()));
    clone->SetDirectory(nullptr);
    return clone;
  };
  g_eff_default = load("default");
  g_eff_pid1 = load("pid1");
  g_eff_pid2 = load("pid2");
  g_eff_pid3 = load("pid3");
  file->Close();
}

std::string EfficiencySetupForCut(const std::string& cut) {
  return cut == "default" || cut == "pid1" || cut == "pid2" || cut == "pid3" ? cut : "default";
}

int FindMixBin(double value, const std::vector<double>& bins) {
  for (size_t i = 0; i + 1 < bins.size(); ++i) {
    if (value >= bins[i] && value < bins[i + 1])
      return static_cast<int>(i);
  }
  return -1;
}

float MixBinCenter(const std::vector<double>& bins, int index) {
  return static_cast<float>(0.5 * (bins[index] + bins[index + 1]));
}

std::unique_ptr<THnD> MakeTHnD(const std::vector<StrVar4Hist>& vars, const TString& title,
                               const TString& tag) {
  TString name = vars[0].fName;
  for (size_t i = 1; i < vars.size(); ++i)
    name += "_" + vars[i].fName;
  if (tag != "")
    name += "_" + tag;

  TString full_title = title;
  full_title += ";";
  std::vector<int> nbins;
  for (const auto& var : vars) {
    TString axis_title = var.fTitle;
    if (var.fUnit != "")
      axis_title += " (" + var.fUnit + ")";
    full_title += axis_title + ";";
    nbins.push_back(var.fNbins);
  }
  auto hist = std::make_unique<THnD>(name, full_title, vars.size(), nbins.data(), nullptr, nullptr);
  for (size_t i = 0; i < vars.size(); ++i)
    hist->SetBinEdges(i, vars[i].fBins.data());
  return hist;
}

struct HistSet {
  std::string cut_name;
  std::string cut_expr;
  std::string efficiency_setup;
  std::unique_ptr<THnD> pair_hist;
  std::unique_ptr<THnD> single_hist;
};

bool PassCut(const HistSet& hist, float e1_nsig_el, float e1_nsig_pi, float e1_nsig_pr,
             float e2_nsig_el, float e2_nsig_pi, float e2_nsig_pr, float ref_ITSChi2NCl,
             float ref_TPCNClsFound, int nITSCluster, float nDcaZ2Dev, float nDcaXY2Dev,
             float ref_pt, float fPosZ) {
  // The config uses a fixed family of scalar variables. Compile each expression
  // once through TFormula to preserve the existing YAML cut surface.
  static std::unordered_map<std::string, std::unique_ptr<TFormula>> formulas;
  auto it = formulas.find(hist.cut_name);
  if (it == formulas.end()) {
    auto formula = std::make_unique<TFormula>(hist.cut_name.c_str(), hist.cut_expr.c_str());
    it = formulas.emplace(hist.cut_name, std::move(formula)).first;
  }
  double values[] = {e1_nsig_el,     e1_nsig_pi,       e1_nsig_pr,
                     e2_nsig_el,     e2_nsig_pi,       e2_nsig_pr,
                     ref_ITSChi2NCl, ref_TPCNClsFound, static_cast<double>(nITSCluster),
                     nDcaZ2Dev,      nDcaXY2Dev,       ref_pt,
                     fPosZ};
  return it->second->EvalPar(values) != 0.;
}

std::string FormulaExpr(const std::string& expr) {
  std::string out = expr;
  const std::vector<std::pair<std::string, std::string>> replacements = {
      {"e1_nsig_el", "x[0]"},     {"e1_nsig_pi", "x[1]"},       {"e1_nsig_pr", "x[2]"},
      {"e2_nsig_el", "x[3]"},     {"e2_nsig_pi", "x[4]"},       {"e2_nsig_pr", "x[5]"},
      {"ref_ITSChi2NCl", "x[6]"}, {"ref_TPCNClsFound", "x[7]"}, {"nITSCluster", "x[8]"},
      {"nDcaZ2Dev", "x[9]"},      {"nDcaXY2Dev", "x[10]"},      {"ref_pt", "x[11]"},
      {"fPosZ", "x[12]"},         {"abs(", "TMath::Abs("}};
  for (const auto& [from, to] : replacements) {
    size_t pos = 0;
    while ((pos = out.find(from, pos)) != std::string::npos) {
      out.replace(pos, from.size(), to);
      pos += to.size();
    }
  }
  return out;
}
} // namespace

vector<pair<ULong64_t, ULong64_t>> MixEvent(unsigned int, const int id, const ULong64_t& event_id) {
  return MixVec<pair<ULong64_t, ULong64_t>, ULong64_t>(
      id, event_id, [](const ULong64_t& a, const ULong64_t& b) { return make_pair(a, b); }, 100);
}

void EventMixingIndexGen(TString path_input_flow, TString path_input_extra,
                         TString path_output_index) {
  ROOT::DisableImplicitMT();
  StrVar4Hist var_posz_mix("fPosZ", "#it{V}_{Z}", "cm", 10, {-10, 10});
  StrVar4Hist var_mult_mix("NumContribCalib", "N_{vtx contrib} Calibrated", "", 10,
                           {0, 5, 8, 11, 14, 18, 23, 28, 36, 48, 300});
  const auto bins_mix_numContrib = var_mult_mix.fBins;
  const auto bins_mix_posZ = var_posz_mix.fBins;

  TChain* tree_flow = MRootIO::OpenChain(path_input_flow.Data(), "O2dqflowvecd");
  TChain* tree_extra = MRootIO::OpenChain(path_input_extra.Data(), "ExtraInfo");
  tree_flow->AddFriend(tree_extra);

  ROOT::RDataFrame rdf(*tree_flow);
  auto df =
      rdf.Define("isntSameBunchPileup", MALICE::IsntSameBunchPileup_NoSlot, {"fSelection"})
          .Define("isntITSROFrameBorder", MALICE::IsntITSROFrameBorder, {"fSelection"})
          .Define("isntTimeFrameBorder", MALICE::IsntTimeFrameBorder, {"fSelection"})
          .Define("isTriggerTVX", MALICE::IsTriggerTVX, {"fSelection"})
          .Filter("isTriggerTVX && isntITSROFrameBorder && isntTimeFrameBorder && "
                  "isntSameBunchPileup && isCBT")
          .Filter("fEta_size>=1", "has Jpsi")
          .Define("EventData", CreateEventData,
                  {"fMultTPC", "fMultTracklets", "fMultNTracksPV", "fMultFT0C", "fNumContrib",
                   "NumContribCalib", "fPosX", "fPosY", "fPosZ", "fSelection", "fHadronicRate",
                   "fPT", "fEta", "fPhi", "fMass", "fSign", "fPTREF", "fEtaREF", "fPhiREF"})
          .Define("isEventGood", [](const EventData& event) { return event.isGood(); },
                  {"EventData"})
          .Filter("isEventGood", "Event is good");

  auto rdf_PartTriggerWithJpsiWithEventWithEventMixing =
      df.Define("IndexMixing_NumContribCalib",
                [bins_mix_numContrib](double numContrib) {
                  int index = -1;
                  for (int i = 0; i < bins_mix_numContrib.size() - 1; i++) {
                    if (numContrib >= bins_mix_numContrib[i] &&
                        numContrib < bins_mix_numContrib[i + 1]) {
                      index = i;
                      break;
                    }
                  }
                  if (index == -1) {
                    return -1; // Invalid index
                  }
                  return int(index);
                },
                {"NumContribCalib"})
          .Filter("IndexMixing_NumContribCalib>=0", "valid NumContribCalib index")
          .Define("IndexMixing_PosZ",
                  [bins_mix_posZ](float posZ) {
                    int index = -1;
                    for (int i = 0; i < bins_mix_posZ.size() - 1; i++) {
                      if (posZ >= bins_mix_posZ[i] && posZ < bins_mix_posZ[i + 1]) {
                        index = i;
                        break;
                      }
                    }
                    if (index == -1) {
                      return -1; // Invalid index
                    }
                    return int(index);
                  },
                  {"fPosZ"})
          .Filter("IndexMixing_PosZ>=0", "valid PosZ index")
          .Define("IndexMixing",
                  [bins_mix_posZ](int indexNumContrib, int indexPosZ) {
                    if (indexNumContrib < 0 || indexPosZ < 0) {
                      return -1; // Invalid index
                    }
                    return int(indexPosZ + indexNumContrib * bins_mix_posZ.size());
                  },
                  {"IndexMixing_NumContribCalib", "IndexMixing_PosZ"})
          .DefineSlot("MixedEvent", MixEvent, {"IndexMixing", "rdfentry_"});
  if (is_interactive())
    ROOT::RDF::Experimental::AddProgressBar(rdf_PartTriggerWithJpsiWithEventWithEventMixing);
  rdf_PartTriggerWithJpsiWithEventWithEventMixing.Snapshot(
      "EventMixing", path_output_index,
      {"MixedEvent", "IndexMixing_NumContribCalib", "IndexMixing_PosZ"});
}

void FillMixedEventHistograms(TString path_input_flow, TString path_input_index,
                              TString path_output_hist, TString path_config,
                              TString only_cut = "") {
  YAML::Node config = YAML::LoadFile(path_config.Data());
  LoadEfficiency(config);
  cout << "start filling" << endl;

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
  while (pairs_reader.Next()) {
    const float mult_value = MixBinCenter(mult_bins, *i_mult);
    const float posz_value = MixBinCenter(posz_bins, *i_posz);
    for (const auto& [event_a, event_b] : *mixed_events) {
      reader_a.SetEntry(event_a);
      reader_b.SetEntry(event_b);
      for (int i_jpsi = 0; i_jpsi < fPT.GetSize(); ++i_jpsi) {
        std::vector<bool> single_filled(hist_sets.size(), false);
        for (int i_ref = 0; i_ref < fPTREF.GetSize(); ++i_ref) {
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
            if (!single_filled[i_hist]) {
              hist.single_hist->Fill(single_values, cut_weight);
              single_filled[i_hist] = true;
            }
          }
        }
      }
    }
  }

  TFile output(path_output_hist, "RECREATE");
  for (auto& hist : hist_sets) {
    hist.pair_hist->Write();
    hist.single_hist->Write();
  }
  output.Close();
}

void JpsiAssoMEPoiEff(TString path_input_flow, TString path_input_extra, TString path_output_hist,
                      TString path_output_index, TString path_config, TString only_cut = "") {
  EventMixingIndexGen(path_input_flow, path_input_extra, path_output_index);
  FillMixedEventHistograms(path_input_flow, path_output_index, path_output_hist, path_config,
                           only_cut);
}

int main(int argc, char** argv) {
  TString path_input_flow = "../input.root";
  TString path_input_extra = "../extra.root";
  TString path_output_hist = "output.root";
  TString path_output_index = "output_mix_index.root";
  TString path_config = "config.yaml";
  TString only_cut = "";
  if (argc > 1)
    path_input_flow = argv[1];
  if (argc > 2)
    path_input_extra = argv[2];
  if (argc > 3)
    path_output_hist = argv[3];
  if (argc > 4)
    path_output_index = argv[4];
  if (argc > 5)
    path_config = argv[5];
  if (argc > 6)
    only_cut = argv[6];
  JpsiAssoMEPoiEff(path_input_flow, path_input_extra, path_output_hist, path_output_index,
                   path_config, only_cut);
  return 0;
}
