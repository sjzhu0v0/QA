#define MRDF
#include "MALICE.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "yaml-cpp/yaml.h"
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TRandom3.h>
#include <TFile.h>
#include <TH2D.h>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

using ROOT::RVec;
using ROOT::VecOps::Take;

namespace {
TH2D* g_eff_default = nullptr;
TH2D* g_eff_pid1 = nullptr;
TH2D* g_eff_pid2 = nullptr;
TH2D* g_eff_pid3 = nullptr;

RVec<int> CountSetBitsVec(const RVec<uint8_t>& values) {
  RVec<int> out;
  out.reserve(values.size());
  for (auto value : values) {
    int count = 0;
    while (value) {
      count += value & 1;
      value >>= 1;
    }
    out.push_back(count);
  }
  return out;
}

float NDca2Dev(float pt, float dca) {
  const double dev_dca = 0.00179344 + 0.000924651 * pow(abs(pt), -1.4062);
  return abs(dca) / dev_dca;
}

RVec<float> NDca2DevVec(const RVec<float>& pt, const RVec<float>& dca) {
  RVec<float> out;
  out.reserve(pt.size());
  for (size_t i = 0; i < pt.size(); ++i)
    out.push_back(NDca2Dev(pt[i], dca[i]));
  return out;
}

RVec<float> AbsDeltaPhiVec(const RVec<float>& phi1, const RVec<float>& phi2) {
  RVec<float> out;
  out.reserve(phi1.size());
  for (size_t i = 0; i < phi1.size(); ++i) {
    float delta = phi1[i] - phi2[i];
    while (delta > M_PI)
      delta -= 2.f * M_PI;
    while (delta < -M_PI)
      delta += 2.f * M_PI;
    out.push_back(abs(delta));
  }
  return out;
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
  if (pt_bin < 1 || pt_bin > hist->GetNbinsX() || eta_bin < 1 ||
      eta_bin > hist->GetNbinsY())
    return 0.f;
  const float eff = hist->GetBinContent(pt_bin, eta_bin);
  return eff > 0.f ? 1.f / eff : 0.f;
}

RVec<float> EfficiencyWeightVec(const RVec<float>& pt, const RVec<float>& eta,
                                const std::string& setup) {
  RVec<float> out;
  out.reserve(pt.size());
  for (size_t i = 0; i < pt.size(); ++i)
    out.push_back(GetEfficiencyWeight(pt[i], eta[i], setup));
  return out;
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
  return cut == "default" || cut == "pid1" || cut == "pid2" || cut == "pid3"
             ? cut
             : "default";
}
} // namespace

void JpsiAssoRawPoiEff(TString path_input_flow,
                       TString path_input_extra,
                       TString path_output,
                       TString path_config,
                       double bootstrap_probability = 0.5,
                       unsigned int bootstrap_seed = 0,
                       TString only_cut = "") {
  YAML::Node config = YAML::LoadFile(path_config.Data());
  LoadEfficiency(config);

  TChain* tree_flow = MRootIO::OpenChain(path_input_flow.Data(), "O2dqflowvecd");
  TChain* tree_extra = MRootIO::OpenChain(path_input_extra.Data(), "ExtraInfo");
  tree_flow->AddFriend(tree_extra);

  ROOT::RDataFrame rdf(*tree_flow);
  auto df_event =
      rdf.Define("isntSameBunchPileup", MALICE::IsntSameBunchPileup_NoSlot, {"fSelection"})
          .Define("isntITSROFrameBorder", MALICE::IsntITSROFrameBorder, {"fSelection"})
          .Define("isntTimeFrameBorder", MALICE::IsntTimeFrameBorder, {"fSelection"})
          .Define("isTriggerTVX", MALICE::IsTriggerTVX, {"fSelection"})
          .Filter("isTriggerTVX && isntITSROFrameBorder && isntTimeFrameBorder && "
                  "isntSameBunchPileup && isCBT")
          .DefineSlot("bootstrap_rand",
                      [bootstrap_seed](unsigned int slot) {
                        thread_local TRandom3 rng(bootstrap_seed + slot);
                        return static_cast<float>(rng.Uniform(0., 1.));
                      })
          .Filter([bootstrap_probability](float r) { return r < bootstrap_probability; },
                  {"bootstrap_rand"}, "event-level bootstrap selection");

  auto df_pairs =
      df_event
          .Define("pair_indices",
                  [](int njpsi, int nref) {
                    std::vector<std::pair<int, int>> out;
                    out.reserve(njpsi * nref);
                    for (int i = 0; i < njpsi; ++i)
                      for (int j = 0; j < nref; ++j)
                        out.emplace_back(i, j);
                    return out;
                  },
                  {"fPT_size", "fPTREF_size"})
          .Define("jpsi_idx",
                  [](const std::vector<std::pair<int, int>>& pairs) {
                    RVec<int> out;
                    out.reserve(pairs.size());
                    for (auto [jpsi, ref] : pairs)
                      out.push_back(jpsi);
                    return out;
                  },
                  {"pair_indices"})
          .Define("ref_idx",
                  [](const std::vector<std::pair<int, int>>& pairs) {
                    RVec<int> out;
                    out.reserve(pairs.size());
                    for (auto [jpsi, ref] : pairs)
                      out.push_back(ref);
                    return out;
                  },
                  {"pair_indices"})
          .Define("jpsi_pt", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPT", "jpsi_idx"})
          .Define("jpsi_eta", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fEta", "jpsi_idx"})
          .Define("jpsi_phi", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPhi", "jpsi_idx"})
          .Define("jpsi_mass", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fMass", "jpsi_idx"})
          .Define("e1_nsig_el", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaEl1", "jpsi_idx"})
          .Define("e1_nsig_pi", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaPi1", "jpsi_idx"})
          .Define("e1_nsig_pr", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaPr1", "jpsi_idx"})
          .Define("e2_nsig_el", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaEl2", "jpsi_idx"})
          .Define("e2_nsig_pi", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaPi2", "jpsi_idx"})
          .Define("e2_nsig_pr", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNSigmaPr2", "jpsi_idx"})
          .Define("ref_pt", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPTREF", "ref_idx"})
          .Define("ref_eta", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fEtaREF", "ref_idx"})
          .Define("ref_phi", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPhiREF", "ref_idx"})
          .Define("ref_ITSChi2NCl",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fITSChi2NCl", "ref_idx"})
          .Define("ref_TPCNClsFound",
                  [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNClsFound", "ref_idx"})
          .Define("ref_itsClusterMap",
                  [](const RVec<uint8_t>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fITSClusterMap", "ref_idx"})
          .Define("ref_dcaxy", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fDcaXY", "ref_idx"})
          .Define("ref_dcaz", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fDcaZ", "ref_idx"})
          .Define("DeltaEta", "jpsi_eta - ref_eta")
          .Define("AbsDeltaPhi", AbsDeltaPhiVec, {"jpsi_phi", "ref_phi"})
          .Define("nITSCluster", CountSetBitsVec, {"ref_itsClusterMap"})
          .Define("nDcaZ2Dev", NDca2DevVec, {"ref_pt", "ref_dcaz"})
          .Define("nDcaXY2Dev", NDca2DevVec, {"ref_pt", "ref_dcaxy"})
          .Define("pair_fPosZ",
                  [](float value, const RVec<int>& idx) { return RVec<float>(idx.size(), value); },
                  {"fPosZ", "jpsi_idx"})
          .Define("pair_NumContribCalib",
                  [](double value, const RVec<int>& idx) {
                    return RVec<float>(idx.size(), static_cast<float>(value));
                  },
                  {"NumContribCalib", "jpsi_idx"})
          .Redefine("fPosZ", "pair_fPosZ")
          .Redefine("NumContribCalib", "pair_NumContribCalib");

  auto var_posz = ParseStrVar4Hist(config["hist_binning"]["fPosZ"]);
  auto var_mult = ParseStrVar4Hist(config["hist_binning"]["NumContribCalibBinned"]);
  auto var_mass = ParseStrVar4Hist(config["hist_binning"]["MassJpsiCandidate"]);
  auto var_pt = ParseStrVar4Hist(config["hist_binning"]["PtJpsiCandidate"]);
  auto var_deta = ParseStrVar4Hist(config["hist_binning"]["DeltaEtaUS"]);
  auto var_dphi = ParseStrVar4Hist(config["hist_binning"]["DeltaPhiUS"]);
  vector<StrVar4Hist> pair_vars = {var_deta, var_dphi, var_posz, var_mass, var_pt, var_mult};
  vector<StrVar4Hist> single_vars = {var_posz, var_mass, var_pt, var_mult};

  vector<std::shared_ptr<THnD>> histograms;
  for (auto it = config["cuts"].begin(); it != config["cuts"].end(); ++it) {
    const std::string cut_name = it->first.as<std::string>();
    if (only_cut != "" && cut_name != only_cut.Data())
      continue;
    const std::string cut_expr = it->second.as<std::string>();
    const std::string mask_name = "pair_mask_" + cut_name;
    auto selected =
        df_pairs.Define(mask_name, cut_expr)
            .Redefine("DeltaEta", "DeltaEta[" + mask_name + "]")
            .Redefine("AbsDeltaPhi", "AbsDeltaPhi[" + mask_name + "]")
            .Redefine("fPosZ", "fPosZ[" + mask_name + "]")
            .Redefine("jpsi_mass", "jpsi_mass[" + mask_name + "]")
            .Redefine("jpsi_pt", "jpsi_pt[" + mask_name + "]")
            .Redefine("jpsi_eta", "jpsi_eta[" + mask_name + "]")
            .Redefine("jpsi_idx", "jpsi_idx[" + mask_name + "]")
            .Redefine("NumContribCalib", "NumContribCalib[" + mask_name + "]")
            .Define("jpsi_eff_weight",
                    [setup = EfficiencySetupForCut(cut_name)](const RVec<float>& pt,
                                                               const RVec<float>& eta) {
                      return EfficiencyWeightVec(pt, eta, setup);
                    },
                    {"jpsi_pt", "jpsi_eta"});
    auto model = GetTHnDModelWithTitle(pair_vars, "", cut_name.c_str());
    auto hist = get<0>(model).GetHistogram();
    selected.Foreach(
        [hist](const RVec<float>& deta, const RVec<float>& dphi, const RVec<float>& posz,
               const RVec<float>& mass, const RVec<float>& pt, const RVec<float>& mult,
               const RVec<float>& weight) {
          double values[6];
          for (size_t i = 0; i < deta.size(); ++i) {
            values[0] = deta[i];
            values[1] = dphi[i];
            values[2] = posz[i];
            values[3] = mass[i];
            values[4] = pt[i];
            values[5] = mult[i];
            hist->Fill(values, weight[i]);
          }
        },
        {"DeltaEta", "AbsDeltaPhi", "fPosZ", "jpsi_mass", "jpsi_pt",
         "NumContribCalib", "jpsi_eff_weight"});
    histograms.push_back(hist);

    auto single_model = GetTHnDModelWithTitle(single_vars, "Single", cut_name.c_str());
    auto single_hist = get<0>(single_model).GetHistogram();
    selected.Foreach(
        [single_hist](const RVec<int>& jpsi_idx, const RVec<float>& posz,
                      const RVec<float>& mass, const RVec<float>& pt, const RVec<float>& mult,
                      const RVec<float>& weight) {
          double values[4];
          for (size_t i = 0; i < mass.size(); ++i) {
            if (i > 0 && jpsi_idx[i] == jpsi_idx[i - 1])
              continue;
            values[0] = posz[i];
            values[1] = mass[i];
            values[2] = pt[i];
            values[3] = mult[i];
            single_hist->Fill(values, weight[i]);
          }
        },
        {"jpsi_idx", "fPosZ", "jpsi_mass", "jpsi_pt", "NumContribCalib", "jpsi_eff_weight"});
    histograms.push_back(single_hist);
  }

  TFile output(path_output, "RECREATE");
  for (auto& hist : histograms)
    hist->Write();
  output.Close();
}

int main(int argc, char** argv) {
  TString path_input_flow = "../input.root";
  TString path_input_extra = "../extra.root";
  TString path_output = "output.root";
  TString path_config = "config.yaml";
  double bootstrap_probability = 0.5;
  unsigned int bootstrap_seed = 0;
  TString only_cut = "";
  if (argc > 1)
    path_input_flow = argv[1];
  if (argc > 2)
    path_input_extra = argv[2];
  if (argc > 3)
    path_output = argv[3];
  if (argc > 4)
    path_config = argv[4];
  if (argc > 5)
    bootstrap_probability = atof(argv[5]);
  if (argc > 6)
    bootstrap_seed = atoi(argv[6]);
  if (argc > 7)
    only_cut = argv[7];
  JpsiAssoRawPoiEff(path_input_flow, path_input_extra, path_output, path_config,
                    bootstrap_probability, bootstrap_seed, only_cut);
  return 0;
}
