#define MRDF
#include "MALICE.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "yaml-cpp/yaml.h"
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/RActionImpl.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TH2D.h>
#include <TProfile3D.h>
#include <TRandom3.h>
#include <cmath>
#include <cstdint>
#include <memory>
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

bool UseEfficiencyCorrection(const YAML::Node& config, bool unit_efficiency) {
  if (unit_efficiency)
    return false;
  const auto eff_config = config["efficiency_correction"];
  return !(eff_config && eff_config["enabled"] &&
           !eff_config["enabled"].as<bool>());
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

bool ParseBoolArg(const char* value) {
  const std::string text(value);
  return text == "1" || text == "true" || text == "True" || text == "TRUE" ||
         text == "yes" || text == "on" || text == "unit_efficiency" ||
         text == "unit-efficiency" || text == "no_eff" || text == "no-eff";
}

std::shared_ptr<TProfile3D> MakeProfile(const StrVar4Hist& mass,
                                        const StrVar4Hist& pt,
                                        const StrVar4Hist& mult,
                                        const std::string& cut_name) {
  const TString name = Form("%s_%s_%s_%s", mass.fName.Data(), pt.fName.Data(),
                            mult.fName.Data(), cut_name.c_str());
  const TString title =
      TString("Mean associated #it{p}_{T};") + mass.fTitle + ";" + pt.fTitle + ";" +
      mult.fTitle;
  auto profile = std::make_shared<TProfile3D>(
      name, title, mass.fNbins, mass.fBins.data(), pt.fNbins, pt.fBins.data(),
      mult.fNbins, mult.fBins.data());
  profile->SetDirectory(nullptr);
  return profile;
}

class Profile3DFillHelper
    : public ROOT::Detail::RDF::RActionImpl<Profile3DFillHelper> {
 public:
  using Result_t = TProfile3D;

  Profile3DFillHelper(std::shared_ptr<TProfile3D> profile, unsigned int nslots)
      : fProfile(std::move(profile)), fSlotProfiles(nslots > 0 ? nslots - 1 : 0) {}
  Profile3DFillHelper(Profile3DFillHelper&&) = default;
  Profile3DFillHelper(const Profile3DFillHelper&) = delete;

  std::shared_ptr<Result_t> GetResultPtr() const { return fProfile; }
  void InitTask(TTreeReader*, unsigned int) {}
  void Initialize() {
    fProfile->Reset();
    fProfile->SetDirectory(nullptr);
    for (auto& slot_profile : fSlotProfiles) {
      slot_profile = std::make_unique<TProfile3D>(*fProfile);
      slot_profile->SetDirectory(nullptr);
    }
  }
  void Exec(unsigned int slot, const RVec<float>& mass, const RVec<float>& pt,
            const RVec<float>& mult, const RVec<float>& ref_pt,
            const RVec<float>& weight) {
    auto* profile = slot == 0 ? fProfile.get() : fSlotProfiles[slot - 1].get();
    for (size_t i = 0; i < mass.size(); ++i)
      profile->Fill(mass[i], pt[i], mult[i], ref_pt[i], weight[i]);
  }
  void Finalize() {
    for (const auto& slot_profile : fSlotProfiles)
      fProfile->Add(slot_profile.get());
  }
  std::string GetActionName() { return "TProfile3D fill"; }

 private:
  std::shared_ptr<TProfile3D> fProfile;
  std::vector<std::unique_ptr<TProfile3D>> fSlotProfiles;
};
} // namespace

void JpsiAssoRawPoiEffProfilePt(TString path_input_flow,
                                TString path_input_extra,
                                TString path_output,
                                TString path_config,
                                double bootstrap_probability = 0.5,
                                unsigned int bootstrap_seed = 0,
                                TString only_cut = "",
                                bool unit_efficiency = false) {
  YAML::Node config = YAML::LoadFile(path_config.Data());
  const bool use_efficiency_correction =
      UseEfficiencyCorrection(config, unit_efficiency);
  if (use_efficiency_correction)
    LoadEfficiency(config);

  ROOT::DisableImplicitMT();

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
          .DefineSlot("bootstrap_rand", [bootstrap_seed](unsigned int slot) {
            thread_local TRandom3 rng(bootstrap_seed + slot);
            return static_cast<float>(rng.Uniform(0., 1.));
          })
          .Filter([bootstrap_probability](float r) { return r < bootstrap_probability; },
                  {"bootstrap_rand"}, "event-level bootstrap selection");

  auto df_pairs =
      df_event
          .Define("pair_indices", [](int njpsi, int nref) {
            std::vector<std::pair<int, int>> out;
            out.reserve(njpsi * nref);
            for (int i = 0; i < njpsi; ++i)
              for (int j = 0; j < nref; ++j)
                out.emplace_back(i, j);
            return out;
          }, {"fPT_size", "fPTREF_size"})
          .Define("jpsi_idx", [](const std::vector<std::pair<int, int>>& pairs) {
            RVec<int> out;
            out.reserve(pairs.size());
            for (auto [jpsi, ref] : pairs)
              out.push_back(jpsi);
            return out;
          }, {"pair_indices"})
          .Define("ref_idx", [](const std::vector<std::pair<int, int>>& pairs) {
            RVec<int> out;
            out.reserve(pairs.size());
            for (auto [jpsi, ref] : pairs)
              out.push_back(ref);
            return out;
          }, {"pair_indices"})
          .Define("jpsi_pt", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPT", "jpsi_idx"})
          .Define("e1_pt", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPt1", "jpsi_idx"})
          .Define("e2_pt", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fPt2", "jpsi_idx"})
          .Define("jpsi_eta", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fEta", "jpsi_idx"})
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
          .Define("ref_ITSChi2NCl", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fITSChi2NCl", "ref_idx"})
          .Define("ref_TPCNClsFound", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fTPCNClsFound", "ref_idx"})
          .Define("ref_itsClusterMap", [](const RVec<uint8_t>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fITSClusterMap", "ref_idx"})
          .Define("ref_dcaxy", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fDcaXY", "ref_idx"})
          .Define("ref_dcaz", [](const RVec<float>& v, const RVec<int>& i) { return Take(v, i); },
                  {"fDcaZ", "ref_idx"})
          .Define("nITSCluster", CountSetBitsVec, {"ref_itsClusterMap"})
          .Define("nDcaZ2Dev", NDca2DevVec, {"ref_pt", "ref_dcaz"})
          .Define("nDcaXY2Dev", NDca2DevVec, {"ref_pt", "ref_dcaxy"})
          .Define("pair_fPosZ", [](float value, const RVec<int>& idx) {
            return RVec<float>(idx.size(), value);
          }, {"fPosZ", "jpsi_idx"})
          .Define("pair_NumContribCalib", [](double value, const RVec<int>& idx) {
            return RVec<float>(idx.size(), static_cast<float>(value));
          }, {"NumContribCalib", "jpsi_idx"})
          .Redefine("fPT", "jpsi_pt")
          .Redefine("fPt1", "e1_pt")
          .Redefine("fPt2", "e2_pt")
          .Redefine("fPosZ", "pair_fPosZ")
          .Redefine("NumContribCalib", "pair_NumContribCalib");

  const auto var_mult = ParseStrVar4Hist(config["hist_binning"]["NumContribCalibBinned"]);
  const auto var_mass = ParseStrVar4Hist(config["hist_binning"]["MassJpsiCandidate"]);
  const auto var_pt = ParseStrVar4Hist(config["hist_binning"]["PtJpsiCandidate"]);

  for (auto it = config["cuts"].begin(); it != config["cuts"].end(); ++it) {
    const std::string cut_name = it->first.as<std::string>();
    if (only_cut != "" && cut_name != only_cut.Data())
      continue;
    const std::string mask_name = "pair_mask_" + cut_name;
    auto selected =
        df_pairs.Define(mask_name, it->second.as<std::string>())
            .Redefine("jpsi_mass", "jpsi_mass[" + mask_name + "]")
            .Redefine("jpsi_pt", "jpsi_pt[" + mask_name + "]")
            .Redefine("jpsi_eta", "jpsi_eta[" + mask_name + "]")
            .Redefine("ref_pt", "ref_pt[" + mask_name + "]")
            .Redefine("NumContribCalib", "NumContribCalib[" + mask_name + "]")
            .Define("jpsi_eff_weight",
                    [use_efficiency_correction,
                     setup = EfficiencySetupForCut(cut_name)](const RVec<float>& pt,
                                                              const RVec<float>& eta) {
                      return use_efficiency_correction ? EfficiencyWeightVec(pt, eta, setup)
                                                       : RVec<float>(pt.size(), 1.f);
                    }, {"jpsi_pt", "jpsi_eta"});

    auto profile = MakeProfile(var_mass, var_pt, var_mult, cut_name);
    auto hist_model = GetTH3DM(var_mass, var_pt, var_mult,
                               Form("Histogram_%s", cut_name.c_str()),
                               "J/#psi candidate yield");
    gRResultHandles.push_back(selected.Histo3D(
        hist_model, "jpsi_mass", "jpsi_pt", "NumContribCalib", "jpsi_eff_weight"));
    auto action = selected.Book<RVec<float>, RVec<float>, RVec<float>, RVec<float>,
                                RVec<float>>(
        Profile3DFillHelper(profile, ROOT::GetThreadPoolSize()),
        {"jpsi_mass", "jpsi_pt", "NumContribCalib", "jpsi_pt", "jpsi_eff_weight"});
    gRResultHandles.push_back(action);
  }

  ROOT::RDF::RunGraphs(gRResultHandles);
  TFile output(path_output, "RECREATE");
  RResultWrite(gRResultHandles);
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
  bool unit_efficiency = false;
  if (argc > 1) path_input_flow = argv[1];
  if (argc > 2) path_input_extra = argv[2];
  if (argc > 3) path_output = argv[3];
  if (argc > 4) path_config = argv[4];
  if (argc > 5) bootstrap_probability = atof(argv[5]);
  if (argc > 6) bootstrap_seed = atoi(argv[6]);
  if (argc > 7) only_cut = argv[7];
  if (argc > 8) unit_efficiency = ParseBoolArg(argv[8]);
  JpsiAssoRawPoiEffProfilePt(path_input_flow, path_input_extra, path_output,
                              path_config, bootstrap_probability, bootstrap_seed,
                              only_cut, unit_efficiency);
  return 0;
}
