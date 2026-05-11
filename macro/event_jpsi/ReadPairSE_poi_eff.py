#!/usr/bin/env python3
import sys
import yaml
import numpy as np
import ROOT
from array import array
from functools import reduce

# Global list to collect histogram handles
gRResultHandles = []
gEfficiencyObjects = []

DEFAULT_EFFICIENCY_FILE = "output/jpsi_pt_eta_2d_efficiency.root"
EFFICIENCY_HISTOGRAMS = {
    "default": "default/jpsi_reconstruction_efficiency_pt_eta_default_low_eff_removed",
    "pid1": "pid1/jpsi_reconstruction_efficiency_pt_eta_pid1_low_eff_removed",
    "pid2": "pid2/jpsi_reconstruction_efficiency_pt_eta_pid2_low_eff_removed",
    "pid3": "pid3/jpsi_reconstruction_efficiency_pt_eta_pid3_low_eff_removed",
}

# Declare efficient functional random number generator in C++
ROOT.gInterpreter.Declare(
    """
#include <cstdint>
#include <string>
#include "TH2D.h"

TH2D* gJpsiEfficiencyHistDefault = nullptr;
TH2D* gJpsiEfficiencyHistPid1 = nullptr;
TH2D* gJpsiEfficiencyHistPid2 = nullptr;
TH2D* gJpsiEfficiencyHistPid3 = nullptr;

void SetJpsiEfficiencyHist(const char* setup, TH2D* hist) {
    const std::string setup_name(setup);
    if (setup_name == "pid1") {
        gJpsiEfficiencyHistPid1 = hist;
    } else if (setup_name == "pid2") {
        gJpsiEfficiencyHistPid2 = hist;
    } else if (setup_name == "pid3") {
        gJpsiEfficiencyHistPid3 = hist;
    } else {
        gJpsiEfficiencyHistDefault = hist;
    }
}

TH2D* GetJpsiEfficiencyHistForCut(const char* cut_name) {
    const std::string name(cut_name);
    if (name == "pid1" && gJpsiEfficiencyHistPid1) {
        return gJpsiEfficiencyHistPid1;
    }
    if (name == "pid2" && gJpsiEfficiencyHistPid2) {
        return gJpsiEfficiencyHistPid2;
    }
    if (name == "pid3" && gJpsiEfficiencyHistPid3) {
        return gJpsiEfficiencyHistPid3;
    }
    return gJpsiEfficiencyHistDefault;
}

float GetJpsiEfficiencyForCut(float pt, float eta, const char* cut_name) {
    TH2D* hist = GetJpsiEfficiencyHistForCut(cut_name);
    if (!hist) {
        return 0.f;
    }
    const int pt_bin = hist->GetXaxis()->FindFixBin(pt);
    const int eta_bin = hist->GetYaxis()->FindFixBin(eta);
    if (
        pt_bin < 1 || pt_bin > hist->GetNbinsX() ||
        eta_bin < 1 || eta_bin > hist->GetNbinsY()
    ) {
        return 0.f;
    }
    const double efficiency = hist->GetBinContent(pt_bin, eta_bin);
    return efficiency > 0. ? static_cast<float>(efficiency) : 0.f;
}

float GetJpsiEfficiencyWeightForCut(float pt, float eta, const char* cut_name) {
    const float efficiency = GetJpsiEfficiencyForCut(pt, eta, cut_name);
    return efficiency > 0.f ? 1.f / efficiency : 0.f;
}

// SplitMix64: fast, high-quality, seedable PRNG (public domain)
inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// Generate deterministic random number in [0,1) from r_event, toyIndex, and global seed
inline double functionalRandom(double r_event, uint64_t toyIndex, uint64_t globalSeed = 0) {
    // Clamp r_event to [0, 1 - ε] to avoid overflow
    if (r_event >= 1.0) r_event = 0.999999999999999; // just below 1.0
    uint64_t base = static_cast<uint64_t>(r_event * (1ULL << 53));
    uint64_t s = base ^ (toyIndex + 0x9e3779b97f4a7c15ULL) ^ globalSeed;
    uint64_t z = splitmix64(s);
    return (z >> 11) * 0x1.0p-53; // use top 53 bits for double precision
}

int countSetBits_uint8(uint8_t x) {
    int count = 0;
    while (x) {
        count += x & 1;
        x >>= 1;
    }
    return count;
}

float nDCA2Dev(float pt, float dca) {
    double dev_dca = 0.00179344 + 0.000924651 * pow(abs(pt), -1.4062);
    return abs(dca) / dev_dca;
}

float GetDeltaPhi(float phi1, float phi2) {
    float delta_phi = phi1 - phi2;
    int n = 0;
    while (delta_phi > 1.5 * M_PI && n < 10) {
        n++;
        delta_phi -= 2 * M_PI;
    }
    while (delta_phi < -0.5 * M_PI && n < 10) {
        n++;
        delta_phi += 2 * M_PI;
    }
    if (n >= 10)
        delta_phi = -999.;
    return delta_phi;
}

"""
)


def RResultWrite(result_handles, output_file):
    """Write histograms, handling duplicate names by appending _0, _1, ..."""
    output_file.cd()
    written_names = {}

    for handle in result_handles:
        try:
            h = handle.GetPtr()
        except Exception as e:
            print(f"Warning: Could not retrieve histogram from handle: {e}")
            continue

        if not hasattr(h, "GetName"):
            continue
        name = h.GetName()

        if name in written_names:
            written_names[name] += 1
            new_name = f"{name}_{written_names[name]}"
            h.SetName(new_name)
            print(f"Renaming duplicate histogram to: {new_name}")
        else:
            written_names[name] = 0

        h.Write()


def get_nested_node(data, path, separator="/"):
    if not path:
        return data
    keys = path.split(separator)
    return reduce(lambda current_dict, key: current_dict[key], keys, data)


class StrVar4Hist:
    # ... __init__ 保持不变 ...
    def __init__(self, name, title, unit, nbins, bins):
        self.fName = name
        self.fTitle = title
        self.fUnit = unit
        self.fNbins = nbins

        if len(bins) != nbins + 1 and len(bins) != 2:
            raise ValueError("bins size is not correct")

        if len(bins) == 2:
            start, stop = bins[0], bins[1]
            self.fBins = [start + i * (stop - start) / nbins for i in range(nbins + 1)]
        else:
            self.fBins = list(bins)

    @classmethod
    def from_config(cls, config, str_subnode=None):
        target_config = get_nested_node(config, str_subnode)

        # ... 提取其他参数 ...
        name = target_config["name"]
        title = target_config["title"]
        unit = target_config["unit"]
        nbins = target_config["nbins"]
        raw_bins = target_config["bins"]

        # 【核心修改】使用 numpy 解析含 pi 的字符串列表
        # dtype='f8' 告诉 numpy 将其转换为 float64
        # 它可以自动识别 'pi', '2*pi', 'pi/2' 等标准数学表达式
        processed_bins = []

        for val in raw_bins:
            if isinstance(val, str):
                # 1. 去除所有空格
                clean_val = val.replace(" ", "")

                # 2. 检查是否包含 'pi'
                if "pi" in clean_val:
                    # 移除 'pi'，只保留系数部分
                    # 例如: "-0.5*pi" -> "-0.5*"
                    # 例如: "1.5pi" -> "1.5"
                    coeff_str = clean_val.replace("pi", "")

                    # 处理末尾可能残留的乘号 (如 "-0.5*" -> "-0.5")
                    coeff_str = coeff_str.strip("*")

                    # 如果系数为空（比如用户只写了 "pi"），默认为 1
                    if not coeff_str:
                        coeff = 1.0
                    else:
                        coeff = float(coeff_str)

                    # 3. 直接相乘
                    processed_bins.append(coeff * ROOT.TMath.Pi())
                else:
                    # 如果不包含 pi，尝试直接转为 float
                    processed_bins.append(float(val))
            else:
                # 如果本来就是数字，直接加入
                processed_bins.append(val)
        return cls(name, title, unit, nbins, processed_bins)


def load_config(path_config: str) -> dict:
    """Load and return configuration from YAML file."""
    with open(path_config, "r") as f:
        return yaml.safe_load(f)


def get_efficiency_correction_config(config: dict) -> tuple:
    """Read efficiency correction settings from config."""
    eff_cfg = config.get("efficiency_correction", {})
    path_efficiency = eff_cfg.get("file", DEFAULT_EFFICIENCY_FILE)
    return path_efficiency


def load_jpsi_efficiency_histograms(config: dict) -> tuple:
    """Load all low-efficiency-removed J/psi efficiency maps."""
    path_efficiency = get_efficiency_correction_config(config)
    file_efficiency = ROOT.TFile.Open(path_efficiency)
    if not file_efficiency or file_efficiency.IsZombie():
        raise RuntimeError(f"Cannot open efficiency file: {path_efficiency}")
    loaded = {}
    for setup, hist_path in EFFICIENCY_HISTOGRAMS.items():
        hist = file_efficiency.Get(hist_path)
        if not hist:
            file_efficiency.Close()
            raise RuntimeError(
                f"Cannot find efficiency histogram '{hist_path}' in {path_efficiency}"
            )
        hist = hist.Clone(f"jpsi_efficiency_correction_map_{setup}")
        hist.SetDirectory(0)
        gEfficiencyObjects.append(hist)
        ROOT.SetJpsiEfficiencyHist(setup, hist)
        loaded[setup] = hist_path
    file_efficiency.Close()
    print(f"Using J/psi efficiency correction file: {path_efficiency}")
    for setup, hist_path in loaded.items():
        print(f"  cut setup {setup}: {hist_path}")
    return path_efficiency, loaded


def open_input_root_file(path_input: str) -> ROOT.TFile:
    """Open input ROOT file and return it."""
    file_flowVecd = ROOT.TFile.Open(path_input)
    if not file_flowVecd or file_flowVecd.IsZombie():
        raise RuntimeError(f"Cannot open input file: {path_input}")
    return file_flowVecd


def get_input_tree(file_flowVecd: ROOT.TFile) -> ROOT.TTree:
    """Get the first TTree from the ROOT file."""
    for key in file_flowVecd.GetListOfKeys():
        if key.GetClassName() == "TTree":
            tree = file_flowVecd.Get(key.GetName())
            if tree:
                return tree
    raise RuntimeError("No TTree found in input file")


def parse_histogram_variables(hist_cfg: dict) -> tuple:
    """Parse histogram variables from configuration."""
    var_fPosZ = StrVar4Hist.from_config(hist_cfg, str_subnode="fPosZ")
    var_NumContribCalibBinned = StrVar4Hist.from_config(
        hist_cfg, str_subnode="NumContribCalibBinned"
    )
    var_MassJpsiCandidate = StrVar4Hist.from_config(
        hist_cfg, str_subnode="MassJpsiCandidate"
    )
    var_PtJpsiCandidate = StrVar4Hist.from_config(
        hist_cfg, str_subnode="PtJpsiCandidate"
    )
    var_DeltaEtaUS = StrVar4Hist.from_config(hist_cfg, str_subnode="DeltaEtaUS")
    var_DeltaPhiUS = StrVar4Hist.from_config(hist_cfg, str_subnode="DeltaPhiUS")

    vec_var = [
        var_DeltaEtaUS,
        var_DeltaPhiUS,
        var_fPosZ,
        var_MassJpsiCandidate,
        var_PtJpsiCandidate,
        var_NumContribCalibBinned,
    ]

    vec_var2 = [
        var_fPosZ,
        var_MassJpsiCandidate,
        var_PtJpsiCandidate,
        var_NumContribCalibBinned,
    ]

    return vec_var, vec_var2


def build_rdf_with_variables(tree_input: ROOT.TTree, toy_index: int) -> ROOT.RDataFrame:
    """Build RDataFrame with all defined columns."""
    rdf_base = ROOT.RDataFrame(tree_input)

    rdf_AllVar = (
        rdf_base.Define("DeltaPhi", "GetDeltaPhi(jpsi_phi, ref_phi)")
        .Define("DeltaEta", "jpsi_eta - ref_eta")
        .Define("nITSCluster", "countSetBits_uint8(ref_itsClusterMap)")
        .Define("nDcaZ2Dev", "nDCA2Dev(ref_pt, ref_dcaz)")
        .Define("nDcaXY2Dev", "nDCA2Dev(ref_pt, ref_dcaxy)")
        .Define("randNew", f"functionalRandom(randTag, {toy_index}ULL)")
    )

    return rdf_AllVar


def get_cuts_from_config(config: dict) -> list:
    """Read cuts from config and return list of (cut_name, cut_expr) tuples."""
    cuts_config = config.get("cuts", {})
    if not cuts_config:
        print("Warning: no 'cuts' section in config. Using default inclusive cut.")
        return [("inclusive", "true")]
    return list(cuts_config.items())


def declare_same_jpsi_filter(cut_name: str) -> None:
    """Declare the isSameJpsi filter function for a given cut."""
    ROOT.gInterpreter.Declare(
        f"""
bool isSameJpsi{cut_name}(float tag)
{{
    static float tag_old = -1.;
    bool aaa = (tag_old == tag);
    tag_old = tag;
    return aaa;
}};
"""
    )


def book_histogram_for_cut(
    rdf_filtered: ROOT.RDataFrame,
    vec_var: list,
    cut_name: str,
    result_handles: list,
) -> ROOT.RDataFrame:
    """Book histogram for a given cut and return filtered RDataFrame."""
    hist_name = "_".join(v.fName for v in vec_var) + "_" + cut_name
    axis_titles = (
        ";".join(
            v.fTitle + (" (" + v.fUnit + ")" if v.fUnit else "") for v in vec_var
        )
        + ";"
    )
    full_title = f"{hist_name};{axis_titles}"

    nbins_list = [v.fNbins for v in vec_var]
    edges_list = [v.fBins for v in vec_var]
    edge_arrays = [array("d", edges) for edges in edges_list]

    thnd_model = ROOT.RDF.THnDModel(
        hist_name, full_title, len(vec_var), nbins_list, edge_arrays
    )

    column_names = [v.fName for v in vec_var]
    hist_handle = rdf_filtered.HistoND(thnd_model, column_names, "jpsi_eff_weight")
    result_handles.append(hist_handle)

    return rdf_filtered


def book_histogram_for_cut2(
    rdf_filtered: ROOT.RDataFrame,
    vec_var2: list,
    cut_name: str,
    result_handles: list,
) -> None:
    """Book histogram for non-same Jpsi pairs for a given cut."""
    rdf_filtered2 = rdf_filtered.Filter("!isSameJpsi")

    hist_name2 = "_".join(v.fName for v in vec_var2) + "_" + cut_name
    axis_titles2 = (
        ";".join(
            v.fTitle + (" (" + v.fUnit + ")" if v.fUnit else "") for v in vec_var2
        )
        + ";"
    )
    full_title2 = f"{hist_name2};{axis_titles2}"

    nbins_list2 = [v.fNbins for v in vec_var2]
    edges_list2 = [v.fBins for v in vec_var2]
    edge_arrays2 = [array("d", edges) for edges in edges_list2]

    thnd_model2 = ROOT.RDF.THnDModel(
        hist_name2, full_title2, len(vec_var2), nbins_list2, edge_arrays2
    )
    column_names2 = [v.fName for v in vec_var2]
    hist_handle2 = rdf_filtered2.HistoND(thnd_model2, column_names2, "jpsi_eff_weight")
    result_handles.append(hist_handle2)


def process_cuts_and_book_histograms(
    rdf_AllVar: ROOT.RDataFrame,
    config: dict,
    vec_var: list,
    vec_var2: list,
    result_handles: list,
) -> None:
    """Process all cuts and book histograms."""
    cut_items = get_cuts_from_config(config)

    for cut_name, cut_expr in cut_items:
        print(f"Applying cut '{cut_name}': {cut_expr}")
        declare_same_jpsi_filter(cut_name)

        rdf_filtered = rdf_AllVar.Filter(cut_expr, cut_name).Define(
            "isSameJpsi", f"isSameJpsi{cut_name}(jpsi_mass)"
        ).Define(
            "jpsi_eff",
            f'GetJpsiEfficiencyForCut(jpsi_pt, jpsi_eta, "{cut_name}")',
        ).Define(
            "jpsi_eff_weight",
            f'GetJpsiEfficiencyWeightForCut(jpsi_pt, jpsi_eta, "{cut_name}")',
        )

        book_histogram_for_cut(rdf_filtered, vec_var, cut_name, result_handles)
        book_histogram_for_cut2(rdf_filtered, vec_var2, cut_name, result_handles)


def EventMixingReadingPair(
    path_input_flowVecd: str, path_output: str, path_config: str, toy_index: int
):
    global gRResultHandles
    gRResultHandles.clear()

    # Load configuration
    config = load_config(path_config)
    efficiency_file, efficiency_histograms = load_jpsi_efficiency_histograms(config)

    # Open input file and get tree
    file_flowVecd = open_input_root_file(path_input_flowVecd)
    tree_input = get_input_tree(file_flowVecd)

    # Parse histogram variables
    hist_cfg = config["hist_binning"]
    vec_var, vec_var2 = parse_histogram_variables(hist_cfg)

    # Build RDataFrame with all variables
    rdf_AllVar = build_rdf_with_variables(tree_input, toy_index)

    # Process cuts and book histograms
    process_cuts_and_book_histograms(rdf_AllVar, config, vec_var, vec_var2, gRResultHandles)

    # Execute all graphs
    ROOT.RDF.RunGraphs(gRResultHandles)

    # Write all results
    output_file = ROOT.TFile(path_output, "RECREATE")
    ROOT.TNamed("jpsi_efficiency_correction_file", efficiency_file).Write()
    ROOT.TNamed(
        "jpsi_efficiency_correction_cut_mapping",
        "cut names default, pid1, pid2, pid3 use matching efficiency; all other cut names use default",
    ).Write()
    for setup, hist_path in efficiency_histograms.items():
        ROOT.TNamed(f"jpsi_efficiency_correction_histogram_{setup}", hist_path).Write()
    ROOT.TNamed(
        "jpsi_efficiency_correction_weight",
        "weight = 1 / efficiency(jpsi_pt, jpsi_eta); invalid or removed bins use weight 0",
    ).Write()
    RResultWrite(gRResultHandles, output_file)
    output_file.Close()
    print(f"Output written to: {path_output}")


def main():
    # Parse command-line arguments: now accept optional toy_index (default=0)
    if len(sys.argv) < 4:
        print(
            "Usage: python script.py <input.root> <output.root> <config.yaml> [toy_index]"
        )
        print("Example: python analysis.py data.root results.root my_config.yaml 42")
        sys.exit(1)

    path_input = sys.argv[1]
    path_output = sys.argv[2]
    path_config = sys.argv[3]
    toy_index = int(sys.argv[4]) if len(sys.argv) >= 5 else 0

    EventMixingReadingPair(path_input, path_output, path_config, toy_index)


if __name__ == "__main__":
    main()
