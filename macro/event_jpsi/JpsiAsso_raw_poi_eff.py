#!/usr/bin/env python3
import argparse
import re
from array import array
from functools import reduce

import ROOT
import yaml


DEFAULT_EFFICIENCY_FILE = "output/jpsi_pt_eta_2d_efficiency.root"
EFFICIENCY_HISTOGRAMS = {
    "default": (
        "default/jpsi_reconstruction_efficiency_pt_eta_default_low_eff_removed",
        "jpsi_reconstruction_efficiency_pt_eta_default_low_eff_removed",
    ),
    "pid1": (
        "pid1/jpsi_reconstruction_efficiency_pt_eta_pid1_low_eff_removed",
        "jpsi_reconstruction_efficiency_pt_eta_pid1_low_eff_removed",
    ),
    "pid2": (
        "pid2/jpsi_reconstruction_efficiency_pt_eta_pid2_low_eff_removed",
        "jpsi_reconstruction_efficiency_pt_eta_pid2_low_eff_removed",
    ),
    "pid3": (
        "pid3/jpsi_reconstruction_efficiency_pt_eta_pid3_low_eff_removed",
        "jpsi_reconstruction_efficiency_pt_eta_pid3_low_eff_removed",
    ),
}
EFFICIENCY_SETUPS = frozenset(EFFICIENCY_HISTOGRAMS)

g_result_handles = []
g_efficiency_objects = []


ROOT.gInterpreter.AddIncludePath("include")
ROOT.gInterpreter.Declare(
    r"""
#include <cmath>
#include <cstdint>
#include <string>
#include "ROOT/RVec.hxx"
#include "TH2D.h"

using ROOT::VecOps::RVec;
using ROOT::VecOps::Take;

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

float GetJpsiEfficiencyForCut(float pt, float eta, const char* cut_name) {
    TH2D* hist = gJpsiEfficiencyHistDefault;
    const std::string name(cut_name);
    if (name == "pid1" && gJpsiEfficiencyHistPid1) hist = gJpsiEfficiencyHistPid1;
    if (name == "pid2" && gJpsiEfficiencyHistPid2) hist = gJpsiEfficiencyHistPid2;
    if (name == "pid3" && gJpsiEfficiencyHistPid3) hist = gJpsiEfficiencyHistPid3;
    if (!hist) return 1.f;
    const int pt_bin = hist->GetXaxis()->FindFixBin(pt);
    const int eta_bin = hist->GetYaxis()->FindFixBin(eta);
    if (
        pt_bin < 1 || pt_bin > hist->GetNbinsX() ||
        eta_bin < 1 || eta_bin > hist->GetNbinsY()
    ) {
        return 0.f;
    }
    const float eff = hist->GetBinContent(pt_bin, eta_bin);
    return eff > 0.f ? eff : 0.f;
}

float GetJpsiEfficiencyWeightForCut(float pt, float eta, const char* cut_name) {
    const float eff = GetJpsiEfficiencyForCut(pt, eta, cut_name);
    return eff > 0.f ? 1.f / eff : 0.f;
}

RVec<float> EfficiencyWeightVec(
    const RVec<float>& jpsi_pt,
    const RVec<float>& jpsi_eta,
    const char* cut_name
) {
    RVec<float> out;
    out.reserve(jpsi_pt.size());
    for (size_t i = 0; i < jpsi_pt.size(); ++i)
        out.push_back(GetJpsiEfficiencyWeightForCut(jpsi_pt[i], jpsi_eta[i], cut_name));
    return out;
}

float nDCA2Dev(float pt, float dca) {
    const double dev_dca = 0.00179344 + 0.000924651 * pow(abs(pt), -1.4062);
    return abs(dca) / dev_dca;
}

RVec<float> NDca2DevVec(const RVec<float>& pt, const RVec<float>& dca) {
    RVec<float> out;
    out.reserve(pt.size());
    for (size_t i = 0; i < pt.size(); ++i)
        out.push_back(nDCA2Dev(pt[i], dca[i]));
    return out;
}

int countSetBits_uchar(unsigned char x) {
    int count = 0;
    while (x) {
        count += x & 1;
        x >>= 1;
    }
    return count;
}

RVec<int> CountSetBitsVec(const RVec<unsigned char>& values) {
    RVec<int> out;
    out.reserve(values.size());
    for (auto value : values)
        out.push_back(countSetBits_uchar(value));
    return out;
}
"""
)

def get_nested_node(data, path, separator="/"):
    if not path:
        return data
    return reduce(lambda current, key: current[key], path.split(separator), data)


class StrVar4Hist:
    def __init__(self, name, title, unit, nbins, bins):
        self.fName = name
        self.fTitle = title
        self.fUnit = unit
        self.fNbins = nbins
        if len(bins) not in (2, nbins + 1):
            raise ValueError(f"invalid bin count for {name}")
        if len(bins) == 2:
            start, stop = bins
            self.fBins = [start + i * (stop - start) / nbins for i in range(nbins + 1)]
        else:
            self.fBins = list(bins)

    @classmethod
    def from_config(cls, config, path):
        node = get_nested_node(config, path)
        bins = []
        for value in node["bins"]:
            if isinstance(value, str) and "pi" in value:
                clean = value.replace(" ", "").replace("pi", "").strip("*")
                bins.append((float(clean) if clean else 1.0) * ROOT.TMath.Pi())
            else:
                bins.append(float(value))
        return cls(node["name"], node["title"], node["unit"], node["nbins"], bins)


def default_histogram_config():
    return {
        "fPosZ": {"name": "fPosZ", "title": "#it{V}_{Z}", "unit": "cm", "nbins": 8, "bins": [-10, 10]},
        "NumContribCalibBinned": {
            "name": "NumContribCalib",
            "title": "N_{vtx contrib} Calibrated",
            "unit": "",
            "nbins": 10,
            "bins": [0, 7, 12, 17, 22, 29, 37, 46, 57, 73, 300],
        },
        "MassJpsiCandidate": {
            "name": "jpsi_mass",
            "title": "M_{ee}",
            "unit": "GeV^{2}/c^{4}",
            "nbins": 45,
            "bins": [1.8, 5.4],
        },
        "PtJpsiCandidate": {"name": "jpsi_pt", "title": "p_{T}", "unit": "GeV/c", "nbins": 10, "bins": [0, 5]},
        "DeltaEtaUS": {
            "name": "DeltaEta",
            "title": "#Delta#eta_{J/#psi, track}",
            "unit": "",
            "nbins": 20,
            "bins": [-2, 2],
        },
        "DeltaPhiUS": {
            "name": "AbsDeltaPhi",
            "title": "|#Delta#phi_{J/#psi, track}|",
            "unit": "",
            "nbins": 30,
            "bins": [0, "pi"],
        },
    }


def load_config(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def merged_histogram_config(config):
    merged = default_histogram_config()
    merged.update(config.get("hist_binning", {}))
    return merged


def parse_histogram_variables(config):
    hist_cfg = merged_histogram_config(config)
    pair_vars = [
        StrVar4Hist.from_config(hist_cfg, "DeltaEtaUS"),
        StrVar4Hist.from_config(hist_cfg, "DeltaPhiUS"),
        StrVar4Hist.from_config(hist_cfg, "fPosZ"),
        StrVar4Hist.from_config(hist_cfg, "MassJpsiCandidate"),
        StrVar4Hist.from_config(hist_cfg, "PtJpsiCandidate"),
        StrVar4Hist.from_config(hist_cfg, "NumContribCalibBinned"),
    ]
    single_vars = [
        StrVar4Hist("fPosZSingle", pair_vars[2].fTitle, pair_vars[2].fUnit, pair_vars[2].fNbins, pair_vars[2].fBins),
        StrVar4Hist("jpsi_massSingle", pair_vars[3].fTitle, pair_vars[3].fUnit, 90, [1.8, 5.4]),
        StrVar4Hist("jpsi_ptSingle", pair_vars[4].fTitle, pair_vars[4].fUnit, pair_vars[4].fNbins, pair_vars[4].fBins),
        StrVar4Hist(
            "NumContribCalibSingle",
            pair_vars[5].fTitle,
            pair_vars[5].fUnit,
            pair_vars[5].fNbins,
            pair_vars[5].fBins,
        ),
    ]
    return pair_vars, single_vars


def load_efficiency_histograms(config, unit_efficiency):
    if unit_efficiency:
        return "unit efficiency", {}
    eff_cfg = config.get("efficiency_correction", {})
    path = eff_cfg.get("file", DEFAULT_EFFICIENCY_FILE)
    root_file = ROOT.TFile.Open(path)
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"cannot open efficiency file: {path}")
    loaded = {}
    for setup, hist_paths in EFFICIENCY_HISTOGRAMS.items():
        hist = None
        hist_path = None
        for candidate in hist_paths:
            hist = root_file.Get(candidate)
            if hist:
                hist_path = candidate
                break
        if not hist:
            raise RuntimeError(f"missing efficiency histogram for '{setup}' in {path}")
        clone = hist.Clone(f"jpsi_efficiency_correction_map_{setup}")
        clone.SetDirectory(0)
        g_efficiency_objects.append(clone)
        ROOT.SetJpsiEfficiencyHist(setup, clone)
        loaded[setup] = hist_path
    root_file.Close()
    return path, loaded


def efficiency_setup_for_cut(cut_name):
    return cut_name if cut_name in EFFICIENCY_SETUPS else "default"


def open_tree(path, tree_name):
    chain = ROOT.TChain(tree_name)
    root_file = ROOT.TFile.Open(path)
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"cannot open input file: {path}")
    if root_file.Get(tree_name):
        chain.Add(f"{path}/{tree_name}")
    for key in root_file.GetListOfKeys():
        if key.GetClassName() == "TDirectoryFile" and key.GetName().startswith("DF_"):
            chain.Add(f"{path}/{key.GetName()}/{tree_name}")
    root_file.Close()
    if chain.GetEntries() == 0:
        raise RuntimeError(f"tree '{tree_name}' not found in {path}")
    return chain


def build_base_rdf(path_flow, path_extra, allow_test_fallback, bootstrap_probability):
    flow_tree = open_tree(path_flow, "O2dqflowvecd")
    if path_extra:
        extra_tree = open_tree(path_extra, "ExtraInfo")
        flow_tree.AddFriend(extra_tree)
        has_extra = True
    elif allow_test_fallback:
        has_extra = False
    else:
        raise RuntimeError(
            "missing --extra-info; accurate rct-style processing requires the ExtraInfo friend tree "
            "from the raw-data preparation step"
        )

    rdf = ROOT.RDataFrame(flow_tree)
    rdf = (
        rdf.Define("isntSameBunchPileup", "bool((fSelection >> 36) & 1)")
        .Define("isntITSROFrameBorder", "bool((fSelection >> 34) & 1)")
        .Define("isntTimeFrameBorder", "bool((fSelection >> 35) & 1)")
        .Define("isTriggerTVX", "bool((fSelection >> 32) & 1)")
    )
    if not has_extra:
        rdf = rdf.Define("isCBT", "true").Define("NumContribCalib", "float(fNumContrib)")
    return (
        rdf.Filter(
            "isTriggerTVX && isntITSROFrameBorder && isntTimeFrameBorder && isntSameBunchPileup && isCBT"
        )
        .Define("bootstrap_rand", "gRandom->Uniform(0., 1.)")
        .Define("Cut_BS", f"bootstrap_rand < {bootstrap_probability}")
        .Filter("Cut_BS", "event-level bootstrap selection")
    )


def define_pair_columns(rdf):
    return (
        rdf.Define("pair_indices", "ROOT::VecOps::Combinations(fPT, fPTREF)")
        .Define("jpsi_idx", "pair_indices[0]")
        .Define("ref_idx", "pair_indices[1]")
        .Define("jpsi_pt", "Take(fPT, jpsi_idx)")
        .Define("jpsi_eta", "Take(fEta, jpsi_idx)")
        .Define("jpsi_phi", "Take(fPhi, jpsi_idx)")
        .Define("jpsi_mass", "Take(fMass, jpsi_idx)")
        .Define("jpsi_sign", "Take(fSign, jpsi_idx)")
        .Define("e1_pt", "Take(fPt1, jpsi_idx)")
        .Define("e1_eta", "Take(fEta1, jpsi_idx)")
        .Define("e1_phi", "Take(fPhi1, jpsi_idx)")
        .Define("e1_sign", "Take(fSign1, jpsi_idx)")
        .Define("e1_ITSChi2NCl", "Take(fITSChi2NCl1, jpsi_idx)")
        .Define("e1_TPCNClsCR", "Take(fTPCNClsCR1, jpsi_idx)")
        .Define("e1_TPCNClsFound", "Take(fTPCNClsFound1, jpsi_idx)")
        .Define("e1_TPCChi2NCl", "Take(fTPCChi2NCl1, jpsi_idx)")
        .Define("e1_TPCSignal", "Take(fTPCSignal1, jpsi_idx)")
        .Define("e1_nsig_el", "Take(fTPCNSigmaEl1, jpsi_idx)")
        .Define("e1_nsig_pi", "Take(fTPCNSigmaPi1, jpsi_idx)")
        .Define("e1_nsig_pr", "Take(fTPCNSigmaPr1, jpsi_idx)")
        .Define("e2_pt", "Take(fPt2, jpsi_idx)")
        .Define("e2_eta", "Take(fEta2, jpsi_idx)")
        .Define("e2_phi", "Take(fPhi2, jpsi_idx)")
        .Define("e2_sign", "Take(fSign2, jpsi_idx)")
        .Define("e2_ITSChi2NCl", "Take(fITSChi2NCl2, jpsi_idx)")
        .Define("e2_TPCNClsCR", "Take(fTPCNClsCR2, jpsi_idx)")
        .Define("e2_TPCNClsFound", "Take(fTPCNClsFound2, jpsi_idx)")
        .Define("e2_TPCChi2NCl", "Take(fTPCChi2NCl2, jpsi_idx)")
        .Define("e2_TPCSignal", "Take(fTPCSignal2, jpsi_idx)")
        .Define("e2_nsig_el", "Take(fTPCNSigmaEl2, jpsi_idx)")
        .Define("e2_nsig_pi", "Take(fTPCNSigmaPi2, jpsi_idx)")
        .Define("e2_nsig_pr", "Take(fTPCNSigmaPr2, jpsi_idx)")
        .Define("ref_pt", "Take(fPTREF, ref_idx)")
        .Define("ref_eta", "Take(fEtaREF, ref_idx)")
        .Define("ref_phi", "Take(fPhiREF, ref_idx)")
        .Define("ref_ITSChi2NCl", "Take(fITSChi2NCl, ref_idx)")
        .Define("ref_TPCNClsCR", "Take(fTPCNClsCR, ref_idx)")
        .Define("ref_TPCNClsFound", "Take(fTPCNClsFound, ref_idx)")
        .Define("ref_TPCChi2NCl", "Take(fTPCChi2NCl, ref_idx)")
        .Define("ref_TPCSignal", "Take(fTPCSignal, ref_idx)")
        .Define("ref_itsClusterMap", "Take(fITSClusterMap, ref_idx)")
        .Define("ref_nsig_el", "Take(fTPCNSigmaEl, ref_idx)")
        .Define("ref_nsig_pi", "Take(fTPCNSigmaPi, ref_idx)")
        .Define("ref_nsig_pr", "Take(fTPCNSigmaPr, ref_idx)")
        .Define("ref_dcaxy", "Take(fDcaXY, ref_idx)")
        .Define("ref_dcaz", "Take(fDcaZ, ref_idx)")
        .Define("DeltaEta", "jpsi_eta - ref_eta")
        .Define("DeltaPhiRaw", "jpsi_phi - ref_phi")
        .Define(
            "DeltaPhiWrapped",
            "ROOT::VecOps::Where("
            "DeltaPhiRaw > float(M_PI), "
            "DeltaPhiRaw - 2.f * float(M_PI), "
            "ROOT::VecOps::Where(DeltaPhiRaw < -float(M_PI), DeltaPhiRaw + 2.f * float(M_PI), DeltaPhiRaw)"
            ")",
        )
        .Define("AbsDeltaPhi", "ROOT::VecOps::Where(DeltaPhiWrapped < 0.f, -DeltaPhiWrapped, DeltaPhiWrapped)")
        .Define("pair_fPosZ", "ROOT::VecOps::RVec<float>(jpsi_idx.size(), fPosZ)")
        .Define("pair_NumContribCalib", "ROOT::VecOps::RVec<float>(jpsi_idx.size(), float(NumContribCalib))")
        .Define("pair_randTag", "ROOT::VecOps::RVec<float>(jpsi_idx.size(), bootstrap_rand)")
        .Define("pair_randNew", "ROOT::VecOps::RVec<float>(jpsi_idx.size(), bootstrap_rand)")
        .Define("nITSCluster", "CountSetBitsVec(ref_itsClusterMap)")
        .Define("nDcaZ2Dev", "NDca2DevVec(ref_pt, ref_dcaz)")
        .Define("nDcaXY2Dev", "NDca2DevVec(ref_pt, ref_dcaxy)")
        .Redefine("fPosZ", "pair_fPosZ")
        .Redefine("NumContribCalib", "pair_NumContribCalib")
        .Define("randNew", "pair_randNew")
    )


EVENT_CUT_COLUMNS = {
    "fPosZ": "pair_fPosZ",
    "NumContribCalib": "pair_NumContribCalib",
    "randTag": "pair_randTag",
    "randNew": "pair_randNew",
}


def vectorize_cut_expression(expr):
    expr = expr.strip()
    if not expr or expr == "true":
        return "ROOT::VecOps::RVec<char>(jpsi_idx.size(), true)"
    for old, new in EVENT_CUT_COLUMNS.items():
        expr = re.sub(rf"\b{old}\b", new, expr)
    return expr


def book_histogram(rdf, variables, cut_name, weight_name=None):
    hist_name = "_".join(variable.fName for variable in variables) + "_" + cut_name
    titles = ";".join(
        variable.fTitle + (f" ({variable.fUnit})" if variable.fUnit else "")
        for variable in variables
    )
    model = ROOT.RDF.THnDModel(
        hist_name,
        f"{hist_name};{titles};",
        len(variables),
        [variable.fNbins for variable in variables],
        [array("d", variable.fBins) for variable in variables],
    )
    columns = [variable.fName for variable in variables]
    if weight_name:
        g_result_handles.append(rdf.HistoND(model, columns, weight_name))
    else:
        g_result_handles.append(rdf.HistoND(model, columns))


def book_cut_histograms(rdf, cut_name, cut_expr, pair_vars, single_vars, unit_efficiency):
    mask_name = f"pair_mask_{cut_name}"
    selected = rdf.Define(mask_name, vectorize_cut_expression(cut_expr))
    for variable in pair_vars:
        selected = selected.Redefine(variable.fName, f"{variable.fName}[{mask_name}]")
    weight_expr = (
        "ROOT::VecOps::RVec<float>(jpsi_pt.size(), 1.f)"
        if unit_efficiency
        else f'::EfficiencyWeightVec(jpsi_pt, jpsi_eta, "{efficiency_setup_for_cut(cut_name)}")'
    )
    selected = (
        selected.Redefine("jpsi_idx", f"jpsi_idx[{mask_name}]")
        .Redefine("ref_idx", f"ref_idx[{mask_name}]")
        .Redefine("jpsi_pt", f"jpsi_pt[{mask_name}]")
        .Redefine("jpsi_eta", f"jpsi_eta[{mask_name}]")
        .Define("jpsi_eff_weight", weight_expr)
    )
    book_histogram(selected, pair_vars, cut_name, None if unit_efficiency else "jpsi_eff_weight")

    singles = (
        selected.Define(f"single_mask_{cut_name}", "ref_idx == 0")
        .Define("fPosZSingle", f"pair_fPosZ[single_mask_{cut_name}]")
        .Define("jpsi_massSingle", f"jpsi_mass[single_mask_{cut_name}]")
        .Define("jpsi_ptSingle", f"jpsi_pt[single_mask_{cut_name}]")
        .Define("NumContribCalibSingle", f"pair_NumContribCalib[single_mask_{cut_name}]")
        .Define("jpsi_eff_weight_single", f"jpsi_eff_weight[single_mask_{cut_name}]")
    )
    book_histogram(singles, single_vars, cut_name, None if unit_efficiency else "jpsi_eff_weight_single")


def write_results(path_output):
    output = ROOT.TFile(path_output, "RECREATE")
    for handle in g_result_handles:
        handle.GetPtr().Write()
    output.Close()


def run(
    path_flow,
    path_output,
    path_config,
    path_extra=None,
    unit_efficiency=False,
    cuts=None,
    allow_test_fallback=False,
    bootstrap_probability=0.5,
    bootstrap_seed=0,
):
    g_result_handles.clear()
    ROOT.gRandom.SetSeed(bootstrap_seed)
    config = load_config(path_config)
    efficiency_file, efficiency_histograms = load_efficiency_histograms(config, unit_efficiency)
    pair_vars, single_vars = parse_histogram_variables(config)
    rdf = define_pair_columns(
        build_base_rdf(path_flow, path_extra, allow_test_fallback, bootstrap_probability)
    )
    cut_items = list((config.get("cuts") or {"default": "true"}).items())
    if cuts:
        wanted = set(cuts)
        cut_items = [(name, expr) for name, expr in cut_items if name in wanted]
    for cut_name, cut_expr in cut_items:
        print(f"Applying cut '{cut_name}': {cut_expr}")
        book_cut_histograms(rdf, cut_name, cut_expr, pair_vars, single_vars, unit_efficiency)
    ROOT.RDF.RunGraphs(g_result_handles)
    write_results(path_output)


def main():
    parser = argparse.ArgumentParser(description="Direct raw-input J/psi association histogrammer")
    parser.add_argument("flow_input", help="ROOT file containing O2dqflowvecd")
    parser.add_argument("output", help="output ROOT file")
    parser.add_argument("config", help="YAML config with cuts, histogram binning, and optional efficiency settings")
    parser.add_argument("--extra-info", help="ROOT file containing ExtraInfo, as used by JpsiAssoPair_rct.cpp")
    parser.add_argument(
        "--bootstrap-probability",
        type=float,
        default=0.5,
        help="event-level bootstrap acceptance probability; default: 0.5",
    )
    parser.add_argument(
        "--bootstrap-seed",
        type=int,
        default=0,
        help="seed for the event-level bootstrap random generator; default: 0",
    )
    parser.add_argument("--unit-efficiency", action="store_true", help="use weight 1 instead of loading efficiency maps")
    parser.add_argument("--cuts", nargs="+", help="only run the named cuts; useful for low-memory tests")
    parser.add_argument(
        "--allow-test-fallback",
        action="store_true",
        help="allow validation without ExtraInfo by using isCBT=true and NumContribCalib=fNumContrib",
    )
    args = parser.parse_args()
    run(
        args.flow_input,
        args.output,
        args.config,
        path_extra=args.extra_info,
        unit_efficiency=args.unit_efficiency,
        cuts=args.cuts,
        allow_test_fallback=args.allow_test_fallback,
        bootstrap_probability=args.bootstrap_probability,
        bootstrap_seed=args.bootstrap_seed,
    )


if __name__ == "__main__":
    main()
