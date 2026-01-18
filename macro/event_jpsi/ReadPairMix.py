#!/usr/bin/env python3
import sys
import yaml
import ROOT
from array import array

# Global list to collect histogram handles
gRResultHandles = []

# Declare efficient functional random number generator in C++
ROOT.gInterpreter.Declare("""
#include <cstdint>
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

#include <vector>
const std::vector<float> vec_mutt = {2.5, 6.5, 9.5, 12.5, 16., 22.5, 25.5, 32., 42., 200.};
const std::vector<float> vec_mutt = {2.5, 6.5, 9.5, 12.5, 16., 22.5, 25.5, 32., 42., 200.};

float MultFromIndex(int index) {
    return vec_mult[index];
}

#include "TRandom3.h"
gRandom->SetSeed(0);
float PosZFromIndex(int index) {
    return (float) (gRandom->Uniform(-10.+2.*index,2.*index));
}

""")

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

        if not hasattr(h, 'GetName'):
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


class StrVar4Hist:
    def __init__(self, name, title, unit, nbins, bins):
        self.fName = name
        self.fTitle = title
        self.fUnit = unit
        self.fNbins = nbins

        if len(bins) != nbins + 1 and len(bins) != 2:
            raise ValueError("bins size is not correct")

        if len(bins) == 2:
            start, stop = bins[0], bins[1]
            self.fBins = [
                start + i * (stop - start) / nbins for i in range(nbins + 1)
            ]
        else:
            self.fBins = list(bins)


def EventMixingReadingPair(path_input_flowVecd: str, path_output: str, path_config: str):
    global gRResultHandles
    gRResultHandles.clear()

    # Open input ROOT file and get TTree
    file_flowVecd = ROOT.TFile.Open(path_input_flowVecd)
    tree_input = None
    for key in file_flowVecd.GetListOfKeys():
        if key.GetClassName() == "TTree":
            tree_input = file_flowVecd.Get(key.GetName())
            break
    if not tree_input:
        raise RuntimeError("No TTree found in input file")

    # Load configuration from specified YAML file
    with open(path_config, "r") as f:
        config = yaml.safe_load(f)

    # Extract binning parameters
    hist_cfg = config["hist_binning"]
    low_edge_deltaPhiToPi = hist_cfg["low_edge_deltaPhiToPi"]
    up_edge_deltaPhiToPi = hist_cfg["up_edge_deltaPhiToPi"]
    n_bins_mass_assoYield = hist_cfg["n_bins_mass_assoYield"]
    delta_eta_cfg = hist_cfg["binning_deltaEta_assoYield"]
    n_bins_deltaEta_assoYield = delta_eta_cfg["n_bins"]
    min_deltaEta_assoYield = delta_eta_cfg["min"]
    max_deltaEta_assoYield = delta_eta_cfg["max"]
    n_bins_deltaPhi_assoYield = hist_cfg["n_bins_deltaPhi_assoYield"]

    # Define histogram axes
    var_fPosZ = StrVar4Hist("fPosZ", "#it{V}_{Z}", "cm", 8, [-10, 10])
    var_NumContribCalibBinned = StrVar4Hist(
        "NumContribCalib", "N_{vtx contrib} Calibrated", "",
        10, [0, 5, 8, 11, 14, 18, 23, 28, 36, 48, 300]
    )
    var_MassJpsiCandidate = StrVar4Hist(
        "jpsi_mass", "M_{ee}", "GeV^{2}/c^{4}",
        n_bins_mass_assoYield, [1.8, 5.4]
    )
    var_PtJpsiCandidate = StrVar4Hist(
        "jpsi_pt", "p_{T}", "GeV/c", 10, [0.0, 5.0]
    )
    var_DeltaEtaUS = StrVar4Hist(
        "DeltaEta", "#Delta#eta_{J/#psi, track}", "",
        n_bins_deltaEta_assoYield, [min_deltaEta_assoYield, max_deltaEta_assoYield]
    )
    var_DeltaPhiUS = StrVar4Hist(
        "DeltaPhi", "#Delta#phi_{J/#psi, track}", "",
        n_bins_deltaPhi_assoYield,
        [low_edge_deltaPhiToPi * ROOT.TMath.Pi(), up_edge_deltaPhiToPi * ROOT.TMath.Pi()]
    )

    vec_var = [
        var_DeltaEtaUS,
        var_DeltaPhiUS,
        var_fPosZ,
        var_MassJpsiCandidate,
        var_PtJpsiCandidate,
        var_NumContribCalibBinned
    ]

    # Build base RDataFrame
    rdf_base = ROOT.RDataFrame(tree_input)

    # Define all needed columns including the random number
    rdf_AllVar = (
        rdf_base.Define("DeltaPhi", "jpsi_phi - ref_phi")
                .Define("DeltaEta", "jpsi_eta - ref_eta")
                .Define("nITSCluster", "countSetBits_uint8(ref_itsClusterMap)")
                .Define("nDcaZ2Dev", "nDCA2Dev(ref_pt, ref_dcaz)")
                .Define("nDcaXY2Dev", "nDCA2Dev(ref_pt, ref_dcaxy)")
                .Define("NumContribCalib", "MultFromIndex(iMult)")
                .Define("fPosZ", "PosZFromIndex(iVtxZ)")
    )

    # Read cuts from config and add random selection to each cut
    cuts_config = config.get("cuts", {})
    if not cuts_config:
        print("Warning: no 'cuts' section in config. Using default inclusive cut.")
        cut_items = [("inclusive", "true")]
    else:
        cut_items = list(cuts_config.items())

    # Book histograms for each cut
    for cut_name, cut_expr in cut_items:
        print(f"Applying cut '{cut_name}': {cut_expr}")
        rdf_filtered = rdf_AllVar.Filter(cut_expr, cut_name)

        hist_name = "_".join(v.fName for v in vec_var) + "_" + cut_name
        axis_titles = ";".join(
            v.fTitle + (" (" + v.fUnit + ")" if v.fUnit else "") for v in vec_var
        )
        full_title = f"{hist_name};{axis_titles}"

        nbins_list = [v.fNbins for v in vec_var]
        edges_list = [v.fBins for v in vec_var]
        edge_arrays = [array('d', edges) for edges in edges_list]

        thnd_model = ROOT.RDF.THnDModel(
            hist_name,
            full_title,
            len(vec_var),
            nbins_list,
            edge_arrays
        )

        column_names = [v.fName for v in vec_var]
        hist_handle = rdf_filtered.HistoND(thnd_model, column_names)
        gRResultHandles.append(hist_handle)

    ROOT.RDF.RunGraphs(gRResultHandles)
    # Write all results
    output_file = ROOT.TFile(path_output, "RECREATE")
    RResultWrite(gRResultHandles, output_file)
    output_file.Close()
    print(f"✅ Output written to: {path_output}")


def main():
    if len(sys.argv) < 4:
        print("Usage: python script.py <input.root> <output.root> <config.yaml>")
        print("Example: python analysis.py data.root results.root my_config.yaml")
        sys.exit(1)

    path_input = sys.argv[1]
    path_output = sys.argv[2]
    path_config = sys.argv[3]

    EventMixingReadingPair(path_input, path_output, path_config)


if __name__ == "__main__":
    main()
