#!/usr/bin/env python3

# Write a macro with ROOT to create and save efficiency histograms
# import necessary libraries
import ROOT
from ROOT import TFile, TH1D, TH2D, TCanvas, TLegend, TStyle, TLine, TEfficiency

"""
  KEY: TH3D     fPosZ_fPt_fEta_mc;1     mc
  KEY: TH3D     fPosZ_fPt_fEta_reco;1   reco
  KEY: TH2D     fPosZ_fNumContrib_reco;1        reco
  KEY: TH2D     fPosZ_fNumContrib_mc;1  mc
  KEY: TH2D     fPt_fNumContrib_reco;1  reco
  KEY: TH2D     fPt_fNumContrib_mc;1    mc
  KEY: TH2D     fEta_fNumContrib_reco;1 reco
  KEY: TH1D     fNumContrib_mc_mult;1   mc_mult, |#eta|<0.9,1GeV/c<p_{T}<5GeV/c
  KEY: TH1D     fNumContrib_reco_mult;1 reco_mult, |#eta|<0.9,1GeV/c<p_{T}<5GeV/c
"""

# Draw pt efficiency histograms and save them as pdf files
def draw_efficiency_histograms_and_compare(input_file_path1, input_file_path2, tag1, tag2):
    # Efficiency: fPt: p_{T} (GeV/c), fEta: #eta, fPosZ: V_{z} (cm), fNumContrib: N_{contrib}
    # Open the input ROOT files
    input_file1 = TFile.Open(input_file_path1, "READ")
    input_file2 = TFile.Open(input_file_path2, "READ")

    # Retrieve the histograms from the input files
    hNumContrib_mc_mult1 = input_file1.Get("fNumContrib_mc_mult")
    hNumContrib_reco_mult1 = input_file1.Get("fNumContrib_reco_mult")
    hNumContrib_mc_mult2 = input_file2.Get("fNumContrib_mc_mult")
    hNumContrib_reco_mult2 = input_file2.Get("fNumContrib_reco_mult")
    # Create efficiency histograms
    eff_NumContrib_mult1 = TEfficiency(hNumContrib_reco_mult1
, hNumContrib_mc_mult1)
    eff_NumContrib_mult2 = TEfficiency(hNumContrib_reco_mult2, hNumContrib_mc_mult2)
    # Set styles for the efficiency histograms
    eff_NumContrib_mult1.SetLineColor(ROOT.kBlue)
    eff_NumContrib_mult1.SetMarkerColor(ROOT.kBlue)
    eff_NumContrib_mult1.SetMarkerStyle(20)
    eff_NumContrib_mult2.SetLineColor(ROOT.kRed)
    eff_NumContrib_mult2.SetMarkerColor(ROOT.kRed)
    eff_NumContrib_mult2.SetMarkerStyle(21)
    # Create a canvas to draw the efficiency histograms
    canvas = TCanvas("canvas", "Efficiency Comparison", 800, 600)
    canvas.SetGrid()
    # Draw the efficiency histograms
    eff_NumContrib_mult1.SetTitle("Reconstruction Efficiency vs N_{contrib}; N_{contrib}; Efficiency")
    eff_NumContrib_mult1.Draw("AP")
    eff_NumContrib_mult2.Draw("same AP")
    # Create a legend
    legend = TLegend(0.6, 0.2, 0.9, 0.4)
    legend.AddEntry(eff_NumContrib_mult1, tag1, "lp")
    legend.AddEntry(eff_NumContrib_mult2, tag2, "lp")
    legend.Draw()
    # Save the canvas as a PDF file
    canvas.SaveAs("Efficiency_Comparison_NumContrib_mult.pdf")

# Main function to execute the drawing and saving of efficiency histograms
if __name__ == "__main__":
    # Define input file paths and tags
    input_file_path1 = "output_file1.root"
    input_file_path2 = "output_file2.root"
    tag1 = "Dataset 1"
    tag2 = "Dataset 2"
    # Call the function to draw and save efficiency histograms
    draw_efficiency_histograms_and_compare(input_file_path1, input_file_path2, tag1, tag2)


