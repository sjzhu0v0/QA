#!/usr/bin/env python3

# Write a macro with ROOT to create and save efficiency histograms
# import necessary libraries
import ROOT
from ROOT import TFile, TH1D, TH2D, TCanvas, TLegend, TStyle, TLine, TEfficiency

def StyleHistCommon(hist_mb):
    hist_mb.GetXaxis().SetLabelSize(0.03)
    hist_mb.GetXaxis().SetTitleSize(0.04)
    hist_mb.GetXaxis().SetTitleOffset(1.2)

    hist_mb.GetYaxis().SetLabelSize(0.03)
    hist_mb.GetYaxis().SetTitleSize(0.04)
    hist_mb.GetYaxis().SetTitleOffset(1.2)

    hist_mb.GetXaxis().SetNdivisions(505)
    hist_mb.GetXaxis().SetTickLength(0.02)

    hist_mb.GetYaxis().SetNdivisions(505)
    hist_mb.GetYaxis().SetTickLength(0.02)

    hist_mb.GetXaxis().SetLabelOffset(0.01)
    hist_mb.GetYaxis().SetLabelOffset(0.01)

    hist_mb.GetXaxis().SetLabelFont(42)
    hist_mb.GetYaxis().SetLabelFont(42)

    hist_mb.GetXaxis().SetTitleFont(42)
    hist_mb.GetYaxis().SetTitleFont(42)

    hist_mb.GetXaxis().CenterTitle()
    hist_mb.GetYaxis().CenterTitle()

    # 如果需要科学计数法，可以打开下面这行
    # hist_mb.GetYaxis().SetLabelFormat("%g")


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

    # Retrieve the histograms from the input files: fPt_fNumContrib_mc and fPt_fNumContrib_reco
    hPtNumContrib_mc1 = TH2.Cast(input_file1.Get("fPt_fNumContrib_mc"))
    hPtNumContrib_reco1 = TH2.Cast(input_file1.Get("fPt_fNumContrib_reco"))
    hPtNumContrib_mc2 = TH2.Cast(input_file2.Get("fPt_fNumContrib_mc"))
    hPtNumContrib_reco2 = TH2.Cast(input_file2.Get("fPt_fNumContrib_reco"))

    # integrate over fNumContrib to get fPt histograms
    hPt_mc1 = hPtNumContrib_mc1.ProjectionX("hPt_mc1")
    hPt_reco1 = hPtNumContrib_reco1.ProjectionX("hPt_reco1")
    hPt_mc2 = hPtNumContrib_mc2.ProjectionX("hPt_mc2")
    hPt_reco2 = hPtNumContrib_reco2.ProjectionX("hPt_reco2")

    # Retrieve Eta histograms
    hEtaNumContrib_mc1 = TH2.Cast(input_file1.Get("fEta_fNumContrib_mc"))
    hEtaNumContrib_reco1 = TH2.Cast(input_file1.Get("fEta_fNumContrib_reco"))
    hEtaNumContrib_mc2 = TH2.Cast(input_file2.Get("fEta_fNumContrib_mc"))
    hEtaNumContrib_reco2 = TH2.Cast(input_file2.Get("fEta_fNumContrib_reco"))
    hEta_mc1 = hEtaNumContrib_mc1.ProjectionX("hEta_mc1")
    hEta_reco1 = hEtaNumContrib_reco1.ProjectionX("hEta_reco1")
    hEta_mc2 = hEtaNumContrib_mc2.ProjectionX("hEta_mc2")
    hEta_reco2 = hEtaNumContrib_reco2.ProjectionX("hEta_reco2")

    # Retrieve Vz histograms
    hVzNumContrib_mc1 = TH2.Cast(input_file1.Get("fPosZ_fNumContrib_mc"))
    hVzNumContrib_reco1 = TH2.Cast(input_file1.Get("fPosZ_fNumContrib_reco"))
    hVzNumContrib_mc2 = TH2.Cast(input_file2.Get("fPosZ_fNumContrib_mc"))
    hVzNumContrib_reco2 = TH2.Cast(input_file2.Get("fPosZ_fNumContrib_reco"))
    hVz_mc1 = hVzNumContrib_mc1.ProjectionX("hVz_mc1")
    hVz_reco1 = hVzNumContrib_reco1.ProjectionX("hVz_reco1")
    hVz_mc2 = hVzNumContrib_mc2.ProjectionX("hVz_mc2")
    hVz_reco2 = hVzNumContrib_reco2.ProjectionX("hVz_reco2")

    # Create efficiency histograms
    effPt1 = TEfficiency("effPt1", "Efficiency vs p_{T} for " + tag1, hPt_mc1.GetNbinsX(), hPt_mc1.GetXaxis().GetXmin(), hPt_mc1.GetXaxis().GetXmax())
    effPt1.SetPassedEvents(hPt_reco1)
    effPt1.SetTotalEvents(hPt_mc1)
    effPt2 = TEfficiency("effPt2", "Efficiency vs p_{T} for " + tag2, hPt_mc2.GetNbinsX(), hPt_mc2.GetXaxis().GetXmin(), hPt_mc2.GetXaxis().GetXmax())
    effPt2.SetPassedEvents(hPt_reco2)
    effPt2.SetTotalEvents(hPt_mc2)

    effEta1 = TEfficiency("effEta1", "Efficiency vs #eta for " + tag1, hEta_mc1.GetNbinsX(), hEta_mc1.GetXaxis().GetXmin(), hEta_mc1.GetXaxis().GetXmax())
    effEta1.SetPassedEvents(hEta_reco1)
    effEta1.SetTotalEvents(hEta_mc1)
    effEta2 = TEfficiency("effEta2", "Efficiency vs #eta for " + tag2, hEta_mc2.GetNbinsX(), hEta_mc2.GetXaxis().GetXmin(), hEta_mc2.GetXaxis().GetXmax())
    effEta2.SetPassedEvents(hEta_reco2)
    effEta2.SetTotalEvents(hEta_mc2)

    effVz1 = TEfficiency("effVz1", "Efficiency vs V_{z} for " + tag1, hVz_mc1.GetNbinsX(), hVz_mc1.GetXaxis().GetXmin(), hVz_mc1.GetXaxis().GetXmax())
    effVz1.SetPassedEvents(hVz_reco1)
    effVz1.SetTotalEvents(hVz_mc1)
    effVz2 = TEfficiency("effVz2", "Efficiency vs V_{z} for " + tag2, hVz_mc2.GetNbinsX(), hVz_mc2.GetXaxis().GetXmin(), hVz_mc2.GetXaxis().GetXmax())
    effVz2.SetPassedEvents(hVz_reco2)
    effVz2.SetTotalEvents(hVz_mc2)

    # Create a canvas to draw the efficiency histograms and compare
    canvas = TCanvas("canvas", "Efficiency Comparison", 1200, 400)
    canvas.Divide(3, 1)
    # Draw pT efficiency
    canvas.cd(1)
    StyleHistCommon(effPt1)
    StyleHistCommon(effPt2)
    effPt1.SetLineColor(ROOT.kRed)
    effPt1.SetMarkerColor(ROOT.kRed)
    effPt1.Draw("AP")
    effPt2.SetLineColor(ROOT.kBlue)
    effPt2.SetMarkerColor(ROOT.kBlue)
    effPt2.Draw("P SAME")
    legend1 = TLegend(0.6, 0.7, 0.9, 0.9)
    legend1.AddEntry(effPt1, tag1, "lp")
    legend1.AddEntry(effPt2, tag2, "lp")
    legend1.Draw()
    # Draw Eta v
    canvas.cd(2)
    StyleHistCommon(effEta1)
    StyleHistCommon(effEta2)
    effEta1.SetLineColor(ROOT.kRed)
    effEta1.SetMarkerColor(ROOT.kRed)
    effEta1.Draw("AP")
    effEta2.SetLineColor(ROOT.kBlue)
    effEta2.SetMarkerColor(ROOT.kBlue)
    effEta2.Draw("P SAME")

    # Draw Vz efficiency
    canvas.cd(3)
    StyleHistCommon(effVz1)
    StyleHistCommon(effVz2)
    effVz1.SetLineColor(ROOT.kRed)
    effVz1.SetMarkerColor(ROOT.kRed)
    effVz1.Draw("AP")
    effVz2.SetLineColor(ROOT.kBlue)
    effVz2.SetMarkerColor(ROOT.kBlue)
    effVz2.Draw("P SAME")
    # Save the canvas as a pdf file
    output_pdf_path = "Efficiency_Comparison_{}_{}.pdf".format(tag1, tag2)
    canvas.SaveAs(output_pdf_path)
    print("Efficiency comparison saved to:", output_pdf_path)

    # Draw efficiency ratio histograms and save them as pdf files
    canvas_ratio = TCanvas("canvas_ratio", "Efficiency Ratio Comparison", 1200, 400)
    canvas_ratio.Divide(3, 1)
    # Draw pT efficiency ratio
    canvas_ratio.cd(1)
    # create ratio histogram
    effPt_ratio = effPt1.Clone("effPt_ratio")
    effPt_ratio.SetTitle("Efficiency Ratio vs p_{T} ({} / {})".
                        format(tag1, tag2))
    for i in range(1, effPt1.GetNumbinsX() + 1):
        eff1 = effPt1.GetEfficiency(i)
        eff2 = effPt2.GetEfficiency(i)
        if eff2 > 0:
            ratio = eff1 / eff2
            effPt_ratio.SetBinContent(i, ratio)
            # calculate error using binomial error propagation
            err1 = effPt1.GetEfficiencyErrorUp(i) if eff1 >= eff2 else effPt1.GetEfficiencyErrorLow(i)
            err2 = effPt2.GetEfficiencyErrorUp(i) if eff2 >= eff1 else effPt2.GetEfficiencyErrorLow(i)
            ratio_err = ratio * ((err1 / eff1) ** 2 + (err2 / eff2) ** 2) ** 0.5
            effPt_ratio.SetBinError(i, ratio_err)
        else:
            effPt_ratio.SetBinContent(i, 0)
            effPt_ratio.SetBinError(i, 0)
    StyleHistCommon(effPt_ratio)
    effPt_ratio.Draw("AP")

    legent2 = TLegend(0.6, 0.7, 0.9, 0.9)
    legent2.AddEntry(effPt_ratio, "{} / {}".format(tag1, tag2), "lp")
    legent2.Draw()

    # Draw Eta efficiency ratio
    canvas_ratio.cd(2)
    effEta_ratio = effEta1.Clone("effEta_ratio")
    effEta_ratio.SetTitle("Efficiency Ratio vs #eta ({} / {})".
                        format(tag1, tag2))
    for i in range(1, effEta1.GetNumbinsX() + 1):
        eff1 = effEta1.GetEfficiency(i)
        eff2 = effEta2.GetEfficiency(i)
        if eff2 > 0:
            ratio = eff1 / eff2
            effEta_ratio.SetBinContent(i, ratio)
            err1 = effEta1.GetEfficiencyErrorUp(i) if eff1 >= eff2 else effEta1.GetEfficiencyErrorLow(i)
            err2 = effEta2.GetEfficiencyErrorUp(i) if eff2 >= eff1 else effEta2.GetEfficiencyErrorLow(i)
            ratio_err = ratio * ((err1 / eff1) ** 2 + (err2 / eff2) ** 2) ** 0.5
            effEta_ratio.SetBinError(i, ratio_err)
        else:
            effEta_ratio.SetBinContent(i, 0)
            effEta_ratio.SetBinError(i, 0)
    StyleHistCommon(effEta_ratio)
    effEta_ratio.Draw("AP")

    # Draw Vz efficiency ratio
    canvas_ratio.cd(3)
    effVz_ratio = effVz1.Clone("effVz_ratio")
    effVz_ratio.SetTitle("Efficiency Ratio vs V_{z} ({} / {})".
                        format(tag1, tag2))
    for i in range(1, effVz1.GetNumbinsX() +
                        1):
            eff1 = effVz1.GetEfficiency(i)
            eff2 = effVz2.GetEfficiency(i)
            if eff2 > 0:
                ratio = eff1 / eff2
                effVz_ratio.SetBinContent(i, ratio)
                err1 = effVz1.GetEfficiencyErrorUp(i) if eff1 >= eff2 else effVz1.GetEfficiencyErrorLow(i)
                err2 = effVz2.GetEfficiencyErrorUp(i) if eff2 >= eff1 else effVz2.GetEfficiencyErrorLow(i)
                ratio_err = ratio * ((err1 / eff1) ** 2 + (err2 / eff2) ** 2) ** 0.5
                effVz_ratio.SetBinError(i, ratio_err)
            else:
                effVz_ratio.SetBinContent(i, 0)
                effVz_ratio.SetBinError(i, 0)   

    StyleHistCommon(effVz_ratio)
    effVz_ratio.Draw("AP")
    # Save the canvas as a pdf file
    output_pdf_ratio_path = "Efficiency_Ratio_Comparison_{}_{}.pdf".format(tag1, tag2)
    canvas_ratio.SaveAs(output_pdf_ratio_path)
    print("Efficiency ratio comparison saved to:", output_pdf_ratio_path)



# Main function to execute the drawing and saving of efficiency histograms
if __name__ == "__main__":
    # Define input file paths and tags
    input_file_path1 = "/lustre/alice/users/szhu/work/ppJpsiFlow/REF_minbias/output/NonUniform/check_af_cluster2.root"
    input_file_path2 = "/lustre/alice/users/szhu/work/ppJpsiFlow/REF_minbias/output/NonUniform/check_an_cluster2.root"
    tag1 = "LHC24af"
    tag2 = "LHC24an"
    # Call the function to draw and save efficiency histograms
    draw_efficiency_histograms_and_compare(input_file_path1, input_file_path2, tag1, tag2)


