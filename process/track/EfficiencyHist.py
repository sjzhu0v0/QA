#!/usr/bin/env python3

import ROOT
from ROOT import TFile, TH1D, TCanvas, TLegend, TEfficiency

def StyleHistCommon(hist):
    hist.GetXaxis().SetLabelSize(0.03)
    hist.GetXaxis().SetTitleSize(0.04)
    hist.GetXaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetLabelSize(0.03)
    hist.GetYaxis().SetTitleSize(0.04)
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetXaxis().SetNdivisions(505)
    hist.GetYaxis().SetNdivisions(505)
    hist.GetXaxis().CenterTitle()
    hist.GetYaxis().CenterTitle()

def style_eff(eff, color):
    eff.SetLineColor(color)
    eff.SetMarkerColor(color)
    eff.SetMarkerStyle(20)

def draw_efficiency_histograms_and_compare(f1, f2, tag1, tag2):

    file1 = TFile.Open(f1)
    file2 = TFile.Open(f2)

    def proj(file, name):
        return file.Get(name).ProjectionX()

    hPt_mc1   = proj(file1, "fPt_fNumContrib_mc")
    hPt_reco1 = proj(file1, "fPt_fNumContrib_reco")
    hPt_mc2   = proj(file2, "fPt_fNumContrib_mc")
    hPt_reco2 = proj(file2, "fPt_fNumContrib_reco")

    hEta_mc1   = proj(file1, "fEta_fNumContrib_mc")
    hEta_reco1 = proj(file1, "fEta_fNumContrib_reco")
    hEta_mc2   = proj(file2, "fEta_fNumContrib_mc")
    hEta_reco2 = proj(file2, "fEta_fNumContrib_reco")

    hVz_mc1   = proj(file1, "fPosZ_fNumContrib_mc")
    hVz_reco1 = proj(file1, "fPosZ_fNumContrib_reco")
    hVz_mc2   = proj(file2, "fPosZ_fNumContrib_mc")
    hVz_reco2 = proj(file2, "fPosZ_fNumContrib_reco")
    
    printf("Creating TEfficiency objects...")
    printf(f"  {tag1}...")
    printf(f"   effPt1")
    effPt1  = TEfficiency(hPt_reco1,  hPt_mc1)
    printf(f"   effPt2")
    effPt2  = TEfficiency(hPt_reco2,  hPt_mc2)
    printf(f"   effEta1")
    effEta1 = TEfficiency(hEta_reco1, hEta_mc1)
    printf(f"   effEta2")
    effEta2 = TEfficiency(hEta_reco2, hEta_mc2)
    printf(f"   effVz1")
    effVz1  = TEfficiency(hVz_reco1,  hVz_mc1)
    printf(f"   effVz2")
    effVz2  = TEfficiency(hVz_reco2,  hVz_mc2)

    style_eff(effPt1,  ROOT.kRed)
    style_eff(effPt2,  ROOT.kBlue)
    style_eff(effEta1, ROOT.kRed)
    style_eff(effEta2, ROOT.kBlue)
    style_eff(effVz1,  ROOT.kRed)
    style_eff(effVz2,  ROOT.kBlue)

    c = TCanvas("c", "", 1200, 400)
    c.Divide(3, 1)

    for i, (e1, e2, xtitle) in enumerate([
        (effPt1,  effPt2,  "p_{T} (GeV/c)"),
        (effEta1, effEta2, "#eta"),
        (effVz1,  effVz2,  "V_{z} (cm)")
    ], start=1):
        c.cd(i)
        e1.Draw("AP")
        ROOT.gPad.Update()
        g = e1.GetPaintedGraph()
        g.GetXaxis().SetTitle(xtitle)
        g.GetYaxis().SetTitle("Efficiency")
        StyleHistCommon(g)
        e2.Draw("P SAME")
        if i == 1:
            leg = TLegend(0.6,0.7,0.9,0.9)
            leg.AddEntry(e1, tag1, "lp")
            leg.AddEntry(e2, tag2, "lp")
            leg.Draw()

    c.SaveAs(f"Efficiency_Comparison_{tag1}_{tag2}.pdf")

    c2 = TCanvas("c2", "", 1200, 400)
    c2.Divide(3, 1)

    def draw_ratio(eff1, eff2, ref_hist, title):
        h = TH1D(title, title,
                 ref_hist.GetNbinsX(),
                 ref_hist.GetXaxis().GetXmin(),
                 ref_hist.GetXaxis().GetXmax())
        for i in range(1, h.GetNbinsX()+1):
            e1 = eff1.GetEfficiency(i)
            e2 = eff2.GetEfficiency(i)
            if e2 > 0 and e1 > 0:
                r = e1 / e2
                err = r * ((eff1.GetEfficiencyErrorUp(i)/e1)**2 +
                           (eff2.GetEfficiencyErrorUp(i)/e2)**2)**0.5
                h.SetBinContent(i, r)
                h.SetBinError(i, err)
        h.GetYaxis().SetTitle("Ratio")
        StyleHistCommon(h)
        h.Draw("E")
        return h

    c2.cd(1); draw_ratio(effPt1,  effPt2,  hPt_mc1,  "p_{T}")
    c2.cd(2); draw_ratio(effEta1, effEta2, hEta_mc1, "#eta")
    c2.cd(3); draw_ratio(effVz1,  effVz2,  hVz_mc1,  "V_{z}")

    c2.SaveAs(f"Efficiency_Ratio_Comparison_{tag1}_{tag2}.pdf")

if __name__ == "__main__":
    draw_efficiency_histograms_and_compare(
        "/lustre/alice/users/szhu/work/ppJpsiFlow/REF_minbias/output/Efficiency/merge/af_cluster3.root",
        "/lustre/alice/users/szhu/work/ppJpsiFlow/REF_minbias/output/Efficiency/merge/an_cluster3.root",
        "LHC24af",
        "LHC24an"
    )

