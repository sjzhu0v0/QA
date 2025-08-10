#include "MALICE.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"
#include "TColor.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TStyle.h"

void TrackInfoMC(TString path_input =
                     "/home/szhu/work/alice/analysis/QA/input/track/"
                     "trackInfoMC_24fd4b_550367_hist.root") {
  TFile *file_input = new TFile(path_input);
  // KEY: TH2D     InversePt_fDcaXY;1      1/p_{T}_DCA_{XY}
  // KEY: TH2D     InversePt_fDcaZ;1       1/p_{T}_DCA_{Z}
  // KEY: TH2D     fPt_DeltaPtOverPt;1
  // p_{T}_(p^{rec}_{T}-p^{MC}_{T})/p^{rec}_{T}

  auto InversePt_fDcaXY =
      MRootIO::GetObjectDiectly<TH2D>(file_input, "InversePt_fDcaXY");
  auto InversePt_fDcaZ =
      MRootIO::GetObjectDiectly<TH2D>(file_input, "InversePt_fDcaZ");
  auto fPt_DeltaPtOverPt =
      MRootIO::GetObjectDiectly<TH2D>(file_input, "fPt_DeltaPtOverPt");

  MRootGraphic::StyleCommon();

  TCanvas *c_InvecrsePt_fDca =
      new TCanvas("c_InversePt_fDcaXY", "c_InversePt_fDcaXY", 800, 400);
  c_InvecrsePt_fDca->Divide(2, 1);
  MRootGraphic::StyleHistCommonHist(InversePt_fDcaXY);
  MRootGraphic::StyleHistCommonHist(InversePt_fDcaZ);

  c_InvecrsePt_fDca->cd(1);
  gPad->SetRightMargin(0.15);
  gPad->SetLogz();
  InversePt_fDcaXY->GetXaxis()->SetMaxDigits(2);
  InversePt_fDcaXY->GetYaxis()->SetMaxDigits(2);
  InversePt_fDcaXY->SetTitle("");
  InversePt_fDcaXY->Draw("colz");
  c_InvecrsePt_fDca->cd(2);
  gPad->SetRightMargin(0.15);
  gPad->SetLogz();
  InversePt_fDcaZ->GetXaxis()->SetMaxDigits(2);
  InversePt_fDcaZ->GetYaxis()->SetMaxDigits(2);
  InversePt_fDcaZ->SetTitle("");
  InversePt_fDcaZ->Draw("colz");
  c_InvecrsePt_fDca->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/track/TrackInfoMC_"
      "InversePt_fDcaXY_Z.pdf");

  // gPublisherCanvas =
  //     new MPublisherCanvas("/home/szhu/work/alice/analysis/QA/plot/track/"
  //                          "TrackInfoMC_InversePt_fDcaXY_Z_diff.pdf",
  //                          3, 3);
  // for (int i = 1; i <= InversePt_fDcaXY->GetXaxis()->GetNbins(); i++) {
  //   auto h1 = (TH1D *)InversePt_fDcaXY->ProjectionY(
  //       Form("h1_fPt_DeltaPtOverPt_%d", i), i, i);
  //   gPublisherCanvas->Draw(h1);
  // }

  // gPublisherCanvas->NextPage();

  // for (int i = 1; i <= InversePt_fDcaZ->GetXaxis()->GetNbins(); i++) {
  //   auto h1 = (TH1D *)InversePt_fDcaZ->ProjectionY(
  //       Form("h1_fPt_DeltaPtOverPt_%d", i), i, i);
  //   gPublisherCanvas->Draw(h1);
  // }
  // gPublisherCanvas->finalize();

  StrVar4Hist var_InversePt("InversePt", "1/p_{T}", "GeV/c", 90, {0.2, 7});
  StrVar4Hist var_DCAxy("fDcaXY", "DCA_{XY}", "cm", 100, {-0.04, 0.04});
  StrVar4Hist var_DCAz("fDcaZ", "DCA_{Z}", "cm", 100, {-0.04, 0.04});
  StrVar4Hist var_Pt("fPt", "p_{T}", "GeV/c", 100, {0.1, 6});
  StrVar4Hist var_DeltaPtOverPt("DeltaPtOverPt",
                                "(p^{rec}_{T}-p^{MC}_{T})/p^{rec}_{T}", "", 100,
                                {-0.3, 0.3});

  MIndexHist indexHistInversePt(var_InversePt, 1, 1);
  MIndexHist indexHistDCAxy(var_DCAxy, 1, 1);
  MIndexHist indexHistDCAz(var_DCAz, 1, 1);
  MIndexHist indexHistPt(var_Pt, 1, 1);
  MIndexHist indexHistDeltaPtOverPt(var_DeltaPtOverPt, 1, 1);

  MHist1D h_InversePt_fDcaXY(indexHistInversePt, "SigmaDCAxy_InversePt");
  MHist1D h_InversePt_fDcaZ(indexHistInversePt, "SigmaDCAz_InversePt");

  gPublisherCanvas =
      new MPublisherCanvas("/home/szhu/work/alice/analysis/QA/plot/track/"
                           "TrackInfoMC_InversePt_fDcaXY_Z_diff.pdf",
                           3, 3);

  for (auto indexHistInversePt : indexHistInversePt) {
    gPublisherCanvas->NewPad()->cd();
    auto h1 = (TH1D *)InversePt_fDcaXY->ProjectionY(
        Form("h1_%d", GenerateUID()), indexHistInversePt, indexHistInversePt);
    h1->Fit("gaus", "Q0", "", -0.04, 0.04);
    double sigma = h1->GetFunction("gaus")->GetParameter(2);
    h1->Fit("gaus", "Q", "", -2 * sigma, 2 * sigma);
    h1->Draw();
    MDouble sigmaDCAxy(h1->GetFunction("gaus")->GetParameter(2),
                       //    h1->GetFunction("gaus")->GetParError(2));
                       0.000002 / h1->GetFunction("gaus")->GetParameter(2));
    h_InversePt_fDcaXY.SetBinInfo(sigmaDCAxy);
  }

  gPublisherCanvas->NextPage();

  for (auto indexHistInversePt : indexHistInversePt) {
    gPublisherCanvas->NewPad()->cd();
    auto h1 = (TH1D *)InversePt_fDcaZ->ProjectionY(
        Form("h1_%d", GenerateUID()), indexHistInversePt, indexHistInversePt);
    h1->Fit("gaus", "Q0", "", -0.04, 0.04);
    double sigma = h1->GetFunction("gaus")->GetParameter(2);
    h1->Fit("gaus", "Q", "", -sigma, sigma);
    h1->Draw();
    MDouble sigmaDCAz(h1->GetFunction("gaus")->GetParameter(2),
                      //   h1->GetFunction("gaus")->GetParError(2));
                      0.000002 / h1->GetFunction("gaus")->GetParameter(2));
    h_InversePt_fDcaZ.SetBinInfo(sigmaDCAz);
  }
  gPublisherCanvas->finalize();
  delete gPublisherCanvas;

  TCanvas *c_sigmaDCA_InversePt =
      new TCanvas("c_sigmaDCA_InversePt", "c_sigmaDCA_InversePt", 400, 400);
  MRootGraphic::StyleHistCommonHist(h_InversePt_fDcaXY.fHisto.get());
  MRootGraphic::StyleHistCommonHist(h_InversePt_fDcaZ.fHisto.get());
  c_sigmaDCA_InversePt->cd();
  h_InversePt_fDcaXY.fHisto->SetLineColor(kRed);
  h_InversePt_fDcaZ.fHisto->SetLineColor(kBlue);
  TF1 *f_SigmaDCAxy =
      new TF1("f_SigmaDCAxy", "[0]+[1]*pow(x,[2])", var_InversePt.GetMin(), 6.);
  f_SigmaDCAxy->SetLineColor(kGreen + 2);
  f_SigmaDCAxy->SetParameter(2, 1);
  f_SigmaDCAxy->SetParameter(1, 0.0025);
  h_InversePt_fDcaXY.fHisto->Fit(f_SigmaDCAxy, "", "", 0.3, 6);

  //  0.000924651 + 0.000924651 * pow(x, 1.4062)

  //   TF1 *f_SigmaDCAz =
  //       new TF1("f_SigmaDCAz", "[0]+[1]*pow(x,[2])",
  //       var_InversePt.GetMin(), 5.);
  //   f_SigmaDCAz->SetParameter(2, 1);
  //   f_SigmaDCAz->SetParameter(1, 0.0025);
  //   h_InversePt_fDcaZ.fHisto->Fit(f_SigmaDCAz, "", "", 0.3, 5);
  h_InversePt_fDcaXY.fHisto->GetYaxis()->SetTitle(
      "#sigma of DCA distributions (cm)");
  h_InversePt_fDcaXY.fHisto->GetYaxis()->SetMaxDigits(2);

  h_InversePt_fDcaXY.fHisto->Draw("PMC");
  f_SigmaDCAxy->Draw("same");
  h_InversePt_fDcaZ.fHisto->Draw("PMC SAME");

  TLegend *legend = new TLegend(0.16, 0.69, 0.46, 0.89);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.04);
  legend->AddEntry(h_InversePt_fDcaXY.fHisto.get(), "DCA_{XY}", "L");
  legend->AddEntry(h_InversePt_fDcaZ.fHisto.get(), "DCA_{Z}", "L");
  legend->AddEntry(f_SigmaDCAxy, "DCA fit", "L");
  legend->Draw();

  c_sigmaDCA_InversePt->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/track/TrackInfoMC_"
      "SigmaDCAxy_DCAz_InversePt.pdf");

  TCanvas *c_fPt_DeltaPtOverPt =
      new TCanvas("c_fPt_DeltaPtOverPt", "c_fPt_DeltaPtOverPt", 400, 400);
  c_fPt_DeltaPtOverPt->cd();
  MRootGraphic::StyleHistCommonHist(fPt_DeltaPtOverPt);
  gPad->SetRightMargin(0.15);
  gPad->SetLogz();
  fPt_DeltaPtOverPt->Draw("colz");
  c_fPt_DeltaPtOverPt->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/track/TrackInfoMC_"
      "fPt_DeltaPtOverPt.pdf");

  gPublisherCanvas =
      new MPublisherCanvas("/home/szhu/work/alice/analysis/QA/plot/track/"
                           "TrackInfoMC_fPt_DeltaPtOverPt_diff.pdf",
                           3, 3);
  MHist1D h_Pt_DeltaPtOverPt(indexHistPt, "SigmaDeltaPtOverPt_fPt");

  // for (int i = 1; i <= fPt_DeltaPtOverPt->GetXaxis()->GetNbins(); i++) {
  //   gPublisherCanvas->NewPad()->cd();
  //   auto h1 = (TH1D *)fPt_DeltaPtOverPt->ProjectionY(
  //       Form("h1_fPt_DeltaPtOverPt_%d", i), i, i);
  //   h1->SetTitle(Form("fPt_DeltaPtOverPt: p_{T}:[%.2f,%.2f];",
  //                     fPt_DeltaPtOverPt->GetXaxis()->GetBinLowEdge(i),
  //                     fPt_DeltaPtOverPt->GetXaxis()->GetBinUpEdge(i)));
  //   h1->Fit("gaus", "Q0", "", -0.5, 0.5);
  //   double sigma = h1->GetFunction("gaus")->GetParameter(2);
  //   h1->Fit("gaus", "Q", "", -2 * sigma, 2 * sigma);
  //   gPad->SetLogy();
  //   h1->Draw();
  // }
  for (auto indexHistPt : indexHistPt) {
    gPublisherCanvas->NewPad()->cd();
    auto h1 = (TH1D *)fPt_DeltaPtOverPt->ProjectionY(
        Form("h1_%d", GenerateUID()), indexHistPt, indexHistPt);
    h1->Fit("gaus", "Q0", "", -0.5, 0.5);
    double sigma = h1->GetFunction("gaus")->GetParameter(2);
    h1->Fit("gaus", "Q", "", -2 * sigma, 2 * sigma);
    gPad->SetLogy();
    h1->Draw();
    // MDouble sigmaDeltaPtOverPt(h1->GetFunction("gaus")->GetParameter(2),
    //                            h1->GetFunction("gaus")->GetParError(2));
    MDouble sigmaDeltaPtOverPt(h1->GetRMS(), h1->GetRMSError());
    h_Pt_DeltaPtOverPt.SetBinInfo(sigmaDeltaPtOverPt);
  }

  gPublisherCanvas->finalize();

  TCanvas *c_sigmaDeltaPtOverPt_fPt = new TCanvas(
      "c_sigmaDeltaPtOverPt_fPt", "c_sigmaDeltaPtOverPt_fPt", 400, 400);
  // MRootGraphic::StyleHistCommonHist(h_Pt_DeltaPtOverPt.fHisto.get());
  c_sigmaDeltaPtOverPt_fPt->cd();
  h_Pt_DeltaPtOverPt.fHisto->GetXaxis()->SetRangeUser(0, 6);
  h_Pt_DeltaPtOverPt.fHisto->Draw("C");
  c_sigmaDeltaPtOverPt_fPt->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/track/TrackInfoMC_"
      "SigmaDeltaPtOverPt_fPt.pdf");
}