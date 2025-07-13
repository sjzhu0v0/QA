#include "MALICE.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "RooAddPdf.h"
#include "RooCrystalBall.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "TApplication.h"
#include "TColor.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TStyle.h"
#include "math.h"

void pt1_2() {
  gROOT->SetBatch(true);
  TFile *file_input = new TFile("/home/szhu/work/alice/analysis/QA/input/jpsi/"
                                "JpsiQA_LHC22pass4_dqfilter.root");
  auto fPosZ_fMass_fPT_NumContribCalib_Binned = MRootIO::GetObjectDiectly<THnD>(
      file_input, "fPosZ_fMass_fPT_NumContribCalib_Binned");

  MHnTool hn_tool(fPosZ_fMass_fPT_NumContribCalib_Binned);
  hn_tool.PrintAllAxis();
  hn_tool.Rebin(0, 200);
  hn_tool.Rebin(3, 2);
  hn_tool.PrintAllAxis();

  auto mass = hn_tool.Project(1, {0, 0, 0});

  MRootGraphic::StyleCommon();
  gPublisherCanvas = new MPublisherCanvas("./mass_pt1_2.pdf", 1, 1, 600, 600);
  int i_pt = 2;
  auto h1_lowMult = hn_tool.Project(1, {0, i_pt, 1});
  auto h1_highMult = hn_tool.Project(1, {0, i_pt, 5});
  h1_lowMult->SetLineColor(kRed);
  h1_highMult->SetLineColor(kBlue);
  double scale_lowMult2highMult =
      h1_highMult->Integral(h1_highMult->FindBin(4.),
                            h1_highMult->FindBin(4.5)) /
      h1_lowMult->Integral(h1_lowMult->FindBin(4.), h1_highMult->FindBin(4.5));

  h1_highMult->Scale(1. / scale_lowMult2highMult);

  double max = std::max(h1_lowMult->GetMaximum(), h1_highMult->GetMaximum());
  h1_lowMult->GetYaxis()->SetRangeUser(0.1, max * 1.2);
  h1_lowMult->SetMaximum(max * 1.2);
  h1_lowMult->GetXaxis()->SetRangeUser(1., 5.);

  auto h1_delta = (TH1D *)h1_lowMult->Clone(Form("h1_delta_pT%d", i_pt));
  h1_delta->Add(h1_highMult, -1);
  h1_delta->SetLineColor(kBlack);
  gPublisherCanvas->Draw(h1_lowMult)->DrawSame(h1_highMult)->DrawSame(h1_delta);
  /////////////////////////////////
  auto h1_highMult2 = hn_tool.Project(1, {0, i_pt, 2});
  double scale_lowMult2highMult2 =
      h1_highMult2->Integral(h1_highMult2->FindBin(4.),
                             h1_highMult->FindBin(4.5)) /
      h1_lowMult->Integral(h1_lowMult->FindBin(4.), h1_highMult->FindBin(4.5));
  h1_highMult2->Scale(1. / scale_lowMult2highMult2);
  auto h1_delta2 = (TH1D *)h1_lowMult->Clone(Form("h1_delta2_pT%d2", i_pt));
  h1_delta2->Add(h1_highMult2, -1);
  h1_delta2->SetLineColor(kGreen);
  gPublisherCanvas->Draw(h1_lowMult)
      ->DrawSame(h1_highMult2)
      ->DrawSame(h1_delta2);
  ///////////////////////////
  double ratio_deltaConvert = h1_delta->GetMaximum() / h1_delta2->GetMaximum();
  auto h_delta2_scaled =
      (TH1D *)h1_delta2->Clone(Form("h1_delta2_scaled_pT%d", i_pt));
  h_delta2_scaled->Scale(ratio_deltaConvert);
  h1_delta->Draw();
  h_delta2_scaled->Draw("same");
  auto h_delta3 = (TH1D *)h1_delta->Clone(Form("h1_delta3_pT%d", i_pt));
  h_delta3->Add(h_delta2_scaled, -1);
  h_delta3->SetLineColor(kMagenta);
  gPublisherCanvas->Draw(h1_delta)->DrawSame(h1_delta2)->DrawSame(h_delta3);
  ///////////////////////////////
  gPublisherCanvas->NewPad()->cd();
  auto lFumula = [](double *x, double *par) {
    double x2 = x[0] - 0.9;
    return par[0] * x2 * x2 / (1 + pow(x2, par[1]));
  };
  TF1 *f1 = new TF1("f1", lFumula, 1., 5., 2);
  f1->SetParameters(1074.42, 6.37887);
  h_delta3->Fit(f1, "", "same", 1., 5.);

  gPublisherCanvas->NewPad()->cd();
  TH1D *hist_jpsi = h1_delta; // 替换为您的直方图数据
  RooRealVar x("x", "Invariant mass [GeV]", 1.6,
               4.5); // 调整范围匹配 J/ψ 质量区间

  // 将 TH1D 转换为 RooFit 可用的数据格式
  RooDataHist data("data", "J/psi ee decay", x, hist_jpsi);

  // --- 信号部分：Crystal Ball 函数 ---
  RooRealVar mean("mean", "Mean of CB", 3.0618e+00, 3.0,
                  3.2); // J/ψ 质量 ~3.097 GeV
  RooRealVar sigma("sigma", "Sigma of CB", 5.0255e-02, 0.01, 0.1);
  RooRealVar alpha("alpha", "Alpha (tail parameter)", 2.3231e-01, 0.1, 3.0);
  RooRealVar n("n", "Power of tail", 6.1057e+00, 0.1, 10.0);
  RooCrystalBall cb("cb", "Crystal Ball PDF", x, mean, sigma, alpha, n);

  // --- 背景部分：自定义函数 ---
  RooRealVar p0("p0", "Background amplitude", 1193.11, 10, 10000);
  RooRealVar p1("p1", "Background exponent", 4.45199, 4, 10);
  // RooRealVar p2("p2", "Background numerator", 2, 0, 4);
  RooGenericPdf bkgPDF("bkgPDF", "Background PDF",
                       "p0 * (x - 0.9)^2 / (1 + (x-0.9)^p1)",
                       RooArgList(x, p0, p1));
  RooRealVar nsig("nsig", "Number of signal events", 3.7184e+04, 1e3, 2e5);
  RooRealVar nbkg("nbkg", "Number of background events", 1.4669e+04, 1e3, 1e5);
  RooAddPdf model("model", "Total PDF", RooArgList(cb, bkgPDF),
                  RooArgList(nsig, nbkg));

  RooArgSet *params = model.getParameters(x);

  // // --- 执行拟合 ---
  // // --- 执行拟合 ---
  // // model.fitTo(data, RooFit::Minimizer("Minuit2"), RooFit::PrintLevel(1));
  model.chi2FitTo(data, RooFit::SumW2Error(true));
  for (RooAbsArg *arg : *params) {
    RooRealVar *var = dynamic_cast<RooRealVar *>(arg);
    if (var) {
      var->removeMin(); // 去掉最小值限制
      var->removeMax(); // 去掉最大值限制
    }
  }
  model.chi2FitTo(data, RooFit::SumW2Error(true));

  // --- 绘制结果 ---
  RooPlot *frame = x.frame();
  data.plotOn(frame);
  model.plotOn(frame);
  model.plotOn(frame, RooFit::Components(cb), RooFit::LineColor(kRed),
               RooFit::LineStyle(kDashed));
  model.plotOn(frame, RooFit::Components(bkgPDF), RooFit::LineColor(kGreen),
               RooFit::LineStyle(kDashed));
  frame->Draw();

  // gPublisherCanvas->NewPad()->cd();
  TF1 *crystalball = new TF1("crystalball", "crystalball", 1.55, 5.);
  // crystalball->SetParameters(1., 3.06213, 0.0491487, 0.234623, 5.12185);
  crystalball->SetParameters(5030., 3.06213, 0.0491487, 0.234623, 5.12185);
  crystalball->SetParNames("alpha", "mean", "n", "sigma");
  gPublisherCanvas->Draw(crystalball);

  auto hist_copy = (TH1D *)h1_lowMult->Clone("hist_copy");
  for (int j = 1; j <= hist_copy->GetNbinsX(); ++j) {
    double content = hist_copy->GetBinContent(j);
    double crystalball_value =
        crystalball->Eval(hist_copy->GetBinCenter(j)) * 0.9;
    hist_copy->SetBinContent(j, content - crystalball_value);
  }
  for (int i = 51; i <= 53; i++) {
    hist_copy->SetBinContent(i, 0);
    hist_copy->SetBinError(i, 0);
  }
  gPublisherCanvas->Draw(hist_copy);
  ////////////////////////////
  TH1D *LogInverse_hist = (TH1D *)hist_copy->Clone("inverse_hist");
  for (int i = 1; i <= LogInverse_hist->GetNbinsX(); ++i) {
    double content = LogInverse_hist->GetBinContent(i);
    if (content > 0) {
      LogInverse_hist->SetBinContent(i, log(1. / content));
      double error = LogInverse_hist->GetBinError(i);
      if (error != 0) {
        LogInverse_hist->SetBinError(i, 1. / (content * content) * error);
      } else {
        LogInverse_hist->SetBinError(i, 0);
      }
    } else {
      LogInverse_hist->SetBinContent(i, 0);
      LogInverse_hist->SetBinError(i, 0);
    }
  }
  LogInverse_hist->GetYaxis()->SetRangeUser(-10, 0);
  gPublisherCanvas->Draw(LogInverse_hist);

  //////////////////////////////////
  //  pol1 fit 3.5-5
  TF1 *f3 = new TF1("f3", "pol1", 3.5, 5.);
  LogInverse_hist->Fit(f3, "0", "same", 3.5, 5.);
  gPublisherCanvas->Draw(LogInverse_hist)->DrawSame(f3);

  // auto LogInverse_hist_lowMass =
  //     (TH1D *)LogInverse_hist->Clone("LogInverse_hist_lowMass");
  // LogInverse_hist_lowMass->GetYaxis()->SetRangeUser(-10., 10);
  // double min_y_LogInverse_hist_lowMass =
  // LogInverse_hist_lowMass->GetMinimum(); for (int i = 1; i <=
  // LogInverse_hist_lowMass->GetNbinsX(); ++i) {
  //   double content = LogInverse_hist_lowMass->GetBinContent(i);
  //   if (content != 0) {
  //     content = content - min_y_LogInverse_hist_lowMass;
  //     // content = 1. / content;
  //     LogInverse_hist_lowMass->SetBinContent(i, content);
  //     cout << "Bin " << i << ": "
  //          << "Content = " << content << ", "
  //          << "Error = " << LogInverse_hist_lowMass->GetBinError(i) << endl;
  //     double error = LogInverse_hist_lowMass->GetBinError(i);
  //     if (error != 0) {
  //       LogInverse_hist_lowMass->SetBinError(i,
  //                                            1. / (content * content) *
  //                                            error);
  //     } else {
  //       LogInverse_hist_lowMass->SetBinError(i, 0);
  //     }
  //   } else {
  //     LogInverse_hist_lowMass->SetBinContent(i, 0);
  //   }
  // }

  // gPublisherCanvas->Draw(LogInverse_hist_lowMass);

  auto lFumula2 = [](double *x, double *par) {
    double term = par[0] / pow(x[0] - par[1], par[2]) + par[3];
    return term;
  };
  TF1 *f4 = new TF1("f4", lFumula2, 1.6, 2.9, 4);
  f4->SetParameters(-6.34629, 1.54991, -0.448589, 0.);
  LogInverse_hist->Fit(f4, "0", "", 1.55, 2.5);
  gPublisherCanvas->DrawSame(f4);
  gPublisherCanvas->Draw(f4);

  gPublisherCanvas->finalize();
}