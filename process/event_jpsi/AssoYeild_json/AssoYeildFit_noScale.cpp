#include "MFit.h"
#include "MHelper.h"
#include "MMath.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"

vector<array<TF1 *, 2>> *gTemplate;

vector<array<TF1 *, 2>> GetTemplate(TString path_file) {
  TF1 *signal0 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/ptInt.root" + ":signal");
  TF1 *bkg0 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/ptInt.root" + ":bkg");
  signal0->SetNormalized(true);
  bkg0->SetNormalized(true);
  array<TF1 *, 2> arr0 = {signal0, bkg0};

  TF1 *signal1 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/pt0_1.root" + ":signal");
  TF1 *bkg1 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/pt0_1.root" + ":bkg");
  signal1->SetNormalized(true);
  bkg1->SetNormalized(true);
  array<TF1 *, 2> arr1 = {signal1, bkg1};

  TF1 *signal2 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/pt1_2.root" + ":signal");
  TF1 *bkg2 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/pt1_2.root" + ":bkg");
  signal2->SetNormalized(true);
  bkg2->SetNormalized(true);
  array<TF1 *, 2> arr2 = {signal2, bkg2};

  TF1 *signal3 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/pt2_3.root" + ":signal");
  TF1 *bkg3 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/pt2_3.root" + ":bkg");
  signal3->SetNormalized(true);
  bkg3->SetNormalized(true);
  array<TF1 *, 2> arr3 = {signal3, bkg3};

  TF1 *signal4 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/pt3_4.root" + ":signal");
  TF1 *bkg4 =
      MRootIO::GetObjectDiectly<TF1>(path_file + "/pt3_4.root" + ":bkg");
  signal4->SetNormalized(true);
  bkg4->SetNormalized(true);
  array<TF1 *, 2> arr4 = {signal4, bkg4};

  return {arr0, arr1, arr2, arr3, arr4};
}

TF1 *ProxyTemplate(int i, int j) {
  int j_temp = j >= 2 ? 1 : j;
  return (*gTemplate)[i][j_temp];
}

funcWithJson(void, AssoYeildFit_noScale)(
    TString path_input = "/home/szhu/work/alice/analysis/QA/test/"
                         "AssoYeildGroupEtagap_NoScale.root",
    TString path_input_mass = "/home/szhu/work/alice/analysis/QA/input/jpsi/"
                              "JpsiMass_LHC24_apass1_DiElectron.root",
    TString path_output = "/home/szhu/work/alice/analysis/QA/test/"
                          "AssoYeilFit_noScale.root",
    TString path_pdf = "/home/szhu/work/alice/analysis/QA/test/"
                       "AssoYeildFit_noScale.pdf") {
  SetUpJson();
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  Configurable<int> config_n_rebin_mass("n_rebin_mass", 3);
  int n_rebin_mass = config_n_rebin_mass.data;

  auto chi2Fit = [](MSignalFit &signal_fit_total) {
    signal_fit_total.FixBkg();
    signal_fit_total.FixSignal();
    signal_fit_total.chi2Fit();
    signal_fit_total.FixBkg(false);
    signal_fit_total.FixSignal();
    signal_fit_total.chi2Fit();
    signal_fit_total.FixBkg(false);
    signal_fit_total.FixSignal(false);
    signal_fit_total.chi2Fit();
    signal_fit_total.RemoveLimit();
    signal_fit_total.chi2Fit();
  };

  auto fit = [](MSignalFit &signal_fit_total) {
    signal_fit_total.FixBkg();
    signal_fit_total.FixSignal();
    signal_fit_total.Fit();
    signal_fit_total.FixBkg(false);
    signal_fit_total.FixSignal();
    signal_fit_total.Fit();
    signal_fit_total.FixBkg(false);
    signal_fit_total.FixSignal(false);
    signal_fit_total.Fit();
    signal_fit_total.RemoveLimit();
    signal_fit_total.Fit();
  };

  auto chi2Fit2 = [](MSignalFit &signal_fit_total) {
    signal_fit_total.FixBkg();
    signal_fit_total.FixSignal();
    signal_fit_total.chi2Fit();
    // signal_fit_total.FixBkg(false);
    // signal_fit_total.FixSignal();
    // signal_fit_total.Fit();
    // signal_fit_total.FixBkg(false);
    signal_fit_total.FixSignal(false);
    signal_fit_total.chi2Fit();
    signal_fit_total.FixBkg(false);
    signal_fit_total.chi2Fit();
  };

  auto chi2Fit3 = [](MSignalFit &signal_fit_total) {
    // signal_fit_total.FixBkg();
    // signal_fit_total.RemoveLimit();
    signal_fit_total.FixSignal();
    signal_fit_total.FixBkg(false);
    signal_fit_total.chi2Fit();
    // signal_fit_total.FixBkg(false);
    // signal_fit_total.FixSignal();
    // signal_fit_total.chi2Fit();
    // signal_fit_total.FixBkg(false);
    // signal_fit_total.FixSignal(false);
    // signal_fit_total.chi2Fit();
    // signal_fit_total.chi2Fit();
  };

  auto chi2Fit4 = [](MSignalFit &signal_fit_total) {
    signal_fit_total.FixBkg();
    signal_fit_total.FixSignal();
    signal_fit_total.chi2Fit();
    signal_fit_total.FixBkg(false);
    signal_fit_total.FixSignal();
    signal_fit_total.chi2Fit();
    // signal_fit_total.FixSignal(false);
    // signal_fit_total.FixBkg();
    // signal_fit_total.chi2Fit();
    // signal_fit_total.FixBkg(false);
    // signal_fit_total.chi2Fit();
    // signal_fit_total.FixBkg(false);
    // signal_fit_total.FixSignal();
    // signal_fit_total.chi2Fit();
    // signal_fit_total.FixBkg(false);
    // signal_fit_total.FixSignal(false);
    // signal_fit_total.chi2Fit();
    // signal_fit_total.RemoveLimit();
    // signal_fit_total.chi2Fit();
  };

  TFile *file_input2 = new TFile("/home/szhu/work/alice/analysis/QA/input/jpsi/"
                                 "JpsiQA_LHC22pass4_dqfilter.root");
  auto vec_arr_template =
      GetTemplate("/home/szhu/work/alice/analysis/QA/output/jpsi/fit_template/"
                  "LHC22pass4_dqfilter");
  gTemplate = &vec_arr_template;
  auto fPosZ_fMass_fPT_NumContribCalib_Binned = MRootIO::GetObjectDiectly<THnD>(
      file_input2, "fPosZ_MassUS_PtUS_NumContribCalib");

  TF1 *gF_signal = new TF1(Form("signal"),
                           "ROOT::Math::crystalball_function(x,[Alpha],[N],["
                           "Sigma],[Mean])",
                           2., 4.);
  gF_signal->SetParameters(0.234623, 3.06213, 5.12185, 0.0491487);
  gF_signal->SetNormalized(true);

  auto mass_total = [gF_signal](double *x, double *par) {
    return par[0] * gF_signal->Eval(x[0]) + par[1] + par[2] * x[0] +
           par[3] * x[0] * x[0] + par[4] * x[0] * x[0] * x[0];
  };

  TFile *file_input = new TFile(path_input);
  TFile *file_input_mass = new TFile(path_input_mass);
  TFile *file_output = new TFile(path_output, "RECREATE");

  struct StrAny_ptV2 {
    const vector<vector<int>> bins = {{1},
                                      {2},
                                      {3},
                                      {4},
                                      {5},
                                      {1, 2},
                                      {2, 3},
                                      {3, 4},
                                      {4, 5},
                                      {5, 6},
                                      {1, 2, 3},
                                      {2, 3, 4},
                                      {3, 4, 5},
                                      {4, 5, 6},
                                      {1, 2, 3, 4},
                                      {2, 3, 4, 5},
                                      {3, 4, 5, 6},
                                      {1, 2, 3, 4, 5},
                                      {6, 7, 8, 9, 10},
                                      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}};
    const vector<int> index_template = {1, 2, 3, 4, 4, 2, 3, 4, 4, 4,
                                        0, 0, 4, 4, 2, 3, 4, 0, 4, 0};
    const int fNbins = bins.size();
    const TString fName = "ptV2";

    vector<int> operator[](int index) { return bins[index]; }
  } strAny_ptV2;

  StrVar4Hist var_fPosZ("PosZUS", "#it{V}_{Z}", "cm", 8, {-10, 10});
  StrVar4Hist var_NumContribCalibBinned(
      "NumContribCalibUS", "N_{vtx contrib} Calibrated", "", 10,
      {0, 23, 31, 37, 43, 48, 54, 61, 69, 81, 297});
  StrVar4Hist var_MassJpsiCandidate("MassUS", "M_{ee}", "GeV^{2}/c^{4}", 90,
                                    {1.8, 5.4});
  StrVar4Hist var_PtJpsiCandidate("PtUS", "p_{T}", "GeV/c", 10, {0., 10.});
  StrVar4Hist var_DeltaEtaUS("DeltaEtaUS", "#Delta#eta_{J/#psi, track}", "", 80,
                             {-4., 4.});
  StrVar4Hist var_DeltaPhiUS("DeltaPhiUS", "#Delta#phi_{J/#psi, track}", "", 10,
                             {-M_PI_2, M_PI + M_PI_2});
  StrVar4Hist var_EtaGap("EtaGap", "#Delta#eta_{gap}", "", 6, {-0.4, 2.});
  StrVar4Hist var_PtV2Jpsi("PtV2Jpsi", "p_{T}", "GeV/c", strAny_ptV2.fNbins,
                           {0., 1.});

  MIndexHist indexHistMass(var_MassJpsiCandidate, 1, n_rebin_mass);
  MIndexHist indexHistPtJpsiCandidate(var_PtJpsiCandidate, 1, 1);
  MIndexHist indexHistDeltaPhiUS(var_DeltaPhiUS, 1, 1);
  MIndexHist indexHistDeltaEtaUS(var_DeltaEtaUS, 1, 2);
  MIndexHist indexHistEtaGap(var_EtaGap, 1, 1);
  MIndexHist indexHistPtV2Jpsi(var_PtV2Jpsi, 1, 1);
  MIndexAny indexAnyPtV2Jpsi(strAny_ptV2, 1);

  MHGroupTool3D *hg3_assoYeild_highMult = new MHGroupTool3D(
      file_input, "MassUS_DeltaEtaUS_DeltaPhiUS_AssoYeild_highMult_ptV2_%d",
      {var_PtV2Jpsi}, {1});
  MHGroupTool3D *hg3_assoYeild_lowMult = new MHGroupTool3D(
      file_input, "MassUS_DeltaEtaUS_DeltaPhiUS_AssoYeild_lowMult_ptV2_%d",
      {var_PtV2Jpsi}, {1});

  auto h_mass =
      (THnD *)file_input_mass->Get("fPosZ_MassUS_PtUS_NumContribCalib");
  MHnTool hnTool_mass(h_mass);
  hnTool_mass.Rebin(0, 200);
  hnTool_mass.Rebin(3, 5);
  hnTool_mass.PrintAllAxis();

  TH2D *mass_pt_lowMult = hnTool_mass.Project(1, 2, {0, 1});
  TH2D *mass_pt_highMult = hnTool_mass.Project(1, 2, {0, 2});

  gDirectory = nullptr;

  MHist2D h2_AssoYeild_lowMult(indexHistDeltaEtaUS, indexHistDeltaPhiUS,
                               "AssoYeild_lowMult");
  MHist2D h2_AssoYeild_highMult(indexHistDeltaEtaUS, indexHistDeltaPhiUS,
                                "AssoYeild_highMult");
  MHist2D h2_AssoYeild_Sub(indexHistDeltaEtaUS, indexHistDeltaPhiUS,
                           "AssoYeild_Sub");
  using MVec1 = MVec<MHist2D, MIndexAny<StrAny_ptV2>>;
  MVec1 h2Vec_AssoYeild_lowMult(indexAnyPtV2Jpsi, h2_AssoYeild_lowMult);
  MVec1 h2Vec_AssoYeild_highMult(indexAnyPtV2Jpsi, h2_AssoYeild_highMult);
  MVec1 h2Vec_AssoYeild_Sub(indexAnyPtV2Jpsi, h2_AssoYeild_Sub);

  MHist1D h1_yeild_lowMult(indexHistPtV2Jpsi, "yield_lowMult");
  MHist1D h1_yeild_highMult(indexHistPtV2Jpsi, "yield_highMult");
  MHist2D h2_yeild_sub(indexHistPtV2Jpsi, indexHistEtaGap, "yield_sub");

  gDirectory = nullptr;
  gPublisherCanvas = new MPublisherCanvas(path_pdf, 3, 1, 600, 600);

  for (auto iPtV2 : indexAnyPtV2Jpsi) {
    if (iPtV2 != 1) {
      gPublisherCanvas->SetCanvasNwNh(3, 1);
    }

    auto bins_pt = indexAnyPtV2Jpsi[iPtV2 - 1];

    TString str_bins_pt = Form("%d p_{T} (GeV/c): ", iPtV2);
    for (auto bin : bins_pt) {
      str_bins_pt += Form("%d_", bin);
    }

    auto mass_highMult = mass_pt_highMult->ProjectionY(
        Form("mass_highMult_ptV2_%d", iPtV2), bins_pt.front(), bins_pt.back());
    MSignalFit signal_fit_highMult(
        Form("JpsiFit_highMult_ptV2_%d", iPtV2),
        ProxyTemplate(strAny_ptV2.index_template[iPtV2 - 1], 0),
        ProxyTemplate(strAny_ptV2.index_template[iPtV2 - 1], 1), 2., 4.);
    signal_fit_highMult.InputData(mass_highMult);
    chi2Fit(signal_fit_highMult);
    signal_fit_highMult >> (gPublisherCanvas->NewPad());

    gPublisherCanvas->AddText(str_bins_pt.Data());

    auto mass_lowMult = mass_pt_lowMult->ProjectionY(
        Form("mass_lowMult_ptV2_%d", iPtV2), bins_pt.front(), bins_pt.back());
    MSignalFit signal_fit_lowMult(
        Form("JpsiFit_lowMult_ptV2_%d", iPtV2),
        ProxyTemplate(strAny_ptV2.index_template[iPtV2 - 1], 0),
        ProxyTemplate(strAny_ptV2.index_template[iPtV2 - 1], 1), 2., 4.);
    // signal_fit_lowMult.CopySignal(signal_fit_highMult);
    // signal_fit_lowMult.CopyBkg(signal_fit_highMult);
    signal_fit_lowMult.InputData(mass_lowMult);
    chi2Fit2(signal_fit_lowMult);
    signal_fit_lowMult >> (gPublisherCanvas->NewPad());

    MSignalFit signal_fit_highMult_new(
        Form("JpsiFit_highMult_ptV2_%d_new", iPtV2),
        ProxyTemplate(strAny_ptV2.index_template[iPtV2 - 1], 0),
        ProxyTemplate(strAny_ptV2.index_template[iPtV2 - 1], 1), 2., 4.);
    signal_fit_highMult_new.InputData(mass_highMult);
    chi2Fit(signal_fit_highMult_new);
    signal_fit_highMult_new.CopySignal(signal_fit_lowMult);
    signal_fit_highMult_new.FixSignal();
    signal_fit_highMult_new.chi2Fit();
    signal_fit_highMult_new >> (gPublisherCanvas->NewPad());

    gPublisherCanvas->SetCanvasNwNh(2, 1);

    auto assoYeild_highMult =
        hg3_assoYeild_highMult->GetHist(vector<int>{iPtV2});
    auto assoYeild_highMult_proj = assoYeild_highMult->ProjectionX();
    auto assoYeild_lowMult = hg3_assoYeild_lowMult->GetHist(vector<int>{iPtV2});
    auto assoYeild_lowMult_proj = assoYeild_lowMult->ProjectionX();

    MSignalFit signal_fit_asso_highMult(
        Form("JpsiFit_asso_highMult_ptV2_%d", iPtV2),
        ProxyTemplate(strAny_ptV2.index_template[iPtV2 - 1], 0),
        ProxyTemplate(strAny_ptV2.index_template[iPtV2 - 1], 1), 2., 4.);
    signal_fit_asso_highMult.InputData(assoYeild_highMult_proj);
    signal_fit_asso_highMult.CopySignal(signal_fit_highMult_new);
    signal_fit_asso_highMult.CopyBkg(signal_fit_highMult_new);
    chi2Fit4(signal_fit_asso_highMult);
    signal_fit_asso_highMult >> (gPublisherCanvas->NewPad());
    StrSignalFit str_fit_result_mass_highMult =
        signal_fit_asso_highMult.getFitResult();
    h1_yeild_highMult.fHisto->SetBinContent(
        iPtV2, str_fit_result_mass_highMult.fNsig[0]);
    h1_yeild_highMult.fHisto->SetBinError(
        iPtV2, str_fit_result_mass_highMult.fNsig[1]);

    MSignalFit signal_fit_asso_lowMult(
        Form("JpsiFit_asso_lowMult_ptV2_%d", iPtV2),
        ProxyTemplate(strAny_ptV2.index_template[iPtV2 - 1], 0),
        ProxyTemplate(strAny_ptV2.index_template[iPtV2 - 1], 1), 2., 4.);
    signal_fit_asso_lowMult.InputData(assoYeild_lowMult_proj);
    signal_fit_asso_lowMult.CopySignal(signal_fit_lowMult);
    signal_fit_asso_lowMult.CopyBkg(signal_fit_asso_highMult);
    chi2Fit4(signal_fit_asso_lowMult);
    signal_fit_asso_lowMult >> (gPublisherCanvas->NewPad());
    StrSignalFit str_fit_result_mass_lowMult =
        signal_fit_asso_lowMult.getFitResult();
    h1_yeild_lowMult.fHisto->SetBinContent(
        iPtV2, str_fit_result_mass_lowMult.fNsig[0]);
    h1_yeild_lowMult.fHisto->SetBinError(iPtV2,
                                         str_fit_result_mass_lowMult.fNsig[1]);

    gPublisherCanvas->SetCanvasNwNh(5, 5);

    for (auto i_deltaEta : indexHistDeltaEtaUS) {
      if (i_deltaEta > 29 || i_deltaEta <= 11)
        continue;
      for (auto i_deltaPhi : indexHistDeltaPhiUS) {
        auto assoYeild_diff_highMult = assoYeild_highMult->ProjectionX(
            Form("assoYeild_diff_highMult_ptV2_%d_dEta_%d_dPhi_%d", iPtV2,
                 i_deltaEta, i_deltaPhi),
            i_deltaEta, i_deltaEta, i_deltaPhi, i_deltaPhi);
        auto assoYeild_diff_lowMult = assoYeild_lowMult->ProjectionX(
            Form("assoYeild_diff_lowMult_ptV2_%d_dEta_%d_dPhi_%d", iPtV2,
                 i_deltaEta, i_deltaPhi),
            i_deltaEta, i_deltaEta, i_deltaPhi, i_deltaPhi);

        MSignalFit signal_fit_asso_diff_highMult(
            Form("JpsiFit_asso_diff_highMult_ptV2_%d_dEta_%d_dPhi_%d", iPtV2,
                 i_deltaEta, i_deltaPhi),
            ProxyTemplate(strAny_ptV2.index_template[iPtV2 - 1], 0),
            ProxyTemplate(strAny_ptV2.index_template[iPtV2 - 1], 1), 2., 4.);
        signal_fit_asso_diff_highMult.InputData(assoYeild_diff_highMult);
        signal_fit_asso_diff_highMult.CopySignal(signal_fit_asso_highMult);
        signal_fit_asso_diff_highMult.CopyBkg(signal_fit_asso_highMult);
        signal_fit_asso_diff_highMult.FixSignal();
        signal_fit_asso_diff_highMult.FixBkg(false);
        signal_fit_asso_diff_highMult.Fit();
        StrSignalFit str_fit_result_asso_diff_highMult =
            signal_fit_asso_diff_highMult.getFitResult();
        h2Vec_AssoYeild_highMult.current().SetBinInfo(
            str_fit_result_asso_diff_highMult.fNsig[0],
            str_fit_result_asso_diff_highMult.fNsig[1]);
        if (i_deltaPhi == 1) {
          signal_fit_asso_diff_highMult >> (gPublisherCanvas->NewPad());
        }
        MSignalFit signal_fit_asso_diff_lowMult(
            Form("JpsiFit_asso_diff_lowMult_ptV2_%d_dEta_%d_dPhi_%d", iPtV2,
                 i_deltaEta, i_deltaPhi),
            ProxyTemplate(strAny_ptV2.index_template[iPtV2 - 1], 0),
            ProxyTemplate(strAny_ptV2.index_template[iPtV2 - 1], 1), 2., 4.);
        signal_fit_asso_diff_lowMult.InputData(assoYeild_diff_lowMult);
        signal_fit_asso_diff_lowMult.CopySignal(signal_fit_asso_lowMult);
        signal_fit_asso_diff_lowMult.CopyBkg(signal_fit_asso_lowMult);
        signal_fit_asso_diff_lowMult.FixSignal();
        signal_fit_asso_diff_lowMult.FixBkg(false);
        signal_fit_asso_diff_lowMult.Fit();
        StrSignalFit str_fit_result_asso_diff_lowMult =
            signal_fit_asso_diff_lowMult.getFitResult();
        h2Vec_AssoYeild_lowMult.current().SetBinInfo(
            str_fit_result_asso_diff_lowMult.fNsig[0],
            str_fit_result_asso_diff_lowMult.fNsig[1]);
        if (i_deltaPhi == 1) {
          signal_fit_asso_diff_lowMult >> (gPublisherCanvas->NewPad());
        }
        assoYeild_diff_highMult->Delete();
        assoYeild_diff_lowMult->Delete();
        signal_fit_asso_diff_highMult.clean();
        signal_fit_asso_diff_lowMult.clean();
      }
    }
    mass_highMult->Delete();
    mass_lowMult->Delete();
    assoYeild_highMult->Delete();
    assoYeild_lowMult->Delete();
    signal_fit_highMult.clean();
    signal_fit_lowMult.clean();
  }

  for (auto iPtV2 : indexAnyPtV2Jpsi) {
    auto assoYeild_lowMult =
        (TH1D *)h2Vec_AssoYeild_lowMult.current().fHisto->Clone();
    auto assoYeild_highMult =
        (TH1D *)h2Vec_AssoYeild_highMult.current().fHisto->Clone();
    double yeild_low = h1_yeild_lowMult.fHisto->GetBinContent(iPtV2);
    double yeild_high = h1_yeild_highMult.fHisto->GetBinContent(iPtV2);

    assoYeild_lowMult->Scale(1. / yeild_low);
    assoYeild_highMult->Scale(1. / yeild_high);

    assoYeild_highMult->Add(assoYeild_lowMult, -1);

    for (auto iEta : indexHistDeltaEtaUS) {
      for (auto iPhi : indexHistDeltaPhiUS) {
        int bin = assoYeild_highMult->GetBin(iEta, iPhi);
        double val = assoYeild_highMult->GetBinContent(bin);
        double err = assoYeild_highMult->GetBinError(bin);
        if (val < 0) {
          val = 0;
        }
        h2Vec_AssoYeild_Sub.current().SetBinInfo(val, 0);
      }
    }
  }

  auto ApplyEtaGap = [](std::shared_ptr<TH2D> h2, double gap_eta) {
    h2->GetXaxis()->SetRangeUser(-1.8, -gap_eta / 2.);
    auto h1_highSubLow_mass1 =
        (TH1D *)h2->ProjectionY(Form("h1_highSubLow_mass1_%d", GenerateUID()));
    h2->GetXaxis()->SetRangeUser(gap_eta / 2., 1.8);
    auto h1_highSubLow_mass =
        (TH1D *)h2->ProjectionY(Form("h1_highSubLow_mass_%d", GenerateUID()));
    h1_highSubLow_mass->Add(h1_highSubLow_mass1);
    return h1_highSubLow_mass;
  };

  for (auto iPtV2 : indexAnyPtV2Jpsi) {
    auto h2_sub = h2Vec_AssoYeild_Sub.current().fHisto;
    for (auto iEtaGap : indexHistEtaGap) {
      double deltaEta = indexHistEtaGap.GetBinUpperEdge();
      auto h1_sub = ApplyEtaGap(h2_sub, deltaEta);
      double value_sub = h1_sub->GetBinContent(iEtaGap);
      h2_yeild_sub.fHisto->SetBinContent(iPtV2, iEtaGap, value_sub);
    }
  }

  gPublisherCanvas->finalize();

  file_output->Write();
  file_output->Close();
}