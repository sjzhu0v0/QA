#include "MALICE.h"
#include "MHead.h"
#include "MRootGraphic.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "fstream"
#include "iostream"

std::vector<std::vector<int>> hierarchicalClustering(const TMatrixD& distanceMatrix,
                                                     double maxDistance) {
  size_t n = distanceMatrix.GetNrows();
  std::vector<std::vector<int>> clusters;

  // 初始化每个直方图为一个簇
  for (size_t i = 0; i < n; ++i) {
    clusters.push_back({static_cast<int>(i)});
  }

  bool merged;
  do {
    merged = false;
    double minMaxDistance = std::numeric_limits<double>::max();
    size_t bestI = 0, bestJ = 0;

    // 遍历所有簇对，寻找可合并的簇
    for (size_t i = 0; i < clusters.size(); ++i) {
      for (size_t j = i + 1; j < clusters.size(); ++j) {
        // 计算合并后的簇内最大距离
        double currentMax = 0.0;
        for (int a : clusters[i]) {
          for (int b : clusters[j]) {
            currentMax = std::max(currentMax, distanceMatrix[a][b]);
          }
        }
        // 检查合并后是否满足阈值
        if (currentMax <= maxDistance && currentMax < minMaxDistance) {
          minMaxDistance = currentMax;
          bestI = i;
          bestJ = j;
          merged = true;
        }
      }
    }

    // 合并最佳簇对
    if (merged) {
      clusters[bestI].insert(clusters[bestI].end(), clusters[bestJ].begin(), clusters[bestJ].end());
      clusters.erase(clusters.begin() + bestJ);
    }
  } while (merged);
  return clusters;
}

double GetStatistic24(vector<int> group) { // total statistic conuting after
                                           // DiElectron selection 5.5213720e+10
  double statistic = 0;
  for (int i = 0; i < group.size(); i++) {
    int run = MALICE_RUN::map_run24_rct[group[i]];
    double statistic_temp = MALICE::EventNumberMinbias(
        run, "/home/szhu/work/alice/analysis/InfoRun/"
             "runinfo_LHC24_pass1_DiElectron.root:bc-selection-task/hCounterTVX");
    if (statistic_temp < 0) {
      cerr << "Run " << run << " not found\n";
      exit(1);
    }
    statistic += statistic_temp;
  }
  return statistic;
}

double GetStatistic24RCT(vector<int> group) { // total statistic conuting after
                                              // DiElectron selection 5.5213720e+10
  double statistic = 0;
  for (int i = 0; i < group.size(); i++) {
    int run = MALICE_RUN::map_run24_rct[group[i]];
    double statistic_temp = MALICE::EventNumberDiElectronRCT(
        run, "/u/szhu/repository/ppJpsiFlow/event/cbt_mean_results.root:h_cbt_good");
    if (statistic_temp < 0) {
      cerr << "Run " << run << " not found\n";
      exit(1);
    }
    statistic += statistic_temp;
  }
  return statistic;
}

double GroupingThreshold(double threshold_distance = 0.01, double threshold_statistic = 1.e10) {
  TMatrixD* ptr_distanceMatrix =
      MRootIO::GetObjectDiectly<TMatrixD>("/home/szhu/work/alice/analysis/QA/output/event/"
                                          "Grouping_24pass1_DiElectron/matrix.root:distanceMatrix");

  TMatrixD& distanceMatrix = *ptr_distanceMatrix;

  std::vector<std::vector<int>> vec_groups =
      hierarchicalClustering(distanceMatrix, threshold_distance);

  std::cout << "Clusters (max distance <= " << threshold_distance << "):\n";
  double low_statistic = 0;
  for (size_t i = 0; i < vec_groups.size(); ++i) {
    std::cout << "Cluster " << i << ": [";
    for (size_t j = 0; j < vec_groups[i].size(); ++j) {
      std::cout << vec_groups[i][j];
      if (j != vec_groups[i].size() - 1)
        std::cout << ", ";
    }
    std::cout << "];";
    double stat = GetStatistic24(vec_groups[i]);
    std::cout << " Statistic: " << stat << "\n";
    if (stat < threshold_statistic) {
      low_statistic += stat;
    }
  }
  cout << "Low statistic: " << low_statistic << "\n";

  std::vector<int> vec_groups_withoutLowStat;
  for (size_t i = 0; i < vec_groups.size(); ++i) {
    double stat = GetStatistic24(vec_groups[i]);
    if (stat > threshold_statistic) {
      vec_groups_withoutLowStat.push_back(i);
    }
  }

  cout << "Clusters (max distance <= " << threshold_distance << "):\n";

  for (size_t i = 0; i < vec_groups_withoutLowStat.size(); ++i) {
    std::cout << "Cluster " << i << ": [";
    for (size_t j = 0; j < vec_groups[vec_groups_withoutLowStat[i]].size(); ++j) {
      // std::cout << vec_groups[vec_groups_withoutLowStat[i]][j];
      std::cout << MALICE_RUN::map_run24_rct[vec_groups[vec_groups_withoutLowStat[i]][j]];
      if (j != vec_groups[vec_groups_withoutLowStat[i]].size() - 1)
        std::cout << ", ";
    }
    std::cout << "];";
    double stat = GetStatistic24(vec_groups[vec_groups_withoutLowStat[i]]);
    std::cout << " Statistic: " << stat << "\n";
  }
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  TH1D* h_distance = new TH1D("h_distance", "Distance", 200, 0, 0.02);
  h_distance->SetTitle(
      Form("Distance matrix: threshold %.4f, %.0f", threshold_distance, low_statistic));
  for (int i = 0; i < distanceMatrix.GetNrows(); i++) {
    for (int j = i + 1; j < distanceMatrix.GetNcols(); j++) {
      h_distance->Fill(distanceMatrix[i][j]);
    }
  }
  h_distance->Scale(1.0 / h_distance->Integral());

  TH1D* h_distance_incluster =
      new TH1D("h_distance_incluster", "Distance in cluster", 200, 0, 0.015);
  for (size_t i = 0; i < vec_groups.size(); ++i) {
    for (size_t j = 0; j < vec_groups[i].size(); ++j) {
      for (size_t k = j + 1; k < vec_groups[i].size(); ++k) {
        int i1 = vec_groups[i][j];
        int i2 = vec_groups[i][k];
        if (distanceMatrix[i1][i2] > threshold_distance)
          cout << "Warning: distance " << distanceMatrix[i1][i2] << " in cluster " << i
               << "Threshold: " << threshold_distance << "\n";
        h_distance_incluster->Fill(distanceMatrix[i1][i2]);
      }
    }
  }
  h_distance_incluster->Scale(1.0 / h_distance_incluster->Integral());
  h_distance_incluster->SetLineColor(kRed);
  h_distance->GetYaxis()->SetRangeUser(0, 0.1);
  h_distance_incluster->GetYaxis()->SetRangeUser(0, 0.1);

  h_distance_incluster->Fit("gaus");
  h_distance->Draw();
  h_distance_incluster->Draw("same");

  TH1D* h1_distance_outcluster =
      new TH1D("h1_distance_outcluster", "Distance out cluster", 200, 0, 0.015);
  for (int i_group = 0; i_group < vec_groups.size(); i_group++) {
    for (int i = 0; i < distanceMatrix.GetNrows(); i++) {
      bool is_in_cluster = false;
      for (int j = 0; j < vec_groups[i_group].size(); j++) {
        if (i == vec_groups[i_group][j]) {
          is_in_cluster = true;
          break;
        }
      }
      if (!is_in_cluster) {
        for (int j = 0; j < vec_groups[i_group].size(); j++) {
          int i1 = i;
          int i2 = vec_groups[i_group][j];
          h1_distance_outcluster->Fill(distanceMatrix[i1][i2]);
        }
      }
    }
  }
  h1_distance_outcluster->Scale(1.0 / h1_distance_outcluster->Integral());
  h1_distance_outcluster->SetLineColor(kGreen);
  h1_distance_outcluster->Draw("same");
  c1->SaveAs(TString::Format("/home/szhu/work/alice/analysis/REF/Grouping/output/plotting/"
                             "distance_%.4f.pdf",
                             threshold_distance));

  // for (int i = 0; i < distanceMatrix.GetNrows(); i++) {
  //   for (int j = i + 1; j < distanceMatrix.GetNcols(); j++) {
  //     if (distanceMatrix[i][j] > 0.01) {
  //       cout << "Histogram " << i << " and " << j << " have distance "
  //            << distanceMatrix[i][j] << "\n";
  //     }
  //   }
  // }

  // for (auto &cluster : cluster_map) {
  //   cout << "Histogram " << cluster.first << " belongs to cluster "
  //        << cluster.second << "\n";
  // }
  // print the cluster map according to the cluster id

  // map<int, vector<int>> cluster_map_reverse;
  // for (auto &cluster : cluster_map) {
  //   cluster_map_reverse[cluster.second].push_back(cluster.first);
  // }
  // for (auto &cluster : cluster_map_reverse) {
  //   cout << "Cluster " << cluster.first << " contains histograms: ";
  //   for (auto &hist : cluster.second) {
  //     cout << hist << ",";
  //   }
  //   cout << "\n";
  //   TH1D *h1 = new TH1D(Form("h1_%d", cluster.first), "Cluster", 100, 0,
  //   0.05); for (auto &hist : cluster.second) {
  //     for (int i = hist + 1; i < distanceMatrix.GetNcols(); i++) {
  //       h1->Fill(distanceMatrix[hist][i]);
  //     }
  //   }
  //   if (cluster.first >= 4)
  //     continue;
  //   int color = GetColor(MColorSpace::RGB_4_0[cluster.first]);
  //   h1->SetLineColor(color);
  //   h1->Scale(1.0 / h1->Integral());
  //   h1->Draw("same");
  // }
  return low_statistic;
}

double GetLowStatistic(std::vector<std::vector<int>> vec_groups, double threshold_statistic) {
  double low_statistic = 0;
  for (size_t i = 0; i < vec_groups.size(); ++i) {
    double stat = GetStatistic24(vec_groups[i]);
    if (stat < threshold_statistic) {
      low_statistic += stat;
    }
  }
  return low_statistic;
}

void Grouping_24pass1_RCT(double min_threshold_distance = 0.003,
                          double max_threshold_distance = 0.009, int n_bins = 60) {
  TMatrixD* ptr_distanceMatrix = MRootIO::GetObjectDiectly<TMatrixD>(
      "/u/szhu/repository/ppJpsiFlow/event/run_grouping/matrix_grouping_24cbt.root:distanceMatrix");

  TMatrixD& distanceMatrix = *ptr_distanceMatrix;

  gPublisherCanvas = new MPublisherCanvas(
      "/u/szhu/repository/ppJpsiFlow/event/run_grouping/plots_grouping_24pass1_rct.pdf", 1, 1);
  // plotting: distance between all the runs
  auto h_distance_all =
      new TH1D("h_distance", "Distance between all the runs; Distance;Counts", 200, 0, 0.02);
  for (int i = 0; i < distanceMatrix.GetNrows(); i++) {
    for (int j = i + 1; j < distanceMatrix.GetNcols(); j++) {
      h_distance_all->Fill(distanceMatrix[i][j]);
    }
  }

  MRootGraphic::StyleCommon();
  MRootGraphic::StyleHistCommon(h_distance_all);
  gPublisherCanvas->Draw(h_distance_all);
  gPad->SetLogy();

  h_distance_all->SaveAs("/u/szhu/repository/ppJpsiFlow/event/run_grouping/"
                         "distance_rct.root");

  // double threshold_distance = 0.01;
  // std::vector<std::vector<int>> vec_groups =
  //     hierarchicalClustering(distanceMatrix, threshold_distance);

  // double threshold_statistic = 1.e10;
  // cout << GetLowStatistic(vec_groups, threshold_statistic) << endl;

  gPublisherCanvas->SetCanvasNwNh(2, 1);

  vector<double> vec_fraction_main;
  vector<double> vec_threshold_distance;

  for (double threshold_distance = min_threshold_distance;
       threshold_distance <= max_threshold_distance;
       threshold_distance += (max_threshold_distance - min_threshold_distance) / (double)n_bins) {
    std::vector<std::vector<int>> vec_groups =
        hierarchicalClustering(distanceMatrix, threshold_distance);

    TH1D* distance_inCluster =
        new TH1D(Form("distance_inCluster_%d", GenerateUID()),
                 Form("Dissimilarity threshold = %.4f", threshold_distance), 200, 0, 0.02);
    distance_inCluster->GetXaxis()->SetTitle("Dissimilarity");
    distance_inCluster->GetYaxis()->SetTitle("N_{run pairs}");
    for (size_t i = 0; i < vec_groups.size(); ++i) {
      for (size_t j = 0; j < vec_groups[i].size(); ++j) {
        for (size_t k = j + 1; k < vec_groups[i].size(); ++k) {
          int i1 = vec_groups[i][j];
          int i2 = vec_groups[i][k];
          distance_inCluster->Fill(distanceMatrix[i1][i2]);
        }
      }
    }
    MRootGraphic::StyleHistCommon(distance_inCluster);
    distance_inCluster->SetLineColor(kRed);
    // distance_inCluster->Scale(1.0 / distance_inCluster->Integral());
    // gPublisherCanvas->DrawClone(distance_inCluster);
    // gPad->SetLogy();
    TH1D* distance_outCluster =
        new TH1D(Form("distance_outCluster_%d", GenerateUID()),
                 Form("Dissimilarity threshold = %.4f", threshold_distance), 200, 0, 0.02);
    distance_outCluster->GetXaxis()->SetTitle("Dissimilarity");
    distance_outCluster->GetYaxis()->SetTitle("N_{run pairs}");

    for (int i_group = 0; i_group < vec_groups.size(); i_group++) {
      for (int i = 0; i < distanceMatrix.GetNrows(); i++) {
        bool is_in_cluster = false;
        for (int j = 0; j < vec_groups[i_group].size(); j++) {
          if (i == vec_groups[i_group][j]) {
            is_in_cluster = true;
            break;
          }
        }
        if (!is_in_cluster) {
          for (int j = 0; j < vec_groups[i_group].size(); j++) {
            int i1 = i;
            int i2 = vec_groups[i_group][j];
            distance_outCluster->Fill(distanceMatrix[i1][i2]);
          }
        }
      }
    }
    MRootGraphic::StyleHistCommon(distance_outCluster);
    distance_outCluster->SetLineColor(kGreen + 2);
    // gPublisherCanvas->DrawClone(distance_outCluster);
    // gPad->SetLogy();
    // distance_inCluster->Scale(1.0 / distance_inCluster->Integral());
    // distance_outCluster->Scale(1.0 / distance_outCluster->Integral());

    gPublisherCanvas->DrawClone(distance_inCluster)->DrawSame(distance_outCluster);
    gPad->SetLogy();
    TLegend* legend = new TLegend(0.6, 0.75, 0.9, 0.88);
    legend->AddEntry(distance_inCluster, "In cluster", "l");
    legend->AddEntry(distance_outCluster, "Out cluster", "l");
    legend->SetTextSize(0.04);
    legend->SetLineColor(0);
    legend->SetFillColor(0);
    legend->Draw();

    vector<double> statistics_cluster;
    for (size_t i = 0; i < vec_groups.size(); ++i) {
      double stat = GetStatistic24RCT(vec_groups[i]);
      statistics_cluster.push_back(stat);
    }

    TH1D* h_statistic = new TH1D(Form("h_statistic_%d", GenerateUID()),
                                 Form("Statistic (max dissimilarity <= %.4f);;Statistic of "
                                      "the cluster",
                                      threshold_distance),
                                 vec_groups.size(), 0, vec_groups.size());
    for (size_t i = 0; i < vec_groups.size(); ++i) {
      h_statistic->SetBinContent(i + 1, statistics_cluster[i]);
      h_statistic->GetXaxis()->SetBinLabel(i + 1, Form("Cluster %zu", i + 1));
    }
    MRootGraphic::StyleHistCommon(h_statistic);
    gPublisherCanvas->DrawClone(h_statistic);
    gPad->SetLogy();
    double fraction_main = h_statistic->GetMaximum() / h_statistic->Integral();
    TLatex* latex = new TLatex();
    latex->SetTextSize(0.04);
    latex->DrawLatexNDC(0.3, 0.85,
                        Form("Proportion of the largest group : %.1f%%", fraction_main * 100));
    vec_fraction_main.push_back(fraction_main);
    vec_threshold_distance.push_back(threshold_distance);
    cout << "integral of statistics: " << h_statistic->Integral()
         << ", maximum: " << h_statistic->GetMaximum()
         << ", Threshold distance: " << threshold_distance
         << ", Proportion of main cluster: " << fraction_main << endl;
    if (abs(threshold_distance - 0.007) <
        (max_threshold_distance - min_threshold_distance) / (double)n_bins) {
      // print the first cluster
      cout << "First cluster: ";
      for (size_t i = 0; i < vec_groups[h_statistic->GetMaximumBin() - 1].size(); ++i) {
        cout << MALICE_RUN::map_run24_rct[vec_groups[h_statistic->GetMaximumBin() - 1][i]];
        if (i != vec_groups[h_statistic->GetMaximumBin() - 1].size() - 1)
          cout << ", ";
      }
      cout << endl;
    }
  }
  TGraph* g_fraction_main = new TGraph(n_bins, &vec_threshold_distance[0], &vec_fraction_main[0]);
  gPublisherCanvas->SetCanvasNwNh(1, 1);
  MRootGraphic::StyleHistCommon(g_fraction_main);
  g_fraction_main->SetTitle("Proportion of the largest group vs threshold "
                            "dissimilarity;Dissimilarity threshold;"
                            "Proportion of the largest group");
  gPublisherCanvas->DrawClone(g_fraction_main);
  gPad->SetGrid();
  g_fraction_main->SaveAs("/u/szhu/repository/ppJpsiFlow/event/run_grouping/"
                          "proportion_main_vs_threshold_distance_fine.root");

  gPublisherCanvas->finalize();
}

void Grouping_24pass1_DiElectron_Detailed(double threshold_distance = 0.0045) {
  TMatrixD* ptr_distanceMatrix =
      MRootIO::GetObjectDiectly<TMatrixD>("/home/szhu/work/alice/analysis/QA/output/event/"
                                          "Grouping_24pass1_DiElectron/matrix.root:distanceMatrix");

  TMatrixD& distanceMatrix = *ptr_distanceMatrix;

  gPublisherCanvas = new MPublisherCanvas("/home/szhu/work/alice/analysis/QA/plot/event/grouping/"
                                          "24pass1_DiElectron_detailed.pdf",
                                          1, 1);
  // plotting: distance between all the runs
  auto h_distance_all =
      new TH1D("h_distance", "Distance between all the runs; Distance;Counts", 200, 0, 0.02);
  for (int i = 0; i < distanceMatrix.GetNrows(); i++) {
    for (int j = i + 1; j < distanceMatrix.GetNcols(); j++) {
      h_distance_all->Fill(distanceMatrix[i][j]);
    }
  }

  MRootGraphic::StyleCommon();
  MRootGraphic::StyleHistCommon(h_distance_all);
  gPublisherCanvas->Draw(h_distance_all);
  gPad->SetLogy();

  h_distance_all->SaveAs("/home/szhu/work/alice/analysis/QA/plot/event/grouping/"
                         "distance_all.root");

  // double threshold_distance = 0.01;
  // std::vector<std::vector<int>> vec_groups =
  //     hierarchicalClustering(distanceMatrix, threshold_distance);

  // double threshold_statistic = 1.e10;
  // cout << GetLowStatistic(vec_groups, threshold_statistic) << endl;

  gPublisherCanvas->SetCanvasNwNh(2, 1);

  vector<double> vec_fraction_main;
  vector<double> vec_threshold_distance;

  std::vector<std::vector<int>> vec_groups =
      hierarchicalClustering(distanceMatrix, threshold_distance);

  TH1D* distance_inCluster =
      new TH1D(Form("distance_inCluster_%d", GenerateUID()),
               Form("Dissimilarity threshold <= %.4f", threshold_distance), 200, 0, 0.02);
  distance_inCluster->GetXaxis()->SetTitle("Dissimilarity");
  distance_inCluster->GetYaxis()->SetTitle("N_{run pairs}");
  for (size_t i = 0; i < vec_groups.size(); ++i) {
    for (size_t j = 0; j < vec_groups[i].size(); ++j) {
      for (size_t k = j + 1; k < vec_groups[i].size(); ++k) {
        int i1 = vec_groups[i][j];
        int i2 = vec_groups[i][k];
        distance_inCluster->Fill(distanceMatrix[i1][i2]);
      }
    }
  }
  MRootGraphic::StyleHistCommon(distance_inCluster);
  distance_inCluster->SetLineColor(kRed);
  // distance_inCluster->Scale(1.0 / distance_inCluster->Integral());
  // gPublisherCanvas->DrawClone(distance_inCluster);
  // gPad->SetLogy();
  TH1D* distance_outCluster =
      new TH1D(Form("distance_outCluster_%d", GenerateUID()),
               Form("Dissimilarity threshold <= %.4f", threshold_distance), 200, 0, 0.02);
  distance_outCluster->GetXaxis()->SetTitle("Dissimilarity");
  distance_outCluster->GetYaxis()->SetTitle("N_{run pairs}");

  for (int i_group = 0; i_group < vec_groups.size(); i_group++) {
    for (int i = 0; i < distanceMatrix.GetNrows(); i++) {
      bool is_in_cluster = false;
      for (int j = 0; j < vec_groups[i_group].size(); j++) {
        if (i == vec_groups[i_group][j]) {
          is_in_cluster = true;
          break;
        }
      }
      if (!is_in_cluster) {
        for (int j = 0; j < vec_groups[i_group].size(); j++) {
          int i1 = i;
          int i2 = vec_groups[i_group][j];
          distance_outCluster->Fill(distanceMatrix[i1][i2]);
        }
      }
    }
  }
  MRootGraphic::StyleHistCommon(distance_outCluster);
  distance_outCluster->SetLineColor(kGreen + 2);
  // gPublisherCanvas->DrawClone(distance_outCluster);
  // gPad->SetLogy();
  // distance_inCluster->Scale(1.0 / distance_inCluster->Integral());
  // distance_outCluster->Scale(1.0 / distance_outCluster->Integral());

  gPublisherCanvas->DrawClone(distance_inCluster)->DrawSame(distance_outCluster);
  gPad->SetLogy();
  TLegend* legend = new TLegend(0.6, 0.75, 0.9, 0.88);
  legend->AddEntry(distance_inCluster, "In cluster", "l");
  legend->AddEntry(distance_outCluster, "Out cluster", "l");
  legend->SetTextSize(0.04);
  legend->SetLineColor(0);
  legend->SetFillColor(0);
  legend->Draw();

  vector<double> statistics_cluster;
  for (size_t i = 0; i < vec_groups.size(); ++i) {
    double stat = GetStatistic24(vec_groups[i]);
    statistics_cluster.push_back(stat);
  }

  TH1D* h_statistic = new TH1D(Form("h_statistic_%d", GenerateUID()),
                               Form("Statistic (max dissimilarity <= %.4f);;Statistic of "
                                    "the cluster",
                                    threshold_distance),
                               vec_groups.size(), 0, vec_groups.size());
  for (size_t i = 0; i < vec_groups.size(); ++i) {
    h_statistic->SetBinContent(i + 1, statistics_cluster[i]);
    h_statistic->GetXaxis()->SetBinLabel(i + 1, Form("Cluster %zu", i + 1));
  }
  MRootGraphic::StyleHistCommon(h_statistic);
  gPad->SetLogy();
  double fraction_main = h_statistic->GetMaximum() / h_statistic->Integral();
  h_statistic->SaveAs("test.root");
  TLatex* latex = new TLatex();
  latex->SetTextSize(0.04);
  latex->DrawLatexNDC(0.3, 0.85,
                      Form("Proportion of the largest group : %.1f%%", fraction_main * 100));
  vec_fraction_main.push_back(fraction_main);
  vec_threshold_distance.push_back(threshold_distance);
  cout << "integral of statistics: " << h_statistic->Integral()
       << ", maximum: " << h_statistic->GetMaximum()
       << ", Threshold distance: " << threshold_distance
       << ", Proportion of main cluster: " << fraction_main << endl;
  // print the first cluster
  cout << "First cluster: ";
  for (size_t i = 0; i < vec_groups[h_statistic->GetMaximumBin() - 1].size(); ++i) {
    cout << MALICE_RUN::map_run24_rct[vec_groups[h_statistic->GetMaximumBin() - 1][i]];
    if (i != vec_groups[h_statistic->GetMaximumBin() - 1].size() - 1)
      cout << ", ";
  }
  cout << endl;

  cout << "Second cluster: ";
  cout << h_statistic->GetMaximumBin() << endl;
  double integral_temp = h_statistic->Integral();
  h_statistic->SetBinContent(h_statistic->GetMaximumBin(), 0);
  cout << "proportion" << h_statistic->GetMaximum() / integral_temp << endl;

  for (size_t i = 0; i < vec_groups[h_statistic->GetMaximumBin() - 1].size(); ++i) {
    550367, 556639, 550369, 550425, 550439, 550375, 550421, 550412, 550424, 554998, 555853, 550417,
        550690, 550707, 550819, 550824, 550843, 550852, 550634, 551926, 551761, 551877, 551875,
        552005, 553588, 553486, 556461, 555543, 554588, 555575, 558330, 557744, 557613, 558410,
        555649, 555759, 555693, 558122, 558284, 555960, 555965, 556542, 556517, 557350, 556497,
        556641, 555881, 555958, 556562, 557021, 557251, 557482, 556640, 557233, 556716, 557321,
        557299, 556491, 556997, 557336, 557339, 555308, 557862, 555451, 555860, 555798, 558150,
        558390, 557926, 555705, 555967, 558182, 558126, 555761, 555676, 557547, 555740, 557726,
        557717, 558387, 558153, 557509, 558155, 558215, 553821, 553903, 552141, 552197, 557659,
        554558, 558273, 551107, 551924, 552201, 551925, 554404, 554613, 554098, 552080, 552381,
        552403, 552177, 552204, 552102, 555259, 555612, 558369, 558383, 555270, 555401, 555723,
        555443, 557749, 555478, 557681, 558437, 557691, 555344, 555591, 555651, 558433, 555431,
        558329, 555883, 556741, 555917, 555933, 556662, 557271, 553555, 553663, 553610, 553633,
        554354, 554569, 554736, 554615, 555370, 558179, 558124, 557074, 555695, 555411, 553187,
        553188, 554208, 553189, 553193, 553530, 554603, 554772, 554791, 555160, 554920, 555150,
        555208, 556767, 556981, 557374, 558327, 557531, 555707, 557012, 556437, 552340, 552401,
        552369, 552400, 556979, 553739, 553862, 554223, 555267, 555850, 557019, 555801, 556939,
        556958, 553512, 554768, 555254, 558288, 553816, 552178, 552203, 552205, 558406, 552283,
        552285, 552304, 556680, 553225, 553536, 555022, 554293, 554808, 554898, 555172, 554873,
        554970, 555156, 554835, 555202, 554837, 555226, 554194, 554322, 554095, 553297, 553299,
        553305, 553785, 553844, 553702, 553880, 553756, 554092, 554526, 555435, 555540, 554701,
        553255, 555121, 555166, 553590, 554538, 554094, 554295, 554261, 554316, 554198, 554564,
        554752, 554880, 555047, 555232, 553655, 554703, 554714, 553700, 554728, 555020, 554774,
        555187, 553660, 554201, 554247, 553294, 554207, 554203, 556816, 554394, 554504, 554408,
        554427, 554524, 554495, 554507, 554413, 556182, 556218, 556370, 556923, 556909, 554494,
        557112, 556210, 556269, 556248, 556889, 556913, 557138, 556372, 557104, 556872, 557149,
        556237, 556412, 557118, 557119, 550630, 550632, 550653, 551391, 551005, 551931, 551398,
        551760, 551989, 551997, 551498, 551890, 551392, 551856, 551921, 551992, 550654, 551468,
        551083, 551943, 551418, 551958, 551365, 551983, 551889, 551982, 551394, 551504, 551066,
        551759, 551843, 551463, 551874, 552103, 552341, 551923, 552176, 552200, 552402, 552384,
        551232, 551260, 551290, 551257, 551427, 550889, 550916, 550997, 551922, 552198, 552139,
        552156, 552140, 550728, 550774, 550781, 550731, 550778, 550756, 550858, 551230, 551013,
        551023, 551780, 551993, 551296, 555482, 553219, 553274, 553807, 555122, 555152, 553253,
        554732, 555722, 555345, 555742, 555790, 557723, 556482, 558422, 555504, 556485, 557291,
        558275, 556734, 558291, 557897, 558244, 558409, 555763, 558354, 552383, 551007, 552179,
        552206, 554462, 556164, 557226, 556284, 553824, 554323, 558221, 555596, 551272, 556907,
        550742, 551008, 552142, 551255, 552382, 554097, 555476, 556664, 555374, 555546, 558217,
        553185, 555023, 557913, 554471, 550760, 552029, 551122, 551127, 551149, 551221, 552353,
        551979, 551977, 551387, 551389, 552138, 550784, 551229, 551027, 551894, 558247, 557876,
        555124, 553250, 555071, 555408, 555900, 557026, 556454, 553825, 554973, 550691, 553819,
        556954, 555073, 555976, 556160,
        556152 cout << MALICE_RUN::map_run24_rct[vec_groups[h_statistic->GetMaximumBin() - 1][i]];
    if (i != vec_groups[h_statistic->GetMaximumBin() - 1].size() - 1)
      cout << ", ";
  }
  cout << endl;
  gPublisherCanvas->finalize();
}
