#include "MALICE.h"
#include "MHead.h"
#include "MRootGraphic.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "fstream"
#include "iostream"

std::vector<std::vector<int>>
hierarchicalClustering(const TMatrixD &distanceMatrix, double maxDistance) {
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
      clusters[bestI].insert(clusters[bestI].end(), clusters[bestJ].begin(),
                             clusters[bestJ].end());
      clusters.erase(clusters.begin() + bestJ);
    }
  } while (merged);
  return clusters;
}

double GetStatistic24(vector<int> group) { // total statistic conuting after
                                           // DiElectron selection 5.5213720e+10
  double statistic = 0;
  for (int i = 0; i < group.size(); i++) {
    int run = MALICE_RUN::map_run24[group[i]];
    double statistic_temp = MALICE::EventNumberMinbias(
        run,
        "/home/szhu/work/alice/analysis/InfoRun/"
        "runinfo_LHC24_pass1_DiElectron.root:bc-selection-task/hCounterTVX");
    if (statistic_temp < 0) {
      cerr << "Run " << run << " not found\n";
      exit(1);
    }
    statistic += statistic_temp;
  }
  return statistic;
}

double GroupingThreshold(double threshold_distance = 0.01,
                         double threshold_statistic = 1.e10) {
  TMatrixD *ptr_distanceMatrix = MRootIO::GetObjectDiectly<TMatrixD>(
      "/home/szhu/work/alice/analysis/QA/output/event/"
      "Grouping_24pass1_DiElectron/matrix.root:distanceMatrix");

  TMatrixD &distanceMatrix = *ptr_distanceMatrix;

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
    for (size_t j = 0; j < vec_groups[vec_groups_withoutLowStat[i]].size();
         ++j) {
      // std::cout << vec_groups[vec_groups_withoutLowStat[i]][j];
      std::cout
          << MALICE_RUN::map_run24[vec_groups[vec_groups_withoutLowStat[i]][j]];
      if (j != vec_groups[vec_groups_withoutLowStat[i]].size() - 1)
        std::cout << ", ";
    }
    std::cout << "];";
    double stat = GetStatistic24(vec_groups[vec_groups_withoutLowStat[i]]);
    std::cout << " Statistic: " << stat << "\n";
  }
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  TH1D *h_distance = new TH1D("h_distance", "Distance", 200, 0, 0.02);
  h_distance->SetTitle(Form("Distance matrix: threshold %.4f, %.0f",
                            threshold_distance, low_statistic));
  for (int i = 0; i < distanceMatrix.GetNrows(); i++) {
    for (int j = i + 1; j < distanceMatrix.GetNcols(); j++) {
      h_distance->Fill(distanceMatrix[i][j]);
    }
  }
  h_distance->Scale(1.0 / h_distance->Integral());

  TH1D *h_distance_incluster =
      new TH1D("h_distance_incluster", "Distance in cluster", 200, 0, 0.015);
  for (size_t i = 0; i < vec_groups.size(); ++i) {
    for (size_t j = 0; j < vec_groups[i].size(); ++j) {
      for (size_t k = j + 1; k < vec_groups[i].size(); ++k) {
        int i1 = vec_groups[i][j];
        int i2 = vec_groups[i][k];
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

  TH1D *h1_distance_outcluster =
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
  c1->SaveAs(TString::Format(
      "/home/szhu/work/alice/analysis/REF/Grouping/output/plotting/"
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

double GetLowStatistic(std::vector<std::vector<int>> vec_groups,
                       double threshold_statistic) {
  double low_statistic = 0;
  for (size_t i = 0; i < vec_groups.size(); ++i) {
    double stat = GetStatistic24(vec_groups[i]);
    if (stat < threshold_statistic) {
      low_statistic += stat;
    }
  }
  return low_statistic;
}

void Grouping_24pass1_DiElectron(double min_threshold_distance = 0.003,
                                 double max_threshold_distance = 0.009,
                                 int n_bins = 60) {
  TMatrixD *ptr_distanceMatrix = MRootIO::GetObjectDiectly<TMatrixD>(
      "/home/szhu/work/alice/analysis/QA/output/event/"
      "Grouping_24pass1_DiElectron/matrix.root:distanceMatrix");

  TMatrixD &distanceMatrix = *ptr_distanceMatrix;

  gPublisherCanvas = new MPublisherCanvas(
      "/home/szhu/work/alice/analysis/QA/plot/event/grouping/"
      "24pass1_DiElectron.pdf",
      1, 1);
  // plotting: distance between all the runs
  auto h_distance_all =
      new TH1D("h_distance", "Distance between all the runs; Distance;Counts",
               200, 0, 0.02);
  for (int i = 0; i < distanceMatrix.GetNrows(); i++) {
    for (int j = i + 1; j < distanceMatrix.GetNcols(); j++) {
      h_distance_all->Fill(distanceMatrix[i][j]);
    }
  }

  MRootGraphic::StyleCommon();
  MRootGraphic::StyleHistCommon(h_distance_all);
  gPublisherCanvas->Draw(h_distance_all);
  gPad->SetLogy();

  h_distance_all->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/grouping/"
      "distance_all.root");

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
       threshold_distance +=
       (max_threshold_distance - min_threshold_distance) / (double)n_bins) {
    std::vector<std::vector<int>> vec_groups =
        hierarchicalClustering(distanceMatrix, threshold_distance);

    TH1D *distance_inCluster = new TH1D(
        Form("distance_inCluster_%d", GenerateUID()),
        Form("Distance in cluster (max distance <= %.4f);Distance;Counts",
             threshold_distance),
        200, 0, 0.02);
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
    TH1D *distance_outCluster = new TH1D(
        Form("distance_outCluster_%d", GenerateUID()),
        Form("Distance out cluster (max distance <= %.4f)", threshold_distance),
        200, 0, 0.02);

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

    gPublisherCanvas->DrawClone(distance_inCluster)
        ->DrawSame(distance_outCluster);
    gPad->SetLogy();
    TLegend *legend = new TLegend(0.6, 0.75, 0.9, 0.88);
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

    TH1D *h_statistic =
        new TH1D(Form("h_statistic_%d", GenerateUID()),
                 Form("Statistic (max distance <= %.4f);;Statistic of "
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
    TLatex *latex = new TLatex();
    latex->SetTextSize(0.04);
    latex->DrawLatexNDC(
        0.3, 0.85,
        Form("Proportion of main cluster: %.1f%%", fraction_main * 100));
    vec_fraction_main.push_back(fraction_main);
    vec_threshold_distance.push_back(threshold_distance);
    cout << "integral of statistics: " << h_statistic->Integral()
         << ", maximum: " << h_statistic->GetMaximum()
         << ", Threshold distance: " << threshold_distance
         << ", Proportion of main cluster: " << fraction_main << endl;
    if (abs(threshold_distance - 0.0049) <
        (max_threshold_distance - min_threshold_distance) / (double)n_bins) {
      // print the first cluster
      cout << "First cluster: ";
      for (size_t i = 0;
           i < vec_groups[h_statistic->GetMaximumBin() - 1].size(); ++i) {
        cout << MALICE_RUN::map_run24[vec_groups[h_statistic->GetMaximumBin() -
                                                 1][i]];
        if (i != vec_groups[h_statistic->GetMaximumBin() - 1].size() - 1)
          cout << ", ";
      }
      cout << endl;
    }
  }
  TGraph *g_fraction_main =
      new TGraph(n_bins, &vec_threshold_distance[0], &vec_fraction_main[0]);
  gPublisherCanvas->SetCanvasNwNh(1, 1);
  MRootGraphic::StyleHistCommon(g_fraction_main);
  g_fraction_main->SetTitle("Proportion of main cluster vs threshold "
                            "distance;Distance threshold;"
                            "Proportion of main cluster");
  gPublisherCanvas->DrawClone(g_fraction_main);
  gPad->SetGrid();
  g_fraction_main->SaveAs(
      "/home/szhu/work/alice/analysis/QA/plot/event/grouping/"
      "proportion_main_vs_threshold_distance_fine.root");

  gPublisherCanvas->finalize();
}

// First cluster: 550367, 550369, 550421, 550425, 550439, 550375, 550417,
// 550412, 550424, 555071, 555958, 550690, 550707, 550819, 550824, 550843,
// 550852, 552383, 550691, 553819, 557233, 555881, 553250, 555150, 555152,
// 550630, 550632, 551027, 551894, 558433, 558273, 555166, 552029, 554873,
// 555254, 550634, 550654, 551083, 551066, 551468, 551394, 551504, 551463,
// 552176, 552200, 552103, 552341, 551759, 551983, 551982, 551843, 551874,
// 551418, 551943, 551958, 551889, 551923, 552384, 552402, 551875, 552139,
// 552156, 552140, 552198, 550728, 550774, 550781, 550731, 550778, 550756,
// 550858, 551398, 551997, 551760, 551989, 551761, 551926, 551877, 551925,
// 550760, 552005, 551122, 551127, 551149, 551221, 552353, 551977, 551391,
// 551392, 551498, 551387, 551979, 551389, 552138, 550784, 550742, 551008,
// 552142, 551005, 551232, 551260, 551427, 551257, 551365, 551290, 551230,
// 553253, 553297, 553299, 553305, 553785, 553844, 553255, 553274, 555411,
// 553739, 553862, 555370, 555917, 555431, 555883, 553486, 555612, 554701,
// 555705, 555801, 555723, 555722, 558122, 558126, 553588, 558291, 558354,
// 558676, 558369, 558390, 555345, 555967, 555374, 555504, 555763, 555443,
// 558627, 558633, 555695, 555546, 558155, 558757, 558179, 558244, 558182,
// 558217, 558656, 555540, 555790, 555408, 555676, 555707, 558752, 558602,
// 555478, 558604, 558329, 558685, 555850, 556182, 554203, 555543, 555742,
// 558482, 555761, 555649, 558409, 555451, 556164, 557350, 556640, 550889,
// 550916, 550997, 551856, 551922, 551921, 554316, 556237, 556734, 557744,
// 556954, 556958, 556269, 555651, 558606, 558744, 558150, 558221, 558288,
// 558215, 558327, 558330, 555900, 555933, 555960, 555965, 555976, 554201,
// 554208, 554408, 554247, 554394, 555121, 555124, 555308, 554920, 555259,
// 554970, 555073, 555156, 555226, 552283, 552285, 552304, 556981, 555591,
// 555596, 555575, 551296, 552382, 556979, 551229, 551255, 555435, 553821,
// 553903, 553185, 555122, 558284, 554524, 551013, 551023, 556662, 558615,
// 553219, 554223, 554615, 554808, 553807, 555160, 555187, 554774, 555798,
// 555860, 555853, 556664, 557681, 557012, 556834, 558275, 558551, 558422,
// 558726, 558449, 558153, 554293, 552080, 552403, 552177, 556767, 557021,
// 556152, 556218, 556160, 556716, 557339, 557613, 557659, 557691, 557717,
// 557926, 556248, 556741, 557723, 557726, 556680, 556816, 556889, 556913,
// 556923, 557019, 558750, 553187, 553188, 553189, 553193, 553555, 553663,
// 555172, 555267, 553225, 553536, 553530, 554703, 554837, 554295, 553590,
// 554323, 554880, 555208, 554998, 555022, 554194, 555482, 554588, 554354,
// 554898, 555020, 554973, 555023, 555202, 555232, 553655, 554198, 554603,
// 554752, 554768, 553700, 554772, 553702, 553880, 553756, 554736, 552102,
// 552178, 552204, 552203, 552205, 556639, 555270, 555759, 557321, 557336,
// 558387, 558124, 553512, 554835, 555344, 558535, 555401, 558247, 554322,
// 552206, 556872, 553294, 554261, 557026, 552340, 556997, 556210, 556641,
// 556907, 556939, 557547, 557374, 558383, 558406, 553610, 553633, 554404,
// 554633, 554791, 554728, 554427, 555476, 555740, 553660, 557271, 557291,
// 557299, 552369, 552381, 552401, 552400, 557226, 557251, 554462, 554558,
// 554494, 554569, 554495, 554504, 554538, 554564, 556437, 556461, 556517,
// 557149, 557119, 556454, 556485, 557112, 557138, 556542, 556491, 557415,
// 557482, 557104, 556482, 556562, 557509, 554526, 557425, 557444, 557481,
// 550653, 551890, 551931, 551780, 551993, 551992, 552179, 554507, 556497,
// 556412, 557531, 551007, 553824, 554413, 555693, 558410, 551272, 557118,
// 553816, 553825, 555047, 551107, 551924, 552201, 554471, 554714, 554207,
// 552141, 556909, 552197, 557749, 558437, 554613 Threshold distance: 0.0067,
// Proportion of main cluster: 0.985052