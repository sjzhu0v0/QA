#include "MALICE.h"
#include "MHead.h"
#include "MRootGraphic.h"
#include "TMatrixD.h"
#include "fstream"
#include "iostream"

// double euclideanDistance(const HistogramFeature& a, const HistogramFeature&
// b) {
//     double distance = 0.0;
//     for (size_t i = 0; i < a.features.size(); ++i) {
//         distance += TMath::Power(a.features[i] - b.features[i], 2);
//     }
//     return TMath::Sqrt(distance);
// }

double euclideanDistance(TH1 *h1, TH1 *h2) {
  double distance = 0.0;
  for (int i = 1; i <= h1->GetNbinsX(); i++) {
    distance += TMath::Power(h1->GetBinContent(i) - h2->GetBinContent(i), 2);
  }
  return sqrt(distance);
}

void calculateAndSaveDistanceMatrix(const std::vector<TH1 *> &histograms,
                                    TString filename) {
  size_t n = histograms.size();
  TMatrixD distanceMatrix(n, n);

  // 计算距离矩阵
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      double distance = euclideanDistance(histograms[i], histograms[j]);
      distanceMatrix[i][j] = distance;
      distanceMatrix[j][i] = distance;
    }
  }

  // 保存距离矩阵到文件
  TFile file(filename, "RECREATE");
  distanceMatrix.Write("distanceMatrix");
  file.Close();
}

TMatrixD loadDistanceMatrix(const std::string &filename) {
  TFile file(filename.c_str(), "READ");
  TMatrixD *distanceMatrix = nullptr;
  file.GetObject("distanceMatrix", distanceMatrix);
  if (!distanceMatrix) {
    std::cerr << "Error: Could not load distance matrix from file!"
              << std::endl;
    exit(1);
  }
  file.Close();
  return *distanceMatrix;
}

map<int, int> ClusteringWithDistanceMatrix(const TMatrixD &distanceMatrix,
                                           double threshold_distance) {
  map<int, int> cluster_map;
  int dim_matrix = distanceMatrix.GetNrows();
  for (int i = 0; i < dim_matrix; i++)
    cluster_map[i] = -1;

  int cluster_id = 0;
  while (true) {
    vector<int> vec_cluster;
    for (int i = 0; i < dim_matrix; i++) {
      if (cluster_map[i] == -1) {
        vec_cluster.push_back(i);
      }
    }
    if (vec_cluster.size() == 0)
      break;
    vector<int> vec_cluster_new;
    vec_cluster_new.push_back(vec_cluster[0]);
    for (int i = 0; i < vec_cluster.size(); i++) {
      int i_row = vec_cluster[i];
      bool doPush = true;
      for (int j = 0; j < vec_cluster_new.size(); j++) {
        int i_col = vec_cluster_new[j];
        if (distanceMatrix[i_row][i_col] > threshold_distance) {
          doPush = false;
        }
      }
      if (doPush) {
        vec_cluster_new.push_back(i_row);
        cluster_map[i_row] = cluster_id;
        break;
      }
    }
    cluster_id++;
  }
  return cluster_map;
}

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

double GetStatistic(vector<int> group) {
  double statistic = 0;
  for (int i = 0; i < group.size(); i++) {
    int run = MALICE_RUN::map_run22[group[i]];
    double statistic_temp = MALICE::EventNumberMinbias(run);
    if (statistic_temp < 0) {
      cerr << "Run " << run << " not found\n";
      exit(1);
    }
    statistic += statistic_temp;
  }
  return statistic;
}

void CreateMatrixRun22(int tag_bin, TString name = "run22") {
  vector<TH1 *> vec_h1;
  for (int i_run = 0; i_run < MALICE_RUN::map_run22.size(); i_run++) {
    int id_run = MALICE_RUN::map_run22[i_run];
    TString path = TString::Format(
        "/home/szhu/work/alice/analysis/REF/Grouping/data/%d.root", id_run);
    TFile *f = new TFile(path);
    TH1 *h1;
    if (tag_bin == -1) {
      TH2 *h2 = (TH2 *)f->Get("Phi_Pt");
      h1 = h2->ProjectionY("h1");
    } else {
      h1 = (TH1 *)f->Get(TString::Format("silceY_Phi_Pt/sliceY_%d", tag_bin));
    }
    h1->SetDirectory(0);
    // cout << path << "\n";
    h1->Scale(1.0 / h1->Integral());
    vec_h1.push_back(h1);
    f->Close();
  }
  TString path_output = TString::Format(
      "/home/szhu/work/alice/analysis/REF/Grouping/output/matrix_%s_%d.root",
      name.Data(), tag_bin);
  calculateAndSaveDistanceMatrix(vec_h1, path_output);
}

void CreateMatrixRun24() {
  vector<TH1 *> vec_h1;
  for (int i_run = 0; i_run < MALICE_RUN::map_run24.size(); i_run++) {
    int id_run = MALICE_RUN::map_run24[i_run];
    TString path = TString::Format("/home/szhu/work/alice/analysis/QA/input/"
                                   "event/phi_24pass1_DiElectro/%d/phi.root",
                                   id_run);
    TFile *f = new TFile(path);
    TH1 *h1 = MRootIO::GetObjectDiectly<TH1>(f, "Phi");
    h1->SetDirectory(0);
    h1->Scale(1.0 / h1->Integral());
    vec_h1.push_back(h1);
    f->Close();
  }
  TString path_output = "/home/szhu/work/alice/analysis/QA/output/event/"
                        "Grouping_24pass1_DiElectron/matrix.root";
  calculateAndSaveDistanceMatrix(vec_h1, path_output);
}

void CreateMatrixRun22All() {
  CreateMatrixRun22(-1);
  for (int i = 2; i <= 15; i++) {
    CreateMatrixRun22(i);
  }
}

double GroupingRun22(int tag_bin, double threshold_distance = 0.01) {
  TString path_matrix = TString::Format(
      "/home/szhu/work/alice/analysis/REF/Grouping/output/matrix_run22_%d.root",
      tag_bin);
  TMatrixD distanceMatrix = loadDistanceMatrix(path_matrix.Data());
  // map<int, int> cluster_map =
  //     ClusteringWithDistanceMatrix(distanceMatrix, threshold_distance);
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
    double stat = GetStatistic(vec_groups[i]);
    std::cout << " Statistic: " << stat << "\n";
    if (stat < 1.e10) {
      low_statistic += stat;
    }
  }
  cout << "Low statistic: " << low_statistic << "\n";

  std::vector<int> vec_groups_withoutLowStat;
  for (size_t i = 0; i < vec_groups.size(); ++i) {
    double stat = GetStatistic(vec_groups[i]);
    if (stat > 1.e10) {
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
          << MALICE_RUN::map_run22[vec_groups[vec_groups_withoutLowStat[i]][j]];
      if (j != vec_groups[vec_groups_withoutLowStat[i]].size() - 1)
        std::cout << ", ";
    }
    std::cout << "];";
    double stat = GetStatistic(vec_groups[vec_groups_withoutLowStat[i]]);
    std::cout << " Statistic: " << stat << "\n";
  }
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  TH1D *h_distance = new TH1D("h_distance", "Distance", 200, 0, 0.015);
  h_distance->SetTitle(Form("Distance matrix for bin %d , threshold %.4f, %.0f",
                            tag_bin, threshold_distance, low_statistic));
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

void Grouping(double threshold_distance = 0.0056) {
  CreateMatrixRun22All();
  double low_stat = GroupingRun22(-1, threshold_distance);
  fstream file;
  file.open("/home/szhu/work/alice/analysis/REF/Grouping/output/"
            "low_statistic.txt",
            ios::app);
  file << threshold_distance << " " << low_stat << "\n";
}