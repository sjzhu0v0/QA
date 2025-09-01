#define MRDF
#include "MALICE.h"
#include "MHead.h"
#include "MHist.h"
#include "MRootIO.h"
#include "TApplication.h"
#include "TLorentzVector.h"

void NJpsiCandidatePerEvent(
    TString path_input =
        "/home/szhu/data/PairFlow/22pass4_highIR/sample/O2dqflowvecd.root",
    TString path_output = "/home/szhu/work/alice/analysis/QA/output/event_jpsi/"
                          "NJpsiCandidatePerEvent.root",
    TString path_output_tree =
        "/home/szhu/work/alice/analysis/QA/output/event_jpsi/"
        "NJpsiCandidatePerEvent.root") {
  TFile *file_event = TFile::Open(path_input);
  TFile *fOutput = new TFile(path_output, "RECREATE");

  // TTree *tree_event = (TTree *)file_event->Get("O2dqflowvecd");
  TTree *tree_event = MRootIO::OpenChain(path_input, "O2dqflowvecd");

  vector<RResultHandle> gRResultHandlesFast;
  ROOT::RDataFrame rdf(*tree_event);

  auto isJPsiCandidate = [](float mass, int sign) {
    return (mass > 1.5f && mass < 5.f && sign == 0);
  };

  /* #region rdf_all definition */
  auto rdf_all =
      rdf.Define("NJpsiCandidata",
                 [isJPsiCandidate](const ROOT::RVec<float> &mass,
                                   const ROOT::RVec<float> &sign) {
                   int nJpsiCandidate = 0;
                   for (size_t i = 0; i < mass.size(); ++i) {
                     if (isJPsiCandidate(mass[i], sign[i]))
                       nJpsiCandidate++;
                   }
                   return nJpsiCandidate;
                 },
                 {"fMass", "fSign"})
          .Define(
              "mass_pair",
              [isJPsiCandidate](
                  const ROOT::RVec<float> &mass, const ROOT::RVec<float> &pt,
                  const ROOT::RVec<float> &eta, const ROOT::RVec<float> &sign,
                  const ROOT::RVec<float> &phi, const ROOT::RVec<float> &pt1,
                  const ROOT::RVec<float> &pt2, const ROOT::RVec<float> &eta1,
                  const ROOT::RVec<float> &eta2, const ROOT::RVec<float> &phi1,
                  const ROOT::RVec<float> &phi2, const ROOT::RVec<int> &sign1,
                  const ROOT::RVec<int> &sign2) {
                ROOT::VecOps::RVec<std::vector<double>> pairs;
                for (size_t i = 0; i < mass.size(); ++i) {
                  for (size_t j = i + 1; j < mass.size(); ++j) {
                    if (sign[i] == 0 && sign[j] == 0) {
                      int val1_sign1 = sign1[i];
                      int val1_sign2 = sign2[i];
                      int val2_sign1 = sign1[j];
                      int val2_sign2 = sign2[j];
                      bool doContainSameDaughter = false;
                      if (val1_sign1 == val2_sign1) {
                        if ((phi1[i] == phi1[j] && eta1[i] == eta1[j] &&
                             pt1[i] == pt1[j]) ||
                            (phi2[i] == phi2[j] && eta2[i] == eta2[j] &&
                             pt2[i] == pt2[j])) {
                          doContainSameDaughter = true;
                        }
                      } else if (val1_sign1 == val2_sign2) {
                        if ((phi1[i] == phi2[j] && eta1[i] == eta2[j] &&
                             pt1[i] == pt2[j]) ||
                            (phi2[i] == phi1[j] && eta2[i] == eta1[j] &&
                             pt2[i] == pt1[j])) {
                          doContainSameDaughter = true;
                        }
                      }
                      if (!doContainSameDaughter) {
                        TLorentzVector p1, p2;
                        p1.SetPtEtaPhiM(pt[i], eta[i], phi[i], 3.0969);
                        p2.SetPtEtaPhiM(pt[j], eta[j], phi[j], 3.0969);
                        double deltaY = p1.Rapidity() - p2.Rapidity();
                        pairs.push_back(
                            {mass[i], mass[j], pt[i], pt[j], deltaY});
                      } /* else {
                        cout << "info1_1:" << sign1[i] << " " << sign2[i] << " "
                             << pt1[i] << " " << eta1[i] << " " << phi1[i]
                             << " " << pt2[i] << " " << eta2[i] << " "
                             << phi2[i] << endl;
                        cout << "info1_2:" << sign1[j] << " " << sign2[j] << " "
                             << pt1[j] << " " << eta1[j] << " " << phi1[j]
                             << " " << pt2[j] << " " << eta2[j] << " "
                             << phi2[j] << endl;
                        cout << endl;
                      } */
                    }
                  }
                }
                return pairs;
              },
              {"fMass", "fPT", "fEta", "fSign", "fPhi", "fPt1", "fPt2", "fEta1",
               "fEta2", "fPhi1", "fPhi2", "fSign1", "fSign2"})
          .Define("mass_pair1",
                  [](const ROOT::VecOps::RVec<std::vector<double>> &pairs) {
                    ROOT::VecOps::RVec<double> mass1;
                    for (const auto &pair : pairs) {
                      mass1.push_back(pair[0]);
                    }
                    return mass1;
                  },
                  {"mass_pair"})
          .Define("mass_pair2",
                  [](const ROOT::VecOps::RVec<std::vector<double>> &pairs) {
                    ROOT::VecOps::RVec<double> mass2;
                    for (const auto &pair : pairs) {
                      mass2.push_back(pair[1]);
                    }
                    return mass2;
                  },
                  {"mass_pair"})
          .Define("pt_pair1",
                  [](const ROOT::VecOps::RVec<std::vector<double>> &pairs) {
                    ROOT::VecOps::RVec<double> mass2;
                    for (const auto &pair : pairs) {
                      mass2.push_back(pair[2]);
                    }
                    return mass2;
                  },
                  {"mass_pair"})
          .Define("pt_pair2",
                  [](const ROOT::VecOps::RVec<std::vector<double>> &pairs) {
                    ROOT::VecOps::RVec<double> mass2;
                    for (const auto &pair : pairs) {
                      mass2.push_back(pair[3]);
                    }
                    return mass2;
                  },
                  {"mass_pair"})
          .Define("deltaY",
                  [](const ROOT::VecOps::RVec<std::vector<double>> &pairs) {
                    ROOT::VecOps::RVec<double> mass2;
                    for (const auto &pair : pairs) {
                      mass2.push_back(pair[4]);
                    }
                    return mass2;
                  },
                  {"mass_pair"})
          .Define("mass_pair_together",
                  [](const ROOT::VecOps::RVec<std::vector<double>> &pairs) {
                    ROOT::VecOps::RVec<double> mass_together;
                    for (const auto &pair : pairs) {
                      mass_together.push_back(pair[0]);
                      mass_together.push_back(pair[1]);
                    }
                    return mass_together;
                  },
                  {"mass_pair"})
          .Define("k_star",
                  [](const ROOT::RVec<float> &mass,
                     const ROOT::RVec<float> &phi, const ROOT::RVec<float> &eta,
                     const ROOT::RVec<float> &pt,
                     const ROOT::RVec<float> &sign) {
                    ROOT::RVec<float> k_star;
                    for (size_t i = 0; i < mass.size(); ++i) {
                      if (mass[i] > 2.5 && mass[i] < 3.2 && sign[i] == 0)
                        k_star.push_back(-1.);
                      for (size_t j = 0; j < i; ++j) {
                        TLorentzVector p1, p2;
                        if (sign[i] == 0 && sign[j] == 0 && mass[i] > 2.5 &&
                            mass[i] < 3.2 && mass[j] > 2.5 && mass[j] < 3.2) {
                          p1.SetPtEtaPhiM(pt[i], eta[i], phi[i], 3.0969);
                          p2.SetPtEtaPhiM(pt[j], eta[j], phi[j], 3.0969);
                          TLorentzVector p_sum = p1 + p2;
                          TLorentzVector p_delta = p1 - p2;
                          p_delta.Boost(-1. * p_sum.BoostVector());
                          // cout << p_delta.Vect().Mag() << endl;
                          k_star.push_back(p_delta.Vect().Mag());
                        }
                      }
                    }
                    return k_star;
                  },
                  {"fMass", "fPhi", "fEta", "fPT", "fSign"});

  // ROOT::RDF::Experimental::AddProgressBar(rdf_all);

  rdf_all.Snapshot(
      "DiJpsi", path_output_tree,
      {"mass_pair1", "mass_pair2", "pt_pair1", "pt_pair2", "deltaY"});
  /* #endregion */

  gRResultHandlesFast.push_back(rdf_all.Histo1D(
      {"NJpsiCandidata", "N_{J/#psi candidates} (2.8<mass<3.2 GeV/c)", 6, -0.5,
       5.5},
      "NJpsiCandidata"));
  gRResultHandlesFast.push_back(rdf_all.Histo1D(
      {"k_star", "k* (GeV); k* (GeV)", 11000, -2, 10}, "k_star"));
  gRResultHandlesFast.push_back(rdf_all.Histo1D(
      {"pair_mass1d", "Mass Pair; Mass (GeV); Counts", 90, 1.8, 5.4},
      "mass_pair_together"));
  gRResultHandlesFast.push_back(rdf_all.Histo2D(
      {"pair_mass2d", "Mass Pair; Mass (GeV/c^{2}); Mass (GeV/c^{2}); Counts",
       90, 1.8, 5.4, 90, 1.8, 5.4},
      "mass_pair1", "mass_pair2"));

  RunGraphs(gRResultHandlesFast);

  fOutput->cd();
  RResultWrite(gRResultHandlesFast);
  fOutput->Close();
}

int main(int argc, char **argv) {
  TString path_input = argv[1];
  TString path_output = argv[2];
  TString path_output_tree = argv[3];
  NJpsiCandidatePerEvent(path_input, path_output, path_output_tree);
}