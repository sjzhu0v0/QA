#include <TFormula.h>
#include <TMath.h>

#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {
struct CutValues {
  double e1_nsig_el = 0.;
  double e1_nsig_pi = 4.;
  double e1_nsig_pr = 4.;
  double e2_nsig_el = 0.;
  double e2_nsig_pi = 4.;
  double e2_nsig_pr = 4.;
  double ref_ITSChi2NCl = 3.;
  double ref_TPCNClsFound = 90.;
  double nITSCluster = 3.;
  double nDcaZ2Dev = 2.;
  double nDcaXY2Dev = 2.;
  double ref_pt = 2.;
  double fPosZ = 0.;
  double fPT = 0.8;
  double fPt1 = 1.2;
  double fPt2 = 1.3;
};

std::string FormulaExpr(const std::string& expr) {
  std::string out = expr;
  const std::vector<std::pair<std::string, std::string>> replacements = {
      {"e1_nsig_el", "x[0]"},     {"e1_nsig_pi", "x[1]"},       {"e1_nsig_pr", "x[2]"},
      {"e2_nsig_el", "x[3]"},     {"e2_nsig_pi", "x[4]"},       {"e2_nsig_pr", "x[5]"},
      {"ref_ITSChi2NCl", "x[6]"}, {"ref_TPCNClsFound", "x[7]"}, {"nITSCluster", "x[8]"},
      {"nDcaZ2Dev", "x[9]"},      {"nDcaXY2Dev", "x[10]"},      {"ref_pt", "x[11]"},
      {"fPosZ", "x[12]"},         {"fPT", "x[13]"},             {"fPt1", "x[14]"},
      {"fPt2", "x[15]"},          {"abs(", "TMath::Abs("}};
  for (const auto& [from, to] : replacements) {
    size_t pos = 0;
    while ((pos = out.find(from, pos)) != std::string::npos) {
      out.replace(pos, from.size(), to);
      pos += to.size();
    }
  }
  return out;
}

bool PassCut(TFormula& formula, const CutValues& v) {
  double values[] = {v.e1_nsig_el,     v.e1_nsig_pi,       v.e1_nsig_pr,
                     v.e2_nsig_el,     v.e2_nsig_pi,       v.e2_nsig_pr,
                     v.ref_ITSChi2NCl, v.ref_TPCNClsFound, v.nITSCluster,
                     v.nDcaZ2Dev,      v.nDcaXY2Dev,       v.ref_pt,
                     v.fPosZ,          v.fPT,              v.fPt1,
                     v.fPt2};
  return formula.EvalPar(values) != 0.;
}

bool CheckCase(TFormula& formula, const char* label, const CutValues& values, bool expected) {
  const bool actual = PassCut(formula, values);
  std::cout << label << ": " << (actual ? "pass" : "fail")
            << " expected " << (expected ? "pass" : "fail") << '\n';
  return actual == expected;
}
} // namespace

int cut_checker_sample() {
  const std::string cut =
      "e1_nsig_pi > 3 && e2_nsig_pi > 3 && "
      "e1_nsig_pr > 3 && e2_nsig_pr > 3 && "
      "e1_nsig_el > -2 && e2_nsig_el > -2 && "
      "e1_nsig_el < 3 && e2_nsig_el < 3 && "
      "ref_ITSChi2NCl < 4 && ref_TPCNClsFound > 80 && "
      "nITSCluster > 2.5 && nDcaZ2Dev < 3 && nDcaXY2Dev < 3 && "
      "ref_pt > 1. && ref_pt < 3 && "
      "(fPt1+fPt2<2.7&&fPT<1 || fPT>1)";

  const std::string formula_expr = FormulaExpr(cut);
  std::cout << "Formula: " << formula_expr << "\n\n";

  TFormula formula("cut_checker_sample", formula_expr.c_str());

  bool ok = true;

  CutValues base;
  ok &= CheckCase(formula, "low-pT J/psi and daughter sum below 2.7", base, true);

  CutValues high_sum_low_jpsi = base;
  high_sum_low_jpsi.fPt1 = 1.5;
  high_sum_low_jpsi.fPt2 = 1.4;
  high_sum_low_jpsi.fPT = 0.8;
  ok &= CheckCase(formula, "low-pT J/psi and daughter sum above 2.7", high_sum_low_jpsi, false);

  CutValues high_jpsi = high_sum_low_jpsi;
  high_jpsi.fPT = 1.2;
  ok &= CheckCase(formula, "high-pT J/psi bypasses daughter-sum cut", high_jpsi, true);

  CutValues failed_ref = base;
  failed_ref.ref_pt = 3.5;
  ok &= CheckCase(formula, "base cut still rejects ref_pt outside range", failed_ref, false);

  return ok ? 0 : 1;
}
