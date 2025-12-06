#include "MALICE.h"
#include "MHead.h"
#include "MRootGraphic.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "fstream"
#include "iostream"

double GetStatistic24(vector<int> group) { // total statistic conuting after
                                           // DiElectron selection 5.5213720e+10
  double statistic = 0;
  for (int i = 0; i < group.size(); i++) {
    // int run = MALICE_RUN::map_run24[group[i]];
    double statistic_temp = MALICE::EventNumberMinbias(
        group[i],
        "/home/szhu/work/alice/analysis/InfoRun/"
        "runinfo_LHC24_pass1_DiElectron.root:bc-selection-task/hCounterTVX");
    if (statistic_temp < 0) {
      cerr << "Run " << group[i] << " not found\n";
      exit(1);
    }
    statistic += statistic_temp;
  }
  return statistic;
}

void getstatistic() {
  vector<int> group = {550367, 550369, 550421, 550425, 550439, 550375,
                       550417, 550412, 550424, 555071, 555958, 550690,
                       550707, 550819, 550824, 550843, 550852, 552383};
  cout << "Statistic: " << GetStatistic24(group) << endl;
}