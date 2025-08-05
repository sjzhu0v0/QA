#include "MALICE.h"
#include "MHead.h"
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
  
}