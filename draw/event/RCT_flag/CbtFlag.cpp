#include "MRootIO.h"
#include "MRootGraphic.h"

void CbtFlag() {
  auto cbt_run = MRootIO::GetObjectSingle<TH1D>("cbt_mean_results_24pass1.root:h_cbt_mean");
  
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c", "c", 1800, 600);
  c->SetRightMargin(0.01);
  c->SetLeftMargin(0.05);
  cbt_run->SetTitle("N_{good events}/N_{total} vs. run id");
  cbt_run->GetYaxis()->SetTitleOffset(0.5);
  cbt_run->GetYaxis()->SetTitle("N_{good events}/N_{total}");
  cbt_run->GetXaxis()->SetTitleSize(0.005);
  cbt_run->Draw();

  for (int i = 1; i < cbt_run->GetXaxis()->GetNbins(); i++) {
    double frac = cbt_run->GetBinContent(i);
    if (frac == 0)
      cout << cbt_run->GetXaxis()->GetBinLabel(i) << endl;
  }

  c->SaveAs("cbt_mean_results_24pass1.pdf");
}
