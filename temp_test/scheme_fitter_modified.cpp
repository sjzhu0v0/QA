#include "MFit.h"
#include "MHelper.h"
#include "MMath.h"
#include "MRootGraphic.h"
#include "MRootIO.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TPad.h"

void scheme_fitter() {
  // Create a canvas to draw on
  TCanvas *canvas = new TCanvas("canvas", "Scheme Fitter Plot", 800, 600);

  TH1D *h1 = new TH1D("h1", "Example Histogram", 20, -5, 5);

  // h1->FillRandom("pol1", 10000);
  h1->Draw();
  h1->GetYaxis()->SetRangeUser(-2,2);

  vector<MDiscreteFunc> vec_func = InitOrthogonalDiscreteFuncs(20, 5);
  for (int i = 0; i < vec_func.size(); i++) {
    TF1 *f1 = new TF1(Form("f%d", i), Form("pol%d", i), -5, 5);
    f1->SetParameters(vec_func[i].GetPars(-5, 0.5).data());
    f1->Draw("same");
  }

  // Save the canvas as a PDF file
  canvas->SaveAs("scheme_fitter_plot.pdf");
}

