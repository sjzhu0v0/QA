#include "MRootIO.h"
#include "RCTSelectionFlags.h"


using namespace o2::aod::rctsel;
using namespace std;

#define MInfo(x) cout << #x ":" << x << endl;

void test(const string& inputPath, const string& outputPath) {
  TChain* tree_input = MRootIO::OpenChain(inputPath.c_str(), "O2rctrawdq");
  MInfo(tree_input->GetEntries());

  RCTFlagsChecker rctChecker{"CBT"};

  MCollision collision;
  collision.init(tree_input);

  TH1D* h_rct = new TH1D("cbt", "cbt pass;cbt value;Counts", 2, -0.5, 1.5);
  for (Long64_t i = 0; i < tree_input->GetEntries(); i++) {
    tree_input->GetEntry(i); // 读取当前事件
    rctChecker(collision) ? h_rct->Fill(1) : h_rct->Fill(0);
  }

  TFile* file_output = new TFile(outputPath.c_str(), "recreate");
  file_output->cd();
  h_rct->Write();
  file_output->Close();
}

int main(int argc, char* argv[]) {
  // 检查命令行参数数量（需传入输入路径、输出路径）
  if (argc != 3) {
    cerr << "Usage: " << argv[0] << " <input_path> <output_path>" << endl;
    return 1;
  }

  string inputPath = argv[1];  // 输入文件路径
  string outputPath = argv[2]; // 输出文件路径

  test(inputPath, outputPath);
  return 0;
}
