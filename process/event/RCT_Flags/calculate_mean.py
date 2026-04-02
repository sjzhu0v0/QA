import ROOT
import os
import re
import argparse

# 保持脚本运行时不弹出图形界面
ROOT.gROOT.SetBatch(True)


def main(input_dir, output_dir):
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)

    # 1. 查找并排序文件
    files = [
        f
        for f in os.listdir(input_dir)
        if f.startswith("output_") and f.endswith(".root")
    ]
    if not files:
        print(f"[-] 在 {input_dir} 未找到 output_*.root 文件")
        return

    # 提取 ID 并排序
    file_list = []
    pattern = re.compile(r"output_(\d+)\.root")

    for filename in files:
        match = pattern.search(filename)
        if match:
            run_id = int(match.group(1))
            file_list.append((filename, run_id))

    file_list.sort(key=lambda x: x[1])

    n_files = len(file_list)
    print(f"[*] 找到 {n_files} 个文件，正在创建直方图...")

    # 2. 创建 TH1D
    hist = ROOT.TH1D(
        "h_cbt_mean", "CBT Mean per Run;Run ID;Mean Value", n_files, 0, n_files
    )
    hist2 = ROOT.TH1D("h_cbt_good", "CBT total;Run ID;N", n_files, 0, n_files)
    hist3 = ROOT.TH1D("h_cbt_all", "CBT good;Run ID;N", n_files, 0, n_files)

    # 3. 填充直方图
    for bin_idx, (filename, run_id) in enumerate(file_list):
        bin_number = bin_idx + 1

        hist.GetXaxis().SetBinLabel(bin_number, str(run_id))
        hist2.GetXaxis().SetBinLabel(bin_number, str(run_id))
        hist3.GetXaxis().SetBinLabel(bin_number, str(run_id))

        file_path = os.path.join(input_dir, filename)

        try:
            f = ROOT.TFile.Open(file_path)
            if not f or f.IsZombie():
                print(f"[-] 无法打开文件: {file_path}")
                continue

            h_cbt = f.Get("cbt")
            if not h_cbt:
                print(f"[-] 文件 {filename} 中未找到 'cbt' 直方图")
                f.Close()
                continue

            mean_val = h_cbt.GetMean()
            total_val = h_cbt.Integral()
            good_val = h_cbt.GetBinContent(2)

            hist.SetBinContent(bin_number, mean_val)
            hist2.SetBinContent(bin_number, total_val)
            hist3.SetBinContent(bin_number, good_val)

            print(f"[*] {filename}: Run {run_id}, Mean = {mean_val:.4f}")

            f.Close()

        except Exception as e:
            print(f"[-] 处理文件 {filename} 时出错: {e}")

    # 4. 输出文件路径
    out_filename = os.path.join(output_dir, "cbt_mean_results.root")

    output_file = ROOT.TFile(out_filename, "RECREATE")
    hist.Write()
    hist2.Write()
    hist3.Write()
    output_file.Close()

    print(f"\n[+] 处理完成！")
    print(f"[+] 结果已保存至: {out_filename}")
    print(f"[+] 直方图包含 {n_files} 个 Run ID 的 Mean 值。")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process ROOT files and extract CBT statistics"
    )

    parser.add_argument(
        "-i",
        "--input",
        default="/lustre/alice/users/szhu/work/ppJpsiFlow/rct_check",
        help="输入目录",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="/lustre/alice/users/tcheng/szhu/work/repository/EventInfo",
        help="输出目录",
    )

    args = parser.parse_args()

    main(args.input, args.output)
