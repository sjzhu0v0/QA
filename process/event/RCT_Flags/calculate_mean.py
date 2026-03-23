import ROOT
import os
import re

# 保持脚本运行时不弹出图形界面
ROOT.gROOT.SetBatch(True)

def main():
    current_dir = os.getcwd()
    
    # 1. 查找并排序文件
    files = [f for f in os.listdir(current_dir) if f.startswith('output_') and f.endswith('.root')]
    if not files:
        print("[-] 未找到 output_*.root 文件")
        return

    # 提取 ID 并排序，确保 Bin 顺序一致
    # 创建 (filename, run_id) 元组列表并按 ID 排序
    file_list = []
    pattern = re.compile(r'output_(\d+)\.root')
    
    for filename in files:
        match = pattern.search(filename)
        if match:
            run_id = int(match.group(1))
            file_list.append((filename, run_id))
    
    # 按 Run ID 升序排序
    file_list.sort(key=lambda x: x[1])
    
    n_files = len(file_list)
    print(f"[*] 找到 {n_files} 个文件，正在创建直方图...")

    # 2. 创建 TH1D
    # Bin 从 1 到 n_files
    hist = ROOT.TH1D("h_cbt_mean", "CBT Mean per Run;Run ID;Mean Value", n_files, 0, n_files)
    hist2 = ROOT.TH1D("h_cbt_good", "CBT good;Run ID;N", n_files, 0, n_files)
    hist3 = ROOT.TH1D("h_cbt_all", "CBT good;Run ID;N", n_files, 0, n_files)
    
    # 准备画布以便后续可能的绘图操作（可选，仅用于生成图像）
    # canvas = ROOT.TCanvas("c1", "Mean Plot", 1200, 800)

    # 3. 填充直方图和设置标签
    for bin_idx, (filename, run_id) in enumerate(file_list):
        # Bin 索引从 1 开始
        bin_number = bin_idx + 1
        
        # 设置 X 轴标签为 Run ID
        hist.GetXaxis().SetBinLabel(bin_number, str(run_id))
        hist2.GetXaxis().SetBinLabel(bin_number, str(run_id))
        hist3.GetXaxis().SetBinLabel(bin_number, str(run_id))
        
        # 打开 ROOT 文件
        try:
            f = ROOT.TFile.Open(filename)
            if not f or f.IsZombie():
                print(f"[-] 无法打开文件: {filename}")
                continue
                
            # 获取直方图
            h_cbt = f.Get("cbt")
            if not h_cbt:
                print(f"[-] 文件 {filename} 中未找到 'cbt' 直方图")
                f.Close()
                continue
            
            # 获取 Mean 值
            mean_val = h_cbt.GetMean()
            total_val = h_cbt.Integral()
            good_val = h_cbt.GetBinContent(2)

            
            # 填充到直方图
            hist.SetBinContent(bin_number, mean_val)
            hist2.SetBinContent(bin_number, total_val)
            hist3.SetBinContent(bin_number, good_val)
            
            print(f"[*] {filename}: Run {run_id}, Mean = {mean_val:.4f}")
            
            f.Close()
            
        except Exception as e:
            print(f"[-] 处理文件 {filename} 时出错: {e}")

    # 4. 保存到输出文件
    out_filename = "cbt_mean_results.root"
    output_file = ROOT.TFile(out_filename, "RECREATE")
    hist.Write()
    hist2.Write()
    hist3.Write()
    output_file.Close()
    
    print(f"\n[+] 处理完成！")
    print(f"[+] 结果已保存至: {out_filename}")
    print(f"[+] 直方图包含 {n_files} 个 Run ID 的 Mean 值。")

if __name__ == "__main__":
    main()
