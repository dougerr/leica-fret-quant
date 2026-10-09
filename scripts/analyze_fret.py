# -*- coding: utf-8 -*-
"""
徕卡共聚焦 FRET 定量（BFP 供体 / 敏化发射），批量分析主脚本。

通道约定（LAS X 导出 0-based）：
    ch03 = 供体 BFP       (405nm 激发 / 450nm 发射)
    ch00 = FRET 敏化发射   (405nm 激发 / 受体发射窗口)
    ch01 = 受体直接激发    (受体表达量)
    ch02 = 明场(弃用)

三个指标（全部输出，供对照）：
    ratio_raw = mean(S_all)/mean(D_all)        整视野均值比(= Trae alpha=0 NFRET 口径，受背景稀释/供体浓度影响)
    ratio_gs  = S_fg / D_fg                    背景扣除前景比(只归一化供体，对供体浓度敏感)
    FRETN     = S_fg / sqrt(D_fg * A_fg)       Xia 归一化(浓度无关)  ★推荐主指标

前景掩膜 = 供体通道 Otsu 阈值 → 细胞区；各通道背景 = mask 外均值。
QC：D_fg < qc_threshold（默认 20）的视野（供体净信号过低）剔除。

样本命名约定（重要）：
    目录名最后一个 '_' 之后为 "基因型-光照-类型-视野"，例：
      20250101-FRET-S1_G1-D-E-1  -> genotype=G1, dark, type=E, position=1
      20250101-FRET-S1_S2-E-2    -> genotype=S2, light, type=E, position=2
      20250101-FRET-S1_G2-D-P-3  -> genotype=G2, dark, type=P, position=3
    规则：'-D' 段出现=黑暗组，无=光照组；'-E' 段=实验组，'-P' 段=阳性对照；
          最后一段纯数字=视野号；第一段=基因型。
    局限：基因型本身不得含孤立字母 'D'/'E'/'P'（会与光照/类型段混淆）。

用法：
    python analyze_fret.py --base "数据目录" [--out-dir 输出目录] [--qc 20]
"""
import os
import math
import glob
import argparse
import numpy as np
import pandas as pd
import tifffile
from scipy.stats import mannwhitneyu


def otsu_threshold(img):
    """纯 numpy Otsu 阈值，避免 skimage 依赖。"""
    img = img.astype(np.float64)
    mx = float(img.max())
    if mx <= 0:
        return 0.0
    hist, edges = np.histogram(img, bins=256, range=(0.0, mx))
    centers = (edges[:-1] + edges[1:]) / 2.0
    total = img.size
    sum_all = float((centers * hist).sum())
    w_bg = sum_bg = 0.0
    best_var, best_t = -1.0, 0.0
    for i in range(256):
        w_bg += hist[i]
        if w_bg == 0:
            continue
        w_fg = total - w_bg
        if w_fg <= 0:
            break
        sum_bg += centers[i] * hist[i]
        m_bg = sum_bg / w_bg
        m_fg = (sum_all - sum_bg) / w_fg
        v = w_bg * w_fg * (m_bg - m_fg) ** 2
        if v > best_var:
            best_var = v
            best_t = centers[i]
    return best_t


def parse_sample(name):
    """从样本目录名解析 (genotype, light, type, position)。见文件头命名约定。"""
    core = name.split("_")[-1]          # e.g. "G1-D-E-1"
    parts = core.split("-")
    genotype = parts[0]
    position = parts[-1] if parts[-1].isdigit() else ""
    light = "dark" if "D" in parts else "light"
    ctype = "E" if "E" in parts else ("P" if "P" in parts else "")
    return genotype, light, ctype, position


def find_channel(subdir, suffix):
    hits = glob.glob(os.path.join(subdir, "*" + suffix))
    return hits[0] if hits else None


def main():
    ap = argparse.ArgumentParser(description="Leica FRET batch quantification")
    ap.add_argument("--base", required=True, help="数据目录（含各样本子目录）")
    ap.add_argument("--out-dir", default=None, help="输出目录（默认=base）")
    ap.add_argument("--qc", type=float, default=20.0, help="D_fg 下限阈值，低于则剔除")
    args = ap.parse_args()

    BASE = os.path.abspath(args.base)
    OUT = os.path.abspath(args.out_dir) if args.out_dir else BASE
    os.makedirs(OUT, exist_ok=True)

    sample_dirs = [os.path.join(BASE, d) for d in sorted(os.listdir(BASE))
                   if os.path.isdir(os.path.join(BASE, d))]

    rows = []
    for d in sample_dirs:
        sample = os.path.basename(d)
        donor_f = find_channel(d, "_ch03.tif")
        fret_f = find_channel(d, "_ch00.tif")
        acc_f = find_channel(d, "_ch01.tif")
        if not (donor_f and fret_f and acc_f):
            print(f"[skip] {sample}: 缺通道文件")
            continue

        D = tifffile.imread(donor_f).astype(np.float64)
        S = tifffile.imread(fret_f).astype(np.float64)
        A = tifffile.imread(acc_f).astype(np.float64)

        thr = otsu_threshold(D)
        mask = D > thr
        if not mask.any():
            continue
        bg_D = float(D[~mask].mean())
        bg_S = float(S[~mask].mean())
        bg_A = float(A[~mask].mean())

        D_fg = float(D[mask].mean() - bg_D)
        S_fg = float(S[mask].mean() - bg_S)
        A_fg = float(A[mask].mean() - bg_A)

        ratio_raw = float(S.mean() / D.mean()) if D.mean() > 0 else np.nan
        ratio_gs = float(S_fg / D_fg) if D_fg > 0 else np.nan
        fretn = float(S_fg / np.sqrt(D_fg * A_fg)) if (D_fg > 0 and A_fg > 0) else np.nan

        g, lt, ct, pos = parse_sample(sample)
        rows.append({
            "sample": sample, "genotype": g, "light": lt, "type": ct, "position": pos,
            "bg_frac": float((~mask).mean()), "bg_D": bg_D, "bg_A": bg_A, "bg_S": bg_S,
            "D_fg": D_fg, "A_fg": A_fg, "S_fg": S_fg,
            "ratio_raw": ratio_raw, "ratio_gs": ratio_gs, "FRETN": fretn,
        })

    df = pd.DataFrame(rows)
    if df.empty:
        print("[error] 未识别到任何样本，请检查 --base 与通道后缀")
        return
    df = df.sort_values(["genotype", "light", "type", "position"]).reset_index(drop=True)

    # QC
    qc_pass = df["D_fg"] >= args.qc
    print(f"总样本 {len(df)}, D_fg<{args.qc:g} 剔除 {int((~qc_pass).sum())} 个")
    for _, r in df[~qc_pass].iterrows():
        print(f"  剔除 {r['sample']}  D_fg={r['D_fg']:.1f}")
    df_qc = df[qc_pass].reset_index(drop=True)
    print(f"QC 后 n={len(df_qc)}\n")

    # 1) 逐视野 CSV（含全部指标 + 原始 100 样本，未做 QC 过滤）
    out_fov = os.path.join(OUT, "fret_fov.csv")
    df.to_csv(out_fov, index=False, encoding="utf-8-sig")
    print("已写:", out_fov)

    # 2) 条件汇总表
    def sem(x):
        x = np.asarray(x, dtype=float)
        return x.std(ddof=1) / math.sqrt(len(x)) if len(x) > 1 else 0.0

    summary = df_qc.groupby(["genotype", "light", "type"]).agg(
        n_fov=("sample", "count"),
        n_pos=("position", "nunique"),
        FRETN_mean=("FRETN", "mean"),
        FRETN_median=("FRETN", "median"),
        FRETN_sem=("FRETN", sem),
        FRETN_std=("FRETN", "std"),
        ratio_raw_mean=("ratio_raw", "mean"),
        ratio_raw_sem=("ratio_raw", sem),
        D_fg_mean=("D_fg", "mean"),
        A_fg_mean=("A_fg", "mean"),
    ).reset_index()
    out_sum = os.path.join(OUT, "fret_condition_summary.csv")
    summary.round(4).to_csv(out_sum, index=False, encoding="utf-8-sig")
    print("已写:", out_sum)

    # 3) dark vs light 对比表（FRETN 主指标 + ratio_raw/ratio_gs 对照）
    cmp_rows = []
    for genotype in sorted(df_qc["genotype"].unique()):
        for ctype in ["E", "P"]:
            dark = df_qc[(df_qc.genotype == genotype) & (df_qc.type == ctype) & (df_qc.light == "dark")]
            light = df_qc[(df_qc.genotype == genotype) & (df_qc.type == ctype) & (df_qc.light == "light")]
            if len(dark) == 0 or len(light) == 0:
                continue
            for metric in ["FRETN", "ratio_raw", "ratio_gs"]:
                d = dark[metric].values.astype(float)
                l = light[metric].values.astype(float)
                p = mannwhitneyu(d, l, alternative="two-sided").pvalue
                fc = l.mean() / d.mean() if d.mean() != 0 else np.nan
                cmp_rows.append({
                    "genotype": genotype, "type": ctype, "metric": metric,
                    "dark_mean": d.mean(), "light_mean": l.mean(),
                    "dark_n": len(d), "light_n": len(l),
                    "fold_light_dark": fc, "mw_p": p,
                })
    cmp = pd.DataFrame(cmp_rows)
    out_cmp = os.path.join(OUT, "fret_dark_vs_light.csv")
    cmp.round(4).to_csv(out_cmp, index=False, encoding="utf-8-sig")
    print("已写:", out_cmp)

    # 终端摘要
    print("\n=== FRETN (Xia 归一化) dark vs light ===")
    f = cmp[cmp.metric == "FRETN"][["genotype", "type", "dark_mean", "light_mean",
                                   "fold_light_dark", "mw_p"]]
    print(f.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
