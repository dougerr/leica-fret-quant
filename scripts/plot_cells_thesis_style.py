#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
论文版式逐细胞 FRET 统计图（灰柱+黑点 / 空心圆散点，论文常见样式）。

输入：analyze_cells_flat.py 产出的 fret_per_cell.csv 与 fret_per_fov.csv
输出：左 FOV 级灰柱+黑点图 + 右 细胞级空心圆散点图

用法：
  python plot_cells_thesis_style.py --per-cell fret_per_cell.csv --per-fov fret_per_fov.csv \
      [--out fret_cells_thesis_style.png] [--metric FRETN_cell]

产物落位：--out 不指定时，PNG 落在 --per-cell 同目录（即原始数据目录），方便原地查找。
"""
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- 论文样式常量（从 PDF 逐像素采样） ----
RED = "#E94A50"      # 红（光照/实验）RGB≈(233,74,80)
BLUE = "#3938D5"     # 蓝（黑暗/对照）RGB≈(57,56,213)
GRAY_BAR = "#96969B" # 柱状图灰填充 RGB≈(150,150,155)
LINE_SD = "#4D4D4D"  # mean±SD 灰线
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "axes.unicode_minus": False,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-cell", required=True)
    ap.add_argument("--per-fov", required=True)
    ap.add_argument("--out", default=None, help="输出 PNG 路径（默认与 --per-cell 同目录，即原始数据目录）")
    ap.add_argument("--metric", default="FRETN_cell")
    ap.add_argument("--fov-metric", default="FRETN_median")
    args = ap.parse_args()
    if not args.out:
        # 默认落回原始数据目录（与 fret_per_cell.csv 同夹），方便原地查找
        args.out = os.path.join(os.path.dirname(os.path.abspath(args.per_cell)),
                                "fret_cells_thesis_style.png")

    dfc = pd.read_csv(args.per_cell, encoding="utf-8-sig")
    dff = pd.read_csv(args.per_fov, encoding="utf-8-sig")
    metric, fov_metric = args.metric, args.fov_metric

    GENOS = sorted(dfc["genotype"].unique())
    TYPES = ["E", "P"]

    # 图宽/字号随组数自适应：2 构建时与原版一致，4+ 构建时横向展宽避免 x 轴标签重叠
    n_fov_groups = len(GENOS) * len(TYPES)          # 左图柱数
    n_cell_ticks = n_fov_groups * 2                 # 右图散点列数
    wL = max(3.6, 0.95 * n_fov_groups)
    wR = max(5.6, 0.85 * n_cell_ticks)   # 右图每列 ~0.85in，保证 "G-E\nlight" 两行标签不重叠
    fs_l = 11 if n_fov_groups <= 5 else 8
    fs_r = 9 if n_cell_ticks <= 10 else 7

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(wL + wR, 4.6),
                                   gridspec_kw={"width_ratios": [wL, wR]})

    # ---- 左：FOV 级灰柱 + 黑点 ----
    xt, xl = [], []
    x = 0.0
    for g in GENOS:
        for t in TYPES:
            v = dff[(dff.genotype == g) & (dff.type == t)][fov_metric].dropna().values
            axL.bar(x, v.mean(), width=0.72, color=GRAY_BAR, edgecolor="black", linewidth=1.0)
            axL.errorbar(x, v.mean(), yerr=v.std(ddof=1) if len(v) > 1 else 0,
                         fmt="none", ecolor="black", elinewidth=1.2, capsize=4, capthick=1.2)
            axL.scatter(np.full(len(v), x), v, s=22, color="black", zorder=5)
            xt.append(x); xl.append(f"{g}-{t}"); x += 1.0
    axL.set_xticks(xt); axL.set_xticklabels(xl, fontsize=fs_l)
    axL.set_ylabel("Ratio  610 nm/450 nm", fontsize=12)
    vmax_l = dff[fov_metric].dropna()
    axL.set_ylim(0, max(0.02, vmax_l.max() * 1.25))
    axL.set_title("FOV-level", fontsize=12)
    axL.spines["top"].set_visible(False); axL.spines["right"].set_visible(False)
    axL.tick_params(labelsize=10)

    # ---- 右：细胞级空心圆散点 + mean±SD 灰线 ----
    rng = np.random.default_rng(11)
    x = 0.0
    xt, xl = [], []
    for g in GENOS:
        for t in TYPES:
            for lt, col in [("light", RED), ("dark", BLUE)]:
                v = dfc[(dfc.genotype == g) & (dfc.type == t) & (dfc.light == lt)][metric].dropna().values
                jit = rng.uniform(-0.22, 0.22, size=len(v))
                axR.scatter(jit + x, v, s=26, facecolors="white", edgecolors=col,
                            linewidths=0.9, alpha=0.95, zorder=4)
                m, sd = v.mean(), v.std(ddof=1)
                axR.hlines(m, x - 0.28, x + 0.28, color=LINE_SD, lw=1.3, zorder=6)
                axR.vlines(x, m - sd, m + sd, color=LINE_SD, lw=1.3, zorder=6)
                axR.hlines(m - sd, x - 0.16, x + 0.16, color=LINE_SD, lw=1.1, zorder=6)
                axR.hlines(m + sd, x - 0.16, x + 0.16, color=LINE_SD, lw=1.1, zorder=6)
                xt.append(x); xl.append(f"{g}-{t}\n{lt}"); x += 1.0
            x += 0.35
    axR.set_xticks(xt); axR.set_xticklabels(xl, fontsize=fs_r)
    axR.set_ylabel("FRETN  (610 nm/450 nm)", fontsize=12)
    axR.set_title("Cell-level", fontsize=12)
    axR.set_xlim(-0.7, x - 0.35)
    # Y 轴上限用 98 分位自适应，防高离群点压缩主体
    ymax_all = dfc[metric].quantile(0.98)
    axR.set_ylim(0, max(0.30, ymax_all * 1.15))
    axR.spines["top"].set_visible(False); axR.spines["right"].set_visible(False)
    axR.tick_params(labelsize=10)

    plt.tight_layout()
    plt.savefig(args.out, bbox_inches="tight", facecolor="white")
    print("Saved:", args.out)


if __name__ == "__main__":
    main()
