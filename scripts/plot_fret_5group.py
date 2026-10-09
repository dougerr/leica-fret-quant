#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FRETN (Xia 归一化) —— 多构建 × 4 条件柱状图 + 右轴 Light/Dark 折线。

左轴：FRETN 柱状图（红=光照、蓝=黑暗，深=E、浅=P），mean±SD + 灰散点。
右轴：E / P 各一条 Light/Dark fold change 点折线，y=1 虚线 = 光暗无差异。

数据源：analyze_fret.py 输出的 fret_fov.csv。
QC：默认剔除 D_fg < 20 的样本（--no-qc 保留全部）。

用法：
    python plot_fret_5group.py --csv fret_fov.csv [--builds S1,S2,S3,S4,S5] [--metric FRETN] [--out 图.png] [--no-qc]

产物落位：--out 不指定时，PNG 落在 --csv 同目录（即原始数据目录），方便原地查找。
"""
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

CONDS = ["E", "D-E", "P", "D-P"]
COL = {
    "E":   "#C0392B",
    "D-E": "#1F4E79",
    "P":   "#E6A9A0",
    "D-P": "#93B4D3",
}
COL_DOT = "#7F8C8D"


def cond(row):
    if row["light"] == "dark" and row["type"] == "E":
        return "D-E"
    if row["light"] == "light" and row["type"] == "E":
        return "E"
    if row["light"] == "dark" and row["type"] == "P":
        return "D-P"
    return "P"


def main():
    ap = argparse.ArgumentParser(description="FRETN 5组×4条件柱状图 + 右轴L/D折线")
    ap.add_argument("--csv", required=True, help="fret_fov.csv 路径")
    ap.add_argument("--builds", default=None, help="逗号分隔构建列表（默认按 CSV 自动排序）")
    ap.add_argument("--metric", default="FRETN", help="指标列名（默认 FRETN）")
    ap.add_argument("--out", default=None, help="输出 PNG 路径（默认与 --csv 同目录）")
    ap.add_argument("--qc", type=float, default=20.0, help="D_fg 下限（默认 20）")
    ap.add_argument("--no-qc", action="store_true", help="不剔除 D_fg<qc 的样本")
    args = ap.parse_args()
    if not args.out:
        # 默认落回原始数据目录（与 fret_fov.csv 同夹），方便原地查找
        args.out = os.path.join(os.path.dirname(os.path.abspath(args.csv)), "fret_5group.png")

    df = pd.read_csv(args.csv, encoding="utf-8-sig")
    if not args.no_qc and "D_fg" in df.columns:
        df = df[df["D_fg"] >= args.qc].copy()

    if args.builds:
        BUILDS = [b.strip() for b in args.builds.split(",")]
    else:
        BUILDS = sorted(df["genotype"].unique().tolist())
    METRIC = args.metric

    df["condition"] = df.apply(cond, axis=1)

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "axes.unicode_minus": False,
        "figure.dpi": 150,
        "savefig.dpi": 300,
    })

    data = {b: {c: df.loc[(df["genotype"] == b) & (df["condition"] == c), METRIC].values
                for c in CONDS} for b in BUILDS}
    mean = {b: np.array([data[b][c].mean() for c in CONDS]) for b in BUILDS}
    sd   = {b: np.array([data[b][c].std(ddof=1) for c in CONDS]) for b in BUILDS}
    n    = {b: np.array([len(data[b][c]) for c in CONDS]) for b in BUILDS}

    fc_E = [data[b]["E"].mean() / data[b]["D-E"].mean() for b in BUILDS]
    fc_P = [data[b]["P"].mean() / data[b]["D-P"].mean() for b in BUILDS]

    barw = 0.58
    group_w = 0.78
    group_gap = 1.55
    xs, centers = [], []
    x = 0.0
    for b in BUILDS:
        xb = x + np.arange(4) * group_w
        xs.append(xb)
        centers.append(xb.mean())
        x = xb[-1] + group_gap

    fig, ax = plt.subplots(figsize=(max(10.0, 2.6 * len(BUILDS)), 5.4))
    rng = np.random.default_rng(7)

    for i, b in enumerate(BUILDS):
        xb = xs[i]
        ax.bar(xb, mean[b], width=barw, color=[COL[c] for c in CONDS],
               edgecolor="black", linewidth=0.7, zorder=3, clip_on=False)
        ax.errorbar(xb, mean[b], yerr=sd[b], fmt="none", ecolor="black",
                    elinewidth=1.1, capsize=3.5, capthick=1.1, zorder=4)
        for j, c in enumerate(CONDS):
            v = data[b][c]
            jitter = rng.uniform(-0.2, 0.2, size=len(v))
            ax.scatter(xb[j] + jitter, v, s=12, color=COL_DOT, alpha=0.5,
                       edgecolor="white", linewidth=0.3, zorder=5, clip_on=False)

    y_abs_max = max((mean[b] + sd[b]).max() for b in BUILDS)
    ax.set_ylim(0, y_abs_max * 1.6)

    # 右轴：E / P 两条 Light/Dark fold change 点折线
    ax2 = ax.twinx()
    xE = [xb[0:2].mean() for xb in xs]
    xP = [xb[2:4].mean() for xb in xs]
    ax2.plot(xE, fc_E, marker="o", markersize=6, color="#1A1A1A",
             linewidth=1.6, linestyle="-", zorder=6, label="E L/D")
    ax2.plot(xP, fc_P, marker="s", markersize=6, color="#7F8C8D",
             linewidth=1.6, linestyle="--", zorder=6, label="P L/D")
    ax2.axhline(1.0, color="#B0B0B0", linestyle=":", linewidth=1.2, zorder=5)
    ax2.set_ylim(0, 1.5)
    ax2.set_yticks(np.arange(0, 1.51, 0.5))
    ax2.set_ylabel("Light / Dark ratio", fontsize=13)
    ax2.tick_params(axis="y", labelsize=11)
    ax2.spines["top"].set_visible(False)
    ax2.spines["left"].set_visible(False)

    ax.set_xticks(centers)
    ax.set_xticklabels(BUILDS, fontsize=14)
    ax.set_xlim(xs[0][0] - 0.9, xs[-1][-1] + 0.9)
    ax.set_ylabel(f"{METRIC}  (Xia normalized FRET)", fontsize=13)
    ax.tick_params(axis="y", labelsize=11)
    ax.tick_params(axis="x", length=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)

    handles = [
        Patch(fc=COL["E"], ec="black", label="Light (E)"),
        Patch(fc=COL["D-E"], ec="black", label="Dark (E)"),
        Patch(fc=COL["P"], ec="black", label="Light (P)"),
        Patch(fc=COL["D-P"], ec="black", label="Dark (P)"),
        Line2D([], [], color="#1A1A1A", marker="o", linestyle="-", label="E  L/D ratio"),
        Line2D([], [], color="#7F8C8D", marker="s", linestyle="--", label="P  L/D ratio"),
    ]
    legend = ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=10, ncol=2,
                       handlelength=1.2, handleheight=1.0, columnspacing=1.0)
    legend.get_frame().set_alpha(0.0)

    plt.tight_layout()
    plt.savefig(args.out, bbox_inches="tight", facecolor="white")
    print("Saved:", args.out)
    print("\n各组 n (E, D-E, P, D-P):")
    for b in BUILDS:
        print(f"{b:6s}", n[b].tolist())
    print("\nLight/Dark fold change:")
    for i, b in enumerate(BUILDS):
        print(f"{b:6s} E(L/D)={fc_E[i]:.3f}   P(L/D)={fc_P[i]:.3f}")


if __name__ == "__main__":
    main()
