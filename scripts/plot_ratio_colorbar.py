# -*- coding: utf-8 -*-
"""
FRET Ratio 伪彩图的色阶图例（LUT colorbar）。

与 per_fov_imaging.py 共用同一套 colormap 与固定标尺（默认 0–0.6），
产出的图例可直接放在 per_fov 伪彩图旁边当 legend 用。

标签用英文 + Times New Roman（出版级，且避免 matplotlib 中文字体缺失导致乱码）。

用法：
  python plot_ratio_colorbar.py [--ratio-max 0.6] [--out fret_ratio_colorbar.png]
      [--dpi 300] [--formats png,tif] [--no-title]
"""
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap

# ★ 必须与 per_fov_imaging.py 的 FRET_CMAP 完全一致（逐像素采样锚点）
FRET_CMAP = LinearSegmentedColormap.from_list("fret", [
    (0.00, (11/255, 226/255, 251/255)),
    (0.25, (10/255, 249/255, 210/255)),
    (0.38, (11/255, 249/255, 113/255)),
    (0.50, (8/255, 249/255, 10/255)),
    (0.62, (55/255, 249/255, 9/255)),
    (0.75, (157/255, 251/255, 13/255)),
    (0.88, (250/255, 192/255, 12/255)),
    (1.00, (245/255, 112/255, 5/255)),
])


def main():
    ap = argparse.ArgumentParser(description="FRET Ratio 色阶图例")
    ap.add_argument("--ratio-max", type=float, default=0.6, help="标尺上限（与伪彩图一致）")
    ap.add_argument("--out", default="fret_ratio_colorbar.png", help="输出路径（默认当前目录）")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--formats", default="png", help="逗号分隔：png,tif")
    ap.add_argument("--step", type=float, default=0.1, help="刻度间隔")
    ap.add_argument("--orientation", default="horizontal", choices=["horizontal", "vertical"],
                    help="horizontal=横版带刻度；vertical=竖版紧凑样式（裸色条，无边框无刻度无文字）")
    ap.add_argument("--ticks", action="store_true",
                    help="竖版时附加刻度数字（原始图版没有，默认关）")
    ap.add_argument("--bg", default="white", choices=["white", "black"],
                    help="背景色：black 用于直接拼到黑底伪彩图旁")
    ap.add_argument("--no-title", action="store_true", help="不画标题")
    ap.add_argument("--no-caption", action="store_true", help="不画底部颜色说明")
    args = ap.parse_args()

    R = args.ratio_max
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "axes.unicode_minus": False,
    })

    bg = "white" if args.bg == "white" else "black"
    fg = "black" if args.bg == "white" else "white"
    low_c = "#0F6E56" if args.bg == "white" else "#5DCAA5"
    high_c = "#993C1D" if args.bg == "white" else "#F0997B"

    if args.orientation == "vertical":
        # 实测样式：裸竖条，宽高比 ≈ 1:6，上=高(橙红) 下=低(青蓝)，
        # 无边框、无刻度线、无数字（--ticks 可加）。
        w = 0.95 if args.ticks else 0.50
        h = 2.20 if args.ticks else 2.10
        fig = plt.figure(figsize=(w, h), dpi=args.dpi, facecolor=bg)
        if args.ticks:
            ax = fig.add_axes([0.10, 0.05, 0.26, 0.90])
        else:
            ax = fig.add_axes([0.22, 0.05, 0.56, 0.90])
        grad = np.linspace(0.0, 1.0, 1024).reshape(-1, 1)
        ax.imshow(grad, aspect="auto", cmap=FRET_CMAP,
                  extent=[0, 1, 0, R], origin="lower", zorder=2)
        ax.axis("off")
        if args.ticks:
            for v in np.arange(0.0, R + 1e-9, args.step):
                ax.plot([1.0, 1.22], [v, v], color=fg, lw=0.9, zorder=3)
                ax.text(1.34, v, f"{v:.1f}", ha="left", va="center",
                        fontsize=9, color=fg)
    else:
        fig = plt.figure(figsize=(6.4, 1.55), dpi=args.dpi, facecolor=bg)
        ax = fig.add_axes([0.045, 0.06, 0.91, 0.80])
        ax.set_xlim(-0.015 * R / 0.6, R * 1.015)
        ax.set_ylim(0.0, 1.30)
        ax.axis("off")

        # --- 色带本体 ---
        y0, y1 = 0.58, 1.06
        grad = np.linspace(0.0, 1.0, 1024).reshape(1, -1)
        ax.imshow(grad, aspect="auto", cmap=FRET_CMAP, extent=[0.0, R, y0, y1],
                  interpolation="bilinear", zorder=2)
        ax.add_patch(Rectangle((0.0, y0), R, y1 - y0, fill=False,
                               edgecolor=fg, linewidth=0.9, zorder=3))

        # --- 刻度线 + 数值 ---
        ticks = np.arange(0.0, R + 1e-9, args.step)
        for v in ticks:
            ax.plot([v, v], [y0 - 0.14, y0], color=fg, lw=0.9, zorder=3)
            ax.text(v, y0 - 0.30, f"{v:.1f}", ha="center", va="center",
                    fontsize=10, color=fg)

        if not args.no_title:
            ax.text(R / 2.0, 1.24, "FRET Ratio  (S/D)", ha="center", va="center",
                    fontsize=12, color=fg)

        if not args.no_caption:
            ax.text(0.0, 0.06, "low FRET", ha="left", va="center",
                    fontsize=10, color=low_c)
            ax.text(R, 0.06, "high FRET", ha="right", va="center",
                    fontsize=10, color=high_c)

    out = os.path.abspath(args.out)
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    formats = [f.strip().lower() for f in args.formats.split(",") if f.strip()]
    for fmt in formats:
        p = out if os.path.splitext(out)[1].lower() == "." + fmt \
            else os.path.splitext(out)[0] + "." + fmt
        fig.savefig(p, dpi=args.dpi, facecolor=bg, bbox_inches="tight")
        print("Saved:", p)
    plt.close(fig)


if __name__ == "__main__":
    main()
