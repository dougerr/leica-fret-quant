# -*- coding: utf-8 -*-
"""
逐细胞 FRET 定量（平铺布局版）。

布局：所有样本的 4 通道 TIF 平铺在同一目录：
    <prefix>_ch00.tif = FRET 敏化发射
    <prefix>_ch01.tif = 受体 (acceptor)
    <prefix>_ch02.tif = 明场(弃用)
    <prefix>_ch03.tif = 供体 BFP
样本命名（两种平铺命名自动识别）：
    A) <基因型>-<D?>-<E|P>_Position<NNN>       如 S1-D-E_Position001
    B) <共同前缀>_<基因型>-<D?>-<E|P>-<视野号>   如 20250101-FRET-S1_2X-D-E-1
    '-D'=黑暗(无=光照)，'-E'=实验组，'-P'=阳性对照，末段数字=视野号。

逐细胞流程：
    1. 供体 Otsu → 掩膜 → 填孔 → 连通域标记
    2. 面积过滤（min_area~max_area），每个连通域=一个细胞
    3. 每细胞在 D/S/A 三通道分别取掩膜内均值并扣各自背景（掩膜外均值）
    4. 指标：ratio_raw=S/D、ratio_gs=S_fg/D_fg、FRETN=S_fg/sqrt(D_fg*A_fg)

用法：
    python fret_cells_flat.py --base "数据目录" [--out-dir 输出目录]
                              [--min-area 300] [--max-area 20000]
"""
import os
import re
import glob
import math
import argparse
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage


def otsu_threshold(img):
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


def parse_sample(prefix):
    """兼容两种平铺命名，返回 (基因型, light|dark, E|P|'', 视野号)。

    A) <基因型>-<D?>-<E|P>_Position<NNN>       如 S1-D-E_Position001
    B) <共同前缀>_<基因型>-<D?>-<E|P>-<视野号>   如 20250101-FRET-S1_2X-D-E-1
    """
    mpos = re.search(r"_Position(\d+)$", prefix, re.IGNORECASE)
    if mpos:                                   # 形式 A
        core = prefix[:mpos.start()]
        position = mpos.group(1)
    else:                                      # 形式 B
        core = prefix.split("_")[-1]
        bits = core.split("-")
        position = bits[-1] if bits[-1].isdigit() else ""
    parts = core.split("-")
    genotype = parts[0]
    # 'D' 必须独立成段才视为 dark，避免误判基因型内的 D
    light = "dark" if "D" in parts[1:] else "light"
    ctype = "E" if "E" in parts[1:] else ("P" if "P" in parts[1:] else "")
    return genotype, light, ctype, position


def segment(donor_bgsub, min_area, max_area):
    mask = donor_bgsub > otsu_threshold(donor_bgsub)
    mask = ndimage.binary_fill_holes(mask)
    labels, n = ndimage.label(mask)
    if n == 0:
        return labels
    areas = ndimage.sum(np.ones_like(labels), labels, np.arange(1, n + 1))
    keep_ids = np.nonzero((areas >= min_area) & (areas <= max_area))[0] + 1
    remap = np.zeros(n + 1, dtype=int)
    for new_id, old in enumerate(keep_ids, start=1):
        remap[old] = new_id
    return remap[labels]


def main():
    ap = argparse.ArgumentParser(description="逐细胞 FRET 定量（平铺布局）")
    ap.add_argument("--base", required=True)
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--min-area", type=int, default=300)
    ap.add_argument("--max-area", type=int, default=20000)
    args = ap.parse_args()

    BASE = os.path.abspath(args.base)
    OUT = os.path.abspath(args.out_dir) if args.out_dir else os.path.join(BASE, "fret_cells_results")
    os.makedirs(OUT, exist_ok=True)

    prefixes = sorted({os.path.basename(f)[:-9]
                       for f in glob.glob(os.path.join(BASE, "*_ch03.tif"))})
    print(f"样本数: {len(prefixes)}   面积过滤: {args.min_area}–{args.max_area} px")

    per_cell, per_fov = [], []
    for prefix in prefixes:
        d_f = os.path.join(BASE, prefix + "_ch03.tif")
        s_f = os.path.join(BASE, prefix + "_ch00.tif")
        a_f = os.path.join(BASE, prefix + "_ch01.tif")
        if not (os.path.exists(d_f) and os.path.exists(s_f) and os.path.exists(a_f)):
            print(f"[skip] {prefix}: 缺通道文件")
            continue

        D = tifffile.imread(d_f).astype(np.float64)
        S = tifffile.imread(s_f).astype(np.float64)
        A = tifffile.imread(a_f).astype(np.float64)

        # 全局背景：整图 Otsu 掩膜外均值（逐细胞扣减在此基础上再做细胞外均值）
        m_all = D > otsu_threshold(D)
        bgD, bgS, bgA = D[~m_all].mean(), S[~m_all].mean(), A[~m_all].mean()
        Dn = np.clip(D - bgD, 0, None)
        Sn = np.clip(S - bgS, 0, None)
        An = np.clip(A - bgA, 0, None)

        labels = segment(Dn, args.min_area, args.max_area)
        ids = np.unique(labels[labels > 0])
        if ids.size == 0:
            print(f"  {prefix}: 未分割出细胞（面积过滤 {args.min_area}–{args.max_area}）")
            continue

        for cid in ids:
            m = labels == cid
            area = int(m.sum())
            d_fg = float(Dn[m].mean())
            s_fg = float(Sn[m].mean())
            a_fg = float(An[m].mean())
            g, lt, ct, pos = parse_sample(prefix)
            per_cell.append({
                "sample": prefix, "cell_id": int(cid), "area": area,
                "genotype": g, "light": lt, "type": ct, "position": pos,
                "D_fg": d_fg, "S_fg": s_fg, "A_fg": a_fg,
                "ratio_raw_cell": s_fg / d_fg if d_fg > 0 else np.nan,
                "ratio_gs_cell": s_fg / d_fg if d_fg > 0 else np.nan,
                "FRETN_cell": s_fg / math.sqrt(d_fg * a_fg) if (d_fg > 0 and a_fg > 0) else np.nan,
            })

        dfv = pd.DataFrame([r for r in per_cell if r["sample"] == prefix])
        per_fov.append({
            "sample": prefix, "genotype": g, "light": lt, "type": ct, "position": pos,
            "n_cells": len(ids),
            "FRETN_median": dfv["FRETN_cell"].median(),
            "FRETN_mean": dfv["FRETN_cell"].mean(),
            "ratio_gs_median": dfv["ratio_gs_cell"].median(),
        })
        print(f"  {prefix}: {len(ids)} cells  FRETN median={dfv['FRETN_cell'].median():.3f}")

    if not per_cell:
        print("[error] 无任何细胞被分割")
        return

    df = pd.DataFrame(per_cell)
    df.to_csv(os.path.join(OUT, "fret_per_cell.csv"), index=False, encoding="utf-8-sig")
    pd.DataFrame(per_fov).to_csv(os.path.join(OUT, "fret_per_fov.csv"), index=False,
                                 encoding="utf-8-sig")
    print(f"\n已写: {os.path.join(OUT, 'fret_per_cell.csv')}  ({len(df)} cells)")
    print(f"已写: {os.path.join(OUT, 'fret_per_fov.csv')}  ({len(per_fov)} FOVs)")

    # 分组统计：细胞级（主） + FOV 级（对照）
    from scipy.stats import mannwhitneyu
    rows = []
    for g in sorted(df["genotype"].unique()):
        for ct in ["E", "P"]:
            d = df[(df.genotype == g) & (df.type == ct) & (df.light == "dark")]["FRETN_cell"].dropna().values
            l = df[(df.genotype == g) & (df.type == ct) & (df.light == "light")]["FRETN_cell"].dropna().values
            if len(d) == 0 or len(l) == 0:
                continue
            p = mannwhitneyu(d, l, alternative="two-sided").pvalue
            rows.append({"genotype": g, "type": ct, "level": "cell",
                         "dark_n": len(d), "light_n": len(l),
                         "dark_mean": d.mean(), "light_mean": l.mean(),
                         "fold_light_dark": l.mean() / d.mean() if d.mean() else np.nan,
                         "mw_p": p})
    for g in sorted(pd.DataFrame(per_fov)["genotype"].unique()):
        for ct in ["E", "P"]:
            d = pd.DataFrame(per_fov)[(pd.DataFrame(per_fov).genotype == g) & (pd.DataFrame(per_fov).type == ct) & (pd.DataFrame(per_fov).light == "dark")]["FRETN_median"].dropna().values
            l = pd.DataFrame(per_fov)[(pd.DataFrame(per_fov).genotype == g) & (pd.DataFrame(per_fov).type == ct) & (pd.DataFrame(per_fov).light == "light")]["FRETN_median"].dropna().values
            if len(d) == 0 or len(l) == 0:
                continue
            p = mannwhitneyu(d, l, alternative="two-sided").pvalue
            rows.append({"genotype": g, "type": ct, "level": "fov",
                         "dark_n": len(d), "light_n": len(l),
                         "dark_mean": d.mean(), "light_mean": l.mean(),
                         "fold_light_dark": l.mean() / d.mean() if d.mean() else np.nan,
                         "mw_p": p})
    cmp = pd.DataFrame(rows)
    out_cmp = os.path.join(OUT, "fret_cells_dark_vs_light.csv")
    cmp.round(4).to_csv(out_cmp, index=False, encoding="utf-8-sig")
    print(f"已写: {out_cmp}")

    print("\n=== 细胞级 FRETN dark vs light ===")
    print(cmp[cmp.level == "cell"].round(4).to_string(index=False))
    print("\n=== FOV 级 FRETN median dark vs light（对照） ===")
    print(cmp[cmp.level == "fov"].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
