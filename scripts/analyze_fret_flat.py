# -*- coding: utf-8 -*-
"""
平铺目录版 FRET 批量分析（analyze_fret.py 的平铺布局变体）。

适用：所有视野所有通道混在同一个文件夹（LAS X 直接导出、未按样本分夹）。
analyze_fret.py 则适配「每个样本一个子目录」布局，二者指标口径一致。

文件命名（两种形式自动识别）：
  A) <基因型>-<D?>-<E|P>_Position<NNN>_ch<CC>.tif
     S1-D-E_Position001_ch00.tif   → 基因型=S1, dark, E, pos=001, ch=FRET
     S1-E_Position001_ch00.tif     → 基因型=S1, light, E, pos=001, ch=FRET
  B) <共同前缀>_<基因型>-<D?>-<E|P>-<视野号>_ch<CC>.tif
     20250101-FRET-S1_2X-D-E-1_ch00.tif → 基因型=2X, dark, E, pos=1, ch=FRET

通道映射（LAS X 0-based）：
  ch00 = FRET 敏化发射   ch01 = 受体直接激发    ch03 = 供体 BFP    ch02 = 明场(弃用)

指标：ratio_raw / ratio_gs / FRETN（Xia 归一化，主指标）
QC：D_fg < 阈值（默认 20）剔除。
产出：fret_fov.csv / fret_condition_summary.csv / fret_dark_vs_light.csv
（与 plot_fret_5group.py 直接对接。）
"""
import os, re, math, argparse
from collections import defaultdict
import numpy as np
import pandas as pd
import tifffile
from scipy.stats import mannwhitneyu


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


# 兼容两类平铺命名（自动识别）：
#   A) <前缀>_Position<NNN>_ch<CC>.tif               如 S1-D-E_Position001_ch00.tif
#   B) <共同前缀>_<基因型>-<D?>-<E|P>-<视野号>_ch<CC>.tif
#                                                     如 20250101-FRET-S1_2X-D-E-1_ch00.tif
FNAME_RE = re.compile(r"^(?P<prefix>.+)_ch(?P<ch>\d{2})\.tif[f]?$", re.IGNORECASE)


def parse_fname(path):
    """返回 dict: prefix / ch / pos，另附 sample（去通道后缀的完整样本名）。"""
    base = os.path.basename(path)
    m = FNAME_RE.match(base)
    if not m:
        return None
    d = m.groupdict()                 # prefix, ch
    raw = d["prefix"]
    mpos = re.search(r"_Position(\d+)$", raw, re.IGNORECASE)
    if mpos:                          # 形式 A
        d["pos"] = mpos.group(1)
        d["prefix"] = raw[:mpos.start()]
        d["sample"] = f"{d['prefix']}_Position{d['pos']}"
    else:                             # 形式 B：末段为 <基因型>-<D?>-<E|P>-<视野号>
        core = raw.split("_")[-1]
        parts = core.split("-")
        d["pos"] = parts[-1] if parts[-1].isdigit() else ""
        d["prefix"] = raw
        d["sample"] = raw             # 形式 B 的 prefix 本身就是完整样本名
    return d


def parse_condition(prefix):
    """prefix 例: 'S1-D-E'（A）或 '20250101-FRET-S1_2X-D-E'（B）。

    取最后一段 '_' 之后的 core 再按 '-' 切分，两种形式通用。
    """
    core = prefix.split("_")[-1]
    parts = core.split("-")
    genotype = parts[0]
    # 'D' 必须独立成段才视为 dark，避免误判基因型内的 D
    light = "dark" if "D" in parts[1:] else "light"
    ctype = "E" if "E" in parts[1:] else ("P" if "P" in parts[1:] else "")
    return genotype, light, ctype


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, help="平铺目录（所有 .tif 在同一文件夹）")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--qc", type=float, default=20.0)
    args = ap.parse_args()

    BASE = os.path.abspath(args.base)
    OUT = os.path.abspath(args.out_dir) if args.out_dir else BASE
    os.makedirs(OUT, exist_ok=True)

    # 收集所有 .tif，按样本分组到通道字典
    groups = defaultdict(dict)  # key=(sample, prefix, pos) → {ch: filepath}
    for fn in sorted(os.listdir(BASE)):
        if not fn.lower().endswith((".tif", ".tiff")):
            continue
        full = os.path.join(BASE, fn)
        info = parse_fname(full)
        if not info:
            print(f"[skip non-matching] {fn}")
            continue
        key = (info["sample"], info["prefix"], info["pos"])
        groups[key][info["ch"]] = full

    rows = []
    for (sample, prefix, pos), chmap in sorted(groups.items()):
        donor_f = chmap.get("03")
        fret_f = chmap.get("00")
        acc_f = chmap.get("01")
        if not (donor_f and fret_f and acc_f):
            missing = [c for c in ("00", "01", "03") if c not in chmap]
            print(f"[skip incomplete] {sample}  缺 ch{missing}")
            continue

        D = tifffile.imread(donor_f).astype(np.float64)
        S = tifffile.imread(fret_f).astype(np.float64)
        A = tifffile.imread(acc_f).astype(np.float64)

        thr = otsu_threshold(D)
        mask = D > thr
        if not mask.any():
            print(f"[skip empty mask] {sample}")
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

        g, lt, ct = parse_condition(prefix)
        rows.append({
            "sample": sample, "prefix": prefix, "genotype": g, "light": lt, "type": ct, "position": pos,
            "bg_frac": float((~mask).mean()), "bg_D": bg_D, "bg_A": bg_A, "bg_S": bg_S,
            "D_fg": D_fg, "A_fg": A_fg, "S_fg": S_fg,
            "ratio_raw": ratio_raw, "ratio_gs": ratio_gs, "FRETN": fretn,
        })

    df = pd.DataFrame(rows)
    if df.empty:
        print("[error] 未识别到任何样本")
        return
    df = df.sort_values(["genotype", "light", "type", "position"]).reset_index(drop=True)

    qc_pass = df["D_fg"] >= args.qc
    print(f"\n总样本 {len(df)}, D_fg<{args.qc:g} 剔除 {int((~qc_pass).sum())} 个")
    for _, r in df[~qc_pass].iterrows():
        print(f"  剔除 {r['sample']}  D_fg={r['D_fg']:.1f}")
    df_qc = df[qc_pass].reset_index(drop=True)
    print(f"QC 后 n={len(df_qc)}\n")

    # 1) 逐视野 CSV
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

    # 3) dark vs light 对比
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
