# -*- coding: utf-8 -*-
"""
每视野伪彩成像图 + 原始灰度图打包（论文同款观感）。

对平铺目录（所有视野所有通道混在一个文件夹）按视野生成：

  <base>/per_fov/<prefix>/
    <prefix>_ch00.tif           ← 原始 FRET 敏化发射（复制）
    <prefix>_ch01.tif           ← 原始 Acceptor（复制）
    <prefix>_ch03.tif           ← 原始 Donor（复制）
    <prefix>_Acceptor.tif       ← 伪彩 受体·红（单页 RGB TIF）
    <prefix>_Donor.tif          ← 伪彩 供体·蓝（单页 RGB TIF）
    <prefix>_FRETRatio.tif      ← FRET Ratio 渐变 固定 0–0.6（单页 RGB TIF）

默认再把整个 per_fov 目录打包为 <base>/per_fov.zip（--no-zip 跳过）。

★ 算法口径（2026-09 定稿，多轮实测确认："终于对了"）——改参数前先留回滚点：

  Ratio 图 = 基准通道灰度（明度/纹理） × FRET 色阶颜色（色相） × brightness

  1) 形态与质感完全一致：Ratio 图与基准通道伪彩图**逐像素共用同一张灰度图**
     （同一 view_norm 输出、同一 sigma、同一 gamma），因此轮廓、细胞内部纹理、
     亮度层次与 Acceptor/Donor 伪彩图完全一样，**只有颜色**反映 FRET 比率。
     ★ 关键教训：绝不要用硬掩膜或阈值化权重（clip((x-lo)/soft)**g）充当明度——
     那会把细胞内部纹理压平成色块，看起来"失真"（v1–v4 全部踩过这个坑）。
  2) 基准通道默认 donor（--base-channel acceptor 可切）：受体亮/供体弱的细胞
     D≈背景时 S/D 被除法放大出假高值（实测虚高至 0.76，正常细胞仅 0.07），
     以 donor 作明度可自动压暗这类不可信区域；如需与受体行完全同貌选 acceptor。
  3) 色阶固定 0–0.6（论文常用标尺，不随视野自适应）；--brightness 控整体明暗
     （默认 0.65，实测 0.8 偏亮）。
  4) Ratio 用原始 S/D（LAS X 口径，不做背景扣除）——与论文图的呈像口径一致。
     定量结论一律走 CSV 的 FRETN，不要从这张伪彩图读数。

命名兼容（两种批次都支持）：
  A) <批次前缀>_<基因型>-<D?>-<E|P>-<视野号>_chCC.tif
  B) <基因型>-<D?>-<E|P>_Position<NNN>_chCC.tif

用法：
  python per_fov_imaging.py --base "数据目录" [--ratio-max 0.6] [--sigma 1.0]
      [--base-channel donor] [--brightness 0.65] [--no-zip] [--preview png]
"""
import os
import re
import shutil
import argparse
import zipfile
from collections import defaultdict
import numpy as np
import tifffile
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter

# 从论文图采样的锚点：0=青 → 1=红
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

ACCEPTOR_RGB = (1.00, 0.18, 0.18)   # 红
DONOR_RGB = (0.28, 0.60, 1.00)      # 蓝

# prefix = 文件名去掉 _chCC.tif 后缀的全部内容（两种命名通用）
FNAME_RE = re.compile(r"^(?P<prefix>.+)_ch(?P<ch>\d{2})\.tif[f]?$", re.IGNORECASE)


def parse_sample(prefix):
    """兼容两种命名，返回 (基因型, light|dark, E|P|'', 视野号)。

    A) <批次前缀>_<基因型>-<D?>-<E|P>-<视野号>
    B) <基因型>-<D?>-<E|P>_Position<NNN>
    """
    mpos = re.search(r"_Position(\d+)$", prefix, re.IGNORECASE)
    if mpos:                                    # 形式 B
        pos = mpos.group(1)
        core = prefix[:mpos.start()]
    else:                                       # 形式 A
        core = prefix.split("_")[-1]
        parts = core.split("-")
        pos = parts[-1] if parts[-1].isdigit() else ""
    parts = core.split("-")
    g = parts[0] if parts else core
    light = "dark" if "D" in parts[1:] else "light"
    ctype = "E" if "E" in parts[1:] else ("P" if "P" in parts[1:] else "")
    return g, light, ctype, pos


def view_norm(x, p_lo=0.5, p_hi=99.8):
    """视野内稳健归一化（分位数拉伸），颜色亮度 ∝ 信号强度。"""
    lo = float(np.percentile(x, p_lo))
    hi = float(np.percentile(x, p_hi))
    if hi <= lo + 1e-6:
        return np.zeros_like(x, dtype=np.float32)
    return np.clip((x - lo) / (hi - lo), 0, 1).astype(np.float32)


def as_rgb(gray, color):
    rgb = np.zeros(gray.shape + (3,), dtype=np.float32)
    for k in range(3):
        rgb[..., k] = gray * color[k]
    return rgb


def write_colorbar(path, ratio_max, dpi=300):
    """在 per_fov 根目录写一张色阶图例（与伪彩图同一 LUT 与标尺）。

    标签用英文 + Times New Roman，避免 matplotlib 中文字体缺失导致乱码。
    独立生成请用 scripts/plot_ratio_colorbar.py（参数更全）。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    R = ratio_max
    fig = plt.figure(figsize=(6.4, 1.55), dpi=dpi, facecolor="white")
    ax = fig.add_axes([0.045, 0.06, 0.91, 0.80])
    ax.set_xlim(-0.015 * R / 0.6, R * 1.015)
    ax.set_ylim(0.0, 1.30)
    ax.axis("off")
    y0, y1 = 0.58, 1.06
    grad = np.linspace(0.0, 1.0, 1024).reshape(1, -1)
    ax.imshow(grad, aspect="auto", cmap=FRET_CMAP, extent=[0.0, R, y0, y1],
              interpolation="bilinear", zorder=2)
    ax.add_patch(Rectangle((0.0, y0), R, y1 - y0, fill=False,
                           edgecolor="black", linewidth=0.9, zorder=3))
    for v in np.arange(0.0, R + 1e-9, 0.1):
        ax.plot([v, v], [y0 - 0.14, y0], color="black", lw=0.9, zorder=3)
        ax.text(v, y0 - 0.30, f"{v:.1f}", ha="center", va="center",
                fontsize=10, color="black")
    ax.text(R / 2.0, 1.24, "FRET Ratio  (S/D)", ha="center", va="center",
            fontsize=12, color="black")
    ax.text(0.0, 0.06, "low FRET", ha="left", va="center",
            fontsize=10, color="#0F6E56")
    ax.text(R, 0.06, "high FRET", ha="right", va="center",
            fontsize=10, color="#993C1D")
    fig.savefig(path, dpi=dpi, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print("colorbar:", path)


def parse_rename(spec):
    """解析 --rename，格式 "OLD=NEW,OLD2=NEW2"；仅替换 TIFF description 里的标注。"""
    if not spec:
        return {}
    out = {}
    for pair in spec.split(","):
        pair = pair.strip()
        if not pair:
            continue
        if "=" not in pair:
            raise SystemExit(f"[error] --rename 格式应为 OLD=NEW，收到: {pair!r}")
        old, new = pair.split("=", 1)
        out[old.strip()] = new.strip()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, help="平铺数据目录（所有 *_chCC.tif 所在文件夹）")
    ap.add_argument("--out-root", default="per_fov", help="相对 base 下的输出根目录")
    ap.add_argument("--zip-name", default="per_fov.zip", help="相对 base 下的 zip 文件名")
    ap.add_argument("--no-zip", action="store_true", help="只生成 per_fov 目录，不打包 zip")
    ap.add_argument("--rename", default="", help='标注替换，如 "S3=S1"（仅改 TIFF description 标注）')
    ap.add_argument("--ratio-max", type=float, default=0.6, help="FRET Ratio 固定标尺上限（论文同款 0.6）")
    ap.add_argument("--sigma", type=float, default=1.0, help="三通道平滑强度（0=不平滑）")
    ap.add_argument("--base-channel", default="donor", choices=["donor", "acceptor"],
                    help="Ratio 图明度基准通道：donor（默认，无除法伪影）/ acceptor（与受体行同貌）")
    ap.add_argument("--disp-gamma", type=float, default=1.0, help="受体/供体显示 gamma（<1 提亮，>1 压暗）")
    ap.add_argument("--w-gamma", type=float, default=1.0, help="基准通道明度 gamma（1.0=原样复制质感）")
    ap.add_argument("--brightness", type=float, default=0.65, help="Ratio 图整体亮度系数")
    ap.add_argument("--preview", default="", help="可选：输出预览图 png 路径（仅供快速查看）")
    ap.add_argument("--colorbar", action="store_true",
                    help="在 per_fov 根目录额外生成 _colorbar.png 色阶图例")
    args = ap.parse_args()

    BASE = os.path.abspath(args.base)
    OUT_ROOT = os.path.join(BASE, args.out_root)
    os.makedirs(OUT_ROOT, exist_ok=True)
    RENAME = parse_rename(args.rename)

    groups = defaultdict(dict)
    for fn in sorted(os.listdir(BASE)):
        if not fn.lower().endswith((".tif", ".tiff")):
            continue
        m = FNAME_RE.match(fn)
        if m:
            groups[m["prefix"]][int(m["ch"])] = os.path.join(BASE, fn)

    n_ok = n_skip = 0
    skipped, preview_items = [], []

    for prefix, chm in sorted(groups.items()):
        if not all(c in chm for c in (0, 1, 3)):
            skipped.append((prefix, "缺通道"))
            n_skip += 1
            continue

        g, light, ctype, pos = parse_sample(prefix)
        gname = RENAME.get(g, g)
        fov_dir = os.path.join(OUT_ROOT, prefix)
        os.makedirs(fov_dir, exist_ok=True)

        # 1) 原始灰度三通道：复制
        for c in (0, 1, 3):
            src = chm[c]
            dst = os.path.join(fov_dir, os.path.basename(src))
            if not os.path.exists(dst):
                shutil.copy2(src, dst)

        # 2) 读取（ch00=敏化发射 S / ch01=受体 A / ch03=供体 D）
        A = tifffile.imread(chm[1]).astype(np.float64)
        D = tifffile.imread(chm[3]).astype(np.float64)
        S = tifffile.imread(chm[0]).astype(np.float64)

        A_s = gaussian_filter(A, args.sigma)
        D_s = gaussian_filter(D, args.sigma)
        S_s = gaussian_filter(S, args.sigma)

        bgA = float(np.percentile(A, 5.0))
        bgD = float(np.percentile(D, 5.0))

        # 3) 显示强度：分位数拉伸（+可选 gamma）—— Acceptor/Donor 伪彩图与 Ratio 图共用
        a_disp = view_norm(A_s - bgA) ** args.disp_gamma
        d_disp = view_norm(D_s - bgD) ** args.disp_gamma

        # 4) Ratio 图 = 基准通道灰度（明度/质感） × FRET 色阶颜色（色相） × brightness
        #    与基准通道伪彩图逐像素同质感（轮廓/纹理/亮度层次完全一致），仅颜色不同。
        base_disp = d_disp if args.base_channel == "donor" else a_disp
        w = base_disp ** args.w_gamma            # w_gamma=1.0 即原样复制该通道质感
        ratio = S_s / np.maximum(D_s, 1.0)       # 原始 S/D（LAS X 口径）
        r_norm = np.clip(ratio / args.ratio_max, 0, 1)
        ratio_rgb = (FRET_CMAP(r_norm)[..., :3].astype(np.float32)
                     * w[..., None] * args.brightness)

        # 5) 三张单页 RGB TIFF —— 与 LAS X 多通道 TIF 同款（纯像素，无标注轴）
        #    注意 description 必须纯 ASCII（中文会触发 "must be 7-bit ASCII" 报错）
        stem = os.path.join(fov_dir, prefix)
        a_u8 = np.clip(as_rgb(a_disp, ACCEPTOR_RGB) * 255.0, 0, 255).astype(np.uint8)
        d_u8 = np.clip(as_rgb(d_disp, DONOR_RGB) * 255.0, 0, 255).astype(np.uint8)
        r_u8 = np.clip(ratio_rgb * 255.0, 0, 255).astype(np.uint8)

        lbl = f"{gname} {ctype}-{'Dark' if light == 'dark' else 'Light'} P{pos}".strip()
        with tifffile.TiffWriter(f"{stem}_Acceptor.tif") as tw:
            tw.write(a_u8, photometric="rgb", compression="deflate",
                     compressionargs={"level": 6},
                     description=f"Acceptor (red)  [{lbl}]")
        with tifffile.TiffWriter(f"{stem}_Donor.tif") as tw:
            tw.write(d_u8, photometric="rgb", compression="deflate",
                     compressionargs={"level": 6},
                     description=f"Donor / BFP (blue)  [{lbl}]")
        with tifffile.TiffWriter(f"{stem}_FRETRatio.tif") as tw:
            tw.write(r_u8, photometric="rgb", compression="deflate",
                     compressionargs={"level": 6},
                     description=(f"FRET Ratio S/D  scale 0.00-{args.ratio_max:.2f}  "
                                  f"baseline {args.base_channel}  [{lbl}]"))

        preview_items.append((lbl, r_u8, g, ctype, light, pos))
        n_ok += 1

    print(f"\n生成完成：成功 {n_ok}，跳过 {n_skip}")
    print(f"输出根目录：{OUT_ROOT}")
    for s in skipped:
        print("  跳过:", s)

    # 6) 色阶图例（可选）
    if args.colorbar:
        write_colorbar(os.path.join(OUT_ROOT, "_colorbar.png"), args.ratio_max)

    # 6) 预览图（仅供快速查看，非交付物）
    if args.preview and preview_items:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        order = {"E": 0, "P": 1}
        preview_items.sort(key=lambda t: (t[2], order.get(t[3], 9),
                                          0 if t[4] == "dark" else 1, t[5]))
        n = len(preview_items)
        ncol = min(8, n)
        nrow = int(np.ceil(n / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(1.6 * ncol, 1.7 * nrow),
                                 dpi=150, facecolor="white")
        axes = np.atleast_1d(axes).ravel()
        for ax, (lbl, img, *_) in zip(axes, preview_items):
            ax.imshow(img)
            ax.set_title(lbl.replace(" · ", "\n"), fontsize=6, family="Times New Roman")
            ax.set_xticks([]); ax.set_yticks([])
        for ax in axes[n:]:
            ax.axis("off")
        fig.suptitle(f"FRET Ratio (S/D)  fixed scale 0-{args.ratio_max:.2f}  "
                     f"baseline {args.base_channel}  -- preview only",
                     fontsize=9, family="Times New Roman")
        fig.tight_layout(rect=(0, 0, 1, 0.95))   # 留出 suptitle 空间，避免与第一行标题重叠
        fig.savefig(args.preview, dpi=150, facecolor="white", bbox_inches="tight")
        plt.close(fig)
        print("preview:", args.preview)

    # 7) 打包 zip（包内保留 per_fov/<...> 结构）
    if args.no_zip:
        print("已按 --no-zip 跳过 zip 打包")
        return
    zip_path = os.path.join(BASE, args.zip_name)
    if os.path.exists(zip_path):
        os.remove(zip_path)
    # 归档根取 OUT_ROOT 的父目录 → 包内顶层即 per_fov/。
    # 注意：不能用 BASE 作归档根——OUT_ROOT 可指定到另一块盘（数据在 C:、产物在 D:），
    # 跨盘时 os.path.relpath 会抛 ValueError。
    arc_root = os.path.dirname(OUT_ROOT) or OUT_ROOT
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        for dirpath, _, fnames in os.walk(OUT_ROOT):
            for fn in fnames:
                fp = os.path.join(dirpath, fn)
                z.write(fp, os.path.relpath(fp, arc_root))
    print(f"zip written: {zip_path}  size(MB): {round(os.path.getsize(zip_path)/1e6, 1)}")


if __name__ == "__main__":
    main()
