# -*- coding: utf-8 -*-
"""
三行组图：受体(红) / 供体(蓝) / FRET Ratio(伪彩)。

版式（逐像素实测标定的版式参数）：
  - 黑底，面板方形、小间隙，组与组之间大间隙
  - 每组上方组标题（带下划线，横跨该组列）
  - 行标签水平文字放在最左
  - Ratio 行右端竖色条（裸条 1:6，上橙下青，无边框无刻度）
  - 受体行每组最后一个面板右下角白色标尺
  - 面板为紧裁剪（donor 掩膜最密窗口），同一视野三行共用同一裁剪框

Ratio 明度算法与 per_fov_imaging.py 定稿口径完全一致：
  Ratio = FRET_LUT(clip(S/D ÷ ratio_max)) × 供体灰度 × brightness
  ★ 明度与 Donor 行共用同一张 view_norm 灰度 → 三行形态逐像素一致。

视野挑选：优先 fret_cells_results/fret_per_fov.csv 的 n_cells（分割出细胞多=质量高），
该文件不存在时现场用 donor Otsu 分割计数。

用法：
  python plot_images_3panel_grid.py --base "数据目录" \
      --builds S1,S2,S3,S4 --type E --n-fov 3 --per-block 4 \
      [--crop 420] [--panel 220] [--brightness 0.65] [--ratio-max 0.6] \
      [--out 组图.png]
"""
import os
import re
import argparse
import numpy as np
import tifffile
from scipy import ndimage
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image, ImageDraw, ImageFont

# 与 per_fov_imaging.py 完全一致的锚点（从论文图采样的锚点）
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
ACCEPTOR_RGB = (1.00, 0.18, 0.18)
DONOR_RGB = (0.28, 0.60, 1.00)

FNAME_RE = re.compile(r"^(?P<prefix>.+)_ch(?P<ch>\d{2})\.tif[f]?$", re.IGNORECASE)


def parse_sample(prefix):
    """兼容 <基因型>-<D?>-<E|P>_Position<NNN> 与 <前缀>_<基因型>-<D?>-<E|P>-<视野>。"""
    m = re.search(r"_Position(\d+)$", prefix, re.IGNORECASE)
    if m:
        core, pos = prefix[:m.start()], m.group(1)
    else:
        core = prefix.split("_")[-1]
        bits = core.split("-")
        pos = bits[-1] if bits[-1].isdigit() else ""
    parts = core.split("-")
    g = parts[0]
    light = "dark" if "D" in parts[1:] else "light"
    ctype = "E" if "E" in parts[1:] else ("P" if "P" in parts[1:] else "")
    return g, light, ctype, pos


def otsu_threshold(img):
    img = img.astype(np.float64)
    mx = float(img.max())
    if mx <= 0:
        return 0.0
    hist, edges = np.histogram(img, bins=256, range=(0.0, mx))
    centers = (edges[:-1] + edges[1:]) / 2.0
    total, sum_all, wb, sb = img.size, float((centers * hist).sum()), 0.0, 0.0
    bv, bt = -1.0, 0.0
    for i in range(256):
        wb += hist[i]
        if wb == 0:
            continue
        wf = total - wb
        if wf <= 0:
            break
        sb += centers[i] * hist[i]
        v = wb * wf * ((sb / wb) - ((sum_all - sb) / wf)) ** 2
        if v > bv:
            bv, bt = v, centers[i]
    return bt


def count_cells(donor):
    mask = donor > otsu_threshold(donor)
    mask = ndimage.binary_fill_holes(mask)
    lab, n = ndimage.label(mask)
    if n == 0:
        return 0
    areas = ndimage.sum(np.ones_like(lab), lab, np.arange(1, n + 1))
    return int(((areas >= 300) & (areas <= 20000)).sum())


def best_crop_window(mask, size, stride=8):
    """donor 掩膜最密窗口 + 窗口内质心精调。返回 (y, x) 左上角。"""
    H, W = mask.shape
    if H <= size and W <= size:
        return 0, 0
    size = min(size, H, W)
    m = mask.astype(np.float64)
    ii = np.zeros((H + 1, W + 1))
    ii[1:, 1:] = m.cumsum(0).cumsum(1)
    best, best_xy = -1.0, (0, 0)
    for y in range(0, H - size + 1, stride):
        for x in range(0, W - size + 1, stride):
            s = ii[y + size, x + size] - ii[y, x + size] - ii[y + size, x] + ii[y, x]
            if s > best:
                best, best_xy = s, (y, x)
    y, x = best_xy
    sub = mask[y:y + size, x:x + size]
    if sub.any():
        ys, xs = np.nonzero(sub)
        ny = int(round(y + ys.mean() - size / 2))
        nx = int(round(x + xs.mean() - size / 2))
        y, x = min(max(0, ny), H - size), min(max(0, nx), W - size)
    return y, x


def view_norm(x, p_lo=2.0, p_hi=99.5):
    lo, hi = float(np.percentile(x, p_lo)), float(np.percentile(x, p_hi))
    if hi <= lo + 1e-6:
        return np.zeros_like(x, dtype=np.float32)
    return np.clip((x - lo) / (hi - lo), 0, 1).astype(np.float32)


def to_u8(rgb):
    return np.clip(rgb * 255.0, 0, 255).astype(np.uint8)


def load_font(size):
    for p in (r"C:\Windows\Fonts\times.ttf", r"C:\Windows\Fonts\timesbd.ttf",
              r"C:\Windows\Fonts\arial.ttf"):
        if os.path.exists(p):
            return ImageFont.truetype(p, size)
    return ImageFont.load_default()


def um_per_px(path, fallback=0.1804):
    try:
        with tifffile.TiffFile(path) as tf:
            tags = tf.pages[0].tags
            if "XResolution" in tags and "ResolutionUnit" in tags:
                num, den = tags["XResolution"].value
                unit = tags["ResolutionUnit"].value
                cm = unit in (3, "centimeter")
                if num and den:
                    return (1.0 / (num / den)) * (10000.0 if cm else 25400.0 / 1000.0)
    except Exception:
        pass
    return fallback


def main():
    ap = argparse.ArgumentParser(description="三行组图")
    ap.add_argument("--base", required=True)
    ap.add_argument("--builds", required=True, help="逗号分隔构建名，如 S1,S2,S3,S4")
    ap.add_argument("--type", default="E", choices=["E", "P"], help="展示的类型（单条件时用）")
    ap.add_argument("--types", default="",
                    help="逗号分隔的并行条件，如 E,P：每张图内每种条件一块，"
                         "排在该构建块之后、S2 同条件块之前；缺省只用 --type")
    ap.add_argument("--conditions", default="light,dark",
                    help="逗号分隔，取 light / dark 的子集（做两两对照时可只取一种）")
    ap.add_argument("--n-fov", type=int, default=3, help="每组视野数")
    ap.add_argument("--per-block", type=int, default=4, help="每个块（横向）放的组数")
    ap.add_argument("--crop", type=int, default=420, help="裁剪窗口边长 px")
    ap.add_argument("--panel", type=int, default=220, help="面板显示边长 px")
    ap.add_argument("--brightness", type=float, default=0.65)
    ap.add_argument("--ratio-max", type=float, default=0.6)
    ap.add_argument("--sigma", type=float, default=1.0)
    ap.add_argument("--scale-bar-um", type=float, default=20.0)
    ap.add_argument("--scale-bar", default="all", choices=["all", "first-group", "none"],
                    help="白色标尺画在哪：all=每个块末面板都画；first-group=只画第一块"
                         "（两两比对时右侧参照组不画）；none=全部不画")
    ap.add_argument("--row-labels", default="Acceptor,Donor / BFP,FRET Ratio")
    ap.add_argument("--bg", default="white", choices=["white", "black"],
                    help="画布背景：论文页面为白底（面板自身仍是黑的）")
    ap.add_argument("--no-labels", action="store_true",
                    help="不画行标签与组标题（组信息由文件名承担），色条仅保留 0.6/0.0 刻度")
    ap.add_argument("--out", default="", help="输出 PNG 路径（默认数据目录下 fret_images_3panel_grid.png）")
    ap.add_argument("--split", action="store_true",
                    help="每组单独出一张图（3 行 × n_fov 列 + 一个色条），文件名加组名后缀")
    ap.add_argument("--per-build", action="store_true", dest="per_build",
                    help="每个构建单独一张图，块顺序由 --layout 指定（不再与参照构建成对）")
    ap.add_argument("--layout", default="light-E,dark-E,light-P,dark-P",
                    help="--per-build 时单张图内的块顺序，逗号分隔的 <light>-<type>，"
                         "light/dark 可用 l/d 简写，如 light-E,dark-E,light-P,dark-P")
    ap.add_argument("--pair-with", default="", dest="pair_with",
                    help="参照构建名（如 S2）：split 模式下，每个非该构建的组与其同光照条件组"
                         "拼成一张两块比对图（左右两组比对），文件名 fret_compare_*")
    ap.add_argument("--save-tif", action="store_true", help="同时输出 TIF")
    args = ap.parse_args()

    BASE = os.path.abspath(args.base)
    OUT = os.path.abspath(args.out) if args.out else os.path.join(BASE, "fret_images_3panel_grid.png")
    os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
    builds = [b.strip() for b in args.builds.split(",") if b.strip()]
    row_labels = [s.strip() for s in args.row_labels.split(",")]

    # ---- 收集样本并挑视野 ----
    ch_files = {}
    for fn in sorted(os.listdir(BASE)):
        m = FNAME_RE.match(fn)
        if m and fn.lower().endswith((".tif", ".tiff")):
            ch_files.setdefault(m["prefix"], {})[int(m["ch"])] = os.path.join(BASE, fn)

    per_fov_csv = os.path.join(BASE, "fret_cells_results", "fret_per_fov.csv")
    ncell_map = {}
    if os.path.exists(per_fov_csv):
        import csv
        with open(per_fov_csv, encoding="utf-8-sig") as f:
            for r in csv.DictReader(f):
                try:
                    ncell_map[r["sample"]] = int(r["n_cells"])
                except Exception:
                    pass

    # 块布局：--per-build 用 --layout 展开成 [(light, type), ...]
    layout_items = []
    if args.per_build:
        for tok in args.layout.split(","):
            tok = tok.strip()
            if "-" not in tok:
                continue
            lt, tp = tok.rsplit("-", 1)
            layout_items.append(("dark" if lt.strip().lower().startswith("d") else "light",
                                 tp.strip().upper()))
    types_list = [t.strip().upper() for t in args.types.split(",") if t.strip()] or [args.type]
    for _, tp in layout_items:
        if tp not in types_list:
            types_list.append(tp)
    multi_t = len(types_list) > 1

    groups = []   # (title, [prefix...])
    group_map = {}  # (build, light, ctype) -> (title, [prefix...])；含 --pair-with 的参照组
    needed = list(builds)
    if args.pair_with and args.pair_with not in needed:
        needed.append(args.pair_with)
    print("=== 视野挑选（按可分割细胞数，括号内为备选池） ===")
    for g in needed:
        for light in ("light", "dark"):
            for t in types_list:
                pool = []
                for prefix, chm in ch_files.items():
                    if not all(c in chm for c in (0, 1, 3)):
                        continue
                    gg, lt, ct, pos = parse_sample(prefix)
                    if gg == g and lt == light and ct == t:
                        n = ncell_map.get(prefix)
                        if n is None:
                            D0 = tifffile.imread(chm[3]).astype(np.float64)
                            n = count_cells(D0)
                        pool.append((n, pos, prefix))
                pool.sort(key=lambda t2: (-t2[0], t2[1]))
                picked = [p for _, _, p in pool[:args.n_fov]]
                picked.sort(key=lambda p: parse_sample(p)[3])
                cn = "Light" if light == "light" else "Dark"
                title = f"{g}-{t}  {cn}" if multi_t else f"{g}  {cn}"
                group_map[(g, light, t)] = (title, picked)
                if g in builds:
                    groups.append((title, picked))
                detail = ", ".join(f"{parse_sample(p)[3]}({n})" for n, _, p in pool)
                print(f"  {g}-{t}-{light:5s}: 选 {['P' + parse_sample(p)[3] for p in picked]}  池: {detail}")

    # ---- 渲染每个视野的三行面板 ----
    P = args.panel
    print("\n=== 渲染面板 ===")
    panels = {}   # prefix -> (acc_u8, don_u8, rat_u8, px_um)
    for _, prefixes in groups:
        for prefix in prefixes:
            chm = ch_files[prefix]
            A = tifffile.imread(chm[1]).astype(np.float64)
            D = tifffile.imread(chm[3]).astype(np.float64)
            S = tifffile.imread(chm[0]).astype(np.float64)
            A_s = ndimage.gaussian_filter(A, args.sigma)
            D_s = ndimage.gaussian_filter(D, args.sigma)
            S_s = ndimage.gaussian_filter(S, args.sigma)
            bgA, bgD = float(np.percentile(A, 5.0)), float(np.percentile(D, 5.0))
            mask = D_s - bgD > otsu_threshold(D_s - bgD)
            mask = ndimage.binary_fill_holes(mask)
            y, x = best_crop_window(mask, args.crop)
            sl = (slice(y, y + args.crop), slice(x, x + args.crop))

            a_disp = view_norm(A_s[sl] - bgA)
            d_disp = view_norm(D_s[sl] - bgD)
            ratio = S_s[sl] / np.maximum(D_s[sl], 1.0)
            r_norm = np.clip(ratio / args.ratio_max, 0, 1)
            rgb_r = FRET_CMAP(r_norm)[..., :3].astype(np.float32) * d_disp[..., None] * args.brightness

            acc = to_u8(a_disp[..., None] * np.array(ACCEPTOR_RGB, dtype=np.float32))
            don = to_u8(d_disp[..., None] * np.array(DONOR_RGB, dtype=np.float32))
            rat = to_u8(rgb_r)
            panels[prefix] = (acc, don, rat, um_per_px(chm[3]))
            print(f"  {prefix}: crop@({y},{x})  细胞px={int(mask.sum())}")

    # ---- 版式常量（论文同款比例：页面白底、面板黑、色条 ≈ 面板高 70% / 宽 12%） ----
    gs = max(3, P // 55)          # 组内面板间隙（论文为细白缝）
    gg_ = max(18, P // 8)         # 组间隙
    if args.no_labels:
        label_w = int(P * 0.06)   # 只留窄边距
        title_h = int(P * 0.04)
    else:
        label_w = int(P * 0.95)   # 左侧行标签区
        title_h = int(P * 0.24)   # 组标题区
    row_gap = max(3, P // 55)
    pad_r = int(P * 0.30)         # 色条 + 0.6/0.0 文字区
    bot = int(P * 0.10)
    margin_x = int(P * 0.06)
    margin_y = int(P * 0.05)
    cb_w = max(8, int(P * 0.12))          # 论文实测 13/107 ≈ 0.12
    cb_h = int(P * 0.70)                  # 论文实测 75/107 ≈ 0.70
    fg = (0, 0, 0) if args.bg == "white" else (255, 255, 255)
    f_title = load_font(max(14, int(P * 0.135)))
    f_label = load_font(max(13, int(P * 0.125)))
    f_cb = load_font(max(11, int(P * 0.15)))

    def block_size(blk):
        gw = args.n_fov * P + (args.n_fov - 1) * gs
        w = (margin_x + label_w + len(blk) * gw + (len(blk) - 1) * gg_
             + pad_r + cb_w + int(P * 0.10))
        h = margin_y + title_h + 3 * P + 2 * row_gap + bot
        return w, h

    def render_block(blk, W, H):
        """blk = [(组标题, [prefix...]), ...]；返回论文样式图（白底、黑面板、细白缝）。"""
        img = Image.new("RGB", (W, H), (255, 255, 255) if args.bg == "white" else (0, 0, 0))
        dr = ImageDraw.Draw(img)
        x0, y_off = margin_x, margin_y
        y_rows = [y_off + title_h + r * (P + row_gap) for r in range(3)]
        x = x0 + label_w
        gw = args.n_fov * P + (args.n_fov - 1) * gs
        for gi, (title, prefixes) in enumerate(blk):
            if not args.no_labels:
                tw = dr.textlength(title, font=f_title)
                dr.text((x + gw / 2 - tw / 2, y_off + title_h * 0.30), title,
                        font=f_title, fill=fg)
                uy = y_off + title_h * 0.86
                dr.line([(x + gw / 2 - tw / 2, uy), (x + gw / 2 + tw / 2, uy)],
                        fill=fg, width=2)
            for fi, prefix in enumerate(prefixes):
                acc, don, rat, px_um = panels[prefix]
                for ri, arr in enumerate((acc, don, rat)):
                    tile = Image.fromarray(arr, "RGB").resize((P, P), Image.LANCZOS)
                    img.paste(tile, (x + fi * (P + gs), y_rows[ri]))
                # 白色标尺：受体行每组最后一个面板右下角（可按 --scale-bar 限制）
                draw_sb = (args.scale_bar == "all") or \
                          (args.scale_bar == "first-group" and gi == 0)
                if fi == len(prefixes) - 1 and draw_sb:
                    L = min(int(args.scale_bar_um / px_um * P / args.crop), P - 20)
                    bx = x + fi * (P + gs) + P - L - int(P * 0.09)
                    by = y_rows[0] + P - int(P * 0.09)
                    dr.rectangle([bx, by, bx + L, by + max(3, P // 40)], fill="white")
            x += gw + gg_
        if not args.no_labels:
            for ri, lab in enumerate(row_labels):
                dr.text((x0, y_rows[ri] + P / 2), lab, font=f_label, fill=fg, anchor="lm")
        # 竖色条：Ratio 行右端（论文样式：无外框，右侧标 0.6 / 0.0）
        cb_x = (x0 + label_w + len(blk) * gw + (len(blk) - 1) * gg_ + int(P * 0.06))
        cb_y = y_rows[2] + (P - cb_h) // 2
        # ★ 数组行数必须等于 cb_h：paste 不会缩放，直接用 512 行会超出画布、裁掉青色端
        grad = np.linspace(1, 0, cb_h)[:, None] * np.ones((1, cb_w))
        img.paste(Image.fromarray(to_u8(FRET_CMAP(grad)[..., :3]), "RGB"), (cb_x, cb_y))
        tx = cb_x + cb_w + int(P * 0.05)
        dr.text((tx, cb_y), f"{args.ratio_max:.1f}", font=f_cb, fill=fg, anchor="lm")
        dr.text((tx, cb_y + cb_h), "0.0", font=f_cb, fill=fg, anchor="lm")
        return img

    stem, ext = os.path.splitext(OUT)

    def save_fig(blk, fname):
        """blk = [(组标题, [prefix...]), ...]，渲染单块并落盘（PNG + 可选 TIF）。"""
        W, H = block_size(blk)
        im = render_block(blk, W, H)
        p = os.path.join(os.path.dirname(OUT), fname + ext)
        im.save(p)
        print(f"  Saved: {p}  {im.size}")
        if args.save_tif:
            tp = os.path.splitext(p)[0] + ".tif"
            im.save(tp)
            print(f"  Saved: {tp}")

    if args.split and args.per_build:
        # 每构建一张：块顺序完全按 --layout（默认 光照-E | 黑暗-E | 光照-P | 黑暗-P），每块 n_fov 个视野
        print(f"\n=== 单构建出图（布局 {args.layout}） ===")
        tag = "_".join(("Light" if lt == "light" else "Dark") + tp
                       for lt, tp in layout_items)
        for g in builds:
            blk = [group_map[(g, lt, tp)] for lt, tp in layout_items]
            save_fig(blk, f"fret_grid_{g}_{tag}")
    elif args.split and args.pair_with:
        # 两两比对：每个非参照构建的条件块 ↔ 参照构建同光照同条件块。
        # 多条件（--types E,P）时一张图 4 块：[本组-E | 本组-P | S2-E | S2-P]
        pw = args.pair_with
        ttag = "-".join(types_list)
        print(f"\n=== 与 {pw} 两两比对出图（{ttag}） ===")
        for g in builds:
            if g == pw:
                continue
            for light in ("light", "dark"):
                cn = "Light" if light == "light" else "Dark"
                blk = [group_map[(g, light, t)] for t in types_list] + \
                      [group_map[(pw, light, t)] for t in types_list]
                if multi_t:
                    fname = f"fret_compare_{g}_{cn}_{ttag}_vs_{pw}"
                else:
                    t0 = types_list[0]
                    fname = f"fret_compare_{g}_{t0}_{cn}_vs_{pw}_{t0}_{cn}"
                save_fig(blk, fname)
    elif args.split:
        # 每组单独一张：3 行 × n_fov 列 + 一个色条（论文左/右半边那种单块）
        print("\n=== 分组出图 ===")
        for title, prefixes in groups:
            slug = re.sub(r"[^0-9A-Za-z]+", "_", title).strip("_")
            save_fig([(title, prefixes)], f"{os.path.basename(stem)}_{slug}")
    else:
        blocks = [groups[i:i + args.per_block] for i in range(0, len(groups), args.per_block)]
        sizes = [block_size(b) for b in blocks]
        W = max(w for w, _ in sizes)
        H = sum(h for _, h in sizes) + 60 * (len(blocks) - 1)
        canvas = Image.new("RGB", (W, H), (0, 0, 0))
        yo = 0
        for b, (_, h) in zip(blocks, sizes):
            canvas.paste(render_block(b, W, h), (0, yo))
            yo += h + 60
        canvas.save(OUT)
        print("\nSaved:", OUT, canvas.size)
        if args.save_tif:
            tif_out = os.path.splitext(OUT)[0] + ".tif"
            canvas.save(tif_out)
            print("Saved:", tif_out)


if __name__ == "__main__":
    main()
