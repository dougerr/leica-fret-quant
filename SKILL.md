---
name: leica-fret-quant
display_name: FRET 定量分析（徕卡共聚焦 · 光暗对照）
display_name_en: FRET Quantification (Leica Confocal, Light/Dark)
description: 批量定量分析徕卡共聚焦 FRET 敏化发射图像（BFP 供体 / 敏化发射双通道，LAS X 导出），支持「每样本一个子目录」与「所有视野通道平铺一夹」两种数据布局；用于光照 vs 黑暗分组的 FOV 级与逐细胞 FRET 定量、条件解析、统计对比、发表级柱状图绘制、按视野生成伪彩成像图打包，以及出版级统计图（灰柱黑点 + 空心圆散点）。核心指标 FRETN = Xia 归一化（S_fg / sqrt(D_fg × A_fg)）。触发词：FRET 定量、FRETN、光暗比、光暗对照、光控 FRET、光遗传学 FRET、敏化发射、徕卡 TIF、逐细胞 FRET、伪彩 FRET 成像。
description_zh: 对徕卡共聚焦导出的 FRET 敏化发射双通道图像做批量定量：解析光照/黑暗分组做光暗对照统计，输出 FOV 级与逐细胞三指标 CSV（含 FRETN 光暗比与 Mann-Whitney U 检验），并绘制发表级柱状图与伪彩成像图。
description_en: Batch quantification of Leica confocal sensitized-emission FRET images with light/dark comparison. Parses light/dark groups, computes FOV-level and per-cell FRET indices (FRETN normalized), exports statistical CSVs, publication-grade bar charts, and per-FOV pseudo-color imaging.
category: data-analysis
version: 1.0.0
author: 豆ger
agent_created: true
---

# Leica FRET 定量

对徕卡共聚焦导出的敏化发射 FRET 图像做批量定量：读取三通道（供体 / 敏化发射 / 受体）、
按命名规则解析实验条件（基因型 / 光照 / 类型 / 视野）、计算三个 FRET 指标、QC、
输出逐视野与分组汇总 CSV，并绘制「多构建 × 4 条件柱状图 + 右轴光暗比折线」；
另支持按视野生成出版级伪彩成像图（伪彩 TIF + 原始灰度 TIF 同夹打包）。

## 何时使用

- 提供徕卡共聚焦 TIF 数据目录，需要计算 FRET 指数、FRETN（Xia 归一化）并做「光照 vs 黑暗」统计对比。
- 需要发表级柱状图（红=光照、蓝=黑暗、深=E/浅=P）。
- 需要每个视野的成像图副本：Acceptor 红、Donor 蓝、FRET Ratio 渐变（绝对标尺 0–0.6），
  与原始灰度 TIF 打包在同一文件夹。
- 需要逐细胞分辨率的定量与统计（`analyze_cells_flat.py`）。
- 需要出版级统计图：灰柱 + 黑点（FOV 级）与空心圆散点 + mean±SD 灰线（细胞级）。

## 通道与指标（先确认，勿猜）

通道后缀严格匹配：`_ch03.tif`=供体 BFP、`_ch00.tif`=FRET 敏化发射、`_ch01.tif`=受体、`_ch02.tif`=明场（弃用）。

三指标全部输出，**推荐 FRETN 为主指标**：
- `ratio_raw` = mean(S)/mean(D)
- `ratio_gs` = S_fg / D_fg
- `FRETN` = S_fg / sqrt(D_fg × A_fg)（浓度无关，消除供体漂白伪影）

## 工作流

> **产物落位约定**：所有产物默认生成在**原始数据目录**下，便于原地查找——
> CSV/PNG 落 `--base`（逐细胞 CSV 落 `<base>/fret_cells_results/`），
> 伪彩成像落 `<base>/per_fov/` 与 `<base>/per_fov.zip`。
> 绘图脚本的 `--out` 也已默认跟随输入 CSV 所在目录。不要另建盘符 / 临时目录存放产物。

1. **确认数据布局与命名**。两种布局，通道后缀规则通用：
   - **平铺布局**：所有视野所有通道混在同一文件夹。`analyze_fret_flat.py` 自动识别两种命名：
     - `<基因型>-<D?>-<E|P>_Position<NNN>_ch<CC>.tif`（如 `S1-D-E_Position001_ch00.tif`）
     - `<共同前缀>_<基因型>-<D?>-<E|P>-<视野>_ch<CC>.tif`（如 `20250101-FRET-S1_2X-D-E-1_ch00.tif`）
   - **子目录布局**：每个样本一个子目录，命名 `<前缀>_<基因型>-<D?>-<E|P>-<视野>`，
     例如 `20250101-FRET-S1_G1-D-E-1`（解析规则见 references/metrics_and_pitfalls.md）。
   若命名不同，改对应脚本的 `parse_sample()` / `parse_condition()`。

2. **运行分析**（在隔离 venv，依赖 numpy/pandas/scipy/tifffile/imagecodecs/matplotlib）：
   ```bash
   # 平铺布局
   python scripts/analyze_fret_flat.py --base "数据目录" [--out-dir 输出目录] [--qc 20]
   # 子目录布局
   python scripts/analyze_fret.py --base "数据目录" [--out-dir 输出目录] [--qc 20]
   ```
   产出 `fret_fov.csv`（逐视野，含全部样本与三指标）、`fret_condition_summary.csv`、
   `fret_dark_vs_light.csv`（dark vs light fold change + Mann-Whitney U p 值）。
   两脚本指标口径一致，CSV 均可直接对接绘图脚本。

3. **绘图**：
   ```bash
   python scripts/plot_fret_5group.py --csv fret_fov.csv [--builds S1,S2,S3,Ctrl] [--metric FRETN] [--no-qc]
   ```
   默认按 CSV 自动识别构建并做 QC（D_fg≥20）；`--no-qc` 保留全部样本。

4. **每视野伪彩成像图打包**（需要成像图 / 按视野交付时）：
   ```bash
   python scripts/per_fov_imaging.py --base "平铺数据目录" [--rename S3=S1] [--no-zip]
   ```
   每视野生成 `per_fov/<prefix>/` 文件夹，含 6 个 TIF：
   3 个原始灰度（`_ch00/_ch01/_ch03` 复制）+ 3 个伪彩单页 RGB TIF
   （`_Acceptor.tif` 红 / `_Donor.tif` 蓝 / `_FRETRatio.tif` 渐变，绝对标尺 0–0.6），
   默认打包 `per_fov.zip`（`--no-zip` 跳过）。
   兼容两种文件名（`<前缀>_<基因型>-<D?>-<E|P>-<视野号>` 与 `_Position<NNN>`），自动识别。
   `--rename OLD=NEW` 只改 TIFF 内标注的显示名（如 S3→S1），不改文件名与数据。
   **定稿算法：Ratio 图 = 基准通道灰度（明度/纹理） × FRET 色阶（色相） × brightness，
   与基准通道伪彩图逐像素共用同一张灰度图——质感形态完全一致，仅颜色不同。**
   基准通道默认 donor（`--base-channel` 可切 acceptor）、`--ratio-max 0.6`、`--sigma 1.0`、
   `--brightness 0.65`。★ 绝不要用硬掩膜或阈值化权重当明度（会把细胞内部压平成色块，v1–v4 均因此失真）。
   算法口径与配色锚点见 references/metrics_and_pitfalls.md「每视野伪彩成像图」一节
   ——**算法已多轮验证定稿，改任何参数前先留回滚点**。

5. **解读**：看 `fret_dark_vs_light.csv` 的 FRETN 行。fold≈1 且 p>0.05 = 无光控响应；
   显著偏离 1 的构建才是有光控响应（结合阴性对照——序列阴性对照与不响应光控的对照构建——
   判断光开关本身是否工作）。

6. **逐细胞定量**（平铺布局，需要单细胞分辨率统计 / 出版级散点图时）：
   ```bash
   python scripts/analyze_cells_flat.py --base "数据目录" [--out-dir 输出目录] [--min-area 300 --max-area 20000]
   ```
   供体 Otsu→填孔→连通域→面积过滤（300–20000 px），每连通域=一个细胞，
   D/S/A 掩膜内均值扣各自背景。产出 `fret_per_cell.csv`（逐细胞）、
   `fret_per_fov.csv`（逐 FOV 中位数）、`fret_cells_dark_vs_light.csv`
   （细胞级 + FOV 级 Mann-Whitney U 对比）。注意细胞级统计存在伪重复，FOV 级列为对照。

7. **出版级统计图**（规范见 references/figure_style_spec.md）：
   ```bash
   # 左 FOV 灰柱+黑点 / 右 细胞空心圆散点（红=光照、蓝=黑暗），mean±SD 灰线
   python scripts/plot_cells_style.py --per-cell fret_per_cell.csv --per-fov fret_per_fov.csv
   ```

## 依赖环境

- Python 需装：`numpy pandas scipy tifffile imagecodecs matplotlib`。
- 徕卡 8-bit LZW TIF 必须有 `imagecodecs`，否则读取报错。
- `per_fov_imaging.py` 仅用 `matplotlib.colors`（colormap 查表），无需 GUI 后端。

## 资源

- `scripts/analyze_fret.py` — 子目录布局批量定量主脚本（Otsu 掩膜 + 三指标 + 条件解析 + 三张 CSV）。
- `scripts/analyze_fret_flat.py` — 平铺布局批量定量（与 analyze_fret.py 同口径，
  自动识别 `_Position<NNN>` 与「基因型在前」两种平铺命名）。
- `scripts/plot_fret_5group.py` — 柱状图 + 右轴光暗比折线绘图脚本。
- `scripts/per_fov_imaging.py` — 每视野伪彩成像图（Acceptor 红 / Donor 蓝 / Ratio 渐变 0–0.6）
  + 原始灰度 TIF 同夹打包 + zip；`--colorbar` 可在 per_fov 根目录附带 `_colorbar.png` 色阶图例。
- `scripts/plot_ratio_colorbar.py` — FRET Ratio 色阶图例（LUT colorbar，0–0.6），
  与伪彩图共用同一 colormap。`--orientation horizontal`（横版带刻度，默认）/
  `vertical`（竖版版式风格：裸色条 1:6，无边框无刻度无文字，`--ticks` 可加数字，
  `--bg black` 用于直接拼黑底伪彩图）。竖版实测宽高比 14:84，上=橙红(高) 下=青(低)。
- `scripts/plot_images_3panel_grid.py` — 三行组图（受体红 / 供体蓝 / Ratio 伪彩；
  `--bg white` 白底页面+黑面板+细白缝，`--bg black` 可选）。组标题带下划线、行标签居左
  （`--no-labels` 全部去掉，组信息由文件名承担）；Ratio 行右端竖色条 ≈ 面板高 70%/宽 12%，
  右侧标 0.6/0.0（黑字白底、无外框）；受体行每组末面板白标尺。
  自动按 `fret_per_fov.csv` 的 n_cells 挑高质量视野；裁剪用 donor 掩膜最密窗口且
  三行共用同一裁剪框；Ratio 明度与 Donor 行共用同一灰度。
  `--split` 每组单独出一张（3 行 × n_fov 列 + 一个色条），否则多组拼成大图（`--per-block` 分块堆叠）。
  `--split --per-build --layout light-E,dark-E,light-P,dark-P` 每构建单独一张，
  块顺序完全按 --layout（逗号分隔 <light>-<type>），每块 n_fov 个视野，
  文件名 `fret_grid_<构建>_<LightE_DarkE_LightP_DarkP>`。
  `--scale-bar all|first-group|none` 控制白色标尺。
  ★ 注意 PIL paste 不缩放：竖色条数组行数必须等于目标高度 cb_h。
- `scripts/analyze_cells_flat.py` — 逐细胞定量（供体 Otsu 分割→面积过滤→逐细胞三指标
  + 细胞级 / FOV 级 dark vs light 对比）。
- `scripts/plot_cells_style.py` — 出版级统计图（FOV 灰柱+黑点 / 细胞空心圆散点）。
- `references/metrics_and_pitfalls.md` — 指标数学定义、通道映射、两种布局命名约定、QC 规则、
  伪彩成像算法口径与配色锚点、历史踩坑清单；需要深入理解口径或排查异常时读取。
- `references/figure_style_spec.md` — 出版级统计图规范（配色常量、柱状图 / 散点图元素、
  坐标轴与字体、Y 轴自适应规则）；绘制出版级图前读取。
