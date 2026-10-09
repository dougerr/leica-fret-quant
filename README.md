# leica-fret-quant

**徕卡共聚焦 FRET 敏化发射图像的批量定量、光暗对照分析与发表级出图** — 一个遵循 [Agent Skills](https://agentskills.io) 开放标准的技能包。

*Batch quantification and publication-grade figure generation for Leica confocal sensitized-emission FRET images,
including light/dark comparison — an [Agent Skills](https://agentskills.io) compatible skill package.*

---

## 这是什么 / What is this

一个"教 AI 怎么分析 FRET 图像"的技能包：把一套可复现的 FRET 定量流程（通道解析 → 掩膜分割 → 三指标计算 → QC → **光照 vs 黑暗统计对比** → 出图）写成 `SKILL.md` 指令 + 可直接执行的 Python 脚本。

装上之后，你只要把徕卡导出的 TIF 目录丢给 AI 并说"跑一下 FRET 定量"，它就会按这套流程做完，输出 CSV 与发表级图。

An Agent Skills package that turns a reproducible FRET quantification pipeline (channel parsing → masking → three indices → QC → statistics → figures) into a `SKILL.md` instruction set plus runnable Python scripts.

## 功能 / Features

| 脚本 / Script | 作用 / Purpose |
|---|---|
| `analyze_fret_flat.py` | 平铺布局批量定量（所有视野所有通道混在一个文件夹） |
| `analyze_fret.py` | 子目录布局批量定量（每样本一个子目录） |
| `analyze_cells_flat.py` | 逐细胞定量（供体 Otsu 分割 → 面积过滤 → 单细胞统计） |
| `plot_fret_5group.py` | 多构建 × 4 条件柱状图 + 右轴光暗比折线 |
| `plot_cells_style.py` | 出版级统计图（FOV 灰柱+黑点 / 细胞空心圆散点） |
| `per_fov_imaging.py` | 每视野伪彩成像图（Acceptor 红 / Donor 蓝 / Ratio 渐变）+ 原始灰度 TIF 打包 |
| `plot_ratio_colorbar.py` | FRET Ratio 色阶图例（0–0.6） |
| `plot_images_3panel_grid.py` | 三行组图（受体 / 供体 / Ratio），用于图版排版 |

## 安装 / Installation

**方式一：WorkBuddy 技能市场**
在 WorkBuddy 的「专家 · 技能 · 连接器 → 技能」里搜索本技能名，点安装。

**方式二：任何支持 Agent Skills 标准的客户端**
（Claude Code、Cursor、Codex CLI、Gemini CLI、GitHub Copilot 等）

```bash
git clone https://github.com/dougerr/leica-fret-quant.git ~/.workbuddy/skills/leica-fret-quant
```

**方式三：手动**
把本仓库整个目录复制到你的 skills 目录下，目录名保持 `leica-fret-quant`。

## 依赖 / Requirements

```bash
pip install -r requirements.txt
```

- Python 3.11+
- `imagecodecs` 是**硬依赖**：徕卡导出的 8-bit LZW 压缩 TIF 没有它读不出来。
- 无 GUI 依赖，服务器/命令行环境可直接跑。

## 数据要求 / Input requirements

- **通道后缀必须严格匹配**（LAS X 0-based 导出）：

  | 后缀 | 通道 |
  |---|---|
  | `_ch00.tif` | FRET 敏化发射（405 nm 激发 / 受体发射窗口） |
  | `_ch01.tif` | 受体直接激发（受体表达量） |
  | `_ch02.tif` | 明场 / 透射光（**弃用**，别当 FRET） |
  | `_ch03.tif` | 供体 BFP（405 nm 激发 / 450 nm 发射） |

- 支持两种文件命名布局，详见 [`references/metrics_and_pitfalls.md`](references/metrics_and_pitfalls.md)。
  命名不同时改脚本里的 `parse_sample()` / `parse_condition()` 即可。
- 8-bit LZW TIF。

## 快速开始 / Quick start

```bash
# 1) 平铺布局批量定量 → fret_fov.csv / fret_condition_summary.csv / fret_dark_vs_light.csv
python scripts/analyze_fret_flat.py --base "<数据目录>"

# 2) 柱状图 + 光暗比折线
python scripts/plot_fret_5group.py --csv "<数据目录>/fret_fov.csv" --metric FRETN

# 3) 逐细胞定量（需要单细胞分辨率时）
python scripts/analyze_cells_flat.py --base "<数据目录>"

# 4) 每视野伪彩成像图打包
python scripts/per_fov_imaging.py --base "<数据目录>"
```

> 所有产物默认落在**原始数据目录**下，便于原地查找。

## 核心指标 / Core indices

设 `D`=供体、`S`=敏化发射、`A`=受体，`_fg` 表示掩膜内均值扣除背景后的净强度：

| 指标 | 公式 | 适用 |
|---|---|---|
| `ratio_raw` | mean(S)/mean(D) | 整视野、不扣背景，作对照 |
| `ratio_gs` | S_fg / D_fg | 扣背景，但只对供体归一化，对供体漂白敏感 |
| **`FRETN`** ★ | S_fg / √(D_fg × A_fg) | **推荐主指标**：对供体+受体做几何平均，浓度无关 |

**为什么要用 FRETN**：`ratio_raw` / `ratio_gs` 只对供体归一化，光照组一旦发生供体光漂白，
光暗比会出现假性升高（实测可造成约 2 倍的假阳性）。FRETN 用 √(D×A) 归一化可消除该伪影。

## 输出 / Outputs

| 文件 | 内容 |
|---|---|
| `fret_fov.csv` | 逐视野，含全部样本与三指标 |
| `fret_condition_summary.csv` | 分组汇总 |
| `fret_dark_vs_light.csv` | 光照 vs 黑暗 fold change + Mann-Whitney U p 值 |
| `fret_cells_results/fret_per_cell.csv` | 逐细胞（逐细胞脚本） |
| `per_fov/` + `per_fov.zip` | 每视野伪彩 TIF + 原始灰度 TIF 打包 |

## 目录结构 / Layout

```
leica-fret-quant/
├── SKILL.md                        # 技能定义与工作流（AI 读这个）
├── scripts/                        # 可执行脚本
├── references/
│   ├── metrics_and_pitfalls.md     # 指标口径、通道映射、命名约定、QC、踩坑清单
│   └── figure_style_spec.md        # 出版级统计图规范
├── requirements.txt
├── CITATION.cff
└── LICENSE
```

## 方法学说明 / Methodology notes

- 指标定义、通道映射、命名解析规则、QC 阈值、伪彩成像算法口径与历史踩坑，
  全部记录在 [`references/metrics_and_pitfalls.md`](references/metrics_and_pitfalls.md)。
- 出版级统计图的配色常量与坐标轴规范（从目标图版逐像素采样）见
  [`references/figure_style_spec.md`](references/figure_style_spec.md)。
- **PR 欢迎**：如果你的命名约定或显微系统与默认不同，改 `parse_*()` 与标尺参数后提 PR。

## 引用 / Citation

见 [`CITATION.cff`](CITATION.cff)。如果在论文中使用本工具，请引用本仓库。

## 更新日志 / Changelog

| 版本 | 日期 | 变化 |
|---|---|---|
| `1.0.1` | 2026-10-09 | 专业表述规范化：全篇「论文版式统计图」统一为「出版级统计图」；样式规范文件更名为 `references/figure_style_spec.md`，统计图脚本更名为 `scripts/plot_cells_style.py`。**功能与 1.0.0 完全一致**（脚本参数、输出文件名、算法、依赖均未改动） |
| `1.0.0` | 2026-10-09 | 首个发布版本 |

## 许可 / License

[MIT](LICENSE) © 2026 豆ger

## 免责声明 / Disclaimer

- 本工具面向**科研用途**，不适用于临床诊断或任何医疗决策。
- 伪彩成像图的参数（绝对标尺 0–0.6、brightness 0.65 等）为特定成像条件标定，
  换显微系统或滤光片组后需重新标定。
- 定量结论以 CSV 为准，**不要从伪彩图读数**。
