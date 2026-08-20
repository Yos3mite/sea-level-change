# Jin et al. Figure 3 复现进度（暂停点）

更新时间：2026-08-20  
工作分支：`feat/figure03-regional-tws`  
状态：代码与运行配置已搭建并通过单元/科学测试；论文数据版尚未执行最终计算与制图。

## 1. 当前目标与约束

目标是建立一套统一的 Figure 3 区域陆地水储量（TWS）处理流程，使同一份代码既能读取论文使用的 CSR、JPL、GSFC Mascon，也能读取本地自制 Level-3 产品，并输出各洲及全球陆地区域对等效全球平均海平面（ESL）的贡献曲线和两次 El Niño 事件的端点变化。

本阶段按以下约束执行：

- 论文 Mascon 路径使用 CSR、JPL、GSFC 三中心产品，三中心先分别处理，再作算术平均。
- 使用 Natural Earth 洲界；排除格陵兰和南极洲。
- 对陆地侧实施 300 km 海岸缓冲剔除。
- 论文 Mascon 若产品已包含 GIA 校正，不再重复施加 GIA；Figure 3 区域 TWS 流程不使用 OBD。
- Xie–Yi 数据只用于补齐原产品缺测月，观测月始终优先；重建月份的使用上限为 2022-12。
- 当前任务已暂停，不合并 `main`，只保存并推送功能分支。

## 2. 已完成内容

### 2.1 设计、计划与数据登记

- `2cecc3a`：统一 Figure 3 区域 TWS 流程设计。
- `17a51d6`：实现计划。
- `47d1d85`：登记 Figure 3 数据源与来源清单。

### 2.2 数据适配与重建月份合并

- `09d9d0a`：建立统一 `MonthlyGridSeries` 数据结构和网格产品适配器。
- 支持数值时间轴及本地 Level-3 的 `YYYY-MM` 字符串时间轴。
- 支持 `mm EWH` 等效水高单位。
- 支持按研究时段切片读取，避免整份多年高分辨率网格全部进入内存。
- `c7c0e40`：按中心合并观测月与 Xie–Yi 重建月；观测值覆盖重建值，并在统一 1° 网格上分别重采样。

### 2.3 区域掩膜、积分和时间处理

- `5d360c4`：生成六大洲掩膜，排除格陵兰/南极洲，实施陆侧 300 km 海岸剔除，使用球面精确格网面积积分。
- 区域海平面贡献公式：

  `ESL_region = sum(EWH_mm * cell_area) / global_ocean_area`

- `d813a3f`：各中心独立执行完整时段月气候态去除、OLS 去趋势、严格居中 3 个月滑动平均，然后计算中心平均和中心间样本离散度。
- 两次事件的变化量采用带符号端点差：
  - 2014–2016 事件：2014-10 至 2015-12。
  - 2023–2024 事件：2023-05 至 2023-12。

### 2.4 输出与配置

- `5da0dfb`：实现可审计的 Figure 3 构建器，可生成图、区域序列、事件指标、掩膜、输入清单和运行元数据等 9 类产物。
- 已新增论文 Mascon 配置：`config/figure03_paper_mascon.json`。
- 已新增本地自制 Level-3 模板：`config/figure03_custom_l3.template.json`。
- 目标网格配置现可使用规则中心网格参数（起点、终点、间距）展开，无须在 JSON 内列出全部经纬度。

## 3. 已登记的输入数据

### 3.1 Xie–Yi 缺测月重建

- 原始压缩包：`C:/Users/Alan/Downloads/mascon solution.zip`
- 文件大小：65,872,990 bytes
- SHA-256：`a9a03d46fa41381855a52a0efa76bf07fc071bd413e04293bced19d3e5c34948`
- 数据 DOI：`10.6084/m9.figshare.25805092.v2`
- 项目内解压目录：`data/external_downloads/xie_yi_2025_mascon_gapfilled/`
- 三个 NetCDF 仅含 35 个重建缺测月份，不是完整连续 Mascon 时序。

### 3.2 Natural Earth 洲界

- 数据版本：Natural Earth Admin-0 Countries 5.1.1。
- 项目内目录：`data/external_downloads/natural_earth/admin_0_countries/ne_10m_admin_0_countries_v5.1.1/`。

### 3.3 论文 Mascon 版本口径

- CSR：RL06.3。
- JPL：RL06.3Mv04 CRI。
- GSFC：论文参考文献口径为 RL06.2；当前拟用官方页面标注的 RL06v2.0 OBP 半度产品，但输入文件仍需最终核验，见下一节。

## 4. 当前阻塞点：GSFC 输入尚未完成可信核验

论文数据版目前不能作为最终结果运行，原因不是计算代码，而是 GSFC 观测产品的版本/文件完整性尚未确认：

- `D:/temp_sealevel_data/gsfc_mascon_rl06v2.0.nc` 以及 Downloads 中同名完整 NetCDF 的内部全局属性均显示 `RL06 v1.0` / `product_version v1.0`，与文件名和网页标注的 v2.0 不一致。
- 项目内 `data/external_downloads/mascon/gsfc.glb_.200204_202603_rl06v2.0_obp-ice6gd_halfdegree.nc` 是一次中断下载留下的约 16 MB 不完整文件，**不得作为有效输入**。
- `D:/temp_sealevel_data/gsfc.glb.200204_202603_RL06v2.0_SLA-ICE6GD.h5` 可能是真正的 v2.0 SLA 产品，但尚未验证内部结构，也尚未建立 HDF5 适配或转换路径。

因此，`config/figure03_paper_mascon.json` 已通过结构契约测试，但其中 GSFC 路径当前只是待替换位置；在解决该输入前不要把论文配置产生的结果视为正式复现结果。

## 5. 测试状态

最后一次 Figure 3 测试结果：`34 passed in 2.54s`。

复现命令（在项目目录执行）：

```powershell
$env:PYTHONPATH="$PWD\.runtime\site-packages;$PWD"
$env:WINDIR='C:\Windows'
$env:MPLCONFIGDIR="$PWD\tmp\matplotlib-config"
$env:MPLBACKEND='Agg'
$figure3Tests = Get-ChildItem -LiteralPath tests/scientific -Filter 'test_figure3_*.py' | ForEach-Object FullName
& 'C:\Users\Alan\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -m pytest $figure3Tests -q
```

既有 Figure 1/2 回归测试此前为 14/14 通过；本暂停点提交前应再次执行 Figure 3 全套测试。

## 6. 本地自制 Level-3 的兼容状态

模板默认读取：

`results/target_b_l2_to_l3_multisource/target_b_custom_l3_multisource_201311_202410.nc`

适配器已支持其字符串月份、`mm EWH`、`valid_month` 和规则网格格式。该产品元数据显示当前未施加 GIA、未施加 OBD。流程为了避免隐式修改，不会额外加入 GIA 或 OBD；因此它可以用于检验统一接口和区域处理方法，但在明确自制产品的 GIA 处理方案前，不应宣称与论文三中心 Mascon 已达到完全同口径比较。

## 7. 恢复任务时的建议顺序

1. 核验 `D:/temp_sealevel_data/gsfc.glb.200204_202603_RL06v2.0_SLA-ICE6GD.h5` 的产品版本、变量、单位、时间轴和 GIA 元数据；或重新获得完整、可验证的官方 GSFC RL06v2.0/RL06.2 观测产品。
2. 计算并登记最终 GSFC 文件的大小与 SHA-256，替换论文配置中的不完整文件路径。
3. 先运行论文 Mascon 配置，逐项检查观测/重建月份优先级、区域掩膜面积、六洲之和与 Total、端点月份和三中心离散度。
4. 对 Figure 3 的 PNG/PDF 做视觉核验，并将区域曲线和论文参考变化量做定量对照。
5. 再运行本地自制 Level-3 模板，输出自制产品版本 Figure 3；结果说明中明确 GIA 口径差异。
6. 运行 Figure 1/2/3 相关回归测试后，再决定是否合并功能分支。

## 8. 暂停时的代码位置

- 主构建器：`pysrc/reference_products/build_figure3_regional_tws.py`
- Figure 3 模块：`pysrc/reference_products/figure3/`
- 论文配置：`config/figure03_paper_mascon.json`
- 自制 Level-3 模板：`config/figure03_custom_l3.template.json`
- 测试：`tests/scientific/test_figure3_*.py`
- 设计：`docs/superpowers/specs/2026-08-20-figure03-regional-tws-design.md`
- 计划：`docs/superpowers/plans/2026-08-20-figure03-regional-tws.md`

该记录的目的，是让下一次恢复任务时可以从 GSFC 数据核验直接继续，而不必重新推导 Figure 3 的数据口径和处理链。
