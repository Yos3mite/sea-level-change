# SaGEA 原版数据清单

更新时间：2026-08-20

## 1. 项目范围与数据边界

- 代码根目录：`D:\AAAA海平面变化\SaGEA-main(L2-L3)原版\SaGEA-main`
- 项目代码和结果只写入上述目录。
- 原始外部数据允许从 `D:\temp_sealevel_data` 读取，但不修改原始文件。
- 本清单区分：原始产品、项目标准化结果、诊断/替代产品。
- 不使用 `reproduce_figure1` 目录、可视化工作流或任何 L2toL3 工作流；本项目只按 `CODEX_PROJECT.md` 执行。

## 2. 原版 SaGEA 数据

| 目录 | 数量 | 用途与状态 |
|---|---:|---|
| `data/auxiliary/` | 3 | 球谐、Love number 等辅助数据 |
| `data/basin_mask/` | 46 | 海洋、流域和矢量掩膜 |
| `data/ddk_data/` | 8 | DDK 滤波矩阵 |
| `data/GIA/` | 4 | ICE-6G-C、ICE-6G-D、Caron 2018/2019 GIA |
| `data/L2_low_degrees/` | 9 | TN-11、TN-13、TN-14 低阶项 |
| `data/L2_SH_products/` | 180 | CSR/GFZ 球谐系数样例 |
| `data/Noah2.1/` | 24 | GLDAS Noah 2.1 辅助数据 |
| `data/topography/` | 13 | 地形和气压辅助数据 |
| `data/vgc_data/` | 3 | VGC 辅助数据 |
| `validation/` | 3 | 原版示例验证数据 |

关键文件已核对：

- `data/GIA/GIA.ICE-6G_D.txt`
- `data/GIA/GIA.Caron_et_al_2018.txt`
- `data/L2_low_degrees/TN-11_C20_SLR_RL06.txt`
- `data/L2_low_degrees/TN-13_GEOC_CSR_RL06.txt`
- `data/L2_low_degrees/TN-14_C30_C20_SLR_GSFC.txt`
- `data/ddk_data/Wbd_2-120.a_1d10p_4`
- `data/basin_mask/SH/Ocean_maskSH.dat`
- `data/basin_mask/Shp/bas200k_shp/bas200k_shp.shp`
- `data/auxiliary/GIF48.gfc`
- `data/auxiliary/LoveNumber.mat`

## 3. Level-2 GSM 覆盖范围

### 3.1 原版目录样例

`data/L2_SH_products/` 保留原版 SaGEA 的 CSR/GFZ/ITSG 样例和项目登记的 ICGEM 文件，用于原版读取接口、低阶项、滤波、球谐积分和掩膜测试。它不自动等同于目标 B 的完整研究窗口输入。

### 3.2 目标 B 本地 CSR DDK1 输入

- 只读目录：`D:\temp_sealevel_data\grace_csr_ddk1`。
- 发现文件数：101 个 `.gfc`，中心为 CSR，包含 GRACE 与 GRACE-FO，文件头记录已使用 DDK1 非各向同性滤波。
- 目标窗口：2013-11 至 2024-10，共 132 个月。
- 月份配准：文件名起止年积日的区间中点，映射为 `YYYY-MM`；重复月份按覆盖区间最长者保留，未发现重复目标月。
- 运行覆盖：101 个有效观测月，31 个缺测月；缺测月在 Level-3 固定轴中保留为 `valid_month=0` 和 NaN，未插值、未补零。
- 由于输入已经是 DDK1-filtered GFC，本次未再次施加 DDK；该输入状态与未经滤波的原始 L2 产品分开记录。
- 每个选中输入文件的 SHA-256、配置 SHA-256 和实际生效参数保存在 `results/target_b_l2_to_l3/target_b_run_manifest.json`。

本次已生成目标 B 的首轮自制 L3，但它是“本地 DDK1、缺测保留、未施加 GIA/OBD/GAX/几何/去混叠/泄漏校正”的阶段性产品，不宣称为完整 132 个月或论文最终 barystatic 结果。完整研究仍缺少 31 个目标窗口 CSR 月份、未经预滤波的原始 L2 对照、配套 GAX/AOD1B 和经登记的泄漏/GIA/OBD 辅助场。

### 3.3 目标 B 首轮输出

- Level-3 NetCDF：`results/target_b_l2_to_l3/target_b_custom_l3_201311_202410.nc`。
- 海洋平均 CSV：`results/target_b_l2_to_l3/target_b_ocean_mean_201311_202410.csv`。
- 月份登记 CSV：`results/target_b_l2_to_l3/target_b_l2_inventory_201311_202410.csv`。
- 运行 manifest：`results/target_b_l2_to_l3/target_b_run_manifest.json`。
- 方法记录：`docs/target_b_l2_to_l3_method.md`。

### 3.4 目标 B 多本地来源覆盖扩展

- 配置：`config/target_b_l2_to_l3_multisource.json`。
- 主来源：`data/L2_SH_products/ICGEM`，128 个 `GSM-2_*.gfc`，优先级 100。
- 只读回退：`D:\temp_sealevel_data\grace_csr_ddk1`，101 个 `.gfc`，优先级 10。
- 合并发现文件：229 个；目标窗口有效月份 108/132，缺失 24 个月。
- 入选来源：项目 ICGEM 107 个月；外部本地回退仅 `2017-06` 1 个月。
- `2015-04` 选择项目内自然月对齐的 `2015091-2015120` 候选；所有竞争候选均登记在 inventory。
- 参考场固定为 2013-11—2024-03 内 101 个有效月，确保与首版共同的 101 个月数值完全一致。
- 输出目录：`results/target_b_l2_to_l3_multisource/`；验证报告：`target_b_multisource_validation_report.md`。
- NetCDF 同时登记 `ocean_mask` 与 `land_mask`；陆地由 `distance_to_land_km == 0` 定义，海岸缓冲区保持为海陆掩膜均为 0，不并入陆地。
- 缺测继续保持 `valid_month=0` 和 NaN；未应用 GIA、OBD、GAX、几何、去混叠、泄漏校正或协方差传播。

### 3.5 目标 B Caron 2018 GIA 敏感性

- 配置：`config/target_b_l2_to_l3_multisource_gia_caron2018.json`。
- GIA：`data/GIA/GIA.Caron_et_al_2018.txt`；SHA-256 `514747ca3eafb8fdf57c4ef6cff58484e65b39bd0d77a1e66dd012d013094d31`。
- 单位与符号：完全归一化无量纲 Stokes 系数 `yr^-1`，从观测系数中相减。
- 一致滤波：GIA 趋势使用 DDK1 矩阵 `Wbd_2-120.a_1d14p_4`；GRACE DDK1 输入不重复过滤。
- 参考历元：固定参考期 101 个有效观测中点平均十进制年 `2019.525821395299`。
- 输出：`results/target_b_l2_to_l3_multisource_gia_caron2018/`；仍为 108 个有效月、24 个缺失月。
- 滤波后 GIA 开阔海洋平均趋势：`-1.2502033517 mm/yr`；GIA 产品相对无 GIA 基线趋势变化 `+1.2502033517 mm/yr`。
- 未加入 OBD、GAX、几何、泄漏或不确定度传播；不替代无 GIA 基线。

### 3.6 目标 B 椭球几何差分敏感性

- 配置：`config/target_b_l2_to_l3_multisource_geometry_ellipsoid.json`；SHA-256 `2d97972b87d044aa789968ad3415d70ae7657c6756b8f7986be828d19117a384`。
- 方法：原版 SaGEA 迭代几何算子的 `input + (Ellipsoid recovery - Sphere recovery)` 控制变量差分；内部 0.5° 网格、4 次迭代。
- 地形依赖：`data/topography/PHISFC_ERA5_invariant.nc`；SHA-256 `6814cb75d5d214665b73f9954bf81595ffc0d017051b3748e20f6bf71354c4fb`。
- Love-number 依赖：`data/auxiliary/LoveNumber.mat`；SHA-256 `462a8c8005e81bc79c10746007e8c7fabd07917a114b78013c7359450a0b8125`。
- 输出：`results/target_b_l2_to_l3_multisource_geometry_ellipsoid/`；108 个有效月、24 个缺失月，时间轴、掩膜和输入来源与无几何基线一致。
- 相对无几何基线：全球场差分 RMS `2.055844 mm`，开阔海洋均值差分 RMS `0.077967 mm`，海洋趋势变化 `+0.0140389 mm/yr`。
- manifest 中 108 个 L2 与 6 个辅助输入共 114 个 SHA-256 均复核匹配；几何算子耗时 `373.1565 s`。
- 验证报告：`target_b_geometry_ellipsoid_validation_report.md`；未加入 GIA、OBD、GAX、泄漏或不确定度传播，不替代无几何基线。

### 3.7 目标 B Forward Modeling 泄漏敏感性

- 实现：`pysrc/level3/leakage.py`；诊断写入：`pysrc/level3/leakage_io.py`。
- DDK1 矩阵：`data/ddk_data/Wbd_2-120.a_1d14p_4`；SHA-256 `b61e429b04965c6b7b436e437bade8e5420d63f5365736070678a095c71b5e81`。
- max50 配置：`config/target_b_l2_to_l3_multisource_leakage_forward.json`；SHA-256 `687639db55eded585118d1774a6c1afc77d1e77e2d27344a20fed6cc04b43bcf`。108 月残差均下降但 0 月达到原双判据，目录只作为非收敛诊断。
- max500 配置：`config/target_b_l2_to_l3_multisource_leakage_forward_max500.json`；SHA-256 `6b71afd261e7c5c7177d5f4c241c886bdaa208054f2f145c69518b2e1062c1d2`。相对 max50 只改变上限、版本和输出命名，不放宽阈值。
- 收敛输出：`results/target_b_l2_to_l3_multisource_leakage_forward_max500/`；108/108 有效月收敛，24 个缺测月保持 NaN，迭代范围 251–467。
- 最终相对残差 `0.001843–0.004334`，最终相对更新 `0.00099698–0.00099999`；Forward Modeling 核心耗时 `362.9400 s`。
- 相对无泄漏基线，海洋均值校正量 RMS `0.772465 mm`，趋势变化 `-0.0021192 mm/yr`；全部逐点恒等式和 113 个输入/辅助哈希复核通过。
- 未加入 GIA、几何、OBD、GAX、去混叠或不确定度传播；不替代无泄漏基线。

### 3.8 目标 B 开阔海洋约束迭代重建敏感性

- 实现：`pysrc/level3/leakage.py` 的 `iterative_ocean_reconstruction`；使用同一 DDK1 矩阵和 300 km 开阔海洋掩膜。
- max500 配置 SHA-256：`f2fdcc9d9bad7e35d3a6d49c397884b735c219b6ffd4a4816ce47c0f0a291765`；108 月残差全部单调下降但 0 月收敛，只保留为诊断。
- max4000 配置 SHA-256：`d9661f42651d06fd2dd1737b5ca8368cf5db3b0ef656d60e0ef2a668d769ef88`；相对 max500 只改变上限、版本和输出名。
- max4000 输出：`results/target_b_l2_to_l3_multisource_leakage_iterative_ocean_max4000/`；108/108 月在 548–3148 次内收敛，中位数 1586，24 个缺测月保持 NaN。
- 113 个输入/辅助哈希、逐点恒等式、历史行数和 CSV 网格积分审计通过；核心耗时 `1233.191 s`。
- 物理稳定性未通过：海洋调整场 RMS `806.547 mm`、极值 `-69.339–110.965 m`，趋势相对基线改变 `-0.495697 mm/yr`。该结果不得进入 Figure 主曲线，只登记为方法失败模式和离散度诊断。
- 仍未传播协方差；两方法差异不能冒充完整观测不确定度。

### 3.9 目标 B 不确定度资格审计

- 108/108 个已选 DDK1 GFC 含 formal sigma 且三角区为有限值；正值范围 `3.204e-14–7.778e-12`。
- 108/108 个文件头同时声明 `errors unchanged (no error propagation applied)`，这些 sigma 不能视为 DDK1 后误差。
- degree-1 与 C20 替换源没有登记相容协方差，完整系数协方差也未登记。
- 实现：`pysrc/level3/uncertainty.py`；报告：`results/target_b_l2_to_l3_multisource/target_b_uncertainty_eligibility_report.md`。
- 决策：不生成伪精确误差场，Level-3 `uncertainty` 继续为 NaN；方法差异只作为敏感性。

## 4. Figure 1 研究窗口

- 研究窗口：2013-11 至 2024-10。
- 总月份数：132。
- 2014-10 至 2015-12：论文事件窗口，15 个月。
- 2023-05 至 2023-12：论文事件窗口，8 个月。
- 论文 Figure 1(a)：去季节后的 GMSL、barystatic、steric 和两者之和。
- 论文 Figure 1(b)：在同一序列上进一步去除线性趋势，并使用 3 个月滑动平均。

## 5. Figure 1 原始及参考产品

| 产品 | 原始文件 | 当前状态 | 论文复现用途 |
|---|---|---|---|
| C3S/CMEMS L4 月 SLA | `data/external_downloads/altimetry/C3S_CMEMS_L4_monthly_SLA_201311_202410.nc`；与 `D:\temp_sealevel_data\cmems_sla_monthly.nc` 内容一致 | 132 个月，0.25° 网格，单位 m | 可作 CMEMS/C3S 等效参考；不能保证是论文精确版本 |
| CMEMS GMSL indicator | `D:\temp_sealevel_data\cmems_gmsl_indicator.nc` | 1999-02 至 2026-01，含全球平均序列及不确定度 | 可作 GMSL 交叉验证，不替代月网格计算 |
| CSR RL06.3 mascon | `data/external_downloads/mascon/CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc` | 已可读取并计算单中心海洋质量参考序列；研究期有缺测 | 可作 CSR 单中心参考；尚未纳入三中心平均 |
| JPL RL06.3Mv04 | 原始：`D:\temp_sealevel_data\jpl_mascon_rl06.3_v04.nc` | NASA/JPL 真实产品；CRI 未应用；研究期 108 个唯一月份；2015-04 有重复记录 | 可作 JPL 参考，但不是论文精确 CRI 版本 |
| JPL RL06.3Mv04 CRI | `C:\Users\Alan\Downloads\GRCTellus.JPL.200204_202605.GLO.RL06.3M.MSCNv04CRI.nc` | 有效 NetCDF4；官方元数据确认 CRI；研究期 109 条记录、108 个唯一月份 | 已纳入 JPL CRI 标准化参考；仍需与论文版本说明逐项核对 |
| GSFC RL06 v2.0 | 原始：`D:\temp_sealevel_data\gsfc_mascon_rl06v2.0.nc`；另有 ASCII 压缩包 | 研究期有 109 条源记录、108 个唯一月份；存在缺测和重复月份 | 可作 GSFC 单中心参考；纳入三中心平均前需处理重复月 |
| IAP 0–2000 m steric | 原始目录：`D:\temp_sealevel_data\IAP_Steric` | 2013-11 至 2023-09，共 119 个月；缺 13 个月 | 可作 IAP 总 steric 网格；不能替代 IAP thermosteric |
| IAPv4 月度温度场（补下载） | `data/external_downloads/IAPv4/` | 已下载 2023-10 至 2024-06 共 9 个月，约 277.6 MB；2024-07 至 2024-10 官方目录返回 404 | 只能作为热成分计算的部分输入；尚无盐度场，不能直接形成 IAP thermosteric |
| SIO steric/OHC | `data/external_downloads/steric/SIO_steric_ohc_time_series_20040116_20260617_v01.nc` | thermosteric、halosteric、totalsteric 均覆盖研究窗口，单位 m | 可直接作 SIO steric 参考分量 |
| SIO/RG 2004–2018 Argo 温盐气候态 | 官方直链已确认；项目目录 `data/external_downloads/Argo_SIO_climatology/` | 温度约 695.5 MB、盐度约 460.9 MB；本机传输停滞，尚未获得完整文件 | 用于独立温盐/steric 交叉验证，不是 Figure 1 主序列 |
| NOAA GMSL | `data/external_downloads/altimetry/NOAA_GMSL_slr_sla_gbl_free_all_66.csv` | 保留多个卫星列；部分卫星覆盖研究窗口 | 可作独立卫星测高交叉验证，不直接拼接为主序列 |
| DMI | `data/external_downloads/climate_indices/DMI_NOAA_PSL_HadISST_monthly.nc` | 研究窗口完整 | 用于 Figure 6 及机制分析 |
| ONI | `data/external_downloads/climate_indices/` 下 NOAA CPC ONI 文件 | 已登记 | 用于 Figure 2；不属于 Figure 1 主输入 |

Figure 1 参考版当前输出：

- 原始缺测诊断图：`results/figures/figure01_reference.png`。
- SSA 填补版：`results/figures/figure01_reference_ssa_gapfilled.png`。
- SSA 填补版指标：`results/tables/figure01_ssa_gapfilled_metrics.csv`。
- SSA 填补版说明：`results/reports/figure01_ssa_gapfilled_notes.json`。
- 论文预处理校正版：`results/figures/figure01_reference_ssa_gapfilled_corrected.png`。
- 论文预处理校正版指标：`results/tables/figure01_ssa_gapfilled_corrected_metrics.csv`。
- 论文预处理校正版说明：`results/reports/figure01_ssa_gapfilled_corrected_notes.json`。
- SSA 质量序列方法记录：`results/reference_baseline/ssa_gap_filling_method.md`。

Figure 1 自制 Level-3 接入版（最终候选 v2）：

- 配置：`config/figure01_custom_l3.json`；构建脚本：`pysrc/reference_products/build_figure1_custom_l3.py`。
- 主图：`results/figure01_custom_l3_20260820_v2/figure01_custom_l3.png` 和同名 PDF。
- 质量输入叠加图：`results/figure01_custom_l3_20260820_v2/figure01_custom_l3_overlay.png` 和同名 PDF。
- 绘图数据：`results/figure01_custom_l3_20260820_v2/figure01_custom_l3_data_201311_202410.csv`。
- 指标、配置快照、方法报告和运行清单位于同一目录；清单登记 3 个输入和全部 8 个非清单输出的 SHA-256。
- 固定月轴为 132 个月；自制 L3 有效 108 个月、缺测 24 个月；3 月中心平滑后可绘 88 个月；不使用 SSA 或插值。
- 自制 L3 趋势为 `0.151992 mm/yr`，Forward 敏感性为 `0.141119 mm/yr`，二者相关系数 `0.995969`、RMSE `0.701636 mm`。
- CSR mascon 参考趋势为 `1.231482 mm/yr`，与自制 L3 存在显著长期差异；该版本是阶段性本地 DDK1 反演，不是 CSR mascon 的等价替代。
- `custom mass + steric` 不含 OBD；物理失稳的 iterative-ocean 方法未进入任何 Figure 1 输入。

## 6. 项目已生成的标准化结果

### 6.1 JPL

- 输出：`results/reference_baseline/jpl_ocean_mass_reference_201311_202410.csv`
- 元数据：`results/reference_baseline/jpl_ocean_mass_reference_metadata.json`
- 脚本：`pysrc/reference_products/standardize_jpl_mascon.py`
- 输入单位：cm 等效水高。
- 输出单位：mm 等效水高。
- 研究期有效唯一月份：108/132。
- 缺测月份保留为 `nan`，未插值。
- 产品明确标记为非 CRI。
- 2015-04 的两条记录按“保留第一条、记录第二条舍弃原因”处理。
- 使用项目 0.25° 300 km 诊断掩膜重映射至 JPL 0.5° 网格，并采用余弦纬度加权。

### 6.1.1 JPL CRI

- 输入：`C:\Users\Alan\Downloads\GRCTellus.JPL.200204_202605.GLO.RL06.3M.MSCNv04CRI.nc`
- 输出：`results/reference_baseline/jpl_cri_ocean_mass_reference_201311_202410.csv`
- 元数据：`results/reference_baseline/jpl_cri_ocean_mass_reference_metadata.json`
- 产品：NASA/JPL RL06.3Mv04-CRI，输入单位 cm，输出单位 mm。
- 研究窗口：132 个月；109 条源记录；去重后 108 个唯一月份。
- 重复记录：`2015-04` 保留第一条；缺测月份保留为 `nan`，未插值。
- 使用同一项目 300 km 诊断掩膜和余弦纬度面积加权，便于与非 CRI 版本及 CSR/GSFC 对比。

### 6.2 GSFC

- 输出：`results/reference_baseline/gsfc_ocean_mass_reference_201311_202410.csv`
- 元数据：`results/reference_baseline/gsfc_ocean_mass_reference_metadata.json`
- 脚本：`pysrc/reference_products/standardize_gsfc_mascon.py`
- 输入单位：cm 等效水高；输出单位：mm 等效水高。
- 研究窗口输出统一为 132 个月；原始源记录 109 条，去重后 108 个唯一月份；重复月份按“保留第一条”处理，其余月份显式保留为 `nan`。
- 缺测月份未插值。
- GSFC 重复月份已处理，并已加入项目 300 km 掩膜在 mascon 中心的最近邻映射；形成三中心平均前仍需统一网格、基准期和不确定度定义。

### 6.4 CSR/JPL/GSFC 诊断对比

- 脚本：`pysrc/reference_products/compare_mass_products.py`
- 指标：`results/reference_baseline/mass_product_metrics_201311_202410.csv`
- 两两对比：`results/reference_baseline/mass_product_pairwise_201311_202410.csv`
- 报告：`results/reference_baseline/mass_product_comparison_201311_202410.md`
- 当前结果显示 JPL CRI 与非 CRI 序列高度接近；CSR 与 JPL/GSFC 存在约 11.5 mm 的平均基准差，不能在未统一基准和掩膜前直接解释为产品质量差异。
- 当前只做诊断对比，尚未计算 CSR/JPL/GSFC 三中心平均。

### 6.3 IAP

- 输出：`results/reference_baseline/iap_steric_2000m_reference_201311_202410.nc`
- 元数据：`results/reference_baseline/iap_steric_2000m_reference_metadata.json`
- 脚本：`pysrc/reference_products/standardize_iap_steric.py`
- 网格：180 × 360，原生 1°。
- `steric_2000m`：IAP `SSL_2000m`，米转毫米，保留原始有限网格。
- `steric_2000m_open_ocean_300km`：应用项目 300 km 诊断掩膜后的变量。
- 已有月份：119/132。
- 缺失月份：2023-10 至 2024-10，共 13 个月；未插值。已另行取得 2023-10 至 2024-06 的 IAPv4 温度场，但它不是 `SSL_2000m` 总 steric 文件，不能填补该总 steric 序列。
- 原始 IAP 在黑海及邻近区域存在异常大有限值；原始变量保留，开阔海洋变量用于全球海洋统计。
- 该文件只有总 steric 网格，没有 IAP thermosteric 独立变量。

### 6.5 SSA 缺测填补

- 方法记录：`results/reference_baseline/ssa_gap_filling_method.md`
- 算法脚本：`pysrc/reference_products/ssa_gap_filling.py`
- 执行脚本：`pysrc/reference_products/run_ssa_gap_filling.py`
- 依据：Xie and Yi (2025) 的 SSA-filling-a / SSA-filling-b；参数来源于 Yi and Sneeuw (2021) 配套公开程序。
- 四条海洋质量序列均由 108/132 个有效月份补为 132/132；24 个原始缺测月份均有 SSA 填充值，原始观测值保留不变。
- `SSA-filling-a` 使用 `M=24, K=10`；`SSA-filling-b` 在 `M=24:12:96`、`K=[1,2,4,6,8,10,12]` 中通过当前窗口的伪缺口交叉验证选参。
- 填补结果另存为 `*_ssa_gapfilled_201311_202410.csv`，不覆盖原始标准化文件；每个文件包含原始值、填补值、填补类别和所用参数。
- 当前是对标准化全球海洋一维序列的 SSA 填补，不是论文原始逐网格 mascon 填补；不能将其等同于论文发布的完整网格产品。

### 6.6 Jin 2026 规格验证模式

- 配置：`config/jin2026.yaml`
- 处理脚本：`pysrc/sea_level_budget/validation.py`
- 输出目录：`results/jin2026_spec_validation/`
- Figure 1：`results/jin2026_spec_validation/figure01_jin2026_spec_validation.png`
- 月度收支表：`results/jin2026_spec_validation/jin2026_budget_timeseries_201311_202410.csv`
- 指标表：`results/jin2026_spec_validation/jin2026_budget_metrics.csv`
- 运行清单：`results/jin2026_spec_validation/run_manifest.json`
- 当前模式使用 C3S 等效 GMSL、SIO total steric 和三条 SSA 填补 mascon；三中心 barystatic 主线取 CSR/JPL/GSFC 中位数。
- 已按规格执行 GMSL `-0.5 mm/yr` 修正、月气候态去季节和 3 个月中心滤波；滤波两端月份标记为无效。
- 当前为 validation mode：完整 IAPv4 TEOS-10 温盐输入、论文精确 CMEMS 文件和 AR(1) 趋势不确定度尚未完成，因此不宣称完整复现论文。

### 6.7 Figure 3 区域陆地水储量输入

- Xie–Yi 数据目录：`data/external_downloads/xie_yi_2025_mascon_gapfilled/`
- 来源：Figshare 数据集 v2，DOI `10.6084/m9.figshare.25805092.v2`；许可为 CC BY 4.0。
- 文件：CSR RL06.2、JPL CRI RL06.1v03、GSFC RL06v2.0 三个 NetCDF；每个文件只包含 35 个 SSA 重建缺测月，单位为 cm 等效水高。
- 使用规则：原始中心 Mascon 的有限观测值始终优先；Xie–Yi 网格仅填补原始缺测月，且不得用于 2022-12 之后。
- Natural Earth 数据目录：`data/external_downloads/natural_earth/admin_0_countries/`。
- Natural Earth 发布物内部版本为 5.1.1；用于构建六大洲、排除格陵兰和南极洲的可复现区域代理，不宣称等同于 Jin 论文未公开的精确大陆边界。
- 两套来源均有 `source_manifest.json`，记录来源版本、许可、字节数和 SHA-256。

## 7. 海洋掩膜审计

- 项目诊断掩膜：`results/reference_baseline/ocean_mask_300km_025deg.nc`
- 诊断元数据：`results/reference_baseline/ocean_mask_300km_025deg.json`
- 来源审计：`results/reference_baseline/ocean_mask_audit.json`
- 文字报告：`results/reference_baseline/ocean_mask_source_audit.md`
- 生成脚本：`pysrc/reference_products/build_ocean_buffer_mask.py`
- 审计脚本：`pysrc/reference_products/audit_ocean_masks.py`

审计结论：

- `D:\temp_sealevel_data\ocean_mask_025deg_300kmbuf.npy` 无坐标、来源和生成元数据，且大量网格不满足项目 300 km 距离条件，不批准替换。
- `D:\temp_sealevel_data\sagea_ocean_mask_1deg_300km.npy` 无来源元数据且为 1° 网格，不适合论文 0.25° 输入。
- `basin_mask.zip/SH/Ocean_maskSH.dat` 与 `aux_data.zip/auxiliary/ocean360_grndline.sh` 内容和 SHA-256 完全相同，是球谐边界辅助文件，不是布尔网格掩膜。
- Natural Earth land 版本 4.1.0 是 110m 陆地多边形，不等同于论文 300 km 海岸缓冲掩膜。
- 项目诊断掩膜可作为可复现的统一处理口径，但不能声称等同论文未公开的精确海岸线实现。

## 8. Figure 1 当前复现状态

当前已有数据可以生成“方法一致的参考重建版”，但不能完整复现论文原 Figure 1。

仍缺失或未完成的关键事项：

1. IAP thermosteric 数据或可复现的 IAP 温度场、盐度场到 thermosteric 的完整计算输入；目前仅补到 9 个月温度场。
2. IAP 2023-10 至 2024-10 的 13 个月数据。
3. CSR/JPL/GSFC 在统一 1° 网格、统一 300 km 掩膜、统一基准期下的三中心平均。
4. GSFC 重复月份和三中心缺测规则。
5. 三中心质量产品的不确定度合成，以及 Figure 1 灰色不确定度带。
6. 论文精确海岸线实现；目前只有项目诊断掩膜。
7. 论文 Figure 1 的完整校正流程：季节项处理、2016 年起 GMSL `-0.5 mm/yr` 修正、3 个月滑动平均、垂直显示偏移和闭合误差带。
8. 论文对应版本的 CMEMS 和 IAP thermosteric 输入文件的版本确认；JPL CRI 已取得，但仍需在论文复现记录中核对版本对应关系。

已完成的阶段性工作：

- 四条质量序列的 24 个缺测月份已按 SSA-filling-a/b 填补，原始观测值保持不变。
- Figure 1 已生成一版使用 SSA 填补 CSR 海洋质量序列的参考图；该图仍不是论文原图的完整复现。
- Figure 1 另已生成按论文原文去季节、并从 2016-01 起对 GMSL 应用 `-0.5 mm/yr` 修正的校正版；该版本仍使用等效 C3S 网格 GMSL 和 CSR 单中心质量序列。

当前 `results/figures/figure01_reference.png` 只能作为早期诊断图，不应标记为论文 Figure 1 的完整复现。

## 9. 其他待补数据

- 目标 B 多来源扩展后仍缺少 24 个目标窗口 CSR 月份，以及未经 DDK1 预滤波的完整原始 Level-2 对照产品。
- GRACE-FO 低阶项、GAX/AOD1B 和配套辅助场。
- ERA5 月降水、蒸发、径流、潜热、感热、短波和长波通量。
- 完整 IAPv4 温度、盐度和 OHC 产品，用于 Figure 5–7。
- 完整 SIO 2004–2018 Argo 温盐网格，用于独立 steric 交叉验证。
- 论文精确使用的 OBD、GIA、缺测重建和质量控制信息。

## 10. 验证记录

科学测试目录：`tests/scientific/`

最近一次全套测试结果：

```text
14 passed in 3.41s
```

该结果验证当前标准化、掩膜审计和参考数据处理脚本的基本行为，不代表已经完成论文 Figure 1 的完整复现。
