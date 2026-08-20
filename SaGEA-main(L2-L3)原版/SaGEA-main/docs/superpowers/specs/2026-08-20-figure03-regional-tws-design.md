# Figure 3 区域 TWS 统一处理与论文复现设计

## 1. 目标

建立一套只依赖配置切换输入产品的 Figure 3 流程。第一版使用 Jin et al. (2026) 采用的 CSR、JPL、GSFC Mascon 与 Xie–Yi (2025) 缺测重建数据，生成论文口径 Figure 3。后续本地自制 Level-3 产品必须通过同一套区域掩膜、质量积分、时间处理、指标和绘图代码生成独立版本，不能另写一套科学定义不同的脚本。

本任务不修改 Figure 1、Figure 2、Level-2 到 Level-3 主处理链，也不调用项目边界之外的代码或输出。

## 2. 论文图定义

Figure 3 固定为四个面板：

1. (a) 2014–2016 El Niño 窗口的 Total、Africa、Asia、Europe、North America、South America、Oceania 区域 TWS 年际序列。
2. (b) 2014-10 到 2015-12 的区域 TWS 端点差值。
3. (c) 2023–2024 El Niño 窗口的同类区域 TWS 年际序列。
4. (d) 2023-05 到 2023-12 的区域 TWS 端点差值。

灰色阴影使用上述两个发展阶段。区域 TWS 增加为正，减少为负。阶段差值统一定义为：

```text
change_mm = value_at_end - value_at_start
```

论文所述“减少 6.37 mm”在指标表中保存为 `change_mm = -6.37`，并可另存 `reduction_magnitude_mm = 6.37`，不得混淆符号。

## 3. 总体架构

```text
论文 Mascon / 自制 Level-3
        ↓
输入适配器：MonthlyGridSeries
        ↓
统一经纬度、单位、月份与来源状态
        ↓
1° 六大陆陆地掩膜与陆侧 300 km 海岸缓冲
        ↓
逐中心区域质量积分和全球海平面等效换算
        ↓
三中心平均或单个自制 L3 区域表
        ↓
去季节、去趋势、3 月中心平均
        ↓
Figure 3、指标、方法报告和运行清单
```

绘图器只读取标准区域表，不读取原始 NetCDF。产品格式、坐标顺序和单位差异全部在适配器层解决。

## 4. 标准输入接口

每个输入适配器输出同一结构：

```text
MonthlyGridSeries
  source_id: str
  months: 一维 YYYY-MM 月份数组
  lat: 一维纬度中心
  lon: 一维经度中心
  ewh_mm: time × lat × lon 等效水高，单位 mm
  valid_month: 一维布尔数组
  month_status: observed / reconstructed / missing
  metadata: 产品、版本、原始单位、改正状态、文件哈希
```

接口必须拒绝重复月份、非单调坐标、无法识别的单位、经纬度维度不匹配和没有来源状态的重建值。

### 4.1 论文 Mascon 适配器

支持：

- CSR 0.25° `lwe_thickness` NetCDF；
- JPL CRI 0.5° `lwe_thickness` NetCDF；
- GSFC 0.5° `lwe_thickness` NetCDF；
- Xie–Yi (2025) CSR/JPL/GSFC 逐网格 SSA 重建产品。

论文主时间轴固定为 2013-11 至 2024-10。原始中心产品中的有效观测月优先保留；只有缺测月可由同中心 Xie–Yi 重建值替换。Xie–Yi 产品截至 2022-12，2023-01 以后只使用本地官方中心产品。

若重建产品与当前中心产品参考基准不同，使用共同观测月估计逐格点加性偏差。偏差对齐只能改变零点，不能改变趋势、季节幅值或年际信号。运行清单保存重叠月份、偏差场统计量和拼接前后 RMS。没有足够共同月份时停止运行。

不同中心分别完成读取、重建月替换、1°重采样和区域积分。不得先把不同原生网格直接平均。最终论文区域序列是三个中心区域序列的逐月算术平均，并同时保存中心间标准差、最小值、最大值和有效中心数。

### 4.2 自制 Level-3 适配器

自制产品通过配置登记以下变量名：

```text
time
latitude
longitude
field
land_mask
valid_month
```

`field` 必须是表面质量或等效水高网格，并提供可验证单位。适配器允许 NetCDF 或项目已登记的 NPZ 回退格式。缺测月默认保留 NaN，不自动调用 SSA、插值或邻月替代。若以后批准重建方案，必须作为新的显式配置和结果版本加入。

自制产品使用与论文版本相同的目标网格、区域掩膜、面积换算和时间处理。这样两版 Figure 3 的差异只来自输入质量场及其已登记的改正，不来自绘图或区域定义。

## 5. 数据下载与登记

Xie–Yi 官方数据从以下 DOI 获取：

```text
https://doi.org/10.6084/m9.figshare.25805092.v2
```

下载目录固定为：

```text
data/external_downloads/xie_yi_2025_mascon_gapfilled/
```

下载流程先保存 Figshare 文件清单，再下载原文件。原文件不改名、不覆盖、不重新压缩。登记内容至少包括 DOI、版本、下载时间、文件名、字节数、发布方校验值和本地 SHA-256。若站点只允许浏览器下载，下载完成后仍由脚本执行同样的文件登记。

Natural Earth 边界文件保存到：

```text
data/external_downloads/natural_earth/admin_0_countries/
```

同样保存来源 URL、版本、许可说明和 SHA-256。

## 6. 区域掩膜

目标网格固定为原项目 1°规则格点：

```text
latitude: -89.5 ... 89.5
longitude: -179.5 ... 179.5
```

Natural Earth 国家多边形按大陆属性合并为 Africa、Asia、Europe、North America、South America、Oceania。格陵兰从 North America 中剥离，南极洲不进入大陆 TWS。海岛按其大陆属性归属，但只有满足格点中心落在区域内且通过海岸缓冲的单元才参与积分。

论文只公开了 Figure S2 的分区图，没有发布精确多边形。因此 Natural Earth 是可复现代理，不得在报告中声称与论文边界逐格点完全相同。

陆侧 300 km 缓冲定义为：距全球海岸线小于 300 km 的陆地格点不参与区域积分。缓冲距离按球面大圆距离计算，不能用经纬度平面距离近似。输出掩膜 NetCDF 必须包含：

```text
region_id
land_mask
distance_to_coast_km
coastal_buffer_excluded
cell_area_m2
```

各区域掩膜互斥，Total 是六大陆有效格点的并集，不单独重新生成边界。

## 7. 区域质量与单位

每个区域的质量异常为：

```text
regional_mass_kg = Σ(EWH_m × water_density_kg_m3 × cell_area_m2)
```

Figure 3 使用全球平均海平面等效毫米：

```text
regional_esl_mm = regional_mass_kg / (water_density_kg_m3 × global_ocean_area_m2) × 1000
```

全球海洋面积及其来源写入配置和运行清单。论文版、自制 L3 版必须使用同一个分母。不能使用区域内 EWH 简单平均代替质量积分，也不能用某一产品的原生 `land_mask` 取代统一大陆掩膜；原生掩膜只用于产品质量控制。

## 8. 时间处理

每个中心先生成 2013-11 至 2024-10 的月度区域原始序列，再执行：

1. 对每个区域、每个中心，以完整研究期月气候态去季节。
2. 对去季节序列拟合并移除普通最小二乘线性趋势。
3. 计算 3 月中心滑动平均，只在三个相邻月均为有限值时输出。
4. 最后计算三中心平均和离散度。

三中心平均放在各中心时间处理之后，便于保留中心差异并防止中心缺测改变季节项。论文版两个事件窗口内部必须连续；若仍有缺月则停止绘图。自制 L3 版允许 NaN，但图、指标和报告必须列出缺测月份及其传播到平滑序列后的月份。

## 9. 配置与输出版本

论文版配置建议命名：

```text
config/figure03_paper_mascon.json
```

自制 L3 配置模板：

```text
config/figure03_custom_l3.template.json
```

论文版输出目录：

```text
results/figure03_paper_mascon_20260820_v1/
```

每次运行至少输出：

```text
figure03_paper_mascon.png
figure03_paper_mascon.pdf
figure03_plotting_data.csv
figure03_regional_by_center.csv
figure03_metrics.csv
continent_masks_1deg.nc
figure03_config_snapshot.json
figure03_method_report.md
figure03_run_manifest.json
```

绘图 CSV 保存图中实际使用值；中心表保存时间处理前后的各中心区域值和 `observed/reconstructed/missing` 状态；指标表同时保存有符号阶段差值、减少量绝对值、论文参考值和差异。

## 10. 论文参考验收值

论文值只作为独立比较，不参与算法调节：

| 阶段 | 区域 | 有符号变化量 (mm) |
|---|---|---:|
| 2014-10—2015-12 | Total | -6.37 |
| 2014-10—2015-12 | Africa | -1.72 |
| 2014-10—2015-12 | North America | +0.82 |
| 2014-10—2015-12 | South America | -3.25 |
| 2023-05—2023-12 | Total | -4.42 |
| 2023-05—2023-12 | Africa | +0.50 |
| 2023-05—2023-12 | North America | -1.30 |
| 2023-05—2023-12 | South America | -3.10 |

Asia、Europe 和 Oceania 在主文中没有全部给出精确文本数值，不从图像反向读取后伪装成论文表值。

## 11. 图形要求

图形结构、颜色和坐标范围尽量贴近论文，但科学数据不为视觉贴合而平移或缩放。

- (a)、(c) 使用相同区域颜色和标记，Total 加粗。
- (b)、(d) 条形颜色与区域时间序列一致。
- 灰色阴影准确覆盖两个发展阶段。
- 所有面板单位为 mm，阶段差值保留正负号。
- 图注说明三中心平均、海岸缓冲、重建月份和 Natural Earth 代理边界。

## 12. 失败条件和警告

以下情况停止论文版运行：

- 输入文件或哈希与配置登记不一致；
- 时间轴重复、倒序或缺少目标月份；
- 单位或正方向不明确；
- 重建值覆盖了有效观测月；
- 三中心区域序列在事件窗口内部仍有缺值；
- 区域掩膜重叠、为空或不满足 Total 恒等式；
- 区域和中心平均的逐点恒等式超出数值容差；
- 重复运行产生不同哈希。

以下情况输出显式警告但允许完成：

- Natural Earth 代理边界与论文 Figure S2 不能证明逐格一致；
- 不同 Mascon 版本拼接存在非零基准偏差；
- 论文参考差值超出配置中的诊断容差；
- 部分月份只有两个中心有效，但事件窗口仍连续。

## 13. 测试与验收

测试采用小型合成网格，不依赖大体积真实数据：

1. 单位换算、经纬度方向和经度范围转换测试。
2. 已知格点面积与区域积分解析值测试。
3. 六大陆掩膜互斥、Total 并集和格陵兰/南极排除测试。
4. 300 km 陆侧缓冲边界测试。
5. 观测月优先、只替换缺测月和拼接基准对齐测试。
6. 三中心平均、有效中心数和离散度测试。
7. 去季节、去趋势、3 月中心平滑顺序测试。
8. 阶段端点差值符号测试。
9. 论文版缺月失败、自制 L3 版保留 NaN 测试。
10. PNG、PDF、CSV、指标、配置、报告和 manifest 完整性测试。

真实数据验收还包括：输入与输出 SHA-256、事件窗口连续性、各中心覆盖月份、重建月份清单、区域面积、逐点恒等式、重复运行一致性，以及 PNG 和渲染 PDF 的人工检查。

## 14. 明确排除项

- 不使用 Figure 1/2 的一维海洋质量序列生成 Figure 3。
- 不把 OBD 加入陆地 TWS。
- 不对已经完成 GIA 改正的 Mascon 再施加 GIA。
- 不把格陵兰和南极冰盖质量并入六大陆 TWS。
- 不用插值、补零或邻月复制替代 Xie–Yi 重建值。
- 不为了命中论文柱高而修改面积、符号、基准或时间端点。
- 本阶段不实现 Figure 4。
