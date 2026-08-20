# Jin et al. (2026) Figure 1/2 and custom GRACE Level-3 data release

本目录汇总了本项目复现 Jin et al. (2026, JGR Oceans) Figure 1 和 Figure 2 时实际使用的标准化月序列，以及使用本地 GRACE/GRACE-FO Level-2 球谐产品生成的自制 Level-3 正式产品。

## 目录

- `paper_reproduction/figure1/`：论文 Figure 1 最终复现的月序列、指标、图件、报告和运行清单。
- `paper_reproduction/figure2/`：论文 Figure 2 最终复现的绘图月序列、事件阶段汇总和图件。
- `paper_reproduction/reference_series/`：Figure 1 参考 GMSL、steric 和 CSR/JPL/GSFC 质量海平面标准化序列。
- `paper_reproduction/climate_indices/`：Figure 2 使用的 NOAA CPC ONI v5 本地文件。
- `custom_l3/product/`：原版 SaGEA 控制链加 Improved MSSA 重建得到的自制 Level-3 NetCDF、球谐系数、海洋平均时间序列、方法报告和验证文件。
- `custom_l3/figure1_v6/`：使用自制 Level-3 替换质量海平面后生成的 Figure 1 v6 数据、指标、图件、配置快照和运行清单。
- `config/`：Figure 2 原始配置副本及指向本发布包 Figure 1 v6 数据的配置。

## 自制 Level-3 口径

- 时间范围：2013-11 至 2024-10，共 132 个月。
- 本地 Level-2 实际观测月：108。
- Improved MSSA 重建月：24；窗口 `M=60`，秩 `rank=14`。
- GIA：Caron et al. (2018)，与 DDK1 状态匹配。
- 海洋掩膜：距陆地至少 300 km 的开阔海洋。
- Figure 1 v6 趋势：GMSL `3.21364551791 mm/yr`；自制质量海平面 `1.38746221764 mm/yr`；steric `1.48777813730 mm/yr`；质量加 steric `2.87524035494 mm/yr`。

## Git LFS

自制 Level-3 NetCDF 和系数 NPZ 通过 Git LFS 保存。克隆后执行：

```bash
git lfs pull
```

## 完整性与范围

`SHA256SUMS.txt` 登记本目录内发布文件的 SHA-256。为控制仓库规模，本包不包含原始 Level-2 `.gfc`、第三方原始 mascon/steric 网格、Python 运行环境或临时缓存；其来源和处理参数记录在各运行清单、配置快照与方法报告中。

Figure 1 v6 原始 CSV 的表头中 `custom_l3_reconstructed_month` 重复两次，且 132 行对应值逐行一致。发布包保留原始 CSV 以匹配运行清单，同时提供后缀为 `_normalized.csv` 的规范化副本；该副本只删除第二个重复字段，不改变任何科学数值，Figure 2 v6 配置默认读取规范化副本。
