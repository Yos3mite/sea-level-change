# 可审计的全球平均海平面主流程

这套流程是新增的干净主流程，不调用、覆盖或修改 `reproduce_figure1` 中的历史复现脚本。正式结果采用标准全球海洋域，不输出 `66°S–66°N` 的 Jin 口径。当前版本处理测高 GMSL、GRACE/GRACE-FO 海洋质量和由同一组球谐系数计算的 OBD；比容项暂缓，因此不生成预算闭合残差。

## 数据与口径

- 测高：`SEALEVEL_GLO_PHY_CLIMATE_L4_MY_008_057` 对应的 C3S/CMEMS 双卫星月平均 SLA，本地文件为 `D:/temp_sealevel_data/cmems_sla_monthly.nc`，时间为 2013-11 至 2024-10。流程直接读取月度产品，不把日数据重新聚合成月数据。
- 测高质量控制：可与 `D:/temp_sealevel_data/cmems_gmsl_indicator.nc` 中的官方面积平均指标比较。官方指标已经包含 GIA 校正、季节调整和六个月滤波，不可当作原始 SLA 面积平均后再次加 GIA。
- GRACE：ICGEM 官方 CSR RL06；GRACE-FO：CSR RL06.3。只选 `BA01`（60 阶）`DDK1` 文件，每个日历月最多保留一个弧段，不插值任务缺月。
- 低阶项与 GIA：TN-11、TN-13 CSR RL06.2、TN-14，以及 Caron et al. (2018) GIA，均取自本地 SaGEA 数据目录。
- 海陆与海岸距离：Climate Data Toolbox 的 `land_mask.mat` 和 `distance2coast.mat`；引用 Greene et al. (2019), DOI `10.1029/2019GC008392`。

## 科学约定

主流程保存两类固定掩膜：

1. `altimetry_global` 是输入期间每个月都有有效 SLA 的标准全球测高域；
2. `budget_common` 在前者基础上保留距海岸至少 300 km 的海洋单元，测高共同域、GRACE 海洋质量和 OBD 使用完全相同的面积权重。

海岸单元不是简单的 0/1 判别，而是从 1/8° CDT 海陆掩膜聚合为 0.25° 海洋面积比例。任何月份如果有效固定权重不足 99.5%，该月记为缺测；分母不会随有效格点动态重归一化。

测高 GIA 标量校正为：

```text
GMSL_GIA(t) = GMSL_raw(t) + 0.30 mm/yr × (t - t0)
```

正号意味着校正后趋势比原始测高趋势高 0.30 mm/yr。代码会拒绝对已经标记为 GIA corrected 的序列重复校正。

GRACE 海洋质量与 OBD 共用一份预处理哈希：低阶项替换、减去 Caron 2018 GIA、去时间均值，并保留 ICGEM 已施加的 DDK1，不再叠加高斯滤波。EWH 转海平面当量使用

```text
rho_freshwater / rho_seawater = 1000 / 1028
```

OBD 采用向上为正；海底下沉为负。负荷 Love 数换算为 `R × h_l / (1 + k_l)`。质量和 OBD 每次只合成一个月并立即面积平均，不在内存中保存全球三维月场。

趋势模型使用实际日历时间，同时拟合截距、线性项、周年正余弦和半年正余弦：

```text
y(t) = b0 + b1(t-t_mean) + a1 sin(2πt) + c1 cos(2πt)
       + a2 sin(4πt) + c2 cos(4πt) + e(t)
```

输出 OLS 标准误和 HAC 标准误；HAC 默认最大滞后 12 个月。缺月不插值。

## 下载与运行

在仓库的 `reproduce_figure1` 目录执行：

```powershell
D:\python\python.exe -m pip install "h5netcdf>=1.3"
D:\python\python.exe scripts\download_icgem_csr.py --destination D:\temp_sealevel_data\grace_csr_ddk1 --start 2013-11 --end 2024-10
D:\python\python.exe run_gmsl_budget.py --config configs\gmsl_budget_main.json
```

下载器只使用 HTTPS 和系统 CA，验证 GFC 头并写 SHA-256 清单；已有有效文件可断点复用。若 ICGEM 返回 429 或出现临时连接中断，下载器进行有限退避重试，不会关闭证书校验。

测试与静态编译：

```powershell
D:\python\python.exe -m pytest tests\gmsl_budget -v
D:\python\python.exe -m compileall -q gmsl_budget run_gmsl_budget.py scripts
```

## 输出

每个配置写入 `output/optimized_budget/<run_id>/`，同一 `run_id` 不允许被不同配置覆盖：

- `monthly_budget.nc`、`monthly_budget.csv`：月序列；
- `trend_summary.csv`：趋势、OLS/HAC 不确定度与时间范围；
- `mask_altimetry_global.nc`、`mask_budget_common.nc`：坐标、海洋比例、固定支持及哈希；
- `config_resolved.json`、`provenance.json`：解析后的配置、输入 SHA-256、软件版本和处理哈希；
- `diagnostics.png`、`run_report.md`：诊断图与科学检查报告。

没有合格比容输入时，NetCDF 中不会出现 `steric` 或 `closure`，也不会用零值占位。

## 私有 GitHub 仓库

远端为 `https://github.com/Yos3mite/sea-level-change.git`。首次上传需先在用户自己的终端完成认证：

```powershell
git credential-manager configure
git ls-remote origin
```

浏览器登录后再推送功能分支。也可使用仅授权该仓库的 fine-grained PAT，权限只需 `Contents: Read and write`。推送前必须先读取远端历史；本流程不使用 force push。

## 当前限制

- 比容项尚未接入，因此当前结果是部分预算而不是闭合预算。
- 当前正式质量项为 CSR 球谐解；验证接口可对照外部 OBD 序列，但外部参考不会替换计算结果。
- 300 km 缓冲可降低陆地泄漏和海岸格点不一致，但不会消除所有 GRACE 泄漏误差；这些限制应随论文结果一并报告。
