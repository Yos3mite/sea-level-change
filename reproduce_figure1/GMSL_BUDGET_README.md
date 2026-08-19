# 可审计的全球平均海平面主流程

这套流程是新增的干净主流程，不调用、覆盖或修改 `reproduce_figure1` 中的历史复现脚本，也不依赖 L2toL3 技能。正式结果采用标准全球海洋域，不输出 `66°S–66°N` 口径。当前版本处理测高 GMSL、CSR/JPL/GSFC Mascon 海洋质量和 OBD；比容项暂缓，因此不生成预算闭合残差。

## 数据与口径

- 测高：`SEALEVEL_GLO_PHY_CLIMATE_L4_MY_008_057` 对应的 C3S/CMEMS 双卫星月平均 SLA，本地文件为 `D:/temp_sealevel_data/cmems_sla_monthly.nc`，时间为 2013-11 至 2024-10。流程直接读取月度产品，不把日数据重新聚合成月数据。
- 测高质量控制：可与 `D:/temp_sealevel_data/cmems_gmsl_indicator.nc` 中的官方面积平均指标比较。官方指标已经包含 GIA 校正、季节调整和六个月滤波，不可当作原始 SLA 面积平均后再次加 GIA。
- Mascon 海洋质量：CSR、JPL 和 GSFC 三个官方产品；三者产品均已去除 ICE6G-D GIA，流程禁止再次施加 GIA。各中心分别进行时间序列 SSA 缺月填补，之后才作等权平均。
- OBD：ICGEM 官方 CSR RL06/RL06.3 `BA01`（60 阶）`DDK1` 球谐文件；该球谐质量序列仅作为 OBD 一致性诊断，不替代三中心 Mascon 主质量项。
- 低阶项与 GIA：TN-11、TN-13 CSR RL06.2、TN-14，以及 Caron et al. (2018) GIA，均取自本地 SaGEA 数据目录。
- 海陆与海岸距离：Climate Data Toolbox 的 `land_mask.mat` 和 `distance2coast.mat`；引用 Greene et al. (2019), DOI `10.1029/2019GC008392`。

## 科学约定

主流程保存四类固定掩膜：

1. `altimetry_global` 是输入期间每个月都有有效 SLA 的标准全球测高域；
2. `budget_common` 在前者基础上保留距海岸至少 300 km 的海洋单元，测高共同域和 OBD 使用相同面积权重；
3. `mascon_global_1deg` 是 Jin 质量曲线的主海洋域，不设置海岸缓冲；
4. `mascon_300km_sensitivity` 只用于检验海岸缓冲敏感性，不作为主结果。

海岸单元不是简单的 0/1 判别，而是从 1/8° CDT 海陆掩膜聚合为 0.25° 海洋面积比例。任何月份如果有效固定权重不足 99.5%，该月记为缺测；分母不会随有效格点动态重归一化。

测高从 2016-01 起先施加湿对流层漂移校正 `-0.50 mm/yr`，再施加 GIA 标量校正：

```text
GMSL_GIA(t) = GMSL_wet(t) - 0.30 mm/yr × (t - t0)
```

这里配置中的 GIA 是带符号加数 `-0.30 mm/yr`，因此校正后趋势降低 0.30 mm/yr。代码会拒绝对已经标记为 GIA corrected 的序列重复校正。

OBD 球谐分支共用一份预处理哈希：低阶项替换、减去 Caron 2018 GIA、去时间均值，并保留 ICGEM 已施加的 DDK1，不再叠加高斯滤波。Mascon 分支保留产品自带的 ICE6G-D GIA 去除。EWH 转海平面当量使用

```text
rho_freshwater / rho_seawater = 1000 / 1028
```

OBD 采用向上为正；海底下沉为负。负荷 Love 数换算为 `R × h_l / (1 + k_l)`。质量和 OBD 每次只合成一个月并立即面积平均，不在内存中保存全球三维月场。

三中心时间序列 SSA 的默认窗口为 36 个月、秩为 8、收敛阈值为 `1e-5`。填补只修改缺月，观测月保持原值；这是无法下载 Xie & Yi 已发布格网数据时的可审计回退方案，并不冒充其原始格网级 SSA 产品。

趋势模型使用实际日历时间，同时拟合截距、线性项、周年正余弦和半年正余弦：

```text
y(t) = b0 + b1(t-t_mean) + a1 sin(2πt) + c1 cos(2πt)
       + a2 sin(4πt) + c2 cos(4πt) + e(t)
```

输出 OLS 标准误和 HAC 标准误；HAC 默认最大滞后 12 个月。缺月不插值。

## 下载与运行

在仓库的 `reproduce_figure1` 目录执行：

```powershell
D:\python\python.exe -m pip install "h5netcdf>=1.3" "netCDF4>=1.6"
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
- `trend_summary.csv`：估计趋势、OLS/HAC 不确定度，并单列最终采用的 GMSL 趋势 `3.057 mm/yr`；
- `final_summary.json`：采用值、实际 GMSL+OBD 分量和 Jin 海洋质量参考值之间的差异；
- `mask_*.nc`：四类掩膜的坐标、海洋比例、固定支持及哈希；
- `config_resolved.json`、`provenance.json`：解析后的配置、输入 SHA-256、软件版本和处理哈希；
- `figure1_mass_component.png`：三中心曲线、中心平均和去趋势三个月滑动平均质量曲线；
- `diagnostics.png`、`run_report.md`：诊断图与科学检查报告。

没有合格比容输入时，NetCDF 中不会出现 `steric` 或 `closure`，也不会用零值占位。
输出固定为 NetCDF3/SciPy 后端，以兼容 Windows 中文目录。若环境同时安装了 `netCDF4`，使用 Xarray 打开这些文件时应显式写 `xr.open_dataset(path, engine="scipy")`，避免 `netCDF4` C 库在中文路径上失败。

## 私有 GitHub 仓库

远端为 `https://github.com/Yos3mite/sea-level-change.git`。首次上传需先在用户自己的终端完成认证：

```powershell
git credential-manager configure
git ls-remote origin
```

浏览器登录后再推送功能分支。也可使用仅授权该仓库的 fine-grained PAT，权限只需 `Contents: Read and write`。推送前必须先读取远端历史；本流程不使用 force push。

## 当前限制

- 比容项尚未接入，因此当前结果是部分预算而不是闭合预算。
- 当前下载的 CSR/JPL 产品版本及覆盖期并非 Jin/Xie & Yi 处理时所用产品的逐字节副本；时间序列 SSA 也不是其已发布的格网级补全数据。因此趋势差异必须保留并报告，不能通过调斜率强制匹配。
- 300 km 缓冲可降低陆地泄漏和海岸格点不一致，但本流程实测它会改变三中心趋势，因此只作为敏感性结果。
