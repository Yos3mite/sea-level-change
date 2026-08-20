# Figure 1 自制 Level-3 方法报告

- 时间范围：2013-11—2024-10。
- 自制 L3 观测月：108/132；MSSA 重建月：24；3 月中心平滑可绘月：132。
- 自制 L3 来源方法：untouched pristine SaGEA control using local CSR DDK1 Level-2, explicit low-degree replacements where available, DDK1-matched Caron et al. (2018) GIA, and the 300 km open-ocean mask; 24 missing Level-2 months reconstructed by direct Improved MSSA with M=60 and rank=14。
- 缺测处理：上游在 Level-2 球谐系数层用 Improved MSSA 重建 24 个月；本绘图不再进行线性插值；重建月不是卫星观测。
- 测高 GMSL GIA 校正：-0.414000 mm/yr；来源：project-registered correction for the fixed local GMSL reference。
- 测高 GMSL 湿对流层校正：-0.500000 mm/yr，自 2016-01 起。
- 显示偏移（仅绘图，不进入指标）：panel (a) {'gmsl': 0.0, 'custom_l3': -10.0, 'steric': -20.0, 'budget_sum': 0.0}；panel (b) {'gmsl': 0.0, 'custom_l3': -5.0, 'steric': -10.0, 'budget_sum': 0.0}。
- 预算恒等式：`custom mass + steric`；OBD 不加入该曲线。
- Forward Modeling 只在叠加图中作为敏感性；物理失稳的 iterative-ocean 方法明确排除。
- 自制 L3 输入原始趋势：1.447873 mm/yr。
- 自制 L3 线性趋势：1.387462 mm/yr；Forward 敏感性：0.141119 mm/yr。
- Forward 相对自制 L3：RMSE 4.191901 mm，相关系数 0.727680。
- CSR mascon 参考线性趋势：1.231482 mm/yr；相对自制 L3 的 RMSE 3.457824 mm，相关系数 0.829136。
- GMSL 趋势：3.213646 mm/yr；自制质量 + steric 趋势：2.875240 mm/yr。
- 自制 L3 与 CSR mascon 的差异仍需结合数据版本、滤波、先验约束和泄漏处理解释，不能视为 CSR mascon 参考的等价替代。
