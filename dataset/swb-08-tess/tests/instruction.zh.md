# AST-TESS-001：活跃恒星周期性暗化事件的耀斑感知验证

## 任务背景

> 目录约定：公开输入固定在 `data/`（容器内为 `/workspace/data`），与平台可见模式 `data/**` 完全一致。所有科学结果必须从该目录独立计算；不得联网获取额外答案，也不得改写输入材料。所有提交只写入 `output/`（容器内为 `/workspace/output`）。

你获得的是 AU Mic 的真实 TESS Sector 1 与 Sector 27 观测派生数据（NASA/MAST 公共产品）。任务要求完成一条可复现的观测天文分析链，并判断周期性暗化信号属于：

- `robust`：周期、逐事件计时、源位置和恢复完备性共同支持；
- `ambiguous`：存在信号，但关键证据相互冲突或不足；
- `artifact`：活动、系统误差或邻近源更能解释信号。

不得以抄录论文数值、只做常规去趋势 + BLS、或只写叙述报告代替分析。所有结果必须从 `data/` 重算，并在 `replay_manifest.json` 中记录命令、环境、随机种子、耗时与输入/输出哈希。

## 必做科研链

1. **时间与质量**：统一 `BTJD = BJD_TDB - 2457000`，保留原始 `QUALITY`，说明拒绝位；分别处理 20 s 与 120 s 曝光，输出模型必须做曝光积分。
2. **活动分离**：利用 32 个 20 s 耀斑窗口建立耀斑判据；至少给出两种不抹平事件的去趋势结果，并说明星斑、耀斑和暗化如何区分。
3. **周期搜索**：执行 BLS 或等价开源方法，同时量化窗口函数；必须比较 `P/2`、`P`、`2P`，按覆盖事件、支持事件和矛盾事件裁决别名。
4. **逐事件计时**：对 5 个事件分别拟合中点、深度、持续时间和误差；至少 300 次 bootstrap 或等价重采样；拟合线性星历、O−C，并做逐事件 leave-one-out。
5. **源位置**：从 5 个差分像事件独立测质心，单位同时为 pixel 和 arcsec；结合 ICRS 邻星目录和 Tmag 通量比估计孔径稀释/背景污染。
6. **恢复与敏感性**：完整运行 `data/config/injection_grid.csv` 的 120 组注入，逐行记录恢复；完整运行 8 个敏感性场景。最终分类必须由前五步证据驱动。

## 天文约定

- 时间列 `time_btjd` 为 TDB 标度；输出中点同时给 BTJD 与 BJD_TDB。
- 天球坐标为 ICRS，角度为 degree；质心偏移为 pixel 与 arcsec，TESS 像素尺度 21 arcsec pixel⁻¹。
- 相对通量无量纲，基线约为 1；原始 flux 列单位为 electron s⁻¹。
- 20 s 与 120 s cadence 不能视为瞬时采样；模型必须进行曝光积分。
- 不确定度、随机种子、阈值、掩膜和异常值处理必须显式记录。

## 要求输出

所有提交仅写入 `output/`，不得修改 `data/`。必须包含且仅包含以下十个文件：

1. `results.json`：最终分类、主参数、六级核心证据（`core_chain` 的 R12–R17）和交叉引用；
2. `event_metrics.csv`：逐事件中点、误差、深度、持续时间、O−C、接受状态；
3. `period_search.csv`：P/2、P、2P 的 BLS、窗口函数、支持/矛盾事件和选择标记；
4. `centroid_audit.csv`：逐事件差分像质心、偏移、S/N、稀释和 on-target 状态；
5. `injection_recovery.csv`：120 个注入逐行输入与恢复结果；
6. `sensitivity.csv`：8 个场景逐行参数、周期/别名和分类稳定性；
7. `analysis_report.md`：中文科学报告，结论先行，含方法、证据、限制与来源；
8. `diagnostic_summary.png`：至少包含时序、别名、O−C、质心、注入恢复、敏感性六个可读面板；
9. `replay_manifest.json`：命令、Python/库版本、随机种子、耗时、输入与输出 SHA-256；
10. `analysis_pipeline.py`：从 `data/` 生成其余九个文件的可运行入口。

### 运行约束

- 普通个人电脑、CPU 路径；总运行时间不得超过 120 分钟。
- 运行时不得依赖网络、私有路径、闭源软件、商业许可证或复杂编译。
- 固定随机种子（`20260819`）；同一环境重复运行的结构化数值结果应在声明容差内一致。
- 禁止更改输入数据、伪造来源或将论文中的答案直接写入输出。
- 推荐使用 NumPy、SciPy、pandas、Astropy 和 matplotlib；均为开源、可通过 pip 安装。

### 结构化输出合同

下列 JSON 是本题完整的机器可读输出合同，只规定文件名、字段、行数和枚举，不含任何 Golden 数值。

```json
{
  "schema_version": "1.0.0",
  "task_id": "AST-TESS-001",
  "required_files": [
    "analysis_pipeline.py",
    "analysis_report.md",
    "centroid_audit.csv",
    "diagnostic_summary.png",
    "event_metrics.csv",
    "injection_recovery.csv",
    "period_search.csv",
    "replay_manifest.json",
    "results.json",
    "sensitivity.csv"
  ],
  "results_required_keys": [
    "task_id",
    "random_seed",
    "time_system",
    "selected_period_days",
    "classification",
    "event_count",
    "accepted_event_count",
    "centroid_event_count",
    "injection_count",
    "injection_recovery_fraction",
    "catalog_dilution_fraction_42arcsec",
    "linear_ephemeris",
    "leave_one_out_periods",
    "core_chain",
    "runtime_seconds"
  ],
  "classification_enum": ["robust", "ambiguous", "artifact"],
  "csv_required_columns": {
    "event_metrics.csv": ["event_id", "epoch", "sector", "midtime_btjd", "midtime_bjd_tdb", "midtime_error_days", "depth", "depth_error", "duration_days", "exposure_seconds", "oc_seconds", "accepted"],
    "period_search.csv": ["alias_label", "period_days", "bls_snr", "bls_t0_btjd", "bls_duration_days", "bls_depth", "window_function_power", "covered_events", "supported_events", "contradictions", "event_consistency", "adjudication_score", "selected"],
    "centroid_audit.csv": ["event_id", "difference_x_pix", "difference_y_pix", "reference_x_pix", "reference_y_pix", "offset_pix", "offset_arcsec", "difference_snr", "catalog_dilution_fraction_42arcsec", "on_target"],
    "injection_recovery.csv": ["injection_id", "depth", "duration_days", "period_days", "phase_offset_days", "recovered_period_days", "recovered_t0_btjd", "recovered_depth", "detection_snr", "phase_error_days", "recovered"],
    "sensitivity.csv": ["scenario_id", "quality_bits", "flare_sigma", "detrend_days", "exposure_integrated", "selected_period_days", "selected_alias", "supported_events", "classification_stable"]
  },
  "row_count_contract": {
    "event_metrics.csv": 5,
    "period_search.csv": 3,
    "centroid_audit.csv": 5,
    "injection_recovery.csv": 120,
    "sensitivity.csv": 8
  },
  "diagnostic_png_minimum": {"width_px": 1200, "height_px": 700}
}
```

## 数据来源与许可

- 观测数据：NASA/MAST 的 TESS SPOC light-curve 与 target-pixel 公共产品（[MAST TESS Archive](https://archive.stsci.edu/missions-and-data/tess)）。
- 科学背景：Gilbert et al. (2022), *Flares, Rotation, and Planets of the AU Mic System from TESS Observations*, [arXiv:2109.03924](https://arxiv.org/abs/2109.03924), [DOI:10.3847/1538-3881/ac23ca](https://doi.org/10.3847/1538-3881/ac23ca)。题包不打包论文全文，所有任务数值均由公开 `data/` 独立重算。
- 报告 `analysis_report.md` 必须明确区分数据来源、科学背景与本次独立分析，并引用 MAST/TESS 与 Gilbert et al. 的 arXiv 或 DOI。
