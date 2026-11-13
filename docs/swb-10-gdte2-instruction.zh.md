# GdTe₂ 界面局域电子调制的自动发现与定量表征

## 任务背景

本题数据来自 GdTe₂ 表面的低温 STM/STS 实验。同一视场包含两个相邻、具有不同周期性电子调制的空间区域。你需要仅依据 `/app/data/` 中的数值数据，建立实空间与倒空间坐标，识别主要调制成分，自动完成区域分区与界面定位，并判断该界面附近是否出现具有能量选择性的额外局域电子调制。不得预设界面坐标、活跃能量或异常谱线编号。

学科标签：condensed matter physics、STM、STS、charge-density wave、Fourier analysis、scientific image analysis。

## 运行环境

- 所有输入数据位于 `/app/data/`（容器内只读）。
- 所有交付结果必须写入 `/app/output/`。
- 环境离线（no-network）：脚本不得访问网络，也不得读取 `/app/data/`、`/app/output/` 之外的题目文件。
- 已预装：Python 3.12、NumPy、SciPy、Matplotlib（含 Pillow 用于图像写出）。

## 坐标与数值约定

- 实空间坐标采用数据文件内的分析坐标：x 随列号增加，y 随存储行号增加，单位 nm。
- FFT 空间频率必须以 cycles/nm 报告；若同时给出 rad/nm，需注明换算 $q_{\mathrm{rad}}=2\pi q_{\mathrm{cycles}}$。
- 界面位置使用 $x_b(y)$ 表示；至少报告 $y=5.00$、$20.00$、$35.00$ nm 三处的 $x_b$。
- q 分量保留 3 位小数；长度保留 2 位小数；bias 保留 1 位小数；无量纲分数保留 2 位小数。
- 不得从文件名、固定像素坐标或人工硬编码 mask 推断答案。

## 任务目标

1. 检查 metadata、单位、像素尺度和 forward/retrace 配准，说明预处理步骤。
2. 计算形貌 FFT，区分一个跨两个区域存在的共同调制分量，以及分别偏向界面两侧的调制分量；给出至少 1 个共同向量、2 个左侧选择性向量和 2 个右侧选择性向量。
3. 通过局域 Fourier/lock-in 等等价的复振幅方法构造两侧调制强度，形成数据驱动的分区指标，并拟合界面曲线 $x_b(y)$。报告定位不确定度。
4. 对所有 bias map 使用同一界面和同一异常指标，自动判定哪些已采样 bias 存在界面局域增强。不要在未采样能量之间外推连续端点。
5. 对检测到的界面信号给出空间 FWHM 或等价宽度；明确宽度定义。
6. 将弯曲界面拉直到界面坐标 $(u,y)$，比较活跃 bias 与对照 bias 的 Fourier 功率，报告主要界面增强特征的 $(q_u,q_y)$、模长和方向。
7. 在 P01/P02/P03 中自动判断哪条谱具有区别于其余两条平均谱的正偏压局域特征，并用固定、可复现的基线和半高宽定义报告峰位、左右边缘和宽度。
8. 将相同判据分别应用于 forward 与 retrace，验证活跃 bias 集合是否一致，并比较两次界面定位。

## 要求输入（位于 `/app/data/`）

- `metadata.json`：测量条件、单位、坐标约定和文件索引。
- `region1_topography.npz`：40 nm × 40 nm 的形貌数据，包含 `forward`、`retrace`、`x_nm`、`y_nm`。
- `map_E001.npz` 至 `map_E004.npz`：同一视场在四个 bias 下的空间微分电导代理信号；每个文件均含 `forward`、`retrace`、`x_nm`、`y_nm`。实际 bias 以 `metadata.json` 为准（依次为 +100.0、+300.0、+400.0、+900.0 mV）。
- `spectra_set_A.csv`：三个中性编号位置 P01、P02、P03 的点谱。包含 bias、同相微分电导代理、正交分量和隧穿电流。

数组单位、像素数、扫描尺寸、setpoint、调制幅度与频率均在 `metadata.json` 中给出。`retrace` 已仅做无损的 x 方向反转，使其列坐标与 `forward` 一致；仍需自行检查并校正残余平移。

## 要求输出（全部写入 `/app/output/`）

- `analysis_report.md`：完整方法、公式、参数、误差来源、物理结论和局限性。
- `summary.json`：机器可评分结果，严格使用下方字段。
- `reciprocal_space.png`
- `domain_segmentation.png`
- `bias_dependence.png`
- `boundary_feature_fft.png`
- `spectroscopy_comparison.png`
- `validation_retrace.png`
- `analysis.py`：从 `/app/data/` 重新生成全部结果的可复现脚本；不得访问网络或 `/app/data/`、`/app/output/` 之外的文件。

`summary.json` 结构：

```json
{
  "schema_version": "1.0",
  "reciprocal_convention": "cycles/nm",
  "common_modulation_vector": {"qx": null, "qy": null, "magnitude": null},
  "left_domain_vectors": [{"qx": null, "qy": null, "magnitude": null}],
  "right_domain_vectors": [{"qx": null, "qy": null, "magnitude": null}],
  "boundary": {
    "sample_y_nm": [5.0, 20.0, 35.0],
    "sample_x_nm": [null, null, null],
    "uncertainty_nm": null
  },
  "boundary_anomaly": {
    "detected": null,
    "active_sampled_biases_mV": [],
    "per_bias_boundary_excess_z": {},
    "active_spatial_fwhm_nm": {}
  },
  "fourier_feature": {
    "q_u_cycles_per_nm": null,
    "q_y_cycles_per_nm": null,
    "magnitude_cycles_per_nm": null,
    "angle_deg_from_u": null
  },
  "spectroscopy": {
    "anomalous_point_id": null,
    "feature_peak_mV": null,
    "left_edge_mV": null,
    "right_edge_mV": null,
    "width_mV": null
  },
  "retrace_validation": {
    "forward_active_biases_mV": [],
    "retrace_active_biases_mV": [],
    "boundary_forward_retrace_rmse_nm": null
  }
}
```

`per_bias_boundary_excess_z` 与 `active_spatial_fwhm_nm` 以 bias（mV，1 位小数）的字符串为键，例如 `"300.0"`。

允许使用与参考流程不同但科学等价的方法；必须在报告中给出定义，使结果可复算。所有判断需同时有数值结果和图像证据。评分只检查 `/app/output/` 中的最终文件、图和 `summary.json`，不依据所用软件、库或代码风格给分。
