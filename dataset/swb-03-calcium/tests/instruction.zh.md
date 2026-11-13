# 活体神经元钙成像视频定量

## 研究背景

某类外周感觉神经元在短时化学刺激前后进行了活体钙成像。活体记录不可避免地带有缓慢视野位移、荧光强度随记录时间衰减以及少量采集异常帧，这些技术性变化会与真实的钙信号混在一起。此外，视频播放帧率与实际采集时间往往并不一致，因此时间轴必须来自实测时间戳。要判断某个细胞是否真的对刺激产生了响应，必须先把运动、背景、周围组织信号和强度漂移逐一扣除。

刺激从 5.5 s 开始，在 6.0 s 结束。记录包含连续荧光视频、每帧实际采集时间、一个无运动参考图像，以及 6 个细胞区域、背景区域和稳定强度参考区域。

## 研究问题

恢复各帧在参考坐标中的位置，识别不可靠帧，并在校正背景、周围组织信号和记录期间强度漂移后，定量比较 6 个细胞的刺激响应。分析必须使用 `timestamps.csv` 中的实际时间，不能用视频播放帧率代替。

## 任务输入

- `data/recording.avi`：连续荧光视频。
- `data/timestamps.csv`：每帧实际采集时间。
- `data/reference_frame.tif`：无运动参考图像。
- `data/regions/`：6 个细胞区域、背景区域和稳定强度参考区域；区域文件与标签对应关系见 `data/regions/region_labels.csv`。
- `data/roi_map.png`：区域位置示意。
- `data/recording_notes.json`：实验室采用的定量口径。该文件中的时间窗、区域统计口径、帧质控阈值和响应判定规则均属于本题数据定义。

## 任务输出

所有结果保存到 `output/`，只提交以下四个文件。

`motion_qc.csv`：列为 `frame,time_s,shift_x_px,shift_y_px,global_mean_ratio,saturated_pixel_fraction,valid`，共 240 行。位移表示将该视频帧校正到参考图像时所施加的 x、y 平移量，`valid` 用 `yes` 或 `no`。

`calcium_traces.csv`：列为 `frame,time_s,cell1,cell2,cell3,cell4,cell5,cell6`，共 240 行。细胞列为校正后的相对荧光响应，单位为百分数。

`summary.csv`：列为 `cell,F0,baseline_sd_percent,peak_dFF0_percent,peak_time_s,latency_half_peak_s,AUC_percent_s,responder`，每个细胞一行，`responder` 用 `yes` 或 `no`。

`analysis.png`：同一张图中给出参考图像及 6 个细胞区域、x/y 位移与被拒绝帧、6 条细胞响应曲线，并标出刺激时间段。
