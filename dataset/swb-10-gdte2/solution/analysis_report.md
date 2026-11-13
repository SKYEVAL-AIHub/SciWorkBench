# 标准答案分析报告

## 1. 结论摘要

提供的 40 nm × 40 nm 视场中存在两个具有不同低 q 调制族的相邻区域。由局域 Fourier 振幅比自动得到的界面从左下向右上延伸。在四个已采样 bias 中，仅 +300 mV 与 +400 mV 同时在 forward 和 retrace 中超过统一的界面增强阈值；+100 mV 与 +900 mV 为对照。点谱 P02 相对 P01/P03 平均谱还存在独立的正偏压局域谱学特征。

本报告中的数值均由 `data/` 重算。论文与 Supporting Information 只用于数量级和数据映射 sanity check，未用于替代 Ground Truth。

## 2. 数据与预处理

- 形貌和谱图：512 × 512 像素，40 nm × 40 nm，像素尺度 0.078125 nm。
- `retrace` 已在导出阶段做无损 x 反转；参考流程再用 phase correlation 的整数像素峰完成残余配准。
- forward/retrace 配准后，四张 map 的平移量依次为 $(1,0)$、$(0,9)$、$(0,9)$、$(1,0)$ 像素，顺序为 $(\Delta y,\Delta x)$。
- 形貌和 map 逐行减去一阶背景，不使用高阶背景、人工界面 mask 或出版图上的边界。
- 所有空间频率使用 cycles/nm：

$$
q_x=\frac{k_x}{N_x\Delta x},\qquad q_y=\frac{k_y}{N_y\Delta y}.
$$

若使用 rad/nm，应整体乘以 $2\pi$。

## 3. 倒空间调制向量

形貌使用二维 Hann window 后做 FFT。为减少界面两侧信号混合，共同分量在全视场拟合；侧选择性分量分别在左、右子区域的 FFT 中，以局部峰强度加权质心精化。结果如下。

| 类型 | $q_x$ (cycles/nm) | $q_y$ (cycles/nm) | $|q|$ (cycles/nm) |
|---|---:|---:|---:|
| 共同分量 C | 0.655 | 0.275 | 0.710 |
| 左侧选择性 L1 | 0.686 | -0.828 | 1.075 |
| 左侧选择性 L2 | 0.221 | 1.007 | 1.031 |
| 右侧选择性 R1 | 0.141 | 0.490 | 0.510 |
| 右侧选择性 R2 | 0.422 | -0.433 | 0.604 |
| 右侧选择性 R3 | 0.286 | 0.981 | 1.022 |

全视场 FFT 间隔为 0.025 cycles/nm；分区子窗口的 $q_x$ 间隔约 0.058-0.060 cycles/nm。因此 rubric 对每个分量使用 ±0.10 cycles/nm，而不按脚本输出的末位小数评分。L/R 标签仅指数据驱动的空间选择性，不要求考生预先赋予论文中的 CDW-II/CDW-III 名称。

## 4. 局域振幅、分区与界面

对每个目标向量 $\mathbf q$，参考流程计算

$$
A_{\mathbf q}(\mathbf r)=G_{\sigma}\ast\left[z(\mathbf r)e^{-i2\pi\mathbf q\cdot\mathbf r}\right],
$$

其中 $G_{\sigma}$ 为实空间 Gaussian，基准 $\sigma=1.7$ nm。对同一侧的分量先除以各自中位数再平均，得到 $A_L$ 与 $A_R$，并构造

$$
\eta(\mathbf r)=\ln\frac{A_R(\mathbf r)+\epsilon}{A_L(\mathbf r)+\epsilon}.
$$

逐行搜索平滑 $\eta$ 的最大正梯度，再以二次曲线稳健拟合：

$$
x_b(y)=ay^2+by+c,
$$

其中

$$
(a,b,c)=(-8.8315\times10^{-4},\ 0.410507,\ 16.266345).
$$

| $y$ (nm) | $x_b$ (nm) |
|---:|---:|
| 5.00 | 18.30 |
| 20.00 | 24.12 |
| 35.00 | 29.55 |

行级梯度点相对拟合的稳健不确定度为 0.70 nm。改变局域幅度和平滑参数后，上述三个位置的范围分别为 17.98-18.68、23.91-24.27、29.43-29.55 nm；forward/retrace 独立界面曲线的 RMSE 为 1.68 nm。综合二者，评分容差取 ±2.00 nm。

## 5. Bias-dependent 界面异常

对每张 map 先减逐行一阶背景，再减 1.8 nm Gaussian 低频分量。对残差平方以 0.7 nm Gaussian 求局域 RMS $A(E,\mathbf r)$。界面带定义为 $|d|\leq3$ nm，bulk 对照定义为 $8\leq|d|\leq15$ nm。统一异常分数为

$$
S(E)=\frac{\operatorname{median}(A)_{|d|\leq3}
-\operatorname{median}(A)_{8\leq|d|\leq15}}
{1.4826\operatorname{MAD}(A)_{8\leq|d|\leq15}}.
$$

判定阈值固定为 $S>1.0$。

| ID | Bias (mV) | 界面中位 RMS (pA) | Bulk 中位 RMS (pA) | $S(E)$ | 检测 | 中心峰 FWHM (nm) |
|---|---:|---:|---:|---:|---|---:|
| E001 | +100 | 0.106 | 0.109 | -0.23 | 否 | 不评分 |
| E002 | +300 | 0.187 | 0.112 | 3.14 | 是 | 8.50 |
| E003 | +400 | 0.119 | 0.093 | 1.75 | 是 | 9.50 |
| E004 | +900 | 0.238 | 0.232 | 0.13 | 否 | 不评分 |

因此 Ground Truth 是“已采样活跃 bias 集合为 {+300, +400} mV”。原始包没有 +200、+500 或 +600 mV map，故不能把论文中的连续能量窗口当作本题数据的精确端点；rubric 只评分已采样 bias 的集合。

空间 FWHM 定义为：对局域 RMS 沿界面平行方向取中位数，得到随 signed distance 的剖面；以 $|d|\geq8$ nm 的中位数为基线，对中心峰使用半高连续区间。参数变化给出的主要范围支持 ±3.00 nm 容差。

## 6. 界面 Fourier 特征

逐行将界面平移到 $u=0$，使用 $-8\leq u\leq8$ nm、$2.5\leq y\leq37.5$ nm 的拉直条带。每个 bias 的 Fourier 幅度用 $0.15\leq|q|\leq1.25$ cycles/nm 环带中位数归一化。定义

$$
R(q_u,q_y)=\ln\frac{\langle P\rangle_{300,400}}
{\langle P\rangle_{100,900}}.
$$

正 $q_y$ 半平面的最强稳定增强为：

- $q_u=0.621$ cycles/nm；
- $q_y=0.513$ cycles/nm；
- $|q|=0.806$ cycles/nm；
- 相对 $+u$ 方向角度 39.55°；
- 局部二阶矩估计的 Fourier FWHM 主/次轴约 0.213/0.149 cycles/nm。

将条带半宽改为 7、8、9 nm 时，$q_u$ 为 0.621-0.665 cycles/nm，$q_y$ 不变。因此分量容差取 ±0.12 cycles/nm，角度允许 ±15°。考生若使用对称的 $(-q_u,-q_y)$ 峰，按同一物理结果评分。

## 7. 点谱证据

对 P01、P02、P03 逐一执行“候选谱减去其余两谱平均”，并在 +80 至 +360 mV 内寻找最大正偏压特征。差谱使用 21 点、三阶 Savitzky-Golay 平滑；基线为 0-80 mV 与 380-500 mV 侧带上的一次拟合。P02 的正偏压峰最高，因此自动选为异常位置。

基线校正后的 P02 差谱结果：

- 峰位：+150.0 mV；
- 左边缘：+127.5 mV；
- 右边缘：+202.5 mV；
- FWHM：75.0 mV；
- 半高区间中点：+165.0 mV；
- 能量步长：2.5 mV；测量调制幅度：4 mV peak。

把平滑窗口改为 17、21、25、31 点时，异常位置始终为 P02，峰位为 147.5-155.0 mV，左右边缘为 125.0-127.5 与 202.5 mV。结合能量步长和调制展宽，采用峰位/边缘 ±15 mV、宽度 ±25 mV 的容差。

该 raw-derived 数值与论文约 +240 mV 中心、约 160 mV 宽度的出版值不一致。根据任务优先级，Ground Truth 采用本包内三条原始点谱与上述固定定义的重算结果，不为贴合论文而调参。可能原因包括点谱批次、位置与出版图所用平均策略不同。

## 8. Forward/retrace 独立验证

在不改变阈值和界面带定义的前提下：

- forward 检测到 {+300, +400} mV；
- retrace 检测到 {+300, +400} mV；
- +300 mV 的 $S$ 为 3.50/3.40；
- +400 mV 的 $S$ 为 1.74/1.71；
- +100 与 +900 mV 在两方向均低于 1.0；
- forward/retrace 界面 RMSE 为 1.68 nm。

原始材料中没有一组能够由 metadata 独立确认且覆盖正偏压活跃窗口的第二条畴壁数据，因此本题没有伪造“第二畴壁复现”，而以独立扫描方向作为验证集。

## 9. 鲁棒性与可复现性

参考脚本固定所有随机性（实际主流程无随机步骤），并测试：

- 局域幅度 Gaussian：1.4、1.7、2.0 nm；
- 界面梯度平滑：0.5、0.7、0.9 nm；
- map 高通尺度：1.5、1.8、2.1 nm；
- 局域 RMS 尺度：0.6、0.7、0.8 nm；
- Fourier 条带半宽：7、8、9 nm；
- 谱线平滑窗口：17、21、25、31 点。

9 组空间参数组合均得到活跃 bias {+300, +400} mV；所有谱线参数均选择 P02。完整试验值见 `robustness_summary.json`，可机读结果见 `reference_results.json` 和 `reference_results.csv`。

## 10. 论文/SI sanity check 与符号核验

论文与 SI 报告界面调制在 +300/+400 mV 明显、在高能对照消失，这与本数据的离散 bias 判断一致。SI S12 caption 将高能图写为 +900 mV，而主文一处文字写成 -900 mV；本包对应 SXM header 明确为 `BIAS = +0.900 V`，因此数据导出与评分均使用 +900 mV。

论文题名：*Emergent One-Dimensional Charge Order at Domain Walls of Competing Charge-Density-Wave Orders*，DOI: 10.1021/acs.nanolett.6c00774。
