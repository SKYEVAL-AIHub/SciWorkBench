# Fe-doped monolayer MoS2 magneto-Raman 参考答案

本文件是批卷用的标准答案与判据说明。`answer/` 不提供参考图片；考生仍须按照 `agent_readme.md` 在 `output/` 中生成 `raman_analysis.png`。图片按是否正确呈现关键结果评分，不要求与某一参考图具有相同配色、版式或绘图风格。

## 1. 样品对照与谱线模型

应首先检查约 200–320 cm⁻¹ 的中低频 helicity-switched Raman 区域，而不是只分析常规高频声子。未掺杂 MoS2 在约 267 cm⁻¹ 表现为近似单一的 orbital–phonon hybridized Raman 模；Fe 掺杂样品在同一区域出现明显非对称线形，需要至少两个相互重叠的谱学分量描述。

公开数据对应论文采用 Gaussian profile 分解得到：

| 样品 | 分量 | Raman shift (cm⁻¹) |
|---|---|---:|
| MoS2 | 原有 hybridized mode | 约 267 |
| Fe:MoS2 | 原有 hybridized mode | **264.92 ± 0.25** |
| Fe:MoS2 | 掺杂诱导的额外 mode | **282.86 ± 0.19** |

因此，核心样品对照结论是：Fe 掺杂后在约 283 cm⁻¹ 出现可重复分辨的额外 Raman 分量，同时原有约 267 cm⁻¹ hybridized mode 发生约 2 cm⁻¹ 的位移。

### 模型选择参考

满分答案必须给出残差、AIC/BIC 或等价模型选择证据。一个独立的参考实现使用 230–310 cm⁻¹、Gaussian 分量加线性背景时，Fe:MoS2 的单分量与双分量模型 AIC 分别约为 2133.4 和 2092.7，即 \(\Delta\mathrm{AIC}\approx40.7\)，明显支持双分量描述。该数值依赖拟合窗口、峰型和背景，不作为硬性目标；评分检查的是是否存在可信的模型选择证据，而不是是否复现这一 AIC 数值。

允许不同合理背景模型、拟合区间和 Gaussian/Lorentzian/Voigt 峰型产生小幅数值差异，但应稳定恢复约 265 cm⁻¹ 与约 283 cm⁻¹ 两个 Fe:MoS2 分量。若参数大量卡在边界或残差呈系统结构，则不能视为成功分峰。

## 2. 磁场与 helicity 依赖

磁场分析必须连续追踪上述两套中低频 hybridized branches。不能在每个磁场点独立拟合后仅按“高频/低频”重新排序，因为这种排序不能保证跨磁场保持同一物理 branch。

正确的场依赖应同时满足：

- 零场 RL 与 LR 对应 branch 在实验不确定度内近似简并；
- 改变磁场方向后，RL 与 LR 的频移方向交换；
- 约 265 cm⁻¹ 和约 283 cm⁻¹ 两套分量都表现出近似线性的场致频移；
- 频率变化是 helicity 相关的，而不是只有强度变化或共同漂移；
- 逐场拟合没有明显 branch switching、参数碰边界或大误差离群点被无条件纳入线性回归。

公开论文对同一数据的 Gaussian branch fit 给出：

| 分量 | helicity | \(\mu_{eff}\) (\(\mu_B\)) | 等价单支斜率 (cm⁻¹/T) |
|---|---|---:|---:|
| 原有 hybridized mode | RL | **2.14 ± 0.07** | **+0.999 ± 0.033** |
| 原有 hybridized mode | LR | **2.29 ± 0.07** | **−1.069 ± 0.033** |
| 掺杂诱导 mode | RL | **2.82 ± 0.09** | **+1.317 ± 0.042** |
| 掺杂诱导 mode | LR | **2.76 ± 0.06** | **−1.289 ± 0.028** |

这里采用单个 helicity branch 的定义

\[
hc\,\Delta\tilde{\nu}=\mu_{eff}B.
\]

若由两 helicity branch 的频差求磁矩，则应使用

\[
hc(\tilde{\nu}_1-\tilde{\nu}_2)=2\mu_{eff}B.
\]

两种定义必须给出一致的 \(\mu_{eff}\)。不能把双支 splitting slope 直接按单支 slope 换算，否则会产生 factor-of-two 错误。

### 数值接受范围与线性拟合质量

公开 CSV 为论文作图数据，独立重拟合时不同背景和约束会带来一定差异，因此评分采用范围而不是单值：

- 原有约 265 cm⁻¹ 分量：\(\mu_{eff}\) 约 **2.0–2.5 \(\mu_B\)**；
- 约 283 cm⁻¹ 新分量：\(\mu_{eff}\) 约 **2.5–3.1 \(\mu_B\)**；
- RL/LR 的斜率方向应相反，并满足磁场反向后的 helicity 交换关系；
- 线性拟合应给出斜率不确定度、拟合优度或残差。由于拟合模型不同，不设置唯一的 \(R^2\) 阈值；若斜率与零在误差内不可区分、残差明显系统弯曲，或核心点由参数边界决定，则不能判为可信 Zeeman-like branch。

把约 450 cm⁻¹ 的普通高频 Raman 分量作为主要 Zeeman branch，并据此得到接近 0 \(\mu_B\) 的结果，不符合这组数据的核心磁场响应。

## 3. 低频 QES 与温度演化

温度序列使用 LL helicity-conserved 通道。对正 Raman shift 的 Stokes 分支，应先通过 Bose occupation 修正构建 Raman response：

\[
\chi''(\omega,T)\propto\frac{I_S(\omega,T)}{1+n_B(\omega,T)}.
\]

随后计算有限窗口的 Kramers–Kronig QES susceptibility proxy：

\[
\chi_{QES}(T)=\frac{2}{\pi}\int_{\omega_{min}}^{\omega_{max}}\frac{\chi''(\omega,T)}{\omega}\,d\omega.
\]

普通面积 \(\int\chi''d\omega\) 不能替代本题要求的 \(\int\chi''/\omega\,d\omega\)。

公开数据中的低频连续散射主要延伸至约 100 cm⁻¹。一个参考实现采用 **5–100 cm⁻¹**，不额外扣除温度依赖常数基线；同时改变为 10–100、5–90 和 10–90 cm⁻¹ 进行 cutoff sensitivity 检查。

主窗口下，11.6 K 的有限窗口 \(\chi_{QES}\) 约为 200 K 的 **3.91 倍**。逐点存在实验散布，但总体表现为**降温时低频 QES 增强**，与公开论文报告的单调低温增强趋势一致。

### 温度模型参考与稳健性

以

\[
\chi(T)=\chi_0+\frac{C}{T-\Theta}
\]

作为 susceptibility-like 参考模型，5–100 cm⁻¹ 实现得到：

- \(\Theta\approx-18.4\pm9.7\) K；
- \(R^2\approx0.929\)；
- AICc \(\approx-28.57\)。

同一数据的简单替代模型参考为：

| 模型 | \(R^2\) | AICc |
|---|---:|---:|
| constant | 0.000 | 2.46 |
| linear in \(T\) | 0.597 | −7.49 |
| linear in \(1/T\) | 0.852 | −21.57 |
| Curie–Weiss-like | **0.929** | **−28.57** |

改变积分窗口后，参考 \(\Theta\) 为：

- 5–100 cm⁻¹：−18.4 ± 9.7 K；
- 10–100 cm⁻¹：−21.0 ± 11.0 K；
- 5–90 cm⁻¹：−18.0 ± 9.5 K；
- 10–90 cm⁻¹：−20.7 ± 10.8 K。

因此，负的 \(\Theta\) 在这一组简单 cutoff 变化下保持同一量级，但由于有限频率窗口、背景选择与低频截止会影响参数，**\(\Theta\) 不设硬性评分单值**。满分要求是完成 Bose correction、\(1/\omega\) Kramers–Kronig 权重、合理温度模型比较和至少一种稳健性检查，并准确陈述参数限制。

背景处理如果导致大量负 spectral weight，或把负值直接裁剪为零而没有进行偏差检查，不能被视为稳健的定量结果。

## 4. 跨数据物理结论与证据边界

三组 Raman 数据共同支持：

1. Fe 掺杂在原有约 265–267 cm⁻¹ hybridized Raman 响应之外引入约 283 cm⁻¹ 的额外谱学分量；
2. 约 265 cm⁻¹ 和约 283 cm⁻¹ 两套分量均具有 helicity-dependent、随磁场近线性变化且在反向磁场下交换的频移；
3. 对应有效 phonon magnetic moment 分别处于约 **2.2 \(\mu_B\)** 和 **2.8 \(\mu_B\)** 的量级；
   这里的有效 phonon magnetic moment 描述的是相关 Raman/phonon hybridized branch 对外磁场的能量响应尺度，不能直接等同为样品的静态宏观磁矩或磁化强度；
4. helicity-conserved 低频 QES 在降温时总体增强，其 susceptibility-like spectral weight 用 Curie–Weiss-like 形式优于简单常数或线性描述，支持显著的低频磁性/自旋涨落响应。

但本题只提供 Raman 数据，因此不能仅据此严格证明某一种零场长程磁有序，也不能把外场 Zeeman-like splitting 或 QES 自动等同于铁磁、反铁磁或其他特定长程磁基态。关于长程磁有序仍需要独立磁性表征。

## 5. 参考结果的使用原则

上述数值用于核验本题公开数据能够恢复出的核心物理结果。不同合理的峰型、背景、拟合窗口和跨磁场约束方式可以造成小幅数值差异，因此评分以谱学对象识别正确、磁场与圆偏振依赖关系自洽、有效磁矩处于正确量级、低频响应物理量定义正确以及结论证据边界合理为主。

参考答案不规定唯一作图样式，也不要求逐点复现单一拟合实现。对于温度模型参数，尤其是有限窗口得到的 Curie–Weiss 参数，应结合频率窗口和背景处理的不确定性解释，不作为脱离处理条件的唯一硬数值。