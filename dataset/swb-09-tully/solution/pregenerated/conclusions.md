# 结论 (conclusions.md)

逐项回答任务目标 4，数值引自 `results.csv`。

## 1. nonadiabatic coupling 是否使 exact quantum 总速率相对 adiabatic 降低？

**是，且为小幅降低，两个温度一致。**

- 143 K：$\log_{10}k_{\mathrm A}=-1.166$，$\log_{10}k_{\mathrm{NA}}=-1.298$，
  即 $k_{\mathrm{NA}}<k_{\mathrm A}$，降低约 0.132 dex（约
  26.2%）。
- 250 K：$\log_{10}k_{\mathrm A}=0.193$，$\log_{10}k_{\mathrm{NA}}=0.113$，
  同样 $k_{\mathrm{NA}}<k_{\mathrm A}$，降低约 0.080 dex（约
  16.8%）。

nonadiabatic coupling 只造成小幅降低，与论文正文结论一致。

## 2. RPSH-BA 是否复现 Fig. 2 的量级与趋势？

**是。**

- 143 K：$\log_{10}k_{\mathrm{BA}}=-1.248$（SE=0.089 dex），
  与 Fig. 2 绿色圆点约 -1.30 相差 0.052 dex，在统计误差内一致；
  且 $k_{\mathrm{BA}}$ 落在 $k_{\mathrm{NA}}$ 附近、略低于 $k_{\mathrm A}$。
- 250 K：$\log_{10}k_{\mathrm{BA}}=0.132$（SE=0.060 dex），
  与 Fig. 2 绿色圆点约 0.14 相差 0.008 dex，一致；
  同样介于 $k_{\mathrm{NA}}$ 与 $k_{\mathrm A}$ 之间。

因此 RPSH-BA 在两个温度点都复现了论文 Fig. 2 的量级与“nonadiabatic 小幅低于
adiabatic”的趋势。

## 3. 温度趋势

从 143 K 到 250 K，三种方法的 $\log_{10}k$ 均显著上升（Arrhenius 行为），
在 $1000/T$ 从 6.99 到 4.00 的区间内三条曲线保持相同的相对次序
$k_{\mathrm{NA}}\lesssim k_{\mathrm{BA}}\lesssim k_{\mathrm A}$（143 K）
与 $k_{\mathrm{NA}}<k_{\mathrm{BA}}<k_{\mathrm A}$（250 K）。
