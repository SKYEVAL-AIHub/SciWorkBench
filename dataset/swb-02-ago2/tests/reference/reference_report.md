# AGO2 任务隐藏标准答案

## 数据与分析单位

输入包含 35 个 guide-level profiles、59429 个唯一 gene symbols 和 15 个 non-targeting controls。主过滤保留 4240 个在三层均有足够 NTC 表达的基因；严格过滤保留 1296 个。whole、nascent 和 preexisting 均为非负整数，且逐元素满足 `whole = nascent + preexisting`。

每个 guide 聚合谱是分析单位。该数据来自一个 pooled HEK293-idCas9 screen，没有独立供体重复。公开 on-target metadata 中只有 AGO2_3，因此 AGO2 guide concordance 无法估计，不能把单细胞数或 15 个 NTC profiles 解释为多个 AGO2 生物学重复。

## 扰动质量与三层效应

AGO2_3 的 AGO2 whole-layer 目标表达相对 NTC 中位数效应为 -2.0387 log2 units，支持该 guide 发生了 target knockdown，但没有第二条可靠 AGO2 guide 作独立复制。AGO2 profile 的总 nascent/whole 比例为 0.185817，NTC 中位数为 0.182380。

在主过滤基因中，AGO2 whole 与 nascent 效应的 Spearman 相关为 0.5370，whole 与 preexisting 为 0.7630，nascent 与 preexisting 为 0.0444。这些数值支持三层联合描述，但不能把 preexisting 层当作独立的 RNA degradation-rate 实验。

基于 15 个 NTC 的经验尾概率最小只能达到 1/16，且 AGO2 只有一个 guide profile。因此完整结果表中的 empirical p/q 只用于相对 NTC 离群程度的描述，不应被解释为具有独立生物学重复的正式差异表达显著性。

## 保持层间构成关系的归因

只分别对 whole、nascent 和 preexisting 做层内归一化再比较相关性，会混合“层内相对丰度变化”和“该层占 whole 的构成比例变化”。参考答案因此同时使用 whole 文库作为共同尺度，并把每个基因的 whole-scale 总变化分解为 nascent 与 preexisting 的对称贡献。该分解是一种可审核的参考实现，不要求参评者逐式复现；其他保持层间构成关系、明确估计量且能够闭合的实现也可接受。

主过滤的 4240 个基因中，共 520 个基因的 common-scale total effect 绝对值不小于 0.5，其中 254 个为正、266 个为负。按绝对贡献判断，186 个为 nascent-dominant，334 个为 preexisting-dominant；254 个大幅正向变化中有 100 个 nascent-dominant，220 个同时收到正向 nascent 和 preexisting 贡献。total effect 与 nascent、preexisting 贡献的 Spearman 相关分别为 0.5841 和 0.7900，最大闭合残差为 8.88e-16。

这些结果支持 mixed contribution，而不是把 whole–nascent 相关直接解释为独立转录证据，或把 preexisting 贡献直接解释为降解速率变化。

## 内部反事实、质量与稳健性

DROSHA、DICER1、XPO5、POLR2B、TAF2、RPL11 和 LRPPRC 的多 guide 数据可用于检查 target knockdown、guide 间一致性和 AGO2 signature 的相似性，但这些 perturbations 不是 AGO2 的生物学重复。不同 profile 聚合细胞数不完全相等，尤其低于 50 个细胞的 profile 应在解释 PCA 和 signature correlation 时标记。

同一个“whole–nascent 相关”判据在不同分析层级给出不同印象：将每个 target 的 guides 取中位数时，AGO2 在 8 个 targets 中排名第 1；将每条 targeting guide 都作为一次同样的 query 时，AGO2_3 只在 20 条 guides 中排名第 7；逐个留出 NTC 并用其余 14 个 NTC 建立参考时，15 个 pseudo-queries 中有 8 个的相关不低于 AGO2，即比例为 0.5333。AGO2 的正向 nascent 大效应比例在 20 条 targeting guides 中仅排名第 17，且 86.7% 的 NTC pseudo-queries 不低于 AGO2。

因此，target 聚合层面的高排名不能单独证明 AGO2-specific transcriptional signature。可接受的替代反事实不必与参考统计量完全相同，但必须把同一判据应用到适当的内部 queries，并说明 guide 层级、target 聚合层级与技术/生物学重复的差别。

主 CPM、严格过滤 CPM 和 NTC-reference poscounts 三种设置均已运行；完整敏感性结果见 `expected_sensitivity_analysis.tsv`。一致方向可提高描述性结论的可信度，但不会弥补 AGO2 单 guide 和无独立生物学重复的设计限制。

## 机制判断

whole 与 nascent 效应呈中等相关，构成保持的分解也确认部分基因有正向 nascent contribution，因此数据与部分转录层贡献相容。但多数大效应 whole 变化为 preexisting-dominant，很多正向变化由两种成分共同贡献；内部 pseudo-query 校准又显示相关性并非 AGO2 独有。经验 q 值没有达到常用显著性水平，长期 CRISPRi 的继发效应仍可能产生相似模式。preexisting 是 `whole - nascent` 的派生量，preexisting 贡献不能直接证明 RNA 降解减慢。

最终判断：该输入支持 AGO2_3 有效降低 AGO2 表达，并提供 AGO2 扰动伴随 nascent RNA 改变的弱描述性证据；它不足以建立 guide 可重复性、基因级正式显著性、直接启动子结合、RNA 半衰期变化或普遍性转录抑制因果机制。关于“AGO2 是普遍性的直接转录抑制因子”的强主张应判为未被本输入证明。
