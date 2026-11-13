# 评分标准：仅根据最终结果评分

总分 100 分。只检查 `output/` 中的机器表、图和报告，不检查实际工具调用。允许采用不同但科学合理的方法；参考数值用于核验，不要求逐位复现。替代方法必须能从提交表格判断其估计量、分析单位、内部反事实和证据边界。

数据结构、样本映射、分析单位、AGO2 target engagement、因果证据边界和最终判断不单独分配过程分，而作为全题成立的基础条件及封顶依据。主要分值集中在以下三个高区分度结果模块。

## 1. 三层依赖、归一化估计量与基因层归因（40 分）

检查完整基因表和报告是否识别 `whole = nascent + preexisting`，并实际处理层内相对归一化与保持构成关系的效应不是同一估计量。高分答案应提供可闭合的基因层成分归因或科学上等价的 compositional/counterfactual 结果，而不是只给三组相关。参考主/严格过滤为 4240/1296；参考共同尺度分解的 520 个大总效应基因中 186 个 nascent-dominant、334 个 preexisting-dominant，254 个大幅正向变化中 100 个 nascent-dominant、220 个两成分均正，total effect 与两贡献相关约 0.584/0.790。合理替代方法按估计量定义、低表达控制、闭合性、机器可追溯性和结论稳定性评分。只做每层各自归一化后以相关性代替成分归因，本项不得高于 20 分；把派生层当作独立实验或直接降解率，本项不得高于 10 分。

## 2. 内部反事实与 AGO2 特异性校准（35 分）

检查机器表是否把同一判断准则应用到 targeting guides、target aggregates 和严格 held-out NTC pseudo-queries，以无 query 泄漏的内部 null/安慰剂评估 AGO2 特异性。参考 whole–nascent 相关在 8 个 target aggregates 中将 AGO2 排第 1，在 20 条 targeting guides 中仅排第 7；15 个 held-out NTC queries 中有 8 个不低于 AGO2。AGO2 正向 nascent 大效应比例在 targeting guides 中排第 17，86.7% 的 NTC queries 不低于 AGO2。可接受其他不造成伪重复的 cross-fitting、嵌套置换或经验 null 设计。只比较 AGO2 对 NTC、只列相关 perturbations，或把 target 聚合排名直接当作特异性证明，本项不得高于 14 分；反事实存在信息泄漏且影响经验显著性解释时不得高于 24 分。

## 3. Profile 质量、稳健性、降维与结果一致性（25 分）

检查 QC、稳健性、PCA 坐标和图是否识别唯一可靠 AGO2 guide `AGO2_3`、AGO2 concordance 不可估计、其他 target whole concordance约 0.064–0.215、TAF2_1/RPL11_1/RPL11_2 低于 50 个聚合细胞，以及 PC1/PC2 方差约 45.4%/8.4%且受到 profile 大小或 nascent fraction 影响。应比较合理的替代过滤或归一化，直接报告其相对主分析的效应一致性，并使完整表、坐标、图和报告方向一致。没有 PCA 机器坐标、替代设置或跨文件核对，本项不得高于 10 分；只提交图或报告而缺少可重算表格，本项不得高于 6 分。

## 总分封顶规则

- 使用模拟数据、编造结果或读取隐藏答案：0 分。
- AGO2/NTC 映射或 guide-level 分析单位错误：总分最高 50 分。
- 制造细胞级伪重复或不存在的 AGO2 guides：总分最高 55 分。
- 仅分析 whole RNA：总分最高 55 分。
- 将 preexisting 视为独立 half-life/degradation-rate 测量：总分最高 60 分。
- 未处理跨层构成依赖且据此作确定性机制归因：总分最高 65 分。
- 没有内部反事实/安慰剂校准却声称 AGO2-specific：总分最高 65 分。
- 把表达关联或 target engagement 直接写成普遍性直接因果机制：总分最高 60 分。
- 只提交图片而没有完整基因结果表：总分最高 55 分。
- 没有机器可追溯的最终判断：总分最高 80 分。
