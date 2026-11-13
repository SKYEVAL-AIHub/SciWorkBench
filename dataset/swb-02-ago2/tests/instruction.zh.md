# AGO2 扰动后的 RNA 变化来自转录还是转录后调控？

## 研究背景

AGO2 是 miRNA-RISC 的核心组分，也被提出可能参与非经典的转录调控。在一项长期 CRISPRi 扰动实验中观察到 AGO2 扰动伴随多种 RNA 表达变化，提出假设：AGO2 是普遍性的直接转录抑制因子，AGO2 敲低会广泛解除转录抑制。

`data/` 提供了从真实 PerturbSci-Kinetics 公共数据整理得到的多层 RNA 表达数据和扰动 metadata。该数据同时包含新生 RNA 与已有 RNA 层，因此原则上可以检验上述机制主张，但也可能受数据质量与分析选择影响。

## 研究问题

请独立评价上述假设。自行检查输入数据、确定合适的分析单位和对照，选择适合该数据的分析方法，回答：

1. AGO2 扰动后的 whole-RNA 变化主要来自新生 RNA、已有 RNA，还是两者共同变化？现有数据能否可靠地区分这些可能？
2. 与其他扰动和 non-targeting controls 相比，AGO2 的变化模式是否具有特异性？是否存在数据质量、分析选择或其他合理解释？
3. 哪些结果是稳定可信的，哪些结论仍然证据不足？本数据是否足以接受上述“普遍性直接转录抑制”主张？

所选方法及证据强度必须与数据的实际结构相匹配；无法可靠估计的内容应明确标记为证据不足。

## 任务输入

`data/` 下全部文件：

- `ago2_perturbsci_guide_pseudobulk.h5ad`：多层 RNA 表达数据（guide 级 pseudobulk）；
- `layer_metadata.tsv`：各 RNA 层的定义；
- `guide_profile_metadata.tsv`：guide 与扰动目标的对应信息；
- `study_context.tsv`：研究背景与设计信息；
- `input_manifest.tsv`、`source_provenance.tsv`：题内文件清单与来源追溯；
- `GSE218566_series_record.txt`、`GSE218566_filelist.txt`、`GSM6752591_record.txt`：公共数据归档记录。

## 任务输出

将全部结果写入 `output/`，包括：

- 可机器读取的完整基因层结果，以及能够复核主要比较、质量和稳健性判断的汇总表；
- 支持主要判断的必要图表；
- 简短分析报告，说明方法选择、关键结果、替代解释和最终结论；
- 可从 `data/` 重现主要表图的分析代码及软件环境记录。

### 格式与合规要求

不得使用已有论文结论替代本次数据分析，不得填造输入中不存在的结果。
