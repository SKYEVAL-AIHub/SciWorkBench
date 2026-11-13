# 任务需求

## 任务背景

Shushkov、Li 和 Tully 在 2012 年提出了 ring-polymer surface hopping 方法，并用一个一维两态模型检验量子统计效应与 nonadiabatic dynamics 对总反应速率的共同影响。论文 Fig. 2 给出了 exact quantum adiabatic、exact quantum nonadiabatic 和 bead approximation ring-polymer surface hopping 三种总速率随温度变化的结果。

本任务要求根据给定的物理模型、论文公式和计算条件，独立完成 143 K 与 250 K 两个温度点的计算，复现 Fig. 2 中上述三组数据。任务重点是从模型势能面出发，完成开放边界量子散射、热平均、ring-polymer 统计取样、surface-hopping 动力学和结果比较，而不是读取或数字化论文图片。

## 任务目标

1. 计算 143 K 与 250 K 下的 `exact_adiabatic` 和 `exact_nonadiabatic` 总反应速率；
2. 计算相同温度下、$P=32$ 的 `rpsh_ba` 总反应速率及其统计标准误；
3. 将六个速率点绘制为 Fig. 2 的局部 Arrhenius 图，横轴使用 $1000/T\ (\mathrm{K^{-1}})$，纵轴使用 $\log_{10}[k/(\mathrm{m\,molecule^{-1}\,s^{-1}})]$；
4. 比较两个温度下三种速率的相对关系，判断 nonadiabatic coupling 是否使 exact quantum 总速率相对 adiabatic approximation 降低，并判断 RPSH-BA 是否复现论文 Fig. 2 的相同趋势；
5. 提供可重新运行的完整代码、数值设置和必要诊断，使所有表格、图和结论能够由同一套最终计算复现。

## 要求输入

- `data/problem.json`：模型参数、物理常数、目标温度、珠子数和单位约定；
- `data/scientific_spec.md`：完成任务所需的论文背景、公式、算法和计算条件；
- `data/output_schema.json`：输出文件的字段和格式要求。

计算条件：

- 模型为一维两态 diabatic Hamiltonian，核质量为 $M=2000$ a.u.；
- 目标温度为 $T=143\ \mathrm K$ 和 $T=250\ \mathrm K$；
- `exact_adiabatic` 只保留 lower adiabatic potential energy surface；
- `exact_nonadiabatic` 使用完整两通道 Hamiltonian；
- `rpsh_ba` 使用论文的 bead approximation，两个温度均取 $P=32$；
- dividing surface 位于 lower adiabatic potential energy surface 的势垒顶，反应物区在其左侧，产物区在其右侧；
- 内部计算采用原子单位，最终总速率同时报告 atomic units 和 $\mathrm{m\,molecule^{-1}\,s^{-1}}$；
- 最终表格中的 $\log_{10}k$ 保留 3 位小数，统计标准误保留 3 位小数。

## 要求输出

在 `output` 文件夹下生成以下文件：

1. `results.csv`：六个目标速率点及统计误差，字段严格遵循 `data/output_schema.json`；
2. `figure2_subset.png`：143 K 与 250 K 的局部 Arrhenius 图，分辨率不低于 150 dpi；
3. `method.md`：计算方法、边界条件、离散设置、热积分设置、取样规模和不确定度估计；
4. `conclusions.md`：逐项回答任务目标 4，并引用 `results.csv` 中的数值；
5. `diagnostics.json`：模型、exact quantum 和 RPSH-BA 的数值诊断；
6. `run_all.py`：从输入文件重新生成全部结果的单一入口；
7. `src/`：完成计算所需的全部自编源代码。
