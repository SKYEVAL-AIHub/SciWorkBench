# 任务：根据 NMR / NOESY / HRMS 谱图数据推断有机分子结构

## 研究背景

本任务涉及三个新合成的有机化合物，编号为 **compound_A、compound_B 和 compound_C**。
其中，compound_C 是 compound_A 经自由基串联环化反应形成的三环产物。
三个化合物均已完成常规表征，相关数据分别存放于：

```
data/compound_A.csv
data/compound_B.csv
data/compound_C.csv
data/compound_C_noesy.png   （compound_C 的相敏 NOESY 谱，用于确定其相对构型）
```

每个 CSV 文件包含该化合物的以下类别数据：

| spectrum | 含义 |
| --- | --- |
| `1H` | 1H NMR 峰表：化学位移（ppm，可为区间）、裂分、偶合常数 J（Hz）、积分（氢原子数） |
| `13C` | 13C NMR 峰表：化学位移（ppm） |
| `MS` | HRMS 实测 m/z；电离方式与离子类型见该行 `note` |

各列含义：`spectrum, shift, multiplicity, j_hz, integration, note`。
`note` 中给出核磁的测试频率与溶剂（均为 500 MHz / 126 MHz，CDCl3）。
注意：1H 的多重峰区间可能覆盖多个化学环境相近的氢；13C 个别信号可能对应两个等价碳（见相应 note）。
HRMS 仅提供实测 m/z 与离子类型，**未提供 calcd 的元素组成**，
分子式需由精确质量结合 NMR 的氢/碳计数推出。

`compound_C_noesy.png` 为 compound_C 的相敏 NOESY 谱（含 1D ¹H 投影）：
图中交叉峰反映质子间的空间邻近关系，结合 1H 归属推断 compound_C 的相对构型。

## 研究问题

1. 对三个化合物分别解析谱图数据，**推断其完整分子结构**。
   compound_A/B 只需报告构造（立体化学不考察）；
   compound_C 需报告**含立体化学的完整结构**（其相对构型由 NOESY 确定）。
2. 将每个结构写成 **SMILES**（compound_C 必须含 `@`/`@@` 立体标注），并用化学信息学程序
   （如 RDKit）由该 SMILES 计算出**标准 InChIKey**，一并提交。
3. 只交一个 CSV 文件，阅卷完全由脚本自动完成，不需要提交报告、解析过程或图片。

## 任务输入

#### 三、约定与立体化学判分口径

- SMILES 允许任意合法写法，不要求规范化形式；判分在规范化之后进行。
- **立体异构只在 compound_C 上考察**：
  - 你提交的 compound_C 立体结构与参考答案**相同**，或与其互为**对映异构体**
    （全部手性中心同时取反），均判为**立体正确**——反应得到外消旋体，NOESY 只能确定
    相对构型、无法区分对映体，因此提交整体镜像构型与参考答案等价。
  - 提交 compound_C 的**任何非对映异构体**（哪怕只错一个手性中心）判为**立体错误**。
  - compound_C **未指定或部分指定**立体化学的 SMILES 同样判为立体错误：
    应报告结构中的全部 6 个手性中心。
- compound_A/B 只考构造：其提交 SMILES 中写或不写立体标注均不影响任何得分点。
- 构造（不考虑立体）与立体是两个独立得分点：compound_C 构造对而立体错，仍可获得构造分。
- InChIKey 必须由你自己写出的 SMILES 用程序计算得到（不得手工编造）。
- CSV 使用 UTF-8 编码、逗号分隔、首行为表头。

## 任务输出

唯一交付物是 `output/structures.csv`，表头顺序固定：

```
compound_id,smiles,inchikey
```

共 3 行数据，`compound_id` 取以下固定值（顺序不限）：

| compound_id | 对应数据文件 |
| --- | --- |
| `compound_A` | `data/compound_A.csv` |
| `compound_B` | `data/compound_B.csv` |
| `compound_C` | `data/compound_C.csv` + `data/compound_C_noesy.png` |

各列含义：

| 列名 | 含义 |
| --- | --- |
| `compound_id` | 上表固定取值之一 |
| `smiles` | 该化合物的结构（compound_C 必须含立体化学标注） |
| `inchikey` | 由同一行 `smiles` 计算出的标准 InChIKey（27 字符，形如 `XXXXXXXXXXXXXX-YYYYYYYYYY-Z`） |
