## 评分标准（过程+结果评分）

允许检查最终文件、源码、运行记录和计算过程。总分 100 分。

### check point1(得分点1)——模型与 adiabatic representation (10分)

- 10分：从 Eq. 57 正确实现 diabatic potential matrix、解析或受控数值导数、连续实 adiabatic gauge、两条 adiabatic potential energy surfaces 和 derivative coupling，并验证本征方程、正交性与 $d_{ij}=-d_{ji}$。
- 6分：主体正确但缺一项连续性或导数检查。
- 2分：势能面可计算，但 coupling 或 gauge 处理不可靠。
- 0分：模型参数、符号或势能公式错误。

### check point2(得分点2)——exact quantum 计算流程 (20分)

- 20分：实际完成 lower-surface 单通道和完整两通道开放边界散射，按入射通道做通量归一和热积分，并对空间边界、离散或能量求积进行真实细化。
- 12分：两种 exact 结果正确，但缺一类收敛或守恒检查。
- 5分：使用经验缩放、WKB 或封闭边界近似，并明确标记局限。
- 0分：未执行 exact quantum 计算或直接填入参考值。

### check point3(得分点3)——RPSH-BA 方法身份 (20分)

- 20分：实际执行 $P=32$ constrained dividing-surface sampling、Eq. 36 bead-averaged electronic propagation、whole-ring hop、backward trial path、渐近 incoming-state reconstruction、Eq. 55–56 reweighting、forward propagation 和 crossing correction。
- 12分：主体为 BA，但一项 rare-event 修正不完整。
- 5分：使用普通单珠 FSSH、centroid approximation 或短时 surrogate，并明确标记失败。
- 0分：未执行 RPSH-BA 或直接填入参考值。

### check point4(得分点4)——两温度独立计算与统计处理 (15分)

- 15分：143 K 与 250 K 均独立运行；两个温度都使用 $P=32$；RPSH-BA 使用独立分块估计标准误，并保留有效历史、零贡献、反应、crossing 与 hop 计数。
- 9分：两温度均运行，但统计分块或一类计数不完整。
- 4分：只完成一个温度，或第二个温度由经验缩放得到。
- 0分：没有可核对的两温度计算记录。

### check point5(得分点5)——数值稳定性与物理诊断 (10分)

- 10分：检查开放通道通量、电子范数、传播幺正性或等价误差、hop 能量守恒、概率预算、权重有限性和热积分尾部。
- 6分：至少五类检查有可观察证据。
- 3分：只有退出码或少量单元测试。
- 0分：无诊断或诊断与源码不符。

### check point6(得分点6)——Fig. 2 六个目标点的准确性 (15分)

- 15分：N2–N5 的六个数值全部达到满分档。
- 10分：至少五个达到满分档，剩余值不低于第二档。
- 5分：至少四个值不低于第二档。
- 0分：少于四个值达到可接受分档。

### check point7(得分点7)——可复现性与输出一致性 (10分)

- 10分：固定随机种子；依赖和环境可安装；`run_all.py` 可从输入重建六行表格、诊断、图和结论；所有产物来自同一最终代码版本，且图表相互一致。
- 6分：可重跑但缺一项次要环境或随机性信息。
- 3分：代码完整但不能单入口重跑全部产物。
- 0分：结果依赖未提供文件、外部专用实现或手工填值。
