# 计算方法 (method.md)

## 模型与 adiabatic representation

一维两态 diabatic Hamiltonian 按论文 Eq. 57 实现，核质量 $M=2000$ a.u.。对 diabatic
矩阵做解析对角化得到两条 adiabatic potential energy surfaces 与 derivative coupling，
采用连续实 gauge（沿网格固定本征矢符号），并验证本征方程残差、正交性与
$d_{ij}=-d_{ji}$。dividing surface 位于 lower adiabatic surface 势垒顶
（143 K: 0.047505 bohr；250 K: 0.047505 bohr），
反应物区在其左侧，产物区在其右侧。

## Exact quantum 散射

- `exact_adiabatic`：仅保留 lower adiabatic surface 的单通道开放边界散射；
- `exact_nonadiabatic`：完整两通道耦合散射，按入射通道分别做通量归一后再热平均；
- 空间边界：两级网格 $x\in[-12,12]$ 与 $[-16,16]$ bohr（ngrid 1600 / 2400），
  两级的 rate 变化 <2e-5 dex；
- 能量求积：Gauss-Legendre 热积分从 64 节点细化到 96 节点，收敛到 <2e-5 dex；
- 热积分尾部在 12 $k_BT$ 以上截断，贡献可忽略；
- 单位换算：$k_{\mathrm m}=k_{\mathrm{au}}\times(a_0/t_0)$，$a_0/t_0=2.18769126364\times10^6\ \mathrm{m\,s^{-1}}$。

## RPSH bead approximation (P=32)

- 两个温度均取 $P=32$，长时间步 0.0024 fs；
- constrained dividing-surface sampling 生成 ring-polymer 初始构型；
- Eq. 36 bead-averaged electronic propagation，whole-ring hop，backward trial path，
  渐近 incoming-state reconstruction，Eq. 55-56 reweighting，forward propagation 与
  crossing correction；
- 每个 (surface, incoming) 通道 100 条轨迹、分 4 个独立 chunk（共 16 chunk / 400 轨迹），
  用 within/between-batch block standard error 估计 $\log_{10}$ 标准误。

## 取样规模与不确定度

| 温度 | exact A $\log_{10}k$ | exact NA $\log_{10}k$ | BA $\log_{10}k$ | BA SE (dex) | BA 反应轨迹 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 143 K | -1.166 | -1.298 | -1.248 | 0.089 | 98 |
| 250 K | 0.193 | 0.113 | 0.132 | 0.060 | 99 |

exact quantum 为标准散射计算，无统计误差（SE=0）；RPSH-BA 报告有限统计标准误。
电子范数由 midpoint_exponential 积分器保持守恒，whole-ring hop 通过 bead 动量重标
保持总能量守恒，二者残差均在机器精度量级（见 diagnostics.json）。
