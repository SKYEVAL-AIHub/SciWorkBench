# 论文方法与本题采用的公式

## 1. 文献来源与使用范围

本文件用中文转述完成本题所需的论文理论、公式、算法和计算细节，不复制论文叙述性原文或原图。

Philip Shushkov, Richard Li, and John C. Tully, “Ring polymer molecular dynamics with surface hopping,” *The Journal of Chemical Physics* **137**, 22A549 (2012), DOI: https://doi.org/10.1063/1.4766449.

本题采用论文的 bead approximation（BA），不采用 centroid approximation（CA）。论文公式中的 bead 数 $N$ 在本题输出中记为 $P$。

## 2. 量子反应速率理论

论文从电子与核的 Hamiltonian 出发：

$$
\hat H(r,R)=\hat T(R)+\hat H_{\mathrm{el}}(r;R).
$$

反应物区和产物区由 $s(R)=0$ 分隔，$s(R)<0$ 为反应物，$s(R)>0$ 为产物。量子速率是 Kubo-transformed flux-side correlation function 的长时间极限除以反应物 partition function：

$$
k_{\mathrm{QM}}=
\frac{1}{Q_A}\lim_{t\rightarrow\infty}C_{fs}(T,t),
$$

$$
C_{fs}(T,t)=
\int_0^\beta\frac{d\lambda}{\beta}
\operatorname{tr}
\left[
e^{-(\beta-\lambda)\hat H}
\hat F
e^{-\lambda\hat H}
e^{i\hat Ht}
\hat h
e^{-i\hat Ht}
\right].
$$

$\hat h=h[s(R)]$ 是产物投影，$\hat F=i[\hat H,\hat h]$ 是 flux operator。

Adiabatic electronic states 满足

$$
\hat H_{\mathrm{el}}(r;R)\phi_i(r;R)
=V_i(R)\phi_i(r;R).
$$

Adiabatic approximation 只保留 ground adiabatic PES，核 Hamiltonian 为

$$
\hat H_{\mathrm{ad}}=\hat T(R)+\hat V_0(R).
$$

对一维气相自由反应物，单位长度的 translational partition function 为

$$
Q_A=\left(\frac{M}{2\pi\beta}\right)^{1/2}.
$$

本题的 `exact_adiabatic` 对应 ground adiabatic PES 上的精确量子散射；`exact_nonadiabatic` 对应完整两态 Hamiltonian 的精确量子散射。

## 3. Classical FSSH

电子波函数在瞬时 adiabatic basis 中展开为

$$
\psi(r,R,t)=\sum_i c_i(t)\phi_i(r;R).
$$

电子振幅满足论文 Eq. 17：

$$
i\dot c_i=V_i(R)c_i
-i\sum_j\dot R\,d_{ij}(R)c_j,
$$

其中 derivative coupling 为论文 Eq. 18：

$$
d_{ij}(R)=
\langle\phi_i(r;R)|\nabla_R|\phi_j(r;R)\rangle.
$$

定义

$$
a_{ij}=c_ic_j^*,
$$

$$
b_{ji}=-2\operatorname{Re}
\left(a_{ji}^*\dot R\,d_{ji}\right).
$$

从 active state $i$ 到 $j$ 的 FSSH probability 是论文 Eq. 19：

$$
g_{ij}=\frac{\Delta t\,b_{ji}}{a_{ii}},
$$

负值设为零。若发生 hop，核从 $V_i$ 切换到 $V_j$，沿 nonadiabatic coupling direction 调整动量以守恒总能量；若可用动量不足，hop 被拒绝。

论文把 classical nonadiabatic rate 写成多 PES transition-state flux 与 state-selected transmission coefficients 的组合。`in`、$i$、`fi` 分别表示远反应物区、dividing surface 和远产物区的电子态。

## 4. Adiabatic RPMD

论文 Eq. 29 的 $P$-bead Hamiltonian 为

$$
H_P^{\mathrm{ad}}=
\sum_{k=1}^{P}\frac{p_k^2}{2M}
+\sum_{k=1}^{P}\frac{M\omega_P^2}{2}(R_{k+1}-R_k)^2
+\sum_{k=1}^{P}V(R_k),
$$

其中

$$
R_{P+1}=R_1,
\qquad
\beta_P=\frac{\beta}{P},
\qquad
\omega_P=\frac{1}{\beta_P}=\frac{P}{\beta}.
$$

centroid position 与 momentum 是

$$
\bar R=\frac1P\sum_{k=1}^{P}R_k,
\qquad
\bar p=\frac1P\sum_{k=1}^{P}p_k.
$$

RPMD flux-side correlation function 在 extended ring-polymer phase space 中使用 $e^{-\beta_PH_P}$ 取样，dividing surface 依赖 $\bar R$。短时间极限给 centroid QTST：

$$
k_{\mathrm{QTST}}=
\frac{1}{Q_A^P}
\int\frac{d\mathbf R\,d\mathbf p}{(2\pi)^P}
e^{-\beta_PH_P}
\delta(\bar R-R^\ddagger)
\frac{\bar p}{M}h(\bar p).
$$

## 5. Ring-polymer surface hopping

### 5.1 Centroid approximation

论文 Eq. 35 的 CA electronic equation 为

$$
i\dot c_i=V_i(\bar R)c_i
-i\sum_j\dot{\bar R}\,d_{ij}(\bar R)c_j.
$$

CA 只在说明 BA 的来源和 $P=1$ 极限时使用；本题最终 RPSH rate 不采用 CA。

### 5.2 Bead approximation

论文把 ring polymer 看作一个有效“molecule”：所有 beads 始终位于同一个 active adiabatic PES，发生电子跃迁时所有 beads 同时切换。

BA electronic equation 是论文 Eq. 36：

$$
i\dot c_i=
\left[\frac1P\sum_{k=1}^{P}V_i(R_k)\right]c_i
-i\sum_j
\left[\frac1P\sum_{k=1}^{P}
\dot R_kd_{ij}(R_k)\right]c_j.
$$

相邻 beads 的 harmonic interaction 不出现在该 electronic equation 中，因为它不影响电子态之间的相对 phase。高温下 ring polymer 收缩，BA 退化为 CA。

核动力学仍由 active surface $i$ 上的 ring-polymer Hamiltonian 产生：

$$
H_P(\mathbf p,\mathbf R;i)=
\sum_{k=1}^{P}\frac{p_k^2}{2M}
+\sum_{k=1}^{P}\frac{M\omega_P^2}{2}(R_{k+1}-R_k)^2
+\sum_{k=1}^{P}V_i(R_k).
$$

BA hop 时沿由各 bead derivative coupling components 组成的多维 nonadiabatic coupling direction 调整 bead momenta，使整个 ring-polymer Hamiltonian 在 hop 前后守恒；无足够动能时拒绝 hop。

论文的 RPSH flux-side correlation function 为

$$
C_{fs}^{\mathrm{RPSH}}(T,t)=
\frac{1}{(2\pi)^P}
\sum_i\int d\mathbf R_0d\mathbf p_0
e^{-\beta_PH_P(\mathbf R_0,\mathbf p_0;i)}
\delta(\bar R_0-R^\ddagger)
\frac{\bar p_0}{M}
h(\bar R_t-R^\ddagger).
$$

rate 是

$$
k_{\mathrm{RPSH}}=
\frac1{Q_A^P}
\lim_{t\rightarrow\infty}
C_{fs}^{\mathrm{RPSH}}(T,t).
$$

state-selected 形式是论文 Eq. 40：

$$
k_{\mathrm{RPSH}}=
\frac1{Q_A^P}
\sum_{\mathrm{in},\mathrm{fi},i}
C_{fs,i}^{P}(T,0^+)
\kappa_{\mathrm{in},\mathrm{fi}}^i.
$$

## 6. MPES-QTST

论文 Eq. 41 定义

$$
k_{\mathrm{MPES-QTST}}=
\frac1{Q_A^P}\sum_i C_{fs,i}^{P}(T,0^+).
$$

为有效计算 rare event，选择反应物渐近区的 centroid surface $\bar R_0$ 和 interaction region 的 variational surface $\bar R_1$，并写成论文 Eq. 42：

$$
k_{\mathrm{MPES-QTST}}=
\sum_i k_i(\bar R_0)
p_i(\bar R_1,\bar R_0).
$$

本一维模型在反应物渐近区有两个简并电子态，因此论文 Eq. 45 为

$$
k_i(\bar R_0)=
\frac12(2\pi M\beta)^{-1/2}.
$$

两 dividing surfaces 的通量比可以化为 centroid probabilities 之比：

$$
p_i(\bar R_1,\bar R_0)=\frac{P_i(1)}{P_i(0)},
$$

$$
P_i(\xi)=
\frac{
\int d\mathbf R\,d\mathbf p\,e^{-\beta_PH_P(\mathbf R,\mathbf p;i)}
\delta(\bar R-\bar R_\xi)
}{
\int d\mathbf R\,d\mathbf p\,e^{-\beta_PH_P(\mathbf R,\mathbf p;i)}
}.
$$

定义论文 Eq. 48 的 potential of mean force：

$$
A_i(\xi)=-\frac1{\beta_P}\ln P_i(\xi),
$$

则论文 Eq. 49 为

$$
p_i(\bar R_1,\bar R_0)=
e^{-\beta_P[A_i(1)-A_i(0)]}.
$$

论文采用 umbrella integration：在 $\bar R_0$ 到 $\bar R_1$ 之间设置相互重叠的 harmonic-restraint windows，在每个 surface、每个温度分别运行 constant-temperature ring-polymer sampling，再由各窗口 centroid distribution 的前两阶矩恢复 mean force。

## 7. State-selected transmission coefficient 与 backward rare-event algorithm

论文 Eq. 50 定义

$$
\kappa_{\mathrm{in},\mathrm{fi}}^i=
\lim_{t\rightarrow\infty}
\frac{
\int d\mathbf R_0d\mathbf p_0\,e^{-\beta_PH_P(\mathbf R_0,\mathbf p_0;i)}
\delta(\bar R_0-R^\ddagger)
\delta_{\mathrm{in},n(-t)}
\frac{\bar p_0}{M}
h(\bar R_t-R^\ddagger)
\delta_{\mathrm{fi},n(t)}
}{
C_{fs,i}^{P}(T,0^+;R^\ddagger)
}.
$$

计算步骤按论文 Sec. VI B：

1. 将 centroid 固定在 interaction region 的 dividing surface；在 surface $i$ 上运行 constant-temperature ring-polymer sampling，消除 centroid motion；
2. 对独立 bead configurations 从 Maxwell distribution 重新抽取 momenta；
3. 以这些 phase points 为起点，先向后、再向前作 constant-energy RPSH propagation；
4. 实时传播时使用 physical bead masses；论文采用 free-ring evolution 与 PES force updates 交替的 multiple-time-step integrator；
5. 论文用 RK4 传播电子振幅，电子 Hamiltonian 在 nuclear long step 的起点、中点和终点更新。

为改善大 energy gap 时的电子传播稳定性，论文分离 adiabatic phases。CA 的 interaction-picture amplitude 是 Eq. 51：

$$
\tilde c_i(t)=c_i(t)
\exp\left[i\int_{-\infty}^{t}V_i(t')dt'\right],
$$

满足 Eq. 52：

$$
\dot{\tilde c}_i=
-\sum_j\dot{\bar R}d_{ij}
\exp\left[-i\int_{-\infty}^{t}(V_j-V_i)dt'\right]
\tilde c_j.
$$

BA 使用 Eq. 36 对应的 bead-averaged $V_i$ 和 $\dot R_kd_{ij}$。

dividing surface 位于 nonadiabatic coupling region，不能预先知道正确的 electronic amplitudes。Backward path 采用论文 Eq. 53 的 trial surface-hopping probability：

$$
f_{ij}=\frac12
\left[1-exp\left(-
\left|\frac1P\sum_k\dot R_kd_{ij}(R_k)\right|
\Delta t\right)
\right].
$$

沿同一 backward nuclear path 同时传播两个线性独立的 electronic solutions，并保存 electronic Hamiltonian。到达 nonadiabatic coupling 消失的反应物渐近区后，对指定 incoming state 施加论文 Eq. 54：

$$
c_{\mathrm{in}}(-\infty)=1,
\qquad
c_{j\ne\mathrm{in}}(-\infty)=0.
$$

满足该渐近条件的解由两个 backward independent solutions 的线性组合唯一确定。由此得到 dividing-surface amplitudes 和整条 backward path 上正确的 amplitudes，再计算 FSSH target probabilities $g_{ij}$。

Backward path 的 reweighting factor 是论文 Eq. 55：

$$
W=\prod_\mu w_\mu.
$$

若 step $\mu$ 尝试 $i\to j$ switch，

$$
w_\mu=\frac{g_{ij}}{f_{ij}}.
$$

若未尝试 switch，论文 Eq. 56 为

$$
w_\mu=
\frac{1-\sum_{j\ne i}g_{ij}}
{1-\sum_{j\ne i}f_{ij}}.
$$

得到 dividing-surface electronic initial conditions 后，forward path 使用正常 FSSH probability。因为同时向后和向前传播，初始 centroid momentum 可只取正 outward components。每条 reactive trajectory 的贡献还要除以它向前穿过 dividing surface 的次数。state-selected transmission coefficient 是平均 reactive flux 与 dividing surface 总 thermal flux 的比值。

## 8. 一维两态模型

论文 Eq. 57 的 diabatic potential matrix elements 为

$$
V_{11}=A e^{-BR^2},
$$

$$
V_{22}=C e^{2E(R-D)}+2C e^{E(R-D)},
$$

$$
V_{12}=F e^{-GR^2}.
$$

参数均为原子单位：

$$
M=2000,
\quad A=5\times10^{-3},
\quad B=2,
\quad C=10^{-5},
$$

$$
D=-2.8,
\quad E=1.2,
\quad F=10^{-3},
\quad G=2.
$$

逐点对角化 diabatic matrix 得到两个 adiabatic PES。ground PES 有 barrier，excited PES 为 repulsive；负 $R$ 渐近区两 adiabatic states 简并，$|R|$ 较大时 nonadiabatic coupling 消失。论文把 dividing surface 放在 ground adiabatic barrier top，反应物渐近区在其左侧，产物渐近区在其右侧。

## 9. 论文 Appendix 的计算细节

### 9.1 Exact quantum rates

论文通过 constant-energy scattering amplitudes 的模平方作数值积分得到 exact rates。Scattering amplitudes 使用扩展到多 PES 的 Feit-Fleck split-operator Fourier-transform wave-packet propagation：

- wave packet 初始化在远反应物区；
- 向前传播至 scattering event 完成；
- 反应物渐近区放置 absorbing potential 消除 reflected component；
- 在产物渐近区投影 transmitted component；
- 对 transmitted fraction 作 Fourier transform，并进行通量归一。

论文计算了 $50$–$5000\ \mathrm K$ 的九个温度点。本题复现其中的 $143\ \mathrm K$ 和 $250\ \mathrm K$ 两个温度点。

### 9.2 Beads、QTST sampling 与时间步

论文在 $T\ge143\ \mathrm K$ 时使用 $P=32$，在 $T\le100\ \mathrm K$ 时使用 $P=128$。因此本题的 $143\ \mathrm K$ 和 $250\ \mathrm K$ 均使用 $P=32$。

MPES-QTST 使用 100 个 umbrella windows：ground PES 的 centroid range 为 $-3a_0$ 到 $3a_0$，excited PES 为 $-3a_0$ 到 $0.5a_0$。Harmonic restraint force constant 为

$$
1.22\ \mathrm{eV}/a_0^2.
$$

每个 window 先平衡 $2.5\times10^5$ steps，再采样 $10^7$ steps。centroid-mode friction 取 umbrella potential frequency 的两倍。论文在 $143\ \mathrm K$ 和 $250\ \mathrm K$ 使用的长时间步均为

$$
2.4\times10^{-3}\ \mathrm{fs}.
$$

### 9.3 Transmission sampling

论文先平衡 $10^6$ steps；在 centroid 固定于 dividing surface 的 ensemble 中每隔 $2.5\times10^3$ steps 抽取 ring-polymer configuration，并从 Maxwell distribution 重新抽取 momenta。每个 PES、每个温度运行 $10^6$ 条 forward/backward trajectories。

在 $T\le250\ \mathrm K$ 时，free-ring-polymer integration 的内部时间步相对长时间步再缩小 10 倍。Classical surface hopping 使用 $P=1$；adiabatic RPMD 在 ground PES 上向前和向后传播。

## 10. 本题的复现范围

本题只复现论文 Fig. 2 中 $143\ \mathrm K$ 和 $250\ \mathrm K$ 的 exact adiabatic、exact nonadiabatic 与 RPSH-BA 总速率点，不要求复现其余温度、CA、adiabatic RPMD、classical surface hopping 或 Figs. 3–5 的 state-selected rates。

珠子数和动力学时间步采用论文 Appendix 对这两个温度给出的设置。论文的生产采样规模为每个 PES、每个温度 $10^6$ 条 forward/backward trajectories；考生应在个人计算机时限内自行选择取样规模、独立分块和收敛检查，并报告 RPSH-BA 统计标准误。RPSH-BA 方法身份、Eq. 36、Eq. 40、Eq. 50–56、whole-ring hop、渐近 incoming-state 初始化和 crossing correction 均保持论文定义。
