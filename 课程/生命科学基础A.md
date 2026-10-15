---
课程名称: 生命科学基础A
学分: 2
总学时: 32
考核方式: 考试
课程属性: 选修
课程性质: 公共基础课程
开课单位: 生命学院
领域: 物质科学（物理·化学·生命）
tags:
  - 课程
  - 领域/物质科学（物理·化学·生命）
---

# 生命科学基础A


### 第1章 生命体系中的量纲分析与标度律

#### 1.1 量纲齐次性原理

定义：物理量 $Q$ 的量纲写作 $[Q]=\mathrm{M}^{\alpha}\mathrm{L}^{\beta}\mathrm{T}^{\gamma}\Theta^{\delta}\mathrm{N}^{\epsilon}\mathrm{J}^{\zeta}$，其中 M 为质量、L 为长度、T 为时间、$\Theta$ 为温度、N 为物质的量、J 为发光强度。

假设：描述生命现象的物理方程在单位制缩放（长度用 $\lambda$、质量用 $\mu$、时间用 $\tau$ 缩放）下形式不变。

公式：设方程 $f(q_1,\dots,q_n)=0$，若各量量纲为 $[q_i]=\prod_k X_k^{a_{ik}}$，则由 $\prod_i \lambda^{a_{ik}}$ 对每个基量纲为零的条件确定无量纲组合。

推导：设物理规律 $q_0 = C\, q_1^{p_1}\cdots q_n^{p_n}$。两边取量纲
$$
[q_0]=\prod_{i=1}^n [q_i]^{p_i}
\Rightarrow
\mathrm{M}^{\alpha_0}\mathrm{L}^{\beta_0}\mathrm{T}^{\gamma_0}
=\prod_i \big(\mathrm{M}^{\alpha_i}\mathrm{L}^{\beta_i}\mathrm{T}^{\gamma_i}\big)^{p_i}.
$$
比较各基量纲指数，得线性方程组
$$
\sum_i \alpha_i p_i=\alpha_0,\quad
\sum_i \beta_i p_i=\beta_0,\quad
\sum_i \gamma_i p_i=\gamma_0.
$$
若未知指数个数 $n$ 大于独立方程个数 $r$，则存在 $n-r$ 个自由参数，对应 $n-r$ 个无量纲群。

物理意义：物理规律中的“形状”由无量纲组合决定，尺度缩放只改变这些组合的数值。几何意义：在 $(q_1,\dots,q_n)$ 空间内，物理关系被约束在由量纲齐次方程组定义的低维流形上。

特例：$n=r$ 时指数唯一确定，如 $T=2\pi\sqrt{l/g}$ 中由 $[l]=\mathrm{L},\ [g]=\mathrm{LT^{-2}}$ 得 $T\propto\sqrt{l/g}$。

应用：估计细胞游动速度、心跳频率与体重的幂律、微流控芯片缩放设计。

#### 1.2 Buckingham π 定理及其证明

定义：若物理问题含 $n$ 个有量纲量，涉及 $r$ 个相互独立的基量纲，则问题等价于 $n-r$ 个无量纲量 $\pi_1,\dots,\pi_{n-r}$ 之间的关系 $\Phi(\pi_1,\dots,\pi_{n-r})=0$。

假设：所研究物理量集合中至少存在 $r$ 个量能张成全部基量纲（构成“量纲基”）。

公式：把 $n$ 个量写成矩阵 $\mathbf{A}=(a_{ik})\in\mathbb{R}^{r\times n}$，其秩为 $r$。取 $r$ 个线性无关列作基变量 $\mathbf{q}_b$，其余为 $\mathbf{q}_s$，则
$$
\pi_j = q_{s,j}\prod_{i=1}^{r} q_{b,i}^{-\,(\mathbf{A}_b^{-1}\mathbf{A}_s)_{ij}},\qquad j=1,\dots,n-r.
$$

推导：设 $\mathbf{q}=(q_1,\dots,q_n)$，对数化得量纲向量 $\boldsymbol{\alpha}_i\in\mathbb{R}^r$，$[\ln q_i]$ 的分量。任意幂次组合 $q_1^{p_1}\cdots q_n^{p_n}$ 无量纲当且仅当
$$
\sum_{i=1}^n p_i\boldsymbol{\alpha}_i=0
\iff \mathbf{A}\mathbf{p}=0.
$$
解空间维数为 $n-\mathrm{rank}(\mathbf{A})=n-r$。取该零空间的一组基 $\mathbf{p}^{(1)},\dots,\mathbf{p}^{(n-r)}$，则
$$
\pi_j=\prod_{i=1}^n q_i^{\,p_i^{(j)}},\qquad j=1,\dots,n-r
$$
为全部独立无量纲量。任何有量纲关系 $f(q_1,\dots,q_n)=0$ 在量纲齐次性下必等价于 $\Phi(\pi_1,\dots,\pi_{n-r})=0$。

物理意义：物理定律不依赖于人为选定的单位，只依赖于无量纲相似参数。几何意义：$\mathbf{A}$ 的零空间即物理规律的“相似子空间”。

特例：单摆 $(n=4:\ l,m,g,T,\ r=3)$，得 $\pi_1=T\sqrt{g/l}$，实验上 $\pi_1=2\pi$。

应用：用 $\mathrm{Re},\ \mathrm{Pe},\ \mathrm{Wi}$ 等数设计血液循环、呼吸道气流与生物反应器。

#### 1.3 扩散时间标度

定义：粒子在无规则热运动下从位置 $x_0$ 出发，其均方位移与时间的关系 $ \langle (x-x_0)^2\rangle = 2Dt$（一维）。

假设：每步位移独立、平稳、等概率左右，步长均方为 $a^2$，每步耗时 $\tau$。

公式：
$$
D=\frac{a^2}{2\tau}\ \text{(一维)},\qquad
t_{\text{diff}}\sim \frac{L^2}{D}.
$$

推导：设第 $k$ 步位移 $s_k$，$\langle s_k\rangle=0$，$\langle s_i s_j\rangle=0\ (i\ne j)$，$\langle s_k^2\rangle=a^2$。则
$$
x_N=\sum_{k=1}^N s_k
\Rightarrow
\langle x_N^2\rangle
=\sum_{i,j}\langle s_i s_j\rangle
=\sum_{k=1}^N a^2 = Na^2.
$$
而 $t=N\tau$，故
$$
\langle x^2\rangle = \frac{a^2}{\tau}\,t \Rightarrow D=\frac{a^2}{2\tau}.
$$
若定义“扩散到尺度 $L$”的时间为 $\langle x^2\rangle=L^2$，则 $t_{\text{diff}}=L^2/(2D)\sim L^2/D$。

物理意义：扩散是“平方根时间”过程，$\sqrt{\langle x^2\rangle}\propto\sqrt{t}$。几何意义：随机游走在空间上以 $t^{1/2}$ 扩展，而非匀速 $t$。

特例：$t\to0$ 时位移为零，符合初始条件；$t\to\infty$ 时均方位移发散，说明无界扩散。

应用：估算氧气在组织内扩散距离、神经元间物质传递、蛋白在细胞质中的到达时间。

#### 1.4 雷诺数与惯性—粘性竞争

定义：$\mathrm{Re}=\dfrac{\rho U L}{\mu}=\dfrac{\text{惯性力}}{\text{粘性力}}$。

假设：流动特征速度 $U$、特征长度 $L$、密度 $\rho$、动力粘度 $\mu$。

公式：由无量纲化 N-S 方程得到惯性项系数 $1/\mathrm{Re}$。

推导：不可压 N-S 方程
$$
\rho(\partial_t\mathbf v+\mathbf v\cdot\nabla\mathbf v)=-\nabla p+\mu\nabla^2\mathbf v.
$$
令 $\mathbf v=U\tilde{\mathbf v},\ \mathbf x=L\tilde{\mathbf x},\ t=(L/U)\tilde t,\ p=\rho U^2\tilde p$，代入并除以 $\rho U^2/L$：
$$
\partial_{\tilde t}\tilde{\mathbf v}+\tilde{\mathbf v}\cdot\nabla_{\tilde x}\tilde{\mathbf v}
=-\nabla_{\tilde x}\tilde p+\frac{\mu}{\rho U L}\nabla^2_{\tilde x}\tilde{\mathbf v}.
$$
故 $\mathrm{Re}=\rho U L/\mu$ 唯一控制解的结构。

物理意义：$\mathrm{Re}\ll1$ 粘性主导（层流、可逆、斯托克斯流）；$\mathrm{Re}\gg1$ 惯性主导（湍流、边界层分离）。几何意义：$\mathrm{Re}$ 是运动方程中非线性对流项与线性粘性扩散项之比。

特例：大肠杆菌游动 $\mathrm{Re}\sim10^{-4}$；主动脉 $\mathrm{Re}\sim10^{3}$；鲸游泳 $\mathrm{Re}\sim10^{8}$。

应用：判定血管、呼吸道、微生物游动与生物反应器搅拌流态。

#### 1.5 布朗运动均方位移与扩散系数

定义：布朗粒子在流体中受随机碰撞，其均方位移 $\langle x^2\rangle$ 随时间线性增长。

假设：过阻尼、无外力、热平衡温度 $T$，随机力满足白噪声。

公式：$\langle x^2\rangle = 2Dt$，$D=\dfrac{k_BT}{6\pi\eta a}$（三维球体）。

推导：过阻尼 Langevin 方程
$$
\gamma\frac{dx}{dt}=\xi(t),\qquad \langle\xi(t)\xi(t')\rangle=2\gamma k_BT\delta(t-t').
$$
积分 $x(t)=x(0)+\frac{1}{\gamma}\int_0^t\xi(t')dt'$，则
$$
\langle x^2\rangle-\langle x\rangle^2
=\frac{1}{\gamma^2}\int_0^t\!\!\int_0^t\langle\xi(t')\xi(t'')\rangle dt'dt''
=\frac{2k_BT}{\gamma}t.
$$
故 $D=k_BT/\gamma$；对半径 $a$ 球体 $\gamma=6\pi\eta a$，得 $D=k_BT/(6\pi\eta a)$。

物理意义：宏观扩散系数由热能与粘性耗散之比决定。几何意义：相空间中的布朗扩散在位置轴上是轴对称的线性增长。

特例：$T\to0$ 时 $D\to0$，无扩散；$a\to\infty$ 时 $D\to0$，大颗粒几乎不扩散。

应用：荧光相关光谱、光镊测量、蛋白扩散系数估算。

#### 1.6 斯托克斯阻力与终端速度

定义：小球在低雷诺数流体中受粘性阻力 $F_d=6\pi\eta a v$。

假设：$\mathrm{Re}\ll1$、不可压牛顿流体、球体、无边界效应。

公式：终端速度 $v_t=\dfrac{2a^2(\rho_p-\rho_f)g}{9\eta}$。

推导：重力与浮力之差 $F_g=\frac{4}{3}\pi a^3(\rho_p-\rho_f)g$，阻力 $F_d=6\pi\eta a v$。匀速时合力为零：
$$
\frac{4}{3}\pi a^3(\rho_p-\rho_f)g=6\pi\eta a v_t
\Rightarrow v_t=\frac{2a^2(\rho_p-\rho_f)g}{9\eta}.
$$
若 $v(0)=0$，由 $m\dot v=F_g-6\pi\eta a v$ 解得
$$
v(t)=v_t\big(1-e^{-t/\tau}\big),\qquad \tau=\frac{m}{6\pi\eta a}=\frac{2\rho_p a^2}{9\eta},
$$
故 $t\gg\tau$ 时趋于终端速度。

物理意义：粘性力与速度线性正比，导致指数趋近终端速度。几何意义：$v_t\propto a^2$ 表明尺寸对沉降速率高度敏感。

特例：$\rho_p<\rho_f$ 时 $v_t<0$，颗粒上浮（脂滴）。$a\to0$ 时 $v_t\to0$，热运动压倒沉降。

应用：血液细胞沉降率、离心分离、微流控颗粒操控。

#### 1.7 代谢标度律（Kleiber 定律）

定义：生物体基础代谢率 $B$ 与体重 $M$ 近似满足 $B\propto M^{3/4}$。

假设：生物体由分形网络输运（West-Brown-Enquist），网络体积 $V\propto M$，末级毛细管尺寸不变。

公式：
$$
B=B_0M^{3/4}.
$$

推导：设网络分 $N$ 级，每级分支比 $n$，末级单元数 $N_c=n^N$，单元质量 $m_c$ 恒定，则
$$
M=N_c m_c\propto n^N.
$$
网络总路径长度为生物体特征规模 $L$，并设第 $k$ 级血管半径 $r_k=r_0\beta^k$、长度 $l_k=l_0\gamma^k$。在“面积守恒”假设下 $n r_k^2=r_{k-1}^2$。若网络体积 $V\propto M$，可写
$$
L=l_0\sum_{k=0}^N \gamma^k\sim l_0\gamma^N,\qquad
M\propto n^N m_c.
$$
单元代谢率 $B_c$ 恒定，总代谢 $B=N_c B_c\propto n^N$。结合 $L\propto\gamma^N$ 与 $M\propto n^N$，并用网络最小化条件（WBE 的 Lagrange 优化）导出
$$
B\propto M^{3/4},\qquad L\propto M^{1/4}.
$$
更简明的标度论证：若每个末级单元供给的质量为常数，而输运路径以分形维度 $d_f$ 展开，则可写 $M\propto L^{d_f}$ 且 $B\propto N_c\propto L^{d_f}$；结合能量耗散最小化给出 $d_f=4$ 的“四维”输运，即 $B\propto M^{3/4}$。

物理意义：代谢率不是与质量成正比（$M^{1}$），而是受输运网络几何的限制呈 $3/4$ 幂律。几何意义：分形输运网络在高维空间有效降维。

特例：单细胞 $M^{3/4}$ 不适用，因其无输运网络；小型恒温动物与大型动物均近似遵守。

应用：药物剂量按体表面积/体重幂律换算、生态能量流估算。

#### 1.8 细胞尺寸的物理限制

定义：细胞大小受扩散、渗透压、细胞骨架与能量供给共同限制。

假设：球形细胞半径 $R$，代谢物扩散系数 $D$，特征消耗速率 $k$。

公式：扩散极限 $R\lesssim\sqrt{D/k}$；渗透压差 $\Delta\pi=RT\Delta c$。

推导：稳态反应-扩散 $\partial c/\partial t=0$ 时
$$
D\nabla^2c-kc=0.
$$
球对称 $c(r)$ 满足 $\frac{D}{r^2}\frac{d}{dr}(r^2c')=kc$，令 $c\sim e^{\lambda r}/r$，得 $\lambda=\sqrt{k/D}$，故浓度随 $r$ 按 $e^{-\sqrt{k/D}\,r}$ 衰减。供给有效尺度即 $\sqrt{D/k}$，大于该尺度中心将缺氧。渗透压由 van't Hoff 公式 $\Pi=cRT$，与膜弹性张力 $\gamma$ 平衡给出 $\Delta p=2\gamma/R$，即
$$
R=\frac{2\gamma}{\Delta c\,RT}.
$$

物理意义：扩散长度决定细胞“不能太大”，渗透压决定细胞“不能太小/太胀”。几何意义：$R=\sqrt{D/k}$ 是浓度场的衰减长度。

特例：$k\to0$ 无消耗时扩散可达任意远；$\Delta c\to0$ 无渗透压差时半径由骨架决定。

应用：解释神经元为何用轴突运输而非纯扩散、红细胞为何无核而小。

### 第2章 生物大分子的热力学基础

#### 2.1 热力学第一定律与焓

定义：内能变化等于吸热与外界做功之和 $\Delta U=Q+W$（化学约定）。

假设：封闭系统、准静态、体积功为主 $W=-p\Delta V$。

公式：
$$
\Delta U=Q-p\Delta V,\qquad H=U+pV,\qquad \Delta H=Q_p.
$$

推导：由定义 $H=U+pV$，取微分
$$
dH=dU+pdV+Vdp.
$$
代入 $dU=\delta Q-p\,dV$ 得
$$
dH=\delta Q+Vdp.
$$
恒压 $dp=0$ 时 $dH=\delta Q$，故等压热等于焓变。

物理意义：生命系统在恒压（大气）下工作，故用焓变描述代谢反应热。几何意义：$H$ 是 $(U,V)$ 到 $(H,p)$ 的 Legendre 变换之一。

特例：理想气体 $H=U+nRT$，等温等压过程 $\Delta H=0$。

应用：量热法测代谢产热、生化反应标准焓 $\Delta H^\circ$。

#### 2.2 熵、第二定律与自由能

定义：熵 $S=k_B\ln\Omega$（Boltzmann），$\Omega$ 为微观状态数。第二定律：孤立系统 $dS\ge0$。

假设：系统与热库接触，温度 $T$ 恒定。

公式：
$$
G=H-TS,\qquad \Delta G=\Delta H-T\Delta S.
$$

推导：Clausius 不等式 $dS\ge \delta Q/T$，可逆过程取等号。恒温恒压自发过程要求
$$
\delta Q\le T\,dS.
$$
由第一定律 $\delta Q=dU+p\,dV$，故
$$
dU+p\,dV-T\,dS\le0.
$$
定义 $G=U+pV-TS$，取微分
$$
dG=dU+p\,dV+V\,dp-T\,dS-S\,dT.
$$
在 $dp=dT=0$ 时 $dG=dU+p\,dV-T\,dS\le0$，即恒温恒压自发过程 $\Delta G\le0$。

物理意义：恒温恒压（活细胞条件）下，反应方向由 $\Delta G$ 符号判定。几何意义：$G$ 是 $U$ 对 $(S,V)\to(T,p)$ 的双重 Legendre 变换（$\Delta$ 形式）。

特例：$\Delta G=0$ 为平衡；$\Delta H<0,\Delta S>0$ 必然自发。

应用：判断 ATP 水解、蛋白折叠、膜形成是否自发。

#### 2.3 化学势与平衡条件

定义：$\mu_i=\left(\dfrac{\partial G}{\partial n_i}\right)_{T,p,n_{j\ne i}}$。

假设：多组分理想稀溶液，组分 $i$ 浓度 $c_i$（活度 $a_i$）。

公式：
$$
\mu_i=\mu_i^\circ+RT\ln a_i,\qquad
\sum_i\nu_i\mu_i=0\ \text{(反应平衡)}.
$$

推导：对理想稀溶液，混合 Gibbs 自由能
$$
G=G^\circ+nRT\sum_i x_i\ln x_i.
$$
对 $n_i$ 求偏导，利用 $\partial(n x_i\ln x_i)/\partial n_i=\ln x_i+1$，得
$$
\mu_i=\mu_i^\circ+RT\ln x_i\approx \mu_i^\circ+RT\ln\frac{c_i}{c^\circ}.
$$
反应 $\sum_i\nu_i A_i=0$ 在平衡时总自由能取极小，约束下 Lagrange 乘子给出 $\sum_i\nu_i\mu_i=0$。对非平衡态，反应亲和势 $\mathcal A=-\sum_i\nu_i\mu_i$ 驱动正向反应。

物理意义：$\mu$ 是物质迁移与化学反应的“驱动力势”。几何意义：$G$-$n$ 曲线斜率即 $\mu$。

特例：$a_i\to0$ 时 $\mu_i\to-\infty$，溶质极稀时化学势极低。

应用：计算膜两侧离子分布、ATP/ADP 比、代谢平衡。

#### 2.4 平衡常数与 van't Hoff 方程

定义：反应 $aA+bB\rightleftharpoons cC+dD$ 的浓度平衡常数 $K=\dfrac{[C]^c[D]^d}{[A]^a[B]^b}$。

假设：理想溶液、恒温恒压、活度以浓度近似。

公式：
$$
\Delta G^\circ=-RT\ln K,\qquad
\frac{d\ln K}{dT}=\frac{\Delta H^\circ}{RT^2}.
$$

推导：平衡条件 $\sum\nu_i\mu_i=0$，
$$
\Delta G^\circ+RT\ln\frac{[C]^c[D]^d}{[A]^a[B]^b}=0
\Rightarrow \Delta G^\circ=-RT\ln K.
$$
对 $T$ 求导，用 Gibbs-Helmholtz 关系 $\Delta G^\circ=\Delta H^\circ-T\Delta S^\circ$ 和 $\partial(\Delta G^\circ/T)/\partial(1/T)=\Delta H^\circ$，得
$$
\frac{d\ln K}{dT}=\frac{\Delta H^\circ}{RT^2}.
$$

物理意义：放热反应升温使 $K$ 减小（Le Chatelier），吸热反应反之。几何意义：$\ln K$ 对 $1/T$ 作图为直线，斜率 $-\Delta H^\circ/R$。

特例：$\Delta H^\circ=0$ 时 $K$ 与温度无关；$\Delta G^\circ=0$ 时 $K=1$。

应用：从量热数据预测生化平衡、酶最适温度分析。

#### 2.5 ATP 水解的自由能与细胞能量货币

定义：$\mathrm{ATP}+\mathrm{H_2O}\rightleftharpoons \mathrm{ADP}+\mathrm{P_i}$，标准自由能 $\Delta G^{\circ\prime}\approx-30.5\ \mathrm{kJ\,mol^{-1}}$（pH 7）。

假设：细胞内非标准态，$[\mathrm{ATP}],[\mathrm{ADP}],[\mathrm{P_i}]$ 已知。

公式：
$$
\Delta G=\Delta G^{\circ\prime}+RT\ln\frac{[\mathrm{ADP}][\mathrm{P_i}]}{[\mathrm{ATP}]}.
$$

推导：由化学势平衡 $\Delta G=\sum\nu_i\mu_i$ 并代入 $\mu_i=\mu_i^{\circ\prime}+RT\ln[\ ]$：
$$
\Delta G=\big(\mu_{\mathrm{ADP}}^\circ+\mu_{\mathrm{P_i}}^\circ-\mu_{\mathrm{ATP}}^\circ-\mu_{\mathrm{H_2O}}^\circ\big)+RT\ln\frac{[\mathrm{ADP}][\mathrm{P_i}]}{[\mathrm{ATP}]},
$$
其中第一项即 $\Delta G^{\circ\prime}$，$[\mathrm{H_2O}]$ 并入常数。

物理意义：$\Delta G$ 是细胞可做功的自由能；其大小随浓度比动态变化。几何意义：在 $([\mathrm{ATP}],[\mathrm{ADP}])$ 平面上 $\Delta G$ 为对数势面。

特例：当 $[\mathrm{ADP}][\mathrm{P_i}]/[\mathrm{ATP}]$ 极大时 $\Delta G$ 可能接近零甚至变正，反应可逆（如 ATP 合成）。

应用：肌肉收缩做功上限、离子泵能量核算、代谢通量控制。

#### 2.6 蛋白质折叠的两态模型

定义：天然态 N 与去折叠态 U 的平衡 $\mathrm N\rightleftharpoons\mathrm U$，折叠自由能 $\Delta G_f=G_U-G_N$。

假设：两态、无稳定中间体、浓度无关。

公式：
$$
\Delta G_f(T)=\Delta H_f-T\Delta S_f,\qquad
K_f=\frac{[U]}{[N]}=e^{-\Delta G_f/RT}.
$$

推导：由平衡常数与自由能关系 $\Delta G_f=-RT\ln K_f$。设热容差 $\Delta C_p$ 恒定，则
$$
\Delta H_f(T)=\Delta H_f(T_m)+\Delta C_p(T-T_m),
$$
$$
\Delta S_f(T)=\Delta S_f(T_m)+\Delta C_p\ln\frac{T}{T_m}.
$$
熔点 $T_m$ 处 $\Delta G_f=0$，故 $\Delta S_f(T_m)=\Delta H_f(T_m)/T_m$。合并得
$$
\Delta G_f(T)=\Delta H_f(T_m)\Big(1-\frac{T}{T_m}\Big)-\Delta C_p\Big[(T_m-T)+T\ln\frac{T}{T_m}\Big].
$$

物理意义：折叠由焓（氢键、范德华）与熵（疏水效应、构象熵）竞争决定。几何意义：$\Delta G_f(T)$ 是开口向上的碗形曲线，$T_m$ 为其零点。

特例：$\Delta C_p\to0$ 时退化为线性 $\Delta G_f=\Delta H_f(1-T/T_m)$；低温变性（冷变性）需 $\Delta C_p>0$ 才能解释。

应用：用化学变性剂滴定曲线拟合 $\Delta G_f$、判断突变稳定性。

#### 2.7 疏水效应与界面自由能

定义：非极性基团在水中的聚集降低体系自由能，称疏水效应。

假设：把非极性分子视为形成水—油界面的球，界面张力 $\gamma$，暴露面积 $A$。

公式：
$$
\Delta G_{\text{hyd}}\approx \gamma\,\Delta A + T\Delta S_{\text{water}}.
$$

推导：把分子从水相移入油相，需破坏“笼形”水结构。单位面积的界面自由能增量近似 $\gamma$。总的 $\Delta G$ 为面积项与水分子熵变项之和。由 $\Delta G=\Delta H-T\Delta S$，疏水水合主要来自 $\Delta S<0$（有序水），故 $\Delta G>0$；当两个疏水面接触，$\Delta A$ 减小使 $\Delta G$ 下降，即 $dG=\gamma\,dA$。

物理意义：疏水效应是蛋白折叠的主要驱动力之一，且以熵为主导。几何意义：$\Delta G$ 正比于溶剂可及表面积变化。

特例：$T$ 升高时 $T\Delta S$ 项增强，疏水效应先增强后减弱。

应用：解释膜自组装、蛋白核心、药物疏水性。

#### 2.8 质量作用定律与反应进度

定义：反应进度 $\xi$，$dn_i=\nu_i\,d\xi$。

假设：理想混合、恒容恒温。

公式：
$$
\frac{d\xi}{dt}=k_+\prod_{j\in\text{react}}c_j^{|\nu_j|}
-k_-\prod_{j\in\text{prod}}c_j^{\nu_j}.
$$

推导：由化学势 $\mu_i=\mu_i^\circ+RT\ln c_i$，反应总自由能
$$
G(\xi)=G^\circ+\sum_i\nu_i\xi\big(\mu_i^\circ+RT\ln c_i\big).
$$
平衡时 $dG/d\xi=0$ 给 $\sum\nu_i\mu_i=0$。远离平衡时，线性不可逆热力学取流 $J\propto\mathcal A=-\sum\nu_i\mu_i$，对小偏离线性化即得正向/逆向速率表达式。对元反应，微观可逆性给出 $k_+/k_-=K$。

物理意义：反应速率的浓度依赖由分子碰撞概率与自由能梯度共同决定。几何意义：$\xi$ 轴上 $G(\xi)$ 的极小点即平衡点。

特例：$k_-=0$ 时不可逆反应；$K=1$ 时正逆速率相等。

应用：酶动力学、代谢网络建模、基因表达。

### 第3章 扩散与输运

#### 3.1 Fick 第一定律的微观推导

定义：扩散通量 $\mathbf J$ 与浓度梯度成正比 $\mathbf J=-D\nabla c$。

假设：粒子独立随机游走，步长 $\ell$，碰撞时间 $\tau$，局部浓度近似线性。

公式：
$$
J_x=-D\frac{\partial c}{\partial x},\qquad D=\frac{\ell^2}{2\tau}.
$$

推导：考虑 $x$ 处截面，每秒从左到右的粒子数为 $\frac12\Gamma c(x-\ell/2)$，从右到左为 $\frac12\Gamma c(x+\ell/2)$，$\Gamma=1/\tau$ 为跳跃频率。净通量
$$
J_x=\frac12\Gamma\big[c(x-\ell/2)-c(x+\ell/2)\big].
$$
Taylor 展开
$$
c(x\pm\ell/2)=c(x)\pm\frac{\ell}{2}c'(x)+\frac{\ell^2}{8}c''(x)+\cdots,
$$
相减得
$$
c(x-\ell/2)-c(x+\ell/2)=-\ell\,c'(x)+O(\ell^3).
$$
故
$$
J_x=-\frac12\Gamma\ell\,c'(x)+O(\ell^3)
=-D\frac{\partial c}{\partial x},
\qquad
D=\frac{\Gamma\ell^2}{2}=\frac{\ell^2}{2\tau}.
$$

物理意义：浓度梯度是扩散的驱动“力”，扩散由随机走动的各向同性导致。几何意义：$\mathbf J$ 与 $-\nabla c$ 平行。

特例：$c$ 均匀时 $J=0$；步长 $\ell\to0$ 但保持 $D$ 有限即连续极限。

应用：膜通透、氧扩散、微流控梯度生成。

#### 3.2 Fick 第二定律（连续性推导）

定义：$\partial_t c=D\nabla^2c$。

假设：无源、无对流、$D$ 常数。

公式：
$$
\frac{\partial c}{\partial t}=D\nabla^2c.
$$

推导：由质量守恒（连续性方程）
$$
\frac{\partial c}{\partial t}+\nabla\cdot\mathbf J=0.
$$
代入 $\mathbf J=-D\nabla c$：
$$
\frac{\partial c}{\partial t}-\nabla\cdot(D\nabla c)=0
\Rightarrow \frac{\partial c}{\partial t}=D\nabla^2c\quad(D\text{ 常数}).
$$

物理意义：浓度的时间变化由曲率（空间不均匀性）驱动。几何意义：$\partial_t c\propto\nabla^2c$，即凹凸处被“抹平”，是典型的扩散平滑方程。

特例：一维稳态 $\nabla^2c=0$ 给线性浓度分布；三维球对称稳态给 $c\propto1/r$（源汇）。

应用：药物释放、氧气输运、组织工程支架设计。

#### 3.3 随机游走与扩散系数（Einstein 关系）

定义：过阻尼粒子在力 $F$ 与热噪声下的漂移-扩散运动。

假设：局部平衡、摩擦系数 $\gamma$、温度 $T$。

公式：
$$
D=\frac{k_BT}{\gamma},\qquad \mathbf J=-D\nabla c+\frac{c}{\gamma}\mathbf F.
$$

推导：粒子在势 $U$ 中，平衡分布 $c\propto e^{-U/k_BT}$，则 $c\mathbf F=-c\nabla U$。平衡时净通量为零：
$$
0=-D\nabla c+\frac{c}{\gamma}\mathbf F
=-D\nabla c-\frac{c}{\gamma}\nabla U.
$$
代入 $\nabla c=-\frac{c}{k_BT}\nabla U$：
$$
0=D\frac{c}{k_BT}\nabla U-\frac{c}{\gamma}\nabla U
\Rightarrow D=\frac{k_BT}{\gamma}.
$$

物理意义：扩散（涨落）与迁移（耗散）由同一摩擦系数联系，是涨落耗散定理的最简形式。几何意义：平衡态通量为零给出 $D$ 与 $\gamma$ 的几何约束。

特例：爱因斯坦关系是 Stokes-Einstein 的基础；带电粒子得 Nernst-Einstein $D=\mu k_BT/q$。

应用：电泳迁移率、离子通道电导、光镊标定。

#### 3.4 反应-扩散方程与 Fisher 波

定义：$\partial_t u=D\nabla^2u+f(u)$，$f(u)$ 为反应项。

假设：一维、$f(u)=ru(1-u/K)$（Logistic 反应）、行波解。

公式：
$$
\frac{\partial u}{\partial t}=D\frac{\partial^2u}{\partial x^2}+ru\Big(1-\frac uK\Big).
$$

推导：令行波 $u(x,t)=U(z)$，$z=x-vt$，则
$$
\frac{\partial u}{\partial t}=-vU',\quad
\frac{\partial^2u}{\partial x^2}=U''.
$$
代入：
$$
D U''+vU'+rU(1-U/K)=0.
$$
在 $U=K$ 前缘 $U\ll K$ 处线性化：$D U''+vU'+rU=0$。设 $U\sim e^{-\lambda z}$，得 $D\lambda^2-v\lambda+r=0$，实根条件 $v^2\ge4Dr$，最小波速
$$
v_{\min}=2\sqrt{Dr}.
$$

物理意义：反应与扩散的竞争形成以有限速度传播的行波。几何意义：在相平面 $(U,U')$ 上，前缘对应异宿轨道，$v_{\min}$ 对应鞍结条件。

特例：$D\to0$ 时 $v_{\min}\to0$，传播受阻；$r\to0$ 同理。

应用：种群扩张、伤口愈合、图灵斑图。

#### 3.5 扩散时间与长度的比较

定义：特征扩散时间 $t_D=L^2/D$；反应时间 $t_R=1/k$。

假设：一级反应、特征尺度 $L$。

公式：
$$
\frac{t_D}{t_R}=\frac{kL^2}{D}=\mathrm{Da}^{-1}\ \text{(Damköhler 数的倒数)}.
$$

推导：无量纲化 $\partial_t c=D\nabla^2c-kc$。令 $\tau=t/t_D,\ \xi=x/L,\ c=c_0\tilde c$：
$$
\frac{1}{t_D}\partial_\tau\tilde c=\frac{D}{L^2}\nabla_\xi^2\tilde c-k\tilde c
\Rightarrow
\partial_\tau\tilde c=\nabla_\xi^2\tilde c-\frac{kL^2}{D}\tilde c.
$$
故无量纲数 $\mathrm{Da}=kL^2/D$（或 $\mathrm{Thiele}$ 模数）唯一控制。

物理意义：$\mathrm{Da}\ll1$ 反应慢，扩散均匀；$\mathrm{Da}\gg1$ 反应快，扩散受限形成梯度。几何意义：$\mathrm{Da}$ 是扩散长度 $\sqrt{D/k}$ 与系统尺度 $L$ 之比的平方。

特例：$\mathrm{Da}=1$ 时二者相当；$\mathrm{Da}\to\infty$ 时扩散极限。

应用：组织氧合、细胞内信号梯度、反应器设计。

#### 3.6 Smoluchowski 扩散限制反应速率

定义：两个分子接近到反应半径 $a$ 即发生反应的速率上限，由扩散供给决定。

假设：稳态、球形吸收边界 $c(a)=0$、远场浓度 $c_\infty$、无对流。

公式：
$$
k=4\pi Da\ (\text{单个靶分子}),\qquad
k_{12}=4\pi(D_1+D_2)(a_1+a_2).
$$

推导：稳态扩散 $\nabla^2c=0$，球对称解
$$
\frac1{r^2}\frac{d}{dr}\Big(r^2\frac{dc}{dr}\Big)=0
\Rightarrow c(r)=c_\infty-\frac{A}{r}.
$$
边界 $c(a)=0$ 给 $A=c_\infty a$。通过半径 $r$ 的球面通量
$$
J=4\pi r^2D\frac{dc}{dr}=4\pi D A=4\pi D a\,c_\infty.
$$
单位浓度下每秒到达靶的分子数即 $k=J/c_\infty=4\pi Da$。两体情形用相对扩散系数 $D_1+D_2$ 与作用距离 $a_1+a_2$ 代替。

物理意义：反应速率不能超过扩散碰撞频率，故酶催化上限由扩散限制。几何意义：$c(r)$ 在靶面为吸收边界、远场趋于 $c_\infty$。

特例：$a\to0$ 或 $D\to0$ 时 $k\to0$；带电吸引可把 $k$ 提高一个数量级以上。

应用：酶-底物结合、蛋白-DNA 搜索、配体-受体结合。

### 第4章 生物膜与跨膜电位

#### 4.1 膜电容模型

定义：细胞膜为脂双层，厚度 $d$，介电常数 $\varepsilon$，构成平行板电容。

假设：膜为均匀绝缘介质、两侧为良导体、面积 $A$。

公式：
$$
C_m=\frac{\varepsilon A}{d}\ \text{(F)},\qquad
Q=C_m V,\qquad I_m=C_m\frac{dV}{dt}.
$$

推导：平行板电容 $C=\varepsilon A/d$。电荷 $Q=CV$，对时间求导得充电电流 $I_m=dQ/dt=C_m\,dV/dt$。

物理意义：膜电位变化需充放电流，决定神经信号时间常数。几何意义：比膜电容 $c_m=C_m/A=\varepsilon/d$，典型 $1\ \mu\mathrm F\,\mathrm{cm^{-2}}$。

特例：稳态 $dV/dt=0$ 时无电容电流；高频时容抗 $1/\omega C$ 减小。

应用：膜片钳、动作电位时间常数 $\tau=R_mC_m$。

#### 4.2 Nernst 方程

定义：单一可通透离子 $i$ 的平衡电位 $E_i$，此时电化学势两侧相等。

假设：膜只对离子 $i$ 通透、理想稀溶液、温度 $T$。

公式：
$$
E_i=\frac{RT}{z_iF}\ln\frac{[i]_{\text{out}}}{[i]_{\text{in}}}.
$$

推导：电化学势
$$
\tilde\mu_i=\mu_i^\circ+RT\ln[i]+z_iF\phi.
$$
平衡时两侧相等：
$$
\mu_i^\circ+RT\ln[i]_{\text{in}}+z_iF\phi_{\text{in}}
=\mu_i^\circ+RT\ln[i]_{\text{out}}+z_iF\phi_{\text{out}}.
$$
整理
$$
z_iF(\phi_{\text{in}}-\phi_{\text{out}})
=RT\ln\frac{[i]_{\text{out}}}{[i]_{\text{in}}}.
$$
定义 $E_i=\phi_{\text{in}}-\phi_{\text{out}}$（膜电位）得
$$
E_i=\frac{RT}{z_iF}\ln\frac{[i]_{\text{out}}}{[i]_{\text{in}}}.
$$

物理意义：离子跨膜浓度梯度与电场力平衡。几何意义：$\ln$ 浓度比与电位差线性对应。

特例：$T=310\ \mathrm K$ 时 $\frac{RT}{F}\ln10\approx61.5\ \mathrm{mV}$，即 $E_i\approx61.5\log_{10}([i]_o/[i]_i)$。

应用：$K^+$、$Na^+$、$Cl^-$ 平衡电位计算，静息电位估算。

#### 4.3 Goldman-Hodgkin-Katz（GHK）电压方程

定义：多离子共同通透时的稳态膜电位。

假设：(1) 膜厚恒定电场；(2) 各离子独立、电扩散服从 Nernst-Planck；(3) 稳态、各离子电流为零。

公式：
$$
V_m=\frac{RT}{F}\ln\frac{P_K[K^+]_o+P_{Na}[Na^+]_o+P_{Cl}[Cl^-]_i}{P_K[K^+]_i+P_{Na}[Na^+]_i+P_{Cl}[Cl^-]_o}.
$$

推导：Nernst-Planck 通量
$$
J_i=-D_i\Big(\frac{dc_i}{dx}+\frac{z_iF}{RT}c_i\frac{d\phi}{dx}\Big).
$$
恒定电场 $d\phi/dx=-V_m/d$。稳态 $J_i$ 恒定，写成
$$
J_i=-D_i\Big(c_i'-\frac{z_iF V_m}{RTd}c_i\Big).
$$
解线性常微分方程，边界 $x=0$ 为内 $c_i(0)=c_i^{in}$，$x=d$ 为外 $c_i(d)=c_i^{out}$，得
$$
J_i=\frac{D_i z_iF V_m}{RTd}\cdot\frac{c_i^{in}-c_i^{out}e^{z_iFV_m/(RT)}}{1-e^{z_iFV_m/(RT)}}.
$$
定义通透系数 $P_i=D_i/d$。总电流 $I=\sum_i z_iFJ_i$。稳态 $I=0$（开路），对一价离子 $K^+,Na^+$（$z=+1$）、$Cl^-$（$z=-1$）整理即得 GHK 公式。

物理意义：膜电位是各离子通透性加权的浓度梯度平均。几何意义：在 $V$ 轴上求总电流零点。

特例：(1) 只有 $K^+$ 通透 $P_K\gg P_{Na},P_{Cl}$，退化为 Nernst $E_K$；(2) 若仅 $Cl^-$ 通透，退化为 $E_{Cl}$ 的等价形式；(3) $V_m=0$ 时无电场项。

应用：静息电位、动作电位平台期、味觉/嗅觉受体电位。

#### 4.4 Donnan 平衡

定义：两室被半透膜隔开，膜允许小离子通透但固定带电荷的大分子不能通透，形成的平衡。

假设：固定电荷浓度 $X$（不可透），可透一价离子 $K^+$、$Cl^-$。

公式：
$$
[K^+]_1[Cl^-]_1=[K^+]_2[Cl^-]_2=K_{eq},
$$
$$
V_{\text{Donnan}}=\frac{RT}{F}\ln\frac{[K^+]_2}{[K^+]_1}.
$$

推导：平衡时 $K^+$ 与 $Cl^-$ 的电化学势分别相等，于是
$$
RT\ln\frac{[K^+]_2}{[K^+]_1}=F(\phi_1-\phi_2),\qquad
RT\ln\frac{[Cl^-]_2}{[Cl^-]_1}=-F(\phi_1-\phi_2).
$$
两式相加消去电位：
$$
\ln\frac{[K^+]_2}{[K^+]_1}+\ln\frac{[Cl^-]_2}{[Cl^-]_1}=0
\Rightarrow [K^+]_1[Cl^-]_1=[K^+]_2[Cl^-]_2.
$$
再由电中性 $[K^+]_1=[Cl^-]_1+X,\ [K^+]_2=[Cl^-]_2$，代入得 $[K^+]_1\big([K^+]_1-X\big)=[K^+]_2^2$。

物理意义：不可透大分子造成跨膜离子不对称与渗透压差。几何意义：Donnan 平衡是电中性平面与等积双曲线交点。

特例：$X=0$ 时两侧浓度相等、电位为零；$X$ 很大时 $[K^+]_1\gg[K^+]_2$，电位近似 $\frac{RT}{F}\ln\frac{X}{[K^+]_2}$。

应用：细胞膨胀压、凝胶离子交换、血浆蛋白与组织液分布。

#### 4.5 电缆方程

定义：沿神经元突起（视为圆柱电缆）膜电位 $V(x,t)$ 的传播方程。

假设：一维、膜为阻容并联、细胞质电阻率 $\rho_i$、膜比电阻 $R_m$、比电容 $C_m$。

公式：
$$
\lambda^2\frac{\partial^2V}{\partial x^2}=\tau_m\frac{\partial V}{\partial t}+V,
$$
$$
\lambda=\sqrt{\frac{R_m d}{4\rho_i}},\qquad \tau_m=R_mC_m.
$$

推导：单位长度膜电流 $i_m=c_m\partial_tV+V/r_m$（$r_m=R_m/(\pi d)$，$c_m=C_m\pi d$）。轴向电流 $i_i=-\frac{1}{r_i}\partial_xV$，单位长度 $r_i=\rho_i/(\pi d^2/4)$。电流守恒：
$$
\partial_x i_i=-i_m.
$$
代入：
$$
\frac{1}{r_i}\partial_x^2V=c_m\partial_tV+\frac{V}{r_m}.
$$
两边乘 $r_m$：
$$
\frac{r_m}{r_i}\partial_x^2V=r_mc_m\partial_tV+V.
$$
定义 $\lambda^2=r_m/r_i=\frac{R_m}{R_i}\frac{d}{4}\cdot$（$R_i=\rho_i$ 单位长定义）与 $\tau_m=R_mC_m$ 即得。

物理意义：$\lambda$ 是电紧张长度（信号衰减到 $1/e$ 的距离），$\tau_m$ 是充放电时间。几何意义：电缆方程是加权扩散-反应方程，$\lambda^2$ 相当于扩散系数，$\tau_m$ 相当于寿命。

特例：稳态 $\partial_tV=0$ 得 $V=V_0e^{-x/\lambda}$；无限长对阶跃电流输入 $V(x,t)$ 为误差函数解。

应用：树突整合、突触电位衰减、动作电位传导速度估算。

#### 4.6 Hodgkin-Huxley 门控变量

定义：$Na^+$、$K^+$ 通道电导由电压依赖门控变量 $m,h,n$ 控制。

假设：门控粒子独立、每个通道需若干粒子同时在位。

公式：
$$
I=C_m\frac{dV}{dt}+\bar g_K n^4(V-E_K)+\bar g_{Na}m^3h(V-E_{Na})+\bar g_l(V-E_l),
$$
$$
\frac{dm}{dt}=\alpha_m(V)(1-m)-\beta_m(V)m,\ \text{等}.
$$

推导：门控变量 $x$ 的一阶动力学可写
$$
\frac{dx}{dt}=\alpha_x(V)(1-x)-\beta_x(V)x
=\frac{x_\infty(V)-x}{\tau_x(V)},
$$
其中 $x_\infty=\alpha_x/(\alpha_x+\beta_x)$，$\tau_x=1/(\alpha_x+\beta_x)$。对 $m^3h$ 形式的通道开放概率，由独立粒子模型 $P_{open}(V)=m^3h$。总膜电流方程由电荷守恒 $I=C_m\dot V+\sum I_{ion}$ 给出。

物理意义：动作电位的产生源于通道电导的电压依赖正反馈（$Na^+$）与失活（$h$）及时延。几何意义：$(m,h,n)$ 三维门控空间的动力学决定 $V$ 的极限环。

特例：$h\to0$ 时 $Na^+$ 完全失活；稳定 $V$ 时门控变量取 $x_\infty$。

应用：神经编码、麻醉机制、心脏节律模型。

### 第5章 分子马达与随机动力学

#### 5.1 Langevin 方程

定义：牛顿方程加随机力的随机微分方程。

假设：粘性摩擦线性、热噪声高斯白噪声、温度 $T$。

公式：
$$
m\frac{dv}{dt}=-\gamma v+F_{\text{ext}}(x)+\xi(t),\qquad
\langle\xi(t)\rangle=0,\quad
\langle\xi(t)\xi(t')\rangle=2\gamma k_BT\delta(t-t').
$$

推导：由涨落耗散定理，系统在无外力时的平衡动能 $\frac12m\langle v^2\rangle=\frac12k_BT$ 要求噪声强度与 $\gamma$ 满足上述关系。验证：对 $F_{\text{ext}}=0$，稳态速度分布为 Boltzmann 分布，可得
$$
\langle v^2\rangle=\frac{k_BT}{m},
$$
而由 Langevin 方程：作 Fourier 变换 $(m i\omega+\gamma)\hat v=\hat\xi$，得功率谱
$$
S_v(\omega)=\frac{2\gamma k_BT}{\gamma^2+m^2\omega^2}.
$$
积分 $\frac{1}{2\pi}\int S_v\,d\omega=\frac{k_BT}{m}$，故噪声强度必须为 $2\gamma k_BT$。

物理意义：随机力与耗散同源，保证能量均分。几何意义：Langevin 轨迹是连续但不可微的布朗路径。

特例：$\gamma\to\infty$ 过阻尼极限（水中小分子）；$t$ 很小时惯性主导，$\langle(x-x_0)^2\rangle\sim v_{th}^2t^2$。

应用：光镊操控、分子马达步进、蛋白构象涨落。

#### 5.2 过阻尼极限与漂移-扩散

定义：当惯性弛豫 $m/\gamma$ 远小于观测时间，忽略惯性项。

假设：$\gamma\,dv/dt\ll\gamma v$，即 $\tau_{in}=m/\gamma\to0$。

公式：
$$
\gamma\frac{dx}{dt}=F(x)+\xi(t),\qquad
\frac{\partial p}{\partial t}=\frac{1}{\gamma}\frac{\partial}{\partial x}\Big[-F(x)p+k_BT\frac{\partial p}{\partial x}\Big].
$$

推导：由 $\gamma v=F+\xi$ 得 $dx/dt=F/\gamma+\xi/\gamma$。写成 Ito 随机微分方程。对应的 Fokker-Planck 方程可用 Kramers-Moyal 展开：对 $\Delta x=\int v\,dt$，
$$
M_1=\frac{F}{\gamma},\qquad
M_2=\frac{2k_BT}{\gamma},
$$
故
$$
\frac{\partial p}{\partial t}=-\frac{\partial}{\partial x}\Big(\frac{F}{\gamma}p\Big)+\frac12\frac{\partial^2}{\partial x^2}\Big(\frac{2k_BT}{\gamma}p\Big)
=-\frac{\partial}{\partial x}\Big(\frac{F p}{\gamma}\Big)+\frac{k_BT}{\gamma}\frac{\partial^2p}{\partial x^2}.
$$
也可写成
$$
\frac{\partial p}{\partial t}=\frac{1}{\gamma}\frac{\partial}{\partial x}\Big(k_BT\frac{\partial p}{\partial x}-Fp\Big).
$$

物理意义：概率流 $J=-\frac{1}{\gamma}\big(k_BT\partial_xp-Fp\big)$。几何意义：概率密度在能量景观中扩散，平衡解 $p\propto e^{-U/k_BT}$。

特例：$F=0$ 时 $\partial_tp=(k_BT/\gamma)\partial_x^2p$，即扩散方程，$D=k_BT/\gamma$。

应用：离子通道、分子马达、蛋白在势阱中的驻留时间。

#### 5.3 涨落耗散定理

定义：系统对微弱外力的响应（耗散）与自发涨落谱相关。

假设：线性响应、平稳、因果。

公式：
$$
S_x(\omega)=\frac{2k_BT}{\omega}\,\mathrm{Im}\,\chi(\omega),
$$
其中 $\chi(\omega)$ 为响应函数。

推导：对线性系统 $\delta x(t)=\int_{-\infty}^t\chi(t-t')f(t')dt'$，其功率谱 $S_x=|\chi(\omega)|^2S_f$。由因果性与平衡条件，涨落谱与响应函数满足上述关系。对粒子位置：$\chi(\omega)=1/(\gamma+m i\omega)$，代入得
$$
S_x(\omega)=\frac{2k_BT}{\omega}\mathrm{Im}\frac{1}{\gamma+im\omega}
=\frac{2k_BT}{\omega}\cdot\frac{m\omega}{\gamma^2+m^2\omega^2}
=\frac{2k_BTm}{\gamma^2+m^2\omega^2}.
$$

物理意义：耗散越强（$\gamma$ 越大），涨落越强，体系能更快回到平衡。几何意义：$\mathrm{Im}\chi$ 的峰位对应系统的耗散共振。

特例：谐振子（光镊）的位移谱在 $\omega_0=\sqrt{k/m}$ 处共振；低频 $S_x\approx2k_BTm/\gamma^2$。

应用：光镊测刚度、原子力显微镜、膜张力测量。

#### 5.4 Fokker-Planck 方程（Kramers-Moyal 展开）

定义：概率密度 $p(x,t)$ 在漂移 $A(x)$ 与扩散 $B(x)$ 下的演化方程。

假设：Markov、连续路径、跳跃短时。

公式：
$$
\frac{\partial p}{\partial t}
=-\frac{\partial}{\partial x}\big[A(x)p\big]+\frac12\frac{\partial^2}{\partial x^2}\big[B(x)p\big].
$$

推导：Chapman-Kolmogorov 方程
$$
p(x,t+\Delta t)=\int p(x-\Delta x,t)\,\phi(\Delta x|x-\Delta x)\,d\Delta x.
$$
对 $\Delta x$ 作 Taylor 且记矩
$$
A(x)=\lim_{\Delta t\to0}\frac{\langle\Delta x\rangle}{\Delta t},\qquad
B(x)=\lim_{\Delta t\to0}\frac{\langle(\Delta x)^2\rangle}{\Delta t}.
$$
保留至二阶：
$$
\frac{\partial p}{\partial t}=-\frac{\partial}{\partial x}(Ap)+\frac12\frac{\partial^2}{\partial x^2}(Bp).
$$

物理意义：漂移与扩散是随机过程的前两阶矩。几何意义：$A$ 驱动概率中心移动，$B$ 使分布展宽。

特例：$B=2D$ 常数、$A=-kx/\gamma$（线性回复力）得 Ornstein-Uhlenbeck 过程，稳态高斯分布。

应用：通道动力学、基因表达噪声、群体遗传漂变。

#### 5.5 Kramers 逃逸率

定义：粒子在势阱中被热噪声驱动的逃逸速率。

假设：高势垒 $\Delta U\gg k_BT$、过阻尼、势阱底部近似抛物。

公式：
$$
k\approx\frac{\omega_0\omega_b}{2\pi\xi}\,e^{-\Delta U/k_BT},
$$
其中 $\omega_0$、$\omega_b$ 为阱底与势垒顶的曲率频率。

推导：在过阻尼 Fokker-Planck 方程中，把势 $U(x)$ 在阱底 $x_a$ 与垒顶 $x_b$ 附近作二次近似：
$$
U\approx U_a+\frac12m\omega_0^2(x-x_a)^2,\qquad
U\approx U_b-\frac12m\omega_b^2(x-x_b)^2.
$$
求准稳态概率流，得 Kramers 公式。核心指数因子来自 Boltzmann 因子 $e^{-\Delta U/k_BT}$，前因子来自两个曲率的比。

物理意义：逃逸时间对势垒高度指数敏感，也称 Arrhenius 律。几何意义：$\Delta U$ 是能量景观纵差，$\omega$ 是阱的横向曲率。

特例：$\Delta U\gg k_BT$ 成立；$\Delta U\sim k_BT$ 时公式失效。过阻尼与低阻尼前因子不同。

应用：蛋白去折叠、配体解离、化学键断裂力谱。

#### 5.6 分子马达效率与失速力

定义：分子马达把化学能（ATP 水解）转化为机械功的效率。

假设：每步消耗 $N$ 个 ATP、步长 $d$、负载力 $F$。

公式：
$$
\eta=\frac{W}{\Delta G}=\frac{F d}{N|\Delta G_{\text{ATP}}|},\qquad
F_{\text{stall}}=\frac{N|\Delta G_{\text{ATP}}|}{d}.
$$

推导：最大功 $W_{\max}=F d$（每步克服负载做功）。输入自由能 $N|\Delta G_{\text{ATP}}|$。由上章的 ATP 水解自由能，令效率不大于 1 得
$$
F\le\frac{N|\Delta G_{\text{ATP}}|}{d}.
$$
取等号即失速力 $F_{\text{stall}}$。考虑 \(\Delta G_{\text{ATP}}\) 的浓度依赖，还可写出
$$
F_{\text{stall}}=\frac{N}{d}\Big(\Delta G^{\circ\prime}+RT\ln\frac{[\mathrm{ADP}][\mathrm{P_i}]}{[\mathrm{ATP}]}\Big).
$$

物理意义：马达在负载增大到失速力时停止。几何意义：在力-速度平面，效率在中等负载取得最大值。

特例：$F\to0$ 时效率趋于零（怠速空耗）；$F=F_{\text{stall}}$ 时速度为零、效率为零。

应用：驱动蛋白运输、肌球蛋白收缩、ATP 合成酶反转。

### 第6章 生物聚合物力学

#### 6.1 自由连接链的构象熵

定义：由 $N$ 个刚性链段（Kuhn 长度 $b$）自由连接构成的柔性链。

假设：链段取向完全独立、无体积排斥、无外力。

公式：
$$
\Omega(\{n_i\})=\frac{N!}{\prod_i n_i!},\qquad
S=k_B\ln\Omega.
$$

推导：设链段有 $m$ 个等概率方向。$N$ 段分布 $n_1,\dots,n_m$ 的微观态数由多项式系数给出。对一维简化模型（每段只能 $\pm b$），选 $n_+$ 段正向、$n_-=N-n_+$ 段负向，端距 $x=(n_+-n_-)b$，则
$$
\Omega=\frac{N!}{n_+!\,n_-!}.
$$
利用 Stirling 公式 $\ln n!\approx n\ln n-n$，
$$
\ln\Omega\approx N\ln N-n_+\ln n_+-n_-\ln n_-.
$$
设 $n_\pm=(N/2)(1\pm x/(Nb))$，则
$$
S(x)=k_B\Big[\ln2-\frac{1}{2}\Big(1+\frac{x}{Nb}\Big)\ln\Big(1+\frac{x}{Nb}\Big)
-\frac{1}{2}\Big(1-\frac{x}{Nb}\Big)\ln\Big(1-\frac{x}{Nb}\Big)\Big]N.
$$

物理意义：链的每一个宏观伸长都对应更少的微观构象，故伸长降低熵。几何意义：$\Omega(x)$ 在 $x=0$ 最大，随 $|x|$ 增大而锐减。

特例：$x\ll Nb$ 时高斯近似 $S\approx S_0-\frac{3k_Bx^2}{2Nb^2}$；$x\to Nb$ 时 $S\to0$。

应用：DNA 拉伸、橡胶弹性、染色质构象。

#### 6.2 熵力的推导

定义：$f=-T\left(\dfrac{\partial S}{\partial x}\right)_U$。

假设：内能 $U$ 与构象无关（理想链）。

公式：
$$
f=\frac{k_BT}{b}\left[\ln\frac{1+x/(Nb)}{1-x/(Nb)}\right]^{-1}\times\cdots
$$
正确表达式由 $x=Nb\,\mathcal L(fb/k_BT)$ 的逆给出。高斯区：
$$
f=\frac{3k_BT}{Nb^2}\,x,
$$
近完全拉伸：
$$
f\approx\frac{k_BT}{2b}\ln\frac{Nb}{Nb-x}.
$$

推导：由 $f=-T\partial S/\partial x$，对一维表达式
$$
S=\frac{k_BN}{2}\Big[-2\ln2+\cdots\Big]
$$
求导。设 $u=x/(Nb)$，
$$
\frac{1}{k_B}\frac{\partial S}{\partial x}
=-\frac{1}{2b}\Big[\ln(1+u)-\ln(1-u)\Big]
+\text{(常数项导数抵消)},
$$
故
$$
f=-T\frac{\partial S}{\partial x}
=\frac{k_BT}{2b}\ln\frac{1+u}{1-u}.
$$
反解 $u=\frac{e^{2fb/(k_BT)}-1}{e^{2fb/(k_BT)}+1}=\tanh\!\Big(\frac{fb}{k_BT}\Big)$，即
$$
x=Nb\tanh\frac{fb}{k_BT}\ \Rightarrow\ f=\frac{k_BT}{b}\,\mathrm{artanh}\frac{x}{Nb}.
$$
高斯区 $\mathrm{artanh}u\approx u$ 给 $f=3k_BTx/(Nb^2)$（用三维对应系数）。

物理意义：拉力来源是构象熵减少而非键能拉伸。几何意义：$f(x)$ 是单调凸增函数，$x\to Nb$ 时发散。

特例：$x/Nb\ll1$ 线性（高斯链）；$x/Nb\to1$ 对数发散，力趋于无穷。

应用：DNA 镊实验、橡胶状态方程、无规卷曲蛋白。

#### 6.3 蠕虫链模型（WLC）与持续性长度

定义：WLC 用连续弯曲能描述半柔性链。

假设：链的弯曲刚度为 $\kappa$、轮廓长度 $L_c$、热扰动产生横向涨落。

公式：
$$
E=\frac{\kappa}{2}\int_0^{L_c}\Big(\frac{\partial\mathbf t}{\partial s}\Big)^2ds,\qquad
L_p=\frac{\kappa}{k_BT},\qquad
\langle\mathbf t(s)\cdot\mathbf t(0)\rangle=e^{-s/L_p}.
$$

推导：定向相关函数 $C(s)=\langle\cos\theta(s)\rangle$。在二维简化下，弯曲角 $\theta(s)$ 满足 $d\theta/ds=\tau(s)$（曲率），能量 $E=\frac\kappa2\int\tau^2ds$。热力学平均
$$
C(s)=e^{-s/(2L_p)}\ \text{(二维)},\quad C(s)=e^{-s/L_p}\ \text{(三维)}.
$$
更严格地，由 $C(s)=\exp\big[-\frac{k_BT}{2\kappa}s\big]\cdot$（维度因子），定义 $L_p=\kappa/k_BT$ 使得 $C(s)=e^{-s/L_p}$（三维）。Marko-Siggia 插值力-伸公式：
$$
f=\frac{k_BT}{L_p}\Big[\frac{1}{4(1-x/L_c)^2}-\frac14+\frac{x}{L_c}\Big].
$$

物理意义：$L_p$ 是链保持方向记忆的特征长度。几何意义：$L_p$ 越大，链越刚，弯曲能越高。

特例：$L_p\to0$ 得柔性高斯链；$L_p\to\infty$ 得刚性杆。$x\ll L_c$ 时 $f\approx\frac{3k_BTx}{2L_pL_c}$，与高斯链一致；$x\to L_c$ 时发散。

应用：双链 DNA（$L_p\approx50\ \mathrm{nm}$）、微管、肌动蛋白。

#### 6.4 橡胶弹性

定义：交联高分子网络的熵弹性。

假设：交联点间链为高斯链、仿射变形（交联点随宏观应变移动）。

公式：
$$
W=\frac{1}{2}n k_BT(\lambda_1^2+\lambda_2^2+\lambda_3^2-3),\qquad
\sigma=n k_BT\Big(\lambda-\frac{1}{\lambda^2}\Big).
$$

推导：单链高斯熵 $S=-3k_Bx^2/(2Nb^2)$，故自由能 $A=-TS=\frac{3k_BT}{2Nb^2}x^2$。单位体积 $n$ 条链，仿射变形 $\mathbf x\to\mathbf F\mathbf x$，则
$$
dA=\frac{3k_BT}{Nb^2}\mathbf x\cdot\delta(\mathbf F\mathbf x)
=\frac{3k_BT}{Nb^2}(\mathbf F^T\mathbf F\mathbf x)\cdot\delta\mathbf x.
$$
对全部链取平均，各向同性初态 $\langle\mathbf x\mathbf x\rangle=\frac13 Nb^2\mathbf I$，得
$$
\frac{\partial A}{\partial\mathbf F}=n k_BT\,\mathbf F,
$$
在不可压 $\lambda_1\lambda_2\lambda_3=1$ 下，单轴拉伸 $\lambda_1=\lambda,\lambda_2=\lambda_3=1/\sqrt\lambda$，工程应力 $\sigma=\partial W/\partial\lambda$（扣去静水压）得
$$
\sigma=n k_BT\Big(\lambda-\frac{1}{\lambda^2}\Big).
$$

物理意义：橡胶弹性几乎纯熵弹性，$\sigma\propto T$。几何意义：熵自由能是应变的各向同性二次型。

特例：$\lambda\to1$ 时 $\sigma\approx3nk_BT(\lambda-1)$ 线性；大拉伸时偏离高斯而硬化（有限伸长）。

应用：弹性蛋白、水凝胶、组织工程支架。

#### 6.5 蛋白质去折叠的 WLC 拟合

定义：用 WLC 力-伸关系拟合单分子力谱中的去折叠台阶。

假设：去折叠态为 WLC、折叠态为刚性核、拉伸准静态。

公式：
$$
f=\frac{k_BT}{L_p}\Big[\frac{1}{4(1-x/L_c)^2}-\frac14+\frac{x}{L_c}\Big].
$$

推导：每次去折叠释放轮廓长度 $\Delta L_c$。把力-伸曲线在恒定力平台处作 WLC 拟合，得 $L_p$ 与 $L_c$。忽略 Marko-Siggia 式中的线性项 $x/L_c$，令 $u=x/L_c$、$\phi=fL_p/(k_BT)$，则
$$
\frac{1}{4(1-u)^2}-\frac14=\phi
\ \Rightarrow\ \frac{1}{4(1-u)^2}=\phi+\frac14
\ \Rightarrow\ 1-u=\frac{1}{\sqrt{1+4\phi}}.
$$
故反解
$$
x=L_c\Big(1-\frac{1}{\sqrt{1+4\phi}}\Big),
$$
在大力 $\phi\gg1$ 时
$$
x\approx L_c\Big(1-\frac{1}{2\sqrt{\phi}}\Big)
=L_c\Big(1-\frac12\sqrt{\frac{k_BT}{L_p f}}\Big).
$$
台阶间距 $\Delta x$ 对应释放的 $\Delta L_c$。

物理意义：力谱平台反映结构域逐个展开。几何意义：$x$-$f$ 曲线上平台与后续 WLC 上升段依次出现。

特例：$f\gg k_BT/L_p$ 时 $x\approx L_c\big(1-\frac12\sqrt{k_BT/(L_pf)}\big)$；小力时接近高斯链。

应用：titin、GFP、DNA 发夹去折叠。

#### 6.6 聚合物网络的网孔尺寸与标度

定义：半稀溶液中高分子链相互缠结，网孔（blob）尺寸 $\xi$。

假设：半稀、无规线团、blob 模型、排除体积指数 $\nu$（良溶剂 $3/5$，$\Theta$ 溶剂 $1/2$）。

公式：
$$
\xi\sim R_F\Big(\frac{c}{c^*}\Big)^{-\nu/(3\nu-1)},
$$
良溶剂 $\xi\propto c^{-3/4}$，$\Theta$ 溶剂 $\xi\propto c^{-1}$。

推导：每个 blob 含 $g$ 个单体，尺寸 $\xi\sim b\,g^{\nu}$；半稀条件 $c\sim g/\xi^3$（体积分数）。于是
$$
c\sim \frac{g}{\xi^3}\sim \frac{(\xi/b)^{1/\nu}}{\xi^3}
=b^{-1/\nu}\xi^{\,1/\nu-3}.
$$
解出
$$
\xi\propto c^{1/(1/\nu-3)}=c^{\nu/(1-3\nu)}=c^{-\nu/(3\nu-1)}.
$$
良溶剂 $\nu=3/5$：$\xi\propto c^{-3/4}$；$\Theta$ 溶剂 $\nu=1/2$：$\xi\propto c^{-1}$。网络剪切模量 $G\sim k_BT/\xi^3$。

物理意义：网孔尺寸由浓度通过标度律决定，并控制扩散与力学。几何意义：blob 是链在 $\xi$ 尺度内近似理想的单元。

特例：$c\to c^*$ 时 $\xi\to R_F$；$c\to\infty$ 时 $\xi\to0$（浓溶液）。

应用：细胞外基质网孔、凝胶筛选、黏液渗透性。

### 第7章 细胞与组织的连续介质力学

#### 7.1 小应变张量

定义：位移场 $\mathbf u(\mathbf x)$，小应变张量
$$
\varepsilon_{ij}=\frac12(\partial_i u_j+\partial_j u_i).
$$

假设：位移梯度 $\|\nabla\mathbf u\|\ll1$。

公式：$\varepsilon_{ij}=\frac12(\partial_i u_j+\partial_j u_i)$。

推导：两点 $\mathbf x$、$\mathbf x+d\mathbf x$ 变形后距离平方
$$
d\ell'^2=|\mathbf dx+d\mathbf u|^2
=d\ell^2+2\,du_i\,dx_i+du_i\,du_i.
$$
其中 $du_i=\partial_j u_i\,dx_j$。忽略二阶 $du_idu_i$：
$$
d\ell'^2-d\ell^2=2\partial_j u_i\,dx_jdx_i
=(\partial_j u_i+\partial_i u_j)dx_idx_j
=2\varepsilon_{ij}dx_idx_j.
$$
故 $\varepsilon_{ij}$ 对称。

物理意义：应变描述线元长度与角度的相对变化。几何意义：对称部分为纯变形，反对称部分为刚体旋转。

特例：纯平移 $\mathbf u=\text{const}$ 时 $\varepsilon=0$；纯旋转 $\mathbf u=\boldsymbol\omega\times\mathbf x$ 时对称部分为零。

应用：组织拉伸测试、超声弹性成像、伤口收缩。

#### 7.2 应力张量与平衡方程

定义：Cauchy 应力张量 $\boldsymbol\sigma$，面上牵引力 $\mathbf t=\boldsymbol\sigma\mathbf n$。

假设：连续介质、角动量守恒（无体力偶）。

公式：
$$
\nabla\cdot\boldsymbol\sigma+\mathbf f=0\ \text{(静态)},\qquad
\rho\frac{D\mathbf v}{Dt}=\nabla\cdot\boldsymbol\sigma+\mathbf f.
$$

推导：对体积 $V$ 作动量守恒：
$$
\frac{d}{dt}\int_V\rho v_i\,dV=\int_{\partial V}\sigma_{ij}n_j\,dA+\int_V f_i\,dV.
$$
对面积分用散度定理：
$$
\frac{d}{dt}\int_V\rho v_i\,dV=\int_V(\partial_j\sigma_{ij}+f_i)\,dV.
$$
对任意 $V$ 成立，得被积函数相等：
$$
\rho\frac{Dv_i}{Dt}=\partial_j\sigma_{ij}+f_i.
$$
静态即 $\nabla\cdot\boldsymbol\sigma+\mathbf f=0$。角动量守恒给 $\sigma_{ij}=\sigma_{ji}$。

物理意义：应力是面力密度，其散度是单位体积净力。几何意义：应力张量的特征值为主应力，特征向量为主方向。

特例：静水压 $\boldsymbol\sigma=-p\mathbf I$ 时 $\nabla\cdot\boldsymbol\sigma=-\nabla p$。

应用：组织内应力分布、血管壁张力（Laplace 定律）、细胞牵引力。

#### 7.3 各向同性线弹性本构（Hooke 定律）

定义：应力与应变的线性关系。

假设：小应变、各向同性、等温、线性弹性。

公式：
$$
\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij},
$$
$$
\lambda=\frac{E\nu}{(1+\nu)(1-2\nu)},\qquad \mu=G=\frac{E}{2(1+\nu)}.
$$

推导：各向同性四阶张量只能由 $\delta_{ij}\delta_{kl}$ 与 $\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}$ 线性组合：
$$
C_{ijkl}=\lambda\delta_{ij}\delta_{kl}+\mu(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}).
$$
代入 $\sigma_{ij}=C_{ijkl}\varepsilon_{kl}$ 得
$$
\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij}.
$$
由 $E$、$\nu$ 反解 $\lambda,\mu$：单轴应力 $\sigma_{11}=\sigma$，其他 $\sigma_{ij}=0$。由 $\sigma_{22}=0$ 得 $\varepsilon_{22}=\varepsilon_{33}=-\nu\varepsilon_{11}$，其中 $\nu=\lambda/[2(\lambda+\mu)]$；$\sigma_{11}=E\varepsilon_{11}$，$E=\mu(3\lambda+2\mu)/(\lambda+\mu)$。反解即得 $\lambda,\mu$ 表达式。

物理意义：$\lambda,\mu$ 为 Lamé 常数，$\mu$ 为剪切模量，$E$ 为杨氏模量，$\nu$ 为泊松比。几何意义：体积模量 $K=\lambda+\frac23\mu$ 度量抗压缩。

特例：不可压 $\nu=0.5$ 时 $K\to\infty$；$\nu=0$ 时无横向收缩。

应用：骨/软骨/血管弹性模量、原子力显微镜压痕、超声速度 $c=\sqrt{E/\rho}$。

#### 7.4 弹性平衡方程（Navier-Cauchy / Lamé）

定义：位移场满足的平衡方程。

假设：各向同性线弹性、小应变、无体力或常体力。

公式：
$$
(\lambda+\mu)\nabla(\nabla\cdot\mathbf u)+\mu\nabla^2\mathbf u+\mathbf f=0.
$$

推导：由 $\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij}$，代入平衡 $\partial_j\sigma_{ij}+f_i=0$：
$$
\partial_j\big(\lambda\partial_k u_k\delta_{ij}+2\mu\varepsilon_{ij}\big)+f_i=0.
$$
用 $\partial_j\varepsilon_{ij}=\frac12(\partial_j\partial_i u_j+\partial_j\partial_j u_i)=\frac12\partial_i(\nabla\cdot\mathbf u)+\frac12\nabla^2u_i$，得
$$
\lambda\partial_i(\nabla\cdot\mathbf u)+\mu\partial_i(\nabla\cdot\mathbf u)+\mu\nabla^2u_i+f_i=0,
$$
即向量形式
$$
(\lambda+\mu)\nabla(\nabla\cdot\mathbf u)+\mu\nabla^2\mathbf u+\mathbf f=0.
$$

物理意义：弹性体平衡是位移的椭圆型方程。几何意义：位移场的旋度部分与散度部分解耦。

特例：$\nabla\cdot\mathbf u=0$（不可压）时退化为 $\mu\nabla^2\mathbf u+\mathbf f=0$，即 Stokes 流。

应用：组织变形有限元、血管壁应力、细胞外基质力学。

#### 7.5 Maxwell 粘弹性模型

定义：弹簧与阻尼串联的粘弹性模型。

假设：总应变率 = 弹性 + 粘性，应力在二者中相等。

公式：
$$
\dot\varepsilon=\frac{\dot\sigma}{E}+\frac{\sigma}{\eta}
\ \Longleftrightarrow\
\sigma+\tau\dot\sigma=\eta\dot\varepsilon,\quad \tau=\frac\eta E.
$$

推导：弹簧 $\varepsilon_e=\sigma/E$，阻尼 $\dot\varepsilon_v=\sigma/\eta$。串联 $\varepsilon=\varepsilon_e+\varepsilon_v$，求导
$$
\dot\varepsilon=\dot\varepsilon_e+\dot\varepsilon_v=\frac{\dot\sigma}{E}+\frac{\sigma}{\eta}.
$$
两边乘 $\eta$ 得 $\sigma+\tau\dot\sigma=\eta\dot\varepsilon$。

应力松弛：$\dot\varepsilon=0$，方程 $\sigma+\tau\dot\sigma=0$，解
$$
\sigma(t)=\sigma_0e^{-t/\tau}.
$$
蠕变：应力阶跃 $\sigma_0$，$\dot\varepsilon=\sigma_0/\eta$（弹性瞬时应变后线性增长）。

物理意义：Maxwell 描述应力松弛与流动；有瞬时弹性、无瞬时粘性。几何意义：在 $\sigma$-$\varepsilon$ 平面轨迹是弹簧线+水平流。

特例：$\eta\to\infty$ 纯弹性；$E\to\infty$ 纯粘性。稳态 $\dot\sigma=0$ 时 $\sigma=\eta\dot\varepsilon$。

应用：血栓、黏液、软组织应力松弛。

#### 7.6 Kelvin-Voigt 模型与标准线性固体

定义：弹簧与阻尼并联（Kelvin-Voigt）；三参数标准线性固体为两者组合。

假设：应变相同（并联），或弹簧 + Maxwell 串联。

公式（Kelvin-Voigt）：
$$
\sigma=E\varepsilon+\eta\dot\varepsilon.
$$
标准线性固体：
$$
\sigma+\tau_\sigma\dot\sigma=E_0(\varepsilon+\tau_\varepsilon\dot\varepsilon).
$$

推导 (KV)：并联故 $\sigma=\sigma_e+\sigma_v=E\varepsilon+\eta\dot\varepsilon$。蠕变阶跃 $\sigma_0$：
$$
E\varepsilon+\eta\dot\varepsilon=\sigma_0
\Rightarrow \varepsilon(t)=\frac{\sigma_0}{E}\big(1-e^{-t/\tau}\big),\ \tau=\frac\eta E.
$$
标准线性固体：弹簧 $E_1$ 与 Maxwell（$E_2,\eta$）并联：
$$
\sigma=E_1\varepsilon+\sigma_M,\quad \dot\varepsilon=\frac{\dot\sigma_M}{E_2}+\frac{\sigma_M}{\eta}.
$$
消去 $\sigma_M$ 得 $\sigma+\tau_\sigma\dot\sigma=E_1\varepsilon+\eta(1+E_1/E_2)\dot\varepsilon$，即 $\tau_\sigma=\eta/E_2$，$\tau_\varepsilon=\tau_\sigma(1+E_1/E_2)$。

物理意义：KV 有延迟弹性但无应力松弛；SLS 同时有松弛与蠕变。几何意义：SLS 的松弛函数为 $E_1+E_2e^{-t/\tau}$。

特例：KV 中 $\eta\to0$ 纯弹性；$\eta\to\infty$ 时 $\dot\varepsilon=0$ 刚体。SLS 中 $E_1=0$ 退化为 Maxwell。

应用：软骨、细胞骨架、水凝胶频率响应。

#### 7.7 复模量与储能/损耗模量

定义：正弦激励下 $\sigma^*=G^*\varepsilon^*$，$G^*=G'+iG''$。

假设：线性粘弹性、稳态、小应变。

公式：
$$
G'(\omega)=G_0\frac{\omega^2\tau^2}{1+\omega^2\tau^2},\qquad
G''(\omega)=G_0\frac{\omega\tau}{1+\omega^2\tau^2}
$$
（标准线性固体的典型形式）。

推导：设 $\varepsilon(t)=\varepsilon_0e^{i\omega t}$，$\sigma(t)=\sigma_0e^{i\omega t}$。代入 SLS：
$$
\sigma_0(1+i\omega\tau_\sigma)=E_0\varepsilon_0(1+i\omega\tau_\varepsilon).
$$
故
$$
G^*=\frac{\sigma_0}{\varepsilon_0}=E_0\frac{1+i\omega\tau_\varepsilon}{1+i\omega\tau_\sigma}.
$$
有理化：
$$
G^*=E_0\frac{(1+i\omega\tau_\varepsilon)(1-i\omega\tau_\sigma)}{1+\omega^2\tau_\sigma^2}
=E_0\frac{1+\omega^2\tau_\varepsilon\tau_\sigma+i\omega(\tau_\varepsilon-\tau_\sigma)}{1+\omega^2\tau_\sigma^2}.
$$
故
$$
G'=E_0\frac{1+\omega^2\tau_\varepsilon\tau_\sigma}{1+\omega^2\tau_\sigma^2},\qquad
G''=E_0\frac{\omega(\tau_\varepsilon-\tau_\sigma)}{1+\omega^2\tau_\sigma^2}.
$$

物理意义：$G'$ 储能（弹性），$G''$ 耗能（粘性），损耗角 $\tan\delta=G''/G'$。几何意义：在复平面 $G^*$ 轨迹为半圆（Cole-Cole）。

特例：$\omega\to0$ 时 $G'\to E_0$（松弛模量）；$\omega\to\infty$ 时 $G'\to E_0\tau_\varepsilon/\tau_\sigma$（瞬时模量）。

应用：流变仪测量细胞/组织粘弹性、微流变。

### 第8章 生物流体力学

#### 8.1 连续性方程

定义：质量守恒的微分形式。

假设：连续介质。

公式：
$$
\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf v)=0,\qquad
\nabla\cdot\mathbf v=0\ (\text{不可压}).
$$

推导：对任意控制体 $V$，
$$
\frac{d}{dt}\int_V\rho\,dV=-\oint_{\partial V}\rho\mathbf v\cdot\mathbf n\,dA.
$$
用 Reynolds 输运定理由固定体积：
$$
\int_V\frac{\partial\rho}{\partial t}dV=-\int_V\nabla\cdot(\rho\mathbf v)\,dV.
$$
对任意 $V$ 成立得
$$
\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf v)=0.
$$
不可压 $\rho$ 常数时得 $\nabla\cdot\mathbf v=0$。

物理意义：质量不生不灭。几何意义：速度场无源无汇（不可压），流线不可压缩。

特例：一维稳态变截面管 $\rho A v=\text{const}$（连续方程积分形式）。

应用：血管分叉、呼吸道流量守恒、心脏瓣膜流量。

#### 8.2 Navier-Stokes 方程推导

定义：牛顿流体动量方程。

假设：牛顿本构 $\boldsymbol\tau=\mu(\nabla\mathbf v+\nabla\mathbf v^T)-\frac23\mu(\nabla\cdot\mathbf v)\mathbf I$，常粘度。

公式：
$$
\rho\Big(\frac{\partial\mathbf v}{\partial t}+\mathbf v\cdot\nabla\mathbf v\Big)=-\nabla p+\mu\nabla^2\mathbf v+\mathbf f.
$$

推导：动量守恒
$$
\rho\frac{D\mathbf v}{Dt}=\nabla\cdot\boldsymbol\sigma+\mathbf f,\qquad
\boldsymbol\sigma=-p\mathbf I+\boldsymbol\tau.
$$
对不可压流体 $\nabla\cdot\mathbf v=0$，牛顿本构 $\boldsymbol\tau=\mu(\nabla\mathbf v+\nabla\mathbf v^T)$。于是
$$
(\nabla\cdot\boldsymbol\tau)_i=\mu\partial_j(\partial_jv_i+\partial_iv_j)
=\mu\nabla^2v_i+\mu\partial_i(\nabla\cdot\mathbf v)
=\mu\nabla^2v_i.
$$
代入得 N-S 方程。$D\mathbf v/Dt=\partial_t\mathbf v+\mathbf v\cdot\nabla\mathbf v$ 为物质导数。

物理意义：惯性、压力、粘性与外力平衡。几何意义：非线性对流项使方程非线、可能湍流。

特例：$\mathrm{Re}\to0$ 得 Stokes 方程 $-\nabla p+\mu\nabla^2\mathbf v=0$；无粘得 Euler 方程。

应用：动脉血流、微循环、精子游动、生物反应器。

#### 8.3 Hagen-Poiseuille 流动

定义：圆管中定常层流不可压牛顿流。

假设：轴对称、充分发展（$\partial_zv_z=0$）、无滑移、定常。

公式：
$$
v_z(r)=\frac{\Delta p}{4\mu L}(R^2-r^2),\qquad
Q=\frac{\pi R^4\Delta p}{8\mu L}.
$$

推导：N-S 的 $z$ 分量在充分发展、轴对称下
$$
0=-\frac{dp}{dz}+\mu\frac{1}{r}\frac{d}{dr}\Big(r\frac{dv_z}{dr}\Big).
$$
积分：
$$
\frac{d}{dr}\Big(r\frac{dv_z}{dr}\Big)=\frac{r}{\mu}\frac{dp}{dz}
\Rightarrow r\frac{dv_z}{dr}=\frac{r^2}{2\mu}\frac{dp}{dz}+C_1.
$$
轴心有限性 $C_1=0$：
$$
\frac{dv_z}{dr}=\frac{r}{2\mu}\frac{dp}{dz}
\Rightarrow v_z=\frac{r^2}{4\mu}\frac{dp}{dz}+C_2.
$$
无滑移 $v_z(R)=0$：
$$
v_z(r)=\frac{1}{4\mu}\Big(\frac{dp}{dz}\Big)(r^2-R^2).
$$
对常压梯度 $dp/dz=-\Delta p/L$：
$$
v_z(r)=\frac{\Delta p}{4\mu L}(R^2-r^2).
$$
流量
$$
Q=\int_0^R2\pi r v_z\,dr
=\frac{\pi\Delta p}{2\mu L}\int_0^R(rR^2-r^3)dr
=\frac{\pi\Delta p}{2\mu L}\Big(\frac{R^4}{2}-\frac{R^4}{4}\Big)
=\frac{\pi R^4\Delta p}{8\mu L}.
$$

物理意义：流量强烈依赖半径（$R^4$）。几何意义：抛物线速度剖面，壁面剪切率 $\dot\gamma=4Q/(\pi R^3)$。

特例：$\Delta p=0$ 无流动；$R\to0$ 流量迅速趋于零。

应用：血管阻力、注射针头、呼吸气流（近似）。

#### 8.4 血流阻力与串并联

定义：流体阻力 $R_f=\Delta p/Q$。Poiseuille：$R_f=\frac{8\mu L}{\pi R^4}$。

假设：层流、刚性圆管。

公式：串联 $R=R_1+R_2$；并联 $\frac1R=\sum\frac1{R_i}$。

推导：串联两段，流量相同 $Q$，总压降 $\Delta p=\Delta p_1+\Delta p_2=R_1Q+R_2Q$，故 $R=R_1+R_2$。并联各支路压降相同，总流量 $Q=\sum Q_i=\Delta p\sum1/R_i$，故 $1/R=\sum1/R_i$。

物理意义：血管网络阻力决定血压分布。几何意义：与电路 Kirchhoff 定律同构。

特例：严重狭窄（$R$ 减小）使阻力剧增，因 $R^{-4}$。

应用：体循环总阻力、冠脉血流、微循环。

#### 8.5 斯托克斯流与低雷诺数游泳

定义：$\mathrm{Re}\ll1$ 时惯性忽略，运动可逆。

假设：不可压、定常或准定常、忽略惯性。

公式：
$$
-\nabla p+\mu\nabla^2\mathbf v=0,\qquad \nabla\cdot\mathbf v=0.
$$

推导：N-S 方程两边除以惯性尺度 $\rho U^2/L$，粘性项为 $\frac{\mu U/L^2}{\rho U^2/L}=\frac1{\mathrm{Re}}$。$\mathrm{Re}\to0$ 时惯性项可忽略，得 Stokes 方程。其解具有时间反演对称性（Scallop 定理）：可逆变形不产生净位移。

物理意义：微生物需非往复运动（螺旋、波状）才能推进。几何意义：Stokes 流场由边界瞬时几何唯一决定，无记忆。

特例：球体阻力 $F=6\pi\mu av$；鞭毛螺旋推进 $F\propto v$。

应用：细菌鞭毛、精子、微纳机器人。

#### 8.6 Womersley 数与脉动流

定义：$\alpha=R\sqrt{\omega\rho/\mu}$（Womersley 数），刻画非定常惯性力与粘性力之比。

假设：圆管中正弦压力梯度驱动、层流、轴对称。

公式：
$$
v_z(r,t)=\frac{\Delta p}{\rho L}\frac{1}{i\omega}\Big[1-\frac{J_0(\alpha i^{3/2}r/R)}{J_0(\alpha i^{3/2})}\Big]e^{i\omega t}.
$$

推导：轴向 N-S 方程（充分发展）
$$
\rho\frac{\partial v_z}{\partial t}=-\frac{\partial p}{\partial z}+\mu\frac{1}{r}\frac{\partial}{\partial r}\Big(r\frac{\partial v_z}{\partial r}\Big).
$$
设 $p=p_0e^{i\omega t}$、$v_z=V(r)e^{i\omega t}$：
$$
i\omega\rho V=\frac{\Delta p}{L}+\mu\Big(V''+\frac1rV'\Big).
$$
无量纲化 $x=r/R$，令 $\alpha^2=\omega\rho R^2/\mu$：
$$
V''+\frac1xV'-\alpha^2 i\,V=-\frac{\Delta p R^2}{\mu L}.
$$
零阶 Bessel 方程的解为 $J_0(\alpha i^{3/2}x)$。由轴心有限、壁面无滑移 $V(1)=0$ 得上式。

物理意义：$\alpha\ll1$ 时准稳态（Poiseuille 剖面），$\alpha\gg1$ 时中心速度平坦、近壁剪切集中，相位滞后。几何意义：$\alpha$ 是管半径与 Stokes 层厚度 $\sqrt{\mu/(\rho\omega)}$ 之比。

特例：$\alpha\to0$ 退化为 Poiseuille；$\alpha\to\infty$ 为无粘振荡流。

应用：主动脉脉动流、脑脊液、呼吸振荡流。

### 第9章 血液流变学与循环力学

#### 9.1 非牛顿流体与 Casson 模型

定义：剪切应力与剪切率非线性关系的流体。

假设：血液含红细胞，可近似为 Casson 塑性流体。

公式：
$$
\sqrt{\tau}=\sqrt{\tau_y}+\sqrt{\mu_p\dot\gamma}\quad(\tau>\tau_y),\qquad \dot\gamma=0\quad(\tau\le\tau_y).
$$

推导：悬浮液微观模型假设红细胞形成屈曲链，其解聚应力为 $\tau_y$。对链的统计平均得到 $\sqrt\tau$ 线性关系。写成
$$
\dot\gamma=\begin{cases}\dfrac{1}{\mu_p}(\sqrt\tau-\sqrt{\tau_y})^2,&\tau>\tau_y\\[4pt]0,&\tau\le\tau_y.\end{cases}
$$
管内积分可得 Casson 流的速度剖面与流量（分“栓塞核”与“剪切区”两部分）。

物理意义：$\tau_y$ 为屈服应力，低于它血液不流动。几何意义：$\sqrt\tau$–$\sqrt{\dot\gamma}$ 为直线，斜率 $\sqrt{\mu_p}$。

特例：$\tau_y=0$ 时退化为牛顿流体；高剪切率时 $\dot\gamma\approx\tau/\mu_p$。

应用：小血管/微循环血液流变、血栓形成、动脉粥样硬化。

#### 9.2 Fahraeus–Lindqvist 效应

定义：血液流经细管时表观粘度随管径减小而降低。

假设：红细胞在细管中向轴心迁移，形成近壁无细胞血浆层。

公式：表观粘度 $\mu_{\text{app}}=\dfrac{\pi R^4\Delta p}{8QL}$ 随 $R$ 下降而下降至约 $R\sim5\ \mu\mathrm m$，再下降转为上升（Fahraeus 逆效应）。

推导：设血浆层厚度 $\delta$、血浆粘度 $\mu_p$、核心区粘度 $\mu_c$。两层流量叠加：
$$
Q=\frac{\pi\Delta p}{8L}\Big[\frac{(R-\delta)^4}{\mu_c}+\frac{R^4-(R-\delta)^4}{\mu_p}\Big].
$$
表观粘度由 $Q$ 反解
$$
\mu_{\text{app}}=\frac{\pi R^4\Delta p}{8QL}
=R^4\Big[\frac{(R-\delta)^4}{\mu_c}+\frac{R^4-(R-\delta)^4}{\mu_p}\Big]^{-1}.
$$
当 $\delta$ 相对 $R$ 增大时，$\mu_{\text{app}}$ 下降。

物理意义：血浆层起“润滑”作用，降低细管阻力。几何意义：分层流的速度剖面对应并联电阻网络。

特例：$\delta=0$ 或 $R\gg\delta$ 时无效应，$\mu_{\text{app}}\to\mu_c$；$R$ 接近红细胞尺度时逆效应出现。

应用：微血管阻力、血液稀释、人工器官设计。

#### 9.3 Windkessel 模型

定义：用集总参数（阻力 $R$、顺应性 $C$）描述动脉系统压力-流量关系。

假设：动脉为大容器、顺应性集中、血流阻力集中在外周。

公式（两元件）：
$$
C\frac{dP}{dt}=Q_{\text{in}}-\frac{P}{R}.
$$

推导：动脉容积 $V$，顺应性 $C=dV/dP$。质量守恒：流入 $Q_{in}$ 等于储存速率 $dV/dt$ 加外流 $P/R$：
$$
Q_{in}=\frac{dV}{dt}+\frac PR=C\frac{dP}{dt}+\frac PR.
$$
舒张期 $Q_{in}=0$，解
$$
P(t)=P_0e^{-t/(RC)}.
$$
时间常数 $\tau=RC$ 决定舒张压衰减。

物理意义：动脉弹性储存心脏搏出能量、平滑血流。几何意义：压力指数衰减，$RC$ 为衰减时间。

特例：$C\to0$ 刚性管时 $P=RQ$ 瞬时响应；$R\to\infty$ 无外流，压力持恒。

应用：血压波形分析、动脉硬化诊断（脉搏波速度）、心室-动脉耦合。

#### 9.4 脉搏波速度（Moens–Korteweg）

定义：压力波沿弹性血管的传播速度。

假设：薄壁弹性圆管、不可压流体、管壁线性弹性。

公式：
$$
c=\sqrt{\frac{Eh}{2\rho R}},\qquad
c_0=\sqrt{\frac{K}{\rho}}\ \text{(纯流体声速)}.
$$

推导：局部管段动量与连续性。设管半径扰动 $\delta R$、压力扰动 $\delta p$。管壁弹性
$$
\frac{\delta R}{R}=\frac{R}{Eh}\delta p
\Rightarrow \delta A=2\pi R\delta R=\frac{2\pi R^3}{Eh}\delta p.
$$
连续性 $\partial_t A+\partial_x(Av)=0$ 线性化：
$$
\partial_t\delta A+A_0\partial_x v=0.
$$
动量线性化 $\rho\partial_tv=-\partial_x\delta p$。消去 $v$：对连续性求 $\partial_t$、对动量求 $\partial_x$：
$$
\partial_t^2\delta A=A_0\partial_t\partial_x v\cdot(-1)\cdots
$$
整理得
$$
\partial_t^2\delta p=\frac{A_0}{2\pi R^3/(Eh)}\cdot\frac1\rho\,\partial_x^2\delta p
=\frac{Eh}{2\rho R}\partial_x^2\delta p.
$$
故波速 $c=\sqrt{Eh/(2\rho R)}$。

物理意义：动脉越硬（$E$ 越大）脉搏波越快，可作心血管风险指标。几何意义：波速由流体惯性与管壁弹性竞争决定。

特例：$E\to\infty$ 刚壁时 $c\to\infty$（刚性管极限的准静态压力传递）；考虑流体可压需加 $K$ 项。

应用：cf-PWV 无创测量、动脉僵硬度评估。

#### 9.5 心脏做功与输出量

定义：心脏每搏做功、功率与效率。

假设：左心每搏输出量 $SV$、平均动脉压 $MAP$、心率 $HR$。

公式：
$$
W_{stroke}=SV\cdot MAP,\qquad
P=W_{stroke}\cdot HR,\qquad
\eta=\frac{P}{E_{\text{metabolic}}}.
$$

推导：左心室将体积 $SV$ 的血从低压（接近 0）推入平均压 $MAP$ 的动脉，克服压力做功
$$
W=\int p\,dV\approx MAP\cdot SV.
$$
功率即单位时间功 $P=W\cdot HR$。心肌耗氧对应代谢能，效率约 $20\%$–$25\%$（与骨骼肌相近）。

物理意义：心脏是恒压泵，$MAP\times SV$ 近似每搏功。几何意义：$p$-$V$ 环面积即搏功。

特例：运动时 $SV$、$HR$ 增大，$P$ 线性增加；心衰时效率下降。

应用：血流动力学监测、心功能评估、循环辅助装置。

#### 9.6 Starling 毛细血管滤过方程

定义：毛细血管壁两侧液体交换的净通量。

假设：半透膜、水与溶质耦合、渗透压近似 van't Hoff、反射系数 $\sigma$。

公式：
$$
J_v=L_pS\big[(P_c-P_i)-\sigma(\pi_c-\pi_i)\big],
$$
$$
J_s=P_sS\,\Delta c+(1-\sigma)\bar c\,J_v.
$$

推导：净驱动压为静水压差减去有效渗透压差：
$$
\Delta P_{\text{eff}}=(P_c-P_i)-\sigma(\pi_c-\pi_i).
$$
滤过通量与驱动压成正比，比例系数 $L_pS$（水力传导×面积）。溶质通量含扩散项 $P_s\Delta c$ 与溶剂拖曳项 $(1-\sigma)\bar c J_v$。渗透压由 van't Hoff 律 $\pi=RT\Delta c$ 给出。

物理意义：净滤过在动脉端为正、静脉端为负，完成组织液交换。几何意义：$J_v=0$ 的等滤过线由静水压与渗透压平衡决定。

特例：$\sigma=1$ 半透膜无溶质拖曳；$\sigma=0$ 无渗透效应，$J_v=L_pS\Delta P$。

应用：水肿形成、肾小球滤过、组织工程灌注。

### 第10章 骨骼与软组织的固体力学

#### 10.1 Euler–Bernoulli 梁

定义：细长梁弯曲理论。

假设：平截面保持平面、小挠度、忽略剪切变形。

公式：
$$
EI\frac{d^2w}{dx^2}=M(x),\qquad
EI\frac{d^4w}{dx^4}=q(x).
$$

推导：曲率 $\kappa=d^2w/dx^2$（小挠度）。应变沿厚度 $\varepsilon=-y\kappa=-y\,w''$。应力 $\sigma=E\varepsilon=-Ey\,w''$。截面弯矩
$$
M=\int_A y\sigma\,dA=-Ew''\int_A y^2\,dA=-EI\,w'',
$$
其中 $I=\int y^2dA$。故 $EIw''=-M$（符号约定）。对 $x$ 求导两次并代入剪力关系 $dM/dx=V$、$dV/dx=q$ 得
$$
EIw''''=q.
$$

物理意义：弯曲刚度 $EI$ 越大越难弯。几何意义：中性轴 $y=0$ 无应变，弯矩与曲率成正比。

特例：$q=0$ 时挠度为三次多项式；简支梁得 $w\propto x(L^3-2Lx^2+x^3)$。

应用：骨梁力学、椎骨、内植物、细胞牵引梁。

#### 10.2 骨的功能适应（Wolff 定律）

定义：骨结构随力学载荷重塑以维持适宜应变。

假设：局部应变反馈、稳态应变 $\varepsilon^*$。

公式：
$$
\frac{dm}{dt}=k(\varepsilon-\varepsilon^*),\qquad
\rho\ddot u=\nabla\cdot\boldsymbol\sigma.
$$

推导：骨密度 $m$ 的变化率与应变偏差成正比（一阶反馈）。结合弹性平衡方程，形成自适应弹性问题。稳态时 $\varepsilon=\varepsilon^*$，结构达到最优。若以能量密度度量，$U=\frac12\sigma\varepsilon$ 被组织趋向于均匀化。

物理意义：载荷决定骨小梁取向与密度。几何意义：应变场均匀化是变分问题的欧拉-拉格朗日条件。

特例：失重时 $\varepsilon$ 减小，骨密度下降；高负荷时增强。

应用：骨质疏松、航天骨丢失、康复训练设计。

#### 10.3 Euler 屈曲

定义：细长压杆失稳的临界载荷。

假设：理想直杆、轴向压力 $P$、小横向挠度。

公式：
$$
P_{cr}=\frac{\pi^2EI}{(\beta L)^2}.
$$

推导：轴向压力 $P$ 作用下弯矩 $M=-Pw$。由 $EIw''=-M$：
$$
EIw''=Pw\ \Rightarrow\ w''+k^2w=0,\quad k^2=\frac{P}{EI}.
$$
通解 $w=A\sin kx+B\cos kx$。边界简支 $w(0)=w(L)=0$ 给 $B=0$，$\sin kL=0$，即 $kL=n\pi$。最小临界
$$
P_{cr}=\frac{\pi^2EI}{L^2}.
$$
其他边界条件引入有效长度 $\beta L$。

物理意义：受压长杆存在分岔失稳。几何意义：$P$-$w$ 曲线在 $P_{cr}$ 处出现叉式分岔。

特例：$P<P_{cr}$ 无横向变形；$P>P_{cr}$ 大变形后硬化（后屈曲）。

应用：骨小梁、细胞微管、内植物稳定性。

#### 10.4 Griffith 断裂力学

定义：裂纹扩展的能量判据。

假设：线弹性、脆性材料、平面裂纹。

公式：
$$
G=\frac{\pi\sigma^2 a}{E}\ \text{(平面应力)},\qquad
G_c=\frac{\pi\sigma_c^2 a}{E}.
$$

推导：含中心裂纹 $2a$ 的无限板，能量释放率由应力强度因子 $K_I=\sigma\sqrt{\pi a}$ 给出
$$
G=\frac{K_I^2}{E}=\frac{\pi\sigma^2 a}{E}.
$$
裂纹扩展释放弹性能率 $\propto\sigma^2a/E$，新表面需功 $2\gamma_s$，故临界
$$
\frac{\pi\sigma_c^2a}{E}=2\gamma_s\ \Rightarrow\ \sigma_c=\sqrt{\frac{2E\gamma_s}{\pi a}}.
$$
若定义 $G_c=2\gamma_s$ 则得 $G_c=\pi\sigma_c^2a/E$。

物理意义：强度随裂纹长度平方根反比下降。几何意义：$K_I$ 为裂尖应力场奇异性强度。

特例：$a\to0$ 理论强度；长裂纹易断。

应用：骨断裂、疲劳、血管夹层、假体裂纹。

#### 10.5 软组织非线性弹性（超弹性）

定义：软组织应力-应变高度非线性，常用应变能函数描述。

假设：不可压、各向同性或横向各向同性、超弹性。

公式（Neo-Hookean）：
$$
W=\frac{\mu}{2}(I_1-3),\qquad
\boldsymbol\sigma=-p\mathbf I+\mu\mathbf b,\quad \mathbf b=\mathbf F\mathbf F^T.
$$

推导：各向同性应变能 $W(I_1,I_2,I_3)$，不可压 $I_3=1$。Cauchy 应力
$$
\boldsymbol\sigma=-p\mathbf I+2\frac{\partial W}{\partial I_1}\mathbf b
-2\frac{\partial W}{\partial I_2}\mathbf b^{-1}.
$$
对 Neo-Hookean $W=\frac\mu2(I_1-3)$，$\partial W/\partial I_1=\mu/2$，得 $\boldsymbol\sigma=-p\mathbf I+\mu\mathbf b$。单轴拉伸 $\lambda$，$\lambda_2=\lambda_3=1/\sqrt\lambda$，由 $\sigma_{22}=\sigma_{33}=0$ 定 $p$，得
$$
\sigma_{11}=\mu\Big(\lambda^2-\frac1\lambda\Big).
$$

物理意义：非线性来自几何大变形与纤维重排。几何意义：$\mathbf b$ 描述局部长度平方变化。

特例：$\lambda\to1$ 时 $\sigma\approx3\mu(\lambda-1)$，与线弹性一致；大拉伸时纤维硬化（加纤维项）。

应用：皮肤、血管、韧带、水凝胶。

#### 10.6 骨的疲劳与累积损伤

定义：循环载荷下微裂纹累积导致疲劳破坏。

假设：线弹性断裂、Paris 律裂纹扩展、Miner 线性损伤累积。

公式：
$$
\frac{da}{dN}=C(\Delta K)^m,\qquad
\Delta K=\Delta\sigma\sqrt{\pi a},\qquad
D=\sum_i\frac{n_i}{N_i}.
$$

推导：裂尖应力强度因子幅值 $\Delta K=\Delta\sigma\sqrt{\pi a}$。由 Paris 律分离变量：
$$
\frac{da}{dN}=C(\Delta\sigma)^m(\pi a)^{m/2}
\Rightarrow
N_f=\frac{1}{C(\Delta\sigma)^m\pi^{m/2}}\int_{a_0}^{a_c}a^{-m/2}da.
$$
当 $m\ne2$：
$$
N_f=\frac{2}{C(2-m)(\Delta\sigma)^m\pi^{m/2}}
\big(a_c^{1-m/2}-a_0^{1-m/2}\big).
$$
Miner 法则把变幅载荷折算为各应力级损伤之和，$D=1$ 判据。

物理意义：裂纹扩展速率是幂律，寿命对裂纹尺寸与应力幅高度敏感。几何意义：$\Delta K$ 度量裂尖奇异场强度，$a_c$ 为断裂韧性阈值。

特例：$m=2$ 时积分给对数形式 $N_f\propto\ln(a_c/a_0)$；$a_0\to0$ 时 $N_f\to\infty$（无初始缺陷）。

应用：应力性骨折、人工关节疲劳、骨水泥裂纹。

### 第11章 肌肉力学与运动

#### 11.1 Hill 力-速方程

定义：肌肉收缩力 $F$ 与收缩速度 $v$ 的经验双曲线关系。

假设：稳态、等张力、能量守恒。

公式：
$$
(F+a)(v+b)=(F_0+a)b,\qquad
v=b\frac{F_0-F}{F+a}.
$$

推导：由能量守恒：肌肉释放的化学能率 $E$ 分为维持热 $a v$、缩短热与机械功：
$$
E=av+Fv.
$$
设最大效率附近能量率与力成正比，$\dot E=b(F_0+a)$（Hill 的经验假设），则
$$
av+Fv=b(F_0+a)
\Rightarrow (F+a)(v+b)=(F_0+a)b.
$$
解出
$$
v=\frac{b(F_0-F)}{F+a}.
$$

物理意义：力越大速度越小，最大力时速度为零。几何意义：$F$-$v$ 为等轴双曲线，渐近 $v=b(F_0-F)/(F+a)$。

特例：$F=0$ 时 $v_{\max}=bF_0/a$；$v=0$ 时 $F=F_0$（等长最大力）。

应用：肌肉功率优化、运动生物力学、康复机器人。

#### 11.2 肌肉功率与最优负载

定义：机械功率 $P=Fv$。

假设：Hill 方程、$a/F_0$ 给定。

公式：
$$
P(F)=F\,b\frac{F_0-F}{F+a}.
$$

推导：由 Hill 方程 $v=b(F_0-F)/(F+a)$，故
$$
P=Fv=b\frac{F(F_0-F)}{F+a}.
$$
最大化：令 $dP/dF=0$，
$$
\frac{d}{dF}\frac{F(F_0-F)}{F+a}=0.
$$
计算
$$
\frac{d}{dF}\Big(\frac{F F_0-F^2}{F+a}\Big)
=\frac{(F_0-2F)(F+a)-(F F_0-F^2)}{(F+a)^2}=0.
$$
分子：
$$
(F_0-2F)(F+a)-F(F_0-F)=-2F^2-2aF+aF_0=0
\Rightarrow 2F^2+2aF-aF_0=0.
$$
解得
$$
F^*=-a+\sqrt{a^2+\frac{aF_0}{2}}.
$$
代入得最大功率。

物理意义：存在最优负载使输出功率最大。几何意义：$P(F)$ 为开口向下抛物线状，峰值在 $F^*$。

特例：$a\ll F_0$ 时 $F^*\approx\sqrt{aF_0/2}$；$a$ 很小时效率更高。

应用：自行车/划船最优阻力、外骨骼助力。

#### 11.3 拉格朗日力学与肢体运动

定义：用广义坐标描述生物运动的力学。

假设：完整约束、保守力为主。

公式：
$$
L=T-V,\qquad
\frac{d}{dt}\frac{\partial L}{\partial\dot q_i}-\frac{\partial L}{\partial q_i}=Q_i^{\text{nc}}.
$$

推导：由 Hamilton 原理 $\delta\int L\,dt=0$，对路径变分
$$
\delta\int_{t_1}^{t_2}L(q,\dot q,t)dt=0.
$$
分部积分
$$
\int\Big(\frac{\partial L}{\partial q_i}\delta q_i+\frac{\partial L}{\partial\dot q_i}\delta\dot q_i\Big)dt
=\int\Big(\frac{\partial L}{\partial q_i}-\frac{d}{dt}\frac{\partial L}{\partial\dot q_i}\Big)\delta q_i\,dt
+\Big[\frac{\partial L}{\partial\dot q_i}\delta q_i\Big]=0.
$$
边界 $\delta q_i=0$，得 Euler-Lagrange 方程。非保守力以广义力 $Q_i^{nc}$ 加入右端。

物理意义：生物运动可归结为构型空间中的测地线或最小作用量。几何意义：$L$ 是切丛上的函数，方程为其变分极值。

特例：无势能时 $L=T$，得自由粒子沿最小弧长运动。

应用：步态建模、机器人腿、关节力矩反解。

#### 11.4 单摆步态模型

定义：把腿近似为无质量杆、躯干为质点的倒立摆/单摆。

假设：无外力、小角度、髋为铰链。

公式：
$$
L=\frac12ml^2\dot\theta^2+mgl\cos\theta,\qquad
\ddot\theta+\frac gl\sin\theta=0.
$$

推导：动能 $T=\frac12ml^2\dot\theta^2$，势能 $V=-mgl\cos\theta$，故 $L=T-V$。Euler-Lagrange：
$$
\frac{\partial L}{\partial\dot\theta}=ml^2\dot\theta,\quad
\frac{d}{dt}=ml^2\ddot\theta,\quad
\frac{\partial L}{\partial\theta}=-mgl\sin\theta.
$$
故
$$
ml^2\ddot\theta+mgl\sin\theta=0
\Rightarrow \ddot\theta+\frac gl\sin\theta=0.
$$
小角度 $\sin\theta\approx\theta$：$\ddot\theta+\frac gl\theta=0$，周期 $T=2\pi\sqrt{l/g}$。

物理意义：步行摆动腿类似单摆，步频与腿长平方根反比。几何意义：在 $(\theta,\dot\theta)$ 相平面闭合轨道。

特例：大角度周期含椭圆积分 $T=4\sqrt{l/g}K(\sin(\theta_0/2))$；$\theta_0\to\pi$ 周期发散。

应用：步态周期、外骨骼频率匹配、运动能耗。

#### 11.5 哈密顿力学在生物运动中的应用

定义：$H=\sum_ip_i\dot q_i-L$。

假设：无显含时的保守系统。

公式：
$$
\dot q_i=\frac{\partial H}{\partial p_i},\qquad
\dot p_i=-\frac{\partial H}{\partial q_i}.
$$

推导：由 $H=\sum p_i\dot q_i-L$，取微分
$$
dH=\sum(\dot q_i\,dp_i+p_i\,d\dot q_i)-\sum\frac{\partial L}{\partial q_i}dq_i-\sum\frac{\partial L}{\partial\dot q_i}d\dot q_i.
$$
利用 $p_i=\partial L/\partial\dot q_i$，$p_i\,d\dot q_i$ 与最后一项抵消，得
$$
dH=\sum\dot q_i\,dp_i-\sum\frac{\partial L}{\partial q_i}dq_i.
$$
由 Euler-Lagrange $\partial L/\partial q_i=\dot p_i$，故
$$
dH=\sum\dot q_i\,dp_i-\sum\dot p_i\,dq_i,
$$
比较 $dH=\sum\frac{\partial H}{\partial p_i}dp_i+\sum\frac{\partial H}{\partial q_i}dq_i$ 得 Hamilton 方程。

物理意义：$H$ 通常为总能量，运动是相空间中的哈密顿流。几何意义：Hamilton 方程保持相空间体积（Liouville）。

特例：$L=\frac12m\dot q^2-V$ 时 $H=\frac{p^2}{2m}+V$。

应用：步态稳定性、生物系统能量守恒、控制律设计。

#### 11.6 力-长度关系与功环

定义：肌肉主动力随肌节长度变化，并结合被动力形成长度循环。

假设：主动力来自粗细肌丝重叠；被动力为指数型；循环准静态。

公式：
$$
F(l)=F_0\exp\!\Big[-\frac{(l-l_0)^2}{2w^2}\Big],\qquad
F_p(l)=A\big(e^{\beta(l-l_0)}-1\big),
$$
$$
W_{\text{net}}=\oint\big[F(l)+F_p(l)\big]\,dl.
$$

推导：可形成的横桥数 $N(l)$ 正比于粗细肌丝重叠长度，平台段近似为钟形，故取高斯形主动力；被动弹性来自 titin 等，取指数硬化。长度循环 $l(t)=l_0+\Delta l\sin\omega t$ 时
$$
W_{\text{net}}=\int_0^{2\pi/\omega}(F+F_p)\dot l\,dt,
$$
即力-长度平面内闭合曲线的面积。结合 Hill 力-速关系，可把 $F$ 写成 $F=F(v(l),l)$ 得到完整的功环模型。

物理意义：肌肉在拉长时储存弹性势能、缩短时释放做功，净功环面积为正。几何意义：$F$-$l$ 平面闭合曲线面积即每循环净机械功。

特例：$\Delta l\to0$ 时 $W\to0$；无被动力时功环退化为主动力环。

应用：跳跃/跑步储能、肌肉效率、外骨骼能量回收。

### 第12章 生物传热与传质

#### 12.1 热传导方程

定义：温度场 $T(\mathbf x,t)$ 的演化。

假设：连续介质、Fourier 定律 $\mathbf q=-k\nabla T$、无对流。

公式：
$$
\rho c_p\frac{\partial T}{\partial t}=k\nabla^2T+q,\qquad
\alpha=\frac{k}{\rho c_p}.
$$

推导：能量守恒：单元体积内能变化率等于流入热量加源：
$$
\rho c_p\frac{\partial T}{\partial t}=-\nabla\cdot\mathbf q+q.
$$
代入 Fourier 定律：
$$
\rho c_p\frac{\partial T}{\partial t}=\nabla\cdot(k\nabla T)+q.
$$
$k$ 常数得 $\partial_tT=\alpha\nabla^2T+q/(\rho c_p)$。

物理意义：热扩散系数 $\alpha$ 决定温度场抹平速率。几何意义：与 Fick 第二定律同型，$\alpha$ 为热扩散系数。

特例：稳态一维 $\nabla^2T=0$ 给线性温度分布；球对称得 $T\propto1/r$（点源）。

应用：组织温度分布、热疗、冷疗、冷冻保存。

#### 12.2 Pennes 生物传热方程

定义：计入血流灌注的组织传热方程。

假设：血液在毛细血管网中与组织局部热平衡、灌注项线性。

公式：
$$
\rho c_p\frac{\partial T}{\partial t}
=\nabla\cdot(k\nabla T)+\omega_b\rho_bc_b(T_a-T)+q_m.
$$

推导：在热传导方程中加入血流携带的热量。血液进入组织温度 $T_a$、流出接近组织温 $T$，单位体积血液灌注率 $\omega_b$，则灌注净热 $\omega_b\rho_bc_b(T_a-T)$。加上代谢产热 $q_m$：
$$
\rho c_p\partial_tT=\nabla\cdot(k\nabla T)+\omega_b\rho_bc_b(T_a-T)+q_m.
$$
稳态且忽略传导时
$$
T=T_a+\frac{q_m}{\omega_b\rho_bc_b}.
$$

物理意义：血液灌注是组织散热/加热的重要机制，可用 $T$ 反演灌注率。几何意义：方程是带线性回复项的反应-扩散方程，其衰减长度 $\sqrt{k/(\omega_b\rho_bc_b)}$。

特例：$\omega_b=0$ 退化为纯热传导；灌注很大时组织温度趋近动脉血温。

应用：肿瘤热疗、微波消融、无创体温监测。

#### 12.3 传质边界层与 Sherwood 数

定义：流动中浓度边界层与对流传质。

假设：稳态层流、平板、低质量通量。

公式：
$$
\mathrm{Sh}=\frac{k_mL}{D},\qquad
\mathrm{Sh}=0.664\,\mathrm{Re}^{1/2}\mathrm{Sc}^{1/3}\ (\text{层流平板}),
$$
$$
\mathrm{Sc}=\frac{\nu}{D}.
$$

推导：浓度边界层方程
$$
u\partial_xc+v\partial_yc=D\partial_y^2c.
$$
速度边界层厚度 $\delta\sim x\mathrm{Re}_x^{-1/2}$，浓度边界层 $\delta_c\sim x\mathrm{Re}_x^{-1/2}\mathrm{Sc}^{-1/3}$。传质系数 $k_m\sim D/\delta_c$，故
$$
k_m\sim\frac{D}{x}\mathrm{Re}_x^{1/2}\mathrm{Sc}^{1/3}
\Rightarrow \mathrm{Sh}=0.664\mathrm{Re}^{1/2}\mathrm{Sc}^{1/3}.
$$

物理意义：对流传质由流动与分子扩散共同控制。几何意义：$\mathrm{Sh}$ 是界面浓度梯度（无量纲）。

特例：$\mathrm{Sc}=1$ 时浓度边界层与速度边界层同厚；$\mathrm{Re}\to0$ 时纯扩散。

应用：血液-组织氧交换、透析器、生物膜反应器。

#### 12.4 热穿透深度

定义：周期热源作用下温度振幅衰减到 $1/e$ 的深度。

假设：半无限介质、正弦表面温度、一维、常物性。

公式：
$$
\delta=\sqrt{\frac{2\alpha}{\omega}},\qquad
T(x,t)=T_0e^{-x/\delta}\cos\!\Big(\omega t-\frac{x}{\delta}\Big).
$$

推导：一维热传导 $\partial_tT=\alpha\partial_x^2T$。设 $T=T_0e^{i\omega t+\lambda x}$，代入
$$
i\omega=\alpha\lambda^2\Rightarrow \lambda=\pm\sqrt{\frac{i\omega}{\alpha}}=(1+i)\sqrt{\frac{\omega}{2\alpha}}.
$$
取衰减解 $\lambda=-(1+i)\sqrt{\omega/(2\alpha)}$，实部衰减 $e^{-x/\delta}$，相位滞后 $x/\delta$，其中 $\delta=\sqrt{2\alpha/\omega}$。

物理意义：频率越高，热穿透越浅。几何意义：温度波是幅值指数衰减、相位线性的行波。

特例：$\omega\to0$ 稳态时 $\delta\to\infty$；$\omega\to\infty$ 时 $\delta\to0$。

应用：射频/微波热疗深度、皮肤温度波、热诊断。

#### 12.5 集总参数热模型与 Biot 数

定义：当内部导热远快于表面散热时，物体温度近似均匀。

假设：Biot 数 $\mathrm{Bi}=hL/k\ll1$。

公式：
$$
\rho cV\frac{dT}{dt}=-hA(T-T_\infty)
\Rightarrow T-T_\infty=(T_0-T_\infty)e^{-t/\tau},
$$
$$
\tau=\frac{\rho cV}{hA}.
$$

推导：$\mathrm{Bi}\ll1$ 时体内无温度梯度，能量守恒为
$$
\rho cV\frac{dT}{dt}=-hA(T-T_\infty).
$$
分离变量积分得指数衰减，时间常数 $\tau=\rho cV/(hA)$。与一阶系统同型。

物理意义：小物体快速趋于环境温度。几何意义：$\mathrm{Bi}$ 为内部导热热阻与表面换热热阻之比。

特例：$h\to0$ 绝热时 $\tau\to\infty$；$V/A\to0$ 薄物体迅速平衡。

应用：体温调节、冷冻保存、微珠热疗。

#### 12.6 Kedem–Katchalsky 膜传质

定义：不可逆热力学描述膜两侧体积流与溶质流的耦合。

假设：线性近平衡、Onsager 倒易关系、稀溶液。

公式：
$$
J_v=L_p(\Delta P-\sigma\Delta\pi),\qquad
J_s=P_s\Delta c+(1-\sigma)\bar c\,J_v.
$$

推导：取力 $X_v=\Delta P$、$X_s=-\Delta\pi$（或 $\Delta c$），流线性叠加：
$$
J_v=L_{vv}X_v+L_{vs}X_s,\qquad J_s=L_{sv}X_v+L_{ss}X_s.
$$
Onsager 倒易 $L_{vs}=L_{sv}$。定义过滤系数 $L_p=L_{vv}$、反射系数 $\sigma=-L_{vs}/L_{vv}$ 与溶质通透 $P_s$，并令 $L_{ss}=P_s$，即得两式。

物理意义：$\sigma$ 度量膜对溶质的截留能力，$\sigma=1$ 完全截留。几何意义：$(J_v,J_s)$ 对流（力）的线性响应平面。

特例：$\sigma=1$ 时 $J_s=P_s\Delta c$（纯扩散）；$\sigma=0$ 时 $J_s=P_s\Delta c+\bar c J_v$（完全拖曳）。

应用：透析、反渗透、上皮转运。

### 第13章 酶动力学与代谢网络

#### 13.1 Michaelis–Menten 方程

定义：单底物酶反应速率。

假设：准稳态近似 $d[ES]/dt=0$、$[S]\gg[E]_0$、不可逆产物生成。

公式：
$$
v=\frac{V_{\max}[S]}{K_m+[S]},\qquad
V_{\max}=k_2[E]_0,\quad K_m=\frac{k_{-1}+k_2}{k_1}.
$$

推导：反应
$$
E+S\underset{k_{-1}}{\overset{k_1}{\rightleftharpoons}}ES\overset{k_2}{\to}E+P.
$$
速率
$$
v=k_2[ES],\quad \frac{d[ES]}{dt}=k_1[E][S]-k_{-1}[ES]-k_2[ES].
$$
酶守恒 $[E]=[E]_0-[ES]$。准稳态
$$
0=k_1([E]_0-[ES])[S]-(k_{-1}+k_2)[ES].
$$
解出
$$
[ES]=\frac{k_1[E]_0[S]}{k_{-1}+k_2+k_1[S]}=\frac{[E]_0[S]}{K_m+[S]},
$$
故
$$
v=\frac{k_2[E]_0[S]}{K_m+[S]}=\frac{V_{\max}[S]}{K_m+[S]}.
$$

物理意义：$K_m$ 是速度为半最大时的底物浓度，反映亲和力。几何意义：$v$-$[S]$ 为饱和双曲线，$[S]\to\infty$ 时趋于 $V_{\max}$。

特例：$[S]\ll K_m$ 时 $v\approx(V_{\max}/K_m)[S]$ 一级；$[S]\gg K_m$ 时 $v\approx V_{\max}$ 零级。

应用：酶活性测定、代谢通量控制、药物靶点分析。

#### 13.2 Lineweaver–Burk 线性化

定义：MM 方程的倒数线性形式。

公式：
$$
\frac1v=\frac{K_m}{V_{\max}}\frac1{[S]}+\frac1{V_{\max}}.
$$

推导：由 $v=V_{\max}[S]/(K_m+[S])$，取倒数
$$
\frac1v=\frac{K_m+[S]}{V_{\max}[S]}
=\frac{K_m}{V_{\max}}\frac1{[S]}+\frac{1}{V_{\max}}.
$$
故 $1/v$ 对 $1/[S]$ 为直线，斜率 $K_m/V_{\max}$，截距 $1/V_{\max}$，横截距 $-1/K_m$。

物理意义：便于从线性拟合提取 $K_m,V_{\max}$。几何意义：双倒数图把双曲线映射为直线。

特例：$1/[S]\to0$ 得 $1/V_{\max}$；$1/v=0$ 得 $[S]=-K_m$。

应用：抑制剂类型鉴别、酶参数拟合。

#### 13.3 竞争性与非竞争性抑制

定义：抑制剂作用方式。

公式：
- 竞争性 $K_m^{\text{app}}=K_m(1+[I]/K_i)$，$V_{\max}$ 不变。
- 非竞争性 $V_{\max}^{\text{app}}=V_{\max}/(1+[I]/K_i)$，$K_m$ 不变。

推导（竞争性）：抑制剂 $I$ 与酶结合 $E+I\rightleftharpoons EI$，$K_i=[E][I]/[EI]$。酶守恒 $[E]_0=[E]+[ES]+[EI]$。准稳态
$$
[ES]=\frac{[E][S]}{K_m}.
$$
又 $[EI]=[E][I]/K_i$，代入守恒：
$$
[E]_0=[E]\Big(1+\frac{[S]}{K_m}+\frac{[I]}{K_i}\Big).
$$
故
$$
[ES]=\frac{[E]_0[S]/K_m}{1+[S]/K_m+[I]/K_i}
=\frac{[E]_0[S]}{K_m(1+[I]/K_i)+[S]}.
$$
即 $K_m\to K_m(1+[I]/K_i)$，$V_{\max}$ 不变。非竞争性：$I$ 与 $ES$ 结合，得 $V_{\max}\to V_{\max}/(1+[I]/K_i)$。

物理意义：竞争性可被高底物浓度克服，非竞争性不能。几何意义：Lineweaver-Burk 图上竞争性绕纵截距旋转，非竞争性改变截距。

特例：$[S]\to\infty$ 时竞争性抑制消失；非竞争性始终降低 $V_{\max}$。

应用：药物设计、代谢调控、毒理学。

#### 13.4 Hill 方程与协同性

定义：多亚基酶/受体的协同结合。

假设：$n$ 个结合位点、完全协同（要么全结合要么全空）、Hill 系数 $n_H$。

公式：
$$
Y=\frac{[S]^{n_H}}{K_d+[S]^{n_H}}.
$$

推导：反应 $P+nS\rightleftharpoons PS_n$，平衡常数 $K_d=[P][S]^n/[PS_n]$。饱和分数
$$
Y=\frac{[PS_n]}{[P]+[PS_n]}=\frac{[S]^n/K_d}{1+[S]^n/K_d}
=\frac{[S]^n}{K_d+[S]^n}.
$$
Hill 图 $\ln(Y/(1-Y))=n_H\ln[S]-\ln K_d$，斜率为 $n_H$。

物理意义：$n_H>1$ 正协同（S 形、开关式响应），$n_H<1$ 负协同。几何意义：Hill 图斜率反映表观结合位点数。

特例：$n_H=1$ 得 MM 型双曲线；$n_H\to\infty$ 得阶跃开关。

应用：血红蛋白氧合、基因开关、受体信号。

#### 13.5 代谢通量平衡（通量平衡分析）

定义：稳态代谢网络中通量满足 $\mathbf S\mathbf v=0$。

假设：稳态、化学计量矩阵 $\mathbf S$、通量边界。

公式：
$$
\mathbf S\mathbf v=0,\qquad
\mathbf v_{\min}\le\mathbf v\le\mathbf v_{\max}.
$$

推导：代谢物浓度变化 $d\mathbf c/dt=\mathbf S\mathbf v$。稳态 $d\mathbf c/dt=0$ 给 $\mathbf S\mathbf v=0$。目标常为 $\max\mathbf c^T\mathbf v$（如生物量），约束为化学计量与边界，构成线性规划。

物理意义：稳态下各代谢物生成与消耗平衡。几何意义：可行通量集合是多面体（通量空间），$(\mathbf S\mathbf v=0)$ 是零空间。

特例：$\mathbf v$ 若在零空间且满足边界即稳态通量分布。

应用：代谢工程、药物靶点、通量组分析。

#### 13.6 代谢控制分析

定义：通量控制系数 $C^J_i$ 度量酶 $i$ 对稳态通量 $J$ 的相对影响。

假设：稳态、酶浓度可独立缩放、系统齐次。

公式：
$$
C^J_i=\frac{\partial\ln J}{\partial\ln E_i},\qquad
\sum_i C^J_i=1.
$$

推导：设所有酶浓度按 $E_i\to\lambda E_i$ 缩放，则各反应速率 $v_i\propto E_i$ 也按 $\lambda$ 缩放，稳态通量满足 $J(\lambda E)=\lambda J(E)$（一次齐次）。由 Euler 齐次函数定理
$$
\sum_i E_i\frac{\partial J}{\partial E_i}=J.
$$
两边除以 $J$：
$$
\sum_i\frac{\partial\ln J}{\partial\ln E_i}=1
\Rightarrow \sum_i C^J_i=1.
$$
浓度控制系数类似定义 $C^S_i=\partial\ln S/\partial\ln E_i$，并满足加和为零的约束。

物理意义：代谢通量由整个网络分担控制，单一“限速步”一般不存在。几何意义：控制系数是通量对酶浓度的对数导数向量，其和为 1。

特例：单一酶反应 $C^J_1=1$；多酶网络中控制可分散到多步。

应用：代谢工程靶点、药物作用位点、通量调控。

### 第14章 种群与传染病动力学

#### 14.1 Malthus 增长

定义：无资源限制下的指数增长。

公式：
$$
\frac{dN}{dt}=rN\ \Rightarrow\ N(t)=N_0e^{rt}.
$$

推导：分离变量
$$
\frac{dN}{N}=r\,dt\ \Rightarrow\ \ln N=rt+C.
$$
初值 $N(0)=N_0$ 得 $C=\ln N_0$，故 $N=N_0e^{rt}$。

物理意义：$r>0$ 增长、$r<0$ 衰减。几何意义：半对数坐标下为直线。

特例：$r=0$ 恒常；$r<0$ 指数衰减。

应用：早期种群扩张、细菌培养、疫情初期。

#### 14.2 Logistic 增长

定义：受环境容量 $K$ 限制的增长。

公式：
$$
\frac{dN}{dt}=rN\Big(1-\frac NK\Big),\qquad
N(t)=\frac{K}{1+\big(\frac{K}{N_0}-1\big)e^{-rt}}.
$$

推导：方程
$$
\frac{dN}{N(1-N/K)}=r\,dt.
$$
部分分式
$$
\frac{1}{N(1-N/K)}=\frac1N+\frac{1/K}{1-N/K}.
$$
积分：
$$
\ln N-\ln(1-N/K)=rt+C
\Rightarrow \frac{N}{1-N/K}=Ce^{rt}.
$$
初值给 $C=\frac{N_0}{1-N_0/K}$。解出
$$
N(t)=\frac{K}{1+\big(\frac{K}{N_0}-1\big)e^{-rt}}.
$$

物理意义：种群趋向容量 $K$，增长率为 S 形。几何意义：$N$-$t$ 曲线为 S 形（logistic），拐点在 $N=K/2$。

特例：$N\ll K$ 时指数增长；$N=K$ 时不增长。

应用：资源限制的生长、肿瘤增长、创新扩散。

#### 14.3 Lotka–Volterra 捕食-被捕食模型

定义：两种群相互作用的周期振荡模型。

公式：
$$
\dot x=\alpha x-\beta xy,\qquad
\dot y=\delta xy-\gamma y.
$$

推导：被捕食者 $x$ 增长 $\alpha$，被捕食率 $\beta xy$；捕食者 $y$ 摄入 $\delta xy$，死亡 $\gamma y$。平衡点 $(\gamma/\delta,\alpha/\beta)$。构造
$$
H(x,y)=\delta(x-x^*\ln x)+\beta(y-y^*\ln y),
$$
则
$$
\frac{dH}{dt}=\delta\Big(1-\frac{x^*}{x}\Big)\dot x+\beta\Big(1-\frac{y^*}{y}\Big)\dot y=0,
$$
故轨道为闭曲线，解为周期振荡。线性化给出特征值 $\pm i\sqrt{\alpha\gamma}$，中性稳定中心。

物理意义：捕食者-被捕食者数量周期性震荡，互不消灭。几何意义：相平面闭轨道，$H$ 为守恒量。

特例：施加捕食者-被捕食者之外的影响可移动中心；加入容量限制可得稳定螺旋（阻尼振荡）。

应用：生态周期、渔业模型、宿主-寄生虫。

#### 14.4 SIR 传染病模型

定义：易感-感染-康复三室模型。

假设：均匀混合、康复后免疫、感染率 $\beta$、恢复率 $\gamma$。

公式：
$$
\dot S=-\frac{\beta SI}{N},\quad
\dot I=\frac{\beta SI}{N}-\gamma I,\quad
\dot R=\gamma I.
$$

推导：每单位时间新感染数 $\beta SI/N$（质量作用），康复 $\gamma I$。对 $S,I,R$ 写微分方程。由 $\dot S<0$ 且 $\dot R>0$。初值 $S(0)=N-1,I(0)=1$。

阈值：感染增长条件 $\dot I>0\iff \beta S/N-\gamma>0$，定义
$$
R_0=\frac{\beta}{\gamma},\qquad
S_c=\frac{\gamma}{\beta}N=\frac{N}{R_0}.
$$
若 $S(0)>S_c$，疫情爆发。

物理意义：$R_0$ 为基本传染数。几何意义：$(S,I)$ 相平面上的阈值线。

特例：$R_0<1$ 疫情消退；$R_0>1$ 疫情扩散。

应用：流感/COVID 防控、疫苗覆盖率。

#### 14.5 群体免疫阈值

定义：使 $R_{\text{eff}}=R_0 S/N<1$ 所需免疫比例。

公式：
$$
p_c=1-\frac1{R_0}.
$$

推导：疫苗接种比例 $p$ 后易感比例变为 $S/N=1-p$，有效传染数
$$
R_{\text{eff}}=R_0(1-p).
$$
要求 $R_{\text{eff}}<1$：
$$
R_0(1-p)<1\ \Rightarrow\ p>1-\frac1{R_0}.
$$
故最小群体免疫比例 $p_c=1-1/R_0$。

物理意义：疫苗通过减少易感者阻断传播链。几何意义：$R_{\text{eff}}=1$ 对应临界免疫线。

特例：$R_0=2$ 时 $p_c=0.5$；$R_0\to\infty$ 时 $p_c\to1$。

应用：疫苗策略、公共卫生阈值。

#### 14.6 平衡点稳定性分析

定义：对 $\dot{\mathbf x}=\mathbf f(\mathbf x)$ 在平衡点 $\mathbf x^*$ 附近线性化。

公式：
$$
\dot{\boldsymbol\xi}=\mathbf J\boldsymbol\xi,\qquad
\mathbf J=\frac{\partial\mathbf f}{\partial\mathbf x}\Big|_{\mathbf x^*},\qquad
\lambda_i=\text{eig}(\mathbf J).
$$

推导：令 $\mathbf x=\mathbf x^*+\boldsymbol\xi$，Taylor 展开
$$
\dot{\boldsymbol\xi}=\mathbf f(\mathbf x^*)+\mathbf J\boldsymbol\xi+O(\|\boldsymbol\xi\|^2)=\mathbf J\boldsymbol\xi+O(\|\boldsymbol\xi\|^2).
$$
解为 $\boldsymbol\xi(t)=\sum c_ie^{\lambda_it}\mathbf v_i$。若 $\mathrm{Re}\,\lambda_i<0$ 全部，稳定；存在 $\mathrm{Re}\,\lambda_i>0$ 则不稳定；纯虚根为临界（需非线性项判断）。

物理意义：稳定性由 Jacobian 特征值实部决定。几何意义：特征向量给出扰动指数增长/衰减方向。

特例：一维 $\dot x=f'(x^*)(x-x^*)$，$f'<0$ 稳定。

应用：生态平衡、代谢稳态、流行病无病平衡。

### 第15章 生物系统的稳定性与分岔

#### 15.1 Lyapunov 稳定性

定义：平衡点稳定若任意邻域内初值轨道保持在更小邻域内。

假设：自治系统、平衡点 $\mathbf x^*=0$。

公式：若存在正定函数 $V(\mathbf x)$ 使 $\dot V\le0$，则稳定；若 $\dot V<0$（除原点）则渐近稳定。

推导：$\dot V=\nabla V\cdot\mathbf f(\mathbf x)$。考虑球 $B_\epsilon$。取 $V$ 的等值面。若 $\dot V\le0$，轨道不能穿越向外的 $V$ 等值面，从而被限制在 $B_\epsilon$ 内，故稳定。若 $\dot V<0$，$V$ 严格递减且有下界，必收敛到不变集；若该集仅含原点，则渐近稳定。

物理意义：$V$ 类似“能量”函数，单调下降保证稳定。几何意义：$V$ 的等值面为嵌套闭曲面，轨道向内。

特例：$V$ 可取二次型 $V=\mathbf x^T\mathbf P\mathbf x$（线性系统）。

应用：种群平衡、代谢网络稳定、控制系统。

#### 15.2 Hopf 分岔

定义：平衡点由稳定变为不稳定并产生极限环。

假设：一维参数 $\mu$，一对复特征值 $\lambda(\mu)=\alpha(\mu)\pm i\omega(\mu)$，在 $\mu_c$ 处 $\alpha(\mu_c)=0,\omega\ne0$，横截条件 $\alpha'(\mu_c)\ne0$。

公式：极限环在 $\mu_c$ 附近振幅
$$
A\sim\sqrt{\frac{\mu-\mu_c}{-a}}\quad(\text{超临界}).
$$

推导：中心流形上约化为复方程
$$
\dot z=\alpha(\mu)z+a|z|^2z+O(|z|^4).
$$
设 $z=re^{i\phi}$，则
$$
\dot r=\alpha(\mu)r+a r^3,\quad \dot\phi=\omega(\mu).
$$
非零平衡 $r^2=-\alpha(\mu)/a$。若 $a<0$，当 $\alpha$ 由负变正时 $r>0$ 出现，即超临界 Hopf 分岔，极限环稳定。振幅 $A=\sqrt{-\alpha/a}\sim\sqrt{\mu-\mu_c}$。

物理意义：生物节律、神经振荡、种群周期由此产生。几何意义：相平面中平衡点失稳后长出极限环。

特例：$a>0$ 为亚临界，环不稳定且跳变。

应用：生物钟、心血管振荡、神经元放电。

#### 15.3 Turing 不稳定性

定义：反应-扩散系统中扩散驱动的不均匀斑图。

假设：两组分、均匀稳态稳定（无扩散）、扩散系数不同。

公式：
$$
\partial_t u=D_u\nabla^2u+f(u,v),\qquad
\partial_t v=D_v\nabla^2v+g(u,v).
$$

推导：均匀稳态 $(u^*,v^*)$，Jacobian
$$
\mathbf J=\begin{pmatrix}f_u&f_v\\g_u&g_v\end{pmatrix}.
$$
无扩散稳定条件 $\mathrm{tr}\mathbf J<0,\ \det\mathbf J>0$。加入波数 $k$ 的扰动 $e^{\lambda t+i kx}$，色散关系
$$
\det\begin{pmatrix}\lambda+D_uk^2-f_u&-f_v\\-g_u&\lambda+D_vk^2-g_v\end{pmatrix}=0.
$$
即
$$
\lambda^2+T_k\lambda+D_k=0,
$$
$$
T_k=(D_u+D_v)k^2-\mathrm{tr}\mathbf J,
$$
$$
D_k=D_uD_vk^4-(D_vf_u+D_ug_v)k^2+\det\mathbf J.
$$
Turing 不稳定性要求 $D_k<0$ 对某 $k$，即 $D_vf_u+D_ug_v>0$ 且判别式为正。临界波数
$$
k_c^2=\frac{D_vf_u+D_ug_v}{2D_uD_v}.
$$

物理意义：短程激活、长程抑制的扩散差异产生空间斑图。几何意义：色散曲线 $\lambda(k)$ 在某 $k_c$ 处变正。

特例：$D_u=D_v$ 时不可能（需差分扩散）；无反应项（$f,g$ 单调）不产生。

应用：动物皮毛、体节形成、菌落图案。

#### 15.4 双稳态与开关

定义：系统存在两个稳定平衡点。

假设：一维 $ \dot x=f(x)$，$f$ 呈三次型。

公式：
$$
\dot x=\alpha-x+x^3/3\ \text{(规范型)}.
$$

推导：平衡满足 $f(x)=0$。对 $f=\alpha-x+x^3/3$，$f'(x)=-1+x^2$。平衡稳定性由 $f'$ 符号：$|x|>1$ 稳定，$|x|<1$ 不稳定。三个平衡存在条件 $|\alpha|<2/3$。分岔点 $\alpha=\pm2/3$ 为鞍结分岔，成对平衡点产生/消失。

物理意义：细胞命运决定、基因开关的双稳切换。几何意义：$f(x)$ 与横轴交点，稳定点对应负斜率。

特例：$\alpha=0$ 三个平衡 $0,\pm\sqrt3$；$|\alpha|>2/3$ 单平衡。

应用：细胞分化、记忆、合成生物学开关。

#### 15.5 Kuramoto 同步与临界耦合

定义：大量相位振子通过正弦耦合趋向同步。

假设：平均场耦合、自然频率分布 $g(\omega)$、弱耦合。

公式：
$$
\frac{d\theta_i}{dt}=\omega_i+\frac{K}{N}\sum_{j=1}^N\sin(\theta_j-\theta_i),
$$
$$
re^{i\psi}=\frac1N\sum_{j=1}^N e^{i\theta_j},\qquad K_c=\frac{2}{\pi g(0)}.
$$

推导：以序参量改写
$$
re^{i\psi}=\frac1N\sum_j e^{i\theta_j}
\Rightarrow \sum_j\sin(\theta_j-\theta_i)=\mathrm{Im}\!\Big(e^{-i\theta_i}\sum_j e^{i\theta_j}\Big)=r\sin(\psi-\theta_i).
$$
故 $\dot\theta_i=\omega_i+Kr\sin(\psi-\theta_i)$。设锁相振子满足 $\dot\theta_i=\psi$，由 $\psi=\Omega$ 得 $\sin(\psi-\theta_i)=(\Omega-\omega_i)/(Kr)$，可锁相条件 $|\omega_i-\Omega|\le Kr$。自洽条件给出
$$
r=\int_{|\omega-\Omega|\le Kr}\cos\theta(\omega)g(\omega)\,d\omega.
$$
在 $r\to0^+$ 展开，临界耦合满足 $K_c\int g(\omega)d\omega=2$ 即 $K_c=2/(\pi g(0))$。

物理意义：耦合强度超过阈值后宏观同步涌现。几何意义：序参量 $r$ 是单位圆上相位点的质心模长，$r=1$ 完全同步。

特例：$K<K_c$ 时 $r=0$ 无序；$K\to\infty$ 时全锁相。

应用：生物钟同步、萤火虫闪光、神经放电节律。

#### 15.6 噪声诱导的双稳态跃迁

定义：双稳势中热噪声驱动系统在两个稳态间随机切换。

假设：过阻尼、势垒 $\Delta U\gg D$、$D=k_BT$。

公式：
$$
k_{\text{switch}}\propto e^{-\Delta U/D},\qquad
\langle T_{\text{dwell}}\rangle\propto e^{\Delta U/D}.
$$

推导：一维过阻尼 Fokker-Planck 方程在双势阱 $U(x)$ 中，准稳态概率流可写
$$
J\approx\frac{D\,p(x_a)\,e^{-U(x_b)/D}}{\int_{x_a}^{x_b}e^{-U/D}dx}\times\cdots
$$
用阱底与垒顶的二次近似（曲率频率 $\omega_a,\omega_b$），整理得 Kramers 律：
$$
k=\frac{\omega_a\omega_b}{2\pi\xi}e^{-\Delta U/D},\qquad \Delta U=U(x_b)-U(x_a).
$$
平均驻留时间 $\langle T\rangle=1/k$ 与 $e^{\Delta U/D}$ 成正比。

物理意义：切换速率对势垒高度/噪声强度指数敏感。几何意义：$\Delta U$ 是双阱能量景观的垒高，噪声在景观上做随机翻越。

特例：$D\to0$ 时不切换；$\Delta U=0$ 时两稳态合并。

应用：基因开关切换、细胞命运跃迁、记忆保持。

### 第16章 生物信息与统计物理基础

#### 16.1 Boltzmann 分布

定义：平衡态下微观态 $i$ 的概率。

假设：正则系综、温度 $T$、能量 $E_i$。

公式：
$$
p_i=\frac{e^{-E_i/k_BT}}{Z},\qquad
Z=\sum_ie^{-E_i/k_BT}.
$$

推导：由最大熵原理：在约束 $\sum p_i=1$ 与 $\sum p_iE_i=\langle E\rangle$ 下最大化 $S=-k_B\sum p_i\ln p_i$。构造 Lagrange 函数
$$
\mathcal L=-k_B\sum p_i\ln p_i-\alpha\Big(\sum p_i-1\Big)-\beta\Big(\sum p_iE_i-\langle E\rangle\Big).
$$
对 $p_i$ 求导置零：
$$
-k_B(\ln p_i+1)-\alpha-\beta E_i=0
\Rightarrow p_i=e^{-1-\alpha/k_B}e^{-\beta E_i/k_B}.
$$
归一化给出 $p_i=e^{-\beta E_i}/Z$。由热力学温度关系 $\beta=1/(k_BT)$（通过 $dS/d\langle E\rangle=1/T$ 可证）。

物理意义：高能态指数稀少。几何意义：概率在能量轴上指数衰减。

特例：$T\to0$ 只占据基态；$T\to\infty$ 各态等概率。

应用：离子通道态分布、蛋白构象、配体结合。

#### 16.2 配分函数与自由能

定义：$Z=\sum_ie^{-\beta E_i}$，Helmholtz 自由能 $F=-k_BT\ln Z$。

公式：
$$
F=-k_BT\ln Z,\qquad
\langle E\rangle=-\frac{\partial\ln Z}{\partial\beta}.
$$

推导：由 $p_i=e^{-\beta E_i}/Z$，$\ln Z=\ln e^{-\beta E_i}-\ln p_i$ 平均：
$$
\ln Z=-\beta\langle E\rangle-S/k_B.
$$
故
$$
-k_BT\ln Z=\langle E\rangle-TS=F.
$$
又
$$
-\frac{\partial\ln Z}{\partial\beta}
=-\frac{1}{Z}\sum_i(-E_i)e^{-\beta E_i}
=\sum_iE_ip_i=\langle E\rangle.
$$

物理意义：$Z$ 编码全部热力学信息。几何意义：$\ln Z$ 是能量谱的 Laplace 变换。

特例：单能级 $E$，$Z=e^{-\beta E}$，$F=E$。

应用：生化平衡、受体占据、生物膜相变。

#### 16.3 熵与信息

定义：Shannon 熵 $H=-\sum_ip_i\log_2p_i$（bit）。

假设：离散概率分布、$p_i>0$。

公式：
$$
S=k_B\ln2\cdot H\ \text{(热力学熵与信息熵关系)}.
$$

推导：Boltzmann 熵 $S=k_B\ln\Omega$。当各微观态等概率 $p_i=1/\Omega$，$H=\log_2\Omega$，故 $S=k_B\ln2\,H$。一般分布下，Gibbs 熵 $S=-k_B\sum p_i\ln p_i$ 与 Shannon 熵差常数因子。

物理意义：信息是熵的负值（负熵），测量/擦除信息有热力学代价。几何意义：熵是概率分布的“不确定性”泛函。

特例：确定事件 $H=0$；均匀分布 $H$ 最大。

应用：生物信息传递、最大熵推断、神经编码。

#### 16.4 最大熵分布

定义：在给定约束下熵最大的分布。

假设：约束 $\langle E\rangle$、$\sum p_i=1$。

公式：一般形式
$$
p_i=\frac{1}{Z}e^{-\sum_k\lambda_k f_k(x_i)}.
$$

推导：同 16.1 的 Lagrange 方法，对多个约束 $\langle f_k\rangle$ 引入 $\lambda_k$，最大化 $-\sum p_i\ln p_i$，得
$$
\ln p_i=-1-\sum_k\lambda_kf_k(x_i)
\Rightarrow p_i\propto e^{-\sum_k\lambda_kf_k(x_i)}.
$$
归一化得 $Z$。

物理意义：给定宏观信息时，最大熵分布是“最无偏”的推断。几何意义：约束集的切空间中最速下降熵方向为零。

特例：只约束均值得指数分布；约束均值与方差得高斯分布。

应用：蛋白接触图谱、神经群体解码、网络推断。

#### 16.5 Landauer 原理

定义：擦除一比特信息至少耗散 $k_BT\ln2$ 热量。

假设：双势阱存储一比特、等概率初始、准静态擦除。

公式：
$$
\Delta Q\ge k_BT\ln2.
$$

推导：擦除把双势阱分布由 $p_0=p_1=1/2$ 变为确定态（如全在 $p_0$）。熵变
$$
\Delta S=k_B\ln2 - 0=k_B\ln2\ \text{(系统熵减少)}.
$$
总熵（系统+环境）不得减少，故环境熵增至少 $k_B\ln2$，对应热耗散 $\Delta Q=T\Delta S=k_BT\ln2$。

物理意义：信息处理有基本热力学极限。几何意义：信息的“擦除”对应相空间体积压缩，须由环境补偿。

特例：可逆计算可低于该极限（渐进接近零耗散）。

应用：分子计算、生物信息处理能耗、DNA 存储。

#### 16.6 基因表达噪声

定义：mRNA/蛋白拷贝数在单细胞中的随机涨落。

假设：生灭过程、转录/翻译为泊松事件、无反馈。

公式：
$$
\frac{d\langle n\rangle}{dt}=k-\gamma\langle n\rangle,\qquad
\langle n\rangle=\frac{k}{\gamma},
$$
简单生灭：$\mathrm{Var}(n)=\langle n\rangle$（Fano 因子 $F=1$）；爆发式表达：$F=1+b$。

推导：主方程
$$
\frac{dp_n}{dt}=k p_{n-1}+\gamma(n+1)p_{n+1}-(k+\gamma n)p_n.
$$
乘 $n$ 并对 $n$ 求和（用 $\sum n p_{n-1}=\langle n\rangle+1$、$\sum n(n+1)p_{n+1}=\langle n^2\rangle-\langle n\rangle$）得
$$
\frac{d\langle n\rangle}{dt}=k-\gamma\langle n\rangle.
$$
乘 $n^2$ 求和得
$$
\frac{d\langle n^2\rangle}{dt}=k(2\langle n\rangle+1)+\gamma(\langle n\rangle-2\langle n^2\rangle).
$$
稳态 $d\langle n\rangle/dt=0$ 给 $k=\gamma\langle n\rangle$，代入第二式：
$$
2\gamma\langle n^2\rangle=2k\langle n\rangle+k+\gamma\langle n\rangle
=2\gamma\langle n\rangle^2+2\gamma\langle n\rangle,
$$
故 $\langle n^2\rangle=\langle n\rangle^2+\langle n\rangle$，$\mathrm{Var}=\langle n\rangle$。若转录成簇（平均 $b$ 个 mRNA），方差变为 $\mathrm{Var}=\langle n\rangle(1+b)$，Fano 因子 $F=1+b>1$。

物理意义：小拷贝数基因有显著噪声，爆发式表达放大约束噪声。几何意义：稳态分布为 Poisson（喷发时为负二项），$F$ 度量相对泊松的超散。

特例：$b=0$ 退化为 Poisson；$b\to\infty$ 时噪声发散。

应用：单细胞异质性、合成基因回路、细胞命运决策。

---
⬆ 上级领域：[[物质科学（物理·化学·生命）]]
