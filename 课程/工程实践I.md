---
课程名称: 工程实践I
学分: 1
总学时: 24
考核方式: 考试
课程属性: 必修
课程性质: 实践环节
开课单位: 空天科学与技术学院
领域: 实践与创新
tags:
  - 课程
  - 领域/实践与创新
---

# 工程实践I


### 第1章 量纲分析与相似理论

**知识点1：基本量与量纲的定义**

定义：物理量 $Q$ 用其数值 $\{Q\}$ 与单位 $[Q]$ 的乘积表示，$Q=\{Q\}\,[Q]$。量纲是物理量对基本量的幂次构成，记 $\dim Q$。力学问题取长度、质量、时间、温度四个基本量，$\dim L=\mathsf{L}$、$\dim m=\mathsf{M}$、$\dim t=\mathsf{T}$、$\dim T=\Theta$。

假设：所研究物理定律在任意单位制下形式不变（单位无关性）。

公式：若 $Q$ 由基本量表达，则 $\dim Q=\mathsf{L}^{\alpha}\mathsf{M}^{\beta}\mathsf{T}^{\gamma}\Theta^{\delta}$。

推导：设两组基本单位缩放因子为 $\lambda_L,\lambda_M,\lambda_T,\lambda_\Theta$，则导出量单位按 $\lambda_Q=\lambda_L^{\alpha}\lambda_M^{\beta}\lambda_T^{\gamma}\lambda_\Theta^{\delta}$ 缩放。要求物理规律 $f(Q_1,\dots,Q_n)=0$ 在缩放下保持，幂指数 $(\alpha,\beta,\gamma,\delta)$ 必须与缩放因子无关，故量纲定义为固定指数向量。

物理意义：量纲反映物理量的内在属性，与单位无关。极限情形：无量纲量 $\alpha=\beta=\gamma=\delta=0$，其数值不随单位改变。

应用：判断公式是否可能正确；由 $\dim E=\mathsf{M}\mathsf{L}^2\mathsf{T}^{-2}$ 验证动能公式。

**知识点2：量纲齐次性原理**

定义：物理方程中各项量纲必须相同。

假设：物理规律可由幂函数型关系表达。

公式：若 $\sum_i Q_i=0$，则 $\dim Q_i$ 全等。

推导：设方程含两项 $A$ 与 $B$，$A+B=0$。单位缩放后 $\lambda_A A'+\lambda_B B'=0$ 须对任意 $\lambda$ 成立，令 $\lambda_A\ne\lambda_B$ 则只有 $A=B=0$，矛盾，故 $\lambda_A=\lambda_B$，即 $\dim A=\dim B$。

物理意义：只有同量纲量才能相加减、比较大小。

特例：$e^x$、$\ln x$、$\sin x$ 的宗量必为无量纲，否则 Taylor 展开各项量纲不一致。典型应用：由 $v=at$ 检查 $\mathsf{L}\mathsf{T}^{-1}=\mathsf{L}\mathsf{T}^{-2}\cdot\mathsf{T}$。

**知识点3：白金汉π定理**

定义：若问题涉及 $n$ 个物理量、含 $k$ 个独立基本量，则存在 $m=n-k$ 个独立的无量纲组合。

假设：物理关系 $\Phi(Q_1,\dots,Q_n)=0$ 可由这些量描述。

公式：$\Phi(\pi_1,\dots,\pi_m)=0$，其中 $\pi_j=\prod_i Q_i^{a_{ij}}$。

推导：将量纲矩阵 $A=[a_{ij}]$（$n$ 列、$k$ 行）写出。$\pi=\prod Q_i^{c_i}$ 无量纲要求 $A\mathbf c=0$。该齐次方程解空间维数 $\dim\ker A=n-\operatorname{rank}A=n-k=m$。取解空间一组基 $\mathbf c^{(1)},\dots,\mathbf c^{(m)}$，即得 $\pi_j$。

物理意义：把高维参数空间压缩为低维无量纲空间，变量数从 $n$ 降到 $n-k$。

特例：$k=n$ 时 $m=0$，无自由无量纲组合，问题被完全确定。

应用：单摆周期 $T_p$ 由 $l,m,g$ 决定，$n=4,k=3,m=1$，得 $\pi=T_p\sqrt{g/l}$，故 $T_p=C\sqrt{l/g}$。

**知识点4：相似三定理**

定义：相似第一定理（正定理）——现象相似则相似准则数相等；第二定理（π定理）——物理关系可写成准则数关系；第三定理（逆定理）——单值条件相似且准则数相等则现象相似。

假设：两系统几何、物理条件可对应。

公式：$\dfrac{Q_{1i}}{Q_{1j}}=\text{const}$（同类量相似比），且 $\pi_{1k}=\pi_{2k}$。

推导（正定理）：设原型与模型对应量成比例 $Q_i'=C_iQ_i$。代入同一物理方程，为使方程形式不变，各 $C_i$ 必须满足与量纲方程相同的关系，即 $C_{Q}=\prod C_{\text{base}}^{a}$，从而 $Q_i'/\prod(\text{base}')^{a}=Q_i/\prod(\text{base})^{a}$，即准则数相等。

物理意义：相似包括几何、运动、动力、热相似。

应用：风洞实验设计保证雷诺数、马赫数相等。

**知识点5：雷诺数、弗劳德数、马赫数等常用准则数**

定义：雷诺数 $Re=\rho UL/\mu=UL/\nu$；弗劳德数 $Fr=U/\sqrt{gL}$；马赫数 $Ma=U/c$；欧拉数 $Eu=\Delta p/(\rho U^2)$。

推导（$Re$）：由 N-S 方程对流域 $L$、速度 $U$、密度 $\rho$、黏度 $\mu$ 无量纲化，惯性项 $\rho U^2/L$ 与黏性项 $\mu U/L^2$ 之比为 $\rho UL/\mu$。

物理意义：$Re$ 表征惯性力与黏性力之比；$Fr$ 表征惯性力与重力之比；$Ma$ 表征流速与声速之比（可压缩性）。

特例：$Re\ll1$ 蠕流（斯托克斯流）；$Re\gg1$ 边界层/湍流。应用：管流层湍转捩约 $Re\approx2300$。

**知识点6：方程无量纲化与特征尺度**

定义：用特征量 $L,U,\rho,\Delta T$ 将变量化为无量纲：$x^*=x/L$，$u^*=u/U$，$p^*=p/(\rho U^2)$。

假设：特征尺度选取使得各项有界。

公式：$\partial_t^* \mathbf u^*+(\mathbf u^*\cdot\nabla^*)\mathbf u^*=-\nabla^*p^*+\dfrac1{Re}\nabla^{*2}\mathbf u^*$。

推导：将 $\nabla=(1/L)\nabla^*$、$\partial_t=(U/L)\partial_t^*$ 代入 N-S 方程 $\rho(\partial_t\mathbf u+\mathbf u\cdot\nabla\mathbf u)=-\nabla p+\mu\nabla^2\mathbf u$，两端同乘 $L/(\rho U^2)$ 得上述无量纲式。

物理意义：无量纲方程中参数个数减少，主导项由 $Re$ 决定。

极限：$Re\to\infty$ 黏性项消失为欧拉方程；$Re\to0$ 惯性项消失为斯托克斯方程。

**知识点7：相似准则的物理组合与模型律**

定义：模型律给出原型与模型之间物理量的缩放关系。

推导：由 $Re$ 相等 $\dfrac{U_pL_p}{\nu_p}=\dfrac{U_mL_m}{\nu_m}$ 得速度比 $\dfrac{U_m}{U_p}=\dfrac{L_p}{L_m}\dfrac{\nu_m}{\nu_p}$。由 $Fr$ 相等得 $\dfrac{U_m}{U_p}=\sqrt{\dfrac{L_m}{L_p}}$。二者同时要求时给出一致性条件 $\dfrac{\nu_m}{\nu_p}=\left(\dfrac{L_m}{L_p}\right)^{3/2}$。

物理意义：多数工程无法同时满足全部准则，须按主导相似选取。

应用：水池/水洞按 $Fr$ 缩比；低速风洞按 $Re$ 缩比。

**知识点8：量纲分析在实验设计中的作用**

定义：通过准则数规划实验点、减少试验次数。

公式：待定函数关系 $\pi_0=f(\pi_1,\dots,\pi_{m})$。

推导：若忽略化为 $f\approx C\prod\pi_i^{a_i}$，取对数 $\ln\pi_0=\ln C+\sum a_i\ln\pi_i$，成为线性回归问题，可用最小二乘求 $a_i$。

物理意义：把多变量问题化为少数无量纲变量，使实验规模由 $n$ 维降至 $m$ 维。

应用：摩擦因子 $f=f(Re,\varepsilon/D)$（莫迪图）；确定经验公式指数。

### 第2章 静力学与刚体平衡

**知识点1：力、力矩与力偶**

定义：力 $\mathbf F$ 是矢量；对点 $O$ 的力矩 $\mathbf M_O=\mathbf r\times\mathbf F$；力偶是等值反向不共线两力，其合力为零、合力矩 $\mathbf M=\mathbf r\times\mathbf F$ 与矩心无关。

假设：刚体、力为滑移矢量（作用于刚体上可沿作用线滑移）。

推导（力偶矩与矩心无关）：设力偶为 $(\mathbf F,-\mathbf F)$ 分别作用于 $\mathbf r_1,\mathbf r_2$。对任意点 $O$ 的矩
$$
\mathbf M_O=(\mathbf r_1-\mathbf r_O)\times\mathbf F+(\mathbf r_2-\mathbf r_O)\times(-\mathbf F)=(\mathbf r_1-\mathbf r_2)\times\mathbf F,
$$
与 $\mathbf r_O$ 无关。

物理意义：力偶只产生转动效应。

特例：$\mathbf F\parallel\mathbf r$ 时 $\mathbf M=0$。应用：方向盘力矩。

**知识点2：力系简化与主矢主矩**

定义：力系主矢 $\mathbf R=\sum\mathbf F_i$，对 $O$ 的主矩 $\mathbf M_O=\sum\mathbf r_i\times\mathbf F_i$。

推导（平移定理）：作用于 $\mathbf r_i$ 的 $\mathbf F_i$ 等效为作用于 $O$ 的 $\mathbf F_i$ 加力偶 $\mathbf r_i\times\mathbf F_i$。求和得等效于 $(O,\mathbf R)$ 加力偶 $\mathbf M_O$。

物理意义：任意力系可化为一个力加一个力偶（力螺旋）。

特例：$\mathbf R=0,\mathbf M_O\ne0$ 为纯力偶；$\mathbf R\cdot\mathbf M_O=0$ 可进一步化为单力。应用：结构整体受力分析。

**知识点3：刚体平衡方程**

定义：平衡条件为 $\mathbf R=0,\ \mathbf M_O=0$。

假设：刚体、静力平衡。

公式（平面）：$\sum F_x=0,\ \sum F_y=0,\ \sum M_O=0$。

推导：刚体不发生平动要求 $\mathbf R=0$；不发生转动要求 $\mathbf M_O=0$。二者相互独立，故三维有六个独立标量方程 $\sum F_x=\sum F_y=\sum F_z=0$，$\sum M_x=\sum M_y=\sum M_z=0$。

物理意义：静定问题未知力数等于独立方程数。

特例：汇交力系 $M$ 自动为零；平面力系仅三个独立方程。应用：支座反力求解。

**知识点4：约束与约束力**

定义：约束限制自由度并产生约束反力。

假设：光滑接触、理想约束（反力做功为零）。

推导（二力杆）：杆两端铰接受力、杆自重忽略，平衡要求两端反力共线反向，故沿杆轴线。

物理意义：约束反力方向由约束类型决定。

特例：固定端阻止平动与转动，约束反力含两个分量与一个力偶。应用：梁支座、铰链。

**知识点5：桁架节点法与截面法**

定义：桁架由二力杆铰接，节点受力汇交。

假设：节点为光滑铰、载荷作用于节点、杆自重忽略。

推导（节点法）：对每个节点列 $\sum F_x=0,\sum F_y=0$。设节点数 $j$、杆数 $b$、支座反力数 $r$，静定条件 $b+r=2j$。

物理意义：节点法逐杆求内力；截面法对部分桁架列 $\sum M=0$ 快速求特定杆。

特例：零杆判断规则——两杆节点无载荷则两杆为零杆（若不共线）。应用：桥梁桁架。

**知识点6：摩擦静力学**

定义：静摩擦 $F_s\le\mu_sN$，动摩擦 $F_k=\mu_kN$，$\mu_s\ge\mu_k$。

假设：库仑摩擦，接触面法向力 $N$。

推导（斜面滑块临界）：滑块重 $W$，斜面倾角 $\theta$。沿面平衡 $W\sin\theta=F_s$，法向 $N=W\cos\theta$。临界 $F_s=\mu_sN$ 得 $\tan\theta_c=\mu_s$，摩擦角 $\phi_s=\arctan\mu_s$。

物理意义：$\theta\le\phi_s$ 时自锁。

特例：$\mu_s=0$ 无摩擦。应用：螺纹自锁、制动。

**知识点7：重心与形心**

定义：重心 $\mathbf r_G=\dfrac{\sum m_i\mathbf r_i}{\sum m_i}=\dfrac{\int\mathbf r\,dm}{\int dm}$。

推导：重力合力矩等于各微元重力矩之和，$\left(\sum m_i\mathbf g\right)\times\mathbf r_G=\sum(m_i\mathbf g\times\mathbf r_i)$，消去 $\mathbf g$ 得定义式。

物理意义：重力作用点。

特例：均质物体重心与形心重合，$\mathbf r_c=\dfrac{\int\mathbf r\,dV}{V}$。应用：组合截面形心。

**知识点8：虚位移原理求解静力**

定义：理想约束下刚体平衡的充要条件是所有主动力在任意虚位移上做功之和为零 $\delta W=\sum\mathbf F_i\cdot\delta\mathbf r_i=0$。

假设：理想约束、微小虚位移。

推导：平衡时每个质点 $\mathbf F_i+\mathbf N_i=0$，乘虚位移求和 $\sum(\mathbf F_i+\mathbf N_i)\cdot\delta\mathbf r_i=0$。理想约束 $\sum\mathbf N_i\cdot\delta\mathbf r_i=0$，故 $\sum\mathbf F_i\cdot\delta\mathbf r_i=0$。反之由虚位移任意性得平衡。

物理意义：用功代替矢量平衡，便于多约束系统。

特例：单自由度系统 $\delta W=Q\delta q=0\Rightarrow Q=0$。应用：机构平衡、杠杆。

### 第3章 应力、应变与广义胡克定律

**知识点1：应力张量与平衡微分方程**

定义：面元 $\Delta A$ 上内力 $\Delta\mathbf F$，应力矢量 $\mathbf t^{(n)}=\lim_{\Delta A\to0}\Delta\mathbf F/\Delta A$；应力张量 $\boldsymbol\sigma$ 使 $\mathbf t^{(n)}=\boldsymbol\sigma\cdot\mathbf n$。

假设：连续介质、Cauchy 应力（仅依赖法向 $\mathbf n$）。

推导（应力张量对称性）：对微元四面体关于 $O$ 取矩，体力矩与惯性力矩为高阶小量，面力矩之和为零，得 $\sigma_{ij}=\sigma_{ji}$。

推导（平衡微分方程）：对微元体 $dV=dxdydz$ 沿 $x$ 列力平衡，$x$ 面应力 $\sigma_{xx}$ 及增量，得
$$
\left(\sigma_{xx}+\frac{\partial\sigma_{xx}}{\partial x}dx\right)dydz-\sigma_{xx}dydz+\cdots+f_x\,dxdydz=0,
$$
化简 $\dfrac{\partial\sigma_{xx}}{\partial x}+\dfrac{\partial\sigma_{xy}}{\partial y}+\dfrac{\partial\sigma_{xz}}{\partial z}+f_x=0$，即 $\nabla\cdot\boldsymbol\sigma+\mathbf f=0$。

物理意义：应力张量对称，6 个独立分量；平衡方程联系应力与体力。

特例：无体力 $\nabla\cdot\boldsymbol\sigma=0$。应用：弹性力学基本方程之一。

**知识点2：应变张量与几何方程**

定义：位移场 $\mathbf u(\mathbf x)$，小应变张量 $\varepsilon_{ij}=\dfrac12(u_{i,j}+u_{j,i})$。

假设：小变形、连续。

推导：变形后微线段 $d\mathbf x'=d\mathbf x+d\mathbf u$，长度平方变化
$$
|d\mathbf x'|^2-|d\mathbf x|^2=2\varepsilon_{ij}dx_idx_j+O(|du|^2)
$$
（格林应变 $\varepsilon_{ij}=\frac12(u_{i,j}+u_{j,i}+u_{k,i}u_{k,j})$），小变形忽略二次项得几何方程。例如 $\varepsilon_{xx}=\partial u/\partial x$，$\gamma_{xy}=\partial u/\partial y+\partial v/\partial x$。

物理意义：正应变度量伸长，切应变度量角变形。

特例：刚体平动 $\mathbf u=\text{const}$ 应变为零；刚体转动 $\mathbf u=\boldsymbol\omega\times\mathbf r$ 应变为零。应用：应变片测量。

**知识点3：主应力与应力不变量**

定义：应力张量特征值 $\sigma$ 满足 $\det(\boldsymbol\sigma-\sigma\mathbf I)=0$。

推导：
$$
\det(\boldsymbol\sigma-\sigma\mathbf I)=-\sigma^3+I_1\sigma^2-I_2\sigma+I_3=0,
$$
其中 $I_1=\operatorname{tr}\boldsymbol\sigma$，$I_2=\frac12[(\operatorname{tr}\boldsymbol\sigma)^2-\operatorname{tr}\boldsymbol\sigma^2]$，$I_3=\det\boldsymbol\sigma$。因特征值不随坐标变换，$I_i$ 为不变量。

物理意义：主应力是剪应力为零方向上的正应力。

特例：二维 $\sigma_{1,2}=\dfrac{\sigma_x+\sigma_y}{2}\pm\sqrt{\left(\dfrac{\sigma_x-\sigma_y}{2}\right)^2+\tau_{xy}^2}$，最大剪应力 $\tau_{\max}=\dfrac{\sigma_1-\sigma_3}{2}$。应用：莫尔圆。

**知识点4：广义胡克定律**

定义：线弹性各向同性材料应力应变关系。

假设：小变形、线弹性、各向同性。

推导（由单向推广）：单向拉伸 $\varepsilon_x=\sigma_x/E$，横向 $\varepsilon_y=\varepsilon_z=-\nu\sigma_x/E$。叠加三向正应力
$$
\varepsilon_x=\frac1E[\sigma_x-\nu(\sigma_y+\sigma_z)],\quad
\gamma_{xy}=\frac{\tau_{xy}}{G}.
$$
由纯剪切与单轴试验约束关系得 $G=\dfrac{E}{2(1+\nu)}$（推导：纯剪 $\tau$ 等效于 $\sigma_1=\tau,\sigma_3=-\tau$，$\varepsilon_1=(1+\nu)\tau/E$，而 $\varepsilon_1=\tau/(2G)$，故 $G=E/[2(1+\nu)]$）。

物理意义：E 表征刚度，ν 表征横向收缩，G 表征剪切刚度。

特例：不可压 $\nu=0.5\Rightarrow G=E/3$。应用：弹性变形计算。

**知识点5：体积应变与体积模量**

定义：体积应变 $\varepsilon_v=\varepsilon_x+\varepsilon_y+\varepsilon_z$。

推导：由广义胡克定律求和
$$
\varepsilon_v=\frac1E[(\sigma_x+\sigma_y+\sigma_z)-2\nu(\sigma_x+\sigma_y+\sigma_z)]=\frac{1-2\nu}{E}(\sigma_x+\sigma_y+\sigma_z).
$$
令平均应力 $\sigma_m=\frac13(\sigma_x+\sigma_y+\sigma_z)$，则 $\varepsilon_v=\dfrac{\sigma_m}{K}$，$K=\dfrac{E}{3(1-2\nu)}$。

物理意义：静水压力只产生体积变化。

特例：$\nu=0.5$ 时 $K\to\infty$ 不可压。应用：橡胶、生物组织。

**知识点6：应变能密度**

定义：单位体积弹性应变能 $w=\dfrac12\sigma_{ij}\varepsilon_{ij}$。

推导：加载过程中外力做功 $\int\sigma\,d\varepsilon$，线弹性 $\sigma=E\varepsilon$ 给出 $\int_0^{\varepsilon}\sigma d\varepsilon=\frac12\sigma\varepsilon$。三维 $w=\frac12\sigma_{ij}\varepsilon_{ij}$，代入胡克定律
$$
w=\frac{1}{2E}\left[\sum\sigma_i^2-2\nu(\sigma_x\sigma_y+\sigma_y\sigma_z+\sigma_z\sigma_x)\right]+\frac{1}{2G}(\tau_{xy}^2+\tau_{yz}^2+\tau_{zx}^2).
$$

物理意义：弹性储能，能量守恒基础。

特例：体积-畸变分解 $w=w_v+w_d$，$w_v=\frac{1-2\nu}{6E}(\sigma_x+\sigma_y+\sigma_z)^2$，$w_d=\frac{1+\nu}{6E}\sum(\sigma_i-\sigma_j)^2$。应用：强度理论。

**知识点7：应变协调方程（相容方程）**

定义：应变分量须满足的附加关系以保证单值连续位移存在。

推导：由几何方程 $\varepsilon_{ij}=\frac12(u_{i,j}+u_{j,i})$，对 $\varepsilon$ 求二阶导消去 $\mathbf u$。二维情形
$$
\frac{\partial^2\varepsilon_x}{\partial y^2}+\frac{\partial^2\varepsilon_y}{\partial x^2}=\frac{\partial^2\gamma_{xy}}{\partial x\partial y}.
$$
推导：$\varepsilon_x=u_{,x}$、$\varepsilon_y=v_{,y}$、$\gamma_{xy}=u_{,y}+v_{,x}$。对 $\varepsilon_x$ 关于 $y$ 二阶导，$\varepsilon_y$ 关于 $x$ 二阶导，求和等于 $\gamma_{xy}$ 混合二阶导。

物理意义：避免出现不连续/裂缝的应变场。

特例：平面问题六个应变分量仅一个协调方程。应用：应力函数法。

**知识点8：应力函数与二维问题**

定义：Airy 应力函数 $\phi$ 使 $\sigma_x=\dfrac{\partial^2\phi}{\partial y^2}$，$\sigma_y=\dfrac{\partial^2\phi}{\partial x^2}$，$\tau_{xy}=-\dfrac{\partial^2\phi}{\partial x\partial y}$。

推导：无体力时先由 $\sigma_{xx,y}+\tau_{xy,x}=0$ 引入 $\phi$ 使自动满足平衡。再代入相容方程（平面应力，体力忽略）得双调和方程
$$
\nabla^4\phi=\frac{\partial^4\phi}{\partial x^4}+2\frac{\partial^4\phi}{\partial x^2\partial y^2}+\frac{\partial^4\phi}{\partial y^4}=0.
$$

物理意义：把求应力问题化为求单一标量函数。

特例：多项式解用于梁弯曲。应用：弹性力学解析解。

### 第4章 杆件基本变形与强度

**知识点1：轴向拉压**

定义：轴力 $N$ 沿杆轴，正应力 $\sigma=N/A$，轴向变形 $\Delta l=Nl/(EA)$。

假设：平截面假设、均匀材料、小变形。

推导：取截面法平衡 $N=P$。应变 $\varepsilon=\Delta l/l$。胡克定律 $\sigma=E\varepsilon$ 得
$$
\frac{N}{A}=E\frac{\Delta l}{l}\Rightarrow\Delta l=\frac{Nl}{EA}.
$$
应变能 $U=\int\frac{N^2}{2EA}dx$。

物理意义：$EA$ 为轴向刚度。

特例：变截面 $\Delta l=\int\frac{N}{EA}dx$。应用：拉杆、柱。

**知识点2：剪切与挤压**

定义：剪应力 $\tau=Q/A$，剪应变 $\gamma=\tau/G$。

假设：均匀分布、线弹性。

推导：连接件受剪，破坏面剪力 $Q$，平均剪应力 $\tau=Q/A_s$。剪切变形 $\Delta=\gamma t=\dfrac{Qt}{GA_s}$。

物理意义：铆钉/螺栓连接强度。

特例：单剪、双剪；挤压 $\sigma_{bs}=P/A_{bs}$。应用：连接件校核。

**知识点3：圆轴扭转**

定义：扭矩 $T$，剪应力 $\tau=\dfrac{T\rho}{I_p}$，扭转角 $\varphi=\dfrac{Tl}{GI_p}$，$I_p=\dfrac{\pi d^4}{32}$。

假设：平截面保持平面、半径直线保持直线（圆轴）。

推导：几何 $\gamma=\rho\dfrac{d\varphi}{dx}$。物理 $\tau=G\gamma=G\rho\dfrac{d\varphi}{dx}$。静力 $T=\int\rho\tau\,dA=G\dfrac{d\varphi}{dx}\int\rho^2dA=GI_p\dfrac{d\varphi}{dx}$。
$$
\therefore\ \frac{d\varphi}{dx}=\frac{T}{GI_p},\quad \tau=\frac{T\rho}{I_p},\quad \tau_{\max}=\frac{T r}{I_p}=\frac{T}{W_p}.
$$

物理意义：$GI_p$ 抗扭刚度，$W_p=\pi d^3/16$。

特例：空心轴 $I_p=\dfrac{\pi(D^4-d^4)}{32}$。应用：传动轴。

**知识点4：梁的弯曲正应力**

定义：$\sigma=\dfrac{My}{I_z}$，中性轴过形心。

假设：平面弯曲、平截面假设、纵向纤维无挤压、线弹性。

推导：几何 $\varepsilon=\dfrac{y}{\rho}$。物理 $\sigma=E\varepsilon=\dfrac{Ey}{\rho}$。静力：轴力为零 $\int\sigma dA=0\Rightarrow\int y\,dA=0$（中性轴过形心）；弯矩 $M=\int y\sigma dA=\dfrac{E}{\rho}\int y^2dA=\dfrac{EI_z}{\rho}$。
$$
\therefore\ \frac1\rho=\frac{M}{EI_z},\quad\sigma=\frac{My}{I_z},\quad\sigma_{\max}=\frac{M}{W_z}.
$$

物理意义：$EI$ 抗弯刚度，$W_z=I_z/y_{\max}$。

特例：矩形 $I_z=bh^3/12$。应用：梁强度校核。

**知识点5：梁的弯曲切应力**

定义：$\tau=\dfrac{QS_z^*}{I_z b}$。

推导：微段 $dx$ 上两截面正应力差产生纵向力 $\int_{A^*}d\sigma\,dA=\dfrac{dM}{I_z}\int_{A^*}y\,dA=\dfrac{dM}{I_z}S_z^*$。该力由侧面剪应力平衡，$\tau b\,dx=\dfrac{dM}{I_z}S_z^*$，又 $dM/dx=Q$，故 $\tau=\dfrac{QS_z^*}{I_zb}$。

物理意义：$S_z^*$ 为面积对中性轴静矩。

特例：矩形截面 $\tau_{\max}=\dfrac32\dfrac{Q}{A}$。应用：短粗梁、腹板剪切。

**知识点6：挠曲线微分方程**

定义：挠度 $w(x)$，转角 $\theta\approx dw/dx$，曲率 $\kappa\approx d^2w/dx^2$。

推导：由 $1/\rho=M/(EI)$，且 $\rho^{-1}=\dfrac{w''}{(1+w'^2)^{3/2}}\approx w''$（小变形），得
$$
EI\,w''=M(x).
$$
积分两次并施加边界条件。剪力与分布载荷关系
$$
EI w^{(4)}=q(x).
$$

物理意义：挠度由弯矩积分得。

特例：简支梁中点集中力 $w_{\max}=\dfrac{Pl^3}{48EI}$。应用：刚度校核。

**知识点7：叠加法与变形比较**

定义：线弹性小变形下多个载荷效应可叠加。

推导：微分方程 $EIw''=\sum M_i$ 线性，解可叠加。超静定问题用变形协调附加方程。

物理意义：分解复杂载荷。

特例：悬臂梁端部集中力 $w=\dfrac{Pl^3}{3EI}$，均布载荷 $w=\dfrac{ql^4}{8EI}$。应用：静不定梁求解。

**知识点8：强度理论与组合变形**

定义：四种经典强度理论——最大拉应力（第一）、最大伸长线应变（第二）、最大剪应力（第三）、畸变能（第四）。

公式：
$$
\sigma_{r1}=\sigma_1,\quad \sigma_{r2}=\sigma_1-\nu(\sigma_2+\sigma_3),
$$
$$
\sigma_{r3}=\sigma_1-\sigma_3,\quad
\sigma_{r4}=\sqrt{\tfrac12[(\sigma_1-\sigma_2)^2+(\sigma_2-\sigma_3)^2+(\sigma_3-\sigma_1)^2]}.
$$

推导（第四理论）：畸变能密度 $w_d=\dfrac{1+\nu}{6E}[(\sigma_1-\sigma_2)^2+(\sigma_2-\sigma_3)^2+(\sigma_3-\sigma_1)^2]$。单向拉伸屈服时 $\sigma_1=\sigma_s,\sigma_2=\sigma_3=0$，$w_d^s=\dfrac{1+\nu}{3E}\sigma_s^2$。令 $w_d=w_d^s$ 得 $\sigma_{r4}$ 表达式。

物理意义：不同材料适用不同理论，脆性用第一/二，塑性用第三/四。

特例：纯剪 $\sigma_1=\tau,\sigma_3=-\tau$，$\sigma_{r3}=2\tau$。应用：弯扭组合 $\sigma_{r3}=\sqrt{\sigma^2+4\tau^2}$，$\sigma_{r4}=\sqrt{\sigma^2+3\tau^2}$。

### 第5章 能量原理与变分法

**知识点1：外力功与应变能**

定义：外力功 $W=\int\mathbf F\cdot d\boldsymbol\delta$，应变能 $U$。

推导：线弹性体上力从 0 线性增至 $F$，$W=\int_0^\delta F'd\delta'=\frac12F\delta$。对杆 $W=\frac12P\Delta l=\frac12P\cdot\dfrac{Pl}{EA}=\dfrac{P^2l}{2EA}=U$。一般
$$
U=\int_V\frac12\sigma_{ij}\varepsilon_{ij}\,dV.
$$

物理意义：弹性体储能，能量守恒。

特例：线弹性 $U=W$；非线弹性不成立。应用：能量法求位移。

**知识点2：虚功原理**

定义：变形体平衡的必要充分条件是对任意协调虚位移，外力虚功等于内力虚功
$$
\delta W_e=\delta W_i,\quad \delta W_e=\int \mathbf f\cdot\delta\mathbf u\,dV+\int\mathbf t\cdot\delta\mathbf u\,dS,\quad \delta W_i=\int_V\sigma_{ij}\delta\varepsilon_{ij}\,dV.
$$

推导：由平衡方程 $\sigma_{ij,j}+f_i=0$ 乘虚位移 $\delta u_i$ 积分并分部积分：
$$
\int_V(\sigma_{ij,j}+f_i)\delta u_i\,dV=0\Rightarrow\int_S\sigma_{ij}n_j\delta u_i\,dS-\int_V\sigma_{ij}\delta u_{i,j}\,dV+\int_Vf_i\delta u_i\,dV=0.
$$
利用 $\sigma_{ij}\delta u_{i,j}=\sigma_{ij}\delta\varepsilon_{ij}$（对称性）与面力边界得 $\delta W_e=\delta W_i$。

物理意义：虚功原理不限于线弹性，是有限元基础。

特例：刚体虚位移 $\delta W_i=0$ 退化为静力学虚位移原理。应用：推导有限元方程。

**知识点3：最小势能原理**

定义：总势能 $\Pi=U-W$，平衡对应于 $\delta\Pi=0$，稳定平衡为极小值。

推导：$U$ 为应变能，$W$ 为给定外力功。对虚位移 $\delta\Pi=\delta U-\delta W=\delta W_i-\delta W_e=0$（虚功原理）。二阶变分 $\delta^2U>0$（正定弹性张量）保证极小。

物理意义：把平衡问题化为泛函极值。

特例：线性弹簧 $\Pi=\frac12kx^2-Fx$，$\delta\Pi=0\Rightarrow kx=F$。应用：瑞利-里兹法。

**知识点4：卡氏定理**

定义：第一定理 $\delta_i=\partial U/\partial F_i$（$U$ 以位移表达则 $\partial U/\partial\delta_i=F_i$）。

推导（第二定理）：设 $U=\frac12\sum_{ij}a_{ij}F_iF_j$，则 $\partial U/\partial F_i=\frac12\sum_j a_{ij}F_j+\frac12\sum_j a_{ji}F_j=\sum_j a_{ij}F_j=\delta_i$（利用对称性 $a_{ij}=a_{ji}$）。含常数力 $F_0$ 时先加虚力后令其为零。

物理意义：由能量对力求位移，对位移求力。

特例：简支梁中点力 $P$，$\delta=\partial U/\partial P$。应用：超静定、桁架位移。

**知识点5：单位载荷法（莫尔积分）**

定义：$\Delta=\int\dfrac{M\bar M}{EI}dx+\int\dfrac{N\bar N}{EA}dx+\int\dfrac{T\bar T}{GI_p}dx$，$\bar M$ 为单位力引起的弯矩。

推导：由卡氏第二定理，真实载荷位移 $\Delta=\partial U/\partial P$。对线性系统 $\partial M/\partial P=\bar M$，
$$
U=\int\frac{M^2}{2EI}dx,\quad \Delta=\frac{\partial U}{\partial P}=\int\frac{M}{EI}\frac{\partial M}{\partial P}dx=\int\frac{M\bar M}{EI}dx.
$$

物理意义：单位力法普适求位移。

特例：桁架 $\Delta=\sum\dfrac{N_i\bar N_iL_i}{EA_i}$。应用：结构位移计算、图形互乘法。

**知识点6：瑞利-里兹法**

定义：设近似解 $\mathbf u=\sum a_i\boldsymbol\phi_i$（满足位移边界），代入势能泛函，由 $\partial\Pi/\partial a_i=0$ 求系数。

推导：$\Pi[\mathbf u]=\int_V w(\varepsilon)dV-\int\mathbf f\cdot\mathbf u dV$。代入试探函数得 $\Pi(a_1,\dots,a_n)=\frac12\mathbf a^T\mathbf K\mathbf a-\mathbf a^T\mathbf F$，极小条件 $\partial\Pi/\partial a_i=0$ 给出 $\mathbf K\mathbf a=\mathbf F$。

物理意义：能量逼近，收敛于真解（势能单调下降）。

特例：一阶近似即得较准解。应用：梁振动频率、板弯曲。

**知识点7：伽辽金法**

定义：残量加权积分为零，权函数取试探函数，$\int_V \mathcal L(\mathbf u_N)\cdot\boldsymbol\phi_i\,dV=0$。

推导：控制方程 $\mathcal L\mathbf u=\mathbf f$。近似 $\mathbf u_N$ 代入残量 $R=\mathcal L\mathbf u_N-\mathbf f\ne0$，要求 $R$ 与基函数正交，即 $\int R\phi_i=0$。分部积分后与虚功/能量法等价。

物理意义：不依赖泛函存在性（可用于非自伴问题）。

特例：与里兹法在自伴问题时等价。应用：流体、传热有限元。

**知识点8：余能与互补虚功原理**

定义：余能 $U^*=\int_V\sigma_{ij}\varepsilon_{ij}dV-U$，线弹性 $U^*=U$。

推导：$\delta U=\int\varepsilon_{ij}\delta\sigma_{ij}dV$（线弹性），互补虚功 $\delta U^*=\int\varepsilon_{ij}\delta\sigma_{ij}dV$，对应给定位移边界的平衡。总余能 $\Pi^*=U^*-W^*$ 取极小对应。

物理意义：与势能原理互补，是对偶变分原理。

特例：线弹性 $U^*=U$。应用：应力法有限元。

### 第6章 质点与质点系动力学

**知识点1：牛顿运动定律**

定义：惯性系中 $\mathbf F=m\mathbf a$。

假设：惯性参考系、质点为常质量。

推导：由动量定义 $\mathbf p=m\mathbf v$，$\mathbf F=\dfrac{d\mathbf p}{dt}=m\dfrac{d\mathbf v}{dt}=m\mathbf a$（$m$ 常数）。

物理意义：力改变动量。

特例：$m$ 可变时用 $\mathbf F=\dot{\mathbf p}$。应用：一般动力学。

**知识点2：质点运动微分方程**

定义：直角坐标 $m\ddot x=F_x$ 等。

推导（自然坐标）：切向 $m\dfrac{dv}{dt}=F_\tau$，法向 $m\dfrac{v^2}{\rho}=F_n$。

物理意义：切向控制速率，法向控制方向。

特例：圆周运动 $a=v^2/R$。应用：轨道、滑块。

**知识点3：动量定理**

定义：$\dfrac{d\mathbf p}{dt}=\mathbf F^{(e)}$。

推导：质点系 $\mathbf p=\sum m_i\mathbf v_i$，$\dot{\mathbf p}=\sum m_i\mathbf a_i=\sum\mathbf F_i=\sum\mathbf F_i^{(e)}+\sum\mathbf F_i^{(i)}$。内力成对抵消 $\sum\mathbf F_i^{(i)}=0$。

物理意义：内力不改变总动量。

特例：$\mathbf F^{(e)}=0$ 动量守恒。应用：碰撞。

**知识点4：角动量定理**

定义：$\dfrac{d\mathbf L_O}{dt}=\mathbf M_O^{(e)}$。

推导：$\mathbf L_O=\sum\mathbf r_i\times m_i\mathbf v_i$，
$$
\dot{\mathbf L}_O=\sum\mathbf r_i\times m_i\mathbf a_i=\sum\mathbf r_i\times\mathbf F_i=\mathbf M_O^{(e)}+\mathbf M_O^{(i)}.
$$
内力矩成对抵消（沿连线），故等于外力矩。

物理意义：外力矩改变角动量。

特例：有心力 $\mathbf M_O=0$ 角动量守恒。应用：开普勒、陀螺。

**知识点5：动能定理**

定义：$dT=\delta W$。

推导：质点 $\mathbf F\cdot d\mathbf r=m\dfrac{d\mathbf v}{dt}\cdot\mathbf vdt=d\left(\frac12mv^2\right)$。求和得 $dT=\sum\mathbf F_i\cdot d\mathbf r_i=\delta W^{(e)}+\delta W^{(i)}$；刚体内部约束功为零。

物理意义：功与动能变化联系。

特例：保守力 $dT=-dV\Rightarrow T+V=E$ 常量。应用：机械能守恒。

**知识点6：达朗贝尔原理**

定义：引入惯性力 $\mathbf F_I=-m\mathbf a$，动力问题化为静力平衡 $\mathbf F+\mathbf F_I=0$。

推导：牛顿方程 $\mathbf F=m\mathbf a\Rightarrow\mathbf F-m\mathbf a=0\Rightarrow\mathbf F+\mathbf F_I=0$。

物理意义：用静力学方法解动力学。

特例：静平衡 $\mathbf a=0$ 退化为静力学。应用：转子动力学、动反力。

**知识点7：变质量系统与火箭方程**

定义：变质量 $\dfrac{d(m\mathbf v)}{dt}=\mathbf F_{ext}$。

推导：设系统质量变化率 $\dot m$，喷出物相对速度 $\mathbf u$。动量守恒：
$$
m\,d\mathbf v=-(\dot m\,dt)\,\mathbf u_{rel}\Rightarrow m\dot{\mathbf v}=-\dot m\,\mathbf u_{rel}.
$$
一维真空无外力，$\mathbf u_{rel}$ 与速度反向大小 $u_e$：
$$
m\frac{dv}{dt}=-\dot m\,u_e\Rightarrow dv=-u_e\frac{dm}{m}\Rightarrow\Delta v=u_e\ln\frac{m_0}{m_f}.
$$

物理意义：齐奥尔科夫斯基公式。

特例：$m_f\to0$ 时 $\Delta v\to\infty$（忽略结构质量）。应用：火箭推进。

**知识点8：有心力与开普勒定律**

定义：力沿质心连线 $\mathbf F=F(r)\hat{\mathbf r}$。

推导（第二定律）：$\mathbf M=0\Rightarrow\mathbf L$ 守恒，面积速度 $\dfrac{dA}{dt}=\dfrac{L}{2m}=$ 常量。
推导（轨道方程）：比内公式
$$
\frac{d^2u}{d\theta^2}+u=-\frac{m}{L^2u^2}F(1/u),\quad u=1/r.
$$
对平方反比引力 $F=-k/r^2$ 得圆锥曲线 $u=\dfrac{mk}{L^2}(1+e\cos\theta)$。

物理意义：行星轨道为椭圆。

特例：圆轨道 $e=0$。应用：卫星轨道设计。

### 第7章 拉格朗日力学

**知识点1：广义坐标与约束**

定义：用 $q_1,\dots,q_n$ 描述系统位形，$n=3N-s$（$s$ 为完整约束数）。

假设：完整约束 $f(\mathbf r_i,t)=0$；理想约束。

推导：$N$ 个质点每个 3 自由度，共 $3N$，减去 $s$ 个约束得独立广义坐标数 $n$。

物理意义：自动消去约束力。

特例：单摆 $n=1$，$q=\theta$。应用：多体系统建模。

**知识点2：虚功与广义力**

定义：$\delta W=\sum\mathbf F_i\cdot\delta\mathbf r_i=\sum_j Q_j\delta q_j$，$Q_j=\sum_i\mathbf F_i\cdot\dfrac{\partial\mathbf r_i}{\partial q_j}$。

推导：$\delta\mathbf r_i=\sum_j\dfrac{\partial\mathbf r_i}{\partial q_j}\delta q_j$，代入 $\delta W$ 得 $Q_j$ 定义。

物理意义：广义力由位形映射确定。

特例：保守力 $Q_j=-\partial V/\partial q_j$。应用：求广义力。

**知识点3：拉格朗日方程的推导**

定义：$\dfrac{d}{dt}\left(\dfrac{\partial T}{\partial\dot q_j}\right)-\dfrac{\partial T}{\partial q_j}=Q_j$。

推导：由达朗贝尔-拉格朗日方程 $\sum_i(\mathbf F_i-m_i\mathbf a_i)\cdot\delta\mathbf r_i=0$。动能项
$$
\sum_i m_i\mathbf a_i\cdot\dfrac{\partial\mathbf r_i}{\partial q_j}
=\sum_i m_i\mathbf v_i\cdot\frac{d}{dt}\frac{\partial\mathbf r_i}{\partial q_j}
=\frac{d}{dt}\sum_i m_i\mathbf v_i\cdot\frac{\partial\mathbf v_i}{\partial\dot q_j}-\sum_i m_i\mathbf v_i\cdot\frac{\partial\mathbf v_i}{\partial q_j}
=\frac{d}{dt}\frac{\partial T}{\partial\dot q_j}-\frac{\partial T}{\partial q_j}.
$$
（用了 $\partial\mathbf r_i/\partial q_j=\partial\mathbf v_i/\partial\dot q_j$。）主动力项为 $Q_j$，故得方程。

物理意义：标量形式，无需内力。

特例：保守系统 $L=T-V$，$\dfrac{d}{dt}\dfrac{\partial L}{\partial\dot q_j}-\dfrac{\partial L}{\partial q_j}=0$。应用：机器人、机构。

**知识点4：循环坐标与动量守恒**

定义：$L$ 不显含 $q_k$ 时 $\dfrac{\partial L}{\partial q_k}=0$，对应广义动量 $p_k=\partial L/\partial\dot q_k$ 守恒。

推导：拉格朗日方程 $\dfrac{d}{dt}p_k=\dfrac{\partial L}{\partial q_k}=0\Rightarrow p_k=$ 常量。

物理意义：对称性对应守恒量（诺特）。

特例：中心力问题中角动量守恒。应用：守恒量化简方程。

**知识点5：诺特定理（简介）**

定义：$L$ 的每种连续对称性对应一个守恒量。

推导：设变换 $q\to q+\epsilon\xi$ 使 $\delta L=\epsilon\dfrac{dF}{dt}$。由欧拉-拉格朗日方程
$$
\delta L=\sum\left(\frac{\partial L}{\partial q}\epsilon\xi+\frac{\partial L}{\partial\dot q}\epsilon\dot\xi\right)=\epsilon\frac{d}{dt}\sum p_k\xi_k,
$$
故 $J=\sum p_k\xi_k-F$ 守恒。

物理意义：时间平移→能量，空间平移→动量，转动→角动量。

特例：自由粒子动量守恒。应用：场论、守恒律。

**知识点6：耗散函数与拉格朗日方程**

定义：瑞利耗散函数 $\mathcal F=\frac12\sum c_{ij}\dot q_i\dot q_j$。

推导：耗散力 $Q_j^{d}=-\dfrac{\partial\mathcal F}{\partial\dot q_j}$，拉格朗日方程变为
$$
\frac{d}{dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}+\frac{\partial\mathcal F}{\partial\dot q_j}=0.
$$

物理意义：黏性阻尼导致能量单调减少。

特例：$\mathcal F=0$ 保守系统。应用：阻尼振动。

**知识点7：非惯性系中的拉格朗日方程**

定义：在非惯性系 $L=T_{rel}-V-U_{inertial}$。

推导：位矢 $\mathbf r=\mathbf R+\mathbf r'$，动能为
$$
T=\frac12m|\dot{\mathbf R}+\boldsymbol\omega\times\mathbf r'+\dot{\mathbf r}'|^2.
$$
展开得惯性力势：平动 $V_{trans}=-m\mathbf a_0\cdot\mathbf r'$；离心 $V_{cf}=-\frac12m|\boldsymbol\omega\times\mathbf r'|^2$；科氏力 $\mathbf F_{cor}=-2m\boldsymbol\omega\times\mathbf v'$ 不能写成势，保留在广义力中。

物理意义：非惯性系引入惯性力。

特例：地球自转科氏力偏转。应用：气象、陀螺。

**知识点8：双摆与耦合系统应用**

定义：双摆广义坐标 $\theta_1,\theta_2$。

推导：位置 $\mathbf r_1=l_1(\sin\theta_1,\cos\theta_1)$，$\mathbf r_2=\mathbf r_1+l_2(\sin\theta_2,\cos\theta_2)$。动能
$$
T=\frac12m_1l_1^2\dot\theta_1^2+\frac12m_2\left[l_1^2\dot\theta_1^2+l_2^2\dot\theta_2^2+2l_1l_2\dot\theta_1\dot\theta_2\cos(\theta_1-\theta_2)\right],
$$
$V=-(m_1+m_2)gl_1\cos\theta_1-m_2gl_2\cos\theta_2$。代入欧拉-拉格朗日方程得耦合非线性方程。

物理意义：混沌典型系统。

特例：小角度线性化得耦合振子。应用：多体动力学。

### 第8章 哈密顿力学

**知识点1：勒让德变换与哈密顿函数**

定义：$H(q,p,t)=\sum_i p_i\dot q_i-L$，$p_i=\partial L/\partial\dot q_i$。

推导（可逆性）：设 $\det\left(\partial^2L/\partial\dot q_i\partial\dot q_j\right)\ne0$，则可解出 $\dot q_i=\dot q_i(q,p,t)$。全微分
$$
dH=\sum\dot q_idp_i-\sum\dot p_idq_i-\frac{\partial L}{\partial t}dt,
$$
故 $H$ 是 $q,p$ 的函数。

物理意义：相空间描述，$H$ 常为能量。

特例：保守、$T$ 为 $\dot q$ 二次齐次时 $H=T+V$。应用：正则方程。

**知识点2：哈密顿正则方程**

定义：$\dot q_i=\dfrac{\partial H}{\partial p_i}$，$\dot p_i=-\dfrac{\partial H}{\partial q_i}$。

推导：由 $H=\sum p_i\dot q_i-L$，
$$
\frac{\partial H}{\partial p_i}=\dot q_i,\quad
\frac{\partial H}{\partial q_i}=-\frac{\partial L}{\partial q_i}=-\frac{d}{dt}\frac{\partial L}{\partial\dot q_i}=-\dot p_i.
$$
（用了拉格朗日方程。）

物理意义：$2n$ 个一阶方程对称。

特例：一维谐振子 $H=p^2/(2m)+\frac12kx^2$，$\dot x=p/m$，$\dot p=-kx$。应用：辛积分、长期演化。

**知识点3：泊松括号**

定义：$\{f,g\}=\sum_i\left(\dfrac{\partial f}{\partial q_i}\dfrac{\partial g}{\partial p_i}-\dfrac{\partial f}{\partial p_i}\dfrac{\partial g}{\partial q_i}\right)$。

推导：$\dfrac{df}{dt}=\sum\left(\dfrac{\partial f}{\partial q}\dot q+\dfrac{\partial f}{\partial p}\dot p\right)+\dfrac{\partial f}{\partial t}=\{f,H\}+\dfrac{\partial f}{\partial t}$。

物理意义：守恒量满足 $\{f,H\}=0$。

特例：$\{q_i,p_j\}=\delta_{ij}$，$\{q_i,q_j\}=\{p_i,p_j\}=0$。应用：量子化（对应原理）。

**知识点4：作用量与哈密顿原理**

定义：作用量 $S=\int_{t_1}^{t_2}L(q,\dot q,t)\,dt$，真实路径使 $\delta S=0$。

推导：$\delta S=\int\left(\dfrac{\partial L}{\partial q}\delta q+\dfrac{\partial L}{\partial\dot q}\delta\dot q\right)dt$。分部积分 $\int\dfrac{\partial L}{\partial\dot q}\delta\dot q\,dt=\left[\dfrac{\partial L}{\partial\dot q}\delta q\right]_{t_1}^{t_2}-\int\dfrac{d}{dt}\dfrac{\partial L}{\partial\dot q}\delta q\,dt$。端点变分为零，故
$$
\delta S=\int\left(\frac{\partial L}{\partial q}-\frac{d}{dt}\frac{\partial L}{\partial\dot q}\right)\delta q\,dt=0,
$$
由 $\delta q$ 任意得欧拉-拉格朗日方程。

物理意义：最小作用量原理是力学最普适表述。

特例：几何光学费马原理。应用：场论、路径积分。

**知识点5：哈密顿-雅可比方程**

定义：寻求正则变换使新哈密顿 $K=0$。

推导：哈密顿主函数 $S(q,t)$ 满足 $\dot q=\partial H/\partial p$，令 $p=\partial S/\partial q$，则
$$
\frac{\partial S}{\partial t}+H\left(q,\frac{\partial S}{\partial q},t\right)=0.
$$
推导：$dS/dt=\sum p_i\dot q_i-L=H$（沿轨道），且 $dS=\sum p_i dq_i-Hdt$，故 $\partial S/\partial t=-H$。

物理意义：一阶偏微分方程，解后积分得运动。

特例：分离变量 $S=W-Et$ 给出定态问题。应用：解析力学、量子力学过渡。

**知识点6：正则变换**

定义：保持正则方程形式的变换 $(q,p)\to(Q,P)$。

推导：要求 $\sum p\,dq-Hdt=\sum P\,dQ-Kdt+dF$。取第一类生成函数 $F_1(q,Q,t)$，则
$$
p=\frac{\partial F_1}{\partial q},\quad P=-\frac{\partial F_1}{\partial Q},\quad K=H+\frac{\partial F_1}{\partial t}.
$$

物理意义：选取使 $H$ 简化的坐标。

特例：恒等变换、谐振子作用量-角变量。应用：可积系统。

**知识点7：相空间与刘维尔定理**

定义：相空间体积 $dV=\prod dq_i\,dp_i$。

推导：相流 $\dot{\mathbf z}=(\partial H/\partial p,-\partial H/\partial q)$。散度
$$
\nabla_{\mathbf z}\cdot\dot{\mathbf z}=\sum\left(\frac{\partial^2H}{\partial q\partial p}-\frac{\partial^2H}{\partial p\partial q}\right)=0.
$$
由连续性方程 $\partial_t\rho+\nabla\cdot(\rho\dot{\mathbf z})=0$ 与散度为零得 $d\rho/dt=0$，相体积守恒。

物理意义：哈密顿系统保相体积。

特例：耗散系统相体积收缩。应用：统计力学。

**知识点8：哈密顿力学应用**

定义：用正则方程解中心力、耦合振子。

推导（谐振子）：$H=p^2/(2m)+\frac12m\omega^2q^2$，$\ddot q=-\omega^2q$，解 $q=A\cos(\omega t+\phi)$。作用量变量 $J=\oint p\,dq=2\pi E/\omega$，$\omega=\partial H/\partial J$。

物理意义：能量与频率关系。

特例：小振动正则化。应用：分子振动、轨道。

### 第9章 连续介质运动学与守恒律

**知识点1：物质描述与空间描述**

定义：拉格朗日（物质）坐标 $\mathbf X$，$\mathbf x=\boldsymbol\chi(\mathbf X,t)$；欧拉（空间）坐标 $\mathbf x$。

推导（物质导数）：$f(\mathbf x(\mathbf X,t),t)$ 对时间求导
$$
\frac{Df}{Dt}=\frac{\partial f}{\partial t}+\frac{\partial f}{\partial x_i}\frac{\partial x_i}{\partial t}=\frac{\partial f}{\partial t}+\mathbf v\cdot\nabla f.
$$

物理意义：随流体质点的变化率。

特例：定常场 $\partial_t f=0$。应用：流体力学。

**知识点2：变形梯度与应变度量**

定义：$\mathbf F=\partial\mathbf x/\partial\mathbf X$，$J=\det\mathbf F$。

推导（格林应变）：$d\mathbf x=\mathbf F\,d\mathbf X$，长度平方差
$$
d\mathbf x\cdot d\mathbf x-d\mathbf X\cdot d\mathbf X=d\mathbf X\cdot(\mathbf F^T\mathbf F-\mathbf I)d\mathbf X=2\,d\mathbf X\cdot\mathbf E\,d\mathbf X,
$$
故 $\mathbf E=\frac12(\mathbf F^T\mathbf F-\mathbf I)$。小变形 $\mathbf F=\mathbf I+\nabla_{\mathbf X}\mathbf u$，得 $\boldsymbol\varepsilon=\frac12(\nabla\mathbf u+\nabla\mathbf u^T)$。

物理意义：$J$ 为体积比。

特例：不可压 $J=1$。应用：有限变形。

**知识点3：小应变张量与变形协调**

定义：$\varepsilon_{ij}=\frac12(u_{i,j}+u_{j,i})$。

推导：与第3章相同，由位移梯度对称化得；协调方程保证位移单值。

物理意义：变形度量。

特例：无旋场仅体积变形。应用：弹性力学。

**知识点4：雷诺输运定理**

定义：$\dfrac{D}{Dt}\int_{V(t)}\phi\,dV=\int_V\left(\dfrac{\partial\phi}{\partial t}+\nabla\cdot(\phi\mathbf v)\right)dV$。

推导：对物质体积，$\dfrac{D}{Dt}\int_V\phi dV=\int_V\left(\dfrac{D\phi}{Dt}+\phi\nabla\cdot\mathbf v\right)dV$。利用 $DJ/Dt=J\nabla\cdot\mathbf v$（Jacobi 公式）。展开 $\dfrac{D\phi}{Dt}=\partial_t\phi+\mathbf v\cdot\nabla\phi$ 再合并得 $\partial_t\phi+\nabla\cdot(\phi\mathbf v)$。

物理意义：把物质积分转化为控制体积分。

特例：$\phi=\rho$ 得连续性方程。应用：守恒律推导。

**知识点5：连续性方程**

定义：$\dfrac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf v)=0$。

推导：质量守恒 $\dfrac{D}{Dt}\int_V\rho dV=0$。由输运定理 $\int_V\left(\partial_t\rho+\nabla\cdot(\rho\mathbf v)\right)dV=0$，$V$ 任意故
$$
\partial_t\rho+\nabla\cdot(\rho\mathbf v)=0.
$$
不可压展开 $\dfrac{D\rho}{Dt}+\rho\nabla\cdot\mathbf v=0$，$\rho$ 常数时 $\nabla\cdot\mathbf v=0$。

物理意义：质量守恒局部形式。

特例：定常流 $\nabla\cdot(\rho\mathbf v)=0$。应用：管流 $A_1v_1=A_2v_2$。

**知识点6：动量守恒与运动方程**

定义：$\rho\dfrac{D\mathbf v}{Dt}=\nabla\cdot\boldsymbol\sigma+\mathbf f$。

推导：动量守恒 $\dfrac{D}{Dt}\int_V\rho\mathbf v\,dV=\int_S\boldsymbol\sigma\cdot\mathbf n\,dS+\int_V\mathbf f\,dV$。由输运定理与散度定理，左端 $=\int_V[\partial_t(\rho\mathbf v)+\nabla\cdot(\rho\mathbf v\otimes\mathbf v)]dV$，右端 $=\int_V(\nabla\cdot\boldsymbol\sigma+\mathbf f)dV$。用连续性方程化简 $\partial_t(\rho v_i)+\partial_j(\rho v_iv_j)=\rho Dv_i/Dt$，得
$$
\rho\frac{Dv_i}{Dt}=\sigma_{ij,j}+f_i.
$$

物理意义：动量守恒。

特例：弹性体静力 $\sigma_{ij,j}+f_i=0$。应用：流体、固体。

**知识点7：能量守恒**

定义：$\rho\dfrac{De}{Dt}=\boldsymbol\sigma:\nabla\mathbf v-\nabla\cdot\mathbf q+r$。

推导：总能量 $\int_V\rho(e+\frac12v^2)dV$。守恒右端为面力功率、热流、热源。结合动量方程消去动能项，得内能方程
$$
\rho\frac{De}{Dt}=\sigma_{ij}\frac{\partial v_i}{\partial x_j}-q_{i,i}+r.
$$
$\boldsymbol\sigma:\nabla\mathbf v$ 为应力功率（可分解为 $\boldsymbol\sigma:\mathbf D$，反对称部分贡献零）。

物理意义：热力学第一定律连续形式。

特例：无热流、绝热 $r=q=0$。应用：传热、可压缩流。

**知识点8：客观性与本构方程**

定义：本构关系须满足客观性（物质标架无差异性）。

推导：转动标架 $\mathbf x^*=\mathbf Q\mathbf x$，要求应力、应变按 $\boldsymbol\sigma^*=\mathbf Q\boldsymbol\sigma\mathbf Q^T$、$\mathbf F^*=\mathbf Q\mathbf F$ 变换。故本构须以客观量（如右柯西-格林张量 $\mathbf C=\mathbf F^T\mathbf F$）表达。

物理意义：排除刚体转动引起的虚假响应。

特例：各向同性弹性 $\mathbf S=2\partial W/\partial\mathbf C$。应用：超弹性、有限元。

### 第10章 流体力学与 Navier-Stokes 方程

**知识点1：流体性质与黏性**

定义：动力黏度 $\mu$，运动黏度 $\nu=\mu/\rho$。牛顿黏性定律 $\tau=\mu\dfrac{du}{dy}$。

推导：平板间距 $h$，上板速度 $U$，$du/dy=U/h$，剪切力 $F=\mu A U/h$，$\tau=F/A=\mu U/h$。

物理意义：黏性由分子动量交换产生。

特例：理想流体 $\mu=0$。应用：润滑、边界层。

**知识点2：静水压力与阿基米德原理**

定义：$p=p_0+\rho gh$，浮力 $F_b=\rho_f gV$。

推导：流体静止 $\nabla p=\rho\mathbf g$。一维 $\dfrac{dp}{dz}=-\rho g\Rightarrow p=p_0-\rho g(z-z_0)=p_0+\rho gh$（$h$ 向下深度）。浮力：浸没体表面压力积分
$$
\mathbf F_b=-\oint p\mathbf n\,dS=-\int_V\nabla p\,dV=-\int_V\rho_f\mathbf g\,dV=-\rho_fVg\hat{\mathbf z}.
$$

物理意义：压力随深度线性增加。

特例：连通器等压面。应用：浮力、水坝。

**知识点3：伯努利方程**

定义：定常、不可压、无黏沿流线 $\dfrac{p}{\rho}+\frac12v^2+gz=$ 常量。

推导：沿流线欧拉方程 $\rho v\dfrac{dv}{ds}=-\dfrac{dp}{ds}-\rho g\dfrac{dz}{ds}$。两端除以 $\rho$ 并沿 $s$ 积分
$$
\frac12v^2+\frac{p}{\rho}+gz=C.
$$

物理意义：机械能守恒。

特例：文丘里管测流量。应用：皮托管、孔板。

**知识点4：欧拉方程**

定义：理想流体 $\rho\dfrac{D\mathbf v}{Dt}=-\nabla p+\mathbf f$。

推导：由动量守恒取 $\sigma_{ij}=-p\delta_{ij}$（理想流体），代入 $\rho Dv_i/Dt=\sigma_{ij,j}+f_i$ 得
$$
\rho\frac{D\mathbf v}{Dt}=-\nabla p+\mathbf f.
$$

物理意义：无黏动量方程。

特例：静力 $D\mathbf v/Dt=0$ 退化为静水。应用：高速流动。

**知识点5：Navier-Stokes 方程推导**

定义：$\rho\dfrac{D\mathbf v}{Dt}=-\nabla p+\mu\nabla^2\mathbf v+\frac13\mu\nabla(\nabla\cdot\mathbf v)+\mathbf f$。

推导：牛顿流体本构
$$
\sigma_{ij}=-p\delta_{ij}+\lambda\delta_{ij}\frac{\partial v_k}{\partial x_k}+2\mu D_{ij},\quad D_{ij}=\frac12(v_{i,j}+v_{j,i}).
$$
代入动量方程 $\rho Dv_i/Dt=\sigma_{ij,j}+f_i$：
$$
\sigma_{ij,j}=-p_{,i}+\lambda(\nabla\cdot\mathbf v)_{,i}+2\mu D_{ij,j}.
$$
对不可压或 Stokes 假设 $\lambda=-\frac23\mu$，且 $D_{ij,j}=\frac12(v_{i,jj}+v_{j,ij})=\frac12\nabla^2v_i+\frac12(\nabla\cdot\mathbf v)_{,i}$，合并得
$$
\rho\frac{D\mathbf v}{Dt}=-\nabla p+\mu\nabla^2\mathbf v+\left(\lambda+\mu\right)\nabla(\nabla\cdot\mathbf v)+\mathbf f,
$$
常写 $\lambda+\mu=\frac13\mu$ 或不可压略去。

物理意义：黏性动量方程，流体力学核心。

特例：$Re\ll1$ 斯托克斯流 $\mu\nabla^2\mathbf v=\nabla p$。应用：管流、边界层。

**知识点6：涡量与涡动力学**

定义：$\boldsymbol\omega=\nabla\times\mathbf v$。

推导（亥姆霍兹方程）：对 N-S 方程取旋度。用 $\mathbf v\cdot\nabla\mathbf v=\nabla(v^2/2)-\mathbf v\times\boldsymbol\omega$，取旋度
$$
\frac{\partial\boldsymbol\omega}{\partial t}+(\mathbf v\cdot\nabla)\boldsymbol\omega=(\boldsymbol\omega\cdot\nabla)\mathbf v+\nu\nabla^2\boldsymbol\omega.
$$
右端第一项为涡拉伸，第二项为黏性扩散。

物理意义：无黏正压流体涡管随流保持（开尔文定理）。

特例：无黏、无涡拉伸时 $\boldsymbol\omega$ 守恒。应用：湍流、涡环。

**知识点7：层流管流（哈根-泊肃叶）**

定义：圆管层流 $u(r)=\dfrac{\Delta p}{4\mu L}(R^2-r^2)$，流量 $Q=\dfrac{\pi R^4\Delta p}{8\mu L}$。

推导：柱坐标 N-S，定常轴对称，$v_r=v_\theta=0$，$v_z=u(r)$，
$$
0=-\frac{dp}{dz}+\mu\frac1r\frac{d}{dr}\left(r\frac{du}{dr}\right).
$$
积分得 $u(r)=\dfrac{1}{4\mu}\dfrac{dp}{dz}r^2+C_1\ln r+C_2$。$r=0$ 有界 $C_1=0$，$u(R)=0$ 定 $C_2$，令 $\dfrac{dp}{dz}=-\dfrac{\Delta p}{L}$ 得表达式。流量 $Q=\int_0^R2\pi r u\,dr=\dfrac{\pi R^4\Delta p}{8\mu L}$。

物理意义：$Q\propto R^4$。

特例：平均速度 $u_{avg}=Q/(\pi R^2)=\dfrac{R^2\Delta p}{8\mu L}=\frac12u_{\max}$。应用：管道、血管。

**知识点8：边界层与量级分析**

定义：边界层厚度 $\delta(x)$，位移厚度 $\delta^*=\int_0^\infty(1-u/U)dy$，动量厚度 $\theta=\int_0^\infty\dfrac{u}{U}(1-\dfrac{u}{U})dy$。

推导（布拉修斯数量级）：沿壁 $u\sim U$，$\partial u/\partial x\sim U/L$，连续性 $\partial v/\partial y\sim U/L\Rightarrow v\sim U\delta/L$。动量方程中惯性项 $\rho u\partial u/\partial x\sim\rho U^2/L$，黏性项 $\mu\partial^2u/\partial y^2\sim\mu U/\delta^2$。平衡 $\rho U^2/L\sim\mu U/\delta^2$ 得
$$
\delta\sim\sqrt{\frac{\nu L}{U}}=\frac{L}{\sqrt{Re_L}}.
$$

物理意义：高 $Re$ 时边界层薄。

特例：层流平板 $\delta\approx5.0x/\sqrt{Re_x}$。应用：摩擦阻力、分离。

**知识点9：流动相似与阻力**

定义：阻力系数 $C_D=\dfrac{F_D}{\frac12\rho U^2A}$。

推导：由量纲分析 $C_D=f(Re,\text{shape})$。层流平板摩擦阻力 $C_f=1.328/\sqrt{Re_L}$（布拉修斯解积分）。

物理意义：无量纲阻力只依赖形状与 $Re$。

特例：湍流 $C_f=0.074/Re_L^{1/5}$。应用：飞行器阻力估算。

### 第11章 振动理论与模态分析

**知识点1：单自由度自由振动**

定义：$m\ddot x+kx=0$，固有频率 $\omega_n=\sqrt{k/m}$。

推导：牛顿定律 $m\ddot x=-kx$。设 $x=e^{rt}$，$mr^2+k=0$，$r=\pm i\omega_n$，故 $x=A\cos\omega_nt+B\sin\omega_nt=C\cos(\omega_nt-\phi)$。

物理意义：简谐振动，振幅由初条件定。

特例：弹簧-质量 $T=2\pi\sqrt{m/k}$。应用：测质量、隔振。

**知识点2：阻尼自由振动**

定义：$m\ddot x+c\dot x+kx=0$，阻尼比 $\zeta=\dfrac{c}{2\sqrt{km}}$。

推导：令 $\omega_n=\sqrt{k/m}$，$2\zeta\omega_n=c/m$，方程 $\ddot x+2\zeta\omega_n\dot x+\omega_n^2x=0$。特征根
$$
r_{1,2}=-\zeta\omega_n\pm\omega_n\sqrt{\zeta^2-1}.
$$
$\zeta<1$ 时 $\omega_d=\omega_n\sqrt{1-\zeta^2}$，$x=e^{-\zeta\omega_nt}(A\cos\omega_dt+B\sin\omega_dt)$。

物理意义：$\zeta$ 控制衰减。

特例：$\zeta=1$ 临界阻尼最快无振荡。应用：减震器、仪表。

**知识点3：受迫振动与共振**

定义：$m\ddot x+c\dot x+kx=F_0\cos\omega t$。

推导：稳态解 $x=X\cos(\omega t-\phi)$，代入得
$$
X=\frac{F_0/k}{\sqrt{(1-r^2)^2+(2\zeta r)^2}},\quad \tan\phi=\frac{2\zeta r}{1-r^2},\quad r=\frac{\omega}{\omega_n}.
$$
推导（放大因子）：设 $x=Xe^{i\omega t}$，$(k-m\omega^2+i c\omega)X=F_0$，取模得 $X$。

物理意义：振幅放大因子 $M=Xk/F_0$。

特例：$\zeta=0$、$r=1$ 共振发散；有阻尼共振峰 $r\approx\sqrt{1-2\zeta^2}$。应用：发电机、地震响应。

**知识点4：多自由度与模态叠加**

定义：$\mathbf M\ddot{\mathbf x}+\mathbf K\mathbf x=\mathbf F$。

推导：设 $\mathbf x=\boldsymbol\phi e^{i\omega t}$，得广义特征值问题 $(\mathbf K-\omega^2\mathbf M)\boldsymbol\phi=0$。$\det(\mathbf K-\omega^2\mathbf M)=0$ 给固有频率。模态正交性：$\boldsymbol\phi_i^T\mathbf M\boldsymbol\phi_j=0$（$i\ne j$），$\boldsymbol\phi_i^T\mathbf K\boldsymbol\phi_j=0$。用模态坐标 $\mathbf x=\sum q_i\boldsymbol\phi_i$ 解耦。

物理意义：系统按模态振动。

特例：两自由度拍频。应用：结构模态测试。

**知识点5：拉格朗日方法建立振动方程**

定义：$L=T-V$，$\dfrac{d}{dt}\dfrac{\partial L}{\partial\dot q}-\dfrac{\partial L}{\partial q}=0$。

推导：小振动 $T=\frac12\dot{\mathbf q}^T\mathbf M\dot{\mathbf q}$，$V=V_0+\frac12\mathbf q^T\mathbf K\mathbf q$（平衡点线性项为零），代入得 $\mathbf M\ddot{\mathbf q}+\mathbf K\mathbf q=0$。

物理意义：自动处理多自由度耦合。

特例：双摆小振动。应用：复杂结构建模。

**知识点6：连续体振动——弦与梁**

定义：弦 $T\dfrac{\partial^2w}{\partial x^2}=\rho A\dfrac{\partial^2w}{\partial t^2}$；梁 $EI\dfrac{\partial^4w}{\partial x^4}+\rho A\dfrac{\partial^2w}{\partial t^2}=0$。

推导（弦）：取微段 $dx$，张力 $T$。横向力
$$
T\sin\theta_{x+dx}-T\sin\theta_x\approx T\left(\frac{\partial w}{\partial x}\Big|_{x+dx}-\frac{\partial w}{\partial x}\Big|_x\right)=T\frac{\partial^2w}{\partial x^2}dx.
$$
牛顿定律 $\rho A\,dx\,\partial_t^2w=T\partial_x^2w\,dx$，得波动方程。梁由弯矩 $M=EIw''$ 与剪力关系 $V=\partial M/\partial x$、$\partial V/\partial x=-\rho A\partial_t^2w$ 得四阶方程。

物理意义：无穷多模态。

特例：弦模态 $\omega_n=\dfrac{n\pi}{L}\sqrt{T/(\rho A)}$；简支梁 $\omega_n=(n\pi/L)^2\sqrt{EI/(\rho A)}$。应用：乐器、桥梁。

**知识点7：模态正交性与展开**

定义：连续体模态满足 $\int_0^L\rho A\phi_m\phi_n dx=0$（$m\ne n$）。

推导：由自伴算子性质，两个不同特征值解作内积并分部积分得正交。

物理意义：任意响应可展开为模态级数。

特例：初始条件定系数。应用：动力响应。

**知识点8：隔振与传递率**

定义：力传递率 $T_r=\dfrac{F_t}{F_0}=\sqrt{\dfrac{1+(2\zeta r)^2}{(1-r^2)^2+(2\zeta r)^2}}$。

推导：由受迫振动支承力 $F_t=kx+c\dot x$，代入稳态解取模得比值。

物理意义：$r>\sqrt2$ 时 $T_r<1$ 隔振有效。

特例：高频 $T_r\to0$。应用：设备隔振、基础设计。

### 第12章 稳定性理论

**知识点1：平衡稳定性定义**

定义：受扰后偏离有限小量、扰动移除后能回到原平衡，则稳定。

假设：系统自治、平衡点存在。

推导（能量）：沿扰动轨道若 $V$ 有下界且 $\dot V\le0$，则 $V$ 单调不增，轨道有界。

物理意义：稳定性是局部性质。

特例：鞍点不稳定。应用：结构、控制。

**知识点2：势能判据**

定义：保守系统平衡点势能极小则稳定，极大或鞍点则不稳定。

推导：在平衡点 $\nabla V=0$。约束下扰动 $\delta\mathbf r$，势能变化 $\Delta V=\frac12\delta\mathbf r^T\mathbf H\delta\mathbf r+O(\delta^3)$，$\mathbf H=\nabla^2V$。若 $\mathbf H$ 正定则 $\Delta V>0$，系统在势阱内往复运动，稳定；若有负特征值则不稳定。

物理意义：势能极小对应稳定平衡。

特例：倒立摆势能极大不稳定。应用：结构屈曲、分子构型。

**知识点3：李雅普诺夫第二方法**

定义：存在正定函数 $V(\mathbf x)$，若 $\dot V\le0$ 则原点稳定；若 $\dot V<0$ 则渐近稳定。

推导：$\dot V=\nabla V\cdot\mathbf f(\mathbf x)\le0$ 保证 $V(t)\le V(0)$，给定 $\epsilon$ 取 $V$ 的等值球内含于 $\epsilon$ 球，即得稳定性。严格负定时 $V\to0\Rightarrow\|\mathbf x\|\to0$。

物理意义：构造“广义能量”判定。

特例：线性系统 $\dot{\mathbf x}=\mathbf A\mathbf x$，$V=\mathbf x^T\mathbf P\mathbf x$，李雅普诺夫方程 $\mathbf A^T\mathbf P+\mathbf P\mathbf A=-\mathbf Q$。应用：控制系统。

**知识点4：线性稳定性与特征值**

定义：$\dot{\mathbf x}=\mathbf A\mathbf x$，若 $\operatorname{Re}\lambda_i<0$ 全部则渐近稳定。

推导：解 $\mathbf x=\sum c_i\mathbf v_ie^{\lambda_it}$。$\operatorname{Re}\lambda_i<0$ 时 $e^{\lambda_it}\to0$；存在 $\operatorname{Re}\lambda_i>0$ 时发散。

物理意义：特征值实部决定稳定性。

特例：纯虚根为临界（中心）。应用：线性化分析。

**知识点5：压杆屈曲（欧拉公式）**

定义：细长压杆临界力 $P_{cr}=\dfrac{\pi^2EI}{(\mu L)^2}$。

推导：压杆微弯 $EIw''=M=-Pw$（两端铰支），$\therefore EIw''+Pw=0$。令 $k^2=P/EI$，$w''+k^2w=0$，解 $w=A\sin kx+B\cos kx$。边界 $w(0)=0\Rightarrow B=0$；$w(L)=0\Rightarrow A\sin kL=0$。非平凡解 $\sin kL=0\Rightarrow kL=n\pi$，$P=\dfrac{n^2\pi^2EI}{L^2}$，取最小 $n=1$：
$$
P_{cr}=\frac{\pi^2EI}{L^2}.
$$
一般支承用有效长度 $\mu L$。

物理意义：分岔点失稳。

特例：一端固定一端自由 $\mu=2$。应用：柱设计。

**知识点6：板屈曲**

定义：薄板临界应力 $\sigma_{cr}=k\dfrac{\pi^2D}{b^2t}$，$D=\dfrac{Et^3}{12(1-\nu^2)}$。

推导：板弯曲方程 $D\nabla^4w=N_xw_{,xx}+2N_{xy}w_{,xy}+N_yw_{,yy}$。取简支矩形板 $w=\sum A_{mn}\sin\frac{m\pi x}{a}\sin\frac{n\pi y}{b}$，代入得特征关系，最小 $N_x$ 即临界载荷。

物理意义：$D$ 为板抗弯刚度。

特例：长板 $k\approx4$。应用：飞机蒙皮。

**知识点7：分岔与后屈曲**

定义：参数越过临界值时平衡解分支出现。

推导（简单分岔）：设 $f(x,\lambda)=0$，$\lambda_c$ 处 $f_x=0$。Taylor 展开 $f\approx f_x\delta x+\frac12f_{xx}\delta x^2+f_\lambda\delta\lambda=0$，解得 $\delta x\approx\pm\sqrt{-2f_\lambda\delta\lambda/f_{xx}}$（叉式分岔）。

物理意义：临界点附近解路径。

特例：跨临界、鞍结分岔。应用：结构非线性。

**知识点8：李雅普诺夫间接法**

定义：由线性化矩阵判定非线性系统稳定。

推导：$\dot{\mathbf x}=\mathbf A\mathbf x+\mathbf g(\mathbf x)$，$\|\mathbf g\|=o(\|\mathbf x\|)$。若 $\operatorname{Re}\lambda_i(\mathbf A)<0$，则存在 $\mathbf P>0$ 使 $V=\mathbf x^T\mathbf P\mathbf x$ 满足 $\dot V\le-\alpha V$，原点渐近稳定。若存在 $\operatorname{Re}\lambda_i>0$，不稳定。

物理意义：非线性稳定性的线性近似。

特例：临界情形需高阶分析。应用：控制、生态。

### 第13章 传热与热应力基础

**知识点1：热传导微分方程**

定义：$\rho c\dfrac{\partial T}{\partial t}=\nabla\cdot(k\nabla T)+q$。

推导：取微元体，能量守恒：流入热流 $-q_{i,i}dV$，内热源 $q\,dV$，内能变化率 $\rho c\,\partial_tT\,dV$。Fourier 定律 $\mathbf q=-k\nabla T$，故
$$
\rho c\frac{\partial T}{\partial t}=-\nabla\cdot\mathbf q+q=\nabla\cdot(k\nabla T)+q.
$$
各向同性常数 $k$：$\dfrac{\partial T}{\partial t}=\alpha\nabla^2T+\dfrac{q}{\rho c}$，$\alpha=k/(\rho c)$。

物理意义：能量守恒+Fourier 定律。

特例：稳态 $\nabla^2T=0$（无源）。应用：散热设计。

**知识点2：稳态一维导热**

定义：平板 $q=-k\dfrac{dT}{dx}=\dfrac{k}{L}(T_1-T_2)$。

推导：稳态无源 $\dfrac{d}{dx}\left(k\dfrac{dT}{dx}\right)=0$，积分 $T=C_1x+C_2$，边界定常，$\dfrac{dT}{dx}=\dfrac{T_2-T_1}{L}$，热流密度 $q=-k\dfrac{T_2-T_1}{L}$。

物理意义：线性温度分布。

特例：多层壁热阻串联 $q=\dfrac{\Delta T}{\sum L_i/k_i}$。应用：保温层。

**知识点3：对流换热（牛顿冷却）**

定义：$q=h(T_s-T_\infty)$。

推导：边界层内流体导热与对流耦合。近壁 $q=-k_f\dfrac{\partial T}{\partial y}\Big|_{y=0}$，定义 $h=\dfrac{-k_f\,\partial_yT|_0}{T_s-T_\infty}$。无量纲化得 $Nu=\dfrac{hL}{k}=f(Re,Pr)$。

物理意义：$h$ 表征对流强度。

特例：层流平板 $Nu=0.332Re^{1/2}Pr^{1/3}$。应用：换热器。

**知识点4：热辐射**

定义：斯特藩-玻尔兹曼 $E=\varepsilon\sigma T^4$，净辐射 $q=\varepsilon\sigma(T_1^4-T_2^4)$。

推导：黑体谱辐射 Planck 定律积分得 $\sigma=2\pi^5k_B^4/(15h^3c^2)$。实际表面乘发射率 $\varepsilon$。

物理意义：辐射与 $T^4$ 成正比。

特例：两个大平行板净交换。应用：高温传热。

**知识点5：集总参数法**

定义：$\dfrac{T-T_\infty}{T_0-T_\infty}=e^{-t/\tau}$，$\tau=\dfrac{\rho cV}{hA_s}$。

推导：能量守恒 $\rho cV\dfrac{dT}{dt}=-hA_s(T-T_\infty)$，分离变量积分得指数律。适用判据 $Bi=\dfrac{hL_c}{k}<0.1$。

物理意义：温度空间均匀。

特例：小金属球冷却。应用：热电偶响应。

**知识点6：热应力**

定义：约束热膨胀产生 $\sigma=E\alpha\Delta T$。

推导：自由热应变 $\varepsilon_{th}=\alpha\Delta T$。若全约束，$总应变$ $\varepsilon=\varepsilon_{mech}+\varepsilon_{th}=0$，故 $\varepsilon_{mech}=-\alpha\Delta T$，$\sigma=E\varepsilon_{mech}=-E\alpha\Delta T$（压应力）。

物理意义：温度变化受约束产生内应力。

特例：两端固定杆升温受压。应用：轨道、管道。

**知识点7：热弹性本构**

定义：$\varepsilon_{ij}=\dfrac{1+\nu}{E}\sigma_{ij}-\dfrac{\nu}{E}\sigma_{kk}\delta_{ij}+\alpha\Delta T\delta_{ij}$。

推导：机械应变由胡克定律给出，总应变加上自由热应变。以柔度形式写并求逆得刚度形式
$$
\sigma_{ij}=\frac{E}{1+\nu}\left(\varepsilon_{ij}+\frac{\nu}{1-2\nu}\varepsilon_{kk}\delta_{ij}\right)-\frac{E\alpha\Delta T}{1-2\nu}\delta_{ij}.
$$

物理意义：热-力耦合本构。

特例：无约束平皿热应力为零。应用：热冲击、焊接残余应力。

**知识点8：热冲击与热扩散率**

定义：热扩散率 $\alpha=k/(\rho c)$，热冲击参数。

推导（半无限体）：初始 $T_0$，表面突降 $T_s$，解
$$
T(x,t)=T_s+(T_0-T_s)\operatorname{erf}\left(\frac{x}{2\sqrt{\alpha t}}\right),
$$
其中 $\operatorname{erf}(z)=\dfrac{2}{\sqrt\pi}\int_0^ze^{-s^2}ds$。由 Fourier 方程相似解得。

物理意义：热渗透深度 $\delta\sim\sqrt{\alpha t}$。

特例：$t\to\infty$ 整体温度趋于表面。应用：淬火、防热。

### 第14章 量测、误差与实验数据处理

**知识点1：误差分类**

定义：系统误差（可修正）、随机误差（统计）、粗大误差（剔除）。

假设：随机误差零均值、独立同分布。

推导：测量值 $x_i=x_{true}+\epsilon_i$。$\mathbb E[\epsilon]=0$ 保证无偏。

物理意义：不同误差处理方式不同。

特例：粗大误差用 $3\sigma$ 准则剔除。应用：实验数据清洗。

**知识点2：随机误差与统计**

定义：样本均值 $\bar x=\frac1n\sum x_i$，样本标准差 $s=\sqrt{\dfrac{1}{n-1}\sum(x_i-\bar x)^2}$。

推导：设 $x_i\sim N(\mu,\sigma^2)$ 独立。$\bar x\sim N(\mu,\sigma^2/n)$，故标准误 $s_{\bar x}=s/\sqrt n$。$t$ 分布置信区间 $\mu=\bar x\pm t_{\alpha/2,n-1}s/\sqrt n$。

物理意义：多次测量降低随机误差。

特例：$n\to\infty$ 标准误趋于零。应用：报告测量结果。

**知识点3：不确定度传播**

定义：$y=f(x_1,\dots,x_n)$，$u_y^2=\sum_i\left(\dfrac{\partial f}{\partial x_i}\right)^2u_{x_i}^2$（独立）。

推导：一阶 Taylor 展开 $\delta y=\sum\dfrac{\partial f}{\partial x_i}\delta x_i$，方差 $\operatorname{Var}(\delta y)=\sum\left(\dfrac{\partial f}{\partial x_i}\right)^2\operatorname{Var}(\delta x_i)+2\sum_{i<j}\dfrac{\partial f}{\partial x_i}\dfrac{\partial f}{\partial x_j}\operatorname{Cov}$，独立时协方差为零。

物理意义：合成不确定度。

特例：$y=x_1/x_2$ 时 $\dfrac{u_y}{y}=\sqrt{(u_1/x_1)^2+(u_2/x_2)^2}$。应用：实验误差估计。

**知识点4：最小二乘拟合**

定义：最小化 $S=\sum_i(y_i-f(x_i;\boldsymbol\beta))^2$。

推导：线性模型 $\mathbf y=\mathbf X\boldsymbol\beta+\boldsymbol\epsilon$。$\nabla_{\boldsymbol\beta}S=-2\mathbf X^T(\mathbf y-\mathbf X\boldsymbol\beta)=0$，得正规方程
$$
\mathbf X^T\mathbf X\boldsymbol\beta=\mathbf X^T\mathbf y\Rightarrow\hat{\boldsymbol\beta}=(\mathbf X^T\mathbf X)^{-1}\mathbf X^T\mathbf y.
$$
协方差 $\operatorname{Cov}(\hat{\boldsymbol\beta})=\sigma^2(\mathbf X^T\mathbf X)^{-1}$。

物理意义：最优线性无偏估计（Gauss-Markov）。

特例：过原点直线斜率 $\hat\beta=\sum x_iy_i/\sum x_i^2$。应用：标定曲线。

**知识点5：线性回归的显著性**

定义：决定系数 $R^2=1-\dfrac{SS_{res}}{SS_{tot}}$。

推导：$SS_{tot}=\sum(y_i-\bar y)^2$，$SS_{reg}=\sum(\hat y_i-\bar y)^2$，$SS_{res}=\sum(y_i-\hat y_i)^2$，且 $SS_{tot}=SS_{reg}+SS_{res}$。$R^2=SS_{reg}/SS_{tot}$。

物理意义：拟合优度。

特例：完全线性 $R^2=1$。应用：判断模型有效性。

**知识点6：传感器原理与标定**

定义：灵敏度 $S=dV/dQ$，分辨率、线性度、迟滞。

推导：标定用已知输入 $Q_i$ 记录输出 $V_i$，最小二乘拟合 $V=SQ+b$。

物理意义：传感器静态特性。

特例：应变片 $S_g=\dfrac{\Delta R/R}{\varepsilon}$。应用：力、温度测量。

**知识点7：信噪比与噪声**

定义：$SNR=10\log_{10}(P_s/P_n)$ dB。

推导：若信号均方 $P_s$、噪声方差 $\sigma_n^2$，$SNR=P_s/\sigma_n^2$。多次平均 $N$ 次使噪声降 $\sqrt N$，$SNR$ 提升 $10\log_{10}N$ dB。

物理意义：信号可检测性。

特例：相干平均。应用：微弱信号测量。

**知识点8：采样与采样定理**

定义：奈奎斯特频率 $f_N=f_s/2$。

推导：连续信号采样为 $x[n]=x(nT_s)$。频域复制 $\hat X(f)=\frac1{T_s}\sum_kX(f-kf_s)$。若 $X(f)$ 带限于 $B$，无混叠条件 $B<f_s/2$。重构 $x(t)=\sum x[n]\operatorname{sinc}(t/T_s-n)$。

物理意义：采样率须大于两倍最高频率。

特例：欠采样混叠。应用：数据采集系统。

### 第15章 空气动力学基础

**知识点1：势流与叠加**

定义：不可压无黏无旋，速度势 $\mathbf v=\nabla\phi$，$\nabla^2\phi=0$。

推导：无旋 $\nabla\times\mathbf v=0\Rightarrow\mathbf v=\nabla\phi$。不可压 $\nabla\cdot\mathbf v=0\Rightarrow\nabla^2\phi=0$。伯努利 $\dfrac{p}{\rho}+\frac12|\nabla\phi|^2=$ 常量。

物理意义：线性方程可叠加。

特例：均匀流+偶极子=绕圆柱流。应用：翼型初步。

**知识点2：库塔-儒可夫斯基定理**

定义：升力 $L=\rho U_\infty\Gamma$。

推导：绕任意封闭曲线的升力 $L'=\rho U_\infty\Gamma$。由伯努利 $p_\infty-p=\frac12\rho(U^2-U_\infty^2)$，升力为压力积分。以圆柱绕流验证：速度 $U_\theta=U_\infty(1+R^2/r^2)+\Gamma/(2\pi r)$ 等。

物理意义：环量产生升力。

特例：$\Gamma=0$ 无升力。应用：机翼设计、马格努斯效应。

**知识点3：薄翼理论**

定义：厚度/弯度小的翼型，扰动速度势满足线性化方程。

推导：小扰动 $\mathbf v=(U_\infty+u',v')$，$u',v'\ll U_\infty$。势方程 $(1-Ma^2)\phi_{xx}+\phi_{yy}=0$。不可压 $Ma=0$ 得 $\nabla^2\phi=0$。边界条件 $v'(x,0)=U_\infty\left(\dfrac{dy_c}{dx}\right)$（弯度）。升力系数
$$
C_L=2\pi(\alpha-\alpha_{L0}),
$$
推导：环量由后缘库塔条件 $\gamma(0)=0$ 确定，$\gamma(\theta)=2U_\infty\left(A_0\dfrac{1+\cos\theta}{\sin\theta}+\sum A_n\sin n\theta\right)$，$C_L=\pi(2A_0+A_1)$ 展开得。

物理意义：升力随攻角线性增加。

特例：对称翼 $\alpha_{L0}=0$。应用：升力线斜率。

**知识点4：边界层分离**

定义：逆压梯度 $\partial p/\partial x>0$ 使近壁动能不足，$\partial u/\partial y|_{y=0}=0$ 处分离。

推导：壁面处动量方程 $0=-\dfrac{1}{\rho}\dfrac{\partial p}{\partial x}+\nu\dfrac{\partial^2u}{\partial y^2}\Big|_0$，故 $\dfrac{\partial^2u}{\partial y^2}\Big|_0=\dfrac{1}{\mu}\dfrac{\partial p}{\partial x}$。逆压梯度时壁面速度剖面凹，直至 $\partial u/\partial y|_0=0$（分离点）。

物理意义：分离导致压差阻力与失速。

特例：顺压梯度不易分离。应用：翼型失速、扩压器。

**知识点5：可压缩流与马赫数**

定义：声速 $a=\sqrt{\gamma RT}=\sqrt{\gamma p/\rho}$，$Ma=v/a$。

推导：等熵扰动 $d\rho/\rho=dp/(\gamma p)$，连续性 $d\rho/\rho=-dv/v$，动量 $dp=-\rho v\,dv$。合并 $dp/d\rho = \gamma p/\rho$ 与 $a^2=dp/d\rho$，得 $a=\sqrt{\gamma p/\rho}=\sqrt{\gamma RT}$。

物理意义：$Ma<1$ 亚声速，$Ma>1$ 超声速。

特例：喷管临界 $Ma=1$ 时 $A=A^*$。应用：高速飞行。

**知识点6：激波关系（正激波）**

定义：兰金-许贡纽关系。

推导：控制体积守恒：连续 $\rho_1v_1=\rho_2v_2$，动量 $p_1+\rho_1v_1^2=p_2+\rho_2v_2^2$，能量 $h_1+\frac12v_1^2=h_2+\frac12v_2^2$。对理想气体可推出
$$
Ma_2^2=\frac{1+\frac{\gamma-1}{2}Ma_1^2}{\gamma Ma_1^2-\frac{\gamma-1}{2}},\quad
\frac{p_2}{p_1}=\frac{2\gamma Ma_1^2-(\gamma-1)}{\gamma+1}.
$$

物理意义：超声速气流经激波减速、增压、升温、熵增。

特例：$Ma_1\to1$ 时激波退化为马赫波。应用：进气道、激波阻力。

**知识点7：升力线理论**

定义：有限翼展，环量沿展向 $\Gamma(y)$。

推导：下洗速度 $w(y)=\dfrac{1}{4\pi}\int_{-b/2}^{b/2}\dfrac{d\Gamma/dy'}{y-y'}dy'$。有效攻角 $\alpha_{eff}=\alpha-\alpha_i$，$\alpha_i=w/U_\infty$。升力 $L'(y)=\rho U_\infty\Gamma(y)$。椭圆分布 $\Gamma=\Gamma_0\sqrt{1-(2y/b)^2}$，诱导阻力
$$
D_i=\frac{L^2}{\frac12\rho U_\infty^2\pi b^2},\quad C_{D_i}=\frac{C_L^2}{\pi\,AR}.
$$

物理意义：翼尖涡产生诱导阻力。

特例：椭圆环量分布诱导阻力最小。应用：机翼平面形状设计。

**知识点8：气动阻力组成**

定义：$C_D=C_{D0}+C_{D_i}$。

推导：$C_{D0}$ 由边界层摩擦与压差贡献，$C_{D_i}=C_L^2/(\pi AR e)$。总升阻比 $\dfrac{L}{D}=\dfrac{C_L}{C_{D0}+C_L^2/(\pi AR e)}$，最大时 $\dfrac{d}{dC_L}(L/D)=0$ 得 $C_L=\sqrt{\pi AR e\,C_{D0}}$。

物理意义：存在最佳升力系数。

特例：$C_{D0}=0$ 时升阻比随 $C_L$ 增大。应用：飞机性能优化。

### 第16章 有限元与计算力学基础

**知识点1：微分方程弱形式**

定义：将强形式乘以试函数并分部积分。

推导：一维问题 $-\dfrac{d}{dx}\left(a\dfrac{du}{dx}\right)=f$，两端乘 $v$ 并在 $[0,L]$ 积分
$$
-\int_0^L v\frac{d}{dx}\left(a\frac{du}{dx}\right)dx=\int_0^L fv\,dx.
$$
分部积分
$$
\int_0^L a\frac{du}{dx}\frac{dv}{dx}dx-\left[av\frac{du}{dx}\right]_0^L=\int_0^L fv\,dx.
$$
取 $v$ 满足齐次边界，得弱形式
$$
\int_0^L a u'v'\,dx=\int_0^L fv\,dx.
$$

物理意义：降低对解的光滑性要求。

特例：$a$ 常数、$f$ 常数解析验证。应用：有限元。

**知识点2：伽辽金有限元与单元刚度矩阵**

定义：$u\approx\sum N_j u_j$，$v=N_i$。

推导：代入弱形式 $\sum_j u_j\int_0^L a N_i'N_j'\,dx=\int_0^L fN_i\,dx$，写成 $\mathbf K\mathbf u=\mathbf F$，$K_{ij}=\int aN_i'N_j'dx$。对两节点线性单元 $N_1=1-\xi$，$N_2=\xi$（$\xi=x/h$），$N_1'=-1/h$，$N_2'=1/h$，
$$
\mathbf K^e=\frac{a}{h}\begin{pmatrix}1&-1\\-1&1\end{pmatrix},\quad
\mathbf F^e=\frac{fh}{2}\begin{pmatrix}1\\1\end{pmatrix}.
$$

物理意义：局部节点平衡。

特例：常应变杆单元。应用：结构分析。

**知识点3：杆单元与梁单元刚度**

定义（杆）：$K^e=\dfrac{EA}{L}\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$。推导同前，$a=EA$。

定义（梁单元）：节点自由度 $(w_1,\theta_1,w_2,\theta_2)$，Hermite 插值
$$
N_1=1-3\xi^2+2\xi^3,\ N_2=L(\xi-2\xi^2+\xi^3),\ N_3=3\xi^2-2\xi^3,\ N_4=L(-\xi^2+\xi^3).
$$
由 $K_{ij}=\int EI N_i''N_j''dx$ 得
$$
\mathbf K^e=\frac{EI}{L^3}\begin{pmatrix}12&6L&-12&6L\\6L&4L^2&-6L&2L^2\\-12&-6L&12&-6L\\6L&2L^2&-6L&4L^2\end{pmatrix}.
$$

物理意义：梁单元含转动自由度。

特例：固支梁静力约化。应用：框架结构。

**知识点4：单元组装与边界条件**

定义：全局矩阵由单元矩阵按自由度映射叠加。

推导：总势能 $\Pi=\sum_e\Pi_e$，对全局自由度求极值，组装规则 $K_{IJ}=\sum_e K^{(e)}_{ij}$（映射 $i\to I,j\to J$）。施加本质边界条件（置大数或消行消列）。

物理意义：整体刚度方程为 $\mathbf K\mathbf U=\mathbf F$。

特例：节点平衡。应用：大型结构。

**知识点5：等参元与数值积分**

定义：几何与场用同一形函数插值，$x=\sum N_i(\xi)x_i$，$u=\sum N_i(\xi)u_i$。

推导：单元积分换元 $dx=J\,d\xi$，$J=dx/d\xi=\sum N_i'(\xi)x_i$。刚度 $K_{ij}=\int_{-1}^1 B_i^T D B_j J\,d\xi$，用高斯求积
$$
\int_{-1}^1 f(\xi)d\xi\approx\sum_{k=1}^n w_kf(\xi_k).
$$
$n$ 点高斯精确积分 $2n-1$ 次多项式。

物理意义：适应曲边单元。

特例：2 点高斯精确到三次。应用：二维/三维单元。

**知识点6：收敛性与完备性**

定义：单元须满足完备性与协调性。

推导：Patch test 要求拼接单元能精确再现常应变状态；完备性要求形函数含常数项与线性项 $\sum N_i=1$，$\sum N_i x_i=x$。协调性要求跨单元场连续。

物理意义：保证网格加密收敛到真解。

特例：非协调单元需通过 patch test。应用：单元选择。

**知识点7：动力学有限元**

定义：$\mathbf M\ddot{\mathbf U}+\mathbf C\dot{\mathbf U}+\mathbf K\mathbf U=\mathbf F$。

推导：动能 $T=\frac12\dot{\mathbf U}^T\mathbf M\dot{\mathbf U}$，$\mathbf M=\int\rho\mathbf N^T\mathbf N\,dV$（一致质量）。由 Hamilton 原理 $\delta\int(T-V+W)dt=0$ 得运动方程。模态分析 $(\mathbf K-\omega^2\mathbf M)\boldsymbol\phi=0$。

物理意义：结构动力学数值求解。

特例：集中质量矩阵对角化。应用：模态、瞬态响应。

**知识点8：显式与隐式时间积分**

定义：显式中心差分，隐式 Newmark-β。

推导（中心差分）：$\ddot{\mathbf U}^n=\dfrac{\mathbf U^{n+1}-2\mathbf U^n+\mathbf U^{n-1}}{\Delta t^2}$，$\dot{\mathbf U}^n=\dfrac{\mathbf U^{n+1}-\mathbf U^{n-1}}{2\Delta t}$。代入运动方程解 $\mathbf U^{n+1}$，条件稳定，$\Delta t\le\Delta t_{cr}\approx\dfrac2{\omega_{\max}}$。隐式 Newmark：$\mathbf U^{n+1}$ 需解线性方程组，无条件稳定（$\beta\ge1/4,\gamma=1/2$）。

物理意义：显式适合冲击，隐式适合长期响应。

特例：稳定极限由最小单元尺寸控制。应用：碰撞、地震分析。

**知识点9：数值误差与网格收敛**

定义：离散误差 $\|\mathbf u_h-\mathbf u\|\to0$（$h\to0$）。

推导：对线性单元，能量范数误差 $O(h)$，$L^2$ 范数 $O(h^2)$（Aubin-Nitsche）。加密网格使误差按幂律下降。

物理意义：网格无关性验证。

特例：奇点处收敛降阶。应用：结果可信度评估。

**知识点10：有限元在工程实践I中的流程**

定义：前处理—求解—后处理。

推导：流程为几何建模→材料/单元赋值→网格→边界与载荷→组装求解→后处理。每一步误差累积由验证与确认（V&V）控制。

物理意义：工程仿真规范化。

特例：简单悬臂梁与解析解对比验证。应用：综合实验项目。

---
⬆ 上级领域：[[实践与创新]]
