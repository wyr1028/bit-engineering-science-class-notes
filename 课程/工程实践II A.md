---
课程名称: 工程实践II A
学分: 1
总学时: 24
考核方式: 考查
课程属性: 必修
课程性质: 实践环节
开课单位: 空天科学与技术学院
领域: 实践与创新
tags:
  - 课程
  - 领域/实践与创新
---

# 工程实践II A


### 第1章 力学建模、量纲分析与相似准则

**知识点1：连续介质假设与尺度分离**
定义：将离散分子构成的物质视为充满空间、物理量（密度、速度、应力）连续可微的连续体，称连续介质。
假设：所研究问题的特征长度 $L$ 远大于分子平均自由程 $\lambda$，即 $\lambda/L\to 0$；且远大于分子间距。对空气在常温常压下 $\lambda\sim 7\times10^{-8}\,\mathrm m$。
判据：克努森数 $\mathrm{Kn}=\lambda/L$。当 $\mathrm{Kn}<0.01$ 时可视为连续介质。
推导连续性条件的几何意义：取体积元 $\Delta V$，平均密度 $\rho_{\text{avg}}=\Delta m/\Delta V$；令 $\Delta V$ 缩小到特征体积 $V^\*$，若在 $V^\*> \lambda^3$ 的范围内 $\rho_{\text{avg}}$ 趋于稳定极限 $\rho$，则
$$\rho=\lim_{\Delta V\to V^\*}\frac{\Delta m}{\Delta V}=\frac{\mathrm dm}{\mathrm dV}.$$
$\because$ 分子数 $N\sim n\Delta V$，相对涨落 $\sim 1/\sqrt N\to 0$ 当 $\Delta V\gg\lambda^3$，$\therefore$ 极限稳定。
意义：连续场量可写成空间坐标函数，才能使用微积分与场论。
极限特例：$\mathrm{Kn}\to\infty$（真空稀薄气体）连续假设失效，需用分子动理论或 DSMC。
应用：气动外形设计、结构应力分析的建模前提。

**知识点2：量纲与基本量纲**
定义：物理量关于基本量（长度 $\mathsf L$、质量 $\mathsf M$、时间 $\mathsf T$、温度 $\mathsf\Theta$）的幂次表达称量纲。
公式：导出量量纲可写为 $[q]=\mathsf M^a\mathsf L^b\mathsf T^c\mathsf\Theta^d$。
推导示例：由 $F=ma$，
$$[F]=[m]\,[a]=\mathsf M\cdot\mathsf L\mathsf T^{-2}=\mathsf{MLT}^{-2}.$$
由 $W=Fs$，$[W]=\mathsf{MLT}^{-2}\cdot\mathsf L=\mathsf{ML}^2\mathsf T^{-2}$。
由应力 $\sigma=F/A$，$[\sigma]=\mathsf{MLT}^{-2}/\mathsf L^2=\mathsf{ML}^{-1}\mathsf T^{-2}$。
意义：量纲是物理量本质属性，与单位制无关。
极限特例：无量纲量量纲为 $\mathsf M^0\mathsf L^0\mathsf T^0=1$。
应用：检查公式、换算单位、导出经验关系。

**知识点3：量纲齐次性原理**
定义：任何完备物理方程的等号两端、以及加减各项，必须具有相同量纲。
推导：设物理规律 $f(q_1,\dots,q_n)=0$。若将基本单位缩小 $\alpha,\beta,\gamma$ 倍，则每个 $q_i$ 数值改变，而方程形式不变。由于物理规律与单位选择无关，
$$f(\alpha^{a_1}\beta^{b_1}\gamma^{c_1}q_1,\dots)=0\ \Longleftrightarrow\ f(q_1,\dots)=0.$$
$\therefore$ 各项量纲指数须平衡。
意义：可用来初判公式对错；不能给出完整系数，但能定出幂次组合。
极限特例：超越函数（$\sin,\ln,\exp$）的自变量必须无量纲；$\because$ 可按泰勒展开为幂级数，非无量纲项无法相加。
应用：判断 $\sigma=E\varepsilon$ 成立而对 $\sigma=E\varepsilon^2$ 不成立（后者量纲为 $\mathsf{ML}^{-1}\mathsf T^{-2}$ 之外的量纲）。

**知识点4：Buckingham Π 定理**
定义：若 $n$ 个物理量满足量纲齐次关系，其中含 $k$ 个基本量纲，则存在 $n-k$ 个独立无量纲组合 $\Pi_1,\dots,\Pi_{n-k}$，且
$$f(\Pi_1,\Pi_2,\dots,\Pi_{n-k})=0.$$
推导：取 $k$ 个“基本变量”作为重复变量（其量纲矩阵秩为 $k$）；其余每个变量用重复变量的幂次消除量纲：
$$\because [q_j]=\prod_{i=1}^{k}[b_i]^{m_{ij}},\qquad \Pi_j=q_j\prod_{i=1}^{k}b_i^{-m_{ij}},$$
$$\therefore [\Pi_j]=1\ (j=k+1,\dots,n).$$
列出量纲指数线性方程求解 $m_{ij}$，得到无量纲数群。
意义：可把多变量问题降维，如把一个含 5 个变量的问题化为含 2 个无量纲群的关系。
极限特例：若 $n=k$，则 $\Pi$ 数为 0，规律唯一。
应用：$P=f(\rho,V,D,\mu)$ 化为 $P/\rho V^2D^2=\phi(\mathrm{Re})$。

**知识点5：常见相似准则的推导**
定义：相似准则为控制流动/变形相似的无量纲数。
推导：
（1）雷诺数：惯性力与黏性力之比
$$\because F_i\sim\rho V^2L^2,\quad F_\mu\sim\mu VL,\qquad \therefore \mathrm{Re}=\frac{F_i}{F_\mu}=\frac{\rho VL}{\mu}=\frac{VL}{\nu}.$$
（2）弗劳德数：惯性力与重力之比
$$\mathrm{Fr}=\frac{V}{\sqrt{gL}},\qquad \because F_g\sim\rho gL^3.$$
（3）马赫数：惯性力与弹性（可压缩）力之比，用体积模量 $K=\rho c^2$，
$$\mathrm{Ma}=\frac{V}{c}.$$
（4）柯西数：惯性力与弹性力之比，$\mathrm{Ca}=\rho V^2/E$。
（5）斯托克斯数：颗粒松弛时间与流场特征时间之比，$\mathrm{St}=\tau_p/t_f=\rho_pd_p^2V/(18\mu L)$。
意义：模型与原型同名无量纲数相等才相似。
极限特例：$\mathrm{Re}\to0$ 黏性主导（斯托克斯流）；$\mathrm{Re}\to\infty$ 惯性主导（理想流）。
应用：风洞缩比试验、水洞流场显示。

**知识点6：无量纲化与特征尺度选取**
定义：用特征量（$L,V,\rho$）将控制方程中各量化为无量纲量，如 $\hat x=x/L,\ \hat u=u/V,\ \hat t=tV/L$。
推导（以不可压 Navier–Stokes 为例）：方程
$$\rho\!\left(\frac{\partial\mathbf u}{\partial t}+\mathbf u\!\cdot\!\nabla\mathbf u\right)=-\nabla p+\mu\nabla^2\mathbf u.$$
代入 $x=L\hat x,\ \mathbf u=V\hat{\mathbf u},\ t=(L/V)\hat t,\ p=\rho V^2\hat p$：
$$\rho\frac{V^2}{L}\!\left(\frac{\partial\hat{\mathbf u}}{\partial\hat t}+\hat{\mathbf u}\!\cdot\!\hat\nabla\hat{\mathbf u}\right)=-\frac{\rho V^2}{L}\hat\nabla\hat p+\mu\frac{V}{L^2}\hat\nabla^2\hat{\mathbf u}.$$
两端同乘 $L/(\rho V^2)$：
$$\frac{\partial\hat{\mathbf u}}{\partial\hat t}+\hat{\mathbf u}\!\cdot\!\hat\nabla\hat{\mathbf u}=-\hat\nabla\hat p+\frac{1}{\mathrm{Re}}\hat\nabla^2\hat{\mathbf u}.$$
意义：显式暴露控制参数 $1/\mathrm{Re}$，便于数值与实验比较。
极限特例：$\mathrm{Re}\gg1$ 忽略黏性项；$\mathrm{Re}\ll1$ 忽略惯性项。
应用：判断缩比模型能否同时满足多个相似数（常发生冲突）。

**知识点7：模型相似设计（风洞/水洞缩比）**
定义：使模型与原型所有相关的 $\Pi$ 数相等，以复现物理现象。
推导：设原型 $p$、模型 $m$，同种介质需 $\mathrm{Re}_p=\mathrm{Re}_m$：
$$\frac{V_pL_p}{\nu_p}=\frac{V_mL_m}{\nu_m}\ \Rightarrow\ V_m=V_p\frac{L_p}{L_m}\frac{\nu_m}{\nu_p}.$$
若缩比 $\lambda=L_m/L_p<1$，则风速需放大 $1/\lambda$。若同时要求 $\mathrm{Ma}$ 相等，则 $V_m=V_p$，与上式矛盾，$\therefore$ 需变密度风洞或不同介质。
意义：相似设计是在多个无量纲数不能同时满足时的折中。
极限特例：$\lambda\to1$ 为全尺寸试验。
应用：飞行器测力、桥梁气弹模型、船模阻力试验（常用 $\mathrm{Fr}$ 相似）。

**知识点8：量纲分析在结构振动中的应用**
定义：用 $\Pi$ 定理导出固有频率的标度律。
推导：单自由度振子 $\omega_n=\sqrt{k/m}$。若更一般地由 $k=EA/L$，则
$$\because \omega_n=\sqrt{\frac{EA/L}{m}}=\sqrt{\frac{E}{\rho}}\frac{1}{L},\qquad \therefore \frac{\omega_n L}{c}=1,\ c=\sqrt{E/\rho}.$$
对几何相似结构，$\omega_n\propto(1/L)\sqrt{E/\rho}$。若 $L\to\lambda L$，则固有频率按 $1/\lambda$ 缩放。
意义：波速 $c=\sqrt{E/\rho}$ 为材料特征速度；时间尺度 $L/c$。
极限特例：相同材料不同尺寸，频率严格成反比于尺寸。
应用：缩比模型模态试验预测全机频率；地震模拟振动台设计。

### 第2章 质点与质点系动力学

**知识点1：牛顿运动定律与惯性系**
定义：惯性系中，不受外力时质点保持静止或匀速直线运动；受力时
$$\mathbf F=m\mathbf a=m\frac{\mathrm d\mathbf v}{\mathrm dt}=\frac{\mathrm d\mathbf p}{\mathrm dt},\qquad \mathbf p=m\mathbf v.$$
假设：质量 $m$ 为常数、参考系为惯性系。
推导：由动量定义 $\mathbf p=m\mathbf v$，若质量恒定，
$$\frac{\mathrm d\mathbf p}{\mathrm dt}=m\frac{\mathrm d\mathbf v}{\mathrm dt}=m\mathbf a.$$
当 $m=m(t)$（变质量）时：
$$\mathbf F=\frac{\mathrm d(m\mathbf v)}{\mathrm dt}=m\frac{\mathrm d\mathbf v}{\mathrm dt}+\mathbf v\frac{\mathrm dm}{\mathrm dt}.$$
意义：力是改变动量（而非维持运动）的原因。
极限特例：$\mathbf F=0\Rightarrow\mathbf p=$ 常量（动量守恒）。
应用：牛顿第二定律建立质点运动微分方程。

**知识点2：动量定理与守恒**
定义：质点系总动量 $\mathbf P=\sum_i m_i\mathbf v_i$；冲量 $\mathbf I=\int_{t_1}^{t_2}\mathbf F_{\text{外}}\,\mathrm dt$。
推导：对第 $i$ 个质点，内力成对出现 $\mathbf F_{ij}=-\mathbf F_{ji}$，求和后
$$\because \sum_i\mathbf F_{ij}=0,\qquad \therefore \frac{\mathrm d\mathbf P}{\mathrm dt}=\sum_i\mathbf F_i^{\text{外}}.$$
积分得动量定理：
$$\Delta\mathbf P=\int_{t_1}^{t_2}\mathbf F^{\text{外}}\,\mathrm dt=\mathbf I.$$
意义：内力不改总动量；只有外力冲量改变总动量。
极限特例：$\mathbf F^{\text{外}}=0\Rightarrow\mathbf P=$ 常量。
应用：碰撞、爆炸、火箭反冲。

**知识点3：角动量定理与守恒**
定义：对原点角动量 $\mathbf L=\sum_i\mathbf r_i\times m_i\mathbf v_i$；力矩 $\mathbf M=\sum_i\mathbf r_i\times\mathbf F_i^{\text{外}}$。
推导：
$$\frac{\mathrm d\mathbf L}{\mathrm dt}=\sum_i\!\left(\frac{\mathrm d\mathbf r_i}{\mathrm dt}\times m_i\mathbf v_i+\mathbf r_i\times m_i\frac{\mathrm d\mathbf v_i}{\mathrm dt}\right).$$
$\because \dfrac{\mathrm d\mathbf r_i}{\mathrm dt}\times m_i\mathbf v_i=\mathbf v_i\times m_i\mathbf v_i=0$，
$$\therefore \frac{\mathrm d\mathbf L}{\mathrm dt}=\sum_i\mathbf r_i\times\mathbf F_i^{\text{外}}+\sum_{i<j}\big[(\mathbf r_i-\mathbf r_j)\times\mathbf F_{ij}\big].$$
又 $\because$ 内力沿 $ij$ 连线（中心力），$(\mathbf r_i-\mathbf r_j)\parallel\mathbf F_{ij}\Rightarrow$ 内力矩之和为 0，$\therefore$
$$\frac{\mathrm d\mathbf L}{\mathrm dt}=\mathbf M^{\text{外}}.$$
意义：内力矩不改变总角动量。
极限特例：$\mathbf M^{\text{外}}=0\Rightarrow\mathbf L=$ 常量。
应用：开普勒第二定律、离心机、陀螺运动。

**知识点4：动能定理与功**
定义：动能 $T=\tfrac12 m v^2$；功 $W=\int_1^2\mathbf F\!\cdot\!\mathrm d\mathbf r$。
推导：由 $\mathbf F=m\dot{\mathbf v}$，两边点乘 $\mathrm d\mathbf r=\mathbf v\,\mathrm dt$：
$$\mathbf F\!\cdot\!\mathrm d\mathbf r=m\dot{\mathbf v}\!\cdot\!\mathbf v\,\mathrm dt=\frac{m}{2}\frac{\mathrm d(\mathbf v\!\cdot\!\mathbf v)}{\mathrm dt}\mathrm dt=\mathrm d\!\left(\frac12 mv^2\right).$$
积分：
$$W_{12}=\int_1^2\mathbf F\!\cdot\!\mathrm d\mathbf r=\tfrac12mv_2^2-\tfrac12mv_1^2=T_2-T_1.$$
意义：合力做功等于动能增量。
极限特例：$W=0$ 时速率不变（如向心力做功用）。
应用：求变力做功、能量守恒建立。

**知识点5：势能与保守力**
定义：若 $\oint\mathbf F\!\cdot\!\mathrm d\mathbf r=0$（等价 $\nabla\times\mathbf F=0$），则存在势能 $V$ 使 $\mathbf F=-\nabla V$。
推导：由 $\mathbf F=-\nabla V$，
$$W_{12}=\int_1^2-\nabla V\!\cdot\!\mathrm d\mathbf r=-\int_1^2\mathrm dV=V_1-V_2.$$
$\therefore$ 保守力做功只与端点有关。
常见势能：重力 $V=mgz$；弹性 $V=\tfrac12kx^2$；万有引力 $V=-GMm/r$。
意义：保守系可引入标量势，简化动力学。
极限特例：摩擦力非保守，$\oint\mathbf F\!\cdot\!\mathrm d\mathbf r\ne0$，无（单值）势能。
应用：机械能守恒 $T+V=E$。

**知识点6：质心运动定理**
定义：质心 $\mathbf r_C=\dfrac{\sum m_i\mathbf r_i}{\sum m_i}$，总质量 $M=\sum m_i$。
推导：
$$\because \mathbf P=\sum m_i\dot{\mathbf r}_i=M\dot{\mathbf r}_C,\qquad \therefore \frac{\mathrm d\mathbf P}{\mathrm dt}=M\ddot{\mathbf r}_C=\mathbf F^{\text{外}}.$$
又
$$\sum m_i\mathbf r_i\times\ddot{\mathbf r}_i=\frac{\mathrm d}{\mathrm dt}\!\left(\sum m_i\mathbf r_i\times\dot{\mathbf r}_i\right)-\underbrace{\sum m_i\dot{\mathbf r}_i\times\dot{\mathbf r}_i}_{0},$$
取质心为矩心即得质心系角动量定理。
意义：质点系平动可等效为质量集中在质心的质点。
极限特例：刚体/$n$ 体系统均适用。
应用：火箭整体轨迹、碰撞问题。

**知识点7：变质量系统与火箭方程**
定义：系统质量随时间变化，如火箭喷出燃气。取 $t$ 时系统质量 $m$，速度 $v$；$t+\mathrm dt$ 喷出质量 $-\mathrm dm$（$\mathrm dm<0$）。
假设：外场力 $\mathbf F_{\text{外}}$，喷流相对火箭速度 $\mathbf u$。
推导：由动量守恒（系统含火箭+已喷气体），
$$m\mathbf v+\mathrm dm\,(\mathbf v-\mathbf u)=(m+\mathrm dm)(\mathbf v+\mathrm d\mathbf v).$$
展开略去二阶小量 $\mathrm dm\,\mathrm d\mathbf v$：
$$m\mathrm d\mathbf v+\mathbf u\,\mathrm dm=0\ \Rightarrow\ m\,\mathrm d\mathbf v=-\mathbf u\,\mathrm dm.$$
有外力时
$$m\frac{\mathrm d\mathbf v}{\mathrm dt}=\mathbf F_{\text{外}}-\mathbf u\frac{\mathrm dm}{\mathrm dt}.$$
无外力、$\mathbf u$ 与 $\mathbf v$ 反平行时积分：
$$v-v_0=u\ln\frac{m_0}{m}.$$
意义：质量比是航天速度增量的核心。
极限特例：$m\to0$ 时 $\Delta v\to\infty$（实际受燃料限制）。
应用：齐奥尔科夫斯基公式、多级火箭。

**知识点8：两体问题与约化质量**
定义：两质点只在相互引力下运动。
推导：质心运动无外力 $\Rightarrow$ 匀速。相对位矢 $\mathbf r=\mathbf r_1-\mathbf r_2$，由牛顿方程
$$m_1\ddot{\mathbf r}_1=\mathbf F_{12},\quad m_2\ddot{\mathbf r}_2=-\mathbf F_{12}.$$
两式分别除以 $m_1,m_2$ 相减：
$$\ddot{\mathbf r}=\mathbf F_{12}\!\left(\frac{1}{m_1}+\frac{1}{m_2}\right)=\frac{\mathbf F_{12}}{\mu},\qquad \mu=\frac{m_1m_2}{m_1+m_2}.$$
对引力 $\mathbf F_{12}=-\dfrac{Gm_1m_2}{r^2}\hat{\mathbf r}$，得
$$\mu\ddot{\mathbf r}=-\frac{Gm_1m_2}{r^2}\hat{\mathbf r}\ \Rightarrow\ \ddot{\mathbf r}=-\frac{G(m_1+m_2)}{r^2}\hat{\mathbf r}.$$
意义：两体问题化为单体在等效中心场中运动。
极限特例：$m_2\gg m_1\Rightarrow\mu\approx m_1$，等效为小质点绕固定中心。
应用：行星轨道、双星系统。

### 第3章 刚体静力学、约束与虚功原理

**知识点1：力系简化与主矢主矩**
定义：空间力系 $\{\mathbf F_i\}$ 对参考点 $O$ 的主矢 $\mathbf R=\sum\mathbf F_i$，主矩 $\mathbf M_O=\sum\mathbf r_i\times\mathbf F_i$。
推导：将每个力移到点 $O$ 并附加力偶 $\mathbf r_i\times\mathbf F_i$，
$$\therefore \text{力系}\equiv\ \{\mathbf R;\ \mathbf M_O\}.$$
换矩心 $O'$：$\mathbf r_i=\mathbf r_{O'}+\mathbf r_i'$，
$$\mathbf M_O=\sum(\mathbf r_{O'}+\mathbf r_i')\times\mathbf F_i=\mathbf M_{O'}+\mathbf r_{O'}\times\mathbf R.$$
意义：主矢与矩心无关，主矩随矩心变化。
极限特例：$\mathbf R=0,\ \mathbf M_O=0$ 为平衡；$\mathbf R=0,\ \mathbf M_O\ne0$ 为纯力偶。
应用：桁架、支座反力简化。

**知识点2：刚体平衡方程**
定义：刚体平衡条件为 $\mathbf R=0,\ \mathbf M_O=0$。
推导（平面力系）：分解到坐标轴
$$\sum F_x=0,\qquad \sum F_y=0,\qquad \sum M_O=0.$$
对空间力系为六式：
$$\sum F_x=\sum F_y=\sum F_z=0,\qquad \sum M_x=\sum M_y=\sum M_z=0.$$
意义：三（六）个独立方程。
极限特例：静定结构未知反力数等于方程数；否则静不定。
应用：梁、框架约束反力求解。

**知识点3：约束、自由度与广义坐标**
定义：限制物体运动的条件为约束；独立坐标数为自由度数 $s$。
推导：$N$ 个质点无约束时有 $3N$ 个坐标，若有 $r$ 个完整约束 $\phi_k(\mathbf r_1,\dots,t)=0$，则
$$s=3N-r.$$
完整约束可用来消去 $r$ 个坐标，取 $q_1,\dots,q_s$ 为广义坐标。
意义：广义坐标可任意选取，不必是长度。
极限特例：刚体空间运动 $s=6$（3 平动 + 3 转动）。
应用：机械臂、机构运动学自由度判定。

**知识点4：虚位移与虚功原理**
定义：在给定瞬时、约束允许的无穷小位移 $\delta\mathbf r_i$ 称虚位移；虚功 $\delta W=\sum\mathbf F_i\!\cdot\!\delta\mathbf r_i$。
推导：理想约束下约束力做功之和为零，
$$\sum\mathbf N_i\!\cdot\!\delta\mathbf r_i=0.$$
质点系平衡 $\mathbf F_i+\mathbf N_i=0$，点乘 $\delta\mathbf r_i$ 求和：
$$\sum\mathbf F_i\!\cdot\!\delta\mathbf r_i+\underbrace{\sum\mathbf N_i\!\cdot\!\delta\mathbf r_i}_{0}=0\ \Rightarrow\ \delta W=\sum\mathbf F_i\!\cdot\!\delta\mathbf r_i=0.$$
反过来，若对任意虚位移 $\delta W=0$，由 $\delta\mathbf r_i$ 独立任意，$\therefore\mathbf F_i+\mathbf N_i=0$，即平衡。
意义：虚功原理把静力学问题化为（广义）力的平衡，可消去理想约束力。
极限特例：有摩擦力时约束非理想，须把摩擦力计入主动力。
应用：机构、桁架、静力求解。

**知识点5：静定与静不定判定**
定义：未知力数 $u$ 与独立平衡方程数 $e$ 比较：$u=e$ 静定，$u>e$ 静不定，次数 $u-e$。
推导（平面桁架）：节点数 $j$、杆数 $b$、支座反力 $r$，方程 $2j$ 个，未知 $b+r$ 个，
$$n_{\text{ind}}=b+r-2j.$$
意义：静不定须引入变形协调或能量法。
极限特例：$n_{\text{ind}}=0$ 静定；$n_{\text{ind}}>0$ 超静定。
应用：判定梁、刚架的可解性。

**知识点6：摩擦与极限平衡**
定义：最大静摩擦 $F_f^{\max}=\mu_s N$；滑动摩擦 $F_f=\mu_k N$，$\mu_k<\mu_s$。
推导：斜面滑块开始下滑的临界角 $\alpha_c$：平衡时
沿斜面 $mg\sin\alpha=F_f$，法向 $N=mg\cos\alpha$，
临界 $F_f=\mu_s N$：
$$\therefore \tan\alpha_c=\mu_s.$$
意义：摩擦角 $\phi=\arctan\mu_s$ 表征材料自锁。
极限特例：$\alpha<\alpha_c$ 静止；$\alpha>\alpha_c$ 加速下滑。
应用：螺纹自锁、边坡稳定。

**知识点7：桁架分析（节点法与截面法）**
定义：节点法对每个节点列二力平衡；截面法用假想截面截断桁架列平衡。
推导：节点 $i$ 处二力平衡
$$\sum F_x=0,\quad\sum F_y=0.$$
截面法对切开的隔离体列
$$\sum M=0\ \Rightarrow\ \text{求特定杆内力}.$$
意义：零力杆判定可先消元（两杆节点无外力则均为零）。
极限特例：静定桁架内力唯一。
应用：桥梁、塔架。

**知识点8：势能驻值原理（静力平衡）**
定义：保守系势能 $V(\mathbf q)$，平衡条件 $\delta V=0$。
推导：虚功 $\delta W=-\delta V=0$，对任意 $\delta\mathbf q$，
$$\therefore \frac{\partial V}{\partial q_i}=0,\quad i=1,\dots,s.$$
二阶判据：若 Hessian $\big(\partial^2V/\partial q_i\partial q_j\big)$ 正定，则势能极小、稳定；若负定，则不稳定；若有零特征值需更高阶分析。
意义：静力平衡 = 势能驻值。
极限特例：线性弹簧 $V=\tfrac12kx^2$，$V'=kx=0\Rightarrow x=0$ 稳定。
应用：结构稳定、机构平衡位形。

### 第4章 拉格朗日力学

**知识点1：广义坐标与位形空间**
定义：用 $s$ 个独立广义坐标 $q_1,\dots,q_s$ 描述系统位形；$\{\mathbf q\}$ 构成 $s$ 维位形空间。
推导：质点位置写成
$$\mathbf r_i=\mathbf r_i(q_1,\dots,q_s,t).$$
虚位移
$$\delta\mathbf r_i=\sum_{k=1}^{s}\frac{\partial\mathbf r_i}{\partial q_k}\delta q_k.$$
意义：约束通过坐标选取自动满足。
极限特例：刚体空间运动 $s=6$。
应用：多体系统建模。

**知识点2：达朗贝尔原理**
定义：对每个质点引入惯性力 $-m_i\ddot{\mathbf r}_i$，则动态问题可写为“静力”平衡：
$$\mathbf F_i+\mathbf N_i-m_i\ddot{\mathbf r}_i=0.$$
推导：牛顿第二定律 $m_i\ddot{\mathbf r}_i=\mathbf F_i+\mathbf N_i$ 移项得
$$\mathbf F_i+\mathbf N_i-m_i\ddot{\mathbf r}_i=0.$$
点乘 $\delta\mathbf r_i$ 求和，理想约束消去：
$$\sum_i(\mathbf F_i-m_i\ddot{\mathbf r}_i)\!\cdot\!\delta\mathbf r_i=0.$$
意义：把动力学化为虚功（达朗贝尔–拉格朗日原理）。
极限特例：静力学时惯性力为零。
应用：推导拉格朗日方程。

**知识点3：拉格朗日方程推导**
定义：动能 $T=\tfrac12\sum m_i\dot{\mathbf r}_i^2$；广义力 $Q_k=\sum_i\mathbf F_i\!\cdot\!\partial\mathbf r_i/\partial q_k$。
推导：由达朗贝尔原理
$$\sum_i m_i\ddot{\mathbf r}_i\!\cdot\!\delta\mathbf r_i=\sum_k Q_k\,\delta q_k.$$
左边代入 $\delta\mathbf r_i=\sum_k(\partial\mathbf r_i/\partial q_k)\delta q_k$：
$$\sum_k\!\left[\sum_i m_i\ddot{\mathbf r}_i\!\cdot\!\frac{\partial\mathbf r_i}{\partial q_k}\right]\delta q_k=\sum_k Q_k\delta q_k.$$
证明两个恒等式：
$$(\mathrm i)\ \frac{\partial\dot{\mathbf r}_i}{\partial\dot q_k}=\frac{\partial\mathbf r_i}{\partial q_k};\qquad
(\mathrm{ii})\ \frac{\mathrm d}{\mathrm dt}\frac{\partial\mathbf r_i}{\partial q_k}=\frac{\partial\dot{\mathbf r}_i}{\partial q_k}.$$
于是
$$\sum_i m_i\ddot{\mathbf r}_i\!\cdot\!\frac{\partial\mathbf r_i}{\partial q_k}=\frac{\mathrm d}{\mathrm dt}\!\left(\frac{\partial T}{\partial\dot q_k}\right)-\frac{\partial T}{\partial q_k}.$$
$\because\delta q_k$ 独立，$\therefore$
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial T}{\partial\dot q_k}-\frac{\partial T}{\partial q_k}=Q_k.$$
对保守力 $Q_k=-\partial V/\partial q_k$，令 $L=T-V$：
$$\boxed{\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_k}-\frac{\partial L}{\partial q_k}=0.}$$
意义：$s$ 个二阶常微分方程，自动消去约束力。
极限特例：$V=0$ 时 $\dfrac{\mathrm d}{\mathrm dt}(m\dot q)=0$（动量守恒）。
应用：多体动力学、机器人、天体轨道。

**知识点4：非保守力与耗散函数**
定义：含耗散力 $Q_k^{\text{d}}=-\partial\mathcal D/\partial\dot q_k$，瑞利耗散函数 $\mathcal D=\tfrac12\sum c_{ij}\dot q_i\dot q_j$。
推导：线性阻尼 $\mathbf F_i=-c_i\dot{\mathbf r}_i$，则
$$Q_k^{\text{d}}=\sum_i(-c_i\dot{\mathbf r}_i)\!\cdot\!\frac{\partial\mathbf r_i}{\partial q_k}=-\frac{\partial}{\partial\dot q_k}\!\left(\frac12\sum_i c_i\dot{\mathbf r}_i^2\right)=-\frac{\partial\mathcal D}{\partial\dot q_k}.$$
拉格朗日方程推广：
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_k}-\frac{\partial L}{\partial q_k}+\frac{\partial\mathcal D}{\partial\dot q_k}=Q_k^{\text{其余}}.$$
意义：阻尼力可由标量函数导出。
极限特例：$\mathcal D=0$ 回到保守系。
应用：振动衰减、结构阻尼建模。

**知识点5：循环坐标与守恒律**
定义：若 $L$ 不显含某坐标 $q_j$（$\partial L/\partial q_j=0$），称 $q_j$ 为循环坐标。
推导：拉格朗日方程
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}-\underbrace{\frac{\partial L}{\partial q_j}}_{0}=0\ \Rightarrow\ \frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}=0.$$
$$\therefore p_j\equiv\frac{\partial L}{\partial\dot q_j}=\text{常量}.$$
意义：坐标不显含 → 共轭动量守恒。
极限特例：无力场中 $L=\tfrac12m\dot x^2$，$p_x=m\dot x$ 守恒。
应用：中心力问题角动量守恒；轴对称问题。

**知识点6：中心力场与开普勒问题**
定义：$L=\tfrac12m(\dot r^2+r^2\dot\theta^2)-V(r)$。
推导：$\theta$ 为循环坐标，
$$p_\theta=mr^2\dot\theta=l=\text{常量}.$$
径向方程
$$m\ddot r-mr\dot\theta^2=-\frac{\mathrm dV}{\mathrm dr}\ \Rightarrow\ m\ddot r-\frac{l^2}{mr^3}=F(r).$$
对引力 $V=-GMm/r$，得
$$m\ddot r=-\frac{GMm}{r^2}+\frac{l^2}{mr^3}.$$
有效势 $V_{\text{eff}}=V(r)+\dfrac{l^2}{2mr^2}$，轨道为圆锥曲线
$$r=\frac{p}{1+e\cos\theta},\qquad p=\frac{l^2}{GMm^2}.$$
意义：平方反比力场必有闭合椭圆轨道。
极限特例：$e=0$ 圆轨道；$e<1$ 椭圆；$e=1$ 抛物；$e>1$ 双曲。
应用：卫星轨道、行星运动。

**知识点7：拉格朗日乘子与约束**
定义：对非完整或不便消去的约束 $\phi_j(\mathbf q,t)=0$，引入乘子 $\lambda_j$。
推导：约束力不出功于虚位移，$\delta\phi_j=\sum_k\partial\phi_j/\partial q_k\,\delta q_k=0$。把 $\sum_j\lambda_j\delta\phi_j=0$ 加到虚功方程：
$$\sum_k\!\left[\frac{\mathrm d}{\mathrm dt}\frac{\partial T}{\partial\dot q_k}-\frac{\partial T}{\partial q_k}-Q_k+\sum_j\lambda_j\frac{\partial\phi_j}{\partial q_k}\right]\delta q_k=0.$$
$\because\delta q_k$ 不再独立，$\therefore$ 各系数等于 $\sum_j\lambda_j\partial\phi_j/\partial q_k$（不直接为零），联立约束方程求解 $q_k,\lambda_j$。
意义：$\lambda_j$ 与约束力同构，可求约束反力。
极限特例：无约束时 $\lambda_j=0$。
应用：含刚性连接的机械系统、滚动约束。

**知识点8：小振动与简正模**
定义：稳定平衡点附近 $L\approx\tfrac12\dot{\mathbf q}^{\mathsf T}M\dot{\mathbf q}-\tfrac12\mathbf q^{\mathsf T}K\mathbf q$。
推导：设平衡 $\mathbf q_0$，令 $\mathbf q=\mathbf q_0+\boldsymbol\eta$，展开 $V$：
$$V(\mathbf q_0+\boldsymbol\eta)=V(\mathbf q_0)+\underbrace{\nabla V(\mathbf q_0)}_{0}\!\cdot\!\boldsymbol\eta+\tfrac12\boldsymbol\eta^{\mathsf T}\frac{\partial^2V}{\partial q_i\partial q_j}\boldsymbol\eta+\cdots$$
取 $K_{ij}=\partial^2V/\partial q_i\partial q_j$，$M_{ij}=\partial^2T/\partial\dot q_i\partial\dot q_j$（在 $\mathbf q_0$）。拉格朗日方程给出
$$M\ddot{\boldsymbol\eta}+K\boldsymbol\eta=0.$$
设 $\boldsymbol\eta=\mathbf a e^{\mathrm i\omega t}$：
$$(K-\omega^2M)\mathbf a=0\ \Rightarrow\ \det(K-\omega^2M)=0.$$
意义：特征值 $\omega_r^2$ 为固有频率，特征矢量为简正模。
极限特例：单自由度为 $\omega=\sqrt{k/m}$。
应用：模态分析、隔振设计。

### 第5章 哈密顿力学与正则方程

**知识点1：勒让德变换**
定义：由 $L(\mathbf q,\dot{\mathbf q},t)$ 定义共轭动量 $p_k=\partial L/\partial\dot q_k$；哈密顿量
$$H(\mathbf q,\mathbf p,t)=\sum_k p_k\dot q_k-L.$$
推导：取 $L$ 的全微分
$$\mathrm dL=\sum_k\!\left(\frac{\partial L}{\partial q_k}\mathrm dq_k+\frac{\partial L}{\partial\dot q_k}\mathrm d\dot q_k\right)+\frac{\partial L}{\partial t}\mathrm dt=\sum_k(\dot p_k\,\mathrm dq_k+p_k\,\mathrm d\dot q_k)+\frac{\partial L}{\partial t}\mathrm dt.$$
由 $H=\sum p_k\dot q_k-L$，
$$\mathrm dH=\sum_k(p_k\,\mathrm d\dot q_k+\dot q_k\,\mathrm dp_k)-\mathrm dL=\sum_k(\dot q_k\,\mathrm dp_k-\dot p_k\,\mathrm dq_k)-\frac{\partial L}{\partial t}\mathrm dt.$$
意义：$(q,p)$ 为独立变量，$L$ 的 Lagrange 描述转为 Hamilton 描述。
极限特例：$L$ 不显含 $t$ 时 $H$ 守恒。
应用：相空间分析、量子力学对应。

**知识点2：哈密顿正则方程**
定义：正则方程
$$\dot q_k=\frac{\partial H}{\partial p_k},\qquad \dot p_k=-\frac{\partial H}{\partial q_k}.$$
推导：由 $\mathrm dH$ 表达式比较系数即得 $\dot q_k=\partial H/\partial p_k,\ \dot p_k=-\partial H/\partial q_k$。又
$$\frac{\partial H}{\partial t}=-\frac{\partial L}{\partial t}.$$
意义：$2s$ 个一阶方程替代 $s$ 个二阶方程。
极限特例：$H=p^2/(2m)+V(q)$ 给出 $\dot q=p/m,\ \dot p=-V'(q)$。
应用：行星轨道数值积分、混沌分析。

**知识点3：相空间与刘维尔定理**
定义：以 $(\mathbf q,\mathbf p)$ 为坐标的 $2s$ 维空间为相空间；相流速度 $\mathbf v=(\dot{\mathbf q},\dot{\mathbf p})$。
推导：相空间“流体”的散度
$$\nabla\!\cdot\!\mathbf v=\sum_k\!\left(\frac{\partial\dot q_k}{\partial q_k}+\frac{\partial\dot p_k}{\partial p_k}\right)=\sum_k\!\left(\frac{\partial^2H}{\partial q_k\partial p_k}-\frac{\partial^2H}{\partial p_k\partial q_k}\right)=0.$$
$\because$ 散度为零，$\therefore$ 相空间体积流不可压缩，密度 $\rho$ 满足连续性方程
$$\frac{\partial\rho}{\partial t}+\nabla\!\cdot\!(\rho\mathbf v)=0\ \Rightarrow\ \frac{\mathrm d\rho}{\mathrm dt}=0.$$
意义：相体积守恒（刘维尔定理）。
极限特例：保守系相轨迹不收缩、不扩张。
应用：统计力学、长时间数值稳定性评估。

**知识点4：泊松括号**
定义：$\{f,g\}=\sum_k\!\left(\dfrac{\partial f}{\partial q_k}\dfrac{\partial g}{\partial p_k}-\dfrac{\partial f}{\partial p_k}\dfrac{\partial g}{\partial q_k}\right)$。
推导：对任意 $f(\mathbf q,\mathbf p,t)$，
$$\frac{\mathrm df}{\mathrm dt}=\sum_k\!\left(\frac{\partial f}{\partial q_k}\dot q_k+\frac{\partial f}{\partial p_k}\dot p_k\right)+\frac{\partial f}{\partial t}=\sum_k\!\left(\frac{\partial f}{\partial q_k}\frac{\partial H}{\partial p_k}-\frac{\partial f}{\partial p_k}\frac{\partial H}{\partial q_k}\right)+\frac{\partial f}{\partial t}.$$
$$\therefore \frac{\mathrm df}{\mathrm dt}=\{f,H\}+\frac{\partial f}{\partial t}.$$
基本括号：$\{q_i,q_j\}=0,\ \{p_i,p_j\}=0,\ \{q_i,p_j\}=\delta_{ij}$。
意义：守恒量满足 $\{f,H\}=0$。
极限特例：$\{f,H\}=0,\ \partial f/\partial t=0\Rightarrow f$ 为运动积分。
应用：识别守恒量、正则量子化对应。

**知识点5：正则变换**
定义：$(q,p)\to(Q,P)$ 保持正则方程形式；存在母函数。
推导：要求 $\sum p\,\mathrm dq-H\,\mathrm dt=\sum P\,\mathrm dQ-K\,\mathrm dt+\mathrm dF$。取 $F_1(q,Q,t)$：
$$\mathrm dF_1=\sum p\,\mathrm dq-\sum P\,\mathrm dQ+(K-H)\mathrm dt.$$
$$\therefore p_k=\frac{\partial F_1}{\partial q_k},\qquad P_k=-\frac{\partial F_1}{\partial Q_k},\qquad K=H+\frac{\partial F_1}{\partial t}.$$
意义：可选取使 $H$ 简化的坐标。
极限特例：恒等变换 $F_1=\sum q_kQ_k$。
应用：哈密顿–雅可比理论、摄动法。

**知识点6：哈密顿–雅可比方程**
定义：求母函数 $S(\mathbf q,\boldsymbol\alpha,t)$ 使新哈密顿 $K=0$。
推导：由 $F_1=S$，$K=H+\partial S/\partial t=0$，且 $p_k=\partial S/\partial q_k$：
$$H\!\left(\mathbf q,\frac{\partial S}{\partial\mathbf q},t\right)+\frac{\partial S}{\partial t}=0.$$
对保守系 $H=E$，分离变量 $S=W(\mathbf q)-\alpha_1 t$：
$$H\!\left(\mathbf q,\frac{\partial W}{\partial\mathbf q}\right)=\alpha_1.$$
意义：一阶偏微分方程，其解给出轨道。
极限特例：可分离变量时化为常微分方程。
应用：解析求解中心力、谐振子。

**知识点7：作用量变量与绝热不变量**
定义：对周期运动，$J_k=\oint p_k\,\mathrm dq_k$ 为作用量变量；角变量 $w_k$。
推导：对一维周期运动，$H=H(J)$，正则方程
$$\dot w=\frac{\partial H}{\partial J}\equiv\omega(J),\qquad \dot J=-\frac{\partial H}{\partial w}=0.$$
$\therefore J=\text{常量},\ w=\omega t+w_0$。
绝热不变量：若参数缓慢变化，$J$ 近似守恒。证明要点：一个周期内 $H$ 的相对变化 $\sim(\dot\lambda/\omega)/\lambda\ll1$。
意义：频率 $\omega=\partial H/\partial J$。
极限特例：谐振子 $J=E/\omega$。
应用：回旋加速器、等离子体磁矩守恒、天体力学。

**知识点8：对称性与守恒量（诺特定理特例）**
定义：若 $L$ 对坐标 $q_j$ 平移不变，则共轭动量守恒；对时间平移不变，则能量守恒；对旋转不变，则角动量守恒。
推导：空间平移 $q_j\to q_j+\epsilon$：
$$0=\frac{\partial L}{\partial q_j}=\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}\ \Rightarrow\ p_j=\text{常量}.$$
时间平移：$\dfrac{\mathrm dL}{\mathrm dt}=\sum\dot p_k\dot q_k+\sum p_k\ddot q_k$，用拉格朗日方程消去 $\dot p_k$ 得 $\dfrac{\mathrm dH}{\mathrm dt}=-\dfrac{\partial L}{\partial t}=0$。
意义：每一种连续对称性对应一个守恒量。
极限特例：自由粒子动量守恒；中心力角动量守恒。
应用：构造守恒格式、识别守恒律。

### 第6章 变分原理与能量法

**知识点1：泛函与变分**
定义：泛函 $I[y]=\int_{x_1}^{x_2}F(x,y,y')\,\mathrm dx$；变分 $\delta I$ 为函数 $y\to y+\epsilon\eta$ 下的一阶增量。
推导：令 $\phi(\epsilon)=I[y+\epsilon\eta]$，展开
$$\phi(\epsilon)=\int F(x,y+\epsilon\eta,y'+\epsilon\eta')\,\mathrm dx.$$
$$\phi'(0)=\int\!\left(\frac{\partial F}{\partial y}\eta+\frac{\partial F}{\partial y'}\eta'\right)\mathrm dx\equiv\delta I.$$
意义：变分是函数空间上的“微分”。
极限特例：$\delta I=0$ 对应驻值。
应用：能量原理、有限元。

**知识点2：欧拉–拉格朗日方程推导**
定义：使 $I$ 取驻值的 $y(x)$ 满足的微分方程。
推导：对 $\delta I$ 分部积分，
$$\int_{x_1}^{x_2}\frac{\partial F}{\partial y'}\eta'\,\mathrm dx=\left[\frac{\partial F}{\partial y'}\eta\right]_{x_1}^{x_2}-\int_{x_1}^{x_2}\frac{\mathrm d}{\mathrm dx}\!\left(\frac{\partial F}{\partial y'}\right)\eta\,\mathrm dx.$$
代入边界固定 $\eta(x_1)=\eta(x_2)=0$：
$$\delta I=\int_{x_1}^{x_2}\!\left[\frac{\partial F}{\partial y}-\frac{\mathrm d}{\mathrm dx}\frac{\partial F}{\partial y'}\right]\eta\,\mathrm dx=0.$$
$\because\eta$ 任意，$\therefore$
$$\boxed{\frac{\partial F}{\partial y}-\frac{\mathrm d}{\mathrm dx}\frac{\partial F}{\partial y'}=0.}$$
意义：变分问题化为微分方程。
极限特例：$F=\sqrt{1+y'^2}$ 给直线（最短路径）。
应用：测地线、最速降线、结构能量法。

**知识点3：哈密顿原理**
定义：真实运动使作用量 $S=\int_{t_1}^{t_2}L(\mathbf q,\dot{\mathbf q},t)\,\mathrm dt$ 取驻值，$\delta S=0$。
推导：对每个自由度 $q_k$ 独立变分，
$$\delta S=\int_{t_1}^{t_2}\sum_k\!\left(\frac{\partial L}{\partial q_k}\delta q_k+\frac{\partial L}{\partial\dot q_k}\delta\dot q_k\right)\mathrm dt.$$
分部积分（端点固定）：
$$\delta S=\int_{t_1}^{t_2}\sum_k\!\left[\frac{\partial L}{\partial q_k}-\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_k}\right]\delta q_k\,\mathrm dt=0.$$
$\because\delta q_k$ 独立，$\therefore$ 得拉格朗日方程。
意义：整个力学可由一个变分原理导出。
极限特例：非保守力须加虚功项 $\int(\mathbf F^{\text{nc}}\!\cdot\!\delta\mathbf r)\mathrm dt$。
应用：连续体动力学、场论。

**知识点4：最小势能原理**
定义：线弹性体总势能 $\Pi=U-W$（$U$ 应变能，$W$ 外力功），平衡时 $\delta\Pi=0$；对稳定平衡为极小。
推导：设位移场 $\mathbf u$，应变能
$$U=\int_\Omega\tfrac12\boldsymbol\sigma\!:\!\boldsymbol\varepsilon\,\mathrm dV.$$
外力功 $W=\int_\Omega\mathbf f\!\cdot\!\mathbf u\,\mathrm dV+\int_\Gamma\mathbf t\!\cdot\!\mathbf u\,\mathrm dS$。平衡要求 $\delta\Pi=0$：
$$\delta U=\int_\Omega\boldsymbol\sigma\!:\!\delta\boldsymbol\varepsilon\,\mathrm dV,\qquad \delta W=\int_\Omega\mathbf f\!\cdot\!\delta\mathbf u\,\mathrm dV+\int_\Gamma\mathbf t\!\cdot\!\delta\mathbf u\,\mathrm dS.$$
由 $\delta\Pi=0$ 与虚功原理等价，得到平衡方程与力的边界条件。
意义：位移法有限元的基础。
极限特例：线性小变形。
应用：Ritz 法求近似解。

**知识点5：余能原理**
定义：余能 $\Pi^*=U^*-W^*$，其中 $U^*=\int_\Omega\tfrac12\boldsymbol\sigma\!:\!\boldsymbol\varepsilon\,\mathrm dV$（线弹性下 $U^*=U$）。
推导：对静力容许应力场 $\boldsymbol\sigma$ 与 $W^*=\int_\Gamma\mathbf u\!\cdot\!\mathbf t\,\mathrm dS$，$\delta\Pi^*=0$ 给出协调方程与位移边界。
意义：与最小势能原理互为对偶；用于应力法。
极限特例：线弹性时 $U^*=U$。
应用：边界元、应力杂交元。

**知识点6：约束变分与自然边界条件**
定义：若端点 $y(x_i)$ 不固定，则变分产生自然边界条件。
推导：由 $\delta I$ 分部积分时边界项不消失：
$$\delta I=\int\!\left[\cdots\right]\eta\,\mathrm dx+\left[\frac{\partial F}{\partial y'}\eta\right]_{x_1}^{x_2}.$$
$\because$ 内部驻值要求 Euler 方程成立，$\therefore$ 剩下
$$\left.\frac{\partial F}{\partial y'}\right|_{x_1,x_2}=0.$$
意义：自然边界由泛函本身决定。
极限特例：自由端给力边界条件。
应用：梁的自由端、变分有限元边界处理。

**知识点7：瑞利–里兹法**
定义：取试探解 $y\approx\sum_{n=1}^{N}c_n\phi_n(x)$，$\phi_n$ 满足几何边界条件。
推导：代入泛函 $I[c_1,\dots,c_N]$，驻值条件
$$\frac{\partial I}{\partial c_n}=0,\quad n=1,\dots,N.$$
对线弹性问题 $I=\tfrac12\mathbf c^{\mathsf T}K\mathbf c-\mathbf c^{\mathsf T}\mathbf F$，
$$\therefore K\mathbf c=\mathbf F.$$
意义：变分近似解；$N$ 增大解收敛于真解。
极限特例：$N=1$ 单模近似。
应用：结构固有频率、临界载荷估算。

**知识点8：伽辽金法与加权残量**
定义：近似解余量 $R=\mathcal L[\tilde u]-f$，要求加权积分
$$\int_\Omega w_n R\,\mathrm dV=0.$$
推导：伽辽金取权 $w_n=\phi_n$（与试探函数同）：
$$\int_\Omega\phi_n\big(\mathcal L[\tilde u]-f\big)\mathrm dV=0.$$
对 $\mathcal L$ 为自伴算子，分部积分后与 Ritz 法一致。
意义：弱形式有限元的出发点。
极限特例：$R=0$ 精确解。
应用：FEM、谱方法。

### 第7章 应力、应变张量与平衡方程

**知识点1：应力张量与柯西公式**
定义：受力体内某点法向 $\mathbf n$ 的截面上的面力矢量 $\mathbf t^{(\mathbf n)}$；应力张量 $\boldsymbol\sigma$ 使
$$\mathbf t^{(\mathbf n)}=\boldsymbol\sigma\!\cdot\!\mathbf n,\qquad t_i^{(\mathbf n)}=\sigma_{ij}n_j.$$
推导（四面体）：取小四面体，斜面法向 $\mathbf n$，面积 $\Delta A$，三坐标面面积 $\Delta A_i=n_i\Delta A$。由动量平衡（体积力为高阶小量）：
$$\mathbf t^{(\mathbf n)}\Delta A-\sum_i\mathbf t^{(-i)}\Delta A_i+\cdots=0,$$
$\because\mathbf t^{(-i)}=-\mathbf t^{(i)}$，$\therefore t_j^{(\mathbf n)}=n_i\sigma_{ij}$。
意义：一点应力状态由二阶张量完全确定。
极限特例：$\mathbf n=\hat e_i$ 得 $\sigma_{ij}$ 分量。
应用：强度校核、有限元应力恢复。

**知识点2：主应力与莫尔圆**
定义：主应力为 $\boldsymbol\sigma$ 的特征值；主方向满足 $\boldsymbol\sigma\!\cdot\!\mathbf n=\sigma\mathbf n$。
推导：$(\boldsymbol\sigma-\sigma\mathbf I)\mathbf n=0$ 有非零解
$$\det(\boldsymbol\sigma-\sigma\mathbf I)=0,$$
展开得特征方程
$$\sigma^3-I_1\sigma^2+I_2\sigma-I_3=0,$$
$$I_1=\operatorname{tr}\boldsymbol\sigma,\quad I_2=\tfrac12\big[(\operatorname{tr}\boldsymbol\sigma)^2-\operatorname{tr}\boldsymbol\sigma^2\big],\quad I_3=\det\boldsymbol\sigma.$$
二维时
$$\sigma_{1,2}=\frac{\sigma_x+\sigma_y}{2}\pm\sqrt{\left(\frac{\sigma_x-\sigma_y}{2}\right)^2+\tau_{xy}^2}.$$
莫尔圆圆心 $\dfrac{\sigma_x+\sigma_y}{2}$，半径 $R=\sqrt{\left(\frac{\sigma_x-\sigma_y}{2}\right)^2+\tau_{xy}^2}$。
意义：不变量 $I_1,I_2,I_3$ 与坐标无关。
极限特例：纯剪 $\sigma_1=-\sigma_3=|\tau|$。
应用：实验应力分析、屈服判据。

**知识点3：平衡微分方程**
定义：应力场须满足
$$\sigma_{ij,j}+f_i=0.$$
推导：对微元体 $\mathrm dV$，面力在 $i$ 方向合力为 $\sigma_{ij,j}\mathrm dV$，体积力 $f_i\mathrm dV$，动量平衡（静力）：
$$\sigma_{ij,j}\mathrm dV+f_i\mathrm dV=0\ \Rightarrow\ \sigma_{ij,j}+f_i=0.$$
动态时
$$\sigma_{ij,j}+f_i=\rho\ddot u_i.$$
展开：
$$\frac{\partial\sigma_x}{\partial x}+\frac{\partial\tau_{xy}}{\partial y}+\frac{\partial\tau_{xz}}{\partial z}+f_x=0,\ \dots$$
意义：静力平衡的局部形式。
极限特例：无体力时为 $\sigma_{ij,j}=0$。
应用：有限元/差分求解基础。

**知识点4：应变张量与几何方程**
定义：小变形下应变张量
$$\varepsilon_{ij}=\tfrac12(u_{i,j}+u_{j,i}).$$
推导：两邻近点 $\mathbf r$ 与 $\mathbf r+\mathrm d\mathbf r$，变形后位移差 $\mathrm du_i=u_{i,j}\mathrm dx_j$。线元平方变化
$$\mathrm ds'^2-\mathrm ds^2=2\varepsilon_{ij}\mathrm dx_i\mathrm dx_j.$$
由
$$\mathrm ds'^2=(u_{i,j}\mathrm dx_j+\mathrm dx_i)(u_{i,k}\mathrm dx_k+\mathrm dx_i),$$
略去高阶项得
$$\therefore \varepsilon_{ij}=\tfrac12(u_{i,j}+u_{j,i}).$$
意义：对称张量，描述长度与角度变化。
极限特例：刚体转动应变为零（反对称部分为转动）。
应用：应变花测量、DIC 计算。

**知识点5：应变协调方程**
定义：由位移单值连续要求的应变间约束。
推导：对几何方程求二阶导，交换求导次序：
$$\varepsilon_{ij,kl}+\varepsilon_{kl,ij}-\varepsilon_{ik,jl}-\varepsilon_{jl,ik}=0.$$
三维共 6 个独立协调方程。二维：
$$\frac{\partial^2\varepsilon_x}{\partial y^2}+\frac{\partial^2\varepsilon_y}{\partial x^2}=\frac{\partial^2\gamma_{xy}}{\partial x\partial y}.$$
意义：任意给定应变场未必对应单值位移。
极限特例：二维问题仅 1 个协调方程。
应用：应力函数法、光弹性相容性检验。

**知识点6：应变花与主应变**
定义：$45^\circ$ 三向应变花测 $\varepsilon_0,\varepsilon_{45},\varepsilon_{90}$。
推导：由应变变换
$$\varepsilon_\theta=\varepsilon_x\cos^2\theta+\varepsilon_y\sin^2\theta+\gamma_{xy}\sin\theta\cos\theta.$$
代入 $\theta=0^\circ,45^\circ,90^\circ$：
$$\varepsilon_0=\varepsilon_x,\qquad\varepsilon_{90}=\varepsilon_y,$$
$$\varepsilon_{45}=\frac{\varepsilon_x+\varepsilon_y}{2}+\frac{\gamma_{xy}}{2}\ \Rightarrow\ \gamma_{xy}=2\varepsilon_{45}-\varepsilon_0-\varepsilon_{90}.$$
主应变
$$\varepsilon_{1,2}=\frac{\varepsilon_0+\varepsilon_{90}}{2}\pm\frac{1}{\sqrt2}\sqrt{(\varepsilon_0-\varepsilon_{45})^2+(\varepsilon_{45}-\varepsilon_{90})^2}.$$
意义：由表面应变反推主应变。
极限特例：单向应力状态 $\varepsilon_{45}=(\varepsilon_0+\varepsilon_{90})/2$。
应用：应变片实验、残余应力测量。

**知识点7：柱坐标与轴对称问题**
定义：轴对称时 $\partial/\partial\theta=0$，$u_\theta=0$。
推导：柱坐标几何方程
$$\varepsilon_r=\frac{\partial u_r}{\partial r},\quad\varepsilon_\theta=\frac{u_r}{r},\quad\varepsilon_z=\frac{\partial u_z}{\partial z},\quad\gamma_{rz}=\frac{\partial u_r}{\partial z}+\frac{\partial u_z}{\partial r}.$$
平衡方程（无体力）
$$\frac{\partial\sigma_r}{\partial r}+\frac{\partial\tau_{rz}}{\partial z}+\frac{\sigma_r-\sigma_\theta}{r}=0.$$
意义：厚壁筒、圆盘问题。
极限特例：平面应变时 $\varepsilon_z=0$。
应用：压力容器、旋转轮盘。

**知识点8：圣维南原理**
定义：作用于弹性体一小块面积上的力系，可用其静力等效的另一力系替代，在远离作用区（距离 $\gg$ 作用区尺寸）应力分布近似相同。
推导（量级估计）：局部力系自平衡（主矢主矩为零），其应力场随距离 $r$ 衰减，量级 $\sim1/r^2$（或 $1/r^3$）；对有限体，影响按 $\exp(-kr/d)$ 衰减，$d$ 为作用区尺寸。
意义：边界力分布的细节只影响局部。
极限特例：严格证明对半无限体成立。
应用：有限元边界条件简化、梁端加载。

### 第8章 本构关系与强度理论

**知识点1：广义胡克定律**
定义：线弹性小变形下
$$\varepsilon_{ij}=\frac{1+\nu}{E}\sigma_{ij}-\frac{\nu}{E}\sigma_{kk}\delta_{ij},$$
或
$$\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2G\varepsilon_{ij}.$$
推导：各向同性弹性张量
$$C_{ijkl}=\lambda\delta_{ij}\delta_{kl}+G(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}),\qquad \sigma_{ij}=C_{ijkl}\varepsilon_{kl}.$$
对 $\sigma_{kk}$：
$$\sigma_{kk}=(3\lambda+2G)\varepsilon_{kk}\ \Rightarrow\ \lambda=\frac{E\nu}{(1+\nu)(1-2\nu)},\ G=\frac{E}{2(1+\nu)}.$$
逆解得柔度形式。
意义：两个独立弹性常数（$E,\nu$）即可描述各向同性材料。
极限特例：$\nu=0.5$ 不可压，$\lambda,G$ 中 $K\to\infty$。
应用：结构变形计算。

**知识点2：弹性常数之间的关系**
定义：$E$（杨氏模量）、$\nu$（泊松比）、$G$（剪切模量）、$K$（体积模量）、$\lambda$（拉梅常数）。
推导：单轴拉伸 $\sigma_x=\sigma$，其余为零，则
$$\varepsilon_x=\frac{\sigma}{E},\quad\varepsilon_y=\varepsilon_z=-\frac{\nu\sigma}{E}.$$
体积应变
$$\frac{\Delta V}{V}=\varepsilon_x+\varepsilon_y+\varepsilon_z=\frac{1-2\nu}{E}\sigma\ \Rightarrow\ K=\frac{E}{3(1-2\nu)}.$$
纯剪可得
$$G=\frac{E}{2(1+\nu)}.$$
意义：常数两两独立，其余可导出。
极限特例：$\nu=1/3$ 时 $K=E$（近似金属）。
应用：材料参数测试。

**知识点3：热弹性本构**
定义：温度变化 $\Delta T=T-T_0$ 引起热应变 $\alpha\Delta T\delta_{ij}$。
推导：总应变 = 弹性应变 + 热应变：
$$\varepsilon_{ij}=\frac{1+\nu}{E}\sigma_{ij}-\frac{\nu}{E}\sigma_{kk}\delta_{ij}+\alpha\Delta T\,\delta_{ij}.$$
逆解：
$$\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2G\varepsilon_{ij}-(3\lambda+2G)\alpha\Delta T\,\delta_{ij}.$$
$\because 3\lambda+2G=3K$，$\therefore$
$$\sigma_{ij}=2G\varepsilon_{ij}+\lambda\varepsilon_{kk}\delta_{ij}-3K\alpha\Delta T\delta_{ij}.$$
意义：自由膨胀不产生应力。
极限特例：全约束一维杆 $\sigma=-E\alpha\Delta T$。
应用：热应力、热防护结构。

**知识点4：正交各向异性本构与工程常数**
定义：三个材料主轴方向，应力–应变关系
$$\begin{pmatrix}\varepsilon_1\\\varepsilon_2\\\varepsilon_3\\\gamma_{23}\\\gamma_{31}\\\gamma_{12}\end{pmatrix}
=\begin{pmatrix}
1/E_1&-\nu_{21}/E_2&-\nu_{31}/E_3&0&0&0\\
-\nu_{12}/E_1&1/E_2&-\nu_{32}/E_3&0&0&0\\
-\nu_{13}/E_1&-\nu_{23}/E_2&1/E_3&0&0&0\\
0&0&0&1/G_{23}&0&0\\
0&0&0&0&1/G_{31}&0\\
0&0&0&0&0&1/G_{12}
\end{pmatrix}
\begin{pmatrix}\sigma_1\\\sigma_2\\\sigma_3\\\tau_{23}\\\tau_{31}\\\tau_{12}\end{pmatrix}.$$
推导：由柔度矩阵的对称性 $\nu_{ij}/E_i=\nu_{ji}/E_j$（能量互等），9 个独立常数（3 个 $E$、3 个 $G$、3 个 $\nu$）。
意义：复合材料层合板分析基础。
极限特例：横观各向同性（纤维单向板）减为 5 个常数。
应用：碳纤维层合结构。

**知识点5：屈服准则（Tresca 与 von Mises）**
定义：Tresca 最大剪应力准则 $\tau_{\max}=\dfrac{\sigma_1-\sigma_3}{2}=k$；von Mises 准则 $J_2=k^2$。
推导：偏应力 $s_{ij}=\sigma_{ij}-\tfrac13\sigma_{kk}\delta_{ij}$，第二不变量
$$J_2=\tfrac12 s_{ij}s_{ij}=\tfrac16\big[(\sigma_1-\sigma_2)^2+(\sigma_2-\sigma_3)^2+(\sigma_3-\sigma_1)^2\big].$$
单轴拉伸屈服 $\sigma_1=\sigma_y$，其余为零：
$$J_2=\tfrac13\sigma_y^2\ \Rightarrow\ \sigma_{\text{vM}}=\sqrt{3J_2}=\sqrt{\tfrac12\big[(\sigma_1-\sigma_2)^2+\cdots\big]}=\sigma_y.$$
Tresca 单轴时 $k=\sigma_y/2$。
意义：两准则在 $\pi$ 平面分别为外接圆与内切六边形。
极限特例：静水压力不影响屈服（$J_2$ 仅含偏量）。
应用：金属强度校核、有限元塑性。

**知识点6：黏弹性本构（Maxwell 与 Kelvin–Voigt）**
定义：Maxwell 模型（弹簧与黏壶串联）：
$$\dot\varepsilon=\frac{\dot\sigma}{E}+\frac{\sigma}{\eta}.$$
Kelvin–Voigt 模型（并联）：
$$\sigma=E\varepsilon+\eta\dot\varepsilon.$$
推导（Maxwell 松弛）：恒应变 $\dot\varepsilon=0$，
$$\frac{\dot\sigma}{E}+\frac{\sigma}{\eta}=0\ \Rightarrow\ \frac{\mathrm d\sigma}{\sigma}=-\frac{E}{\eta}\mathrm dt,$$
$$\therefore \sigma(t)=\sigma_0 e^{-t/\tau},\qquad \tau=\frac{\eta}{E}.$$
（Kelvin–Voigt 蠕变）：恒应力 $\sigma_0$，
$$\dot\varepsilon+\frac{E}{\eta}\varepsilon=\frac{\sigma_0}{\eta}\ \Rightarrow\ \varepsilon(t)=\frac{\sigma_0}{E}\big(1-e^{-t/\tau}\big).$$
意义：分别描述应力松弛与蠕变。
极限特例：$t\to\infty$ Maxwell 应力趋于零；K–V 趋于稳态。
应用：聚合物、生物组织、阻尼材料。

**知识点7：经典强度理论与失效判据**
定义：四种经典强度理论。
推导：
（1）最大拉应力：$\sigma_1\le[\sigma]$。
（2）最大伸长线应变：$\varepsilon_1\le[\varepsilon]$，即 $\sigma_1-\nu(\sigma_2+\sigma_3)\le[\sigma]$。
（3）最大剪应力：$\sigma_1-\sigma_3\le[\sigma]$。
（4）畸变能（von Mises）：$\sigma_{\text{vM}}\le[\sigma]$。
其中 $[\sigma]=\sigma_u/n$，$n$ 为安全系数。
意义：脆性材料用（1）（2），塑性材料用（3）（4）。
极限特例：三向等拉时（1）判危险，剪应力理论失效。
应用：轴、压力容器设计。

**知识点8：本构参数实验标定**
定义：由单轴拉伸、纯剪、体积压缩试验确定 $E,\nu,G,K$。
推导：单轴 $\sigma$–$\varepsilon$ 曲线初始斜率 $E=\dfrac{\Delta\sigma}{\Delta\varepsilon}$；横向应变比 $\nu=-\dfrac{\Delta\varepsilon_y}{\Delta\varepsilon_x}$。扭转试验 $\tau=G\gamma$，$G=\dfrac{T L}{J\phi}$（$J$ 极惯性矩，$\phi$ 扭转角）。体积模量 $K=-\dfrac{\Delta p}{\Delta V/V}$。
意义：参数一致性检验 $G=\dfrac{E}{2(1+\nu)}$。
极限特例：各向同性材料三法测得的 $E,G,\nu$ 应自洽。
应用：材料数据库、仿真输入。

### 第9章 弹性力学基本方程与能量原理

**知识点1：弹性力学三类基本方程与边界条件**
定义：15 个未知量（6 应力 + 6 应变 + 3 位移）对应 15 个方程：3 平衡 + 6 几何 + 6 本构。
边界：
$$\text{位移边界}\ \Gamma_u:\ u_i=\bar u_i;\qquad \text{力边界}\ \Gamma_t:\ \sigma_{ij}n_j=\bar t_i.$$
推导：未知数计数
$$\underbrace{3}_{\text{平衡}}+\underbrace{6}_{\text{几何(含协调)}}+\underbrace{6}_{\text{本构}}=15.$$
加上 $\Gamma_u\cup\Gamma_t=\Gamma,\ \Gamma_u\cap\Gamma_t=\varnothing$，问题适定。
意义：弹性力学边值问题。
极限特例：纯位移法或纯应力法须消元。
应用：解析与数值求解。

**知识点2：位移法（Navier 方程）**
定义：以位移为基本未知量
$$G\nabla^2 u_i+(\lambda+G)u_{j,ji}+f_i=0.$$
推导：本构 $\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2G\varepsilon_{ij}$ 代入平衡 $\sigma_{ij,j}+f_i=0$：
$$\lambda\varepsilon_{kk,i}+2G\varepsilon_{ij,j}+f_i=0.$$
由 $\varepsilon_{kk,i}=u_{k,ki}$，$\varepsilon_{ij,j}=\tfrac12(u_{i,jj}+u_{j,ij})$，得
$$\lambda u_{k,ki}+G(u_{i,jj}+u_{j,ij})+f_i=0,$$
即
$$G u_{i,jj}+(\lambda+G)u_{j,ji}+f_i=0.$$
矢量形式（$\mathbf u$）：
$$G\nabla^2\mathbf u+(\lambda+G)\nabla(\nabla\!\cdot\!\mathbf u)+\mathbf f=0.$$
意义：三个二阶 PDE，位移边界直接可用。
极限特例：无体力、$\nabla\!\cdot\!\mathbf u=0$（不可压）时 $\nabla^2\mathbf u=0$。
应用：FEM 位移元。

**知识点3：应力法（Beltrami–Michell 方程）**
定义：以应力为基本未知量，需满足平衡 + 协调。
推导：由协调方程 $\varepsilon_{ij,kl}+\varepsilon_{kl,ij}-\varepsilon_{ik,jl}-\varepsilon_{jl,ik}=0$ 代入本构并利用平衡消去体力，得
$$\nabla^2\sigma_{ij}+\frac{1}{1+\nu}\Theta_{,ij}=-\frac{\nu}{1-\nu}f_{k,k}\delta_{ij}-(f_{i,j}+f_{j,i}),$$
其中 $\Theta=\sigma_{kk}$。无体力时
$$\nabla^2\sigma_{ij}+\frac{1}{1+\nu}\Theta_{,ij}=0,\qquad \nabla^2\Theta=0.$$
意义：应力函数法的基础。
极限特例：$\Theta$ 为调和函数。
应用：平面问题 Airy 函数。

**知识点4：应变能密度**
定义：$\mathcal U_0=\displaystyle\int_0^{\varepsilon_{ij}}\sigma_{ij}\,\mathrm d\varepsilon_{ij}$；线弹性
$$\mathcal U_0=\tfrac12\sigma_{ij}\varepsilon_{ij}=\tfrac12 C_{ijkl}\varepsilon_{ij}\varepsilon_{kl}.$$
推导：一维弹簧 $U=\int_0^\varepsilon E\varepsilon\,\mathrm d\varepsilon=\tfrac12E\varepsilon^2$。三维按分量叠加（线弹性无耦合路径依赖）：
$$\mathcal U_0=\int_0^{\varepsilon}\sum\sigma_{ij}\mathrm d\varepsilon_{ij}=\tfrac12\sigma_{ij}\varepsilon_{ij}.$$
分解：
$$\mathcal U_0=\underbrace{\tfrac12 K\varepsilon_{kk}^2}_{\text{体积}}+\underbrace{G\,\varepsilon_{ij}'\varepsilon_{ij}'}_{\text{畸变}}.$$
意义：能量是正定二次型（稳定性）。
极限特例：$K\to\infty$ 时体积变为零。
应用：能量法、断裂力学 $J$ 积分。

**知识点5：卡氏定理**
定义：线弹性体位移 $u_i=\partial U/\partial P_i$，其中 $U$ 为余能（线弹性下等于应变能）。
推导：总应变能 $U=U(P_1,\dots,P_n)$，外力功 $W=\sum\tfrac12P_i u_i$。由能量守恒 $U=W$。虚位移 $\delta u_i$，
$$\delta U=\sum_i\frac{\partial U}{\partial P_i}\delta P_i.$$
但更直接：对力 $P_i$ 的广义位移
$$u_i=\frac{\partial U}{\partial P_i}.$$
证明：令 $P_i$ 增 $\delta P_i$，做功增量 $\tfrac12\delta P_i\delta u_i$ 为二阶小量，故 $\mathrm dU=u_i\,\delta P_i$。
意义：求任意点位移的简便方法。
极限特例：须保证 $P_i$ 处有对应位移。
应用：桁架、梁挠度计算。

**知识点6：互等定理（Betti 与 Maxwell）**
定义：Betti 互等：第一组力在第二组位移上做功等于第二组力在第一组位移上做功。
推导：线弹性应变能
$$U=\tfrac12\sum_i P_i u_i.$$
考虑两种加载状态 (1) 与 (2)，先加 (1) 后加 (2)，总功
$$W=\tfrac12\sum P_i^{(1)}u_i^{(1)}+\sum P_i^{(1)}u_i^{(2)}+\tfrac12\sum P_i^{(2)}u_i^{(2)}.$$
反向加载得
$$\sum P_i^{(1)}u_i^{(2)}=\sum P_i^{(2)}u_i^{(1)}.$$
Maxwell 互等（单位力）：$\delta_{ij}=\delta_{ji}$。
意义：柔度矩阵对称。
极限特例：多力系统适用。
应用：影响线、有限元刚度矩阵对称性。

**知识点7：唯一性定理**
定义：给定体力、面力、位移边界，线弹性解唯一。
推导（反证）：设两组解 $\mathbf u^{(1)},\mathbf u^{(2)}$，差 $\Delta\mathbf u=\mathbf u^{(1)}-\mathbf u^{(2)}$ 对应零体力、零面力、零位移边界。其应变能
$$U=\int_\Omega\tfrac12\Delta\sigma_{ij}\Delta\varepsilon_{ij}\,\mathrm dV=\int_\Gamma\Delta t_i\Delta u_i\,\mathrm dS-\int_\Omega(\Delta\sigma_{ij,j})\Delta u_i\,\mathrm dV=0.$$
$\because U$ 正定，$\therefore\Delta\varepsilon_{ij}=0$，刚体位移由边界零位移消除，$\therefore\Delta\mathbf u=0$。
意义：解的存在性之外，唯一性保证数值收敛到正确解。
极限特例：理想塑性、不稳定材料可能失去唯一性。
应用：保证 FEM 解可信。

**知识点8：平面问题的 Airy 应力函数**
定义：无体力平面问题，取 $\phi(x,y)$ 使
$$\sigma_x=\frac{\partial^2\phi}{\partial y^2},\quad\sigma_y=\frac{\partial^2\phi}{\partial x^2},\quad\tau_{xy}=-\frac{\partial^2\phi}{\partial x\partial y}.$$
推导：平衡方程自动满足：
$$\frac{\partial\sigma_x}{\partial x}+\frac{\partial\tau_{xy}}{\partial y}=\frac{\partial^3\phi}{\partial x\partial y^2}-\frac{\partial^3\phi}{\partial x\partial y^2}=0.$$
代入协调方程
$$\left(\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}\right)(\sigma_x+\sigma_y)=0,$$
得双调和方程
$$\nabla^4\phi=\frac{\partial^4\phi}{\partial x^4}+2\frac{\partial^4\phi}{\partial x^2\partial y^2}+\frac{\partial^4\phi}{\partial y^4}=0.$$
意义：平面弹性问题化为求双调和函数。
极限特例：多项式解对应均匀应力场。
应用：孔边应力集中、裂纹尖端场。

### 第10章 梁的弯曲、剪切与扭转

**知识点1：梁的弯曲正应力**
定义：纯弯曲时横截面正应力
$$\sigma=\frac{My}{I},\qquad \sigma_{\max}=\frac{M c}{I}.$$
假设：平截面、小变形、材料线弹性、纵向纤维无挤压。
推导：曲率 $\kappa=1/\rho$，应变沿高度线性
$$\varepsilon=\frac{y}{\rho}.$$
由胡克定律 $\sigma=E\varepsilon=E y/\rho$。由轴向力为零：
$$\int_A\sigma\,\mathrm dA=\frac{E}{\rho}\int_A y\,\mathrm dA=0\ \Rightarrow\ \bar y=0\ (\text{中性轴过形心}).$$
由弯矩定义：
$$M=\int_A\sigma y\,\mathrm dA=\frac{E}{\rho}\int_A y^2\mathrm dA=\frac{EI}{\rho}\ \Rightarrow\ \frac{1}{\rho}=\frac{M}{EI}.$$
代入 $\sigma=E y/\rho$ 得 $\sigma=My/I$。
意义：$I$ 为截面二次矩，中性轴在形心。
极限特例：$y=0$ 应力为零；$y=c$ 最大。
应用：梁强度校核。

**知识点2：挠曲线微分方程**
定义：挠度 $w(x)$ 满足
$$EI\frac{\mathrm d^2w}{\mathrm dx^2}=M(x).$$
推导：由曲率 $\dfrac{1}{\rho}=\dfrac{\mathrm d^2w/\mathrm dx^2}{[1+(w')^2]^{3/2}}\approx w''$（小变形 $w'\ll1$）。结合 $1/\rho=M/EI$：
$$w''=\frac{M}{EI}.$$
又 $M=EIw''$，剪力 $V=\dfrac{\mathrm dM}{\mathrm dx}=EIw'''$，分布载荷 $q=\dfrac{\mathrm dV}{\mathrm dx}=-EIw''''$：
$$EI\frac{\mathrm d^4w}{\mathrm dx^4}=q(x).$$
意义：四阶 ODE，四个边界条件。
极限特例：纯弯 $q=0\Rightarrow w$ 为二次多项式。
应用：梁挠度求解。

**知识点3：弯曲切应力与剪流**
定义：矩形截面切应力
$$\tau=\frac{VQ}{Ib},\qquad Q=\int_{y}^{c}\tilde y\,\mathrm dA.$$
推导：取长 $\mathrm dx$ 微段，考察 $y$ 以上部分，正应力差合力
$$\Delta F=\int_A\frac{\Delta M\,\tilde y}{I}\mathrm dA=\frac{V\,\mathrm dx}{I}Q.$$
平衡 $\tau b\,\mathrm dx=\Delta F$：
$$\therefore \tau=\frac{VQ}{Ib}.$$
剪流 $q_s=\tau t=\dfrac{VQ}{I}$。
意义：$Q$ 为部分面积对中性轴一次矩。
极限特例：矩形截面 $\tau_{\max}=\dfrac32\dfrac{V}{A}$（中性轴处）。
应用：工字梁、组合梁连接件设计。

**知识点4：弯曲中心**
定义：横向载荷通过弯曲中心时梁只弯不扭。
推导：开口薄壁截面剪流对形心力矩
$$M_s=\oint q_s\,\rho\,\mathrm ds=\oint\frac{VQ}{I}\rho\,\mathrm ds.$$
令对形心合成的等效横向力作用点使 $M_s=0$，即弯曲中心坐标
$$y_c=\frac{\int q_s\,z\,\mathrm ds}{V},\qquad z_c=-\frac{\int q_s\,y\,\mathrm ds}{V}.$$
意义：开口薄壁梁的弯曲中心一般不在形心。
极限特例：双轴对称截面弯曲中心与形心重合。
应用：机翼盒段、薄壁梁加载点选择。

**知识点5：自由扭转（圣维南扭转）**
定义：圆轴扭转
$$\tau=\frac{T\rho}{J},\qquad \phi=\frac{TL}{GJ},\qquad \tau_{\max}=\frac{Tc}{J}.$$
推导：平面假设，扭转角沿轴线性，剪应变
$$\gamma=\rho\frac{\mathrm d\phi}{\mathrm dx}=\frac{\rho\phi}{L}.$$
$\tau=G\gamma=G\rho\,\mathrm d\phi/\mathrm dx$。扭矩
$$T=\int_A\tau\rho\,\mathrm dA=G\frac{\mathrm d\phi}{\mathrm dx}\int_A\rho^2\mathrm dA=GJ\frac{\mathrm d\phi}{\mathrm dx}.$$
$\therefore \dfrac{\mathrm d\phi}{\mathrm dx}=\dfrac{T}{GJ}$，$\tau=\dfrac{T\rho}{J}$。
意义：$J$ 极惯性矩；非圆截面须用普朗特应力函数。
极限特例：圆截面 $\tau_{\max}$ 在外缘。
应用：传动轴、扭力杆。

**知识点6：开口与闭口薄壁杆的扭转**
定义：闭口薄壁 $J_t=\dfrac{4A_m^2}{\oint\dfrac{\mathrm ds}{t}}$，剪流 $q=\dfrac{T}{2A_m}$；开口薄壁 $J_t=\dfrac13\int t^3\mathrm ds$。
推导（闭口 Bredt）：由剪流 $q=\tau t$ 沿壁为常数。扭矩
$$T=\oint q\,\rho\,\mathrm ds=q\oint\rho\,\mathrm ds=2q A_m\ \Rightarrow\ q=\frac{T}{2A_m}.$$
扭转角能法：$U=\int\dfrac{\tau^2}{2G}\mathrm dV=\dfrac{T^2L}{8GA_m^2}\oint\dfrac{\mathrm ds}{t}$，又 $U=\dfrac{T^2L}{2GJ_t}$，$\therefore J_t=\dfrac{4A_m^2}{\oint\mathrm ds/t}$。
开口截面用矩形条带近似，$J_t\approx\frac13\int t^3\mathrm ds$。
意义：闭口截面抗扭远高于开口。
极限特例：同一材料闭口 $J_t$ 比开口大 $O((R/t)^2)$ 倍。
应用：机翼盒段（多闭室）抗扭。

**知识点7：组合变形与叠加**
定义：拉（压）+ 弯 + 扭同时作用。
推导：线弹性小变形可叠加。危险点正应力
$$\sigma=\frac{N}{A}+\frac{M_y z}{I_y}+\frac{M_z y}{I_z},$$
剪应力
$$\tau=\frac{T\rho}{J}+\frac{V Q}{I b}.$$
危险点用强度理论：
$$\sigma_{\text{vM}}=\sqrt{\sigma^2+3\tau^2}\le[\sigma].$$
意义：多维载荷下的强度评估。
极限特例：单向拉弯 $\tau=0$。
应用：轴类、框架节点。

**知识点8：超静定梁与能量法**
定义：多余约束使未知反力数大于平衡方程数。
推导：用单位载荷法（莫尔积分）
$$\Delta=\int\frac{M\bar M}{EI}\mathrm dx+\int\frac{N\bar N}{EA}\mathrm dx+\int\frac{kV\bar V}{GA}\mathrm dx.$$
协调条件（多余约束处位移为零）给出补充方程。
意义：力法/位移法求解超静定。
极限特例：简支梁为静定；固支梁一次超静定。
应用：连续梁、框架。

### 第11章 板壳理论与薄壁结构

**知识点1：薄板 Kirchhoff 假设**
定义：小挠度薄板（厚度 $h$ 远小于面内尺寸 $a$）：
（1）中面法线变形后仍为直线且垂直于中面（直法线）；
（2）厚度方向正应力 $\sigma_z\approx0$；
（3）中面无面内位移（薄膜忽略）。
推导：由（1）位移场
$$u=-z\frac{\partial w}{\partial x},\qquad v=-z\frac{\partial w}{\partial y},\qquad w=w(x,y).$$
应变
$$\varepsilon_x=-z\frac{\partial^2w}{\partial x^2},\quad\varepsilon_y=-z\frac{\partial^2w}{\partial y^2},\quad\gamma_{xy}=-2z\frac{\partial^2w}{\partial x\partial y}.$$
意义：全部应变由挠度 $w$ 的曲率表示。
极限特例：$w$ 常数时无应变（刚体平动）。
应用：薄板弯曲分析。

**知识点2：薄板平衡微分方程**
定义：薄板弯曲方程
$$D\nabla^4w=q(x,y),\qquad D=\frac{Eh^3}{12(1-\nu^2)}.$$
推导：弯矩由应力沿厚度积分：
$$M_x=\int_{-h/2}^{h/2}z\sigma_x\,\mathrm dz=-D\!\left(\frac{\partial^2w}{\partial x^2}+\nu\frac{\partial^2w}{\partial y^2}\right),$$
$$M_y=-D\!\left(\frac{\partial^2w}{\partial y^2}+\nu\frac{\partial^2w}{\partial x^2}\right),\qquad M_{xy}=-D(1-\nu)\frac{\partial^2w}{\partial x\partial y}.$$
剪力 $Q_x=\partial M_x/\partial x+\partial M_{xy}/\partial y$。代入 $z$ 方向平衡
$$\frac{\partial Q_x}{\partial x}+\frac{\partial Q_y}{\partial y}+q=0,$$
消去剪力得 $D\nabla^4w=q$。
意义：板弯曲为四阶双调和型方程。
极限特例：$q=0$ 为双调和方程。
应用：板挠度、屈曲。

**知识点3：圆板轴对称弯曲**
定义：轴对称时 $w=w(r)$，方程
$$D\!\left(\frac{\mathrm d^2}{\mathrm dr^2}+\frac1r\frac{\mathrm d}{\mathrm dr}\right)^2w=q.$$
推导：拉普拉斯算子 $\nabla^2=\dfrac{\mathrm d^2}{\mathrm dr^2}+\dfrac1r\dfrac{\mathrm d}{\mathrm dr}=\dfrac1r\dfrac{\mathrm d}{\mathrm dr}\!\left(r\dfrac{\mathrm d}{\mathrm dr}\right)$。积分四次得通解
$$w(r)=\frac{qr^4}{64D}+C_1r^2\ln r+C_2r^2+C_3\ln r+C_4.$$
对实心板省略奇异项 $\ln r,\ r^2\ln r$：
$$w=\frac{qr^4}{64D}+C_2r^2+C_4.$$
边界：周边固支 $w(a)=w'(a)=0$，周边简支 $w(a)=0,\ M_r(a)=0$。
意义：中心挠度可解析。
极限特例：周边固支圆板中心挠度 $w_{\max}=\dfrac{qa^4}{64D}$。
应用：压力容器封头、膜片。

**知识点4：板边界条件**
定义：常见边界：
简支 $w=0,\ M_n=0$；固支 $w=0,\ \partial w/\partial n=0$；自由 $M_n=0,\ V_n=0$（Kirchhoff 有效剪力）。
推导（有效剪力）：自由边须把扭矩 $M_{nt}$ 等效为分布剪力
$$V_n=Q_n+\frac{\partial M_{nt}}{\partial s}.$$
由扭矩沿边界等效：
$$\therefore V_n=\left(\frac{\partial M_n}{\partial n}+\frac{\partial M_{nt}}{\partial s}\right)+\frac{\partial M_{nt}}{\partial s}\ \text{的恰当组合}.$$
意义：板弯曲四个边界条件（每边两个）。
极限特例：简支板可用双正弦级数。
应用：板壳有限元边界施加。

**知识点5：加筋板与有效宽度**
定义：蒙皮受纵向压力时，与桁条共同承力，蒙皮有效宽度 $b_{\text{eff}}$。
推导：等刚度等效
$$b_{\text{eff}} h E=\int_{-b/2}^{b/2}E\,\mathrm dx\ \Rightarrow\ b_{\text{eff}}=\int_{-b/2}^{b/2}\frac{\sigma(x)}{\sigma_{\max}}\mathrm dx.$$
由板屈曲应力分布近似 $\dfrac{\sigma}{\sigma_0}=\dfrac{1}{1+(x/\lambda)^2}$，得
$$b_{\text{eff}}\approx0.9\lambda.$$
意义：蒙皮屈曲后仍能传力（张力场）。
极限特例：厚蒙皮 $b_{\text{eff}}\to b$。
应用：航空薄壁加筋结构。

**知识点6：壳体无矩（薄膜）理论**
定义：忽略弯矩、剪力，壳体内只有面内薄膜力 $N_\phi,N_\theta$。
推导：旋转壳无矩平衡（法向）
$$\frac{N_\phi}{R_1}+\frac{N_\theta}{R_2}=p.$$
沿子午向
$$\frac{\mathrm d}{\mathrm d\phi}(N_\phi r)+N_\theta\frac{\mathrm dr}{\mathrm d\phi}\cdot(\cdots)=0.$$
由整体 $z$ 向平衡得
$$\frac{N_\phi}{R_1}+\frac{N_\theta}{R_2}=p\ (\text{拉普拉斯方程}).$$
意义：薄膜解简单但只在大范围光滑区域有效。
极限特例：球壳 $N_\phi=N_\theta=pR/2$。
应用：压力容器、穹顶。

**知识点7：旋转壳薄膜解与储罐**
定义：圆柱壳 $R_1=\infty,R_2=R$：
$$N_\theta=pR,\qquad N_x=-\frac{\text{轴向合力}}{2\pi R}.$$
推导：由拉普拉斯方程 $\dfrac{N_x}{\infty}+\dfrac{N_\theta}{R}=p\Rightarrow N_\theta=pR$；轴向由端盖压力平衡得 $N_x=pR/2$（封闭筒体）。
应力
$$\sigma_\theta=\frac{pR}{t},\qquad\sigma_x=\frac{pR}{2t}.$$
意义：环向应力是轴向两倍。
极限特例：开口筒体 $N_x=0$。
应用：气瓶、火箭贮箱。

**知识点8：壳体稳定性简介**
定义：薄壳受压力易屈曲，临界应力远低于材料屈服。
推导（经典球壳）：经典理论 $\sigma_{cr}=\dfrac{E}{\sqrt{3(1-\nu^2)}}\dfrac{t}{R}$；实际因初始缺陷下降
$$\sigma_{cr}^{\text{real}}\approx\gamma\,\sigma_{cr},\qquad \gamma\sim0.1\!-\!0.3.$$
$\because$ 缺陷引发局部失稳并扩展，$\therefore$ 用 Knockdown 因子 $\gamma$。
意义：壳体稳定性是航空结构关键。
极限特例：加筋、加压可提高稳定。
应用：火箭壳段、机身壁板。

### 第12章 稳定性、屈曲与分岔

**知识点1：平衡稳定性概念**
定义：受扰后能否回到原平衡；Lyapunov 稳定、渐近稳定、不稳定。
推导（能量判据）：保守系平衡 $\delta V=0$。若 $\delta^2V>0$ 对任意非零扰动，则势能极小、稳定；若存在方向 $\delta^2V<0$，则可释放势能、不稳定。
$$\because V(\mathbf q_0+\boldsymbol\eta)\approx V_0+\tfrac12\boldsymbol\eta^{\mathsf T}[\partial^2V]\boldsymbol\eta,$$
$$\therefore \text{稳定}\Leftrightarrow[\partial^2V]\ \text{正定}.$$
意义：稳定由势能二阶变分判定。
极限特例：$\det[\partial^2V]=0$ 为临界（分岔点）。
应用：压杆、壳、机构。

**知识点2：欧拉压杆屈曲**
定义：两端铰支等截面压杆临界载荷
$$P_{cr}=\frac{\pi^2EI}{L^2}.$$
推导：设挠度 $y(x)$，弯矩 $M=-Py$（$P$ 为压力）。由 $EIy''=M$：
$$EIy''+Py=0\ \Rightarrow\ y''+k^2y=0,\quad k^2=\frac{P}{EI}.$$
通解
$$y=A\sin kx+B\cos kx.$$
边界 $y(0)=0\Rightarrow B=0$；$y(L)=0\Rightarrow A\sin kL=0$。非零解要求
$$\sin kL=0\ \Rightarrow\ kL=n\pi.$$
最小非零载荷（$n=1$）
$$P_{cr}=\frac{\pi^2EI}{L^2}.$$
意义：屈曲是特征值问题；中性平衡出现。
极限特例：$P<P_{cr}$ 直线稳定；$P=P_{cr}$ 任意幅值。
应用：柱、撑杆。

**知识点3：不同边界条件的临界载荷**
定义：有效长度 $L_e=\mu L$，$P_{cr}=\pi^2EI/L_e^2$。
推导：对一般边界，令 $y=A\sin kx+B\cos kx+Cx+D$，代入四个边界条件，非零解 $\det=0$ 给出最小 $k$。结果 $\mu$：
两端铰支 $1$；两端固支 $0.5$；一端固支一端自由 $2$；一端固支一端铰支 $0.7$。
例（两端固支）：$y(0)=y'(0)=0,\ y(L)=y'(L)=0$ 得 $kL=2\pi$，$\therefore P_{cr}=4\pi^2EI/L^2$。
意义：约束越强临界载荷越高。
极限特例：悬臂柱最弱（$\mu=2$）。
应用：支架、机身框。

**知识点4：能量法求屈曲**
定义：用势能二阶变分或瑞利–里兹求 $P_{cr}$。
推导：屈曲时总势能
$$\Pi=\frac12\int EI(y'')^2\mathrm dx-\frac{P}{2}\int (y')^2\mathrm dx.$$
取 $y\approx A\phi(x)$，则
$$\Pi=\frac{A^2}{2}\!\left[\int EI(\phi'')^2\mathrm dx-P\int(\phi')^2\mathrm dx\right].$$
失稳条件（$\Pi$ 由正变负）：
$$P_{cr}=\frac{\int EI(\phi'')^2\mathrm dx}{\int(\phi')^2\mathrm dx}.$$
取精确屈曲模态可给精确解；近似模态给上界。
意义：能量法求临界载荷便于复杂结构。
极限特例：取正弦模态得欧拉解。
应用：变截面柱、加筋板。

**知识点5：初始缺陷与折减**
定义：实际杆有初弯曲 $y_0$，载荷–挠度非线性。
推导：设 $y_0=a\sin(\pi x/L)$，总挠度满足
$$EIy''+P(y+y_0)=0.$$
代入正弦展开：
$$A\!\left(\frac{\pi^2EI}{L^2}-P\right)=P a\ \Rightarrow\ A=\frac{a}{P_{cr}/P-1}.$$
最大总挠度 $=a+A$，弯矩 $M_{\max}=P(a+A)$。当 $P\to P_{cr}$ 时挠度发散（渐近）。
意义：缺陷使屈曲呈渐进失稳而非突变。
极限特例：$a\to0$ 回欧拉。
应用：实测缺陷折减（Perry 公式）。

**知识点6：板与壳的屈曲**
定义：四边简支矩形板单向受压临界应力
$$\sigma_{cr}=\frac{\pi^2D}{b^2 h}\!\left(m\frac{b}{a}+\frac{n^2 a}{m b}\right)^2.$$
推导：设 $w=\sin\dfrac{m\pi x}{a}\sin\dfrac{n\pi y}{b}$，代入板屈曲方程
$$D\nabla^4w+N_x\frac{\partial^2w}{\partial x^2}=0$$
得
$$N_x=\frac{\pi^2D}{b^2}\!\left(m\frac{b}{a}+\frac{n^2a}{mb}\right)^2.$$
取 $m,n$ 最小值得临界。
意义：板屈曲应力随 $a/b$ 振荡趋于 $4\pi^2D/(b^2h)$。
极限特例：$a/b$ 较大时 $k_{\min}=4$。
应用：壁板、腹板。

**知识点7：颤振与动态稳定性**
定义：气动力与弹性耦合导致自激振动，超过临界速度振幅发散。
推导（二元机翼简化）：俯仰 $\alpha$ 与沉浮 $h$ 方程
$$m\ddot h+S_\alpha\ddot\alpha+K_h h=-L,\qquad S_\alpha\ddot h+I_\alpha\ddot\alpha+K_\alpha\alpha=M.$$
气动力 $L=\rho V^2 c\,C_L(\alpha,\dot\alpha),\ M=\rho V^2c^2 C_M$。设解 $e^{\mathrm i\omega t}$，得复特征方程
$$\det([K]-[M]\omega^2-\mathrm i\omega[C]_{\text{aero}}-\rho V^2[G])=0.$$
当特征根实部由负变正（$\mathrm{Re}\,\lambda=0$）时颤振发生，对应临界速度 $V_f$。
意义：颤振是气动弹性失稳。
极限特例：低速时阻尼使振动衰减。
应用：机翼、桥梁（塔科马）。

**知识点8：分岔与突变**
定义：参数变化使平衡解数目/稳定性改变的点称分岔点。
推导（一维）：平衡方程 $f(x,\lambda)=0$。设分支，隐函数定理条件 $\partial f/\partial x\ne0$ 时解唯一；当
$$f=0,\qquad \frac{\partial f}{\partial x}=0$$
同时成立为分岔（极限点）条件。鞍结分岔：$f=x^2-\lambda=0\Rightarrow x=\pm\sqrt\lambda$；叉式分岔：$f=x^3-\lambda x=0\Rightarrow x=0,\pm\sqrt\lambda$。
意义：屈曲是叉式分岔的力学体现。
极限特例：$\lambda<0$ 无实分岔。
应用：结构失稳、非线性动力学。

### 第13章 流体力学基础与 Navier–Stokes 方程

**知识点1：连续方程**
定义：质量守恒
$$\frac{\partial\rho}{\partial t}+\nabla\!\cdot\!(\rho\mathbf u)=0.$$
推导：取控制体 $V$，质量变化率 = 净流入
$$\frac{\mathrm d}{\mathrm dt}\int_V\rho\,\mathrm dV=-\oint_S\rho\mathbf u\!\cdot\!\mathbf n\,\mathrm dS.$$
对固定控制体交换微分积分，用高斯定理：
$$\int_V\!\left[\frac{\partial\rho}{\partial t}+\nabla\!\cdot\!(\rho\mathbf u)\right]\mathrm dV=0.$$
$\because V$ 任意，$\therefore$ 得连续方程。
不可压 $\rho=$ 常量时
$$\nabla\!\cdot\!\mathbf u=0.$$
意义：质量守恒的局部形式。
极限特例：一维定常 $\rho uA=$ 常量。
应用：管道流量、CFD。

**知识点2：动量方程与应力张量**
定义：流体微元动量方程
$$\rho\frac{\mathrm D\mathbf u}{\mathrm Dt}=\nabla\!\cdot\!\boldsymbol\sigma+\rho\mathbf f.$$
推导：对控制体动量守恒
$$\frac{\mathrm d}{\mathrm dt}\int_V\rho\mathbf u\,\mathrm dV=-\oint_S\rho\mathbf u(\mathbf u\!\cdot\!\mathbf n)\mathrm dS+\oint_S\boldsymbol\sigma\!\cdot\!\mathbf n\,\mathrm dS+\int_V\rho\mathbf f\,\mathrm dV.$$
整理得 Cauchy 动量方程
$$\rho\frac{\partial\mathbf u}{\partial t}+\rho(\mathbf u\!\cdot\!\nabla)\mathbf u=\nabla\!\cdot\!\boldsymbol\sigma+\rho\mathbf f.$$
其中物质导数
$$\frac{\mathrm D}{\mathrm Dt}=\frac{\partial}{\partial t}+\mathbf u\!\cdot\!\nabla.$$
意义：$\nabla\!\cdot\!\boldsymbol\sigma$ 为面力散度。
极限特例：无黏 $\boldsymbol\sigma=-p\mathbf I$ 得欧拉方程。
应用：流固耦合、湍流建模。

**知识点3：牛顿流体本构与 Navier–Stokes 方程**
定义：牛顿流体
$$\boldsymbol\sigma=-p\mathbf I+\boldsymbol\tau,\qquad \tau_{ij}=\mu\!\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\right)+\lambda'\frac{\partial u_k}{\partial x_k}\delta_{ij}.$$
推导：Stokes 假设 $\lambda'=-\tfrac23\mu$（体积黏度为零）。代入 Cauchy 方程：
$$\rho\frac{\mathrm D\mathbf u}{\mathrm Dt}=-\nabla p+\mu\nabla^2\mathbf u+\big(\lambda'+\mu\big)\nabla(\nabla\!\cdot\!\mathbf u)+\rho\mathbf f.$$
不可压 $\nabla\!\cdot\!\mathbf u=0$：
$$\boxed{\rho\!\left(\frac{\partial\mathbf u}{\partial t}+\mathbf u\!\cdot\!\nabla\mathbf u\right)=-\nabla p+\mu\nabla^2\mathbf u+\rho\mathbf f.}$$
意义：黏性流体运动的基本方程。
极限特例：$\mu=0$ 欧拉方程；$\mathrm{Re}\ll1$ 斯托克斯方程。
应用：管流、边界层、CFD。

**知识点4：伯努利方程**
定义：理想不可压定常流沿流线
$$\frac{p}{\rho}+\frac12u^2+gz=\text{常量}.$$
推导：欧拉方程沿流线 $\mathrm d\mathbf r=\mathbf u\,\mathrm dt$。定常、无黏、不可压：
$$\rho(\mathbf u\!\cdot\!\nabla)\mathbf u=-\nabla p-\rho g\hat z.$$
点乘 $\mathbf u$，利用 $(\mathbf u\!\cdot\!\nabla)\mathbf u\!\cdot\!\mathbf u=\tfrac12\mathbf u\!\cdot\!\nabla(u^2)$：
$$\frac12\rho\,\mathbf u\!\cdot\!\nabla(u^2)=-\mathbf u\!\cdot\!\nabla p-\rho g u_z.$$
沿流线积分
$$\int\!\left[\frac12\rho\,\mathrm d(u^2)+\mathrm dp+\rho g\,\mathrm dz\right]=0\ \Rightarrow\ \frac{p}{\rho}+\frac12u^2+gz=\text{常量}.$$
意义：机械能沿流线守恒。
极限特例：静止流体 $u=0\Rightarrow p+\rho gz=$ 常量（静压）。
应用：文丘里、机翼升力（无黏近似）。

**知识点5：涡量方程与亥姆霍兹定理**
定义：涡量 $\boldsymbol\omega=\nabla\times\mathbf u$。
推导：对不可压 N–S 取旋度：
$$\frac{\partial\boldsymbol\omega}{\partial t}=\nabla\times(\mathbf u\times\boldsymbol\omega)+\nu\nabla^2\boldsymbol\omega.$$
利用矢量恒等式
$$\nabla\times(\mathbf u\times\boldsymbol\omega)=(\boldsymbol\omega\!\cdot\!\nabla)\mathbf u-(\mathbf u\!\cdot\!\nabla)\boldsymbol\omega,$$
得
$$\frac{\mathrm D\boldsymbol\omega}{\mathrm Dt}=(\boldsymbol\omega\!\cdot\!\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega.$$
意义：涡管随流体拉伸/倾斜；无黏、正压、有势体力时涡量守恒（亥姆霍兹）。
极限特例：$\nu=0$ 且初始无旋 $\Rightarrow$ 永远无旋（拉格朗日定理）。
应用：翼尖涡、湍流。

**知识点6：层流边界层（Blasius）**
定义：平板边界层方程
$$u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=\nu\frac{\partial^2u}{\partial y^2},\qquad \frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0.$$
假设：$\mathrm{Re}\gg1$，边界层厚度 $\delta\ll L$，压力沿法向不变。
推导（量级）：$u\sim U,\ x\sim L,\ y\sim\delta$，连续方程给 $v\sim U\delta/L$。惯性项 $\sim U^2/L$，黏性项 $\sim\nu U/\delta^2$，平衡
$$\frac{U^2}{L}\sim\frac{\nu U}{\delta^2}\ \Rightarrow\ \delta\sim\sqrt{\frac{\nu L}{U}}=\frac{L}{\sqrt{\mathrm{Re}}}.$$
引入相似变量 $\eta=y\sqrt{U/(\nu x)}$ 与流函数 $\psi=\sqrt{\nu Ux}f(\eta)$，得 Blasius 方程
$$2f'''+ff''=0.$$
边界 $f(0)=f'(0)=0,\ f'(\infty)=1$。数值解 $f''(0)=0.332$。
壁面剪应力
$$\tau_w=\mu\!\left.\frac{\partial u}{\partial y}\right|_{y=0}=0.332\rho U^2\sqrt{\frac{\nu}{Ux}}.$$
意义：边界层是黏性主导薄层。
极限特例：后缘 $\delta$ 增厚，可能分离。
应用：摩阻、换热。

**知识点7：湍流与雷诺平均**
定义：速度分解 $\mathbf u=\bar{\mathbf u}+\mathbf u'$，$\overline{\mathbf u'}=0$。
推导：对 N–S 取时间平均。不可压时连续方程
$$\frac{\partial\bar u_i}{\partial x_i}=0.$$
动量方程平均后出现雷诺应力
$$\rho\frac{\partial\bar u_i}{\partial t}+\rho\bar u_j\frac{\partial\bar u_i}{\partial x_j}=-\frac{\partial\bar p}{\partial x_i}+\frac{\partial}{\partial x_j}\!\left(\mu\frac{\partial\bar u_i}{\partial x_j}-\rho\overline{u_i'u_j'}\right).$$
湍动能方程可由脉动方程乘 $u_i'$ 平均得到。
意义：雷诺应力需湍流模型封闭（$k$–$\varepsilon$、SST）。
极限特例：层流 $\overline{u_i'u_j'}=0$。
应用：湍流 CFD、风载荷。

**知识点8：可压缩流动与激波**
定义：马赫数 $\mathrm{Ma}=u/c$，$c=\sqrt{\gamma RT}$，$\gamma=c_p/c_v$。
推导（等熵关系）：由 $p/\rho^\gamma=$ 常量与状态方程 $p=\rho RT$，
$$\frac{T_0}{T}=1+\frac{\gamma-1}{2}\mathrm{Ma}^2,$$
$$\frac{p_0}{p}=\left(1+\frac{\gamma-1}{2}\mathrm{Ma}^2\right)^{\gamma/(\gamma-1)},\qquad \frac{\rho_0}{\rho}=\left(1+\frac{\gamma-1}{2}\mathrm{Ma}^2\right)^{1/(\gamma-1)}.$$
正激波关系（Rankine–Hugoniot）：质量、动量、能量守恒联立得
$$\frac{p_2}{p_1}=\frac{2\gamma\mathrm{Ma}_1^2-(\gamma-1)}{\gamma+1},\qquad
\mathrm{Ma}_2^2=\frac{(\gamma-1)\mathrm{Ma}_1^2+2}{2\gamma\mathrm{Ma}_1^2-(\gamma-1)}.$$
意义：$\mathrm{Ma}>1$ 出现激波、膨胀波。
极限特例：$\mathrm{Ma}\to0$ 回到不可压。
应用：超声速进气道、再入。

### 第14章 实验力学与测量原理

**知识点1：电阻应变片原理**
定义：应变片电阻 $R=\rho_e L/A$，应变引起 $R$ 变化。
推导：
$$\frac{\mathrm dR}{R}=\frac{\mathrm d\rho_e}{\rho_e}+\frac{\mathrm dL}{L}-\frac{\mathrm dA}{A}.$$
对圆丝半径 $r$，$A=\pi r^2$，$\dfrac{\mathrm dA}{A}=2\dfrac{\mathrm dr}{r}=-2\nu\dfrac{\mathrm dL}{L}$，
$$\therefore \frac{\mathrm dR/R}{\varepsilon}=\underbrace{\frac{\mathrm d\rho_e/\rho_e}{\varepsilon}}_{1}+2\nu.$$
对金属，压阻项近似 $0.4$，$1+2\nu\approx1.6$，敏感系数
$$GF=\frac{\Delta R/R}{\varepsilon}\approx2.$$
意义：通过电阻变化测应变。
极限特例：半导体应变片压阻项很大，$GF\sim100$。
应用：结构应变测量。

**知识点2：惠斯通电桥与温度补偿**
定义：电桥输出
$$U_o=U_s\!\left(\frac{R_1}{R_1+R_2}-\frac{R_4}{R_3+R_4}\right).$$
推导：等臂 $R_1=R_2=R_3=R_4=R$，单臂工作 $R_1\to R+\Delta R$：
$$U_o=U_s\!\left(\frac{R+\Delta R}{2R+\Delta R}-\frac12\right)\approx U_s\frac{\Delta R}{4R}=\frac{U_s}{4}GF\,\varepsilon.$$
温度补偿：用工作片与补偿片接成半桥/全桥，使温度引起同向电阻变化相消。
意义：电桥把微小电阻变化转为电压。
极限特例：$R$ 变化线性仅当 $\Delta R\ll R$。
应用：应变仪、传感器。

**知识点3：光弹性法**
定义：透明双折射材料在应力下折射率随主应力变化（应力光学定律）
$$n_1-n_2=C(\sigma_1-\sigma_2).$$
推导：光通过厚度 $h$ 产生相位差
$$\delta=\frac{2\pi h}{\lambda}(n_1-n_2)=\frac{2\pi h C}{\lambda}(\sigma_1-\sigma_2).$$
等色线条纹 $N=\delta/(2\pi)$：
$$\sigma_1-\sigma_2=\frac{N\lambda}{hC}=\frac{N f_\sigma}{h},$$
$f_\sigma=\lambda/C$ 为材料条纹值。
意义：直接显示主应力差等值线。
极限特例：$N=0$ 表示主应力相等。
应用：应力集中系数、接触应力。

**知识点4：数字图像相关（DIC）**
定义：通过匹配变形前后子区灰度求位移场。
推导：设参考子区 $f(x,y)$、变形子区 $g(x',y')$，归一化互相关
$$C=\frac{\sum[f(\xi,\eta)-\bar f][g(\xi+u,\eta+v)-\bar g]}{\sqrt{\sum[f-\bar f]^2\sum[g-\bar g]^2}}.$$
位移 $u,v$ 使 $C$ 最大；应变由位移梯度中心差分
$$\varepsilon_{xx}=\frac{\partial u}{\partial x}\approx\frac{u(i+1)-u(i-1)}{2\Delta x}.$$
意义：非接触全场测量，可后处理得应变。
极限特例：亚像素插值提高精度至 $0.01$ 像素。
应用：断裂、高温变形。

**知识点5：模态测试与频响函数**
定义：频响函数 $H(\omega)=X(\omega)/F(\omega)$。
推导：多自由度系统 $M\ddot{\mathbf x}+C\dot{\mathbf x}+K\mathbf x=\mathbf f$，傅里叶变换：
$$\big[-\omega^2M+\mathrm i\omega C+K\big]\mathbf X=\mathbf F.$$
$$\therefore H_{ij}=\frac{X_i}{F_j}=\sum_{r=1}^{N}\frac{\phi_{ir}\phi_{jr}}{m_r(\omega_r^2-\omega^2+2\mathrm i\zeta_r\omega_r\omega)}.$$
峰值对应 $\omega\approx\omega_r$。半功率带宽
$$\zeta_r=\frac{\omega_2-\omega_1}{2\omega_r}.$$
意义：由 $H$ 识别固有频率、阻尼、振型。
极限特例：单自由度 $H=1/(k-m\omega^2+\mathrm i c\omega)$。
应用：模态试验、模型修正。

**知识点6：振动传感器原理**
定义：压电加速度计输出电荷 $q=d_{ij}F$，$F=m a$。
推导：
$$q=d_{33}ma\ \Rightarrow\ V_o=\frac{q}{C}=K_a a,\qquad K_a=\frac{d_{33}m}{C}.$$
频率范围受安装共振限制，可用弹簧–质量模型
$$\frac{a_{\text{meas}}}{a_{\text{base}}}=\frac{\omega_n^2}{\sqrt{(\omega_n^2-\omega^2)^2+(2\zeta\omega_n\omega)^2}}.$$
意义：在 $\omega\ll\omega_n$ 时响应平坦。
极限特例：$\omega/\omega_n<0.2$ 时误差 $<2\%$。
应用：振动台、冲击测试。

**知识点7：误差与不确定度**
定义：测量误差 = 真值 − 测量值；不确定度合成
$$u_c^2=\sum_i\!\left(\frac{\partial f}{\partial x_i}\right)^2u_{x_i}^2+2\sum_{i<j}\frac{\partial f}{\partial x_i}\frac{\partial f}{\partial x_j}\operatorname{cov}(x_i,x_j).$$
推导：对 $y=f(x_1,\dots,x_n)$ 在均值处做一阶泰勒展开：
$$\Delta y\approx\sum_i\frac{\partial f}{\partial x_i}\Delta x_i.$$
取方差并假设独立：
$$\therefore u_c^2=\sum_i\!\left(\frac{\partial f}{\partial x_i}\right)^2u_{x_i}^2.$$
扩展不确定度 $U=k u_c$（$k=2$ 约 95%）。
意义：A 类（统计）与 B 类（非统计）合成。
极限特例：线性函数 $y=\sum a_ix_i\Rightarrow u_c^2=\sum a_i^2u_i^2$。
应用：实验报告、标定证书。

**知识点8：数据拟合与模型验证**
定义：最小二乘拟合 $\min\sum_i[y_i-f(x_i;\boldsymbol\beta)]^2$。
推导：对 $\boldsymbol\beta$ 求极小
$$\frac{\partial}{\partial\beta_j}\sum_i[y_i-f_i]^2=-2\sum_i[y_i-f_i]\frac{\partial f_i}{\partial\beta_j}=0.$$
线性模型 $f=\mathbf x^{\mathsf T}\boldsymbol\beta$ 得正规方程
$$\therefore \boldsymbol\beta=(X^{\mathsf T}X)^{-1}X^{\mathsf T}\mathbf y.$$
可决系数
$$R^2=1-\frac{\sum(y_i-\hat y_i)^2}{\sum(y_i-\bar y)^2}.$$
意义：参数辨识与模型验证。
极限特例：残差正态时 LS 等价 MLE。
应用：本构参数反演、仿真对标。

### 第15章 计算力学：有限元、有限差分与数值方法

**知识点1：微分方程的弱形式**
定义：将强形式两边乘试探函数并积分，分部积分降低导数阶。
推导（一维泊松）：强形式
$$-\frac{\mathrm d}{\mathrm dx}\!\left(k\frac{\mathrm du}{\mathrm dx}\right)=f,\quad u(0)=0,\ k u'(L)=\bar t.$$
乘 $v$（$v(0)=0$）积分：
$$\int_0^L v\!\left[-\frac{\mathrm d}{\mathrm dx}\!\left(k u'\right)\right]\mathrm dx=\int_0^L fv\,\mathrm dx.$$
分部积分
$$-\!\left[vu'k\right]_0^L+\int_0^L k u'v'\,\mathrm dx=\int_0^L fv\,\mathrm dx.$$
$$\therefore \int_0^L k u'v'\,\mathrm dx=\int_0^L fv\,\mathrm dx+k\bar t\,v(L).$$
意义：只要求 $u\in H^1$，降低连续性要求。
极限特例：自然边界自动进入弱形式。
应用：FEM 基础。

**知识点2：伽辽金加权残量与有限元离散**
定义：$u_h=\sum_j N_j u_j$，取权 $v=N_i$。
推导：代入弱形式
$$\sum_j\!\left[\int_0^L k N_i'N_j'\,\mathrm dx\right]u_j=\int_0^L fN_i\,\mathrm dx+k\bar t\,N_i(L).$$
$$\therefore K_{ij}=\int_0^L kN_i'N_j'\,\mathrm dx,\qquad F_i=\int_0^L fN_i\,\mathrm dx+k\bar t\,N_i(L).$$
意义：刚度矩阵 $K$ 对称、稀疏、正定。
极限特例：$N_j$ 线性形函数得三对角 $K$。
应用：结构分析、热传导。

**知识点3：单元插值与形函数**
定义：一维两节点线性单元，节点坐标 $x_1,x_2$，长度 $h=x_2-x_1$。
推导：设 $u(\xi)=a_0+a_1\xi$，$\xi\in[0,1]$。由 $u(0)=u_1,u(1)=u_2$：
$$a_0=u_1,\quad a_1=u_2-u_1.$$
$$\therefore u(\xi)=(1-\xi)u_1+\xi u_2=N_1u_1+N_2u_2,$$
$$N_1=1-\xi,\qquad N_2=\xi.$$
单元刚度
$$K^e=\frac{k}{h}\begin{pmatrix}1&-1\\-1&1\end{pmatrix}.$$
意义：形函数满足 $N_i(x_j)=\delta_{ij}$。
极限特例：$h\to0$ 收敛。
应用：二维三角/四边形单元推广。

**知识点4：刚度矩阵与整体组装**
定义：整体 $K=\sum_e A_e^{\mathsf T}K^eA_e$，$A_e$ 为布尔（组装）矩阵。
推导：单元内能
$$U_e=\tfrac12\mathbf u_e^{\mathsf T}K^e\mathbf u_e.$$
全局位移 $\mathbf u$ 与单元位移 $\mathbf u_e=A_e\mathbf u$：
$$U=\sum_e\tfrac12(A_e\mathbf u)^{\mathsf T}K^e(A_e\mathbf u)=\tfrac12\mathbf u^{\mathsf T}\!\left(\sum_e A_e^{\mathsf T}K^eA_e\right)\mathbf u.$$
$$\therefore K=\sum_e A_e^{\mathsf T}K^eA_e.$$
意义：组装体现位移协调。
极限特例：节点重合处刚度叠加。
应用：大型结构 FEM。

**知识点5：收敛性与误差估计**
定义：$C^0$ 线性单元在能量范数下误差 $O(h)$，$L^2$ 范数下 $O(h^2)$。
推导：插值误差
$$\|u-u_h\|_{L^2}\le Ch^2\|u''\|_{L^2},\qquad |u-u_h|_{H^1}\le Ch\|u''\|_{L^2}.$$
由 Céa 引理
$$\|u-u_h\|_E=\min_{v_h}\|u-v_h\|_E,$$
结合插值估计得能量误差 $O(h)$。$p$ 次单元误差 $O(h^p)$。
意义：网格加密提高精度。
极限特例：奇异点使收敛阶下降。
应用：网格收敛性研究。

**知识点6：显式与隐式时间积分**
定义：半离散方程 $M\ddot{\mathbf u}+C\dot{\mathbf u}+K\mathbf u=\mathbf F$。
推导：
（中心差分显式）：
$$\ddot{\mathbf u}^n\approx\frac{\mathbf u^{n+1}-2\mathbf u^n+\mathbf u^{n-1}}{\Delta t^2},\quad\dot{\mathbf u}^n\approx\frac{\mathbf u^{n+1}-\mathbf u^{n-1}}{2\Delta t}.$$
代入得 $\mathbf u^{n+1}$ 显式；稳定性条件
$$\Delta t\le\Delta t_{cr}=\frac{2}{\omega_{\max}}\quad(\text{无阻尼}).$$
（Newmark 隐式）：取 $\gamma=1/2,\beta=1/4$ 无条件稳定。
意义：显式适合冲击/波传播，隐式适合准静态。
极限特例：$\Delta t>\Delta t_{cr}$ 显式发散。
应用：碰撞、爆炸仿真。

**知识点7：几何非线性与牛顿–拉夫森**
定义：残量 $\mathbf R(\mathbf u)=\mathbf F_{\text{ext}}-\mathbf F_{\text{int}}(\mathbf u)=0$。
推导：泰勒展开
$$\mathbf R(\mathbf u+\Delta\mathbf u)\approx\mathbf R(\mathbf u)+\frac{\partial\mathbf R}{\partial\mathbf u}\Delta\mathbf u.$$
切线刚度
$$K_T=-\frac{\partial\mathbf R}{\partial\mathbf u}=\frac{\partial\mathbf F_{\text{int}}}{\partial\mathbf u}.$$
迭代
$$\therefore K_T\Delta\mathbf u=\mathbf R(\mathbf u),\qquad \mathbf u\leftarrow\mathbf u+\Delta\mathbf u.$$
意义：二次收敛（接近解时）。
极限特例：线性问题一步收敛。
应用：大变形、屈曲后路径。

**知识点8：有限差分与 CFL 条件**
定义：空间导数用差商
$$u_{i}'\approx\frac{u_{i+1}-u_{i-1}}{2\Delta x},\qquad u_i''\approx\frac{u_{i+1}-2u_i+u_{i-1}}{\Delta x^2}.$$
推导（一维对流）：$u_t+au_x=0$，显式迎风
$$u_i^{n+1}=u_i^n-\frac{a\Delta t}{\Delta x}(u_i^n-u_{i-1}^n)\quad(a>0).$$
von Neumann 稳定性分析设 $u_i^n=\xi^n e^{\mathrm i k i\Delta x}$，放大因子
$$G=1-\lambda(1-e^{-\mathrm i k\Delta x}),\qquad \lambda=\frac{a\Delta t}{\Delta x}.$$
$|G|\le1$ 要求
$$0\le\lambda\le1\ \Rightarrow\ \Delta t\le\frac{\Delta x}{|a|}\quad(\text{CFL 条件}).$$
意义：显式格式时间步受空间步长限制。
极限特例：$\lambda=1$ 精确输运。
应用：CFD、波传播。

### 第16章 综合实践：振动、冲击、疲劳与结构优化

**知识点1：单自由度振动**
定义：$m\ddot x+c\dot x+kx=F(t)$。
推导：无阻尼自由振动 $m\ddot x+kx=0$，设 $x=e^{\lambda t}$：
$$\lambda^2+\omega_n^2=0,\quad\omega_n=\sqrt{\frac km}\ \Rightarrow\ x=A\cos\omega_nt+B\sin\omega_nt.$$
有阻尼，令 $\zeta=\dfrac{c}{2\sqrt{mk}}$：
$$\lambda=-\zeta\omega_n\pm\mathrm i\omega_n\sqrt{1-\zeta^2},$$
$$x=e^{-\zeta\omega_nt}\big(A\cos\omega_dt+B\sin\omega_dt\big),\quad\omega_d=\omega_n\sqrt{1-\zeta^2}.$$
意义：$\zeta<1$ 欠阻尼振荡衰减；$\zeta=1$ 临界；$\zeta>1$ 过阻尼。
极限特例：$\zeta=0$ 等幅；$\zeta\ge1$ 不振荡。
应用：隔振、阻尼器。

**知识点2：多自由度系统与模态叠加**
定义：$M\ddot{\mathbf x}+K\mathbf x=\mathbf F$。
推导：广义特征问题 $K\boldsymbol\phi=\omega^2M\boldsymbol\phi$，$\boldsymbol\phi_r$ 满足 $M$ 正交
$$\boldsymbol\phi_r^{\mathsf T}M\boldsymbol\phi_s=m_r\delta_{rs},\qquad \boldsymbol\phi_r^{\mathsf T}K\boldsymbol\phi_s=k_r\delta_{rs}=m_r\omega_r^2\delta_{rs}.$$
令 $\mathbf x=\sum_r\boldsymbol\phi_r q_r$，代入并前乘 $\boldsymbol\phi_s^{\mathsf T}$：
$$m_s\ddot q_s+m_s\omega_s^2q_s=\boldsymbol\phi_s^{\mathsf T}\mathbf F\ \Rightarrow\ \ddot q_s+\omega_s^2q_s=\frac{\boldsymbol\phi_s^{\mathsf T}\mathbf F}{m_s}.$$
$$\therefore q_s(t)=\frac{1}{m_s\omega_s}\int_0^t(\boldsymbol\phi_s^{\mathsf T}\mathbf F)\sin\omega_s(t-\tau)\,\mathrm d\tau.$$
意义：解耦为多个单自由度。
极限特例：无阻尼振型叠加。
应用：地震响应、振动控制。

**知识点3：连续系统振动**
定义：杆纵向振动方程
$$\frac{\partial}{\partial x}\!\left(EA\frac{\partial u}{\partial x}\right)=\rho A\frac{\partial^2u}{\partial t^2}.$$
推导（微元）：取长 $\mathrm dx$ 微元，轴向力差 $N(x+\mathrm dx)-N(x)=\dfrac{\partial N}{\partial x}\mathrm dx$，$N=EA\partial u/\partial x$，惯性力 $\rho A\,\mathrm dx\,\ddot u$，平衡：
$$\frac{\partial}{\partial x}\!\left(EA\frac{\partial u}{\partial x}\right)\mathrm dx=\rho A\frac{\partial^2u}{\partial t^2}\mathrm dx.$$
等截面时
$$c^2\frac{\partial^2u}{\partial x^2}=\frac{\partial^2u}{\partial t^2},\qquad c=\sqrt{E/\rho}.$$
分离变量 $u=X(x)T(t)$ 得
$$X''+k^2X=0,\quad T''+\omega^2T=0,\quad\omega=ck.$$
边界给出特征频率 $\omega_n$。
意义：连续体有无穷多固有频率。
极限特例：两端固定杆 $\omega_n=\dfrac{n\pi}{L}\sqrt{E/\rho}$。
应用：杆、梁、板振动。

**知识点4：冲击响应谱**
定义：单自由度系统对基座加速度冲击 $\ddot y(t)$ 的最大响应随固有频率变化曲线。
推导：方程
$$\ddot x+2\zeta\omega_n\dot x+\omega_n^2x=-\ddot y(t).$$
相对位移 $z=x-y$ 满足
$$\ddot z+2\zeta\omega_n\dot z+\omega_n^2z=-\ddot y.$$
Duhamel 积分
$$z(t)=-\frac{1}{\omega_d}\int_0^t\ddot y(\tau)e^{-\zeta\omega_n(t-\tau)}\sin\omega_d(t-\tau)\,\mathrm d\tau.$$
SRS 定义为
$$\text{SRS}(\omega_n)=\max_t|\omega_n^2 z(t)|\quad(\text{伪加速度}).$$
意义：描述冲击对结构的最坏影响。
极限特例：$\omega_n\to\infty$ 响应趋于零；低频趋于基座位移。
应用：航天器冲击试验规范。

**知识点5：疲劳 S–N 曲线与损伤累积**
定义：应力幅 $S$ 与失效循环数 $N$ 关系（Basquin）
$$S^aN=C\ \Rightarrow\ \lg N=\lg C-a\lg S.$$
推导：对幂律取对数即得线性 S–N。Miner 线性损伤累积：
$$D=\sum_i\frac{n_i}{N_i},\qquad D\ge1\ \text{估算失效}.$$
意义：变幅载荷下按循环比累加。
极限特例：恒幅时 $D=n/N$。
应用：机翼疲劳寿命、螺栓。

**知识点6：断裂力学简介**
定义：应力强度因子 $K_I=\beta\sigma\sqrt{\pi a}$。
推导（Griffith 能量）：裂纹扩展释放应变能
$$U_e=\frac{\pi\sigma^2a^2}{E'}\ (\text{平面应力 }E'=E,\ \text{平面应变 }E'=E/(1-\nu^2)),$$
新增表面能 $U_s=2\gamma a$。总能量
$$\Pi=\frac{\pi\sigma^2a^2}{E'}+2\gamma a.$$
失稳条件 $\dfrac{\partial\Pi}{\partial a}=0$：
$$\frac{2\pi\sigma^2a}{E'}+2\gamma=0\ \Rightarrow\ \sigma_f=\sqrt{\frac{E'\gamma}{\pi a}}.$$
关系
$$G=\frac{K^2}{E'},\qquad K=\sigma\sqrt{\pi a}.$$
意义：含裂纹结构剩余强度。
极限特例：$a\to0$ 强度趋于理论强度（受其他机制限制）。
应用：损伤容限设计。

**知识点7：结构优化与设计**
定义：$\min f(\mathbf x)$ s.t. $g_j(\mathbf x)\le0,\ h_k(\mathbf x)=0$。
推导（Lagrange）：构造
$$\mathcal L=f+\sum_j\mu_jg_j+\sum_k\lambda_kh_k.$$
KKT 条件
$$\nabla f+\sum_j\mu_j\nabla g_j+\sum_k\lambda_k\nabla h_k=0,\quad\mu_j\ge0,\ \mu_jg_j=0.$$
意义：约束极值的必要条件。
极限特例：无约束时 $\nabla f=0$。
应用：轻量化、拓扑优化。

**知识点8：实验–仿真联合验证**
定义：用实验数据校正仿真模型参数，使残差最小。
推导：设仿真输出 $\mathbf y_{\text{sim}}(\boldsymbol\theta)$，实验 $\mathbf y_{\text{exp}}$，目标
$$\min_{\boldsymbol\theta}\ J=\big(\mathbf y_{\text{sim}}-\mathbf y_{\text{exp}}\big)^{\mathsf T}W\big(\mathbf y_{\text{sim}}-\mathbf y_{\text{exp}}\big).$$
线性化 $\mathbf y_{\text{sim}}\approx\mathbf y_0+S\,\delta\boldsymbol\theta$（$S$ 灵敏度矩阵），正规方程
$$\therefore \delta\boldsymbol\theta=(S^{\mathsf T}WS)^{-1}S^{\mathsf T}W\!\left(\mathbf y_{\text{exp}}-\mathbf y_0\right).$$
意义：模型修正、参数反演、可信度评估。
极限特例：$S$ 列满秩时参数可辨识。
应用：模态修正、材料参数反演。

---
⬆ 上级领域：[[实践与创新]]
