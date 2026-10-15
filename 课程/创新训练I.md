---
课程名称: 创新训练I
学分: 1
总学时: 32
考核方式: 考查
课程属性: 必修
课程性质: 专业课
开课单位: 空天科学与技术学院
领域: 实践与创新
tags:
  - 课程
  - 领域/实践与创新
---

# 创新训练I


### 第1章 量纲分析与相似理论

**知识点1　物理量的量纲与基本量纲**

定义：在力学与航空航天问题中，取长度 $L$、质量 $M$、时间 $T$ 为基本量纲，任一物理量 $Q$ 的量纲记作
$$[Q]=L^{a}M^{b}T^{c},\qquad a,b,c\in\mathbb R.$$

推导运算规则：设 $Q_1=L^{a_1}M^{b_1}T^{c_1}$，$Q_2=L^{a_2}M^{b_2}T^{c_2}$。由量纲的可乘性
$$
\begin{aligned}
[Q_1Q_2]&=[Q_1][Q_2]=L^{a_1+a_2}M^{b_1+b_2}T^{c_1+c_2},\\
[Q_1/Q_2]&=L^{a_1-a_2}M^{b_1-b_2}T^{c_1-c_2},\\
[Q_1^n]&=L^{na_1}M^{nb_1}T^{nc_1}.
\end{aligned}
$$

物理/几何意义：量纲刻画物理量对基本尺度的依赖结构，是物理方程在单位变换下的不变特征。纯数、角度、应变、相似准则为无量纲量 $[Q]=L^0M^0T^0$。

极限与特例：速度 $[v]=LT^{-1}$，加速度 $[a]=LT^{-2}$，力 $[F]=MLT^{-2}$，能量 $[E]=ML^2T^{-2}$，压强 $[\sigma]=ML^{-1}T^{-2}$，动力黏度 $[\mu]=ML^{-1}T^{-1}$。

典型应用：判断经验公式是否可能成立，构造无量纲变量以简化问题。

**知识点2　量纲齐次性原理**

定义：任何有物理意义的方程，其相加各项必须具有相同量纲。

假设：定律形式不随基本单位的人为选取而改变。

推导：设有物理方程 $\sum_i f_i(Q_1,\dots,Q_n)=0$。作单位变换 $L'=\lambda_L L,\ M'=\lambda_M M,\ T'=\lambda_T T$，则量纲为 $L^{a_i}M^{b_i}T^{c_i}$ 的量数值变为原值的 $\lambda_L^{-a_i}\lambda_M^{-b_i}\lambda_T^{-c_i}$ 倍。要求方程形式不变：
$$\sum_i \lambda_L^{-a_i}\lambda_M^{-b_i}\lambda_T^{-c_i} f_i=0,\qquad \forall \lambda_L,\lambda_M,\lambda_T>0.$$
各 $f_i$ 相互独立，$\therefore$ 只能各项指数相同，$(a_i,b_i,c_i)=\text{const}$。

物理意义：量纲齐次性是物理定律的普适约束，可用于排除错误公式、确定指数。

特例：$s=v_0t+\tfrac12at^2$ 各项均为 $L$；$E=mc^2$ 两边均为 $ML^2T^{-2}$。

典型应用：校核实验拟合式，识别隐含的无量纲组合。

**知识点3　瑞利量纲法**

定义：设被求量 $Q$ 为自变量 $Q_1,\dots,Q_n$ 的幂次积
$$Q=k\,Q_1^{x_1}Q_2^{x_2}\cdots Q_n^{x_n},$$
由量纲齐次性解出指数 $x_i$。

推导：写量纲方程
$$[Q]=\prod_{i=1}^{n}[Q_i]^{x_i}.$$
对第 $j$ 个基本量纲取指数，得线性方程组
$$\sum_{i=1}^{n}a_{ji}x_i=a_j,\qquad j=1,\dots,m,$$
其中 $a_{ji}$ 为 $Q_i$ 对第 $j$ 个基本量纲的指数，$a_j$ 为 $Q$ 对应指数。写成矩阵形式
$$A\mathbf x=\mathbf a,\qquad A=(a_{ji})\in\mathbb R^{m\times n}.$$
$\because$ 方程数 $m$，未知数 $n$，$\therefore$ 当 $n=m$ 且 $A$ 满秩时解唯一。

几何意义：可用 ${[Q_i]}$ 的对数作坐标，量纲齐次性成为线性约束的超平面。

特例：单摆周期 $T_p$ 依赖 $l,m,g$，设 $T_p=kl^{x}m^{y}g^{z}$，则
$$T=L^{x}M^{y}(LT^{-2})^{z}=L^{x+z}M^{y}T^{-2z}.$$
比较得 $x+z=0,\ y=0,\ -2z=1$，$\Rightarrow z=-\tfrac12,\ x=\tfrac12,\ y=0$，$\therefore T_p=k\sqrt{l/g}$。

典型应用：估算周期、阻力、频率等对参数的依赖关系。

**知识点4　白金汉 $\pi$ 定理**

定义：若问题含 $n$ 个有量纲物理量、涉及 $m$ 个基本量纲，量纲矩阵 $A\in\mathbb R^{m\times n}$ 的秩为 $r$，则存在 $n-r$ 个相互独立的无量纲组合 $\pi_1,\dots,\pi_{n-r}$，且物理关系可化为
$$F(\pi_1,\dots,\pi_{n-r})=0.$$

推导：寻找指数向量 $\mathbf x=(x_1,\dots,x_n)^T$，使
$$\prod_{i=1}^{n}Q_i^{x_i}\ \text{无量纲}\quad\Longleftrightarrow\quad A\mathbf x=0.$$
这是齐次线性方程组。$\because \operatorname{rank}A=r$，由秩—零化度定理，
$$\dim\ker A=n-r.$$
取 $\ker A$ 的一组基 $\mathbf x^{(1)},\dots,\mathbf x^{(n-r)}$，令
$$\pi_k=\prod_{i=1}^{n}Q_i^{x_i^{(k)}},\qquad k=1,\dots,n-r,$$
则 $[\pi_k]=1$，且诸 $\pi_k$ 相互独立。$\therefore$ 结论成立。

物理/几何意义：$\pi$ 定理把 $n$ 维参数空间压缩到 $n-r$ 维无量纲空间，实验点数与维度大幅降低。

极限与特例：$n=r$ 时无独立无量纲量，关系退化为常数。

典型应用：圆管压降 $\Delta p=f(D,L,\rho,\mu,v)$，$n=6,m=3,r=3$，得
$$\frac{\Delta p}{\rho v^2}=\Phi\!\left(\frac{L}{D},\frac{\rho vD}{\mu}\right).$$

**知识点5　常见相似准则数**

定义：由 $\pi$ 定理得到的表征两种物理效应之比的无量纲数。

推导（以雷诺数为例）：不可压 N–S 方程
$$\rho\frac{D\mathbf v}{Dt}=-\nabla p+\mu\nabla^2\mathbf v.$$
取无量纲化 $\mathbf x=L\mathbf x^*,\ \mathbf v=V\mathbf v^*,\ t=\frac{L}{V}t^*,\ p=\rho V^2p^*$，代入并整理：
$$
\begin{aligned}
\rho\frac{V^2}{L}\frac{D\mathbf v^*}{Dt^*}&=-\frac{\rho V^2}{L}\nabla^*p^*+\mu\frac{V}{L^2}\nabla^{*2}\mathbf v^*,\\
\frac{D\mathbf v^*}{Dt^*}&=-\nabla^*p^*+\frac{\mu}{\rho VL}\nabla^{*2}\mathbf v^*.
\end{aligned}
$$
$\therefore$ 自然出现雷诺数
$$\mathrm{Re}=\frac{\rho VL}{\mu}=\frac{VL}{\nu}.$$

诸准则及其含义：
$$\mathrm{Fr}=\frac{v}{\sqrt{gL}}\ (\text{惯性/重力}),\quad \mathrm{Ma}=\frac{v}{c}\ (\text{流速/声速}),$$
$$\mathrm{St}=\frac{fL}{v}\ (\text{非定常惯性/对流惯性}),\quad \mathrm{Eu}=\frac{\Delta p}{\rho v^2}\ (\text{压力/惯性}).$$

物理意义：准则相等是流动动力相似的充要条件。

极限与特例：$\mathrm{Re}\to0$ 为斯托克斯（蠕流）流；$\mathrm{Re}\to\infty$ 为无粘流；$\mathrm{Ma}<0.3$ 近似不可压。

典型应用：风洞缩比试验、水洞与飞行器模型设计。

**知识点6　相似三定理与部分相似**

定义：第一相似定理：相似系统同名相似准则相等；第二相似定理：物理关系可表为准则间函数关系；第三相似定理：单值条件相似且准则相等则系统完全相似。

推导思路：设原型 $p$ 与模型 $m$ 的量满足尺度变换
$$x_m=\lambda_L x_p,\quad v_m=\lambda_v v_p,\quad t_m=\lambda_t t_p.$$
运动相似要求 $\lambda_v=\lambda_L/\lambda_t$；动力相似要求力比 $\lambda_F=\lambda_\rho\lambda_L^2\lambda_v^2$。若所有准则相等，诸系统的无量纲方程 $F(\pi)=0$ 完全相同，$\therefore$ 无量纲解相同，即相似。

物理意义：模型试验结果可反演到原型。

极限与特例：同时满足 $\mathrm{Re}$ 与 $\mathrm{Fr}$ 需 $\lambda_\nu=\lambda_L^{3/2}$，通常难以实现，须采用部分相似、畸变模型或松弛准则。

典型应用：飞机缩比模型、桥梁气弹模型、船舶阻力试验。

### 第2章 矢量与张量代数及场论

**知识点1　指标记法与爱因斯坦求和约定**

定义：三维空间中矢量 $\mathbf a=a_i\mathbf e_i$（$i=1,2,3$）。约定同一项中重复出现两次的指标表示对该指标从 $1$ 到 $3$ 求和，称为哑标；只出现一次的指标为自由标。

假设：笛卡尔标准正交基满足 $\mathbf e_i\cdot\mathbf e_j=\delta_{ij}$，$\mathbf e_i\times\mathbf e_j=\epsilon_{ijk}\mathbf e_k$。

推导：
$$\mathbf a\cdot\mathbf b=(a_i\mathbf e_i)\cdot(b_j\mathbf e_j)=a_ib_j\delta_{ij}=a_ib_i,$$
$$(\mathbf a\times\mathbf b)_i=\epsilon_{ijk}a_jb_k,$$
$$\mathbf a\cdot(\mathbf b\times\mathbf c)=\epsilon_{ijk}a_ib_jc_k.$$
证 $\epsilon$–$\delta$ 恒等式：
$$\epsilon_{ijk}\epsilon_{ilm}=\delta_{jl}\delta_{km}-\delta_{jm}\delta_{kl}.$$
由之推导矢量三重积：
$$
\begin{aligned}
[\mathbf a\times(\mathbf b\times\mathbf c)]_i
&=\epsilon_{ijk}a_j(\mathbf b\times\mathbf c)_k
=\epsilon_{ijk}\epsilon_{klm}a_jb_lc_m\\
&=\epsilon_{ijk}\epsilon_{lmk}a_jb_lc_m
=(\delta_{il}\delta_{jm}-\delta_{im}\delta_{jl})a_jb_lc_m\\
&=a_jb_ic_j-a_jb_jc_i
=b_i(\mathbf a\cdot\mathbf c)-c_i(\mathbf a\cdot\mathbf b).
\end{aligned}
$$
$\therefore \mathbf a\times(\mathbf b\times\mathbf c)=\mathbf b(\mathbf a\cdot\mathbf c)-\mathbf c(\mathbf a\cdot\mathbf b)$。

物理意义：指标记法把矢量恒等式化为纯代数恒等式，便于机器推导与张量推广。

极限与特例：$\epsilon_{ijk}\epsilon_{ijl}=2\delta_{kl}$（对 $k,l$ 缩并）。

典型应用：推导刚体运动学、电磁学与流体力学矢量恒等式。

**知识点2　张量的定义与坐标变换**

定义：$n$ 阶张量 $T_{i_1\cdots i_n}$ 在正交变换 $x_i'=Q_{ij}x_j$ 下满足
$$T'_{i_1\cdots i_n}=Q_{i_1j_1}\cdots Q_{i_nj_n}T_{j_1\cdots j_n},$$
其中 $Q$ 为常正交矩阵。

推导 $Q$ 的正交性：长度不变要求
$$x_i'x_i'=Q_{ij}x_jQ_{ik}x_k=x_jx_k\delta_{jk},$$
对任意 $\mathbf x$ 成立，$\therefore Q_{ij}Q_{ik}=\delta_{jk}$，即 $Q^TQ=I$。

物理意义：张量是不依赖观察坐标系的客观量；标量是 0 阶、矢量是 1 阶、矩阵是 2 阶张量。

极限与特例：单位二阶张量 $\delta_{ij}$；置换张量 $\epsilon_{ijk}$ 为三阶赝张量（在镜像变换下多一个负号）。

典型应用：应力、应变、惯性张量、各向异性材料本构、坐标旋转。

**知识点3　张量的代数运算与对称分解**

定义：张量积 $(A\otimes B)_{ijkl}=A_{ij}B_{kl}$；缩并如 $C_{ik}=A_{ijk}$（对某对指标求和）；对称化与反对称化。

推导对称—反对称分解：设 $A_{ij}$ 任意。令
$$S_{ij}=\frac12(A_{ij}+A_{ji}),\qquad W_{ij}=\frac12(A_{ij}-A_{ji}),$$
则显然 $A_{ij}=S_{ij}+W_{ij}$，且 $S_{ij}=S_{ji}$，$W_{ij}=-W_{ji}$。唯一性：设另有 $A=S'+W'$，则转置 $A^T=S'-W'$，两式相加得 $S'=\frac12(A+A^T)=S$，相减得 $W'=W$。

二次型性质：对任意 $\mathbf x$
$$A_{ij}x_ix_j=\frac12(A_{ij}+A_{ji})x_ix_j=S_{ij}x_ix_j,$$
$\because W_{ij}x_ix_j=0$（交换哑标 $i,j$ 得 $W_{ij}x_ix_j=-W_{ij}x_ix_j$）。

物理意义：对称部分贡献变形/做功，反对称部分对应转动。实对称张量可正交对角化，主值实、主轴正交。

极限与特例：反对称二阶张量与轴矢量对应 $W_{ij}=-\epsilon_{ijk}\omega_k$，$\omega_k=-\frac12\epsilon_{kij}W_{ij}$。

典型应用：应力偏量/静水分量、应变/转动张量、涡量。

**知识点4　高斯散度定理**

定义：对光滑矢量场 $\mathbf F$ 与有界区域 $\Omega$（边界 $\partial\Omega$，外法线 $\mathbf n$），
$$\int_\Omega \nabla\cdot\mathbf F\,dV=\oint_{\partial\Omega}\mathbf F\cdot\mathbf n\,dS.$$

推导：先证一维基本情形
$$\int_a^b\frac{\partial F_x}{\partial x}\,dx=F_x(b)-F_x(a).$$
对长方体区域 $[0,L_x]\times[0,L_y]\times[0,L_z]$，
$$
\begin{aligned}
\int_\Omega\frac{\partial F_x}{\partial x}\,dV
&=\int_0^{L_y}\!\!\int_0^{L_z}\!\big[F_x(L_x)-F_x(0)\big]\,dy\,dz\\
&=\oint_{\partial\Omega}F_xn_x\,dS,
\end{aligned}
$$
因为在 $x=L_x$ 面 $n_x=1$、在 $x=0$ 面 $n_x=-1$。对 $F_y,F_z$ 同理相加：
$$\int_\Omega\Big(\frac{\partial F_x}{\partial x}+\frac{\partial F_y}{\partial y}+\frac{\partial F_z}{\partial z}\Big)dV=\oint_{\partial\Omega}\mathbf F\cdot\mathbf n\,dS.$$
对一般光滑区域，用 $C^1$ 坐标变换将小立方体拼接并取极限即可。

物理/几何意义：散度的体积分等于通量的面积分，是散度作为“源强度密度”的积分表述。

极限与特例：$\mathbf F=\mathbf r$，$\nabla\cdot\mathbf r=3$，得 $3V=\oint_{\partial\Omega}\mathbf r\cdot\mathbf n\,dS$。

典型应用：守恒律积分—微分形式转换、有限体积法、格林恒等式。

**知识点5　斯托克斯定理**

定义：对光滑曲面 $S$（边界为闭曲线 $\partial S$，$\mathbf n$ 为法向），
$$\int_S(\nabla\times\mathbf F)\cdot\mathbf n\,dS=\oint_{\partial S}\mathbf F\cdot d\mathbf l.$$

推导：先取 $xy$ 平面内区域 $S$，$\mathbf n=\mathbf e_z$，$\mathbf F=(F_x,F_y,0)$。则
$$(\nabla\times\mathbf F)\cdot\mathbf e_z=\frac{\partial F_y}{\partial x}-\frac{\partial F_x}{\partial y}.$$
由格林定理
$$
\begin{aligned}
\int_S(\nabla\times\mathbf F)\cdot\mathbf n\,dS
&=\int_S\Big(\frac{\partial F_y}{\partial x}-\frac{\partial F_x}{\partial y}\Big)dx\,dy\\
&=\oint_{\partial S}(F_x\,dx+F_y\,dy)
=\oint_{\partial S}\mathbf F\cdot d\mathbf l.
\end{aligned}
$$
对一般空间曲面，将其分割为小片，每片法向近似常值，相邻内部边界方向相反积分相消，取极限得证。

物理意义：旋度的通量等于沿边界的环量。

极限与特例：$\nabla\times\mathbf F=0$ 的场环量为零，$\mathbf F=\nabla\phi$（保守场）。

典型应用：开尔文环量定理、安培环路定理、斯托克斯流函数。

**知识点6　正交曲线坐标中的微分算子**

定义：设位置 $\mathbf r=\mathbf r(q_1,q_2,q_3)$，拉梅系数 $h_i=|\partial\mathbf r/\partial q_i|$，单位基 $\mathbf e_i=\dfrac1{h_i}\dfrac{\partial\mathbf r}{\partial q_i}$。

推导梯度：由全微分
$$d\phi=\sum_i\frac{\partial\phi}{\partial q_i}dq_i=\nabla\phi\cdot d\mathbf r,\qquad d\mathbf r=\sum_i h_i\,dq_i\,\mathbf e_i,$$
比较各分量得
$$\nabla\phi=\sum_i\frac1{h_i}\frac{\partial\phi}{\partial q_i}\mathbf e_i.$$

推导散度：考虑微元体积 $dV=h_1h_2h_3\,dq_1dq_2dq_3$，由高斯定理
$$\nabla\cdot\mathbf F=\lim_{V\to0}\frac1V\oint_{\partial V}\mathbf F\cdot\mathbf n\,dS.$$
沿 $q_1$ 方向两对面净通量为
$$\frac{\partial}{\partial q_1}(F_1h_2h_3)\,dq_1dq_2dq_3,$$
依次相加并除以 $dV$：
$$\nabla\cdot\mathbf F=\frac1{h_1h_2h_3}\left[\frac{\partial(F_1h_2h_3)}{\partial q_1}+\frac{\partial(F_2h_3h_1)}{\partial q_2}+\frac{\partial(F_3h_1h_2)}{\partial q_3}\right].$$

物理/几何意义：在贴合问题对称性的坐标下，方程形式最简。

极限与特例：柱坐标 $(r,\theta,z)$：$h_r=1,\ h_\theta=r,\ h_z=1$；球坐标 $(r,\theta,\varphi)$：$h_r=1,\ h_\theta=r,\ h_\varphi=r\sin\theta$。球坐标拉普拉斯算子
$$\nabla^2\phi=\frac1{r^2}\frac{\partial}{\partial r}\!\Big(r^2\frac{\partial\phi}{\partial r}\Big)+\frac1{r^2\sin\theta}\frac{\partial}{\partial\theta}\!\Big(\sin\theta\frac{\partial\phi}{\partial\theta}\Big)+\frac1{r^2\sin^2\theta}\frac{\partial^2\phi}{\partial\varphi^2}.$$

典型应用：圆管流、球对称场、圆柱壳与球壳问题。

### 第3章 质点与质点系动力学

**知识点1　牛顿运动定律与惯性系**

定义：存在惯性参考系，质点在其中不受力时保持静止或匀速直线运动；受力时
$$\mathbf F=\frac{d\mathbf p}{dt}=\frac{d}{dt}(m\mathbf v).$$

推导：当质量 $m$ 为常数，
$$\mathbf F=\frac{d(m\mathbf v)}{dt}=m\frac{d\mathbf v}{dt}=m\mathbf a.$$
伽利略变换：设 $S'$ 相对 $S$ 以常速 $\mathbf V$ 平动，$x'=x-Vt,\ t'=t$。则
$$\mathbf v'=\frac{d\mathbf x'}{dt'}=\mathbf v-\mathbf V,\qquad \mathbf a'=\frac{d\mathbf v'}{dt'}=\mathbf a.$$
$\therefore$ 加速度在惯性系间不变，$\mathbf F=m\mathbf a$ 形式不变。

物理意义：力是改变动量（而非维持运动）的原因；加速度是绝对量。

极限与特例：$m$ 随速度变化时须用 $\mathbf F=d\mathbf p/dt$；相对论下 $\mathbf p=\gamma m_0\mathbf v$。

典型应用：轨道运动、刚体平动、实验班力学建模的出发点。

**知识点2　质点动量定理与冲量**

定义：冲量 $\mathbf J=\displaystyle\int_{t_1}^{t_2}\mathbf F\,dt$，动量 $\mathbf p=m\mathbf v$。

推导：由牛顿第二定律
$$\mathbf F=\frac{d\mathbf p}{dt}\ \Rightarrow\ \int_{t_1}^{t_2}\mathbf F\,dt=\int_{\mathbf p_1}^{\mathbf p_2}d\mathbf p=\mathbf p_2-\mathbf p_1.$$
$\therefore \mathbf J=m\mathbf v_2-m\mathbf v_1$。

物理意义：冲量等于动量增量；力—时间曲线下的面积决定动量变化。

极限与特例：碰撞时间趋零时力趋于 $\delta$ 函数，冲量有限；动量守恒当 $\mathbf F=0$。

典型应用：碰撞、爆炸、缓冲设计、着陆冲击。

**知识点3　质点的角动量定理**

定义：对固定点 $O$，角动量 $\mathbf L=\mathbf r\times m\mathbf v$，力矩 $\mathbf M=\mathbf r\times\mathbf F$。

推导：
$$\frac{d\mathbf L}{dt}=\frac{d\mathbf r}{dt}\times m\mathbf v+\mathbf r\times\frac{d(m\mathbf v)}{dt}.$$
$\because \dfrac{d\mathbf r}{dt}=\mathbf v\parallel\mathbf v$，$\therefore \mathbf v\times m\mathbf v=0$。故
$$\frac{d\mathbf L}{dt}=\mathbf r\times\mathbf F=\mathbf M.$$

物理/几何意义：力矩是角动量变化率；中心力场中 $\mathbf M=0$，角动量守恒，运动在过原点的平面内。

极限与特例：开普勒第二定律即 $\mathbf L$ 守恒的几何表述（面积速度 $=\frac{L}{2m}$）。

典型应用：行星轨道、陀螺进动、有心力问题。

**知识点4　功与动能定理**

定义：功 $W=\displaystyle\int_{\mathbf r_1}^{\mathbf r_2}\mathbf F\cdot d\mathbf r$，动能 $T=\tfrac12mv^2$。

推导：由牛顿定律
$$
\begin{aligned}
W&=\int \mathbf F\cdot d\mathbf r=\int m\frac{d\mathbf v}{dt}\cdot\mathbf v\,dt
=\int m\,\mathbf v\cdot d\mathbf v\\
&=m\int v\,dv=\frac12mv_2^2-\frac12mv_1^2=T_2-T_1.
\end{aligned}
$$
$\therefore W=\Delta T$。

物理意义：合外力做功等于动能增量；切向力做功，法向力不做功。

极限与特例：$W=0$ 时速率不变（约束力、洛伦兹力）；保守力做功与路径无关。

典型应用：能量法求速度、保守系统快速求解。

**知识点5　势能与机械能守恒**

定义：保守力满足 $\oint\mathbf F\cdot d\mathbf r=0$，可引入势能 $V$ 使
$$\mathbf F=-\nabla V.$$

推导：保守力 $\Rightarrow$ 存在 $V(\mathbf r)$ 使得 $dV=-\mathbf F\cdot d\mathbf r$。由动能定理
$$dT=\mathbf F\cdot d\mathbf r=-dV\ \Rightarrow\ d(T+V)=0,$$
$$\therefore\ E=T+V=\text{const}.$$
验证负梯度：在笛卡尔坐标
$$\mathbf F\cdot d\mathbf r=F_xdx+F_ydy+F_zdz=-\Big(\frac{\partial V}{\partial x}dx+\frac{\partial V}{\partial y}dy+\frac{\partial V}{\partial z}dz\Big),$$
比较分量即得 $\mathbf F=-\nabla V$。

物理/几何意义：$V$ 的等值面为等势面，$\mathbf F$ 垂直于等势面指向 $V$ 下降方向。

极限与特例：重力 $V=mgz$；弹簧 $V=\tfrac12kx^2$；有心力 $V(r)=-\dfrac{GMm}{r}$。

典型应用：单摆、抛体、守恒量快速求解与稳定性分析。

**知识点6　质点系动力学一般定理**

定义：质点系质量 $M=\sum_im_i$，质心 $\mathbf r_C=\dfrac1M\sum_im_i\mathbf r_i$。

推导质心运动定理：对每个质点 $\mathbf F_i+\mathbf F_i^{\text{int}}=m_i\mathbf a_i$，求和
$$\sum_i\mathbf F_i^{\text{ext}}+\sum_i\mathbf F_i^{\text{int}}=\sum_im_i\mathbf a_i.$$
$\because$ 内力的矢量和为零（牛顿第三定律成对抵消），且 $\sum_im_i\mathbf a_i=M\mathbf a_C$，$\therefore$
$$\mathbf F^{\text{ext}}=M\mathbf a_C.$$
对内力的矩同样成对抵消（以及力共线），得
$$\frac{d\mathbf L_C}{dt}=\mathbf M^{\text{ext}}.$$
柯尼希定理：设 $\mathbf v_i=\mathbf v_C+\mathbf v_i'$，
$$
\begin{aligned}
T&=\frac12\sum_im_i(\mathbf v_C+\mathbf v_i')\cdot(\mathbf v_C+\mathbf v_i')\\
&=\frac12Mv_C^2+\mathbf v_C\cdot\frac{d}{dt}\sum_im_i\mathbf r_i'+T'\\
&=\frac12Mv_C^2+T',
\end{aligned}
$$
因 $\sum m_i\mathbf r_i'=0$ 其导数为零。

物理/几何意义：质点系运动可分解为质心平动与相对运动。

极限与特例：刚体即 $T'=\tfrac12 I\omega^2$ 的质点系；碰撞中 $\mathbf F^{\text{ext}}=0$ 时质心动量守恒。

典型应用：火箭、爆轰、多体系统、刚体动力学。

**知识点7　变质量系统与火箭方程**

定义：主体质量为 $m(t)$，以相对速度 $\mathbf u$（喷射物相对主体的速度）排出物质。

推导：取时刻 $t$ 系统（主体 $+$ 即将排出质量）用动量守恒。时刻 $t$ 动量 $m\mathbf v$；时刻 $t+dt$ 主体动量 $(m+dm)(\mathbf v+d\mathbf v)$，排出物动量 $(-dm)(\mathbf v+\mathbf u)$，其中 $dm<0$。外力冲量 $\mathbf F^{\text{ext}}dt$。由动量定理
$$\mathbf F^{\text{ext}}dt=(m+dm)(\mathbf v+d\mathbf v)+(-dm)(\mathbf v+\mathbf u)-m\mathbf v.$$
展开并略去二阶小量 $dm\,d\mathbf v$：
$$\mathbf F^{\text{ext}}dt=m\,d\mathbf v+dm\,(\mathbf u-\mathbf v)+dm\,\mathbf v=m\,d\mathbf v+dm\,\mathbf u.$$
$\therefore$
$$m\frac{d\mathbf v}{dt}=\mathbf F^{\text{ext}}-\mathbf u\frac{dm}{dt}.$$
无外力、沿 $\mathbf u$ 方向积分得齐奥尔科夫斯基公式
$$\Delta v=\int_{m_0}^{m_f}\frac{-u\,dm}{m}=u\ln\frac{m_0}{m_f}.$$

物理意义：喷出物带来反冲推力 $-\mathbf u\,dm/dt$。

极限与特例：$m_f/m_0\to1$ 时 $\Delta v\to0$；多级火箭通过逐级抛离提高 $\Delta v$。

典型应用：运载火箭、变质量落体、喷气推进。

### 第4章 拉格朗日力学

**知识点1　广义坐标与完整约束**

定义：用 $s$ 个广义坐标 $q_1,\dots,q_s$ 描述系统位形。完整约束为可积的位形约束
$$f_\alpha(\mathbf r_1,\dots,\mathbf r_N,t)=0.$$
自由度数 $s=3N-c$（$N$ 个质点、$c$ 个独立完整约束）。

推导：约束将 $3N$ 维位形空间限制为 $3N-c$ 维流形，可用独立参数 $q_j$ 参数化：
$$\mathbf r_i=\mathbf r_i(q_1,\dots,q_s,t).$$
于是
$$\delta\mathbf r_i=\sum_{j=1}^{s}\frac{\partial\mathbf r_i}{\partial q_j}\delta q_j,\qquad
\dot{\mathbf r}_i=\sum_j\frac{\partial\mathbf r_i}{\partial q_j}\dot q_j+\frac{\partial\mathbf r_i}{\partial t}.$$

物理/几何意义：广义坐标自动消去约束，无需引入约束力。

极限与特例：非完整约束（如滚动）不可积，须用乘子或拉格朗日—达朗贝尔方程。

典型应用：多连杆、机构、刚体、机器人运动学建模。

**知识点2　虚功原理与理想约束**

定义：虚位移 $\delta\mathbf r_i$ 为满足约束的无穷小、瞬间的坐标变分。理想约束指约束力做的虚功之和为零
$$\sum_i\mathbf N_i\cdot\delta\mathbf r_i=0.$$

推导：对静平衡系统，各质点 $\mathbf F_i+\mathbf N_i=0$，取虚位移点乘并求和
$$\sum_i(\mathbf F_i+\mathbf N_i)\cdot\delta\mathbf r_i=0.$$
由理想约束条件，$\sum_i\mathbf N_i\cdot\delta\mathbf r_i=0$，$\therefore$
$$\sum_i\mathbf F_i\cdot\delta\mathbf r_i=0.$$
用广义力 $Q_j=\sum_i\mathbf F_i\cdot\dfrac{\partial\mathbf r_i}{\partial q_j}$，得
$$\sum_{j=1}^{s}Q_j\,\delta q_j=0.$$
$\because \delta q_j$ 独立，$\therefore Q_j=0$。

物理意义：理想约束力不做功，可从平衡方程中消去。

极限与特例：刚体、无滑动接触、理想铰链均为理想约束。

典型应用：静力学平衡、结构分析。

**知识点3　达朗贝尔原理与拉格朗日方程**

定义：将动力学问题化为“惯性力”参与下的静力平衡，即对每个质点
$$\mathbf F_i-m_i\mathbf a_i=0.$$

推导：由虚功原理作用于动力学
$$\sum_i(\mathbf F_i-m_i\mathbf a_i)\cdot\delta\mathbf r_i=0.$$
将 $\delta\mathbf r_i=\sum_j\dfrac{\partial\mathbf r_i}{\partial q_j}\delta q_j$ 代入并交换求和：
$$\sum_j\Big[\sum_i\mathbf F_i\cdot\frac{\partial\mathbf r_i}{\partial q_j}-\sum_i m_i\mathbf a_i\cdot\frac{\partial\mathbf r_i}{\partial q_j}\Big]\delta q_j=0.$$
$\because \delta q_j$ 独立，括号为零。记广义力 $Q_j=\sum_i\mathbf F_i\cdot\dfrac{\partial\mathbf r_i}{\partial q_j}$。对惯性项利用
$$\frac{\partial\mathbf r_i}{\partial q_j}=\frac{\partial\dot{\mathbf r}_i}{\partial\dot q_j},\qquad
\frac{d}{dt}\frac{\partial\mathbf r_i}{\partial q_j}=\frac{\partial\dot{\mathbf r}_i}{\partial q_j},$$
得
$$
\begin{aligned}
\sum_i m_i\mathbf a_i\cdot\frac{\partial\mathbf r_i}{\partial q_j}
&=\sum_i m_i\dot{\mathbf v}_i\cdot\frac{\partial\mathbf v_i}{\partial\dot q_j}
=\frac{d}{dt}\sum_i m_i\mathbf v_i\cdot\frac{\partial\mathbf v_i}{\partial\dot q_j}-\sum_i m_i\mathbf v_i\cdot\frac{\partial\mathbf v_i}{\partial q_j}\\
&=\frac{d}{dt}\frac{\partial T}{\partial\dot q_j}-\frac{\partial T}{\partial q_j}.
\end{aligned}
$$
$\therefore$
$$\frac{d}{dt}\frac{\partial T}{\partial\dot q_j}-\frac{\partial T}{\partial q_j}=Q_j.$$
若有势能 $V$ 且 $Q_j=-\dfrac{\partial V}{\partial q_j}$，令 $L=T-V$，即得拉格朗日方程
$$\boxed{\ \frac{d}{dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}=0\ },\qquad j=1,\dots,s.$$

物理/几何意义：运动方程由动能与势能构造，自动满足约束，形式与坐标选取无关。

极限与特例：$V$ 不显含速度时 $L=T-V$；含速度依赖势（电磁）时须用广义势。

典型应用：复杂约束系统、机器人、航天器建模。

**知识点4　拉格朗日函数的循环积分与能量守恒**

定义：若 $\dfrac{\partial L}{\partial q_j}=0$，称 $q_j$ 为循环坐标，其共轭动量
$$p_j=\frac{\partial L}{\partial\dot q_j}.$$

推导（循环积分）：拉格朗日方程
$$\frac{d}{dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}=0.$$
当 $\dfrac{\partial L}{\partial q_j}=0$，得 $\dfrac{d}{dt}p_j=0$，$\therefore p_j=\text{const}$。

推导（能量守恒）：当 $\dfrac{\partial L}{\partial t}=0$，求 $L$ 的全时间导数
$$\frac{dL}{dt}=\sum_j\Big(\frac{\partial L}{\partial q_j}\dot q_j+\frac{\partial L}{\partial\dot q_j}\ddot q_j\Big).$$
由拉格朗日方程 $\dfrac{\partial L}{\partial q_j}=\dfrac{d}{dt}p_j$，代入：
$$
\begin{aligned}
\frac{dL}{dt}&=\sum_j\Big[\dot q_j\frac{dp_j}{dt}+p_j\ddot q_j\Big]
=\frac{d}{dt}\sum_j p_j\dot q_j,\\
\Rightarrow\ \frac{d}{dt}\Big(\sum_j p_j\dot q_j-L\Big)&=0.
\end{aligned}
$$
定义广义能量（雅可比积分）
$$E=\sum_j p_j\dot q_j-L=\text{const}.$$
当 $T=\tfrac12\sum M_{ij}\dot q_i\dot q_j$ 为速度二次齐次、$V=V(q)$ 时，由欧拉齐次定理 $\sum_j p_j\dot q_j=2T$，$\therefore E=T+V$。

物理/几何意义：时间平移对称 $\Rightarrow$ 能量守恒；坐标平移/旋转对称 $\Rightarrow$ 动量/角动量守恒。

极限与特例：非惯性系中 $L$ 可能显含时间，机械能不守恒。

典型应用：中心力、陀螺、守恒量降维。

**知识点5　诺特定理（连续对称与守恒律）**

定义：若作用量 $S=\int L\,dt$ 在单参数连续变换 $q_j\to q_j+\epsilon\,\eta_j(q,t)$ 下不变，则存在守恒量。

推导：设变换的生成元 $\eta_j$。作用量变分
$$\delta S=\int\Big[\sum_j\Big(\frac{\partial L}{\partial q_j}\eta_j+\frac{\partial L}{\partial\dot q_j}\dot\eta_j\Big)\Big]dt.$$
在运动方程成立时，对时间积分项分部积分：
$$
\begin{aligned}
\delta S&=\int\sum_j\Big[\frac{\partial L}{\partial q_j}-\frac{d}{dt}\frac{\partial L}{\partial\dot q_j}\Big]\eta_j\,dt
+\Big[\sum_j p_j\eta_j\Big]_{t_1}^{t_2}.
\end{aligned}
$$
若变换为对称，$\delta S=0$ 对任意端点及运动成立，故
$$\sum_j p_j\eta_j=\text{const}.$$

对应关系：时间平移 $\eta_t=1$ $\Rightarrow$ 能量 $E$；空间平移 $\eta_{\mathbf r}=\mathbf a$ $\Rightarrow$ 动量 $\mathbf p$；旋转 $\boldsymbol\eta=\boldsymbol\omega\times\mathbf r$ $\Rightarrow$ 角动量 $\mathbf L$。

物理/几何意义：守恒律源于系统的连续对称性，是分析力学最深刻的结构。

极限与特例：显式破缺（外力、阻尼）时对应量不守恒。

典型应用：从对称性快速得到守恒量，指导降维与数值守恒格式。

**知识点6　小振动与线性化**

定义：在稳定平衡位形 $q_j=q_j^0$ 附近作小扰动，线性化拉格朗日方程得到
$$M\ddot{\mathbf q}+K\mathbf q=0.$$

推导：平衡点满足 $\dfrac{\partial V}{\partial q_j}\big|_{q^0}=0$。将 $V,M$ 在 $q^0$ 处展开
$$
\begin{aligned}
V(q^0+\mathbf q)&=V(q^0)+\frac12\sum_{i,j}K_{ij}q_iq_j+\cdots,\quad
K_{ij}=\frac{\partial^2 V}{\partial q_i\partial q_j}\Big|_{q^0},\\
T&=\frac12\sum_{i,j}M_{ij}\dot q_i\dot q_j+\cdots,\quad
M_{ij}=M_{ji}>0.
\end{aligned}
$$
代入 $L=T-V$ 并保留二次项，拉格朗日方程给出
$$\sum_j M_{ij}\ddot q_j+\sum_j K_{ij}q_j=0.$$
设简谐解 $q_j=A_je^{i\omega t}$，得
$$(K-\omega^2M)\mathbf A=0.$$
非平凡解存在要求
$$\det(K-\omega^2M)=0.$$
$\therefore$ 得到 $s$ 个固有频率 $\omega_k$ 与振型 $\mathbf A^{(k)}$，且 $\mathbf A^{(i)T}M\mathbf A^{(j)}=0$（$i\neq j$），$\mathbf A^{(i)T}K\mathbf A^{(j)}=0$。

物理/几何意义：小振动是二次型广义特征值问题；$M$ 正定、$K$ 正定对应稳定平衡。

极限与特例：$K$ 非正定 $\Rightarrow$ 不稳定（频率虚数）；单自由度 $\omega=\sqrt{K/M}$。

典型应用：结构模态分析、分子振动、多体系统线性稳定性。

**知识点7　拉格朗日乘子与约束力**

定义：对不便于消去的约束，引入乘子 $\lambda_\alpha$。

推导：设约束 $f_\alpha(q,t)=0$。对广义坐标变分
$$\sum_j\frac{\partial f_\alpha}{\partial q_j}\delta q_j=0.$$
将约束力做的虚功写成 $\sum_\alpha\lambda_\alpha\dfrac{\partial f_\alpha}{\partial q_j}\delta q_j$。达朗贝尔原理变为
$$\sum_j\Big[\frac{d}{dt}\frac{\partial T}{\partial\dot q_j}-\frac{\partial T}{\partial q_j}-Q_j-\sum_\alpha\lambda_\alpha\frac{\partial f_\alpha}{\partial q_j}\Big]\delta q_j=0.$$
此时 $\delta q_j$ 不再独立，但 $\lambda_\alpha$ 的选取可抵消约束方向的组合，最终
$$\frac{d}{dt}\frac{\partial T}{\partial\dot q_j}-\frac{\partial T}{\partial q_j}=Q_j+\sum_\alpha\lambda_\alpha\frac{\partial f_\alpha}{\partial q_j},\qquad f_\alpha=0.$$
未知量为 $q_j$ 与 $\lambda_\alpha$，共 $s+c$ 个，方程数亦相同。

物理/几何意义：$\lambda_\alpha$ 即约束力（如绳张力、铰链反力、接触力）。

极限与特例：完整约束可消去，乘子法便于求约束力。

典型应用：机构动力学、接触/滚动问题、多体系统仿真。

### 第5章 哈密顿力学与正则方程

**知识点1　勒让德变换与哈密顿量**

定义：由拉格朗日函数 $L(q,\dot q,t)$ 定义共轭动量 $p_j=\dfrac{\partial L}{\partial\dot q_j}$，哈密顿量
$$H(q,p,t)=\sum_j p_j\dot q_j-L.$$

推导：由 $p_j=\partial L/\partial\dot q_j$，若 $\partial^2L/\partial\dot q_i\partial\dot q_j$ 非奇异，可反解 $\dot q_j=\dot q_j(q,p,t)$。对 $H$ 求微分：
$$
dH=\sum_j\Big(p_j\,d\dot q_j+\dot q_j\,dp_j-\frac{\partial L}{\partial q_j}dq_j-\frac{\partial L}{\partial\dot q_j}d\dot q_j\Big)-\frac{\partial L}{\partial t}dt.
$$
$\because p_j=\partial L/\partial\dot q_j$，$\therefore$ $d\dot q_j$ 项相消：
$$
dH=\sum_j\Big(\dot q_j\,dp_j-\frac{\partial L}{\partial q_j}dq_j\Big)-\frac{\partial L}{\partial t}dt.
$$
另一方面 $dH=\sum_j\big(\partial H/\partial p_j\,dp_j+\partial H/\partial q_j\,dq_j\big)+\partial H/\partial t\,dt$，比较得
$$\dot q_j=\frac{\partial H}{\partial p_j},\qquad \frac{\partial H}{\partial q_j}=-\frac{\partial L}{\partial q_j}.$$

物理/几何意义：勒让德变换把速度变量换成动量变量，$H$ 与 $L$ 互为对偶，构成 $(q,p)$ 相空间上的函数。

极限与特例：$L=T-V$ 且 $V=V(q)$ 时 $H=T+V=E$；含速度相关势时 $H\neq T+V$。

典型应用：保守系统、相空间分析、辛几何与几何数值积分。

**知识点2　哈密顿正则方程**

定义：在相空间 $(q,p)$ 中，运动由一阶方程组描述。

推导：由拉格朗日方程 $\dfrac{d}{dt}\dfrac{\partial L}{\partial\dot q_j}=\dfrac{\partial L}{\partial q_j}$，注意 $p_j=\partial L/\partial\dot q_j$，得
$$\dot p_j=\frac{\partial L}{\partial q_j}.$$
又由勒让德变换结果 $\dfrac{\partial H}{\partial q_j}=-\dfrac{\partial L}{\partial q_j}$，故
$$\dot p_j=-\frac{\partial H}{\partial q_j}.$$
结合 $\dot q_j=\dfrac{\partial H}{\partial p_j}$，得正则方程
$$\boxed{\ \dot q_j=\frac{\partial H}{\partial p_j},\qquad \dot p_j=-\frac{\partial H}{\partial q_j}\ }.$$

能量演化：
$$\frac{dH}{dt}=\sum_j\Big(\frac{\partial H}{\partial q_j}\dot q_j+\frac{\partial H}{\partial p_j}\dot p_j\Big)+\frac{\partial H}{\partial t}=\frac{\partial H}{\partial t}.$$
$\therefore$ 若 $H$ 不显含时间，则 $H=\text{const}$。

物理/几何意义：正则方程是相空间中的哈密顿向量场，其流保持辛形式 $\omega=\sum_j dq_j\wedge dp_j$。

极限与特例：$H=\dfrac{p^2}{2m}+V(q)$ 时 $\dot q=p/m$，$\dot p=-V'(q)$。

典型应用：轨道摄动、粒子加速器、辛积分器设计。

**知识点3　泊松括号与守恒量**

定义：对相空间函数 $f(q,p,t),g(q,p,t)$，
$$\{f,g\}=\sum_j\Big(\frac{\partial f}{\partial q_j}\frac{\partial g}{\partial p_j}-\frac{\partial f}{\partial p_j}\frac{\partial g}{\partial q_j}\Big).$$

推导运动方程：
$$
\begin{aligned}
\frac{df}{dt}&=\sum_j\Big(\frac{\partial f}{\partial q_j}\dot q_j+\frac{\partial f}{\partial p_j}\dot p_j\Big)+\frac{\partial f}{\partial t}\\
&=\sum_j\Big(\frac{\partial f}{\partial q_j}\frac{\partial H}{\partial p_j}-\frac{\partial f}{\partial p_j}\frac{\partial H}{\partial q_j}\Big)+\frac{\partial f}{\partial t}\\
&=\{f,H\}+\frac{\partial f}{\partial t}.
\end{aligned}
$$
$\therefore$ $f$ 为守恒量当且仅当 $\{f,H\}+\partial f/\partial t=0$。

泊松括号性质：反对称 $\{f,g\}=-\{g,f\}$；双线性；莱布尼茨律 $\{f,gh\}=g\{f,h\}+\{f,g\}h$；雅可比恒等式
$$\{f,\{g,h\}\}+\{g,\{h,f\}\}+\{h,\{f,g\}\}=0.$$
基本括号 $\{q_i,q_j\}=0,\{p_i,p_j\}=0,\{q_i,p_j\}=\delta_{ij}$。

物理/几何意义：泊松括号对应量子对易子的经典极限 $\{,\}\to\frac{1}{i\hbar}[\,] $。

极限与特例：$\{H,H\}=0$ 即能量守恒。

典型应用：寻找守恒量、正则微扰、约化系统。

**知识点4　正则变换**

定义：保持正则方程形式不变的变换 $(q,p)\to(Q,P)$，即存在新哈密顿 $K(Q,P,t)$ 使 $\dot Q_j=\partial K/\partial P_j,\ \dot P_j=-\partial K/\partial Q_j$。

推导：正则方程来自变分原理 $\delta\int(\sum p_j\dot q_j-H)dt=0$。要求新变量同形式 $\delta\int(\sum P_j\dot Q_j-K)dt=0$，则两被积函数之差必须是某个 $F$ 的全微分：
$$\sum_j p_j\dot q_j-H=\sum_j P_j\dot Q_j-K+\frac{dF}{dt}.$$
取 $F=F_1(q,Q,t)$，则
$$\frac{dF_1}{dt}=\sum_j\frac{\partial F_1}{\partial q_j}\dot q_j+\sum_j\frac{\partial F_1}{\partial Q_j}\dot Q_j+\frac{\partial F_1}{\partial t}.$$
比较 $\dot q_j,\dot Q_j$ 系数得
$$p_j=\frac{\partial F_1}{\partial q_j},\qquad P_j=-\frac{\partial F_1}{\partial Q_j},\qquad K=H+\frac{\partial F_1}{\partial t}.$$
由第一式反解 $Q=Q(q,p,t)$ 并代入即得变换。

物理/几何意义：正则变换是相空间保辛同胚，可用来简化 $H$（如作用量—角变量）。

极限与特例：$F_1=\sum q_jP_j$ 给恒等变换；谐振子用作用量—角变量化为常量 $H$。

典型应用：摄动理论、可积系统、哈密顿—雅可比方法。

**知识点5　哈密顿—雅可比方程**

定义：寻求正则变换使新哈密顿 $K=0$，从而 $Q,P$ 为常数。由 $P_j=\partial S/\partial q_j$ 与 $K=H+\partial S/\partial t=0$ 得
$$H\Big(q,\frac{\partial S}{\partial q},t\Big)+\frac{\partial S}{\partial t}=0.$$

推导：取第二类生成函数 $S(q,P,t)$，则 $p_j=\partial S/\partial q_j$。要求 $K=H+\partial S/\partial t=0$，即得一阶偏微分方程。若 $H$ 不显含时间，可分离
$$S=W(q)-Et,\qquad H\Big(q,\frac{\partial W}{\partial q}\Big)=E.$$
解得 $W(q,\alpha)$ 后，$\beta_j=\partial W/\partial\alpha_j$ 为常数，反解 $q(t)$。

物理/几何意义：$S$ 是作用量，其等值面为波前，特征线即粒子轨道，是经典力学与波动力学的桥梁。

极限与特例：单自由度，$W=\int\sqrt{2m(E-V)}\,dq$。

典型应用：可积系统求解析解、几何光学类比。

**知识点6　刘维尔定理与相体积守恒**

定义：哈密顿系统的相空间流保持相体积（更一般地保持辛体积）不变。

推导：相空间中概率密度 $\rho(q,p,t)$ 满足连续性方程
$$\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\,\mathbf V)=0,$$
其中 $\mathbf V=(\dot q,\dot p)$。利用正则方程
$$\nabla\cdot\mathbf V=\sum_j\Big(\frac{\partial\dot q_j}{\partial q_j}+\frac{\partial\dot p_j}{\partial p_j}\Big)=\sum_j\Big(\frac{\partial^2H}{\partial q_j\partial p_j}-\frac{\partial^2H}{\partial p_j\partial q_j}\Big)=0.$$
$\therefore$ 连续性方程化为
$$\frac{\partial\rho}{\partial t}+\sum_j\Big(\dot q_j\frac{\partial\rho}{\partial q_j}+\dot p_j\frac{\partial\rho}{\partial p_j}\Big)=0,\qquad \frac{d\rho}{dt}=0.$$
即相流不可压缩，相体积守恒。

物理/几何意义：相空间“流体”保持体积，是统计力学系综理论的基础。

极限与特例：无耗散系统成立；有阻尼时 $\nabla\cdot\mathbf V\neq0$，相体积收缩。

典型应用：系综平均、混沌与遍历性、辛数值格式的守恒性。

**知识点7　最小作用量原理**

定义：物理系统的真实运动使作用量
$$S=\int_{t_1}^{t_2}L(q,\dot q,t)\,dt$$
取驻值（通常为极小）。

推导：设 $q(t)\to q(t)+\delta q(t)$，端点固定 $\delta q(t_1)=\delta q(t_2)=0$。
$$
\begin{aligned}
\delta S&=\int_{t_1}^{t_2}\sum_j\Big(\frac{\partial L}{\partial q_j}\delta q_j+\frac{\partial L}{\partial\dot q_j}\delta\dot q_j\Big)dt\\
&=\int_{t_1}^{t_2}\sum_j\frac{\partial L}{\partial q_j}\delta q_j\,dt+\sum_j\Big[\frac{\partial L}{\partial\dot q_j}\delta q_j\Big]_{t_1}^{t_2}-\int_{t_1}^{t_2}\sum_j\frac{d}{dt}\Big(\frac{\partial L}{\partial\dot q_j}\Big)\delta q_j\,dt\\
&=\int_{t_1}^{t_2}\sum_j\Big[\frac{\partial L}{\partial q_j}-\frac{d}{dt}\frac{\partial L}{\partial\dot q_j}\Big]\delta q_j\,dt.
\end{aligned}
$$
$\because \delta q_j$ 任意，$\therefore$ 由 $\delta S=0$ 得欧拉—拉格朗日方程
$$\frac{d}{dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}=0.$$

物理/几何意义：真实路径是位形空间中作用量泛函的驻值路径，体现自然界的变分原理。

极限与特例：短时高动量路径作用量较大，故通常为极小。

典型应用：场的拉格朗日密度、相对论与量子力学的路径积分。

### 第6章 刚体动力学

**知识点1　刚体运动学与角速度**

定义：刚体上任意两点距离不变。速度场
$$\mathbf v(\mathbf r)=\mathbf v_O+\boldsymbol\omega\times(\mathbf r-\mathbf r_O).$$

推导：由 $|\mathbf r_i-\mathbf r_j|^2=\text{const}$ 微分得 $(\mathbf r_i-\mathbf r_j)\cdot(\mathbf v_i-\mathbf v_j)=0$。设固定点 $O$，一般解为 $\mathbf v=\boldsymbol\omega\times\mathbf r$（$\boldsymbol\omega$ 可待定，反对称张量对应轴矢量）：
$$W_{ij}=\frac{\partial v_i}{\partial x_j}\ \text{反对称}\ \Rightarrow\ W_{ij}=-\epsilon_{ijk}\omega_k,\quad v_i=W_{ij}x_j=\epsilon_{ijk}\omega_jx_k=(\boldsymbol\omega\times\mathbf r)_i.$$
$\therefore$ 成立。

几何意义：瞬时角速度轴，速度由转动分量决定。

极限与特例：纯平动 $\boldsymbol\omega=0$；对任意参考点 $\mathbf v_P=\mathbf v_O+\boldsymbol\omega\times\mathbf{OP}$。

典型应用：刚体位姿运动学、姿态运动学建模。

**知识点2　惯性张量与角动量**

定义：对刚体质量分布，惯性张量
$$I_{ij}=\int_V\rho(\mathbf r)\big(\delta_{ij}r^2-r_ir_j\big)dV,$$
角动量 $\mathbf L=\int\mathbf r\times(\boldsymbol\omega\times\mathbf r)\,dm$，动能 $T=\tfrac12\boldsymbol\omega\cdot\mathbf L$。

推导：
$$
\begin{aligned}
\mathbf L&=\int\mathbf r\times(\boldsymbol\omega\times\mathbf r)\,dm
=\int\big[r^2\boldsymbol\omega-(\boldsymbol\omega\cdot\mathbf r)\mathbf r\big]dm\\
&=\int\big[\delta_{ij}r^2\boldsymbol\omega_j-r_ir_j\boldsymbol\omega_j\big]dm\,\mathbf e_i
=I_{ij}\omega_j\,\mathbf e_i=\mathbf I\boldsymbol\omega.
\end{aligned}
$$
动能：
$$T=\frac12\int\mathbf v\cdot\mathbf v\,dm=\frac12\int(\boldsymbol\omega\times\mathbf r)\cdot(\boldsymbol\omega\times\mathbf r)dm=\frac12\boldsymbol\omega\cdot\int\mathbf r\times(\boldsymbol\omega\times\mathbf r)dm=\frac12\boldsymbol\omega\cdot\mathbf L.$$

物理/几何意义：惯性张量对称正定，可在主轴系对角化 $\mathbf I=\operatorname{diag}(I_1,I_2,I_3)$，$\mathbf L$ 一般不与 $\boldsymbol\omega$ 平行。

极限与特例：薄板绕法向 $I_z=I_x+I_y$（垂直轴定理）；长方体
$$I_x=\tfrac1{12}m(b^2+c^2),\ I_y=\tfrac1{12}m(a^2+c^2),\ I_z=\tfrac1{12}m(a^2+b^2).$$

典型应用：航天器姿态、机器人、陀螺。

**知识点3　欧拉方程**

定义：在刚体主轴系中，角动量定理给出
$$I_1\dot\omega_1-(I_2-I_3)\omega_2\omega_3=M_1,$$
$$I_2\dot\omega_2-(I_3-I_1)\omega_3\omega_1=M_2,$$
$$I_3\dot\omega_3-(I_1-I_2)\omega_1\omega_2=M_3.$$

推导：在固定系 $\dot{\mathbf L}=\mathbf M$。设主轴系以 $\boldsymbol\omega$ 转动，转动系中导数满足
$$\frac{d\mathbf L}{dt}\Big|_{\text{inertial}}=\frac{d\mathbf L}{dt}\Big|_{\text{body}}+\boldsymbol\omega\times\mathbf L=\mathbf M.$$
在主轴系 $\mathbf L=(I_1\omega_1,I_2\omega_2,I_3\omega_3)$，计算叉积分量：
$$(\boldsymbol\omega\times\mathbf L)_1=\omega_2L_3-\omega_3L_2=(I_3-I_2)\omega_2\omega_3,$$
故
$$I_1\dot\omega_1+(I_3-I_2)\omega_2\omega_3=M_1,$$
整理即欧拉方程。

物理/几何意义：陀螺项源于坐标系转动，导致耦合非线性。

极限与特例：无力矩 $\mathbf M=0$，自由刚体可绕最大/最小惯量主轴稳定旋转，中间轴不稳定（网球拍定理）。

典型应用：卫星自旋稳定、陀螺仪、网球拍翻转。

**知识点4　欧拉角与姿态运动学**

定义：欧拉角 $(\varphi,\theta,\psi)$（进动、章动、自转）将体轴系绕 $z,x',z''$ 依次旋转到目标姿态。

推导角速度：三次旋转角速度分别为 $\dot\varphi\mathbf e_z$、$\dot\theta\mathbf e_{x'}$、$\dot\psi\mathbf e_{z''}$，投影到体轴系得
$$
\begin{aligned}
\omega_1&=\dot\varphi\sin\theta\sin\psi+\dot\theta\cos\psi,\\
\omega_2&=\dot\varphi\sin\theta\cos\psi-\dot\theta\sin\psi,\\
\omega_3&=\dot\varphi\cos\theta+\dot\psi.
\end{aligned}
$$
证明思路：单位矢量变换
$$\mathbf e_z=\sin\theta\sin\psi\,\mathbf e_1+\sin\theta\cos\psi\,\mathbf e_2+\cos\theta\,\mathbf e_3,$$
$$\mathbf e_{x'}=\cos\psi\,\mathbf e_1-\sin\psi\,\mathbf e_2,$$
代入 $\boldsymbol\omega=\dot\varphi\mathbf e_z+\dot\theta\mathbf e_{x'}+\dot\psi\mathbf e_3$ 即得上式。

几何意义：欧拉角给出刚体姿态的最小三维参数化（存在 $\theta=0$ 的万向节锁奇点）。

极限与特例：$\theta=0$ 时 $\varphi,\psi$ 不可分（奇点）；对称陀螺常用 $\omega_3=$ 常数。

典型应用：卫星姿态解算、飞行器飞控、机器人学。

**知识点5　对称陀螺的进动**

定义：轴对称刚体 $I_1=I_2=I_\perp$，$I_3$，重力矩作用下绕定点转动。

推导（自由对称陀螺）：欧拉方程
$$I_\perp\dot\omega_1-(I_\perp-I_3)\omega_2\omega_3=0,$$
$$I_\perp\dot\omega_2-(I_3-I_\perp)\omega_3\omega_1=0,$$
$$I_3\dot\omega_3=0\ \Rightarrow\ \omega_3=\text{const}.$$
令 $\Omega=\dfrac{I_3-I_\perp}{I_\perp}\omega_3$，则
$$\dot\omega_1=-\Omega\omega_2,\qquad \dot\omega_2=\Omega\omega_1.$$
解得 $\omega_1=A\cos(\Omega t+\phi)$，$\omega_2=A\sin(\Omega t+\phi)$，$\therefore$ 横向角速度在体轴系中匀速旋转，即在空间系中产生进动。

物理意义：陀螺的定向性与进动性；进动角速度由 $I_3,I_\perp,\omega_3$ 决定。

极限与特例：$I_3=I_\perp$ 球陀螺 $\Omega=0$，$\boldsymbol\omega$ 恒定。

典型应用：陀螺罗盘、惯性导航、卫星自旋稳定。

**知识点6　刚体定点转动稳定性**

定义：无力矩刚体绕主轴旋转，研究小扰动的稳定性。

推导：设 $\boldsymbol\omega=(\epsilon_1,\Omega+\epsilon_2,\epsilon_3)$ 绕主轴 2。线性化欧拉方程
$$I_1\dot\epsilon_1-(I_2-I_3)\Omega\epsilon_3=0,$$
$$I_3\dot\epsilon_3-(I_1-I_2)\Omega\epsilon_1=0.$$
求导消元：
$$\dot\epsilon_1=\frac{I_2-I_3}{I_1}\Omega\epsilon_3\ \Rightarrow\ \ddot\epsilon_1=\frac{I_2-I_3}{I_1}\Omega\dot\epsilon_3
=\frac{(I_2-I_3)(I_1-I_2)}{I_1I_3}\Omega^2\epsilon_1.$$
$\therefore$ 稳定条件为系数为负：
$$(I_2-I_3)(I_1-I_2)<0.$$
即 $I_2$ 不是中间值：绕最大惯量轴（$I_2>I_1,I_3$）或最小惯量轴（$I_2<I_1,I_3$）稳定，绕中间轴不稳定。

物理意义：转动动能在主轴方向分布决定线性稳定性，即中间轴定理（网球拍定理）。

极限与特例：$I_1=I_3$ 时绕 2 轴中性稳定。

典型应用：卫星姿态稳定设计、空间系绳、刚体翻滚。

**知识点7　动静法与刚体平衡**

定义：刚体平衡需合力与合力矩为零：
$$\sum\mathbf F_i=0,\qquad \sum\mathbf r_i\times\mathbf F_i=0.$$

推导：对刚体，由质点系一般定理 $\mathbf F^{\text{ext}}=M\mathbf a_C$，平衡要求 $\mathbf a_C=0$ 即 $\mathbf F^{\text{ext}}=0$；角动量定理要求 $\mathbf M_O^{\text{ext}}=d\mathbf L_O/dt=0$，即 $\mathbf M^{\text{ext}}=0$。对任意参考点，若 $\mathbf F^{\text{ext}}=0$，则 $\mathbf M_O$ 与参考点无关：
$$\mathbf M_{O'}=\sum(\mathbf r_i-\mathbf{OO'})\times\mathbf F_i=\mathbf M_O-\mathbf{OO'}\times\mathbf F^{\text{ext}}=-\mathbf{OO'}\times\mathbf F^{\text{ext}}=\mathbf M_O.$$

物理/几何意义：平衡的力多边形闭合与力矩闭合。

极限与特例：静定与静不定问题；摩擦平衡的静摩擦锥。

典型应用：桁架、支架、起重与着陆器支腿静力分析。

### 第7章 应力与应变分析

**知识点1　应力张量与柯西公式**

定义：物体内一点处法向为 $\mathbf n$ 的面上的面力密度
$$\mathbf t^{(\mathbf n)}=\lim_{\Delta S\to0}\frac{\Delta\mathbf F}{\Delta S}.$$ 存在二阶张量 $\boldsymbol\sigma$ 使 $\mathbf t^{(\mathbf n)}=\boldsymbol\sigma\cdot\mathbf n$，即 $t_i^{(\mathbf n)}=\sigma_{ij}n_j$。

推导（柯西四面体）：取一四面体，三个面平行坐标面，斜面法向 $\mathbf n=(n_x,n_y,n_z)$。设斜面面积 $\Delta S$，其余三面面积 $n_x\Delta S,n_y\Delta S,n_z\Delta S$。由动量平衡（体积力与惯性力为高阶小量）：
$$\mathbf t^{(\mathbf n)}\Delta S-\mathbf t^{(x)}n_x\Delta S-\mathbf t^{(y)}n_y\Delta S-\mathbf t^{(z)}n_z\Delta S=0.$$
$\therefore \mathbf t^{(\mathbf n)}=\mathbf t^{(x)}n_x+\mathbf t^{(y)}n_y+\mathbf t^{(z)}n_z$。记 $(\mathbf t^{(j)})_i=\sigma_{ij}$，即得 $t_i=\sigma_{ij}n_j$。

物理/几何意义：应力是面上单位面积的作用力，完整描述需 9 个分量；由力矩平衡 $\sigma_{ij}=\sigma_{ji}$，独立分量 6 个。

极限与特例：主应力面法向为 $\mathbf n$ 且 $\mathbf t\parallel\mathbf n$，即 $(\boldsymbol\sigma-\lambda\mathbf I)\mathbf n=0$。

典型应用：强度校核、有限元、断裂力学。

**知识点2　应力平衡方程**

定义：静止连续体内部每点的应力满足
$$\sigma_{ij,j}+f_i=0,$$
即 $\nabla\cdot\boldsymbol\sigma+\mathbf f=0$。

推导：取任意体积 $V$，总面力 $\oint_{\partial V}\boldsymbol\sigma\cdot\mathbf n\,dS$ 与体积力 $\int_V\mathbf f\,dV$ 平衡：
$$\oint_{\partial V}\sigma_{ij}n_j\,dS+\int_Vf_i\,dV=0.$$
由高斯散度定理
$$\oint_{\partial V}\sigma_{ij}n_j\,dS=\int_V\frac{\partial\sigma_{ij}}{\partial x_j}dV.$$
$\therefore \int_V(\sigma_{ij,j}+f_i)dV=0$ 对任意 $V$ 成立，$\therefore$ 被积函数为零。

物理意义：应力散度平衡体积力；是静力学基本方程。

极限与特例：无体积力时 $\sigma_{ij,j}=0$；平面问题 $\partial_x\sigma_{xx}+\partial_y\sigma_{xy}=0$ 等。

典型应用：结构应力分析、边界层与薄壁结构。

**知识点3　应变张量与几何方程**

定义：小变形位移场 $u_i$，应变张量
$$\varepsilon_{ij}=\frac12(u_{i,j}+u_{j,i}).$$

推导：取邻近两点 $d\mathbf r$，变形后线元 $d\mathbf r'=d\mathbf r+d\mathbf u$，其中 $du_i=u_{i,j}dx_j$。长度平方变化
$$
\begin{aligned}
|d\mathbf r'|^2-|d\mathbf r|^2&=2\varepsilon_{ij}dx_idx_j+O(|\nabla\mathbf u|^2),\\
2\varepsilon_{ij}dx_idx_j&=(u_{i,j}+u_{j,i})dx_idx_j.
\end{aligned}
$$
\therefore 仅对称部分贡献长度变化，$\varepsilon_{ij}=\tfrac12(u_{i,j}+u_{j,i})$。旋转部分
$$\omega_{ij}=\tfrac12(u_{i,j}-u_{j,i})$$
不改变长度，对应刚体转动。

物理/几何意义：$\varepsilon_{ii}$ 为轴向伸长率；$\varepsilon_{ij}\ (i\neq j)$ 为工程切应变的一半。体积应变 $\varepsilon_{kk}=u_{k,k}=\nabla\cdot\mathbf u$。

极限与特例：协调方程保证应变可由单值位移导出：
$$\varepsilon_{ij,kl}+\varepsilon_{kl,ij}-\varepsilon_{ik,jl}-\varepsilon_{jl,ik}=0.$$

典型应用：应变片测量、有限元位移法、热变形。

**知识点4　主应力与应力不变量**

定义：主应力为特征方程
$$\det(\sigma_{ij}-\sigma\delta_{ij})=0$$
的根，即
$$\sigma^3-I_1\sigma^2+I_2\sigma-I_3=0,$$
其中
$$I_1=\sigma_{kk},\quad I_2=\frac12(\sigma_{ii}\sigma_{jj}-\sigma_{ij}\sigma_{ji}),\quad I_3=\det(\sigma_{ij}).$$

推导：主轴方向要求 $\sigma_{ij}n_j=\sigma n_i$，即 $(\sigma_{ij}-\sigma\delta_{ij})n_j=0$。非零解要求系数矩阵奇异。$\because \boldsymbol\sigma$ 对称，$\therefore$ 特征值为实、特征向量正交。$I_k$ 为相似不变量，因为在正交变换下特征多项式不变。

物理/几何意义：主应力是法向应力的极值，主方向相互正交；应力状态可由 $I_1,I_2,I_3$ 完全确定（等效应力基础）。

极限与特例：三向等应力 $\sigma_{ij}=p\delta_{ij}$，所有方向皆主方向；von Mises 应力
$$\sigma_v=\sqrt{\tfrac12[(\sigma_1-\sigma_2)^2+(\sigma_2-\sigma_3)^2+(\sigma_3-\sigma_1)^2]}.$$

典型应用：强度理论、屈服判据、实验应力分析。

**知识点5　应力偏量与静水分量**

定义：
$$\sigma_{ij}=-p\delta_{ij}+s_{ij},\qquad p=-\frac13\sigma_{kk},\qquad s_{ij}=\sigma_{ij}-\frac13\sigma_{kk}\delta_{ij}.$$

推导：将张量分解为各向同性（球形）部分与偏斜部分。验证偏量迹为零：
$$s_{kk}=\sigma_{kk}-\frac13\sigma_{kk}\cdot3=0.$$
且 $s_{ij}$ 与 $\delta_{ij}$ 正交：$s_{ij}\delta_{ij}=s_{kk}=0$（指标缩并），$\therefore$ 分解唯一。

物理/几何意义：$p$ 描述体积压缩（静水压），$s_{ij}$ 描述形状畸变，是金属塑性与蠕变的主要驱动力。

极限与特例：纯剪切 $\sigma_{kk}=0$，$p=0$，$s=\sigma$。

典型应用：J2 塑性理论、偏应力张量不变量
$$J_2=\tfrac12s_{ij}s_{ij},\qquad \sigma_v=\sqrt{3J_2}.$$

**知识点6　莫尔圆**

定义：对给定法向 $\mathbf n$，正应力 $\sigma_n=\sigma_{ij}n_in_j$，切应力 $\tau=\sqrt{|\mathbf t|^2-\sigma_n^2}$。莫尔圆给出 $\sigma_n$-$\tau$ 平面临界区域。

推导（平面应力简化）：设主应力 $\sigma_1\ge\sigma_2\ge\sigma_3$，斜面上
$$\sigma_{ij}=\operatorname{diag}(\sigma_1,\sigma_2,\sigma_3),\quad n=(\cos\alpha,\cos\beta,\cos\gamma).$$
则
$$\sigma_n=\sigma_1n_1^2+\sigma_2n_2^2+\sigma_3n_3^2,\qquad \tau^2+|\sigma|^2=\sigma_1^2n_1^2+\sigma_2^2n_2^2+\sigma_3^2n_3^2.$$
利用 $n_1^2+n_2^2+n_3^2=1$ 消元可得三个圆
$$(\sigma_n-C_i)^2+\tau^2=R_i^2,\quad C_i=\frac{\sigma_j+\sigma_k}{2},\ n_i=0.$$
对二维平面应力，取与主向夹角 $\theta$：
$$\sigma_{\theta}=\frac{\sigma_x+\sigma_y}{2}+\frac{\sigma_x-\sigma_y}{2}\cos2\theta+\tau_{xy}\sin2\theta,$$
$$\tau_{\theta}=-\frac{\sigma_x-\sigma_y}{2}\sin2\theta+\tau_{xy}\cos2\theta.$$

物理/几何意义：可能应力状态落在三圆之间的区域内；最大切应力为 $(\sigma_1-\sigma_3)/2$。

极限与特例：单向拉伸 $\sigma_1=\sigma,\ \sigma_2=\sigma_3=0$，最大切应力 $\sigma/2$ 在 $45^\circ$ 面。

典型应用：实验数据处理、脆性材料破坏面预估。

**知识点7　体积应变与应力功率**

定义：线弹性时体积应变 $\varepsilon_{kk}$ 与静水压关系，及单位体积应变能（应变能密度）。

推导：小变形下体积相对变化率为 $\dot\varepsilon_{kk}=v_{i,i}=\nabla\cdot\mathbf v$。应力做功的功率密度
$$\dot W=\sigma_{ij}\dot\varepsilon_{ij}=\sigma_{ij}\dot u_{i,j}=\frac{\partial}{\partial x_j}(\sigma_{ij}\dot u_i)-\sigma_{ij,j}\dot u_i.$$
在平衡方程 $\sigma_{ij,j}+f_i=0$ 下，第二项为 $f_i\dot u_i$，配合高斯定理得外力功率等于应力内功率加体力功率。将 $\sigma=\tfrac13\sigma_{kk}\delta+s$，$\varepsilon=\tfrac13\varepsilon_{kk}\delta+e$（$e$ 偏应变），则
$$\dot W=\frac13\sigma_{kk}\dot\varepsilon_{kk}+s_{ij}\dot e_{ij}.$$
$\therefore$ 功率分解为体积变形与形状畸变两部分。

物理意义：体积改变与畸变在能量上解耦（各向同性弹性）。

极限与特例：不可压缩 $\varepsilon_{kk}=0$ 时仅畸变能。

典型应用：应力功共轭对、有限元虚功方程、塑性耗散。

### 第8章 弹性本构关系与边值问题

**知识点1　广义胡克定律**

定义：线弹性各向异性材料
$$\sigma_{ij}=C_{ijkl}\varepsilon_{kl},\qquad \varepsilon_{ij}=S_{ijkl}\sigma_{kl},$$
其中 $C$ 为四阶弹性张量、$S$ 为柔度张量。

推导对称性：由 $\sigma_{ij}=\sigma_{ji}$ 得 $C_{ijkl}=C_{jikl}$；由 $\varepsilon_{kl}=\varepsilon_{lk}$ 得 $C_{ijkl}=C_{ijlk}$。由应变能密度 $w=\tfrac12C_{ijkl}\varepsilon_{ij}\varepsilon_{kl}$，
$$\sigma_{ij}=\frac{\partial w}{\partial\varepsilon_{ij}}=C_{ijkl}\varepsilon_{kl}.$$
对 $i\leftrightarrow k,\ j\leftrightarrow l$ 交换二阶偏导，得 $C_{ijkl}=C_{klij}$。$\therefore$ 一般各向异性独立常数 $21$ 个（三斜），单斜 $13$，正交 $9$，横观各向同性 $5$，各向同性 $2$。

物理意义：本构关系反映材料内部微观结构的对称性。

极限与特例：各向同性只有两个独立常数 $\lambda,\mu$。

典型应用：复合材料层合板、晶体弹性。

**知识点2　各向同性广义胡克定律**

定义：
$$\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij},$$
或
$$\varepsilon_{ij}=\frac{1+\nu}{E}\sigma_{ij}-\frac{\nu}{E}\sigma_{kk}\delta_{ij}.$$

推导：各向同性四阶张量只能是 $\delta_{ij}\delta_{kl}$ 与 $\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}$ 的线性组合：
$$C_{ijkl}=\lambda\delta_{ij}\delta_{kl}+\mu(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}).$$
代入 $\sigma_{ij}=C_{ijkl}\varepsilon_{kl}$：
$$\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+\mu(\varepsilon_{ij}+\varepsilon_{ji})=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij}.$$
对之求迹：
$$\sigma_{kk}=3\lambda\varepsilon_{kk}+2\mu\varepsilon_{kk}=(3\lambda+2\mu)\varepsilon_{kk}.$$
再反解 $\varepsilon_{ij}$ 可得柔度形式，其中
$$\lambda=\frac{E\nu}{(1+\nu)(1-2\nu)},\qquad \mu=G=\frac{E}{2(1+\nu)}.$$

物理意义：$\lambda,\mu$ 为拉梅常数；$E$ 杨氏模量、$\nu$ 泊松比、$G$ 剪切模量。

极限与特例：$\nu\to1/2$ 时 $\lambda\to\infty$ 但 $\mu$ 有限，接近不可压缩。

典型应用：金属弹性、橡胶近似、有限元材料输入。

**知识点3　弹性力学基本方程组**

定义：三维线弹性静力学由平衡、几何、本构及边界条件组成。

推导汇总：
平衡方程 $\sigma_{ij,j}+f_i=0$（3 个）；
几何方程 $\varepsilon_{ij}=\tfrac12(u_{i,j}+u_{j,i})$（6 个）；
本构方程 $\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij}$（6 个）；
边界 $u_i=\bar u_i$ 在 $S_u$，$\sigma_{ij}n_j=\bar t_i$ 在 $S_t$。
未知量 $u_i$（3）、$\varepsilon_{ij}$（6）、$\sigma_{ij}$（6）共 15 个，方程数 15，适定。

物理意义：三类方程分别代表动量守恒、几何相容、材料响应。

极限与特例：动力学问题加惯性项 $\rho\ddot u_i$，成为弹性动力学。

典型应用：结构仿真、接触与热应力。

**知识点4　位移法：拉梅—纳维方程**

定义：以位移为基本未知量，
$$(\lambda+\mu)u_{j,ji}+\mu u_{i,jj}+f_i=\rho\ddot u_i.$$

推导：将本构代入平衡方程（动力学形式 $\sigma_{ij,j}+f_i=\rho\ddot u_i$）：
$$\sigma_{ij,j}=\lambda\varepsilon_{kk,j}\delta_{ij}+2\mu\varepsilon_{ij,j}=\lambda\varepsilon_{kk,i}+\mu(u_{i,jj}+u_{j,ij}).$$
$\therefore$
$$\lambda u_{j,ji}+\mu u_{i,jj}+\mu u_{j,ij}+f_i=\rho\ddot u_i.$$
合并得 $(\lambda+\mu)u_{j,ji}+\mu u_{i,jj}+f_i=\rho\ddot u_i$。矢量形式
$$(\lambda+\mu)\nabla(\nabla\cdot\mathbf u)+\mu\nabla^2\mathbf u+\mathbf f=\rho\ddot{\mathbf u}.$$

物理意义：用纳维算子表示弹性波传播。

极限与特例：静力问题删除惯性项；对无旋部分得纵波速度
$$c_p=\sqrt{(\lambda+2\mu)/\rho},\qquad c_s=\sqrt{\mu/\rho}.$$

典型应用：弹性波、无损检测、地震波。

**知识点5　应力法：贝蒂—米歇尔斯方程**

定义：以应力为基本未知量，需满足平衡与协调：
$$\nabla^2\sigma_{ij}+\frac{1}{1+\nu}\sigma_{kk,ij}=-\Big(\frac{\nu}{1-\nu}\delta_{ij}f_{k,k}+f_{i,j}+f_{j,i}\Big).$$

推导：将本构反解 $\varepsilon_{ij}=\frac{1}{E}[(1+\nu)\sigma_{ij}-\nu\sigma_{kk}\delta_{ij}]$ 代入协调方程
$$\varepsilon_{ij,kl}+\varepsilon_{kl,ij}-\varepsilon_{ik,jl}-\varepsilon_{jl,ik}=0.$$
缩并指标并利用平衡方程消去应力偏导，即得 Beltrami–Michell 方程。无体力时为
$$\nabla^2\sigma_{ij}+\frac{1}{1+\nu}\sigma_{kk,ij}=0.$$
对之求迹得 $\nabla^2\sigma_{kk}=0$，$\therefore$ 静水应力调和。

物理意义：应力法在边界以面力给定时自然适用。

极限与特例：二维问题化为艾里应力函数的双调和方程。

典型应用：接触、孔洞应力集中。

**知识点6　平面问题与艾里应力函数**

定义：平面应力/平面应变下，引入艾里函数 $\phi(x,y)$：
$$\sigma_{xx}=\frac{\partial^2\phi}{\partial y^2},\quad \sigma_{yy}=\frac{\partial^2\phi}{\partial x^2},\quad \sigma_{xy}=-\frac{\partial^2\phi}{\partial x\partial y}.$$

推导：无体力平衡方程
$$\frac{\partial\sigma_{xx}}{\partial x}+\frac{\partial\sigma_{xy}}{\partial y}=0,\qquad \frac{\partial\sigma_{xy}}{\partial x}+\frac{\partial\sigma_{yy}}{\partial y}=0.$$
第一式可写 $\partial_x(\sigma_{xx})+\partial_y(\sigma_{xy})=0$，令 $\sigma_{xx}=\partial_y\psi$，则 $\sigma_{xy}=-\partial_x\psi$。第二式 $\partial_x(\sigma_{xy})+\partial_y(\sigma_{yy})=0$，令 $\sigma_{yy}=\partial_x\chi$，比较得 $\chi=\partial_y\phi$ 与 $\psi=\partial_y\phi$，即得上面三式。将其代入协调方程得
$$\nabla^4\phi=\Big(\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}\Big)^2\phi=0.$$

物理意义：自动满足平衡，问题归为求双调和函数。

极限与特例：多项式解可表示均匀应力；Williams 解给出裂纹尖端 $r^{-1/2}$ 奇性。

典型应用：孔边应力集中、断裂力学、接触问题。

**知识点7　圣维南原理**

定义：作用于弹性体小区域内的自平衡力系，对远处应力的影响可忽略。

推导（量级估计）：在边界局部 $d\times d$ 上作用自平衡面力，合力与合力矩为零。对距离 $r\gg d$ 处，由量纲分析应力衰减
$$\sigma\sim\frac{Fd^2}{r^4}\cdot\text{（偶极子/四极子）},$$
而合力的应力衰减为 $F/r^2$。自平衡使低阶多极矩消失，$\therefore$ 实际衰减更快。更严格地，可由圣维南扭率与衰减估计给出指数衰减。

物理意义：边界条件的等效替换对远离边界区域影响小。

极限与特例：靠近加载区的应力分布不同，须使用精确解。

典型应用：杆件端部加载等效、结构简化建模。

### 第9章 能量原理与变分方法

**知识点1　应变能与余应变能**

定义：线弹性体应变能
$$U=\frac12\int_V\sigma_{ij}\varepsilon_{ij}\,dV=\frac12\int_V C_{ijkl}\varepsilon_{ij}\varepsilon_{kl}\,dV.$$
余应变能
$$U^*=\int_V\sigma_{ij}\varepsilon_{ij}\,dV-U.$$

推导：外力在弹性变形中做功转换为应变能。对单位体积，应力从 0 增加到 $\sigma_{ij}$，做功为 $\int_0^{\varepsilon_{ij}}\sigma_{ij}\,d\varepsilon_{ij}$；线性时等于 $\frac12\sigma_{ij}\varepsilon_{ij}$。由本构 $\sigma_{ij}=C_{ijkl}\varepsilon_{kl}$，所以
$$w=\frac12C_{ijkl}\varepsilon_{ij}\varepsilon_{kl}=\frac12\sigma_{ij}\varepsilon_{ij}.$$
对之求导验证共轭 $\sigma_{ij}=\partial w/\partial\varepsilon_{ij}$，$\varepsilon_{ij}=\partial w/\partial\sigma_{ij}$。对非线性弹性，$U+U^*=\int\sigma:d\varepsilon$。

物理/几何意义：在应力—应变图中，应变能是曲线下方面积，余能是右侧/上方面积。

极限与特例：线弹性两者相等 $U=U^*$；非线性时不等。

典型应用：能量法求位移、卡氏定理、断裂能量释放率。

**知识点2　虚功原理与虚功方程**

定义：对任意运动允许（满足位移边界与协调）的虚位移 $\delta u_i$，平衡等价于
$$\int_V\sigma_{ij}\delta\varepsilon_{ij}\,dV=\int_Vf_i\delta u_i\,dV+\int_{S_t}\bar t_i\delta u_i\,dS.$$

推导：由平衡方程 $\sigma_{ij,j}+f_i=0$，两边乘以 $\delta u_i$ 并在 $V$ 上积分：
$$\int_V(\sigma_{ij,j}+f_i)\delta u_i\,dV=0.$$
对第一项分部积分
$$\int_V\sigma_{ij,j}\delta u_i\,dV=\int_{S_t+\Sigma}\sigma_{ij}\delta u_i n_j\,dS-\int_V\sigma_{ij}\delta u_{i,j}\,dV.$$
在位移已知边界 $S_u$ 上 $\delta u_i=0$，$\therefore$ 面积分只在 $S_t$。又
$$\sigma_{ij}\delta u_{i,j}=\sigma_{ij}\delta\varepsilon_{ij},$$
因为 $\sigma_{ij}$ 对称且 $\delta\varepsilon_{ij}$ 为对称部分。$\therefore$
$$\int_V\sigma_{ij}\delta\varepsilon_{ij}\,dV-\int_{S_t}\bar t_i\delta u_i\,dS-\int_Vf_i\delta u_i\,dV=0.$$

物理/几何意义：外力虚功等于内应力虚功。

极限与特例：刚体虚位移使两边为零，故不提供信息。

典型应用：有限元弱形式、夹层结构分析。

**知识点3　最小势能原理**

定义：总势能
$$\Pi=U-W_{\text{ext}}=\frac12\int_VC_{ijkl}\varepsilon_{ij}\varepsilon_{kl}\,dV-\Big(\int_Vf_iu_i\,dV+\int_{S_t}\bar t_iu_i\,dS\Big),$$
线弹性体的真实位移使 $\Pi$ 取极小。

推导：设 $u$ 为真实解，$u+\delta u$ 为运动允许，$\delta u|_{S_u}=0$。展开
$$\Pi(u+\delta u)-\Pi(u)=\delta\Pi+\frac12\delta^2\Pi.$$
一次变分
$$\delta\Pi=\int_V\sigma_{ij}\delta\varepsilon_{ij}\,dV-\int_Vf_i\delta u_i\,dV-\int_{S_t}\bar t_i\delta u_i\,dS=0,$$
正是虚功方程。二次变分
$$\delta^2\Pi=\int_V\delta\varepsilon_{ij}C_{ijkl}\delta\varepsilon_{kl}\,dV\ge0,$$
$\because C$ 正定。$\therefore$ 驻点即极小点。

物理/几何意义：弹性平衡是势能泛函的最小点，可用于近似求解。

极限与特例：不稳定（几何非线性软化）时可能为极大/鞍点。

典型应用：瑞利—里兹法、有限元。

**知识点4　最小余能原理**

定义：总余能
$$\Pi^*=U^*-W^*=\frac12\int_V\sigma_{ij}S_{ijkl}\sigma_{kl}\,dV-\int_{S_u}\sigma_{ij}n_j\bar u_i\,dS,$$
真实应力场使 $\Pi^*$ 取极小。

推导：设静力允许应力场 $\sigma$，在 $V$ 内满足 $\sigma_{ij,j}+f_i=0$，在 $S_t$ 上满足 $\sigma_{ij}n_j=\bar t_i$。由余虚功原理
$$\delta\Pi^*=\int_V\varepsilon_{ij}\delta\sigma_{ij}\,dV-\int_{S_u}\bar u_i\delta(\sigma_{ij}n_j)\,dS=0,$$
这正是位移边界条件（通过协调方程）的弱形式。二次变分非负。

物理/几何意义：应力场近似（如力法、平衡元）的变分基础。

极限与特例：线弹性时与最小势能互为对偶。

典型应用：杆系内力法、下限塑性分析。

**知识点5　瑞利—里兹法**

定义：取试探解
$$u_i\approx\sum_{m=1}^{N}a_m\phi_i^{(m)}(\mathbf x),$$
其中 $\phi^{(m)}$ 满足几何边界条件，代入 $\Pi$ 得关于 $a_m$ 的二次函数，由 $\partial\Pi/\partial a_m=0$ 求解。

推导：代入应变能
$$U=\frac12\sum_{m,n}a_ma_nK_{mn},\qquad K_{mn}=\int_VC_{ijkl}\phi_{i,j}^{(m)}\phi_{k,l}^{(n)}\,dV,$$
外力功 $W=\sum_m a_mF_m$。$\therefore$
$$\Pi=\frac12\mathbf a^T K\mathbf a-\mathbf F^T\mathbf a.$$
由 $\partial\Pi/\partial\mathbf a=0$：
$$K\mathbf a=\mathbf F,\qquad a_m=\ldots$$

物理/几何意义：将无穷维变分问题投影到有限维试函数空间，是有限元的前身。

极限与特例：试探解越接近真实解，能量越小（上界性质）；$N$ 增大收敛。

典型应用：梁板振动、屈曲、桁架。

**知识点6　伽辽金法**

定义：若微分方程 $\mathcal L(u)=0$，取 $u\approx\sum a_m\phi_m$，要求残差与试函数正交：
$$\int_V\mathcal L(u_n)\,\phi_m\,dV=0,\qquad m=1,\dots,N.$$

推导：从加权残差角度，令残差 $R=\mathcal L(u_n)$。在最小二乘意义下选择权函数 $w_m=\partial R/\partial a_m$ 得最小二乘；若取 $w_m=\phi_m$ 即伽辽金。对线性弹性平衡 $\mathcal L(u)=\sigma_{ij,j}+f_i=0$：
$$\int_V(\sigma_{ij,j}+f_i)\phi_m\,dV=0.$$
分部积分并与虚功方程比较，$\therefore$ 伽辽金法与虚功/最小势能等价（对于自伴算子）。

物理/几何意义：试函数无需满足自然边界，但必须满足本质边界；自动包含自然边界条件。

极限与特例：非自伴（对流扩散）时伽辽金可能不稳定，需用 Petrov–Galerkin/SUPG。

典型应用：有限元弱形式、流体/传热。

**知识点7　卡氏定理与单位载荷法**

定义：线弹性体应变能对某广义力的偏导数等于对应广义位移：
$$\Delta_j=\frac{\partial U}{\partial F_j}.$$
更一般地，对广义位移用单位载荷法：
$$\Delta=\sum\int\frac{M\bar M}{EI}dx+\sum\int\frac{N\bar N}{EA}dx+\sum\int\frac{T\bar T}{GJ}dx.$$

推导：由 $\Pi=U-W$，对广义力 $F_j$ 求导并用 $\partial U/\partial\varepsilon_{ij}=\sigma_{ij}$ 与链式法则：
$$\frac{\partial U}{\partial F_j}=\int_V\sigma_{ij}\frac{\partial\varepsilon_{ij}}{\partial F_j}\,dV.$$
对线性弹性 $\varepsilon=C^{-1}\sigma$，由 $\partial W/\partial u_j=F_j$ 及 $\Pi$ 的极值条件得
$$\frac{\partial\Pi}{\partial F_j}=\frac{\partial U}{\partial F_j}-\Delta_j=0\ \Rightarrow\ \Delta_j=\frac{\partial U}{\partial F_j}.$$
对单位载荷法：将 $F=F_0+\bar F$ 代入 $\Delta=\sum\partial U/\partial\bar F_j|_{\bar F=0}$，交叉项给出 $M\bar M/(EI)$ 等。

物理/几何意义：能量对力的梯度给出共轭位移。

极限与特例：非线性弹性须用 $\Delta_j=\partial U^*/\partial F_j$ 或补余能量。

典型应用：桁架/梁挠度快速求解、结构试验。

### 第10章 流体运动学与守恒方程

**知识点1　拉格朗日描述与欧拉描述**

定义：拉格朗日追踪流体质点 $\mathbf x=\mathbf x(\mathbf X,t)$（$\mathbf X$ 为初始位置）；欧拉观察空间场 $\mathbf v=\mathbf v(\mathbf x,t)$。

推导物质导数：拉格朗日物理量 $f=f(\mathbf x(\mathbf X,t),t)$，链式法则
$$\frac{df}{dt}=\frac{\partial f}{\partial t}+\frac{\partial f}{\partial x_i}\frac{dx_i}{dt}=\frac{\partial f}{\partial t}+v_i\frac{\partial f}{\partial x_i}=\frac{\partial f}{\partial t}+(\mathbf v\cdot\nabla)f.$$
加速度
$$\mathbf a=\frac{D\mathbf v}{Dt}=\frac{\partial\mathbf v}{\partial t}+(\mathbf v\cdot\nabla)\mathbf v.$$

物理/几何意义：物质导数 = 当地变化 + 对流变化；流场线由 $\mathbf v$ 积分得流线、迹线。

极限与特例：定常 $\partial_t=0$ 时物质导数=对流项；流动快时对流主导。

典型应用：CFD、粒子追踪、流场测量。

**知识点2　雷诺输运定理**

定义：对流动中物质体 $V_m(t)$ 内的物理量 $\int_{V_m}\phi\,dV$，
$$\frac{d}{dt}\int_{V_m}\phi\,dV=\int_{V_m}\frac{\partial\phi}{\partial t}dV+\oint_{\partial V_m}\phi\mathbf v\cdot\mathbf n\,dS.$$

推导：先写拉格朗日体积积分。用雅可比 $J=\det(\partial x_i/\partial X_j)$，$dV=J\,dV_0$，且
$$\frac{dJ}{dt}=J\,\nabla\cdot\mathbf v.$$
于是
$$\frac{d}{dt}\int_{V_m}\phi\,dV=\int_{V_0}\Big(\frac{d\phi}{dt}J+\phi\frac{dJ}{dt}\Big)dV_0
=\int_{V_m}\Big(\frac{d\phi}{dt}+\phi\nabla\cdot\mathbf v\Big)dV.$$
又 $\dfrac{d\phi}{dt}=\partial_t\phi+(\mathbf v\cdot\nabla)\phi$，代入并合并：
$$\frac{d\phi}{dt}+\phi\nabla\cdot\mathbf v=\frac{\partial\phi}{\partial t}+\nabla\cdot(\phi\mathbf v).$$
$\therefore$ 得雷诺输运定理。

物理/几何意义：物质体内总量的变化等于当地变化加边界通量。

极限与特例：$\phi=1$ 给 $dV_m/dt=\int\nabla\cdot\mathbf v\,dV$。

典型应用：推导连续、动量、能量方程。

**知识点3　连续性方程**

定义：质量守恒。

推导：取物质体，质量 $m=\int_{V_m}\rho\,dV$ 不随时间变化，$dm/dt=0$。由雷诺输运定理（$\phi=\rho$）：
$$0=\frac{d}{dt}\int_{V_m}\rho\,dV=\int_{V_m}\Big(\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf v)\Big)dV.$$
对任意物质体成立，$\therefore$
$$\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf v)=0.$$
不可压 $\rho=$ 常数时
$$\nabla\cdot\mathbf v=0.$$
展开形式
$$\frac{\partial\rho}{\partial t}+v_i\frac{\partial\rho}{\partial x_i}+\rho\frac{\partial v_i}{\partial x_i}=0\ \Leftrightarrow\ \frac{D\rho}{Dt}+\rho\nabla\cdot\mathbf v=0.$$

物理意义：密度的物质变化由体积膨胀/压缩引起。

极限与特例：定常可压 $\nabla\cdot(\rho\mathbf v)=0$；一维定常 $\rho vA=$ const。

典型应用：管道、喷管、气象与 CFD。

**知识点4　动量方程**

定义：牛顿第二定律对连续体。

推导：对物质体，合外力（面力+体力）等于动量变化率：
$$\int_{V_m}\mathbf f\,dV+\oint_{\partial V_m}\boldsymbol\sigma\cdot\mathbf n\,dS=\frac{d}{dt}\int_{V_m}\rho\mathbf v\,dV.$$
由雷诺输运定理（$\phi=\rho\mathbf v$）与高斯定理：
$$\text{右边}=\int_{V_m}\Big(\frac{\partial(\rho v_i)}{\partial t}+\frac{\partial(\rho v_iv_j)}{\partial x_j}\Big)dV,$$
面力 $=\int\partial_j\sigma_{ij}\,dV$。合并并对任意体积取零：
$$\frac{\partial(\rho v_i)}{\partial t}+\frac{\partial(\rho v_iv_j)}{\partial x_j}=f_i+\frac{\partial\sigma_{ij}}{\partial x_j}.$$
展开左侧并用连续性方程
$$\rho\frac{Dv_i}{Dt}=\frac{\partial\sigma_{ij}}{\partial x_j}+f_i.$$
矢量形式
$$\rho\frac{D\mathbf v}{Dt}=\nabla\cdot\boldsymbol\sigma+\mathbf f.$$

物理意义：流体微团动量变化由应力散度与体力驱动。

极限与特例：无粘（欧拉）$\sigma=-p\mathbf I$，得 $\rho D\mathbf v/Dt=-\nabla p+\mathbf f$。

典型应用：管道受力、叶轮机械、洪水波。

**知识点5　能量方程**

定义：热力学能与动能守恒。

推导：对物质体，总能量 $E=\int\rho(e+\tfrac12v^2)dV$。外力功率加面力功率加传导热：
$$\frac{dE}{dt}=\int\mathbf f\cdot\mathbf v\,dV+\oint_{\partial V}\mathbf v\cdot(\boldsymbol\sigma\cdot\mathbf n)\,dS-\oint_{\partial V}\mathbf q\cdot\mathbf n\,dS.$$
由雷诺输运与高斯定理化简，并用动量方程消去动能项，可得内能方程
$$\rho\frac{De}{Dt}=\sigma_{ij}\frac{\partial v_i}{\partial x_j}-\nabla\cdot\mathbf q.$$
对牛顿流体 $\sigma_{ij}=-p\delta_{ij}+\tau_{ij}$，$\tau_{ij}\partial_jv_i=\Phi$ 为黏性耗散。温度形式
$$\rho c_v\frac{DT}{Dt}=\Phi-\nabla\cdot\mathbf q.$$

物理意义：机械功与黏性耗散转化为内能，热传导再分配温度。

极限与特例：无粘绝热时 $Ds/Dt=0$（等熵）。

典型应用：高超音速气动加热、燃烧、传热。

**知识点6　伯努利方程**

定义：沿流线/无旋流动的理想不可压定常流
$$\frac{p}{\rho}+\frac12v^2+gz=\text{const}.$$

推导：沿流线的欧拉方程（无粘、定常、不可压）：
$$\rho\frac{D\mathbf v}{Dt}=\rho(\mathbf v\cdot\nabla)\mathbf v=-\nabla p-\rho g\nabla z.$$
利用矢量恒等式
$$(\mathbf v\cdot\nabla)\mathbf v=\nabla\Big(\frac{v^2}{2}\Big)-\mathbf v\times(\nabla\times\mathbf v).$$
沿流线方向 $d\mathbf l=\mathbf v\,dt$，叉积项点乘 $\mathbf v$ 为零，$\therefore$
$$\rho\,d\Big(\frac{v^2}{2}\Big)=-dp-\rho g\,dz.$$
除以 $\rho$：
$$d\Big(\frac{p}{\rho}+\frac{v^2}{2}+gz\Big)=0.$$
$\therefore$ 沿流线常数。无旋时 $\mathbf v\times\boldsymbol\omega=0$，全流场同一常数。

物理意义：压强能、动能、位能之和守恒，是能量守恒的一维形式。

极限与特例：滞止点 $v=0$，$p_0=p+\tfrac12\rho v^2$；可压时用等熵关系。

典型应用：文丘里管、机翼吸力、空速管。

**知识点7　涡量方程与开尔文环量定理**

定义：涡量 $\boldsymbol\omega=\nabla\times\mathbf v$。

推导涡量方程：对动量方程取旋度。无体力、不可压：
$$\frac{\partial\mathbf v}{\partial t}+(\mathbf v\cdot\nabla)\mathbf v=-\frac1\rho\nabla p+\nu\nabla^2\mathbf v.$$
取旋度，利用 $\nabla\times(\mathbf v\cdot\nabla)\mathbf v=(\mathbf v\cdot\nabla)\boldsymbol\omega-(\boldsymbol\omega\cdot\nabla)\mathbf v+\boldsymbol\omega\nabla\cdot\mathbf v$（不可压末项为零），以及 $\nabla\times\nabla p=0$：
$$\frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf v+\nu\nabla^2\boldsymbol\omega.$$

开尔文环量定理：环量 $\Gamma=\oint\mathbf v\cdot d\mathbf l$，由斯托克斯 $\Gamma=\int\boldsymbol\omega\cdot\mathbf n\,dS$。对无粘正压且体力有势的流动，$D\Gamma/Dt=0$。证明：
$$\frac{D\Gamma}{Dt}=\oint\frac{D\mathbf v}{Dt}\cdot d\mathbf l+\underbrace{\oint\mathbf v\cdot d\Big(\frac{D\mathbf l}{Dt}\Big)}_{=0}.$$
第二项因 $D(d\mathbf l)/Dt=d\mathbf v$ 且 $\oint\mathbf v\cdot d\mathbf v=0$。代入欧拉方程
$$\oint\Big(-\frac{\nabla p}{\rho}-\nabla\Psi\Big)\cdot d\mathbf l=0$$
（正压时 $\nabla p/\rho=\nabla\int dp/\rho$ 为全微分）。$\therefore$ 环量守恒。

物理意义：无粘涡旋不生成也不消亡；黏性使涡量扩散耗散。

极限与特例：二维无粘 $D\omega/Dt=0$（涡量冻结）。

典型应用：机翼尾涡、涡动力学、升力产生机理。

### 第11章 纳维—斯托克斯方程与本构

**知识点1　牛顿流体本构关系**

定义：偏应力张量 $\tau_{ij}$ 与变形率张量 $D_{ij}=\frac12(v_{i,j}+v_{j,i})$ 线性且各向同性：
$$\tau_{ij}=2\mu D_{ij}+\lambda(\nabla\cdot\mathbf v)\delta_{ij}.$$

推导：各向同性张量表示定理：$\tau$ 为 $D$ 与 $\delta$ 的各向同性线性张量函数，只能是
$$\tau_{ij}=\alpha D_{ij}+\beta D_{kk}\delta_{ij}.$$
令 $\alpha=2\mu,\ \beta=\lambda$。斯托克斯假设：平均正应力等于热力学压强，即 $\frac13\tau_{ii}=0$（无体积粘性），则
$$\frac13(2\mu+3\lambda)D_{kk}=0\ \Rightarrow\ \lambda=-\frac23\mu.$$
$\therefore$
$$\tau_{ij}=2\mu\Big(D_{ij}-\frac13D_{kk}\delta_{ij}\Big)=2\mu\Big(D_{ij}-\frac13\nabla\cdot\mathbf v\,\delta_{ij}\Big).$$

物理意义：黏性应力抵抗变形速率；不可压时 $\tau=2\mu D$。

极限与特例：$\mu=0$ 为理想流体；非牛顿流体本构为非线性。

典型应用：水、空气、机油流动；湍流模型。

**知识点2　应力分量与应变率分量的对应**

定义：笛卡尔坐标中
$$\tau_{xx}=2\mu\frac{\partial u}{\partial x}-\frac23\mu\nabla\cdot\mathbf v,\quad \tau_{xy}=\mu\Big(\frac{\partial u}{\partial y}+\frac{\partial v}{\partial x}\Big),$$
等。

推导：由 $\tau_{ij}=2\mu(D_{ij}-\frac13D_{kk}\delta_{ij})$ 直接展开。对 $i=j=x$：
$$D_{xx}=\partial u/\partial x,\quad D_{kk}=\partial_ku_k=\nabla\cdot\mathbf v,$$
$\therefore \tau_{xx}=2\mu(\partial_xu-\frac13\nabla\cdot\mathbf v)$。对 $i=x,j=y$：
$$D_{xy}=\frac12(\partial_yu+\partial_xv),\ \delta_{xy}=0,\ \tau_{xy}=2\mu D_{xy}=\mu(\partial_yu+\partial_xv).$$

物理意义：法向黏应力含体积黏度修正；切应力与速度梯度张量对称部分成正比。

极限与特例：不可压 $\nabla\cdot\mathbf v=0$，$\tau_{xx}=2\mu\partial_xu$。

典型应用：管流、边界层、润滑。

**知识点3　纳维—斯托克斯方程推导**

定义：牛顿流体动量方程
$$\rho\frac{D\mathbf v}{Dt}=-\nabla p+\nabla\cdot\boldsymbol\tau+\mathbf f.$$

推导：动量方程 $\rho Dv_i/Dt=\partial_j\sigma_{ij}+f_i$。牛顿流体总应力
$$\sigma_{ij}=-p\delta_{ij}+2\mu D_{ij}+\lambda(\nabla\cdot\mathbf v)\delta_{ij}.$$
逐项
$$
\partial_j\sigma_{ij}=-\partial_ip+\mu\partial_j(\partial_jv_i+\partial_iv_j)+\lambda\partial_i(\partial_jv_j).
$$
当 $\mu,\lambda$ 为常数时
$$\partial_j\sigma_{ij}=-\partial_ip+\mu\nabla^2v_i+\mu\partial_i(\partial_jv_j)+\lambda\partial_i(\partial_jv_j)
=-\partial_ip+\mu\nabla^2v_i+(\mu+\lambda)\partial_i(\nabla\cdot\mathbf v).$$
$\therefore$ 矢量形式
$$\boxed{\ \rho\frac{D\mathbf v}{Dt}=-\nabla p+\mu\nabla^2\mathbf v+(\mu+\lambda)\nabla(\nabla\cdot\mathbf v)+\mathbf f\ }.$$
不可压时 $\nabla\cdot\mathbf v=0$：
$$\rho\Big(\frac{\partial\mathbf v}{\partial t}+(\mathbf v\cdot\nabla)\mathbf v\Big)=-\nabla p+\mu\nabla^2\mathbf v+\mathbf f.$$

物理意义：非线性对流项与黏性扩散项共同决定流动；非线性能产生涡与湍流。

极限与特例：$\mu=0$ 欧拉方程；$\mathrm{Re}\to0$ 斯托克斯方程；$\mathrm{Re}\to\infty$ 边界层/无粘。

典型应用：CFD 求解、雷诺数依赖、湍流建模。

**知识点4　无量纲化与雷诺数**

定义：特征长度 $L$、速度 $U$、压力 $\rho U^2$。

推导：令 $\mathbf x=L\mathbf x^*,\ \mathbf v=U\mathbf v^*,\ t=\frac LU t^*,\ p=\rho U^2p^*$。代入不可压 N–S：
$$\frac{\rho U^2}{L}\frac{D\mathbf v^*}{Dt^*}=-\frac{\rho U^2}{L}\nabla^*p^*+\frac{\mu U}{L^2}\nabla^{*2}\mathbf v^*.$$
除以 $\rho U^2/L$：
$$\frac{D\mathbf v^*}{Dt^*}=-\nabla^*p^*+\frac1{\mathrm{Re}}\nabla^{*2}\mathbf v^*,\qquad \mathrm{Re}=\frac{\rho UL}{\mu}.$$

物理意义：$\mathrm{Re}$ 表示惯性力与黏性力量级之比。

极限与特例：$\mathrm{Re}\ll1$ 蠕流可忽略对流；$\mathrm{Re}\gg1$ 黏性限于边界层。

典型应用：管流临界雷诺数约 2300；平板转捩 $\mathrm{Re}_x\sim5\times10^5$。

**知识点5　层流精确解**

定义：泊肃叶流（圆管）与库埃特流（平板间）。

推导（泊肃叶）：定常、层流、不可压、圆管 $r$ 方向。N–S 轴向：
$$0=-\frac{dp}{dx}+\mu\frac1r\frac{d}{dr}\Big(r\frac{du}{dr}\Big),$$
边界 $u(r=R)=0$，中心 $du/dr|_{r=0}=0$。两次积分：
$$r\frac{du}{dr}=\frac{1}{2\mu}\frac{dp}{dx}r^2+C_1,$$
$C_1=0$ 由中心条件。再积：
$$u(r)=\frac{1}{4\mu}\frac{dp}{dx}(r^2-R^2).$$
令 $dp/dx=-\Delta p/L$：
$$u(r)=\frac{\Delta p}{4\mu L}(R^2-r^2).$$
流量
$$Q=2\pi\int_0^Ru(r)r\,dr=\frac{\pi R^4\Delta p}{8\mu L}.$$

物理意义：抛物线速度剖面，流量与半径四次方成正比（泊肃叶定律）。

极限与特例：平板库埃特 $u(y)=\frac{G}{2\mu}(h y-y^2)+U\frac yh$。

典型应用：血液流动、润滑、微流控。

**知识点6　边界层方程**

定义：高 $\mathrm{Re}$ 下黏性效应集中在壁面薄层。

推导量级：设 $x\sim L$，$y\sim\delta\ll L$，$u\sim U$。由连续性 $\partial_xu+\partial_yv=0$ 得 $v\sim U\delta/L$。对 $x$ 动量方程：
$$u\partial_xu+v\partial_yu=-\frac1\rho\partial_xp+\nu(\partial_{xx}u+\partial_{yy}u).$$
量级：$u\partial_xu\sim U^2/L$；$\nu\partial_{yy}u\sim\nu U/\delta^2$。两项同阶得
$$\frac{\nu U}{\delta^2}\sim\frac{U^2}{L}\ \Rightarrow\ \delta\sim\sqrt{\frac{\nu L}{U}}=\frac{L}{\sqrt{\mathrm{Re}}}.$$
忽略 $\partial_{xx}u$（相对 $\partial_{yy}u$ 为 $O(\delta^2/L^2)$），得普朗特边界层方程
$$u\partial_xu+v\partial_yu=-\frac1\rho\frac{dp}{dx}+\nu\partial_{yy}u,\qquad \partial_xu+\partial_yv=0.$$

物理意义：黏性项与惯性项同阶，压力沿法向近似不变。

极限与特例：布拉修斯平板解 $\delta/x\approx5/\sqrt{\mathrm{Re}_x}$；分离时 $\partial_yu|_w=0$。

典型应用：摩阻/传热估算、气动设计。

**知识点7　雷诺平均与湍流应力**

定义：将速度分解为平均与脉动 $v_i=\bar v_i+v_i'$，$\overline{v_i'}=0$。

推导：对 N–S 取时间平均。对流项
$$\overline{v_j\partial_jv_i}=\overline{(\bar v_j+v_j')\partial_j(\bar v_i+v_i')}=\bar v_j\partial_j\bar v_i+\overline{v_j'\partial_jv_i'}.$$
$\therefore$
$$\rho\bar v_j\partial_j\bar v_i=-\partial_i\bar p+\mu\nabla^2\bar v_i-\rho\partial_j\overline{v_i'v_j'}.$$
新增项 $-\rho\overline{v_i'v_j'}$ 为雷诺应力，湍流黏性假设
$$-\rho\overline{v_i'v_j'}\approx\mu_t\Big(\partial_j\bar v_i+\partial_i\bar v_j-\frac23\partial_k\bar v_k\delta_{ij}\Big)-\frac23\rho k\delta_{ij}.$$

物理意义：脉动对平均流动产生附加动量输运。

极限与特例：层流脉动为零；湍流黏性需模型（$k$–$\varepsilon$ 等）。

典型应用：CFD 湍流模拟、大气边界层。

### 第12章 振动理论与模态分析

**知识点1　单自由度无阻尼自由振动**

定义：质量—弹簧系统 $m\ddot x+kx=0$。

推导：牛顿第二定律或能量守恒。设 $x=Ae^{\lambda t}$，$m\lambda^2+k=0$，$\lambda=\pm i\omega_n$，$\omega_n=\sqrt{k/m}$。通解
$$x(t)=A\cos\omega_nt+B\sin\omega_nt=C\cos(\omega_nt-\phi),$$
其中 $C=\sqrt{A^2+B^2}$，周期 $T=2\pi/\omega_n$。

物理/几何意义：简谐振动是相空间中椭圆。

极限与特例：$\omega_n=0$（$k=0$）退化为匀速直线运动。

典型应用：弹簧振子、单摆小振动、结构固有频率估计。

**知识点2　阻尼自由振动**

定义：$m\ddot x+c\dot x+kx=0$，阻尼比 $\zeta=\dfrac{c}{2\sqrt{km}}$。

推导：设 $x=e^{\lambda t}$，特征方程 $m\lambda^2+c\lambda+k=0$，
$$\lambda_{1,2}=-\frac{c}{2m}\pm\sqrt{\Big(\frac{c}{2m}\Big)^2-\frac km}=-\zeta\omega_n\pm\omega_n\sqrt{\zeta^2-1}.$$
- $\zeta<1$ 欠阻尼：$\lambda=-\zeta\omega_n\pm i\omega_d$，$\omega_d=\omega_n\sqrt{1-\zeta^2}$，解
$$x=e^{-\zeta\omega_nt}\big(A\cos\omega_dt+B\sin\omega_dt\big).$$
- $\zeta=1$ 临界阻尼：$x=(A+Bt)e^{-\omega_nt}$。
- $\zeta>1$ 过阻尼：两个实指数衰减。
对数减缩 $\delta=\ln\frac{x_n}{x_{n+1}}=\frac{2\pi\zeta}{\sqrt{1-\zeta^2}}$。

物理意义：阻尼决定衰减速率与振荡与否。

极限与特例：$\zeta=1$ 无振荡且最快回到平衡。

典型应用：减震器、阻尼材料、模态参数识别。

**知识点3　受迫振动与共振**

定义：$m\ddot x+c\dot x+kx=F_0\cos\omega t$。

推导：设稳态解 $x=X\cos(\omega t-\phi)$，代入得
$$X=\frac{F_0/k}{\sqrt{(1-r^2)^2+(2\zeta r)^2}},\qquad \tan\phi=\frac{2\zeta r}{1-r^2},\qquad r=\frac{\omega}{\omega_n}.$$
推导：将 $x=Xe^{i\omega t}$，$F=F_0e^{i\omega t}$ 代入 $(-m\omega^2+ic\omega+k)X=F_0$，
$$X=\frac{F_0}{k-m\omega^2+ic\omega}=\frac{F_0/k}{1-r^2+2i\zeta r}.$$
取模与幅角即得幅频与相频。

物理意义：$r\to1$ 时振幅峰值；$\zeta\to0$ 时峰值发散（共振）。

极限与特例：低频 $r\ll1$，$X\approx F_0/k$；高频 $r\gg1$，$X\approx F_0/(m\omega^2)$。

典型应用：旋转机械临界转速、地震响应、振动台。

**知识点4　多自由度系统与固有频率**

定义：$M\ddot{\mathbf q}+K\mathbf q=\mathbf F$。

推导：设模态解 $\mathbf q=\mathbf A e^{i\omega t}$，得
$$(K-\omega^2M)\mathbf A=0.$$
$\therefore$ 特征值问题 $K\mathbf A=\omega^2M\mathbf A$。$M$ 对称正定、$K$ 对称，$\therefore$ 广义特征值实非负，特征向量 $M$ 正交：
$$\mathbf A^{(i)T}M\mathbf A^{(j)}=0,\quad \mathbf A^{(i)T}K\mathbf A^{(j)}=0,\ i\neq j.$$
可采用质量归一化 $\mathbf A^{(i)T}M\mathbf A^{(i)}=1$，则 $\mathbf A^{(i)T}K\mathbf A^{(i)}=\omega_i^2$。

物理/几何意义：多自由度小振动可解耦为独立的单自由度模态。

极限与特例：重根对应相同固有频率的多个模态；刚性模态 $\omega=0$。

典型应用：结构模态分析、模态叠加法。

**知识点5　模态叠加法**

定义：将响应展开在模态基上 $\mathbf q(t)=\sum_i\eta_i(t)\mathbf A^{(i)}$。

推导：代入 $M\ddot{\mathbf q}+K\mathbf q=\mathbf F$：
$$\sum_i\ddot\eta_iM\mathbf A^{(i)}+\sum_i\eta_iK\mathbf A^{(i)}=\mathbf F.$$
左乘 $\mathbf A^{(j)T}$，利用正交性得
$$\ddot\eta_j+\omega_j^2\eta_j=\mathbf A^{(j)T}\mathbf F\equiv f_j(t).$$
$\therefore$ 每个模态解耦，解为杜哈梅积分
$$\eta_j(t)=\frac1{\omega_j}\int_0^tf_j(\tau)\sin\omega_j(t-\tau)\,d\tau.$$

物理/几何意义：用少数低阶模态近似即可捕捉主要响应。

极限与特例：无阻尼自由振动 $f_j=0$，各模态独立振荡。

典型应用：地震响应、风机叶片、模型降阶。

**知识点6　连续系统振动与弦/梁方程**

定义：弦横向振动 $w(x,t)$，张力 $T$，单位长度质量 $\rho A$。

推导：微元 $dx$ 各受张力，横向合力
$$T\sin\theta|_{x+dx}-T\sin\theta|_x\approx T\frac{\partial^2w}{\partial x^2}dx.$$
由牛顿定律 $\rho A\,dx\,\partial_{tt}w=T\partial_{xx}w\,dx$，$\therefore$
$$\frac{\partial^2w}{\partial t^2}=c^2\frac{\partial^2w}{\partial x^2},\qquad c=\sqrt{T/(\rho A)}.$$
分离变量 $w=X(x)T(t)$，$T''/T=c^2X''/X=-\omega^2$，得 $X''+(\omega/c)^2X=0$。两端固定边界，得特征频率
$$\omega_n=\frac{n\pi}{cL}c=\frac{n\pi}{L}\sqrt{\frac{T}{\rho A}},\qquad n=1,2,\dots$$

物理/几何意义：连续系统有无穷多个模态，频率谱不连续。

极限与特例：悬臂梁 $\omega_n=(\beta_nL)^2\sqrt{EI/(\rho AL^4)}$，$\beta_nL$ 为梁特征值。

典型应用：弦乐器、拉索、太阳能帆板动力学。

**知识点7　隔振与动力吸振器**

定义：振源与设备间安装弹簧阻尼，传递率
$$T_r=\frac{|F_{\text{trans}}|}{|F_{\text{excite}}|}=\sqrt{\frac{1+(2\zeta r)^2}{(1-r^2)^2+(2\zeta r)^2}}.$$

推导：基础简谐位移或力激励下，系统位移幅值 $X$，传递力 $F_t=kX+c\dot X$。由受迫响应 $X=F_0/k\cdot1/(1-r^2+2i\zeta r)$，
$$\frac{F_t}{F_0}=\frac{kX+c(i\omega X)}{F_0}=\frac{1+2i\zeta r}{1-r^2+2i\zeta r}.$$
取模即 $T_r$。$T_r<1$ 要求 $r>\sqrt2$。

物理/几何意义：高频隔振，低频放大；阻尼降低共振峰但高频效果变差。

极限与特例：$r\gg1$，$T_r\approx1/(r^2)$。

典型应用：设备减震、汽车悬架、航天器振动隔离。

### 第13章 稳定性理论基础

**知识点1　平衡与轨道稳定性定义**

定义：设自治系统 $\dot{\mathbf x}=\mathbf f(\mathbf x)$，$\mathbf f(\mathbf x_e)=0$。李雅普诺夫意义下稳定：$\forall\varepsilon>0\ \exists\delta>0$，$\|\mathbf x(0)-\mathbf x_e\|<\delta\Rightarrow\|\mathbf x(t)-\mathbf x_e\|<\varepsilon,\ \forall t\ge0$；渐近稳定额外要求 $t\to\infty$ 时收敛。

推导（线性化稳定性）：设 $\mathbf x=\mathbf x_e+\boldsymbol\xi$，将 $\mathbf f$ 泰勒展开 $\mathbf f(\mathbf x_e)=\mathbf f(\mathbf x_e)+J\boldsymbol\xi+O(\|\boldsymbol\xi\|^2)$，$J=\partial\mathbf f/\partial\mathbf x|_{\mathbf x_e}$。$\therefore$ 线性系统 $\dot{\boldsymbol\xi}=J\boldsymbol\xi$。若 $J$ 全部特征值实部为负，则渐近稳定；若有正实部，则不稳定；临界情形需非线性分析。

物理/几何意义：稳定的平衡点邻近轨道不远离；相图上是稳定的结点/焦点/中心。

极限与特例：哈密顿系统无阻尼时特征值纯虚，需非线性分析。

典型应用：控制、结构屈曲、生态/动力学系统。

**知识点2　李雅普诺夫直接法**

定义：若存在连续可微函数 $V(\mathbf x)$（李雅普诺夫函数）在 $\mathbf x_e$ 的邻域内正定，且 $\dot V=\nabla V\cdot\mathbf f(\mathbf x)\le0$，则平衡稳定；$\dot V<0$（除 $\mathbf x_e$）则渐近稳定。

推导：$\dot V\le0$ 表明能量不增，$\therefore$ 不离开自初始的等值面内。严格地，对任意 $\varepsilon$，取 $B_\varepsilon$，因 $V$ 连续且正定，存在 $c>0$ 使 $\{V<c\}\subset B_\varepsilon$；由 $\dot V\le0$ 轨迹保持在 $\{V\le V(\mathbf x(0))\}$ 内，取 $\delta$ 充分小使 $B_\delta\subset\{V<c\}$ 即得稳定。渐近稳定性再用 LaSalle 不变集定理：轨迹趋于 $\{\dot V=0\}$ 的最大不变集。

物理/几何意义：能量函数单调不增，提供全局或局部稳定性。

极限与特例：保守系统能量 $V$ 即自然候选（中心点稳定但非渐近）。

典型应用：非线性控制、机器人、电源系统。

**知识点3　欧拉压杆屈曲**

定义：细长压杆临界载荷 $P_{cr}$。

推导：小挠度梁挠曲微分方程
$$EI\frac{d^2w}{dx^2}=-Mw=-Pw.$$
即 $w''+k^2w=0$，$k^2=P/(EI)$。简支边界 $w(0)=w(L)=0$：
$$w=A\sin kx,\qquad \sin kL=0\ \Rightarrow\ kL=n\pi.$$
$\therefore$
$$P_{cr}=\frac{n^2\pi^2EI}{L^2},\qquad n=1,2,\dots$$
最低阶 $n=1$：$P_{cr}=\pi^2EI/L^2$。

物理/几何意义：屈曲是几何非线性失稳；临界载荷与抗弯刚度成正、与长度平方成反比。

极限与特例：不同端部约束用有效长度 $\mu L$：两端简支 $\mu=1$，两端固定 $\mu=0.5$，悬臂 $\mu=2$。

典型应用：受压杆件、薄壁结构、蜂窝板。

**知识点4　铁木辛柯—里兹法求临界载荷**

定义：用能量法近似。

推导：压杆总势能
$$\Pi=\frac12\int_0^LEI(w'')^2dx-\frac P2\int_0^L(w')^2dx.$$
取 $w=\sum a_n\phi_n$，代入 $\delta\Pi=0$ 得广义特征值 $K\mathbf a=P\,G\mathbf a$，其中
$$K_{mn}=\int EI\phi_m''\phi_n''\,dx,\qquad G_{mn}=\int\phi_m'\phi_n'\,dx.$$
最小特征值 $P_{cr}$。

物理/几何意义：屈曲是刚度竞争（弯曲能 vs 轴向载荷做功）。

极限与特例：$\phi=\sin(\pi x/L)$ 一步得精确解。

典型应用：复杂边界、变截面杆。

**知识点5　分岔与极限点**

定义：系统参数变化时平衡解数目或稳定性发生变化，称为分岔。

推导（鞍结分岔）：一维范式
$$\dot x=\mu-x^2.$$
平衡 $x=\pm\sqrt\mu$。$\mu<0$ 无实平衡；$\mu>0$ 两个平衡。稳定性：$f'(x)=-2x$，$\therefore x=+\sqrt\mu$ 稳定，$x=-\sqrt\mu$ 不稳定。$\therefore \mu=0$ 处鞍结分岔。

推导（叉式分岔）：$\dot x=\mu x-x^3$，$x=0$ 及 $x=\pm\sqrt\mu$。$f'(x)=\mu-3x^2$。$\mu<0$ 时 $x=0$ 稳定；$\mu>0$ 时 $x=0$ 不稳定，$\pm\sqrt\mu$ 稳定。

物理/几何意义：分岔图上稳定性交换，对应结构屈曲/失稳的跳变。

极限与特例：霍普夫分岔产生极限环。

典型应用：结构后屈曲、颤振、非线性振动。

**知识点6　线性稳定性与特征值判据**

定义：线性系统 $\dot{\mathbf x}=A\mathbf x$，解 $\mathbf x(t)=e^{At}\mathbf x(0)$。

推导：若 $A$ 可对角化 $A=V\Lambda V^{-1}$，则 $e^{At}=Ve^{\Lambda t}V^{-1}$，
$$\mathbf x(t)=\sum_ic_i\mathbf v_ie^{\lambda_it}.$$
$\therefore$ 稳定当且仅当 $\mathrm{Re}\,\lambda_i\le0$，且零实部特征值对应约当块为 $1\times1$（否则多项式增长）。渐近稳定当且仅当 $\mathrm{Re}\,\lambda_i<0$。
Routh–Hurwitz：对特征多项式 $a_0s^n+\dots+a_n$，稳定当且仅当 Hurwitz 行列式全正。

物理/几何意义：模态衰减/增长速率由特征值实部决定，频率由虚部决定。

极限与特例：二阶系统 $s^2+2\zeta\omega_ns+\omega_n^2$，稳定当 $\zeta,\omega_n>0$。

典型应用：控制系统、转子动力、飞行稳定性。

**知识点7　结构后屈曲与 imperfection 敏感性**

定义：屈曲后载荷—挠度关系。

推导（Koiter 展开）：将总势能展为振幅 $a$ 与载荷参数 $\lambda$：
$$\Pi(a,\lambda)=\frac12K_2a^2+\frac13K_3a^3+\frac14K_4a^4-\lambda\frac12Pa^2+\cdots$$
平衡 $\partial\Pi/\partial a=0$：
$$K_2a+K_3a^2+K_4a^3-\lambda Pa=0.$$
$K_3\neq0$ 时 $a$ 与 $\lambda$ 呈不对称分岔（限点），初始缺陷 $\bar a$ 使临界载荷降低：
$$\frac{\Delta\lambda}{\lambda_{cr}}\propto-\big(\bar a\big)^{2/3}.$$
$\therefore$ 缺陷敏感结构（如轴压圆柱壳）临界载荷可大幅降低。

物理/几何意义：缺陷使完美分岔退化为极限点，稳定性差。

极限与特例：$K_3=0$ 且 $K_4>0$ 超临界稳定；$K_4<0$ 亚临界不稳定。

典型应用：轴压圆柱壳、加强筋板、薄壁结构设计 knockdown 因子。

### 第14章 空气动力学基础

**知识点1　势流与流函数**

定义：不可压无旋流，$\mathbf v=\nabla\phi$，$\nabla\cdot\mathbf v=0$，得拉普拉斯方程
$$\nabla^2\phi=0.$$
二维不可压可引入流函数 $\psi$：
$$u=\frac{\partial\psi}{\partial y},\qquad v=-\frac{\partial\psi}{\partial x}.$$

推导：无旋 $\nabla\times\mathbf v=0$ 在单连通区域保证 $\phi$ 存在；不可压代入得 $\nabla^2\phi=0$。对 $\psi$：连续性 $u_x+v_y=0$ 是 $(-v,u)$ 无散的充分条件，$\therefore$ 存在 $\psi$ 使 $u=\partial_y\psi,\ v=-\partial_x\psi$。流线为 $\psi=$ const；等势线 $\phi=$ const 与流线正交。

物理/几何意义：复势 $W(z)=\phi+i\psi$，$dW/dz=u-iv$。

极限与特例：圆柱绕流复势 $W=U(z+R^2/z)$。

典型应用：翼型低速绕流、保角映射。

**知识点2　伯努利方程与压力系数**

定义：不可压势流任意点
$$p+\frac12\rho v^2=p_\infty+\frac12\rho U_\infty^2,$$
压力系数
$$C_p=\frac{p-p_\infty}{\tfrac12\rho U_\infty^2}=1-\frac{v^2}{U_\infty^2}.$$

推导：无旋时伯努利常数全流场相同（见第10章），$\therefore$ 任意两点成立。代入 $C_p$ 定义即得。

物理意义：速度大处压力小；$C_p$ 无量纲且 $C_p\in[-?,1]$，滞止 $C_p=1$。

极限与特例：不可压极限 $\mathrm{Ma}\to0$；可压修正 $C_p=C_{p,0}/\sqrt{1-\mathrm{Ma}_\infty^2}$（Prandtl–Glauert）。

典型应用：翼型压力分布、风洞测量。

**知识点3　库塔—儒可夫斯基升力定理**

定义：二维不可压势流中单位展长升力
$$L'=\rho U_\infty\Gamma,$$
方向为来流旋转 $90^\circ$（儒可夫斯基方向）。

推导：绕圆柱（带环量）复势
$$W(z)=U_\infty\Big(z+\frac{R^2}{z}\Big)-\frac{i\Gamma}{2\pi}\ln z.$$
柱面速度由 $dW/dz$ 求得，$v_\theta(U_\infty,\Gamma)$。由伯努利 $p$，升力
$$L'=-\oint p\sin\theta\,R\,d\theta.$$
代入 $v_\theta=-\big(2U_\infty\sin\theta+\Gamma/(2\pi R)\big)$，积分中与 $\Gamma$ 无关项为零，交叉项给出
$$L'=\rho U_\infty\Gamma.$$

物理/几何意义：环量是升力的根源；库塔条件确定环量使后缘平滑。

极限与特例：$\Gamma=0$ 无升力（达朗贝尔佯谬）；有限翼展需考虑诱导阻力。

典型应用：机翼升力、旋翼、帆船。

**知识点4　薄翼理论**

定义：小迎角、薄翼型，中弧线 $y_c(x)$，弦长 $c$。

推导：将翼型用涡分布 $\gamma(x)$ 代替，儒可夫斯基升力。边界条件（流动相切于中弧线）
$$\frac12\int_0^c\frac{\gamma(\xi)}{\pi(x-\xi)}d\xi=U_\infty\Big(\alpha-\frac{dy_c}{dx}\Big).$$
解为 Glauert 级数。升力系数
$$C_L=2\pi\Big(\alpha+\frac1\pi\int_0^\pi\frac{dy_c}{dx}(1-\cos\theta)\,d\theta\Big),$$
其中 $x=\frac c2(1-\cos\theta)$。对平板 $C_L=2\pi\alpha$，升力线斜率 $dC_L/d\alpha=2\pi$（弧度）。力矩系数关于前缘
$$C_{M,LE}=-\frac{\pi\alpha}{2}.$$

物理/几何意义：升力与迎角线性，焦点在 $1/4$ 弦点。

极限与特例：$\alpha$ 过大发生失速（黏性分离）；可压修正 $C_L/\sqrt{1-\mathrm{Ma}^2}$。

典型应用：机翼初步设计、气动载荷。

**知识点5　可压缩流与马赫数**

定义：$\mathrm{Ma}=v/c$，声速 $c=\sqrt{\gamma RT}=\sqrt{\gamma p/\rho}$。

推导声速：一维小扰动，连续与动量线性化：
$$\frac{\partial\rho'}{\partial t}+\rho_0\frac{\partial u'}{\partial x}=0,\qquad \rho_0\frac{\partial u'}{\partial t}+\frac{\partial p'}{\partial x}=0.$$
等熵 $p'=c^2\rho'$，$c^2=(\partial p/\partial\rho)_s$。消元得波动方程 $\partial_{tt}\rho'=c^2\partial_{xx}\rho'$，$\therefore c=\sqrt{(\partial p/\partial\rho)_s}=\sqrt{\gamma RT}$。

等熵关系
$$\frac{p_0}{p}=\Big(1+\frac{\gamma-1}{2}\mathrm{Ma}^2\Big)^{\frac\gamma{\gamma-1}},\qquad \frac{T_0}{T}=1+\frac{\gamma-1}{2}\mathrm{Ma}^2.$$

物理意义：$\mathrm{Ma}<1$ 亚声、$=1$ 声速、$>1$ 超声；超声时扰动限于马赫锥。

极限与特例：不可压 $\mathrm{Ma}\to0$；高超 $\mathrm{Ma}\gg1$ 时温度剧升。

典型应用：喷管、高超声速飞、激波。

**知识点6　正激波关系**

定义：正激波前后跳跃关系，控制体为激波。

推导：取定常一维绝热，连续/动量/能量：
$$\rho_1u_1=\rho_2u_2,$$
$$p_1+\rho_1u_1^2=p_2+\rho_2u_2^2,$$
$$h_1+\frac12u_1^2=h_2+\frac12u_2^2.$$
对完全气体 $h=c_pT$，等熵能量关系。由连续与动量得
$$u_2-u_1=-\frac{p_2-p_1}{\rho_1u_1}=-\frac{p_2-p_1}{\rho_2u_2}.$$
联立与状态方程解得
$$\mathrm{Ma}_2^2=\frac{1+\frac{\gamma-1}{2}\mathrm{Ma}_1^2}{\gamma\mathrm{Ma}_1^2-\frac{\gamma-1}{2}},$$
$$\frac{p_2}{p_1}=1+\frac{2\gamma}{\gamma+1}(\mathrm{Ma}_1^2-1),$$
$$\frac{\rho_2}{\rho_1}=\frac{(\gamma+1)\mathrm{Ma}_1^2}{2+(\gamma-1)\mathrm{Ma}_1^2}.$$
总压损失
$$\frac{p_{02}}{p_{01}}=\Big[\frac{(\gamma+1)\mathrm{Ma}_1^2}{2+(\gamma-1)\mathrm{Ma}_1^2}\Big]^{\gamma/(\gamma-1)}\Big[\frac{\gamma+1}{2\gamma\mathrm{Ma}_1^2-(\gamma-1)}\Big]^{1/(\gamma-1)}<1.$$

物理意义：激波是熵增的不可逆压缩过程；亚声波后为亚声。

极限与特例：$\mathrm{Ma}_1\to1$ 弱激波，$\mathrm{Ma}_1\to\infty$ 强激波，密度比极限 $(\gamma+1)/(\gamma-1)$。

典型应用：超声进气道、再入、爆炸波。

**知识点7　边界层与摩擦阻力**

定义：壁面剪切应力 $\tau_w=\mu(\partial u/\partial y)_{y=0}$，局部摩擦系数 $c_f=\tau_w/(\tfrac12\rho U^2)$。

推导（布拉修斯平板）：引入流函数与相似变量 $\eta=y\sqrt{U/(\nu x)}$，$\psi=\sqrt{\nu Ux}f(\eta)$。代入边界层方程得布拉修斯方程
$$f'''+\frac12ff''=0,\quad f(0)=f'(0)=0,\ f'(\infty)=1.$$
数值解得 $f''(0)=0.332$。壁面剪应力
$$\tau_w=\mu U\sqrt{\frac U{\nu x}}f''(0),\qquad c_f=\frac{0.664}{\sqrt{\mathrm{Re}_x}}.$$

物理/几何意义：层流摩擦系数随 $x^{-1/2}$ 衰减，边界层增长 $\delta/x\sim5/\sqrt{\mathrm{Re}_x}$。

极限与特例：湍流平板 $c_f\sim0.059\mathrm{Re}_x^{-1/5}$；转捩后阻力增大。

典型应用：气动阻力估算、表面摩阻控制。

### 第15章 航天器轨道与姿态动力学

**知识点1　二体问题与开普勒轨道**

定义：两质点引力 $\mathbf F=-\dfrac{GMm}{r^2}\hat{\mathbf r}$，约化质量 $\mu=G(M+m)\approx GM$。

推导运动方程：设 $\mathbf r=\mathbf r_2-\mathbf r_1$，由牛顿引力
$$\ddot{\mathbf r}=-\frac{\mu}{r^3}\mathbf r.$$
取 $\mathbf r\times$ 两边：$\mathbf r\times\ddot{\mathbf r}=0\Rightarrow \frac{d}{dt}(\mathbf r\times\dot{\mathbf r})=0$，$\therefore \mathbf h=\mathbf r\times\dot{\mathbf r}=$ const（角动量），运动在垂直于 $\mathbf h$ 的平面内。再求 Laplace–Runge–Lenz 矢量
$$\mathbf e=\frac{\dot{\mathbf r}\times\mathbf h}{\mu}-\frac{\mathbf r}{r},$$
可证 $d\mathbf e/dt=0$。$\therefore \mathbf r\cdot\mathbf e=r\frac{h^2}{\mu}-r$，得轨道方程
$$r=\frac{h^2/\mu}{1+e\cos\theta}.$$

物理/几何意义：轨道为圆锥曲线，$e=0$ 圆、$0<e<1$ 椭圆、$e=1$ 抛物、$e>1$ 双曲。

极限与特例：$e=0$ 圆轨道；近日点 $\theta=0$，$r_p=\frac{h^2/\mu}{1+e}$。

典型应用：卫星轨道设计、行星运动。

**知识点2　开普勒三定律**

定义：行星轨道为椭圆，太阳在一焦点；面积速度恒定；周期平方正比于长半轴立方。

推导：第一定律即轨道方程 $r=p/(1+e\cos\theta)$，$p=h^2/\mu$。第二定律：面积速度
$$\frac{dA}{dt}=\frac12r^2\dot\theta=\frac12h=\text{const}.$$
第三定律：椭圆面积 $A=\pi a^2\sqrt{1-e^2}=\pi ab$，周期
$$T=\frac{A}{h/2}=\frac{2\pi ab}{h}.$$
又 $b=a\sqrt{1-e^2}$，$p=a(1-e^2)=h^2/\mu$，$\therefore h^2=\mu a(1-e^2)$，$b^2=a^2(1-e^2)=h^2a/\mu$，$\therefore$
$$T^2=\frac{4\pi^2a^2b^2}{h^2}=\frac{4\pi^2a^3}{\mu}.$$

物理/几何意义：万有引力与向心运动的自洽结果。

极限与特例：$e\to0$，$T^2=4\pi^2a^3/\mu$ 圆。

典型应用：卫星周期、同步轨道半径估算。

**知识点3　轨道能量与活力公式**

定义：比机械能 $\varepsilon=\frac12v^2-\frac\mu r$。

推导：由能量守恒 $\varepsilon=$ const。在近日点，$v_p=h/r_p$，$r_p=p/(1+e)=a(1-e)$，$h=\sqrt{\mu p}$。
$$\varepsilon=\frac12\frac{h^2}{r_p^2}-\frac\mu{r_p}=\frac{\mu p}{2r_p^2}-\frac\mu{r_p}=\frac{\mu}{2r_p}\Big(\frac p{r_p}-2\Big).$$
由 $p=a(1-e^2)$，$r_p=a(1-e)$，得
$$\varepsilon=-\frac{\mu}{2a}.$$
$\therefore$ 活力公式（vis-viva）
$$v^2=\mu\Big(\frac2r-\frac1a\Big).$$

物理/几何意义：$a$ 决定能量，$e$ 决定形状；圆轨道 $v=\sqrt{\mu/r}$。

极限与特例：$\varepsilon=0$ 抛物线逃逸，$v_{esc}=\sqrt{2\mu/r}$。

典型应用：速度增量预算、变轨。

**知识点4　霍曼转移轨道**

定义：共面圆轨道间的双脉冲最优（近似）转移。

推导：初始圆轨道半径 $r_1$，目标 $r_2$。转移椭圆近日点 $r_1$、远日点 $r_2$，$a_t=(r_1+r_2)/2$。由活力公式：
$$v_1=\sqrt{\frac\mu{r_1}},\quad v_{p}=\sqrt{\mu\Big(\frac2{r_1}-\frac1{a_t}\Big)},$$
$$\Delta v_1=v_p-v_1,\qquad v_{a}=\sqrt{\mu\Big(\frac2{r_2}-\frac1{a_t}\Big)},\quad v_2=\sqrt{\frac\mu{r_2}},\qquad \Delta v_2=v_2-v_a.$$
总速度增量 $\Delta v=\Delta v_1+\Delta v_2$。转移时间
$$t=\frac{T_t}{2}=\pi\sqrt{\frac{a_t^3}{\mu}}.$$

物理/几何意义：两次切向脉冲把圆轨道变为与目标相切的椭圆。

极限与特例：$r_2\gg r_1$，$\Delta v\approx(\sqrt2-1)v_1+\ldots$；双椭圆转移在比很大时更省。

典型应用：地球同步轨道入轨、行星际转移。

**知识点5　刚体姿态动力学**

定义：体轴系角速度 $\boldsymbol\omega$，惯量张量 $\mathbf I$，欧拉方程
$$\mathbf I\dot{\boldsymbol\omega}+\boldsymbol\omega\times(\mathbf I\boldsymbol\omega)=\mathbf M.$$

推导：角动量 $\mathbf L=\mathbf I\boldsymbol\omega$。惯性系 $\dot{\mathbf L}=\mathbf M$，体轴系转动
$$\dot{\mathbf L}\big|_i=\dot{\mathbf L}\big|_b+\boldsymbol\omega\times\mathbf L=\mathbf M.$$
$\therefore \mathbf I\dot{\boldsymbol\omega}+\boldsymbol\omega\times(\mathbf I\boldsymbol\omega)=\mathbf M$。主轴系展开为三个欧拉方程。动力学与运动学（欧拉角/四元数）联立。

物理/几何意义：陀螺项引起耦合，惯量分布决定稳定性。

极限与特例：自旋稳定绕最小惯量轴；$\mathbf M=0$ 时角动量守恒。

典型应用：卫星姿态控制、动量轮、重力梯度。

**知识点6　重力梯度稳定**

定义：在中心引力场中，细长航天器因两端引力差异产生恢复力矩。

推导：设两质量 $m$ 位于 $r\pm l\cos\theta$ 方向不同距离，引力差
$$\Delta F=GMm\Big(\frac1{(r-l)^2}-\frac1{(r+l)^2}\Big)\approx\frac{4GMm l}{r^3}.$$
力矩 $N\approx\Delta F\,l\sin\theta\approx\frac{4GMm l^2}{r^3}\sin\theta$。对均匀杆惯量 $I\approx2ml^2$，线性化
$$I\ddot\theta=-\frac{6GMm l^2}{r^3}\theta\ \Rightarrow\ \omega_{lib}^2=\frac{3GM}{r^3}=\frac{3\mu}{r^3}.$$

物理/几何意义：长轴指地是稳定平衡，短轴指地不稳定；天平动周期约为轨道周期。

极限与特例：$\omega_{lib}$ 为轨道角速度的 $\sqrt3$ 倍。

典型应用：无源姿态稳定、空间系绳。

**知识点7　姿态控制基础**

定义：用执行机构（动量轮、喷气、磁力矩器）施加控制力矩。

推导（PID 控制一轴简化）：单轴 $I\ddot\theta=N_c$，PID
$$N_c=-\big(K_p\theta+K_d\dot\theta+K_i\!\int\theta\,dt\big).$$
闭环
$$I\ddot\theta+K_d\dot\theta+K_p\theta+K_i\!\int\theta\,dt=0.$$
对谐振动 LTI 系统，特征方程 $Is^3+K_ds^2+K_ps+K_i=0$，由 Routh–Hurwitz 得稳定性条件 $K_dK_p>IK_i$。

物理/几何意义：比例稳定姿态、微分提供阻尼、积分消除常值偏差。

极限与特例：无反作用飞轮时用喷气力矩；磁力矩器受地磁场限制。

典型应用：卫星三轴稳定、太阳能帆板指向。

### 第16章 测量不确定度与实验数据处理

**知识点1　误差与测量不确定度**

定义：真值 $x_0$，测量值 $x$，误差 $\Delta x=x-x_0$。不确定度 $u(x)$ 表征分散性。A 类不确定度
$$u_A=\frac{s}{\sqrt n},\qquad s=\sqrt{\frac1{n-1}\sum_{i=1}^n(x_i-\bar x)^2}.$$
B 类不确定度 $u_B=\Delta/\sqrt3$（均匀分布）。

推导 $s$ 的无偏性：设 $x_i$ 独立同分布，$\bar x=\frac1n\sum x_i$。则
$$\mathbb E\!\left[\sum(x_i-\bar x)^2\right]=\sum\mathbb E[(x_i-\mu)^2]-n\mathbb E[(\bar x-\mu)^2]=(n-1)\sigma^2.$$
$\therefore \mathbb E[s^2]=\sigma^2$，无偏。$\mathrm{s.e.}(\bar x)=\sigma/\sqrt n$，以其估计得 $u_A=s/\sqrt n$。

物理/几何意义：多次测量平均值的分散性随 $\sqrt n$ 减小。

极限与特例：$n=1$ 不能给 A 类；$n$ 小用 t 分布修正。

典型应用：实验报告、校准。

**知识点2　随机误差与统计分布**

定义：随机误差服从一定分布，常见正态
$$p(x)=\frac1{\sqrt{2\pi}\sigma}\exp\!\Big(-\frac{(x-\mu)^2}{2\sigma^2}\Big).$$

推导（中心极限定理要点）：设 $x_i$ 独立同分布，均值 $\mu$、方差 $\sigma^2$，则
$$\frac{\bar x-\mu}{\sigma/\sqrt n}\xrightarrow{d}\mathcal N(0,1).$$
证明思路：特征函数 $\varphi(t)=\mathbb E[e^{itx}]$，对标准化和取对数并展开到二阶，高阶项消失，得标准正态特征函数 $e^{-t^2/2}$。

物理/几何意义：大量独立小误差叠加趋正态。

极限与特例：有限样本用 t 分布 $\bar x\pm t_{\alpha/2,n-1}s/\sqrt n$。

典型应用：置信区间、异常值判定（$3\sigma$）。

**知识点3　不确定度传播律**

定义：间接量 $y=f(x_1,\dots,x_m)$，各 $x_i$ 不相关时
$$u_c^2(y)=\sum_{i=1}^m\Big(\frac{\partial f}{\partial x_i}\Big)^2u^2(x_i).$$

推导：对 $f$ 作一阶泰勒展开 $y\approx f(\boldsymbol\mu)+\sum_i\frac{\partial f}{\partial x_i}(x_i-\mu_i)$。取方差：
$$\mathrm{Var}(y)=\sum_{i,j}\frac{\partial f}{\partial x_i}\frac{\partial f}{\partial x_j}\mathrm{Cov}(x_i,x_j).$$
不相关时协方差为零，$\therefore$ 得传播律。相关时保留协方差项。

物理/几何意义：灵敏度系数 $\partial f/\partial x_i$ 放大各输入不确定度。

极限与特例：$y=x_1\pm x_2$，$u_c=\sqrt{u_1^2+u_2^2}$；$y=x_1x_2$，相对不确定度平方相加。

典型应用：实验误差分析、工程公差。

**知识点4　最小二乘拟合**

定义：对模型 $y=f(x;\boldsymbol\theta)$，最小化
$$S(\boldsymbol\theta)=\sum_{i=1}^n\big(y_i-f(x_i;\boldsymbol\theta)\big)^2.$$

推导（线性最小二乘）：$f=\mathbf a(\mathbf x)^T\boldsymbol\theta$。极值条件
$$\frac{\partial S}{\partial\boldsymbol\theta}=2\sum_i\big(y_i-\mathbf a_i^T\boldsymbol\theta\big)(-\mathbf a_i)=0.$$
$\therefore$ 正规方程
$$\Big(\sum_i\mathbf a_i\mathbf a_i^T\Big)\boldsymbol\theta=\sum_i\mathbf a_iy_i\ \Leftrightarrow\ A^TA\boldsymbol\theta=A^T\mathbf y.$$
解 $\boldsymbol\theta=(A^TA)^{-1}A^T\mathbf y$。若 $\mathbf y$ 误差同方差且独立，$\mathrm{Cov}(\hat{\boldsymbol\theta})=\sigma^2(A^TA)^{-1}$。

物理/几何意义：最小二乘是 $\mathbf y$ 在 $A$ 列空间上的正交投影。

极限与特例：非线性用高斯—牛顿/Levenberg–Marquardt 迭代；过拟合需正则化。

典型应用：标定曲线、参数辨识（模态、气动）。

**知识点5　回归显著性与相关系数**

定义：相关系数
$$r=\frac{\sum(x_i-\bar x)(y_i-\bar y)}{\sqrt{\sum(x_i-\bar x)^2\sum(y_i-\bar y)^2}}.$$

推导：由柯西—施瓦茨不等式
$$\Big(\sum a_ib_i\Big)^2\le\sum a_i^2\sum b_i^2,$$
取 $a_i=x_i-\bar x,\ b_i=y_i-\bar y$，得 $r^2\le1$。$|r|=1$ 当且仅当点共线。线性拟合中
$$r^2=1-\frac{S_{\text{res}}}{S_{\text{tot}}},$$
$\therefore$ $r^2$ 为被解释的方差比例。

物理/几何意义：$r$ 度量线性相关强度，不反映因果关系。

极限与特例：$r=0$ 线性无关但可能有非线性关系。

典型应用：实验数据相关性分析、模型评价。

**知识点6　量纲分析与相似实验设计**

定义：用 $\pi$ 定理设计缩比实验。

推导：设原型与模型满足几何相似 $L_m=\lambda L_p$。要求 $\mathrm{Re}$ 相等：
$$\frac{V_pL_p}{\nu_p}=\frac{V_mL_m}{\nu_m}\ \Rightarrow\ V_m=V_p\frac{L_p}{L_m}\frac{\nu_m}{\nu_p}=\frac{V_p\nu_m}{\lambda\nu_p}.$$
若同一介质 $\nu_m=\nu_p$，则 $V_m=V_p/\lambda$，试验速度须提高。若同时要求 $\mathrm{Fr}$ 相等，则需
$$\frac{V_m}{\sqrt{gL_m}}=\frac{V_p}{\sqrt{gL_p}}\ \Rightarrow\ V_m=V_p\sqrt\lambda,$$
与上式矛盾，$\therefore$ 完全相似不可同时满足。

物理/几何意义：相似准则约束模型尺寸、速度、介质选择。

极限与特例：部分相似允许松弛次要准则；取 $\lambda_\nu=\lambda^{3/2}$ 可同时满足。

典型应用：风洞/水洞缩比、疲劳试验。

**知识点7　创新训练中的实验方案与数据检验**

定义：从假设到验证的闭环：建模—量纲分析—试验设计—数据采集—不确定度评定—模型修正。

推导（信噪比与试验量）：设效应量 $\Delta$，噪声标准差 $\sigma$，单次信噪比 $SNR=\Delta/\sigma$。$n$ 次平均后 $\sigma_{\bar x}=\sigma/\sqrt n$，要求检测置信度对应
$$n\ge\Big(\frac{z_{\alpha/2}+z_\beta}{\Delta/\sigma}\Big)^2.$$
推导：检验统计量 $Z=\bar x/(\sigma/\sqrt n)$。在显著水平 $\alpha$、功效 $1-\beta$ 下
$$P\big(Z>z_{\alpha/2}\mid \Delta\big)=1-\beta\ \Rightarrow\ \frac{\Delta\sqrt n}{\sigma}-z_{\alpha/2}=z_\beta,$$
$\therefore$ 解出 $n$。

物理/几何意义：试验重复次数由效应/噪声比与显著水平、功效共同决定。

极限与特例：$\Delta/\sigma$ 越小所需 $n$ 越大（平方增长）。

典型应用：项目实验设计、A/B 测试、风洞/模态试验规划。

---
⬆ 上级领域：[[实践与创新]]
