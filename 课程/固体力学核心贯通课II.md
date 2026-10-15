---
课程名称: 固体力学核心贯通课II
学分: 4
总学时: 64
考核方式: 考试
课程属性: 必修
课程性质: 专业基础课程
开课单位: 空天科学与技术学院
领域: 力学
tags:
  - 课程
  - 领域/力学
---

# 固体力学核心贯通课II


### 第1章 张量代数与场论基础

**1.1 指标记号与 Einstein 求和约定**

约定：在单项式中，同一字母下标（或上标）重复出现且只出现两次时，表示对该指标在空间维数 $n=3$ 内求和；重复指标称为哑标，可以任意改名；非重复指标称为自由指标，等式两端自由指标必须一致。

$$\because a_i b_i \equiv \sum_{i=1}^{3} a_i b_i,\qquad \delta_{ii}=\delta_{11}+\delta_{22}+\delta_{33}=3.$$

Kronecker 符号 $\delta_{ij}$ 定义为

$$\delta_{ij}=\begin{cases}1,& i=j,\\ 0,& i\neq j,\end{cases}\qquad \delta_{ij}a_j=a_i.$$

置换符号（Levi-Civita）$\varepsilon_{ijk}$ 定义为

$$\varepsilon_{ijk}=\begin{cases}+1,& (i,j,k)\ \text{为偶排列},\\ -1,& (i,j,k)\ \text{为奇排列},\\ 0,& \text{有重复指标}.\end{cases}$$

**1.2 置换符号与 Kronecker 符号的恒等式**

$$\because \varepsilon_{ijk}\varepsilon_{lmn}=\begin{vmatrix}\delta_{il}&\delta_{im}&\delta_{in}\\ \delta_{jl}&\delta_{jm}&\delta_{jn}\\ \delta_{kl}&\delta_{km}&\delta_{kn}\end{vmatrix},$$

将 $k=n$ 缩并得

$$\therefore \varepsilon_{ijk}\varepsilon_{lmk}=\delta_{il}\delta_{jm}-\delta_{im}\delta_{jl}.$$

再令 $j=m$ 缩并：

$$\because \varepsilon_{ijk}\varepsilon_{ljk}=2\delta_{il},\qquad \varepsilon_{ijk}\varepsilon_{ijk}=6.$$

这些恒等式是把向量叉乘、旋度写成指标形式的工具：

$$(\mathbf{a}\times\mathbf{b})_i=\varepsilon_{ijk}a_jb_k,\qquad (\nabla\times\mathbf{u})_i=\varepsilon_{ijk}\partial_j u_k.$$

**1.3 坐标变换与张量定义**

设旧坐标 $x_i$ 与新坐标 $x_i'$ 间满足 $x_i'=a_{ij}x_j$，其中 $a_{ij}=\partial x_i'/\partial x_j$ 为正交矩阵，满足

$$a_{ik}a_{jk}=\delta_{ij},\qquad a_{ki}a_{kj}=\delta_{ij}.$$

零阶张量（标量）满足 $\phi'(\mathbf{x}')=\phi(\mathbf{x})$；一阶张量（向量）满足 $v_i'=a_{ij}v_j$；二阶张量满足

$$T_{ij}'=a_{ik}a_{jl}T_{kl}.$$

$n$ 阶张量按 $n$ 个正交变换矩阵相乘变换。由变换式与正交性可验证张量加法、外积、缩并均保持张量性质：

$$\because (T_{ij}'+S_{ij}')=a_{ik}a_{jl}(T_{kl}+S_{kl})\Rightarrow\text{和仍为二阶张量};$$

$$\because T_{ij}'S_{j}'=a_{ik}a_{jl}T_{kl}\,a_{jm}S_m=a_{ik}\delta_{lm}T_{kl}S_m=a_{ik}(T_{kl}S_l),$$

故缩并后为一阶张量。

**1.4 商法则（判别张量的充分条件）**

若对任意一阶张量 $v_j$ 都有 $T_{ij}v_j=w_i$，且 $w_i$ 为一阶张量，则 $T_{ij}$ 为二阶张量。

$$\because w_i'=a_{ik}w_k=a_{ik}T_{kj}v_j=a_{ik}T_{kj}a_{lj}v_l',$$

另一方面 $w_i'=T_{ij}'v_j'$，由 $v'$ 任意性得

$$\therefore T_{il}'=a_{ik}a_{lj}T_{kj}.$$

**1.5 主值与不变量**

二阶张量 $T$ 的特征方程为

$$\det(T_{ij}-\lambda\delta_{ij})=0\Rightarrow -\lambda^3+I_1\lambda^2-I_2\lambda+I_3=0,$$

其中

$$I_1=T_{ii},\qquad I_2=\tfrac12\!\left(T_{ii}T_{jj}-T_{ij}T_{ji}\right),\qquad I_3=\det(T_{ij}).$$

在坐标变换下 $\operatorname{tr}T=T_{ii}'=a_{ik}a_{ik}T_{kl}$，而 $a_{ik}a_{ik}=\delta_{kl}$，故

$$\therefore T_{ii}'=T_{kk},$$

不变量是张量内蕴量。对称张量三个主值均为实数，且存在相互正交的主方向；将 $T$ 在主轴系下对角化得 $T=\sum_{i=1}^{3}\lambda_i\,\mathbf{n}^{(i)}\otimes\mathbf{n}^{(i)}$。

**1.6 张量的分解：球量与偏量**

任何二阶张量可唯一分解

$$T_{ij}=\underbrace{\tfrac13T_{kk}\delta_{ij}}_{\text{球量(各向同性)}}+\underbrace{\left(T_{ij}-\tfrac13T_{kk}\delta_{ij}\right)}_{\text{偏量}\ \operatorname{dev}T}.$$

$$\because \operatorname{tr}(\operatorname{dev}T)=\left(T_{ii}-\tfrac13T_{kk}\delta_{ii}\right)=T_{ii}-T_{ii}=0,$$

故球量只改变体积、偏量只改变形状（在应力/应变中即体积变形与畸变）。偏张量第二、第三不变量记 $J_2,J_3$，是塑性屈服判据的核心量。

**1.7 场论微分算子与积分定理**

梯度、散度、旋度在指标记号下为

$$(\nabla\phi)_i=\partial_i\phi,\qquad \nabla\!\cdot\!\mathbf{v}=\partial_i v_i,\qquad (\nabla\times\mathbf{v})_i=\varepsilon_{ijk}\partial_j v_k.$$

散度定理（Gauss）与 Stokes 定理：

$$\int_V \partial_i v_i\,dV=\oint_{\partial V}v_i n_i\,dS,\qquad \oint_C v_i\,dx_i=\int_S \varepsilon_{ijk}\partial_j v_k\,n_i\,dS.$$

由 Gauss 定理可把局部平衡微分方程与整体平衡联系：对任意体积取 $v_i=\sigma_{ij}\phi$，

$$\because \int_V\partial_j(\sigma_{ij}\phi)\,dV=\oint_{\partial V}\sigma_{ij}\phi n_j\,dS,$$

$$\therefore \int_V\partial_j\sigma_{ij}\,\phi\,dV=\oint_{\partial V}t_i\phi\,dS-\int_V\sigma_{ij}\partial_j\phi\,dV,$$

这正是虚功原理的场论基础。

**1.8 各向同性张量与四阶张量**

若张量在任意正交变换下不变，则为各向同性张量。二阶各向同性张量必为 $\alpha\delta_{ij}$。四阶弹性张量若各向同性，必可写成

$$C_{ijkl}=\lambda\delta_{ij}\delta_{kl}+\mu\left(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}\right),$$

其中 $\lambda,\mu$ 为两个独立常数。$\because$ 四阶各向同性张量空间维数为 2，基底为 $\delta_{ij}\delta_{kl}$ 与 $\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}$，由对称性 $C_{ijkl}=C_{klij}=C_{ijlk}$ 线性组合即得该式。

**1.9 张量微分的两个关键恒等式**

对对称二阶张量 $\sigma$ 与向量 $\mathbf{u}$：

$$\partial_j(\sigma_{ij}u_i)=\sigma_{ij,j}u_i+\sigma_{ij}u_{i,j},$$

对反对称部分 $\omega$ 有 $\sigma_{ij}\omega_{ij}=0$。$\because$ $\sigma$ 对称、$\omega$ 反对称，缩并 $T_{ij}A_{ij}=0$。此性质用于证明应力张量的对称部分与应变能只依赖对称部分，以及小变形中转动不产生应变能。

**1.10 典型应用与量纲**

应力 $\sigma_{ij}$ 量纲 $\mathrm{Pa}=\mathrm{N/m^2}$；应变无量纲；弹性常数 $E,\lambda,\mu,K$ 量纲同应力；能量密度量纲 $\mathrm{J/m^3}$。指标记法的价值在于：使力学守恒律的局部形式统一写成 $\partial_j(\cdot)+\text{源}=0$，并把坐标变换、材料各向异性、数值实现统一到同一代数框架。

---

### 第2章 应力理论

**2.1 面力（traction）与 Cauchy 引理**

设体内任意点处法向为 $\mathbf{n}$ 的面元 $dS$ 上作用面力 $\Delta\mathbf{F}$，定义面力矢量

$$\mathbf{t}^{(n)}=\lim_{\Delta S\to0}\frac{\Delta\mathbf{F}}{\Delta S}=\sigma_{ij}n_j\,\mathbf{e}_i.$$

Cauchy 引理断言：同一点两个外法向相反的面上，$\mathbf{t}^{(-n)}=-\mathbf{t}^{(n)}$，由动量守恒取面元趋于零的极限直接得到。

**2.2 应力张量与 Cauchy 四面体推导**

取一点处由三个坐标面与一斜面 $ABC$（外法向 $\mathbf{n}$）围成的四面体，高为 $h$。面力平衡（忽略体力对 $h$ 的高阶项）：

$$\because \mathbf{t}^{(n)}\,dS+\mathbf{t}^{(-x_1)}\,dS_1+\mathbf{t}^{(-x_2)}\,dS_2+\mathbf{t}^{(-x_3)}\,dS_3=\mathbf{0},$$

而投影关系 $dS_i=n_i\,dS$，$\mathbf{t}^{(-x_i)}=-\mathbf{t}^{(x_i)}=-\sigma_{ji}\mathbf{e}_j$。代入得

$$\therefore \mathbf{t}^{(n)}=\sigma_{ji}n_j\mathbf{e}_i,\qquad t_i^{(n)}=\sigma_{ji}n_j.$$

将 $\sigma_{ij}$ 排列为矩阵，其第 $i$ 行即法向为 $x_i$ 的面上的面力分量，故 $\sigma_{ij}$ 为二阶应力张量。

**2.3 应力平衡微分方程**

对体内任意微元体作动量平衡。体力 $f_i$ 与面力作用，由 Gauss 散度定理：

$$\because 0=\int_V f_i\,dV+\oint_{\partial V}\sigma_{ij}n_j\,dS=\int_V\!\left(f_i+\partial_j\sigma_{ij}\right)dV,$$

由体积任意性得

$$\therefore \partial_j\sigma_{ij}+f_i=0\quad(\text{静态}),\qquad \partial_j\sigma_{ij}+f_i=\rho\ddot u_i\quad(\text{动态}).$$

写成向量形式 $\nabla\!\cdot\!\boldsymbol{\sigma}+\mathbf{f}=\rho\ddot{\mathbf{u}}$。分量展开即

$$\frac{\partial\sigma_{xx}}{\partial x}+\frac{\partial\sigma_{xy}}{\partial y}+\frac{\partial\sigma_{xz}}{\partial z}+f_x=0,\quad\cdots$$

**2.4 应力张量的对称性**

对微元体取关于 $x_3$ 轴的矩平衡。成对剪力构成力偶，力偶矩等于应力乘面积乘力臂：

$$\because \left(\sigma_{12}\,dy\,dz\right)dx-\left(\sigma_{21}\,dx\,dz\right)dy=0,$$

$$\therefore \sigma_{12}=\sigma_{21},\qquad \sigma_{ij}=\sigma_{ji}.$$

物理意义：应力张量只有 6 个独立分量；对称性来源于无体力矩假设，若材料存在体力偶（如极性介质）则应力不再对称。

**2.5 坐标变换与应力转轴公式**

新老坐标系之间 $\sigma_{ij}'=a_{ik}a_{jl}\sigma_{kl}$。对平面应力，设绕 $z$ 轴转 $\theta$：

$$\sigma_{x'x'}=\sigma_{xx}\cos^2\theta+\sigma_{yy}\sin^2\theta+2\sigma_{xy}\sin\theta\cos\theta,$$

$$\sigma_{y'y'}=\sigma_{xx}\sin^2\theta+\sigma_{yy}\cos^2\theta-2\sigma_{xy}\sin\theta\cos\theta,$$

$$\sigma_{x'y'}=\left(\sigma_{yy}-\sigma_{xx}\right)\sin\theta\cos\theta+\sigma_{xy}\left(\cos^2\theta-\sin^2\theta\right).$$

用半角可化为

$$\sigma_{x'x'}=\frac{\sigma_{xx}+\sigma_{yy}}{2}+\frac{\sigma_{xx}-\sigma_{yy}}{2}\cos2\theta+\sigma_{xy}\sin2\theta,$$

$$\sigma_{x'y'}=-\frac{\sigma_{xx}-\sigma_{yy}}{2}\sin2\theta+\sigma_{xy}\cos2\theta.$$

**2.6 主应力与应力不变量**

主平面上剪力为零，$t_i=\sigma n_i$，即

$$\left(\sigma_{ij}-\sigma\delta_{ij}\right)n_j=0.$$

非零解条件为 $\det(\sigma_{ij}-\sigma\delta_{ij})=0$，展开得特征方程

$$-\sigma^3+I_1\sigma^2-I_2\sigma+I_3=0,$$

其中

$$I_1=\sigma_{ii},\quad I_2=\tfrac12\left(\sigma_{ii}\sigma_{jj}-\sigma_{ij}\sigma_{ji}\right),\quad I_3=\det\sigma_{ij}.$$

由于对称，三个主应力 $\sigma_1\ge\sigma_2\ge\sigma_3$ 为实数，主方向正交。物理意义：主应力是面力法向分量的极值，$I_1,I_2,I_3$ 与坐标无关。

**2.7 最大剪应力与 Mohr 圆**

在主应力系下取任意斜面 $\mathbf{n}=(\alpha,\beta,\gamma)$，则

$$t_n=\sigma_1\alpha^2+\sigma_2\beta^2+\sigma_3\gamma^2,\qquad |\mathbf{t}|^2=\sigma_1^2\alpha^2+\sigma_2^2\beta^2+\sigma_3^2\gamma^2,$$

剪切分量平方为 $\tau^2=|\mathbf{t}|^2-t_n^2$。令 $\alpha^2+\beta^2+\gamma^2=1$，把 $\tau^2$ 写成 $\alpha^2$ 等线性式，得三个圆方程：

$$\tau^2+\left(t_n-\frac{\sigma_2+\sigma_3}{2}\right)^2=\left(\frac{\sigma_2-\sigma_3}{2}\right)^2+\alpha^2(\sigma_1-\sigma_2)(\sigma_1-\sigma_3),$$

其余两式由轮换得到。故 $(\!t_n,\tau\!)$ 点落在三个 Mohr 圆围成的区域。最大剪应力为

$$\therefore \tau_{\max}=\frac{\sigma_1-\sigma_3}{2},$$

其作用面法向平分 $\sigma_1$ 与 $\sigma_3$ 方向。

**2.8 偏应力、von Mises 等效应力与八面体剪应力**

分解

$$s_{ij}=\sigma_{ij}-\frac{\sigma_{kk}}{3}\delta_{ij},\qquad J_k\ \text{为}\ s_{ij}\ \text{的不变量}.$$

第二不变量

$$J_2=\frac12 s_{ij}s_{ij}=\frac16\left[(\sigma_1-\sigma_2)^2+(\sigma_2-\sigma_3)^2+(\sigma_3-\sigma_1)^2\right].$$

$\because$ 在主应力系下 $s_i=\sigma_i-\bar\sigma$，$\sum s_i=0$，故 $\sum s_i^2=2J_2$，展开即得上式。von Mises 等效应力

$$\sigma_{\mathrm{eq}}=\sqrt{3J_2}=\sqrt{\tfrac12\left[(\sigma_1-\sigma_2)^2+(\sigma_2-\sigma_3)^2+(\sigma_3-\sigma_1)^2\right]}.$$

八面体面法向 $\mathbf{n}=(\pm1,\pm1,\pm1)/\sqrt3$，其剪应力

$$\tau_{\mathrm{oct}}=\sqrt{\tfrac23 J_2}=\frac{\sqrt2}{3}\sigma_{\mathrm{eq}}.$$

物理意义：$J_2$ 度量畸变能，$\sigma_{\mathrm{eq}}$ 把复杂应力等效为单轴应力。

**2.9 应力边界条件**

在表面 $\partial V$ 上给定外法向 $\mathbf{n}$ 与面力 $\bar{\mathbf{t}}$：

$$t_i^{(n)}=\sigma_{ij}n_j=\bar t_i\quad\text{于}\ \partial V_{\!t}.$$

位移边界条件为 $u_i=\bar u_i$ 于 $\partial V_{\!u}$，$\partial V_{\!u}\cup\partial V_{\!t}=\partial V$，$\partial V_{\!u}\cap\partial V_{\!t}=\varnothing$（混合边值问题）。法向与切向分解：

$$\sigma_{nn}=n_i\sigma_{ij}n_j,\qquad t_{n s}=s_i\sigma_{ij}n_j,$$

其中 $\mathbf{s}\perp\mathbf{n}$。

**2.10 应力圆典型应用：薄壁圆柱压力容器**

筒壁受内压 $p$、内径 $D$、壁厚 $t$。轴向取截面平衡：

$$\because \sigma_x(2\pi a t)=p\pi a^2\Rightarrow \sigma_x=\frac{pa}{2t}.$$

周向取半筒平衡：

$$\because \sigma_\theta(2tL)=p(2aL)\Rightarrow \sigma_\theta=\frac{pa}{t}.$$

径向应力 $\sigma_r$ 由内壁 $-p$ 到外壁 $0$。最大剪应力

$$\therefore \tau_{\max}=\frac{\sigma_\theta-\sigma_r}{2}\approx\frac{pa}{2t}.$$

这解释了压力容器沿纵向开裂而非环向开裂。

---

### 第3章 应变与变形几何

**3.1 变形梯度与位移梯度**

物质点位置映射 $\mathbf{x}=\boldsymbol{\chi}(\mathbf{X},t)$，变形梯度

$$F_{iJ}=\frac{\partial x_i}{\partial X_J},\qquad dx_i=F_{iJ}\,dX_J.$$

位移 $\mathbf{u}=\mathbf{x}-\mathbf{X}$，故

$$F_{iJ}=\delta_{iJ}+\frac{\partial u_i}{\partial X_J}=I+\nabla_X\mathbf{u}.$$

小变形时物质坐标与空间坐标近似不分：$F\approx I+\nabla\mathbf{u}$，$|\nabla\mathbf{u}|\ll1$。

**3.2 Green–Lagrange 应变张量**

比较变形前后微元长度平方：

$$\because ds^2-dS^2=dx_i\,dx_i-dX_I\,dX_I=F_{iJ}F_{iK}\,dX_J\,dX_K-\delta_{JK}\,dX_J\,dX_K,$$

$$\therefore ds^2-dS^2=2E_{JK}\,dX_J\,dX_K,\qquad E_{JK}=\tfrac12\left(F_{iJ}F_{iK}-\delta_{JK}\right).$$

展开 $F=I+\nabla_X\mathbf{u}$：

$$E_{JK}=\tfrac12\left(\partial_J u_K+\partial_K u_J+\partial_J u_i\,\partial_K u_i\right).$$

物理意义：$E$ 为物质坐标下的有限应变，对刚体转动即使 $ds^2=dS^2$ 也给出 $E=0$ 的正确结果（客观性）。

**3.3 小应变张量与转动张量**

当 $|\nabla\mathbf{u}|\ll1$ 时略去二次项，并把物质/空间导数合并：

$$\varepsilon_{ij}=\tfrac12\left(\partial_i u_j+\partial_j u_i\right)=\tfrac12\left(u_{i,j}+u_{j,i}\right).$$

同时定义转动张量

$$\omega_{ij}=\tfrac12\left(u_{i,j}-u_{j,i}\right),\qquad \omega_{ij}=-\omega_{ji},$$

即 $\nabla\mathbf{u}=\boldsymbol{\varepsilon}+\boldsymbol{\omega}$。物理意义：$\boldsymbol{\varepsilon}$ 描述纯变形，$\boldsymbol{\omega}$ 描述刚体转动；转动不产生应变，$\because \varepsilon_{ij}\omega_{ij}=0$。

**3.4 正应变与剪应变的几何解释**

沿 $x$ 方向线元相对伸长

$$\varepsilon_{xx}=\frac{ds-dS}{dS}=\frac{\partial u}{\partial x},$$

两条原正交线元夹角改变量

$$\gamma_{xy}=2\varepsilon_{xy}=\frac{\partial u}{\partial y}+\frac{\partial v}{\partial x}.$$

体积应变

$$\varepsilon_V=\varepsilon_{ii}=\frac{dV-dV_0}{dV_0}=\partial_i u_i=\nabla\!\cdot\!\mathbf{u}.$$

$\because$ 小变形下 Jacobian $J=\det F\approx1+\operatorname{tr}\nabla\mathbf{u}=1+\varepsilon_{ii}$。

**3.5 应变协调方程（Saint-Venant）**

给定 6 个应变求 3 个位移，方程超定，必须满足协调条件。直接对 $\varepsilon_{ij}$ 求导消去 $\mathbf{u}$：

$$\because \varepsilon_{ij,kl}=\tfrac12\left(u_{i,jkl}+u_{j,ikl}\right),$$

轮换指标并相加减：

$$\varepsilon_{ij,kl}+\varepsilon_{kl,ij}-\varepsilon_{ik,jl}-\varepsilon_{jl,ik}=\tfrac12\left(u_{i,jkl}+u_{j,ikl}+u_{k,ijl}+u_{l,kij}-u_{i,klj}-u_{k,ijl}-u_{j,kli}-u_{l,ikj}\right)=0,$$

$$\therefore \varepsilon_{ij,kl}+\varepsilon_{kl,ij}-\varepsilon_{ik,jl}-\varepsilon_{jl,ik}=0.$$

独立方程数为 6 个（其余由对称性重复）。平面问题退化为

$$\frac{\partial^2\varepsilon_{xx}}{\partial y^2}+\frac{\partial^2\varepsilon_{yy}}{\partial x^2}=2\frac{\partial^2\varepsilon_{xy}}{\partial x\partial y}.$$

**3.6 主应变与应变不变量**

应变张量对称，故

$$\det(\varepsilon_{ij}-\varepsilon\delta_{ij})=0\Rightarrow -\varepsilon^3+I_1\varepsilon^2-I_2\varepsilon+I_3=0,$$

主应变 $\varepsilon_1\ge\varepsilon_2\ge\varepsilon_3$。体积应变

$$\varepsilon_V=I_1=\varepsilon_1+\varepsilon_2+\varepsilon_3,$$

偏应变张量 $e_{ij}=\varepsilon_{ij}-\tfrac13\varepsilon_{kk}\delta_{ij}$，$e_{ii}=0$，其第二不变量度量畸变：

$$J_2^{e}=\tfrac12 e_{ij}e_{ij}=\tfrac16\left[(\varepsilon_1-\varepsilon_2)^2+(\varepsilon_2-\varepsilon_3)^2+(\varepsilon_3-\varepsilon_1)^2\right].$$

**3.7 应变转轴与应变花**

平面应变转轴同应力转轴公式（张量同构）：

$$\varepsilon_{x'x'}=\frac{\varepsilon_{xx}+\varepsilon_{yy}}{2}+\frac{\varepsilon_{xx}-\varepsilon_{yy}}{2}\cos2\theta+\varepsilon_{xy}\sin2\theta,$$

$$\varepsilon_{x'y'}=-\frac{\varepsilon_{xx}-\varepsilon_{yy}}{2}\sin2\theta+\varepsilon_{xy}\cos2\theta.$$

工程中 $0^\circ/45^\circ/90^\circ$ 三向应变花由三片应变片读数联立求解 $\varepsilon_{xx},\varepsilon_{yy},\gamma_{xy}$；主应变方向

$$\tan2\theta_p=\frac{\gamma_{xy}}{\varepsilon_{xx}-\varepsilon_{yy}}.$$

**3.8 柱坐标与球坐标下的几何方程**

柱坐标 $(r,\theta,z)$ 位移 $(u_r,u_\theta,u_z)$：

$$\varepsilon_{rr}=\frac{\partial u_r}{\partial r},\qquad \varepsilon_{\theta\theta}=\frac1r\frac{\partial u_\theta}{\partial\theta}+\frac{u_r}{r},\qquad \varepsilon_{zz}=\frac{\partial u_z}{\partial z},$$

$$\gamma_{r\theta}=\frac1r\frac{\partial u_r}{\partial\theta}+\frac{\partial u_\theta}{\partial r}-\frac{u_\theta}{r},\qquad \gamma_{\theta z}=\frac{\partial u_\theta}{\partial z}+\frac1r\frac{\partial u_z}{\partial\theta},\qquad \gamma_{zr}=\frac{\partial u_z}{\partial r}+\frac{\partial u_r}{\partial z}.$$

$\theta\theta$ 分量中的 $u_r/r$ 来源于弧长随半径的变化，是曲率效应。

**3.9 刚体位移与应变的零空间**

若 $\varepsilon_{ij}=0$，则

$$\because u_{i,j}+u_{j,i}=0\Rightarrow u_{i,jk}=0\Rightarrow u_i=a_i+\omega_{ij}x_j,$$

$$\therefore \mathbf{u}=\mathbf{a}+\boldsymbol{\omega}\times\mathbf{x},$$

即只能为平动加刚体转动 $\boldsymbol{\omega}$，共 6 个自由度。这说明应变确定位移需附加约束，也解释了协调方程的必要性。

**3.10 典型应用：纯弯曲梁的应变分布**

纯弯梁平截面假设下中性轴曲率 $1/R$，离中性轴 $y$ 处纤维长度改变

$$\because ds=\left(R+y\right)d\theta,\quad dS=R\,d\theta\Rightarrow \varepsilon_x=\frac{ds-dS}{dS}=\frac{y}{R},$$

$$\therefore \varepsilon_x=-\,y\,\kappa,\qquad \kappa=\frac{d^2 w}{dx^2}.$$

结合 $\varepsilon_x=\sigma_x/E$ 与截面弯矩积分，可推出 $M=EI\kappa$。此即本章几何方程在梁理论中的应用。

---

### 第4章 运动学与守恒定律

**4.1 物质导数与 Reynolds 输运定理**

物质导数

$$\frac{D\phi}{Dt}=\frac{\partial\phi}{\partial t}+v_i\,\partial_i\phi=\frac{\partial\phi}{\partial t}+\mathbf{v}\!\cdot\!\nabla\phi.$$

对任意物理量 $A$ 的体积分：

$$\because \frac{d}{dt}\int_{V(t)}A\,dV=\int_{V}\left[\frac{\partial A}{\partial t}+\partial_i(Av_i)\right]dV,$$

$$\therefore \frac{D}{Dt}\int_{V(t)}A\,dV=\int_{V(t)}\left(\frac{DA}{Dt}+A\,\partial_i v_i\right)dV.$$

这就是 Reynolds 输运定理，用于把体积分的物质导数化为被积函数的物质导数。

**4.2 连续性方程（质量守恒）**

取 $A=\rho$ 并令质量守恒：

$$\because 0=\frac{D}{Dt}\int_V\rho\,dV=\int_V\left(\frac{D\rho}{Dt}+\rho\,\partial_i v_i\right)dV,$$

$$\therefore \frac{\partial\rho}{\partial t}+\partial_i(\rho v_i)=0\quad\Longleftrightarrow\quad \frac{D\rho}{Dt}+\rho\,\nabla\!\cdot\!\mathbf{v}=0.$$

不可压时 $D\rho/Dt=0$，故

$$\nabla\!\cdot\!\mathbf{v}=v_{i,i}=0.$$

**4.3 动量守恒与 Cauchy 运动方程**

对任意控制体：

$$\because \frac{d}{dt}\int_V\rho v_i\,dV=\int_V f_i\,dV+\oint_{\partial V}t_i\,dS=\int_V\!\left(f_i+\partial_j\sigma_{ij}\right)dV,$$

用连续性方程化简左端：

$$\frac{d}{dt}\int_V\rho v_i\,dV=\int_V\rho\frac{Dv_i}{Dt}\,dV,$$

$$\therefore \partial_j\sigma_{ij}+f_i=\rho\frac{Dv_i}{Dt}=:\rho a_i.$$

**4.4 角动量守恒**

对任意点取矩：

$$\because \frac{d}{dt}\int_V\rho\,\varepsilon_{ijk}x_jv_k\,dV=\int_V\varepsilon_{ijk}x_jf_k\,dV+\oint_{\partial V}\varepsilon_{ijk}x_jt_k\,dS,$$

展开面力项并利用散度定理，得 $\varepsilon_{ijk}\sigma_{jk}=0$，即

$$\therefore \sigma_{ij}=\sigma_{ji}.$$

此给出对称性的另一种整体证明，与微元力矩平衡一致。

**4.5 机械能守恒与能量平衡**

对动量方程点乘速度并积分：

$$\because v_i\left(\partial_j\sigma_{ij}+f_i-\rho\dot v_i\right)=0,$$

$$\therefore \int_V\rho v_i\dot v_i\,dV=\int_V v_i f_i\,dV+\oint_{\partial V}v_i t_i\,dS-\int_V\sigma_{ij}\partial_j v_i\,dV.$$

右端最后一项为应力功率。对弹性材料 $\sigma_{ij}\dot\varepsilon_{ij}=\dot W$ 为应变能密度率，故

$$\frac{D}{Dt}\left(K+U\right)=\dot Q_{\text{外}}+\dot Q_{\text{热}},$$

其中 $K=\tfrac12\int\rho|\mathbf{v}|^2 dV$，$U=\int W\,dV$。

**4.6 热力学第一、第二定律与 Clausius–Duhem 不等式**

局部能量平衡

$$\rho\dot e=\sigma_{ij}\dot\varepsilon_{ij}-\partial_i q_i+\rho r,$$

其中 $e$ 比内能、$\mathbf{q}$ 热流、$r$ 热源。第二定律（Clausius–Duhem）

$$\rho\dot s+\partial_i\!\left(\frac{q_i}{T}\right)-\frac{\rho r}{T}\ge0.$$

引入 Helmholtz 自由能 $\psi=e-Ts$，在等温、小变形下合并两式得

$$\sigma_{ij}\dot\varepsilon_{ij}-\rho\dot\psi\ge0,$$

这是推导本构方程（第5章）的核心不等式。

**4.7 Navier–Stokes 方程推导**

对 Newton 流体，应力张量取

$$\sigma_{ij}=-p\delta_{ij}+2\mu D_{ij}+\lambda D_{kk}\delta_{ij},\qquad D_{ij}=\tfrac12(v_{i,j}+v_{j,i}).$$

代入动量方程 $v$：

$$\because \partial_j\sigma_{ij}=-p_{,i}+\mu\left(\partial_j v_{i,j}+\partial_j v_{j,i}\right)+\lambda\partial_i(\partial_k v_k),$$

$$\therefore \rho\frac{Dv_i}{Dt}=f_i-p_{,i}+\mu\nabla^2v_i+(\lambda+\mu)\partial_i(\partial_k v_k).$$

不可压 $\partial_kv_k=0$ 时

$$\rho\!\left(\frac{\partial\mathbf{v}}{\partial t}+\mathbf{v}\!\cdot\!\nabla\mathbf{v}\right)=-\nabla p+\mu\nabla^2\mathbf{v}+\mathbf{f}.$$

物理意义：左端惯性，右端压力、黏性耗散与体积力。

**4.8 涡量与黏性耗散**

取旋度得涡量方程（不可压）：

$$\because \nabla\times(\mathbf{v}\!\cdot\!\nabla\mathbf{v})=\mathbf{v}\!\cdot\!\nabla\boldsymbol{\omega}-\boldsymbol{\omega}\!\cdot\!\nabla\mathbf{v},$$

$$\therefore \frac{D\boldsymbol{\omega}}{Dt}=\boldsymbol{\omega}\!\cdot\!\nabla\mathbf{v}+\nu\nabla^2\boldsymbol{\omega},\qquad \nu=\frac{\mu}{\rho}.$$

第一项是涡线拉伸，第二项是黏性扩散。黏性耗散率

$$\Phi=2\mu D_{ij}D_{ij}\ge0,$$

$\because D_{ij}$ 对称且 $\mu>0$。

**4.9 分析力学形式：Lagrange 方程**

广义坐标 $q^\alpha$，拉格朗日函数 $L=T-V$，作用量

$$S=\int_{t_1}^{t_2}L\left(q,\dot q,t\right)dt.$$

对 $S$ 取驻值 $\delta S=0$：

$$\because \delta S=\int\left(\frac{\partial L}{\partial q^\alpha}\delta q^\alpha+\frac{\partial L}{\partial\dot q^\alpha}\delta\dot q^\alpha\right)dt=\int\left(\frac{\partial L}{\partial q^\alpha}-\frac{d}{dt}\frac{\partial L}{\partial\dot q^\alpha}\right)\delta q^\alpha\,dt+\left[\frac{\partial L}{\partial\dot q^\alpha}\delta q^\alpha\right]_{t_1}^{t_2},$$

由端点变分为零、$\delta q^\alpha$ 任意：

$$\therefore \frac{d}{dt}\frac{\partial L}{\partial\dot q^\alpha}-\frac{\partial L}{\partial q^\alpha}=0.$$

**4.10 Hamilton 原理与正则方程**

定义广义动量 $p_\alpha=\partial L/\partial\dot q^\alpha$，Hamilton 量

$$H=p_\alpha\dot q^\alpha-L.$$

对 $H$ 取全微分并用 Lagrange 方程：

$$\because dH=\dot q^\alpha dp_\alpha-\dot p_\alpha dq^\alpha,\qquad \dot p_\alpha=\frac{\partial L}{\partial q^\alpha},$$

$$\therefore \dot q^\alpha=\frac{\partial H}{\partial p_\alpha},\qquad \dot p_\alpha=-\frac{\partial H}{\partial q^\alpha}.$$

物理意义：$H$ 通常是总能量，正则方程把 $2n$ 个一阶方程写成辛结构。弹性动力学、结构振动、波传播均可用此框架导出。

---

### 第5章 本构方程与线弹性

**5.1 广义 Hooke 定律与四阶弹性张量**

小变形、等温、线弹性下应力与应变成线性齐次关系：

$$\sigma_{ij}=C_{ijkl}\,\varepsilon_{kl},$$

$C_{ijkl}$ 为四阶弹性刚度张量。由 $\sigma$ 的对称性得 $C_{ijkl}=C_{jikl}$；由 $\varepsilon$ 的对称性得 $C_{ijkl}=C_{ijlk}$。Voigt 记法把 $(11,22,33,23,31,12)$ 缩并为 $(1,2,3,4,5,6)$：

$$\{\sigma\}=[C]\{\varepsilon\},\qquad [C]\ \text{为}\ 6\times6.$$

**5.2 弹性张量的主对称性与独立常数**

若存在应变能密度 $W=\tfrac12 C_{ijkl}\varepsilon_{ij}\varepsilon_{kl}$，则

$$\because \sigma_{ij}=\frac{\partial W}{\partial\varepsilon_{ij}}=C_{ijkl}\varepsilon_{kl},$$

对 $W$ 的二阶混合偏导可交换：

$$\therefore C_{ijkl}=C_{klij}.$$

故 $[C]$ 对称，独立常数至多 $21$ 个（三斜晶系）。这一结论把最一般的各向异性弹性体从 $36$ 个常数降到 $21$ 个。

**5.3 各向异性到各向同性的常数递减**

由晶体对称性逐步化简：

$$\underbrace{21}_{\text{三斜}}\to\underbrace{13}_{\text{单斜}}\to\underbrace{9}_{\text{正交}}\to\underbrace{7}_{\text{六方/横观各向同性}}\to\underbrace{5}_{\text{立方}}\to\underbrace{2}_{\text{各向同性}}.$$

各向同性时

$$C_{ijkl}=\lambda\delta_{ij}\delta_{kl}+\mu\left(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}\right),$$

$$\therefore \sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij}.$$

其中 $\lambda,\mu$ 为 Lamé 常数，$\mu$ 即剪切模量 $G$。

**5.4 工程常数与 Lamé 常数的关系**

单轴拉伸 $\sigma_1=\sigma$，其余为零。由 $\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij}$ 先定泊松比：

$$\because \sigma_2=\lambda\varepsilon_{kk}+2\mu\varepsilon_2=0,\quad \varepsilon_{kk}=\varepsilon_1+2\varepsilon_2,$$

$$\therefore \lambda\varepsilon_1+2(\lambda+\mu)\varepsilon_2=0\Rightarrow \varepsilon_2=-\frac{\lambda}{2(\lambda+\mu)}\varepsilon_1=-\nu\varepsilon_1,\qquad \nu=\frac{\lambda}{2(\lambda+\mu)}.$$

再代入轴向分量，并用 $\varepsilon_{kk}=(1-2\nu)\varepsilon_1=\dfrac{\mu}{\lambda+\mu}\varepsilon_1$：

$$\because \sigma_1=\lambda(1-2\nu)\varepsilon_1+2\mu\varepsilon_1=\left[\frac{\lambda\mu}{\lambda+\mu}+2\mu\right]\varepsilon_1=\frac{\mu(3\lambda+2\mu)}{\lambda+\mu}\varepsilon_1,$$

$$\therefore E=\frac{\sigma_1}{\varepsilon_1}=\frac{\mu(3\lambda+2\mu)}{\lambda+\mu},\qquad G=\mu=\frac{E}{2(1+\nu)},$$

$$\lambda=\frac{E\nu}{(1+\nu)(1-2\nu)},\qquad K=\lambda+\frac{2\mu}{3}=\frac{E}{3(1-2\nu)}.$$

以上关系可两两反解，构成 $E,\nu,G,K,\lambda$ 之间的换算表。

**5.5 应变能密度与正定性**

$$W=\tfrac12\sigma_{ij}\varepsilon_{ij}=\tfrac12\lambda\left(\varepsilon_{kk}\right)^2+\mu\,\varepsilon_{ij}\varepsilon_{ij}.$$

分成体积与畸变两部分：

$$W=\underbrace{\tfrac12 K\left(\varepsilon_{kk}\right)^2}_{\text{体积能}}+\underbrace{\mu\, e_{ij}e_{ij}}_{\text{畸变能}},\qquad e_{ij}=\varepsilon_{ij}-\tfrac13\varepsilon_{kk}\delta_{ij}.$$

$\because \lambda+\tfrac23\mu=K$。热力学稳定性要求 $W\ge0$ 且仅当 $\varepsilon=0$ 取零，于是

$$K>0,\qquad \mu>0\Rightarrow -1<\nu<\tfrac12.$$

**5.6 热弹性本构**

自由能取 $\psi=\tfrac12\varepsilon_{ij}C_{ijkl}\varepsilon_{kl}-\beta_{ij}\varepsilon_{ij}\Delta T-\tfrac{c}{2T_0}(\Delta T)^2$。由 $\sigma_{ij}=\partial\psi/\partial\varepsilon_{ij}$：

$$\therefore \sigma_{ij}=C_{ijkl}\varepsilon_{kl}-\beta_{ij}\Delta T.$$

各向同性时 $\beta_{ij}=\left(3\lambda+2\mu\right)\alpha\delta_{ij}$，故

$$\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij}-(3\lambda+2\mu)\alpha\Delta T\,\delta_{ij}.$$

自由膨胀 $\sigma=0$ 给出 $\varepsilon_{kk}=3\alpha\Delta T$，即体积热膨胀系数为 $3\alpha$。

**5.7 柔度与工程形式**

$$[S]=[C]^{-1},\qquad \varepsilon_{ij}=S_{ijkl}\sigma_{kl}.$$

各向同性展开为

$$\varepsilon_{x}=\frac1E\left[\sigma_x-\nu(\sigma_y+\sigma_z)\right]+\alpha\Delta T,\quad \gamma_{xy}=\frac{\sigma_{xy}}{G},\ \dots$$

平面应力 ($\sigma_z=0$)：

$$\varepsilon_x=\frac1E(\sigma_x-\nu\sigma_y),\quad \varepsilon_y=\frac1E(\sigma_y-\nu\sigma_x),\quad \gamma_{xy}=\frac{\sigma_{xy}}{G}.$$

平面应变 ($\varepsilon_z=0$)：

$$\varepsilon_x=\frac{1-\nu^2}{E}\left(\sigma_x-\frac{\nu}{1-\nu}\sigma_y\right),\quad \sigma_z=\nu(\sigma_x+\sigma_y).$$

**5.8 弹性常数的物理界限与极限**

$$\because G=\frac{E}{2(1+\nu)}>0\Rightarrow\nu>-1;\qquad K=\frac{E}{3(1-2\nu)}>0\Rightarrow\nu<\tfrac12.$$

$\nu=0$ 为无侧向收缩材料；$\nu=\tfrac12$ 为不可压橡胶极限（$K\to\infty$）；$\nu=-1$ 为极端 auxetic 上限。$E\to0$ 对应软材料，$G/E=1/[2(1+\nu)]$。

**5.9 本构方程的客观性与材料对称性**

材料对称性由不变群 $G$ 描述：$C_{ijkl}=Q_{ip}Q_{jq}Q_{kr}Q_{ls}C_{pqrs}$ 对所有 $Q\in G$ 成立。各向同性群为全正交群；横观各向同性群为绕某轴旋转群。$\because$ 自由能必须在该群下不变，故 $C$ 是相应不变张量，独立常数由群表示论决定。

**5.10 典型应用：三向应力下体积模量与波速**

静水压 $p$：$\sigma_{ij}=-p\delta_{ij}$，则 $\varepsilon_{kk}=-p/K$。纵波与横波波速（各向同性弹性体）

$$c_L=\sqrt{\frac{\lambda+2\mu}{\rho}}=\sqrt{\frac{E(1-\nu)}{\rho(1+\nu)(1-2\nu)}},\qquad c_T=\sqrt{\frac{\mu}{\rho}}=\sqrt{\frac{E}{2\rho(1+\nu)}}.$$

$\because$ 由第4章波动方程与弹性本构代入可得，$c_L>c_T$ 恒成立。

---

### 第6章 弹性力学边值问题与变分原理

**6.1 Navier 位移方程**

把本构代入平衡方程：

$$\because \partial_j\left(\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij}\right)+f_i=0,\qquad \varepsilon_{ij}=\tfrac12(u_{i,j}+u_{j,i}),$$

$$\therefore \mu\,u_{i,jj}+(\lambda+\mu)\,u_{j,ji}+f_i=0.$$

向量形式

$$\mu\nabla^2\mathbf{u}+(\lambda+\mu)\nabla(\nabla\!\cdot\!\mathbf{u})+\mathbf{f}=\mathbf{0},$$

用 $G,\nu$ 表示为

$$G\nabla^2\mathbf{u}+\frac{G}{1-2\nu}\nabla(\nabla\!\cdot\!\mathbf{u})+\mathbf{f}=\mathbf{0}.$$

**6.2 Beltrami–Michell 应力协调方程**

无体力时，把本构代入协调方程 $\varepsilon_{ij,kl}+\varepsilon_{kl,ij}-\varepsilon_{ik,jl}-\varepsilon_{jl,ik}=0$，利用 $\varepsilon_{ij}=\frac{1+\nu}{E}\sigma_{ij}-\frac{\nu}{E}\sigma_{kk}\delta_{ij}$ 并用平衡方程消去体力项：

$$\therefore \nabla^2\sigma_{ij}+\frac{1}{1+\nu}\left(\sigma_{kk}\right)_{,ij}=0.$$

缩并 $i=j$ 得 $\nabla^2\sigma_{kk}=0$，即第一应力不变量调和。该方程组与平衡方程共同构成应力解法基础。

**6.3 弹性力学边值问题的提法**

求 $\mathbf{u}\in V$ 满足

$$\mu u_{i,jj}+(\lambda+\mu)u_{j,ji}+f_i=0\quad\text{in}\ V,$$

$$u_i=\bar u_i\ \text{on}\ \partial V_u,\qquad \sigma_{ij}n_j=\bar t_i\ \text{on}\ \partial V_t.$$

这属于椭圆型偏微分方程，需 $u_i\in C^2$ 且边界光滑以保证解存在。

**6.4 解的唯一性（Kirchhoff）**

设同一问题有两个解 $\mathbf{u}^{(1)},\mathbf{u}^{(2)}$，令差 $\tilde{\mathbf{u}}=\mathbf{u}^{(1)}-\mathbf{u}^{(2)}$，则 $\tilde{\mathbf{u}}$ 满足齐次方程与齐次边界条件。应变能为

$$\because U=\tfrac12\int_V\sigma_{ij}\tilde\varepsilon_{ij}\,dV=\tfrac12\int_V\left(\sigma_{ij}\tilde u_{i}\right)_{,j}dV-\tfrac12\int_V\sigma_{ij,j}\tilde u_i\,dV,$$

用散度定理与齐次边界：

$$U=\tfrac12\oint_{\partial V}\tilde t_i\tilde u_i\,dS-\tfrac12\int_V f_i^{h}\tilde u_i\,dV=0,$$

而 $U\ge0$ 且 $U=0\Leftrightarrow\tilde\varepsilon_{ij}=0$，故 $\tilde{\mathbf{u}}$ 为刚体位移；在给定位移约束消除刚体位移后

$$\therefore \mathbf{u}^{(1)}=\mathbf{u}^{(2)}.$$

**6.5 虚功原理**

对任意运动学容许的虚位移 $\delta u_i$（满足 $\delta u_i=0$ on $\partial V_u$），有

$$\because \int_V\left(\sigma_{ij,j}+f_i\right)\delta u_i\,dV=0\Rightarrow \int_V\sigma_{ij,j}\delta u_i\,dV+\int_Vf_i\delta u_i\,dV=0,$$

用散度定理：

$$\therefore \int_V\sigma_{ij}\delta\varepsilon_{ij}\,dV=\int_Vf_i\delta u_i\,dV+\oint_{\partial V_t}\bar t_i\delta u_i\,dS.$$

即内虚功等于外虚功。这是所有能量法与有限元法的出发点。

**6.6 最小势能原理**

定义总势能

$$\Pi=\int_V\left(W-f_iu_i\right)dV-\oint_{\partial V_t}\bar t_iu_i\,dS.$$

取一阶变分并用 $\delta W=\sigma_{ij}\delta\varepsilon_{ij}$：

$$\because \delta\Pi=\int_V\sigma_{ij}\delta\varepsilon_{ij}\,dV-\int_Vf_i\delta u_i\,dV-\oint_{\partial V_t}\bar t_i\delta u_i\,dS,$$

由虚功原理 $\delta\Pi=0$。进一步对二阶变分

$$\because \delta^2\Pi=\int_V\delta\sigma_{ij}\delta\varepsilon_{ij}\,dV=\int_V C_{ijkl}\delta\varepsilon_{ij}\delta\varepsilon_{kl}\,dV>0\ \ (\delta\varepsilon\neq0),$$

$$\therefore \Pi\ \text{在真实解处取严格极小}.$$

**6.7 余能与最小余能原理**

余能

$$\Pi_c=\int_V W_c(\sigma)\,dV-\oint_{\partial V_u}\bar u_i\,t_i\,dS,\qquad W_c=\tfrac12 S_{ijkl}\sigma_{ij}\sigma_{kl}.$$

对静力容许的应力场取变分（$\delta\sigma_{ij}n_j=0$ on $\partial V_t$，$\delta\sigma_{ij,j}=0$）：

$$\because \delta\Pi_c=\int_V\varepsilon_{ij}\delta\sigma_{ij}\,dV-\oint_{\partial V_u}\bar u_i\delta t_i\,dS=0,$$

$$\therefore \Pi_c\ \text{取极小}.$$

势能原理给出位移上界（刚度过高），余能原理给出柔度上界，二者共同界定真解。

**6.8 Ritz 法与 Rayleigh–Ritz 法**

设近似解 $u_i=\sum_{m=1}^{N}a_{im}\phi_m(\mathbf{x})$，$\phi_m$ 满足位移边界。代入

$$\Pi\approx\Pi(a_{im}),\qquad \frac{\partial\Pi}{\partial a_{im}}=0\ \ (m=1,\dots,N),$$

得到线性代数方程组 $[K]\{a\}=\{F\}$。$\because$ $\Pi$ 为 $a$ 的二次函数，驻值条件即线性方程组。Rayleigh 法用于振动：$\omega^2=\min\dfrac{\langle u,\mathcal{L}u\rangle}{\langle u,u\rangle}$，即 Rayleigh 商。

**6.9 Hamilton 原理用于弹性动力学**

$$S=\int_{t_1}^{t_2}\left(T-U+W_{\text{ext}}\right)dt.$$

其中 $T=\tfrac12\int_V\rho\dot u_i\dot u_i\,dV$，$U=\int_VW\,dV$。令 $\delta S=0$：

$$\because \delta T=-\int_V\rho\ddot u_i\delta u_i\,dV,\qquad \delta U=\int_V\sigma_{ij}\delta\varepsilon_{ij}\,dV,$$

$$\therefore \int_{t_1}^{t_2}\!\int_V\left(\sigma_{ij,j}+f_i-\rho\ddot u_i\right)\delta u_i\,dVdt+\text{边界项}=0,$$

回到运动方程。物理意义：动力学控制方程是作用量驻值条件。

**6.10 典型应用：用 Ritz 法求简支梁中点挠度**

取 $w(x)=a\sin(\pi x/L)$，$W=\tfrac12EI(w'')^2$，外力 $\delta P$ 作功 $P a$（中点位移为 $a$）。则

$$\Pi=\frac{EI}{2}\int_0^L\left(\frac{\pi^2}{L^2}\right)^2a^2\sin^2\frac{\pi x}{L}dx-Pa=\frac{EI\pi^4}{4L^3}a^2-Pa,$$

$$\because \frac{\partial\Pi}{\partial a}=0\Rightarrow a=\frac{2PL^3}{EI\pi^4}=0.02053\,\frac{PL^3}{EI},$$

与精确解 $PL^3/(48EI)=0.02083PL^3/EI$ 相差不足 $1.5\%$，说明能量法收敛快。

---

### 第7章 能量法

**7.1 应变能与余能**

对线弹性体

$$U=\frac12\int_V\sigma_{ij}\varepsilon_{ij}\,dV=\int_VW\,dV,$$

余能

$$U_c=\frac12\int_V\sigma_{ij}\varepsilon_{ij}\,dV=U\quad(\text{线弹性}).$$

$\because$ 线性本构下 $\sigma$-$\varepsilon$ 图与 $\varepsilon$-$\sigma$ 图同面积，故 $U=U_c$。非线性时一般 $U\neq U_c$。

**7.2 Clapeyron 定理**

对全部外力（比例加载）作功

$$W=\frac12\sum_iP_i\Delta_i,$$

而能量平衡 $W=U$：

$$\therefore \sum_iP_i\Delta_i=2U.$$

即外力功等于应变能的两倍，这是卡氏定理与单位载荷法的前置结论。

**7.3 卡氏第一定理（Castigliano I）**

设 $U=U(\Delta_1,\dots,\Delta_n)$，外力做功 $W=\sum P_i\Delta_i$，总势能 $\Pi=U-\sum P_i\Delta_i$。由 $\partial\Pi/\partial\Delta_i=0$：

$$\because \frac{\partial\Pi}{\partial\Delta_i}=\frac{\partial U}{\partial\Delta_i}-P_i=0,$$

$$\therefore P_i=\frac{\partial U}{\partial\Delta_i}.$$

物理意义：广义力等于应变能对相应广义位移的偏导。

**7.4 卡氏第二定理（Castigliano II）**

以 $U=U(P_1,\dots,P_n)$ 表述，由

$$\because dU=\sum_i\frac{\partial U}{\partial P_i}dP_i=\sum_i\Delta_i\,dP_i,$$

比较系数得

$$\therefore \Delta_i=\frac{\partial U}{\partial P_i}.$$

对线性弹性同样成立；对非线性弹性需用余能 $U_c$：$\Delta_i=\partial U_c/\partial P_i$。物理意义：位移等于应变能对相应载荷的偏导，是结构位移计算最常用的方法。

**7.5 单位载荷法与 Mohr 积分**

把 $\Delta_i$ 写成对广义力的偏导并分布于梁：

$$\because U=\int_0^L\frac{M^2}{2EI}dx\Rightarrow \Delta=\frac{\partial U}{\partial P}=\int_0^L\frac{M}{EI}\frac{\partial M}{\partial P}dx,$$

而 $\partial M/\partial P=m$ 为单位载荷引起的弯矩：

$$\therefore \Delta=\int_0^L\frac{M\,m}{EI}\,dx.$$

一般形式

$$\Delta=\sum\int\frac{N n}{EA}dx+\sum\int\frac{M m}{EI}dx+\sum\int\frac{T t}{GJ}dx+\sum\int k\frac{V v}{GA}dx.$$

**7.6 Betti 互等定理**

对线弹性体，两组载荷 $(P_i,\Delta_i)$ 与 $(P_j,\Delta_j)$。先加第一组再加第二组所作总功与次序无关：

$$\because W=\tfrac12\sum P_i\Delta_i^{(1)}+\tfrac12\sum P_j\Delta_j^{(2)}+\sum P_i\Delta_i^{(j)},$$

交换次序并令两式相等：

$$\therefore \sum_iP_i\Delta_i^{(j)}=\sum_jP_j\Delta_j^{(i)}.$$

即第一组力在第二组力位移上作的功等于第二组力在第一组力位移上作的功。

**7.7 Maxwell 互等定理与柔度矩阵对称性**

在 Betti 定理中取 $P_i=P_j=1$（单位载荷）：

$$\therefore \delta_{ij}=\delta_{ji},$$

$\delta_{ij}$ 为 $j$ 处单位力在 $i$ 处引起的位移。其矩阵形式 $[\delta]=[\delta]^T$ 即柔度矩阵对称，也是刚度矩阵对称 $[K]=[K]^T$ 的能量根源。

**7.8 静不定结构的最小功与力法**

多余约束力 $X_1,\dots,X_m$。协调条件为多余约束处位移为零：

$$\delta_i=\sum_j\delta_{ij}X_j+\delta_{i0}=0,$$

其中 $\delta_{i0}$ 为外载引起的位移。此即力法正则方程。矩阵形式

$$[\delta]\{X\}=-\{\delta_0\},$$

由 Maxwell 定理 $[\delta]$ 对称，可直接求解。

**7.9 势能法求静不定：最小势能**

对静不定结构，用 Ritz 法对总势能取极小自动满足协调。设多余位移参数 $a_j$：

$$\frac{\partial\Pi}{\partial a_j}=0\Rightarrow [K]\{a\}=\{F\},$$

无需事先判断静定基，这是有限元法与位移法的统一之处。

**7.10 典型应用：悬臂梁端部挠度（卡氏定理）**

端部集中力 $P$，$M(x)=-Px$，则

$$U=\int_0^L\frac{P^2x^2}{2EI}dx=\frac{P^2L^3}{6EI},$$

$$\therefore \delta=\frac{\partial U}{\partial P}=\frac{PL^3}{3EI},$$

与直接积分 $EIw''=-Px$ 所得一致。用于验证能量法正确性。

---

### 第8章 平面问题

**8.1 平面应力与平面应变**

平面应力：$\sigma_z=\sigma_{xz}=\sigma_{yz}=0$，适用于薄板；平面应变：$\varepsilon_z=\gamma_{xz}=\gamma_{yz}=0$，适用于长柱体。二者均为二维问题但本构不同：

$$\text{平面应力}:\ \varepsilon_z=-\frac{\nu}{E}(\sigma_x+\sigma_y);\qquad \text{平面应变}:\ \sigma_z=\nu(\sigma_x+\sigma_y).$$

令 $E'=E,\ \nu'=\nu$（平面应力）；$E'=\dfrac{E}{1-\nu^2},\ \nu'=\dfrac{\nu}{1-\nu}$（平面应变），则两者本构形式统一。

**8.2 Airy 应力函数**

无体力时平衡方程

$$\partial_x\sigma_x+\partial_y\sigma_{xy}=0,\qquad \partial_x\sigma_{xy}+\partial_y\sigma_y=0.$$

令

$$\sigma_x=\frac{\partial^2\phi}{\partial y^2},\qquad \sigma_y=\frac{\partial^2\phi}{\partial x^2},\qquad \sigma_{xy}=-\frac{\partial^2\phi}{\partial x\partial y},$$

$$\because \partial_x\sigma_x+\partial_y\sigma_{xy}=\phi_{,xyy}-\phi_{,xyy}=0,\quad \partial_x\sigma_{xy}+\partial_y\sigma_y=-\phi_{,xxy}+\phi_{,xxy}=0,$$

自动满足平衡。$\phi$ 称为 Airy 应力函数。

**8.3 双调和方程推导**

把应力函数代入协调方程 $\varepsilon_{x,yy}+\varepsilon_{y,xx}=2\varepsilon_{xy,xy}$。用平面应力本构

$$\therefore \nabla^4\phi=\frac{\partial^4\phi}{\partial x^4}+2\frac{\partial^4\phi}{\partial x^2\partial y^2}+\frac{\partial^4\phi}{\partial y^4}=0.$$

推导：$\varepsilon_x=(\phi_{,yy}-\nu\phi_{,xx})/E$，$\varepsilon_y=(\phi_{,xx}-\nu\phi_{,yy})/E$，$\varepsilon_{xy}=-(1+\nu)\phi_{,xy}/E$，代入协调并经整理即得双调和方程。平面应变只需换 $E,\nu$，方程形式不变。

**8.4 极坐标下的 Airy 函数与双调和方程**

$$\sigma_r=\frac1r\frac{\partial\phi}{\partial r}+\frac1{r^2}\frac{\partial^2\phi}{\partial\theta^2},\quad \sigma_\theta=\frac{\partial^2\phi}{\partial r^2},\quad \sigma_{r\theta}=-\frac{\partial}{\partial r}\!\left(\frac1r\frac{\partial\phi}{\partial\theta}\right),$$

$$\left(\frac{\partial^2}{\partial r^2}+\frac1r\frac{\partial}{\partial r}+\frac1{r^2}\frac{\partial^2}{\partial\theta^2}\right)^2\phi=0.$$

对轴对称问题 $\phi=\phi(r)$，方程退化为

$$\left(\frac{d^2}{dr^2}+\frac1r\frac{d}{dr}\right)^2\phi=0\Rightarrow\phi=A\ln r+Br^2\ln r+Cr^2+D.$$

**8.5 Kirsch 圆孔应力集中解**

无限大板开半径 $a$ 的圆孔，远场 $\sigma_x=\sigma$，$\sigma_y=0$。利用分离变量，双调和函数取

$$\phi=\phi_0(r)+\phi_2(r)\cos2\theta,\qquad \phi_0=A_0r^2+B_0\ln r,\qquad \phi_2=A_2r^2+B_2r^4+C_2r^{-2}+D_2.$$

$\because$ $r^m\cos n\theta$ 满足 $\nabla^4\phi=0$ 当 $(m^2-n^2)\left[(m-2)^2-n^2\right]=0$，$n=2$ 时 $m\in\{2,-2,4,0\}$。远场要求 $r^4$ 项消失（$B_2=0$），远端 $\sigma_r\to\dfrac{\sigma}{2}(1+\cos2\theta)$ 定出 $A_0=\dfrac{\sigma}{4}$、$A_2=-\dfrac{\sigma}{4}$，再代入孔边条件 $\sigma_r(a)=0$、$\sigma_{r\theta}(a)=0$ 定出 $B_0,C_2,D_2$，得

$$\sigma_r=\frac{\sigma}{2}\!\left(1-\frac{a^2}{r^2}\right)+\frac{\sigma}{2}\!\left(1-\frac{4a^2}{r^2}+\frac{3a^4}{r^4}\right)\cos2\theta,$$

$$\sigma_\theta=\frac{\sigma}{2}\!\left(1+\frac{a^2}{r^2}\right)-\frac{\sigma}{2}\!\left(1+\frac{3a^4}{r^4}\right)\cos2\theta,$$

$$\sigma_{r\theta}=-\frac{\sigma}{2}\!\left(1+\frac{2a^2}{r^2}-\frac{3a^4}{r^4}\right)\sin2\theta.$$

在孔边 $r=a,\ \theta=\pm\pi/2$：

$$\therefore \sigma_\theta(\theta=\pi/2)=3\sigma,\qquad K_t=3.$$

即圆孔应力集中系数为 3，与孔尺寸无关。

**8.6 应力集中与孔边塑性**

$K_t=\sigma_{\max}/\sigma_{\text{nom}}$。对圆孔 $K_t=3$；椭圆孔 $K_t=1+2a/b$，当 $b\to0$（裂纹）时无限增大，是断裂力学的动机。孔边进入塑性后应力重分布，Neuber 公式

$$\sigma_{\max}\varepsilon_{\max}=K_t^2\sigma_{\text{nom}}\varepsilon_{\text{nom}}.$$

**8.7 楔与半平面问题**

对楔（张角 $2\alpha$），取 $\phi=r^2 f(\theta)$ 或 $r^{\lambda+1}$ 形式；半无限平面受集中力采用 Flamant 解：

$$\sigma_r=-\frac{2P}{\pi}\frac{\cos\theta}{r},\qquad \sigma_\theta=\sigma_{r\theta}=0,$$

$$\because \sigma_r\ \text{沿径向分布，}\ \int_{-\pi/2}^{\pi/2}\sigma_r\cos\theta\,r\,d\theta=-P.$$

该解是接触力学 Boussinesq 解的基础。

**8.8 应力函数解法的一般步骤**

(i) 选择 $\phi$ 使其形式满足边界对称性；(ii) 由 $\sigma=\mathcal{D}\phi$ 求应力；(iii) 代入 $\nabla^4\phi=0$ 确定待定函数；(iv) 用边界条件定常数；(v) 由本构求应变、积分求位移并检查多值性。对多连通域，位移单值性要求回路积分

$$\oint du_i=0,$$

否则需对数项修正。

**8.9 典型应用：厚壁圆筒（Lamé 解）**

内外压 $p_i,p_o$，轴对称。由 Lamé 解

$$\sigma_r=\frac{A}{r^2}+B,\qquad \sigma_\theta=-\frac{A}{r^2}+B,$$

边界 $\sigma_r(a)=-p_i,\ \sigma_r(b)=-p_o$ 定常数：

$$A=\frac{a^2b^2(p_o-p_i)}{b^2-a^2},\qquad B=\frac{p_ia^2-p_ob^2}{b^2-a^2},$$

$$\therefore \sigma_\theta-\sigma_r=\frac{2A}{r^2}.$$

最大剪应力在内壁，这决定了厚壁筒的屈服内压。

**8.10 典型应用：无限大板中裂纹的 I 型场（预告）**

对半无限裂纹，Airy 函数可取含 $\sqrt{r}\cos(\theta/2)$ 项（Westergaard），由此在第14章给出 $K_I=\sigma\sqrt{\pi a}$。本章说明平面 Airy 方法可推广到含奇异性的裂纹问题。

---

### 第9章 空间问题与典型解

**9.1 柱体扭转：半逆解法与 Prandtl 应力函数**

等截面柱体绕 $z$ 轴扭转，位移取半逆设

$$u=-\theta z\,y,\qquad v=\theta z\,x,\qquad w=w(x,y),$$

$\theta$ 为单位长度扭转角。代入几何方程得 $\varepsilon_x=\varepsilon_y=\varepsilon_z=\gamma_{xy}=0$，仅

$$\gamma_{xz}=\frac{\partial w}{\partial x}-\theta y,\qquad \gamma_{yz}=\frac{\partial w}{\partial y}+\theta x.$$

由本构：$\sigma_{xz}=G(\partial_xw-\theta y)$，$\sigma_{yz}=G(\partial_yw+\theta x)$，其余应力为零。

**9.2 Prandtl 应力函数与 Poisson 方程**

令

$$\sigma_{xz}=\frac{\partial\phi}{\partial y},\qquad \sigma_{yz}=-\frac{\partial\phi}{\partial x},$$

则平衡方程 $\partial_x\sigma_{xz}+\partial_y\sigma_{yz}=0$ 自动满足。由协调（或应变协调）得

$$\frac{\partial^2\phi}{\partial x^2}+\frac{\partial^2\phi}{\partial y^2}=-2G\theta.$$

$\because$ 对 $\gamma_{xz},\gamma_{yz}$ 消去 $w$：$\partial_y(\gamma_{xz})-\partial_x(\gamma_{yz})=-2\theta$，代入本构即得。

**9.3 边界条件与扭矩公式**

侧表面无面力要求 $\phi=\text{const}$ on $\partial A$，取 $\phi=0$。扭矩

$$\because T=\int_A\left(x\sigma_{yz}-y\sigma_{xz}\right)dA=\int_A\left(-x\phi_{,x}-y\phi_{,y}\right)dA=-\int_A\left(x\phi\right)_{,x}+\left(y\phi\right)_{,y}-2\phi\,dA,$$

用 $\phi=0$ 于边界与散度定理：

$$\therefore T=2\int_A\phi\,dA.$$

**9.4 薄膜比拟**

Poisson 方程 $\nabla^2\phi=-2G\theta$ 与受均匀压力薄膜方程 $\nabla^2 z=-p/S$ 同型。因此：

$$T=2\int_A\phi\,dA\ \longleftrightarrow\ 2\times\text{薄膜体积};\qquad \tau_{\max}\propto\ \text{薄膜最大坡度}.$$

薄膜比拟把扭转问题化为可直观想象的问题：孔洞（空心截面）对应薄膜在该处为常数，从而剪应力为零。

**9.5 椭圆截面与矩形截面扭转**

椭圆 $\dfrac{x^2}{a^2}+\dfrac{y^2}{b^2}=1$，取 $\phi=C\left(\dfrac{x^2}{a^2}+\dfrac{y^2}{b^2}-1\right)$，代入 $\nabla^2\phi=-2G\theta$：

$$C\left(\frac2{a^2}+\frac2{b^2}\right)=-2G\theta\Rightarrow C=-G\theta\frac{a^2b^2}{a^2+b^2},$$

$$T=2\int_A\phi\,dA=\frac{\pi a^3b^3}{a^2+b^2}G\theta,\qquad \tau_{\max}=\frac{2T}{\pi ab^2}\ \text{于短轴端点}.$$

矩形截面用级数解，$\tau_{\max}$ 在长边中点，角点为零。

**9.6 空心薄壁截面：Bredt 公式**

设壁厚 $\delta$，中线围成面积 $A_m$。由薄膜比拟或平衡可得

$$\tau=\frac{T}{2A_m\delta},\qquad \theta=\frac{T}{4A_m^2G}\oint\frac{ds}{\delta},$$

$$\therefore GJ=\frac{4A_m^2}{\oint ds/\delta}.$$

物理意义：闭口薄壁截面扭转刚度远高于开口截面，因剪力流形成闭环。

**9.7 空间轴对称问题的 Love 位移函数**

轴对称位移 $(u_r,0,u_z)$，引入 Love 函数 $\chi(r,z)$：

$$u_r=-\frac{1}{2G}\frac{\partial^2\chi}{\partial r\partial z},\qquad u_z=\frac{1}{2G}\left[2(1-\nu)\nabla^2\chi-\frac{\partial^2\chi}{\partial z^2}\right],$$

控制方程为双调和 $\nabla^4\chi=0$。适用于圆盘、圆柱、半空间受轴对称载荷。

**9.8 Kelvin 集中力解（全空间）**

无限弹性体原点作用集中力 $P\mathbf{e}_z$，位移

$$u_i=\frac{F}{16\pi G(1-\nu)}\frac{1}{r}\left[(3-4\nu)\delta_{iz}+\frac{x_i x_z}{r^2}\right],\qquad r=\sqrt{x^2+y^2+z^2}.$$

$\because$ 控制方程 $(\lambda+\mu)u_{j,ji}+\mu u_{i,jj}+F\delta(\mathbf{x})\delta_{iz}=0$，将 $1/r$ 型奇性解代入并用 $r_{,i}=x_i/r$、$\nabla^2(1/r)=-4\pi\delta(\mathbf{x})$ 配平奇异项即可验证。

该 $1/r$ 奇性即边界元法的基本解；对集中力矩、双力偶均有对应基本解，构成弹性力学 Green 函数体系。

**9.9 应力波速与微结构**

各向同性弹性体中纵波与横波速度由第5章给出。对于含孔洞或夹杂的复合介质，低频有效模量由 Hashin–Shtrikman 界或 Mori–Tanaka 估计。$\because$ 能量法 + 变分可证两相复合材料体积模量满足

$$K_1+\frac{\phi_2}{\dfrac1{K_2-K_1}+\dfrac{3\phi_1}{3K_1+4\mu_1}}\le K\le\cdots,$$

即 Hashin–Shtrikman 上下界。

**9.10 典型应用：扭转轴的截面形状比较**

同为面积 $A$ 的圆、方、矩形截面，$J/A^2$ 越大扭转刚度越高。$\because$ 薄膜比拟中体积/面积比在圆形时最大，故

$$\therefore J_{\text{圆}}>J_{\text{方}}>J_{\text{狭长矩形}}.$$

开口薄壁截面 $J\approx\tfrac13\sum bt^3$，远小于闭口截面，是航空结构采用闭口盒段的原因。

---

### 第10章 薄板与薄壳

**10.1 Kirchhoff–Love 假设与位移场**

薄板（厚度 $h$，挠度 $w$ 满足 $w\ll h\ll$ 面内尺度）假设：直法线保持直线且垂直中面，厚度无伸缩。位移场

$$u=-z\frac{\partial w}{\partial x},\qquad v=-z\frac{\partial w}{\partial y},\qquad w=w(x,y).$$

由此得曲率与扭率

$$\varepsilon_x=-z\frac{\partial^2w}{\partial x^2},\qquad \varepsilon_y=-z\frac{\partial^2w}{\partial y^2},\qquad \gamma_{xy}=-2z\frac{\partial^2w}{\partial x\partial y}.$$

**10.2 板的内力弯矩与扭矩**

$$\because M_x=\int_{-h/2}^{h/2}\sigma_x z\,dz,\quad M_y=\int_{-h/2}^{h/2}\sigma_yz\,dz,\quad M_{xy}=\int_{-h/2}^{h/2}\sigma_{xy}z\,dz,$$

用平面应力本构与曲率关系积分得

$$M_x=-D\!\left(\frac{\partial^2w}{\partial x^2}+\nu\frac{\partial^2w}{\partial y^2}\right),\quad M_y=-D\!\left(\frac{\partial^2w}{\partial y^2}+\nu\frac{\partial^2w}{\partial x^2}\right),$$

$$M_{xy}=-D(1-\nu)\frac{\partial^2w}{\partial x\partial y},\qquad D=\frac{Eh^3}{12(1-\nu^2)}.$$

$D$ 为板的抗弯刚度。

**10.3 板的平衡方程与双调和方程**

对微元板作力矩与横向力平衡：

$$\because \frac{\partial Q_x}{\partial x}+\frac{\partial Q_y}{\partial y}+q=0,\quad Q_x=\frac{\partial M_x}{\partial x}+\frac{\partial M_{xy}}{\partial y},\quad Q_y=\frac{\partial M_{xy}}{\partial x}+\frac{\partial M_y}{\partial y},$$

消去 $Q_x,Q_y$：

$$\therefore D\nabla^4w=q(x,y).$$

无横向载荷 $q=0$ 时 $\nabla^4w=0$，与平面双调和方程同型。

**10.4 边界条件**

简支边 $x=0$：$w=0$，$M_x=0$（即 $\partial^2w/\partial x^2+\nu\partial^2w/\partial y^2=0$）。固支边：$w=0,\ \partial w/\partial x=0$。自由边 Kirchhoff 等效剪力条件：

$$V_x=Q_x+\frac{\partial M_{xy}}{\partial y}=0,\qquad M_x=0.$$

自由边只有两个条件，这是 Kirchhoff 板理论的著名特征。

**10.5 Navier 解（简支矩形板）**

四边简支 $a\times b$ 板，受分布载荷 $q$ 展开为双重 Fourier 级数：

$$w=\sum_{m,n}A_{mn}\sin\frac{m\pi x}{a}\sin\frac{n\pi y}{b},$$

$$q=\sum_{m,n}q_{mn}\sin\frac{m\pi x}{a}\sin\frac{n\pi y}{b}.$$

代入 $D\nabla^4w=q$：

$$\because \nabla^4\left(\sin\frac{m\pi x}{a}\sin\frac{n\pi y}{b}\right)=\pi^4\!\left(\frac{m^2}{a^2}+\frac{n^2}{b^2}\right)^2\sin\frac{m\pi x}{a}\sin\frac{n\pi y}{b},$$

$$\therefore A_{mn}=\frac{q_{mn}}{D\pi^4\left(\dfrac{m^2}{a^2}+\dfrac{n^2}{b^2}\right)^2}.$$

**10.6 均布载荷下简支板的最大挠度**

对均布 $q_0$，$q_{mn}=\dfrac{16q_0}{\pi^2mn}$（$m,n$ 奇），故

$$w_{\max}=\frac{16q_0}{\pi^6D}\sum_{m,n\ \text{odd}}\frac{1}{mn\left(\dfrac{m^2}{a^2}+\dfrac{n^2}{b^2}\right)^2}.$$

方形板 $a=b$ 时 $w_{\max}=0.00406\,q_0a^4/D$，与级数首项 $0.00416\,q_0a^4/D$ 接近。

**10.7 Lévy 解（对边简支）**

设对边 $y=0,b$ 简支，$w=\sum_m Y_m(y)\sin\dfrac{m\pi x}{a}$，代入得 $Y_m$ 的四阶常微分方程：

$$Y_m''''-2\left(\frac{m\pi}{a}\right)^2Y_m''+\left(\frac{m\pi}{a}\right)^4Y_m=\frac{q_m(y)}{D},$$

解为双曲函数加特解，再由 $y=0,b$ 的边界确定常数。适用于一边自由、一边固支等组合。

**10.8 圆板的轴对称弯曲**

圆板轴对称，$w=w(r)$。极坐标下

$$D\!\left(\frac{d^2}{dr^2}+\frac1r\frac{d}{dr}\right)^2w=q(r),$$

即

$$D\frac1r\frac{d}{dr}\!\left[r\frac{d}{dr}\!\left(\frac1r\frac{d}{dr}\!\left(r\frac{dw}{dr}\right)\right)\right]=q(r).$$

弯矩

$$M_r=-D\!\left(\frac{d^2w}{dr^2}+\frac{\nu}{r}\frac{dw}{dr}\right),\qquad M_\theta=-D\!\left(\frac1r\frac{dw}{dr}+\nu\frac{d^2w}{dr^2}\right).$$

**10.9 薄壳的无矩理论与有矩理论**

无矩（薄膜）理论假设弯矩、扭矩、剪力为零，只保留薄膜力 $N_\theta,N_\phi$。回转壳平衡：

$$\frac{d}{d\phi}\left(rN_\phi\right)-r_1N_\theta\cos\phi+r_1r\,p_r=0,$$

$$\because \text{由法向力平衡，} \frac{N_\phi}{r_1}+\frac{N_\theta}{r_2}=p_n.$$

有矩理论保留弯矩与横向剪力，控制方程为 8 阶（轴对称回转壳为 4 阶），边界效应只在壳边附近衰减，衰减长度 $\sim\sqrt{rh}$。

**10.10 典型应用：受内压薄壁球壳与圆柱壳**

球壳半径 $R$、厚 $t$、内压 $p$：薄膜解

$$\because N_\theta=N_\phi=\frac{pR}{2}\Rightarrow\sigma=\frac{pR}{2t}.$$

圆柱壳：环向 $\sigma_\theta=pR/t$，轴向 $\sigma_x=pR/(2t)$，与第2章一致。边界处有矩效应产生附加弯曲应力，最长衰减长度 $\sqrt{Rt}$，是压力容器开孔补强设计的依据。

---

### 第11章 稳定性与屈曲

**11.1 平衡稳定性的能量判据**

总势能 $\Pi(\mathbf{u},\lambda)$，$\lambda$ 为载荷参数。平衡解满足 $\delta\Pi=0$。稳定性由二阶变分判定：

$$\delta^2\Pi>0\ \text{稳定};\qquad \delta^2\Pi=0\ \text{临界};\qquad \delta^2\Pi<0\ \text{不稳定}.$$

$\because$ 势能极小对应稳定平衡，临界时 Hessian 奇异，即

$$\det\!\left(\frac{\partial^2\Pi}{\partial q_i\partial q_j}\right)=0.$$

**11.2 分叉屈曲与极值屈曲**

分叉屈曲：基本平衡路径与相邻路径在临界点相交，临界条件为线性化特征值问题。极值屈曲：路径达到极限载荷后软化，无分叉。缺陷对极值屈曲敏感，对分叉屈曲次之（Koiter 理论）。

**11.3 压杆屈曲（Euler 公式）**

细长压杆，设失稳时挠度 $w(x)$，弯矩 $M=-Pw$，挠曲线

$$EI\frac{d^2w}{dx^2}=-Pw\Rightarrow\frac{d^2w}{dx^2}+k^2w=0,\qquad k^2=\frac{P}{EI}.$$

通解 $w=A\sin kx+B\cos kx$。两端铰支 $w(0)=w(L)=0$：

$$\because B=0,\ A\sin kL=0\Rightarrow\sin kL=0\Rightarrow kL=n\pi,$$

$$\therefore P_{\mathrm{cr}}=\frac{n^2\pi^2EI}{L^2},\qquad n=1\ \text{为最低阶}.$$

对应波形 $w=A\sin(n\pi x/L)$。物理意义：屈曲是特征值问题，临界力与抗弯刚度成正比、与长度平方成反比。

**11.4 不同约束下的有效长度**

统一写

$$P_{\mathrm{cr}}=\frac{\pi^2EI}{(\mu L)^2},$$

其中 $\mu$：两端铰支 $1$，一端固支一端自由 $2$，两端固支 $0.5$，一端固支一端铰支 $0.7$。$\because$ 通过比较各约束下特征方程的最小根可得。

**11.5 大挠度与临界力的高阶修正**

当挠度不小时，挠曲线微分方程

$$EI\frac{d\theta}{ds}=-P\sin\theta,$$

积分得第一类完全椭圆积分，载荷–挠度关系为椭圆积分比。临界附近

$$P=P_{\mathrm{cr}}\!\left(1+\frac{1}{8}\theta_0^2+\cdots\right),$$

说明后屈曲路径轻微上升（稳定后屈曲）。

**11.6 梁柱方程与放大因子**

同时有横向载荷 $q$ 与轴力 $P$：$EIw''''+Pw''=q$。对两端简支、正弦载荷 $q=q_0\sin(\pi x/L)$：

$$w=\frac{q_0}{EI\left(\dfrac{\pi}{L}\right)^4-P\left(\dfrac{\pi}{L}\right)^2}\sin\frac{\pi x}{L}=\frac{w_0}{1-\dfrac{P}{P_{\mathrm{cr}}}},$$

$$\therefore \text{放大因子}\ AF=\frac{1}{1-P/P_{\mathrm{cr}}}.$$

这是二阶效应（P-$\Delta$ 效应）的精确来源，用于框架稳定分析。

**11.7 Rayleigh 商与能量法求临界力**

设试函数 $w$，弯曲应变能

$$U=\frac12\int_0^L EI\left(w''\right)^2dx,\qquad \text{轴向功}\ T=-\frac12\int_0^L P\left(w'\right)^2dx.$$

临界时 $U+T$ 的二阶变分零：

$$P_{\mathrm{cr}}=\frac{\displaystyle\int_0^LEI(w'')^2dx}{\displaystyle\int_0^L(w')^2dx}.$$

$\because$ 对任意试函数该商是 $P_{\mathrm{cr}}$ 的上界，取最小即得。用于复杂边界与变截面杆。

**11.8 板的屈曲**

四边简支矩形板受面内单向压 $N_x$，设 $w=\sum A_{mn}\sin\dfrac{m\pi x}{a}\sin\dfrac{n\pi y}{b}$，屈曲条件为

$$\therefore N_x=\frac{D\pi^2\left(\dfrac{m^2}{a^2}+\dfrac{n^2}{b^2}\right)^2}{\dfrac{m^2}{a^2}}.$$

对 $n=1$，对 $m$ 取极小得临界

$$\left(N_x\right)_{\mathrm{cr}}=\frac{4\pi^2D}{b^2},\qquad \sigma_{\mathrm{cr}}=\frac{4\pi^2D}{b^2h}=\frac{4\pi^2E}{12(1-\nu^2)}\left(\frac{h}{b}\right)^2=k\frac{\pi^2E}{12(1-\nu^2)}\left(\frac{h}{b}\right)^2,\quad k=4.$$

即经典的板屈曲应力公式。$\because$ 对长板取 $n=1$，对 $m$ 优化得最小值恰对应 $k=4$；短边约束、应力梯度会使 $k$ 改变。

**11.9 圆柱壳轴压屈曲与缺陷敏感性**

薄膜应力 $\sigma_x=pR/(2t)$ 或轴向压 $\sigma$。经典临界应力

$$\sigma_{\mathrm{cr}}=\frac{E}{\sqrt{3(1-\nu^2)}}\frac{t}{R},$$

远高于实验值，实验折减因子可达 $0.2$。$\because$ 后屈曲路径不稳定、缺陷敏感，Koiter 理论给出对初始缺陷幅值 $\bar w$ 的强依赖。

**11.10 典型应用：压杆设计的三种失效模式**

短粗杆强度失效 $\sigma=P/A>\sigma_y$；中长杆弹塑性屈曲由切线模量 $E_t$ 控制：

$$P_{\mathrm{cr}}=\frac{\pi^2E_tI}{(\mu L)^2}\quad(\text{切线模量理论}),$$

细长杆弹性 Euler 屈曲。转点由 $P_{\mathrm{cr}}=\sigma_yA$ 确定。这是航空杆件稳定性设计的核心判据。

---

### 第12章 塑性本构与极限分析

**12.1 屈服现象与屈服条件**

屈服：材料由弹性进入塑性，应力点落在屈服面上

$$f\!\left(\sigma_{ij}\right)=0.$$

对金属，静水压力基本不影响屈服，故 $f$ 只依赖偏应力 $s_{ij}$，即 $f(J_2,J_3)=0$。$\because$ 塑性体积不可压，$\varepsilon_{kk}^p=0$。

**12.2 Tresca 屈服准则**

最大剪应力达到临界值：

$$\tau_{\max}=\frac{\sigma_1-\sigma_3}{2}=k=\frac{\sigma_y}{2},$$

$$\therefore \max\!\left\{|\sigma_1-\sigma_2|,|\sigma_2-\sigma_3|,|\sigma_3-\sigma_1|\right\}=\sigma_y.$$

在主应力空间中为六棱柱，$\pi$ 平面上为正六边形。

**12.3 von Mises 屈服准则**

由 $J_2$ 判据

$$\therefore \sqrt{3J_2}=\sigma_y,\qquad J_2=\tfrac12s_{ij}s_{ij}.$$

在主应力空间为圆柱面。$\because$ Mises 圆通过 Tresca 六边形顶点，两准则在纯剪时相差 $\sqrt3/2$（Mises 纯剪屈服 $\tau_y=\sigma_y/\sqrt3$，Tresca 为 $\sigma_y/2$）。

**12.4 硬化模型与后继屈服面**

等向硬化：$f=J_2-\kappa^2(\bar\varepsilon^p)=0$，$\kappa$ 由单轴曲线确定，$\bar\varepsilon^p=\int\sqrt{\tfrac23 d\varepsilon_{ij}^p d\varepsilon_{ij}^p}$。随动硬化：$f=\tfrac12(s_{ij}-\alpha_{ij})(s_{ij}-\alpha_{ij})-\kappa^2=0$，$\alpha_{ij}$ 为背应力，用于循环塑性（Bauschinger 效应）。

**12.5 正交流动法则与 Drucker 公设**

Drucker 公设要求塑性功率非负且塑性应变增量沿屈服面外法向。对稳定材料：

$$d\varepsilon_{ij}^p=d\lambda\frac{\partial f}{\partial\sigma_{ij}},\qquad d\lambda\ge0.$$

对 von Mises：$\partial f/\partial\sigma_{ij}=3s_{ij}/(2\sigma_{\mathrm{eq}})$，故

$$\therefore d\varepsilon_{ij}^p=d\lambda\,s_{ij},\qquad d\varepsilon_{kk}^p=0.$$

即 Prandtl–Reuss 流动法则，塑性应变增量与偏应力同向。

**12.6 一致性条件**

加载过程中应力点保持在屈服面上：$df=0$。对等向硬化

$$\because df=\frac{\partial f}{\partial\sigma_{ij}}d\sigma_{ij}+\frac{\partial f}{\partial\kappa}\frac{\partial\kappa}{\partial\bar\varepsilon^p}d\bar\varepsilon^p=0,$$

$$\therefore d\lambda=\frac{\dfrac{\partial f}{\partial\sigma_{ij}}\,C_{ijkl}\,d\varepsilon_{kl}}{H+\dfrac{\partial f}{\partial\sigma_{mn}}C_{mnpq}\dfrac{\partial f}{\partial\sigma_{pq}}},$$

其中 $H$ 为硬化模量。由此得弹塑性切线模量

$$C_{ijkl}^{ep}=C_{ijkl}-\frac{C_{ijmn}\dfrac{\partial f}{\partial\sigma_{mn}}\dfrac{\partial f}{\partial\sigma_{pq}}C_{pqkl}}{H+\dfrac{\partial f}{\partial\sigma_{rs}}C_{rstu}\dfrac{\partial f}{\partial\sigma_{tu}}}.$$

**12.7 加卸载准则**

$$df=0\ \text{且}\ \frac{\partial f}{\partial\sigma_{ij}}d\sigma_{ij}>0:\ \text{加载（塑性）},$$

$$df<0\ \text{或}\ \frac{\partial f}{\partial\sigma_{ij}}d\sigma_{ij}<0:\ \text{卸载（弹性）}.$$

物理意义：只有应力增量指向屈服面外侧才产生新的塑性流动。

**12.8 理想刚塑性极限分析定理**

设外载由载荷乘子 $\lambda$ 描述，$P=\lambda P_0$。

下限定理由静力容许应力场（满足平衡、屈服 $f\le0$、力的边界）给出 $P^-$，满足 $P_{\text{极限}}\ge P^-$；上限定理由运动容许速度场（满足速度边界、几何协调）由内功率等于外功率给出 $P^+$，满足 $P_{\text{极限}}\le P^+$。$\because$ 虚功率原理：

$$\int_V\sigma_{ij}\dot\varepsilon_{ij}\,dV=\int_{\partial V}t_i v_i\,dS+\int_Vf_iv_i\,dV.$$

**12.9 塑形极限载荷典型解**

两端铰支梁受跨中集中力 $P$，形成塑性铰时截面弯矩达到 $M_p=\sigma_yW_p$。静力法：$M_p=\dfrac{P_uL}{4}\Rightarrow P_u=\dfrac{4M_p}{L}$。机动法：设机构角速度 $\dot\theta$，跨中速度 $\dfrac{L}{2}\dot\theta$，内功率 $2M_p\dot\theta$，外功率 $P_u\dfrac{L}{2}\dot\theta$，令相等得同一结果。二者一致说明限界收敛。

**12.10 安定性与循环载荷**

给定循环载荷区间 $[\lambda^-,\lambda^+]$，若存在不随时间累积塑性的应力状态，则结构安定。Melan 定理：存在自平衡残余应力场 $\rho_{ij}$ 使 $\lambda^-\sigma^{E}+\rho$ 与 $\lambda^+\sigma^{E}+\rho$ 均满足屈服，则安定。$\because$ 由弹性包络与屈服面相交性判定。这是航空发动机热端部件寿命设计基础。

---

### 第13章 黏弹性与流变模型

**13.1 黏弹性现象与本构特征**

黏弹性材料同时具有弹性储能与黏性耗散，响应依赖时间。基本实验：蠕变（恒定应力下应变随时间增加）、应力松弛（恒定应变下应力衰减）、滞回（循环加载能量耗散）。线性黏弹性本构可写为

$$\sigma(t)=\int_{-\infty}^{t}G(t-\tau)\dot\varepsilon(\tau)\,d\tau=\int_{-\infty}^{t}G(t-\tau)\,d\varepsilon(\tau).$$

**13.2 Maxwell 模型（串联）**

弹性弹簧 $E$ 与黏壶 $\eta$ 串联，总应力相等，应变率相加：

$$\because \dot\varepsilon=\dot\varepsilon_e+\dot\varepsilon_v=\frac{\dot\sigma}{E}+\frac{\sigma}{\eta},$$

$$\therefore \dot\sigma+\frac{E}{\eta}\sigma=E\dot\varepsilon.$$

应力松弛（$\varepsilon=\varepsilon_0$ 恒定，$\dot\varepsilon=0$）：

$$\sigma(t)=\sigma_0e^{-t/\tau},\qquad \tau=\frac{\eta}{E}.$$

物理意义：Maxwell 描述松弛与流动，不能描述稳态蠕变约束下的弹性回复。

**13.3 Kelvin–Voigt 模型（并联）**

弹簧与黏壶并联，应变相等，应力相加：

$$\therefore \sigma=E\varepsilon+\eta\dot\varepsilon.$$

蠕变（$\sigma=\sigma_0$ 恒定）：

$$\varepsilon(t)=\frac{\sigma_0}{E}\left(1-e^{-t/\tau}\right),\qquad \tau=\frac{\eta}{E}.$$

物理意义：描述延迟弹性（可回复蠕变），但无瞬时弹性响应（初始应变率为有限而非无穷）。

**13.4 三元件与标准线性固体**

标准线性固体（弹簧 $E_1$ 与 Kelvin 元件串联，或 Maxwell 与弹簧并联）本构为

$$\sigma+\tau_\sigma\dot\sigma=E_0\varepsilon+E_1\tau_\sigma\dot\varepsilon,$$

其蠕变柔量

$$J(t)=\frac1{E_0}+\frac1{E_1}\left(1-e^{-t/\tau}\right),$$

松弛模量

$$G(t)=E_\infty+\left(E_0-E_\infty\right)e^{-t/\tau}.$$

**13.5 Boltzmann 叠加原理**

线性黏弹性材料对历时应变的总响应等于各应变增量响应之和：

$$\therefore \sigma(t)=\int_{-\infty}^{t}G(t-\tau)\,d\varepsilon(\tau)=G(0)\varepsilon(t)+\int_0^t\dot G(t-\tau)\varepsilon(\tau)\,d\tau.$$

物理意义：叠加原理把任意载荷历史化为增量响应的卷积，是线性黏弹性的核心。

**13.6 蠕变柔量与松弛模量的关系**

对单位阶跃应变，$\sigma(t)=G(t)$；对单位阶跃应力，$\varepsilon(t)=J(t)$。二者通过卷积互逆：

$$\therefore \int_0^t G(t-\tau)J(\tau)\,d\tau=t.$$

Laplace 变换下乘积关系

$$\bar G(s)\bar J(s)=\frac{1}{s^2},$$

$\because$ 卷积定理与阶跃函数变换可得。

**13.7 对应原理（弹性–黏弹性类比）**

对边界条件为阶跃或时间的慢变函数，黏弹性解可由弹性解经 Laplace 变换得到：把弹性常数 $E,\nu$ 替换为 $s\bar E(s),\bar\nu(s)$，求得变换域解后再逆变换。$\because$ 黏弹性场方程在 Laplace 域与弹性场方程同构。

**13.8 复模量与动态力学分析**

稳态循环 $\varepsilon=\varepsilon_0e^{i\omega t}$，$\sigma=\sigma_0e^{i(\omega t+\delta)}$，定义

$$E^*(\omega)=\frac{\sigma_0}{\varepsilon_0}e^{i\delta}=E'(\omega)+iE''(\omega),$$

$$\tan\delta=\frac{E''}{E'},$$

$E'$ 储能模量、$E''$ 损耗模量，$\tan\delta$ 为损耗因子。$\because$ 每周期单位体积耗能

$$W_d=\oint\sigma\,d\varepsilon=\pi E''\varepsilon_0^2.$$

**13.9 时温等效原理**

高分子等热流变简单材料满足

$$E(t,T)=E\!\left(\frac{t}{a_T(T)},T_0\right),$$

移位因子 $a_T$ 常由 WLF 方程给出：

$$\log_{10}a_T=-\frac{C_1(T-T_0)}{C_2+(T-T_0)}.$$

物理意义：温度升高等效于时间缩短，可把蠕变数据平移成主曲线。

**13.10 典型应用：黏弹性梁与阻尼结构**

做动载荷下梁的弯曲，用对应原理把 $E$ 换为 $E^*(\omega)$，得复刚度：

$$D^*(\omega)=\frac{E^*(\omega)h^3}{12\left(1-\nu^2\right)}.$$

固有频率由 $\det([K^*]-m\omega^2)=0$ 求出，实部定频率、虚部定模态阻尼。这是航空结构减振与复合材料阻尼设计的基础。

---

### 第14章 断裂力学基础

**14.1 Griffith 能量判据**

含裂纹长 $a$ 的板，释放应变能与新增表面能满足

$$\frac{d\Pi}{da}=\frac{dW_s}{da}.$$

对无限大薄板中心裂纹，远场拉应力 $\sigma$，Griffith 用 Inglis 解给出裂纹释放能

$$U_a=-\frac{\pi\sigma^2a^2}{E'},$$

表面能 $W_s=4a\gamma$（两个裂纹面各两面）。临界条件

$$\because \frac{d}{da}\left(4a\gamma-\frac{\pi\sigma^2a^2}{E'}\right)=0,$$

$$\therefore \sigma_c=\sqrt{\frac{2E'\gamma}{\pi a}},\qquad G_c=2\gamma.$$

物理意义：裂纹扩展由能量释放率与表面能竞争决定。

**14.2 能量释放率 G**

定义

$$G=-\frac{d\Pi}{da}\quad(\text{恒位移或恒载荷对应不同符号约定}),$$

线弹性下

$$\therefore G=\frac{K^2}{E'},\qquad E'=\begin{cases}E,&\text{平面应力},\\ \dfrac{E}{1-\nu^2},&\text{平面应变}.\end{cases}$$

$\because$ Irwin 用裂纹闭合功积分证明 $G=K^2/E'$。

**14.3 裂纹尖端应力场与应力强度因子**

Westergaard 函数 $Z(z)$，对 I 型裂纹

$$Z(z)=\frac{\sigma z}{\sqrt{z^2-a^2}}.$$

令 $z=a+re^{i\theta}$，展开得应力

$$\sigma_{xx}=\frac{K_I}{\sqrt{2\pi r}}\cos\frac\theta2\!\left(1-\sin\frac\theta2\sin\frac{3\theta}2\right),$$

$$\sigma_{yy}=\frac{K_I}{\sqrt{2\pi r}}\cos\frac\theta2\!\left(1+\sin\frac\theta2\sin\frac{3\theta}2\right),$$

$$\sigma_{xy}=\frac{K_I}{\sqrt{2\pi r}}\cos\frac\theta2\sin\frac\theta2\cos\frac{3\theta}2,$$

$$K_I=\sigma\sqrt{\pi a}\quad(\text{无限大板中心裂纹}).$$

其余两型 $K_{II},K_{III}$ 形式类似。

**14.4 K 与 G 的关系推导**

用裂纹面位移与裂纹闭合功：

$$\because G=\lim_{\Delta a\to0}\frac{1}{2\Delta a}\int_0^{\Delta a}\sigma_{yy}(r)\,\delta v(\Delta a-r)\,dr,$$

裂纹张开位移

$$v(r)=\frac{K_I}{G'}\sqrt{\frac{r}{2\pi}}\times2,\quad G'=\frac E{2(1+\nu)},\ \delta v=v(\Delta a-r),$$

$$\therefore G=\frac{K_I^2}{E'},\qquad E'=\frac E{1-\nu^2}\ (\text{平面应变}),\ E'=E\ (\text{平面应力}).$$

**14.5 混合模式与断裂判据**

$K_I,K_{II},K_{III}$ 分别对应张开、滑开、撕开。最大周向应力判据：

$$\theta_0=2\arctan\frac14\!\left[\frac{K_I}{K_{II}}\pm\sqrt{\left(\frac{K_I}{K_{II}}\right)^2+8}\right],$$

临界条件 $K_{I\theta}(\theta_0)=K_{Ic}$。能量判据 $G=G_c$，$G=\dfrac{K_I^2+K_{II}^2}{E'}+\dfrac{K_{III}^2}{2G}$。

**14.6 J 积分：定义与路径无关性**

对二维问题

$$J=\oint_\Gamma\left(W\,dy-t_i\frac{\partial u_i}{\partial x}\,ds\right),\qquad W=\int_0^{\varepsilon}\sigma_{ij}\,d\varepsilon_{ij}.$$

路径无关性证明：对闭合回路 $\Gamma=\Gamma_1+\Gamma_2+\Gamma^++\Gamma^-$，由散度定理

$$\because \oint_\Gamma\left(Wn_x-t_i u_{i,x}\right)ds=\int_A\left[\frac{\partial W}{\partial x}-\partial_j\!\left(\sigma_{ij}u_{i,x}\right)\right]dA=0,$$

对弹性材料 $\partial_xW=\sigma_{ij}\partial_x\varepsilon_{ij}=\sigma_{ij}\partial_j u_{i,x}$，故被积函数为零，回路积分为零，裂纹面无牵引时 $J$ 与路径无关。对小范围屈服，$J=G$。

**14.7 裂纹尖端场 HRR 与 J 控制**

幂硬化材料裂纹尖端应力

$$\sigma_{ij}=\sigma_0\!\left(\frac{J}{\sigma_0\varepsilon_0I_nr}\right)^{\frac{1}{n+1}}\tilde\sigma_{ij}(\theta),\qquad \varepsilon_{ij}=\varepsilon_0\!\left(\frac{J}{\sigma_0\varepsilon_0I_nr}\right)^{\frac{n}{n+1}}\tilde\varepsilon_{ij}(\theta).$$

$n$ 为幂硬化指数，$I_n$ 为数值常数。J 积分是弹塑性断裂的控制参量，但有 $J$ 控制区尺寸限制。

**14.8 断裂韧性测试与尺寸效应**

$K_{Ic}$ 需满足平面应变与小范围屈服：

$$a,\ B,\ (W-a)\ge2.5\left(\frac{K_{Ic}}{\sigma_y}\right)^2.$$

$\because$ 塑性区尺寸 $r_p\approx\dfrac{1}{2\pi}\left(\dfrac{K_I}{\sigma_y}\right)^2$，试样尺寸须远大于 $r_p$ 方能测得平面应变值。

**14.9 疲劳裂纹扩展：Paris 定律**

$$\therefore \frac{da}{dN}=C\left(\Delta K\right)^m,$$

$C,m$ 为材料常数，$\Delta K=K_{\max}-K_{\min}$。存在门槛值 $\Delta K_{th}$ 与快速失稳区。积分得寿命

$$N=\int_{a_0}^{a_c}\frac{da}{C\left(\Delta K(a)\right)^m}.$$

对中心裂纹 $\Delta K=\Delta\sigma\sqrt{\pi a}$，积分给出 $N\!-\!a$ 关系。

**14.10 内聚区与断裂过程区**

Dugdale–Barenblatt 内聚区模型：裂纹尖端存在长度 $c$ 的塑性/内聚区，区内外加应力与内聚力平衡。Dugdale 解

$$c=\frac{\pi K^2}{8\sigma_y^2},\qquad \delta_t=\frac{K^2}{E'\sigma_y},$$

裂纹张开口位移与断裂能

$$\therefore \Gamma=\int_0^{\delta_t}\sigma(\delta)\,d\delta=G_c.$$

这统一了线弹性断裂与内聚断裂框架，用于界面与复合材料分层。

---

### 第15章 结构动力学与波动

**15.1 单自由度系统**

$$m\ddot x+c\dot x+kx=F(t).$$

无阻尼自由振动 $\ddot x+\omega_n^2x=0$，$\omega_n=\sqrt{k/m}$，解 $x=A\cos(\omega_nt-\phi)$。阻尼比 $\zeta=c/(2\sqrt{km})$，欠阻尼解

$$x=e^{-\zeta\omega_nt}\left(C_1\cos\omega_dt+C_2\sin\omega_dt\right),\qquad \omega_d=\omega_n\sqrt{1-\zeta^2}.$$

**15.2 频率响应与共振**

谐激励 $F=F_0e^{i\omega t}$，稳态 $x=Xe^{i\omega t}$：

$$X=\frac{F_0}{k-m\omega^2+ic\omega}=\frac{F_0/k}{1-\beta^2+i2\zeta\beta},\qquad \beta=\frac{\omega}{\omega_n},$$

幅值

$$\therefore |X|=\frac{F_0/k}{\sqrt{\left(1-\beta^2\right)^2+\left(2\zeta\beta\right)^2}}.$$

共振峰 $\beta\approx\sqrt{1-2\zeta^2}$，峰值放大 $\approx1/(2\zeta)$。

**15.3 多自由度与模态分析**

$$[M]\{\ddot x\}+[K]\{x\}=\{F\}.$$

设 $\{x\}=\{\phi\}e^{i\omega t}$，得广义特征值问题

$$\therefore \left([K]-\omega^2[M]\right)\{\phi\}=0,\qquad \det\!\left([K]-\omega^2[M]\right)=0.$$

模态正交性：

$$\{\phi^{(r)}\}^T[M]\{\phi^{(s)}\}=\delta_{rs}m_r,\qquad \{\phi^{(r)}\}^T[K]\{\phi^{(s)}\}=\delta_{rs}k_r,$$

$\because$ 对称矩阵广义特征向量按 $M$ 正交。用模态坐标 $x=\Phi q$ 解耦：

$$m_r\ddot q_r+k_rq_r=\phi^{(r)T}F.$$

**15.4 连续系统：杆的纵向振动**

一维杆，微元平衡

$$\because \frac{\partial N}{\partial x}=\rho A\frac{\partial^2u}{\partial t^2},\qquad N=EA\frac{\partial u}{\partial x},$$

$$\therefore \frac{\partial^2u}{\partial t^2}=c^2\frac{\partial^2u}{\partial x^2},\qquad c=\sqrt{\frac E\rho}.$$

分离变量 $u=X(x)T(t)$ 得

$$X''+k^2X=0,\qquad T''+\omega^2T=0,\qquad \omega=ck.$$

固定-自由杆 $k_n=(2n-1)\pi/(2L)$，$\omega_n=\dfrac{(2n-1)\pi}{2L}\sqrt{E/\rho}$。

**15.5 梁的横向振动（Euler–Bernoulli）**

$$\because \frac{\partial V}{\partial x}=\rho A\frac{\partial^2w}{\partial t^2},\quad M=-EI\frac{\partial^2w}{\partial x^2},\quad V=\frac{\partial M}{\partial x},$$

$$\therefore EI\frac{\partial^4w}{\partial x^4}+\rho A\frac{\partial^2w}{\partial t^2}=0.$$

设 $w=W(x)e^{i\omega t}$：

$$W''''-\beta^4W=0,\qquad \beta^4=\frac{\rho A\omega^2}{EI}.$$

简支梁边界 $W=W''=0$ 给出 $W_n=\sin\dfrac{n\pi x}{L}$，$\omega_n=\left(\dfrac{n\pi}{L}\right)^2\sqrt{\dfrac{EI}{\rho A}}$。

**15.6 波动方程与色散**

一维波动方程 $\partial_t^2u=c^2\partial_x^2u$，d'Alembert 解

$$u(x,t)=f(x-ct)+g(x+ct),$$

$$\because \partial_t^2u-c^2\partial_x^2u=\left(\partial_t-c\partial_x\right)\left(\partial_t+c\partial_x\right)u=0,$$

通解即左右行波叠加。Euler 梁弯曲波满足

$$\omega=\sqrt{\frac{EI}{\rho A}}\,\beta^2,$$

相位速度 $c_p=\omega/\beta=\sqrt{EI/(\rho A)}\,\beta$ 随波数变化，故弯曲波是色散波；杆纵波无色散。

**15.7 弹性体中的体波与表面波**

由弹性动力学方程 $\rho\ddot u_i=\partial_j\sigma_{ij}$ 代入本构，设平面波 $\mathbf{u}=\mathbf{A}e^{i(\mathbf{k}\cdot\mathbf{x}-\omega t)}$，得 Christoffel 方程

$$\left(C_{ijkl}k_jk_l-\rho\omega^2\delta_{ik}\right)A_k=0.$$

各向同性时

$$\therefore c_L=\sqrt{\frac{\lambda+2\mu}{\rho}},\qquad c_T=\sqrt{\frac{\mu}{\rho}}.$$

Rayleigh 表面波速度满足

$$\left(2-\frac{c_R^2}{c_T^2}\right)^2=4\sqrt{1-\frac{c_R^2}{c_L^2}}\sqrt{1-\frac{c_R^2}{c_T^2}},\qquad c_R\approx0.92c_T.$$

**15.8 模态叠加与瞬态响应**

$$x(t)=\sum_r\phi^{(r)}q_r(t),\qquad q_r(t)=\frac{1}{m_r\omega_{dr}}\int_0^tF_r(\tau)e^{-\zeta_r\omega_r(t-\tau)}\sin\omega_{dr}(t-\tau)\,d\tau.$$

$\because$ 每个模态是单自由度杜哈梅积分。这是结构地震响应谱法和冲响应分析的依据。

**15.9 能量法与 Rayleigh 商求固有频率**

$$R(w)=\frac{\displaystyle\int_0^L EI(w'')^2dx}{\displaystyle\int_0^L\rho A w^2dx}.$$

对真实模态 $R=\omega^2$，对任意试函数 $R$ 给出 $\omega_1^2$ 的上界。$\because$ 由 Rayleigh 原理：

$$\therefore \omega_1^2=\min_w R(w).$$

Rayleigh–Ritz 用多个试函数即得前几阶频率。

**15.10 典型应用：旋转梁与离心刚化**

高速旋转叶片受离心力产生几何刚度，横向振动方程

$$\rho A\ddot w-\frac{\partial}{\partial x}\!\left(T(x)\frac{\partial w}{\partial x}\right)+EI\frac{\partial^4w}{\partial x^4}=0,$$

其中 $T(x)=\int_x^L\rho A\Omega^2r\,dr$。$\because$ 离心张力提供正几何刚度，故随转速升高固有频率增大，可能穿越共振区。这是航空发动机叶片动力学核心问题。

---

### 第16章 非线性连续介质与有限变形

**16.1 有限变形度量**

变形梯度 $F_{iJ}=\partial x_i/\partial X_J$，右 Cauchy–Green 张量

$$C=F^TF,\qquad C_{IJ}=F_{iI}F_{iJ},$$

左 Cauchy–Green 张量

$$b=FF^T,\qquad b_{ij}=F_{iJ}F_{jJ}.$$

Green–Lagrange 应变

$$E=\tfrac12\left(C-I\right),$$

Almansi 应变

$$e=\tfrac12\left(I-b^{-1}\right).$$

$\because$ $ds^2=d\mathbf{x}\!\cdot\!d\mathbf{x}=d\mathbf{X}^T C\,d\mathbf{X}$，$dS^2=d\mathbf{X}\!\cdot\!d\mathbf{X}$，故 $E$ 直接度量长度平方变化。

**16.2 极分解与转动客观性**

$$F=RU=VR,$$

$R$ 正交（转动），$U,V$ 对称正定（右/左伸长）。$\because$ $C=F^TF=U^2$，$b=FF^T=V^2$，$R=FU^{-1}=V^{-1}F$。

在叠加刚体运动 $\mathbf{x}^*=Q\mathbf{x}+\mathbf{c}$ 下，$F^*=QF$，$C^*=C$，故 $E$ 不变，体现物质客观性。

**16.3 应力度量：Cauchy、PK1、PK2**

Cauchy 应力 $\boldsymbol{\sigma}$（现时构形）。第一类 Piola–Kirchhoff 应力 $P$：

$$P_{iJ}=J\sigma_{ij}F_{Jj}^{-1},\qquad J=\det F,$$

$$\because d\mathbf{f}=\boldsymbol{\sigma}\,d\mathbf{a}=\boldsymbol{\sigma}J F^{-T}d\mathbf{A}=P\,d\mathbf{A}.$$

第二类 Piola–Kirchhoff 应力 $S$：

$$S=F^{-1}P=JF^{-1}\boldsymbol{\sigma}F^{-T},\qquad \boldsymbol{\sigma}=\frac1J F S F^T.$$

$S$ 在参考构形对称，与 $E$ 能量共轭。

**16.4 参考构形下的守恒律**

质量守恒：

$$\rho_0=\rho J.$$

动量守恒：

$$\frac{\partial P_{iJ}}{\partial X_J}+f_{0i}=\rho_0\ddot x_i.$$

$\because$ 用 Nanson 公式 $da=JF^{-T}dA$ 与 $\rho J=\rho_0$ 变换即得。能量共轭对

$$\dot W=S_{IJ}\dot E_{IJ},$$

是超弹性本构的基础。

**16.5 超弹性（Green 弹性）**

存在应变能函数 $\Psi(E)$，使

$$S_{IJ}=\frac{\partial\Psi}{\partial E_{IJ}},\qquad C_{IJKL}=\frac{\partial^2\Psi}{\partial E_{IJ}\partial E_{KL}},$$

$C_{IJKL}$ 具有主对称性。$\because$ 由热力学第二定律 $\dot W=S:\dot E$，且 $\Psi$ 为势函数。各向同性超弹性 $\Psi=\Psi(I_1,I_2,I_3)$，$I_k$ 为 $C$ 的不变量。

**16.6 Neo-Hookean 与 Mooney–Rivlin 模型**

不可压 Neo-Hookean：

$$\Psi=\frac{\mu}{2}\left(I_1-3\right),$$

$$\because \boldsymbol{\sigma}=-pI+\mu b,\qquad J=1.$$

Mooney–Rivlin：

$$\Psi=C_{10}\left(I_1-3\right)+C_{01}\left(I_2-3\right),$$

$$\therefore \boldsymbol{\sigma}=-pI+2\!\left(C_{10}+C_{01}I_1\right)b-2C_{01}b^2.$$

用于橡胶与生物软组织大变形。

**16.7 线性化与弹性张量**

小变形时 $F\approx I+\nabla\mathbf{u}$，$E\approx\boldsymbol{\varepsilon}$，$S\approx\boldsymbol{\sigma}$。由 $C_{ijkl}=\left.\dfrac{\partial^2\Psi}{\partial E_{ij}\partial E_{kl}}\right|_{E=0}$，各向同性超弹性线性化得

$$\therefore C_{ijkl}=\lambda\delta_{ij}\delta_{kl}+\mu\left(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}\right),$$

回到第5章线弹性。物理意义：线弹性是有限变形理论的几何与材料双重线性化。

**16.8 客观率与应力率**

有限变形塑性需用客观应力率。Jaumann 率

$$\overset{\triangle}{\boldsymbol{\sigma}}=\dot{\boldsymbol{\sigma}}-\mathbf{W}\boldsymbol{\sigma}+\boldsymbol{\sigma}\mathbf{W},\qquad \mathbf{W}=\tfrac12\!\left(\mathbf{L}-\mathbf{L}^T\right),$$

$$\because \overset{\triangle}{\boldsymbol{\sigma}}=Q\!\left[\overset{\triangle}{\boldsymbol{\sigma}}\right]Q^T\ \text{在刚体转动下客观},$$

其中 $\mathbf{L}=\dot F F^{-1}$ 为速度梯度。Truesdell 率、Green–Naghdi 率亦常用。

**16.9 有限变形下的稳定性与分叉**

增量平衡方程

$$\int_V\left(\dot S:\delta\dot E+\boldsymbol{\sigma}:\delta\!\left(\delta\mathbf{u}\cdot\nabla\right)\!\cdot\nabla\right)dV=\text{外载率功}.$$

更规范地，以参考构形上的虚功率原理写出：

$$\because\ \int_V\left(S_{IJ}\,\delta\dot E_{IJ}+\dot S_{IJ}\,\delta E_{IJ}\right)dV_0=\int_{\partial V}t_i\,\delta v_i\,dA+\int_V f_i\,\delta v_i\,dV,$$

$$\therefore\ \int_V\left(\dot S_{IJ}\,\delta E_{IJ}+S_{IJ}\,\delta\dot E_{IJ}\right)dV_0=\delta\dot W_{\text{ext}}.$$

当切线刚度算子的二阶变分失去正定性，结构发生分叉或极限失稳。$\because$ 令 $\delta^2\Pi=0$ 得临界条件。几何非线性（大位移）与材料非线性（塑性、损伤）耦合，是航空航天柔性结构与薄壁结构分析的核心。

**16.10 典型应用：不可压橡胶块的单轴压缩**

不可压 Neo-Hookean 单轴压缩，沿 $x_1$ 方向伸长比 $\lambda$，不可压 $\lambda_1\lambda_2\lambda_3=1$，取横向等比 $\lambda_2=\lambda_3=\lambda^{-1/2}$。Cauchy 应力

$$\sigma_{11}=\mu\left(\lambda^2-\lambda^{-1}\right),\qquad \sigma_{22}=\sigma_{33}=0,$$

$$\because \boldsymbol{\sigma}=-pI+\mu b,\ b_{11}=\lambda^2,\ b_{22}=b_{33}=\lambda^{-1},\ p=\mu\lambda^{-1}.$$

$\therefore \sigma=\mu\left(\lambda^2-\lambda^{-1}\right)$。小 $\lambda-1$ 展开发散得 $\sigma\approx3\mu(\lambda-1)$，即小应变杨氏模量 $E=3\mu$，与线弹性不可压关系一致。

---
⬆ 上级领域：[[力学]]
