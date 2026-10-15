---
课程名称: 固体力学核心贯通课I
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

# 固体力学核心贯通课I


### 第1章 张量代数、不变量与场论基础

#### 知识点1. 指标记号、求和约定与Kronecker符号

**定义**：在三维欧氏空间中，矢量 $\mathbf a$ 的分量记为 $a_i\ (i=1,2,3)$。为使书写简洁，规定：当一单项式中某一哑指标（dummy index）重复出现且恰好两次时，表示对该指标在 $1,2,3$ 上求和，称 Einstein 求和约定，即
$$a_i b_i \equiv \sum_{i=1}^{3} a_i b_i.$$
Kronecker 符号定义为
$$\delta_{ij}=\begin{cases}1,& i=j,\\0,& i\ne j.\end{cases}$$

**假设**：讨论限于三维 Euclidean 空间、正交直线坐标系的可微场。

**公式与推导**：
(1) 分量代换性质。$\because \delta_{ij}$ 只在 $i=j$ 时为 $1$，对任意分量 $a_j$ 有
$$\delta_{ij}a_j=\sum_{j=1}^3\delta_{ij}a_j=a_i.$$
$\therefore$ $\delta_{ij}$ 起"换标"作用。

(2) $\delta_{ii}=3$。$\because$
$$\delta_{ii}=\sum_{i=1}^3\delta_{ii}=1+1+1=3.$$

(3) 正交性。设单位基 $\mathbf e_i$，则
$$\mathbf e_i\cdot\mathbf e_j=\delta_{ij}.$$

**几何意义**：$\delta_{ij}$ 是恒等算子的分量，代表正交归一基的内积结构。

**极限/特例**：当维数为 $n$ 时 $\delta_{ii}=n$；二维 $\delta_{ii}=2$。

**典型应用**：用 $\delta_{ij}$ 与 $\epsilon_{ijk}$ 书写各向同性张量，如各向同性四阶张量 $C_{ijkl}=\lambda\delta_{ij}\delta_{kl}+\mu(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk})$。

#### 知识点2. 置换符号与矢量叉积恒等式

**定义**：置换符号 $\epsilon_{ijk}$ 定义为
$$\epsilon_{ijk}=\begin{cases}+1,&(i,j,k)\text{ 为 }(1,2,3)\text{ 的偶排列},\\-1,&(i,j,k)\text{ 为奇排列},\\0,&\text{有重复指标}.\end{cases}$$

**公式与推导**：
(1) 叉积分量式。设 $\mathbf c=\mathbf a\times\mathbf b$，则 $\mathbf e_i\times\mathbf e_j=\epsilon_{ijk}\mathbf e_k$，于是
$$\mathbf c=\mathbf a\times\mathbf b=a_i\mathbf e_i\times b_j\mathbf e_j=a_ib_j\,\epsilon_{ijk}\mathbf e_k,$$
$\therefore c_k=\epsilon_{kij}a_ib_j.$

(2) 核心恒等式（$\epsilon$-$\delta$ 恒等式）：
$$\epsilon_{ijk}\epsilon_{ilm}=\delta_{jl}\delta_{km}-\delta_{jm}\delta_{kl}.$$

**完整证明**：利用行列式表示
$$\epsilon_{ijk}=\begin{vmatrix}\delta_{i1}&\delta_{i2}&\delta_{i3}\\\delta_{j1}&\delta_{j2}&\delta_{j3}\\\delta_{k1}&\delta_{k2}&\delta_{k3}\end{vmatrix}.$$
于是
$$\epsilon_{ijk}\epsilon_{ilm}=\begin{vmatrix}\delta_{i1}&\delta_{i2}&\delta_{i3}\\\delta_{j1}&\delta_{j2}&\delta_{j3}\\\delta_{k1}&\delta_{k2}&\delta_{k3}\end{vmatrix}\begin{vmatrix}\delta_{i1}&\delta_{i2}&\delta_{i3}\\\delta_{l1}&\delta_{l2}&\delta_{l3}\\\delta_{m1}&\delta_{m2}&\delta_{m3}\end{vmatrix}.$$
对行指标 $i$ 求和等价于"对三行逐行做内积"。由行列式乘法与对 $i$ 的求和，得
$$\epsilon_{ijk}\epsilon_{ilm}=\begin{vmatrix}\delta_{jl}&\delta_{jm}\\\delta_{kl}&\delta_{km}\end{vmatrix}=\delta_{jl}\delta_{km}-\delta_{jm}\delta_{kl}.$$
更严格地，利用 $\epsilon_{ijk}\epsilon_{ilm}=\det(\delta_{i\cdot}\text{行组合})$，展开 Cauchy–Binet：
$$\epsilon_{ijk}\epsilon_{ilm}=\delta_{jl}\delta_{km}-\delta_{jm}\delta_{kl}.$$
另一常用形式（对第二个取和）：
$$\epsilon_{ijk}\epsilon_{ijl}=2\delta_{kl},\qquad \epsilon_{ijk}\epsilon_{ijk}=6.$$
**推导**：令上式 $m\to j$（即第二个 $\epsilon$ 的第三个指标与第一个的第二个指标相同）：
$$\epsilon_{ijk}\epsilon_{ilj}=\delta_{jl}\delta_{kj}-\delta_{jj}\delta_{kl}=\delta_{kl}-3\delta_{kl}=-2\delta_{kl}.$$
注意 $\epsilon_{ilj}=-\epsilon_{ijl}$，故 $\epsilon_{ijk}\epsilon_{ijl}=2\delta_{kl}$。再令 $l\to k$：$\epsilon_{ijk}\epsilon_{ijk}=2\delta_{kk}=6$。

(3) 标量三重积 $\mathbf a\cdot(\mathbf b\times\mathbf c)=\epsilon_{ijk}a_ib_jc_k$，即 $[\mathbf a,\mathbf b,\mathbf c]$，几何上为平行六面体体积。

(4) 矢量三重积：
$$\mathbf a\times(\mathbf b\times\mathbf c)=\mathbf b(\mathbf a\cdot\mathbf c)-\mathbf c(\mathbf a\cdot\mathbf b).$$
**证明**：$(\mathbf a\times(\mathbf b\times\mathbf c))_i=\epsilon_{ijk}a_j(\mathbf b\times\mathbf c)_k=\epsilon_{ijk}a_j\epsilon_{klm}b_lc_m=\epsilon_{ijk}\epsilon_{klm}a_jb_lc_m$。用 $\epsilon_{ijk}=\epsilon_{kij}$ 及 $\epsilon_{kij}\epsilon_{klm}=\delta_{il}\delta_{jm}-\delta_{im}\delta_{jl}$，得
$$=\delta_{il}\delta_{jm}a_jb_lc_m-\delta_{im}\delta_{jl}a_jb_lc_m=a_jb_ic_j-a_jb_jc_i,$$
即分量为 $b_i(\mathbf a\cdot\mathbf c)-c_i(\mathbf a\cdot\mathbf b)$。

**几何意义**：$\epsilon$ 编码了右手系的定向与反对称性。

**极限/特例**：共面时 $\epsilon_{ijk}a_ib_jc_k=0$。

**典型应用**：涡量 $\boldsymbol\omega=\nabla\times\mathbf v$、角速度、刚体转动运动学。

#### 知识点3. 二阶张量、运算与坐标变换

**定义**：二阶张量 $\mathbf T$ 是将矢量线性映射为矢量的线性算子，$\mathbf T=T_{ij}\mathbf e_i\otimes\mathbf e_j$，作用为 $\mathbf T\mathbf a=T_{ij}a_j\mathbf e_i$。分量满足坐标变换
$$T'_{ij}=Q_{ik}Q_{jl}T_{kl},\qquad Q_{ik}Q_{jl}\delta_{kl}=\delta_{ij},$$
其中 $Q_{ij}=\mathbf e'_i\cdot\mathbf e_j$ 为正交矩阵。

**公式与推导**：
(1) 正交矩阵性质。$\because \mathbf e'_i\cdot\mathbf e'_j=\delta_{ij}$，而 $\mathbf e'_i=Q_{ik}\mathbf e_k$，$\therefore Q_{ik}Q_{jl}\delta_{kl}=\delta_{ij}$，即 $Q Q^{\mathrm T}=I$，$Q^{\mathrm T}=Q^{-1}$，$\det Q=\pm1$（正常转动取 $+1$）。

(2) 张量不变量。迹、行列式、Frobenius 范数在正交变换下不变。以迹为例：
$$T'_{ii}=Q_{ik}Q_{il}T_{kl}=\delta_{kl}T_{kl}=T_{kk}.$$

(3) 张量积与缩并。$(\mathbf a\otimes\mathbf b)_{ij}=a_ib_j$；双点积 $\mathbf A:\mathbf B=A_{ij}B_{ij}$。

**几何意义**：二阶张量可分解为"拉伸+旋转"；对称部分对应纯变形，反对称部分对应刚体转动。

**极限/特例**：$\mathbf T$ 对称时 $T_{ij}=T_{ji}$，其独立分量由 $9$ 个减为 $6$ 个；反对称时 $T_{ij}=-T_{ji}$，可表示为轴矢量 $W_{ij}=-\epsilon_{ijk}\omega_k$。

**典型应用**：应变张量、应力张量、惯性张量、变形梯度。

#### 知识点4. 商法则与张量识别

**定义（商法则）**：若对任意矢量 $a_j$，量 $p_i=T_{ij}a_j$ 均为矢量，则 $T_{ij}$ 为二阶张量。

**证明**：设 $p$ 为矢量，$\therefore p'_i=Q_{ik}p_k$。又 $p'_i=T'_{ij}a'_j=T'_{ij}Q_{jl}a_l$。而 $p_i=T_{ij}a_j\Rightarrow Q_{ik}p_k=Q_{ik}T_{kj}a_j$。比较任意 $a_l$ 的系数：
$$T'_{ij}Q_{jl}=Q_{ik}T_{kl}\ \Rightarrow\ T'_{ij}=Q_{ik}T_{kl}Q_{jl}=Q_{ik}Q_{jl}T_{kl},$$
满足二阶张量变换律。$\blacksquare$

**几何意义**：商法则给出"由物理关系识别张量阶数"的判据，无需直接验证变换律。

**极限/特例**：若 $p_i=T_{ijkl}a_{jkl}$ 且左边为矢量，则四阶。

**典型应用**：由 $\mathbf t=\boldsymbol\sigma\mathbf n$ 断定 $\boldsymbol\sigma$ 为二阶张量（Cauchy 定理）。

#### 知识点5. 特征值、主不变量与Cayley–Hamilton定理

**定义**：若 $\mathbf T\mathbf n=\lambda\mathbf n$（$\mathbf n\ne\mathbf 0$），称 $\lambda$ 为特征值（主值），$\mathbf n$ 为特征矢量（主方向）。特征方程为
$$\det(T_{ij}-\lambda\delta_{ij})=0.$$

**公式与推导**：
(1) 特征方程展开。三维下
$$\det(\mathbf T-\lambda\mathbf I)=-\lambda^3+I_1\lambda^2-I_2\lambda+I_3=0,$$
其中主不变量
$$I_1=\operatorname{tr}\mathbf T=T_{ii},$$
$$I_2=\tfrac12\big[(T_{ii})^2-T_{ij}T_{ji}\big],\qquad I_3=\det\mathbf T.$$
**推导**：直接展开 $\det(T_{ij}-\lambda\delta_{ij})$。对 $3\times3$，常数项 $\det T$，一次项来自去掉一行一列的主子式之和。可证明 $I_2=\frac12(T_{ii}T_{jj}-T_{ij}T_{ji})$，$I_1=T_{ii}$。

(2) Cayley–Hamilton 定理：
$$\mathbf T^3-I_1\mathbf T^2+I_2\mathbf T-I_3\mathbf I=\mathbf 0.$$
即张量满足自身特征方程（对每个特征方向显然，按连续性/密实性推广到全空间）。

**几何意义**：$I_1,I_2,I_3$ 是张量的坐标不变量，刻画体积、形状等内禀属性。

**极限/特例**：对称张量特征值全为实、特征方向正交。对偏张量 $s_{ij}=\sigma_{ij}-p\delta_{ij}$，$I_1=0$。

**典型应用**：主应力 $J_1=\sigma_{ii}$，Mises 等效应力依赖 $J_2=\frac12 s_{ij}s_{ij}$。

#### 知识点6. 张量分解：对称/反对称、球张量/偏张量

**定义与公式**：
(1) 任意二阶张量分解为
$$T_{ij}=\underbrace{\tfrac12(T_{ij}+T_{ji})}_{\text{对称 }T_{(ij)}}+\underbrace{\tfrac12(T_{ij}-T_{ji})}_{\text{反对称 }T_{[ij]}}.$$
**证明**：令 $S_{ij}=\frac12(T_{ij}+T_{ji})$，$A_{ij}=\frac12(T_{ij}-T_{ji})$，则 $S_{ij}=S_{ji}$，$A_{ij}=-A_{ji}$，且 $S+A=T$。$\blacksquare$

(2) 应力张量的球–偏分解：
$$\sigma_{ij}=s_{ij}+p\delta_{ij},\qquad p=\tfrac13\sigma_{kk},$$
其中 $p$ 为平均应力（球张量），$s_{ij}$ 为应力偏量，$s_{ii}=0$。
**推导**：$\because s_{ii}=\sigma_{ii}-\frac13\sigma_{kk}\delta_{ii}=\sigma_{ii}-\sigma_{kk}=0$。

**几何意义**：球张量改变体积（静水压），偏张量改变形状（剪切/畸变）。

**极限/特例**：不可压缩时体积应变为零，变形全由偏张量承担。

**典型应用**：J2 塑性理论仅用 $s_{ij}$；应力波分解为膨胀波与畸变波。

#### 知识点7. 场论：梯度、散度、旋度与积分定理

**定义**：标量场 $\phi$、矢量场 $\mathbf u$ 的微分算子
$$(\nabla\phi)_i=\partial_i\phi,\qquad \nabla\cdot\mathbf u=\partial_iu_i,\qquad (\nabla\times\mathbf u)_i=\epsilon_{ijk}\partial_ju_k.$$

**公式与推导**：
(1) 散度定理（Gauss）：
$$\int_V \partial_i u_i\,dV=\oint_{\partial V}u_in_i\,dS.$$
**证明思路**：以一维微元 $\partial_i u_i\,dV$ 沿 $x_i$ 方向累积为两端面通量之差，三维推广逐分量求和。等价 $\int_V\nabla\cdot\mathbf u\,dV=\oint\mathbf u\cdot\mathbf n\,dS$。

(2) 分部积分（Green 恒等式）：
$$\int_V (\partial_i u_i)v\,dV=\oint_{\partial V}u_in_i v\,dS-\int_V u_i\partial_i v\,dV.$$
由乘积法则 $\partial_i(u_iv)=(\partial_iu_i)v+u_i\partial_iv$ 加散度定理即得。

(3) 常用恒等式，例如
$$\nabla\cdot(\nabla\times\mathbf u)=0,\qquad \nabla\times(\nabla\phi)=\mathbf 0.$$
**证明**：$\partial_i\epsilon_{ijk}\partial_ju_k=\epsilon_{ijk}\partial_i\partial_ju_k$，而 $\partial_i\partial_j$ 对称、$\epsilon$ 反对称，缩并为 $0$。同理 $\epsilon_{ijk}\partial_j\partial_k\phi=0$。

(4) 向量 Laplacian 恒等式：
$$\nabla\times(\nabla\times\mathbf u)=\nabla(\nabla\cdot\mathbf u)-\nabla^2\mathbf u.$$
**证明**：$(\nabla\times(\nabla\times\mathbf u))_i=\epsilon_{ijk}\partial_j\epsilon_{klm}\partial_lu_m=\epsilon_{ijk}\epsilon_{klm}\partial_j\partial_lu_m=(\delta_{il}\delta_{jm}-\delta_{im}\delta_{jl})\partial_j\partial_lu_m=\partial_i\partial_ju_j-\partial_j\partial_ju_i=\partial_i(\nabla\cdot\mathbf u)-\nabla^2u_i.$

**几何意义**：散度为源密度，旋度为环量密度，梯度为最大变化方向。

**极限/特例**：无旋场可写为势梯度；无源场可写为旋度。

**典型应用**：连续介质守恒律、Navier–Stokes 的矢量形式、弹性动力学。

#### 知识点8. 正交曲线坐标：柱坐标与球坐标

**定义**：设正交曲线坐标 $(q_1,q_2,q_3)$，Lamé 系数
$$h_i=\left|\frac{\partial\mathbf r}{\partial q_i}\right|,\qquad ds^2=h_1^2dq_1^2+h_2^2dq_2^2+h_3^2dq_3^2.$$

**公式与推导**：
(1) 柱坐标 $(r,\theta,z)$：$h_r=1,h_\theta=r,h_z=1$，梯度、散度、Laplacian
$$\nabla\phi=\left(\partial_r\phi,\tfrac1r\partial_\theta\phi,\partial_z\phi\right),$$
$$\nabla\cdot\mathbf u=\frac1r\partial_r(ru_r)+\frac1r\partial_\theta u_\theta+\partial_zu_z,$$
$$\nabla^2\phi=\frac1r\partial_r(r\partial_r\phi)+\frac1{r^2}\partial_\theta^2\phi+\partial_z^2\phi.$$
**推导**：微元弧长 $ds^2=dr^2+r^2d\theta^2+dz^2$。散度由面通量定义，径向两面净通量
$$\frac{(ru_r)|_{r+dr}-(ru_r)|_r\,d\theta dz}{r\,dr\,d\theta dz}=\frac1r\partial_r(ru_r),$$
其余方向类似相加得证。

(2) 球坐标 $(r,\theta,\varphi)$：$h_r=1,h_\theta=r,h_\varphi=r\sin\theta$，Laplacian
$$\nabla^2\phi=\frac1{r^2}\partial_r(r^2\partial_r\phi)+\frac1{r^2\sin\theta}\partial_\theta(\sin\theta\,\partial_\theta\phi)+\frac1{r^2\sin^2\theta}\partial_\varphi^2\phi.$$
**推导**：由 $ds^2=dr^2+r^2d\theta^2+r^2\sin^2\theta\,d\varphi^2$，按 $\nabla^2\phi=\frac1{h_1h_2h_3}\partial_i\!\left(\frac{h_1h_2h_3}{h_i^2}\partial_i\phi\right)$ 代入 $h_i$ 即得。

**几何意义**：Lamé 系数刻画坐标网格沿各方向的伸缩；拉普拉斯算子是坐标不变量。

**极限/特例**：球对称时 $\nabla^2\phi=\frac1{r^2}\partial_r(r^2\partial_r\phi)$。

**典型应用**：厚壁圆筒（柱坐标）、球孔扩张（球坐标）、Kirsch 解、Boussinesq 问题。

### 第2章 变形几何、运动学与应变协调

#### 知识点1. 构型、物质坐标与变形梯度

**定义**：连续体在参考构型 $\mathcal B_0$ 中质点用物质坐标 $\mathbf X$ 标记，在现时构型 $\mathcal B$ 中位置为 $\mathbf x=\boldsymbol\varphi(\mathbf X,t)$。变形梯度
$$F_{iJ}=\frac{\partial x_i}{\partial X_J},\qquad \mathbf F=\frac{\partial\mathbf x}{\partial\mathbf X}.$$

**假设**：连续、单值、足够光滑；$\det\mathbf F=J>0$（保持定向、无物质穿透）。

**公式与推导**：
(1) 功共轭线元变换：$dx_i=F_{iJ}\,dX_J$。
(2) 体积变换：面元 $d\mathbf a$ 与体元
$$dV=J\,dV_0,\qquad J=\det\mathbf F.$$
**推导**：$\because d\mathbf x_1,d\mathbf x_2,d\mathbf x_3$ 张成的体积 $dV=[d\mathbf x_1,d\mathbf x_2,d\mathbf x_3]$，而 $d\mathbf x_i=\mathbf F d\mathbf X_i$，由行列式多重线性 $\therefore dV=(\det\mathbf F)[d\mathbf X_1,d\mathbf X_2,d\mathbf X_3]=J\,dV_0$。
(3) Nanson 公式（面元）：
$$n_i\,da=J\,F^{-T}_{iJ}N_J\,dA.$$
**推导**：由 $dV=\mathbf n\cdot d\mathbf x\,da$ 与 $dV=J\,dV_0$，结合线元映射并用 $F^{-T}$ 变换面积法向，即得。

**物理/几何意义**：$J$ 为局部体积比；$J=1$ 表示等容（不可压缩）。

**极限/特例**：小变形 $F_{iJ}=\delta_{iJ}+u_{i,J}$，$J\approx 1+\partial_iu_i$。

**典型应用**：有限变形理论、超弹性本构、塑性不可压缩约束。

#### 知识点2. 应变张量：Green–Lagrange与Almansi

**定义**：右 Cauchy–Green 张量 $C_{IJ}=F_{kI}F_{kJ}$，Green–Lagrange 应变
$$E_{IJ}=\tfrac12(C_{IJ}-\delta_{IJ})=\tfrac12(F_{kI}F_{kJ}-\delta_{IJ}).$$
左 Cauchy–Green 张量 $b_{ij}=F_{iK}F_{jK}$，Almansi 应变 $e_{ij}=\frac12(\delta_{ij}-b^{-1}_{ij})$。

**公式与推导**：线元长度平方
$$ds^2-dS^2=(F_{iJ}dX_J)(F_{iK}dX_K)-dX_IdX_I=(C_{IJ}-\delta_{IJ})dX_IdX_J=2E_{IJ}dX_IdX_J.$$
$\therefore E$ 度量长度平方变化，且以物质坐标表达，故客观（对刚体转动不变，因 $C$ 不随空间转动改变）。

**物理/几何意义**：$E_{IJ}$ 是"物质型"应变，无刚性转动污染。

**极限/特例**：小变形下 $E_{IJ}\to\varepsilon_{ij}=\frac12(u_{i,j}+u_{j,i})$。

**典型应用**：橡胶等大变形、几何非线性有限元。

#### 知识点3. 小应变张量与几何意义

**定义**：小位移梯度 $u_{i,j}$，小应变张量
$$\varepsilon_{ij}=\tfrac12(u_{i,j}+u_{j,i}),\qquad \omega_{ij}=\tfrac12(u_{i,j}-u_{j,i}).$$

**公式与推导**：
(1) 正应变：考虑沿 $x$ 的微元 $dX$，伸长 $\Delta=u_x(x+dx)-u_x(x)=\partial_xu_x\,dx$，
$$\varepsilon_{xx}=\lim_{dx\to0}\frac{\Delta}{dx}=\partial_xu_x.$$
(2) 剪应变：两正交微元夹角变化
$$\gamma_{xy}=2\varepsilon_{xy}=\partial_xu_y+\partial_yu_x.$$
**推导**：$x$ 向微元转角度 $\alpha\approx u_{y,x}$，$y$ 向微元转角度 $\beta\approx u_{x,y}$，夹角改变为 $\alpha+\beta$。
(3) 转动：$\omega_x=\frac12(\partial_yu_z-\partial_zu_y)$ 等。

**物理/几何意义**：对称部分纯变形，反对称部分刚体转动。

**极限/特例**：$\varepsilon_{xy}=\frac12\gamma_{xy}$；工程剪应变 $\gamma_{xy}$ 是张量剪应变的 2 倍。

**典型应用**：应变片测量、线性弹性理论。

#### 知识点4. 应变协调方程（Saint–Venant）

**定义**：给定对称应变场 $\varepsilon_{ij}$，存在单值位移场 $u_i$ 使其满足 $\varepsilon_{ij}=\frac12(u_{i,j}+u_{j,i})$ 的充要条件，称协调方程。

**公式与推导**：
对 $\varepsilon_{ij}$ 二阶求导并交换指标：
$$\partial_k\partial_l\varepsilon_{ij}+\partial_i\partial_j\varepsilon_{kl}-\partial_j\partial_l\varepsilon_{ik}-\partial_i\partial_k\varepsilon_{jl}=0.$$
**完整推导**：由 $\varepsilon_{ij}=\frac12(u_{i,j}+u_{j,i})$。求组合
$$\varepsilon_{ij,kl}+\varepsilon_{kl,ij}-\varepsilon_{ik,jl}-\varepsilon_{jl,ik}:$$
$$\varepsilon_{ij,kl}=\tfrac12(u_{i,jkl}+u_{j,ikl}),\quad \varepsilon_{kl,ij}=\tfrac12(u_{k,lij}+u_{l,kij}),$$
$$\varepsilon_{ik,jl}=\tfrac12(u_{i,kjl}+u_{k,ijl}),\quad \varepsilon_{jl,ik}=\tfrac12(u_{j,lik}+u_{l,jik}).$$
由于 $u$ 光滑，偏导可交换，上式右边配对抵消为零。$\therefore$ 协调条件成立。三维中独立方程由 $81$ 个减为 $6$ 个（Beltrami 形式）。

**物理/几何意义**：应变场必须"来自"一个连续位移场，否则会出现位错/裂缝。

**极限/特例**：平面问题协调方程只有一个：
$$\frac{\partial^2\varepsilon_{xx}}{\partial y^2}+\frac{\partial^2\varepsilon_{yy}}{\partial x^2}=2\frac{\partial^2\varepsilon_{xy}}{\partial x\partial y}.$$

**典型应用**：应力函数法（须使应力对应的应变满足协调），光弹性。

#### 知识点5. 速度梯度、应变率与旋率

**定义**：速度场 $\mathbf v(\mathbf x,t)$，速度梯度 $L_{ij}=\partial v_i/\partial x_j$，变形率张量 $D_{ij}=\frac12(L_{ij}+L_{ji})$，旋率张量 $W_{ij}=\frac12(L_{ij}-L_{ji})$。

**公式与推导**：对两相邻质点
$$dv_i=L_{ij}dx_j=D_{ij}dx_j+W_{ij}dx_j.$$
$D$ 引起长度变化与剪切，$W$ 对应刚体转动，$\omega_i=-\frac12\epsilon_{ijk}W_{jk}$ 为角速度。

**物理/几何意义**：$D$ 是速率型应变（用于流体与塑性），$W$ 是涡量的一半。

**极限/特例**：刚体运动 $D=0$；不可压缩 $D_{ii}=0$。

**典型应用**：Navier–Stokes 中黏性应力正比于 $D$；大变形塑性流动。

#### 知识点6. 极分解与客观性

**定义（极分解）**：任何 $\det\mathbf F>0$ 的 $\mathbf F$ 可唯一分解为
$$\mathbf F=\mathbf R\mathbf U=\mathbf V\mathbf R,$$
其中 $\mathbf R$ 正交（$\det\mathbf R=1$），$\mathbf U,\mathbf V$ 对称正定（右/左伸长张量）。

**证明**：令 $\mathbf C=\mathbf F^{\mathrm T}\mathbf F$，对称正定，存在唯一对称正定 $\mathbf U=\sqrt{\mathbf C}$。令 $\mathbf R=\mathbf F\mathbf U^{-1}$，则
$$\mathbf R^{\mathrm T}\mathbf R=(\mathbf U^{-1})^{\mathrm T}\mathbf F^{\mathrm T}\mathbf F\mathbf U^{-1}=\mathbf U^{-1}\mathbf C\mathbf U^{-1}=\mathbf I,$$
故 $\mathbf R$ 正交，$\det\mathbf R=1$（由连续性且 $\det\mathbf F>0$）。$\mathbf V=\mathbf R\mathbf U\mathbf R^{\mathrm T}$ 为对称正定。$\blacksquare$

**物理/几何意义**：先纯拉伸 $\mathbf U$ 再转动 $\mathbf R$（或先转动后拉伸）；$\mathbf U,\mathbf V$ 与 $\mathbf C,\mathbf b$ 同特征值。

**极限/特例**：小变形 $\mathbf U\approx\mathbf I+\boldsymbol\varepsilon$，$\mathbf R\approx\mathbf I+\boldsymbol\omega$。

**典型应用**：大变形本构的客观性、有限元共旋格式。

#### 知识点7. 体积应变与不可压缩条件

**定义**：小变形体积应变
$$\varepsilon_V=\varepsilon_{ii}=\varepsilon_{xx}+\varepsilon_{yy}+\varepsilon_{zz}=\mathrm{div}\,\mathbf u.$$
有限变形体积比 $J=\det\mathbf F$。

**公式与推导**：小变形 $\det\mathbf F\approx\det(\mathbf I+\nabla\mathbf u)\approx 1+\operatorname{tr}(\nabla\mathbf u)=1+\varepsilon_{ii}$，故体积应变 $\frac{dV-dV_0}{dV_0}\approx\varepsilon_{ii}$。
不可压缩条件：
$$\varepsilon_{ii}=0\ (\text{小变形}),\qquad J=1\ (\text{有限变形}),\qquad D_{ii}=0\ (\text{速率型}).$$

**物理/几何意义**：反映体积变化，与球张量/静水压共轭。

**极限/特例**：各向同性线弹性中 $\varepsilon_{ii}=\sigma_{kk}/(3K)$，$K$ 为体积模量。

**典型应用**：橡胶/金属塑性体积不可压、热膨胀修正 $\varepsilon_{ii}=3\alpha\Delta T+\sigma_{kk}/(3K)$。

#### 知识点8. 应变张量的主值与应变莫尔圆

**定义**：$\varepsilon_{ij}$ 的特征值 $\varepsilon_1\ge\varepsilon_2\ge\varepsilon_3$ 为主应变，特征方向为主应变轴。

**公式与推导**：特征方程
$$\det(\varepsilon_{ij}-\lambda\delta_{ij})=0\Rightarrow -\lambda^3+J_1^\varepsilon\lambda^2-J_2^\varepsilon\lambda+J_3^\varepsilon=0.$$
在主轴系中 $\varepsilon_{ij}=\mathrm{diag}(\varepsilon_1,\varepsilon_2,\varepsilon_3)$。斜截面法向 $\mathbf n$，正应变 $\varepsilon=\varepsilon_in_i^2$，剪应变平方
$$\gamma^2/4=\tfrac12\big[(\varepsilon_1-\varepsilon_2)^2n_1^2n_2^2+(\varepsilon_2-\varepsilon_3)^2n_2^2n_3^2+(\varepsilon_3-\varepsilon_1)^2n_3^2n_1^2\big].$$
由 $\sum n_i^2=1$ 可导出在 $(\varepsilon,\gamma/2)$ 平面上的三圆（莫尔圆）。

**物理/几何意义**：$\varepsilon_1$ 最大伸长，$\varepsilon_3$ 最大压缩；最大剪应变 $\gamma_{\max}=\varepsilon_1-\varepsilon_3$。

**极限/特例**：体积不变 $\varepsilon_1+\varepsilon_2+\varepsilon_3=0$。

**典型应用**：应变分析、屈服准则、光弹性。

### 第3章 应力理论与平衡方程

#### 知识点1. Cauchy应力矢量与Cauchy四面体定理

**定义**：过点 $\mathbf x$、外法向 $\mathbf n$ 的截面上的面力密度（单位面积力）为应力矢量
$$\mathbf t^{(\mathbf n)}=\lim_{\Delta S\to0}\frac{\Delta\mathbf F}{\Delta S}.$$
Cauchy 定理断言 $\mathbf t^{(\mathbf n)}=\boldsymbol\sigma\mathbf n$，即 $t_i=\sigma_{ji}n_j$，$\boldsymbol\sigma$ 为二阶应力张量。

**证明（四面体）**：取微小四面体，三个面与坐标面平行，斜法向 $\mathbf n$。设斜面积 $\Delta S$，坐标面积 $\Delta S_i=\Delta S\,n_i$。由动量定理：
$$\rho\Delta V\,\mathbf a=\mathbf t^{(\mathbf n)}\Delta S-\sum_i\mathbf t^{(-\mathbf e_i)}\Delta S_i+\rho\Delta V\,\mathbf f.$$
令体积 $\Delta V\sim\frac13 h\Delta S\to0$（$h\to0$），动量和体力项为高阶小量。用 $\mathbf t^{(-\mathbf e_i)}=-\mathbf t^{(\mathbf e_i)}$（作用反作用），得
$$\mathbf t^{(\mathbf n)}=\sum_i n_i\,\mathbf t^{(\mathbf e_i)}.$$
定义 $\mathbf t^{(\mathbf e_j)}=\sigma_{ij}\mathbf e_i$（即向量 $t_i^{(\mathbf e_j)}=\sigma_{ij}$），则
$$t_i^{(\mathbf n)}=n_j\sigma_{ij}=\sigma_{ij}n_j\ \Rightarrow\ \mathbf t=\boldsymbol\sigma^{\mathrm T}\mathbf n.$$
习惯将 $\sigma_{ij}$ 的第一个指标记为力方向、第二个为面法向，则 $\mathbf t=\boldsymbol\sigma\mathbf n$。$\blacksquare$

**物理/几何意义**：一点的应力状态完全由九个分量确定；斜截面应力是张量与法向的线性作用。

**极限/特例**：法向 $\mathbf e_j$ 上 $\mathbf t^{(\mathbf e_j)}=\sigma_{ij}\mathbf e_i$。

**典型应用**：所有应力分析与强度计算的基础。

#### 知识点2. 应力张量的对称性

**公式**：$\sigma_{ij}=\sigma_{ji}$。

**证明（角动量守恒）**：对微元取矩。应力偶矩净值为 $\epsilon_{ijk}\sigma_{jl}\,dV$（对体力矩、惯性力矩为高阶小量）。由角动量定理，对微元 $\int\epsilon_{ijk}x_j t_k\,dS+\int\epsilon_{ijk}x_j\rho f_k\,dV=\frac{d}{dt}\int\epsilon_{ijk}x_j\rho v_k\,dV$，化简到点即
$$\epsilon_{ijk}\sigma_{jk}=0.$$
由 $\epsilon$ 反对称，得 $\sigma_{jk}=\sigma_{kj}$。$\blacksquare$

**物理/几何意义**：无体力偶时有互等剪应力；应力张量对称，独立分量 $6$ 个。

**极限/特例**：存在体力偶（极性介质）时不成立。

**典型应用**：主应力分析、本构对称性。

#### 知识点3. 平衡微分方程与运动方程

**公式**：静力平衡
$$\sigma_{ji,j}+f_i=0,$$
运动方程
$$\sigma_{ji,j}+f_i=\rho\ddot u_i,$$
其中 $f_i$ 为体力的分量密度。

**推导**：对任意控制体 $V$，动量平衡
$$\int_V\rho\ddot u_i\,dV=\int_{\partial V}\sigma_{ji}n_j\,dS+\int_V f_i\,dV.$$
由散度定理 $\int_{\partial V}\sigma_{ji}n_j\,dS=\int_V\sigma_{ji,j}\,dV$，得
$$\int_V(\rho\ddot u_i-\sigma_{ji,j}-f_i)\,dV=0.$$
由 $V$ 任意且被积函数连续，得 $\sigma_{ji,j}+f_i=\rho\ddot u_i$。$\blacksquare$

**物理/几何意义**：内力散度与体力/惯性力平衡。

**极限/特例**：静止 $\ddot u=0$ 得平衡方程；忽略体力时 $\nabla\cdot\boldsymbol\sigma=\rho\ddot{\mathbf u}$。

**典型应用**：弹性动力学、静力学边值问题。

#### 知识点4. 主应力、不变量与应力莫尔圆

**定义**：满足 $\sigma_{ji}n_j=\sigma n_i$（即 $(\sigma_{ij}-\sigma\delta_{ij})n_j=0$）的非零 $\mathbf n$ 给出主应力 $\sigma$，特征方程
$$\det(\sigma_{ij}-\sigma\delta_{ij})=0\Rightarrow-\sigma^3+J_1\sigma^2-J_2\sigma+J_3=0.$$

**公式与推导**：
$$J_1=\sigma_{ii},\quad J_2=\tfrac12(\sigma_{ii}\sigma_{jj}-\sigma_{ij}\sigma_{ji}),\quad J_3=\det(\sigma_{ij}).$$
主轴系中 $\sigma=\mathrm{diag}(\sigma_1,\sigma_2,\sigma_3)$。斜截面正应力 $\sigma_n=\sigma_in_i^2$，剪应力
$$\tau^2=\sigma_i^2n_i^2-(\sigma_in_i^2)^2.$$
令三主应力给定，$n_i^2$ 的变化在 $(\sigma_n,\tau)$ 平面给出三个圆（莫尔圆）：
$$\left(\sigma_n-\frac{\sigma_2+\sigma_3}{2}\right)^2+\tau^2=\left(\frac{\sigma_2-\sigma_3}{2}\right)^2\ (\text{圆A})+\cdots$$

**证明**：由 $n_1^2+n_2^2+n_3^2=1$ 及 $\sigma_n=\sum\sigma_in_i^2$，视 $n_i^2$ 为未知量解线性方程组，非负条件 $n_i^2\ge0$ 给出三圆区域。

**物理/几何意义**：最大剪应力 $\tau_{\max}=\frac{\sigma_1-\sigma_3}{2}$，位于与 $\sigma_1,\sigma_3$ 主方向成 $45^\circ$ 的截面。

**极限/特例**：静水压 $\sigma_1=\sigma_2=\sigma_3$ 时莫尔圆退化为一点。

**典型应用**：Tresca 屈服准则、岩石/混凝土强度。

#### 知识点5. 应力偏量、von Mises等效应力

**定义**：
$$s_{ij}=\sigma_{ij}-p\delta_{ij},\qquad p=\tfrac13\sigma_{kk}.$$
偏应力不变量
$$J_2=\tfrac12 s_{ij}s_{ij},\qquad J_3=\tfrac13 s_{ij}s_{jk}s_{ki}.$$

**公式与推导**：
(1) $J_2$ 用主应力表示：
$$J_2=\tfrac16\big[(\sigma_1-\sigma_2)^2+(\sigma_2-\sigma_3)^2+(\sigma_3-\sigma_1)^2\big].$$
**证明**：$\because s_{ii}=0$，有 $s_{ij}s_{ij}=(\sigma_{ij}-p\delta_{ij})(\sigma_{ij}-p\delta_{ij})=\sigma_{ij}\sigma_{ij}-3p^2$。而 $\sigma_{ij}\sigma_{ij}=\sigma_1^2+\sigma_2^2+\sigma_3^2$，$p=(\sigma_1+\sigma_2+\sigma_3)/3$。代入 $J_2=\frac12(\sigma_1^2+\sigma_2^2+\sigma_3^2)-\frac16(\sigma_1+\sigma_2+\sigma_3)^2=\frac16[(\sigma_1-\sigma_2)^2+(\sigma_2-\sigma_3)^2+(\sigma_3-\sigma_1)^2]$。
(2) von Mises 等效应力
$$\sigma_{\mathrm{eq}}=\sqrt{3J_2}=\sqrt{\tfrac12[(\sigma_1-\sigma_2)^2+(\sigma_2-\sigma_3)^2+(\sigma_3-\sigma_1)^2]}.$$
单轴拉伸 $\sigma_1=\sigma_0,\sigma_2=\sigma_3=0$ 时 $\sigma_{\mathrm{eq}}=\sigma_0$，故可作屈服量度。

**物理/几何意义**：$J_2$ 表征畸变能；$\sigma_{\mathrm{eq}}$ 把复杂应力等效为单轴应力。

**极限/特例**：纯剪切 $\sigma_1=-\sigma_3=\tau$ 时 $\sigma_{\mathrm{eq}}=\sqrt3\tau$。

**典型应用**：金属屈服（Mises 准则）、疲劳等效应力。

#### 知识点6. 应力圆与坐标变换（转轴公式）

**公式**：转轴 $Q$ 下
$$\sigma'_{ij}=Q_{ik}Q_{jl}\sigma_{kl}.$$
平面应力转轴角 $\theta$：
$$\sigma_{x'}=\frac{\sigma_x+\sigma_y}{2}+\frac{\sigma_x-\sigma_y}{2}\cos2\theta+\tau_{xy}\sin2\theta,$$
$$\tau_{x'y'}=-\frac{\sigma_x-\sigma_y}{2}\sin2\theta+\tau_{xy}\cos2\theta.$$
**推导**：由 $Q(\theta)$ 代入，用倍角公式化简。第二式对 $\theta$ 取极值给出 $\tan2\theta_p=\frac{2\tau_{xy}}{\sigma_x-\sigma_y}$ 确定主方向。

**物理/几何意义**：应力随观察面旋转按正弦/余弦规律变化，可作莫尔圆图解。

**典型应用**：应变片花、主应力方向测定。

#### 知识点7. 应力张量的分解与八面体应力

**定义**：八面体面法向 $\mathbf n=(1,1,1)/\sqrt3$，其上正应力
$$\sigma_{\mathrm{oct}}=\tfrac13\sigma_{kk}=p,$$
剪应力
$$\tau_{\mathrm{oct}}=\sqrt{\tfrac23J_2}=\tfrac13\sqrt{(\sigma_1-\sigma_2)^2+(\sigma_2-\sigma_3)^2+(\sigma_3-\sigma_1)^2}.$$
**推导**：$\sigma_n=\sigma_in_i^2=\frac13\sigma_{ii}=p$。$\tau_{\mathrm{oct}}^2=\sigma_i^2n_i^2-\sigma_n^2=\frac13\sigma_i^2-\frac19(\sigma_{ii})^2=\frac23 J_2$。

**物理/几何意义**：八面体剪应力与 $J_2$ 成正比，是塑性理论的自然量。

**典型应用**：Mises 屈服、等效塑性应变。

#### 知识点8. 有限变形中的应力张量

**定义**：设参考构型面元 $N\,dA$，力 $d\mathbf P$。
(1) 第一类 Piola–Kirchhoff 应力 $P_{iJ}$：$dP_i=P_{iJ}N_J\,dA$，两点张量。
(2) 第二类 Piola–Kirchhoff 应力 $S_{IJ}$：由 $P_{iJ}=F_{iK}S_{KJ}$ 定义。
(3) Cauchy 应力 $\sigma_{ij}$：$dP_i=\sigma_{ij}n_j\,da$。

**关系推导**：由 Nanson 公式 $n_j\,da=J F^{-T}_{jK}N_K\,dA$，代入 Cauchy 应力：
$$dP_i=\sigma_{ij}\,n_j\,da=\sigma_{ij} J F^{-T}_{jK}N_K\,dA\ \Rightarrow\ P_{iK}=J\sigma_{ij}F^{-T}_{jK}=J F^{-T}_{Kj}\sigma_{ij}.$$
又 $P_{iJ}=F_{iK}S_{KJ}$，故
$$S_{IJ}=J F^{-1}_{Ii}F^{-1}_{Jj}\sigma_{ij},\qquad \sigma_{ij}=J^{-1}F_{iI}F_{jJ}S_{IJ}.$$

**物理/几何意义**：$P$ 混合参考面与现时力（非对称）；$S$ 完全参考型且对称，与格林应变功共轭。

**极限/特例**：小变形下三者趋于同一 $\sigma_{ij}$。

**典型应用**：有限元 Total/Updated Lagrangian 格式、超弹性。

### 第4章 连续介质守恒定律与Navier–Stokes方程

#### 知识点1. 雷诺输运定理

**定义**：对随物质运动的控制体，标量/张量场 $\phi$ 的体积分变化率
$$\frac{d}{dt}\int_{V(t)}\phi\,dV=\int_V\frac{\partial\phi}{\partial t}\,dV+\oint_{\partial V}\phi\,\mathbf v\cdot\mathbf n\,dS.$$

**证明**：用物质坐标变换 $\int_{V(t)}\phi(\mathbf x,t)\,dV=\int_{V_0}\phi(\boldsymbol\varphi(\mathbf X,t),t)J\,dV_0$，对 $t$ 求导：
$$\frac{d}{dt}\int_{V_0}\phi J\,dV_0=\int_{V_0}\left(\frac{D\phi}{Dt}J+\phi\frac{DJ}{Dt}\right)dV_0.$$
由 $\frac{DJ}{Dt}=J\,\nabla\cdot\mathbf v$（Euler 展开式），得
$$=\int_{V_0}\left(\frac{D\phi}{Dt}+\phi\nabla\cdot\mathbf v\right)J\,dV_0=\int_{V(t)}\left(\frac{D\phi}{Dt}+\phi\nabla\cdot\mathbf v\right)dV.$$
又 $\frac{D\phi}{Dt}+\phi\nabla\cdot\mathbf v=\frac{\partial\phi}{\partial t}+\nabla\cdot(\phi\mathbf v)$，再由散度定理得所证。$\blacksquare$

**物理/几何意义**：体积分变化率 = 局部（固定体积）变化 + 表面对流通量。

**极限/特例**：定常不可压时 $\frac{d}{dt}\int\phi\,dV=\oint\phi\mathbf v\cdot\mathbf n\,dS$。

**典型应用**：质量、动量、能量守恒方程的统一起点。

#### 知识点2. 质量守恒与连续性方程

**公式**：
$$\frac{D\rho}{Dt}+\rho\nabla\cdot\mathbf v=0,\qquad \frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf v)=0.$$

**推导**：质量守恒 $\frac{d}{dt}\int_V\rho\,dV=0$。用雷诺输运定理：
$$0=\int_V\left(\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf v)\right)dV\ \Rightarrow\ \frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf v)=0.$$
展开 $\nabla\cdot(\rho\mathbf v)=\mathbf v\cdot\nabla\rho+\rho\nabla\cdot\mathbf v$，得物质导数形式
$$\frac{D\rho}{Dt}+\rho\nabla\cdot\mathbf v=0,\quad \frac{D(\cdot)}{Dt}=\frac{\partial(\cdot)}{\partial t}+\mathbf v\cdot\nabla(\cdot).$$

**物理/几何意义**：密度变化与体积膨胀/压缩对应。

**极限/特例**：不可压缩 $\frac{D\rho}{Dt}=0\Rightarrow\nabla\cdot\mathbf v=0$。固体常用 $\rho J=\rho_0$。

**典型应用**：流体、多孔介质、爆炸波。

#### 知识点3. 动量守恒与Cauchy运动方程

**公式**：
$$\nabla\cdot\boldsymbol\sigma+\mathbf f=\rho\frac{D\mathbf v}{Dt}\quad(\text{分量 } \sigma_{ji,j}+f_i=\rho\dot v_i).$$

**推导**：对控制体 $V$ 用动量守恒
$$\frac{d}{dt}\int_V\rho v_i\,dV=\oint_{\partial V}\sigma_{ji}n_j\,dS+\int_V f_i\,dV.$$
左边由雷诺输运并用连续性方程化简为 $\int_V\rho\frac{Dv_i}{Dt}\,dV$：
$$\frac{d}{dt}\int_V\rho v_i\,dV=\int_V\left[\frac{\partial(\rho v_i)}{\partial t}+\nabla\cdot(\rho v_i\mathbf v)\right]dV=\int_V\rho\frac{Dv_i}{Dt}\,dV.$$
右边用散度定理，得 $\int_V(\rho Dv_i/Dt-\sigma_{ji,j}-f_i)dV=0$，由 $V$ 任意得运动方程。$\blacksquare$

**物理/几何意义**：牛顿第二定律在连续体中的局部形式。

**极限/特例**：静力平衡、固体动力学的统一形式。

**典型应用**：结构动力学、流体力学、地球物理。

#### 知识点4. 能量守恒与热力学第一定律

**公式**：控制体能量守恒
$$\frac{d}{dt}\int_V\left(\tfrac12\rho v_i^2+\rho e\right)dV=\oint_{\partial V}v_i\sigma_{ji}n_j\,dS+\int_V f_iv_i\,dV-\oint_{\partial V}q_in_i\,dS+\int_V r\,dV,$$
其中 $e$ 比内能，$q_i$ 热流，$r$ 体热源。

**推导（局部形式）**：用雷诺输运、散度定理与运动方程消去动能项，得
$$\rho\frac{De}{Dt}=\sigma_{ij}D_{ij}-\nabla\cdot\mathbf q+r,$$
其中 $D_{ij}$ 为变形率张量。
**关键步骤**：动能方程
$$\frac{D}{Dt}\left(\tfrac12\rho v^2\right)=v_i(\sigma_{ji,j}+f_i).$$
应力做功项 $\oint v_i\sigma_{ji}n_j=\int[(v_i\sigma_{ji})_{,j}]=\int(v_{i,j}\sigma_{ji}+v_i\sigma_{ji,j})$。两式相减并用 $\sigma_{ij}D_{ij}$ 得内能方程。$\blacksquare$

**物理/几何意义**：内能增加 = 应力做功率（机械耗散/储存）+ 热源 - 热流散度。

**极限/特例**：弹性储能时 $\sigma_{ij}D_{ij}=\dot W$ 为应变能率；黏性时含耗散。

**典型应用**：热弹性、塑性耗散、Navier–Stokes 能量方程。

#### 知识点5. 熵不等式与Clausius–Duhem不等式

**定义**：Clausius–Duhem 不等式（等温或一般热力学第二定律局部形式）
$$\rho\frac{D\eta}{Dt}\ge \frac{r}{T}-\nabla\cdot\left(\frac{\mathbf q}{T}\right)\ \Rightarrow\ \rho T\dot\eta-\sigma_{ij}D_{ij}+\frac{q_i}{T}\partial_iT\ge0.$$

**推导**：由 $\rho\dot e=\sigma_{ij}D_{ij}-\nabla\cdot\mathbf q+r$ 及 $\psi=e-T\eta$（Helmholtz 自由能），代入熵不等式，整理得耗散不等式
$$\sigma_{ij}D_{ij}-\rho(\dot\psi+\eta\dot T)-\frac{q_i}{T}\partial_iT\ge0.$$

**物理/几何意义**：限制本构关系，保证热力学一致性（不产生负熵）。

**极限/特例**：弹性材料在等温下 $\sigma_{ij}=\rho\partial\psi/\partial\varepsilon_{ij}$；Fourier 导热 $q_i=-k\partial_iT$ 使 $-q_i\partial_iT/T\ge0$ 要求 $k\ge0$。

**典型应用**：本构公理、耗散型材料建模。

#### 知识点6. Newton流体本构与Navier–Stokes方程

**定义（牛顿流体）**：偏应力线性正比于变形率，
$$\sigma_{ij}=-p\delta_{ij}+2\mu D_{ij}+\lambda D_{kk}\delta_{ij},$$
其中 $\mu$ 动力黏度，$\lambda$ 第二黏度。

**Stokes 假设**：$\lambda+\frac23\mu=0$（体积黏性为零），于是
$$\sigma_{ij}=-p\delta_{ij}+2\mu\left(D_{ij}-\tfrac13D_{kk}\delta_{ij}\right).$$

**推导 Navier–Stokes**：将本构代入运动方程 $\sigma_{ji,j}+f_i=\rho Dv_i/Dt$：
$$\partial_j\left[-p\delta_{ij}+\mu(\partial_jv_i+\partial_iv_j)+\lambda\partial_kv_k\delta_{ij}\right]+f_i=\rho\frac{Dv_i}{Dt}.$$
若 $\mu,\lambda$ 为常数，
$$-p_{,i}+\mu\nabla^2v_i+\mu\partial_i(\partial_kv_k)+\lambda\partial_i(\partial_kv_k)+f_i=\rho\frac{Dv_i}{Dt},$$
即
$$\rho\frac{D\mathbf v}{Dt}=-\nabla p+\mu\nabla^2\mathbf v+(\mu+\lambda)\nabla(\nabla\cdot\mathbf v)+\mathbf f.$$
用 Stokes 假设 $\mu+\lambda=\frac13\mu$... 更常用不可压 $\nabla\cdot\mathbf v=0$：
$$\rho\frac{D\mathbf v}{Dt}=-\nabla p+\mu\nabla^2\mathbf v+\mathbf f,\qquad \nabla\cdot\mathbf v=0.$$

**物理/几何意义**：惯性力、压力梯度、黏性扩散、体力平衡。

**极限/特例**：$\mu=0$ 得 Euler 方程；$Re\to0$ 得 Stokes 流；层流管流得 Hagen–Poiseuille 抛物线速度。

**典型应用**：空气动力学、润滑、血流；同时是固体动量方程的特例。

#### 知识点7. 固体场方程体系与初边值条件

**基本方程（小变形线弹性静力）**：
$$\begin{aligned}\text{平衡：}&\quad \sigma_{ji,j}+f_i=0,\\ \text{几何：}&\quad \varepsilon_{ij}=\tfrac12(u_{i,j}+u_{j,i}),\\ \text{本构：}&\quad \sigma_{ij}=C_{ijkl}\varepsilon_{kl}.\end{aligned}$$

**边界条件**：
$$u_i=\bar u_i\ (\partial B_u),\qquad \sigma_{ji}n_j=\bar t_i\ (\partial B_t),\qquad \partial B_u\cup\partial B_t=\partial B,\ \partial B_u\cap\partial B_t=\varnothing.$$

**推导（位移法）**：将几何代入本构再代入平衡，得位移形式
$$C_{ijkl}u_{k,lj}+f_i=0.$$
各向同性时 $C_{ijkl}=\lambda\delta_{ij}\delta_{kl}+\mu(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk})$：
$$\mu u_{i,jj}+(\lambda+\mu)u_{j,ji}+f_i=0,$$
即 Lamé–Navier 方程 $(\lambda+\mu)\nabla(\nabla\cdot\mathbf u)+\mu\nabla^2\mathbf u+\mathbf f=0$。

**物理/几何意义**：三类场量的联立封闭。

**极限/特例**：平面应力/平面应变约化；轴对称。

**典型应用**：所有弹性静力学解析与数值求解。

#### 知识点8. 物质导数、对流与随体守恒的等价形式

**定义**：物质导数 $\frac{D\phi}{Dt}=\frac{\partial\phi}{\partial t}+v_i\partial_i\phi$。

**推导与应用**：
(1) $\frac{DJ}{Dt}=J\nabla\cdot\mathbf v$。
**证明**：$J=\det F$，$d(\det F)=(\det F)F^{-T}_{Jj}dF_{jJ}$，故 $\frac{DJ}{Dt}=J F^{-T}_{Jj}\partial_Jv_j=J\partial_jv_j$。
(2) 质量守恒 $\rho J=\rho_0$。
**证明**：由 $\frac{D}{Dt}(\rho J)=J\frac{D\rho}{Dt}+\rho J\nabla\cdot\mathbf v=J(\frac{D\rho}{Dt}+\rho\nabla\cdot\mathbf v)=0$。
(3) 张量随体形式（如应力）：$\dot\sigma_{ij}=\frac{D\sigma_{ij}}{Dt}$ 不是客观张量率，客观率需减去旋率（Jaumann、Green–Naghdi）。

**物理/几何意义**：区分 Lagrange 与 Euler 描述，客观性要求。

**极限/特例**：定常 $\partial_t=0$ 时 $D/Dt=\mathbf v\cdot\nabla$。

**典型应用**：大变形本构、流固耦合、Abaqus/ANSYS 中 Lagrangian/Eulerian 算法。

### 第5章 弹性本构关系与材料对称性

#### 知识点1. 广义胡克定律与弹性张量

**定义**：线性弹性本构
$$\sigma_{ij}=C_{ijkl}\varepsilon_{kl},\qquad \varepsilon_{ij}=S_{ijkl}\sigma_{kl},$$
$C_{ijkl}$ 为弹性刚度张量，$S_{ijkl}=C^{-1}_{ijkl}$ 为柔度张量。

**假设**：小变形、等温、应力应变为线性、无初始应力。

**推导对称性**：
(1) 由 $\sigma_{ij}=\sigma_{ji}$ 得 $C_{ijkl}=C_{jikl}$；
(2) 由 $\varepsilon_{kl}=\varepsilon_{lk}$ 得 $C_{ijkl}=C_{ijlk}$；
(3) 超弹性（应变能函数存在）：$\sigma_{ij}=\partial W/\partial\varepsilon_{ij}$，$W=\frac12C_{ijkl}\varepsilon_{ij}\varepsilon_{kl}$，故 $C_{ijkl}=C_{klij}$。
$\therefore$ 一般各向异性材料独立常数由 $81$ 降为 $21$。

**物理/几何意义**：$C$ 是四阶张量，联系应变与应力，编码材料刚度各向异性。

**极限/特例**：各向同性独立常数减为 $2$。

**典型应用**：复合材料层合板、木材、晶体。

#### 知识点2. 各向同性弹性本构与Lamé常数

**公式**：
$$\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij},$$
其中 $\lambda,\mu$ 为 Lamé 常数（$\mu$ 即剪切模量 $G$）。

**推导（由一般形式）**：各向同性四阶张量必为
$$C_{ijkl}=\lambda\delta_{ij}\delta_{kl}+\mu(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}).$$
代入 $\sigma_{ij}=C_{ijkl}\varepsilon_{kl}$：
$$\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij}.$$
**验证对称性**：$C_{ijkl}=C_{klij}$、$C_{ijkl}=C_{jikl}$ 均满足。独立常数 $\lambda,\mu$ 两者。

**物理/几何意义**：$\lambda$ 联系体积与正应力，$\mu$ 联系剪切与剪应变。

**极限/特例**：$\lambda=\mu$ 时某些弹性解析解简化；$\mu=0$ 为无剪切弹性（不稳定）。

**典型应用**：金属、玻璃、各向同性材料的统一描述。

#### 知识点3. 弹性常数关系（E, ν, G, K）

**定义**：工程弹性常数：杨氏模量 $E$、泊松比 $\nu$、剪切模量 $G$、体积模量 $K$。

**推导**：
(1) 单轴应力 $\sigma_{11}=\sigma,\ \sigma_{22}=\sigma_{33}=0$：
$$\varepsilon_{22}=\varepsilon_{33}=-\nu\frac{\sigma}{E},\qquad \varepsilon_{11}=\frac{\sigma}{E}.$$
由 $\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij}$ 反解得
$$E=\frac{\mu(3\lambda+2\mu)}{\lambda+\mu},\qquad \nu=\frac{\lambda}{2(\lambda+\mu)}.$$
**推导**：单轴条件代入本构并解 $\varepsilon_{11},\varepsilon_{22}$ 即得。
(2) 纯剪切得 $G=\mu$。
(3) 静水压 $\sigma_{ij}=-p\delta_{ij}$：$\varepsilon_{kk}=-3p/(3\lambda+2\mu)$，故
$$K=\frac{3\lambda+2\mu}{3}=\frac{E}{3(1-2\nu)}.$$
(4) 反解：
$$\mu=G=\frac{E}{2(1+\nu)},\qquad \lambda=\frac{E\nu}{(1+\nu)(1-2\nu)}.$$

**物理/几何意义**：弹性常数之间存在约束，独立只有 $2$ 个。

**极限/特例**：各向同性稳定性要求 $-1<\nu<1/2$；$\nu=1/2$ 不可压（$K\to\infty$）；$\nu=-1$ 时 $G\to\infty$ 但 $K\to E/9$。

**典型应用**：由拉伸实验测定 $E,\nu$，扭转变形测定 $G$。

#### 知识点4. 应变能密度与正定性

**定义**：线弹性应变能密度
$$W=\tfrac12\sigma_{ij}\varepsilon_{ij}=\tfrac12C_{ijkl}\varepsilon_{ij}\varepsilon_{kl}.$$

**公式与推导**：
(1) $\sigma_{ij}=\partial W/\partial\varepsilon_{ij}$。
**证明**：$\frac{\partial W}{\partial\varepsilon_{ij}}=\frac12C_{klmn}\frac{\partial(\varepsilon_{kl}\varepsilon_{mn})}{\partial\varepsilon_{ij}}=\frac12C_{klmn}(\delta_{ki}\delta_{lj}\varepsilon_{mn}+\varepsilon_{kl}\delta_{mi}\delta_{nj})=\frac12(C_{ijmn}\varepsilon_{mn}+C_{klij}\varepsilon_{kl})=C_{ijkl}\varepsilon_{kl}=\sigma_{ij}$（用了 $C_{ijkl}=C_{klij}$）。
(2) 各向同性分解：
$$W=\tfrac12\lambda(\varepsilon_{kk})^2+\mu\varepsilon_{ij}\varepsilon_{ij}=\tfrac12K(\varepsilon_{kk})^2+2\mu\,e_{ij}e_{ij},\quad e_{ij}=\varepsilon_{ij}-\tfrac13\varepsilon_{kk}\delta_{ij}.$$
**证明**：$\varepsilon_{ij}=e_{ij}+\frac13\varepsilon_{kk}\delta_{ij}$，代入并用 $e_{ii}=0$，交叉项 $\lambda\varepsilon_{kk}\delta_{ij}e_{ij}=0$，得体积项 $\frac12(\lambda+\frac23\mu)(\varepsilon_{kk})^2=\frac12K(\varepsilon_{kk})^2$ 与畸变项 $2\mu e_{ij}e_{ij}$。
(3) 正定性：$W\ge0$ 且仅 $\varepsilon=0$ 时为零 $\Rightarrow K>0,\mu>0$。

**物理/几何意义**：$W$ 分体积能（球）与畸变能（偏）；正定性是材料稳定性的必要条件。

**极限/特例**：单轴 $\frac12E\varepsilon^2$；纯剪 $\frac12G\gamma^2$。

**典型应用**：能量法、Rice $J$ 积分、有限元刚度矩阵。

#### 知识点5. 各向异性与材料对称群

**定义**：材料对称是指存在正交变换 $Q$ 使 $C_{ijkl}=Q_{ip}Q_{jq}Q_{kr}Q_{ls}C_{pqrs}$。

**推导与分类**：
(1) 单斜（一个对称面）：$13$ 独立常数；
(2) 正交各向异性（三个正交对称面）：$9$ 独立常数，
$$C=\begin{pmatrix}C_{11}&C_{12}&C_{13}&0&0&0\\C_{12}&C_{22}&C_{23}&0&0&0\\C_{13}&C_{23}&C_{33}&0&0&0\\0&0&0&C_{44}&0&0\\0&0&0&0&C_{55}&0\\0&0&0&0&0&C_{66}\end{pmatrix}.$$
(3) 横观各向同性（一个对称轴）：$5$ 独立常数；
(4) 各向同性：$2$ 独立常数。
**推导纲要**：令 $Q$ 为关于某平面/轴旋转，代入变换式并令不变，逐一消去分量。Voigt 记法 $\{11,22,33,23,31,12\}\to\{1,\dots,6\}$ 把 $C$ 写成 $6\times6$。

**物理/几何意义**：对称性越高，独立常数越少；正交各向异性在复合材料、木材中常见。

**极限/特例**：横观各向同性沿轴向为各向同性。

**典型应用**：层合板、木材、岩石、晶体弹性。

#### 知识点6. 热弹性本构与Duhamel–Neumann

**公式**：温度变化 $T-T_0=\Delta T$，热应变 $\varepsilon^{th}_{ij}=\alpha\Delta T\delta_{ij}$，本构
$$\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij}-(3\lambda+2\mu)\alpha\Delta T\,\delta_{ij}.$$

**推导**：总应变 $\varepsilon_{ij}=\varepsilon^{el}_{ij}+\alpha\Delta T\delta_{ij}$，弹性部分 $\sigma_{ij}=\lambda\varepsilon^{el}_{kk}\delta_{ij}+2\mu\varepsilon^{el}_{ij}$，代入得
$$\sigma_{ij}=\lambda(\varepsilon_{kk}-3\alpha\Delta T)\delta_{ij}+2\mu(\varepsilon_{ij}-\alpha\Delta T\delta_{ij})=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij}-(3\lambda+2\mu)\alpha\Delta T\delta_{ij}.$$
无约束自由膨胀时 $\sigma_{ij}=0\Rightarrow\varepsilon_{ij}=\alpha\Delta T\delta_{ij}$。

**物理/几何意义**：温度引起初应变，受约束时产生热应力。

**极限/特例**：全约束 $\varepsilon_{ij}=0$ 时 $\sigma_{ij}=-3K\alpha\Delta T\delta_{ij}$（三向等压）。

**典型应用**：热障涂层、焊接残余应力、航空发动机热端部件。

#### 知识点7. 超弹性（Hyperelasticity）与本构客观性

**定义**：存在应变能函数 $W(\mathbf C)$（$\mathbf C=\mathbf F^{\mathrm T}\mathbf F$）使
$$\mathbf S=2\frac{\partial W}{\partial\mathbf C},\qquad S_{IJ}=2\frac{\partial W}{\partial C_{IJ}}.$$

**推导**：由功共轭 $\dot W=S_{IJ}\dot E_{IJ}=\frac12S_{IJ}\dot C_{IJ}$，对任意速度成立，故 $S_{IJ}=2\partial W/\partial C_{IJ}$。转换为 Cauchy 应力：
$$\sigma_{ij}=J^{-1}F_{iI}F_{jJ}S_{IJ}=J^{-1}F_{iI}F_{jJ}2\frac{\partial W}{\partial C_{IJ}}.$$
客观性要求 $W$ 只通过 $\mathbf C$（或主伸长）依赖变形，故 $W(\mathbf C)=W(\mathbf R^{\mathrm T}\mathbf C\mathbf R)$。各向同性时 $W=W(I_1,I_2,I_3)$，$I_1=\operatorname{tr}\mathbf C,I_2=\frac12[(tr\mathbf C)^2-tr\mathbf C^2],I_3=\det\mathbf C$。
Neo-Hookean：$W=\frac{\mu}{2}(I_1-3)-\mu\ln J+\frac{\lambda}{2}(\ln J)^2$。Mooney–Rivlin：$W=C_{10}(I_1-3)+C_{01}(I_2-3)$。

**物理/几何意义**：弹性变形可逆、无耗散，应力由自由能梯度给出。

**极限/特例**：小变形展开即得线弹性 $C_{ijkl}$。

**典型应用**：橡胶、生物软组织、封装材料。

### 第6章 弹性力学边值问题、唯一性与基本解

#### 知识点1. 边值问题提法与适定性

**定义**：求位移场 $\mathbf u$ 满足
$$\begin{cases}\sigma_{ji,j}+f_i=0,& \mathbf x\in B,\\ \sigma_{ij}=C_{ijkl}\varepsilon_{kl},& \mathbf x\in B,\\ \varepsilon_{ij}=\frac12(u_{i,j}+u_{j,i}),& \mathbf x\in B,\\ u_i=\bar u_i\ (\partial B_u),\quad \sigma_{ji}n_j=\bar t_i\ (\partial B_t).\end{cases}$$

**推导封闭性**：未知量 $u_i(3)+\sigma_{ij}(6)+\varepsilon_{ij}(6)=15$；方程 $3+6+6=15$，封闭。

**物理/几何意义**：位移边界与力边界互补（不相交、覆盖整个边界）。

**极限/特例**：全位移边界为 Dirichlet 问题；全力边界为 Neumann 问题（需满足整体平衡）。

**典型应用**：结构强度、接触、约束。

#### 知识点2. 位移法：Lamé–Navier方程

**公式**：各向同性无体力
$$(\lambda+\mu)\nabla(\nabla\cdot\mathbf u)+\mu\nabla^2\mathbf u+\mathbf f=\mathbf 0.$$
**推导**：由 $\sigma_{ij}=\lambda u_{k,k}\delta_{ij}+\mu(u_{i,j}+u_{j,i})$ 代入 $\sigma_{ij,j}+f_i=0$：
$$\lambda u_{k,ki}+\mu(u_{i,jj}+u_{j,ij})+f_i=0\Rightarrow(\lambda+\mu)u_{j,ji}+\mu u_{i,jj}+f_i=0.$$

**物理/几何意义**：用位移作基本未知量，适合给定位移边界。

**极限/特例**：Helmholtz 分解 $\mathbf u=\nabla\phi+\nabla\times\boldsymbol\psi$ 可解耦为 $\nabla^2\phi\propto\nabla\cdot\mathbf u$ 与 $\nabla^2\boldsymbol\psi\propto\nabla\times\mathbf u$。

**典型应用**：弹性静力学有限元、Boussinesq/Navier 问题。

#### 知识点3. 应力法：Beltrami–Michell方程

**公式**：无体力时，应力满足
$$\nabla^2\sigma_{ij}+\frac{1}{1+\nu}\sigma_{kk,ij}=0.$$
**推导**：由协调方程 $\varepsilon_{ij,kl}+\varepsilon_{kl,ij}-\varepsilon_{ik,jl}-\varepsilon_{jl,ik}=0$ 与柔度 $\varepsilon_{ij}=\frac{1+\nu}{E}\sigma_{ij}-\frac{\nu}{E}\sigma_{kk}\delta_{ij}$ 代入。利用平衡方程消去双力项：
对 $\nabla^2\sigma_{ij}$ 展开，并用 $\sigma_{ij,jj}=0$、$\sigma_{kk,jj}=0$（由平衡取散度）得 Beltrami 形式。更完整地：
$$\nabla^2\sigma_{ij}+\frac{1}{1+\nu}\sigma_{kk,ij}=0.$$

**物理/几何意义**：以应力为未知量，需额外满足协调；适合力边界。

**极限/特例**：平面问题化为 $\nabla^4(U)=0$（Airy）。

**典型应用**：解析弹性力学、光弹性。

#### 知识点4. 叠加原理

**陈述**：线弹性边值问题解对载荷线性叠加：若 $\mathbf u^{(1)}$ 对应 $(\mathbf f^{(1)},\bar{\mathbf t}^{(1)},\bar{\mathbf u}^{(1)})$，$\mathbf u^{(2)}$ 对应 $(\mathbf f^{(2)},\bar{\mathbf t}^{(2)},\bar{\mathbf u}^{(2)})$，则 $\mathbf u^{(1)}+\mathbf u^{(2)}$ 对应载荷之和。

**证明**：控制方程、本构、边界条件均线性，故线性组合仍满足。$\blacksquare$

**物理/几何意义**：复杂载荷可分解为简单载荷之和。

**极限/特例**：几何非线性、塑性、接触问题不适用。

**典型应用**：影响线、Green 函数、模态叠加。

#### 知识点5. 唯一性定理（Kirchhoff）

**陈述**：在给定体力、位移边界 $\partial B_u$、面力边界 $\partial B_t$ 下，线弹性静力解唯一（相差刚体位移仅在纯力边界时出现）。

**证明（能量法）**：设两组解 $\mathbf u^{(1)},\mathbf u^{(2)}$，差 $\Delta\mathbf u=\mathbf u^{(1)}-\mathbf u^{(2)}$ 对应零体力、零面力、零给定位移。其应变能为
$$U=\frac12\int_B\Delta\varepsilon_{ij}C_{ijkl}\Delta\varepsilon_{kl}\,dV\ge0.$$
由虚功/能量守恒，$U=0$。因 $C$ 正定，$\Delta\varepsilon_{ij}=0$，故 $\Delta\mathbf u$ 为刚体位移。在 $\partial B_u$ 上为零，故刚体位移为零。$\therefore$ 唯一。$\blacksquare$

**物理/几何意义**：弹性问题适定性，是数值方法与叠加原理的基础。

**极限/特例**：纯力边界解唯一到刚体位移；非保守/不稳定材料可能非唯一。

**典型应用**：证明圣维南原理的合理性、极限分析。

#### 知识点6. 圣维南原理及能量估计

**陈述**：作用于弹性体局部表面、静力等效的两组载荷，在远离该区域处引起的应力差别可忽略（随距离指数衰减）。

**推导（衰减估计思路）**：对半无限带/柱，将位移展开为特征模态，远离端部时高阶模态 $\propto e^{-\lambda_n x}$ 指数衰减，最低非零自平衡模态衰减最慢。Toupin 给出能量衰减不等式
$$E(x)\le E(0)e^{-2kx},\quad k>0.$$
**证明纲要**：以应变能 $E(x)=\int_{x}^{\infty}\int_A C_{ijkl}\varepsilon_{ij}\varepsilon_{kl}\,dA\,d\xi$ 与截面能量流不等式，用 Poincaré 型不等式得 $\frac{dE}{dx}\le-2kE$。

**物理/几何意义**：局部作用影响局部，工程上可用等效力系简化。

**极限/特例**：薄壁、长细结构衰减长度大，需谨慎。

**典型应用**：梁端部加载、局部支承简化、有限元加载。

#### 知识点7. 厚壁圆筒（Lamé解）

**问题**：内半径 $a$、外半径 $b$，内压 $p_i$、外压 $p_o$，平面应变或轴对称。

**公式与推导**：轴对称平衡
$$\frac{d\sigma_r}{dr}+\frac{\sigma_r-\sigma_\theta}{r}=0.$$
几何 $\varepsilon_r=du/dr,\ \varepsilon_\theta=u/r$。用位移法（$\sigma_r=\lambda(\varepsilon_r+\varepsilon_\theta)+2\mu\varepsilon_r$ 等）代入平衡式，整理得
$$\frac{d}{dr}\left[\frac1r\frac{d(ru)}{dr}\right]=0\Rightarrow u=Ar+\frac{B}{r},$$
于是 $\varepsilon_r=A-B/r^2,\ \varepsilon_\theta=A+B/r^2$。注意到 $\varepsilon_r+\varepsilon_\theta=2A$，
$$\sigma_r=\lambda(2A)+2\mu(A-B/r^2)=\frac{E}{1-2\nu}A-\frac{E}{1+\nu}\frac{B}{r^2},$$
$$\sigma_\theta=\lambda(2A)+2\mu(A+B/r^2)=\frac{E}{1-2\nu}A+\frac{E}{1+\nu}\frac{B}{r^2}.$$
标准结果（无体力）：
$$\sigma_r=\frac{a^2p_i-b^2p_o}{b^2-a^2}-\frac{(p_i-p_o)a^2b^2}{r^2(b^2-a^2)},$$
$$\sigma_\theta=\frac{a^2p_i-b^2p_o}{b^2-a^2}+\frac{(p_i-p_o)a^2b^2}{r^2(b^2-a^2)}.$$
**推导细节**：由边界 $\sigma_r(a)=-p_i,\ \sigma_r(b)=-p_o$ 解出
$$\frac{E}{1-2\nu}A=\frac{a^2p_i-b^2p_o}{b^2-a^2},\qquad \frac{E}{1+\nu}B=\frac{(p_i-p_o)a^2b^2}{b^2-a^2},$$
回代即得上面两式。
轴向应力（闭端）$\sigma_z=\frac{a^2p_i-b^2p_o}{b^2-a^2}$（平面应变时另有 $\sigma_z=\nu(\sigma_r+\sigma_\theta)$）。

**物理/几何意义**：内壁 $\sigma_\theta$ 最大，$r=a$ 处 $\sigma_\theta=p_i\frac{b^2+a^2}{b^2-a^2}$（$p_o=0$）。

**极限/特例**：薄壁 $b\to a$ 得 $\sigma_\theta\approx p_ir/t$；$a\to0$ 得实心轴。

**典型应用**：压力容器、炮筒、过盈配合。

#### 知识点8. 圆柱扭转与翘曲

**问题**：圆轴受端扭矩 $M_t$，Saint–Venant 半逆解 $\mathbf u=(u,v,w)$，设 $w=\theta\psi(x,y)$ 为翘曲，$u=-\theta yz,\ v=\theta xz$。

**公式与推导**：
(1) 应变 $\varepsilon_{xz}=\frac12(\partial_xw+\partial_zu)=\frac{\theta}{2}(\partial_x\psi-y),\ \varepsilon_{yz}=\frac{1}{2}(\partial_y\psi+x)$。
(2) 应力 $\sigma_{xz}=G\theta(\partial_x\psi-y),\ \sigma_{yz}=G\theta(\partial_y\psi+x)$。
(3) 平衡 $\partial_x\sigma_{xz}+\partial_y\sigma_{yz}=0\Rightarrow\nabla^2\psi=0$（翘曲调和）。
(4) 边界无面力 $\Rightarrow\frac{\partial\psi}{\partial n}=y n_x-x n_y$。
(5) 扭矩 $M_t=\int(x\sigma_{yz}-y\sigma_{xz})dA=G\theta J$，$J$ 为极惯性矩（圆截面翘曲 $\psi=0$）。
圆形截面 $\psi=0$，$\tau=\frac{M_t r}{J}$，$J=\frac{\pi R^4}{2}$。
非圆截面翘曲非零，需解 Poisson 方程。

**物理/几何意义**：圆形截面不翘曲，非圆截面翘曲；自由扭转刚度 $GJ$。

**极限/特例**：薄壁闭口管 $\tau=\frac{M_t}{2A_m t}$（Bredt 公式）。

**典型应用**：轴类零件、机翼盒段抗扭。

### 第7章 平面问题与复变函数方法

#### 知识点1. 平面应力与平面应变

**定义**：
(1) 平面应力：薄板，$\sigma_{zz}=\sigma_{xz}=\sigma_{yz}=0$，仅 $\sigma_{xx},\sigma_{yy},\sigma_{xy}$ 非零。
(2) 平面应变：长柱，$\varepsilon_{zz}=\varepsilon_{xz}=\varepsilon_{yz}=0$，$w=0$，$u=u(x,y),v=v(x,y)$。

**公式与推导**：
平面应力：
$$\varepsilon_{xx}=\frac1E(\sigma_{xx}-\nu\sigma_{yy}),\quad \varepsilon_{yy}=\frac1E(\sigma_{yy}-\nu\sigma_{xx}),\quad \gamma_{xy}=\frac{\tau_{xy}}{G},$$
$$\varepsilon_{zz}=-\frac{\nu}{E}(\sigma_{xx}+\sigma_{yy})\ (\ne0).$$
平面应变：$\varepsilon_{zz}=0$ 得 $\sigma_{zz}=\nu(\sigma_{xx}+\sigma_{yy})$，且
$$\varepsilon_{xx}=\frac{1+\nu}{E}[(1-\nu)\sigma_{xx}-\nu\sigma_{yy}],\quad \gamma_{xy}=\frac{\tau_{xy}}{G}.$$
**转换**：平面应变弹性常数替换 $E\to\frac{E}{1-\nu^2},\ \nu\to\frac{\nu}{1-\nu}$。
**推导**：由三维本构令对应分量为零并消去 $\sigma_{zz}$。

**物理/几何意义**：平面应力厚度方向自由，平面应变厚度方向受约束。

**极限/特例**：薄板平面应力，厚截面平面应变。

**典型应用**：板、长坝、隧道横截面。

#### 知识点2. Airy应力函数与双调和方程

**定义**：无体力平面问题取
$$\sigma_{xx}=\frac{\partial^2\Phi}{\partial y^2},\quad \sigma_{yy}=\frac{\partial^2\Phi}{\partial x^2},\quad \sigma_{xy}=-\frac{\partial^2\Phi}{\partial x\partial y}.$$

**推导**：
(1) 自动满足平衡方程：
$$\partial_x\sigma_{xx}+\partial_y\sigma_{xy}=\partial_x\partial_y^2\Phi-\partial_x\partial_y^2\Phi=0;\quad \partial_x\sigma_{xy}+\partial_y\sigma_{yy}=-\partial_x^2\partial_y\Phi+\partial_x^2\partial_y\Phi=0.$$
(2) 协调方程（无体力、常 $\nu$）：
$$\nabla^4\Phi=\frac{\partial^4\Phi}{\partial x^4}+2\frac{\partial^4\Phi}{\partial x^2\partial y^2}+\frac{\partial^4\Phi}{\partial y^4}=0.$$
**推导**：将应力代入协调
$$\frac{\partial^2\varepsilon_{xx}}{\partial y^2}+\frac{\partial^2\varepsilon_{yy}}{\partial x^2}=2\frac{\partial^2\varepsilon_{xy}}{\partial x\partial y},$$
用平面应力/应变本构 $\varepsilon_{xx}\propto(\sigma_{xx}-\nu\sigma_{yy})$ 等，代入后交叉项抵消，得双调和方程 $\nabla^4\Phi=0$；有体力势 $V$（$\mathbf f=-\nabla V$）时平面应力下右端为 $-(1-\nu)\nabla^2V$。

**物理/几何意义**：Airy 函数使平衡自动满足，问题归结为求双调和函数。

**极限/特例**：多项式 $\Phi$ 给出均匀应力/纯弯/纯剪解。

**典型应用**：梁的应力、孔口、缺口、裂纹尖端场。

#### 知识点3. 平面问题的极坐标解答

**公式**：极坐标 Airy 函数 $\Phi(r,\theta)$，应力
$$\sigma_{rr}=\frac1r\frac{\partial\Phi}{\partial r}+\frac1{r^2}\frac{\partial^2\Phi}{\partial\theta^2},\quad \sigma_{\theta\theta}=\frac{\partial^2\Phi}{\partial r^2},\quad \sigma_{r\theta}=-\frac{\partial}{\partial r}\left(\frac1r\frac{\partial\Phi}{\partial\theta}\right).$$
**推导**：由柱坐标平衡方程
$$\frac{\partial\sigma_{rr}}{\partial r}+\frac1r\frac{\partial\sigma_{r\theta}}{\partial\theta}+\frac{\sigma_{rr}-\sigma_{\theta\theta}}{r}=0,\quad \frac{\partial\sigma_{r\theta}}{\partial r}+\frac1r\frac{\partial\sigma_{\theta\theta}}{\partial\theta}+\frac{2\sigma_{r\theta}}{r}=0,$$
引入 $\Phi$ 使前两式恒成立（Morera）。双调和方程
$$\left(\frac{\partial^2}{\partial r^2}+\frac1r\frac{\partial}{\partial r}+\frac1{r^2}\frac{\partial^2}{\partial\theta^2}\right)^2\Phi=0.$$

**物理/几何意义**：适合圆孔、圆盘、楔、轴对称问题。

**极限/特例**：轴对称 $\Phi=\Phi(r)$ 时 $\nabla^4\Phi=(\frac{d^2}{dr^2}+\frac1r\frac d{dr})^2\Phi=0$，通解 $\Phi=A\ln r+Br^2\ln r+Cr^2+D$。

**典型应用**：厚壁圆筒、圆盘旋转、Kirsch 圆孔。

#### 知识点4. Kirsch圆孔应力集中解

**问题**：无限大板，远场单轴拉伸 $\sigma_{xx}=\sigma_0$，中心圆孔半径 $a$。

**公式**：Michell 解
$$\sigma_{rr}=\frac{\sigma_0}{2}\left(1-\frac{a^2}{r^2}\right)+\frac{\sigma_0}{2}\left(1-4\frac{a^2}{r^2}+3\frac{a^4}{r^4}\right)\cos2\theta,$$
$$\sigma_{\theta\theta}=\frac{\sigma_0}{2}\left(1+\frac{a^2}{r^2}\right)-\frac{\sigma_0}{2}\left(1+3\frac{a^4}{r^4}\right)\cos2\theta,$$
$$\sigma_{r\theta}=-\frac{\sigma_0}{2}\left(1+2\frac{a^2}{r^2}-3\frac{a^4}{r^4}\right)\sin2\theta.$$

**推导**：设 $\Phi=\Phi_0(r)+\Phi_2(r)\cos2\theta$，代入双调和，解 Euler 型常微分方程，由远场条件与孔边 $\sigma_{rr}(a)=\sigma_{r\theta}(a)=0$ 定常数。

**物理/几何意义**：孔边 $\theta=\pm\pi/2$ 处 $\sigma_{\theta\theta}=3\sigma_0$，应力集中系数 $K_t=3$。$r\gg a$ 恢复远场。

**极限/特例**：$a\to0$ 无孔；$r=a$ 最大环向应力 $3\sigma_0$。

**典型应用**：孔边疲劳、结构开口补强。

#### 知识点5. 复变函数方法（Kolosov–Muskhelishvili）

**定义**：用两个解析函数 $\varphi(z),\psi(z)$（$z=x+iy$）表示
$$\sigma_{xx}+\sigma_{yy}=2[\varphi'(z)+\overline{\varphi'(z)}]=4\,\mathrm{Re}\,\varphi'(z),$$
$$\sigma_{yy}-\sigma_{xx}+2i\sigma_{xy}=2[\bar z\varphi''(z)+\psi'(z)],$$
$$2\mu(u_x+iu_y)=\kappa\varphi(z)-z\overline{\varphi'(z)}-\overline{\psi(z)},$$
其中 $\kappa=3-4\nu$（平面应变），$\kappa=(3-\nu)/(1+\nu)$（平面应力）。

**推导纲要**：由双调和函数可表示为
$$\Phi(x,y)=\mathrm{Re}[\bar z\varphi(z)+\chi(z)],\ \chi'(z)=\psi(z).$$
代入应力函数与位移-应变关系并利用解析函数性质（Cauchy–Riemann）即得 Kolosov 公式。

**物理/几何意义**：复势将边值问题转化为解析函数边值问题，便于保角映射处理复杂边界。

**极限/特例**：孔口/裂纹/夹杂问题可求出封闭解。

**典型应用**：裂纹尖端 Williams/Westergaard 展开、夹杂 Eshelby、孔口。

#### 知识点6. 半逆解法与多项式应力函数

**定义**：假定应力分量（或 $\Phi$）的某种函数形式，代入方程与边界定系数。

**推导示例**：
(1) 纯弯梁：$\Phi=Ay^3$，
$$\sigma_{xx}=\partial_y^2\Phi=6Ay,\ \sigma_{yy}=0,\ \sigma_{xy}=0.$$
由弯矩 $M=\int_A\sigma_{xx}y\,dA=6A I$，得 $A=\frac{M}{6I}$，$\sigma_{xx}=\frac{My}{I}$。
(2) 纯剪：$\Phi=Bxy$，$\sigma_{xy}=-B$ 为常数，$\sigma_{xx}=\sigma_{yy}=0$。
(3) 悬臂梁受端部力：$\Phi=Axy^3+Bxy$ 等。

**物理/几何意义**：工程常用；Airy 函数的多项式对应均匀/线性/二次应力场。

**极限/特例**：高次多项式可满足复杂边界。

**典型应用**：梁、槽、齿根应力分析。

#### 知识点7. 复变函数中的边值问题与映射法

**公式**：保角映射 $z=\omega(\zeta)$，复势变换为
$$\varphi(z)=\varphi[\omega(\zeta)]=\Phi(\zeta),\quad \psi(z)=\Psi(\zeta).$$
应力边界条件化为
$$\varphi(z)+z\overline{\varphi'(z)}+\overline{\psi(z)}=i\int(\bar t_x+i\bar t_y)\,ds+C.$$

**推导**：由应力矢量与复势关系，结合面力平衡积分即得。对含孔区域 $\omega(\zeta)=R(\zeta+m/\zeta)$（椭圆孔），在单位圆 $\zeta=e^{i\theta}$ 上展开边界条件，比较 Laurent 级数系数求复势。

**物理/几何意义**：把复杂区域映为圆/上半平面，再解析求解。

**极限/特例**：$m=0$ 圆孔（Kirsch）；$m\to1$ 裂纹。

**典型应用**：孔口、裂纹、夹杂的解析解。

#### 知识点8. 平面问题通解与Michell结构

**陈述**：齐次双调和函数 $\Phi$ 的一般解可由两个解析函数 $\varphi(z),\chi(z)$ 表示为 $\Phi=\mathrm{Re}[\bar z\,\varphi(z)+\chi(z)]$，即前面复势形式的实值化。Michell 给出极坐标通解：
$$\sigma_{rr}=\frac1r\frac{\partial\Phi}{\partial r}+\frac1{r^2}\frac{\partial^2\Phi}{\partial\theta^2},\ \dots$$
一般解
$$\Phi=a_0\ln r+b_0r^2+c_0r^2\ln r+d_0+\sum_{n}\Phi_n(r)\cos n\theta+\cdots$$

**推导**：极坐标双调和方程分离变量 $\Phi=R(r)\Theta(\theta)$，径向方程
$$r^4R''''+2r^3R'''-(\cdots)R''+\cdots=0,$$
解为 $R=r^n,r^{2-n},r^{-n},r^{n+2}$ 等，组合得通解。

**物理/几何意义**：任何平面弹性解可由此通解叠加。

**极限/特例**：轴对称与谐波项分别对应不同载荷模式。

**典型应用**：圆盘、旋转盘、楔、孔群。

### 第8章 能量原理与变分法

#### 知识点1. 虚功原理（虚位移原理）

**陈述**：变形体平衡的充要条件是对任意虚位移 $\delta u_i$（满足位移边界条件），内力虚功等于外力虚功：
$$\int_B\sigma_{ij}\delta\varepsilon_{ij}\,dV=\int_B f_i\delta u_i\,dV+\int_{\partial B_t}\bar t_i\delta u_i\,dS.$$

**推导**：由平衡 $\sigma_{ij,j}+f_i=0$，乘以 $\delta u_i$ 并在 $B$ 积分：
$$\int_B(\sigma_{ij,j}+f_i)\delta u_i\,dV=0.$$
分部积分
$$\int_B\sigma_{ij,j}\delta u_i\,dV=\int_{\partial B}\sigma_{ij}n_j\delta u_i\,dS-\int_B\sigma_{ij}\delta u_{i,j}\,dV.$$
用 $\delta\varepsilon_{ij}=\frac12(\delta u_{i,j}+\delta u_{j,i})$ 及 $\sigma_{ij}$ 对称，$\sigma_{ij}\delta u_{i,j}=\sigma_{ij}\delta\varepsilon_{ij}$；在 $\partial B_u$ 上 $\delta u_i=0$，在 $\partial B_t$ 上 $\sigma_{ij}n_j=\bar t_i$，代入即得虚功方程。反之由任意 $\delta u_i$ 可推回平衡与力边界。$\blacksquare$

**物理/几何意义**：平衡的弱形式，是有限元与能量法的基石。

**极限/特例**：线弹性时左端 $\delta U$，右端 $\delta W$。

**典型应用**：FEM、Ritz、Galerkin、单位载荷法。

#### 知识点2. 最小势能原理

**定义**：总势能
$$\Pi=\underbrace{\frac12\int_B C_{ijkl}\varepsilon_{ij}\varepsilon_{kl}\,dV}_{U}-\underbrace{\left(\int_Bf_iu_i\,dV+\int_{\partial B_t}\bar t_iu_i\,dS\right)}_{W}.$$

**推导**：对 $\mathbf u$ 取一阶变分并令 $\delta\Pi=0$：
$$\delta\Pi=\int_B\sigma_{ij}\delta\varepsilon_{ij}\,dV-\int_Bf_i\delta u_i\,dV-\int_{\partial B_t}\bar t_i\delta u_i\,dS=0,$$
即虚功原理，故平衡。二阶变分
$$\delta^2\Pi=\int_B C_{ijkl}\delta\varepsilon_{ij}\delta\varepsilon_{kl}\,dV\ge0$$
（$C$ 正定），故势能取极小。$\therefore$ 平衡解使 $\Pi$ 取最小值。$\blacksquare$

**物理/几何意义**：真实位移使总势能最小（稳定平衡）。

**极限/特例**：纯力边界时刚体位移不改变 $\Pi$。

**典型应用**：近似解法、稳定性判据、接触。

#### 知识点3. 最小余能原理与Castigliano定理

**定义**：余能 $U^*=\int_B W^*\,dV-\int_{\partial B_u}\bar u_i t_i\,dS$，线弹性时 $W^*=\frac12S_{ijkl}\sigma_{ij}\sigma_{kl}=W$。在满足平衡与力边界的应力中，真实应力使余能取极小（或驻值）。

**推导**：对静力许可应力场取变分，令 $\delta U^*=0$，得位移协调条件。

**Castigliano 定理**：
$$\delta U=\sum_k F_k\,\delta\Delta_k,\qquad \Delta_k=\frac{\partial U}{\partial F_k}\quad(\text{线性弹性}).$$
广义形式：$\Delta_k=\int_B\varepsilon_{ij}\frac{\partial\sigma_{ij}}{\partial F_k}\,dV$。
第一定理：$\delta\Delta_k=\partial U/\partial F_k$；第二定理：$\Delta_k=\partial U/\partial F_k$（线性弹性）。
**推导**：$U=\frac12\int\sigma_{ij}\varepsilon_{ij}$，$\frac{\partial U}{\partial F_k}=\frac12\int(\frac{\partial\sigma_{ij}}{\partial F_k}\varepsilon_{ij}+\sigma_{ij}\frac{\partial\varepsilon_{ij}}{\partial F_k})=\int\varepsilon_{ij}\frac{\partial\sigma_{ij}}{\partial F_k}$（用对称性）。而 $\int\varepsilon_{ij}\frac{\partial\sigma_{ij}}{\partial F_k}=\int u_{i,j}\frac{\partial\sigma_{ij}}{\partial F_k}$，分部积分并用平衡与边界，得 $=\Delta_k$。$\blacksquare$

**物理/几何意义**：位移可由能量对力求导得到。

**极限/特例**：求超静定反力用最小余能。

**典型应用**：桁架、框架位移计算。

#### 知识点4. 功互等定理（Betti与Maxwell）

**Betti 定理**：同一弹性体两组载荷状态 $(1),(2)$：
$$\int_{\partial B}t_i^{(1)}u_i^{(2)}\,dS+\int_Bf_i^{(1)}u_i^{(2)}\,dV=\int_{\partial B}t_i^{(2)}u_i^{(1)}\,dS+\int_Bf_i^{(2)}u_i^{(1)}\,dV.$$

**证明**：由 $C_{ijkl}$ 对称及散度定理，$\int\sigma_{ij}^{(1)}\varepsilon_{ij}^{(2)}=\int\sigma_{ij}^{(2)}\varepsilon_{ij}^{(1)}$，展开两边积分即得。$\blacksquare$

**Maxwell 互等**：单位力引起的位移矩阵对称，$\delta_{ij}=\delta_{ji}$。

**物理/几何意义**：影响系数对称，是结构力学力法基础。

**极限/特例**：线性弹性成立；非线性不成立。

**典型应用**：影响线、柔度矩阵、位移法。

#### 知识点5. 变分法基础与Euler–Lagrange方程

**定义**：泛函 $I[y]=\int_{x_1}^{x_2}F(x,y,y')\,dx$ 的极值条件 $\delta I=0$。

**推导**：设 $y_\epsilon=y+\epsilon\eta$，$\eta(x_1)=\eta(x_2)=0$，
$$\frac{dI}{d\epsilon}\Big|_{\epsilon=0}=\int_{x_1}^{x_2}\left(\frac{\partial F}{\partial y}\eta+\frac{\partial F}{\partial y'}\eta'\right)dx=0.$$
分部积分 $\int\frac{\partial F}{\partial y'}\eta'\,dx=[\frac{\partial F}{\partial y'}\eta]-\int\frac{d}{dx}(\frac{\partial F}{\partial y'})\eta\,dx$，边界项为零，故
$$\int\left(\frac{\partial F}{\partial y}-\frac{d}{dx}\frac{\partial F}{\partial y'}\right)\eta\,dx=0.$$
由 $\eta$ 任意得 Euler–Lagrange 方程
$$\frac{\partial F}{\partial y}-\frac{d}{dx}\left(\frac{\partial F}{\partial y'}\right)=0.$$
多维/高阶：$\frac{\partial F}{\partial y}-\partial_i\frac{\partial F}{\partial y_{,i}}=0$，$\frac{\partial F}{\partial u}-\partial_x\frac{\partial F}{\partial u_{xx}}-\partial_y\frac{\partial F}{\partial u_{yy}}-\partial_x\frac{\partial F}{\partial u_{xy}}=0$。

**物理/几何意义**：极值/驻值条件给出控制微分方程。

**极限/特例**：悬链线、最速降线、最小势能。

**典型应用**：梁板方程、场论、最优控制。

#### 知识点6. Ritz法与Galerkin法

**Ritz 法**：取试探解 $u_i\approx\sum_m a_{im}\phi_{im}$（满足位移边界），代入泛函，令
$$\frac{\partial\Pi}{\partial a_{im}}=0\ \Rightarrow\ \mathbf K\mathbf a=\mathbf F,$$
求解系数。

**推导**：$\Pi(\mathbf a)=\frac12\mathbf a^{\mathrm T}\mathbf K\mathbf a-\mathbf F^{\mathrm T}\mathbf a$，$\nabla_{\mathbf a}\Pi=\mathbf K\mathbf a-\mathbf F=0$。

**Galerkin 法（加权残量）**：将近似解代入微分方程得残量 $R$，要求与试探函数正交，检验函数取 $\phi_i$：
$$\int_B R(\mathbf u_N)\phi_i\,dV=0.$$

**物理/几何意义**：把无限维问题投影到有限维子空间。

**极限/特例**：Ritz 用于自伴（保守）问题；Galerkin 更一般。

**典型应用**：有限元、谱方法、加权残量法。

#### 知识点7. 能量形式的稳定性判据

**定义**：平衡态 $\mathbf u_0$ 稳定的充要条件是总势能二阶变分非负：
$$\delta^2\Pi\ge0\ \forall\ \delta\mathbf u.$$

**推导**：设 $u=u_0+\epsilon\delta u$，展开
$$\Pi(u)=\Pi(u_0)+\epsilon\,\delta\Pi+\tfrac12\epsilon^2\delta^2\Pi+O(\epsilon^3).$$
平衡 $\delta\Pi=0$；稳定极小要求 $\delta^2\Pi\ge0$。失稳时存在 $\delta u$ 使 $\delta^2\Pi<0$，$\epsilon$ 方向能量下降。$\blacksquare$

**物理/几何意义**：稳定的能量判据。对线性化问题即最小特征值非负。

**极限/特例**：屈曲临界点为 $\delta^2\Pi=0$。

**典型应用**：Euler 压杆、板屈曲、分叉分析。

#### 知识点8. 单位载荷法与Mohr积分

**公式**：线弹性结构位移
$$\Delta=\int_B\sigma_{ij}\tilde\varepsilon_{ij}\,dV\ \xrightarrow{\text{梁}}\ \Delta=\int_0^l\frac{M\tilde M}{EI}dx+\int_0^l\frac{N\tilde N}{EA}dx+\int_0^l k\frac{Q\tilde Q}{GA}dx,$$
其中 $\tilde M,\tilde N,\tilde Q$ 为单位载荷引起的内力。

**推导**：由虚功原理，对真实位移取单位载荷虚力系：
$$1\cdot\Delta=\int\sigma_{ij}\tilde\varepsilon_{ij}\,dV.$$
对梁用 $\sigma=M y/I,\tilde\sigma=\tilde M y/I$，积分得 $\int M\tilde M/(EI)dx$。

**物理/几何意义**：求位移的通用方法。

**极限/特例**：桁架用 $\sum N\tilde N l/(EA)$。

**典型应用**：刚架位移、超静定、温度位移。

### 第9章 分析力学：拉格朗日、哈密顿与连续体推广

#### 知识点1. 广义坐标、约束与自由度

**定义**：用 $n$ 个广义坐标 $q_1,\dots,q_n$ 描述系统位形。完整约束 $f_\alpha(\mathbf r,t)=0$ 可消去；非完整约束含速度不可积。

**公式与推导**：笛卡尔坐标 $\mathbf r_\alpha=\mathbf r_\alpha(q_1,\dots,q_n,t)$，速度
$$\dot{\mathbf r}_\alpha=\sum_{k}\frac{\partial\mathbf r_\alpha}{\partial q_k}\dot q_k+\frac{\partial\mathbf r_\alpha}{\partial t}.$$
自由度 $=3N-s$（$s$ 为完整约束数）。

**物理/几何意义**：广义坐标把约束"内化"，减少独立变量。

**极限/特例**：定常约束时 $\partial\mathbf r/\partial t=0$。

**典型应用**：多体动力学、机构、连续体离散化。

#### 知识点2. 达朗贝尔原理与虚功

**陈述**：对理想约束，主动力与惯性力的虚功之和为零：
$$\sum_\alpha(\mathbf F_\alpha-m_\alpha\ddot{\mathbf r}_\alpha)\cdot\delta\mathbf r_\alpha=0.$$

**推导**：牛顿方程 $\mathbf F_\alpha-m_\alpha\ddot{\mathbf r}_\alpha=\mathbf R_\alpha$（约束力）。对理想约束，约束力虚功 $\sum\mathbf R_\alpha\cdot\delta\mathbf r_\alpha=0$，两边点乘 $\delta\mathbf r_\alpha$ 求和即得。$\blacksquare$

**物理/几何意义**：动力学问题的静力学化，是拉格朗日方程的出发点。

**极限/特例**：静力学 $\ddot{\mathbf r}=0$ 退回虚功原理。

**典型应用**：复杂约束系统的动力学。

#### 知识点3. 拉格朗日方程

**公式**：
$$\frac{d}{dt}\left(\frac{\partial L}{\partial\dot q_k}\right)-\frac{\partial L}{\partial q_k}=Q_k^{nc},\qquad L=T-V.$$

**推导**：由达朗贝尔原理，将 $\delta\mathbf r_\alpha=\sum_k\frac{\partial\mathbf r_\alpha}{\partial q_k}\delta q_k$ 代入：
$$\sum_k\left[\sum_\alpha\mathbf F_\alpha\cdot\frac{\partial\mathbf r_\alpha}{\partial q_k}-\sum_\alpha m_\alpha\ddot{\mathbf r}_\alpha\cdot\frac{\partial\mathbf r_\alpha}{\partial q_k}\right]\delta q_k=0.$$
定义广义力 $Q_k=\sum_\alpha\mathbf F_\alpha\cdot\partial\mathbf r_\alpha/\partial q_k$。对第二项，利用恒等式
$$\dot{\mathbf r}_\alpha\cdot\frac{\partial\mathbf r_\alpha}{\partial q_k}=\frac{d}{dt}\left(\mathbf r_\alpha\cdot\frac{\partial\mathbf r_\alpha}{\partial q_k}\right)-\mathbf r_\alpha\cdot\frac{d}{dt}\frac{\partial\mathbf r_\alpha}{\partial q_k},$$
以及 $\frac{\partial\dot{\mathbf r}_\alpha}{\partial\dot q_k}=\frac{\partial\mathbf r_\alpha}{\partial q_k}$、$\frac{d}{dt}\frac{\partial\mathbf r_\alpha}{\partial q_k}=\frac{\partial\dot{\mathbf r}_\alpha}{\partial q_k}$，得
$$\sum_\alpha m_\alpha\ddot{\mathbf r}_\alpha\cdot\frac{\partial\mathbf r_\alpha}{\partial q_k}=\frac{d}{dt}\frac{\partial T}{\partial\dot q_k}-\frac{\partial T}{\partial q_k}.$$
若力保守 $Q_k=-\partial V/\partial q_k$，则
$$\frac{d}{dt}\frac{\partial T}{\partial\dot q_k}-\frac{\partial T}{\partial q_k}+\frac{\partial V}{\partial q_k}=0\ \Rightarrow\ \frac{d}{dt}\frac{\partial L}{\partial\dot q_k}-\frac{\partial L}{\partial q_k}=0.$$
$\blacksquare$

**物理/几何意义**：能量形式的最小作用量原理的欧拉-拉格朗日方程。

**极限/特例**：$L=T-V$ 对保守系统；耗散函数 $R$ 给出 $Q_k^{nc}=-\partial R/\partial\dot q_k$。

**典型应用**：振动、刚体、机器人、连续体。

#### 知识点4. 哈密顿原理（最小作用量）

**陈述**：在给定初末时刻位形间，真实运动使作用量
$$S=\int_{t_1}^{t_2}L(q,\dot q,t)\,dt$$
取驻值，$\delta S=0$。

**推导**：对 $q_k(t)+\epsilon\eta_k(t)$，$\eta_k(t_1)=\eta_k(t_2)=0$：
$$\delta S=\int_{t_1}^{t_2}\left(\frac{\partial L}{\partial q_k}\eta_k+\frac{\partial L}{\partial\dot q_k}\dot\eta_k\right)dt=\int_{t_1}^{t_2}\left(\frac{\partial L}{\partial q_k}-\frac{d}{dt}\frac{\partial L}{\partial\dot q_k}\right)\eta_k\,dt+\left[\frac{\partial L}{\partial\dot q_k}\eta_k\right]_{t_1}^{t_2}=0,$$
边界项为零，由 $\eta_k$ 任意得拉格朗日方程。$\blacksquare$

**物理/几何意义**：自然界的"最经济"原理，统一经典力学与场论。

**极限/特例**：相对论、电磁场、连续介质均由作用量导出。

**典型应用**：场论、弹性动力学变分、数值积分。

#### 知识点5. 哈密顿正则方程与相空间

**定义**：广义动量 $p_k=\partial L/\partial\dot q_k$，哈密顿函数
$$H=\sum_kp_k\dot q_k-L.$$

**推导**：
$$\delta H=\sum_k\left(\dot q_k\delta p_k+p_k\delta\dot q_k-\frac{\partial L}{\partial q_k}\delta q_k-\frac{\partial L}{\partial\dot q_k}\delta\dot q_k\right).$$
由 $p_k=\partial L/\partial\dot q_k$，含 $\delta\dot q_k$ 的项抵消；又 $\partial L/\partial q_k=\dot p_k$，故
$$dH=\sum_k(\dot q_k\,dp_k-\dot p_k\,dq_k-\frac{\partial L}{\partial t}dt).$$
而 $H=H(q,p,t)$，$dH=\sum(\frac{\partial H}{\partial q_k}dq_k+\frac{\partial H}{\partial p_k}dp_k)+\frac{\partial H}{\partial t}dt$。比较得哈密顿正则方程
$$\dot q_k=\frac{\partial H}{\partial p_k},\qquad \dot p_k=-\frac{\partial H}{\partial q_k},\qquad \frac{\partial H}{\partial t}=-\frac{\partial L}{\partial t}.$$

**物理/几何意义**：相空间 $(q,p)$ 中的辛流；$H$ 通常为总能量。

**极限/特例**：保守定常系统 $H=E=$ 常数。

**典型应用**：统计力学、量子力学对应、辛积分。

#### 知识点6. 泊松括号与守恒律

**定义**：$\{f,g\}=\sum_k\left(\frac{\partial f}{\partial q_k}\frac{\partial g}{\partial p_k}-\frac{\partial f}{\partial p_k}\frac{\partial g}{\partial q_k}\right)$。

**推导**：任意力学量 $f(q,p,t)$ 的时间演化
$$\dot f=\sum_k\left(\frac{\partial f}{\partial q_k}\dot q_k+\frac{\partial f}{\partial p_k}\dot p_k\right)+\frac{\partial f}{\partial t}=\{f,H\}+\frac{\partial f}{\partial t}.$$
$\therefore$ $f$ 守恒 $\iff \{f,H\}+\partial_tf=0$。若 $f$ 不显含 $t$，守恒 $\iff\{f,H\}=0$。

**物理/几何意义**：守恒量对应与哈密顿量的对合性；Noether 定理的相空间表述。

**极限/特例**：$\{q_i,p_j\}=\delta_{ij}$，$\{q_i,q_j\}=\{p_i,p_j\}=0$。

**典型应用**：可积系统、守恒量、量子化。

#### 知识点7. 诺特定理（对称与守恒）

**陈述**：作用量在连续对称变换下不变，则存在守恒量。

**推导（时间平移）**：若 $L$ 不显含 $t$，则
$$\frac{dL}{dt}=\sum_k\left(\frac{\partial L}{\partial q_k}\dot q_k+\frac{\partial L}{\partial\dot q_k}\ddot q_k\right)=\sum_k\left(\dot p_k\dot q_k+p_k\ddot q_k\right)=\frac{d}{dt}\sum_kp_k\dot q_k,$$
$\therefore\frac{d}{dt}(\sum p_k\dot q_k-L)=0$，即 $H$ 守恒。
（空间平移）$L$ 在 $q_k\to q_k+\epsilon$ 下不变 $\Rightarrow\sum_k\partial L/\partial q_k=0\Rightarrow\sum\dot p_k=0$，动量守恒。
（转动）角动量守恒。$\blacksquare$

**物理/几何意义**：对称性 $\leftrightarrow$ 守恒律（Noether）。

**极限/特例**：场论中能量-动量张量守恒 $T^{\mu\nu}_{,\nu}=0$。

**典型应用**：连续介质守恒律、粒子物理。

#### 知识点8. 连续体的拉格朗日-哈密顿描述

**定义**：位移场 $u_i(\mathbf X,t)$，拉格朗日密度
$$\mathcal L=\tfrac12\rho_0\dot u_i\dot u_i-W(\nabla\mathbf u)+\rho_0f_iu_i.$$

**推导**：作用量 $S=\int\mathcal L\,dV_0\,dt$。Euler–Lagrange：
$$\frac{\partial\mathcal L}{\partial u_i}-\partial_t\frac{\partial\mathcal L}{\partial\dot u_i}-\partial_J\frac{\partial\mathcal L}{\partial u_{i,J}}=0.$$
代入得
$$\rho_0f_i+\frac{\partial P_{iJ}}{\partial X_J}-\rho_0\ddot u_i=0,\qquad P_{iJ}=-\frac{\partial W}{\partial F_{iJ}},$$
即动量方程（参考构型）。正则动量 $\pi_i=\partial\mathcal L/\partial\dot u_i=\rho_0\dot u_i$，哈密顿密度 $\mathcal H=\pi_i\dot u_i-\mathcal L$。

**物理/几何意义**：连续体是无穷自由度拉格朗日系统，场方程由变分导出。

**极限/特例**：弹性动力学、流体变分（Clebsch 表示）。

**典型应用**：非线性弹性波、变分有限元、多尺度。

### 第10章 梁、板、壳的弯曲理论

#### 知识点1. Euler–Bernoulli梁

**假设**：直梁小变形；平截面假设（变形后仍垂直中性轴）；忽略剪切变形。

**公式与推导**：
(1) 几何：曲率 $\kappa\approx w''$，应变 $\varepsilon_x=-z w''$（$z$ 从中性轴起）。
(2) 本构与内力：
$$M=\int_A\sigma_xz\,dA=\int_A E\varepsilon_xz\,dA=-Ew''\int_Az^2\,dA=-EIw'',$$
其中 $I=\int_A z^2dA$。
(3) 平衡：微段 $\frac{dQ}{dx}=-q$，$\frac{dM}{dx}=Q$，合得
$$\frac{d^2}{dx^2}\left(EI\frac{d^2w}{dx^2}\right)=q(x).$$
**推导**：对微段竖向力平衡 $Q-(Q+dQ)-q\,dx=0\Rightarrow dQ/dx=-q$；力矩平衡 $dM/dx=Q$；代入 $M=-EIw''$。

**物理/几何意义**：弯曲主导，$\sigma_x=-Mz/I$。

**极限/特例**：简支梁均布载荷 $w=\frac{q x(l^3-2lx^2+x^3)}{24EI}$。

**典型应用**：细长梁、桁架、建筑结构。

#### 知识点2. Timoshenko梁（考虑剪切）

**假设**：横截面转动 $\phi$ 独立于挠度斜率，剪应变 $\gamma=\phi-w'$。

**公式与推导**：
$$M=EI\phi',\qquad Q=\kappa GA(\phi-w'),$$
平衡
$$\frac{d}{dx}\left[\kappa GA(w'-\phi)\right]+q=0,\qquad \frac{d}{dx}(EI\phi')+\kappa GA(w'-\phi)=0.$$
**推导**：位移场 $u_x=-z\phi(x),u_z=w(x)$，$\varepsilon_x=-z\phi'$，$\gamma_{xz}=w'-\phi$。内力由本构积分，微段平衡给出。消去 $\phi$：
$$EIw''''=q-\frac{EI}{\kappa GA}q''.$$

**物理/几何意义**：适合中厚/短梁，剪切变形不可忽略。

**极限/特例**：$\kappa GA\to\infty$ 退回 Euler–Bernoulli。

**典型应用**：复合材料梁、高频振动、层合梁。

#### 知识点3. Kirchhoff薄板

**假设**：中面法线保持直线且垂直中面（直法线假设），忽略横向剪切与法向应力。

**公式与推导**：
(1) 位移 $u=-z\partial_xw,\ v=-z\partial_yw,\ w=w(x,y)$。
(2) 应变 $\varepsilon_x=-z w_{,xx},\ \varepsilon_y=-z w_{,yy},\ \gamma_{xy}=-2z w_{,xy}$。
(3) 弯矩/扭矩
$$M_x=-D(w_{,xx}+\nu w_{,yy}),\quad M_y=-D(w_{,yy}+\nu w_{,xx}),\quad M_{xy}=-D(1-\nu)w_{,xy},$$
其中弯曲刚度 $D=\frac{Eh^3}{12(1-\nu^2)}$。
(4) 平衡（Love–Kirchhoff）
$$D\nabla^4w=q(x,y).$$
**推导**：由 $\sigma$ 沿厚度积分成弯矩，代入板微元力矩/横向力平衡：
$$\frac{\partial^2M_x}{\partial x^2}+2\frac{\partial^2M_{xy}}{\partial x\partial y}+\frac{\partial^2M_y}{\partial y^2}=-q,$$
代入弯矩表达式即得 $D\nabla^4w=q$。

**物理/几何意义**：薄板以弯曲为主，挠度满足双调和方程。

**极限/特例**：四边简支矩形板 Navier 解（双三角级数）；轴对称圆板 $D\nabla^4w=q$，$\nabla^2=\frac1r\frac d{dr}(r\frac d{dr})$。

**典型应用**：机翼蒙皮、压力容器封头、建筑楼板。

#### 知识点4. Mindlin–Reissner中厚板

**假设**：中面法线保持直线但不再垂直中面，转动 $\phi_x,\phi_y$ 独立。

**公式与推导**：
$$M_x=D(\phi_{x,x}+\nu\phi_{y,y}),\quad M_{xy}=\frac{D(1-\nu)}{2}(\phi_{x,y}+\phi_{y,x}),$$
$$Q_x=\kappa Gh(\phi_x+w_{,x}),\quad Q_y=\kappa Gh(\phi_y+w_{,y}).$$
平衡
$$Q_{x,x}+Q_{y,y}+q=0,\qquad M_{x,x}+M_{xy,y}-Q_x=0,\qquad M_{xy,x}+M_{y,y}-Q_y=0.$$
**推导**：位移 $u=z\phi_x,v=z\phi_y$，剪应变 $\gamma_{xz}=\phi_x+w_{,x}$ 等；内力积分与微元平衡。

**物理/几何意义**：考虑横向剪切，适合厚板/复合材料；需剪切修正系数 $\kappa$（常取 $5/6$）。

**极限/特例**：$h\to0$ 时需注意剪切闭锁。

**典型应用**：厚板、层合板、Reissner–Mindlin 有限元。

#### 知识点5. 薄壳薄膜理论

**定义**：壳中面为曲面，薄膜状态忽略弯矩，只考虑膜力 $N_{ij}$。

**公式**：球壳内压 $p$、半径 $R$、壁厚 $t$：
$$\sigma=\frac{pR}{2t}\quad(\text{环向/经向}).$$
圆柱壳内压：环向 $\sigma_\theta=\frac{pR}{t}$，轴向 $\sigma_z=\frac{pR}{2t}$。
**推导**：取半球/半圆柱，压力合力与膜力平衡。球壳：$p\pi R^2=N\cdot2\pi R$，$N=\sigma t$，$\sigma=pR/(2t)$。圆柱环向：$p\,2R\,L=2\sigma_\theta tL\Rightarrow\sigma_\theta=pR/t$。

**物理/几何意义**：曲率使膜力承担横向载荷。

**极限/特例**：薄壳失稳临界压 $p_{cr}$ 由稳定性理论给出。

**典型应用**：压力容器、储罐、机身蒙皮。

#### 知识点6. 变分推导梁板方程

**推导（梁）**：势能
$$\Pi=\int_0^l\left(\frac12EI(w'')^2-qw\right)dx.$$
变分
$$\delta\Pi=\int_0^l(EIw''\delta w''-q\delta w)dx=0.$$
两次分部积分：
$$\int EIw''\delta w''\,dx=[EIw''\delta w']-\int(EIw'')'\delta w'\,dx=[EIw''\delta w'- (EIw'')'\delta w]+\int(EIw'')''\delta w\,dx.$$
边界项给出自然边界条件（弯矩、剪力），内部得 $(EIw'')''=q$。$\blacksquare$

**推导（板）**：
$$\Pi=\iint\left[\frac{D}{2}(w_{,xx}+w_{,yy})^2-(\cdots)\right]dxdy.$$
变分并分部积分得 $D\nabla^4w=q$，自然边界为 Kirchhoff 剪力 $V_n=Q_n+\partial_sM_{ns}$。

**物理/几何意义**：能量泛函给出控制方程与自然边界条件。

**极限/特例**：Ritz 近似可求近似解。

**典型应用**：有限元板单元、稳定性分析。

#### 知识点7. 梁的剪切中心与约束扭转

**定义**：开口薄壁梁扭转含自由扭转 $GJ$ 与约束扭转（翘曲受约束）$EC_w$。

**公式**：
$$M_t=GJ\theta'-EC_w\theta''',\qquad EC_w\theta''''-GJ\theta''=m.$$
**推导**：扭转角 $\theta$，扇性坐标 $\omega$，翘曲位移 $w_s=-\omega\theta'$，双力矩 $B=-EC_w\theta''$，由微元扭矩平衡 $B'=M_t-GJ\theta'$ 消去即得方程。

**物理/几何意义**：开口截面翘曲大，约束扭转产生正应力。

**极限/特例**：自由扭转型（长梁）退化为 $M_t=GJ\theta'$。

**典型应用**：工字梁、薄壁结构。

#### 知识点8. 梁柱效应与轴向-横向耦合

**定义**：轴向压力 $N$ 下梁的横向平衡含几何非线性项。

**公式**：
$$EIw''''+Nw''=q.$$
**推导**：微段力矩平衡多出 $N w$ 项：$\frac{d^2M}{dx^2}+N\frac{d^2w}{dx^2}=q$，代入 $M=-EIw''$。

**物理/几何意义**：压力软化结构，导致失稳。

**极限/特例**：$q=0$ 时本征值问题给临界压力（第11章）。

**典型应用**：压杆、板、框架稳定性。

### 第11章 稳定性理论与屈曲

#### 知识点1. 稳定性概念（Lyapunov与能量）

**定义**：平衡态 $\mathbf u_0$ 稳定：任给 $\epsilon>0$，存在 $\delta>0$，当初始扰动 $\|\mathbf u(0)-\mathbf u_0\|<\delta$ 时，$\|\mathbf u(t)-\mathbf u_0\|<\epsilon$ 对所有 $t>0$。渐近稳定再加 $\to0$。

**能量判据**：保守系统平衡稳定 $\iff$ 势能取严格局部极小，即 $\delta^2\Pi\ge0$。

**推导**：由 $\Pi=\Pi_0+\frac12\delta^2\Pi+O(\epsilon^3)$，若 $\delta^2\Pi>0$，扰动有界（能量守恒）；若存在 $\delta^2\Pi<0$ 的方向，则该方向能量无下界而失稳。$\blacksquare$

**物理/几何意义**：势能曲面极小点为稳定平衡。

**极限/特例**：中性平衡 $\delta^2\Pi=0$（临界）。

**典型应用**：屈曲、相变、结构失稳。

#### 知识点2. 分叉、极限点与缺陷敏感性

**定义**：载荷-位移路径在临界点分叉（pitchfork/transcritical）为分叉屈曲；在极限点（snap-through）为极值屈曲。

**推导（单自由度模型）**：势能 $\Pi=\frac12k u^2-\lambda f(u)$，平衡 $\Pi'=0$。临界条件 $\Pi''=0$。对对称分叉，展开 $\Pi=\frac12a u^2+\frac14bu^4-\lambda cu^2$，平衡 $u(a-2\lambda c+bu^2)=0$，分叉点 $\lambda_c=a/(2c)$。

**物理/几何意义**：临界点后出现多个平衡解，结构可能跳跃。

**极限/特例**：初始缺陷使分叉变为极值点，临界载荷下降（Koiter）。

**典型应用**：薄壳、拱、压力容器。

#### 知识点3. Euler压杆临界力

**问题**：两端铰支压杆，长度 $l$，弯曲刚度 $EI$，轴压 $P$。

**公式**：
$$P_{cr}=\frac{\pi^2EI}{l^2}.$$
**推导**：微弯后平衡 $EIw''''+Pw''=0$。令 $\lambda^2=P/(EI)$：
$$w''''+\lambda^2w''=0.$$
通解 $w=A\sin\lambda x+B\cos\lambda x+Cx+D$。边界 $w(0)=w''(0)=0$，$w(l)=w''(l)=0$。由 $w''=0$：$B=0,\ D=0$；$w(l)=A\sin\lambda l=0$。非平凡解要求
$$\sin\lambda l=0\Rightarrow\lambda l=n\pi\Rightarrow P_n=\frac{n^2\pi^2EI}{l^2}.$$
最小 $n=1$ 给 $P_{cr}=\pi^2EI/l^2$。$\blacksquare$

**物理/几何意义**：临界力与 $EI$ 成正比、与 $l^2$ 成反比；失稳模态为半正弦。

**极限/特例**：其他支承用有效长度 $\mu l$，$P_{cr}=\pi^2EI/(\mu l)^2$：固-自由 $\mu=2$，固-铰 $\mu\approx0.7$，固-固 $\mu=0.5$。

**典型应用**：柱、桁架压杆、细长构件。

#### 知识点4. 梁柱方程与放大系数

**公式**：横向载荷+轴压 $P$：
$$EIw''''+Pw''=q.$$
简支梁均布 $q$ 解
$$w_{\max}\approx\frac{5ql^4}{384EI}\cdot\frac{1}{1-P/P_{cr}},\qquad k=\sqrt{P/EI},\ P_{cr}=\frac{\pi^2EI}{l^2}.$$
**推导**：令 $w_p$ 为特解，通解加齐次项，边界定常数。轴向力使挠度放大，放大系数
$$AF=\frac{1}{1-P/P_{cr}}.$$

**物理/几何意义**：轴压放大横向变形（P-$\Delta$ 效应）。

**极限/特例**：$P\to0$ 退回普通弯曲；$P\to P_{cr}$ 挠度发散。

**典型应用**：框架二阶效应、稳定校核。

#### 知识点5. 板的屈曲

**公式**：四边简支矩形板（$a\times b$）面内双轴或单轴受压，临界应力
$$\sigma_{cr}=k\frac{\pi^2D}{b^2h},\qquad D=\frac{Eh^3}{12(1-\nu^2)}.$$
单轴简支 $k=\left(\frac{mb}{a}+\frac{a}{mb}\right)^2$，最小 $k=4$（$a/b$ 为整数比）。
**推导**：设 $w=\sum A_{mn}\sin\frac{m\pi x}{a}\sin\frac{n\pi y}{b}$，代入平衡
$$D\nabla^4w+N_xw_{,xx}+N_yw_{,yy}=0,$$
得 $N_x$ 特征值，取最小。

**物理/几何意义**：板屈曲形成半波模态；临界应力随宽厚比平方下降。

**极限/特例**：长板 $k\to4$；加劲板提高临界值。

**典型应用**：机翼壁板、船体、加筋板。

#### 知识点6. 能量法求临界载荷（Rayleigh–Ritz）

**方法**：设满足边界条件的挠曲模态 $w=a\phi(x)$，令 $\delta^2\Pi=0$ 即弯曲应变能等于外力做功：
$$U_b=\frac12\int EI(w'')^2dx,\qquad W=\frac12P\int(w')^2dx.$$
临界条件 $U_b=W(P_{cr})$：
$$P_{cr}=\frac{\int EI(w'')^2dx}{\int(w')^2dx}.$$
**推导**：由 $EIw''''+Pw''=0$ 乘 $w$ 积分：$\int EI(w'')^2=\int P(w')^2$（分部积分，边界项为零），故上式。对给定模态得近似临界力，且 Ritz 给上界。

**物理/几何意义**：弯曲刚度与几何刚度的竞争。

**极限/特例**：取真实模态得精确值；取近似模态得偏大值（保守）。

**典型应用**：复杂边界压杆/板、加筋结构。

#### 知识点7. 动力稳定性与参激振动

**方程（Mathieu）**：
$$\ddot q+(\delta+2\epsilon\cos t)q=0.$$
轴向周期载荷下的梁：
$$EIw''''+(P_0+P_1\cos\Omega t)w''+\rho A\ddot w=0.$$
**推导**：由梁柱方程加惯性项，设 $w=\phi(x)q(t)$ 分离变量，得 Mathieu–Hill 方程。不稳定区（共振舌）由 Floquet 理论确定。

**物理/几何意义**：周期轴力可引发参激共振，频率为固有频率 $2/n$ 附近。

**极限/特例**：$\epsilon=0$ 退化为简谐振动；无阻尼时无穷多舌。

**典型应用**：直升机旋翼、管道、周期载荷结构。

#### 知识点8. 后屈曲与Koiter渐近

**公式**：对缺陷 $\bar\lambda$ 敏感结构，Koiter 展开
$$\lambda_c/\lambda_{cl}=1-C(\bar\lambda/\rho)^{2/3}+\cdots$$
（缺陷敏感，指数 $2/3$）。对分叉结构，后屈曲路径
$$\lambda=\lambda_{cl}+a\,\xi^2+\cdots,\quad \xi\text{ 为模态幅值}.$$
**推导**：将势能在临界模态方向展开到高阶，加入缺陷项，求平衡路径的渐近解。$a>0$ 为稳定后屈曲，$a<0$ 为不稳定。

**物理/几何意义**：缺陷使临界载荷显著降低（薄壳尤为敏感，可达 $1/5$）。

**极限/特例**：完善结构 $C=0$ 时分叉。

**典型应用**：薄壳屈曲设计、Koiter 理论。

### 第12章 弹性动力学与波传播

#### 知识点1. 弹性动力学运动方程

**公式**：各向同性无体力
$$\rho\ddot u_i=(\lambda+\mu)u_{j,ji}+\mu u_{i,jj}.$$

**推导**：由动量方程 $\sigma_{ji,j}=\rho\ddot u_i$ 与各向同性本构即得，见第6章。

**物理/几何意义**：弹性波由惯性项与刚度项平衡控制。

**极限/特例**：静力 $\ddot u=0$ 退回平衡。

**典型应用**：地震波、超声检测、冲击。

#### 知识点2. Helmholtz分解与波方程解耦

**定义**：$\mathbf u=\nabla\phi+\nabla\times\boldsymbol\psi$，规范 $\nabla\cdot\boldsymbol\psi=0$。

**推导**：代入运动方程
$$\rho(\nabla\ddot\phi+\nabla\times\ddot{\boldsymbol\psi})=(\lambda+2\mu)\nabla(\nabla^2\phi)+\mu\nabla\times(\nabla^2\boldsymbol\psi).$$
使用 $\nabla(\nabla\cdot\mathbf u)=\nabla(\nabla^2\phi)$（因 $\nabla\cdot(\nabla\times\boldsymbol\psi)=0$）与 $\nabla^2\mathbf u=\nabla(\nabla^2\phi)+\nabla\times(\nabla^2\boldsymbol\psi)$，分离得
$$\nabla^2\phi=\frac{1}{c_p^2}\ddot\phi,\qquad \nabla^2\boldsymbol\psi=\frac{1}{c_s^2}\ddot{\boldsymbol\psi},$$
$$c_p=\sqrt{\frac{\lambda+2\mu}{\rho}},\qquad c_s=\sqrt{\frac{\mu}{\rho}}.$$

**物理/几何意义**：$c_p$ 为膨胀（纵）波速度，$c_s$ 为畸变（横）波速度；$c_p>c_s$。

**极限/特例**：不可压 $\lambda\to\infty$，$c_p\to\infty$；$\mu=0$ 无横波。

**典型应用**：地震 P/S 波、无损检测。

#### 知识点3. 一维杆波与梁波

**公式**：
(1) 杆纵波：$\rho\ddot u=E u''$，$c=\sqrt{E/\rho}$。
**推导**：一维本构 $\sigma=E u'$，运动 $\sigma'=\rho\ddot u$。
(2) Euler 梁弯曲波：$w''''+\frac{\rho A}{EI}\ddot w=0$，简谐解 $w=e^{i(kx-\omega t)}$，色散
$$\omega=\sqrt{\frac{EI}{\rho A}}k^2\ (\text{弥散}).$$
(3) Timoshenko 梁波含两支（弯曲/剪切）。

**物理/几何意义**：杆波非弥散，梁波弥散（相速度随频率变）。

**极限/特例**：细长梁低阶。

**典型应用**：冲击杆（Hopkinson）、梁振动、超声。

#### 知识点4. 平面波与Christoffel方程

**定义**：平面波 $\mathbf u=\mathbf A\,e^{i(\mathbf k\cdot\mathbf x-\omega t)}$。

**推导**：代入各向同性运动方程：
$$-\rho\omega^2A_i+(\lambda+\mu)(\mathbf k\cdot\mathbf A)k_i+\mu k^2A_i=0.$$
取 $\mathbf k=k\mathbf n$，设纵波 $\mathbf A\parallel\mathbf n$：
$$[(\lambda+2\mu)k^2-\rho\omega^2]A=0\Rightarrow c_p=\omega/k=\sqrt{(\lambda+2\mu)/\rho}.$$
横波 $\mathbf A\perp\mathbf n$：
$$[\mu k^2-\rho\omega^2]A=0\Rightarrow c_s=\sqrt{\mu/\rho}.$$
一般各向异性得 Christoffel 方程
$$|C_{ijkl}n_jn_l-\rho c^2\delta_{ik}|=0.$$

**物理/几何意义**：波速由弹性常数与密度确定；横波两偏振。

**极限/特例**：各向同性横波两种偏振简并。

**典型应用**：晶体声学、地壳波速反演。

#### 知识点5. 界面的反射与折射

**条件**：位移连续、面力连续（$\llbracket u_i\rrbracket=0,\ \llbracket\sigma_{ij}n_j\rrbracket=0$）。

**Snell 定律**：$\frac{\sin\theta_i}{c_1}=\frac{\sin\theta_r}{c_1}=\frac{\sin\theta_t}{c_2}$。
**推导**：波沿界面相速度必须匹配，即 $k_x$ 连续。

**反射系数**：垂直入射时（P 波）
$$R=\frac{\rho_2c_2-\rho_1c_1}{\rho_2c_2+\rho_1c_1},\qquad T=\frac{2\rho_2c_2}{\rho_2c_2+\rho_1c_1},$$
其中 $\rho c$ 为声阻抗。
**推导**：由位移连续 $A_i+A_r=A_t$ 与应力连续 $\rho_1c_1(A_i-A_r)=\rho_2c_2A_t$ 联立。

**物理/几何意义**：阻抗失配决定反射；斜入射有波型转换。

**极限/特例**：自由面 $R=-1$；匹配阻抗 $R=0$。

**典型应用**：无损检测、地震勘探、超声换能器。

#### 知识点6. Rayleigh表面波

**定义**：半空间自由表面附近的表面波，位移随深度指数衰减。

**公式**：波速 $c_R$ 满足
$$\left(2-\frac{c_R^2}{c_s^2}\right)^2=4\sqrt{1-\frac{c_R^2}{c_p^2}}\sqrt{1-\frac{c_R^2}{c_s^2}}.$$
近似 $c_R\approx\frac{0.87+1.12\nu}{1+\nu}c_s$。
**推导**：取 $\phi=Ae^{-qz}e^{i(kx-\omega t)}$、$\psi=Be^{-sz}e^{i(kx-\omega t)}$，由自由面 $\sigma_{zz}=\sigma_{xz}=0$ 得特征方程，解出 $c_R$。

**物理/几何意义**：能量集中于表面，振幅随深度衰减。

**极限/特例**：$c_R<c_s<c_p$；二维地震波主要成分。

**典型应用**：地震破坏、表面波无损检测、SAW 器件。

#### 知识点7. 梁与板的自由振动

**公式**：简支梁固有频率
$$\omega_n=(n\pi)^2\sqrt{\frac{EI}{\rho A l^4}},\qquad n=1,2,\dots$$
**推导**：$w=\phi(x)e^{i\omega t}$，$EI\phi''''-\rho A\omega^2\phi=0$，边界条件给特征值。$\phi=\sin\frac{n\pi x}{l}$，代入得 $\omega_n$。
板（四边简支）：
$$\omega_{mn}=\pi^2\left(\frac{m^2}{a^2}+\frac{n^2}{b^2}\right)\sqrt{\frac{D}{\rho h}}.$$

**物理/几何意义**：离散模态与频率；模态正交。

**极限/特例**：瑞利法 $\omega^2=U_b/W_k$。

**典型应用**：振动控制、模态分析、共振避让。

#### 知识点8. 波的弥散与群速度

**定义**：相速度 $c_p=\omega/k$，群速度 $c_g=d\omega/dk$。

**推导**：两相近波叠加
$$u=e^{i(kx-\omega t)}+e^{i((k+dk)x-(\omega+d\omega)t)}=2\cos\left(\frac{dk\,x-d\omega\,t}{2}\right)e^{i(kx-\omega t)},$$
包络速度即 $c_g=d\omega/dk$；能量传播速度等于群速度。

**物理/几何意义**：弥散介质中不同频率波速不同，波形展宽。

**极限/特例**：非弥散介质 $c_g=c_p$；梁波 $c_g=2c_p$（Euler 梁）。

**典型应用**：波包、地震波、光纤、声子晶体。

### 第13章 塑性力学基础

#### 知识点1. 屈服准则：Tresca与von Mises

**定义**：屈服准则给出材料由弹性进入塑性的应力条件 $f(\sigma_{ij})=0$。

**Tresca（最大剪应力）**：
$$\tau_{\max}=\frac{\sigma_1-\sigma_3}{2}=k\ \Rightarrow\ \sigma_1-\sigma_3=\sigma_Y.$$
**von Mises（J2）**：
$$J_2=k^2,\qquad \sigma_{\mathrm{eq}}=\sqrt{3J_2}=\sigma_Y.$$
单轴下 $\sigma_{\mathrm{eq}}=\sigma_1$，故 $\sigma_Y$ 为单轴屈服应力。纯剪切下 Tresca 给 $k=\sigma_Y/2$，Mises 给 $k=\sigma_Y/\sqrt3$。

**推导**：Tresca 由最大剪应力达到临界；Mises 由畸变能达到临界，$J_2=\frac16[(\sigma_1-\sigma_2)^2+\cdots]$。在 $\pi$ 平面（$\sigma_1+\sigma_2+\sigma_3=0$）上，Tresca 为正六边形，Mises 为其外接圆。

**物理/几何意义**：静水压不影响金属屈服（压力不敏感性）；两准则相差约 $15\%$。

**极限/特例**：$\nu=1/2$ 时 Mises 屈服等价于畸变能密度达限。

**典型应用**：金属成形、压力容器、疲劳。

#### 知识点2. 加载条件、硬化与一致性条件

**定义**：屈服面 $f(\sigma_{ij},q)=0$，$q$ 为硬化参量（等效塑性应变 $\bar\varepsilon^p$ 或塑性功 $W^p$）。加载/中性/卸载判据
$$f=0,\ \frac{\partial f}{\partial\sigma_{ij}}\dot\sigma_{ij}\begin{cases}>0&\text{加载}\\=0&\text{中性}\\<0&\text{卸载}\end{cases}$$

**一致性条件**：
$$\dot f=\frac{\partial f}{\partial\sigma_{ij}}\dot\sigma_{ij}+\frac{\partial f}{\partial q}\dot q=0\quad(\text{塑性加载}).$$

**推导**：塑性加载时应力点保持在屈服面上，故 $\dot f=0$。对等向硬化 $q=\bar\varepsilon^p$，$\dot{\bar\varepsilon}^p=\sqrt{\frac23\dot\varepsilon^p_{ij}\dot\varepsilon^p_{ij}}$。

**物理/几何意义**：屈服面随塑性变形扩展（等向）或移动（随动）。

**极限/特例**：理想塑性 $\partial f/\partial q=0$。

**典型应用**：金属循环塑性、Bauschinger 效应。

#### 知识点3. 塑性流动理论（Prandtl–Reuss）

**公式**：总应变率 $\dot\varepsilon_{ij}=\dot\varepsilon^e_{ij}+\dot\varepsilon^p_{ij}$，
$$\dot\varepsilon^e_{ij}=\frac{1+\nu}{E}\dot\sigma_{ij}-\frac{\nu}{E}\dot\sigma_{kk}\delta_{ij},$$
$$\dot\varepsilon^p_{ij}=\dot\lambda\frac{\partial f}{\partial\sigma_{ij}}=\dot\lambda\,s_{ij}\quad(\text{Mises }f=J_2).$$

**推导**：关联流动法则 $\dot\varepsilon^p_{ij}=\dot\lambda\partial f/\partial\sigma_{ij}$（正交流动）。对 Mises $f=\frac12s_{ij}s_{ij}-k^2$，$\partial f/\partial\sigma_{ij}=s_{ij}$，故 $\dot\varepsilon^p_{ij}=\dot\lambda s_{ij}$。$\dot\lambda$ 由一致性条件确定：
$$\dot\lambda=\frac{s_{ij}\dot s_{ij}}{2k^2H'}\ \text{等}.$$
具体地对 $f=J_2-k^2(\bar\varepsilon^p)$，$s_{ij}\dot s_{ij}$ 由弹性关系表示，$\dot\lambda=\frac{s_{ij}\dot\sigma_{ij}}{s_{ij}s_{ij}H'}$（$H'$ 为塑性模量）。

**物理/几何意义**：塑性应变率垂直于屈服面（法向），体积不可压（$s_{ii}=0$）。

**极限/特例**：理想塑性 $H'=0$ 得 Lévy–Mises $\dot\varepsilon^p_{ij}=\dot\lambda s_{ij}$。

**典型应用**：金属成形、弹塑性有限元。

#### 知识点4. 增量弹塑性本构与算法

**公式**：弹塑性切线刚度
$$\dot\sigma_{ij}=C^{ep}_{ijkl}\dot\varepsilon_{kl},$$
$$C^{ep}_{ijkl}=C_{ijkl}-\frac{C_{ijmn}\frac{\partial f}{\partial\sigma_{mn}}\frac{\partial f}{\partial\sigma_{pq}}C_{pqkl}}{\frac{\partial f}{\partial\sigma_{rs}}C_{rsuv}\frac{\partial f}{\partial\sigma_{uv}}+H'}.$$

**推导**：由 $\dot\sigma=C(\dot\varepsilon-\dot\varepsilon^p)$，$\dot\varepsilon^p=\dot\lambda\partial f/\partial\sigma$，一致性条件联立消去 $\dot\lambda$，即得。

**物理/几何意义**：塑性加载时刚度降低；卸载用弹性刚度。

**极限/特例**：理想塑性 $H'=0$。

**典型应用**：径向返回算法（return mapping）、ABAQUS/ANSYS 塑性。

#### 知识点5. 形变理论（Hencky）与全量关系

**公式**：
$$\varepsilon^p_{ij}=\frac{3}{2}\frac{\bar\varepsilon^p}{\sigma_{\mathrm{eq}}}s_{ij},$$
$$\bar\varepsilon^p=\frac{\sqrt2}{3}\sqrt{(\varepsilon^p_1-\varepsilon^p_2)^2+(\varepsilon^p_2-\varepsilon^p_3)^2+(\varepsilon^p_3-\varepsilon^p_1)^2}.$$
**推导**：假设塑性应变张量与应力偏量成比例（比例加载），由等价关系定义 $\bar\varepsilon^p$ 使单轴时 $\bar\varepsilon^p=\varepsilon^p$。

**物理/几何意义**：仅适用于简单（比例）加载，但形式简单。

**极限/特例**：卸载须用增量理论（否则卸载路径错误）。

**典型应用**：塑性理论解析解、平板塑性弯曲。

#### 知识点6. 极限分析：上限与下限定理

**下限定理**：任何满足平衡、力边界与 $f(\sigma)\le0$ 的静力许可应力场，其对应载荷 $\lambda^-$ 满足 $\lambda^-\le\lambda_p$。
**上限定理**：任何满足位移边界与流动法则的运动许可速度场，其耗散功率 $\int\sigma_Y|\dot\varepsilon^p|dV$ 与外功率之比 $\lambda^+\ge\lambda_p$。
$$\lambda_p=\max\lambda^-=\min\lambda^+.$$

**推导纲要**：由虚功率原理与塑性耗散不等式 $(\sigma_{ij}-\sigma^*_{ij})\dot\varepsilon^p_{ij}\ge0$（稳定材料）证明上下限夹逼。

**物理/几何意义**：真极限载荷夹在静力许可与运动许可之间。

**极限/特例**：理想塑性、小变形、载荷比例加载。

**典型应用**：压力容器爆破压力、梁板极限承载。

#### 知识点7. 弹塑性弯曲与极限弯矩

**公式**：矩形截面 $b\times h$，理想弹塑性（屈服 $\sigma_Y$）。
(1) 弹性极限弯矩 $M_e=\sigma_Y bh^2/6$。
(2) 塑性极限弯矩 $M_p=\sigma_Y bh^2/4$，形状系数 $\frac{M_p}{M_e}=1.5$。
**推导**：塑性极限时截面全塑性，应力分布 $\pm\sigma_Y$：
$$M_p=\int_{-h/2}^{h/2}\sigma_Y|z|b\,dz=2\sigma_Yb\frac{(h/2)^2}{2}=\frac{\sigma_Ybh^2}{4}.$$
弹塑性阶段中性轴位置随加载移动（轴对称时保持）。

**物理/几何意义**：形成塑性铰，结构可继续承载但变形增大。

**极限/特例**：理想塑性；硬化材料有渐进破坏。

**典型应用**：钢框架塑性设计、钢筋混凝土。

#### 知识点8. 塑性势与Drucker公设

**Drucker 公设**：稳定材料在应力循环中附加应力做功非负：
$$\oint(\sigma_{ij}-\sigma^0_{ij})\dot\varepsilon^p_{ij}\,dt\ge0.$$

**推论**：屈服面外凸；塑性应变率沿外法向（关联流动）
$$\dot\varepsilon^p_{ij}=\dot\lambda\frac{\partial f}{\partial\sigma_{ij}},\qquad \dot\lambda\ge0.$$
**推导**：由公设对任意 $\sigma^0$ 在屈服面内成立，导出 $(\sigma_{ij}-\sigma^0_{ij})\dot\varepsilon^p_{ij}\ge0$，取极限得外法向。

**物理/几何意义**：保证塑性功非负、解唯一稳定。

**极限/特例**：非关联材料（土、岩）不满足。

**典型应用**：岩土塑性、本构公理。

### 第14章 粘弹性与流变模型

#### 知识点1. Maxwell模型

**构造**：弹簧（$E$）与黏壶（$\eta$）串联。

**方程**：
$$\dot\varepsilon=\frac{\dot\sigma}{E}+\frac{\sigma}{\eta}.$$
**推导**：串联应力相同，总应变率 $\dot\varepsilon=\dot\varepsilon_e+\dot\varepsilon_v=\dot\sigma/E+\sigma/\eta$。
蠕变（$\sigma=\sigma_0$）：$\varepsilon=\sigma_0/E+\sigma_0t/\eta$（稳态流动）。
松弛（$\varepsilon=\varepsilon_0$）：$\sigma=E\varepsilon_0e^{-t/\tau}$，$\tau=\eta/E$。

**物理/几何意义**：瞬时弹性+黏性流动，无长期平衡。

**极限/特例**：$\eta\to\infty$ 弹性；$E\to\infty$ 牛顿流体。

**典型应用**：应力松弛、聚合物。

#### 知识点2. Kelvin–Voigt模型

**构造**：弹簧与黏壶并联。

**方程**：
$$\sigma=E\varepsilon+\eta\dot\varepsilon.$$
**推导**：并联应变相同，应力相加。
蠕变：$\varepsilon=\frac{\sigma_0}{E}(1-e^{-t/\tau})$，$\tau=\eta/E$（延迟弹性）。
松弛：$t\to$ 快速时无瞬时应变。

**物理/几何意义**：延迟弹性，无瞬时弹性。

**极限/特例**：$\eta\to0$ 弹性；$E\to0$ 纯黏性。

**典型应用**：蠕变回复、阻尼器。

#### 知识点3. 标准线性固体（三参量）

**构造**：Maxwell 与弹簧并联，或 Kelvin 与弹簧串联。

**方程**：
$$\sigma+p_1\dot\sigma=q_0\varepsilon+q_1\dot\varepsilon,$$
$$p_1=\frac{\eta}{E_1+E_2},\ q_0=\frac{E_1E_2}{E_1+E_2},\ q_1=\frac{E_1\eta}{E_1+E_2}.$$
**推导**：对弹簧并联 Maxwell（$E_1$ 弹簧 $+\ E_2,\eta$ Maxwell）：$\sigma=\sigma_{s}+E_1\varepsilon$，$\dot\varepsilon=\dot\sigma_s/E_2+\sigma_s/\eta$，消去 $\sigma_s$ 得。

**物理/几何意义**：兼具瞬时弹性、延迟弹性与松弛，有平衡模量。

**极限/特例**：$E_1\to0$ 退回 Maxwell；$\eta\to\infty$ 弹性。

**典型应用**：玻璃态聚合物、生物组织。

#### 知识点4. 玻尔兹曼叠加原理与折积

**陈述**：线性粘弹性材料的应变是应力历史的线性叠加
$$\varepsilon(t)=\int_{-\infty}^{t}J(t-\tau)\dot\sigma(\tau)\,d\tau=J(0)\sigma(t)+\int_0^t\dot J(t-\tau)\sigma(\tau)\,d\tau,$$
其中 $J(t)$ 为蠕变柔量。对阶跃应力 $\sigma=\sigma_0H(t)$，
$$\varepsilon(t)=\int_0^tJ(t-\tau)\sigma_0\delta(\tau)d\tau=J(t)\sigma_0.$$

**推导**：把应力历史分解为阶跃之和，线性响应叠加，取极限得折积。

**物理/几何意义**：记忆效应，满足线性叠加。

**极限/特例**：$J(t)=1/E$ 弹性。

**典型应用**：蠕变/松弛数据拟合、对应原理。

#### 知识点5. 复模量与损耗因子

**定义**：谐波激励 $\varepsilon=\varepsilon_0e^{i\omega t}$，稳态 $\sigma=\sigma_0e^{i(\omega t+\delta)}$，
$$E^*(\omega)=E'(\omega)+iE''(\omega),\qquad \tan\delta=\frac{E''}{E'}.$$

**推导**：对 $p_0\sigma+p_1\dot\sigma=q_0\varepsilon+q_1\dot\varepsilon$ 代入谐波，得
$$E^*=\frac{q_0+i\omega q_1}{1+i\omega p_1}.$$
实部 $E'$ 储能模量，虚部 $E''$ 损耗模量。

**物理/几何意义**：滞回环面积 $\Delta W=\pi E''\varepsilon_0^2$ 为单位体积耗散。

**极限/特例**：弹性 $E''=0$；纯黏性 $E'\propto\omega$。

**典型应用**：阻尼材料、DMA 测试、减振。

#### 知识点6. 对应原理

**陈述**：弹性解中弹性常数用复模量替换，即得对应粘弹性谐稳态解。若弹性准静态解含 $E,\nu$，则粘弹性谐响应解取相同形式但用 $E^*(\omega),\nu^*(\omega)$。

**推导**：粘弹性本构在频域为代数关系（$E^*$），与弹性本构形式一致，边值问题在频域是同一代数方程，故解同构。

**物理/几何意义**：用弹性解析解直接得粘弹性频域解。

**极限/特例**：Laplace 域对应原理用 $\bar E(s)$。

**典型应用**：粘弹性梁/板的动力响应。

#### 知识点7. 时温等效与WLF方程

**公式**：时间-温度叠加
$$E(t,T)=E(t/a_T,T_0),\qquad \log a_T=-\frac{C_1(T-T_0)}{C_2+(T-T_0)}.$$
**推导**：自由体积理论：黏度随温度按 Doolittle 关系，取参考温度展开得 WLF。$C_1=17.44,C_2=51.6$（$T_g$ 参考）。

**物理/几何意义**：高温等效短时，主曲线覆盖宽时间域。

**极限/特例**：$T$ 远离 $T_g$ 用 Arrhenius。

**典型应用**：聚合物长期性能预测。

#### 知识点8. 拟弹性与蠕变柔量/松弛模量关系

**公式**：
$$\int_0^t E(t-\tau)J(\tau)\,d\tau=t\quad(\text{折积关系}),\qquad \hat E(s)\hat J(s)=1/s^2.$$
**推导**：对应原理在 Laplace 域的乘积关系，卷积定理得时域折积。

**物理/几何意义**：松弛模量与蠕变柔量互逆（折积意义下）。

**极限/特例**：标准固体可用 Prony 级数表示。

**典型应用**：参数识别、有限元粘弹性。

### 第15章 断裂力学基础

#### 知识点1. 应力强度因子与裂纹尖端场

**定义**：I/II/III 型裂纹尖端应力场
$$\sigma_{ij}=\frac{K}{\sqrt{2\pi r}}f_{ij}(\theta).$$
I 型：
$$\sigma_{xx}=\frac{K_I}{\sqrt{2\pi r}}\cos\frac\theta2\left(1-\sin\frac\theta2\sin\frac{3\theta}{2}\right),$$
$$\sigma_{yy}=\frac{K_I}{\sqrt{2\pi r}}\cos\frac\theta2\left(1+\sin\frac\theta2\sin\frac{3\theta}{2}\right),$$
$$\sigma_{xy}=\frac{K_I}{\sqrt{2\pi r}}\cos\frac\theta2\sin\frac\theta2\cos\frac{3\theta}{2}.$$

**推导（Westergaard）**：取 Westergaard 应力函数
$$Z_I(z)=\frac{K_I}{\sqrt{2\pi z}},$$
由 $\sigma_{xx}=\mathrm{Re}Z-y\,\mathrm{Im}Z'$ 等，取 $z=re^{i\theta}$ 展开即得。

**物理/几何意义**：$K$ 度量尖端奇异性强度，$\sigma\sim r^{-1/2}$ 奇异。

**极限/特例**：无限板中心裂纹 $K_I=\sigma\sqrt{\pi a}$；单边裂纹 $K_I=1.12\sigma\sqrt{\pi a}$。

**典型应用**：断裂判据 $K_I=K_{IC}$、剩余寿命。

#### 知识点2. Griffith能量释放率

**定义**：裂纹扩展能量释放率
$$G=-\frac{\partial\Pi}{\partial a}.$$

**Griffith 准则**：$G\ge G_c=2\gamma_s$（表面能）。
无限板远处拉应力 $\sigma$、中心穿透裂纹长 $2a$：
$$\Pi=-\frac{\sigma^2\pi a^2}{E'}B+2(2a)\gamma_sB,$$
**推导**：含裂纹的应变能降低 $\Delta U=\frac{\sigma^2\pi a^2}{E'}$（Inglis），表面能增 $4a\gamma_sB$。令 $\partial\Pi/\partial a=0$：
$$-\frac{\partial}{\partial a}\left(\frac{\sigma^2\pi a^2}{E'}\right)+4\gamma_s=0\Rightarrow\frac{\pi\sigma^2 a}{E'}=2\gamma_s,$$
即 $\sigma_c=\sqrt{\frac{2E'\gamma_s}{\pi a}}$，其中 $E'=E$（平面应力）或 $E/(1-\nu^2)$（平面应变）。

**物理/几何意义**：裂纹扩展由能量释放驱动、表面能抵抗。

**极限/特例**：$a\to0$ 时 $\sigma_c\to\infty$（理论强度）。Orowan 修正加塑性耗散。

**典型应用**：脆性材料断裂、玻璃、陶瓷。

#### 知识点3. Irwin关系（G–K）

**公式**：
$$G=\frac{K^2}{E'}\quad(\text{单一型}),$$
$$G=\frac{K_I^2+K_{II}^2}{E'}+\frac{K_{III}^2}{2\mu}.$$

**推导**：由裂纹闭合积分（Irwin）
$$G=\lim_{\delta\to0}\frac{1}{2\delta}\int_0^\delta\sigma_{yy}(r)\,\Delta u_y(\delta-r)\,dr,$$
其中 $\Delta u_y=\frac{4K_I}{E'}\sqrt{\frac{\delta-r}{2\pi}}$，$\sigma_{yy}=\frac{K_I}{\sqrt{2\pi r}}$，积分得 $G=K_I^2/E'$。

**物理/几何意义**：能量释放率与应力强度因子等价。

**极限/特例**：平面应力 $E'=E$；平面应变 $E'=E/(1-\nu^2)$。

**典型应用**：把能量法与 $K$ 判据统一。

#### 知识点4. Westergaard方法与裂纹尖端位移

**公式**：I 型裂纹位移
$$u_y=\frac{K_I}{2\mu}\sqrt{\frac{r}{2\pi}}\sin\frac\theta2\left(\kappa+1-2\cos^2\frac\theta2\right),$$
$$u_x=\frac{K_I}{2\mu}\sqrt{\frac{r}{2\pi}}\cos\frac\theta2\left(\kappa-1+2\sin^2\frac\theta2\right).$$
裂纹面张开位移（$\theta=\pi$）：精确解 $\Delta u_y=\frac{2K_I}{E'\sqrt{\pi a}}\sqrt{a^2-x^2}$，近尖端 $r=a-x\to0$ 时 $\Delta u_y\approx\frac{4K_I}{E'}\sqrt{\frac{r}{2\pi}}$。

**推导**：由复势与 Kolosov 公式，取 Westergaard 函数积分。

**物理/几何意义**：裂尖位移场，用于 COD 判据。

**极限/特例**：塑性区修正（Irwin）$r_p=\frac{1}{2\pi}\frac{K^2}{\sigma_Y^2}$。

**典型应用**：COD、弹塑性断裂。

#### 知识点5. J积分

**定义**：
$$J=\int_\Gamma\left(W n_x-\sigma_{ij}n_j\frac{\partial u_i}{\partial x}\right)d\Gamma,$$
$\Gamma$ 为绕裂尖逆时针回路，$W$ 为应变能密度。

**推导（路径无关）**：由散度定理与平衡、无体力，证明对任意两条环绕裂尖的路径 $J$ 相等。同时 $J=G$。

**物理/几何意义**：J 是能量释放率，适用于弹塑性（HRR 场）。

**极限/特例**：线弹性 $J=G=K^2/E'$。

**典型应用**：弹塑性断裂、核电管道。

#### 知识点6. 复合型断裂准则

**最大周向应力准则**：裂纹沿 $\theta_0$ 扩展，满足
$$\left.\frac{\partial\sigma_{\theta\theta}}{\partial\theta}\right|_{\theta_0}=0\Rightarrow K_I\sin\theta_0+K_{II}(3\cos\theta_0-1)=0,$$
$$\sigma_{\theta\theta}(\theta_0)=\frac{1}{\sqrt{2\pi r}}\cos\frac{\theta_0}{2}\left(K_I\cos^2\frac{\theta_0}{2}-\frac32K_{II}\sin\theta_0\right)=K_{IC}.$$

**推导**：取裂尖周向应力极大方向为扩展方向（周向应力准则）。

**最小应变能密度因子 S 准则**：$S=a_{11}K_I^2+2a_{12}K_IK_{II}+a_{22}K_{II}^2$，取 $dS/d\theta=0$ 定方向，$S_{\min}=S_c$ 破裂。

**物理/几何意义**：混合型裂纹不沿原裂纹面扩展。

**极限/特例**：纯 I 型 $\theta_0=0$；纯 II 型 $\theta_0\approx-70.5^\circ$。

**典型应用**：疲劳裂纹转折、复合载荷。

#### 知识点7. 疲劳裂纹扩展（Paris律）

**公式**：
$$\frac{da}{dN}=C(\Delta K)^m,\qquad \Delta K=K_{\max}-K_{\min}.$$

**推导**：由实验结果，双对数 $\log(da/dN)$ vs $\log\Delta K$ 呈线性。门槛值 $\Delta K_{th}$ 以下不扩展；$K_{\max}\to K_{IC}$ 加速。

**物理/几何意义**：每循环裂纹扩展量随应力强度因子幅幂律增长。

**极限/特例**：$R$ 比影响；Forman 公式含 $K_{IC}$。

**典型应用**：剩余寿命积分
$$N_f=\int_{a_0}^{a_c}\frac{da}{C(\Delta K)^m}.$$

#### 知识点8. 断裂韧性与尺寸效应

**定义**：$K_{IC}$ 为平面应变断裂韧性。厚度足够大时 $K_C$ 趋于 $K_{IC}$。

**公式**：塑性区 $r_y=\frac{1}{2\pi}\left(\frac{K_I}{\sigma_Y}\right)^2$（平面应力），平面应变 $r_y=\frac{1}{6\pi}\left(\frac{K_I}{\sigma_Y}\right)^2$。有效性条件
$$B,a,W-a\ge2.5\left(\frac{K_{IC}}{\sigma_Y}\right)^2.$$

**物理/几何意义**：只有小范围屈服时 $K$ 才有效。

**极限/特例**：薄板 $K_C$ 高于 $K_{IC}$。

**典型应用**：材料选型、损伤容限设计。

### 第16章 有限元法与数值方法初步

#### 知识点1. 加权残量与弱形式

**定义**：微分方程 $L(u)=f$，残量 $R=L(u_N)-f$，要求
$$\int_B W_i R\,dV=0.$$
Galerkin 取 $W_i=N_i$。

**弱形式推导（一维）**：对 $-(EAu')'=f$，乘检验函数 $v$ 积分：
$$\int_0^l EAu'v'\,dx=\int_0^l fv\,dx+[EAu'v]_0^l.$$
分部积分降低导数阶，自然边界条件进入右端。

**物理/几何意义**：弱形式允许低阶试探函数，便于分段多项式。

**极限/特例**：Galerkin 与 Ritz 对自伴问题等价。

**典型应用**：FEM、谱元法。

#### 知识点2. 单元刚度矩阵与插值

**定义**：一维杆单元，节点位移 $\mathbf d=(u_1,u_2)$，线性形函数
$$N_1=1-\frac{x}{l},\qquad N_2=\frac{x}{l}.$$

**推导**：
$$u=N_1u_1+N_2u_2,\qquad \varepsilon=\frac{du}{dx}=\left[-\frac1l,\frac1l\right]\mathbf d=\mathbf B\mathbf d.$$
单元刚度
$$\mathbf k=\int_0^l \mathbf B^{\mathrm T}EA\mathbf B\,dx=\frac{EA}{l}\begin{pmatrix}1&-1\\-1&1\end{pmatrix}.$$
**推导**：$\mathbf B$ 为常数，积分得。载荷向量 $\mathbf f=\int\mathbf N^{\mathrm T}q\,dx$。

**物理/几何意义**：离散化把连续体变为弹簧网络。

**极限/特例**：高阶单元用二次/三次形函数。

**典型应用**：桁架、框架、连续体。

#### 知识点3. 总体组装与边界条件

**组装**：$\mathbf K=\sum_e\mathbf A_e^{\mathrm T}\mathbf k_e\mathbf A_e$，$\mathbf K\mathbf d=\mathbf F$。
**推导**：系统的总势能是单元势能之和，对公共节点位移求极值，节点贡献相加。

**边界条件施行**：置零行列或罚函数/乘子法。

**物理/几何意义**：稀疏带状矩阵，节点自由度耦合。

**极限/特例**：对称结构可缩聚。

**典型应用**：所有商业 FEM 软件。

#### 知识点4. 二维/三维单元与等参元

**定义**：四边形 4 节点双线性，$N_i=\frac14(1+\xi\xi_i)(1+\eta\eta_i)$。等参映射
$$x=\sum N_i(\xi,\eta)x_i,\qquad u=\sum N_i(\xi,\eta)u_i.$$
Jacobian
$$J=\begin{vmatrix}\frac{\partial x}{\partial\xi}&\frac{\partial y}{\partial\xi}\\\frac{\partial x}{\partial\eta}&\frac{\partial y}{\partial\eta}\end{vmatrix},\qquad \mathbf B=\mathbf J^{-1}\begin{pmatrix}\partial_\xi N\\\partial_\eta N\end{pmatrix}.$$
刚度
$$\mathbf k=\int_{-1}^1\int_{-1}^1\mathbf B^{\mathrm T}\mathbf D\mathbf B\,|J|\,d\xi\,d\eta\quad(\text{Gauss 积分}).$$

**推导**：由链式法则把 $(\xi,\eta)$ 导数转成 $(x,y)$ 导数。

**物理/几何意义**：等参元可表示曲边，坐标与位移同插值。

**极限/特例**：四节点单元出现剪切闭锁需减缩积分。

**典型应用**：通用结构分析。

#### 知识点5. 收敛性、完备性与Patch Test

**收敛条件**：(1) 完备性（能表示刚体位移与常应变）；(2) 协调性（$C^0$ 连续）；(3) 分片试验通过。

**推导**：由有限元势能 $\Pi_N\ge\Pi$ 的单调性与误差估计
$$\|u-u_N\|_1\le Ch^p\|u\|_{p+1},$$
其中 $h$ 单元尺寸、$p$ 完备阶。

**物理/几何意义**：网格加密解收敛到真解。

**极限/特例**：畸变网格降低精度。

**典型应用**：网格收敛性验证。

#### 知识点6. 动力学有限元与模态分析

**方程**：
$$\mathbf M\ddot{\mathbf d}+\mathbf C\dot{\mathbf d}+\mathbf K\mathbf d=\mathbf F(t),$$
其中 $\mathbf M=\int\rho\mathbf N^{\mathrm T}\mathbf N\,dV$。
**模态**：$\mathbf K\boldsymbol\phi=\omega^2\mathbf M\boldsymbol\phi$。
**推导**：由动能 $T=\frac12\dot{\mathbf d}^{\mathrm T}\mathbf M\dot{\mathbf d}$、势能 $\frac12\mathbf d^{\mathrm T}\mathbf K\mathbf d$，拉格朗日方程得运动方程。

**物理/几何意义**：广义特征值给固有频率与模态；模态叠加解耦。

**极限/特例**：Rayleigh 商 $\omega^2=\frac{\boldsymbol\phi^{\mathrm T}\mathbf K\boldsymbol\phi}{\boldsymbol\phi^{\mathrm T}\mathbf M\boldsymbol\phi}$。

**典型应用**：振动、地震响应、模态试验。

#### 知识点7. 时间积分与数值稳定性

**Newmark-β法**：
$$\mathbf d_{n+1}=\mathbf d_n+\Delta t\dot{\mathbf d}_n+\Delta t^2\left[\left(\tfrac12-\beta\right)\ddot{\mathbf d}_n+\beta\ddot{\mathbf d}_{n+1}\right],$$
$$\dot{\mathbf d}_{n+1}=\dot{\mathbf d}_n+\Delta t\left[(1-\gamma)\ddot{\mathbf d}_n+\gamma\ddot{\mathbf d}_{n+1}\right].$$
**稳定性**：$\beta\ge1/4,\ \gamma=1/2$ 无条件稳定；$\beta=0,\gamma=1/2$ 中心差分条件稳定（$\Delta t\le2/\omega_{\max}$）。
**推导**：将加速度假设为常数/线性变化，积分得递推式；对线性系统分析放大矩阵谱半径 $\rho\le1$ 得稳定条件。

**物理/几何意义**：隐式稳定但每步需求解；显式适合冲击波。

**典型应用**：碰撞、地震、波传播。

#### 知识点8. 变分 crime：不可压缩与体积闭锁

**问题**：近不可压（$\nu\to0.5$）时位移法有限元体积闭锁。

**方法**：混合 u-p 列式，用 Lagrange 乘子/罚函数引入 $\nabla\cdot\mathbf u=0$ 或 $p=K\nabla\cdot\mathbf u$。Hu–Washizu 三场变分。

**推导**：势能泛函加约束项
$$\Pi_{u,p}=\int_B\left(\mu\,e_{ij}e_{ij}-\frac{1}{2K}p^2\right)dV+\int_B p\,\nabla\cdot\mathbf u\,dV,$$
其中 $e_{ij}=\varepsilon_{ij}-\frac13\varepsilon_{kk}\delta_{ij}$。对 $p$ 取驻值 $\delta_p\Pi=0$ 得 $-p/K+\nabla\cdot\mathbf u=0$，即 $p=K\nabla\cdot\mathbf u$。$p$ 与 $\mathbf u$ 独立插值可消除体积闭锁。

**物理/几何意义**：位移插值不能同时满足常应变与不可压约束。

**极限/特例**：B-bar、F-bar 方法。

**典型应用**：橡胶、塑性、生物软组织有限元。

---
⬆ 上级领域：[[力学]]
