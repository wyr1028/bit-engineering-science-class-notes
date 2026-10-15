---
课程名称: 人工智能与计算科学A
学分: 2
总学时: 32
考核方式: 考试
课程属性: 必修
课程性质: 公共基础课程
开课单位: 计算机学院
领域: 计算与人工智能
tags:
  - 课程
  - 领域/计算与人工智能
---

# 人工智能与计算科学A


### 第1章 计算科学的数学基础与误差分析

#### 知识点1.1 浮点数系统与机器精度

**定义** 设基数为 $\beta$，有效位数为 $p$，规格化浮点数可表示为
$$x = \pm \left(1 + \sum_{i=1}^{p-1} d_i \beta^{-i}\right)\beta^{e},\quad 0\le d_i\le \beta-1.$$
**假设** 采用 IEEE-754 双精度：$\beta=2,\ p=53$，指数范围决定溢出与下溢。

**公式与推导** 相邻浮点数间距为机器 epsilon
$$\text{eps} = \beta^{1-p}.$$
$\because$ 尾数最低位权重为 $\beta^{-(p-1)}$，规格化后间隔恰为 $\beta^{-(p-1)}$，
$$\therefore \text{eps}=\beta^{1-p}=2^{-52}\approx 2.22\times10^{-16}.$$

**物理/几何意义** eps 表示“1 与下一个可表示数”的间距，是相对舍入误差的上界。

**极限/特例** 单精度 $p=24$，$\text{eps}=2^{-23}\approx1.19\times10^{-7}$。

**典型应用** 判断两实数相等须用 $|x-y|\le \text{eps}\cdot\max(|x|,|y|)$，而非 $x=y$。

#### 知识点1.2 误差来源、绝对误差与相对误差

**定义** 设精确值 $x$，近似值 $\hat x$。绝对误差 $e_a=\hat x-x$；相对误差 $e_r=e_a/x$（$x\ne0$）。

**假设** 计算过程满足浮点运算模型
$$\text{fl}(a\star b)=(a\star b)(1+\delta),\quad |\delta|\le u,\ \star\in\{+,-,\times,\div\},$$
其中单位舍入误差 $u=\text{eps}/2$。

**公式与推导** 加减法误差放大：
$$\because \text{fl}(a+b)=a+b+(\eta_a+\eta_b),\ |\eta_a|\le u|a|,\ |\eta_b|\le u|b|,$$
$$\therefore \frac{|e_a(a+b)|}{|a+b|}\le u\frac{|a|+|b|}{|a+b|}.$$
当 $a\approx -b$ 时右端分母趋零，产生**相消误差**。

**物理/几何意义** 有效数字因相消而损失，几何上相当于用差很小的两个大量做减法。

**极限/特例** $a=1+\epsilon,\ b=-1$，则 $a+b=\epsilon$，相对误差可由 $u$ 放大到 $u/\epsilon$。

**典型应用** 求解二次方程 $ax^2+bx+c=0$，当 $b^2\gg 4ac$ 时用 $x_1=- \frac{2c}{b\pm\sqrt{b^2-4ac}}$ 避免相消。

#### 知识点1.3 截断误差与 Taylor 展开

**定义** 用有限项近似无穷级数产生的误差称为截断误差。

**假设** $f\in C^{n+1}[a,b]$。

**公式与推导** 带 Lagrange 余项的 Taylor 公式
$$f(x+h)=\sum_{k=0}^{n}\frac{f^{(k)}(x)}{k!}h^k+\frac{f^{(n+1)}(\xi)}{(n+1)!}h^{n+1},\ \xi\in(x,x+h).$$
$\because$ 令 $R_n(h)=f(x+h)-P_n(h)$，反复罗尔定理可证存在 $\xi$。
$$\therefore |R_n|\le \frac{\max|f^{(n+1)}|}{(n+1)!}|h|^{n+1}=O(h^{n+1}).$$

**物理/几何意义** 用 $n$ 次多项式在一点附近拟合光滑函数，误差由 $(n+1)$ 阶导数控制。

**极限/特例** $n=1$：线性近似 $f(x+h)\approx f(x)+f'(x)h$，误差 $O(h^2)$。

**典型应用** 数值微分、误差阶估计、构造差分格式。

#### 知识点1.4 收敛阶与相容性

**定义** 若误差 $E(h)=C h^p+O(h^{p+1})$，则方法具有 $p$ 阶收敛。

**假设** 离散化参数 $h\to0$。

**公式与推导** 由 $E(h)=Ch^p$ 得
$$\because \frac{E(h)}{E(h/2)}=\frac{Ch^p}{C(h/2)^p}=2^p,\quad \therefore p=\log_2\frac{E(h)}{E(h/2)}.$$

**物理/几何意义** 网格加密一倍，误差按 $2^{-p}$ 缩小。

**极限/特例** $p=1$ 线性收敛，$p=2$ 二次收敛。

**典型应用** 用两套网格结果估计实际收敛阶（网格收敛性验证）。

#### 知识点1.5 数值稳定性与条件数

**定义** 算法对输入扰动不敏感称为数值稳定；问题本身对扰动敏感性用条件数度量。

**假设** 计算映射 $y=f(x)$ 可微。

**公式与推导** 条件数
$$\kappa = \left|\frac{\mathrm{d}y/y}{\mathrm{d}x/x}\right|=\left|\frac{x f'(x)}{f(x)}\right|.$$
$\because$ $dy=f'(x)dx$，
$$\therefore \frac{dy}{y}=\frac{xf'(x)}{f(x)}\frac{dx}{x}.$$

**物理/几何意义** $\kappa$ 为相对误差放大因子。

**极限/特例** $\kappa\gg1$ 为病态问题，几何上函数曲线陡峭。

**典型应用** 线性方程组条件数 $\kappa(A)=\|A\|\|A^{-1}\|$。

#### 知识点1.6 前向误差与后向误差

**定义** 前向误差为计算解与真解之差；后向误差为“使计算解成为精确解”所需的最小扰动。

**假设** 映射 $f$，计算解 $\hat y=f(\hat x)$。

**公式与推导** 后向误差 $\eta$ 满足 $f(x+\Delta x)=\hat y$，相对后向误差 $\|\Delta x\|/\|x\|$。
$$\because \hat y-y=f'(x)\Delta x,\quad \therefore \frac{|\hat y-y|}{|y|}\le\kappa\frac{\|\Delta x\|}{\|x\|}.$$

**物理/几何意义** 后向稳定算法把小后向误差转为小前向误差。

**极限/特例** 若后向误差 $O(u)$ 且 $\kappa$ 有界，则前向误差 $O(u)$。

**典型应用** 高斯消元的后向稳定性分析。

#### 知识点1.7 量纲分析与无量纲化

**定义** 用基本量纲 $[M],[L],[T]$ 表示物理量。

**假设** 物理规律与单位无关（Buckingham $\Pi$ 定理）。

**公式与推导** 若 $n$ 个变量含 $k$ 个独立量纲，则存在 $n-k$ 个无量纲数。以单摆为例 $\theta,T,g,L$，
$$T=C\sqrt{L/g}\ \Rightarrow\ \Pi=\frac{T}{\sqrt{L/g}}=const.$$

**物理/几何意义** 无量纲数刻画相似性，便于数值实验。

**极限/特例** Reynolds 数 $Re=\rho U L/\mu$ 表征惯性/黏性比。

**典型应用** 风洞与数值模拟相似准则。

#### 知识点1.8 病态矩阵与条件数估计

**定义** 对方程 $Ax=b$，定义 $\kappa(A)=\|A\|\|A^{-1}\|$。

**假设** $A$ 非奇异。

**公式与推导** 设 $b$ 扰动 $\delta b$，解扰动 $\delta x$：
$$\because A\delta x=\delta b,\quad \therefore \delta x=A^{-1}\delta b,\ \|\delta x\|\le\|A^{-1}\|\|\delta b\|.$$
又 $\|b\|\le\|A\|\|x\|$，
$$\therefore \frac{\|\delta x\|}{\|x\|}\le \kappa(A)\frac{\|\delta b\|}{\|b\|}.$$
若 $A$ 亦扰动，可得 $\frac{\|\delta x\|}{\|x\|}\le\frac{\kappa(A)}{1-\kappa(A)\|\delta A\|/\|A\|}\left(\frac{\|\delta A\|}{\|A\|}+\frac{\|\delta b\|}{\|b\|}\right)$。

**物理/几何意义** 谱条件数等于最大与最小奇异值之比，几何上是超椭球最扁方向之比。

**极限/特例** 正交矩阵 $\kappa_2=1$。

**典型应用** 判断迭代解法是否可靠，预处理降条件数。

### 第2章 插值与函数逼近

#### 知识点2.1 多项式插值的存在唯一性

**定义** 给定 $n+1$ 个节点 $(x_i,y_i)$，求次数不超过 $n$ 的多项式 $P_n$ 满足 $P_n(x_i)=y_i$。

**假设** $x_i$ 互异。

**公式与推导** 设 $P_n(x)=\sum_{k=0}^n a_k x^k$，代入得 Vandermonde 方程
$$Va=y,\quad V_{ij}=x_i^{j},\ i,j=0,\dots,n.$$
$\because \det V=\prod_{0\le i<j\le n}(x_j-x_i)\ne0$，
$$\therefore a \text{ 唯一存在}.$$

**物理/几何意义** $n+1$ 个条件确定 $n+1$ 个自由度。

**极限/特例** 节点重合则退化，导数条件需 Hermite 插值。

**典型应用** 数据拟合、函数表构造。

#### 知识点2.2 Lagrange 插值公式

**定义** 基函数
$$\ell_i(x)=\prod_{\substack{j=0\\j\ne i}}^n\frac{x-x_j}{x_i-x_j}.$$

**假设** 节点互异。

**公式与推导** 由 $\ell_i(x_k)=\delta_{ik}$ 得
$$P_n(x)=\sum_{i=0}^n y_i\ell_i(x).$$
$\because$ 每个 $\ell_i$ 次数 $n$ 且满足插值条件，
$$\therefore P_n(x_i)=\sum_j y_j\delta_{ij}=y_i.$$

**物理/几何意义** 用“开关”型基函数线性组合穿过所有点。

**极限/特例** $n=1$ 时 $P_1=y_0\frac{x-x_1}{x_0-x_1}+y_1\frac{x-x_0}{x_1-x_0}$，即线性插值。

**典型应用** 符号推导、误差分析。

#### 知识点2.3 插值余项定理

**定义** 误差 $R_n(x)=f(x)-P_n(x)$。

**假设** $f\in C^{n+1}$，节点互异。

**公式与推导** 固定 $x\ne x_i$，构造辅助函数
$$\phi(t)=f(t)-P_n(t)-\frac{f(x)-P_n(x)}{W(x)}W(t),\quad W(t)=\prod_{i=0}^n(t-x_i).$$
$\because \phi$ 在 $n+2$ 个点为零，反复罗尔定理得存在 $\xi$ 使 $\phi^{(n+1)}(\xi)=0$，且 $W^{(n+1)}=(n+1)!$，$P_n^{(n+1)}=0$，
$$\therefore f(x)-P_n(x)=\frac{f^{(n+1)}(\xi)}{(n+1)!}W(x).$$

**物理/几何意义** 误差由节点分布 $W(x)$ 与高阶导数共同决定。

**极限/特例** 等距节点 $n$ 大时 $W$ 在端点附近巨大，导致 Runge 现象。

**典型应用** 误差估计与节点优化（Chebyshev 节点）。

#### 知识点2.4 Newton 差商插值

**定义** 差商
$$f[x_i]=f(x_i),\quad f[x_i,\dots,x_{i+k}]=\frac{f[x_{i+1},\dots,x_{i+k}]-f[x_i,\dots,x_{i+k-1}]}{x_{i+k}-x_i}.$$

**假设** 节点互异。

**公式与推导** 由插值唯一性可将 $P_n$ 写成
$$P_n(x)=f[x_0]+\sum_{k=1}^n f[x_0,\dots,x_k]\prod_{j=0}^{k-1}(x-x_j).$$
$\because$ 该式逐点满足插值条件且次数 $\le n$，由唯一性即所求。

**物理/几何意义** 逐次加密节点、复用已有计算。

**极限/特例** 差商与导数关系 $f[x_0,\dots,x_n]=\frac{f^{(n)}(\xi)}{n!}$。

**典型应用** 自适应插值、数值微分。

#### 知识点2.5 Hermite 插值与重节点差商

**定义** 同时匹配函数值与导数值的插值。

**假设** $f\in C^{2n+2}$。

**公式与推导** 余项为
$$f(x)-H_{2n+1}(x)=\frac{f^{(2n+2)}(\xi)}{(2n+2)!}\prod_{i=0}^n(x-x_i)^2.$$
$\because$ 与 Lagrange 余项同法构造，$W(t)=\prod(x_i-t)^2$，其 $(2n+2)$ 阶导为 $(2n+2)!$。

**物理/几何意义** 曲线不仅过点，还与函数相切。

**极限/特例** 两点 Hermite 即三次 Hermite，用于有限元形函数。

**典型应用** 梁单元挠度插值。

#### 知识点2.6 Chebyshev 节点与最小最大误差

**定义** $T_n(x)=\cos(n\arccos x)$，$x\in[-1,1]$。

**假设** $f$ 光滑。

**公式与推导** 由 $T_n(\cos\theta)=\cos n\theta$ 及三角恒等式
$$T_{n+1}=2xT_n-T_{n-1}.$$
$\because$ $T_n$ 的极值在 $\theta=k\pi/n$，交替取 $\pm1$，
$$\therefore \min_{\text{节点}}\max_{x\in[-1,1]}|W(x)|=\frac{1}{2^{n}} \text{（取 } T_n \text{ 零点为节点）}.$$

**物理/几何意义** 等幅振荡使误差在区间上均匀分布。

**极限/特例** Chebyshev 零点 $x_k=\cos\frac{(2k+1)\pi}{2n+2}$。

**典型应用** 谱方法、数值积分节点。

#### 知识点2.7 分段插值与三次样条

**定义** 在每个 $[x_i,x_{i+1}]$ 上为三次多项式，且在节点处 $S,S',S''$ 连续。

**假设** 节点递增。

**公式与推导** 记 $M_i=S''(x_i)$，由线性插值
$$S''(x)=\frac{x_{i+1}-x}{h_i}M_i+\frac{x-x_i}{h_i}M_{i+1}.$$
两次积分并利用 $S(x_i)=y_i$ 得
$$S(x)=\frac{(x_{i+1}-x)^3}{6h_i}M_i+\frac{(x-x_i)^3}{6h_i}M_{i+1}+C_i(x-x_i)+D_i(x_{i+1}-x).$$
强制 $S'$ 连续得到三对角方程
$$\frac{h_{i-1}}{6}M_{i-1}+\frac{h_{i-1}+h_i}{3}M_i+\frac{h_i}{6}M_{i+1}=\frac{y_{i+1}-y_i}{h_i}-\frac{y_i-y_{i-1}}{h_{i-1}}.$$

**物理/几何意义** 相当于细梁的曲率场，能抵抗 Runge 振荡。

**极限/特例** 自然边界 $M_0=M_n=0$。

**典型应用** 曲线插值、计算几何。

#### 知识点2.8 最小二乘逼近与正交投影

**定义** 求 $p\in\text{span}\{\phi_j\}$ 使 $\sum_i (y_i-p(x_i))^2$ 最小。

**假设** 基函数线性无关。

**公式与推导** 令 $p=\sum_j c_j\phi_j$，目标 $J(c)=\|A c-y\|_2^2$，$A_{ij}=\phi_j(x_i)$。
$$\because \nabla_c J=2A^\top (Ac-y)=0,\quad \therefore A^\top A c=A^\top y.$$
这就是法方程，几何上残差正交于列空间。

**物理/几何意义** $y$ 在 $\text{col}(A)$ 上的正交投影。

**极限/特例** 若 $A$ 列正交，$c_j=\frac{\langle y,\phi_j\rangle}{\langle\phi_j,\phi_j\rangle}$。

**典型应用** 数据拟合、回归。

#### 知识点2.9 正交多项式与逼近误差

**定义** 关于权 $w$ 正交：$\int_a^b w(x)\phi_m\phi_n dx=\gamma_n\delta_{mn}$。

**假设** $w>0$。

**公式与推导** 由 Gram-Schmidt 递推
$$\phi_{n+1}=(x-\alpha_n)\phi_n-\beta_n\phi_{n-1},$$
$$\alpha_n=\frac{\langle x\phi_n,\phi_n\rangle}{\langle\phi_n,\phi_n\rangle},\quad \beta_n=\frac{\langle\phi_n,\phi_n\rangle}{\langle\phi_{n-1},\phi_{n-1}\rangle}.$$
$\because$ 令 $\phi_{n+1}-x\phi_n$ 在低维空间正交投影，
$$\therefore \text{递推成立}.$$

**物理/几何意义** 在函数空间做正交基，投影即最佳逼近。

**极限/特例** Legendre（$w=1$）、Chebyshev（$w=(1-x^2)^{-1/2}$）。

**典型应用** 谱展开、Gauss 求积。

### 第3章 数值积分与数值微分

#### 知识点3.1 插值型求积公式

**定义** 用插值多项式积分近似 $\int_a^b f\,dx$。

**假设** $f$ 足够光滑。

**公式与推导** 由 $P_n\approx f$，
$$\int_a^b f\,dx\approx \sum_i w_i f(x_i),\quad w_i=\int_a^b \ell_i(x)\,dx.$$
$\because$ 代入 Lagrange 插值并交换求和与积分，
$$\therefore w_i=\int_a^b\prod_{j\ne i}\frac{x-x_j}{x_i-x_j}dx.$$

**物理/几何意义** 用多项式下方面积近似曲线下方面积。

**极限/特例** 梯形 $w_0=w_1=(b-a)/2$。

**典型应用** 构造 Newton-Cotes 族。

#### 知识点3.2 梯形公式与误差

**定义** $\int_a^b f\,dx\approx \frac{b-a}{2}[f(a)+f(b)]$。

**假设** $f\in C^2[a,b]$。

**公式与推导** 由线性插值余项
$$f(x)=P_1(x)+\frac{f''(\xi)}{2}(x-a)(x-b).$$
积分 $(x-a)(x-b)$ 在 $[a,b]$ 上积分为 $-\frac{(b-a)^3}{6}$，由积分中值定理
$$\therefore \int_a^b f\,dx-\frac{b-a}{2}[f(a)+f(b)]=-\frac{(b-a)^3}{12}f''(\eta).$$

**物理/几何意义** 误差正比于曲率与区间三次方。

**极限/特例** $f$ 线性时精确。

**典型应用** 复合梯形法。

#### 知识点3.3 Simpson 公式与误差

**定义** $\int_a^b f\,dx\approx \frac{b-a}{6}\left[f(a)+4f\!\left(\frac{a+b}{2}\right)+f(b)\right]$。

**假设** $f\in C^4[a,b]$。

**公式与推导** 对二次插值积分得该式。为求误差，用三次 Hermite 展开，在对称点上 $(x-a)(x-m)^2(x-b)$ 的积分为零，故三阶项消失，得
$$\int_a^b f\,dx-S=-\frac{(b-a)^5}{2880}f^{(4)}(\eta).$$

**物理/几何意义** 用抛物线下方面积近似，误差由四阶导控制。

**极限/特例** 三次以下多项式精确。

**典型应用** 工程积分标准格式。

#### 知识点3.4 复合求积与收敛阶

**定义** 将 $[a,b]$ 等分 $N$ 段，每段用低阶公式。

**假设** $h=(b-a)/N$。

**公式与推导** 复合梯形
$$T_N=h\left[\frac{f_0+f_N}{2}+\sum_{i=1}^{N-1}f_i\right],\quad E_T=-\frac{b-a}{12}h^2 f''(\eta).$$
复合 Simpson（$N$ 偶）
$$S_N=\frac{h}{3}\left[f_0+f_N+4\!\!\sum_{i\ \text{odd}}\!\!f_i+2\!\!\sum_{i\ \text{even}}\!\!f_i\right],\quad E_S=-\frac{b-a}{180}h^4 f^{(4)}(\eta).$$
$\because$ 误差对子区间求和，$N\propto h^{-1}$，
$$\therefore E_T=O(h^2),\ E_S=O(h^4).$$

**物理/几何意义** 加密网格快速降低误差。

**极限/特例** Romberg 外推可进一步提升阶。

**典型应用** 数值积分库。

#### 知识点3.5 Gauss-Legendre 求积

**定义** 选节点与权重使 $\int_{-1}^1 f\,dx\approx\sum_{i=1}^n w_i f(x_i)$ 对次数 $\le 2n-1$ 的多项式精确。

**假设** 节点为 Legendre 多项式 $P_n$ 的零点。

**公式与推导** 有 $2n$ 个自由参数（节点+权重），可满足 $2n$ 个矩条件（次数 $0$ 到 $2n-1$）。设 $f=P_n q+r$，$\deg r<n$。因 $P_n\perp$ 低次多项式，
$$\int_{-1}^1 P_n q\,dx=0,\quad \int_{-1}^1 f\,dx=\int_{-1}^1 r\,dx=\sum w_i r(x_i)=\sum w_i f(x_i).$$
$\because$ $P_n(x_i)=0$，
$$\therefore \text{代数精度 } 2n-1.$$

**物理/几何意义** 最优节点分布，精度翻倍。

**极限/特例** $n=1$ 为中点公式，$n=2$ 为两点 Gauss。

**典型应用** 有限元单元积分。

#### 知识点3.6 数值微分与差分公式

**定义** 用差商近似导数。

**假设** $f$ 光滑。

**公式与推导** 由 Taylor：
$$f(x+h)=f(x)+hf'(x)+\frac{h^2}{2}f''(x)+\frac{h^3}{6}f'''(x)+\cdots$$
$$\Rightarrow\ f'(x)=\frac{f(x+h)-f(x)}{h}-\frac{h}{2}f''(x)+O(h^2),$$
$$\therefore \text{前向差商误差 }O(h).$$
中心差商：
$$f(x+h)-f(x-h)=2hf'(x)+\frac{h^3}{3}f'''(x)+O(h^5),$$
$$\therefore f'(x)=\frac{f(x+h)-f(x-h)}{2h}-\frac{h^2}{6}f'''(\xi).$$
二阶：
$$f''(x)=\frac{f(x+h)-2f(x)+f(x-h)}{h^2}+O(h^2).$$

**物理/几何意义** 差商是切线斜率/曲率的离散逼近。

**极限/特例** $h\to0$ 时趋近真导数，但舍入误差 $\sim u/h$ 使存在最优 $h\sim u^{1/3}$。

**典型应用** 自动微分验证、PDE 离散。

#### 知识点3.7 舍入误差与最优步长

**定义** 数值微分总误差 = 截断误差 + 舍入误差。

**假设** 函数值有绝对误差 $\epsilon$。

**公式与推导** 中心差商
$$E(h)\approx \frac{h^2}{6}M+\frac{\epsilon}{h}.$$
$$\because \frac{dE}{dh}=\frac{hM}{3}-\frac{\epsilon}{h^2}=0,\quad \therefore h^*=\left(\frac{3\epsilon}{M}\right)^{1/3},\ E_{\min}=O(\epsilon^{2/3}).$$

**物理/几何意义** 步长过小被噪声淹没。

**极限/特例** 前向差商 $h^*=(\epsilon/M)^{1/2}$。

**典型应用** 实验数据求导。

### 第4章 线性方程组的数值解法

#### 知识点4.1 高斯消元与 LU 分解

**定义** 用初等行变换化 $A$ 为上三角。

**假设** 主元非零（或使用部分主元）。

**公式与推导** 第 $k$ 步消元乘子 $m_{ik}=a_{ik}^{(k)}/a_{kk}^{(k)}$，
$$a_{ij}^{(k+1)}=a_{ij}^{(k)}-m_{ik}a_{kj}^{(k)}.$$
将乘子存入 $L$ 下三角，
$$A=LU,\quad L=\begin{pmatrix}1&&\\ m_{21}&1&\\ \vdots&&\ddots\end{pmatrix},\ U \text{ 上三角}.$$
$\because$ 消元等价于左乘单位下三角矩阵的逆，
$$\therefore A=LU.$$

**物理/几何意义** 把一般线性系统化为两个三角系统。

**极限/特例** 对称正定可用 Cholesky $A=LL^\top$。

**典型应用** 每次解多个右端项。

#### 知识点4.2 部分主元与增长因子

**定义** 选当前列绝对值最大元素作主元并交换行。

**假设** $A$ 非奇异。

**公式与推导** 交换矩阵 $P$，则 $PA=LU$。消元中元素满足
$$|a_{ij}^{(k+1)}|\le 2\max|a_{ij}^{(k)}|,$$
$$\therefore \text{增长因子 }\rho\le 2^{n-1}.$$

**物理/几何意义** 抑制舍入放大，几何上避免小主元除法。

**极限/特例** 若 $\rho=O(1)$，后向误差小。

**典型应用** 通用稠密求解器。

#### 知识点4.3 迭代法：Jacobi 与 Gauss-Seidel

**定义** 将 $A=D+L+U$，迭代求解。

**假设** 对角线非零。

**公式与推导** Jacobi：$Dx^{(k+1)}=-(L+U)x^{(k)}+b$，迭代矩阵 $G_J=-D^{-1}(L+U)$。
Gauss-Seidel：$(D+L)x^{(k+1)}=-Ux^{(k)}+b$，$G_{GS}=-(D+L)^{-1}U$。
收敛充要条件
$$\rho(G)<1.$$
$\because$ 误差 $e^{(k+1)}=Ge^{(k)}$，$\therefore \|e^{(k)}\|\le\rho(G)^k\|e^{(0)}\|$。

**物理/几何意义** 逐步松弛，谱半径决定收敛速度。

**极限/特例** 严格对角占优保证收敛。

**典型应用** 大型稀疏系统。

#### 知识点4.4 谱半径与收敛速度

**定义** 渐进收敛因子 $R=-\ln\rho(G)$（每步减少的位数）。

**假设** $\rho(G)<1$。

**公式与推导** 使 $\rho^k\le\epsilon$ 所需步数
$$k\ge \frac{\ln(1/\epsilon)}{-\ln\rho(G)}=\frac{\ln(1/\epsilon)}{R}.$$

**物理/几何意义** $R$ 越大收敛越快。

**极限/特例** SOR 可加速，最优松弛 $\omega_{opt}=\frac{2}{1+\sqrt{1-\rho_J^2}}$。

**典型应用** 迭代求解器配置。

#### 知识点4.5 幂法求主特征值

**定义** 迭代 $v^{(k+1)}=Av^{(k)}/\|Av^{(k)}\|$。

**假设** $A$ 有唯一模最大特征值 $\lambda_1$，$|\lambda_1|>|\lambda_2|\ge\cdots$。

**公式与推导** 将 $v^{(0)}=\sum_i c_i u_i$，
$$A^k v^{(0)}=\sum_i c_i\lambda_i^k u_i=\lambda_1^k\left(c_1u_1+\sum_{i>1}c_i\left(\frac{\lambda_i}{\lambda_1}\right)^k u_i\right).$$
$$\therefore \lambda_1=\lim_{k\to\infty}\frac{(Av^{(k)})_j}{v^{(k)}_j},\quad \text{收敛比 }|\lambda_2/\lambda_1|.$$

**物理/几何意义** 反复作用放大主方向。

**极限/特例** 存在谱隙则快速收敛。

**典型应用** PageRank、模态分析。

#### 知识点4.6 奇异值分解（SVD）

**定义** 任意 $A\in\mathbb{R}^{m\times n}$ 可写 $A=U\Sigma V^\top$。

**假设** $U,V$ 正交，$\Sigma$ 对角含奇异值。

**公式与推导** $A^\top A=V\Sigma^2 V^\top$ 为对称正定谱分解，$\sigma_i=\sqrt{\lambda_i(A^\top A)}$。
$$\because \|Ax\|^2=x^\top A^\top A x,\quad \therefore \sigma_{\max}=\|A\|_2,\ \kappa_2=\sigma_{\max}/\sigma_{\min}.$$

**物理/几何意义** 把单位球映射为超椭球，轴长为奇异值。

**极限/特例** 低秩逼近截断 SVD 为最佳。

**典型应用** 降维、最小二乘稳定解。

#### 知识点4.7 最小二乘的正交化方法

**定义** 用 QR 解 $\min\|Ax-b\|_2$。

**假设** $A$ 列满秩。

**公式与推导** $A=QR$，$Q$ 列正交，$R$ 上三角。
$$\because \|Ax-b\|^2=\|QRx-b\|^2=\|Rx-Q^\top b\|^2+\|(I-QQ^\top)b\|^2,$$
$$\therefore Rx=Q^\top b \text{ 给出最小二乘解}.$$
法方程 $A^\top A=R^\top R$，故 $\kappa(A^\top A)=\kappa(R)^2$，QR 较直接法方程更稳定。

**物理/几何意义** 投影到 $Q$ 的列空间。

**极限/特例** Householder/Givens 实现。

**典型应用** 稳健回归。

### 第5章 常微分方程数值解

#### 知识点5.1 初值问题与 Euler 法

**定义** 求解 $y'=f(t,y),\ y(t_0)=y_0$。

**假设** $f$ 满足 Lipschitz 条件。

**公式与推导** 由向前 Taylor
$$y(t_{n+1})=y(t_n)+hf(t_n,y_n)+\frac{h^2}{2}y''(\xi).$$
$$\therefore y_{n+1}=y_n+hf(t_n,y_n),\quad E_n=O(h^2)\ \text{局部},\ O(h)\ \text{整体}.$$

**物理/几何意义** 沿切线方向前进。

**极限/特例** 隐式 Euler $y_{n+1}=y_n+hf(t_{n+1},y_{n+1})$ 更适合刚性。

**典型应用** 入门与稳定性分析。

#### 知识点5.2 局部与整体截断误差

**定义** 局部误差为一步误差；整体误差为 $N$ 步累积。

**假设** 方法一致。

**公式与推导** 设局部误差 $C h^{p+1}$，全局误差 $e_n=y(t_n)-y_n$ 满足
$$e_{n+1}=e_n+h[f(t_n,y(t_n))-f(t_n,y_n)]+Ch^{p+1}.$$
由 Lipschitz，$|e_{n+1}|\le (1+hL)|e_n|+Ch^{p+1}$，
$$\therefore |e_n|\le \frac{C h^p}{L}(e^{LT}-1)=O(h^p).$$

**物理/几何意义** 每步误差按指数放大，整体阶比局部低一阶。

**极限/特例** $p=1$ Euler，$p=4$ 经典 RK4。

**典型应用** 收敛性证明。

#### 知识点5.3 Runge-Kutta 方法

**定义** 用若干斜率加权平均
$$y_{n+1}=y_n+h\sum_{i=1}^s b_i k_i,\quad k_i=f\!\left(t_n+c_ih,\ y_n+h\sum_j a_{ij}k_j\right).$$

**假设** 系数满足阶条件。

**公式与推导** 以二阶为例，展开
$$k_1=f,\quad k_2=f+h c_2 f_t+h a_{21} f f_y+O(h^2).$$
与 Taylor 解
$$y(t_{n+1})=y_n+hf+\frac{h^2}{2}(f_t+ff_y)+O(h^3)$$
比较得
$$b_1+b_2=1,\quad b_2c_2=\frac12,\quad b_2a_{21}=\frac12.$$
取 $b_2=1/2,\ c_2=1,\ a_{21}=1$ 得改进 Euler。

**物理/几何意义** 多次采样斜率提高精度。

**极限/特例** RK4 每步 4 次函数求值，全局 $O(h^4)$。

**典型应用** 通用 ODE 求解。

#### 知识点5.4 绝对稳定性与 A-稳定

**定义** 对测试方程 $y'=\lambda y$，放大因子 $R(z)$，$z=\lambda h$。若 $\forall \text{Re}(z)<0,\ |R(z)|\le1$，称 A-稳定。

**假设** $\lambda$ 可为复数。

**公式与推导** Euler：$R=1+z$，$|1+z|\le1$ 只在单位圆，故非 A-稳定。
隐式 Euler：$R=1/(1-z)$，$|R|=\frac{1}{|1-z|}$，
$$\because \text{Re}(z)<0\Rightarrow |1-z|>1,\quad \therefore |R|<1.$$
梯形法 $R=\frac{1+z/2}{1-z/2}$ 亦 A-稳定。

**物理/几何意义** 刚性方程要求大步长仍稳定。

**极限/特例** RK4 稳定域有限。

**典型应用** 化学反应、结构动力学。

#### 知识点5.5 线性多步法

**定义** 利用多个历史点
$$\sum_{j=0}^k \alpha_j y_{n+j}=h\sum_{j=0}^k \beta_j f_{n+j}.$$

**假设** 系数满足相容性 $\sum\alpha_j=0,\ \sum\alpha_j j=\sum\beta_j$。

**公式与推导** 在 $t_n$ 处 Taylor 展开两边，消去低阶项得到阶条件。Adams-Bashforth（显式）$\beta_0=0$；Adams-Moulton（隐式）含 $\beta_0$。
$$\because \text{特征多项式 }\rho(\xi)=\sum\alpha_j\xi^j,\ \sigma(\xi)=\sum\beta_j\xi^j,$$
$$\therefore \text{零点在单位圆且单根}\Rightarrow \text{零稳定}.$$

**物理/几何意义** 复用历史函数值，减少每步计算。

**极限/特例** 预测-校正组合。

**典型应用** 长时间积分。

#### 知识点5.6 刚性方程与稳定域

**定义** 特征值实部相差悬殊的线性系统。

**假设** 最大与最小模特征值比 $\gg1$。

**公式与推导** 显式法稳定要求
$$h\le \frac{C}{|\lambda_{\max}|},$$
而精度只要求 $h\sim 1/|\lambda_{\min}|$，两者矛盾。

**物理/几何意义** 快模态限制步长但无需精度。

**极限/特例** 隐式法可用大 $h$。

**典型应用** 含快慢时间尺度系统。

#### 知识点5.7 辛积分与几何积分思想

**定义** 保持 Hamilton 系统辛结构的积分。

**假设** 系统 $\dot q=\partial H/\partial p,\ \dot p=-\partial H/\partial q$。

**公式与推导** Störmer-Verlet
$$p_{n+1/2}=p_n-\frac h2\nabla_qV(q_n),\quad q_{n+1}=q_n+hM^{-1}p_{n+1/2},\quad p_{n+1}=p_{n+1/2}-\frac h2\nabla_qV(q_{n+1}).$$
其 Jacobian 满足 $J^\top\Omega J=\Omega$，
$$\therefore \text{保持面积/辛形式，能量长期有界}.$$

**物理/几何意义** 相空间面积守恒，数值轨道稳定。

**极限/特例** 谐振子下能量不漂移。

**典型应用** 分子动力学、天体力学。

### 第6章 偏微分方程数值解

#### 知识点6.1 有限差分与网格

**定义** 用差商代替偏导数。

**假设** 网格 $(x_i,t_n)$，步长 $\Delta x,\Delta t$。

**公式与推导** 由 Taylor
$$\frac{\partial u}{\partial x}\Big|_i=\frac{u_{i+1}-u_{i-1}}{2\Delta x}+O(\Delta x^2),$$
$$\frac{\partial^2 u}{\partial x^2}\Big|_i=\frac{u_{i+1}-2u_i+u_{i-1}}{\Delta x^2}+O(\Delta x^2).$$

**物理/几何意义** 用邻近点线性组合逼近局部微分。

**极限/特例** 非均匀网格需变系数格式。

**典型应用** 热传导、波动方程。

#### 知识点6.2 热传导方程显式格式与 CFL

**定义** $u_t=\alpha u_{xx}$。

**假设** 显式 FTCS：$u_i^{n+1}=u_i^n+r(u_{i+1}^n-2u_i^n+u_{i-1}^n)$，$r=\alpha\Delta t/\Delta x^2$。

**公式与推导** 用 von Neumann 分析，令 $u_i^n=\xi^n e^{ik i\Delta x}$，代入得
$$\xi=1-2r(1-\cos k\Delta x)=1-4r\sin^2\frac{k\Delta x}{2}.$$
稳定性 $|\xi|\le1$ 对所有 $k$ 成立要求
$$0\le r\le \frac12.$$
$$\therefore \Delta t\le \frac{\Delta x^2}{2\alpha}.$$

**物理/几何意义** 信息传播速度受网格限制（CFL 条件）。

**极限/特例** $r=1/2$ 时 $k=\pi/\Delta x$ 模态 $\xi=0$。

**典型应用** 显式时间推进。

#### 知识点6.3 隐式格式与无条件稳定

**定义** 后向 Euler：$\frac{u_i^{n+1}-u_i^n}{\Delta t}=\alpha\frac{u_{i+1}^{n+1}-2u_i^{n+1}+u_{i-1}^{n+1}}{\Delta x^2}$。

**假设** 每步解三对角系统。

**公式与推导** 令 $r=\alpha\Delta t/\Delta x^2$，
$$-r u_{i-1}^{n+1}+(1+2r)u_i^{n+1}-r u_{i+1}^{n+1}=u_i^n.$$
von Neumann：$\xi=\frac{1}{1+4r\sin^2(k\Delta x/2)}$，
$$\therefore |\xi|\le1\ \forall r>0,\ \text{无条件稳定}.$$

**物理/几何意义** 隐式平均未来信息，允许大步长。

**极限/特例** 一阶精度 $O(\Delta t)$。

**典型应用** 刚性扩散。

#### 知识点6.4 Crank-Nicolson 格式

**定义** 对空间用 $n$ 与 $n+1$ 两层平均
$$\frac{u_i^{n+1}-u_i^n}{\Delta t}=\frac{\alpha}{2}\left(\delta_{xx}u_i^{n+1}+\delta_{xx}u_i^n\right).$$

**假设** 均匀网格。

**公式与推导** von Neumann
$$\xi=\frac{1-2r\sin^2(k\Delta x/2)}{1+2r\sin^2(k\Delta x/2)},\quad r=\frac{\alpha\Delta t}{\Delta x^2}.$$
$$\because \text{分子绝对值}\le\text{分母},\quad \therefore |\xi|\le1\ \text{无条件稳定}.$$
截断误差 $O(\Delta t^2+\Delta x^2)$。

**物理/几何意义** 时间中心差分，二阶精度。

**极限/特例** 对不连续初值有振荡。

**典型应用** 高精度扩散求解。

#### 知识点6.5 Laplace 方程与五点格式

**定义** $\nabla^2 u=0$。

**假设** 矩形网格。

**公式与推导** 五点格式
$$u_{i+1,j}+u_{i-1,j}+u_{i,j+1}+u_{i,j-1}-4u_{i,j}=0.$$
其截断误差为 $\frac{\Delta x^2}{12}(u_{xxxx}+u_{yyyy})+O(\Delta x^4)$。
极值原理说明离散解满足最大值原理。

**物理/几何意义** 每点等于邻域平均。

**极限/特例** 调和函数光滑。

**典型应用** 势流、稳态传热。

#### 知识点6.6 波动方程与 CFL 条件

**定义** $u_{tt}=c^2u_{xx}$。

**假设** 显式中心格式
$$u_i^{n+1}=2u_i^n-u_i^{n-1}+\lambda^2(u_{i+1}^n-2u_i^n+u_{i-1}^n),\quad \lambda=\frac{c\Delta t}{\Delta x}.$$
**公式与推导** von Neumann 令 $u_i^n=\xi^n e^{iki\Delta x}$，得
$$\xi^2-\left(2-4\lambda^2\sin^2\frac{k\Delta x}{2}\right)\xi+1=0.$$
稳定要求 $|\xi|\le1$ 即
$$\lambda\le1\quad\Leftrightarrow\quad \Delta t\le\frac{\Delta x}{c}.$$

**物理/几何意义** 数值信息传播须快于物理波速。

**极限/特例** $\lambda=1$ 精确传播（1D）。

**典型应用** 声学、弹性波。

#### 知识点6.7 有限体积与守恒律

**定义** 对守恒律 $u_t+\nabla\cdot F=0$ 在控制体积上积分。

**假设** 网格单元 $\Omega_i$。

**公式与推导** 由 Gauss 散度定理
$$\frac{d}{dt}\int_{\Omega_i}u\,dV=-\oint_{\partial\Omega_i}F\cdot n\,dS.$$
$$\therefore \frac{u_i^{n+1}-u_i^n}{\Delta t}=-\frac{F_{i+1/2}-F_{i-1/2}}{\Delta x}.$$

**物理/几何意义** 通量差决定守恒量变化。

**极限/特例** 迎风格式保证单调性。

**典型应用** 流体、激波捕捉。

#### 知识点6.8 双曲守恒律与熵条件

**定义** $u_t+f(u)_x=0$。

**假设** $f$ 凸。

**公式与推导** 特征速度 $f'(u)$，特征线 $dx/dt=f'(u)$。特征线相交产生激波，Rankine-Hugoniot
$$s=\frac{f(u_R)-f(u_L)}{u_R-u_L}.$$
熵条件筛除膨胀激波。

**物理/几何意义** 非线性导致波形陡化。

**极限/特例** Burgers 方程 $f=u^2/2$。

**典型应用** CFD 激波。

### 第7章 变分原理与力学能量法

#### 知识点7.1 泛函与变分

**定义** 泛函 $J[y]=\int_a^b F(x,y,y')\,dx$。

**假设** $y\in C^2$，端点固定。

**公式与推导** 设变分 $\delta y$，$\delta y(a)=\delta y(b)=0$。一阶变分
$$\delta J=\int_a^b\left(\frac{\partial F}{\partial y}\delta y+\frac{\partial F}{\partial y'}\delta y'\right)dx.$$
分部积分
$$\int_a^b F_{y'}\delta y'\,dx=\Big[F_{y'}\delta y\Big]_a^b-\int_a^b\frac{d}{dx}F_{y'}\delta y\,dx.$$
由边界项为零及 $\delta y$ 任意，
$$\therefore \frac{\partial F}{\partial y}-\frac{d}{dx}\frac{\partial F}{\partial y'}=0.$$

**物理/几何意义** 泛函驻值对应物理平衡/运动。

**极限/特例** 多点/多函数推广。

**典型应用** 推导欧拉-拉格朗日方程。

#### 知识点7.2 Euler-Lagrange 方程

**定义** 使 $J$ 取驻值的曲线满足上方程。

**假设** $F_{y'y'}\ne0$。

**公式与推导** 展开
$$\frac{d}{dx}F_{y'}=F_{y'x}+F_{y'y}y'+F_{y'y'}y'',$$
$$\therefore F_{y'y'}y''+F_{y'y}y'+F_{y'x}-F_y=0.$$
当 $F$ 不显含 $x$ 时有 Beltrami 恒等式
$$F-y'F_{y'}=const.$$

**物理/几何意义** 二阶非线性 ODE，边界值问题。

**极限/特例** 最速降线、测地线。

**典型应用** 变分法基础。

#### 知识点7.3 最小势能原理

**定义** 弹性体总势能 $\Pi=U-W$，其中 $U$ 为应变能，$W$ 为外力功。

**假设** 小变形、线性弹性、保守力。

**公式与推导** 对位移变分 $\delta u$，
$$\delta\Pi=\int_\Omega \sigma:\delta\epsilon\,dV-\int_\Omega f\cdot\delta u\,dV-\int_{\Gamma_t}\bar t\cdot\delta u\,dS=0.$$
由 $\sigma:\delta\epsilon=\sigma:\nabla_s\delta u$ 及散度定理
$$\int_\Omega\sigma:\delta\epsilon\,dV=\int_{\partial\Omega}(\sigma n)\cdot\delta u\,dS-\int_\Omega(\nabla\cdot\sigma)\cdot\delta u\,dV.$$
$$\therefore \int_\Omega(\nabla\cdot\sigma+f)\cdot\delta u\,dV+\int_{\Gamma_t}(t-\bar t)\cdot\delta u\,dS=0.$$
由 $\delta u$ 任意得平衡方程与力边界条件。

**物理/几何意义** 真实位移使总势能最小。

**极限/特例** 稳定平衡为极小。

**典型应用** FEM 理论基础。

#### 知识点7.4 虚功原理

**定义** 平衡时外力虚功等于内力虚功
$$\int_\Omega\sigma:\delta\epsilon\,dV=\int_\Omega f\cdot\delta u\,dV+\int_{\Gamma_t}\bar t\cdot\delta u\,dS.$$

**假设** 任意协调虚位移。

**公式与推导** 由最小势能 $\delta\Pi=0$ 直接得到；反之由虚功原理结合散度定理得平衡方程，二者等价。

**物理/几何意义** 功的守恒表述，不依赖本构。

**极限/特例** 动力问题加惯性项 $-\int\rho\ddot u\cdot\delta u\,dV$。

**典型应用** 梁、板壳理论。

#### 知识点7.5 加权残值与 Galerkin 法

**定义** 设近似解 $u_h=\sum_j c_j\phi_j$，残差 $R=L u_h-f$，令 $\int w_i R=0$。

**假设** 试探/权函数空间。

**公式与推导** Galerkin 取 $w_i=\phi_i$，
$$\int_\Omega \phi_i L u_h\,dV=\int_\Omega \phi_i f\,dV.$$
对含二阶导算子，分部积分降低阶数，自然引入弱形式。

**物理/几何意义** 残差在试探空间正交投影为零。

**极限/特例** 最小二乘、配点法是不同权函数选择。

**典型应用** FEM、谱方法。

#### 知识点7.6 哈密顿原理与最小作用量

**定义** 作用量 $S=\int_{t_1}^{t_2}L(q,\dot q,t)\,dt$。

**假设** 端点 $\delta q(t_1)=\delta q(t_2)=0$。

**公式与推导** $\delta S=0$ 给出
$$\frac{d}{dt}\frac{\partial L}{\partial\dot q_i}-\frac{\partial L}{\partial q_i}=0,$$
其中 $L=T-V$。推导同 Euler-Lagrange。

**物理/几何意义** 真实运动使作用量取驻值。

**极限/特例** 相对论、场论推广。

**典型应用** 分析力学统一框架。

#### 知识点7.7 诺特定理与守恒量

**定义** 作用量对连续对称性不变对应守恒量。

**假设** $L$ 在变换下不变。

**公式与推导** 设 $q_i\to q_i+\epsilon K_i$，
$$\delta L=\epsilon\left(\frac{\partial L}{\partial q_i}K_i+\frac{\partial L}{\partial\dot q_i}\dot K_i\right)=\epsilon\frac{d}{dt}\left(\frac{\partial L}{\partial\dot q_i}K_i\right).$$
由 E-L 方程及对称性 $\delta L=0$，
$$\therefore \frac{d}{dt}\left(\frac{\partial L}{\partial\dot q_i}K_i\right)=0.$$

**物理/几何意义** 时间平移→能量，空间平移→动量，旋转→角动量。

**极限/特例** 守恒量降低系统维数。

**典型应用** 积分与约化。

#### 知识点7.8 场论的变分与 Noether 流

**定义** 场 $\phi(x^\mu)$，作用量 $S=\int\mathcal L(\phi,\partial_\mu\phi)d^4x$。

**假设** 拉格朗日密度 $\mathcal L$。

**公式与推导** $\delta S=0$ 给出
$$\partial_\mu\frac{\partial\mathcal L}{\partial(\partial_\mu\phi)}-\frac{\partial\mathcal L}{\partial\phi}=0.$$
Noether 流 $j^\mu=\frac{\partial\mathcal L}{\partial(\partial_\mu\phi)}\delta\phi-\mathcal L\delta x^\mu$ 满足 $\partial_\mu j^\mu=0$。

**物理/几何意义** 连续性方程形式。

**极限/特例** 电磁场、弹性动力学。

**典型应用** 连续介质场论。

### 第8章 有限元方法

#### 知识点8.1 弱形式与函数空间

**定义** 对 $-\frac{d}{dx}(a\frac{du}{dx})=f$，乘以试探函数 $v$ 并积分。

**假设** $u\in H^1,\ v\in H_0^1$。

**公式与推导** 分部积分
$$\int_0^L a u'v'\,dx=\int_0^L fv\,dx+\Big[a u'v\Big]_0^L.$$
由 $v(0)=v(L)=0$，得弱形式
$$a(u,v):=\int_0^L a u'v'\,dx=\int_0^L fv\,dx=:\ell(v).$$

**物理/几何意义** 降低导数阶，容许更弱的解。

**极限/特例** 自然边界自动满足。

**典型应用** FEM 起点。

#### 知识点8.2 Galerkin 离散与单元刚度矩阵

**定义** $u_h=\sum_j N_j u_j$，取 $v=N_i$。

**假设** 局部支撑基函数。

**公式与推导** 代入得
$$\sum_j\left(\int a N_i'N_j'\,dx\right)u_j=\int fN_i\,dx,$$
即 $Ku=f$。一维两节点线性单元，$N_1=1-\xi,\ N_2=\xi,\ \xi=x/h$，
$$K^{(e)}=\frac{a}{h}\begin{pmatrix}1&-1\\-1&1\end{pmatrix},\quad f^{(e)}=\frac{fh}{2}\begin{pmatrix}1\\1\end{pmatrix}.$$

**物理/几何意义** 单元矩阵表示局部能量。

**极限/特例** 高阶单元更精确。

**典型应用** 结构分析。

#### 知识点8.3 组装与全局矩阵

**定义** 将单元贡献叠加到全局系统。

**假设** 节点编号映射。

**公式与推导** 
$$K=\sum_e A^{(e)\top}K^{(e)}A^{(e)},\quad f=\sum_e A^{(e)\top}f^{(e)},$$
其中 $A^{(e)}$ 为自由度选择矩阵。
$\because$ 总能量为单元能量之和，
$$\therefore \text{叠加即组装}.$$

**物理/几何意义** 局部到全局的映射。

**极限/特例** 稀疏带状。

**典型应用** 大型结构。

#### 知识点8.4 一维杆单元完整推导

**定义** 轴向位移 $u(x)$，应变 $\epsilon=u'$，应力 $\sigma=E\epsilon$。

**假设** 小变形、线性弹性、等截面 $A$。

**公式与推导** 总势能
$$\Pi=\frac12\int_0^L EA(u')^2dx-\int_0^L bu\,dx-F_1u_1-F_2u_2.$$
离散 $u=N_1u_1+N_2u_2$，$u'=\frac{u_2-u_1}{h}$，
$$\therefore \Pi=\frac{EA}{2h}(u_2-u_1)^2-\frac{bh}{2}(u_1+u_2)-F_1u_1-F_2u_2.$$
由 $\partial\Pi/\partial u_i=0$ 得刚度矩阵。

**物理/几何意义** 杆的轴向刚度 $EA/h$。

**极限/特例** 多点单元。

**典型应用** 桁架、杆系。

#### 知识点8.5 二维三角形单元与形函数

**定义** 面积坐标 $L_1,L_2,L_3$，$L_1+L_2+L_3=1$。

**假设** 线性位移场。

**公式与推导** $u=\sum N_i u_i$，$N_i=L_i$。应变-位移矩阵
$$B=\frac{1}{2A}\begin{pmatrix}y_2-y_3&0&y_3-y_1&0&y_1-y_2&0\\0&x_3-x_2&0&x_1-x_3&0&x_2-x_1\\x_3-x_2&y_2-y_3&x_1-x_3&y_3-y_1&x_2-x_1&y_1-y_2\end{pmatrix}.$$
单元刚度
$$K^{(e)}=\int_{A_e}B^\top D B\,t\,dA=A_e t B^\top D B.$$

**物理/几何意义** 常应变单元。

**极限/特例** 面积 $A_e\to0$ 时退化。

**典型应用** 平面应力/应变。

#### 知识点8.6 等参单元与 Jacobian

**定义** 用同一形函数插值坐标与位移。

**假设** 映射 $x(\xi,\eta)=\sum N_i(\xi,\eta)x_i$ 可逆。

**公式与推导** 积分换元
$$\int_{\Omega_e}f\,dx\,dy=\int_{-1}^1\int_{-1}^1 f\,|J|\,d\xi\,d\eta,$$
$$J=\begin{pmatrix}\partial x/\partial\xi&\partial y/\partial\xi\\\partial x/\partial\eta&\partial y/\partial\eta\end{pmatrix}.$$
导数变换
$$\begin{pmatrix}\partial N/\partial x\\\partial N/\partial y\end{pmatrix}=J^{-1}\begin{pmatrix}\partial N/\partial\xi\\\partial N/\partial\eta\end{pmatrix}.$$

**物理/几何意义** 曲边单元映射为标准单元。

**极限/特例** $|J|>0$ 保证不翻转。

**典型应用** 复杂几何。

#### 知识点8.7 有限元误差估计

**定义** 能量范数 $\|u\|_E^2=\int a(u')^2dx$。

**假设** 精确解 $u\in H^{k+1}$。

**公式与推导** Céa 引理
$$\|u-u_h\|_E\le C\inf_{v_h}\|u-v_h\|_E.$$
对 $P_k$ 单元插值误差
$$\|u-\Pi_h u\|_{H^s}\le C h^{k+1-s}|u|_{H^{k+1}},$$
$$\therefore \|u-u_h\|_E=O(h^k),\ \|u-u_h\|_{L^2}=O(h^{k+1}).$$

**物理/几何意义** 网格加密按阶收敛。

**极限/特例** 奇点会降低收敛阶。

**典型应用** 自适应网格。

#### 知识点8.8 时间相关有限元：半离散

**定义** 对 $u_t=\nabla\cdot(a\nabla u)+f$ 空间 FEM。

**假设** 时间连续。

**公式与推导** 弱形式得
$$M\dot U+KU=F,$$
其中 $M_{ij}=\int N_iN_jdx$ 为质量矩阵。对时间用 $\theta$ 法
$$\left(\frac{M}{\Delta t}+\theta K\right)U^{n+1}=\left(\frac{M}{\Delta t}-(1-\theta)K\right)U^n+\theta F^{n+1}+(1-\theta)F^n.$$
$\theta=1$ 后向 Euler，$\theta=1/2$ Crank-Nicolson。

**物理/几何意义** 质量矩阵刻画惯性。

**极限/特例** 集中质量可显式。

**典型应用** 瞬态热/结构。

#### 知识点8.9 混合有限元与鞍点问题

**定义** 同时求解多个场，如 Stokes 的 $(u,p)$。

**假设** 满足 inf-sup 条件。

**公式与推导** 弱形式
$$\begin{cases}\int\nabla u:\nabla v-\int p\nabla\cdot v=\int fv,\\ \int q\nabla\cdot u=0.\end{cases}$$
离散系统为
$$\begin{pmatrix}A&B^\top\\B&0\end{pmatrix}\begin{pmatrix}U\\P\end{pmatrix}=\begin{pmatrix}F\\0\end{pmatrix},$$
稳定性需 $\inf_q\sup_v\frac{b(v,q)}{\|v\|\|q\|}\ge\beta>0$。

**物理/几何意义** 压力-速度耦合。

**极限/特例** Taylor-Hood 单元满足条件。

**典型应用** 不可压流动。

### 第9章 拉格朗日力学与哈密顿力学

#### 知识点9.1 广义坐标与约束

**定义** 用 $n$ 个独立广义坐标 $q_i$ 描述 $3N$ 质点系统，自由度 $n=3N-k$。

**假设** 完整约束 $\phi_\alpha(r,t)=0$。

**公式与推导** 虚位移满足
$$\sum_i\frac{\partial\phi_\alpha}{\partial r_i}\cdot\delta r_i=0.$$
理想约束虚功为零，从而可用广义坐标消除约束力。

**物理/几何意义** 约束流形上的坐标。

**极限/特例** 非完整约束需乘子。

**典型应用** 机构、多体系统。

#### 知识点9.2 拉格朗日方程的推导

**定义** $L=T-V$。

**假设** 完整、理想约束、保守力。

**公式与推导** 从达朗贝尔原理
$$\sum_i(m_i\ddot r_i-F_i)\cdot\delta r_i=0.$$
利用 $\dot r_i=\sum_j\frac{\partial r_i}{\partial q_j}\dot q_j$ 及恒等式
$$\frac{d}{dt}\frac{\partial r_i}{\partial q_j}=\frac{\partial\dot r_i}{\partial q_j},\quad \frac{\partial\dot r_i}{\partial\dot q_j}=\frac{\partial r_i}{\partial q_j},$$
可证
$$\sum_i m_i\ddot r_i\cdot\frac{\partial r_i}{\partial q_j}=\frac{d}{dt}\frac{\partial T}{\partial\dot q_j}-\frac{\partial T}{\partial q_j}.$$
$$\therefore \frac{d}{dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}=Q_j^{nc}.$$

**物理/几何意义** 广义力平衡。

**极限/特例** 直角坐标回到牛顿方程。

**典型应用** 复杂系统建模。

#### 知识点9.3 对称性与守恒量

**定义** 循环坐标：$L$ 不显含 $q_j$。

**假设** $L$ 正则。

**公式与推导** 若 $\partial L/\partial q_j=0$，则
$$\frac{d}{dt}\frac{\partial L}{\partial\dot q_j}=0\ \Rightarrow\ p_j=\frac{\partial L}{\partial\dot q_j}=const.$$
能量 $E=\sum_j\dot q_jp_j-L$ 守恒，若 $L$ 不显含 $t$：
$$\frac{dE}{dt}=-\frac{\partial L}{\partial t}=0.$$

**物理/几何意义** 对称性$\to$守恒律。

**极限/特例** 中心力场角动量守恒。

**典型应用** 约化问题。

#### 知识点9.4 哈密顿正则方程

**定义** $H(q,p,t)=\sum_i p_i\dot q_i-L$。

**假设** 变换 $(q,\dot q)\to(q,p)$ 可逆（Legendre 变换）。

**公式与推导** 全微分
$$dH=\sum_i\left(\dot q_i dp_i+p_id\dot q_i-\frac{\partial L}{\partial q_i}dq_i-\frac{\partial L}{\partial\dot q_i}d\dot q_i\right)-\frac{\partial L}{\partial t}dt.$$
由 $p_i=\partial L/\partial\dot q_i$ 消去 $d\dot q_i$，
$$\therefore \dot q_i=\frac{\partial H}{\partial p_i},\quad \dot p_i=-\frac{\partial H}{\partial q_i},\quad \frac{\partial H}{\partial t}=-\frac{\partial L}{\partial t}.$$

**物理/几何意义** 相空间一阶系统，流保持辛结构。

**极限/特例** 可积系统。

**典型应用** 辛积分、统计力学。

#### 知识点9.5 泊松括号与守恒量

**定义** $\{f,g\}=\sum_i\left(\frac{\partial f}{\partial q_i}\frac{\partial g}{\partial p_i}-\frac{\partial f}{\partial p_i}\frac{\partial g}{\partial q_i}\right)$。

**假设** 正则变量。

**公式与推导** 沿轨道
$$\frac{df}{dt}=\{f,H\}+\frac{\partial f}{\partial t}.$$
$\because$ 代入正则方程，
$$\therefore f \text{ 守恒}\Leftrightarrow \{f,H\}=0.$$
基本括号 $\{q_i,q_j\}=0,\ \{p_i,p_j\}=0,\ \{q_i,p_j\}=\delta_{ij}$。

**物理/几何意义** 括号描述流对函数的作用。

**极限/特例** 角动量分量括号给出李代数。

**典型应用** 量子化、可积性。

#### 知识点9.6 正则变换与生成函数

**定义** 保持正则方程形式的变换 $(q,p)\to(Q,P)$。

**假设** 存在生成函数 $F$。

**公式与推导** 要求 $p\,dq-H\,dt=P\,dQ-K\,dt+dF$。取 $F_2(q,P)$，
$$p=\frac{\partial F_2}{\partial q},\quad Q=\frac{\partial F_2}{\partial P},\quad K=H+\frac{\partial F_2}{\partial t}.$$

**物理/几何意义** 相空间保辛坐标变换。

**极限/特例** 恒等变换 $F_2=qP$。

**典型应用** Hamilton-Jacobi 理论。

#### 知识点9.7 Hamilton-Jacobi 方程

**定义** 找生成函数 $S(q,t)$ 使新哈密顿 $K=0$。

**假设** 可分离系统。

**公式与推导** 由 $p=\partial S/\partial q$ 及 $K=H+\partial S/\partial t=0$，
$$\therefore \frac{\partial S}{\partial t}+H\!\left(q,\frac{\partial S}{\partial q},t\right)=0.$$
若 $H=E$ 守恒，$S=W(q)-Et$，得定态方程 $H(q,\partial W/\partial q)=E$。

**物理/几何意义** 波前演化，连接经典与量子。

**极限/特例** 可积系统分离变量。

**典型应用** 解析力学、WKB 近似。

#### 知识点9.8 刚体动力学欧拉方程

**定义** 主惯量 $I_1,I_2,I_3$，角速度 $\omega$。

**假设** 无外力矩。

**公式与推导** 拉格朗日 $L=\frac12(I_1\omega_1^2+I_2\omega_2^2+I_3\omega_3^2)$，角速度与欧拉角关系代入 E-L 方程，得
$$I_1\dot\omega_1=(I_2-I_3)\omega_2\omega_3,$$
$$I_2\dot\omega_2=(I_3-I_1)\omega_3\omega_1,\quad I_3\dot\omega_3=(I_1-I_2)\omega_1\omega_2.$$

**物理/几何意义** 角动量守恒下的进动/章动。

**极限/特例** 对称陀螺 $I_1=I_2$。

**典型应用** 卫星姿态。

### 第10章 辛几何与几何数值积分

#### 知识点10.1 辛结构与哈密顿流

**定义** 相空间 $z=(q,p)$，辛形式 $\omega=\sum_i dq_i\wedge dp_i$，矩阵 $\Omega=\begin{pmatrix}0&I\\-I&0\end{pmatrix}$。

**假设** 正则系统。

**公式与推导** 运动方程 $\dot z=\Omega\nabla H$。流映射 $\Phi_t$ 满足
$$\left(\frac{\partial\Phi_t}{\partial z}\right)^\top\Omega\left(\frac{\partial\Phi_t}{\partial z}\right)=\Omega,$$
即保辛。由 $\frac{d}{dt}(M^\top\Omega M)=0$ 且初值 $M(0)=I$ 得证，其中 $M=\partial\Phi_t/\partial z$。

**物理/几何意义** 相面积/体积守恒（Liouville）。

**极限/特例** 一自由度面积守恒。

**典型应用** 长期积分。

#### 知识点10.2 辛 Euler 与隐式中点

**定义** 辛 Euler：
$$p_{n+1}=p_n-h\nabla_qH(q_n,p_{n+1}),\quad q_{n+1}=q_n+h\nabla_pH(q_n,p_{n+1}).$$

**假设** $H$ 可微。

**公式与推导** 其 Jacobian 满足辛条件，因此保持辛结构但不显式。

**物理/几何意义** 一阶辛格式。

**极限/特例** 对线性系统能量有界振荡。

**典型应用** 需要守恒的长时间模拟。

#### 知识点10.3 Störmer-Verlet 格式与分裂

**定义** 将 $H=T(p)+V(q)$ 分裂，交替精确求解。

**假设** $T,V$ 分别可解。

**公式与推导** 
$$\Phi_h=\Phi_V^{h/2}\circ\Phi_T^{h}\circ\Phi_V^{h/2}+O(h^3).$$
展开得
$$p_{n+1/2}=p_n-\frac h2\nabla V(q_n),\quad q_{n+1}=q_n+h\nabla T(p_{n+1/2}),\quad p_{n+1}=p_{n+1/2}-\frac h2\nabla V(q_{n+1}).$$
每步均为辛映射，复合仍辛。

**物理/几何意义** 二阶精度、时间可逆、能量长期有界。

**极限/特例** 谐振子能量精确有界。

**典型应用** 分子动力学。

#### 知识点10.4 后向误差与影子哈密顿量

**定义** 辛积分的数值轨道精确对应某个修正哈密顿 $\tilde H$ 的轨道。

**假设** 格式为辛且可逆。

**公式与推导** 存在 $\tilde H=H+h^2H_2+\cdots$ 使
$$\Phi_h=\Phi_h^{\tilde H}.$$
$$\therefore \text{数值能量 }\tilde H \text{ 守恒，原能量有界振荡}.$$

**物理/几何意义** 长期稳定性来源。

**极限/特例** 振荡幅度 $O(h^2)$。

**典型应用** 守恒性证明。

#### 知识点10.5 变分积分器

**定义** 离散作用量 $S_d=\sum_k L_d(q_k,q_{k+1})$，令 $\delta S_d=0$。

**假设** $L_d$ 逼近 $\int Ldt$。

**公式与推导** 离散 E-L：
$$D_2L_d(q_{k-1},q_k)+D_1L_d(q_k,q_{k+1})=0.$$
对 $L_d=hL(\frac{q_k+q_{k+1}}2,\frac{q_{k+1}-q_k}{h})$ 得 Verlet。

**物理/几何意义** 保持离散变分结构，自动辛。

**极限/特例** 高阶变分积分器。

**典型应用** 约束系统。

#### 知识点10.6 保体积与李群积分

**定义** 对 $\dot z=f(z)$，若 $\nabla\cdot f=0$ 则流保体积。

**假设** 可分/李群结构。

**公式与推导** 分裂法每一子流保体积，复合保体积。

**物理/几何意义** 不可压流、哈密顿系统。

**极限/特例** 辛$\Rightarrow$保体积。

**典型应用** 不可压流体、刚体。

### 第11章 连续介质力学与 Navier-Stokes 方程

#### 知识点11.1 变形梯度与应变张量

**定义** 变形梯度 $F=\partial x/\partial X$。

**假设** 连续、可微、$\det F>0$。

**公式与推导** 右 Cauchy-Green 张量 $C=F^\top F$，Green-Lagrange 应变
$$E=\frac12(C-I)=\frac12(F^\top F-I).$$
小变形 $F=I+\nabla u$，则
$$E\approx\frac12(\nabla u+\nabla u^\top)=\epsilon.$$

**物理/几何意义** 长度平方变化 $ds^2-dS^2=2dX^\top E\,dX$。

**极限/特例** 刚体运动 $E=0$。

**典型应用** 大变形/小变形理论。

#### 知识点11.2 应力张量与平衡方程

**定义** Cauchy 应力 $t=n\cdot\sigma$。

**假设** 连续、动量守恒。

**公式与推导** 对控制体用动量定理
$$\int_\Omega\rho\dot v\,dV=\int_{\partial\Omega}t\,dS+\int_\Omega\rho b\,dV.$$
由散度定理 $\int_{\partial\Omega}\sigma n\,dS=\int_\Omega\nabla\cdot\sigma\,dV$，
$$\therefore \rho\dot v=\nabla\cdot\sigma+\rho b.$$
角动量守恒给出 $\sigma=\sigma^\top$。

**物理/几何意义** 牛顿第二定律的连续形式。

**极限/特例** 静力平衡 $\nabla\cdot\sigma+\rho b=0$。

**典型应用** 结构/流体。

#### 知识点11.3 线性弹性本构与 Lamé 参数

**定义** 广义 Hooke 定律 $\sigma=\mathbb C:\epsilon$。

**假设** 各向同性、线弹性。

**公式与推导** 由各向同性张量表示
$$\sigma=2\mu\epsilon+\lambda\text{tr}(\epsilon)I.$$
反解
$$\epsilon=\frac{1+\nu}{E}\sigma-\frac{\nu}{E}\text{tr}(\sigma)I,\quad \mu=\frac{E}{2(1+\nu)},\ \lambda=\frac{E\nu}{(1+\nu)(1-2\nu)}.$$

**物理/几何意义** 拉压与剪切解耦。

**极限/特例** $\nu\to1/2$ 不可压。

**典型应用** 弹性力学。

#### 知识点11.4 Navier-Stokes 方程推导

**定义** 牛顿流体本构
$$\sigma=-pI+2\mu D+\lambda_v(\nabla\cdot v)I,\quad D=\frac12(\nabla v+\nabla v^\top).$$

**假设** Stokes 假设 $\lambda_v+2\mu/3=0$，不可压 $\nabla\cdot v=0$。

**公式与推导** 代入动量方程
$$\rho\left(\frac{\partial v}{\partial t}+v\cdot\nabla v\right)=\nabla\cdot\sigma+\rho b.$$
对不可压牛顿流体
$$\nabla\cdot(2\mu D)=2\mu\nabla\cdot D=\mu\nabla^2v+\mu\nabla(\nabla\cdot v)=\mu\nabla^2v.$$
$$\therefore \rho\left(\frac{\partial v}{\partial t}+v\cdot\nabla v\right)=-\nabla p+\mu\nabla^2v+\rho b.$$

**物理/几何意义** 惯性 = 压力梯度 + 黏性扩散 + 体积力。

**极限/特例** $\mu=0$ Euler 方程。

**典型应用** CFD、湍流。

#### 知识点11.5 无量纲化与 Reynolds 数

**定义** 用特征速度 $U$、长度 $L$ 无量纲化。

**假设** 常物性。

**公式与推导** 令 $v=U\tilde v,\ x=L\tilde x,\ t=(L/U)\tilde t,\ p=\rho U^2\tilde p$，
$$\frac{\partial \tilde v}{\partial\tilde t}+\tilde v\cdot\tilde\nabla\tilde v=-\tilde\nabla\tilde p+\frac{1}{Re}\tilde\nabla^2\tilde v+\tilde b,$$
$$Re=\frac{\rho UL}{\mu}.$$

**物理/几何意义** 惯性力/黏性力之比。

**极限/特例** $Re\ll1$ Stokes 流，$Re\gg1$ 边界层/湍流。

**典型应用** 相似实验。

#### 知识点11.6 涡量与涡量方程

**定义** $\omega=\nabla\times v$。

**假设** 不可压均匀流体。

**公式与推导** 对 N-S 取旋度，利用 $\nabla\times(\nabla p)=0$、$\nabla\times(\nabla^2v)=\nabla^2\omega$，
$$\frac{\partial\omega}{\partial t}+\nabla\times(\omega\times v)=\nu\nabla^2\omega.$$
利用 $\nabla\times(\omega\times v)=v\cdot\nabla\omega-\omega\cdot\nabla v$，
$$\therefore \frac{D\omega}{Dt}=\omega\cdot\nabla v+\nu\nabla^2\omega.$$

**物理/几何意义** 拉伸涡管放大涡量，黏性扩散。

**极限/特例** 二维 $\omega\cdot\nabla v=0$。

**典型应用** 旋涡动力学。

#### 知识点11.7 能量方程与耗散

**定义** 内能 $e$，热流 $q=-k\nabla T$。

**假设** Fourier 定律。

**公式与推导** 能量守恒
$$\rho\frac{De}{Dt}=\sigma:\nabla v-\nabla\cdot q+\rho r.$$
对牛顿流体 $\sigma:\nabla v=-p\nabla\cdot v+\Phi$，耗散函数
$$\Phi=2\mu D:D+\lambda_v(\nabla\cdot v)^2\ge0.$$

**物理/几何意义** 机械能不可逆转为热。

**极限/特例** 不可压 $\Phi=2\mu D:D$。

**典型应用** 传热/热力学耦合。

#### 知识点11.8 边界层方程

**定义** 高 $Re$ 下近壁薄层。

**假设** 厚度 $\delta\ll L$，$\delta\sim L/\sqrt{Re}$。

**公式与推导** 量级分析保留主导项得
$$u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=-\frac{1}{\rho}\frac{dp}{dx}+\nu\frac{\partial^2 u}{\partial y^2},\quad \frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0.$$
边界条件 $u=v=0$（壁），$u\to U_e$（外缘）。

**物理/几何意义** 黏性仅作用于薄层。

**极限/特例** Blasius 相似解 $f'''+\frac12ff''=0$。

**典型应用** 阻力预测。

#### 知识点11.9 湍流与 Reynolds 平均

**定义** $v=\bar v+v',\ p=\bar p+p'$。

**假设** 统计平均。

**公式与推导** 代入 N-S 取平均，非线性项生成 Reynolds 应力
$$\rho\frac{D\bar v}{Dt}=-\nabla\bar p+\mu\nabla^2\bar v-\rho\nabla\cdot\overline{v'v'}.$$
湍动能方程含产生项 $-\rho\overline{v'v'}:\nabla\bar v$ 与耗散 $\rho\nu\overline{(\nabla v')^2}$。

**物理/几何意义** 湍流脉动输运动量。

**极限/特例** 需湍流模型闭合。

**典型应用** 工程湍流计算。

### 第12章 本构关系与材料模型

#### 知识点12.1 弹性本构的热力学约束

**定义** 应变能密度 $W(\epsilon)$，$\sigma=\partial W/\partial\epsilon$。

**假设** 超弹性、等温。

**公式与推导** 由热力学第二定律，Clausius-Duhem 不等式
$$\sigma:\dot\epsilon-\dot W\ge0.$$
取 $\dot W=\frac{\partial W}{\partial\epsilon}:\dot\epsilon$，
$$\therefore \left(\sigma-\frac{\partial W}{\partial\epsilon}\right):\dot\epsilon=0\ \Rightarrow\ \sigma=\frac{\partial W}{\partial\epsilon}.$$

**物理/几何意义** 应力由势函数导出。

**极限/特例** 线弹性 $W=\frac12\epsilon:\mathbb C:\epsilon$。

**典型应用** 材料稳定性（凸性）。

#### 知识点12.2 各向同性弹性张量

**定义** 四阶张量 $\mathbb C$ 满足各向同性。

**假设** 旋转不变。

**公式与推导** 由各向同性张量基
$$\mathbb C_{ijkl}=\lambda\delta_{ij}\delta_{kl}+\mu(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}).$$
$$\therefore \sigma=2\mu\epsilon+\lambda\text{tr}\epsilon I.$$

**物理/几何意义** 体积模量 $K=\lambda+\frac23\mu$，剪切模量 $\mu$。

**极限/特例** $\mu=0$ 无剪切刚度。

**典型应用** 弹性波速 $c_p,c_s$。

#### 知识点12.3 弹塑性本构与屈服准则

**定义** 屈服面 $f(\sigma,\kappa)\le0$。

**假设** 关联流动法则。

**公式与推导** 塑性应变率
$$\dot\epsilon^p=\dot\lambda\frac{\partial f}{\partial\sigma},$$
一致性条件 $\dot f=0$：
$$\dot f=\frac{\partial f}{\partial\sigma}:\dot\sigma+\frac{\partial f}{\partial\kappa}\dot\kappa=0,$$
给出塑性乘子 $\dot\lambda$。von Mises 屈服
$$f=\sqrt{\frac32 s:s}-\sigma_y,\quad s=\sigma-\frac13\text{tr}\sigma I.$$

**物理/几何意义** 应力点保持在屈服面上。

**极限/特例** 理想塑性 $\partial f/\partial\kappa=0$。

**典型应用** 金属成形。

#### 知识点12.4 粘弹性：Maxwell 与 Kelvin-Voigt

**定义** 弹簧与黏壶组合。

**假设** 线性粘弹性。

**公式与推导** Maxwell：$\dot\epsilon=\frac{\dot\sigma}{E}+\frac{\sigma}{\eta}$，松弛
$$\sigma(t)=\sigma_0e^{-Et/\eta},\quad \tau=\eta/E.$$
Kelvin-Voigt：$\sigma=E\epsilon+\eta\dot\epsilon$，蠕变
$$\epsilon(t)=\frac{\sigma_0}{E}\left(1-e^{-Et/\eta}\right).$$

**物理/几何意义** 松弛与蠕变。

**极限/特例** 标准线性固体可同时描述。

**典型应用** 聚合物、生物组织。

#### 知识点12.5 塑性增量与一致性切线

**定义** 返回映射算法。

**假设** 小应变、关联塑性。

**公式与推导** 弹性预测
$$\sigma^{tr}=\mathbb C:(\epsilon^{n+1}-\epsilon^{p,n}).$$
若 $f(\sigma^{tr})>0$，修正
$$\sigma^{n+1}=\sigma^{tr}-\Delta\lambda\mathbb C:\frac{\partial f}{\partial\sigma}.$$
由一致性条件解 $\Delta\lambda$。

**物理/几何意义** 先弹性试探再塑性投影。

**极限/特例** 径向返回用于 von Mises。

**典型应用** 非线性有限元。

#### 知识点12.6 超弹性大变形模型

**定义** 应变能 $W(\lambda_1,\lambda_2,\lambda_3)$。

**假设** 各向同性、不可压。

**公式与推导** 由 $\sigma_i=\frac{1}{J}\lambda_i\frac{\partial W}{\partial\lambda_i}$，不可压 $\lambda_1\lambda_2\lambda_3=1$。Neo-Hookean
$$W=\frac{\mu}{2}(I_1-3),\quad I_1=\lambda_1^2+\lambda_2^2+\lambda_3^2,$$
$$\therefore \sigma=2\mu B-pI,\quad B=FF^\top.$$

**物理/几何意义** 熵弹性主导橡胶。

**极限/特例** Mooney-Rivlin 增加 $I_2$ 项。

**典型应用** 橡胶、软材料。

#### 知识点12.7 损伤与断裂

**定义** 损伤变量 $d\in[0,1]$，有效应力 $\tilde\sigma=\sigma/(1-d)$。

**假设** 各向同性损伤。

**公式与推导** 应变能
$$W=\frac12(1-d)\epsilon:E\epsilon,$$
$$\therefore \sigma=(1-d)E\epsilon,\quad Y=-\frac{\partial W}{\partial d}=\frac12\epsilon:E\epsilon.$$
损伤演化 $\dot d=\dot\lambda\frac{\partial F}{\partial Y}$。

**物理/几何意义** 刚度退化。

**极限/特例** $d=1$ 完全破坏。

**典型应用** 混凝土、复合材料。

### 第13章 稳定性与分岔

#### 知识点13.1 平衡点与线性稳定性

**定义** 系统 $\dot x=f(x)$，平衡点 $f(x^*)=0$。

**假设** $f$ 可微。

**公式与推导** 令 $x=x^*+\delta x$，
$$\dot{\delta x}=J\delta x+O(\|\delta x\|^2),\quad J=\frac{\partial f}{\partial x}\Big|_{x^*}.$$
$$\therefore \text{Re}\lambda_i(J)<0\ \forall i\Rightarrow \text{渐近稳定}.$$

**物理/几何意义** 小扰动指数衰减/增长。

**极限/特例** 纯虚特征值需非线性判断。

**典型应用** 结构屈曲、流稳。

#### 知识点13.2 Lyapunov 稳定性与函数

**定义** 若存在 $V(x)>0$ 且 $\dot V\le0$，则稳定；$\dot V<0$ 渐近稳定。

**假设** $V$ 正定、连续可微。

**公式与推导** 沿轨道
$$\dot V=\nabla V\cdot f(x).$$
$$\because V \text{ 单调不增下方有界},\quad \therefore \text{收敛}.$$
LaSalle 不变原理可加强结论。

**物理/几何意义** 能量函数单调下降。

**极限/特例** 线性系统 $V=x^\top Px$，$A^\top P+PA=-Q$。

**典型应用** 非线性稳定分析。

#### 知识点13.3 分岔类型

**定义** 参数变化时平衡点数目/稳定性改变。

**假设** 单参数族。

**公式与推导** 鞍结分岔 $\dot x=\mu-x^2$，平衡 $x=\pm\sqrt\mu$，$\mu=0$ 合并。
跨临界 $\dot x=\mu x-x^2$。
叉式 $\dot x=\mu x-x^3$，$\mu>0$ 三个平衡。

**物理/几何意义** 解的拓扑改变。

**极限/特例** Hopf 分岔产生极限环。

**典型应用** 屈曲、颤振。

#### 知识点13.4 结构屈曲与特征值问题

**定义** 压杆平衡分支。

**假设** 小变形、轴向压力 $P$。

**公式与推导** 梁方程
$$EIw''''+Pw''=0,$$
设 $w\sim\sin(n\pi x/L)$ 得
$$P_{cr}=\frac{n^2\pi^2EI}{L^2}.$$
最低模态 $n=1$。

**物理/几何意义** 压力达临界值失去稳定性。

**极限/特例** 边界条件改变 $P_{cr}$。

**典型应用** 结构设计。

#### 知识点13.5 线性稳定性与 Orr-Sommerfeld

**定义** 平行剪切流扰动。

**假设** 模态扰动 $v'\sim e^{i(kx-\omega t)}$。

**公式与推导** 线性化 N-S 得 Orr-Sommerfeld
$$\left[(\partial_y^2-k^2)^2-i k Re(U-c)(\partial_y^2-k^2)+ikReU''\right]\hat v=0.$$
$c=\omega/k$ 为相速度，$c_i>0$ 不稳定。

**物理/几何意义** 黏性与惯性竞争。

**极限/特例** 无黏 Rayleigh 方程。

**典型应用** 边界层转捩。

### 第14章 机器学习基础

#### 知识点14.1 线性回归与最小二乘

**定义** 模型 $y=Xw+\epsilon$。

**假设** $\epsilon$ 零均值、同方差。

**公式与推导** 损失
$$J(w)=\frac12\|Xw-y\|^2.$$
$$\because \nabla_wJ=X^\top(Xw-y)=0,\quad \therefore w=(X^\top X)^{-1}X^\top y.$$
在概率视角，若 $\epsilon\sim\mathcal N(0,\sigma^2)$，
$$p(y|X,w)\propto\exp\left(-\frac{\|Xw-y\|^2}{2\sigma^2}\right),$$
极大似然等价最小二乘。

**物理/几何意义** 正交投影。

**极限/特例** 岭回归 $w=(X^\top X+\lambda I)^{-1}X^\top y$ 处理共线性。

**典型应用** 拟合。

#### 知识点14.2 逻辑回归与交叉熵

**定义** $p(y=1|x)=\sigma(w^\top x)=\frac{1}{1+e^{-w^\top x}}$。

**假设** 伯努利似然。

**公式与推导** 对数似然
$$\ell(w)=\sum_i[y_i\ln\sigma_i+(1-y_i)\ln(1-\sigma_i)].$$
梯度
$$\because \frac{\partial\sigma}{\partial z}=\sigma(1-\sigma),\quad \therefore \nabla_w\ell=X^\top(y-\sigma).$$

**物理/几何意义** 线性决策边界，概率输出。

**极限/特例** 系数可加正则。

**典型应用** 分类。

#### 知识点14.3 感知机与多层网络

**定义** 单层 $\hat y=\phi(Wx+b)$。

**假设** 激活 $\phi$ 非线性。

**公式与推导** 多层前向
$$z^{(l)}=W^{(l)}a^{(l-1)}+b^{(l)},\quad a^{(l)}=\phi(z^{(l)}).$$
万能逼近定理：紧集上连续函数可被单隐层网络一致逼近。

**物理/几何意义** 逐层特征变换。

**极限/特例** 无隐层退化为线性。

**典型应用** 非线性回归。

#### 知识点14.4 损失函数与极大似然

**定义** 回归用 MSE，分类用交叉熵。

**假设** 独立同分布样本。

**公式与推导** 高斯噪声下负对数似然
$$-\ln p=-\frac{N}{2}\ln(2\pi\sigma^2)+\frac{1}{2\sigma^2}\sum(y_i-\hat y_i)^2,$$
$$\therefore \text{最小化 MSE}.$$
交叉熵来自伯努利/Categorical 分布。

**物理/几何意义** 损失是负对数似然。

**极限/特例** 稳健回归用 Huber。

**典型应用** 训练目标设计。

#### 知识点14.5 过拟合、偏差-方差分解

**定义** 泛化误差期望
$$\mathbb E[(y-\hat f)^2]=\text{Bias}^2+\text{Var}+\sigma^2.$$

**假设** $y=f+\epsilon$。

**公式与推导** 令 $\bar f=\mathbb E[\hat f]$，
$$\because \mathbb E[(y-\hat f)^2]=\mathbb E[(f-\bar f)^2]+(\bar f-f)^2+\sigma^2,$$
$$\therefore \text{分解成立}.$$

**物理/几何意义** 模型复杂度权衡。

**极限/特例** 正则化/早停降低方差。

**典型应用** 模型选择。

#### 知识点14.6 卷积神经网络基础

**定义** 卷积层 $(f*k)_{ij}=\sum_{m,n}f_{i+m,j+n}k_{mn}$。

**假设** 局部连接、权值共享。

**公式与推导** 输出尺寸
$$H_{out}=\left\lfloor\frac{H_{in}+2P-K}{S}\right\rfloor+1.$$
平移等变形来自卷积与平移算子交换。

**物理/几何意义** 抽取局部模式。

**极限/特例** 池化降采样。

**典型应用** 图像/场数据。

#### 知识点14.7 循环网络与时间序列

**定义** $h_t=\phi(W_hh_{t-1}+W_xx_t+b)$。

**假设** 序列依赖。

**公式与推导** 展开时间后梯度含连乘 $\prod_t\frac{\partial h_t}{\partial h_{t-1}}$，导致梯度消失/爆炸。LSTM 用门控：
$$f_t=\sigma(W_f[h_{t-1},x_t]),\ i_t=\sigma(\cdot),\ o_t=\sigma(\cdot),$$
$$c_t=f_t\odot c_{t-1}+i_t\odot\tanh(\cdot).$$

**物理/几何意义** 记忆单元缓解长期依赖。

**极限/特例** 梯度裁剪。

**典型应用** 动力系统预测。

### 第15章 反向传播与自动微分

#### 知识点15.1 计算图与链式法则

**定义** 将复合函数表示为有向无环图。

**假设** 每个节点可微。

**公式与推导** 设 $y=f(g(x))$，
$$\frac{dy}{dx}=\frac{dy}{dg}\frac{dg}{dx}.$$
对多变量，节点 $i$ 接收上游梯度 $\bar a_i=\partial L/\partial a_i$，向下游传播
$$\bar a_j=\sum_{i\in\text{children}(j)}\bar a_i\frac{\partial a_i}{\partial a_j}.$$

**物理/几何意义** 局部导数沿图回传。

**极限/特例** 反向模式当输出少时高效。

**典型应用** 深度学习框架。

#### 知识点15.2 反向传播完整推导

**定义** 损失 $L$，网络 $a^{(l)}=\phi(z^{(l)}),\ z^{(l)}=W^{(l)}a^{(l-1)}+b^{(l)}$。

**假设** 标量损失。

**公式与推导** 定义误差 $\delta^{(l)}=\partial L/\partial z^{(l)}$。输出层
$$\delta^{(L)}=\nabla_{a^{(L)}}L\odot\phi'(z^{(L)}).$$
回传
$$\because z^{(l+1)}=W^{(l+1)}a^{(l)}+b^{(l+1)},\quad \therefore \delta^{(l)}=(W^{(l+1)\top}\delta^{(l+1)})\odot\phi'(z^{(l)}).$$
参数梯度
$$\frac{\partial L}{\partial W^{(l)}}=\delta^{(l)}a^{(l-1)\top},\quad \frac{\partial L}{\partial b^{(l)}}=\delta^{(l)}.$$

**物理/几何意义** 误差逐层反向加权。

**极限/特例** 激活饱和导致梯度消失。

**典型应用** 训练。

#### 知识点15.3 前向模式与反向模式自动微分

**定义** 前向模式计算 $(\dot x,\dot y)$；反向模式计算伴随。

**假设** 计算图可微。

**公式与推导** 前向：
$$\dot y=\frac{\partial f}{\partial x}\dot x,$$
一次传播一个输入方向。反向：
$$\bar x=\bar y\frac{\partial y}{\partial x},$$
一次传播一个输出方向。
$$\therefore \text{输入维大输出维小时用反向}.$$

**物理/几何意义** Jacobian-向量积与向量-Jacobian 积。

**极限/特例** 混合模式。

**典型应用** 梯度计算、敏感性分析。

#### 知识点15.4 通过 ODE/PDE 求解器微分

**定义** 离散求解器为前向图，直接反传。

**假设** 存储中间状态。

**公式与推导** 伴随法：对 $\dot z=f(z,t,\theta)$，定义伴随 $\lambda=\partial L/\partial z$，
$$\dot\lambda=-\lambda^\top\frac{\partial f}{\partial z},\quad \frac{dL}{d\theta}=-\int\lambda^\top\frac{\partial f}{\partial\theta}dt.$$
$$\therefore \text{无需存储全部时间步}.$$

**物理/几何意义** 连续伴随。

**极限/特例** 检查点折中内存与计算。

**典型应用** Neural ODE、PINN。

#### 知识点15.5 梯度检查与数值验证

**定义** 用有限差分校验解析梯度。

**假设** 函数可微。

**公式与推导** 中心差商
$$g_i\approx\frac{L(\theta+\epsilon e_i)-L(\theta-\epsilon e_i)}{2\epsilon},$$
相对误差
$$\text{rel}=\frac{|g_i-\hat g_i|}{\max(|g_i|,|\hat g_i|)}.$$

**物理/几何意义** 验证实现正确性。

**极限/特例** $\epsilon\sim10^{-6}$。

**典型应用** 调试。

### 第16章 优化算法

#### 知识点16.1 凸集、凸函数与最优性

**定义** $f$ 凸：$f(\lambda x+(1-\lambda)y)\le\lambda f(x)+(1-\lambda)f(y)$。

**假设** 可微。

**公式与推导** 一阶条件
$$f(y)\ge f(x)+\nabla f(x)^\top(y-x).$$
$$\therefore x^* \text{ 全局最优}\Leftrightarrow \nabla f(x^*)=0.$$
二阶条件 $\nabla^2f\succeq0$。

**物理/几何意义** 局部信息全局有效。

**极限/特例** 强凸保证唯一解。

**典型应用** 机器学习训练。

#### 知识点16.2 梯度下降与收敛性

**定义** $\theta_{k+1}=\theta_k-\eta\nabla f(\theta_k)$。

**假设** $f$ 为 $L$-光滑且 $\mu$-强凸。

**公式与推导** 由光滑性
$$f(\theta_{k+1})\le f(\theta_k)-\eta\|\nabla f\|^2+\frac{L\eta^2}{2}\|\nabla f\|^2.$$
取 $\eta=1/L$ 得单调下降。强凸下
$$\|\theta_k-\theta^*\|\le\left(1-\frac{\mu}{L}\right)^k\|\theta_0-\theta^*\|,$$
$$\therefore O\left(\frac{L}{\mu}\ln\frac1\epsilon\right)\ \text{次迭代}.$$

**物理/几何意义** 沿负梯度下降，条件数控制速度。

**极限/特例** $\eta>2/L$ 发散。

**典型应用** 大规模训练。

#### 知识点16.3 随机梯度下降（SGD）

**定义** 用无偏估计 $g_k=\nabla f_{i_k}(\theta_k)$，$\mathbb E[g_k]=\nabla f$。

**假设** 有界方差 $\mathbb E\|g-\nabla f\|^2\le\sigma^2$。

**公式与推导** 
$$\mathbb E[f(\theta_{k+1})]\le f(\theta_k)-\eta\|\nabla f\|^2+\frac{L\eta^2}{2}\sigma^2.$$
取 $\eta\propto1/\sqrt k$ 得 $O(1/\sqrt k)$ 收敛。

**物理/几何意义** 噪声换速度。

**极限/特例** mini-batch 降方差。

**典型应用** 深度学习。

#### 知识点16.4 动量与 Nesterov 加速

**定义** 动量：$v_{k+1}=\beta v_k+\nabla f$，$\theta_{k+1}=\theta_k-\eta v_{k+1}$。

**假设** 凸光滑。

**公式与推导** 重球法在 $\beta=1-\sqrt{\mu/L}$ 下达到 $O(\sqrt{L/\mu}\ln(1/\epsilon))$，优于梯度下降。
Nesterov 加速梯度满足
$$f(\theta_k)-f^*\le\frac{2L\|\theta_0-\theta^*\|^2}{(k+1)^2}.$$

**物理/几何意义** 惯性抑制振荡、提前观测。

**极限/特例** 非凸下为启发式。

**典型应用** 加速训练。

#### 知识点16.5 牛顿法与拟牛顿法

**定义** $\theta_{k+1}=\theta_k-H^{-1}\nabla f$。

**假设** $H=\nabla^2f$ 正定。

**公式与推导** 由二阶 Taylor
$$f(\theta+\Delta)\approx f(\theta)+\nabla f^\top\Delta+\frac12\Delta^\top H\Delta,$$
令导数为零得 $\Delta=-H^{-1}\nabla f$。
拟牛顿 BFGS 用 $B_k$ 逼近 $H$，满足割线方程
$$B_{k+1}s_k=y_k,\quad s_k=\theta_{k+1}-\theta_k,\ y_k=\nabla f_{k+1}-\nabla f_k.$$

**物理/几何意义** 二次模型步，曲率信息加速。

**极限/特例** $H$ 非正定时需修正。

**典型应用** 中小规模、二阶优化。

#### 知识点16.6 Adam 与自适应学习率

**定义** 
$$m_k=\beta_1m_{k-1}+(1-\beta_1)g_k,\quad v_k=\beta_2v_{k-1}+(1-\beta_2)g_k^2.$$

**假设** 指数滑动平均。

**公式与推导** 偏差修正
$$\hat m_k=\frac{m_k}{1-\beta_1^k},\quad \hat v_k=\frac{v_k}{1-\beta_2^k},$$
$$\therefore \theta_{k+1}=\theta_k-\eta\frac{\hat m_k}{\sqrt{\hat v_k}+\epsilon}.$$

**物理/几何意义** 每个参数步长按梯度尺度自适应。

**极限/特例** 非凸下收敛到驻点。

**典型应用** 深度网络默认优化器。

#### 知识点16.7 约束优化与拉格朗日乘子

**定义** $\min f(x)$ s.t. $g_i(x)=0,\ h_j(x)\le0$。

**假设** 约束规范成立。

**公式与推导** 构造
$$\mathcal L(x,\lambda,\mu)=f+\sum\lambda_ig_i+\sum\mu_jh_j.$$
KKT 条件
$$\nabla_x\mathcal L=0,\quad g_i=0,\quad h_j\le0,\quad \mu_j\ge0,\quad \mu_jh_j=0.$$
$\because$ 可行下降方向须不违反约束，
$$\therefore \text{KKT 为最优性必要条件}.$$

**物理/几何意义** 梯度被约束梯度线性组合抵消。

**极限/特例** 凸问题 KKT 充分。

**典型应用** 支持向量机、最优控制。

### 第17章 科学机器学习与物理信息方法

#### 知识点17.1 物理信息神经网络（PINN）

**定义** 用网络 $u_\theta(x,t)$ 近似 PDE 解。

**假设** PDE $N[u]=0$，初边值条件。

**公式与推导** 损失
$$L=\frac1{N_f}\sum_i\|N[u_\theta](x_i,t_i)\|^2+\frac1{N_b}\sum_j\|B[u_\theta]\|^2.$$
通过自动微分计算 $u_t,u_x,u_{xx}$，梯度下降优化 $\theta$。

**物理/几何意义** 软约束满足物理规律。

**极限/特例** 刚性/多尺度问题训练困难。

**典型应用** 正反问题、参数辨识。

#### 知识点17.2 本征正交分解（POD）

**定义** 对快照矩阵 $X=[x_1,\dots,x_m]$ 求最优低维基。

**假设** 能量意义下最优。

**公式与推导** 最小化投影误差
$$\min_{\Phi^\top\Phi=I}\sum_i\|x_i-\Phi\Phi^\top x_i\|^2.$$
由 SVD $X=U\Sigma V^\top$，最优 $\Phi$ 为前 $r$ 个左奇异向量，
$$\therefore \Phi=U_r,\quad x\approx\Phi a.$$

**物理/几何意义** 提取主导模态。

**极限/特例** 能量占比 $\sum_{i\le r}\sigma_i^2/\sum\sigma_i^2$。

**典型应用** 降阶模型。

#### 知识点17.3 Galerkin 降阶模型

**定义** $u\approx\Phi a(t)$ 代入 PDE 并投影。

**假设** 基 $\Phi$ 固定。

**公式与推导** 
$$\dot a=\Phi^\top f(\Phi a).$$
对线性系统 $\dot u=Au$，
$$\dot a=\Phi^\top A\Phi\,a=:A_r a.$$
若 $\Phi$ 为 POD 基，$A_r$ 低维。

**物理/几何意义** 快照子空间动力学。

**极限/特例** 非线性项需超缩减。

**典型应用** 实时仿真。

#### 知识点17.4 算子学习与 DeepONet

**定义** 学习映射 $G:u\mapsto G(u)$。

**假设** 有监督数据对。

**公式与推导** DeepONet
$$G(u)(y)\approx\sum_{k=1}^p b_k(u)t_k(y),$$
其中 branch 网络编码输入函数，trunk 网络编码查询点。损失
$$L=\sum_i\|G(u_i)(y_i)-G_i\|^2.$$

**物理/几何意义** 在函数空间学习算子。

**极限/特例** 万能算子逼近定理。

**典型应用** 参数化 PDE 求解。

#### 知识点17.5 神经算子与 Fourier 神经算子

**定义** 用积分算子
$$(\mathcal K v)(x)=\int\kappa(x,y)v(y)dy.$$

**假设** 核平移不变 $\kappa(x,y)=\kappa(x-y)$。

**公式与推导** 卷积定理
$$(\mathcal K v)=\mathcal F^{-1}[\mathcal F[\kappa]\cdot\mathcal F[v]].$$
FNO 在频域参数化
$$v_{l+1}=\sigma\left(Wv_l+\mathcal F^{-1}(R_l\cdot\mathcal F[v_l])\right).$$

**物理/几何意义** 谱卷积学习全局耦合。

**极限/特例** 截断高频。

**典型应用** 湍流、天气。

#### 知识点17.6 高斯过程回归与不确定性

**定义** $f\sim\mathcal{GP}(m,k)$。

**假设** 观测 $y=f(x)+\epsilon,\ \epsilon\sim\mathcal N(0,\sigma_n^2)$。

**公式与推导** 联合高斯
$$\begin{pmatrix}y\\f_*\end{pmatrix}\sim\mathcal N\left(0,\begin{pmatrix}K+\sigma_n^2I&K_*\\K_*^\top&K_{**}\end{pmatrix}\right).$$
条件分布
$$f_*|y\sim\mathcal N\!\left(K_*^\top(K+\sigma_n^2I)^{-1}y,\ K_{**}-K_*^\top(K+\sigma_n^2I)^{-1}K_*\right).$$

**物理/几何意义** 核函数编码先验光滑性，方差给出置信区间。

**极限/特例** 核宽度控制光滑度。

**典型应用** 代理模型、贝叶斯优化。

#### 知识点17.7 贝叶斯优化

**定义** 用代理模型指导昂贵函数寻优。

**假设** 目标 $f$ 评估代价高。

**公式与推导** 采集函数 EI：
$$\text{EI}(x)=\mathbb E[\max(f(x)-f^+,0)].$$
对高斯代理，有闭式
$$\text{EI}=( \mu-f^+)\Phi(z)+\sigma\phi(z),\quad z=\frac{\mu-f^+}{\sigma}.$$

**物理/几何意义** 权衡探索与利用。

**极限/特例** UCB、PI 为替代。

**典型应用** 实验设计、超参调优。

#### 知识点17.8 数据驱动本构与发现

**定义** 用网络/稀疏回归发现本构或控制方程。

**假设** 候选库 $\Theta$。

**公式与推导** SINDy：
$$\dot X=\Theta(X)\Xi,$$
稀疏回归
$$\min_\Xi\|\dot X-\Theta\Xi\|_2^2+\lambda\|\Xi\|_1.$$
物理约束网络则把热力学不等式加入损失。

**物理/几何意义** 从数据提炼可解释规律。

**极限/特例** 噪声大时稀疏识别退化。

**典型应用** 材料模型、流体本构。

---
⬆ 上级领域：[[计算与人工智能]]
