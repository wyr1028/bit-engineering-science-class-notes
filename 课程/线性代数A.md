---
课程名称: 线性代数A
学分: 4
总学时: 64
考核方式: 考试
课程属性: 必修
课程性质: 公共基础课程
开课单位: 数学与统计学院
领域: 数学
tags:
  - 课程
  - 领域/数学
---

# 线性代数A


### 第1章 线性方程组与矩阵表示

#### 1.1 线性方程组与增广矩阵

**定义** 含 $m$ 个方程、$n$ 个未知量 $x_1,\dots,x_n$ 的线性方程组为
$$
\begin{cases}
a_{11}x_1+a_{12}x_2+\cdots+a_{1n}x_n=b_1,\\
a_{21}x_1+a_{22}x_2+\cdots+a_{2n}x_n=b_2,\\
\qquad\vdots\\
a_{m1}x_1+a_{m2}x_2+\cdots+a_{mn}x_n=b_m.
\end{cases}
$$
记 $A=(a_{ij})\in\mathbb{R}^{m\times n}$，$\mathbf{x}=(x_1,\dots,x_n)^{\mathsf T}$，$\mathbf{b}=(b_1,\dots,b_m)^{\mathsf T}$，则方程组等价于矩阵方程
$$
A\mathbf{x}=\mathbf{b}.
$$
**假设** $a_{ij},b_i\in\mathbb{R}$（或一般域 $\mathbb{F}$），未知量为一次且无乘积项。

**推导**（矩阵乘法的展开） 由 $A\mathbf{x}$ 的第 $i$ 个分量为 $\sum_{j=1}^{n}a_{ij}x_j$，令其等于 $b_i$，即得原方程组，$\therefore$ 两表示一一对应。

**意义** 把“解方程”转化为“矩阵列的线性组合”：$A\mathbf{x}=\sum_{j=1}^{n}x_j\mathbf{a}_j=\mathbf{b}$，其中 $\mathbf{a}_j$ 为 $A$ 的第 $j$ 列。解存在当且仅当 $\mathbf{b}$ 落在诸列张成的空间中。

**特例** $n=1$ 时退化为 $a_{i1}x_1=b_i$，相容条件为所有 $b_i/a_{i1}$ 相等；$A$ 为方阵且可逆时解唯一为 $\mathbf{x}=A^{-1}\mathbf{b}$。

**应用** 桁架结构节点受力平衡：设第 $k$ 个节点受各杆内力 $x_j$ 沿方向余弦 $(c_{kx},c_{ky})$ 的分量之和与外载 $b_k$ 平衡，则
$$
\sum_{j\in\mathcal{N}(k)}c_{kj}x_j=b_k,
$$
即 $A\mathbf{x}=\mathbf{b}$；节点数即方程数，杆数即未知量数。

#### 1.2 高斯消元法与初等行变换

**定义** 对增广矩阵 $(A\mid\mathbf{b})$ 施行三种初等行变换：(1) 交换两行；(2) 某行乘非零常数；(3) 把某行的 $k$ 倍加到另一行。

**假设** 变换仅在行上进行，不改变未知量排列。

**公式** 三种变换均对应左乘一个可逆初等矩阵 $E$，故 $(A\mid\mathbf{b})\mapsto(E A\mid E\mathbf{b})$。

**推导** 设原方程 $A\mathbf{x}=\mathbf{b}$。左乘可逆 $E$ 得 $EA\mathbf{x}=E\mathbf{b}$，即新方程。反之由 $E^{-1}$ 可还原，$\therefore$ 解集不变。这证明初等行变换是同解变换。

**意义** 消元的几何意义是把约束超平面作可逆线性重组，不改变其交（解集）。

**特例** 若某步出现 $0=c\ne0$，说明变换后第 $i$ 行为矛盾方程，原方程组无解（这等价于 $\operatorname{rank}A<\operatorname{rank}(A\mid\mathbf{b})$）。

**应用** 电路网络分析：对节点写 KCL、对回路写 KVL 得到线性方程组，用消元求各支路电流。

#### 1.3 行阶梯形与行最简形

**定义** 矩阵称为行阶梯形（REF），若每行首个非零元素（主元）所在列严格位于上一行主元之右，且全零行在底部。进一步若每个主元为 $1$ 且其所在列其余元素为 $0$，则称为行最简形（RREF）。

**假设** 只做行变换。

**推导** 高斯消元：先用第 $1$ 行消去第 $1$ 列下方元素，再对子矩阵递归。设第 $j$ 步主元为 $a^{(j)}_{jj}\ne0$，则
$$
\text{第 }i\text{行}\leftarrow\text{第 }i\text{行}-\frac{a^{(j)}_{ij}}{a^{(j)}_{jj}}\times\text{第 }j\text{行},\quad i>j,
$$
使 $a^{(j+1)}_{ij}=0$。回代由最底非零行向上依次求出主元变量，$\therefore$ 得到唯一 RREF（RREF 唯一性定理）。

**意义** RREF 给出解的参数化：主元列对应基本变量，非主元列对应自由变量。

**特例** 每列都有主元且 $m=n$ 时 REF 为上三角，回代即得唯一解。

**应用** 结构力学中把刚度方程化为 REF 以判断机构（机制）自由度与约束冗余。

#### 1.4 齐次方程组与解空间

**定义** $A\mathbf{x}=\mathbf{0}$ 称为齐次方程组，其解集记为 $\mathcal{N}(A)$（零空间）。

**推导** 任取 $\mathbf{u},\mathbf{v}\in\mathcal{N}(A)$，$\alpha,\beta\in\mathbb{R}$，
$$
A(\alpha\mathbf{u}+\beta\mathbf{v})=\alpha A\mathbf{u}+\beta A\mathbf{v}=\mathbf{0},
$$
$\therefore$ $\mathcal{N}(A)$ 对加法与数乘封闭，是子空间。其维数 $\dim\mathcal{N}(A)=n-\operatorname{rank}A$（秩-零化度定理，见 §6）。

**意义** 齐次方程的解是“自平衡”的内力/位移模式：无外载也能存在的状态。

**特例** $\operatorname{rank}A=n$ 时只有零解；$\operatorname{rank}A<n$ 时有非零解且解空间至少一条直线。

**应用** 桁架零特征位移（刚体模态）即 $\mathcal{N}(K)$ 中的向量，对应整体平移/转动，是刚度矩阵奇异的原因。

#### 1.5 非齐次方程组解的结构

**定理** 若 $A\mathbf{x}=\mathbf{b}$ 有特解 $\mathbf{x}_p$，则通解为
$$
\mathbf{x}=\mathbf{x}_p+\mathbf{x}_h,\qquad \mathbf{x}_h\in\mathcal{N}(A).
$$
**推导** 设 $\mathbf{x}$ 为任一解，$\mathbf{x}_p$ 为特解。则
$$
A(\mathbf{x}-\mathbf{x}_p)=A\mathbf{x}-A\mathbf{x}_p=\mathbf{b}-\mathbf{b}=\mathbf{0},
$$
$\therefore$ $\mathbf{x}-\mathbf{x}_p=\mathbf{x}_h\in\mathcal{N}(A)$，即 $\mathbf{x}=\mathbf{x}_p+\mathbf{x}_h$。反之 $\mathbf{x}_p+\mathbf{x}_h$ 满足方程。

**相容判据（Rouché–Capelli）** $A\mathbf{x}=\mathbf{b}$ 有解 $\iff \operatorname{rank}(A)=\operatorname{rank}(A\mid\mathbf{b})$。

**推导** 有解 $\iff\mathbf{b}\in\operatorname{span}\{\mathbf{a}_1,\dots,\mathbf{a}_n\}$ $\iff$ 加入 $\mathbf{b}$ 不增加列空间的维数 $\iff \operatorname{rank}(A\mid\mathbf{b})=\operatorname{rank}A$。

**意义** 解集是仿射子空间：若 $\mathcal{N}(A)=\{\mathbf{0}\}$ 则为单点，否则为过 $\mathbf{x}_p$ 的平行子空间。

**特例** $m=n$ 且 $\det A\ne0$ 时解唯一 $\mathbf{x}=A^{-1}\mathbf{b}$。

**应用** 超静定结构（方程数多于独立未知量）中相容性由几何协调方程保证，$\operatorname{rank}$ 判据给出自平衡预应力状态的存在条件。

### 第2章 矩阵代数与分块矩阵

#### 2.1 矩阵乘法及其结合律

**定义** $A=(a_{ij})\in\mathbb{R}^{m\times n}$，$B=(b_{jk})\in\mathbb{R}^{n\times p}$，乘积 $C=AB\in\mathbb{R}^{m\times p}$，其中
$$
c_{ik}=\sum_{j=1}^{n}a_{ij}b_{jk}.
$$
**假设** 维数匹配：左矩阵列数等于右矩阵行数。

**推导**（结合律） 设 $A\in\mathbb{R}^{m\times n},B\in\mathbb{R}^{n\times p},C\in\mathbb{R}^{p\times q}$，则
$$
\big((AB)C\big)_{il}=\sum_{k=1}^{p}(AB)_{ik}c_{kl}=\sum_{k=1}^{p}\sum_{j=1}^{n}a_{ij}b_{jk}c_{kl},
$$
$$
\big(A(BC)\big)_{il}=\sum_{j=1}^{n}a_{ij}(BC)_{jl}=\sum_{j=1}^{n}\sum_{k=1}^{p}a_{ij}b_{jk}c_{kl}.
$$
两者为同一有限和的重排，$\therefore (AB)C=A(BC)$。

**意义** 矩阵是线性变换的表示，乘法对应变换复合，结合律即复合的结合性；一般 $AB\ne BA$（非交换），$AB=0$ 不蕴含 $A=0$ 或 $B=0$（零因子）。

**特例** 标量矩阵 $\lambda I$ 与一切同阶矩阵交换；两个上三角阵之积仍为上三角。

**应用** 连续坐标变换：$\mathbf{x}'=B\mathbf{x}$，$\mathbf{x}''=A\mathbf{x}'$，则 $\mathbf{x}''=(AB)\mathbf{x}$，复合矩阵即乘积。

#### 2.2 转置、对称与反对称

**定义** $(A^{\mathsf T})_{ij}=a_{ji}$。若 $A^{\mathsf T}=A$ 称对称；$A^{\mathsf T}=-A$ 称反对称。

**公式** $(AB)^{\mathsf T}=B^{\mathsf T}A^{\mathsf T}$，$(A+B)^{\mathsf T}=A^{\mathsf T}+B^{\mathsf T}$，$(\alpha A)^{\mathsf T}=\alpha A^{\mathsf T}$。

**推导** 逐元素：
$$
\big((AB)^{\mathsf T}\big)_{ij}=(AB)_{ji}=\sum_k a_{jk}b_{ki}=\sum_k (B^{\mathsf T})_{ik}(A^{\mathsf T})_{kj}=(B^{\mathsf T}A^{\mathsf T})_{ij}.
$$
$\therefore (AB)^{\mathsf T}=B^{\mathsf T}A^{\mathsf T}$。

**分解定理** 任意方阵 $A$ 可唯一分解为
$$
A=\frac{A+A^{\mathsf T}}{2}+\frac{A-A^{\mathsf T}}{2}=S+K,
$$
$S$ 对称、$K$ 反对称。证明：显然右端两项分别对称与反对称且和为 $A$；若 $A=S_1+K_1=S_2+K_2$，则 $S_1-S_2=K_2-K_1$ 既对称又反对称，故为零阵，唯一性成立。

**意义** 应力张量分解为静水压（球量 $\sigma_m I$）与偏应力（对称无迹）；速度梯度张量分解为对称应变率与反对称涡量。

**特例** 反对称 $3\times3$ 矩阵可写成
$$
K=\begin{pmatrix}0&-k_3&k_2\\ k_3&0&-k_1\\ -k_2&k_1&0\end{pmatrix}=[\mathbf{k}]_\times,\qquad K\mathbf{x}=\mathbf{k}\times\mathbf{x},
$$
对应叉乘。

**应用** 刚体角速度 $\boldsymbol{\omega}$ 与速度场 $\mathbf{v}=\boldsymbol{\omega}\times\mathbf{r}=[\boldsymbol{\omega}]_\times\mathbf{r}$。

#### 2.3 分块矩阵与分块乘法

**定义** 将矩阵按行列分成若干子块，如
$$
A=\begin{pmatrix}A_{11}&A_{12}\\ A_{21}&A_{22}\end{pmatrix},\qquad
B=\begin{pmatrix}B_{11}&B_{12}\\ B_{21}&B_{22}\end{pmatrix}.
$$
**假设** 各子块维数使相应乘积有意义。

**公式**
$$
AB=\begin{pmatrix}A_{11}B_{11}+A_{12}B_{21}&A_{11}B_{12}+A_{12}B_{22}\\ A_{21}B_{11}+A_{22}B_{21}&A_{21}B_{12}+A_{22}B_{22}\end{pmatrix}.
$$
**推导** $(AB)_{ij}=\sum_k a_{ik}b_{kj}$，按块求和即得；与普通乘法规则形式相同，因为块内元素乘法为标量乘法，块间加法逐元素对应。

**意义** 分块是“矩阵的矩阵”，便于表示复合系统、分块求解与大矩阵计算。

**特例** 分块对角阵 $\operatorname{diag}(A_1,\dots,A_k)$ 的行列式等于 $\prod_i\det A_i$，逆为 $\operatorname{diag}(A_1^{-1},\dots,A_k^{-1})$。

**应用** 多自由度系统质量矩阵常为分块对角（部件独立），刚度矩阵含耦合块。

#### 2.4 迹及其性质

**定义** 方阵 $A=(a_{ij})_{n\times n}$ 的迹 $\operatorname{tr}A=\sum_{i=1}^{n}a_{ii}$。

**公式** $\operatorname{tr}(AB)=\operatorname{tr}(BA)$，$\operatorname{tr}(A+B)=\operatorname{tr}A+\operatorname{tr}B$，$\operatorname{tr}(A^{\mathsf T})=\operatorname{tr}A$。

**推导**
$$
\operatorname{tr}(AB)=\sum_i (AB)_{ii}=\sum_i\sum_j a_{ij}b_{ji}=\sum_j\sum_i b_{ji}a_{ij}=\sum_j (BA)_{jj}=\operatorname{tr}(BA).
$$
另外 $\operatorname{tr}(P^{-1}AP)=\operatorname{tr}(APP^{-1})=\operatorname{tr}A$，$\therefore$ 迹是相似不变量。

**意义** 迹等于特征值之和（§8），反映矩阵的“散度型”总量。

**特例** 反对称阵迹为零；$\operatorname{tr}(I_n)=n$。

**应用** 应力张量第一不变量 $I_1=\operatorname{tr}\boldsymbol{\sigma}=3\sigma_m$，与坐标无关。

#### 2.5 矩阵多项式与幂

**定义** 对多项式 $p(t)=\sum_{k=0}^{d}c_kt^k$，定义 $p(A)=\sum_{k=0}^{d}c_kA^k$，$A^0=I$。

**公式** 若 $A=PJP^{-1}$，则 $A^k=PJ^kP^{-1}$，$p(A)=Pp(J)P^{-1}$。

**推导** $A^2=(PJP^{-1})(PJP^{-1})=PJ^2P^{-1}$，归纳得 $A^k=PJ^kP^{-1}$，线性组合即得 $p(A)=Pp(J)P^{-1}$。

**意义** 多项式演算可搬到相似标准形，简化计算。

**特例** 对幂零阵 $N$（$N^k=0$），$p(N)$ 为有限项之和。

**应用** 离散动力系统 $\mathbf{x}_{k}=A^k\mathbf{x}_0$ 的稳定性由 $A$ 的特征值模长决定。

#### 2.6 应用：应力张量、转动惯量与坐标变换

**应力张量** 柯西应力 $\boldsymbol{\sigma}$ 为对称二阶张量，面力
$$
\mathbf{t}=\boldsymbol{\sigma}\mathbf{n},
$$
$\mathbf{n}$ 为单位外法向。$(\boldsymbol{\sigma}\mathbf{n})_i=\sum_j\sigma_{ij}n_j$ 即矩阵-向量乘积。

**转动惯量张量** 对连续体
$$
\mathbf{J}=\int_V\rho\big(\|\mathbf{r}\|^2I-\mathbf{r}\mathbf{r}^{\mathsf T}\big)\,dV,
$$
其中 $\mathbf{r}\mathbf{r}^{\mathsf T}$ 为并矢（秩-1 矩阵），$\mathbf{J}$ 对称。角动量 $\mathbf{L}=\mathbf{J}\boldsymbol{\omega}$。

**坐标变换** 设 $\mathbf{x}'=Q\mathbf{x}$，$Q$ 正交，则应力按二阶张量变换
$$
\boldsymbol{\sigma}'=Q\boldsymbol{\sigma}Q^{\mathsf T}.
$$
**推导** 要求面力变换一致：$\mathbf{t}'=Q\mathbf{t}$ 且 $\mathbf{n}'=Q\mathbf{n}$，则 $\boldsymbol{\sigma}'\mathbf{n}'=Q\boldsymbol{\sigma}\mathbf{n}=Q\boldsymbol{\sigma}Q^{\mathsf T}\mathbf{n}'$，$\therefore \boldsymbol{\sigma}'=Q\boldsymbol{\sigma}Q^{\mathsf T}$。

**意义** 张量分量随坐标旋转按合同变换改变，特征值（主应力）不变。

**特例** 平面应力
$$
\boldsymbol{\sigma}=\begin{pmatrix}\sigma_x&\tau_{xy}\\\tau_{xy}&\sigma_y\end{pmatrix},
$$
其主应力为 $\sigma_{1,2}=\frac{\sigma_x+\sigma_y}{2}\pm\sqrt{\big(\frac{\sigma_x-\sigma_y}{2}\big)^2+\tau_{xy}^2}$，将在 §11 用特征值重新导出。

#### 2.7 应用：线性系统与叠加原理

**线性叠加** 系统 $A(\alpha\mathbf{x}_1+\beta\mathbf{x}_2)=\alpha A\mathbf{x}_1+\beta A\mathbf{x}_2$。力学小变形线弹性问题满足叠加原理：多组载荷的响应等于各自响应之和。

**应用** 梁的多个集中力挠度叠加；热-力耦合中温度载荷与机械载荷解相加。这是有限元法把总体载荷向量按单元组装的理论基础。

### 第3章 行列式

#### 3.1 排列、逆序数与行列式定义

**定义** $n$ 阶排列 $\sigma=(\sigma_1,\dots,\sigma_n)$ 是 $1,\dots,n$ 的一个重排。逆序数 $\tau(\sigma)$ 为满足 $i<j$ 且 $\sigma_i>\sigma_j$ 的对数；$\operatorname{sgn}\sigma=(-1)^{\tau(\sigma)}$。$n$ 阶行列式
$$
\det A=\sum_{\sigma\in S_n}\operatorname{sgn}(\sigma)\,a_{1\sigma_1}a_{2\sigma_2}\cdots a_{n\sigma_n}.
$$
**假设** $A$ 为 $n\times n$ 方阵。

**推导**（与递归展开一致） 按定义验证 $n=1$：$\det(a_{11})=a_{11}$。设按第一行展开公式 $\det A=\sum_j(-1)^{1+j}a_{1j}M_{1j}$，其中 $M_{1j}$ 为余子式。把每个排列按 $\sigma_1=j$ 分组，对每组固定 $a_{1j}$ 后剩余列的排列恰是 $n-1$ 阶排列，符号满足 $\operatorname{sgn}\sigma=(-1)^{1+j}\operatorname{sgn}\sigma'$，$\therefore$ 两种定义等价。

**意义** 行列式是列向量构成的平行多面体的有向体积，也是线性变换的体积伸缩因子。

**特例** 上（下）三角阵 $\det=\prod_i a_{ii}$；$I_n$ 的 $\det=1$。

**应用** 判断线性方程组是否有唯一解、判断基是否为右手系。

#### 3.2 行列式的基本性质

**公式** (1) $\det A^{\mathsf T}=\det A$；(2) 交换两行（列）变号；(3) 某行乘 $k$ 使 $\det$ 乘 $k$；(4) 两行相同或成比例则 $\det=0$；(5) 某行加上另一行的倍数，$\det$ 不变；(6) $\det(AB)=\det A\det B$；(7) $\det(kA)=k^n\det A$。

**推导**（乘积性） 按定义展开
$$
\det(AB)=\sum_{\sigma}\operatorname{sgn}(\sigma)\prod_i\Big(\sum_k a_{ik}b_{k\sigma_i}\Big).
$$
把乘积中每项用排列 $\pi$ 选取列指标，重组得
$$
=\sum_{\pi}\Big(\sum_{\sigma}\operatorname{sgn}(\sigma)\operatorname{sgn}(\pi)\,a_{1\pi_1}\cdots a_{n\pi_n}\Big)b_{\pi_1\sigma(1)}\cdots
=\det A\cdot\det B,
$$
其中用到 $\operatorname{sgn}(\sigma)\operatorname{sgn}(\pi)=\operatorname{sgn}(\sigma\pi)$ 与 $\sum_\sigma \operatorname{sgn}(\sigma\pi)=\sum_\sigma \operatorname{sgn}\sigma$（左乘一固定置换只是重排）。

**意义** $\det(AB)=\det A\det B$ 表明体积伸缩因子可乘。由 $\det(AA^{-1})=1$ 得 $\det(A^{-1})=1/\det A$。

**特例** 若 $A$ 奇异则 $\det A=0$，反之 $\det A\ne0\Rightarrow A$ 可逆。

**应用** 雅可比行列式 $J=\det(\partial x_i/\partial \xi_j)$ 度量坐标变换的局部体积比，出现在虚功原理换元与 N-S 张力项中。

#### 3.3 余子式与 Laplace 展开

**定义** 删除第 $i$ 行第 $j$ 列所得 $n-1$ 阶行列式记为 $M_{ij}$（余子式），代数余子式 $A_{ij}=(-1)^{i+j}M_{ij}$。

**公式（展开定理）**
$$
\det A=\sum_{j=1}^{n}a_{ij}A_{ij}\quad(\text{按第 }i\text{ 行}),\qquad
\det A=\sum_{i=1}^{n}a_{ij}A_{ij}\quad(\text{按第 }j\text{ 列}),
$$
$$
\sum_{j=1}^{n}a_{kj}A_{ij}=0\quad(k\ne i).
$$
**推导** 按定义把 $a_{ij}A_{ij}$ 中属于第 $i$ 行的因子提出，剩余为 $n-1$ 阶行列式；对 $k\ne i$ 的交替式，相当于把第 $i$ 行换成第 $k$ 行后求行列式，出现两行相同，故为零。

**意义** 展开定理把 $n$ 阶行列式化为 $n$ 个 $n-1$ 阶行列式的带号和，是递归计算与理论推导的桥梁。

**特例** 上三角阵按第一列展开立即得 $\det=\prod a_{ii}$。

**应用** 稀疏矩阵（如刚度矩阵）按非零元素展开可大幅降低计算量。

#### 3.4 范德蒙德行列式

**公式** 对 $x_1,\dots,x_n$，
$$
V_n=\begin{vmatrix}
1&x_1&\cdots&x_1^{n-1}\\
1&x_2&\cdots&x_2^{n-1}\\
\vdots&\vdots&&\vdots\\
1&x_n&\cdots&x_n^{n-1}
\end{vmatrix}
=\prod_{1\le i<j\le n}(x_j-x_i).
$$
**推导** 从第 $n$ 列起依次用 $C_k\leftarrow C_k-x_1C_{k-1}$ 消去第一行非对角元，得
$$
V_n=\prod_{i=2}^{n}(x_i-x_1)\cdot V_{n-1}(x_2,\dots,x_n).
$$
递归展开即得 $\prod_{i<j}(x_j-x_i)$。

**意义** 节点互异时 $V_n\ne0$，说明多项式插值、幂基线性无关。

**特例** $n=2$：$V_2=x_2-x_1$。

**应用** 结构动力学中特征频率的辨识、离散数据的 Vandermonde 系统求解（见最小二乘）。

#### 3.5 克拉默法则

**定理** 设 $A$ 可逆，$A\mathbf{x}=\mathbf{b}$ 的唯一解为
$$
x_j=\frac{\det A_j}{\det A},\qquad A_j=[\mathbf{a}_1,\dots,\mathbf{a}_{j-1},\mathbf{b},\mathbf{a}_{j+1},\dots,\mathbf{a}_n].
$$
**推导** 由 $A\mathbf{x}=\mathbf{b}$ 得 $\mathbf{x}=A^{-1}\mathbf{b}$。用伴随矩阵表示 $A^{-1}=\frac{1}{\det A}\operatorname{adj}A$（§4），故 $x_j=\frac{1}{\det A}\sum_i A_{ij}b_i$。而 $\sum_i A_{ij}b_i$ 恰是把 $A$ 第 $j$ 列换成 $\mathbf{b}$ 后按第 $j$ 列 Laplace 展开的结果，即 $\det A_j$，$\therefore x_j=\det A_j/\det A$。

**意义** 把解写成行列式之比，理论价值高于计算价值。

**特例** 齐次方程 $A\mathbf{x}=\mathbf{0}$ 有非零解 $\iff \det A=0$。

**应用** 2×2 桁架、电路网络的显式解；主应力公式的代数来源。

#### 3.6 行列式的几何意义与变换

**几何意义** 设 $\mathbf{a}_1,\dots,\mathbf{a}_n$ 为列向量，则 $|\det A|$ 等于它们张成的平行多面体体积；$\operatorname{sgn}\det A$ 表示定向。

**推导** 体积函数 $V(\mathbf{a}_1,\dots,\mathbf{a}_n)$ 满足：对每个变量线性、交换两变量变号、单位立方体体积为 $1$。这些恰是行列式的刻画性质，由唯一性得 $V=\det A$。

**变换公式** 若 $\mathbf{y}=A\mathbf{x}$，则 $d\mathbf{y}=|\det A|\,d\mathbf{x}$；对一般非线性映射 $\mathbf{y}=\mathbf{f}(\mathbf{x})$，
$$
dV_y=\Big|\det\frac{\partial \mathbf{y}}{\partial \mathbf{x}}\Big|\,dV_x,\qquad \frac{\partial \mathbf{y}}{\partial \mathbf{x}}=\Big(\frac{\partial y_i}{\partial x_j}\Big).
$$
**极限特例** 当 $\det A\to0$，至少一个方向被压成零体积，变换降维。

**应用** 变形梯度 $\mathbf{F}=\partial\mathbf{x}/\partial\mathbf{X}$，$J=\det\mathbf{F}$ 为体积比；不可压缩条件 $J=1$；连续性方程 $\dot\rho+\rho\nabla\cdot\mathbf{v}=0$ 可由 $dJ/dt=J\,\nabla\cdot\mathbf{v}$ 推出。

### 第4章 逆矩阵与初等变换

#### 4.1 逆矩阵的定义与唯一性

**定义** 方阵 $A$ 可逆，若存在 $B$ 使 $AB=BA=I$，记 $B=A^{-1}$。

**推导**（唯一性） 设 $AB=BA=I$ 且 $AC=CA=I$，则
$$
B=BI=B(AC)=(BA)C=IC=C,
$$
$\therefore$ 逆唯一。

**公式** $(AB)^{-1}=B^{-1}A^{-1}$；$(A^{\mathsf T})^{-1}=(A^{-1})^{\mathsf T}$；$(\alpha A)^{-1}=\alpha^{-1}A^{-1}$。

**推导** $(AB)(B^{-1}A^{-1})=A(BB^{-1})A^{-1}=AA^{-1}=I$；$A^{\mathsf T}(A^{-1})^{\mathsf T}=(A^{-1}A)^{\mathsf T}=I$。

**意义** 可逆等价于线性变换是双射，等价于 $\det A\ne0$，等价于列（行）线性无关。

**特例** 对角阵 $D=\operatorname{diag}(d_1,\dots,d_n)$ 可逆 $\iff d_i\ne0$，$D^{-1}=\operatorname{diag}(d_1^{-1},\dots,d_n^{-1})$。

**应用** 求 $A^{-1}\mathbf{b}$ 即解线性系统；坐标反变换。

#### 4.2 伴随矩阵

**定义** $\operatorname{adj}A=(A_{ji})$，即代数余子式矩阵的转置。

**公式** $A\,\operatorname{adj}A=(\det A)I$，从而 $A^{-1}=\frac{1}{\det A}\operatorname{adj}A$。

**推导** $(A\operatorname{adj}A)_{ij}=\sum_k a_{ik}A_{jk}$。当 $i=j$ 时为 $\det A$；当 $i\ne j$ 时为第 $j$ 行替换第 $i$ 行后的行列式，两行相同故为 $0$。$\therefore A\operatorname{adj}A=(\det A)I$。

**意义** 给出逆的显式，是克拉默法则的来源。

**特例** $2\times2$：$\begin{pmatrix}a&b\\c&d\end{pmatrix}^{-1}=\dfrac{1}{ad-bc}\begin{pmatrix}d&-b\\-c&a\end{pmatrix}$。

**应用** 解析求逆用于小系统与公式推导；大系统用数值分解。

#### 4.3 初等矩阵与求逆

**定义** 对 $I$ 施行一次初等行变换所得矩阵称为初等矩阵。

**公式** 若 $E$ 由 $I$ 经行变换得到，则对任意 $A$，$EA$ 等于对 $A$ 施行同一行变换。三种初等矩阵的行列式分别为 $-1,\ k,\ 1$。

**推导** 只用矩阵乘法验证：交换 $I$ 的第 $i,j$ 行得 $E$，则 $(EA)$ 的第 $k$ 行是 $A$ 的相应行重排，与行变换一致。乘 $k$ 与倍加同理。

**求逆算法** 若 $A$ 可逆，则 $A\sim I$，等价地存在初等矩阵 $E_1,\dots,E_s$ 使 $E_s\cdots E_1A=I$，于是
$$
A^{-1}=E_s\cdots E_1 I.
$$
对增广矩阵 $(A\mid I)$ 行变换化左块为 $I$，右块即为 $A^{-1}$。

**意义** 初等矩阵刻画可逆变换可由“剪切、缩放、置换”复合而成。

**特例** 倍加变换不改变行列式，故可用于求 $\det$ 而不变符号。

**应用** 解多个右端项 $AX=B$ 时一次求逆或一次 LU 分解。

#### 4.4 LU 分解

**定义** 若 $A=LU$，$L$ 为单位下三角，$U$ 为上三角，称为 LU 分解。

**推导** 高斯消元每步消去用初等矩阵 $E_k$（倍加），$E_m\cdots E_1 A=U$，故 $A=E_1^{-1}\cdots E_m^{-1}U=LU$，其中 $L=E_1^{-1}\cdots E_m^{-1}$ 为单位下三角且乘子恰为消元系数。

**意义** 把 $A\mathbf{x}=\mathbf{b}$ 拆成 $L\mathbf{y}=\mathbf{b}$（前代）与 $U\mathbf{x}=\mathbf{y}$（回代），复杂度 $O(n^3/3)$。

**特例** 对称正定矩阵可用 Cholesky 分解 $A=LL^{\mathsf T}$，$L$ 下三角。

**应用** 有限元刚度方程 $K\mathbf{u}=\mathbf{f}$ 反复求解。

#### 4.5 分块求逆与 Schur 补

**公式** 设 $A=\begin{pmatrix}A_{11}&A_{12}\\A_{21}&A_{22}\end{pmatrix}$，若 $A_{11}$ 可逆，令 Schur 补 $S=A_{22}-A_{21}A_{11}^{-1}A_{12}$，则 $A$ 可逆 $\iff S$ 可逆，且
$$
A^{-1}=\begin{pmatrix}A_{11}^{-1}+A_{11}^{-1}A_{12}S^{-1}A_{21}A_{11}^{-1}&-A_{11}^{-1}A_{12}S^{-1}\\ -S^{-1}A_{21}A_{11}^{-1}&S^{-1}\end{pmatrix}.
$$
**推导** 做块消元
$$
\begin{pmatrix}I&0\\-A_{21}A_{11}^{-1}&I\end{pmatrix}A=\begin{pmatrix}A_{11}&A_{12}\\0&S\end{pmatrix},
$$
两边求逆并整理即得。

**意义** Schur 补是不在 $A_{11}$ 中体现的“有效算子”，其奇异性对应子系统的失稳。

**特例** 块对角时 $S=A_{22}$，逆为块逆。

**应用** 子结构法（superelement）把大结构划分为子域，通过 Schur 补凝聚内部自由度到界面。

### 第5章 向量空间与线性相关性

#### 5.1 向量空间与子空间

**定义** 域 $\mathbb{F}$ 上向量空间 $V$ 为集合，带加法 $+:V\times V\to V$ 与数乘 $\cdot:\mathbb{F}\times V\to V$，满足八条公理（加法交换/结合/零元/负元，数乘分配/结合/单位）。

**定义** 非空 $W\subseteq V$ 为子空间，若对任意 $\mathbf{u},\mathbf{v}\in W$，$\alpha,\beta\in\mathbb{F}$，$\alpha\mathbf{u}+\beta\mathbf{v}\in W$。

**推导** 验证若干子空间：$\mathcal{N}(A)$、列空间 $\mathcal{C}(A)$、行空间 $\mathcal{R}(A)$、特征子空间 $E_\lambda$、解空间均在相应封闭性下构成子空间。

**意义** 子空间是“线性闭合”的集合，维数是自由度计数。

**特例** $\{\mathbf{0}\}$ 与 $V$ 是平凡子空间；两个子空间之交仍是子空间，并一般不是。

**应用** 力学的运动许可位移场构成函数空间的子空间；应变场属于对称张量子空间。

#### 5.2 线性组合、张成与生成集

**定义** $\operatorname{span}\{\mathbf{v}_1,\dots,\mathbf{v}_k\}=\{\sum_{i=1}^{k}\alpha_i\mathbf{v}_i:\alpha_i\in\mathbb{F}\}$。

**推导** 该集合对线性组合封闭：$\sum\beta_j\mathbf{w}_j=\sum\beta_j\sum_i\alpha_{ji}\mathbf{v}_i=\sum_i(\sum_j\beta_j\alpha_{ji})\mathbf{v}_i$，$\therefore$ 为子空间。

**意义** 张成空间是“可达状态”的全体；矩阵列空间 $\mathcal{C}(A)=\operatorname{span}\{\mathbf{a}_1,\dots,\mathbf{a}_n\}$ 即 $A\mathbf{x}$ 的取值范围。

**特例** 单个非零向量张成一条过原点的直线；两个不共线向量张成平面。

**应用** 桁架节点的位移若可表示为若干基模态的线性组合，则整体响应由模态叠加给出。

#### 5.3 线性无关与基

**定义** 向量组 $\mathbf{v}_1,\dots,\mathbf{v}_k$ 线性无关，若
$$
\sum_{i=1}^{k}\alpha_i\mathbf{v}_i=\mathbf{0}\ \Rightarrow\ \alpha_1=\cdots=\alpha_k=0.
$$
否则线性相关。$V$ 的一组基是线性无关的生成集。

**推导**（基表示唯一性） 设 $\mathbf{v}=\sum\alpha_i\mathbf{v}_i=\sum\beta_i\mathbf{v}_i$，则 $\sum(\alpha_i-\beta_i)\mathbf{v}_i=\mathbf{0}$，由无关性得 $\alpha_i=\beta_i$，$\therefore$ 坐标唯一。

**意义** 基是“最省”的坐标架；维数 $\dim V$ 为任一基的向量个数（基的势不变）。

**特例** $\mathbb{R}^n$ 的标准基 $\mathbf{e}_1,\dots,\mathbf{e}_n$；多项式空间 $\mathcal{P}_n$ 的基 $1,t,\dots,t^n$，维数 $n+1$。

**应用** 有限元形函数即单元位移空间的基；模态坐标是物理坐标在模态基下的表示。

#### 5.4 极大无关组与替换定理

**定理（Steinitz 替换）** 若 $W=\operatorname{span}\{\mathbf{w}_1,\dots,\mathbf{w}_m\}$ 且 $\mathbf{v}_1,\dots,\mathbf{v}_k\in W$ 线性无关，则 $k\le m$，且可用诸 $\mathbf{v}_i$ 替换等量个 $\mathbf{w}_j$ 仍张成 $W$。

**推导** 由 $\mathbf{v}_1\ne\mathbf{0}$ 写成 $\mathbf{w}$ 的线性组合，必含某 $\mathbf{w}_{j_1}$ 的系数非零，解出 $\mathbf{w}_{j_1}$ 并由其余 $\mathbf{w}$ 与 $\mathbf{v}_1$ 表示；逐步替换，每次增加一个 $\mathbf{v}$，最终得 $k$ 个。若 $k>m$，替换完后 $\mathbf{v}_{m+1}$ 可由前 $m$ 个表示，与无关性矛盾，$\therefore k\le m$。

**意义** 无关组长度不超过任一生成集长度；由此定义维数良定。

**特例** $\mathbb{R}^n$ 中任意 $n+1$ 个向量必相关。

**应用** 判断超静定次数与独立约束数。

#### 5.5 坐标与基变换

**定义** 设 $\mathcal{B}=\{\mathbf{b}_1,\dots,\mathbf{b}_n\}$ 为 $V$ 的基，$\mathbf{v}=\sum x_i\mathbf{b}_i$，则 $[\mathbf{v}]_{\mathcal{B}}=(x_1,\dots,x_n)^{\mathsf T}$。

**公式** 设两组基，过渡矩阵 $P$ 满足 $\mathbf{b}'_j=\sum_i P_{ij}\mathbf{b}_i$，即 $\mathcal{B}'=\mathcal{B}P$，则
$$
[\mathbf{v}]_{\mathcal{B}}=P\,[\mathbf{v}]_{\mathcal{B}'},\qquad [\mathbf{v}]_{\mathcal{B}'}=P^{-1}[\mathbf{v}]_{\mathcal{B}}.
$$
**推导** $\mathbf{v}=\mathcal{B}'[\mathbf{v}]_{\mathcal{B}'}=\mathcal{B}P[\mathbf{v}]_{\mathcal{B}'}$，与 $\mathbf{v}=\mathcal{B}[\mathbf{v}]_{\mathcal{B}}$ 比较坐标唯一性，$\therefore [\mathbf{v}]_{\mathcal{B}}=P[\mathbf{v}]_{\mathcal{B}'}$。

**意义** 换基不改变向量本身，只改坐标；$P$ 的列是新基在旧基下的坐标。

**特例** 标准基到新基的过渡矩阵列即新基向量。

**应用** 结构分析中整体坐标与局部单元坐标的转换矩阵 $\mathbf{T}$。

#### 5.6 应用：函数空间与应变/应力空间

**函数空间** 满足齐次边界条件的可微函数构成无限维向量空间，内积 $\langle f,g\rangle=\int fg$。微分算子 $L$ 把该空间映到自身，线性性 $L(\alpha f+\beta g)=\alpha Lf+\beta Lg$。

**应变空间** 小应变张量 $\boldsymbol{\varepsilon}=\frac12(\nabla\mathbf{u}+\nabla\mathbf{u}^{\mathsf T})$ 属于对称张量的 6 维子空间（三维情形），其分量 $\varepsilon_{ij}$ 构成坐标。

**应用** 应变能
$$
U=\frac12\int \boldsymbol{\varepsilon}^{\mathsf T}\mathbf{C}\,\boldsymbol{\varepsilon}\,dV
$$
是应变空间上的二次型，$\mathbf{C}$ 为弹性张量（正定对称，见 §11）。

### 第6章 矩阵的秩与线性方程组解理论

#### 6.1 秩的定义与计算

**定义** 矩阵 $A$ 的秩 $\operatorname{rank}A$ 为非零子式的最高阶数，也等于行阶梯形中主元个数，等于列空间维数。

**等价刻画** $\operatorname{rank}A=\dim\mathcal{C}(A)=\dim\mathcal{R}(A)$（行秩=列秩）。

**推导**（行秩=列秩） 对 $A$ 作行变换化 REF，主元所在列对应 $\mathcal{C}(A)$ 的一组基，共有 $r$ 个；同时 REF 的非零行给出 $\mathcal{R}(A)$ 的一组基，也为 $r$ 个，$\therefore$ 行秩 $=$ 列秩 $=r$。

**意义** 秩度量矩阵的“信息量”或线性变换像的维数。

**特例** $\operatorname{rank}A\le\min(m,n)$；满秩时取等。

**应用** 判断结构是否静定、系统是否可控/可观。

#### 6.2 秩的不等式

**公式** $\operatorname{rank}(A+B)\le\operatorname{rank}A+\operatorname{rank}B$；$\operatorname{rank}(AB)\le\min(\operatorname{rank}A,\operatorname{rank}B)$；Sylvester 秩不等式
$$
\operatorname{rank}A+\operatorname{rank}B-n\le\operatorname{rank}(AB)\le\min(\operatorname{rank}A,\operatorname{rank}B),\quad A,B\in\mathbb{R}^{n\times n}.
$$
**推导**（第一个） $\mathcal{C}(A+B)\subseteq\mathcal{C}(A)+\mathcal{C}(B)$，故维数不等式成立。（乘积下界）考虑线性映射 $B:\mathbb{R}^n\to\mathbb{R}^n$ 与 $A$，$\mathcal{C}(AB)=A(\mathcal{C}(B))$，$\dim A(\mathcal{C}(B))=\dim\mathcal{C}(B)-\dim(\mathcal{C}(B)\cap\mathcal{N}(A))$，再用 $\dim(\mathcal{C}(B)\cap\mathcal{N}(A))\le\dim\mathcal{N}(A)=n-\operatorname{rank}(A)$，$\therefore \operatorname{rank}(AB)\ge\operatorname{rank}B-(n-\operatorname{rank}A)=\operatorname{rank}A+\operatorname{rank}B-n$。

**意义** 复合不能增加信息，且可能因核与像重叠而损失。

**特例** $A,B$ 均可逆时 $\operatorname{rank}(AB)=n$。

**应用** 串联系统的能观测性损失；控制中 PBH 判据。

#### 6.3 秩-零化度定理

**定理** 设 $A\in\mathbb{R}^{m\times n}$，则
$$
\dim\mathcal{N}(A)+\operatorname{rank}A=n.
$$
**推导** 把 $A$ 行变换到 RREF，主元列数 $=r=\operatorname{rank}A$，自由变量个数 $=n-r$ 即零空间维数，$\therefore \dim\mathcal{N}(A)=n-r$。

**意义** 定义域维数 $=$ 像维数 $+$ 核维数，是线性映射的“守恒律”。

**特例** 方阵可逆时 $\dim\mathcal{N}(A)=0$。

**应用** 自由度计算：独立广义坐标数 $=$ 总坐标数 $-$ 约束数。

#### 6.4 解的存在唯一性与通解结构

**定理** $A\mathbf{x}=\mathbf{b}$：
(1) 无解 $\iff \operatorname{rank}A<\operatorname{rank}(A\mid\mathbf{b})$；
(2) 唯一解 $\iff \operatorname{rank}A=\operatorname{rank}(A\mid\mathbf{b})=n$；
(3) 无穷多解 $\iff \operatorname{rank}A=\operatorname{rank}(A\mid\mathbf{b})<n$，通解维数 $n-\operatorname{rank}A$。

**推导** 由 §1.5 的相容判据与秩-零化度定理直接得到。

**意义** 用秩一次判定系统的适定性。

**特例** 方阵唯一解 $\iff\det A\ne0$。

**应用** 结构力学：$K\mathbf{u}=\mathbf{f}$ 中 $K$ 奇异（存在刚体模态）时需施加足够约束使 $\operatorname{rank}K=n$。

#### 6.5 广义逆（Moore–Penrose）

**定义** $A^+\in\mathbb{R}^{n\times m}$ 满足
$$
AA^+A=A,\quad A^+AA^+=A^+,\quad (AA^+)^{\mathsf T}=AA^+,\quad (A^+A)^{\mathsf T}=A^+A.
$$
**公式** $A\mathbf{x}=\mathbf{b}$ 的最小范数最小二乘解为 $\mathbf{x}=A^+\mathbf{b}$。

**推导** 若 $A=U\Sigma V^{\mathsf T}$ 为 SVD（§12），取 $A^+=V\Sigma^+U^{\mathsf T}$，其中 $\Sigma^+$ 把非零奇异值取倒数再转置。逐条验证四个条件。

**意义** 对不相容或欠定系统给出唯一“最优”解。

**特例** $A$ 可逆时 $A^+=A^{-1}$；$A$ 满列秩时 $A^+=(A^{\mathsf T}A)^{-1}A^{\mathsf T}$。

**应用** 病态逆问题、参数辨识、机器人逆运动学。

#### 6.6 应用：静定与超静定判断

设三维桁架有 $j$ 个节点、$m$ 根杆、$r$ 个支座约束，平衡方程矩阵 $A\in\mathbb{R}^{3j\times m}$，未知杆力 $\mathbf{x}\in\mathbb{R}^m$。

**判据** 静定 $\iff \operatorname{rank}A=m$ 且 $3j=m+r$；超静定 $\iff m>3j-r$（存在自平衡力态 $\mathcal{N}(A)\ne\{\mathbf{0}\}$）。

**推导** 平衡方程 $A\mathbf{x}=\mathbf{f}$，支座反力作为附加未知量。$\mathcal{N}(A)$ 维数为 $m-\operatorname{rank}A$，即自平衡内力模式数，就是超静定次数。

**应用** 桁架设计、桥梁力学分析。

### 第7章 线性变换

#### 7.1 线性变换的定义与例子

**定义** 映射 $T:V\to W$ 线性，若
$$
T(\alpha\mathbf{u}+\beta\mathbf{v})=\alpha T(\mathbf{u})+\beta T(\mathbf{v})\quad\forall\mathbf{u},\mathbf{v},\ \alpha,\beta.
$$
**推导** 线性性等价于 $T(\mathbf{0})=\mathbf{0}$ 且保持加法与数乘。由 $T(\mathbf{0})=T(0\cdot\mathbf{0})=0\,T(\mathbf{0})=\mathbf{0}$。

**意义** 线性变换是向量空间之间保持结构的映射，矩阵是其坐标表示。

**特例** 恒等 $I$、零映射、伸缩、旋转、投影、微分 $D=\frac{d}{dt}$、积分都是线性变换。

**应用** 应力变换 $\boldsymbol{\sigma}\mapsto Q\boldsymbol{\sigma}Q^{\mathsf T}$ 是张量空间上的线性变换。

#### 7.2 核与像、维数定理

**定义** $\mathcal{N}(T)=\{\mathbf{v}:T\mathbf{v}=\mathbf{0}\}$，$\mathcal{R}(T)=\{T\mathbf{v}:\mathbf{v}\in V\}$，均为子空间。

**定理** $\dim\mathcal{N}(T)+\dim\mathcal{R}(T)=\dim V$。

**推导** 取 $\mathcal{N}(T)$ 的基 $\mathbf{n}_1,\dots,\mathbf{n}_k$，扩充为 $V$ 的基 $\mathbf{n}_1,\dots,\mathbf{n}_k,\mathbf{v}_1,\dots,\mathbf{v}_{n-k}$。证明 $T\mathbf{v}_1,\dots,T\mathbf{v}_{n-k}$ 是 $\mathcal{R}(T)$ 的基：张成性由任意 $\mathbf{v}$ 分解可得；无关性：若 $\sum c_iT\mathbf{v}_i=\mathbf{0}$，则 $T(\sum c_i\mathbf{v}_i)=\mathbf{0}$，故 $\sum c_i\mathbf{v}_i\in\mathcal{N}(T)$，可写成 $\sum d_j\mathbf{n}_j$，由基无关性得所有 $c_i=d_j=0$。$\therefore \dim\mathcal{R}(T)=n-k$。

**意义** 核是“丢失的信息”，像是“保留的信息”，维数守恒。

**特例** 微分算子 $D:\mathcal{P}_n\to\mathcal{P}_{n-1}$ 的核为常数（维数 1），像为 $\mathcal{P}_{n-1}$（维数 $n$），总和 $n+1=\dim\mathcal{P}_n$。

**应用** 结构自由度与约束的匹配分析。

#### 7.3 线性变换的矩阵表示

**定义** 取 $V$ 的基 $\mathcal{B}=\{\mathbf{b}_j\}$、$W$ 的基 $\mathcal{C}=\{\mathbf{c}_i\}$。$T$ 在基下的矩阵 $[T]_{\mathcal{C}\leftarrow\mathcal{B}}$ 的第 $j$ 列为 $[T\mathbf{b}_j]_{\mathcal{C}}$，即
$$
T\mathbf{b}_j=\sum_i A_{ij}\mathbf{c}_i.
$$
**公式** $[T\mathbf{v}]_{\mathcal{C}}=A[\mathbf{v}]_{\mathcal{B}}$。

**推导** $\mathbf{v}=\sum_j x_j\mathbf{b}_j$，则 $T\mathbf{v}=\sum_j x_jT\mathbf{b}_j=\sum_j x_j\sum_i A_{ij}\mathbf{c}_i=\sum_i\big(\sum_j A_{ij}x_j\big)\mathbf{c}_i$，$\therefore$ 坐标按 $A\mathbf{x}$ 变换。

**意义** 选定基后，线性变换与矩阵一一对应；变换复合对应矩阵乘积。

**特例** 恒等变换的矩阵为单位阵。

**应用** 局部单元坐标系下的刚度矩阵到整体坐标的装配。

#### 7.4 基变换与相似

**定理** 设 $V$ 的两组基 $\mathcal{B},\mathcal{B}'$ 过渡矩阵为 $P$（$[\mathbf{v}]_{\mathcal{B}}=P[\mathbf{v}]_{\mathcal{B}'}$），同一线性变换的矩阵满足
$$
A'=P^{-1}AP.
$$
**推导** $[T\mathbf{v}]_{\mathcal{B}'}=P^{-1}[T\mathbf{v}]_{\mathcal{B}}=P^{-1}A[\mathbf{v}]_{\mathcal{B}}=P^{-1}AP[\mathbf{v}]_{\mathcal{B}'}$，$\therefore A'=P^{-1}AP$。

**意义** 相似矩阵描述同一变换在不同基下的表示，特征值、迹、行列式、秩为相似不变量。

**特例** 若 $A'$ 为对角阵，则找到了由特征向量组成的基。

**应用** 模态分析中把耦合的运动方程解耦为对角形式。

#### 7.5 不变子空间

**定义** $W\subseteq V$ 为 $T$-不变子空间，若 $T(W)\subseteq W$。

**公式** 若 $W$ 不变，则存在基使 $[T]$ 为分块上三角
$$
[T]=\begin{pmatrix}A_{11}&A_{12}\\0&A_{22}\end{pmatrix}.
$$
**推导** 取 $W$ 的基并扩充为 $V$ 的基，$T$ 把 $W$ 中向量仍映到 $W$，故表示矩阵左下块为零。

**意义** 不变子空间是系统可解耦的代数结构。

**特例** 特征子空间 $E_\lambda$ 是一维不变子空间；若 $V$ 分解为不变子空间直和，则矩阵分块对角。

**应用** 结构对称性导致模态按对称/反对称分组，形成不变子空间。

#### 7.6 应用：坐标变换与张量

**坐标旋转** 设 $Q$ 正交，$\mathbf{x}'=Q\mathbf{x}$，向量按 $[\mathbf{v}']=Q[\mathbf{v}]$ 变换；二阶张量按 $[\boldsymbol{\sigma}']=Q[\boldsymbol{\sigma}]Q^{\mathsf T}$。

**推导** 张量是不依赖坐标的多线性对象，其分量变换由基变换 $P=Q^{-1}=Q^{\mathsf T}$ 得 $\sigma'=P^{-1}\sigma P=Q\sigma Q^{\mathsf T}$。

**意义** 特征值（主值）与不变量在相似/正交相似下保持。

**应用** 复合材料铺层方向变换、应变花分析。

### 第8章 特征值与特征向量

#### 8.1 特征值与特征多项式

**定义** 对 $A\in\mathbb{C}^{n\times n}$，若存在 $\mathbf{v}\ne\mathbf{0}$ 与 $\lambda$ 使
$$
A\mathbf{v}=\lambda\mathbf{v},
$$
则 $\lambda$ 为特征值，$\mathbf{v}$ 为特征向量。

**公式** $\lambda$ 是特征值 $\iff \det(\lambda I-A)=0$。特征多项式
$$
p_A(\lambda)=\det(\lambda I-A)=\lambda^n-c_1\lambda^{n-1}+\cdots+(-1)^n c_n,
$$
其中 $c_1=\operatorname{tr}A$，$c_n=\det A$。

**推导** $A\mathbf{v}=\lambda\mathbf{v}\iff(\lambda I-A)\mathbf{v}=\mathbf{0}$ 有非零解 $\iff \lambda I-A$ 奇异 $\iff\det(\lambda I-A)=0$。展开 det 得首一多项式，由 Vieta 定理，根之和 $=c_1=\operatorname{tr}A$，根之积 $=(-1)^n c_n=\det A$。

**意义** 特征值是线性变换在不变方向上的伸缩率；特征向量是仅被伸缩的方向。

**特例** 三角阵特征值即对角元；$A$ 与 $A^{\mathsf T}$ 有相同特征多项式。

**应用** 结构固有频率由 $\det(K-\omega^2M)=0$ 给出（广义特征值问题）。

#### 8.2 特征子空间与相似不变性

**定义** $E_\lambda=\mathcal{N}(\lambda I-A)=\{\mathbf{v}:A\mathbf{v}=\lambda\mathbf{v}\}\cup\{\mathbf{0}\}$。

**推导** $E_\lambda$ 是子空间（核空间）；不同特征值的特征向量线性无关：设 $\mathbf{v}_1,\dots,\mathbf{v}_k$ 对应互异 $\lambda_i$，用 $A$ 反复作用并消元，可归纳证明其无关。

**相似不变性** 若 $A'=P^{-1}AP$，则
$$
\det(\lambda I-A')=\det\big(P^{-1}(\lambda I-A)P\big)=\det(\lambda I-A),
$$
故特征多项式、特征值、迹、行列式均相似不变。

**意义** 特征值刻画变换的内禀性质。

**特例** 实矩阵的复特征值成共轭对出现。

**应用** 稳定性判据：所有特征值实部 $<0$ 则渐近稳定。

#### 8.3 特征值的代数与几何重数

**定义** $\lambda_0$ 的代数重数 $m_a$ 为特征多项式根的重数；几何重数 $m_g=\dim E_{\lambda_0}$。

**定理** $1\le m_g\le m_a$。

**推导** $m_g\ge1$ 显然。取 $E_\lambda$ 的基扩充为全空间基，则 $A$ 在该基下为 $\begin{pmatrix}\lambda I&B\\0&C\end{pmatrix}$，故
$$
\det(\mu I-A)=(\mu-\lambda)^{m_g}\det(\mu I-C),
$$
$\therefore m_a\ge m_g$。

**意义** 代数重数是“应该”有的特征向量个数，几何重数是“真正”有的；差异导致不可对角化。

**特例** $m_g=m_a$ 对所有特征值成立 $\iff A$ 可对角化。

**应用** 重根模态对应的耦合振动；缺陷系统出现若尔当链。

#### 8.4 凯莱–哈密顿定理

**定理** 方阵满足其特征多项式：$p_A(A)=0$。

**推导**（借助伴随矩阵） 对任意 $\lambda$，$(\lambda I-A)\operatorname{adj}(\lambda I-A)=\det(\lambda I-A)I=p_A(\lambda)I$。伴随矩阵元素是 $\lambda$ 的至多 $n-1$ 次多项式，设
$$
\operatorname{adj}(\lambda I-A)=\sum_{k=0}^{n-1}\lambda^k B_k,\qquad p_A(\lambda)=\sum_{k=0}^{n}a_k\lambda^k\ (a_n=1).
$$
比较 $(\lambda I-A)\sum B_k\lambda^k=\sum a_k\lambda^k I$ 两边 $\lambda$ 的各次幂系数：
$$
B_{n-1}=I,\quad B_{k-1}-AB_k=a_kI,\quad -AB_0=a_0I.
$$
依次左乘 $A^k$ 并相加，中间项相消，得 $\sum_{k=0}^{n}a_kA^k=0$，即 $p_A(A)=0$。

**意义** 高次幂可由低次幂线性表示，从而 $A^{-1}$、$A^k$ 可化为 $A$ 的多项式。

**特例** 若 $A$ 可逆，由 $p_A(A)=0$ 得 $A^{-1}=-\frac{1}{a_0}(A^{n-1}+a_{n-1}A^{n-2}+\cdots+a_1I)$。

**应用** 快速计算矩阵幂、离散系统演化。

#### 8.5 特征值与迹、行列式、秩

**公式** $\operatorname{tr}A=\sum_i\lambda_i$，$\det A=\prod_i\lambda_i$，$\operatorname{rank}A\ge$ 非零特征值个数；$A$ 可逆 $\iff$ 无零特征值。

**推导** 由 $p_A(\lambda)=\prod_i(\lambda-\lambda_i)$ 比较系数即得。

**意义** 用少量标量快速估计整体性质。

**特例** 幂零矩阵所有特征值为 0。

**应用** 估计结构刚度矩阵的奇异性（零特征值对应刚体模态）。

#### 8.6 应用：主应力与固有频率

**主应力** 求 $\boldsymbol{\sigma}\mathbf{n}=\sigma\mathbf{n}$，即 $\det(\boldsymbol{\sigma}-\sigma I)=0$。平面应力时
$$
\sigma_{1,2}=\frac{\sigma_x+\sigma_y}{2}\pm\sqrt{\Big(\frac{\sigma_x-\sigma_y}{2}\Big)^2+\tau_{xy}^2}.
$$
**固有频率** 单自由度弹簧-质量：$m\ddot u+ku=0\Rightarrow \omega=\sqrt{k/m}$。多自由度：
$$
M\ddot{\mathbf{u}}+K\mathbf{u}=\mathbf{0},\qquad \det(K-\omega^2M)=0.
$$
**推广** 设 $\mathbf{u}=\boldsymbol{\phi}e^{i\omega t}$ 代入得 $(K-\omega^2M)\boldsymbol{\phi}=\mathbf{0}$，非零 $\boldsymbol{\phi}$ 存在要求 $K-\omega^2M$ 奇异，$\therefore$ 广义特征值问题 $K\boldsymbol{\phi}=\omega^2M\boldsymbol{\phi}$，$\omega^2$ 为广义特征值。

**意义** 主应力与主频率均为不变量。

**特例** $M=I$ 时退化为标准特征值问题。

**应用** 机翼颤振、航天器模态试验。

### 第9章 相似对角化与若尔当标准形

#### 9.1 可对角化条件

**定理** $A$ 可对角化 $\iff$ 存在 $n$ 个线性无关的特征向量 $\iff$ 每个特征值的几何重数等于代数重数。

**推导** $A=PDP^{-1}\iff AP=PD$。设 $P=[\mathbf{v}_1,\dots,\mathbf{v}_n]$，$D=\operatorname{diag}(\lambda_1,\dots,\lambda_n)$，则 $A\mathbf{v}_i=\lambda_i\mathbf{v}_i$，$\therefore P$ 的列必须是 $A$ 的特征向量；可逆即线性无关。重数条件由 §8.3 的 $m_g\le m_a$ 与 $\sum m_g=n$ 给出。

**意义** 对角化把耦合系统解耦；对角元即各方向伸缩率。

**特例** $n$ 个互异特征值必可对角化；实对称矩阵必可正交对角化（§11）。

**应用** 振型叠加法：$\mathbf{u}=P\boldsymbol{\eta}$ 使 $M,K$ 同时对角化（当 $M=I$ 或已做质量归一化）。

#### 9.2 相似标准形与若尔当块

**定义** 若尔当块
$$
J_k(\lambda)=\begin{pmatrix}\lambda&1&&\\&\lambda&\ddots&\\&&\ddots&1\\&&&\lambda\end{pmatrix}_{k\times k}.
$$
**定理（若尔当标准形）** 任意复方阵 $A$ 相似于分块对角阵 $J=\operatorname{diag}(J_{k_1}(\lambda_1),\dots)$，除去块排列外唯一。

**构造** 对广义特征向量链 $\mathbf{v}_1,\dots,\mathbf{v}_k$ 满足
$$
(A-\lambda I)\mathbf{v}_1=\mathbf{0},\quad (A-\lambda I)\mathbf{v}_{i+1}=\mathbf{v}_i,
$$
则 $A\mathbf{v}_{i}=\lambda\mathbf{v}_i+\mathbf{v}_{i-1}$（约定 $\mathbf{v}_0=0$），在这些向量构成的基下 $A$ 表现为 $J_k(\lambda)$。

**意义** 若不可对角化，若尔当形是最接近对角的形式，刻画特征值的“缺陷”结构。

**特例** 若所有 $k_i=1$ 则为对角阵。

**应用** 具有重根且阻尼耦合的振动系统；线性化动力系统的非半单分岔。

#### 9.3 最小多项式

**定义** 使 $m(A)=0$ 的最低次首一多项式 $m_A$ 为最小多项式。

**公式** $m_A$ 的根恰为 $A$ 的特征值（不重复）；$A$ 可对角化 $\iff m_A$ 无重根。

**推导** 由凯莱–哈密顿，$m_A\mid p_A$。若 $\lambda$ 为特征值、$\mathbf{v}$ 为特征向量，$m_A(A)\mathbf{v}=m_A(\lambda)\mathbf{v}=\mathbf{0}$，$\therefore m_A(\lambda)=0$。若 $m_A$ 无重根且整除 $p_A$，则每个特征值的若尔当块必为 $1\times1$（否则 $(A-\lambda I)^2$ 项需要），$\therefore$ 可对角化；反之显然。

**意义** 最小多项式给出使 $A$ 满足的最短代数关系。

**特例** 幂零阵 $N$ 的最小多项式为 $\lambda^k$，$k$ 为最大若尔当块尺寸。

**应用** 矩阵函数的定义域（§9.5）。

#### 9.4 矩阵函数

**定义** 若 $f$ 在 $A$ 的谱上解析，定义
$$
f(A)=\sum_{k=0}^{\infty}c_kA^k,\qquad f(t)=\sum c_kt^k,
$$
收敛由谱半径保证。

**公式（若尔当形式）** 对 $J_k(\lambda)$，用 Taylor 展开
$$
f(J_k(\lambda))=\begin{pmatrix}
f(\lambda)&f'(\lambda)&\frac{f''(\lambda)}{2!}&\cdots\\
&f(\lambda)&f'(\lambda)&\ddots\\
&&\ddots&\ddots\\
&&&f(\lambda)
\end{pmatrix}.
$$
**推导** $J_k(\lambda)=\lambda I+N$，$N^k=0$，故
$$
f(\lambda I+N)=\sum_{r=0}^{\infty}\frac{f^{(r)}(\lambda)}{r!}N^r=\sum_{r=0}^{k-1}\frac{f^{(r)}(\lambda)}{r!}N^r,
$$
$N^r$ 恰为第 $r$ 条上对角线为 1 的矩阵，$\therefore$ 得所示上三角 Toeplitz 结构。

**意义** 把标量函数解析延拓到矩阵。

**特例** $e^{At}$、$\sqrt{A}$、$\log A$。

**应用** 线性 ODE 解 $\mathbf{x}(t)=e^{At}\mathbf{x}_0$；材料本构的分数阶算子。

#### 9.5 应用：线性动力系统解耦

**系统** $\dot{\mathbf{x}}=A\mathbf{x}$。若 $A=PDP^{-1}$，令 $\mathbf{x}=P\mathbf{y}$，则
$$
P\dot{\mathbf{y}}=AP\mathbf{y}\ \Rightarrow\ \dot{\mathbf{y}}=P^{-1}AP\mathbf{y}=D\mathbf{y},
$$
即 $\dot y_i=\lambda_iy_i$，$y_i(t)=y_i(0)e^{\lambda_it}$，$\therefore \mathbf{x}(t)=Pe^{Dt}P^{-1}\mathbf{x}_0=e^{At}\mathbf{x}_0$。

**稳定性** $\operatorname{Re}\lambda_i<0\ \forall i\Rightarrow$ 渐近稳定；某 $\operatorname{Re}\lambda_i>0\Rightarrow$ 不稳定；纯虚且半单 $\Rightarrow$ 中性稳定。

**应用** 航天器姿态动力学线性化、结构主动控制。

### 第10章 内积空间与正交性

#### 10.1 内积与范数

**定义** 实向量空间内积 $\langle\cdot,\cdot\rangle:V\times V\to\mathbb{R}$ 满足对称、双线性、正定。范数 $\|\mathbf{v}\|=\sqrt{\langle\mathbf{v},\mathbf{v}\rangle}$。

**公式（柯西–施瓦茨）** $|\langle\mathbf{u},\mathbf{v}\rangle|\le\|\mathbf{u}\|\,\|\mathbf{v}\|$，等号 $\iff\mathbf{u},\mathbf{v}$ 线性相关。

**推导** 对任意 $t\in\mathbb{R}$，$0\le\|\mathbf{u}+t\mathbf{v}\|^2=\|\mathbf{u}\|^2+2t\langle\mathbf{u},\mathbf{v}\rangle+t^2\|\mathbf{v}\|^2$。若 $\mathbf{v}\ne\mathbf{0}$，取 $t^*=-\langle\mathbf{u},\mathbf{v}\rangle/\|\mathbf{v}\|^2$ 使二次式最小，得
$$
0\le\|\mathbf{u}\|^2-\frac{\langle\mathbf{u},\mathbf{v}\rangle^2}{\|\mathbf{v}\|^2}\ \Rightarrow\ \langle\mathbf{u},\mathbf{v}\rangle^2\le\|\mathbf{u}\|^2\|\mathbf{v}\|^2.
$$
开方即得；等号 $\iff\mathbf{u}+t^*\mathbf{v}=\mathbf{0}$。

**意义** 定义夹角 $\cos\theta=\frac{\langle\mathbf{u},\mathbf{v}\rangle}{\|\mathbf{u}\|\|\mathbf{v}\|}$。

**特例** $\mathbb{R}^n$ 标准内积；函数空间 $\langle f,g\rangle=\int_a^b fg\,dx$。

**应用** 能量范数 $\|\mathbf{u}\|_K=\sqrt{\mathbf{u}^{\mathsf T}K\mathbf{u}}$；功 $\int\mathbf{f}\cdot d\mathbf{u}$ 为内积。

#### 10.2 正交与正交补

**定义** $\mathbf{u}\perp\mathbf{v}\iff\langle\mathbf{u},\mathbf{v}\rangle=0$；$W^\perp=\{\mathbf{v}:\langle\mathbf{v},\mathbf{w}\rangle=0\ \forall\mathbf{w}\in W\}$。

**定理** $\dim W+\dim W^\perp=\dim V$，且 $V=W\oplus W^\perp$。

**推导** 取 $W$ 的标准正交基 $\mathbf{q}_1,\dots,\mathbf{q}_k$，扩充为 $V$ 的标准正交基 $\mathbf{q}_1,\dots,\mathbf{q}_n$，则 $W^\perp=\operatorname{span}\{\mathbf{q}_{k+1},\dots,\mathbf{q}_n\}$，$\therefore$ 维数互补且直和分解成立。

**意义** 正交投影是把向量分为“在内”与“在外”两部分。

**特例** $\mathbb{R}^3$ 中过原点平面 $W$ 的正交补是法线。

**应用** 应力分解为球量与偏量；信号分解到正交子空间。

#### 10.3 Gram–Schmidt 正交化

**算法** 给定线性无关 $\mathbf{a}_1,\dots,\mathbf{a}_k$：
$$
\mathbf{q}_1=\frac{\mathbf{a}_1}{\|\mathbf{a}_1\|},\qquad
\mathbf{v}_i=\mathbf{a}_i-\sum_{j=1}^{i-1}\langle\mathbf{a}_i,\mathbf{q}_j\rangle\mathbf{q}_j,\qquad
\mathbf{q}_i=\frac{\mathbf{v}_i}{\|\mathbf{v}_i\|}.
$$
**推导** 归纳证明 $\langle\mathbf{v}_i,\mathbf{q}_m\rangle=\langle\mathbf{a}_i,\mathbf{q}_m\rangle-\sum_{j}\langle\mathbf{a}_i,\mathbf{q}_j\rangle\delta_{jm}=0$（$m<i$），且 $\mathbf{v}_i\ne\mathbf{0}$（否则 $\mathbf{a}_i$ 由前 $i-1$ 个线性表示，矛盾），$\therefore$ $\mathbf{q}_i$ 规范正交。

**意义** 构造标准正交基；数值实现为 QR 分解。

**特例** 两个向量：$\mathbf{q}_2\propto\mathbf{a}_2-\frac{\langle\mathbf{a}_2,\mathbf{a}_1\rangle}{\langle\mathbf{a}_1,\mathbf{a}_1\rangle}\mathbf{a}_1$。

**应用** 模态正交化、多项式正交基（Legendre、Chebyshev）。

#### 10.4 QR 分解

**公式** 列满秩 $A\in\mathbb{R}^{m\times n}$ 可分解为 $A=QR$，$Q$ 列标准正交，$R$ 上三角可逆。

**推导** 对 $A$ 的列作 Gram–Schmidt：$\mathbf{a}_j=\sum_{i=1}^{j}r_{ij}\mathbf{q}_i$，其中 $r_{ij}=\langle\mathbf{a}_j,\mathbf{q}_i\rangle$（$i<j$），$r_{jj}=\|\mathbf{v}_j\|>0$。写成矩阵即 $A=QR$。

**意义** QR 把列空间的正交基与坐标耦合分开。

**特例** 方阵可逆时 $Q$ 正交，$R$ 上三角。

**应用** 最小二乘、特征值 QR 迭代。

#### 10.5 正交投影与最小二乘

**定义** 到子空间 $W$ 的投影 $P_W$，满足 $P_W^2=P_W$，$P_W^{\mathsf T}=P_W$。

**公式** 若 $Q$（$m\times k$，列标准正交）为 $W$ 的基，则
$$
P_W=QQ^{\mathsf T}.
$$
**推导** 对任意 $\mathbf{x}$，$QQ^{\mathsf T}\mathbf{x}\in W$，且 $\mathbf{x}-QQ^{\mathsf T}\mathbf{x}\perp W$，因为对每个基列 $\mathbf{q}_j$，$\mathbf{q}_j^{\mathsf T}(\mathbf{x}-QQ^{\mathsf T}\mathbf{x})=\mathbf{q}_j^{\mathsf T}\mathbf{x}-\mathbf{q}_j^{\mathsf T}\mathbf{x}=0$，$\therefore$ 是正交投影。

**最小二乘** 解 $\min_{\mathbf{x}}\|A\mathbf{x}-\mathbf{b}\|$ 等价于法方程
$$
A^{\mathsf T}A\mathbf{x}=A^{\mathsf T}\mathbf{b}.
$$
**推导** 残差 $\mathbf{r}=\mathbf{b}-A\mathbf{x}$ 在最优时垂直于 $\mathcal{C}(A)$，即 $A^{\mathsf T}\mathbf{r}=\mathbf{0}$，$\therefore A^{\mathsf T}A\mathbf{x}=A^{\mathsf T}\mathbf{b}$。

**意义** 无解时求最优近似。

**特例** $A$ 列满秩时 $\mathbf{x}=(A^{\mathsf T}A)^{-1}A^{\mathsf T}\mathbf{b}$。

**应用** 数据拟合、参数反演、有限元投影。

#### 10.6 应用：函数正交基与能量投影

**函数内积** $\langle u,v\rangle=\int_\Omega uv\,d\Omega$。Fourier 级数：$u=\sum_k c_k\phi_k$，$\phi_k$ 正交，
$$
c_k=\frac{\langle u,\phi_k\rangle}{\langle\phi_k,\phi_k\rangle}.
$$
**推导** 两边与 $\phi_j$ 内积，正交性消去其余项，$\therefore$ 得系数。

**应用** Galerkin/加权残值法：把偏微分方程残差投影到有限维子空间使其正交，得到有限元方程。

### 第11章 实对称矩阵与二次型

#### 11.1 二次型与矩阵表示

**定义** 二次型 $Q(\mathbf{x})=\sum_{i,j}a_{ij}x_ix_j=\mathbf{x}^{\mathsf T}A\mathbf{x}$，可取 $A$ 对称（因 $x_ix_j$ 与 $x_jx_i$ 合并）。

**推导** 任意 $A$ 的二次型只依赖对称部分：$\mathbf{x}^{\mathsf T}A\mathbf{x}=\mathbf{x}^{\mathsf T}\frac{A+A^{\mathsf T}}{2}\mathbf{x}$，$\therefore$ 可设 $A$ 对称。

**意义** 二次型是能量、度量、方差的代数表达。

**特例** $\|\mathbf{x}\|^2=\mathbf{x}^{\mathsf T}\mathbf{x}$。

**应用** 应变能 $\frac12\mathbf{u}^{\mathsf T}K\mathbf{u}$、动能 $\frac12\dot{\mathbf{u}}^{\mathsf T}M\dot{\mathbf{u}}$。

#### 11.2 实对称矩阵的谱定理

**定理** 实对称矩阵 $A$ 的特征值全实，且存在正交矩阵 $Q$ 使
$$
A=Q\Lambda Q^{\mathsf T},\qquad \Lambda=\operatorname{diag}(\lambda_1,\dots,\lambda_n).
$$
**推导** (1) 特征值实：设 $A\mathbf{v}=\lambda\mathbf{v}$，$\mathbf{v}\ne\mathbf{0}$，则
$$
\lambda\|\mathbf{v}\|^2=\mathbf{v}^*A\mathbf{v}=\overline{\mathbf{v}^{\mathsf T}A\mathbf{v}}=\overline{\lambda}\|\mathbf{v}\|^2,
$$
$\therefore \lambda=\overline{\lambda}\in\mathbb{R}$。(2) 不同特征值的特征向量正交：$A\mathbf{u}=\lambda\mathbf{u},A\mathbf{v}=\mu\mathbf{v},\lambda\ne\mu$，
$$
(\lambda-\mu)\mathbf{u}^{\mathsf T}\mathbf{v}=(A\mathbf{u})^{\mathsf T}\mathbf{v}-\mathbf{u}^{\mathsf T}A\mathbf{v}=0\Rightarrow\mathbf{u}\perp\mathbf{v}.
$$
(3) 用数学归纳法或 Schur 分解：存在单位特征向量 $\mathbf{q}_1$，将 $A$ 限制在 $\mathbf{q}_1^\perp$（不变子空间）上仍是实对称，递归得到 $n$ 个标准正交特征向量，拼成 $Q$，$\therefore A=Q\Lambda Q^{\mathsf T}$。

**意义** 实对称矩阵总可正交对角化，主方向互相垂直。

**特例** $n=2$ 时正交矩阵 $Q=\begin{pmatrix}\cos\theta&-\sin\theta\\ \sin\theta&\cos\theta\end{pmatrix}$，$\theta$ 为主方向角。

**应用** 主轴定理：主应力、主应变、主惯性轴。

#### 11.3 正定性与判据

**定义** $A$ 正定（$A\succ0$）若 $\mathbf{x}^{\mathsf T}A\mathbf{x}>0\ \forall\mathbf{x}\ne\mathbf{0}$；半正定 $\succeq0$ 同理非负。

**判据** 实对称 $A\succ0$ $\iff$ 所有特征值 $>0$ $\iff$ 所有顺序主子式 $>0$（Sylvester）$\iff$ 存在可逆 $R$ 使 $A=R^{\mathsf T}R$（Cholesky）。

**推导**（特征值判据） 由谱定理 $\mathbf{x}^{\mathsf T}A\mathbf{x}=\sum_i\lambda_i y_i^2$（$\mathbf{y}=Q^{\mathsf T}\mathbf{x}$），$y$ 可取任意向量，故正定 $\iff\lambda_i>0$。（$R^{\mathsf T}R$）若 $A=R^{\mathsf T}R$ 且 $R$ 可逆，则 $\mathbf{x}^{\mathsf T}A\mathbf{x}=\|R\mathbf{x}\|^2>0$。反之由 Cholesky 分解。

**意义** 正定是能量函数有下界、存在最小点的代数条件。

**特例** 单位阵正定；$K$ 正定表示结构无刚体模态。

**应用** 刚度矩阵半正定、质量矩阵正定；弹性张量正定保证能量凸性。

#### 11.4 合同变换与惯性定律

**定义** 若存在可逆 $C$ 使 $B=C^{\mathsf T}AC$，称 $A,B$ 合同。

**定理（Sylvester 惯性定律）** 二次型经合同变换，正、负、零特征值（惯性指数）个数不变。

**推导** 由谱定理 $A=Q\Lambda Q^{\mathsf T}$，令 $C=Q\Lambda^{-1/2}$ 可将 $A$ 合同于 $\operatorname{diag}(I_p,-I_q,0)$，$p,q$ 为正负惯性指数。若两种对角化给出不同 $(p,q)$，在交点处由正定性产生矛盾（对正部零化后 $x^{\mathsf T}Ax>0$ 与 $<0$ 同时成立），$\therefore$ 唯一。

**意义** 二次型可按惯性分类，与坐标选择无关。

**特例** 秩 $r=p+q$。

**应用** 应力张量符号分类（拉/压/混合）、能量正负定判定稳定性。

#### 11.5 瑞利商与极值刻画

**定义** 瑞利商
$$
R(\mathbf{x})=\frac{\mathbf{x}^{\mathsf T}A\mathbf{x}}{\mathbf{x}^{\mathsf T}\mathbf{x}},\qquad\mathbf{x}\ne\mathbf{0}.
$$
**定理** $A$ 对称，则 $\lambda_{\min}\le R(\mathbf{x})\le\lambda_{\max}$，且
$$
\lambda_{\max}=\max_{\mathbf{x}\ne\mathbf{0}}R(\mathbf{x}),\qquad\lambda_{\min}=\min_{\mathbf{x}\ne\mathbf{0}}R(\mathbf{x}).
$$
**推导** 令 $\mathbf{x}=Q\mathbf{y}$，则 $R=\frac{\sum\lambda_iy_i^2}{\sum y_i^2}$，是 $\lambda_i$ 的加权平均，$\therefore$ 介于最小与最大特征值之间；取 $\mathbf{y}=\mathbf{e}_i$ 得等号。

**意义** 特征值是二次型在单位球上的极值，梯度 $\nabla R$ 在特征向量处为零。

**推广** $\lambda_k=\max\{\min R(\mathbf{x}):\mathbf{x}\in U,\dim U=k\}$（Courant–Fischer）。

**应用** 结构稳定性的能量判据、主应力估计、条件数 $\kappa=\lambda_{\max}/\lambda_{\min}$。

#### 11.6 应用：应变能与稳定性

**应变能** 线弹性
$$
U=\frac12\boldsymbol{\varepsilon}^{\mathsf T}\mathbf{C}\boldsymbol{\varepsilon}\ge0,
$$
$\mathbf{C}$ 对称半正定；$\mathbf{C}\succ0$ 保证材料稳定（Drucker 稳定性）。

**势能极小** 平衡 $K\mathbf{u}=\mathbf{f}$ 是势能
$$
\Pi(\mathbf{u})=\frac12\mathbf{u}^{\mathsf T}K\mathbf{u}-\mathbf{f}^{\mathsf T}\mathbf{u}
$$
的驻点：$\nabla\Pi=K\mathbf{u}-\mathbf{f}=\mathbf{0}$。$K\succ0$ 保证唯一极小。**推导** 对 $\Pi$ 求梯度并用 $K$ 对称。

**稳定性** 保守系统平衡稳定 $\iff$ 势能二阶变分正定，即切线刚度阵 $K_T\succ0$。当最小特征值过零时失稳（分岔）。

**应用** 欧拉压杆屈曲：$P_{cr}=\lambda_{\min}$ 对应广义问题 $K\phi=P_{cr}K_G\phi$。

### 第12章 奇异值分解

#### 12.1 Gram 矩阵与奇异值

**定义** 非零 $A\in\mathbb{R}^{m\times n}$ 的奇异值为 $A^{\mathsf T}A$（或 $AA^{\mathsf T}$）特征值的非负平方根 $\sigma_i=\sqrt{\lambda_i(A^{\mathsf T}A)}$。

**公式** $A^{\mathsf T}A$ 对称半正定，特征值 $\lambda_1\ge\cdots\ge\lambda_r>0=\cdots=0$，$\sigma_i=\sqrt{\lambda_i}$。

**推导** 对任意 $\mathbf{x}$，$\mathbf{x}^{\mathsf T}A^{\mathsf T}A\mathbf{x}=\|A\mathbf{x}\|^2\ge0$，$\therefore$ 半正定；特征值非负，$\therefore$ 平方根有意义。

**意义** 奇异值度量 $A$ 在各正交方向上的伸缩量。

**特例** 方阵可逆时 $\sigma_i>0$；对称矩阵的奇异值 $=|\lambda_i|$。

**应用** 条件数 $\kappa_2(A)=\sigma_{\max}/\sigma_{\min}$。

#### 12.2 奇异值分解定理

**定理** 任意 $A\in\mathbb{R}^{m\times n}$ 可分解为
$$
A=U\Sigma V^{\mathsf T},\qquad U\in\mathbb{R}^{m\times m},\ V\in\mathbb{R}^{n\times n}\text{ 正交},
$$
$\Sigma$ 为对角（奇异值降序排列）。

**推导** $A^{\mathsf T}A$ 对称半正定，由谱定理有标准正交特征向量 $\mathbf{v}_1,\dots,\mathbf{v}_n$，$A^{\mathsf T}A\mathbf{v}_i=\sigma_i^2\mathbf{v}_i$。对 $\sigma_i>0$ 定义
$$
\mathbf{u}_i=\frac{A\mathbf{v}_i}{\sigma_i},
$$
则 $\|\mathbf{u}_i\|=1$ 且
$$
\mathbf{u}_i^{\mathsf T}\mathbf{u}_j=\frac{\mathbf{v}_i^{\mathsf T}A^{\mathsf T}A\mathbf{v}_j}{\sigma_i\sigma_j}=\frac{\sigma_j^2\mathbf{v}_i^{\mathsf T}\mathbf{v}_j}{\sigma_i\sigma_j}=\delta_{ij}.
$$
把 $\{\mathbf{u}_i\}$ 扩充为 $\mathbb{R}^m$ 的标准正交基。于是对任意 $\mathbf{v}_i$，$A\mathbf{v}_i=\sigma_i\mathbf{u}_i$（$\sigma_i=0$ 时两边为零），写成矩阵即 $AV=U\Sigma$，$\therefore A=U\Sigma V^{\mathsf T}$。

**意义** 任意线性变换可分解为：右正交旋转 $\to$ 沿坐标轴伸缩 $\to$ 左正交旋转。

**特例** $A$ 对称正定时 $U=V$，SVD 退化为特征分解。

**应用** 图像压缩、本征正交分解（POD）、降阶模型。

#### 12.3 低秩近似（Eckart–Young）

**定理** 在 $\operatorname{rank}B\le k$ 的矩阵中，$\|A-A_k\|_2$ 的最小值由
$$
A_k=\sum_{i=1}^{k}\sigma_i\mathbf{u}_i\mathbf{v}_i^{\mathsf T}
$$
达到，最小值为 $\sigma_{k+1}$；Frobenius 范数意义下最小值为 $\sqrt{\sum_{i>k}\sigma_i^2}$。

**推导** 由 $A=U\Sigma V^{\mathsf T}$，$\|A-B\|_2=\|U\Sigma V^{\mathsf T}-B\|_2=\|\Sigma-U^{\mathsf T}BV\|_2$。对任意秩 $\le k$ 的 $C=U^{\mathsf T}BV$，其前 $k$ 行零化后第 $k+1$ 个奇异值下界由 $\sigma_{k+1}$ 给出（用秩不等式），取 $A_k$ 时达到。

**意义** 用少数主方向最优逼近原矩阵。

**特例** $k=r=\operatorname{rank}A$ 时精确。

**应用** 主成分分析、图像压缩、数据去噪。

#### 12.4 伪逆与最小二乘

**公式** 设 $A=U\Sigma V^{\mathsf T}$，$A^+=V\Sigma^+U^{\mathsf T}$，其中 $(\Sigma^+)_{ii}=1/\sigma_i$（$\sigma_i>0$），否则 $0$。

**定理** $\mathbf{x}=A^+\mathbf{b}$ 是所有最小二乘解中范数最小者。

**推导** 令 $\mathbf{c}=U^{\mathsf T}\mathbf{b}$，$\mathbf{z}=V^{\mathsf T}\mathbf{x}$，则 $\|A\mathbf{x}-\mathbf{b}\|^2=\|\Sigma\mathbf{z}-\mathbf{c}\|^2$。对 $\sigma_i>0$ 取 $z_i=c_i/\sigma_i$ 使该项为零；对 $\sigma_i=0$ 需 $c_i=0$ 才相容，否则残差固定，取 $z_i=0$ 使 $\|\mathbf{z}\|$ 最小，$\therefore \mathbf{x}=V\Sigma^+U^{\mathsf T}\mathbf{b}$。

**应用** 病态系统正则化、全最小二乘。

#### 12.5 应用：降阶与主成分

**POD/模态降阶** 设快照矩阵 $X=[\mathbf{x}_1,\dots,\mathbf{x}_N]$，其 SVD 的 $\mathbf{u}_i$ 为最优正交模态基，投影得降阶方程。

**主成分** 数据中心化后协方差 $C=\frac{1}{N}XX^{\mathsf T}$，其特征分解等价于 $X$ 的 SVD；主方向为 $\mathbf{u}_i$。

**意义** 用低维子空间捕捉主要能量。

**特例** 秩 $k$ 近似保留能量比 $\sum_{i\le k}\sigma_i^2/\sum_i\sigma_i^2$。

**应用** 流场降阶、结构模型降阶、气动弹性分析。

### 第13章 矩阵微分方程与矩阵指数

#### 13.1 矩阵指数

**定义** $e^{At}=\sum_{k=0}^{\infty}\frac{A^kt^k}{k!}$。

**公式（性质）** $e^{A(t+s)}=e^{At}e^{As}$（同一 $A$），$\frac{d}{dt}e^{At}=Ae^{At}=e^{At}A$，$(e^{At})^{-1}=e^{-At}$。

**推导**（导数） 逐项求导
$$
\frac{d}{dt}e^{At}=\sum_{k=1}^{\infty}\frac{A^kt^{k-1}}{(k-1)!}=A\sum_{k=0}^{\infty}\frac{A^kt^k}{k!}=Ae^{At}.
$$
（可交换性）$AB=BA\Rightarrow e^{A+B}=e^Ae^B$，由 Cauchy 乘积与二项式定理。

**意义** 矩阵指数是标量指数的自然推广，是线性系统的时间演化算子。

**特例** 对角阵 $e^{\Lambda t}=\operatorname{diag}(e^{\lambda_it})$；幂零 $N$，$e^{Nt}=\sum_{k=0}^{n-1}\frac{N^kt^k}{k!}$。

**应用** 求解 $\dot{\mathbf{x}}=A\mathbf{x}$。

#### 13.2 线性常微分方程组的解

**公式** 初值问题 $\dot{\mathbf{x}}=A\mathbf{x}+\mathbf{f}(t)$，$\mathbf{x}(0)=\mathbf{x}_0$ 的解
$$
\mathbf{x}(t)=e^{At}\mathbf{x}_0+\int_0^t e^{A(t-s)}\mathbf{f}(s)\,ds.
$$
**推导** 齐次：$\mathbf{x}=e^{At}\mathbf{c}$，代入初值得 $\mathbf{c}=\mathbf{x}_0$。非齐次用常数变易 $\mathbf{x}=e^{At}\mathbf{c}(t)$，则
$$
\dot{\mathbf{x}}=Ae^{At}\mathbf{c}+e^{At}\dot{\mathbf{c}}=A\mathbf{x}+\mathbf{f}\Rightarrow e^{At}\dot{\mathbf{c}}=\mathbf{f}\Rightarrow \mathbf{c}(t)=\mathbf{x}_0+\int_0^t e^{-As}\mathbf{f}(s)\,ds,
$$
$\therefore \mathbf{x}=e^{At}\mathbf{x}_0+\int_0^t e^{A(t-s)}\mathbf{f}(s)\,ds$。

**意义** 卷积形式的叠加积分。

**特例** 若 $A$ 可对角化，$e^{At}=Pe^{\Lambda t}P^{-1}$。

**应用** 受迫振动、控制系统阶跃/脉冲响应。

#### 13.3 稳定性与李雅普诺夫

**定义** 平衡点 $\mathbf{x}=\mathbf{0}$ 渐近稳定，若所有 $\operatorname{Re}\lambda_i(A)<0$。

**李雅普诺夫方程** 若存在 $P\succ0$ 使
$$
A^{\mathsf T}P+PA=-Q\prec0,
$$
则系统渐近稳定。

**推导** 取 $V(\mathbf{x})=\mathbf{x}^{\mathsf T}P\mathbf{x}$，沿轨迹
$$
\dot V=\dot{\mathbf{x}}^{\mathsf T}P\mathbf{x}+\mathbf{x}^{\mathsf T}P\dot{\mathbf{x}}=\mathbf{x}^{\mathsf T}(A^{\mathsf T}P+PA)\mathbf{x}=-\mathbf{x}^{\mathsf T}Q\mathbf{x}<0,
$$
$\therefore V$ 单调递减到零。

**意义** 用二次型能量函数判定稳定性，无需显式求特征值。

**特例** 对任意 $Q\succ0$，$P=\int_0^\infty e^{A^{\mathsf T}t}Qe^{At}\,dt$ 是解（当稳定时收敛）。

**应用** 结构主动控制、飞控系统。

#### 13.4 二阶系统的矩阵形式

**无阻尼振动** $M\ddot{\mathbf{u}}+K\mathbf{u}=\mathbf{f}$。令 $\mathbf{y}=\begin{pmatrix}\mathbf{u}\\\dot{\mathbf{u}}\end{pmatrix}$，则
$$
\dot{\mathbf{y}}=\underbrace{\begin{pmatrix}0&I\\-M^{-1}K&0\end{pmatrix}}_{\mathcal{A}}\mathbf{y}+\begin{pmatrix}0\\M^{-1}\mathbf{f}\end{pmatrix}.
$$
**有阻尼** 加入 $C\dot{\mathbf{u}}$：
$$
\dot{\mathbf{y}}=\begin{pmatrix}0&I\\-M^{-1}K&-M^{-1}C\end{pmatrix}\mathbf{y}+\cdots.
$$
**特征方程** $\det(\lambda^2M+\lambda C+K)=0$。

**推导** 把二阶方程降阶为状态空间；$\mathcal{A}$ 的特征值 $\lambda$ 满足 $(M\lambda^2+C\lambda+K)\boldsymbol{\phi}=\mathbf{0}$。

**意义** 二阶系统可完全用矩阵谱分析。

**特例** 无阻尼时 $\lambda=\pm i\omega$，$\omega^2$ 为广义特征值。

**应用** 模态叠加、颤振分析。

#### 13.5 应用：控制系统的可控性与可观性

**可控性矩阵** $\mathcal{C}=[B,AB,A^2B,\dots,A^{n-1}B]$。系统可控 $\iff\operatorname{rank}\mathcal{C}=n$。

**可观性矩阵** $\mathcal{O}=[C^{\mathsf T},A^{\mathsf T}C^{\mathsf T},\dots,(A^{\mathsf T})^{n-1}C^{\mathsf T}]^{\mathsf T}$，可观 $\iff\operatorname{rank}\mathcal{O}=n$。

**推导（可控性）** 解 $e^{-At}B=e^{-At}Be^{-At}e^{At}$... 由 $\mathbf{x}(t)=e^{At}\mathbf{x}_0+\int_0^t e^{A(t-s)}B\mathbf{u}(s)ds$，把 $e^{A(t-s)}B$ 按凯莱–哈密顿展开成 $B,AB,\dots,A^{n-1}B$ 的线性组合，$\therefore$ 可达集维数 $=\operatorname{rank}\mathcal{C}$。

**应用** 航天器推力器配置、主动振动控制作动器布点。

### 第14章 力学与航空航天中的线性代数方法

#### 14.1 应力张量与主应力正交分解

**胡克定律** 各向同性线弹性
$$
\boldsymbol{\sigma}=\lambda\operatorname{tr}(\boldsymbol{\varepsilon})I+2\mu\boldsymbol{\varepsilon},\qquad \lambda=\frac{E\nu}{(1+\nu)(1-2\nu)},\ \mu=\frac{E}{2(1+\nu)}.
$$
**推导** 由应变能密度 $W=\frac12\lambda(\operatorname{tr}\boldsymbol{\varepsilon})^2+\mu\,\boldsymbol{\varepsilon}:\boldsymbol{\varepsilon}$ 对 $\boldsymbol{\varepsilon}$ 求导即得（对称张量求导）。

**主应力** 解 $\det(\boldsymbol{\sigma}-\sigma I)=0$ 得 $\sigma_1\ge\sigma_2\ge\sigma_3$，对应正交单位向量 $\mathbf{n}_i$，$\therefore\boldsymbol{\sigma}=Q\operatorname{diag}(\sigma_1,\sigma_2,\sigma_3)Q^{\mathsf T}$。

**应力不变量** $I_1=\operatorname{tr}\boldsymbol{\sigma}$，$I_2=\frac12[(\operatorname{tr}\boldsymbol{\sigma})^2-\operatorname{tr}(\boldsymbol{\sigma}^2)]$，$I_3=\det\boldsymbol{\sigma}$。

**意义** 主方向互相垂直，最大剪应力 $\tau_{\max}=(\sigma_1-\sigma_3)/2$。

**应用** 屈服准则 von Mises：$J_2=\frac12\operatorname{tr}(\mathbf{s}^2)=\text{const}$，$\mathbf{s}=\boldsymbol{\sigma}-\frac13I_1I$ 为偏应力。

#### 14.2 应变能与有限元刚度

**单元应变能** 位移插值 $\mathbf{u}=N\mathbf{d}$，应变 $\boldsymbol{\varepsilon}=B\mathbf{d}$，则
$$
U_e=\frac12\int_{V_e}\boldsymbol{\varepsilon}^{\mathsf T}\mathbf{C}\boldsymbol{\varepsilon}\,dV=\frac12\mathbf{d}^{\mathsf T}\Big(\int_{V_e}B^{\mathsf T}\mathbf{C}B\,dV\Big)\mathbf{d}=\frac12\mathbf{d}^{\mathsf T}K_e\mathbf{d},
$$
其中单元刚度阵 $K_e=\int_{V_e}B^{\mathsf T}\mathbf{C}B\,dV$ 对称半正定。

**总体方程** 组装 $K=\sum_e A_e^{\mathsf T}K_eA_e$，得 $K\mathbf{u}=\mathbf{f}$。

**推导** 势能 $\Pi=\sum_e U_e-\mathbf{f}^{\mathsf T}\mathbf{u}$，对 $\mathbf{u}$ 取驻值 $\nabla\Pi=0$ 得平衡方程。

**意义** 刚度阵由能量二次型导出，正定性对应弹性稳定。

**特例** 刚体模态使 $K\mathbf{1}=\mathbf{0}$，故 $K$ 仅半正定。

**应用** 飞机机翼、机身有限元建模。

#### 14.3 拉格朗日方程的线性化

**拉格朗日方程** 设广义坐标 $\mathbf{q}$，
$$
\frac{d}{dt}\Big(\frac{\partial L}{\partial\dot{\mathbf{q}}}\Big)-\frac{\partial L}{\partial\mathbf{q}}=Q,\qquad L=T-V.
$$
**动能** $T=\frac12\dot{\mathbf{q}}^{\mathsf T}M(\mathbf{q})\dot{\mathbf{q}}$，**势能** $V=V(\mathbf{q})$。

**线性化** 在平衡 $\mathbf{q}_0$（$\partial V/\partial\mathbf{q}|_0=Q$）处设 $\mathbf{q}=\mathbf{q}_0+\boldsymbol{\eta}$，展开到二阶：
$$
M\ddot{\boldsymbol{\eta}}+K\boldsymbol{\eta}=\mathbf{f},\qquad K=\frac{\partial^2V}{\partial\mathbf{q}\partial\mathbf{q}^{\mathsf T}}\Big|_{\mathbf{q}_0},\ M=M(\mathbf{q}_0).
$$
**推导** 计算
$$
\frac{\partial L}{\partial\dot{\mathbf{q}}}=M\dot{\mathbf{q}},\qquad \frac{d}{dt}=M\ddot{\mathbf{q}}+\dot M\dot{\mathbf{q}},
$$
$$
\frac{\partial L}{\partial\mathbf{q}}=\frac12\dot{\mathbf{q}}^{\mathsf T}\frac{\partial M}{\partial\mathbf{q}}\dot{\mathbf{q}}-\frac{\partial V}{\partial\mathbf{q}}.
$$
在平衡附近 $\dot M\dot{\mathbf{q}}$ 与速度二次项为高阶小量，保留一阶得 $M\ddot{\boldsymbol{\eta}}+\frac{\partial^2V}{\partial\mathbf{q}^2}\boldsymbol{\eta}=\mathbf{0}$。

**意义** 小振动问题完全由对称正定质量阵与对称刚度阵描述。

**特例** 单摆线性化 $\ddot\theta+\frac{g}{l}\theta=0$。

**应用** 卫星柔性附件振动、直升机旋翼动力学。

#### 14.4 哈密顿正则方程的辛结构

**哈密顿方程** 广义动量 $\mathbf{p}=\partial L/\partial\dot{\mathbf{q}}$，哈密顿量 $H=\mathbf{p}^{\mathsf T}\dot{\mathbf{q}}-L$，则
$$
\dot{\mathbf{q}}=\frac{\partial H}{\partial\mathbf{p}},\qquad \dot{\mathbf{p}}=-\frac{\partial H}{\partial\mathbf{q}}.
$$
**矩阵形式** 令 $\mathbf{z}=(\mathbf{q},\mathbf{p})^{\mathsf T}$，
$$
\dot{\mathbf{z}}=J\nabla H,\qquad J=\begin{pmatrix}0&I\\-I&0\end{pmatrix},\qquad J^{\mathsf T}=-J,\ J^2=-I.
$$
**推导** 由全微分 $dH=\dot{\mathbf{q}}^{\mathsf T}d\mathbf{p}-\frac{\partial L}{\partial\mathbf{q}}^{\mathsf T}d\mathbf{q}$，与 $\dot{\mathbf{p}}=\partial L/\partial\mathbf{q}$ 比较得正则方程。

**线性化** $H=\frac12\mathbf{z}^{\mathsf T}H_{zz}\mathbf{z}$，得 $\dot{\mathbf{z}}=JH_{zz}\mathbf{z}$，其特征值成对 $\pm i\omega$（辛谱）。

**意义** 哈密顿系统保持相空间体积（$\operatorname{tr}J=0$），特征值对称。

**特例** 简谐振子 $H=\frac{p^2}{2m}+\frac12kq^2$。

**应用** 航天器轨道摄动、辛数值积分。

#### 14.5 刚体惯量与欧拉方程

**角动量** $\mathbf{L}=J\boldsymbol{\omega}$，$J$ 对称正定。主轴系下 $J=\operatorname{diag}(J_1,J_2,J_3)$。

**欧拉方程**
$$
J_1\dot\omega_1-(J_2-J_3)\omega_2\omega_3=M_1,
$$
$$
J_2\dot\omega_2-(J_3-J_1)\omega_3\omega_1=M_2,
$$
$$
J_3\dot\omega_3-(J_1-J_2)\omega_1\omega_2=M_3.
$$
**推导** 由 $\dot{\mathbf{L}}=\mathbf{M}$，在体坐标下 $\dot{\mathbf{L}}=\dot{\mathbf{L}}_b+\boldsymbol{\omega}\times\mathbf{L}$，展开得各分量式。

**线性化与稳定性** 绕主轴 $\omega=(0,0,\Omega)$ 的平衡，令 $\boldsymbol{\omega}=(\delta\omega_1,\delta\omega_2,\Omega+\delta\omega_3)$，保留一阶得
$$
J_1\delta\dot\omega_1-(J_2-J_3)\Omega\,\delta\omega_2=0,\qquad J_2\delta\dot\omega_2-(J_3-J_1)\Omega\,\delta\omega_1=0.
$$
特征方程给出 $\lambda^2=\frac{(J_2-J_3)(J_3-J_1)\Omega^2}{J_1J_2}$。$\therefore$ 绕最大或最小惯量主轴旋转稳定，绕中间轴不稳定（网球拍定理）。

**应用** 卫星姿态稳定（自旋稳定、重力梯度稳定）。

#### 14.6 结构稳定性与特征值屈曲

**广义特征值问题** 线性屈曲
$$
K\boldsymbol{\phi}=\lambda K_G\boldsymbol{\phi},\qquad \lambda=P_{cr},
$$
$K$ 为弹性刚度阵，$K_G$ 为几何刚度阵。

**推导** 总势能二阶变分 $\delta^2\Pi=\frac12\boldsymbol{\phi}^{\mathsf T}(K-PK_G)\boldsymbol{\phi}$；失稳对应存在非零 $\boldsymbol{\phi}$ 使二阶变分为零，即 $(K-PK_G)\boldsymbol{\phi}=\mathbf{0}$，$\therefore$ 为广义特征值问题。最小正特征值即临界载荷。

**意义** 屈曲是刚度矩阵出现零特征值的分岔。

**特例** 欧拉压杆 $P_{cr}=\pi^2EI/(\mu L)^2$。

**应用** 机翼壁板屈曲、火箭壳体轴压稳定性。

#### 14.7 流体力学方程的线性化与稳定性

**Navier–Stokes 方程** 不可压缩
$$
\frac{\partial\mathbf{v}}{\partial t}+(\mathbf{v}\cdot\nabla)\mathbf{v}=-\frac1\rho\nabla p+\nu\nabla^2\mathbf{v},\qquad\nabla\cdot\mathbf{v}=0.
$$
**线性化** 设 $\mathbf{v}=\mathbf{U}+\mathbf{v}'$，$\mathbf{U}$ 为定常基流，忽略 $\mathbf{v}'\cdot\nabla\mathbf{v}'$：
$$
\frac{\partial\mathbf{v}'}{\partial t}+(\mathbf{U}\cdot\nabla)\mathbf{v}'+(\mathbf{v}'\cdot\nabla)\mathbf{U}=-\frac1\rho\nabla p'+\nu\nabla^2\mathbf{v}'.
$$
**模态分析** 令 $\mathbf{v}'=\hat{\mathbf{v}}e^{i(\mathbf{k}\cdot\mathbf{x}-\omega t)}$，得线性特征值问题
$$
\mathcal{L}\hat{\mathbf{v}}=\omega\hat{\mathbf{v}},
$$
其中 $\mathcal{L}$ 为线性算子。增长模态 $\operatorname{Im}\omega>0$ 对应失稳。

**推导** 把扰动代入并用连续性条件消去压力（投影到无散子空间），得到作用在无散向量场上（或其离散矩阵）的线性算子。

**意义** 失稳分析归结为线性算子/矩阵的谱。

**特例** 平行剪切流的 Orr–Sommerfeld 方程。

**应用** 边界层转捩预测、气动弹性颤振。

#### 14.8 正交分解与降阶模型

**本征正交分解（POD）** 采样矩阵 $X=[\mathbf{x}(t_1),\dots,\mathbf{x}(t_N)]$，SVD $X=U\Sigma V^{\mathsf T}$，取前 $k$ 个左奇异向量 $\Phi=[\mathbf{u}_1,\dots,\mathbf{u}_k]$ 为模态基。

**Galerkin 投影** 令 $\mathbf{x}\approx\Phi\mathbf{a}$，代入 $\dot{\mathbf{x}}=f(\mathbf{x})$，投影
$$
\dot{\mathbf{a}}=\Phi^{\mathsf T}f(\Phi\mathbf{a}),
$$
得到 $k$ 维降阶模型。

**推导** 要求残差与基正交：$\Phi^{\mathsf T}(\dot{\mathbf{x}}-f)=0$，$\therefore$ 得投影方程。

**意义** 利用 SVD 最优低秩性把高维动力学压到低维。

**特例** 线性系统 $f=A\mathbf{x}$ 得 $\dot{\mathbf{a}}=\Phi^{\mathsf T}A\Phi\,\mathbf{a}$。

**应用** 气动弹性降阶、流场快速仿真、数字孪生。

#### 14.9 结构模态正交性与叠加

**模态正交性** 对第 $r,s$ 阶模态 $\boldsymbol{\phi}_r,\boldsymbol{\phi}_s$（对应 $\omega_r^2,\omega_s^2$）：
$$
\boldsymbol{\phi}_r^{\mathsf T}M\boldsymbol{\phi}_s=0,\qquad\boldsymbol{\phi}_r^{\mathsf T}K\boldsymbol{\phi}_s=0\quad(r\ne s).
$$
**推导** 由 $K\boldsymbol{\phi}_r=\omega_r^2M\boldsymbol{\phi}_r$ 与 $K\boldsymbol{\phi}_s=\omega_s^2M\boldsymbol{\phi}_s$，左乘 $\boldsymbol{\phi}_s^{\mathsf T}$、$\boldsymbol{\phi}_r^{\mathsf T}$ 并利用 $M,K$ 对称相减：
$$
(\omega_r^2-\omega_s^2)\boldsymbol{\phi}_r^{\mathsf T}M\boldsymbol{\phi}_s=0\Rightarrow\boldsymbol{\phi}_r^{\mathsf T}M\boldsymbol{\phi}_s=0,
$$
进而 $\boldsymbol{\phi}_r^{\mathsf T}K\boldsymbol{\phi}_s=0$。

**模态叠加** 令 $\mathbf{u}=\Phi\boldsymbol{\eta}$，方程变为解耦
$$
\ddot\eta_r+\omega_r^2\eta_r=\frac{\boldsymbol{\phi}_r^{\mathsf T}\mathbf{f}}{\boldsymbol{\phi}_r^{\mathsf T}M\boldsymbol{\phi}_r}.
$$
**质量归一化** 常取 $\boldsymbol{\phi}_r^{\mathsf T}M\boldsymbol{\phi}_r=1$。

**意义** 用广义正交基把耦合系统化为独立振子。

**应用** 飞机地面共振试验、航天器模态综合。

#### 14.10 最小二乘与参数辨识

**问题** 已知输入输出数据，辨识参数 $\boldsymbol{\theta}$。由模型 $y_i=\mathbf{x}_i^{\mathsf T}\boldsymbol{\theta}+\epsilon_i$，堆叠
$$
\mathbf{y}=X\boldsymbol{\theta}+\boldsymbol{\epsilon},\qquad\hat{\boldsymbol{\theta}}=(X^{\mathsf T}X)^{-1}X^{\mathsf T}\mathbf{y}.
$$
**推导** 最小化 $\|\mathbf{y}-X\boldsymbol{\theta}\|^2$，法方程由 §10.5 给出；$X$ 列满秩时系数阵可逆。

**意义** 大量测量数据的统计平均，等价于投影到试验数据张成的空间。

**特例** 加权最小二乘 $\hat{\boldsymbol{\theta}}=(X^{\mathsf T}WX)^{-1}X^{\mathsf T}W\mathbf{y}$，$W\succ0$。

**应用** 气动参数辨识、结构损伤识别、系统辨识。

#### 14.11 条件数与数值稳定性

**条件数** $\kappa_2(A)=\|A\|_2\|A^{-1}\|_2=\sigma_{\max}/\sigma_{\min}$。

**误差放大** 解 $A\mathbf{x}=\mathbf{b}$ 的扰动满足
$$
\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|}\le\kappa_2(A)\frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}.
$$
**推导** 由 $\delta\mathbf{x}=A^{-1}\delta\mathbf{b}$，$\|\delta\mathbf{x}\|\le\|A^{-1}\|\|\delta\mathbf{b}\|$，又 $\|\mathbf{b}\|\le\|A\|\|\mathbf{x}\|$，两式相除即得。

**意义** 条件数衡量求解对误差的敏感度；$\kappa$ 大则病态。

**特例** 正交矩阵 $\kappa_2=1$（最佳条件）。

**应用** 有限元网格病态、正则化（Tikhonov，利用 SVD 滤波因子 $\frac{\sigma_i}{\sigma_i^2+\alpha}$）。

#### 14.12 综合：从连续介质到矩阵方程

**离散化链条** 连续场 $\to$ 试探/检验函数空间（基的选取）$\to$ 刚度/质量矩阵 $\to$ 特征值或线性系统 $\to$ 物理量。

**核心对应** 内积 $\leftrightarrow$ 虚功/能量；正交 $\leftrightarrow$ 模态解耦；特征值 $\leftrightarrow$ 固有频率/主应力/临界载荷；秩 $\leftrightarrow$ 独立约束/可控性；SVD $\leftrightarrow$ 降阶与最优逼近；条件数 $\leftrightarrow$ 数值可靠性与稳定性。

**统一原理** 线性代数把连续介质的场方程投影到有限维子空间，用矩阵的谱与分解回答稳定性、响应、可观测与最优问题。极限情形（网格加密、模态截断增加）收敛到连续算子的谱与解。

---
⬆ 上级领域：[[数学]]
