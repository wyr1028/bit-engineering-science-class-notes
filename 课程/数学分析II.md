---
课程名称: 数学分析II
学分: 6
总学时: 96
考核方式: 考试
课程属性: 必修
课程性质: 公共基础课程
开课单位: 数学与统计学院
领域: 数学
tags:
  - 课程
  - 领域/数学
---

# 数学分析II


### 第1章 欧氏空间与多元函数的极限与连续

**1.1 欧氏空间中的点集结构**
定义 $\mathbb{R}^n=\{(x_1,\dots,x_n):x_i\in\mathbb{R}\}$，内积 $\langle x,y\rangle=\sum_{i=1}^n x_iy_i$，范数 $\|x\|=\sqrt{\langle x,x\rangle}$。
约定：以下点集均指 $\mathbb{R}^n$ 的子集。
由内积可导出 Cauchy–Schwarz 不等式
$$\langle x,y\rangle^2\le \|x\|^2\|y\|^2.$$
推导：对任意 $t\in\mathbb{R}$，
$$\because 0\le \|x+ty\|^2=\|x\|^2+2t\langle x,y\rangle+t^2\|y\|^2,$$
上式右端为 $t$ 的二次多项式且恒非负，故判别式
$$\Delta=4\langle x,y\rangle^2-4\|x\|^2\|y\|^2\le 0\ \Rightarrow\ \langle x,y\rangle^2\le\|x\|^2\|y\|^2.$$
再取平方根得 $|\langle x,y\rangle|\le\|x\|\|y\|$，从而三角不等式
$$\because \|x+y\|^2=\|x\|^2+2\langle x,y\rangle+\|y\|^2\le(\|x\|+\|y\|)^2\ \therefore\ \|x+y\|\le\|x\|+\|y\|.$$
几何意义：$\|x-y\|$ 是两点距离，上述不等式说明距离满足三角不等式。

**1.2 邻域、开集、闭集、聚点**
定义：$U(x_0,r)=\{x:\|x-x_0\|<r\}$ 为 $x_0$ 的 $r$-邻域。若 $A$ 中每点都有某邻域含于 $A$，称 $A$ 为开集。若 $\mathbb{R}^n\setminus A$ 为开集，称 $A$ 为闭集。若 $x_0$ 的任一邻域都含 $A\setminus\{x_0\}$ 的点，称 $x_0$ 为 $A$ 的聚点。
推导性质：开集任意并仍为开集；有限交仍为开集。验证：
$$\because x\in\bigcap_{i=1}^m U_i,\ \forall i\ \exists r_i>0:U(x,r_i)\subset U_i,\ \therefore r=\min_i r_i>0,\ U(x,r)\subset\bigcap_i U_i.$$
说明：无限交不一定开，如 $\bigcap_{k=1}^\infty(-1/k,1/k)=\{0\}$ 非开集。
物理意义：连续介质的“内部”“边界”即点集的内部与边界。

**1.3 多元函数的二重极限**
定义：$f:D\subset\mathbb{R}^n\to\mathbb{R}$，$x_0$ 为 $D$ 的聚点，若 $\forall\varepsilon>0,\exists\delta>0$，当 $0<\|x-x_0\|<\delta$ 且 $x\in D$ 时 $|f(x)-L|<\varepsilon$，记 $\lim_{x\to x_0}f(x)=L$。
与累次极限的区别：二重极限要求沿任意路径趋于 $x_0$ 时极限相同；累次极限 $\lim_{x\to a}\lim_{y\to b}f(x,y)$ 是先固定一变量。
定理：若二重极限与两个累次极限都存在，则三者相等。证明思路：
$$\because \forall\varepsilon>0\ \exists\delta>0:\|(x,y)-(a,b)\|<\delta\Rightarrow|f(x,y)-L|<\varepsilon,$$
令 $y\to b$（保持 $|x-a|<\delta$）得 $|\lim_{y\to b}f(x,y)-L|\le\varepsilon$，再令 $x\to a$ 得累次极限等于 $L$。
反例：$f=\dfrac{xy}{x^2+y^2}$，沿 $y=kx$ 有 $f\to k/(1+k^2)$，与 $k$ 有关，故二重极限不存在。
应用：判断函数可微、连续，需先确认二重极限。

**1.4 多元函数的连续性**
定义：$f$ 在 $x_0$ 连续 $\iff\lim_{x\to x_0}f(x)=f(x_0)$。
运算性质：连续函数四则运算（分母非零）仍连续；连续函数的复合仍连续。
推导复合连续性：设 $g$ 在 $x_0$ 连续、$f$ 在 $u_0=g(x_0)$ 连续。$\forall\varepsilon>0,\exists\eta>0:\|u-u_0\|<\eta\Rightarrow\|f(u)-f(u_0)\|<\varepsilon$；对 $\eta$ 取 $\delta>0$ 使 $\|x-x_0\|<\delta\Rightarrow\|g(x)-g(x_0)\|<\eta$，则复合连续。
物理意义：小扰动引起小响应，是“适定”概念的雏形。

**1.5 紧集上连续函数的性质**
定理（有界性与最值性）：若 $K\subset\mathbb{R}^n$ 有界闭（紧），$f:K\to\mathbb{R}$ 连续，则 $f$ 在 $K$ 上有界且取到最大最小值。
证明（最值）：令 $M=\sup_K f$，取 $\{x_k\}\subset K$ 使 $f(x_k)\to M$。由 $K$ 紧，存在子列 $x_{k_j}\to x^*\in K$。由连续性 $f(x^*)=\lim_j f(x_{k_j})=M<\infty$，故最大值可取。
定理（一致连续性）：紧集上连续函数一致连续。证明：若不，则 $\exists\varepsilon_0>0$ 及 $x_k,y_k\in K$，$\|x_k-y_k\|<1/k$ 但 $|f(x_k)-f(y_k)|\ge\varepsilon_0$；取子列 $x_{k_j}\to x^*$，则 $y_{k_j}\to x^*$，由连续性两端趋于 $f(x^*)$，矛盾。
应用：能量泛函在容许集上极小点的存在性论证。

**1.6 连续映射与同胚初步**
定义：$F:D\subset\mathbb{R}^n\to\mathbb{R}^m$ 连续指各分量连续。若 $F$ 为双射且 $F,F^{-1}$ 均连续，称同胚。
几何意义：同胚保持开集、紧性、连通性，是拓扑分类的基本工具。
特例：线性映射 $F(x)=Ax$ 连续，$\|Ax\|\le\|A\|\|x\|$，其中 $\|A\|=\sup_{\|x\|=1}\|Ax\|$。

**1.7 连通性与区域**
定义：若 $D$ 不能分解为两个非空不相交开集之并，称 $D$ 连通。连通开集称为区域。若区域内任两点可用区域内折线连接，称道路连通。
定理：$\mathbb{R}^n$ 中开集连通 $\iff$ 道路连通。
物理意义：区域是积分与场论的基本定义域；保证积分与路径无关讨论有意义。

**1.8 典型应用：极限存在的夹逼判定**
例：证明 $\displaystyle\lim_{(x,y)\to(0,0)}\frac{x^2y^2}{x^2+y^2}=0$。
$$\because 0\le\frac{x^2y^2}{x^2+y^2}\le\frac{x^2y^2}{2|xy|}=\frac{|xy|}{2}\to0,\ \therefore \text{由夹逼定理极限为 }0.$$
量纲：若 $x,y$ 为长度，则 $x^2y^2/(x^2+y^2)$ 量纲为 $\mathrm{L}^2$，与 $|xy|$ 一致。

### 第2章 偏导数、全微分与泰勒展开

**2.1 偏导数**
定义：$f:D\subset\mathbb{R}^n\to\mathbb{R}$，在 $x_0$ 处对 $x_i$ 的偏导数
$$\frac{\partial f}{\partial x_i}(x_0)=\lim_{h\to0}\frac{f(x_0+he_i)-f(x_0)}{h}.$$
几何意义：沿坐标方向切线的斜率。
量纲：$\partial f/\partial x_i$ 的量纲为 $[f]/[x_i]$。

**2.2 高阶偏导数与混合偏导相等定理（Schwarz）**
定义：$\partial_{ij}f=\dfrac{\partial}{\partial x_j}\Big(\dfrac{\partial f}{\partial x_i}\Big)$。
定理：若 $\partial_{ij}f,\partial_{ji}f$ 在 $x_0$ 的某邻域存在且在 $x_0$ 连续，则
$$\partial_{ij}f(x_0)=\partial_{ji}f(x_0).$$
推导：对二元函数，构造
$$\Delta=f(x_0+h,y_0+k)-f(x_0+h,y_0)-f(x_0,y_0+k)+f(x_0,y_0).$$
令 $\varphi(x)=f(x,y_0+k)-f(x,y_0)$，由一元中值定理
$$\Delta=\varphi(x_0+h)-\varphi(x_0)=h\,\varphi'(x_0+\theta h)=h\big[f_x(x_0+\theta h,y_0+k)-f_x(x_0+\theta h,y_0)\big].$$
再对 $y$ 用中值定理：
$$\Delta=hk\,f_{xy}(x_0+\theta h,y_0+\theta' k).$$
同理令 $\psi(y)=f(x_0+h,y)-f(x_0,y)$ 得 $\Delta=hk\,f_{yx}(x_0+\theta_1 h,y_0+\theta_1' k)$。令 $h,k\to0$，由连续性得 $f_{xy}(x_0,y_0)=f_{yx}(x_0,y_0)$。
意义：$C^2$ 函数二阶导矩阵（Hessian）对称，是变分与稳定性的基础。

**2.3 全微分**
定义：$f$ 在 $x_0$ 可微，若存在线性映射 $A$ 使
$$f(x_0+h)=f(x_0)+Ah+o(\|h\|),\quad \|h\|\to0.$$
且 $A$ 的矩阵即 Jacobi 矩阵 $Df(x_0)=(\partial_j f_i)$。
推导：若可微则各偏导存在且 $A$ 唯一。
$$\because f(x_0+he_i)-f(x_0)=A(he_i)+o(|h|)=h\,A e_i+o(|h|),\ \therefore \partial_i f(x_0)=A e_i.$$
几何意义：切平面近似 $f(x_0+h)\approx f(x_0)+\nabla f(x_0)\cdot h$。

**2.4 可微的充分条件**
定理：若 $f$ 的偏导数在 $x_0$ 邻域存在且在 $x_0$ 连续，则 $f$ 在 $x_0$ 可微。
推导：沿路径逐段估计
$$f(x_0+h)-f(x_0)=\Big[f(x_0+h)-\prod\text{逐步改变}\Big]=\sum_{i=1}^n\big[f(\cdots)\big],$$
对第 $i$ 段用一元中值定理
$$\because f(\dots,h_i,\dots)-f(\dots,0,\dots)=h_i\,\partial_i f(\xi_i),\quad \xi_i\text{ 在 }x_0\text{ 与}x_0+h\text{ 之间},$$
$$\therefore f(x_0+h)-f(x_0)-\sum_i h_i\partial_i f(x_0)=\sum_i h_i\big[\partial_i f(\xi_i)-\partial_i f(x_0)\big]=o(\|h\|),$$
由偏导连续性得上式成立，故可微。
反例：$f=\sqrt{x^2+y^2}$ 在原点连续、偏导存在但不可微（$f(h,0)/\|h\|=1\not\to0$）。

**2.5 复合函数链式法则**
设 $z=f(u,v)$ 可微，$u=u(x,y),v=v(x,y)$ 可微，则
$$\frac{\partial z}{\partial x}=\frac{\partial f}{\partial u}\frac{\partial u}{\partial x}+\frac{\partial f}{\partial v}\frac{\partial v}{\partial x}.$$
推导：由可微性
$$\Delta z=f_u\Delta u+f_v\Delta v+o(\sqrt{\Delta u^2+\Delta v^2}),$$
而 $\Delta u=u_x\Delta x+u_y\Delta y+o(\|(\Delta x,\Delta y)\|)$，代入并注意 $\sqrt{\Delta u^2+\Delta v^2}=O(\|(\Delta x,\Delta y)\|)$，比较 $\Delta x$ 系数即得。
特例（极坐标）：$z=f(x,y),\ x=r\cos\theta,\ y=r\sin\theta$，
$$f_r=f_x\cos\theta+f_y\sin\theta,\quad f_\theta=-r f_x\sin\theta+r f_y\cos\theta.$$
应用：拉格朗日力学中广义坐标与直角坐标之间的导数变换。

**2.6 方向导数**
定义：沿单位向量 $l$ 的方向导数
$$\frac{\partial f}{\partial l}(x_0)=\lim_{t\to0^+}\frac{f(x_0+tl)-f(x_0)}{t}.$$
定理：若 $f$ 可微，则 $\dfrac{\partial f}{\partial l}=\nabla f\cdot l$。
推导：$f(x_0+tl)-f(x_0)=\nabla f(x_0)\cdot(tl)+o(t)$，除以 $t$ 令 $t\to0^+$。
几何意义：方向导数为梯度在该方向的投影；当 $l=\nabla f/\|\nabla f\|$ 时取最大增长方向，最大值 $\|\nabla f\|$。

**2.7 梯度与等值面**
定义：$\nabla f=(\partial_1 f,\dots,\partial_n f)$。
性质：$\nabla f$ 垂直于等值面 $f=c$。推导：设曲线 $\gamma(t)$ 在等值面上，则 $f(\gamma(t))=c$，求导得 $\nabla f(\gamma(t))\cdot\gamma'(t)=0$，而 $\gamma'$ 为任意切向，故 $\nabla f\perp$ 切平面。
物理意义：温度梯度指向升温最快方向；势场 $V$ 的力 $\mathbf F=-\nabla V$。
量纲：$[\nabla f]=[f]/[x]$。

**2.8 多元泰勒公式**
定理：$f$ 在 $x_0$ 邻域有 $m+1$ 阶连续偏导，则
$$f(x_0+h)=\sum_{|\alpha|\le m}\frac{D^\alpha f(x_0)}{\alpha!}h^\alpha+R_m,$$
$$R_m=\sum_{|\alpha|=m+1}\frac{D^\alpha f(x_0+\theta h)}{\alpha!}h^\alpha,\quad 0<\theta<1,$$
其中 $\alpha=(\alpha_1,\dots,\alpha_n)$，$\alpha!=\prod\alpha_i!$，$h^\alpha=\prod h_i^{\alpha_i}$。
推导：令 $\phi(t)=f(x_0+th)$，由一元泰勒
$$\phi(t)=\sum_{k=0}^m\frac{\phi^{(k)}(0)}{k!}t^k+\frac{\phi^{(m+1)}(\theta t)}{(m+1)!}t^{m+1},$$
而 $\phi^{(k)}(0)=\dfrac{d^k}{dt^k}f(x_0+th)\big|_{t=0}=\sum_{|\alpha|=k}\dfrac{k!}{\alpha!}D^\alpha f(x_0)h^\alpha$，代入 $t=1$ 即得。
二元二阶特例：
$$f(x_0+h)=f(x_0)+\nabla f\cdot h+\tfrac12 h^\top H_f(x_0)h+o(\|h\|^2).$$
应用：极值判别、误差传播、摄动分析。

**2.9 误差传播与切平面近似**
设 $y=f(x_1,\dots,x_n)$，各 $x_i$ 有误差 $\Delta x_i$，一阶近似
$$\Delta y\approx\sum_i\frac{\partial f}{\partial x_i}\Delta x_i.$$
若误差独立随机，则平方和
$$\sigma_y^2\approx\sum_i\Big(\frac{\partial f}{\partial x_i}\Big)^2\sigma_{x_i}^2.$$
推导：由全微分 $\Delta y=\sum_i f_{x_i}\Delta x_i+o$，取方差并注意协方差项在独立假设下为零。
应用：实验测量、结构不确定性量化。

### 第3章 隐函数定理与逆映射定理

**3.1 隐函数存在唯一性定理**
定理：设 $F(x,y)$ 在 $(x_0,y_0)$ 邻域连续可微，$F(x_0,y_0)=0$，且 $F_y(x_0,y_0)\ne0$，则存在 $x_0$ 邻域 $U$ 及唯一 $y=\varphi(x)$ 使 $\varphi(x_0)=y_0$，$F(x,\varphi(x))=0$，且 $\varphi$ 可微。
推导（存在性用压缩映射）：由 $F_y\ne0$，不妨 $F_y>0$。取闭矩形 $U\times V$，使 $F_y\ge m>0$。对固定 $x$，$F(x,\cdot)$ 严格增，结合 $F(x_0,y_0)=0$ 与连续性，$\exists\,y_1<y_0<y_2$ 使
$$F(x_0,y_1)<0<F(x_0,y_2).$$
由连续性，$x$ 充分接近 $x_0$ 时仍有 $F(x,y_1)<0<F(x,y_2)$；由中值定理存在唯一 $y=\varphi(x)\in(y_1,y_2)$ 使 $F=0$。
可微性与导数：由
$$0=F(x+\Delta x,\varphi(x+\Delta x))-F(x,\varphi(x))=F_x\Delta x+F_y\Delta\varphi+o(\cdot),$$
$$\therefore \varphi'(x)=-\frac{F_x(x,\varphi(x))}{F_y(x,\varphi(x))}.$$
几何意义：曲线 $F=0$ 在 $F_y\ne0$ 处可局部表为函数图像。

**3.2 多方程隐函数组**
设 $F:\mathbb{R}^{n+m}\to\mathbb{R}^m$，$F(x_0,y_0)=0$，$\det F_y(x_0,y_0)\ne0$，则存在 $y=\varphi(x)$ 使 $F(x,\varphi(x))=0$，且
$$\varphi'(x)=-\big[F_y\big]^{-1}F_x.$$
推导：对 $F(x,\varphi(x))=0$ 求导得 $F_x+F_y\varphi'=0$，由 $\det F_y\ne0$ 可解 $\varphi'=-F_y^{-1}F_x$。
应用：约束方程组降维，力学中完整约束的广义坐标选取。

**3.3 逆映射定理**
定理：设 $f:U\subset\mathbb{R}^n\to\mathbb{R}^n$ 连续可微，$x_0\in U$，$\det Df(x_0)\ne0$，则存在 $x_0$ 邻域 $U_0$ 与 $f(x_0)$ 邻域 $V_0$，使 $f:U_0\to V_0$ 为 $C^1$ 双射，且
$$Df^{-1}(y_0)=\big[Df(x_0)\big]^{-1}.$$
推导：令 $g(x)=x-Df(x_0)^{-1}(f(x)-y_0)$。因 $Dg(x_0)=I-Df(x_0)^{-1}Df(x_0)=0$，且 $g$ 连续可微，在 $x_0$ 小邻域内 $\|Dg\|\le1/2$，故 $g$ 为压缩映射，存在唯一不动点 $x$ 满足 $g(x)=x$，即 $f(x)=y_0$。对任意 $y$ 靠近 $y_0$ 同理。
几何意义：局部坐标变换（微分同胚），雅可比行列式给出体积伸缩率。

**3.4 雅可比行列式**
定义：$\dfrac{\partial(u_1,\dots,u_n)}{\partial(x_1,\dots,x_n)}=\det\Big(\dfrac{\partial u_i}{\partial x_j}\Big)$。
链式相乘：
$$\frac{\partial(u)}{\partial(x)}=\frac{\partial(u)}{\partial(v)}\cdot\frac{\partial(v)}{\partial(x)}.$$
推导：由 $Du=D(u)\circ Dv$ 取行列式得 $\det Du=\det(Du)\det(Dv)$。
几何意义：$n$ 维体积元之比 $dV_u=|\det J|\,dV_x$。
应用：坐标变换积分、连续体运动学变形梯度 $F=\partial x/\partial X$，$J=\det F$。

**3.5 变量代换公式**
设变换 $x=\Phi(u)$ 为 $C^1$ 微分同胚，则
$$\int_{\Phi(D)}f(x)\,dx=\int_D f(\Phi(u))\,|\det D\Phi(u)|\,du.$$
推导（无穷小体积）：局部线性化
$$\Delta x=\Phi(u+\Delta u)-\Phi(u)=D\Phi(u)\Delta u+o(\|\Delta u\|),$$
平行体体积 $|\det D\Phi|\,|\Delta u|$，积分求和取极限。
特例：极坐标 $\det J=r$，$dxdy=r\,drd\theta$。

**3.6 反函数高阶导与内蕴几何**
设 $y=f(x)$，则 $f^{-1}$ 的二阶导
$$(f^{-1})''(y)=-\frac{f''(x)}{[f'(x)]^3},\quad y=f(x).$$
推导：由 $(f^{-1})'(f(x))=1/f'(x)$，两边对 $x$ 求导
$$(f^{-1})''(f(x))f'(x)=-\frac{f''(x)}{[f'(x)]^2}\ \Rightarrow\ (f^{-1})''(f(x))=-\frac{f''(x)}{[f'(x)]^3}.$$
应用：参数曲线曲率 $\kappa=\dfrac{|x'y''-y'x''|}{(x'^2+y'^2)^{3/2}}$。

### 第4章 多元函数极值与拉格朗日乘数法

**4.1 无条件极值的必要条件**
定义：$f$ 在 $x_0$ 取局部极小，若 $\exists r>0$ 使 $\|x-x_0\|<r\Rightarrow f(x)\ge f(x_0)$。
定理：若 $f$ 在 $x_0$ 可微且取局部极值，则 $\nabla f(x_0)=0$。
推导：取 $g_i(t)=f(x_0+te_i)$，在 $t=0$ 取极值，故 $g_i'(0)=0$，即 $\partial_i f(x_0)=0$。
几何意义：极值点处切平面水平；驻点可能是鞍点。

**4.2 二阶充分条件**
定理：若 $\nabla f(x_0)=0$ 且 Hessian $H_f(x_0)$ 正定，则 $x_0$ 为严格局部极小；负定则为严格局部极大；不定则为鞍点。
推导：由泰勒
$$f(x_0+h)-f(x_0)=\tfrac12 h^\top H h+o(\|h\|^2).$$
正定时最小特征值 $\lambda_{\min}>0$，故 $h^\top Hh\ge\lambda_{\min}\|h\|^2$，当 $\|h\|$ 充分小余项不超过 $\tfrac14\lambda_{\min}\|h\|^2$，故差为正。
二元判别：$A=f_{xx},B=f_{xy},C=f_{yy}$，$\Delta=AC-B^2$；$\Delta>0,A>0$ 极小；$\Delta>0,A<0$ 极大；$\Delta<0$ 鞍点；$\Delta=0$ 需高阶判断。
应用：弹性势能极小对应稳定平衡，$H$ 即刚度矩阵。

**4.3 条件极值与拉格朗日乘数法**
问题：$\min f(x)$ s.t. $g(x)=0$，$g:\mathbb{R}^n\to\mathbb{R}^m$。
定理：若 $x_0$ 为正则约束（$\nabla g_i(x_0)$ 线性无关）下的极值，则存在 $\lambda\in\mathbb{R}^m$ 使
$$\nabla f(x_0)=\sum_{j=1}^m\lambda_j\nabla g_j(x_0).$$
推导：约束曲面 $M=\{g=0\}$ 在 $x_0$ 的切空间 $T=\{\xi:Dg(x_0)\xi=0\}$。若 $x_0$ 为约束极值，则对任意 admissible 变分曲线 $\gamma(t)\subset M$，$\gamma(0)=x_0$，
$$\because \frac{d}{dt}f(\gamma(t))\big|_{t=0}=\nabla f(x_0)\cdot\gamma'(0)=0\ \forall\gamma'(0)\in T,$$
$$\therefore \nabla f(x_0)\perp T=\mathrm{span}\{\nabla g_j(x_0)\}^\perp,$$
由线性代数，$\nabla f(x_0)\in\mathrm{span}\{\nabla g_j(x_0)\}$，即存在 $\lambda_j$。
几何意义：极值点处目标函数与约束曲面法向共线。
应用：最小势能约束平衡、最大熵、最优控制必要条件。

**4.4 拉格朗日函数与求解**
构造 $L(x,\lambda)=f(x)-\sum_j\lambda_j g_j(x)$，令 $\nabla_xL=0,\ \nabla_\lambda L=0$。
例：求原点到曲面 $g(x,y,z)=0$ 最短距离，即 $\min(x^2+y^2+z^2)$ s.t. $g=0$，条件为 $2x=\lambda g_x$ 等。
特例：三维中两约束、单约束的叉乘形式 $\nabla f=\lambda\nabla g$。

**4.5 不等式约束的KKT条件**
问题 $\min f$ s.t. $g_j(x)\le0,\ h_k(x)=0$。若约束规范成立，则存在 $\mu_j\ge0,\lambda_k$ 使
$$\nabla f+\sum_j\mu_j\nabla g_j+\sum_k\lambda_k\nabla h_k=0,\quad \mu_jg_j=0.$$
推导（互补松弛）：若 $g_j(x_0)<0$（非活跃），则局部可视为无约束，故 $\mu_j=0$；活跃约束按等式处理。
应用：接触力学、塑性屈服（Kuhn–Tucker 条件对应屈服面正交流动法则）。

**4.6 全局极值与紧性**
方法：闭有界区域 $D$ 上 $f$ 连续必取到全局最值；候选点 = 内部驻点 + 边界约束极值点 + 不可微点。
推导依据：第1.5节紧集最值定理。
应用：结构优化、参数反演。

**4.7 最小势能原理的数学表述**
设弹性系统总势能 $\Pi(u)=\int_\Omega W(\nabla u)\,dV-\int_\Omega b\cdot u\,dV-\int_{\Gamma_t}t\cdot u\,dS$。
平衡满足 $\delta\Pi=0$。推导：对任意 admissible 变分 $\delta u$（满足位移边界条件），
$$\delta\Pi=\int_\Omega \sigma:\nabla\delta u\,dV-\int_\Omega b\cdot\delta u\,dV-\int_{\Gamma_t}t\cdot\delta u\,dS.$$
由散度定理与平衡方程 $\nabla\cdot\sigma+b=0$，$\sigma n=t$ 知 $\delta\Pi=0$ 对所有 $\delta u$ 成立，反之亦然。
极限/特例：线性弹性时 $\Pi$ 二次，二阶变分正定对应唯一稳定平衡。

**4.8 典型应用：结构压杆临界载荷**
欧拉压杆势能
$$\Pi=\frac12\int_0^L EI(w'')^2dx-\frac P2\int_0^L (w')^2dx.$$
变分 $\delta\Pi=0$ 得 $EIw''''+Pw''=0$，边界 $w(0)=w(L)=0$，$w''(0)=w''(L)=0$。令 $w=\sin(n\pi x/L)$ 得
$$P_{cr}=n^2\pi^2\frac{EI}{L^2},$$
最小临界载荷 $n=1$。量纲验证：$EI/L^2$ 为力，正确。

### 第5章 重积分及其坐标变换

**5.1 二重积分的定义**
定义：$f$ 在可求面积区域 $D$ 上有界，分割 $D=\cup\Delta D_i$，$\|\Delta\|=\max\mathrm{diam}\Delta D_i$，取 $\xi_i\in\Delta D_i$，
$$\iint_D f\,dA=\lim_{\|\Delta\|\to0}\sum_i f(\xi_i)\Delta A_i.$$
可积的必要条件：有界。充分条件：连续或分片连续有界。
推导（Darboux）：用上和 $U(f,\Delta)$、下和 $L(f,\Delta)$，可积 $\iff\lim(U-L)=0$。

**5.2 Fubini 定理与累次积分**
定理：若 $f$ 在矩形 $D=[a,b]\times[c,d]$ 可积且 $\int_c^d f(x,y)dy$ 对每个 $x$ 存在，则
$$\iint_D f=\int_a^b\Big(\int_c^d f(x,y)\,dy\Big)dx.$$
推导思路：作分割，先在 $y$ 方向合并求和；由可积性上下和夹逼。
应用：将高维积分降为一维。

**5.3 极坐标变换**
$$\iint_D f(x,y)\,dxdy=\iint_{D'} f(r\cos\theta,r\sin\theta)\,r\,drd\theta.$$
推导：变换 $(x,y)=(r\cos\theta,r\sin\theta)$，Jacobi 矩阵
$$J=\begin{pmatrix}\cos\theta&-r\sin\theta\\ \sin\theta&r\cos\theta\end{pmatrix},\quad \det J=r\cos^2\theta+r\sin^2\theta=r.$$
故 $dxdy=|\det J|\,drd\theta=r\,drd\theta$。
几何意义：面积微元为扇形环，$r\,dr\,d\theta$。

**5.4 三重积分的柱坐标与球坐标**
柱坐标 $(r,\theta,z)$：$x=r\cos\theta,y=r\sin\theta,z=z$，
$$\det J=r\ \Rightarrow\ dV=r\,dr\,d\theta\,dz.$$
球坐标 $(r,\theta,\varphi)$：$x=r\sin\varphi\cos\theta,\ y=r\sin\varphi\sin\theta,\ z=r\cos\varphi$，其中 $\varphi\in[0,\pi]$ 为极角。
$$J=\begin{pmatrix}\sin\varphi\cos\theta&r\cos\varphi\cos\theta&-r\sin\varphi\sin\theta\\ \sin\varphi\sin\theta&r\cos\varphi\sin\theta&r\sin\varphi\cos\theta\\ \cos\varphi&-r\sin\varphi&0\end{pmatrix}.$$
计算行列式：
$$\det J=r^2\sin\varphi\ \Rightarrow\ dV=r^2\sin\varphi\,dr\,d\varphi\,d\theta.$$
推导（几何法）：纬度 $\varphi$ 处环带面积 $\sim r\sin\varphi\,d\theta\cdot r\,d\varphi$，径向厚度 $dr$，故 $dV=r^2\sin\varphi\,dr\,d\varphi\,d\theta$。
应用：球对称势、流体球、转动惯量。

**5.5 重积分的物理应用**
质量：$m=\iiint_\Omega\rho\,dV$。
质心：$\bar x=\dfrac1m\iiint x\rho\,dV$，其余同理。
转动惯量：$I_z=\iiint(x^2+y^2)\rho\,dV$。
平行轴定理推导：设 $\bar x=0$，则
$$I_z=\iiint[(x')^2+(y')^2]\,dm=I_{z_c}+md^2,$$
其中 $d$ 为轴间距，交叉项因 $\iiint x\,dm=0$ 消失。
量纲：$[m]=\mathrm{M}$，$[I]=\mathrm{M\,L^2}$。

**5.6 反常重积分与高斯积分**
定义：无界区域上 $\iint_{\mathbb{R}^2}e^{-x^2-y^2}dA=\lim_{R\to\infty}\iint_{B_R}$。
计算：
$$\because \iint_{\mathbb{R}^2}e^{-x^2-y^2}dxdy=\int_0^{2\pi}\!\!\int_0^\infty e^{-r^2}r\,dr\,d\theta=\int_0^{2\pi}\frac12 d\theta=\pi,$$
$$\therefore \int_{-\infty}^\infty e^{-x^2}dx=\sqrt\pi.$$
极坐标下逐元分离：
$$\Big(\int_{-\infty}^\infty e^{-x^2}dx\Big)^2=\iint_{\mathbb{R}^2}e^{-(x^2+y^2)}dxdy=\pi.$$
应用：概率正态化、热核、扩散方程基本解。

**5.7 重积分的变量代换一般公式**
$$\iint_D f(x,y)dxdy=\iint_{D'}f(x(u,v),y(u,v))|J(u,v)|\,du\,dv,$$
$$J=\frac{\partial(x,y)}{\partial(u,v)}\ne0.$$
推导：局部线性映射体积伸缩因子 $|J|$，细分取极限。
特例：椭圆 $\frac{x^2}{a^2}+\frac{y^2}{b^2}\le1$ 经 $x=ar\cos\theta,y=br\sin\theta$，$|J|=abr$，面积 $\pi ab$。

**5.8 典型应用：均匀球的转动惯量**
取球 $x^2+y^2+z^2\le R^2$，密度 $\rho$。
$$I_z=\iiint(x^2+y^2)\rho\,dV.$$
球坐标：
$$I_z=\rho\int_0^{2\pi}\!\!\int_0^\pi\!\!\int_0^R r^2\sin^2\varphi\cdot r^2\sin\varphi\,dr\,d\varphi\,d\theta=\rho\cdot2\pi\cdot\frac{R^5}{5}\int_0^\pi\sin^3\varphi\,d\varphi.$$
$$\int_0^\pi\sin^3\varphi\,d\varphi=\int_{-1}^{1}(1-u^2)du=\frac43.$$
$$\therefore I_z=\rho\cdot2\pi\cdot\frac{R^5}{5}\cdot\frac43=\frac{8\pi}{15}\rho R^5=\frac25 MR^2,\quad M=\frac43\pi R^3\rho.$$
量纲：$\mathrm{M\,L^2}$。

### 第6章 曲线积分与Green公式

**6.1 第一型曲线积分**
定义：对弧长微元 $ds$，
$$\int_L f\,ds=\lim\sum_i f(\xi_i)\Delta s_i.$$
参数表示 $\gamma(t)=(x(t),y(t)),\ t\in[a,b]$：
$$\int_L f\,ds=\int_a^b f(x(t),y(t))\sqrt{x'^2+y'^2}\,dt.$$
推导：弧长微元 $ds=\sqrt{dx^2+dy^2}=\sqrt{x'^2+y'^2}\,dt$。
应用：曲线质量、质心。

**6.2 第二型曲线积分（功）**
$$\int_L P\,dx+Q\,dy=\lim\sum_i\big[P(\xi_i)\Delta x_i+Q(\xi_i)\Delta y_i\big]=\int_a^b\big[P x'+Q y'\big]dt.$$
推导：$dx=x'dt,dy=y'dt$ 代入。
物理意义：力 $\mathbf F=(P,Q)$ 沿曲线做功 $W=\int_L \mathbf F\cdot d\mathbf r$。
量纲：$[W]=\mathrm{M\,L^2T^{-2}}$。

**6.3 Green 公式**
定理：$D$ 为有界闭区域，边界 $\partial D$ 分段光滑正向，$P,Q\in C^1(D)$，则
$$\oint_{\partial D}P\,dx+Q\,dy=\iint_D\Big(\frac{\partial Q}{\partial x}-\frac{\partial P}{\partial y}\Big)dxdy.$$
推导（矩形区域）：设 $D=[a,b]\times[c,d]$。
$$\because \iint_D\frac{\partial Q}{\partial x}dxdy=\int_c^d\!\!\int_a^b Q_x\,dx\,dy=\int_c^d[Q(b,y)-Q(a,y)]dy,$$
而沿正向边界
$$\oint_{\partial D}Q\,dy=\int_c^dQ(b,y)dy+\int_d^cQ(a,y)dy=\int_c^d[Q(b,y)-Q(a,y)]dy,$$
$$\therefore \iint_D Q_x\,dxdy=\oint_{\partial D}Q\,dy.$$
同理 $\iint_D(-P_y)\,dxdy=\oint_{\partial D}P\,dx$。相加得 Green 公式。
一般区域经分割与逼近。
物理意义：二维环量—旋度关系；几何意义：面积 $A=\frac12\oint(x\,dy-y\,dx)$。

**6.4 曲线积分与路径无关**
定理：在单连通区域 $D$ 上，以下等价：
$$\oint_C P\,dx+Q\,dy=0\ \forall\text{闭曲线}\iff \int_L\text{ 与路径无关}\iff \frac{\partial P}{\partial y}=\frac{\partial Q}{\partial x}\ \text{on }D.$$
推导：由 Green，$\oint_C=\iint(Q_x-P_y)$，故若 $Q_x=P_y$ 则环量为零。反之若某点 $Q_x\ne P_y$，取小邻域内闭曲线环量非零，矛盾。
几何意义：向量场无旋，存在势函数。

**6.5 势函数与全微分**
若 $P\,dx+Q\,dy=du$，即 $P=u_x,Q=u_y$，则
$$\int_{(x_0,y_0)}^{(x,y)}P\,dx+Q\,dy=u(x,y)-u(x_0,y_0).$$
推导：链式法则 $\dfrac{du}{dt}=u_x x'+u_y y'=Px'+Qy'$，积分即得。
求法：$u(x,y)=\int_{x_0}^xP(s,y_0)ds+\int_{y_0}^yQ(x,t)dt$。
应用：保守力做功与路径无关，势能 $V=-u$。

**6.6 典型应用：环量与面积**
由 Green 令 $P=-y,Q=x$：
$$A=\frac12\oint_{\partial D}(x\,dy-y\,dx).$$
对椭圆 $x=a\cos t,y=b\sin t$：
$$A=\frac12\int_0^{2\pi}(ab\cos^2t+ab\sin^2t)dt=\pi ab.$$
与解析几何一致。

### 第7章 曲面积分与场论三大公式

**7.1 曲面参数化与第一型曲面积分**
设曲面 $S:\mathbf r(u,v)=(x,y,z)(u,v)\in D$。
面积元：
$$dS=\|\mathbf r_u\times\mathbf r_v\|\,du\,dv.$$
推导：参数网格的切向量 $\mathbf r_u du,\mathbf r_v dv$ 张成平行四边形，面积为其叉乘模。
$$\int_S f\,dS=\iint_D f(\mathbf r(u,v))\|\mathbf r_u\times\mathbf r_v\|\,du\,dv.$$
量纲：$[dS]=\mathrm{L^2}$。

**7.2 第二型曲面积分（通量）**
$$\iint_S P\,dydz+Q\,dzdx+R\,dxdy=\iint_S\mathbf F\cdot\mathbf n\,dS.$$
方向由 $\mathbf n=\dfrac{\mathbf r_u\times\mathbf r_v}{\|\mathbf r_u\times\mathbf r_v\|}$ 给出。
微元关系：$dydz=n_x\,dS,\ dzdx=n_y\,dS,\ dxdy=n_z\,dS$。
推导：面积向量 $d\mathbf S=\mathbf r_u\times\mathbf r_v\,du\,dv$，其分量为 $d\mathbf S=(dydz,dzdx,dxdy)$。
物理意义：单位时间通过曲面的体积流量 $\Phi=\iint_S\mathbf v\cdot\mathbf n\,dS$。

**7.3 高斯（散度）公式**
定理：$\Omega$ 有界闭区域，边界 $\partial\Omega$ 分片光滑取外法向，$\mathbf F\in C^1(\Omega)$，则
$$\oiint_{\partial\Omega}\mathbf F\cdot d\mathbf S=\iiint_\Omega\nabla\cdot\mathbf F\,dV.$$
推导（长方体）：设 $\Omega=[a,b]\times[c,d]\times[e,f]$，$\mathbf F=(P,Q,R)$。
$$\because \iiint_\Omega R_z\,dV=\int_a^b\!\!\int_c^d[R(x,y,f)-R(x,y,e)]dxdy,$$
$$\oiint_{\partial\Omega}R\,dxdy=\int_a^b\!\!\int_c^d R(x,y,f)dxdy-\int_a^b\!\!\int_c^dR(x,y,e)dxdy,$$
两式相等；对 $P,Q$ 同理，相加即得。一般区域通过分割逼近。
物理意义：通量等于源密度体积分。极限特例：缩为一点得
$$\nabla\cdot\mathbf F=\lim_{V\to0}\frac{1}{V}\oiint_{\partial V}\mathbf F\cdot d\mathbf S.$$

**7.4 斯托克斯公式**
定理：$S$ 为分片光滑定向曲面，边界 $\partial S$ 取诱导正向，$\mathbf F\in C^1$，则
$$\oint_{\partial S}\mathbf F\cdot d\mathbf r=\iint_S(\nabla\times\mathbf F)\cdot d\mathbf S.$$
推导（平直情形）：设 $S$ 为 $z$ 平面区域 $D$，取 $\mathbf F=(P(x,y),Q(x,y),0)$，则 $\nabla\times\mathbf F=(0,0,Q_x-P_y)$，$d\mathbf S=(0,0,dxdy)$，
$$\iint_S(\nabla\times\mathbf F)\cdot d\mathbf S=\iint_D(Q_x-P_y)dxdy,$$
由 Green 公式等于 $\oint_{\partial D}Pdx+Qdy=\oint_{\partial S}\mathbf F\cdot d\mathbf r$。
一般曲面通过参数化与拉回。
物理意义：环量等于旋度通量。极限特例得旋度定义
$$\mathbf n\cdot(\nabla\times\mathbf F)=\lim_{A\to0}\frac{1}{A}\oint_{\partial A}\mathbf F\cdot d\mathbf r.$$

**7.5 散度与旋度的场论恒等式**
定义：$\nabla\cdot\mathbf F=\partial_xF_x+\partial_yF_y+\partial_zF_z$，$\nabla\times\mathbf F=(\partial_yF_z-\partial_zF_y,\dots)$。
恒等式：
$$\nabla\cdot(\nabla\times\mathbf F)=0,\qquad \nabla\times(\nabla f)=\mathbf 0.$$
推导：
$$\because \nabla\cdot(\nabla\times\mathbf F)=\partial_x(\partial_yF_z-\partial_zF_y)+\cdots=\partial_{xy}F_z-\partial_{xz}F_y+\cdots=0,$$
混合偏导相消（$C^2$）。
应用：涡量场无散 $\nabla\cdot\boldsymbol\omega=0$；无旋场为梯度场。

**7.6 连续性方程的推导**
设控制体 $\Omega(t)$，质量守恒
$$\frac{d}{dt}\iiint_{\Omega(t)}\rho\,dV=0.$$
由雷诺输运定理（详见第18章）
$$\frac{d}{dt}\iiint_{\Omega(t)}\rho\,dV=\iiint_{\Omega(t)}\Big(\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf v)\Big)dV=0.$$
由于 $\Omega(t)$ 任意，
$$\boxed{\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf v)=0}.$$
不可压时 $\rho$ 常数，得 $\nabla\cdot\mathbf v=0$。
量纲：$[\partial_t\rho]=\mathrm{M\,L^{-3}T^{-1}}$，$[\nabla\cdot(\rho\mathbf v)]$ 同。

**7.7 典型应用：Maxwell 方程组积分形式**
由 Gauss 与 Stokes 可统一写出
$$\oiint_S\mathbf D\cdot d\mathbf S=\iiint_V\rho_e\,dV,\quad \oint_C\mathbf E\cdot d\mathbf r=-\frac{d}{dt}\iint_S\mathbf B\cdot d\mathbf S.$$
推导：对微分形式两边积分并用高斯/斯托克斯公式。
意义：场论公式是电磁、流体、弹性统一积分定理的语言。

### 第8章 含参变量积分

**8.1 定义与连续性**
定义：$F(y)=\int_a^b f(x,y)\,dx,\ y\in[c,d]$。
定理：若 $f$ 在 $[a,b]\times[c,d]$ 连续，则 $F$ 连续。
推导：$f$ 在紧集一致连续，$\forall\varepsilon,\exists\delta:\ |y-y_0|<\delta\Rightarrow|f(x,y)-f(x,y_0)|<\varepsilon/(b-a)$，
$$|F(y)-F(y_0)|\le\int_a^b|f(x,y)-f(x,y_0)|dx<\varepsilon.$$

**8.2 可微性与莱布尼茨公式**
定理：若 $f,f_y$ 连续，则
$$F'(y)=\int_a^b f_y(x,y)\,dx.$$
推导：由中值定理
$$\frac{F(y+h)-F(y)}{h}=\int_a^b\frac{f(x,y+h)-f(x,y)}{h}dx=\int_a^b f_y(x,y+\theta(x)h)dx,$$
由 $f_y$ 一致连续，令 $h\to0$ 得 $F'(y)=\int f_y\,dx$。
变限推广：
$$\frac{d}{dt}\int_{a(t)}^{b(t)}f(x,t)\,dx=f(b(t),t)b'(t)-f(a(t),t)a'(t)+\int_{a(t)}^{b(t)}f_t(x,t)\,dx.$$
推导：令 $G(u,v,t)=\int_u^vf(x,t)dx$，用链式法则 $\frac{d}{dt}G(a,b,t)=G_u a'+G_v b'+G_t$，其中 $G_u=-f(a,t),G_v=f(b,t)$。

**8.3 积分次序交换**
定理：$f$ 在矩形上连续，则
$$\int_a^b\!\!\int_c^d f\,dydx=\int_c^d\!\!\int_a^b f\,dxdy.$$
推导：令 $F(y)=\int_a^b f\,dx$，由 Fubini 亦可直接由二重积分两种次序相等推出。
应用：简化难算积分，如
$$\int_0^1\frac{x^b-x^a}{\ln x}dx=\int_0^1\!\!\int_a^b x^y\,dy\,dx=\int_a^b\!\!\int_0^1x^y\,dx\,dy=\ln\frac{b+1}{a+1}.$$

**8.4 Euler 积分**
Gamma 函数 $\Gamma(s)=\int_0^\infty x^{s-1}e^{-x}dx,\ s>0$。
递推：分部积分
$$\Gamma(s+1)=\int_0^\infty x^se^{-x}dx=\big[-x^se^{-x}\big]_0^\infty+s\int_0^\infty x^{s-1}e^{-x}dx=s\Gamma(s).$$
$s=n\in\mathbb{N}$ 时 $\Gamma(n+1)=n!$。
Beta 函数 $B(p,q)=\int_0^1t^{p-1}(1-t)^{q-1}dt$，满足
$$B(p,q)=\frac{\Gamma(p)\Gamma(q)}{\Gamma(p+q)}.$$
推导：$\Gamma(p)\Gamma(q)=\int_0^\infty\!\!\int_0^\infty u^{p-1}v^{q-1}e^{-(u+v)}dudv$，令 $u=x^2,v=y^2$ 再转极坐标，或令 $u=st,v=s(1-t)$ 作变量代换，得 $B$。
量纲：$\Gamma$ 无量纲（$x$ 已无量纲化）。

**8.5 含参积分的物理应用**
热传导基本解：
$$u(x,t)=\frac{1}{\sqrt{4\pi kt}}\int_{-\infty}^\infty e^{-(x-\xi)^2/(4kt)}u_0(\xi)\,d\xi.$$
验证：由第5.6节高斯积分，$\int u\,dx=\int u_0\,d\xi$ 质量守恒。当 $t\to0^+$，核趋于 $\delta(x-\xi)$。

**8.6 一致有界与积分号下取极限**
若 $f(x,y)\to g(x)$ 当 $y\to y_0$ 对 $x$ 一致，则
$$\lim_{y\to y_0}\int_a^b f(x,y)dx=\int_a^b g(x)dx.$$
推导：$\big|\int(f-g)\big|\le(b-a)\sup_x|f-g|\to0$。
应用：由含参极限求定积分。

### 第9章 反常积分与含参反常积分

**9.1 无穷限反常积分**
定义：$\int_a^\infty f\,dx=\lim_{R\to\infty}\int_a^R f\,dx$。
比较判别：若 $0\le f\le g$，$\int g$ 收敛则 $\int f$ 收敛。
Cauchy 准则：$\forall\varepsilon,\exists M:\ A_2>A_1>M\Rightarrow\big|\int_{A_1}^{A_2}f\big|<\varepsilon$。
$p$-积分：$\int_1^\infty x^{-p}dx$ 当 $p>1$ 收敛，$p\le1$ 发散。推导：
$$p\ne1:\int_1^R x^{-p}dx=\frac{R^{1-p}-1}{1-p}\to\begin{cases}\frac{1}{p-1},&p>1\\\infty,&p<1\end{cases}.$$

**9.2 无界函数反常积分**
定义：$\int_a^b f\,dx=\lim_{\eta\to0^+}\int_a^{b-\eta}f\,dx$（$b$ 为奇点）。
$a$ 附近 $f\sim(x-a)^{-p}$，收敛 $\iff p<1$。推导同上。

**9.3 含参反常积分的一致收敛**
Weierstrass 判别：若 $|f(x,y)|\le M(x)$（与 $y$ 无关）且 $\int_a^\infty M(x)dx$ 收敛，则 $\int_a^\infty f(x,y)dx$ 在 $y\in Y$ 上一致收敛。
推导：由 $\big|\int_{A_1}^{A_2}f\,dx\big|\le\int_{A_1}^{A_2}M\,dx$，对 $\varepsilon$ 取 $M$ 使尾积分 $<\varepsilon$，与 $y$ 无关。
Dirichlet 判别：$\int_a^A f$ 对 $A,y$ 一致有界且 $g(x,y)$ 关于 $x$ 单调、一致趋于 $0$，则 $\int fg$ 一致收敛。

**9.4 一致收敛积分的性质**
若 $f,f_y$ 连续且 $\int f$ 一致收敛、$\int f_y$ 一致收敛，则
$$F'(y)=\int_a^\infty f_y(x,y)dx.$$
推导：用有限截断 $F_R$ 的莱布尼茨公式，再令 $R\to\infty$，由一致收敛交换极限。

**9.5 积分号下求导计算**
例：$F(a)=\int_0^\infty e^{-ax}\frac{\sin x}{x}dx$。求导
$$F'(a)=-\int_0^\infty e^{-ax}\sin x\,dx=-\frac{1}{a^2+1}.$$
再由 $F(\infty)=0$ 得
$$F(a)=\frac{\pi}{2}-\arctan a.$$
令 $a\to0^+$ 得 Dirichlet 积分
$$\int_0^\infty\frac{\sin x}{x}dx=\frac{\pi}{2}.$$
推导中用到的 $\int_0^\infty e^{-ax}\sin x\,dx$ 由两次分部积分或 $\mathrm{Im}\int e^{(-a+i)x}dx=1/(a^2+1)$。

**9.6 典型的振荡积分**
Fresnel 积分：$\int_0^\infty\cos(x^2)dx=\int_0^\infty\sin(x^2)dx=\frac12\sqrt{\frac\pi2}$。
推导：令 $t=x^2$，则
$$\int_0^\infty e^{ix^2}dx=\frac12\int_0^\infty t^{-1/2}e^{it}dt,$$
绕复平面扇形围道，利用高斯积分 $\int_0^\infty t^{-1/2}e^{-t}dt=\Gamma(1/2)=\sqrt\pi$ 与留数，得 $\sqrt\pi e^{i\pi/4}$，取实虚部。
应用：衍射光学、波的色散。

### 第10章 数项级数

**10.1 收敛定义与Cauchy准则**
定义：$S_n=\sum_{k=1}^n a_k$，若 $S_n$ 收敛称级数 $\sum a_k$ 收敛。
Cauchy 准则：$\sum a_k$ 收敛 $\iff\forall\varepsilon>0,\exists N:\ n>m\ge N\Rightarrow\big|\sum_{k=m+1}^n a_k\big|<\varepsilon$。
推论：收敛的必要条件 $\lim a_k=0$。推导：取 $n=m+1$，则 $|a_{m+1}|<\varepsilon$。

**10.2 正项级数判别法**
比较判别：$0\le a_n\le b_n$，$\sum b_n$ 收敛 $\Rightarrow\sum a_n$ 收敛。
比值判别（d'Alembert）：若 $\lim\dfrac{a_{n+1}}{a_n}=q$，$q<1$ 收敛，$q>1$ 发散。
推导：$q<1$ 时取 $r\in(q,1)$，$n$ 大时 $a_{n+1}\le r a_n$，故 $a_n\le a_Nr^{n-N}$，与几何级数比较。
根值判别（Cauchy）：$\lim\sqrt[n]{a_n}=q$，结论同上。推导：$q<1$ 时 $a_n\le r^n$。
积分判别：$f$ 非负递减，$\sum f(n)$ 与 $\int_1^\infty f(x)dx$ 同敛散。推导：
$$\int_1^{N+1}f(x)dx\le\sum_{n=1}^N f(n)\le f(1)+\int_1^Nf(x)dx.$$

**10.3 交错级数与Leibniz判别**
定理：$a_n\downarrow0$，则 $\sum(-1)^{n-1}a_n$ 收敛，且余和 $|R_N|\le a_{N+1}$。
推导：偶子列 $S_{2k}$ 递增有上界 $a_1$，奇子列 $S_{2k+1}$ 递减有下界，二者同趋于 $S$；余和由交错结构夹逼。
例：$\sum(-1)^{n-1}/n=\ln2$。

**10.4 绝对收敛与条件收敛**
定义：$\sum|a_n|$ 收敛称绝对收敛，$\sum a_n$ 收敛但 $\sum|a_n|$ 发散称条件收敛。
定理：绝对收敛 $\Rightarrow$ 收敛。推导：由 $\big|\sum_{k=m}^na_k\big|\le\sum|a_k|$ 及 Cauchy 准则。
Riemann 重排定理：条件收敛级数可经重排收敛到任意实数或发散。直观：正负项各发散且趋于零，可按目标部分和插入项逼近。
应用：幂级数边界收敛性、Fourier 级数条件收敛。

**10.5 级数运算**
线性：收敛级数可逐项线性组合。
乘积（Cauchy）：$\sum a_n,\sum b_n$ 绝对收敛时其 Cauchy 乘积 $\sum c_n,\ c_n=\sum_{k=0}^na_kb_{n-k}$ 绝对收敛且
$$\sum c_n=\Big(\sum a_n\Big)\Big(\sum b_n\Big).$$
推导：由绝对收敛可重排二重和，分组对角求和即得。
应用：指数函数加法公式 $e^{x+y}=e^xe^y$ 的级数证明。

**10.6 典型应用：调和级数与Euler常数**
$\sum1/n$ 发散（积分判别，$\int_1^Ndx/x=\ln N\to\infty$）。
部分和 $\gamma=\lim_{N\to\infty}(\sum_{n=1}^N1/n-\ln N)$ 存在且 $0<\gamma<1$。推导：$a_n=\sum_{k=1}^n1/k-\ln n$ 单调递减有下界，由单调有界定理收敛。

### 第11章 函数项级数与一致收敛

**11.1 逐点收敛与一致收敛**
定义：$\sum u_n(x)$ 在 $E$ 上逐点收敛于 $S(x)$，若对每个 $x$，部分和 $S_n(x)\to S(x)$。
一致收敛：$\sup_{x\in E}|S_n(x)-S(x)|\to0$。
区别：一致收敛误差界与 $x$ 无关。反例 $u_n=x^n,\ x\in[0,1)$，逐点收敛到 $0$ 但不一致（$x\to1$ 时误差趋 $1$）。

**11.2 Weierstrass 判别法（M-判别法）**
定理：$|u_n(x)|\le M_n\ \forall x\in E$ 且 $\sum M_n$ 收敛，则 $\sum u_n$ 一致收敛。
推导：由 Cauchy 准则
$$\Big|\sum_{k=m}^n u_k(x)\Big|\le\sum_{k=m}^nM_k<\varepsilon,$$
右端与 $x$ 无关，故一致 Cauchy。
应用：幂级数在收敛区间内闭一致收敛。

**11.3 和函数的连续、可积、可导**
连续性：$u_n$ 连续且 $\sum u_n$ 一致收敛，则 $S$ 连续。推导：
$$|S(x)-S(x_0)|\le|S(x)-S_N(x)|+|S_N(x)-S_N(x_0)|+|S_N(x_0)-S(x_0)|,$$
取 $N$ 使首末项 $<\varepsilon/3$，中间用 $S_N$ 连续。
逐项积分：$u_n$ 连续、一致收敛，则 $\int_a^x S=\sum\int_a^x u_n$。
逐项求导：$u_n\in C^1$，$\sum u_n$ 在一点收敛，$\sum u_n'$ 一致收敛，则
$$S'(x)=\sum u_n'(x).$$
推导：对 $\sum u_n'$ 逐项积分得 $C^1$ 函数，与 $S$ 相差常数；由一点收敛定常数。
量纲：逐项求导不改变量纲，但需一致收敛保证交换极限合法。

**11.4 Dini 定理**
定理：$u_n$ 非负连续，$\sum u_n$ 逐点收敛到连续 $S$，则 $\sum u_n$ 一致收敛。
推导：余项 $R_n(x)=S(x)-S_n(x)$ 连续、非负、逐点递减到 $0$；由 Dini 定理（紧集上单调递减连续函数列一致收敛）得一致收敛。
应用：判定 Fourier 级数一致收敛。

**11.5 典型应用：级数定义函数的可微性**
$S(x)=\sum_{n=1}^\infty\frac{\sin nx}{n^3}$。因 $|u_n|\le n^{-3}$，M-判别一致收敛，$S$ 连续。又 $u_n'=\frac{\cos nx}{n^2}$，$|u_n'|\le n^{-2}$，一致收敛，故
$$S'(x)=\sum_{n=1}^\infty\frac{\cos nx}{n^2}.$$
量纲：各项无量纲。

### 第12章 幂级数

**12.1 收敛半径**
定义：$\sum a_n(x-x_0)^n$。Cauchy–Hadamard 公式
$$R=\frac{1}{\limsup_{n\to\infty}|a_n|^{1/n}}.$$
推导：根值判别
$$\limsup\sqrt[n]{|a_n(x-x_0)^n|}=|x-x_0|\limsup|a_n|^{1/n},$$
其 $<1$ 即收敛，故 $|x-x_0|<R$ 收敛、$>R$ 发散。
比值形式：若 $\lim|a_{n+1}/a_n|=\rho$，则 $R=1/\rho$。
收敛域：区间 $(x_0-R,x_0+R)$ 端点需单独判定。

**12.2 和函数的解析性质**
在收敛区间内，幂级数可逐项求导、逐项积分，且收敛半径不变：
$$S'(x)=\sum_{n=1}^\infty n a_n(x-x_0)^{n-1},\quad \int_{x_0}^xS(t)dt=\sum_{n=0}^\infty\frac{a_n}{n+1}(x-x_0)^{n+1}.$$
推导：在 $|x-x_0|\le r<R$ 上，$|na_nr^{n-1}|$ 构成收敛级数（因 $\limsup(na_n)^{1/n}=1/R$），由 M-判别一致收敛，故可逐项求导/积分。
应用：求 $\sum n x^n=\frac{x}{(1-x)^2}$（$|x|<1$）。

**12.3 Taylor 级数与余项**
定理：$f$ 在 $x_0$ 邻域无穷可微，则
$$f(x)=\sum_{n=0}^\infty\frac{f^{(n)}(x_0)}{n!}(x-x_0)^n\iff R_n(x)\to0.$$
$$R_n(x)=\frac{1}{n!}\int_{x_0}^x f^{(n+1)}(t)(x-t)^n dt.$$
推导（积分余项）：由 $f(x)=f(x_0)+\int_{x_0}^xf'(t)dt$，反复分部积分
$$\int_{x_0}^xf'(t)dt=\big[f'(t)(x-t)\big]_{x_0}^x+\int_{x_0}^xf''(t)(x-t)dt=\cdots.$$
Lagrange 余项由中值定理得 $R_n=\dfrac{f^{(n+1)}(\xi)}{(n+1)!}(x-x_0)^{n+1}$。
反例：$f=e^{-1/x^2},f(0)=0$ 的所有导数为零但 $f\ne0$，说明 $C^\infty$ 不保证解析。

**12.4 初等函数展开**
$$e^x=\sum_{n=0}^\infty\frac{x^n}{n!},\quad \sin x=\sum_{k=0}^\infty\frac{(-1)^kx^{2k+1}}{(2k+1)!},\quad \cos x=\sum_{k=0}^\infty\frac{(-1)^kx^{2k}}{(2k)!},$$
$$\ln(1+x)=\sum_{n=1}^\infty\frac{(-1)^{n-1}x^n}{n},\quad (1+x)^\alpha=\sum_{n=0}^\infty\binom{\alpha}{n}x^n.$$
推导（$e^x$）：$f^{(n)}=e^x$，$f^{(n)}(0)=1$，余项 $|R_n|\le\dfrac{e^{|x|}}{(n+1)!}|x|^{n+1}\to0$。
量纲：展开变量须无量纲；物理中先作特征量归一化。

**12.5 幂级数解微分方程**
例：$y''+y=0$。设 $y=\sum a_nx^n$，代入
$$\sum_{n=2}^\infty n(n-1)a_nx^{n-2}+\sum_{n=0}^\infty a_nx^n=0.$$
整理同幂次
$$\sum_{k=0}^\infty\big[(k+2)(k+1)a_{k+2}+a_k\big]x^k=0\ \Rightarrow\ a_{k+2}=-\frac{a_k}{(k+1)(k+2)}.$$
得 $a_{2k}=\frac{(-1)^ka_0}{(2k)!},a_{2k+1}=\frac{(-1)^ka_1}{(2k+1)!}$，故 $y=a_0\cos x+a_1\sin x$。
应用：力学振动、量子力学级数解（Frobenius 方法）。

**12.6 典型应用：近似计算与误差**
用交错级数估计误差：$\sin x\approx x-\frac{x^3}{6}$，误差 $\le\frac{|x|^5}{120}$。推导：Leibniz 余项界。
量纲检查：$x$ 无量纲（弧度）。

### 第13章 傅里叶级数

**13.1 三角函数系的正交性**
在 $[-\pi,\pi]$ 上
$$\int_{-\pi}^\pi\cos mx\cos nx\,dx=\pi\delta_{mn}\ (m,n\ge1),\quad \int_{-\pi}^\pi\sin mx\sin nx\,dx=\pi\delta_{mn},$$
$$\int_{-\pi}^\pi\sin mx\cos nx\,dx=0,\quad \int_{-\pi}^\pi1\cdot\cos nx\,dx=\int_{-\pi}^\pi1\cdot\sin nx\,dx=0.$$
推导：用积化和差
$$\cos mx\cos nx=\tfrac12[\cos(m-n)x+\cos(m+n)x],$$
当 $m\ne n$ 时积分均为零；$m=n$ 时第一项积分 $2\pi$，第二项为零，故为 $\pi$。

**13.2 傅里叶系数的推导**
设 $f(x)=\frac{a_0}{2}+\sum_{n=1}^\infty(a_n\cos nx+b_n\sin nx)$ 且可逐项积分。
两边乘以 $\cos kx$ 后积分，由正交性仅 $a_k$ 项存活：
$$\because \int_{-\pi}^\pi f(x)\cos kx\,dx=a_k\int_{-\pi}^\pi\cos^2kx\,dx=\pi a_k,$$
$$\therefore a_k=\frac1\pi\int_{-\pi}^\pi f(x)\cos kx\,dx.$$
同理
$$b_k=\frac1\pi\int_{-\pi}^\pi f(x)\sin kx\,dx,\quad a_0=\frac1\pi\int_{-\pi}^\pi f(x)dx.$$
周期 $2l$ 推广：$\dfrac1l\int_{-l}^l$，频率 $n\pi/l$。
量纲：$[a_n]=[f]$（系数与原函数同量纲）。

**13.3 收敛定理（Dirichlet）**
定理：$f$ 以 $2\pi$ 为周期，在 $[-\pi,\pi]$ 上分片光滑（至多有限个第一类间断与极值），则其 Fourier 级数处处收敛于
$$\frac{f(x^+)+f(x^-)}{2}.$$
推导要点：部分和用 Dirichlet 核表示
$$S_N(x)=\frac1{2\pi}\int_{-\pi}^\pi f(x-t)D_N(t)dt,\quad D_N(t)=\frac{\sin(N+\frac12)t}{\sin\frac t2},$$
利用 $\frac1{2\pi}\int_{-\pi}^\pi D_N=\frac12$ 及 Riemann–Lebesgue 引理，振荡积分趋于零，只剩跳跃平均值。
特例：连续点收敛到 $f(x)$；跳跃点收敛到中点。

**13.4 Parseval 等式**
$$\frac{a_0^2}{2}+\sum_{n=1}^\infty(a_n^2+b_n^2)=\frac1\pi\int_{-\pi}^\pi f^2(x)\,dx.$$
推导：对 $f$ 与自身内积：将级数代入 $\frac1\pi\int f^2$，由正交性交叉项消失，得 $\frac{a_0^2}{2}+\sum(a_n^2+b_n^2)$。
物理意义：能量守恒——时域能量等于各谐波能量之和。
特例：$f=x$ 在 $[-\pi,\pi]$ 上给 $\sum_{n=1}^\infty\frac1{n^2}=\frac{\pi^2}{6}$。

**13.5 正弦级数与余弦级数**
奇延拓得正弦级数 $b_n=\frac2\pi\int_0^\pi f(x)\sin nx\,dx$；偶延拓得余弦级数 $a_n=\frac2\pi\int_0^\pi f(x)\cos nx\,dx$。
推导：奇函数 $a_n=0$，偶函数 $b_n=0$。
应用：热传导（两端固定温度的杆）用正弦级数；Neumann 边界用余弦级数。

**13.6 复指数形式**
$$f(x)=\sum_{n=-\infty}^\infty c_ne^{inx},\quad c_n=\frac1{2\pi}\int_{-\pi}^\pi f(x)e^{-inx}dx.$$
与三角形式关系：$c_0=a_0/2,\ c_n=(a_n-ib_n)/2,\ c_{-n}=(a_n+ib_n)/2$。推导：由 Euler 公式 $e^{inx}=\cos nx+i\sin nx$ 展开比较。
应用：信号频谱、量子力学动量表象。

**13.7 典型应用：一维热传导分离变量**
热方程 $u_t=ku_{xx}$，边界 $u(0,t)=u(L,t)=0$，初值 $u(x,0)=u_0(x)$。设 $u=X(x)T(t)$，得
$$\frac{T'}{kT}=\frac{X''}{X}=-\lambda\ \Rightarrow\ X''+\lambda X=0,\ T'=-\lambda kT.$$
$X_n=\sin\frac{n\pi x}{L},\lambda_n=\frac{n^2\pi^2}{L^2}$，$T_n=e^{-\lambda_nkt}$。
$$u(x,t)=\sum_{n=1}^\infty b_n\sin\frac{n\pi x}{L}e^{-n^2\pi^2kt/L^2},$$
$$b_n=\frac2L\int_0^Lu_0(x)\sin\frac{n\pi x}{L}dx.$$
量纲：$\lambda_n$ 量纲 $\mathrm{L^{-2}}$，$k\lambda_n t$ 无量纲。

### 第14章 常微分方程与力学应用

**14.1 一阶变量分离与线性方程**
可分离：$\dfrac{dy}{dx}=f(x)g(y)$，则 $\int\dfrac{dy}{g(y)}=\int f(x)dx+C$。
一阶线性 $\dfrac{dy}{dx}+P(x)y=Q(x)$。积分因子 $\mu=e^{\int P dx}$：
$$\because (\mu y)'=\mu y'+\mu Py=\mu(y'+Py)=\mu Q,\ \therefore y=e^{-\int Pdx}\Big(\int Qe^{\int Pdx}dx+C\Big).$$
应用：阻尼运动 $m\dot v=-kv+mg$。

**14.2 恰当方程与积分因子**
$M(x,y)dx+N(x,y)dy=0$ 恰当 $\iff M_y=N_x$，此时存在 $F$ 使 $F_x=M,F_y=N$，解为 $F=C$。
推导：由 $F_x=M$ 积分 $F=\int M dx+\phi(y)$，令 $F_y=N$ 定 $\phi'$。
量纲：$M\,dx$ 与 $N\,dy$ 量纲须一致。

**14.3 二阶常系数线性齐次方程**
$y''+py'+qy=0$。设 $y=e^{rx}$ 得特征方程
$$r^2+pr+q=0.$$
实根 $r_1\ne r_2$：$y=C_1e^{r_1x}+C_2e^{r_2x}$。
重根 $r$：$y=(C_1+C_2x)e^{rx}$。推导：设 $y=ue^{rx}$，代入得 $u''=0$，故 $u$ 线性。
复根 $\alpha\pm i\beta$：$y=e^{\alpha x}(C_1\cos\beta x+C_2\sin\beta x)$。
应用：阻尼振子 $m\ddot x+c\dot x+kx=0$。
$$\ddot x+2\zeta\omega_0\dot x+\omega_0^2x=0,\quad \omega_0=\sqrt{k/m},\ \zeta=\frac{c}{2\sqrt{km}}.$$
$\zeta<1$ 欠阻尼，$\zeta=1$ 临界，$\zeta>1$ 过阻尼。量纲：$[\omega_0]=\mathrm{T^{-1}}$，$\zeta$ 无量纲。

**14.4 常数变易法**
对 $y''+py'+qy=f(x)$，将齐次解 $C_1y_1+C_2y_2$ 中常数变为函数，令
$$y=u_1y_1+u_2y_2,\quad u_1'y_1+u_2'y_2=0.$$
推导：代入使 $y''$ 不出现 $u_i''$，附加条件消去 $u'$ 项，得
$$u_1'y_1'+u_2'y_2'=f,$$
解得
$$u_1'=-\frac{y_2f}{W},\quad u_2'=\frac{y_1f}{W},\quad W=y_1y_2'-y_1'y_2.$$
$$y_p=-y_1\int\frac{y_2f}{W}dx+y_2\int\frac{y_1f}{W}dx.$$
量纲：$[W]=[y]^2[x]^{-1}$，$[y_p]=[y]$。

**14.5 线性微分方程组与矩阵指数**
$\dot{\mathbf x}=A\mathbf x$，$\mathbf x(0)=\mathbf x_0$，解
$$\mathbf x(t)=e^{At}\mathbf x_0,\quad e^{At}=\sum_{k=0}^\infty\frac{A^kt^k}{k!}.$$
推导：逐项求导
$$\frac{d}{dt}e^{At}=\sum_{k=1}^\infty\frac{A^kt^{k-1}}{(k-1)!}=Ae^{At}.$$
特征分解：$A=V\Lambda V^{-1}\Rightarrow e^{At}=Ve^{\Lambda t}V^{-1}$。
应用：耦合振子、刚体姿态运动学。

**14.6 典型应用：受迫振动与共振**
$\ddot x+\omega_0^2x=F\cos\omega t$。稳态特解
$$x_p=\frac{F}{\omega_0^2-\omega^2}\cos\omega t.$$
推导：设 $x_p=A\cos\omega t$，代入得 $A=\frac{F}{\omega_0^2-\omega^2}$。当 $\omega\to\omega_0$ 共振振幅无界；有阻尼时
$$A(\omega)=\frac{F}{\sqrt{(\omega_0^2-\omega^2)^2+4\zeta^2\omega_0^2\omega^2}},$$
在 $\omega_r=\omega_0\sqrt{1-2\zeta^2}$ 处最大。

### 第15章 变分法基础

**15.1 泛函与变分**
定义：泛函 $J[y]=\int_{x_1}^{x_2}L(x,y,y')dx$。变分 $\delta y$ 为 admissible 扰动。
定义：若 $J[y+\delta y]-J[y]=\delta J+o(\|\delta y\|)$，称 $\delta J$ 为第一变分。
假设：$L\in C^2$，$y\in C^2$，端点固定 $\delta y(x_1)=\delta y(x_2)=0$。

**15.2 Euler–Lagrange 方程推导**
$$J[y+\delta y]-J[y]=\int_{x_1}^{x_2}\big[L_y\delta y+L_{y'}\delta y'\big]dx+o(\|\delta y\|).$$
对第二项分部积分
$$\int_{x_1}^{x_2}L_{y'}\delta y'\,dx=\big[L_{y'}\delta y\big]_{x_1}^{x_2}-\int_{x_1}^{x_2}\frac{d}{dx}L_{y'}\,\delta y\,dx.$$
端点消失，故
$$\delta J=\int_{x_1}^{x_2}\Big(L_y-\frac{d}{dx}L_{y'}\Big)\delta y\,dx.$$
由于 $\delta y$ 任意（基本引理），得
$$\boxed{\frac{\partial L}{\partial y}-\frac{d}{dx}\frac{\partial L}{\partial y'}=0}.$$
展开：$L_y-L_{xy'}-L_{yy'}y'-L_{y'y'}y''=0$，二阶 ODE。
几何/物理意义：真实路径使作用量取驻值；等价的 Newton 方程。

**15.3 高阶与多变量泛函**
含 $y''$：$\dfrac{\partial L}{\partial y}-\dfrac{d}{dx}\dfrac{\partial L}{\partial y'}+\dfrac{d^2}{dx^2}\dfrac{\partial L}{\partial y''}=0$。推导：两次分部积分。
多自变量 $u(x_1,\dots,x_n)$：$L=L(x,u,\nabla u)$，
$$\frac{\partial L}{\partial u}-\sum_i\frac{\partial}{\partial x_i}\frac{\partial L}{\partial u_{x_i}}=0=\frac{\partial L}{\partial u}-\nabla\cdot\frac{\partial L}{\partial\nabla u}.$$
推导：对 $\int L_{u_{x_i}}\delta u_{x_i}$ 用散度定理，边界项为零。
应用：拉普拉斯方程、极小曲面、弹性板。

**15.4 约束变分**
约束 $\int G(x,y,y')dx=C$，引入乘子
$$J^*=\int(L-\lambda G)dx,$$
Euler–Lagrange 方程用于 $L^*=L-\lambda G$。推导：约束流形上取极值，Lagrange 乘子法（函数空间版本）。
应用：等周问题、定长悬链线。

**15.5 Beltrami 恒等式与首次积分**
若 $L$ 不显含 $x$，则
$$L-y'L_{y'}=C.$$
推导：
$$\frac{d}{dx}(L-y'L_{y'})=L_xy'+L_yy'+L_{y'}y''-\frac{d}{dx}(y'L_{y'})=y'(L_y-\frac{d}{dx}L_{y'})=0.$$
应用：最速降线、悬链线。
最速降线：$L=\sqrt{(1+y'^2)/(2gy)}$，得 $y(1+y'^2)=C$，解为摆线
$$x=\frac C2(\theta-\sin\theta),\quad y=\frac C2(1-\cos\theta).$$
悬链线：$L=y\sqrt{1+y'^2}$，得 $y=C\cosh\frac{x-a}{C}$。

**15.6 典型应用：最小作用量原理**
力学系统 $L=T-V$，作用量 $S=\int_{t_1}^{t_2}L(q,\dot q,t)dt$。$\delta S=0$ 给 Lagrange 方程（见第16章）。
量纲：$[S]=\mathrm{M\,L^2T^{-1}}$（作用量），与 $\hbar$ 同量纲，是量子经典对应的桥梁。

### 第16章 拉格朗日力学与哈密顿力学

**16.1 广义坐标与完整约束**
含 $N$ 个质点、$s$ 个完整约束，自由度 $n=3N-s$。广义坐标 $q_1,\dots,q_n$。
虚功原理：
$$\sum_{i=1}^N(\mathbf F_i-m_i\ddot{\mathbf r}_i)\cdot\delta\mathbf r_i=0.$$
理想约束下约束力虚功为零。推导：D'Alembert 原理。

**16.2 拉格朗日方程推导**
由虚功原理将 $\delta\mathbf r_i=\sum_j\frac{\partial\mathbf r_i}{\partial q_j}\delta q_j$ 代入，并利用
$$\dot{\mathbf r}_i=\sum_j\frac{\partial\mathbf r_i}{\partial q_j}\dot q_j+\frac{\partial\mathbf r_i}{\partial t}.$$
推导关键恒等式：
$$\frac{\partial\dot{\mathbf r}_i}{\partial\dot q_j}=\frac{\partial\mathbf r_i}{\partial q_j},\quad \frac{d}{dt}\frac{\partial\mathbf r_i}{\partial q_j}=\frac{\partial\dot{\mathbf r}_i}{\partial q_j}.$$
动能 $T=\frac12\sum m_i\dot{\mathbf r}_i^2$，则
$$\sum_im_i\ddot{\mathbf r}_i\cdot\frac{\partial\mathbf r_i}{\partial q_j}=\frac{d}{dt}\frac{\partial T}{\partial\dot q_j}-\frac{\partial T}{\partial q_j}.$$
广义力 $Q_j=\sum_i\mathbf F_i\cdot\frac{\partial\mathbf r_i}{\partial q_j}=-\frac{\partial V}{\partial q_j}$（保守力），故
$$\boxed{\frac{d}{dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}=0},\quad L=T-V.$$
物理意义：广义动量 $p_j=\partial L/\partial\dot q_j$，方程即广义力动量定理。

**16.3 哈密顿正则方程**
定义广义动量 $p_j=\dfrac{\partial L}{\partial\dot q_j}$。假设 Legendre 变换可逆，即 Hessian $\partial^2L/\partial\dot q_i\partial\dot q_j$ 非奇异。
Hamiltonian：
$$H(q,p,t)=\sum_jp_j\dot q_j-L(q,\dot q,t),$$
其中 $\dot q=\dot q(q,p,t)$ 由 $p_j=\partial L/\partial\dot q_j$ 反解。
推导：
$$dH=\sum_j\dot q_jdp_j+p_jd\dot q_j-\sum_j\frac{\partial L}{\partial q_j}dq_j-\sum_j\frac{\partial L}{\partial\dot q_j}d\dot q_j-\frac{\partial L}{\partial t}dt.$$
由 $p_j=\partial L/\partial\dot q_j$ 消去 $d\dot q_j$，并用 Lagrange 方程 $\partial L/\partial q_j=\dot p_j$：
$$dH=\sum_j\Big(\dot q_jdp_j-\dot p_jdq_j\Big)-\frac{\partial L}{\partial t}dt.$$
另一方面 $dH=\sum_j(\partial_{q_j}H\,dq_j+\partial_{p_j}H\,dp_j)+\partial_tH\,dt$，比较得
$$\boxed{\dot q_j=\frac{\partial H}{\partial p_j},\quad \dot p_j=-\frac{\partial H}{\partial q_j},\quad \frac{\partial H}{\partial t}=-\frac{\partial L}{\partial t}}.$$
物理意义：相空间 $(q,p)$ 中 Hamilton 向量场，流保持辛结构。
量纲：$[p_j]=[L][\dot q]^{-1}=\mathrm{M\,L^2T^{-1}}$（动量或角动量），$[H]=\mathrm{M\,L^2T^{-2}}$（能量）。

**16.4 守恒律与Noether定理（初步）**
若 $L$ 不显含 $q_j$（循环坐标），则 $p_j$ 守恒：由 Lagrange 方程 $\dot p_j=\partial L/\partial q_j=0$。
若 $L$ 不显含 $t$，则 $H=\sum p_j\dot q_j-L$ 守恒：
$$\frac{dH}{dt}=\sum_j\Big(\dot q_j\dot p_j+p_j\ddot q_j-\frac{\partial L}{\partial q_j}\dot q_j-\frac{\partial L}{\partial\dot q_j}\ddot q_j\Big)-\frac{\partial L}{\partial t}=-\frac{\partial L}{\partial t}=0.$$
Noether 定理：每一连续对称性对应一个守恒量。时间平移→能量；空间平移→动量；旋转→角动量。
推导思路：对单参数变换 $q\to q+\epsilon\xi$ 使 $L$ 不变，Euler-Lagrange 给出 $\frac{d}{dt}(L_{\dot q}\xi)=0$。

**16.5 泊松括号与运动方程**
$$\{f,g\}=\sum_j\Big(\frac{\partial f}{\partial q_j}\frac{\partial g}{\partial p_j}-\frac{\partial f}{\partial p_j}\frac{\partial g}{\partial q_j}\Big).$$
运动方程：$\dot f=\{f,H\}+\dfrac{\partial f}{\partial t}$。
推导：由 Hamilton 方程
$$\dot f=\sum_j(f_{q_j}\dot q_j+f_{p_j}\dot p_j)+f_t=\sum_j(f_{q_j}H_{p_j}-f_{p_j}H_{q_j})+f_t=\{f,H\}+f_t.$$
性质：反对称 $\{f,g\}=-\{g,f\}$，Jacobi 恒等式，Leibniz 法则。
量子对应：$\{f,g\}\to\frac{1}{i\hbar}[\hat f,\hat g]$。

**16.6 典型应用：二体问题与中心力**
相对坐标 $\mathbf r$，约化质量 $\mu=m_1m_2/(m_1+m_2)$，
$$L=\frac12\mu\dot{\mathbf r}^2-V(r).$$
球坐标下角动量 $p_\theta=\mu r^2\dot\theta$ 守恒（循环坐标 $\theta$），得
$$r^2\dot\theta=h=\mathrm{const},$$
即 Kepler 第二定律。有效势 $V_{\rm eff}=V(r)+\dfrac{\mu h^2}{2r^2}$，稳定圆轨道由 $V_{\rm eff}'=0$ 给出。

### 第17章 稳定性、线性化与分岔初步

**17.1 平衡点与Lyapunov稳定性**
系统 $\dot{\mathbf x}=\mathbf f(\mathbf x)$，平衡点 $\mathbf f(\mathbf x_0)=0$。
定义：Lyapunov 稳定：$\forall\varepsilon>0,\exists\delta>0:\|\mathbf x(0)-\mathbf x_0\|<\delta\Rightarrow\|\mathbf x(t)-\mathbf x_0\|<\varepsilon\ \forall t\ge0$。渐近稳定：稳定且 $\mathbf x(t)\to\mathbf x_0$。
意义：结构受小扰动是否返回平衡。

**17.2 线性化与特征值判据**
令 $\mathbf x=\mathbf x_0+\mathbf u$，展开
$$\dot{\mathbf u}=A\mathbf u+O(\|\mathbf u\|^2),\quad A=D\mathbf f(\mathbf x_0).$$
定理（Lyapunov 间接法）：若 $A$ 的全部特征值有负实部，则 $\mathbf x_0$ 渐近稳定；若存在正实部特征值，则不稳定。
推导：解 $\mathbf u(t)=e^{At}\mathbf u(0)$，由 Jordan 形 $\|e^{At}\|\le Ce^{\alpha t}$，$\alpha=\max\mathrm{Re}\lambda_i<0$。
临界情形 $\mathrm{Re}\lambda=0$ 需非线性项判断。

**17.3 Lyapunov函数法**
定理：若存在 $V(\mathbf x)$ 在 $\mathbf x_0$ 邻域正定、$V(\mathbf x_0)=0$、且 $\dot V=\nabla V\cdot\mathbf f\le0$，则稳定；若 $\dot V<0$（$\mathbf x\ne\mathbf x_0$）则渐近稳定。
推导：$\dot V\le0\Rightarrow V(\mathbf x(t))\le V(\mathbf x(0))$；由正定，$\|\mathbf x\|$ 被 $V$ 控制，故不能远离。
特例：机械能 $E=T+V$，阻尼时 $\dot E=-c\dot q^2\le0$，故阻尼系统稳定。

**17.4 结构屈曲的特征值问题**
压杆平衡：$EIw''''+Pw''=0$，边界 $w(0)=w(L)=0$，$w''(0)=w''(L)=0$。
设 $w=e^{rx}$ 得 $r^2(EIr^2+P)=0$，非平凡解 $w=A\sin kx+B\cos kx+Cx+D$。边界给出
$$\sin kL=0\ \Rightarrow\ kL=n\pi,\quad P_n=\frac{n^2\pi^2EI}{L^2}.$$
意义：$P_1$ 为分岔点，平衡解由直线分岔为弯曲；$n$ 为模态阶数。
量纲：$EI$ 为 $\mathrm{M\,L^3T^{-2}}$，除以 $L^2$ 为力。

**17.5 分岔与中心流形初步**
在临界参数处，平衡点数目或稳定性改变。Pitchfork 分岔范式
$$\dot x=\mu x-x^3.$$
平衡 $x=0$ 与 $x=\pm\sqrt\mu$（$\mu>0$）。稳定性：$\mu<0$ 时 $x=0$ 稳定；$\mu>0$ 时 $x=0$ 不稳定而 $\pm\sqrt\mu$ 稳定。
推导：线性化 $\lambda=\mu$，故 $\mu=0$ 为分岔点；非线性项决定超临界/亚临界。

**17.6 典型应用：单摆稳定性**
$\ddot\theta+\frac gL\sin\theta=0$，状态 $(\theta,\dot\theta)$，$\mathbf f=(\dot\theta,-\frac gL\sin\theta)$。
平衡 $(0,0)$ 与 $(\pi,0)$。
$$A=\begin{pmatrix}0&1\\-\frac gL\cos\theta&0\end{pmatrix}.$$
在 $(0,0)$：$\lambda^2=-g/L$ 纯虚，中心；能量法证明稳定。
在 $(\pi,0)$：$\lambda^2=+g/L$，一正一负，鞍点，不稳定。
量纲：$[g/L]=\mathrm{T^{-2}}$，$\lambda$ 为 $\mathrm{T^{-1}}$。

### 第18章 连续介质力学守恒律与Navier–Stokes方程

**18.1 连续介质假设与雷诺输运定理**
假设：存在尺度 $\ell$ 使宏观量可定义且连续可微；场量 $\rho(\mathbf x,t),\mathbf v(\mathbf x,t)$。
雷诺输运定理：对随流体运动的控制体 $\Omega(t)$，
$$\frac{d}{dt}\iiint_{\Omega(t)}F(\mathbf x,t)\,dV=\iiint_{\Omega(t)}\Big(\frac{\partial F}{\partial t}+\nabla\cdot(F\mathbf v)\Big)dV.$$
推导：将 $\Omega(t)$ 拉回到参考构形 $\Omega_0$，$dV=J\,dV_0$，$\dot J=J\nabla\cdot\mathbf v$（见下）。
$$\frac{d}{dt}\int_{\Omega_0}FJ\,dV_0=\int_{\Omega_0}\big(\dot FJ+F\dot J\big)dV_0=\int_{\Omega_0}\big(\dot F+F\nabla\cdot\mathbf v\big)J\,dV_0.$$
物质导数 $\dot F=\partial_tF+\mathbf v\cdot\nabla F$，代入得
$$\dot F+F\nabla\cdot\mathbf v=\partial_tF+\nabla\cdot(F\mathbf v).$$
行列式导数引理推导：$\dot J=J\,\mathrm{tr}(F^{-1}\dot F)=J\nabla\cdot\mathbf v$，其中 $F=\partial\mathbf x/\partial\mathbf X$ 为变形梯度。

**18.2 质量守恒（连续性方程）**
取 $F=\rho$，质量守恒 $\frac{d}{dt}\int_{\Omega(t)}\rho\,dV=0$。由雷诺定理
$$\int_{\Omega(t)}\Big(\partial_t\rho+\nabla\cdot(\rho\mathbf v)\Big)dV=0\ \forall\Omega(t)\ \Rightarrow\ \partial_t\rho+\nabla\cdot(\rho\mathbf v)=0.$$
展开：$\partial_t\rho+\mathbf v\cdot\nabla\rho+\rho\nabla\cdot\mathbf v=0$，即 $\dot\rho+\rho\nabla\cdot\mathbf v=0$。
不可压：$\rho$ 常数，$\nabla\cdot\mathbf v=0$。
量纲：$[\rho]=\mathrm{M\,L^{-3}}$，$[\rho\mathbf v]=\mathrm{M\,L^{-2}T^{-1}}$，$[\nabla\cdot(\rho\mathbf v)]=\mathrm{M\,L^{-3}T^{-1}}$。

**18.3 动量守恒与Cauchy应力方程**
动量变化率等于表面力与体积力之和：
$$\frac{d}{dt}\iiint_{\Omega(t)}\rho\mathbf v\,dV=\oiint_{\partial\Omega(t)}\boldsymbol\sigma\cdot\mathbf n\,dS+\iiint_{\Omega(t)}\rho\mathbf b\,dV.$$
左边由雷诺定理
$$\frac{d}{dt}\int\rho\mathbf v\,dV=\int\Big(\partial_t(\rho\mathbf v)+\nabla\cdot(\rho\mathbf v\otimes\mathbf v)\Big)dV.$$
右边由散度定理
$$\oiint_{\partial\Omega}\boldsymbol\sigma\cdot\mathbf n\,dS=\iiint_\Omega\nabla\cdot\boldsymbol\sigma\,dV.$$
利用连续性方程可将左边化为 $\int\rho\dot{\mathbf v}\,dV$：
$$\partial_t(\rho\mathbf v)+\nabla\cdot(\rho\mathbf v\otimes\mathbf v)=\rho\partial_t\mathbf v+\mathbf v[\partial_t\rho+\nabla\cdot(\rho\mathbf v)]+\rho\mathbf v\cdot\nabla\mathbf v=\rho\dot{\mathbf v}.$$
故得 Cauchy 动量方程
$$\boxed{\rho\dot{\mathbf v}=\nabla\cdot\boldsymbol\sigma+\rho\mathbf b},\quad \dot{\mathbf v}=\partial_t\mathbf v+(\mathbf v\cdot\nabla)\mathbf v.$$
量纲检验：$[\rho\dot{\mathbf v}]=\mathrm{M\,L^{-2}T^{-2}}$，$[\nabla\cdot\boldsymbol\sigma]$ 同，$[\boldsymbol\sigma]=\mathrm{M\,L^{-1}T^{-2}}=\mathrm{Pa}$。

**18.4 角动量守恒与应力对称性**
角动量方程
$$\frac{d}{dt}\int\Omega\times(\rho\mathbf v)\,dV=\int\mathbf r\times(\rho\mathbf b)dV+\oint\mathbf r\times(\boldsymbol\sigma\mathbf n)dS.$$
推导：对 $\oint r_i\sigma_{jk}n_k$ 用散度定理，动量方程代入，得
$$\sigma_{ij}=\sigma_{ji}.$$
推导细节：$\varepsilon_{ijk}\sigma_{jk}=0\Rightarrow$ 反对称部分为零。
物理意义：无体力偶时应力张量对称，独立分量从 9 降为 6。

**18.5 本构方程与Navier–Stokes方程**
牛顿流体假设：应力为压力项与粘性项之和，粘性项线性依赖于速度梯度：
$$\boldsymbol\sigma=-p\mathbf I+\boldsymbol\tau,\quad \boldsymbol\tau=2\mu\mathbf D+\lambda(\nabla\cdot\mathbf v)\mathbf I,$$
$$\mathbf D=\frac12(\nabla\mathbf v+\nabla\mathbf v^\top),\quad D_{ij}=\frac12(\partial_jv_i+\partial_iv_j).$$
假设：各向同性、Galilean 不变、Stokes 假设 $3\lambda+2\mu=0$，则 $\lambda=-2\mu/3$；不可压时 $\nabla\cdot\mathbf v=0$，$\boldsymbol\tau=2\mu\mathbf D$。
代入 Cauchy 方程（不可压、$\mu$ 常数）：
$$\rho(\partial_t\mathbf v+(\mathbf v\cdot\nabla)\mathbf v)=-\nabla p+\mu\nabla\cdot(\nabla\mathbf v+\nabla\mathbf v^\top)+\rho\mathbf b.$$
因 $\nabla\cdot(\nabla\mathbf v^\top)=\nabla(\nabla\cdot\mathbf v)=0$（不可压），且 $\nabla\cdot\nabla\mathbf v=\Delta\mathbf v$，得
$$\boxed{\rho\Big(\frac{\partial\mathbf v}{\partial t}+(\mathbf v\cdot\nabla)\mathbf v\Big)=-\nabla p+\mu\Delta\mathbf v+\rho\mathbf b,\quad \nabla\cdot\mathbf v=0}.$$
量纲：$[\mu]=\mathrm{M\,L^{-1}T^{-1}}=\mathrm{Pa\cdot s}$，$[\nu]=\mu/\rho=\mathrm{L^2T^{-1}}$。

**18.6 能量方程、耗散与Bernoulli方程**
动能定理：对不可压 N-S 点乘 $\mathbf v$ 并积分。先看局部：
$$\rho\mathbf v\cdot\dot{\mathbf v}=\mathbf v\cdot\nabla\cdot\boldsymbol\sigma+\rho\mathbf v\cdot\mathbf b.$$
利用 $\mathbf v\cdot(\mathbf v\cdot\nabla)\mathbf v=\frac12\mathbf v\cdot\nabla|\mathbf v|^2$ 及
$$\mathbf v\cdot\nabla\cdot\boldsymbol\sigma=\nabla\cdot(\boldsymbol\sigma^\top\mathbf v)-\boldsymbol\sigma:\nabla\mathbf v,$$
$$\boldsymbol\sigma:\nabla\mathbf v=-p\nabla\cdot\mathbf v+2\mu\mathbf D:\mathbf D=-p\nabla\cdot\mathbf v+2\mu\|\mathbf D\|^2,$$
得
$$\frac{\partial}{\partial t}\Big(\frac12\rho|\mathbf v|^2\Big)+\nabla\cdot\Big[\mathbf v\Big(\frac12\rho|\mathbf v|^2+p\Big)-2\mu\mathbf D\mathbf v\Big]=-2\mu\|\mathbf D\|^2+\rho\mathbf v\cdot\mathbf b.$$
量纲：$2\mu\|\mathbf D\|^2$ 为 $\mathrm{M\,L^{-1}T^{-3}}$（单位体积耗散率）。
理想无粘、定常、不可压、无旋、沿流线时，
$$\frac12\rho|\mathbf v|^2+p+\rho gz=\mathrm{const},$$
即 Bernoulli 方程。推导：对 $\frac{\mathbf v^2}{2}+\frac p\rho+gz$ 沿流线求导为零。

**18.7 无量纲化与雷诺数**
取特征长度 $L$、速度 $U$、时间 $L/U$、压力 $\rho U^2$。令 $\mathbf v=U\mathbf v^*,\ \mathbf x=L\mathbf x^*,\ t=(L/U)t^*,\ p=\rho U^2p^*$，代入 N-S：
$$\rho\frac{U^2}{L}\Big(\partial_{t^*}\mathbf v^*+(\mathbf v^*\cdot\nabla^*)\mathbf v^*\Big)=-\frac{\rho U^2}{L}\nabla^*p^*+\frac{\mu U}{L^2}\Delta^*\mathbf v^*,$$
$$\Rightarrow\ \partial_{t^*}\mathbf v^*+(\mathbf v^*\cdot\nabla^*)\mathbf v^*=-\nabla^*p^*+\frac{1}{Re}\Delta^*\mathbf v^*,\quad Re=\frac{\rho UL}{\mu}.$$
物理意义：$Re$ 为惯性力与粘性力之比；$Re\ll1$ 层流/Stokes，$Re\gg1$ 边界层与湍流。
量纲：$Re$ 无量纲。

**18.8 典型应用：管道Poiseuille流**
圆管半径 $R$，定常层流 $\mathbf v=(0,0,u(r))$，压差驱动 $\partial_zp=-G$。
N-S 化为
$$-G+\mu\frac1r\frac{d}{dr}\Big(r\frac{du}{dr}\Big)=0.$$
积分：
$$r\frac{du}{dr}=-\frac{G}{2\mu}r^2+C_1\ \Rightarrow\ \frac{du}{dr}=-\frac{G}{2\mu}r+\frac{C_1}{r},$$
有界性要求 $C_1=0$；$u(R)=0$ 给
$$u(r)=\frac{G}{4\mu}(R^2-r^2).$$
体积流量
$$Q=\int_0^R2\pi ru(r)dr=\frac{\pi GR^4}{8\mu},$$
平均速度 $\bar u=GR^2/(8\mu)$，壁面切应力 $\tau_w=GR/2$。
量纲：$G=\mathrm{Pa\,m^{-1}}$，$GR^4/\mu=\mathrm{m^4\,s^{-1}}$，故 $Q$ 为 $\mathrm{m^3\,s^{-1}}$，正确。

---
⬆ 上级领域：[[数学]]
