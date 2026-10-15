---
课程名称: 数学分析I
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

# 数学分析I


### 第1章 实数系的完备性与函数结构

**1.1 实数域与有序域公理**

定义：实数集 $\mathbb{R}$ 是含加、乘两种运算的集合，满足域公理、全序公理，以及完备性（确界）公理。

假设：$\mathbb{R}$ 上定义了全序 $\le$，且与运算相容：$a\le b\Rightarrow a+c\le b+c$，$a\le b,\,c\ge0\Rightarrow ac\le bc$。

完备性公理：任意非空有上界的子集 $S\subseteq\mathbb{R}$ 在 $\mathbb{R}$ 中存在上确界 $\sup S$。

由完备性推 Archimedes 性质：$\forall x\in\mathbb{R},\,\exists n\in\mathbb{N}$ 使 $n>x$。

推导（反证）：若 $\mathbb{N}$ 有上界，设 $\alpha=\sup\mathbb{N}$，则 $\alpha-1$ 不是上界，存在 $n\in\mathbb{N}$ 使 $n>\alpha-1$，从而 $n+1>\alpha$，与 $\alpha$ 为上界矛盾。$\therefore$ 无上界。

几何意义：实数轴无“空隙”，这是极限存在、连续函数介值性、微积分基本定理的根基。

特例：有理数集 $\mathbb{Q}$ 满足有序域公理但不完备，如 $S=\{q\in\mathbb{Q}:q^2<2\}$ 在 $\mathbb{Q}$ 中无上确界。

典型应用：力学中一切连续运动的位置函数 $x(t)$ 都取值于完备的 $\mathbb{R}$；能量泛函极小点的存在性依赖确界原理。

**1.2 绝对值与三角不等式**

定义：$|x|=\begin{cases}x,&x\ge0\\-x,&x<0\end{cases}$。

三角不等式：$\forall a,b\in\mathbb{R}$，$|a+b|\le|a|+|b|$。

推导：由 $-|a|\le a\le|a|$，$-|b|\le b\le|b|$，两式相加得
$$-(|a|+|b|)\le a+b\le |a|+|b|\ \Rightarrow\ |a+b|\le|a|+|b|.$$
反向不等式：$|a-b|\ge\big||a|-|b|\big|$，由 $|a|=|(a-b)+b|\le|a-b|+|b|$ 得 $|a|-|b|\le|a-b|$，交换 $a,b$ 同理。

几何意义：折线路径不短于直线，是误差估计与收敛证明的工具。

典型应用：数值积分误差 $|I-I_n|\le\sum w_i|x_i-x_{i-1}|$ 的估计。

**1.3 确界、上下确界与唯一性**

定义：$\sup S$ 是最小上界；$\inf S$ 是最大下界。

推导（唯一性）：设 $s_1,s_2$ 均为 $\sup S$。$s_2$ 是上界 $\Rightarrow s_1\le s_2$；$s_1$ 是上界 $\Rightarrow s_2\le s_1$。$\therefore s_1=s_2$。

公式：$\sup(S+T)=\sup S+\sup T$（$S,T$ 非空有上界）。

几何意义：确界是“逼近可达的最小/最大界限”，是极限语言 $\varepsilon$ 叙述的母语。

特例：$S=(0,1)$ 时 $\sup S=1\notin S$，说明确界未必属于集合。

典型应用：最小势能原理中实际位移使总势能取极小，极小值存在即用确界。

**1.4 函数、复合与反函数**

定义：映射 $f:D\to\mathbb{R}$，$D\subseteq\mathbb{R}$。复合 $(g\circ f)(x)=g(f(x))$。

假设：$f:D\to E$ 为双射，则存在反函数 $f^{-1}:E\to D$，满足 $f^{-1}(f(x))=x,\ f(f^{-1}(y))=y$。

推导（反函数导数关系预告）：设 $y=f(x)$，$\Delta y=f(x+\Delta x)-f(x)$，则
$$\frac{\Delta x}{\Delta y}=\frac{1}{\Delta y/\Delta x}\ \Rightarrow\ (f^{-1})'(y)=\frac{1}{f'(x)}\quad(f'(x)\ne0).$$

几何意义：反函数图像关于直线 $y=x$ 对称。

特例：$y=\sin x$ 在 $[-\pi/2,\pi/2]$ 上双射，$\arcsin$ 的导数由上式得 $1/\cos y=1/\sqrt{1-x^2}$。

典型应用：由 $v(t)$ 求 $t(v)$，再由 $x(t(v))$ 描述相轨迹。

**1.5 初等函数与双曲函数恒等式**

定义：$\sinh x=\dfrac{e^x-e^{-x}}{2}$，$\cosh x=\dfrac{e^x+e^{-x}}{2}$。

公式与推导：
$$\cosh^2x-\sinh^2x=\frac{(e^x+e^{-x})^2-(e^x-e^{-x})^2}{4}=\frac{4}{4}=1.$$
$$\sinh(x+y)=\sinh x\cosh y+\cosh x\sinh y,$$
由定义展开 $e^{x+y}-e^{-(x+y)}$ 即得。导数：
$$(\sinh x)'=\cosh x,\qquad(\cosh x)'=\sinh x.$$

几何意义：$(\cosh t,\sinh t)$ 是双曲线 $u^2-v^2=1$ 的参数化，类比单位圆。

极限：$\sinh x\sim x$，$\cosh x\sim 1+\frac{x^2}{2}$（$x\to0$）。

典型应用：悬链线 $y=a\cosh(x/a)$ 是重力下均质柔索的形状，由后面变分法导出。

**1.6 参数方程、极坐标与周期运动**

定义：曲线 $\begin{cases}x=x(t)\\y=y(t)\end{cases}$；极坐标 $x=r\cos\theta,\ y=r\sin\theta$。

公式：若 $x'(t)\ne0$，则 $\dfrac{dy}{dx}=\dfrac{dy/dt}{dx/dt}=\dfrac{y'(t)}{x'(t)}$。

推导：由链式法则 $dy=dy/dt\,dt$，$dx=dx/dt\,dt$，相除即得。

几何意义：参数 $t$ 可为时间，$(x(t),y(t))$ 为质点在平面上的轨迹。

特例：圆周 $x=R\cos\omega t,\,y=R\sin\omega t$，则 $dy/dx=-\cot\omega t$。

典型应用：简谐振动 $x=A\cos(\omega t+\varphi)$ 的相图是椭圆，周期 $T=2\pi/\omega$。

### 第2章 数列极限与收敛判据

**2.1 数列极限的 $\varepsilon$-$N$ 定义**

定义：$\displaystyle\lim_{n\to\infty}a_n=a\iff\forall\varepsilon>0,\ \exists N\in\mathbb{N},\ \forall n>N:\ |a_n-a|<\varepsilon$。

假设：$a_n$ 为实数列。

推导（极限唯一性）：设 $a_n\to a$ 且 $a_n\to b$，$a\ne b$。取 $\varepsilon=|a-b|/2>0$，存在 $N$ 使 $n>N$ 时 $|a_n-a|<\varepsilon$ 且 $|a_n-b|<\varepsilon$，于是
$$|a-b|=|(a-a_n)+(a_n-b)|\le|a_n-a|+|a_n-b|<\varepsilon+\varepsilon=|a-b|,$$
矛盾。$\therefore a=b$。

推论：收敛数列有界。

几何意义：充分靠后的项全部落入 $a$ 的 $\varepsilon$-邻域。

极限或特例：常数列极限为其本身；$a_n=n$ 发散。

典型应用：迭代格式 $x_{n+1}=g(x_n)$ 的收敛性分析（如求解非线性方程）。

**2.2 单调有界定理**

定理：递增且有上界的数列收敛，且极限为上确界；递减且有下界的数列收敛。

推导：设 $\{a_n\}$ 递增有上界，$L=\sup\{a_n\}$。由 $\sup$ 定义，$\forall\varepsilon>0$，$L-\varepsilon$ 不是上界，故 $\exists N$ 使 $a_N>L-\varepsilon$。因递增，$n>N\Rightarrow a_n\ge a_N>L-\varepsilon$，又 $a_n\le L<L+\varepsilon$，故 $|a_n-L|<\varepsilon$。$\therefore a_n\to L$。

几何意义：单调逼近必收敛，即使不能显式求极限。

特例：$a_1=\sqrt2,\ a_{n+1}=\sqrt{2+a_n}$ 递增有上界 2，极限满足 $L=\sqrt{2+L}\Rightarrow L=2$。

典型应用：证明 $\big(1+1/n\big)^n$ 收敛（下面 2.7）。

**2.3 夹逼定理**

定理：若 $x_n\le z_n\le y_n$ 且 $x_n\to a,\ y_n\to a$，则 $z_n\to a$。

推导：$\forall\varepsilon>0$，取 $N_1,N_2$ 使 $n>N_1$ 时 $|x_n-a|<\varepsilon$，$n>N_2$ 时 $|y_n-a|<\varepsilon$。令 $N=\max(N_1,N_2)$，则
$$a-\varepsilon<x_n\le z_n\le y_n<a+\varepsilon\ \Rightarrow\ |z_n-a|<\varepsilon.$$

几何意义：被两条收敛到同一值的曲线“夹住”的曲线同收敛。

特例：$\lim_{n\to\infty}\dfrac{n}{2^n}=0$，用二项式 $2^n=(1+1)^n\ge n(n-1)/2$ 夹逼。

典型应用：估计数值算法的误差上界。

**2.4 Cauchy 收敛准则**

定理：$\{a_n\}$ 收敛 $\iff$ $\forall\varepsilon>0,\ \exists N,\ \forall m,n>N:\ |a_n-a_m|<\varepsilon$。

推导（必要性）：设 $a_n\to a$，取 $N$ 使 $n>N$ 时 $|a_n-a|<\varepsilon/2$，则 $m,n>N$ 时
$$|a_n-a_m|\le|a_n-a|+|a-a_m|<\varepsilon.$$
充分性依赖实数完备性（Bolzano–Weierstrass）。

几何意义：收敛的本性是项与项越来越靠近，无需预知极限。

特例：$a_n=\sum_{k=1}^n 1/k$ 发散，因 $a_{2n}-a_n\ge 1/2$ 不满足 Cauchy。

典型应用：判定级数收敛的理论基础。

**2.5 Stolz 定理**

定理（$\infty/\infty$ 型）：若 $y_n$ 严格递增且 $y_n\to+\infty$，且 $\dfrac{x_n-x_{n-1}}{y_n-y_{n-1}}\to L$，则 $\dfrac{x_n}{y_n}\to L$。

推导思路：令 $b_n=y_n-y_{n-1}>0$，$a_n=x_n-x_{n-1}$。对 $n$ 求和得 $x_n-x_0=\sum a_k$，$y_n-y_0=\sum b_k$。将 $a_k$ 用 $L b_k$ 加误差 $\varepsilon_k b_k$ 表示，求和后 $\sum\varepsilon_k b_k / \sum b_k$ 加权平均可任意小，即得。

特例：$\lim \dfrac{1+2+\cdots+n}{n^2}=\lim\dfrac{n}{2n-1}=\dfrac12$。

典型应用：求离散平均、对数增长等极限。

**2.6 两个基本极限与等价无穷小**

公式：$\displaystyle\lim_{x\to0}\frac{\sin x}{x}=1,\qquad \lim_{x\to\infty}\Big(1+\frac1x\Big)^x=e$。

$\sin x/x$ 推导（几何法）：单位圆中，$0<x<\pi/2$ 时面积满足
$$\frac12\sin x<\frac12 x<\frac12\tan x\ \Rightarrow\ \cos x<\frac{\sin x}{x}<1,$$
由 $\cos x\to1$ 及夹逼即得。

几何意义：小角近似 $\sin x\approx\tan x\approx x$，是振动线性化的根据。

特例：$\lim_{x\to0}(1-\cos x)/x^2=1/2$，由 $1-\cos x=2\sin^2(x/2)$ 得。

典型应用：单摆方程 $\ddot\theta+\frac{g}{L}\sin\theta=0$ 在小角下线性化为 $\ddot\theta+\frac{g}{L}\theta=0$。

**2.7 自然常数 $e$ 的推导**

定义：$e=\displaystyle\lim_{n\to\infty}\Big(1+\frac1n\Big)^n$。

推导：
$$\Big(1+\frac1n\Big)^n=\sum_{k=0}^n\binom nk\frac1{n^k}=\sum_{k=0}^n\frac{1}{k!}\prod_{j=1}^{k-1}\Big(1-\frac{j}{n}\Big).$$
当 $n$ 增大，每个因子 $1-j/n$ 增大，且项数增多，故数列递增。又有
$$\Big(1+\frac1n\Big)^n\le\sum_{k=0}^n\frac1{k!}\le1+1+\sum_{k=2}^n\frac1{2^{k-1}}<3,$$
故有上界，由单调有界定理收敛，记为 $e$。

数值：$e=2.718281828\ldots$。

典型应用：连续复利 $A=A_0e^{rt}$，一阶动力学 $e^{-\lambda t}$。

### 第3章 函数极限、连续与闭区间性质

**3.1 函数极限的 $\varepsilon$-$\delta$ 定义**

定义：$\displaystyle\lim_{x\to x_0}f(x)=A\iff\forall\varepsilon>0,\ \exists\delta>0,\ 0<|x-x_0|<\delta\Rightarrow|f(x)-A|<\varepsilon$。

推导（局部有界性）：设极限为 $A$，取 $\varepsilon=1$，存在 $\delta_1>0$ 使 $0<|x-x_0|<\delta_1$ 时 $|f(x)-A|<1$，故 $|f(x)|<|A|+1$。

几何意义：无论给出多窄的水平带 $(A-\varepsilon,A+\varepsilon)$，总能在 $x_0$ 附近找到去心邻域使图像落入其中。

特例：$\lim_{x\to0}x\sin(1/x)=0$ 由夹逼；$\lim_{x\to0}\sin(1/x)$ 不存在。

典型应用：导数定义为差商的函数极限。

**3.2 无穷小与无穷大阶的比较**

定义：$x\to x_0$ 时 $\alpha,\beta\to0$。若 $\alpha/\beta\to0$，记 $\alpha=o(\beta)$；若 $\alpha/\beta\to1$，记 $\alpha\sim\beta$；若 $\alpha/\beta$ 有界，记 $\alpha=O(\beta)$。

公式：$o(\beta)$ 运算封闭：$o(\beta)+o(\beta)=o(\beta)$，$o(\beta)\cdot O(\gamma)=o(\beta\gamma)$。

推导：设 $\alpha_1/\beta\to0,\alpha_2/\beta\to0$，则 $(\alpha_1+\alpha_2)/\beta\to0$。

几何意义：阶描述收敛速度，$O$ 上界 $o$ 严格更快。

特例：$x\to0$ 时 $\sin x=x+o(x)$，$1-\cos x=O(x^2)$。

典型应用：数值格式的局部截断误差，如 Euler 法 $x(t+h)=x(t)+h f+O(h^2)$。

**3.3 连续性、间断点分类**

定义：$f$ 在 $x_0$ 连续 $\iff \lim_{x\to x_0}f(x)=f(x_0)$。

间断点分类：可去（左右极限存在相等 $\ne f(x_0)$）；跳跃（左右极限存在不等）；第二类（至少一侧极限不存在）。

推导（连续函数四则运算）：若 $f,g$ 在 $x_0$ 连续，则 $fg$ 连续。由 $f(x)g(x)-f(x_0)g(x_0)=f(x)(g(x)-g(x_0))+g(x_0)(f(x)-f(x_0))$，利用局部有界与极限运算得。

特例：$f(x)=x$、常数连续；由乘积得多项式连续；再由商得有理函数在分母非零处连续。

典型应用：力学中若力为位置的连续函数，则解对初值连续（定性理论起点）。

**3.4 闭区间上连续函数的四大定理**

定理（有界性、最值性、介值性、一致连续性）：$f\in C[a,b]$。

推导（介值性）：设 $f(a)<\lambda<f(b)$，令 $S=\{x\in[a,b]:f(x)<\lambda\}$，$c=\sup S$。由连续性与确界性质可证 $f(c)=\lambda$。

推导（零点定理特例）：$\lambda=0$ 即得。

几何意义：连续曲线从下到上必穿过每条中间水平线。

特例：奇次多项式必有实根，因 $x\to\pm\infty$ 时符号相反。

典型应用：证明单摆方程存在周期解、方程 $x-e^{-x}=0$ 有根。

**3.5 一致连续与 Cantor 定理**

定义：$f$ 在 $I$ 一致连续 $\iff\forall\varepsilon>0,\exists\delta>0,\forall x,y\in I:|x-y|<\delta\Rightarrow|f(x)-f(y)|<\varepsilon$。

Cantor 定理：$f\in C[a,b]\Rightarrow f$ 在 $[a,b]$ 一致连续。

推导（反证）：若不一致连续，存在 $\varepsilon_0>0$ 和 $x_n,y_n$ 使 $|x_n-y_n|<1/n$ 而 $|f(x_n)-f(y_n)|\ge\varepsilon_0$。由有界性取子列 $x_{n_k}\to c$，则 $y_{n_k}\to c$，与连续矛盾。

特例：$f(x)=1/x$ 在 $(0,1)$ 上连续但不一致连续。

典型应用：保证积分和存在、函数序列一致收敛。

### 第4章 导数、微分与求导法则

**4.1 导数定义与几何意义**

定义：$f'(x_0)=\displaystyle\lim_{\Delta x\to0}\frac{f(x_0+\Delta x)-f(x_0)}{\Delta x}=\lim_{x\to x_0}\frac{f(x)-f(x_0)}{x-x_0}$。

推导（可导必连续）：
$$f(x)-f(x_0)=\frac{f(x)-f(x_0)}{x-x_0}(x-x_0)\to f'(x_0)\cdot0=0,$$
故 $\lim_{x\to x_0}f(x)=f(x_0)$。

几何意义：切线斜率；物理意义：瞬时变化率（速度、加速度、电流、应力率）。

特例：$|x|$ 在 $0$ 连续不可导；$f'(x_0)$ 存在 $\iff$ 左右导数相等。

典型应用：$x(t)$ 的导数 $v(t)$，$v(t)$ 的导数 $a(t)$。

**4.2 求导四则运算与链式法则**

公式：
$$(uv)'=u'v+uv',\qquad\Big(\frac uv\Big)'=\frac{u'v-uv'}{v^2},\qquad \big(f(g(x))\big)'=f'(g(x))g'(x).$$

推导（乘积）：$\Delta(uv)=u(x+\Delta x)v(x+\Delta x)-u(x)v(x)$。加減项 $u(x+\Delta x)v(x)$：
$$\frac{\Delta(uv)}{\Delta x}=v(x)\frac{\Delta u}{\Delta x}+u(x+\Delta x)\frac{\Delta v}{\Delta x}.$$
令 $\Delta x\to0$，用 $u$ 连续 ($u(x+\Delta x)\to u(x)$) 得结论。

推导（链式）：令 $\Delta y=f(g(x+\Delta x))-f(g(x))=f'(g(x))\Delta g+o(\Delta g)$，两边除以 $\Delta x$ 令 $\Delta x\to0$ 得 $dy/dx=f'(g)g'$。

特例：$(\sin x)'=\cos x$ 由和差化积与 $\sin x/x\to1$。

典型应用：复合运动 $x(\theta(t))$ 的速度。

**4.3 微分、线性主部与一阶近似**

定义：$f$ 在 $x$ 可微 $\iff \Delta f=A\Delta x+o(\Delta x)$，且 $A=f'(x)$，记 $df=f'(x)dx$。

推导（可导与可微等价）：若可导，$\Delta f=\frac{\Delta f}{\Delta x}\Delta x$，令 $\frac{\Delta f}{\Delta x}=f'(x)+\alpha$，$\alpha\to0$，则 $\Delta f=f'(x)\Delta x+\alpha\Delta x$。反之由 $o(\Delta x)/\Delta x\to0$ 得差商趋于 $A$。

几何意义：微分是切线增量，$f(x+\Delta x)\approx f(x)+f'(x)\Delta x$。

特例：$\sqrt{1+x}\approx1+x/2$，$\sin x\approx x$。

典型应用：误差传播 $\delta f\approx\sum|\partial f/\partial x_i|\delta x_i$。

**4.4 高阶导数与 Leibniz 公式**

定义：$f^{(n)}=(f^{(n-1)})'$。

Leibniz 公式：
$$(uv)^{(n)}=\sum_{k=0}^n\binom nk u^{(k)}v^{(n-k)}.$$

推导（归纳）：$n=0$ 显然。设对 $n$ 成立，求导一次：
$$(uv)^{(n+1)}=\sum_{k=0}^n\binom nk\big(u^{(k+1)}v^{(n-k)}+u^{(k)}v^{(n-k+1)}\big),$$
重组并用 $\binom nk+\binom n{k-1}=\binom{n+1}k$ 得。

特例：$(x^2e^x)^{(n)}=e^x\big[x^2+2nx+n(n-1)\big]$。

典型应用：$n$ 阶导数的递推格式、微分方程降阶。

**4.5 隐函数与参数方程求导**

隐函数定理（一元）：$F(x,y)=0$ 在 $(x_0,y_0)$ 附近确定 $y=y(x)$，若 $F_y\ne0$，则
$$\frac{dy}{dx}=-\frac{F_x}{F_y}.$$
推导：对 $F(x,y(x))=0$ 关于 $x$ 求全导：$F_x+F_y y'=0$，解出。

参数方程：$\dfrac{dy}{dx}=\dfrac{y_t}{x_t}$；二阶
$$\frac{d^2y}{dx^2}=\frac{1}{x_t}\frac{d}{dt}\Big(\frac{y_t}{x_t}\Big)=\frac{y_{tt}x_t-y_tx_{tt}}{x_t^3}.$$

特例：椭圆 $x=a\cos t,y=b\sin t$，$dy/dx=-(b/a)\tan t$。

典型应用：由约束 $g(x,y)=0$ 表示的机构运动分析。

**4.6 相关变化率与运动学**

定义：若 $x,y$ 均为 $t$ 的函数且满足约束，求 $\dot y$ 由 $\dot x$ 表示。

推导（链式）：对约束求导，$F_x\dot x+F_y\dot y=0\Rightarrow\dot y=-F_x\dot x/F_y$。

特例：圆面积 $A=\pi r^2$，$\dot A=2\pi r\dot r$；体积 $V=\frac43\pi r^3$，$\dot V=4\pi r^2\dot r$。

典型应用：气球充气、雷达测速、机械连杆。

### 第5章 微分中值定理与 Taylor 展开

**5.1 Fermat 引理**

定理：若 $f$ 在 $x_0$ 取局部极值且可导，则 $f'(x_0)=0$。

推导：设极大，当 $x>x_0$ 时差商 $\le0$，令 $x\downarrow x_0$ 得 $f'(x_0)\le0$；当 $x<x_0$ 时差商 $\ge0$，令 $x\uparrow x_0$ 得 $f'(x_0)\ge0$。$\therefore f'(x_0)=0$。

几何意义：可导极值点切线水平。

典型应用：变分法中泛函极值的一阶必要条件。

**5.2 Rolle 定理**

定理：$f\in C[a,b]$，$(a,b)$ 可导，$f(a)=f(b)$，则 $\exists\xi\in(a,b):f'(\xi)=0$。

推导：由最值定理，$f$ 取最大 $M$、最小 $m$。若 $M=m$，$f$ 常值，任意点导数为 0。否则至少一个最值点在内部 $\xi$，由 Fermat 引理 $f'(\xi)=0$。

典型应用：方程根的计数、证明唯一性。

**5.3 Lagrange 中值定理**

定理：$f\in C[a,b]$，$(a,b)$ 可导，则 $\exists\xi$：$f(b)-f(a)=f'(\xi)(b-a)$。

推导：构造辅助函数
$$F(x)=f(x)-f(a)-\frac{f(b)-f(a)}{b-a}(x-a).$$
则 $F(a)=F(b)=0$，由 Rolle 定理 $\exists\xi:F'(\xi)=0$。而 $F'(x)=f'(x)-\frac{f(b)-f(a)}{b-a}$，即得。

几何意义：曲线上存在一点切线平行于端点弦。

特例：$|f'|\le M\Rightarrow|f(b)-f(a)|\le M|b-a|$（Lipschitz 估计）。

典型应用：证明不等式、误差估计、解唯一性。

**5.4 Cauchy 中值定理与 L'Hospital 法则**

Cauchy 定理：$f,g\in C[a,b]$，$(a,b)$ 可导，$g'\ne0$，则 $\exists\xi$：$\dfrac{f(b)-f(a)}{g(b)-g(a)}=\dfrac{f'(\xi)}{g'(\xi)}$。

推导：令 $F(x)=f(x)-f(a)-\frac{f(b)-f(a)}{g(b)-g(a)}(g(x)-g(a))$，用 Rolle。

L'Hospital（$0/0$）：若 $f,g\to0$，在 $x_0$ 去心邻域可导 $g'\ne0$，且 $f'/g'\to L$，则 $f/g\to L$。推导由 Cauchy 定理在 $[x_0,x]$ 上用，令 $x\to x_0$ 得。

特例：$\lim_{x\to0}\frac{\sin x-x}{x^3}=-\frac16$。

典型应用：求极限、分析数值格式精度。

**5.5 Taylor 公式（Peano 余项）**

定理：$f$ 在 $x_0$ 有 $n$ 阶导数，则
$$f(x)=\sum_{k=0}^n\frac{f^{(k)}(x_0)}{k!}(x-x_0)^k+o\big((x-x_0)^n\big).$$

推导：令 $P_n$ 为 Taylor 多项式，$R_n=f-P_n$。由 $R_n^{(j)}(x_0)=0,\ j=0,\dots,n$。对 $R_n(x)/ (x-x_0)^n$ 反复用 L'Hospital，或由导数定义逐步得 $R_n=o((x-x_0)^n)$。

几何意义：用 $n$ 次多项式在 $x_0$ 局部最佳逼近函数。

特例：$e^x=1+x+x^2/2+o(x^2)$。

典型应用：小参数展开、摄动法。

**5.6 Taylor 公式（Lagrange 余项）**

定理：$f$ 在 $[a,b]$ 有 $n+1$ 阶导数，则
$$f(x)=\sum_{k=0}^n\frac{f^{(k)}(x_0)}{k!}(x-x_0)^k+\frac{f^{(n+1)}(\xi)}{(n+1)!}(x-x_0)^{n+1},\quad\xi\text{ 介于 }x_0,x.$$

推导：设 $R_n(x)=f(x)-P_n(x)$。构造
$$\Phi(t)=f(x)-\sum_{k=0}^n\frac{f^{(k)}(t)}{k!}(x-t)^k,$$
则 $\Phi(x_0)=R_n(x)$，$\Phi(x)=0$。对 $\Phi$ 用 Rolle 型论证（或构造 $\Psi(t)=\Phi(t)-(x-t)^{n+1}/(x-x_0)^{n+1}\Phi(x_0)$），得 $\Phi'(t)$ 在某 $\xi$ 为零，计算 $\Phi'(t)=-\frac{f^{(n+1)}(t)}{n!}(x-t)^n$，解得余项。

典型应用：收敛性证明、误差上界 $|R_n|\le\frac{M}{(n+1)!}|x-x_0|^{n+1}$。

**5.7 Taylor 展开的物理近似**

公式：$\sin x=x-\frac{x^3}{6}+O(x^5)$，$\cos x=1-\frac{x^2}{2}+\frac{x^4}{24}+O(x^6)$，$(1+x)^\alpha=1+\alpha x+\frac{\alpha(\alpha-1)}2x^2+O(x^3)$。

推导（$\sin$）：由 $f^{(k)}(0)$ 周期性取值 $0,1,0,-1,\dots$，代入 Taylor 公式。

物理意义：单摆周期 $T=2\pi\sqrt{L/g}\big(1+\frac1{16}\theta_0^2+\cdots\big)$，含非线性修正。

特例：相对论能量 $E=mc^2/\sqrt{1-v^2/c^2}\approx mc^2+\frac12mv^2+\cdots$。

典型应用：有限元刚度矩阵、几何非线性分析。

### 第6章 导数的应用：极值、凸性与曲率

**6.1 单调性与极值判别**

定理：$f'$ 在区间上 $>0$（$<0$）则 $f$ 递增（递减）。

推导：对 $x_1<x_2$ 用 Lagrange 中值定理，$f(x_2)-f(x_1)=f'(\xi)(x_2-x_1)$，符号由 $f'$ 决定。

极值判别：若 $f'(x_0)=0$ 且 $f'$ 在 $x_0$ 两侧由正变负，则极大；反向则极小。二阶判别：$f'(x_0)=0,\ f''(x_0)<0$ 极大，$f''(x_0)>0$ 极小。

推导（二阶）：Taylor $f(x_0+h)=f(x_0)+\frac12f''(x_0)h^2+o(h^2)$，符号由 $f''$ 决定。

典型应用：最优化、能量极小状态判定。

**6.2 凸函数与 Jensen 不等式**

定义：$f$ 凸 $\iff f(\lambda x+(1-\lambda)y)\le\lambda f(x)+(1-\lambda)f(y)$。

判别：$f$ 二阶可导且 $f''\ge0$ 则凸。

推导（一阶条件）：凸 $\iff f(y)\ge f(x)+f'(x)(y-x)$。由割线斜率单调，令 $y\to x$ 得证。

Jensen 不等式：凸函数 $f$ 下，$f\big(\sum\lambda_i x_i\big)\le\sum\lambda_i f(x_i)$，$\lambda_i\ge0,\sum\lambda_i=1$。用归纳与凸定义证明。

典型应用：熵不等式、屈曲能量、凸优化。

**6.3 曲率与挠率（平面曲线）**

定义：曲率 $\kappa=\big|\frac{d\theta}{ds}\big|$，$\theta$ 为切向角，$s$ 为弧长。

公式：$y=f(x)$ 时
$$\kappa=\frac{|y''|}{\big(1+y'^2\big)^{3/2}}.$$
推导：$\tan\theta=y'$，$\sec^2\theta\,d\theta=y''dx$，$ds=\sqrt{1+y'^2}dx$，故 $\frac{d\theta}{ds}=\frac{y''}{(1+y'^2)^{3/2}}$。

参数式：$\kappa=\dfrac{|x'y''-y'x''|}{(x'^2+y'^2)^{3/2}}$。

典型应用：梁弯曲曲率与弯矩 $M=EI\kappa$，力学本构的几何来源。

**6.4 渐近线与函数作图**

定义：$y=kx+b$ 为斜渐近线，若 $f(x)-(kx+b)\to0$。

公式：$k=\lim_{x\to\infty}f(x)/x$，$b=\lim_{x\to\infty}[f(x)-kx]$。

推导：由 $\frac{f(x)-kx-b}{x}\to0$ 得 $k=\lim f/x$，再代入原式得 $b$。

典型应用：绘制相图、稳定性分析中的 nullcline。

**6.5 用导数证明不等式**

方法：构造 $F(x)=$ 左 $-$ 右，证 $F\ge0$。

推导示例：证 $\ln(1+x)<x\ (x>0)$。令 $F(x)=x-\ln(1+x)$，$F(0)=0$，$F'(x)=1-\frac1{1+x}=\frac{x}{1+x}>0$，故 $F>0$。

特例：$e^x\ge1+x$，$\sin x<x\ (x>0)$。

典型应用：误差界、收敛性判据。

### 第7章 不定积分

**7.1 原函数与不定积分**

定义：$F'=f$ 则 $F$ 为 $f$ 的原函数，$\int f\,dx=F+C$。

定理：连续函数有原函数；同一函数的任意两原函数相差常数。

推导：若 $F_1'=F_2'=f$，则 $(F_1-F_2)'=0$，由中值定理 $F_1-F_2\equiv C$。

典型应用：由速度求位移公式。

**7.2 换元积分法**

公式：$\int f(g(x))g'(x)\,dx=\int f(u)\,du,\quad u=g(x)$。

推导：设 $\int f(u)du=F(u)$，则 $F(g(x))'=f(g(x))g'(x)$，故 $\int f(g(x))g'(x)dx=F(g(x))+C$。

特例：$\int \frac{2x}{1+x^2}dx=\ln(1+x^2)+C$。

典型应用：力学换元（如 $v\,dv=a\,dx$）。

**7.3 分部积分法**

公式：$\int u\,dv=uv-\int v\,du$。

推导：$(uv)'=u'v+uv'$ 积分得 $uv=\int u'v+\int uv'$，移项。

特例：$\int x e^x dx=xe^x-e^x+C$，$\int\ln x\,dx=x\ln x-x+C$。

典型应用：功-能积分、拉普拉斯变换。

**7.4 有理函数积分与部分分式**

定义：$\int \frac{P(x)}{Q(x)}dx$，$P,Q$ 多项式。

推导（部分分式）：将 $Q$ 分解为线性与二次不可约因子，设
$$\frac{P(x)}{Q(x)}=\sum_i\frac{A_i}{x-\alpha_i}+\sum_j\frac{B_jx+C_j}{x^2+px+q},$$
比较系数求 $A_i,B_j,C_j$。再逐项积分：$\int\frac{dx}{x-a}=\ln|x-a|$；$\int\frac{2x+p}{x^2+px+q}=\ln(x^2+px+q)$；$\int\frac{dx}{x^2+px+q}$ 配方为 $\arctan$。

典型应用：传递函数部分分式分解。

**7.5 三角有理函数积分**

方法：万能代换 $t=\tan(x/2)$，则
$$\sin x=\frac{2t}{1+t^2},\quad\cos x=\frac{1-t^2}{1+t^2},\quad dx=\frac{2\,dt}{1+t^2}.$$
推导：由 $\cos x=\frac{1-\tan^2(x/2)}{1+\tan^2(x/2)}$ 等半角公式得。

典型应用：含 $\sin,\cos$ 的积分，如振荡响应。

**7.6 简单无理函数积分**

类型：$\int R(x,\sqrt{ax^2+bx+c})dx$，用三角代换或 Euler 代换。

推导示例：$\int\sqrt{a^2-x^2}dx$ 令 $x=a\sin t$，$dx=a\cos t\,dt$，
$$\int a\cos t\cdot a\cos t\,dt=a^2\int\frac{1+\cos2t}{2}dt=\frac{a^2}{2}t+\frac{a^2}{4}\sin2t+C.$$
回代 $t=\arcsin(x/a)$ 得。

典型应用：圆弧、圆面积、摆线。

### 第8章 定积分与微积分基本定理

**8.1 Riemann 积分定义**

定义：分割 $a=x_0<\cdots<x_n=b$，$\Delta x_i=x_i-x_{i-1}$，$\lambda=\max\Delta x_i$，取 $\xi_i\in[x_{i-1},x_i]$，积分和
$$S=\sum_{i=1}^n f(\xi_i)\Delta x_i.$$
若 $\lambda\to0$ 时极限存在且与分割、取点无关，记 $I=\int_a^b f(x)dx$。

几何意义：曲边梯形面积的代数值。

特例：$f$ 连续或单调有界则可积。

典型应用：由分布力求合力、由密度求质量、由流率求总量。

**8.2 Darboux 上、下和与可积条件**

定义：$M_i=\sup_{[x_{i-1},x_i]}f$，$m_i=\inf f$，上、下和 $U=\sum M_i\Delta x_i$，$L=\sum m_i\Delta x_i$。

定理：$f$ 可积 $\iff$ $\forall\varepsilon>0$，存在分割使 $U-L<\varepsilon$。

推导：由 $L\le S\le U$ 及确界性质，可积当且仅当上下和之差可任意小。

特例：Dirichlet 函数不可积（$U-L\equiv1$）。

典型应用：可积性判断、误差控制。

**8.3 变上限积分与微积分基本定理**

定义：$\Phi(x)=\int_a^x f(t)dt$。

定理（第一基本定理）：$f$ 连续 $\Rightarrow\Phi'(x)=f(x)$。

推导：$\Phi(x+\Delta x)-\Phi(x)=\int_x^{x+\Delta x}f(t)dt$。由积分中值定理 $\exists\xi\in[x,x+\Delta x]$：
$$\frac{\Phi(x+\Delta x)-\Phi(x)}{\Delta x}=f(\xi).$$
令 $\Delta x\to0$，$\xi\to x$，由连续性得 $\Phi'(x)=f(x)$。

几何意义：面积函数的变化率等于当前高度。

典型应用：引入原函数，连接微分与积分。

**8.4 Newton–Leibniz 公式**

定理：$f\in C[a,b]$，$F'=f$，则 $\int_a^b f(x)dx=F(b)-F(a)$。

推导：由 $\Phi'(x)=f(x)$ 及 $\Phi(a)=0$，$\Phi$ 与 $F$ 差常数：$F(x)=\Phi(x)+C$。于是 $F(b)-F(a)=\Phi(b)-\Phi(a)=\int_a^b f$。

特例：$\int_0^\pi\sin x\,dx=2$。

典型应用：一切可求原函数的定积分计算。

**8.5 定积分的换元与分部**

换元：$\int_a^b f(g(x))g'(x)dx=\int_{g(a)}^{g(b)}f(u)du$。

分部：$\int_a^b u\,dv=[uv]_a^b-\int_a^b v\,du$。

推导（换元）：由链式法则与 N–L 公式，$(F(g(x)))'=f(g(x))g'(x)$，积分得。

特例：$\int_0^1 xe^x dx=1$。

典型应用：积分变换、能量积分。

**8.6 积分中值定理与平均值**

定理：$f\in C[a,b]$，则 $\exists\xi:\int_a^b f=f(\xi)(b-a)$。

推导：由最值 $m(b-a)\le\int f\le M(b-a)$，除以 $(b-a)$ 得 $m\le$ 平均值 $\le M$，由介值定理取到。

推广（带权）：$\int_a^b f(x)g(x)dx=f(\xi)\int_a^b g(x)dx$，$g$ 不变号。

典型应用：等效集中力、均方根、信号平均。

### 第9章 定积分的几何与力学应用

**9.1 平面面积与曲线弧长**

面积：$A=\int_a^b|f(x)-g(x)|dx$。

弧长：$s=\int_a^b\sqrt{1+f'(x)^2}\,dx$。

推导（弧长）：分割后每段近似直线长 $\sqrt{\Delta x_i^2+\Delta y_i^2}=\sqrt{1+(f'(\xi_i))^2}\Delta x_i$，求和取极限。

参数式：$s=\int_\alpha^\beta\sqrt{x'(t)^2+y'(t)^2}\,dt$。

典型应用：轨迹长度、传动带长度。

**9.2 旋转体体积与旋转面面积**

体积（圆盘法）：$V=\int_a^b\pi f(x)^2dx$。

推导：每个薄片近似圆柱 $\Delta V=\pi f(\xi_i)^2\Delta x_i$。

壳层法：$V=\int 2\pi x f(x)dx$。

旋转面面积：$S=\int_a^b2\pi|f(x)|\sqrt{1+f'^2}dx$。

特例：球体 $V=\frac43\pi R^3$，球面 $S=4\pi R^2$。

**9.3 质心、转动惯量与压力**

质心：$\bar x=\dfrac{\int_a^b x\rho(x)dx}{\int_a^b\rho(x)dx}$。

推导：力矩 $\int x\,dm$ 除以总质量，$dm=\rho dx$。

转动惯量：$I=\int r^2dm$。薄板绕 $y$ 轴 $I_y=\int x^2\rho f(x)dx$。

静水压力：$P=\int\rho g\,h(x)w(x)dx$。

典型应用：梁截面惯性矩、闸门受力。

**9.4 变力做功、能量与引力**

做功：$W=\int_a^b F(x)dx$。

推导：$dW=F(x)dx$，沿路径积分。

弹性能：$W=\int_0^\delta kx\,dx=\frac12k\delta^2$。

引力势能：由 $F=-GMm/r^2$，
$$U(r)-U(\infty)=-\int_\infty^r\frac{GMm}{s^2}ds=-\frac{GMm}{r}.$$

典型应用：能量法解静力问题、逃逸速度 $v_e=\sqrt{2GM/R}$。

**9.5 功能定理与连续体能量**

功-能定理：$W=\Delta T=\frac12mv_2^2-\frac12mv_1^2$。

推导：$W=\int F\,dx=\int m\frac{dv}{dt}dx=\int mv\,dv=\frac12mv^2\big|_0^t$。

连续体动能：$T=\int_V\frac12\rho v^2dV$。

典型应用：能量守恒、冲击问题。

### 第10章 反常积分与特殊函数

**10.1 无穷限反常积分**

定义：$\int_a^{+\infty}f=\lim_{b\to\infty}\int_a^b f$。

$p$-积分判据：$\int_1^\infty\frac{dx}{x^p}$ 收敛 $\iff p>1$。

推导：
$$\int_1^b x^{-p}dx=\begin{cases}\frac{b^{1-p}-1}{1-p},&p\ne1\\\ln b,&p=1\end{cases}$$
令 $b\to\infty$：$p>1$ 时 $b^{1-p}\to0$，收敛于 $1/(p-1)$；$p\le1$ 发散。

典型应用：幂律分布、引力势尾项。

**10.2 无界函数反常积分**

定义：$\int_a^b f$，$f$ 在 $b$ 附近无界，$\int_a^b f=\lim_{\eta\to0^+}\int_a^{b-\eta}f$。

判据：$\int_0^1\frac{dx}{x^p}$ 收敛 $\iff p<1$。

推导：$\int_\eta^1x^{-p}dx=\frac{1-\eta^{1-p}}{1-p}$（$p\ne1$），$\eta^{1-p}\to0$ 当 $p<1$。

典型应用：边裂纹应力奇异场 $r^{-1/2}$。

**10.3 比较判别法与绝对收敛**

比较法：若 $0\le f\le g$，$\int g$ 收敛则 $\int f$ 收敛；$\int f$ 发散则 $\int g$ 发散。

极限形式：若 $f/g\to L\in(0,\infty)$，则同敛散。

绝对收敛：$\int|f|<\infty\Rightarrow\int f$ 收敛。

推导：由 $0\le f^+\le|f|,0\le f^-\le|f|$ 比较得。

典型应用：Fourier 积分的收敛性。

**10.4 Gamma 与 Beta 函数**

定义：$\Gamma(s)=\int_0^\infty x^{s-1}e^{-x}dx\ (s>0)$。

递推：$\Gamma(s+1)=s\Gamma(s)$。

推导：分部积分
$$\Gamma(s+1)=\int_0^\infty x^se^{-x}dx=\big[-x^se^{-x}\big]_0^\infty+s\int_0^\infty x^{s-1}e^{-x}dx=s\Gamma(s).$$
特例：$\Gamma(1)=1\Rightarrow\Gamma(n+1)=n!$，$\Gamma(1/2)=\sqrt\pi$。

Beta 函数：$B(p,q)=\int_0^1t^{p-1}(1-t)^{q-1}dt=\frac{\Gamma(p)\Gamma(q)}{\Gamma(p+q)}$。

典型应用：统计分布、分数阶微积分。

### 第11章 数项级数

**11.1 收敛定义与 Cauchy 准则**

定义：$\sum a_n$ 收敛 $\iff S_n=\sum_{k=1}^na_k$ 收敛。

必要条件：收敛 $\Rightarrow a_n\to0$。

Cauchy 准则：$\sum a_n$ 收敛 $\iff\forall\varepsilon>0,\exists N,\forall n>m>N:\big|\sum_{k=m+1}^na_k\big|<\varepsilon$。

典型应用：错误格式的发散性诊断。

**11.2 正项级数判别法**

比较、比值、根值、积分判别：

比值：$\lim\frac{a_{n+1}}{a_n}=r<1$ 收敛，$r>1$ 发散。

推导：$a_{n+1}\le q a_n\ (q<1)$，则 $a_n\le a_N q^{n-N}$，与几何级数比较。

根值：$\lim a_n^{1/n}=r<1$ 收敛。

积分判别：$f$ 正、减，$\sum f(n)$ 与 $\int_1^\infty f(x)dx$ 同敛散。推导由
$$\int_n^{n+1}f(x)dx\le f(n)\le\int_{n-1}^nf(x)dx.$$

$p$-级数：$\sum1/n^p$ 收敛 $\iff p>1$。

典型应用：确定迭代法收敛阶。

**11.3 交错级数与 Leibniz 判别**

定理：$a_n\downarrow0$，则 $\sum(-1)^{n-1}a_n$ 收敛，余项 $|R_n|\le a_{n+1}$。

推导：部分和 $S_{2m}$ 递增有上界 $a_1$，$S_{2m+1}$ 递减有下界，且差 $S_{2m+1}-S_{2m}=a_{2m+1}\to0$，故两子列同极限。

特例：$\sum(-1)^{n-1}/n=\ln2$。

**11.4 绝对收敛与条件收敛**

定义：$\sum|a_n|$ 收敛为绝对收敛；$\sum a_n$ 收敛而 $\sum|a_n|$ 发散为条件收敛。

定理：绝对收敛 $\Rightarrow$ 收敛；绝对收敛级数可任意重排。

推导：由 Cauchy 准则 $|\sum a_k|\le\sum|a_k|$。

Riemann 重排定理：条件收敛级数可重排为任意和。

典型应用：Fourier 级数的绝对/条件收敛。

**11.5 级数与积分的联系、数值估计**

方法：用积分估计级数余项。

推导：$\sum_{k=n+1}^\infty f(k)\le\int_n^\infty f(x)dx\le\sum_{k=n}^\infty f(k)$。

特例：$\sum1/n^2=\pi^2/6$。

典型应用：截断误差控制。

### 第12章 函数项级数与一致收敛

**12.1 逐点收敛与一致收敛**

定义：$f_n\to f$ 逐点：固定 $x$，$f_n(x)\to f(x)$。一致：$\sup_{x\in I}|f_n(x)-f(x)|\to0$。

推导（一致收敛的 Cauchy 准则）：$f_n$ 一致收敛 $\iff\forall\varepsilon,\exists N,\forall n,m>N,\forall x:|f_n(x)-f_m(x)|<\varepsilon$。

特例：$f_n(x)=x^n$ 在 $[0,1)$ 逐点收敛不连续极限，不一致。

典型应用：保证极限运算与积分、导数交换。

**12.2 Weierstrass M-判别法**

定理：$|f_n(x)|\le M_n$ 且 $\sum M_n<\infty$，则 $\sum f_n$ 一致收敛。

推导：由 Cauchy 准则，$\big|\sum_{k=m+1}^nf_k(x)\big|\le\sum_{k=m+1}^nM_k<\varepsilon$，与 $x$ 无关。

典型应用：幂级数、Fourier 级数的一致收敛。

**12.3 一致收敛级数的分析性质**

定理：$f_n\in C(I)$ 且 $f_n\rightrightarrows f$，则 $f\in C(I)$；$\int_a^b f=\lim\int_a^b f_n$；若 $f_n'$ 一致收敛且 $f_n$ 逐点收敛，则 $f'=\lim f_n'$。

推导（连续）：$|f(x)-f(x_0)|\le|f(x)-f_n(x)|+|f_n(x)-f_n(x_0)|+|f_n(x_0)-f(x_0)|$，选 $n$ 大使前两项小。

典型应用：幂级数逐项求导、积分。

**12.4 一致收敛与极限交换、Dini 定理**

Dini 定理：$f_n$ 连续、$f_n\to f$ 连续、$f_n$ 单调，则在紧集上一致收敛。

推导：利用 $\sup|f_n-f|$ 定义反证，若不一致则存在 $x_k$ 子列使差 $\ge\varepsilon_0$，取收敛子列与极限连续矛盾。

典型应用：单调逼近格式。

### 第13章 幂级数与 Taylor 级数

**13.1 收敛半径与 Abel 定理**

定义：$\sum a_n(x-x_0)^n$，收敛半径 $R$ 使 $|x-x_0|<R$ 收敛，$>R$ 发散。

公式：$R=\big(\limsup|a_n|^{1/n}\big)^{-1}$ 或 $R=\lim|a_n/a_{n+1}|$。

推导（根值）：令 $L=\limsup|a_n|^{1/n}$，则 $|a_n(x-x_0)^n|^{1/n}\to L|x-x_0|$，由根值判别：$L|x-x_0|<1$ 收敛，得 $R=1/L$。

Abel 定理：在收敛半径内绝对收敛，在 $|x-x_0|<r<R$ 一致收敛。

典型应用：解线性常微分方程的级数解。

**13.2 幂级数运算与逐项微积分**

定理：$|x-x_0|<R$ 内，$f(x)=\sum a_n(x-x_0)^n$ 可逐项求导、积分，且收敛半径不变。

公式：$f'(x)=\sum n a_n(x-x_0)^{n-1}$，$\int f=\sum\frac{a_n}{n+1}(x-x_0)^{n+1}$。

推导：由 Abel 定理在 $|x-x_0|\le r<R$ 一致收敛及第 12 章性质。

典型应用：求导积分法求和函数，如 $\sum nx^n=x/(1-x)^2$。

**13.3 初等函数的幂级数展开**

公式：
$$e^x=\sum_{n=0}^\infty\frac{x^n}{n!},\quad \ln(1+x)=\sum_{n=1}^\infty\frac{(-1)^{n-1}x^n}{n}\ (|x|<1),$$
$$\frac1{1-x}=\sum_{n=0}^\infty x^n,\quad \sin x=\sum_{n=0}^\infty\frac{(-1)^n x^{2n+1}}{(2n+1)!}.$$

推导（几何级数）：由 $(1-x)\sum_{k=0}^nx^k=1-x^{n+1}$，令 $n\to\infty$ 当 $|x|<1$。

典型应用：数值计算、微扰展开。

**13.4 Taylor 级数收敛性与余项**

定理：$f$ 在 $x_0$ 邻域无穷可导，$f(x)=\sum\frac{f^{(n)}(x_0)}{n!}(x-x_0)^n\iff R_n\to0$。

推导：由 Taylor 公式 $f=P_n+R_n$，级数收敛当且仅当余项趋于 0。

特例：$f(x)=e^{-1/x^2}$ 在 0 处 Taylor 级数为 0 但不等于 $f$。

典型应用：判断展开的可用性。

**13.5 幂级数求解微分方程**

方法：设 $y=\sum a_nx^n$，代入方程比较系数。

推导示例：$y'-y=0$，则 $\sum na_nx^{n-1}-\sum a_nx^n=0$，即 $\sum(n+1)a_{n+1}x^n=\sum a_nx^n$，得 $(n+1)a_{n+1}=a_n\Rightarrow a_n=a_0/n!$，故 $y=a_0e^x$。

典型应用：Airy 方程、Bessel 方程、量子力学级数解。

### 第14章 Fourier 级数与振动分析

**14.1 三角级数与正交性**

定义：$f(x)\sim\frac{a_0}{2}+\sum_{n=1}^\infty(a_n\cos nx+b_n\sin nx)$。

正交性：
$$\int_{-\pi}^\pi\cos mx\cos nx\,dx=\pi\delta_{mn},\quad \int_{-\pi}^\pi\sin mx\sin nx\,dx=\pi\delta_{mn},\quad \int_{-\pi}^\pi\cos mx\sin nx\,dx=0.$$
推导：用积化和差，如 $\cos mx\cos nx=\frac12[\cos(m+n)x+\cos(m-n)x]$，$m\ne n$ 时积分为 0。

典型应用：振动模态的正交性、能量解耦。

**14.2 Fourier 系数公式**

公式：
$$a_n=\frac1\pi\int_{-\pi}^\pi f(x)\cos nx\,dx,\quad b_n=\frac1\pi\int_{-\pi}^\pi f(x)\sin nx\,dx.$$
推导：对级数两边乘 $\cos nx$ 并逐项积分，由正交性只剩 $n$ 项，得 $a_n$，同法得 $b_n$。

典型应用：频谱分析、信号分解。

**14.3 收敛定理与 Dirichlet 条件**

定理：$f$ 在 $[-\pi,\pi]$ 分段光滑，则 Fourier 级数在连续点收敛于 $f(x)$，在跳跃点收敛于 $\frac{f(x^-)+f(x^+)}2$。

推导：将部分和写成 Dirichlet 积分 $S_N(x)=\frac1{2\pi}\int_{-\pi}^\pi f(x+t)D_N(t)dt$，利用核的性质与 Riemann-Lebesgue 引理得。

典型应用：方波展开 $\frac4\pi\sum_{k=1}^\infty\frac{\sin(2k-1)x}{2k-1}$。

**14.4 Parseval 等式与能量**

公式：
$$\frac1\pi\int_{-\pi}^\pi f^2dx=\frac{a_0^2}{2}+\sum_{n=1}^\infty(a_n^2+b_n^2).$$
推导：对 $f^2$ 积分并代入级数，由正交性交叉项消失。

物理意义：时域能量等于各频率分量能量之和。

典型应用：信号功率谱、振动能量分布。

**14.5 半区间展开与复数形式**

奇/偶延拓得正弦/余弦级数。

复数形式：$f(x)=\sum_{n=-\infty}^\infty c_ne^{inx}$，$c_n=\frac1{2\pi}\int_{-\pi}^\pi f(x)e^{-inx}dx$。

推导：由 Euler 公式 $e^{inx}=\cos nx+i\sin nx$，合并 $a_n,b_n$ 得 $c_n=(a_n-ib_n)/2$ 等。

典型应用：Fourier 变换、频响函数。

**14.6 受迫振动与共振**

方程：$m\ddot x+c\dot x+kx=F_0\cos\omega t$。

稳态解：$x_p=A\cos(\omega t-\varphi)$，$A=\frac{F_0}{\sqrt{(k-m\omega^2)^2+(c\omega)^2}}$。

推导：设特解代入，比较 $\cos,\sin$ 系数，解得
$$A=\frac{F_0}{\sqrt{(k-m\omega^2)^2+(c\omega)^2}},\quad\tan\varphi=\frac{c\omega}{k-m\omega^2}.$$
当 $\omega=\sqrt{k/m}$ 且 $c$ 小时共振。

典型应用：结构振动、颤振分析。

### 第15章 多元函数微分学

**15.1 多元极限与连续**

定义：$\lim_{(x,y)\to(x_0,y_0)}f=L$：$\forall\varepsilon,\exists\delta$，$0<\sqrt{(x-x_0)^2+(y-y_0)^2}<\delta\Rightarrow|f-L|<\varepsilon$。

推导（极限不存在）：沿不同路径极限不同则不存在。例 $f=\frac{xy}{x^2+y^2}$，沿 $y=x$ 极限 $1/2$，沿 $y=0$ 极限 0。

典型应用：判断函数光滑性。

**15.2 偏导数与全微分**

定义：$f_x=\lim_{h\to0}\frac{f(x+h,y)-f(x,y)}h$。

全微分：$df=f_xdx+f_ydy$ 且 $\Delta f-f_x\Delta x-f_y\Delta y=o(\sqrt{\Delta x^2+\Delta y^2})$。

推导（可微 $\Rightarrow$ 偏导存在）：沿 $y$ 固定，$\Delta f=f_x\Delta x+o(\Delta x)$，故偏导存在且等于 $f_x$。

特例：$f_x,f_y$ 存在不一定可微；$f_x,f_y$ 连续 $\Rightarrow$ 可微。

典型应用：线性化、误差分析、形变梯度。

**15.3 复合函数链式法则**

公式：$z=f(u,v),\ u=u(x,y),v=v(x,y)$，则
$$\frac{\partial z}{\partial x}=\frac{\partial z}{\partial u}\frac{\partial u}{\partial x}+\frac{\partial z}{\partial v}\frac{\partial v}{\partial x}.$$
推导：由全微分 $dz=z_udu+z_vdv$，而 $du=u_xdx+u_ydy,\ dv=v_xdx+v_ydy$，代入比较 $dx,dy$ 系数。

典型应用：多元运动学、材料坐标变换。

**15.4 隐函数定理与雅可比**

定理：$F(x,y,z)=0$ 确定 $z=z(x,y)$，$F_z\ne0$，则
$$\frac{\partial z}{\partial x}=-\frac{F_x}{F_z},\quad\frac{\partial z}{\partial y}=-\frac{F_y}{F_z}.$$
推导：对 $F(x,y,z(x,y))=0$ 分别对 $x,y$ 求偏导。

雅可比行列式：$\frac{\partial(u,v)}{\partial(x,y)}=\begin{vmatrix}u_x&u_y\\v_x&v_y\end{vmatrix}$，用于变量替换。

典型应用：约束力学、坐标变换。

**15.5 方向导数与梯度**

定义：$\frac{\partial f}{\partial\boldsymbol l}=\nabla f\cdot\boldsymbol l$。

公式：$\nabla f=(f_x,f_y,f_z)$；$\frac{\partial f}{\partial\boldsymbol l}=\nabla f\cdot\boldsymbol l=|\nabla f|\cos\theta$。

推导：沿单位向量 $\boldsymbol l$ 的参数化，用链式法则得。

几何意义：梯度指向最速上升方向，模为最大方向导数。

典型应用：热流 $q=-k\nabla T$，Fourier 定律。

**15.6 多元 Taylor 公式与极值**

Taylor：$f(\boldsymbol x_0+\boldsymbol h)=f(\boldsymbol x_0)+\nabla f\cdot\boldsymbol h+\frac12\boldsymbol h^TH\boldsymbol h+o(|\boldsymbol h|^2)$，$H=(f_{x_ix_j})$。

极值判别：$\nabla f=0$。$H$ 正定 $\Rightarrow$ 极小；负定 $\Rightarrow$ 极大；不定 $\Rightarrow$ 鞍点。

推导：二阶 Taylor，符号由二次型决定。

典型应用：势能极值、稳定性。

**15.7 Lagrange 乘子法**

问题：$\min f(\boldsymbol x)$ s.t. $g_i(\boldsymbol x)=0$。

公式：$\nabla f=\sum_i\lambda_i\nabla g_i$，即 $\nabla(f-\sum\lambda_i g_i)=0$。

推导：约束曲面切空间内任意方向 $\delta\boldsymbol x$ 满足 $\nabla g_i\cdot\delta\boldsymbol x=0$；极值要求 $\nabla f\cdot\delta\boldsymbol x=0$，故 $\nabla f$ 属于 $\nabla g_i$ 张成空间。

典型应用：最小势能原理、能量法、屈曲临界载荷。

### 第16章 重积分、场论与变分力学

**16.1 二重积分与累次积分**

定义：$\iint_D f\,dA=\lim\sum f(\xi_i,\eta_i)\Delta A_i$。

Fubini 定理：$D=\{(x,y):a\le x\le b,\ \varphi_1(x)\le y\le\varphi_2(x)\}$，则
$$\iint_Df\,dA=\int_a^b\!\!\int_{\varphi_1(x)}^{\varphi_2(x)}f(x,y)\,dy\,dx.$$
推导思路：用细网格将区域划分为小块，先沿 $y$ 向累加再沿 $x$ 向，Riemann 和重排。

典型应用：截面法求体积、质量、惯性矩。

**16.2 三重积分与柱、球坐标**

公式：$dV=dxdydz$；柱坐标 $dV=r\,drd\theta dz$；球坐标 $dV=\rho^2\sin\varphi\,d\rho d\varphi d\theta$。

推导（球坐标 Jacobian）：
$$\frac{\partial(x,y,z)}{\partial(\rho,\varphi,\theta)}=\rho^2\sin\varphi,\quad dV=\rho^2\sin\varphi\,d\rho\,d\varphi\,d\theta.$$

典型应用：球对称引力场、转动惯量。

**16.3 曲线积分与 Green 公式**

第一类（弧长）：$\int_L f\,ds$；第二类（坐标）：$\int_L Pdx+Qdy$。

Green 公式：$L=\partial D$ 正向，
$$\oint_L P\,dx+Q\,dy=\iint_D\Big(\frac{\partial Q}{\partial x}-\frac{\partial P}{\partial y}\Big)dA.$$
推导：对 $D$ 作矩形分割，用微积分基本定理：
$$\int_a^b\int_{y_1(x)}^{y_2(x)}\frac{\partial P}{\partial y}dy\,dx=\int_a^b[P(x,y_2)-P(x,y_1)]dx,$$
与边界积分对应，两式相加得证。

典型应用：平面流的环量、面积公式 $A=\frac12\oint x\,dy-y\,dx$。

**16.4 曲面积分、Gauss 与 Stokes 公式**

Gauss（散度）定理：
$$\oiint_{\partial V}\boldsymbol F\cdot\boldsymbol n\,dS=\iiint_V\nabla\cdot\boldsymbol F\,dV.$$
Stokes（旋度）定理：
$$\oint_{\partial S}\boldsymbol F\cdot d\boldsymbol r=\iint_S(\nabla\times\boldsymbol F)\cdot\boldsymbol n\,dS.$$

推导（Gauss 一维情形）：对 $\partial F_x/\partial x$ 项沿 $x$ 积分，用微积分基本定理化为两端面通量，三项相加。

典型应用：守恒律、Navier–Stokes 推导。

**16.5 变分法与 Euler–Lagrange 方程**

问题：$J[y]=\int_{x_1}^{x_2}f(x,y,y')dx$ 取极值，$y(x_1)=y_1,y(x_2)=y_2$。

推导：设 $y_\varepsilon=y+\varepsilon\eta$，$\eta(x_1)=\eta(x_2)=0$。
$$\frac{dJ}{d\varepsilon}\Big|_{\varepsilon=0}=\int_{x_1}^{x_2}\Big(\frac{\partial f}{\partial y}\eta+\frac{\partial f}{\partial y'}\eta'\Big)dx.$$
对第二项分部积分：
$$\int_{x_1}^{x_2}\frac{\partial f}{\partial y'}\eta'\,dx=\Big[\frac{\partial f}{\partial y'}\eta\Big]_{x_1}^{x_2}-\int_{x_1}^{x_2}\frac{d}{dx}\Big(\frac{\partial f}{\partial y'}\Big)\eta\,dx=-\int_{x_1}^{x_2}\frac{d}{dx}\Big(\frac{\partial f}{\partial y'}\Big)\eta\,dx.$$
故
$$\delta J=\int_{x_1}^{x_2}\Big[\frac{\partial f}{\partial y}-\frac{d}{dx}\Big(\frac{\partial f}{\partial y'}\Big)\Big]\eta\,dx=0.$$
由 $\eta$ 任意（变分法基本引理）得 Euler–Lagrange 方程
$$\boxed{\ \frac{\partial f}{\partial y}-\frac{d}{dx}\Big(\frac{\partial f}{\partial y'}\Big)=0\ }$$

特例（最速降线/悬链线/短程线）均可由上式导出。

典型应用：力学中的 Hamilton 原理。

**16.6 力学的变分原理：Lagrange 与 Hamilton**

Lagrange 量：$L=T-V$。作用量 $S=\int_{t_1}^{t_2}L(q,\dot q,t)dt$。Hamilton 原理：真实运动使 $S$ 取驻值，$\delta S=0$。

推导：对每个广义坐标 $q_i$，Euler–Lagrange 方程给出
$$\frac{d}{dt}\Big(\frac{\partial L}{\partial\dot q_i}\Big)-\frac{\partial L}{\partial q_i}=0.$$

牛顿例（一维）：$T=\frac12m\dot x^2$，$V=V(x)$，$L=\frac12m\dot x^2-V(x)$。则 $\partial L/\partial\dot x=m\dot x$，$\partial L/\partial x=-V'(x)$，代入得 $m\ddot x=-V'(x)=F$，即牛顿第二定律。

推广 N 质点：$L=\sum\frac12m_i\dot{\boldsymbol r}_i^2-\sum V$，E–L 方程给出 $m_i\ddot{\boldsymbol r}_i=-\nabla_iV$。

**16.7 广义动量、Hamilton 方程与守恒律**

广义动量：$p_i=\dfrac{\partial L}{\partial\dot q_i}$。

Hamilton 量（Legendre 变换）：$H=\sum_ip_i\dot q_i-L(q,\dot q,t)$。

推导：对 $H$ 取全微分
$$dH=\sum_i\dot q_idp_i+\sum_ip_id\dot q_i-\sum_i\frac{\partial L}{\partial q_i}dq_i-\sum_i\frac{\partial L}{\partial\dot q_i}d\dot q_i-\frac{\partial L}{\partial t}dt.$$
由 $p_i=\partial L/\partial\dot q_i$ 消去 $d\dot q_i$ 项，并用 E–L 方程 $\dot p_i=\partial L/\partial q_i$，得
$$dH=\sum_i\dot q_idp_i-\sum_i\dot p_idq_i-\frac{\partial L}{\partial t}dt,$$
故
$$\dot q_i=\frac{\partial H}{\partial p_i},\qquad \dot p_i=-\frac{\partial H}{\partial q_i},\qquad \frac{\partial H}{\partial t}=-\frac{\partial L}{\partial t}.$$

守恒律：若 $\partial L/\partial q_i=0$，则 $p_i$ 守恒；若 $\partial L/\partial t=0$，则 $H=T+V$ 守恒（能量守恒）。

典型应用：刚体动力学、轨道力学、规范场。

**16.8 连续介质与 Navier–Stokes 方程推导**

思想：对控制体 $V$ 用 Newton 第二定律，再用 Reynolds 输运与 Gauss 定理化为局部方程。

Reynolds 输运：
$$\frac{d}{dt}\int_{V(t)}\rho\boldsymbol v\,dV=\int_V\Big[\frac{\partial(\rho\boldsymbol v)}{\partial t}+\nabla\cdot(\rho\boldsymbol v\otimes\boldsymbol v)\Big]dV.$$
动量守恒：
$$\frac{d}{dt}\int_V\rho\boldsymbol v\,dV=\int_V\rho\boldsymbol f\,dV+\oint_{\partial V}\boldsymbol\sigma\cdot\boldsymbol n\,dS.$$
由 Gauss 定理 $\oint\boldsymbol\sigma\cdot\boldsymbol n\,dS=\int_V\nabla\cdot\boldsymbol\sigma\,dV$，得
$$\frac{\partial(\rho\boldsymbol v)}{\partial t}+\nabla\cdot(\rho\boldsymbol v\otimes\boldsymbol v)=\rho\boldsymbol f+\nabla\cdot\boldsymbol\sigma.$$
不可压牛顿流体本构：
$$\boldsymbol\sigma=-p\boldsymbol I+2\mu\boldsymbol D+\lambda(\nabla\cdot\boldsymbol v)\boldsymbol I,\qquad \boldsymbol D=\frac12\big(\nabla\boldsymbol v+\nabla\boldsymbol v^T\big).$$
对不可压 $\nabla\cdot\boldsymbol v=0$，$\lambda$ 项消失，代入得
$$\rho\Big(\frac{\partial\boldsymbol v}{\partial t}+\boldsymbol v\cdot\nabla\boldsymbol v\Big)=-\nabla p+\mu\nabla^2\boldsymbol v+\rho\boldsymbol f.$$

量纲：$[\rho]=M/L^3$，$[p]=M/(LT^2)$，$[\mu]=M/(LT)$。

特例（静止流体）：$\nabla p=\rho\boldsymbol f$，重力下 $p=p_0+\rho gh$。

**16.9 线弹性本构与能量法**

应变张量：$\boldsymbol\varepsilon=\frac12(\nabla\boldsymbol u+\nabla\boldsymbol u^T)$。

本构（各向同性线弹性）：
$$\boldsymbol\sigma=\lambda\,\mathrm{tr}(\boldsymbol\varepsilon)\boldsymbol I+2\mu\boldsymbol\varepsilon,$$
$$\sigma_{ij}=\frac{E}{1+\nu}\Big(\varepsilon_{ij}+\frac{\nu}{1-2\nu}\varepsilon_{kk}\delta_{ij}\Big).$$
推导：由应变能密度 $W=\frac12\lambda(\varepsilon_{kk})^2+\mu\varepsilon_{ij}\varepsilon_{ij}$，对 $\varepsilon_{ij}$ 求导得 $\sigma_{ij}=\partial W/\partial\varepsilon_{ij}$。

应变能：$U=\int_VW\,dV$。最小势能原理：$\Pi=U-W_{\text{ext}}$ 取极小，等价于平衡方程，可由变分推导。

典型应用：有限元刚度矩阵由 $W$ 双线性化得 $K_{ij}=\int B_i^TDB_jdV$。

**16.10 Lyapunov 稳定性与线性化**

定义：平衡点 $x_e$，Lyapunov 稳定：$\forall\varepsilon>0,\exists\delta>0,\|x(0)-x_e\|<\delta\Rightarrow\|x(t)-x_e\|<\varepsilon,\forall t\ge0$。

Lyapunov 定理：若存在正定 $V(x)$ 使 $\dot V=\nabla V\cdot f(x)\le0$，则稳定；若 $\dot V<0\ (x\ne x_e)$，则渐近稳定。

推导：由 $\dot V\le0$，$V(x(t))\le V(x(0))$。选定 $V$ 的水平集 $\{V\le c\}\subseteq B_\varepsilon$，令 $\delta$ 对应最小水平集，则轨道不越出。

线性化：$\dot x=Ax$，$A=Df(x_e)$。特征值实部全负 $\Rightarrow$ 渐近稳定；存在正实部 $\Rightarrow$ 不稳定。

特例：单摆平衡 $(0,0)$ 线性化矩阵
$$A=\begin{pmatrix}0&1\\-g/L&0\end{pmatrix},\quad \lambda=\pm i\sqrt{g/L},$$
临界情形，非线性项决定稳定性（中心）。

典型应用：结构稳定性、控制镇定、屈曲后行为。

---
⬆ 上级领域：[[数学]]
