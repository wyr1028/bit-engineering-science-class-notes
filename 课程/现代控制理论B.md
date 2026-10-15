---
课程名称: 现代控制理论B
学分: 2
总学时: 32
考核方式: 考试
课程属性: 选修
课程性质: 专业课
开课单位: 空天科学与技术学院
领域: 电子·控制·系统
tags:
  - 课程
  - 领域/电子·控制·系统
---

# 现代控制理论B


### 第1章 动态系统的状态空间描述与线性化

**知识点1.1　状态、状态变量与状态空间的定义**

定义：系统在时刻 $t_0$ 的状态，是指能够完全刻画系统未来演化所需的、与输入历史无关的最小信息集合；这些变量组成状态向量 $x(t)\in\mathbb{R}^n$，其张成的空间称为状态空间。

假设：系统是确定性的、因果的；状态变量数目 $n$ 等于系统独立储能元件的数目（力学系统中为独立广义坐标与广义速度的数目的两倍）。

公式：一般非线性系统写为
$$\dot x=f(x,u,t),\qquad y=g(x,u,t).$$

推导（状态的存在性）：设系统可由 $n$ 阶微分方程描述
$$\frac{d^n y}{dt^n}=\phi\!\left(y,\dot y,\dots,y^{(n-1)},u,t\right).$$
取
$$x_1=y,\ x_2=\dot y,\ \dots,\ x_n=y^{(n-1)},$$
则
$$\dot x_1=x_2,\ \dot x_2=x_3,\ \dots,\ \dot x_{n-1}=x_n,\ \dot x_n=\phi(x_1,\dots,x_n,u,t).$$
写成向量形式即为 $\dot x=f(x,u,t)$，故状态变量集合使 $t_0$ 以后的行为只由 $(x(t_0),u_{[t_0,t]})$ 决定。

物理/几何意义：$x(t)$ 是状态空间中的一条轨线，$f$ 是状态空间上的向量场，轨线处处与向量场相切。

极限/特例：线性时不变系统 $f=Ax+Bu$，$g=Cx+Du$，向量场不随时间变化，轨线由线性叠加决定。

典型应用：卫星姿态、电路网络、机器人关节、飞行器纵向运动的统一描述。

**知识点1.2　由高阶微分方程建立状态空间（相变量法）**

定义：相变量法是把 $n$ 阶标量微分方程化为 $n$ 个一阶方程的标准步骤。

假设：微分方程最高阶项系数可归一；初始条件 $y(0),\dot y(0),\dots,y^{(n-1)}(0)$ 已知。

公式：对
$$y^{(n)}+a_{n-1}y^{(n-1)}+\cdots+a_1\dot y+a_0y=b_0u,$$
取相变量 $x_i=y^{(i-1)}$，得
$$\dot x=Ax+Bu,\quad A=\begin{bmatrix}0&1&0&\cdots&0\\0&0&1&\cdots&0\\\vdots&&&\ddots&\vdots\\0&0&0&\cdots&1\\-a_0&-a_1&-a_2&\cdots&-a_{n-1}\end{bmatrix},\quad B=\begin{bmatrix}0\\0\\\vdots\\0\\b_0\end{bmatrix}.$$

推导：由 $x_i=y^{(i-1)}$ 直接求导 $\dot x_i=y^{(i)}=x_{i+1}$（$i<n$）；对 $i=n$ 用原方程解出 $y^{(n)}-b_0u=-a_0y-\cdots-a_{n-1}y^{(n-1)}$，即 $\dot x_n=-a_0x_1-\cdots-a_{n-1}x_n+b_0u$。

物理/几何意义：$A$ 称为友矩阵（companion matrix），其特征多项式恰为原方程特征多项式 $\lambda^n+a_{n-1}\lambda^{n-1}+\cdots+a_0$。

极限/特例：$a_0=0$ 时系统含积分环节，$A$ 奇异，存在零特征值。

典型应用：单自由度振动 $\ddot y+2\zeta\omega_n\dot y+\omega_n^2y=\omega_n^2u$ 的标准状态化。

**知识点1.3　传递函数与能控标准型的相互转换**

定义：单输入单输出系统传递函数 $G(s)=C(sI-A)^{-1}B+D$；能控标准型是使 $B$ 为坐标单位向量、$A$ 为友矩阵的实现。

假设：传递函数为真分式 $G(s)=\dfrac{b_{n-1}s^{n-1}+\cdots+b_0}{s^n+a_{n-1}s^{n-1}+\cdots+a_0}$。

公式：能控标准型
$$\dot x=\begin{bmatrix}0&1&&\\ \vdots&&\ddots&\\ 0&&&1\\ -a_0&-a_1&\cdots&-a_{n-1}\end{bmatrix}x+\begin{bmatrix}0\\\vdots\\0\\1\end{bmatrix}u,\quad y=[\,b_0\ b_1\ \cdots\ b_{n-1}\,]x.$$

推导：由 $X_1(s)=\dfrac{1}{s^n+a_{n-1}s^{n-1}+\cdots+a_0}U(s)$ 出发，设中间变量 $V$，则 $V=U/(D(s))$，$Y=(N(s))V$。取相变量 $x_1=v,\dots,x_n=v^{(n-1)}$ 得 $A$、$B$；又 $Y= b_0 v+b_1\dot v+\cdots+b_{n-1}v^{(n-1)}$，代入 $x_i$ 得 $C=[b_0\ \cdots\ b_{n-1}]$。

物理/几何意义：能控标准型把传递函数的分子系数直接放在输出矩阵 $C$ 中，分母系数放在 $A$ 的最后一行。

极限/特例：$D\neq0$（$N(s)$ 阶数与 $D(s)$ 相同）时需先做多项式除法 $G(s)=D+N'(s)/D(s)$。

典型应用：用 MATLAB `tf2ss`/`ss2tf` 在传递函数与状态空间之间转换。

**知识点1.4　非线性系统平衡点与雅可比线性化**

定义：平衡点 $x_e$ 满足 $f(x_e,u_e)=0$；雅可比线性化是在平衡点附近用一阶泰勒展开近似原系统。

假设：$f,g$ 在 $(x_e,u_e)$ 附近连续可微。

公式：令 $\delta x=x-x_e,\ \delta u=u-u_e$，则
$$\dot{\delta x}=A\delta x+B\delta u,\quad \delta y=C\delta x+D\delta u,$$
其中
$$A=\left.\frac{\partial f}{\partial x}\right|_{(x_e,u_e)},\quad B=\left.\frac{\partial f}{\partial u}\right|_{(x_e,u_e)},\quad C=\left.\frac{\partial g}{\partial x}\right|_{(x_e,u_e)},\quad D=\left.\frac{\partial g}{\partial u}\right|_{(x_e,u_e)}.$$

推导：写
$$f(x_e+\delta x,u_e+\delta u)=f(x_e,u_e)+\frac{\partial f}{\partial x}\delta x+\frac{\partial f}{\partial u}\delta u+O(\|\delta x\|^2,\|\delta u\|^2).$$
由 $f(x_e,u_e)=0$ 并对时间求导：
$$\dot{\delta x}=\dot x=f(x,u)=A\delta x+B\delta u+O(\|\cdot\|^2),$$
略去高阶项即得线性化模型。

物理/几何意义：线性化用平衡点处切平面上的向量场代替原向量场，刻画小扰动运动。

极限/特例：当某特征值实部为零（临界情形）时，线性化不能判定原非线性系统的稳定性，须用中心流形定理。

典型应用：倒立摆、飞行器配平点小扰动线性模型、电力系统小信号分析。

**知识点1.5　状态变换下的系统等价性**

定义：对非奇异矩阵 $T$，令 $\bar x=Tx$，所得系统 $(\bar A,\bar B,\bar C,\bar D)$ 与原系统代数等价。

假设：$T$ 可逆，且不改变输入输出。

公式：
$$\dot{\bar x}=TAT^{-1}\bar x+TBu,\quad y=CT^{-1}\bar x+Du.$$

推导：由 $x=T^{-1}\bar x$，$\dot x=T^{-1}\dot{\bar x}$，代入 $\dot x=Ax+Bu$：
$$T^{-1}\dot{\bar x}=AT^{-1}\bar x+Bu\ \Rightarrow\ \dot{\bar x}=TAT^{-1}\bar x+TBu,$$
输出 $y=CT^{-1}\bar x+Du$。

物理/几何意义：状态空间基底的更换；特征值、传递函数、能控性、能观性均不变。

极限/特例：$T$ 为模态矩阵时，$\bar A$ 为对角或约当标准型。

典型应用：把系统化为约当型以解耦模态、便于极点配置与观测器设计。

**知识点1.6　线性系统解与状态转移矩阵**

定义：状态转移矩阵 $\Phi(t,t_0)$ 是把初始状态映射到当前状态的算子。

假设：$A(t)$ 与 $A(\tau)$ 可交换（对时变系统）或 $A$ 为常阵。

公式（LTI）：$\Phi(t)=e^{At}$，解为
$$x(t)=e^{A(t-t_0)}x(t_0)+\int_{t_0}^{t}e^{A(t-\tau)}Bu(\tau)\,d\tau.$$

推导（常数变易法）：齐次解 $x_h=e^{A(t-t_0)}x_0$。设非齐次解 $x(t)=e^{A(t-t_0)}c(t)$，求导
$$\dot x=Ae^{A(t-t_0)}c+e^{A(t-t_0)}\dot c=Ae^{A(t-t_0)}c+Bu,$$
因为 $\dot x=Ax+Bu=Ae^{A(t-t_0)}c+Bu$，故 $e^{A(t-t_0)}\dot c=Bu$，即
$$\dot c=e^{-A(t-t_0)}Bu,\qquad c(t)=x_0+\int_{t_0}^{t}e^{-A(\tau-t_0)}Bu(\tau)d\tau,$$
代回得结论。

物理/几何意义：齐次项为初值的“自由漂移”，卷积项为输入的“受迫响应”。

极限/特例：$u\equiv0$ 时退化为齐次解；$A$ 稳定时 $e^{At}\to0$。

典型应用：系统响应仿真、采样离散化。

**知识点1.7　矩阵指数的性质与计算**

定义：$e^{At}=\sum_{k=0}^{\infty}\dfrac{A^k t^k}{k!}$。

假设：级数对任意有限 $t$ 绝对收敛。

公式与性质：
$$e^{A\cdot0}=I,\quad \frac{d}{dt}e^{At}=Ae^{At}=e^{At}A,\quad e^{A(t_1+t_2)}=e^{At_1}e^{At_2},\quad (e^{At})^{-1}=e^{-At}.$$
若 $AB=BA$，则 $e^{A+B}=e^Ae^B$；一般 $e^{A+B}\neq e^Ae^B$。

推导（求导性质）：逐项求导
$$\frac{d}{dt}\sum_{k=0}^\infty\frac{A^kt^k}{k!}=\sum_{k=1}^\infty\frac{A^k t^{k-1}}{(k-1)!}=A\sum_{k=0}^\infty\frac{A^kt^k}{k!}=Ae^{At}.$$

凯莱–哈密顿法：$n$ 阶矩阵满足特征方程 $A^n+a_{n-1}A^{n-1}+\cdots+a_0I=0$，故
$$e^{At}=\sum_{k=0}^{n-1}\alpha_k(t)A^k,$$
其中 $\alpha_k$ 由特征值 $\lambda_i$ 满足 $e^{\lambda_i t}=\sum_{k=0}^{n-1}\alpha_k(t)\lambda_i^k$ 确定。

极限/特例：$A$ 对角化 $A=V\Lambda V^{-1}$ 时 $e^{At}=Ve^{\Lambda t}V^{-1}$。

典型应用：解析求解线性系统、判别稳定性、采样周期离散化。

### 第2章 分析力学基础：拉格朗日方程

**知识点2.1　广义坐标、约束与自由度**

定义：用 $n$ 个独立参数 $q_1,\dots,q_n$ 确定系统位形，$q_i$ 称广义坐标；独立参数个数 $n$ 为自由度。

假设：约束为完整约束（可积分为 $f(r_1,\dots,r_N,t)=0$）或非完整约束；本节重点完整约束。

公式：对 $N$ 个质点、$k$ 个完整约束，自由度 $n=3N-k$。位置表示为
$$r_i=r_i(q_1,\dots,q_n,t).$$

推导（速度）：链式法则
$$v_i=\dot r_i=\sum_{j=1}^{n}\frac{\partial r_i}{\partial q_j}\dot q_j+\frac{\partial r_i}{\partial t}.$$

物理/几何意义：广义坐标把受约束运动限制在 $n$ 维位形流形上，约束力不做虚功因而被自动消去。

极限/特例：$n=1$ 单自由度；非完整约束（如冰刀、滚动盘）不能用此直接消元，需用拉格朗日乘子。

典型应用：机械臂、多体航天器、车辆动力学。

**知识点2.2　虚位移、虚功与达朗贝尔原理**

定义：虚位移 $\delta r_i$ 是在固定时刻满足约束的无穷小位移；达朗贝尔把惯性力 $-m_i\ddot r_i$ 视为附加的“有效力”。

假设：理想约束，即约束力在任意虚位移上做功之和为零 $\sum_i R_i\cdot\delta r_i=0$。

公式（达朗贝尔–拉格朗日原理）：
$$\sum_{i}(F_i-m_i\ddot r_i)\cdot\delta r_i=0.$$

推导：对每个质点写牛顿第二定律 $F_i+R_i=m_i\ddot r_i$，即 $F_i-m_i\ddot r_i=-R_i$。两边点乘 $\delta r_i$ 并求和：
$$\sum_i(F_i-m_i\ddot r_i)\cdot\delta r_i=-\sum_i R_i\cdot\delta r_i=0,$$
由理想约束假设得结论。

物理/几何意义：把动力学问题转化为在约束流形上的“平衡”问题，约束力被消去。

极限/特例：无约束时 $\delta r_i$ 独立，退化为牛顿方程。

典型应用：推导受约束刚体、连杆系统动力学。

**知识点2.3　拉格朗日方程的完整推导**

定义：拉格朗日函数 $L=T-V$，其中 $T$ 为动能、$V$ 为势能（广义力可由 $V$ 或 $Q_j$ 表示）。

假设：约束理想、完整；$L$ 对 $q,\dot q$ 足够光滑。

公式：
$$\frac{d}{dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}=Q_j,\qquad j=1,\dots,n,$$
其中 $Q_j=\sum_i F_i\cdot\dfrac{\partial r_i}{\partial q_j}$ 为非势广义力。

推导：由达朗贝尔原理，把 $\delta r_i=\sum_j\frac{\partial r_i}{\partial q_j}\delta q_j$ 代入：
$$\sum_j\left[\sum_i F_i\cdot\frac{\partial r_i}{\partial q_j}-\sum_i m_i\ddot r_i\cdot\frac{\partial r_i}{\partial q_j}\right]\delta q_j=0.$$
因 $\delta q_j$ 独立，括号为零，即
$$\sum_i m_i\ddot r_i\cdot\frac{\partial r_i}{\partial q_j}=Q_j.\qquad(\ast)$$
利用 $\dfrac{\partial v_i}{\partial\dot q_j}=\dfrac{\partial r_i}{\partial q_j}$ 及
$$\frac{d}{dt}\frac{\partial r_i}{\partial q_j}=\frac{\partial v_i}{\partial q_j},$$
有
$$\sum_i m_i\ddot r_i\cdot\frac{\partial r_i}{\partial q_j}=\sum_i m_i\left[\frac{d}{dt}\!\left(v_i\cdot\frac{\partial v_i}{\partial\dot q_j}\right)-v_i\cdot\frac{\partial v_i}{\partial q_j}\right]
=\frac{d}{dt}\frac{\partial T}{\partial\dot q_j}-\frac{\partial T}{\partial q_j}.$$
代入 $(\ast)$ 得
$$\frac{d}{dt}\frac{\partial T}{\partial\dot q_j}-\frac{\partial T}{\partial q_j}=Q_j.$$
若力均来自势能 $Q_j=-\partial V/\partial q_j$，则
$$\frac{d}{dt}\frac{\partial (T-V)}{\partial\dot q_j}-\frac{\partial (T-V)}{\partial q_j}=0,$$
即欧拉–拉格朗日方程。

物理/几何意义：拉格朗日方程是约束流形位形空间上的二阶微分方程，$T-V$ 反映运动的能量结构。

极限/特例：$L$ 不显含某 $q_j$（循环坐标）时对应广义动量守恒 $\dfrac{d}{dt}\dfrac{\partial L}{\partial\dot q_j}=0$。

典型应用：双摆、机械臂、卫星太阳翼展开动力学建模。

**知识点2.4　动能、惯性矩阵与广义力**

定义：动能一般形式 $T=\dfrac12\dot q^TM(q)\dot q$，$M(q)$ 为对称正定惯性矩阵。

假设：质点系总动能由平动与转动组成。

公式：对第 $i$ 刚体
$$T_i=\frac12 m_i v_{ci}^Tv_{ci}+\frac12\omega_i^TI_i\omega_i,$$
求和得 $T=\frac12\sum_i\dot q_i^TM_i\dot q_i$。

推导（展开拉格朗日方程）：对 $T=\frac12\dot q^TM(q)\dot q$，
$$\frac{\partial T}{\partial\dot q}=M(q)\dot q,\qquad \frac{d}{dt}\frac{\partial T}{\partial\dot q}=M\ddot q+\dot M\dot q,$$
$$\frac{\partial T}{\partial q_j}=\frac12\dot q^T\frac{\partial M}{\partial q_j}\dot q.$$
由 $\dot M=\sum_k\frac{\partial M}{\partial q_k}\dot q_k$，得
$$M\ddot q+\dot M\dot q-\frac12\begin{bmatrix}\dot q^T\frac{\partial M}{\partial q_1}\dot q\\ \vdots\\ \dot q^T\frac{\partial M}{\partial q_n}\dot q\end{bmatrix}=Q.$$
用克里斯托费尔符号 $\Gamma_{jk}^i=\frac12\sum_l M^{il}\!\left(\partial_jM_{lk}+\partial_kM_{lj}-\partial_lM_{jk}\right)$，可写成
$$\ddot q^i+\sum_{j,k}\Gamma_{jk}^i\dot q^j\dot q^k=(M^{-1}Q)^i,$$
即位形流形上的测地线型方程加外力。

物理/几何意义：惯性矩阵是对称正定的“质量度规”，克里斯托费尔符号表示坐标弯曲带来的惯性力。

极限/特例：$M$ 为常数时 $\dot M=0$，方程退化为 $M\ddot q=Q$。

典型应用：多连杆机器人前向动力学、航天器多体动力学递推算法。

**知识点2.5　保守系统的能量守恒与诺特定理**

定义：若 $L$ 不显含时间 $t$，则哈密顿量 $H=\sum_j\dot q_j\dfrac{\partial L}{\partial\dot q_j}-L$ 守恒。

假设：$L(q,\dot q)$ 对 $t$ 无显式依赖。

公式：
$$\frac{dH}{dt}=0,\qquad H=\sum_j\dot q_jp_j-L,\quad p_j=\frac{\partial L}{\partial\dot q_j}.$$

推导：
$$\frac{dH}{dt}=\sum_j\left(\dot q_j\frac{d}{dt}\frac{\partial L}{\partial\dot q_j}+\ddot q_j\frac{\partial L}{\partial\dot q_j}\right)-\frac{dL}{dt}.$$
由欧拉–拉格朗日方程 $\frac{d}{dt}\frac{\partial L}{\partial\dot q_j}=\frac{\partial L}{\partial q_j}$，故
$$\frac{dH}{dt}=\sum_j\left(\dot q_j\frac{\partial L}{\partial q_j}+\ddot q_j\frac{\partial L}{\partial\dot q_j}\right)-\left(\sum_j\frac{\partial L}{\partial q_j}\dot q_j+\sum_j\frac{\partial L}{\partial\dot q_j}\ddot q_j+\frac{\partial L}{\partial t}\right)=-\frac{\partial L}{\partial t}.$$
故 $\partial L/\partial t=0\Rightarrow dH/dt=0$。

物理/几何意义：哈密顿量是系统总机械能（在自然系统中），时间平移对称性对应能量守恒（诺特定理）。

极限/特例：$L$ 显含 $t$（如时变外场）则能量不守恒，变化率等于 $-\partial L/\partial t$。

典型应用：判定无阻尼机械系统能量守恒、构造李雅普诺夫函数。

**知识点2.6　拉格朗日方程与牛顿方程的等价性**

定义：两种方法给出同一运动，只是描述框架不同。

假设：完整约束、理想约束。

公式：牛顿方程 $m_i\ddot r_i=F_i+R_i$；拉格朗日方程由知识点2.3给出。

推导：由 $r_i=r_i(q,t)$ 两边求二阶导
$$\ddot r_i=\sum_j\frac{\partial r_i}{\partial q_j}\ddot q_j+\sum_{j,k}\frac{\partial^2 r_i}{\partial q_j\partial q_k}\dot q_j\dot q_k+2\sum_j\frac{\partial^2r_i}{\partial q_j\partial t}\dot q_j+\frac{\partial^2r_i}{\partial t^2}.$$
把牛顿方程两端左乘 $\frac{\partial r_i}{\partial q_j}$ 并求和，约束力项因理想约束消去，惯性项经知识点2.3的恒等变形化为拉格朗日项，故二者等价。

物理/几何意义：牛顿法在笛卡尔空间并显式含约束力；拉格朗日法在广义坐标空间且自动消去约束力。

极限/特例：$n=3N$、$q_i$ 取笛卡尔坐标时拉格朗日方程退化为牛顿方程。

典型应用：复杂约束系统优先用拉格朗日法，简单单质点用牛顿法。

### 第3章 哈密顿力学与正则方程

**知识点3.1　勒让德变换与哈密顿函数**

定义：勒让德变换把以 $(q,\dot q)$ 为自变量的 $L$ 换成以 $(q,p)$ 为自变量的 $H$，其中 $p_j=\partial L/\partial\dot q_j$。

假设：映射 $\dot q\mapsto p=\partial L/\partial\dot q$ 可逆，即赫森矩阵 $\partial^2L/\partial\dot q\partial\dot q$ 非奇异。

公式：
$$H(q,p,t)=\sum_j p_j\dot q_j-L(q,\dot q,t),\qquad \dot q_j=\dot q_j(q,p,t).$$

推导：全微分
$$dL=\sum_j\frac{\partial L}{\partial q_j}dq_j+\sum_j\frac{\partial L}{\partial\dot q_j}d\dot q_j+\frac{\partial L}{\partial t}dt
=\sum_j\dot p_jdq_j+\sum_jp_jd\dot q_j+\frac{\partial L}{\partial t}dt.$$
又 $dH=\sum_j\dot q_jdp_j+\sum_jp_jd\dot q_j-dL$，代入得
$$dH=\sum_j\dot q_jdp_j-\sum_j\dot p_jdq_j-\frac{\partial L}{\partial t}dt,$$
故 $\dot q_j=\partial H/\partial p_j$，$\dot p_j=-\partial H/\partial q_j$。

物理/几何意义：勒让德变换实现从切丛 $(q,\dot q)$ 到余切丛 $(q,p)$ 的相空间描述。

极限/特例：$L=T-V$ 且 $T$ 为速度二次型时 $H=T+V$。

典型应用：相空间分析、辛几何、量子力学对应原理。

**知识点3.2　哈密顿正则方程的推导**

定义：正则方程是相空间中的一阶对称方程组。

假设：同上，$H$ 可微。

公式：
$$\dot q_j=\frac{\partial H}{\partial p_j},\qquad \dot p_j=-\frac{\partial H}{\partial q_j},\qquad j=1,\dots,n.$$

推导：从 $H=\sum_k p_k\dot q_k-L$ 出发，对 $q_j$ 求偏导
$$\frac{\partial H}{\partial q_j}=p_j\frac{\partial\dot q_j}{\partial q_j}-\frac{\partial L}{\partial q_j}-\frac{\partial L}{\partial\dot q_j}\frac{\partial\dot q_j}{\partial q_j}=-\frac{\partial L}{\partial q_j}.$$
由欧拉–拉格朗日方程 $\frac{\partial L}{\partial q_j}=\frac{d}{dt}\frac{\partial L}{\partial\dot q_j}=\dot p_j$，故 $\dot p_j=-\partial H/\partial q_j$。对 $p_j$ 求偏导
$$\frac{\partial H}{\partial p_j}=\dot q_j+p_j\frac{\partial\dot q_j}{\partial p_j}-\frac{\partial L}{\partial\dot q_j}\frac{\partial\dot q_j}{\partial p_j}=\dot q_j,$$
其中用了 $p_j=\partial L/\partial\dot q_j$。证毕。

物理/几何意义：正则方程把 $2n$ 维相空间中的流写成辛形式 $\dot z=J\nabla H$，$J=\begin{bmatrix}0&I\\-I&0\end{bmatrix}$。

极限/特例：$H$ 不显含 $t$ 时能量守恒；$H$ 不显含某坐标时 $p_j$ 守恒。

典型应用：天体力学、最优控制中的哈密顿系统。

**知识点3.3　泊松括号与守恒量**

定义：对相空间函数 $F(q,p),G(q,p)$，
$$\{F,G\}=\sum_j\left(\frac{\partial F}{\partial q_j}\frac{\partial G}{\partial p_j}-\frac{\partial F}{\partial p_j}\frac{\partial G}{\partial q_j}\right).$$

假设：$F,G$ 可微。

公式：$\dfrac{dF}{dt}=\{F,H\}+\dfrac{\partial F}{\partial t}$；若 $\{F,H\}=0$ 且 $\partial F/\partial t=0$，则 $F$ 守恒。

推导：由正则方程
$$\frac{dF}{dt}=\sum_j\left(\frac{\partial F}{\partial q_j}\dot q_j+\frac{\partial F}{\partial p_j}\dot p_j\right)+\frac{\partial F}{\partial t}
=\sum_j\left(\frac{\partial F}{\partial q_j}\frac{\partial H}{\partial p_j}-\frac{\partial F}{\partial p_j}\frac{\partial H}{\partial q_j}\right)+\frac{\partial F}{\partial t},$$
即得。

物理/几何意义：泊松括号是相空间函数的李代数结构，守恒量对应与哈密顿量对易的量。

极限/特例：$\{q_i,q_j\}=0$，$\{p_i,p_j\}=0$，$\{q_i,p_j\}=\delta_{ij}$。

典型应用：寻找运动积分、判定可积性。

**知识点3.4　正则变换**

定义：保持正则方程形式不变的相空间坐标变换 $(q,p)\to(Q,P)$ 称为正则变换。

假设：变换可逆且足够光滑。

公式：存在生成函数 $S$ 使
$$p_j=\frac{\partial S}{\partial q_j},\quad Q_j=\frac{\partial S}{\partial P_j},\quad K=H+\frac{\partial S}{\partial t}.$$

推导：要求变换后仍满足 $\dot Q=\partial K/\partial P$、$\dot P=-\partial K/\partial Q$，等价于存在 $S$ 使
$$\sum_j p_jdq_j-Hdt=\sum_jP_jdQ_j-Kdt+dS.$$
若 $S=S(q,P,t)$，则 $dS=\sum_j\frac{\partial S}{\partial q_j}dq_j+\sum_j\frac{\partial S}{\partial P_j}dP_j+\frac{\partial S}{\partial t}dt$，比较系数即得。

物理/几何意义：正则变换保持相空间体积（刘维尔定理）与辛结构。

极限/特例：$S=\sum q_jP_j$ 给出恒等变换。

典型应用：用哈密顿–雅可比方法求解可积系统。

**知识点3.5　哈密顿–雅可比方程**

定义：寻找生成函数 $S(q,t)$ 使新哈密顿量 $K=0$，从而新坐标守恒。

假设：$H$ 可微，$S$ 二次可微。

公式：
$$\frac{\partial S}{\partial t}+H\!\left(q,\frac{\partial S}{\partial q},t\right)=0.$$

推导：取 $S(q,P,t)$，由知识点3.4，$p_j=\partial S/\partial q_j$，且 $K=H+\partial S/\partial t$。要求 $K=0$ 得
$$\frac{\partial S}{\partial t}+H(q,\partial S/\partial q,t)=0.$$
解得 $S(q,\alpha,t)$，其中 $\alpha$ 为积分常数（新动量），再由 $\beta_j=\partial S/\partial\alpha_j$ 给出轨线。

物理/几何意义：把 $2n$ 维正则方程化为一个一阶偏微分方程，是经典力学与量子力学（WKB）的桥梁。

极限/特例：$S=-Et+W(q)$（$H$ 不显含 $t$）时得定态方程 $H(q,\partial W/\partial q)=E$。

典型应用：可积系统求解；与最优控制的 HJB 方程同源。

**知识点3.6　哈密顿力学与最优控制的对应**

定义：最优控制的哈密顿量与经典力学哈密顿量在形式上一致，协态对应广义动量。

假设：最优控制问题有可微值函数。

公式：对系统 $\dot x=f(x,u,t)$、代价 $\int L(x,u,t)dt$，定义
$$\mathcal H(x,u,\lambda,t)=L(x,u,t)+\lambda^Tf(x,u,t).$$
协态方程 $\dot\lambda=-\partial\mathcal H/\partial x$，最优性 $\partial\mathcal H/\partial u=0$。

推导：把控制 $u$ 视为“速度”，把 $\lambda$ 视为动量，则 $\mathcal H=L+\lambda^Tf$ 与 $H=p^T\dot q-L$ 结构相同（符号相反）。对 $\lambda$ 求偏导给状态方程 $\dot x=\partial\mathcal H/\partial\lambda$；对 $x$ 求偏导并对照正则方程 $\dot p=-\partial H/\partial q$ 给 $\dot\lambda=-\partial\mathcal H/\partial x$。

物理/几何意义：最优轨线是扩展相空间中的哈密顿轨线，庞特里亚金极小值原理即“哈密顿量取极值”。

极限/特例：$L=\frac12u^2$、$f=u$ 时 $\mathcal H=\frac12u^2+\lambda u$，最优 $u=-\lambda$，退化为经典质点。

典型应用：轨迹优化、火箭制导。

### 第4章 连续介质与结构动力学建模（N–S、本构、能量法）

**知识点4.1　连续体运动学：应变与速度梯度**

定义：位移场 $u(x,t)$，小应变张量 $\varepsilon=\frac12(\nabla u+\nabla u^T)$；速度梯度 $L=\nabla v$。

假设：连续介质、小变形（线性应变）。

公式：
$$\varepsilon_{ij}=\frac12\!\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\right),\qquad D=\frac12(L+L^T),\quad W=\frac12(L-L^T).$$

推导：两邻近物质点 $dx$ 变形后变为 $dx'=dx+du$，$du=\nabla u\,dx$。长度平方变化
$$|dx'|^2-|dx|^2=dx^T(\nabla u+\nabla u^T+(\nabla u)^T\nabla u)dx.$$
小变形下忽略二阶项，得 $2\,dx^T\varepsilon\,dx$，故 $\varepsilon$ 对称部分。速度场同理分解为应变率 $D$ 与旋转率 $W$。

物理/几何意义：$\varepsilon$ 度量长度与角度的相对变化；$W$ 表示刚体转动。

极限/特例：大变形需用格林–拉格朗日应变 $E=\frac12(F^TF-I)$。

典型应用：结构应力分析、流体应变率。

**知识点4.2　应力张量与柯西动量方程**

定义：柯西应力张量 $\sigma$ 通过 $t=\sigma n$ 给出面力；动量守恒给出柯西运动方程。

假设：连续介质、无体力偶。

公式：
$$\rho\frac{Dv}{Dt}=\nabla\cdot\sigma+f,\qquad \frac{D}{Dt}=\frac{\partial}{\partial t}+v\cdot\nabla.$$

推导：对任意控制体 $\Omega$，
$$\frac{d}{dt}\int_\Omega\rho v\,dV=\int_{\partial\Omega}\sigma n\,dS+\int_\Omega f\,dV.$$
由雷诺输运定理，左端 $=\int_\Omega\rho\frac{Dv}{Dt}dV$；由散度定理 $\int_{\partial\Omega}\sigma n\,dS=\int_\Omega\nabla\cdot\sigma\,dV$。因 $\Omega$ 任意，得
$$\rho\frac{Dv}{Dt}=\nabla\cdot\sigma+f.$$

物理/几何意义：牛顿第二定律的连续介质形式，$\nabla\cdot\sigma$ 为内力合力密度。

极限/特例：静力平衡时 $Dv/Dt=0$ 得 $\nabla\cdot\sigma+f=0$。

典型应用：流体与固体动力学统一框架。

**知识点4.3　本构方程：线弹性与牛顿粘性**

定义：本构方程联系应力与运动学量，刻画材料性质。

假设：小变形线弹性（胡克定律）；牛顿流体（应力与应变率线性）。

公式：
$$\text{固体：}\ \sigma_{ij}=C_{ijkl}\varepsilon_{kl};\qquad \text{各向同性：}\ \sigma=2\mu\varepsilon+\lambda\,(\mathrm{tr}\varepsilon)I.$$
$$\text{流体：}\ \sigma=-pI+\lambda(\nabla\cdot v)I+2\mu D,\qquad D=\frac12(\nabla v+\nabla v^T).$$

推导（各向同性弹性）：设 $\sigma$ 是 $\varepsilon$ 的线性、各向同性函数，由表示定理
$$\sigma=\alpha I+\beta\varepsilon+\gamma\varepsilon^2+\cdots$$
线性项给出 $\sigma=\alpha\,\mathrm{tr}(\varepsilon)I+2\mu\varepsilon$。定义 $\lambda=\alpha$，得 $\sigma=\lambda(\mathrm{tr}\varepsilon)I+2\mu\varepsilon$。流体本构同理由“各向同性张量函数 + 静止时 $\sigma=-pI$”得到，取 $\varepsilon$ 换成 $D$。

物理/几何意义：$\lambda,\mu$ 为拉梅常数，$\mu$ 为剪切模量；流体中 $\mu$ 为动力粘度。

极限/特例：不可压缩流体 $\nabla\cdot v=0$，$\sigma=-pI+2\mu D$；不可压缩弹性 $\mathrm{tr}\varepsilon=0$ 只由 $\mu$ 决定。

典型应用：结构有限元、空气动力学、生物力学。

**知识点4.4　纳维–斯托克斯方程的推导**

定义：把牛顿粘性本构代入柯西动量方程得到的流体运动方程。

假设：牛顿流体、各向同性、常粘度、可能可压缩。

公式：
$$\rho\!\left(\frac{\partial v}{\partial t}+v\cdot\nabla v\right)=-\nabla p+\lambda\nabla(\nabla\cdot v)+2\nabla\cdot(\mu D)+f.$$
常粘度、$\lambda,\mu$ 常数：
$$\rho\!\left(\frac{\partial v}{\partial t}+v\cdot\nabla v\right)=-\nabla p+(\lambda+\mu)\nabla(\nabla\cdot v)+\mu\nabla^2v+f.$$
不可压缩：
$$\rho\!\left(\frac{\partial v}{\partial t}+v\cdot\nabla v\right)=-\nabla p+\mu\nabla^2v+f,\qquad \nabla\cdot v=0.$$

推导：由 $\sigma=-pI+\lambda(\nabla\cdot v)I+2\mu D$，
$$\nabla\cdot\sigma=-\nabla p+\nabla(\lambda\nabla\cdot v)+2\nabla\cdot(\mu D).$$
利用恒等式 $\nabla\cdot(\nabla v)= \nabla(\nabla\cdot v)$、$\nabla\cdot(\nabla v^T)=\nabla(\nabla\cdot v)$，故
$$2\nabla\cdot D=\nabla^2v+\nabla(\nabla\cdot v).$$
代入柯西方程得结论。

物理/几何意义：左端为惯性，右端为压力梯度、粘性扩散与体积粘性。$\mu\nabla^2v$ 使动量扩散。

极限/特例：$\mu=0$ 得欧拉方程；$\mathrm{Re}=UL/\nu\gg1$ 边界层近似；$\mathrm{Re}\ll1$ 得斯托克斯方程。

典型应用：CFD、气动加热、湍流建模。

**知识点4.5　能量法：虚功原理与最小势能原理**

定义：虚功原理指平衡时外力虚功等于内力虚功；最小势能原理指稳定平衡时总势能取局部极小。

假设：保守系统、小变形、线性弹性。

公式：
$$\delta W_{\text{ext}}=\delta U,\qquad \Pi=U-W_{\text{ext}},\qquad \delta\Pi=0,\quad \delta^2\Pi>0.$$

推导：由虚功原理 $\int\sigma:\delta\varepsilon\,dV=\int f\cdot\delta u\,dV+\int t\cdot\delta u\,dS$。对线弹性 $\sigma=C:\varepsilon$，
$$\delta U=\int\varepsilon:C:\delta\varepsilon\,dV=\delta\!\left(\frac12\int\varepsilon:C:\varepsilon\,dV\right).$$
故平衡方程等价于 $\delta(U-W_{\text{ext}})=0$。稳定平衡要求二阶变分正定。

物理/几何意义：势能泛函的驻点给出平衡位形，极值类型判别稳定性。

极限/特例：线性问题能量泛函二次，可化为线性代数方程组。

典型应用：结构位移法、有限元、稳定屈曲分析。

**知识点4.6　瑞利–里兹法与伽辽金法**

定义：用有限个试函数的线性组合近似解，把泛函驻值问题化为代数方程。

假设：试函数满足几何边界条件；瑞利–里兹要求泛函存在，伽辽金要求残差加权正交。

公式：设 $u\approx\sum_{i=1}^N c_i\phi_i$，瑞利–里兹：
$$\frac{\partial \Pi}{\partial c_i}=0\ \Rightarrow\ Kc=f,\qquad K_{ij}=\int \phi_i\,\mathcal L\,\phi_j\,dV.$$
伽辽金：残差 $R=\mathcal L u-f$，要求
$$\int\phi_i R\,dV=0.$$

推导（瑞利–里兹）：$\Pi(c)=\frac12c^TKc-c^Tf$，令 $\partial\Pi/\partial c=Kc-f=0$，其中 $K$ 由应变能的二次型系数得到。伽辽金由弱形式
$$\int\phi_i(\mathcal Lu-f)dV=0$$
分部积分回到变分方程。

物理/几何意义：把无穷维变分问题投影到有限维子空间，是最小二乘意义下的最优近似。

极限/特例：$N\to\infty$ 在完备基下收敛到精确解；对称正定问题瑞利–里兹与伽辽金等价。

典型应用：板壳振动、流固耦合近似求解。

**知识点4.7　模态叠加与结构动力学**

定义：线性结构运动方程 $M\ddot q+C\dot q+Kq=F$ 经模态变换解耦。

假设：$M,K$ 对称正定，$C$ 为比例阻尼 $C=\alpha M+\beta K$。

公式：无阻尼特征问题 $K\phi_i=\omega_i^2M\phi_i$，模态矩阵 $\Phi=[\phi_1,\dots,\phi_n]$ 使
$$\Phi^TM\Phi=I,\qquad \Phi^TK\Phi=\mathrm{diag}(\omega_i^2).$$
令 $q=\Phi\eta$，得
$$\ddot\eta_i+2\zeta_i\omega_i\dot\eta_i+\omega_i^2\eta_i=\phi_i^TF.$$

推导：代入 $q=\Phi\eta$，左乘 $\Phi^T$：
$$\Phi^TM\Phi\ddot\eta+\Phi^TC\Phi\dot\eta+\Phi^TK\Phi\eta=\Phi^TF.$$
由 $M,K$ 正交性及比例阻尼 $\Phi^TC\Phi=\alpha I+\beta\,\mathrm{diag}(\omega_i^2)$，各方程解耦。第 $i$ 式除以质量归一化后得所示形式。

物理/几何意义：模态是结构的固有振型，任意响应为各模态响应之和。

极限/特例：截断低阶模态可降阶；共振时 $\omega=\omega_i$ 响应放大 $1/(2\zeta_i)$。

典型应用：航天器柔性附件、桥梁风振、模态试验。

### 第5章 线性系统的运动分析与稳定性

**知识点5.1　状态转移矩阵及其性质**

定义：$\Phi(t,t_0)$ 是使 $x(t)=\Phi(t,t_0)x(t_0)$ 成立的 $n\times n$ 矩阵。

假设：$\dot x=A(t)x$ 的解存在唯一。

公式与性质：
$$\frac{\partial}{\partial t}\Phi(t,t_0)=A(t)\Phi(t,t_0),\quad \Phi(t_0,t_0)=I,$$
$$\Phi(t_2,t_0)=\Phi(t_2,t_1)\Phi(t_1,t_0),\quad \Phi^{-1}(t,t_0)=\Phi(t_0,t).$$

推导（半群性质）：$x(t_2)=\Phi(t_2,t_1)x(t_1)=\Phi(t_2,t_1)\Phi(t_1,t_0)x(t_0)$，又 $x(t_2)=\Phi(t_2,t_0)x(t_0)$，对任意 $x(t_0)$ 成立，故 $\Phi(t_2,t_0)=\Phi(t_2,t_1)\Phi(t_1,t_0)$。

物理/几何意义：$\Phi$ 是状态空间上的线性流映射，满足群（半群）结构。

极限/特例：LTI 时 $\Phi(t,t_0)=e^{A(t-t_0)}$；$A$ 常数且可交换时结论简单。

典型应用：离散化、能控性 Gramian 构造。

**知识点5.2　矩阵指数的计算方法**

定义：$e^{At}$ 的计算方法包括级数、拉普拉斯反变换、对角化/约当、凯莱–哈密顿。

假设：$A$ 为常阵。

公式：
$$e^{At}=\mathcal L^{-1}\left[(sI-A)^{-1}\right],\qquad e^{At}=Ve^{\Lambda t}V^{-1}.$$

推导（拉普拉斯法）：对 $\dot x=Ax$ 取拉普拉斯变换 $sX(s)-x(0)=AX(s)$，则 $X(s)=(sI-A)^{-1}x(0)$，取反变换与 $x(t)=e^{At}x(0)$ 对比即得。对角化法：$A=V\Lambda V^{-1}$，则
$$e^{At}=\sum_k\frac{(V\Lambda V^{-1})^kt^k}{k!}=V\left(\sum_k\frac{\Lambda^kt^k}{k!}\right)V^{-1}=Ve^{\Lambda t}V^{-1}.$$

物理/几何意义：$(sI-A)^{-1}$ 是预解矩阵，极点即特征值，决定模态衰减率。

极限/特例：重根时用约当块，$e^{At}$ 含 $t^ke^{\lambda t}$ 项。

典型应用：解析仿真、稳定性判别、控制器设计。

**知识点5.3　模态分解与自由响应**

定义：把状态响应分解为各特征模态的线性组合。

假设：$A$ 可对角化，特征值 $\lambda_i$、左右特征向量 $v_i,w_i$。

公式：
$$x(t)=\sum_{i=1}^n e^{\lambda_i t}v_iw_i^Tx(0).$$

推导：设 $x=\sum_i z_iv_i$，由 $A v_i=\lambda_iv_i$ 得 $\dot z_i=\lambda_iz_i$，$z_i(t)=e^{\lambda_it}z_i(0)$。由双正交性 $w_i^Tv_j=\delta_{ij}$ 得 $z_i(0)=w_i^Tx(0)$，故结论。

物理/几何意义：各模态独立演化，$\mathrm{Re}\lambda_i<0$ 时衰减，$\mathrm{Re}\lambda_i>0$ 时发散。

极限/特例：复特征值成对出现对应振荡模态 $\lambda=-\sigma\pm j\omega$。

典型应用：振荡模态分析、参与因子计算。

**知识点5.4　线性系统稳定性判据**

定义：平衡点 $x=0$ 渐近稳定指数值解趋于零；系统稳定指所有特征值实部非正且临界特征值对应简单根。

假设：LTI 系统 $\dot x=Ax$。

公式：渐近稳定 $\iff \mathrm{Re}\lambda_i(A)<0,\ \forall i$；李雅普诺夫方程 $A^TP+PA=-Q$ 对任意 $Q=Q^T>0$ 有唯一正定解 $P$。

推导（充分性）：取 $V(x)=x^TPx$，$P>0$。沿轨线
$$\dot V=x^T(A^TP+PA)x=-x^TQx<0.$$
故 $V$ 为正定、$\dot V$ 负定，系统渐近稳定。

物理/几何意义：特征值位于左半平面；$V$ 是能量型函数，$\dot V<0$ 表示能量耗散。

极限/特例：虚轴上单重根为临界稳定（不渐近）；重根则不稳定。

典型应用：Routh–Hurwitz 判据、根轨迹、李雅普诺夫方程数值求解。

**知识点5.5　连续系统离散化**

定义：在采样周期 $T$ 下把连续系统化为离散状态方程。

假设：输入在采样间隔内零阶保持（ZOH）。

公式：
$$x_{k+1}=A_dx_k+B_du_k,\quad A_d=e^{AT},\quad B_d=\int_0^Te^{A\tau}\,d\tau\,B.$$

推导：$x(t_{k+1})=e^{AT}x(t_k)+\int_{t_k}^{t_{k+1}}e^{A(t_{k+1}-\tau)}Bu(\tau)d\tau$。ZOH 下 $u(\tau)=u_k$，令 $s=t_{k+1}-\tau$，则
$$\int_{t_k}^{t_{k+1}}e^{A(t_{k+1}-\tau)}d\tau=\int_0^Te^{As}ds,$$
得 $A_d,B_d$。当 $A$ 可逆时 $B_d=A^{-1}(e^{AT}-I)B$。

物理/几何意义：离散化保持采样时刻状态一致，是数字控制的基础。

极限/特例：$T\to0$ 时 $A_d\to I+AT$、$B_d\to BT$；离散特征值 $\mu_i=e^{\lambda_iT}$。

典型应用：数字控制器实现、模型预测控制。

**知识点5.6　脉冲响应与卷积**

定义：单位脉冲响应 $g(t)=Ce^{At}B+D\delta(t)$，零状态响应为卷积。

假设：初始状态为零、系统线性时不变。

公式：
$$y(t)=\int_0^tg(t-\tau)u(\tau)d\tau+Ce^{At}x(0).$$

推导：零初值下由解公式
$$y(t)=C\int_0^te^{A(t-\tau)}Bu(\tau)d\tau+Du(t)=\int_0^tCe^{A(t-\tau)}Bu(\tau)d\tau+Du(t).$$
含 $D$ 项与 $\delta$ 卷积即 $Du$，故 $g(t)=Ce^{At}B+D\delta(t)$。

物理/几何意义：卷积表示系统对输入的线性叠加记忆。

极限/特例：传递函数 $G(s)=C(sI-A)^{-1}B+D=\mathcal L[g(t)]$。

典型应用：系统辨识、滤波器设计。

### 第6章 能控性

**知识点6.1　能控性的定义**

定义：状态能控指存在有限时间 $t_1>t_0$ 和容许控制 $u$，使系统从任意初态 $x_0$ 转移到原点；状态能达指从原点转移到任意 $x_1$。LTI 系统两者等价。

假设：系统 $\dot x=Ax+Bu$，$A,B$ 为常阵。

公式：能控 $\iff$ 对任意 $x_0$ 和 $t_1>t_0$，存在 $u$ 使
$$x(t_1)=e^{A(t_1-t_0)}x_0+\int_{t_0}^{t_1}e^{A(t_1-\tau)}Bu(\tau)d\tau=0.$$

推导（化归）：令目标为零，即要求
$$\int_{t_0}^{t_1}e^{A(t_1-\tau)}Bu(\tau)d\tau=-e^{A(t_1-t_0)}x_0.$$
这是关于 $u$ 的算子方程，若对任意右端有解则系统能控。

物理/几何意义：能控性刻画输入对状态各模态的“激励能力”。

极限/特例：单输入时能控即能影响所有模态；不可控模态对输入无响应。

典型应用：判断系统可否用状态反馈任意配置极点。

**知识点6.2　格拉姆矩阵判据**

定义：能控性 Gramian
$$W_c(t_0,t_1)=\int_{t_0}^{t_1}e^{A(t_0-\tau)}BB^Te^{A^T(t_0-\tau)}d\tau.$$

假设：$A,B$ 连续。

公式：系统能控 $\iff W_c(t_0,t_1)$ 非奇异。

推导（充分性）：若 $W_c$ 非奇异，取
$$u(t)=-B^Te^{A^T(t_0-t)}W_c^{-1}(t_0,t_1)e^{A(t_1-t_0)}x_0,$$
代入状态方程可验证 $x(t_1)=0$。计算
$$\int_{t_0}^{t_1}e^{A(t_1-\tau)}Bu(\tau)d\tau=-e^{A(t_1-t_0)}\!\left(\int_{t_0}^{t_1}e^{A(t_0-\tau)}BB^Te^{A^T(t_0-\tau)}d\tau\right)W_c^{-1}e^{A(t_1-t_0)}x_0
=-e^{A(t_1-t_0)}x_0,$$
故 $x(t_1)=0$。
必要性：若奇异，存在非零 $x_0$ 使 $W_cx_0=0$，则 $x_0^TW_cx_0=\int\|B^Te^{A^T(t_0-\tau)}x_0\|^2d\tau=0$，即 $B^Te^{A^T(t_0-t)}x_0\equiv0$。此时对任意 $u$，
$$x_0^Te^{A^T(t_0-t_1)}x(t_1)=\|x_0\|^2+\int(\cdots)u\,d\tau=\|x_0\|^2\neq0,$$
与 $x(t_1)=0$ 矛盾。

物理/几何意义：$W_c$ 是输入能“到达”的状态集合的椭球度量，奇异方向为不可控方向。

极限/特例：$A$ 稳定时稳态 Gramian 满足李雅普诺夫方程 $AW_c+W_cA^T=-BB^T$。

典型应用：可控性度量、平衡截断模型降阶。

**知识点6.3　卡尔曼秩判据**

定义：用 $(A,B)$ 构造能控性矩阵 $\mathcal C=[B\ AB\ \cdots\ A^{n-1}B]$。

假设：$A\in\mathbb R^{n\times n}$，$B\in\mathbb R^{n\times m}$。

公式：系统能控 $\iff \mathrm{rank}\,\mathcal C=n$。

推导（必要性）：若 $\mathrm{rank}\,\mathcal C<n$，存在非零 $q$ 使 $q^T\mathcal C=0$，即 $q^TB=0,q^TAB=0,\dots,q^TA^{n-1}B=0$。由凯莱–哈密顿，$A^k$ 可表为 $I,\dots,A^{n-1}$ 的线性组合，故 $q^TA^kB=0$ 对所有 $k\ge0$，从而 $q^Te^{At}B=0$。于是对任意 $u$，$q^T[x(t_1)-e^{A(t_1-t_0)}x_0]=0$。若系统能控则取 $x_0=q$、目标为零，得 $0=q^Te^{A(t_1-t_0)}q\neq q^Tq$ 矛盾，故必满秩。
充分性：$\mathrm{rank}\,\mathcal C=n$ 时对任意目标可用 $u(t)=\sum$ 脉冲构造，见知识点6.2。

物理/几何意义：$\mathcal C$ 的列张成的子空间为能控子空间。

极限/特例：$m=1$ 时 $\mathcal C$ 为方阵，满秩即行列式非零。

典型应用：软件 `ctrb`、系统能控性快速判别。

**知识点6.4　PBH 判据与模态能控性**

定义：PBH 判据用特征值判别能控性。

假设：$A,B$ 为常阵。

公式：$(A,B)$ 能控 $\iff$ 对所有 $\lambda\in\mathbb C$，$\mathrm{rank}[\lambda I-A\ \ B]=n$。等价地，左特征向量 $w^TA=\lambda w^T$ 满足 $w^TB\neq0$。

推导：若存在 $\lambda$ 与非零 $w$ 使 $w^TA=\lambda w^T$ 且 $w^TB=0$，则 $w^TA^kB=\lambda^kw^TB=0$，故 $\mathrm{rank}\,\mathcal C<n$，不可控。反之，若 $\mathrm{rank}\,\mathcal C<n$，存在 $w$ 使 $w^TA^kB=0$（$k=0,\dots,n-1$），取 $w$ 为 $A^T$ 的特征向量（$A^Tw=\lambda w$）即可，则 $w^TB=0$。

物理/几何意义：每个模态 $\lambda$ 必须能被输入“看见”，否则该模态不可控。

极限/特例：不可控模态对应 $w^TB=0$ 的特征方向。

典型应用：识别不可控模态、设计执行器位置。

**知识点6.5　能控标准型与极点配置的能控性条件**

定义：单输入能控系统可通过变换化为能控标准型。

假设：$(A,B)$ 能控，$B$ 为单列。

公式：存在 $T$ 使
$$\bar A=TAT^{-1}=\begin{bmatrix}0&1&&\\ \vdots&&\ddots&\\ 0&&&1\\ -a_0&-a_1&\cdots&-a_{n-1}\end{bmatrix},\quad \bar B=TB=\begin{bmatrix}0\\\vdots\\0\\1\end{bmatrix}.$$

推导：能控性矩阵 $\mathcal C=[B\ AB\ \cdots\ A^{n-1}B]$ 满秩，其特征多项式系数 $a_i$ 由 $\det(\lambda I-A)=\lambda^n+a_{n-1}\lambda^{n-1}+\cdots+a_0$ 给出。取 $T=\mathcal C^{-1}\mathcal C_c$（$\mathcal C_c$ 为能控标准型对应的能控性矩阵），则直接验证 $\bar A,\bar B$ 具标准型。

物理/几何意义：能控标准型把系统化为“积分链 + 反馈”，极点配置问题化为多项式的根配置。

极限/特例：多输入系统需化为龙伯格能控标准型。

典型应用：Ackermann 公式极点配置。

**知识点6.6　能控性指数与能达性**

定义：能控性指数 $\mu$ 是使 $\mathrm{rank}[B\ AB\ \cdots\ A^{k}B]=n$ 的最小 $k+1$；能达性关注从原点可达的状态集。

假设：一般多输入系统。

公式：能控性指数满足 $\mu\ge n/m$（$m$ 为输入维数）且 $\sum\mu_i=n$。

推导（估计）：每一列 $B$ 有 $m$ 个输入，每次乘 $A$ 至少新增一个线性无关方向；总维数 $n$，故需 $\lceil n/m\rceil$ 次。严格证明用卡尔曼分解的能控子空间维数。

物理/几何意义：$\mu$ 反映输入到状态的“传递深度”。

极限/特例：$m=n$ 且 $B$ 满秩时 $\mu=1$。

典型应用：最小拍控制、降阶建模。

### 第7章 能观性

**知识点7.1　能观性与能检性的定义**

定义：系统能观指由有限时间段内的输出 $y(t)$ 与已知输入 $u(t)$ 可唯一确定初始状态 $x_0$；若不可观部分渐近稳定则称能检测。

假设：$\dot x=Ax+Bu,\ y=Cx+Du$。

公式：能观 $\iff$ 对任意 $t_1>t_0$，
$$y(t)=Ce^{A(t-t_0)}x_0+\int_{t_0}^tCe^{A(t-\tau)}Bu(\tau)d\tau+Du(t)$$
能唯一确定 $x_0$。

推导：已知 $u$ 可算出输出中的已知项，剩下
$$\bar y(t)=Ce^{A(t-t_0)}x_0.$$
若存在两个初态 $x_0\neq x_0'$ 产生相同 $\bar y$，则 $Ce^{A(t-t_0)}(x_0-x_0')=0$，即存在非零 $z$ 使 $Ce^{A(t-t_0)}z=0$，系统不可观。

物理/几何意义：能观性刻画由测量反推内部状态的能力。

极限/特例：能观且稳定系统可构造观测器；不能检的对偶是不可稳。

典型应用：状态估计、故障诊断。

**知识点7.2　能观性秩判据与对偶原理**

定义：对偶系统为 $(A^T,C^T)$；能观性矩阵 $\mathcal O=\begin{bmatrix}C\\CA\\\vdots\\CA^{n-1}\end{bmatrix}$。

假设：$A,C$ 为常阵。

公式：能观 $\iff \mathrm{rank}\,\mathcal O=n$。

推导（对偶）：$(A,C)$ 能观 $\iff (A^T,C^T)$ 能控。由能控秩判据，$(A^T,C^T)$ 的能控性矩阵为
$$[C^T\ A^TC^T\ \cdots\ (A^T)^{n-1}C^T]=\mathcal O^T.$$
故 $\mathrm{rank}\,\mathcal O^T=n\iff\mathrm{rank}\,\mathcal O=n$。

物理/几何意义：能观性是能控性在转置系统上的镜像，观测器设计等价于对偶系统状态反馈设计。

极限/特例：PBH：能观 $\iff\mathrm{rank}\begin{bmatrix}\lambda I-A\\C\end{bmatrix}=n$，等价于右特征向量 $v$ 满足 $Cv\neq0$。

典型应用：观测器极点配置、对偶极点配置。

**知识点7.3　能观标准型**

定义：单输出能观系统可化为能观标准型，$C$ 为单位向量。

假设：$(A,C)$ 能观，$C$ 为单行。

公式：
$$\bar A=\begin{bmatrix}0&0&\cdots&0&-a_0\\ 1&0&\cdots&0&-a_1\\ 0&1&\cdots&0&-a_2\\ \vdots&&\ddots&&\vdots\\ 0&0&\cdots&1&-a_{n-1}\end{bmatrix},\quad \bar C=[0\ \cdots\ 0\ 1].$$

推导：由对偶原理，$(A^T,C^T)$ 能控，其能控标准型转置即得能观标准型。能观性矩阵满秩保证变换存在。

物理/几何意义：能观标准型使输出直接等于最后一个状态，观测器设计简化。

极限/特例：与能控标准型互为转置关系。

典型应用：观测器增益计算、系统辨识。

**知识点7.4　卡尔曼能观性分解与能控性分解**

定义：通过坐标变换把状态空间分解为能控/不能控、能观/不能观子空间。

假设：$(A,B,C)$ 一般系统。

公式：能控性分解
$$\bar A=\begin{bmatrix}A_{11}&A_{12}\\0&A_{22}\end{bmatrix},\quad \bar B=\begin{bmatrix}B_1\\0\end{bmatrix},\quad \bar C=[C_1\ C_2].$$
能观性分解类似。

推导：设 $\mathcal C$ 的秩为 $r$，取 $\mathcal C$ 中 $r$ 个线性无关列并补足为基构成变换矩阵 $T^{-1}$。在此基下，$A$ 作用保持能控子空间不变，故分块下三角；$B$ 的相应分量为零。

物理/几何意义：系统可分解为四个子空间：能控能观、能控不能观、不能控能观、不能控不能观。

极限/特例：系统最小实现的阶数等于能控能观子空间维数。

典型应用：模型降阶、最小实现。

**知识点7.5　最小实现**

定义：给定传递函数 $G(s)$，阶数最小且能控能观的实现称为最小实现。

假设：$G(s)$ 为真有理函数。

公式：最小实现阶数等于 $G(s)$ 的麦克米伦阶（SISO 时为分子分母互质后的分母阶数）。

推导：任一实现经卡尔曼分解后，只有能控能观子空间对 $G(s)$ 有贡献，且 $G(s)=C_1(sI-A_{11})^{-1}B_1$。因此最小阶数等于能控能观子空间维数。又 $G(s)$ 的互质分母次数是该子空间的维数下界，故得结论。

物理/几何意义：消去传递函数中的零极点对消，得到系统“本质”的动态阶数。

极限/特例：$G(s)$ 有零极点对消时，非最小实现阶数大于最小阶数。

典型应用：系统辨识、控制器降阶。

**知识点7.6　能观性与卡尔曼滤波的可观测条件**

定义：随机系统中，能观性是滤波估计误差协方差有界、滤波器可收敛的必要条件。

假设：线性高斯系统，过程噪声与量测噪声统计已知。

公式：能观性矩阵满秩是滤波误差协方差黎卡提方程存在正定稳态解的条件之一。

推导：卡尔曼滤波误差协方差满足黎卡提方程；若存在不可观且不稳定的模态，该模态误差协方差发散。对不可观模态，量测不含其信息，滤波增益对该方向为零，故其估计误差仅由过程噪声驱动，若该模态不稳定则发散。

物理/几何意义：能观性保证量测足以抑制所有不稳定模态的估计误差。

极限/特例：可检测（不可观模态稳定）时滤波器仍可收敛。

典型应用：导航、目标跟踪。

### 第8章 李雅普诺夫稳定性理论

**知识点8.1　稳定性的定义（李雅普诺夫意义）**

定义：平衡点 $x_e=0$ 稳定指对任意 $\varepsilon>0$ 存在 $\delta>0$，使 $\|x(0)\|<\delta\Rightarrow\|x(t)\|<\varepsilon,\forall t\ge0$；渐近稳定指稳定且存在 $\delta$ 使 $\|x(0)\|<\delta\Rightarrow x(t)\to0$；指数稳定指 $\|x(t)\|\le\alpha\|x(0)\|e^{-\beta t}$。

假设：$\dot x=f(x)$，$f(0)=0$，$f$ 局部利普希茨。

公式：李雅普诺夫稳定 $\iff\forall\varepsilon>0,\exists\delta>0:\|x_0\|<\delta\Rightarrow\sup_{t\ge0}\|x(t)\|<\varepsilon$。

推导：由解对初值的连续依赖性结合一致有界性可证稳定性等价于轨线不逃逸出任意给定邻域。渐近稳定进一步要求吸引性。

物理/几何意义：稳定性是不“发散”，渐近稳定是最终回到平衡点。

极限/特例：稳定但不渐近（如无阻尼振子 $\ddot x+x=0$）；指数稳定 $\Rightarrow$ 渐近稳定。

典型应用：控制系统的首要性能指标。

**知识点8.2　李雅普诺夫第二法（直接法）**

定义：不求解微分方程，通过构造标量函数 $V(x)$ 的定号性判别稳定性。

假设：$V:\mathcal D\to\mathbb R$ 连续可微，$V(0)=0$。

公式：
- $V>0$，$\dot V\le0$ $\Rightarrow$ 稳定；
- $V>0$，$\dot V<0$ $\Rightarrow$ 渐近稳定；
- $V>0$，$\dot V<0$，且 $V$ 径向无界 $\Rightarrow$ 全局渐近稳定。

推导（渐近稳定）：$\dot V=\nabla V^Tf<0$ 使 $V$ 沿轨线严格递减，故 $V(x(t))<V(x_0)$。取 $c>0$ 使 $\Omega_c=\{V(x)\le c\}$ 含于吸引域。对任意 $x_0\in\Omega_c$，轨线始终在 $\Omega_c$ 内且 $V$ 单调下降有下界 $0$，故收敛。设极限 $V^\ast>0$，则轨线趋于紧集 $\{V=V^\ast\}$，其上 $\dot V<0$ 与 $V$ 收敛矛盾，故 $V^\ast=0$，即 $x(t)\to0$。

物理/几何意义：$V$ 是广义能量，$\dot V<0$ 表示能量持续耗散；$V$ 的等值面（李雅普诺夫面）被轨线由外向内穿越。

极限/特例：$\dot V\le0$（半负定）时只能判稳定，需 LaSalle 不变集原理加强为渐近稳定。

典型应用：非线性系统镇定、自适应控制稳定性分析。

**知识点8.3　线性系统的李雅普诺夫方程**

定义：对线性系统 $\dot x=Ax$，取二次型 $V=x^TPx$。

假设：$P=P^T>0$。

公式：
$$A^TP+PA=-Q,\qquad Q=Q^T>0.$$
系统渐近稳定 $\iff$ 对任意 $Q>0$ 存在唯一 $P>0$ 满足该方程。

推导：$\dot V=x^T(A^TP+PA)x=-x^TQx$。若 $P$ 存在，则 $\dot V<0$，由知识点8.2 得渐近稳定。反之，若 $A$ 稳定，定义
$$P=\int_0^\infty e^{A^Tt}Qe^{At}\,dt,$$
此积分收敛（因 $e^{At}$ 指数衰减），且
$$A^TP+PA=\int_0^\infty\frac{d}{dt}(e^{A^Tt}Qe^{At})dt=e^{A^Tt}Qe^{At}\Big|_0^\infty=-Q,$$
且 $P>0$。

物理/几何意义：$P$ 度量状态“能量”，$Q$ 度量耗散率；线性稳定性可完全由代数方程判定。

极限/特例：$A$ 有虚轴特征值时方程可能无正定解；$Q$ 可选为 $C^TC$ 用于能观性分析。

典型应用：二次型调节器、稳定性数值判据、可控 Gramian 求解。

**知识点8.4　克拉索夫斯基法与变量梯度法**

定义：克拉索夫斯基法取 $V=(f(x))^T P(f(x))$；变量梯度法直接构造梯度 $\nabla V$。

假设：$f$ 可微、$f(0)=0$。

公式（克拉索夫斯基）：若存在 $P>0$ 使雅可比 $J(x)=\partial f/\partial x$ 满足
$$J^TP+PJ<0,$$
则 $V=f^TPf$ 正定且 $\dot V=f^T(J^TP+PJ)f<0$，系统渐近稳定。

推导：$\dot f=J\dot x=Jf$，故
$$\dot V=\dot f^TPf+f^TP\dot f=f^TJ^TPf+f^TPJf=f^T(J^TP+PJ)f.$$
证毕。

变量梯度法：设 $\nabla V=g(x)$，要求 $\partial g_i/\partial x_j=\partial g_j/\partial x_i$（旋度为零）以保证 $V$ 存在，再积分得 $V$。

物理/几何意义：克拉索夫斯基法把非线性稳定性化归为雅可比矩阵的负定性。

极限/特例：$f=Ax$ 时即线性李雅普诺夫方程。

典型应用：非线性电路、非线性机械系统。

**知识点8.5　拉萨尔不变集原理**

定义：若 $V$ 正定、$\dot V\le0$，则轨线收敛到最大不变集 $\mathcal M\subset\{x:\dot V(x)=0\}$。

假设：$V$ 径向无界时全局成立；$\Omega$ 为有界闭不变集。

公式：
$$\mathcal M=\text{包含于}\ \{x:\dot V(x)=0\}\ \text{的最大不变集},\qquad x(t)\to\mathcal M.$$

推导：$V$ 有下界且单调不增，故收敛到 $L$。轨线有界，其 $\omega$ 极限集 $\Omega^+$ 非空紧不变。在 $\Omega^+$ 上 $V\equiv L$，故 $\dot V\equiv0$，即 $\Omega^+\subset\{x:\dot V=0\}$。由于 $\Omega^+$ 不变，$\Omega^+\subset\mathcal M$，故 $x(t)\to\mathcal M$。

物理/几何意义：即使 $\dot V$ 只在部分集合为零，轨线也只能收敛到其中不变的部分。

极限/特例：若 $\{x:\dot V=0\}$ 中除原点外无不变集，则原点渐近稳定。

典型应用：自适应控制、机器人无源性分析。

**知识点8.6　全局渐近稳定与吸引域估计**

定义：全局渐近稳定指平衡点稳定且吸引域为全空间；吸引域为所有使 $x(t)\to0$ 的初态集合。

假设：$V$ 径向无界、$\dot V<0$（全局）。

公式：$V$ 径向无界且 $\dot V<0,\ \forall x\neq0$ $\Rightarrow$ 全局渐近稳定。

推导：径向无界保证有界集中 $V$ 等值面所围区域有界，从而任意初态始终于某有界李雅普诺夫面内，轨线有界且 $V$ 单调下降趋于零，故收敛到原点。

吸引域估计：取最大 $c$ 使 $\Omega_c=\{V\le c\}$ 满足 $\dot V<0$（除原点），则 $\Omega_c\subset$ 吸引域。

物理/几何意义：径向无界保证能量面“包裹”整个状态空间，排除有界发散。

极限/特例：$V$ 非径向无界时可能出现初态逃逸，只能得局部结果。

典型应用：非线性控制器全局稳定证明。

**知识点8.7　指数稳定性与比较原理**

定义：指数稳定指存在 $\alpha,\beta>0$ 使 $\|x(t)\|\le\alpha\|x(0)\|e^{-\beta t}$。

假设：$\dot x=f(x,t)$，存在 $V$ 满足
$$c_1\|x\|^2\le V\le c_2\|x\|^2,\quad \dot V\le-c_3\|x\|^2.$$

公式：则 $\|x(t)\|\le\sqrt{c_2/c_1}\,\|x(0)\|e^{-c_3t/(2c_2)}$。

推导：由 $\dot V\le-\dfrac{c_3}{c_2}V$，得 $\dfrac{d}{dt}(Ve^{c_3t/c_2})\le0$，故 $V(t)\le V(0)e^{-c_3t/c_2}$。再用 $c_1\|x\|^2\le V$、$V(0)\le c_2\|x(0)\|^2$，开方得结论。

物理/几何意义：指数稳定给出收敛速率的定量估计。

极限/特例：线性稳定系统必指数稳定；非线性可能只有渐近（多项式）收敛。

典型应用：收敛速度估计、观测器性能分析。

### 第9章 状态反馈与极点配置

**知识点9.1　状态反馈的基本结构**

定义：状态反馈控制律 $u=-Kx+v$，$K\in\mathbb R^{m\times n}$，$v$ 为参考输入。

假设：状态可测；$(A,B)$ 能控或至少能稳。

公式：闭环系统
$$\dot x=(A-BK)x+Bv,\qquad y=(C-DK)x+Dv.$$

推导：把 $u=-Kx+v$ 代入 $\dot x=Ax+Bu$：
$$\dot x=Ax+B(-Kx+v)=(A-BK)x+Bv;$$
代入输出 $y=Cx+Du=(C-DK)x+Dv$。

物理/几何意义：状态反馈改变系统矩阵特征值（极点），从而改变动态响应；不改变输入输出零点结构（见知识点9.5）。

极限/特例：$K=0$ 为开环。

典型应用：全状态反馈控制器、LQR 增益。

**知识点9.2　极点配置定理**

定义：若 $(A,B)$ 能控，则可用 $u=-Kx$ 把闭环特征值配置到任意共轭对称集合。

假设：$(A,B)$ 能控；期望特征值可任选（复共轭成对）。

公式：
$$\det(\lambda I-A+BK)=\prod_{i=1}^n(\lambda-\lambda_i^\ast).$$

推导（单输入）：先化能控标准型 $(\bar A,\bar B)$。设 $\bar K=[k_0\ \cdots\ k_{n-1}]$，则
$$\bar A-\bar B\bar K=\begin{bmatrix}0&1&&\\ \vdots&&\ddots&\\ 0&&&1\\ -a_0-k_0&-a_1-k_1&\cdots&-a_{n-1}-k_{n-1}\end{bmatrix}.$$
其首一特征多项式为 $\lambda^n+(a_{n-1}+k_{n-1})\lambda^{n-1}+\cdots+(a_0+k_0)$。期望多项式 $\prod(\lambda-\lambda_i^\ast)=\lambda^n+\alpha_{n-1}\lambda^{n-1}+\cdots+\alpha_0$。逐项比较得 $k_i=\alpha_i-a_i$。再变换回原坐标 $K=\bar KT$。多输入时需先选输入方向化为单输入可配条件或解线性方程。

物理/几何意义：能控性保证输入能独立影响每个特征方向，从而自由配置模态。

极限/特例：不能控模态不可配；$n$ 个极点需 $n$ 个独立参数时需 $m\ge1$ 且能控。

典型应用：指定超调、调节时间设计控制器。

**知识点9.3　阿克曼公式**

定义：单输入系统极点配置增益的显式公式。

假设：$(A,B)$ 能控，期望特征多项式
$$\alpha(\lambda)=\lambda^n+\alpha_{n-1}\lambda^{n-1}+\cdots+\alpha_0.$$

公式：
$$K=[0\ \cdots\ 0\ 1]\,\mathcal C^{-1}\,\alpha(A),\qquad \alpha(A)=A^n+\alpha_{n-1}A^{n-1}+\cdots+\alpha_0I.$$

推导：在能控标准型下，$\bar K=[\alpha_0-a_0\ \cdots\ \alpha_{n-1}-a_{n-1}]$。由 $\bar A$ 的友矩阵性质，
$$\alpha(\bar A)=\bar A^n+\alpha_{n-1}\bar A^{n-1}+\cdots+\alpha_0I$$
的最后一行恰为 $[a_0-\alpha_0\ \cdots\ a_{n-1}-\alpha_{n-1}]$（因 $\bar A$ 满足 $\det(\lambda I-\bar A)=\lambda^n+a_{n-1}\lambda^{n-1}+\cdots+a_0$）。故
$$\bar K=-[0\ \cdots\ 1]\alpha(\bar A).$$
利用 $T=\mathcal C^{-1}\mathcal C_c$、$\bar K=KT^{-1}$ 及 $[0\cdots1]\mathcal C_c^{-1}\mathcal C=[0\cdots1]\mathcal C^{-1}$（当 $\mathcal C_c=I$ 时）得 $K=[0\cdots1]\mathcal C^{-1}\alpha(A)$。

物理/几何意义：阿克曼公式把极点配置化为矩阵多项式的计算，避免显式求变换矩阵。

极限/特例：仅适用单输入；多输入可用循环设计或直接法。

典型应用：MATLAB `acker`、`place`。

**知识点9.4　状态反馈对能控性与能观性的影响**

定义：状态反馈保持能控性不变，但可能破坏能观性。

假设：$u=-Kx+v$。

公式：$(A,B)$ 能控 $\iff(A-BK,B)$ 能控；能观性 $(C,A)$ 一般不保持。

推导（能控性）：若存在 $z$ 使 $z^T(A-BK)^k B=0$ 对所有 $k$，可反复展开证明 $z^TA=0$ 或 $z^TB=0$，与能控矛盾。更简洁地：$A-BK=A+B(-K)$，能控性在反馈 $u=-Kx+v$ 下是输入变换下的不变量，可用能控子空间的不变性证明。
能观性反例：$A=\begin{bmatrix}0&1\\0&0\end{bmatrix},B=\begin{bmatrix}0\\1\end{bmatrix},C=[1\ 0]$，取 $K=[1\ 1]$ 后能观性矩阵秩下降。

物理/几何意义：反馈只能影响输入通道，不能改变输入对状态的“可达范围”；但会改变零点与输出通道。

极限/特例：能控必然能稳（能配置不稳定模态）。

典型应用：观测器与控制器联合设计需检查能观性。

**知识点9.5　状态反馈与零点**

定义：状态反馈不改变系统的输入输出零点（传递函数零点），只改变极点。

假设：$(A,B,C)$ 可观测足够；零点由 Rosenbrock 系统矩阵定义。

公式：$\det\begin{bmatrix}sI-A&-B\\C&D\end{bmatrix}$ 的零点不因 $A\to A-BK$ 改变。

推导：考察闭环 Rosenbrock 矩阵
$$\begin{bmatrix}sI-A+BK&-B\\C-DK&D\end{bmatrix}.$$
做行变换 $R=[I\ \ -K]$ 左乘前块并对应列变换：
$$\begin{bmatrix}I&0\\-K&I\end{bmatrix}\begin{bmatrix}sI-A&-B\\C&D\end{bmatrix}=\begin{bmatrix}sI-A&-B\\C+KsI-KA&D-KB\end{bmatrix}.$$
进一步可化为与原系统零点相同的多项式（乘以可逆矩阵）。故不变零点保持。

物理/几何意义：反馈移动极点但不移动零点，零极点相对位置决定响应形状。

极限/特例：若零点与配置极点接近，会出现大的响应峰值。

典型应用：非最小相位系统限制控制带宽。

**知识点9.6　镇定与能稳性**

定义：镇定指设计反馈使闭环渐近稳定；能稳指不可控模态渐近稳定。

假设：$(A,B)$ 可能不完全能控。

公式：存在 $K$ 使 $\mathrm{Re}\lambda_i(A-BK)<0\iff$ 所有不可控模态稳定（Popov–Belevitch–Hautus：对 $\mathrm{Re}\lambda\ge0$，$\mathrm{rank}[\lambda I-A\ B]=n$）。

推导：由 PBH 判据，能控模态可任意配置；不可控模态特征值由 $\bar A_{22}$ 决定，反馈无法改变，故可镇定当且仅当 $\mathrm{Re}\lambda(\bar A_{22})<0$。

物理/几何意义：镇定只需求稳定，不需求性能，适用于不稳定但可稳的系统。

极限/特例：能控 $\Rightarrow$ 能稳；能稳不一定能控。

典型应用：不稳定飞行器增稳。

**知识点9.7　参考输入跟踪与伺服补偿**

定义：为使输出跟踪常值参考 $r$，需引入参考增益或积分器。

假设：系统能控能观、参考为常值或阶跃。

公式：$u=-Kx+Nr$，其中 $N$ 满足稳态
$$y_\infty=C(-A+BK)^{-1}BNr=r.$$
引入积分状态 $\dot z=r-y$ 后增广系统
$$\begin{bmatrix}\dot x\\\dot z\end{bmatrix}=\begin{bmatrix}A&0\\-C&0\end{bmatrix}\begin{bmatrix}x\\ z\end{bmatrix}+\begin{bmatrix}B\\0\end{bmatrix}u+\begin{bmatrix}0\\1\end{bmatrix}r.$$

推导：稳态时 $\dot x=0$，$x_\infty=-(A-BK)^{-1}BNr$，输出 $y_\infty=Cx_\infty$，令 $y_\infty=r$ 解 $N$。积分器使误差 $e=r-y$ 趋于零：由 $\dot z=e$，$z$ 为误差积分，闭环稳定时 $e\to0$。

物理/几何意义：积分作用消除稳态误差，实现无差跟踪。

极限/特例：型别系统需更高阶积分器跟踪斜坡。

典型应用：伺服电机、飞行器航向跟踪。

### 第10章 状态观测器

**知识点10.1　全维状态观测器**

定义：用可测输入输出重构状态，形如
$$\dot{\hat x}=A\hat x+Bu+L(y-C\hat x),$$
$L$ 为观测器增益。

假设：$(A,C)$ 能观或能检测。

公式：估计误差 $e=x-\hat x$ 满足
$$\dot e=(A-LC)e.$$

推导：$\dot x=Ax+Bu$，$\dot{\hat x}=A\hat x+Bu+L(Cx-C\hat x)=A\hat x+Bu+LCe$。相减
$$\dot e=\dot x-\dot{\hat x}=Ae-LCe=(A-LC)e.$$

物理/几何意义：观测器复制系统动态，并用输出误差校正；误差动态由 $A-LC$ 决定。

极限/特例：$L=0$ 为开环观测器（对初值误差和扰动敏感）；$A-LC$ 稳定则 $\hat x\to x$。

典型应用：状态反馈需要但状态不可测时。

**知识点10.2　观测器极点配置与对偶性**

定义：通过选择 $L$ 配置 $A-LC$ 特征值。

假设：$(A,C)$ 能观。

公式：$A-LC$ 与 $(A^T-C^TL^T)$ 特征值相同；对偶地令 $K=L^T$ 配置 $(A^T,C^T)$ 的反馈极点。

推导：$\det(\lambda I-A+LC)=\det(\lambda I-A^T+C^TL^T)$（转置不改变行列式）。故观测器极点配置等价于对偶系统 $(A^T,C^T)$ 的状态反馈极点配置，由知识点9.2 能观保证可任意配置。增益 $L=K^T$（对偶增益）。

物理/几何意义：观测器设计是控制器设计的对偶问题。

极限/特例：能观可任意配置；能检测可配置稳定。

典型应用：龙伯格观测器、`place` 对偶调用。

**知识点10.3　降维观测器**

定义：若输出直接给出部分状态，可只估计剩余状态，降低观测器阶数。

假设：$C$ 满行秩，不妨设 $C=[I\ \ 0]$，状态分块 $x=[x_1;x_2]$。

公式：
$$\begin{aligned}\dot x_1&=A_{11}x_1+A_{12}x_2+B_1u\\ \dot x_2&=A_{21}x_1+A_{22}x_2+B_2u\\ y&=x_1.\end{aligned}$$
令 $z=x_2-Ly$，则
$$\dot z=(A_{22}-LA_{12})z+[(A_{21}-LA_{11})+(A_{22}-LA_{12})L]y+(B_2-LB_1)u.$$

推导：把 $x_2=z+Ly$ 代入 $\dot x_2$：
$$\dot z=\dot x_2-L\dot y=A_{21}y+A_{22}(z+Ly)+B_2u-L(A_{11}y+A_{12}(z+Ly)+B_1u).$$
整理得所示方程；误差 $e=x_2-\hat x_2$ 满足 $\dot e=(A_{22}-LA_{12})e$，配置其极点使 $e\to0$，由 $\hat x_2=\hat z+Ly$ 重构。

物理/几何意义：降维观测器阶数为 $n-\mathrm{rank}C$，结构更经济。

极限/特例：$C=I$ 时无需观测器。

典型应用：仅部分状态可测的机械系统。

**知识点10.4　分离原理**

定义：状态反馈增益 $K$ 与观测器增益 $L$ 可独立设计，组合系统稳定性由 $A-BK$ 与 $A-LC$ 共同决定。

假设：$(A,B)$ 能控、$(A,C)$ 能观。

公式：用 $\hat x$ 代替 $x$ 的控制律 $u=-K\hat x$，闭环系统
$$\begin{bmatrix}\dot x\\\dot e\end{bmatrix}=\begin{bmatrix}A-BK&BK\\0&A-LC\end{bmatrix}\begin{bmatrix}x\\ e\end{bmatrix},$$
特征值为 $\lambda(A-BK)\cup\lambda(A-LC)$。

推导：$u=-Kx+Ke$，故 $\dot x=Ax-BKx+BKe=(A-BK)x+BKe$，且 $\dot e=(A-LC)e$。矩阵为分块上三角，特征值即对角块特征值之并。

物理/几何意义：控制器与观测器设计解耦，可独立配置。

极限/特例：理论上成立；实际采样/噪声下需权衡观测器速度。

典型应用：LQG、输出反馈控制器。

**知识点10.5　龙伯格观测器与最优观测器**

定义：龙伯格观测器即全维或降维观测器；卡尔曼滤波是噪声下的最优观测器。

假设：噪声统计已知（$Q,R$）。

公式：卡尔曼增益
$$L=P C^TR^{-1},\qquad \dot P=AP+PA^T+Q-PC^TR^{-1}CP.$$

推导：由正交投影原理或最小化误差协方差 $\mathrm{tr}(P)$，对增益求导 $\partial\mathrm{tr}(P)/\partial L=0$ 得 $L=PC^TR^{-1}$，代入协方差方程得黎卡提方程。详见第14章。

物理/几何意义：龙伯格观测器按极点速度设计，卡尔曼滤波按噪声权衡最优。

极限/特例：$Q\to0$ 时卡尔曼滤波增益变小，趋向开环观测。

典型应用：导航滤波、定位。

**知识点10.6　观测器与干扰观测器**

定义：干扰观测器估计并补偿未知外部干扰或模型不确定性。

假设：干扰慢变或可由外系统建模。

公式：对 $\dot x=Ax+Bu+B_dd$，扩张干扰为状态 $\dot d=0$，构造
$$\dot{\hat x}=A\hat x+Bu+L(y-C\hat x),\quad \dot{\hat d}=L_d(y-C\hat x).$$

推导：增广系统状态 $[x;d]$，$\begin{bmatrix}\dot x\\\dot d\end{bmatrix}=\begin{bmatrix}A&B_d\\0&0\end{bmatrix}\begin{bmatrix}x\\ d\end{bmatrix}+\begin{bmatrix}B\\0\end{bmatrix}u$，$y=[C\ 0][x;d]$。对增广系统设计观测器得扩张状态观测器（ESO）。误差动态由增广系统能观性决定。

物理/几何意义：把干扰视为扩张状态，估计后前馈补偿。

极限/特例：干扰快变时零阶扩张模型失配；可提高阶数。

典型应用：自抗扰控制（ADRC）、运动控制。

### 第11章 最优控制：变分法基础

**知识点11.1　泛函、变分与一阶变分**

定义：泛函 $J[x]=\int_{t_0}^{t_f}L(x,\dot x,t)dt$；变分 $\delta x$ 是函数的无穷小改变；一阶变分 $\delta J$ 是 $J$ 对 $\delta x$ 的线性主部。

假设：$L$ 二次连续可微；端点或自由或固定。

公式：
$$\delta J=\int_{t_0}^{t_f}\left(\frac{\partial L}{\partial x}\delta x+\frac{\partial L}{\partial\dot x}\delta\dot x\right)dt+[\text{端点项}].$$

推导：$J[x+\epsilon\eta]=\int L(x+\epsilon\eta,\dot x+\epsilon\dot\eta,t)dt$，对 $\epsilon$ 求导
$$\left.\frac{d}{d\epsilon}J\right|_{\epsilon=0}=\int\left(\frac{\partial L}{\partial x}\eta+\frac{\partial L}{\partial\dot x}\dot\eta\right)dt.$$
令 $\eta=\delta x$ 即得。

物理/几何意义：变分是函数空间中的方向导数；泛函驻值点是“最优轨线”。

极限/特例：$L$ 不含 $\dot x$ 时可直接逐点极小化。

典型应用：最速降线、最优控制。

**知识点11.2　欧拉–拉格朗日方程（变分法）**

定义：使泛函取极值的函数满足的微分方程。

假设：$L$ 光滑、端点固定 $\delta x(t_0)=\delta x(t_f)=0$。

公式：
$$\frac{\partial L}{\partial x}-\frac{d}{dt}\frac{\partial L}{\partial\dot x}=0.$$

推导：对上式分部积分
$$\int_{t_0}^{t_f}\frac{\partial L}{\partial\dot x}\delta\dot x\,dt=\frac{\partial L}{\partial\dot x}\delta x\Big|_{t_0}^{t_f}-\int_{t_0}^{t_f}\frac{d}{dt}\frac{\partial L}{\partial\dot x}\delta x\,dt,$$
端点固定使边界项为零。故
$$\delta J=\int_{t_0}^{t_f}\left(\frac{\partial L}{\partial x}-\frac{d}{dt}\frac{\partial L}{\partial\dot x}\right)\delta x\,dt=0.$$
由 $\delta x$ 任意且连续，括号恒为零（变分法基本引理）。

物理/几何意义：欧拉–拉格朗日方程是泛函极值的必要条件，与力学拉格朗日方程形式一致（这里 $L$ 是代价密度）。

极限/特例：若 $L$ 不显含 $t$，则 $L-\dot x\,\partial L/\partial\dot x=\mathrm{const}$（贝尔特拉米恒等式）。

典型应用：经典变分问题。

**知识点11.3　横截条件**

定义：当端点自由或终端时刻自由时，边界项给出额外的横截条件。

假设：$t_f$ 自由、$x(t_f)$ 自由。

公式：
$$\left(\frac{\partial L}{\partial\dot x}\delta x+L\,\delta t\right)\Big|_{t_0}^{t_f}=0.$$
常见情形：
$$\frac{\partial L}{\partial\dot x}\Big|_{t_f}=0\ (x(t_f)\ \text{自由}),\qquad \left(L-\dot x\frac{\partial L}{\partial\dot x}\right)\Big|_{t_f}=0\ (t_f\ \text{自由}).$$

推导：若 $t_f$ 自由，终端状态变化 $\delta x_f+\dot x(t_f)\delta t_f$，故 $\delta x(t_f)$ 的独立变分为 $\delta x_f+\dot x\delta t_f$。边界项
$$\frac{\partial L}{\partial\dot x}\delta x\Big|_{t_f}+L\delta t_f=\frac{\partial L}{\partial\dot x}\delta x_f+\left(L+\frac{\partial L}{\partial\dot x}\dot x\right)\delta t_f$$
（变分法中对 $t_f$ 的贡献需仔细处理，标准结果为 $[L\,\delta t+\partial L/\partial\dot x(\delta x-\dot x\delta t)]_{t_f}=0$）。分别令 $\delta x_f=0$ 或 $\delta t_f=0$ 得上述条件。

物理/几何意义：自由端点对应自然边界条件，决定最优轨线如何“击中”目标集。

极限/特例：终点固定在流形 $N(x,t)=0$ 上时用拉格朗日乘子。

典型应用：拦截问题、最优落点。

**知识点11.4　拉格朗日问题、梅耶尔问题与博尔扎问题**

定义：拉格朗日问题 $J=\int L\,dt$；梅耶尔问题 $J=\Phi(x(t_f),t_f)$；博尔扎问题 $J=\Phi+\int L\,dt$（最一般）。

假设：终点条件 $x(t_f,t_f)\in S$。

公式：博尔扎问题通过引入附加状态 $x_{n+1}$，$\dot x_{n+1}=L$，$x_{n+1}(t_0)=0$，化为梅耶尔问题。

推导：令 $\dot x_{n+1}=L(x,u,t)$，则 $x_{n+1}(t_f)=\int L\,dt$，代价变为 $\Phi(x(t_f))+x_{n+1}(t_f)$，即纯终端代价。故三类问题可互相转换。

物理/几何意义：终端代价与运行代价的统一描述。

极限/特例：$\Phi=0$ 为拉格朗日；$L=0$ 为梅耶尔。

典型应用：最优控制通用表述。

**知识点11.5　哈密顿形式与协态方程**

定义：用哈密顿函数把变分问题化为一阶系统。

假设：$L$ 可微，控制问题 $\dot x=f(x,u,t)$，代价 $\int L$。

公式：
$$\mathcal H=L+\lambda^Tf,\qquad \dot x=\frac{\partial\mathcal H}{\partial\lambda},\quad \dot\lambda=-\frac{\partial\mathcal H}{\partial x},\quad \frac{\partial\mathcal H}{\partial u}=0.$$

推导：把动态作为约束引入拉格朗日乘子 $\lambda$：
$$J_a=\int\left[L+\lambda^T(f-\dot x)\right]dt+\Phi.$$
分部积分 $\int\lambda^T\dot x\,dt=\lambda^T\delta x|-\int\dot\lambda^T\delta x\,dt$，对 $x$ 变分得协态方程 $\dot\lambda=-\partial L/\partial x-(\partial f/\partial x)^T\lambda$，对 $u$ 变分得 $\partial L/\partial u+\lambda^T\partial f/\partial u=0$，即 $\partial\mathcal H/\partial u=0$。

物理/几何意义：$\lambda$ 是状态约束的影子价格（边际价值），$\mathcal H$ 是瞬时“扩展代价”。

极限/特例：$\mathcal H$ 不显含 $t$ 时 $\mathcal H=$ 常数。

典型应用：最优轨线数值求解（打靶法）。

**知识点11.6　最小作用量原理与变分法联系**

定义：哈密顿最小作用量原理：真实运动使作用量 $S=\int L\,dt$ 取驻值。

假设：保守系统、完整约束。

公式：$\delta S=0$ 给出 $\dfrac{d}{dt}\dfrac{\partial L}{\partial\dot q}-\dfrac{\partial L}{\partial q}=0$。

推导：与知识点11.2 完全相同，只把 $x$ 换成广义坐标 $q$。作用量的变分为零即欧拉–拉格朗日方程。

物理/几何意义：自然界沿“最省作用量”路径演化，是力学与分析力学的统一原理。

极限/特例：相对论作用量含 $-mc^2\sqrt{1-v^2/c^2}$；量子力学路径积分以作用量为相。

典型应用：场论、最优控制与力学统一。

### 第12章 线性二次型最优控制（LQR）与 HJB

**知识点12.1　LQR 问题提法**

定义：线性系统 $\dot x=Ax+Bu$，二次型性能指标
$$J=\frac12\int_{t_0}^{t_f}(x^TQx+u^TRu)dt+\frac12x^T(t_f)Sx(t_f),$$
$Q\ge0$、$R>0$、$S\ge0$。求最优 $u^\ast$。

假设：$(A,B)$ 能稳，$Q$ 使 $(A,Q^{1/2})$ 能检测（无限时域）。

公式：目标是最小化 $J$。

推导：由二次型结构猜想值函数为二次型，见知识点12.2。有限时域用哈密顿形式；无限时域用 HJB。

物理/几何意义：$Q$ 惩罚状态偏差，$R$ 惩罚控制能量，权衡性能与代价。

极限/特例：$R\to0$ 得最小能量；$Q\to0$ 得最小控制。

典型应用：飞行器增稳、机器人轨迹跟踪。

**知识点12.2　黎卡提方程推导（HJB 方法）**

定义：值函数 $V(x,t)=\min_u J$ 满足 HJB 方程
$$-\frac{\partial V}{\partial t}=\min_u\left[\frac12x^TQx+\frac12u^TRu+\frac{\partial V}{\partial x}^T(Ax+Bu)\right].$$

假设：$V$ 二次可微、LQR 结构下 $V=\frac12x^TP(t)x$。

公式：最小化右端：对 $u$ 求导 $Ru+B^T\frac{\partial V}{\partial x}^T=0$，得 $u^\ast=-R^{-1}B^TP(t)x$。代入并令 $x$ 各二次项系数相等得
$$-\dot P=PA+A^TP-PBR^{-1}B^TP+Q,\qquad P(t_f)=S.$$

推导：设 $V=\frac12x^TPx$，$\partial V/\partial x=Px$。括号内函数
$$\phi(u)=\frac12x^TQx+\frac12u^TRu+x^TP(Ax+Bu).$$
$\partial\phi/\partial u=Ru+B^TPx=0\Rightarrow u^\ast=-R^{-1}B^TPx$。代回：
$$\phi(u^\ast)=\frac12x^TQx+x^TPAx-\frac12x^TPBR^{-1}B^TPx.$$
HJB：$-\frac12x^T\dot Px=\phi(u^\ast)$，对任意 $x$ 成立，比较二次型
$$-\dot P=PA+A^TP-PBR^{-1}B^TP+Q.$$
终端 $V(x,t_f)=\frac12x^TSx$ 给 $P(t_f)=S$。

物理/几何意义：$P(t)$ 是未来代价的“价值矩阵”，非线性矩阵微分方程沿时间反向积分。

极限/特例：无限时域 $t_f\to\infty$、$A-BR^{-1}B^TP$ 稳定时 $\dot P=0$，得代数黎卡提方程。

典型应用：LQR 控制器设计。

**知识点12.3　代数黎卡提方程与最优反馈**

定义：无限时域 LQR 满足代数黎卡提方程（ARE）。

假设：$(A,B)$ 能稳、$(A,Q^{1/2})$ 能检测。

公式：
$$A^TP+PA-PBR^{-1}B^TP+Q=0,\qquad u^\ast=-R^{-1}B^TPx=-Kx.$$

推导：令有限时域 $\dot P=0$ 得 ARE。存在唯一正定稳定解 $P$，使 $A-BR^{-1}B^TP$ 稳定。

物理/几何意义：$K=R^{-1}B^TP$ 即最优状态反馈增益，兼顾性能与代价。

极限/特例：$Q=C^TC$ 时 LQR 与频域回路传函满足 Kalman 不等式 $|1+C(j\omega I-A)^{-1}BK|\ge1$。

典型应用：工业 MPC、飞行控制。

**知识点12.4　LQR 的稳定性**

定义：LQR 闭环自动稳定（在能稳能检测条件下）。

假设：$Q\ge0$、$R>0$、能稳能检测。

公式：$A_{cl}=A-BR^{-1}B^TP$ 稳定。

推导：取 $V=x^TPx$。由 ARE，
$$A_{cl}^TP+PA_{cl}=(A-BR^{-1}B^TP)^TP+P(A-BR^{-1}B^TP)
=A^TP+PA-2PBR^{-1}B^TP=-Q-PBR^{-1}B^TP<0.$$
（对任意 $x\neq0$，$x^T(Q+PBR^{-1}B^TP)x>0$ 由 $Q\ge0$ 与能检测保证。）故 $V$ 是李雅普诺夫函数，$A_{cl}$ 稳定。

物理/几何意义：最优性蕴含稳定性，无需额外设计。

极限/特例：$Q=0$ 时仅当 $(A,B)$ 能控且 $R>0$ 才保证；否则 ARE 可能退化解。

典型应用：稳定性保证的控制器综合。

**知识点12.5　加权矩阵选择与性能权衡**

定义：$Q$、$R$ 的选取决定响应速度与控制幅度。

假设：一般为对角阵，$Q=\mathrm{diag}(q_i)$，$R=\mathrm{diag}(r_j)$。

公式：增大 $Q$ 相对 $R$ 提高状态收敛速度、增大控制增益；反之节省控制。

推导（标量情形）：一阶系统 $\dot x=ax+bu$，$J=\int(qx^2+ru^2)dt$。ARE $2ap-b^2p^2/r+q=0$，解 $p=\frac{ar+\sqrt{a^2r^2+qb^2r}}{b^2}$，$K=pb/r$。可见 $q/r$ 增大时 $K$ 增大。

物理/几何意义：$Q/R$ 是性能与代价的权衡旋钮。

极限/特例：Bryson 规则 $Q=\mathrm{diag}(1/x_{i,\max}^2)$、$R=\mathrm{diag}(1/u_{j,\max}^2)$ 做量纲归一。

典型应用：工程调参。

**知识点12.6　哈密顿–雅可比–贝尔曼方程（一般形式）**

定义：连续时间最优控制的值函数满足的一阶偏微分方程。

假设：系统 $\dot x=f(x,u,t)$，代价 $\int L+\Phi(x(t_f))$，值函数可微。

公式：
$$-\frac{\partial V}{\partial t}=\min_{u}\left[L(x,u,t)+\nabla_xV^Tf(x,u,t)\right],\quad V(x,t_f)=\Phi(x).$$
最优控制 $u^\ast=\arg\min_u[\cdots]$。

推导：由最优性原理，$V(x,t)=\min_u\{\int_t^{t+\Delta}L\,d\tau+V(x(t+\Delta),t+\Delta)\}$。对 $\Delta$ 展开 $V(x(t+\Delta),t+\Delta)\approx V+V_t\Delta+\nabla V^Tf\Delta$，整理令 $\Delta\to0$：
$$0=\min_u[L+\nabla V^Tf+V_t],$$
即 HJB。

物理/几何意义：HJB 是最优控制的“主方程”，值函数是问题的完整解。

极限/特例：LQR 时 HJB 化为黎卡提方程；确定性 HJB 与经典哈密顿–雅可比方程同源。

典型应用：非线性最优控制、微分对策。

### 第13章 庞特里亚金极小值原理

**知识点13.1　极小值原理的陈述**

定义：对系统 $\dot x=f(x,u,t)$、代价 $J=\Phi(x(t_f),t_f)+\int_{t_0}^{t_f}L(x,u,t)dt$、控制约束 $u\in U$，最优控制使哈密顿量逐点取极小。

假设：$f,L,\Phi$ 对 $x$ 连续可微，对 $u$ 连续；$U$ 为闭凸集（或一般集合）。

公式：定义
$$\mathcal H(x,u,\lambda,t)=L(x,u,t)+\lambda^Tf(x,u,t),$$
则最优轨线满足
$$\dot x=\frac{\partial\mathcal H}{\partial\lambda},\qquad \dot\lambda=-\frac{\partial\mathcal H}{\partial x},$$
$$u^\ast(t)=\arg\min_{u\in U}\mathcal H(x^\ast(t),u,\lambda(t),t),$$

横截条件
$$\lambda(t_f)=\frac{\partial\Phi}{\partial x}\Big|_{t_f},\qquad \mathcal H(x^\ast,u^\ast,\lambda,t)\Big|_{t_f}=-\frac{\partial\Phi}{\partial t}\Big|_{t_f}.$$

推导（思路）：在控制变分 $\delta u$ 下，状态变分满足 $\delta\dot x=f_x\delta x+f_u\delta u$。代价一阶变分
$$\delta J=\frac{\partial\Phi}{\partial x}\delta x(t_f)+\int(L_x\delta x+L_u\delta u)dt.$$
引入协态 $\lambda$ 满足 $\dot\lambda=-f_x^T\lambda-L_x^T$、$\lambda(t_f)=\Phi_x^T$，则
$$\delta J=\int(\mathcal H_u)\delta u\,dt+\lambda^T\delta x\Big|_{t_0},$$
端点固定消去边界项。最优要求 $\delta J\ge0$ 对任意容许 $\delta u$。由于 $\delta u$ 受约束不能任意变号，只能要求 $\mathcal H$ 在 $u^\ast$ 处对 $U$ 取极小。

物理/几何意义：协态是约束的价格；PMP 是点wise 最优性条件，$U$ 有界时给出 bang-bang 结构。

极限/特例：$U=\mathbb R^m$ 时 $\partial\mathcal H/\partial u=0$，退化为变分法。

典型应用：时间最优、燃料最优控制。

**知识点13.2　哈密顿量的性质与守恒**

定义：$\mathcal H$ 沿最优轨线的取值性质。

假设：$f,L$ 不显含 $t$、$t_f$ 自由、$U$ 给定。

公式：沿最优轨线
$$\frac{d\mathcal H}{dt}=\frac{\partial\mathcal H}{\partial t},\qquad \text{若}\ \frac{\partial\mathcal H}{\partial t}=0\ \text{则}\ \mathcal H=\mathrm{const}.$$

推导：
$$\frac{d\mathcal H}{dt}=\mathcal H_x^T\dot x+\mathcal H_u^T\dot u+\mathcal H_\lambda^T\dot\lambda+\frac{\partial\mathcal H}{\partial t}.$$
对最优控制 $\mathcal H_u=0$（无约束）或 $\mathcal H$ 关于 $u$ 的极小值不随 $t$ 变化；$\mathcal H_x=-\dot\lambda$、$\mathcal H_\lambda=\dot x$，故前两项相消：
$$\mathcal H_x^T\dot x+\mathcal H_\lambda^T\dot\lambda=-\dot\lambda^T\dot x+\dot x^T\dot\lambda=0,$$
余下 $\partial\mathcal H/\partial t$。

物理/几何意义：$\mathcal H$ 是“最优控制系统的能量”，时间平移对称给出守恒。

极限/特例：终端时刻自由且不显含 $t$ 时 $\mathcal H=0$（与横截条件一致）。

典型应用：时间最优问题的求解。

**知识点13.3　时间最优控制与 bang-bang 原理**

定义：时间最优问题 $J=\int_{t_0}^{t_f}1\,dt$，控制约束 $|u_j|\le1$。

假设：系统可控、$f$ 线性于 $u$。

公式：对 $\dot x=Ax+Bu$，
$$\mathcal H=1+\lambda^T(Ax+Bu)=1+\lambda^TAx+\sum_j\lambda^TB_ju_j.$$
$u_j^\ast=-\mathrm{sgn}(\lambda^TB_j)$，即每个控制分量取边界值，切换由 $\lambda^TB_j$ 的符号决定。

推导：$\mathcal H$ 对 $u_j$ 是线性的，$u_j\in[-1,1]$，故最小值在端点取得：若系数 $\lambda^TB_j>0$ 取 $u_j=-1$，反之为 $+1$，即 $u_j^\ast=-\mathrm{sgn}(\lambda^TB_j)$。当系数为零时出现奇异弧。

物理/几何意义：最优控制始终“满舵”，只在有限个时刻切换，故为 bang-bang 控制。

极限/特例：无约束时可得连续控制；有约束时 bang-bang。

典型应用：快速定位、姿态机动。

**知识点13.4　奇异控制与切换曲面**

定义：若在某时间区间上 $\lambda^TB_j\equiv0$，则 PMP 不能直接给出 $u$，称奇异弧。

假设：奇异弧上 $\mathcal H_u\equiv0$ 且其各阶导数均为零。

公式：由 $\mathcal H_u=\lambda^TB=0$ 及
$$\frac{d}{dt}(\lambda^TB)=\dot\lambda^TB+\lambda^T B\dot{} =(-A^T\lambda-\cdots)^TB+\lambda^T\dot B=0,$$
逐步求导可确定奇异控制。

推导：沿奇异弧，$\lambda^TB=0$ 对 $t$ 反复求导并代入协态方程与状态方程，得到关于 $u$ 的代数方程，解出奇异控制 $u_s$。

物理/几何意义：奇异弧对应控制“恰好平衡”的退化情形，常出现在最优曲线中段。

极限/特例：线性系统时间最优一般无奇异；非线性或能耗型问题常有奇异弧。

典型应用：双积分器问题、机械臂能耗最优。

**知识点13.5　横截条件与自由终端**

定义：终端状态或时刻自由时补充的边界条件。

假设：终端约束 $N(x(t_f),t_f)=0$。

公式：
$$\lambda(t_f)=\frac{\partial\Phi}{\partial x}+\frac{\partial N}{\partial x}\nu,\qquad \mathcal H(t_f)=-\frac{\partial\Phi}{\partial t}-\frac{\partial N}{\partial t}\nu,$$
$\nu$ 为乘子。

推导：把约束经乘子 $\nu$ 加入终端代价 $\Phi+\nu^TN$，对 $x(t_f)$ 变分得第一式；对 $t_f$ 变分得第二式。

物理/几何意义：约束的梯度方向给出协态终端的“法向”分量。

极限/特例：$x(t_f)$ 自由时 $\nu=0$，$\lambda(t_f)=\Phi_x$；$t_f$ 自由时 $\mathcal H(t_f)=-\Phi_t$。

典型应用：最优拦截、轨道转移。

**知识点13.6　极小值原理与 HJB 的关系**

定义：PMP 是 HJB 沿最优轨线的必要条件，HJB 是全局充要条件（可微时）。

假设：值函数 $V$ 可微。

公式：沿最优轨线 $\lambda(t)=\nabla_xV(x^\ast(t),t)$，且
$$\mathcal H(x^\ast,u^\ast,\lambda,t)=-\frac{\partial V}{\partial t}.$$

推导：由 HJB，$u^\ast=\arg\min_u[L+\nabla V^Tf]$，即 $u^\ast$ 使 $\mathcal H(x,u,\nabla V,t)$ 极小，故 $\nabla V$ 扮演协态角色。沿最优轨线对 $x$ 求梯度可得 $\dot\lambda=-\mathcal H_x$。

物理/几何意义：协态是值函数的梯度（边际价值），PMP 是最优性沿单条轨线的局部表述，HJB 是全空间的偏微分方程。

极限/特例：$V$ 不可微时 PMP 仍可用，HJB 需粘性解。

典型应用：数值最优控制算法设计。

### 第14章 卡尔曼滤波与随机系统

**知识点14.1　随机过程与高斯白噪声**

定义：白噪声 $w(t)$ 均值为零、自相关 $E[w(t)w^T(\tau)]=Q(t)\delta(t-\tau)$；高斯过程由均值与协方差完全刻画。

假设：过程噪声与量测噪声为互不相关的高斯白噪声。

公式：随机系统
$$\dot x=Ax+Bu+Gw,\qquad y=Cx+v,$$
$$E[w]=0,\ E[v]=0,\ E[ww^T]=Q\delta,\ E[vv^T]=R\delta,\ E[wv^T]=0.$$

推导（白噪声作为形式导数）：维纳过程 $\beta(t)$ 满足 $E[d\beta d\beta^T]=Qdt$，形式导数 $w=d\beta/dt$ 即为白噪声。随机微分方程 $dx=Axdt+Bu dt+Gd\beta$ 的均方意义解即上述模型。

物理/几何意义：白噪声是理想化的无记忆扰动，用于建模随机力、量测误差。

极限/特例：有色噪声可用成形滤波器白化。

典型应用：导航、控制中的随机建模。

**知识点14.2　卡尔曼滤波递推公式**

定义：线性高斯系统的最小方差状态估计。

假设：初始状态 $x_0\sim\mathcal N(\hat x_0,P_0)$，噪声统计已知。

公式（连续时间）：
$$\dot{\hat x}=A\hat x+Bu+L(y-C\hat x),\qquad L=PC^TR^{-1},$$
$$\dot P=AP+PA^T+GQG^T-PC^TR^{-1}CP.$$

推导（正交投影法）：设估计为 $\hat x$，要求误差 $e=x-\hat x$ 与量测正交。由 $E[(x-\hat x)y^T]=0$ 与先验\更新结构，可推出上述增益与协方差方程。另一途径：由黎卡提方程 $\dot P=AP+PA^T+GQG^T-PC^TR^{-1}CP$，$L=PC^TR^{-1}$ 使 $\mathrm{tr}\dot P$ 最小。

物理/几何意义：卡尔曼滤波在“预测”（$A$ 演化，协方差增大）与“校正”（用量测 $y$）之间最优权衡。

极限/特例：$R\to0$ 完全信任量测，$L$ 大；$Q\to0$ 信任模型，$L$ 小。

典型应用：GPS/INS 组合导航、目标跟踪。

**知识点14.3　离散卡尔曼滤波**

定义：离散系统的递推估计器。

假设：$x_{k+1}=A_dx_k+B_du_k+w_k$，$y_k=C_dx_k+v_k$，$w_k\sim\mathcal N(0,Q_k)$，$v_k\sim\mathcal N(0,R_k)$。

公式：
- 预测：$\hat x_{k|k-1}=A_d\hat x_{k-1|k-1}+B_du_{k-1}$，$P_{k|k-1}=A_dP_{k-1|k-1}A_d^T+Q_{k-1}$；
- 增益：$K_k=P_{k|k-1}C_d^T(C_dP_{k|k-1}C_d^T+R_k)^{-1}$；
- 更新：$\hat x_{k|k}=\hat x_{k|k-1}+K_k(y_k-C_d\hat x_{k|k-1})$，$P_{k|k}=(I-K_kC_d)P_{k|k-1}$。

推导：设后验是高斯，条件密度 $p(x_k|Y_k)\propto p(y_k|x_k)p(x_k|Y_{k-1})$，两高斯相乘仍为高斯。比较指数二次项系数得 $P_{k|k}$ 与 $K_k$；线性项得均值更新。等价于对 $\hat x$ 最小化 $E[\|x_k-\hat x\|^2]$。

物理/几何意义：递推贝叶斯滤波在线性高斯假设下的闭式解。

极限/特例：稳态 Kalman 增益满足离散 ARE；$K_k$ 收敛则滤波器稳态。

典型应用：数字控制、传感器融合。

**知识点14.4　扩展卡尔曼滤波（EKF）**

定义：对非线性系统 $x_{k+1}=f(x_k,u_k)+w_k$、$y_k=h(x_k)+v_k$ 在线性化后套用卡尔曼滤波。

假设：$f,h$ 可微，状态近似高斯。

公式：$\hat x_{k|k-1}=f(\hat x_{k-1|k-1},u_{k-1})$，
$$F_k=\left.\frac{\partial f}{\partial x}\right|_{\hat x_{k-1|k-1}},\quad H_k=\left.\frac{\partial h}{\partial x}\right|_{\hat x_{k|k-1}},$$
$$P_{k|k-1}=F_kP_{k-1|k-1}F_k^T+Q_{k-1},\quad K_k=P_{k|k-1}H_k^T(H_kP_{k|k-1}H_k^T+R_k)^{-1}.$$

推导：对 $f$ 在 $\hat x$ 处一阶泰勒展开 $f(x)\approx f(\hat x)+F(x-\hat x)$，$h$ 同理。代入线性卡尔曼滤波（知识点14.3）即得。

物理/几何意义：用每步局部线性化近似非线性，适用于弱非线性。

极限/特例：强非线性时 EKF 发散，需 UKF/粒子滤波。

典型应用：姿态估计、SLAM。

**知识点14.5　LQG 与分离定理**

定义：线性二次高斯问题是 LQR 与卡尔曼滤波的组合：随机系统下最小化期望二次代价。

假设：噪声高斯、代价二次、能稳能检测。

公式：最优控制
$$u^\ast=-K\hat x,\qquad K=R_u^{-1}B^TP,\quad \hat x\ \text{由卡尔曼滤波给出},$$
其中 $P$ 满足控制 ARE，$P_f$ 满足滤波 ARE。

推导（分离定理）：期望代价 $E[\int(x^TQx+u^TR_uu)dt]$ 中，$x$ 的估计精度与控制器可分离。由确定性等价原理，最优控制等于用 $\hat x$ 代替 $x$ 的 LQR 解。严格证明：代价可分解为 $E[\hat x^TQ\hat x]+E[e^TQe]$，第二项与 $u$ 无关，故最小化第一项即 LQR，得 $u=-K\hat x$。

物理/几何意义：估计与控制可独立设计，闭环特征值为 $A-BK$ 与 $A-LC$ 之并。

极限/特例：确定性等价仅在二次代价、高斯噪声、线性系统下成立。

典型应用：随机最优控制、航空导航控制系统。

**知识点14.6　协方差与信息滤波的等价形式**

定义：信息滤波用信息矩阵 $P^{-1}$ 与信息向量 $P^{-1}\hat x$ 递推。

假设：$R>0$，$P>0$。

公式：信息形式更新
$$P_{k|k}^{-1}=P_{k|k-1}^{-1}+C_d^TR_k^{-1}C_d,\qquad P_{k|k}^{-1}\hat x_{k|k}=P_{k|k-1}^{-1}\hat x_{k|k-1}+C_d^TR_k^{-1}y_k.$$

推导：由后验协方差 $P_{k|k}=(P_{k|k-1}^{-1}+C_d^TR_k^{-1}C_d)^{-1}$（用矩阵求逆引理），两端求逆即第一式。第二式由高斯乘积的指数线性项匹配得到。

物理/几何意义：信息形式在量测维数大于状态维数时计算更高效，且可分布式融合多传感器。

极限/特例：无先验信息时 $P^{-1}\to0$，退化为最小二乘。

典型应用：多传感器融合、分布式估计。

### 第15章 离散时间系统与采样控制

**知识点15.1　离散状态空间模型与 z 变换**

定义：离散系统 $x_{k+1}=A_dx_k+B_du_k$，$y_k=C_dx_k+D_du_k$；z 变换 $X(z)=\sum_{k=0}^\infty x_kz^{-k}$。

假设：采样周期 $T$ 常数、ZOH。

公式：传递函数
$$G(z)=C_d(zI-A_d)^{-1}B_d+D_d=Z\!\left\{ZOH\!\left[G(s)\right]\right\}.$$

推导：对差分方程取 z 变换并利用时移性质 $Z\{x_{k+1}\}=zX(z)-zx_0$，零初值下 $zX=A_dX+B_dU$，故 $X=(zI-A_d)^{-1}B_dU$，输出 $Y=(C_d(zI-A_d)^{-1}B_d+D_d)U$。

物理/几何意义：z 变换把时域差分化为代数方程，极点位置用单位圆判别稳定性。

极限/特例：$z=e^{sT}$ 联系连续与离散极点；$T\to0$ 时 $G(z)\to G(s)$（近似）。

典型应用：数字控制器、DSP 实现。

**知识点15.2　离散系统稳定性判据**

定义：离散系统渐近稳定指 $x_k\to0$，由特征值模长判别。

假设：$x_{k+1}=A_dx_k$。

公式：渐近稳定 $\iff|\lambda_i(A_d)|<1,\ \forall i$；李雅普诺夫方程
$$A_d^TPA_d-P=-Q,\quad Q>0\ \Rightarrow\ P>0.$$

推导：取 $V_k=x_k^TPx_k$，则
$$\Delta V=V_{k+1}-V_k=x_k^T(A_d^TPA_d-P)x_k=-x_k^TQx_k<0.$$
故渐近稳定。反之稳定时 $P=\sum_{k=0}^\infty(A_d^T)^kQ(A_d)^k$ 收敛且满足方程。

物理/几何意义：单位圆内为稳定域，与连续系统左半平面对应。

极限/特例：$|\lambda|=1$ 为临界；$|\lambda|>1$ 发散。

典型应用：数字控制器稳定性判别。

**知识点15.3　离散能控性与能观性**

定义：$x_{k+1}=A_dx_k+B_du_k$ 的能控/能观性。

假设：$A_d,B_d,C_d$ 常阵。

公式：能控 $\iff\mathrm{rank}[B_d\ A_dB_d\ \cdots\ A_d^{n-1}B_d]=n$；能观 $\iff\mathrm{rank}[C_d;C_dA_d;\cdots;C_dA_d^{n-1}]=n$。

推导：$n$ 步状态
$$x_n=A_d^nx_0+\sum_{k=0}^{n-1}A_d^{n-1-k}B_du_k.$$
把 $u_0,\dots,u_{n-1}$ 作为待定量，要求对任意 $x_0$ 能到零，即能控性矩阵满秩；能观性由对偶同理。由凯莱–哈密顿，$n$ 步已足够。

物理/几何意义：与连续系统结构相同，只是用 $A_d$ 代替 $A$。

极限/特例：采样可能破坏能控性（如 $T$ 使某些模态不可辨），出现“采样病态”。

典型应用：数字 LQR、数字极点配置。

**知识点15.4　离散极点配置与死拍控制**

定义：用 $u_k=-K_dx_k$ 配置 $A_d-B_dK_d$ 特征值；死拍控制指在有限步内把状态驱到零。

假设：$(A_d,B_d)$ 能控。

公式：若把闭环特征值全配到 $z=0$，则
$$A_d-B_dK_d\ \text{幂零}\ \Rightarrow x_k=0,\ \forall k\ge n.$$
增益由 $\det(zI-A_d+B_dK_d)=z^n$ 求得。

推导：特征值全为零时，由凯莱–哈密顿 $(A_d-B_dK_d)^n=0$，故
$$x_n=(A_d-B_dK_d)^nx_0=0.$$

物理/几何意义：死拍控制在 $n$ 个采样周期内无误差到达，但控制量可能很大、鲁棒性差。

极限/特例：采样周期越小控制量越大；实际受饱和限制。

典型应用：快速点定位、数控机床。

**知识点15.5　离散 LQR**

定义：离散二次型代价
$$J=\frac12\sum_{k=0}^{N-1}(x_k^TQx_k+u_k^TRu_k)+\frac12x_N^TSx_N.$$

假设：$Q\ge0,R>0$。

公式：最优控制 $u_k=-K_kx_k$，
$$K_k=(R+B_d^TP_{k+1}B_d)^{-1}B_d^TP_{k+1}A_d,$$
$$P_k=Q+A_d^TP_{k+1}A_d-A_d^TP_{k+1}B_d(R+B_d^TP_{k+1}B_d)^{-1}B_d^TP_{k+1}A_d,\quad P_N=S.$$

推导（动态规划）：设值函数 $V_k=\frac12x_k^TP_kx_k$，Bellman 方程
$$V_k=\min_{u_k}\left[\frac12x_k^TQx_k+\frac12u_k^TRu_k+V_{k+1}(A_dx_k+B_du_k)\right].$$
对 $u_k$ 求导：$Ru_k+B_d^TP_{k+1}(A_dx_k+B_du_k)=0$，得增益。代回比较二次项系数得 $P_k$ 递推。

物理/几何意义：离散 ARE 的有限时域形式，是 MPC 的基础。

极限/特例：$N\to\infty$ 得稳态离散 ARE；$R\to0$ 增益增大。

典型应用：模型预测控制、数字最优控制。

**知识点15.6　采样系统与零阶保持器的影响**

定义：连续被控对象经采样—ZOH 构成的采样系统；ZOH 传递函数 $G_{ZOH}(s)=\frac{1-e^{-sT}}{s}$。

假设：采样周期远小于系统时间常数。

公式：
$$G(z)=Z\!\left[\frac{1-e^{-sT}}{s}G_p(s)\right]=(1-z^{-1})Z\!\left[\frac{G_p(s)}{s}\right].$$

推导：ZOH 输出为分段常值，其拉普拉斯变换为 $\frac{1-e^{-sT}}{s}$。采样后
$$Z\{G_{ZOH}(s)G_p(s)\}=(1-z^{-1})Z\!\left\{\frac{G_p(s)}{s}\right\},$$
其中用了 $Z\{e^{-sT}F(s)\}=z^{-1}F(z)$。

物理/几何意义：ZOH 引入半采样周期等效延迟，降低相位裕度。

极限/特例：采样快时 $G(z)\approx G_p(s)$；采样慢时可能失稳。

典型应用：数字控制器采样率选择。

### 第16章 航空航天姿态动力学与控制应用

**知识点16.1　刚体姿态运动学（欧拉角与四元数）**

定义：用欧拉角 $(\phi,\theta,\psi)$ 或四元数 $q=(q_0,q_1,q_2,q_3)$ 描述刚体相对于参考系的姿态。

假设：刚体、参考系正交。

公式（四元数运动学）：
$$\dot q=\frac12\begin{bmatrix}0&-\omega^T\\ \omega&-[\omega\times]\end{bmatrix}q,\qquad [\omega\times]=\begin{bmatrix}0&-\omega_3&\omega_2\\ \omega_3&0&-\omega_1\\ -\omega_2&\omega_1&0\end{bmatrix}.$$

推导：单位四元数 $q=(\cos\frac{\vartheta}{2},\ \sin\frac{\vartheta}{2}\,e)$，表示绕单位轴 $e$ 转 $\vartheta$。由 $\dot q=\frac12q\otimes(0,\omega)$ 的矩阵形式即得。对欧拉角，$ZYX$ 序有
$$\begin{bmatrix}\dot\phi\\\dot\theta\\\dot\psi\end{bmatrix}=\frac{1}{\cos\theta}\begin{bmatrix}\cos\theta&\sin\phi\sin\theta&\cos\phi\sin\theta\\0&\cos\phi\cos\theta&-\sin\phi\cos\theta\\0&\sin\phi&\cos\phi\end{bmatrix}\begin{bmatrix}\omega_1\\\omega_2\\\omega_3\end{bmatrix}.$$

物理/几何意义：四元数无奇点、计算高效；欧拉角直观但存在万向锁。

极限/特例：$\theta=\pm90^\circ$ 欧拉角奇异；四元数全局非奇异。

典型应用：卫星姿态确定、无人机飞控。

**知识点16.2　刚体转动动力学（欧拉方程）**

定义：在体坐标系下，刚体角动量守恒给出欧拉方程。

假设：刚体、主轴坐标系、转动惯量矩阵对角 $J=\mathrm{diag}(J_1,J_2,J_3)$。

公式：
$$J\dot\omega+\omega\times(J\omega)=M,$$
展开
$$J_1\dot\omega_1-(J_2-J_3)\omega_2\omega_3=M_1,$$
$$J_2\dot\omega_2-(J_3-J_1)\omega_3\omega_1=M_2,$$
$$J_3\dot\omega_3-(J_1-J_2)\omega_1\omega_2=M_3.$$

推导：角动量 $H=J\omega$，在惯性系守恒 $\dot H_{in}=M$。变换到体坐标系 $\dot H_{in}=\dot H_b+\omega\times H_b$，故
$$J\dot\omega+\omega\times J\omega=M.$$

物理/几何意义：$\omega\times J\omega$ 为陀螺力矩，是刚体转动的本质非线性项。

极限/特例：$M=0$、$J_1=J_2=J_3$ 时 $\omega$ 常数；自由刚体绕最大/最小惯量轴稳定、绕中间轴不稳定（网球拍定理）。

典型应用：卫星姿态机动、陀螺仪。

**知识点16.3　无力矩刚体的稳定性（中间轴定理）**

定义：自由刚体（$M=0$）绕主轴的转动稳定性。

假设：$J_1<J_2<J_3$，平衡 $\omega=(0,\omega_2,0)$ 附近小扰动。

公式：线性化方程
$$\dot{\tilde\omega}_1=\frac{J_3-J_2}{J_1}\omega_2\tilde\omega_3,\quad \dot{\tilde\omega}_3=\frac{J_2-J_1}{J_3}\omega_2\tilde\omega_1.$$
特征方程
$$\lambda^2=-\frac{(J_3-J_2)(J_2-J_1)}{J_1J_3}\omega_2^2.$$

推导：在 $\omega=(0,\omega_2,0)$ 处对欧拉方程线性化，忽略二阶扰动项，得上述方程组。令 $\tilde\omega_i\propto e^{\lambda t}$ 得特征方程。因 $J_1<J_2<J_3$，$(J_3-J_2)>0$、$(J_2-J_1)>0$，故 $\lambda^2<0$，为纯虚根（临界）。对中间轴 $J_2$ 扰动，实际非线性效应导致不稳定（网球拍定理）。

物理/几何意义：绕最大或最小惯量轴自旋稳定，绕中间轴不稳定。

极限/特例：能量—角动量椭球的几何分析给出同一结论。

典型应用：卫星自旋稳定、空间碎片翻滚。

**知识点16.4　姿态稳定性与李雅普诺夫分析**

定义：用李雅普诺夫函数分析姿态控制闭环稳定性。

假设：姿态误差用四元数 $q_e$ 表示，角速度误差 $\omega_e$。

公式：构造
$$V=\frac12\omega_e^TJ\omega_e+k\left[(q_{e0}-1)^2+q_{ev}^Tq_{ev}\right],$$
控制律 $M=-k_pq_{ev}-k_d\omega_e$。

推导：对 $V$ 求导，利用四元数运动学 $\dot q_{ev}=\frac12(q_{e0}I+[q_{ev}\times])\omega_e$，可得
$$\dot V=\omega_e^T(J\dot\omega_e)+2kq_{ev}^T\dot q_{ev}.$$
代入 $J\dot\omega_e=\cdots+M$ 与 $M=-k_pq_{ev}-k_d\omega_e$，交叉项通过陀螺项反对称消去，最终
$$\dot V=-(k_d-k\cdots)\|\omega_e\|^2\le0.$$
由拉萨尔不变集原理得 $\omega_e\to0$、$q_{ev}\to0$。

物理/几何意义：$V$ 是姿态误差的“能量”，PD 控制使能量单调下降。

极限/特例：$k_d$ 足够大保证 $\dot V\le0$；仅比例项产生常值偏差。

典型应用：卫星姿态控制器设计。

**知识点16.5　柔性航天器与模态振动控制**

定义：带柔性附件的航天器动力学含刚性姿态与弹性模态耦合。

假设：小变形、模态截断。

公式：
$$\begin{aligned}J\ddot\theta+\sum_i c_i\ddot\eta_i&=M,\\ \ddot\eta_i+2\zeta_i\omega_i\dot\eta_i+\omega_i^2\eta_i+c_i\ddot\theta&=0.\end{aligned}$$

推导：用拉格朗日方程（第2、4章），广义坐标取 $\theta$ 与模态坐标 $\eta_i$，动能含耦合项，势能为弹性势能，得耦合方程。

物理/几何意义：姿态机动激发柔性振动，振动又通过耦合 $c_i\ddot\theta$ 反作用姿态。

极限/特例：$c_i=0$ 时解耦；用输入成型（input shaping）抑制残余振动。

典型应用：大型太阳翼卫星、空间机械臂。

**知识点16.6　编队飞行与相对运动（Hill–Clohessy–Wiltshire 方程）**

定义：近圆轨道上两个航天器相对运动的线性化方程。

假设：参考轨道近圆、相对距离远小于轨道半径。

公式：在旋转的 LVLH 坐标系中
$$\ddot x-3n^2x-2n\dot y=a_x,$$
$$\ddot y+2n\dot x=a_y,$$
$$\ddot z+n^2z=a_z,$$
其中 $n$ 为轨道角速度，$x$ 径向、$y$ 沿迹、$z$ 法向。

推导：在旋转参考系中写牛顿第二定律，对地心引力做一阶泰勒展开：引力梯度张量在 LVLH 下为 $\mathrm{diag}(2n^2,-n^2,-n^2)$（相对项），加上离心与科氏力，得上述方程。

物理/几何意义：$x$-$y$ 平面存在漂移与本征振荡（周期轨、绕飞），$z$ 为简谐运动。

极限/特例：$n\to0$ 退化为自由空间二体相对运动；$z$ 运动与面内解耦。

典型应用：卫星编队构型保持、交会对接。

**知识点16.7　轨迹跟踪控制与非线性逆动力学**

定义：用反馈线性化或逆动力学使系统输出跟踪期望轨迹。

假设：系统相对阶等于状态数、零动态稳定（最小相位）。

公式：对 $\dot x=f(x)+g(x)u$、$y=h(x)$，求李导数
$$y^{(r)}=L_f^rh+L_gL_f^{r-1}h\,u.$$
令其等于期望 $y_d^{(r)}$，得
$$u=\frac{1}{L_gL_f^{r-1}h}\left(y_d^{(r)}-L_f^rh\right).$$

推导：对输出反复求导 $r$ 次直到出现 $u$。由 Lie 导数定义
$$L_fh=\frac{\partial h}{\partial x}f,\quad L_f^rh=L_f(L_f^{r-1}h),$$
第 $r$ 阶导含 $u$ 的系数 $L_gL_f^{r-1}h$。若该系数非零，可直接解出 $u$，实现输入—输出线性化。

物理/几何意义：用反馈精确抵消非线性，构造期望的线性误差动态。

极限/特例：零动态不稳定时内态发散；系数为零时相对阶需提高。

典型应用：无人机/飞行器轨迹跟踪、机器人。

**知识点16.8　模型预测控制在航天器上的应用**

定义：在线求解有限时域优化问题，取首步控制实施并滚动更新。

假设：系统模型（线性或线性化）已知、约束为凸。

公式：每步求解
$$\min_{U}\ \sum_{k=0}^{N-1}\left(\|x_k-x_k^{ref}\|_Q^2+\|u_k\|_R^2\right)+\|x_N-x_N^{ref}\|_S^2,$$
$$\text{s.t.}\ x_{k+1}=A_dx_k+B_du_k,\quad x_k\in\mathcal X,\ u_k\in\mathcal U.$$

推导：由离散 LQR（知识点15.5）加约束。把 $x_k$ 用 $x_0$ 与 $U$ 表示（预测方程），代入代价得二次规划，解出 $U^\ast=\{u_0^\ast,\dots\}$，实施 $u_0^\ast$ 后移位重复。

物理/几何意义：预测未来、滚动优化、显式处理约束，适合容错和受限任务。

极限/特例：无约束时退化为 LQR；终端代价/约束保证递归可行与稳定。

典型应用：交会对接、推力受限姿态机动。

**知识点16.9　航天器姿态确定与卡尔曼滤波**

定义：用星敏感器、陀螺、磁强计等量测估计姿态与角速度。

假设：陀螺零偏建模为常值状态，星敏量测为矢量观测。

公式：状态 $x=[q;\ b]$（四元数与陀螺零偏），陀螺模型
$$\dot q=\frac12\Omega(\omega)q,\qquad \dot b=w_b,$$
量测 $y=h(q)+v$。用 EKF（第14章）或 UKF 递推。

推导：陀螺输出 $\omega_m=\omega+b+\eta$，把 $b$ 扩张为状态并用 $\dot b=0$ 建模。对非线性姿态运动学与矢量观测做线性化，套用 EKF 预测—校正。

物理/几何意义：陀螺高频差、星敏低频好，融合后得到无偏、高带宽姿态估计。

极限/特例：纯陀螺积分随时间的漂移被星敏校正；$q$ 需归一化处理。

典型应用：卫星姿态确定、深空探测器导航。

---
⬆ 上级领域：[[电子·控制·系统]]
