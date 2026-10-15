---
课程名称: 动力学与控制核心贯通课I
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

# 动力学与控制核心贯通课I


### 第1章 运动学基础与坐标系描述

#### 知识点1.1 位矢、速度与加速度的定义与微分结构

**定义**：在惯性参考系中取定原点 $O$，质点位置由位置矢量 $\vec r(t)$ 表示，速度与加速度定义为
$$\vec v=\frac{\mathrm d\vec r}{\mathrm dt},\qquad \vec a=\frac{\mathrm d\vec v}{\mathrm dt}=\frac{\mathrm d^2\vec r}{\mathrm dt^2}.$$

**假设**：时间 $t$ 为绝对参量（牛顿绝对时空），参考系为惯性系，$\vec r$ 二次连续可微。

**推导**：由极限定义
$$\vec v=\lim_{\Delta t\to0}\frac{\vec r(t+\Delta t)-\vec r(t)}{\Delta t}.$$
将 $\vec r(t+\Delta t)$ 在 $t$ 处 Taylor 展开：
$$\vec r(t+\Delta t)=\vec r(t)+\dot{\vec r}(t)\Delta t+\frac12\ddot{\vec r}(t)\Delta t^2+O(\Delta t^3).$$
代入并令 $\Delta t\to0$，得 $\vec v=\dot{\vec r}$；再对 $\vec v$ 重复上述过程得 $\vec a=\ddot{\vec r}$。

**物理意义**：速度是位矢对时间的一阶变化率，加速度是位矢的二阶变化率；加速度不为零是受力的运动学标志（配合第二定律）。

**特例**：$\vec a\equiv0$ 时 $\vec r=\vec r_0+\vec v_0t$（匀速直线运动）；$\vec v\equiv0$ 时质点静止。

**应用**：轨迹消参、机构运动分析、数值积分（Verlet、Runge–Kutta）中加速度的离散实现。

#### 知识点1.2 直角坐标与极坐标的运动学变换

**定义**：平面直角坐标 $(x,y)$ 与极坐标 $(r,\theta)$ 的关系为
$$x=r\cos\theta,\qquad y=r\sin\theta.$$

**假设**：$r>0$，$\theta$ 连续可微，$\hat e_r,\hat e_\theta$ 为随体单位正交基。

**推导**：设 $\hat e_r=\cos\theta\,\hat i+\sin\theta\,\hat j$，$\hat e_\theta=-\sin\theta\,\hat i+\cos\theta\,\hat j$。对时间求导
$$\because \dot{\hat e}_r=-\sin\theta\,\dot\theta\,\hat i+\cos\theta\,\dot\theta\,\hat j=\dot\theta\,\hat e_\theta,$$
$$\therefore \dot{\hat e}_\theta=-\cos\theta\,\dot\theta\,\hat i-\sin\theta\,\dot\theta\,\hat j=-\dot\theta\,\hat e_r.$$
又 $\vec r=r\hat e_r$，故
$$\vec v=\dot r\hat e_r+r\dot{\hat e}_r=\dot r\hat e_r+r\dot\theta\hat e_\theta.$$
再求导
$$\vec a=\ddot r\hat e_r+\dot r\dot{\hat e}_r+\frac{\mathrm d}{\mathrm dt}(r\dot\theta)\hat e_\theta+r\dot\theta\dot{\hat e}_\theta$$
$$=(\ddot r-r\dot\theta^2)\hat e_r+(r\ddot\theta+2\dot r\dot\theta)\hat e_\theta.$$

**几何意义**：$-r\dot\theta^2$ 为向心加速度，$2\dot r\dot\theta$ 为科里奥利加速度（极坐标下面积速度导致的耦合项）。

**特例**：圆周运动 $r=R$ 常数时 $\vec a=-R\dot\theta^2\hat e_r+R\ddot\theta\hat e_\theta$。

**应用**：行星运动、旋转机械、极坐标下的中心力问题。

#### 知识点1.3 柱坐标与球坐标的加速度公式

**定义**：柱坐标 $(\rho,\phi,z)$：
$$x=\rho\cos\phi,\ y=\rho\sin\phi,\ z=z.$$
球坐标 $(r,\theta,\phi)$（$\theta$ 为极角，$\phi$ 为方位角）：
$$x=r\sin\theta\cos\phi,\ y=r\sin\theta\sin\phi,\ z=r\cos\theta.$$

**假设**：坐标变换非奇异（$\rho=0$ 或 $r=0$ 除外），单位基正交。

**推导（柱坐标）**：$\hat e_\rho=\cos\phi\,\hat i+\sin\phi\,\hat j$，$\hat e_\phi=-\sin\phi\,\hat i+\cos\phi\,\hat j$，$\hat e_z=\hat k$，且有 $\dot{\hat e}_\rho=\dot\phi\hat e_\phi,\ \dot{\hat e}_\phi=-\dot\phi\hat e_\rho$。则
$$\vec v=\dot\rho\hat e_\rho+\rho\dot\phi\hat e_\phi+\dot z\hat e_z,$$
$$\vec a=(\ddot\rho-\rho\dot\phi^2)\hat e_\rho+(\rho\ddot\phi+2\dot\rho\dot\phi)\hat e_\phi+\ddot z\hat e_z.$$

**推导（球坐标）**：利用 $\dot{\hat e}_r=\dot\theta\hat e_\theta+\dot\phi\sin\theta\hat e_\phi$，$\dot{\hat e}_\theta=-\dot\theta\hat e_r+\dot\phi\cos\theta\hat e_\phi$，$\dot{\hat e}_\phi=-\dot\phi(\sin\theta\hat e_r+\cos\theta\hat e_\theta)$，可得
$$\vec a=(\ddot r-r\dot\theta^2-r\dot\phi^2\sin^2\theta)\hat e_r+(r\ddot\theta+2\dot r\dot\theta-r\dot\phi^2\sin\theta\cos\theta)\hat e_\theta$$
$$+\,(r\ddot\phi\sin\theta+2\dot r\dot\phi\sin\theta+2r\dot\theta\dot\phi\cos\theta)\hat e_\phi.$$

**物理意义**：各坐标方向的加速度由径向二阶项、离心项和科氏项构成。

**特例**：球坐标中 $\theta=\pi/2$ 退化为极坐标。

**应用**：大气轨道动力学、地球物理流体、天线指向控制。

#### 知识点1.4 自然坐标与 Frenet 标架

**定义**：沿轨迹弧长 $s$ 定义切向量 $\hat t$、主法向量 $\hat n$、副法向量 $\hat b=\hat t\times\hat n$，曲率 $\kappa$、挠率 $\tau$ 满足 Frenet–Serret 公式
$$\frac{\mathrm d\hat t}{\mathrm ds}=\kappa\hat n,\quad \frac{\mathrm d\hat n}{\mathrm ds}=-\kappa\hat t+\tau\hat b,\quad \frac{\mathrm d\hat b}{\mathrm ds}=-\tau\hat n.$$

**假设**：轨迹光滑且 $\kappa>0$，$\hat t$ 为单位向量。

**推导**：由 $\hat t=\mathrm d\vec r/\mathrm ds$ 且 $\hat t\cdot\hat t=1$，求导得 $2\hat t\cdot(\mathrm d\hat t/\mathrm ds)=0$，故 $\mathrm d\hat t/\mathrm ds\perp\hat t$，方向即 $\hat n$，大小定义为 $\kappa$。同理对 $\hat n\cdot\hat n=1$ 求导得 $\mathrm d\hat n/\mathrm ds\perp\hat n$，可分解为 $-\kappa\hat t+\tau\hat b$；对 $\hat b=\hat t\times\hat n$ 求导并利用叉乘法则得第三式。

**运动学分解**：
$$\vec v=v\hat t,\qquad \vec a=\dot v\hat t+\frac{v^2}{R}\hat n,\qquad R=1/\kappa.$$

**几何意义**：$\kappa$ 描述轨迹弯曲程度，$\tau$ 描述轨迹偏离平面的扭转程度。

**特例**：直线 $\kappa=0$；圆周 $\tau=0$，$R$ 恒定。

**应用**：铁路曲线设计、粒子加速器、车辆横向动力学。

#### 知识点1.5 刚体运动学、欧拉角与角速度

**定义**：刚体位形由质心位置 $\vec R$ 与姿态矩阵 $\mathbf A$ 描述，任意物质点
$$\vec r_i=\vec R+\mathbf A\,\vec\rho_i',$$
其中 $\vec\rho_i'$ 为体坐标系中的常向量。

**假设**：刚体为连续介质，任意两点距离不变，$\mathbf A$ 为正交矩阵 $\mathbf A^\top\mathbf A=\mathbf I$。

**推导**：正交性求导 $\dot{\mathbf A}^\top\mathbf A+\mathbf A^\top\dot{\mathbf A}=0$，令 $\boldsymbol\Omega=\dot{\mathbf A}\mathbf A^\top$，则 $\boldsymbol\Omega+\boldsymbol\Omega^\top=0$，即 $\boldsymbol\Omega$ 反对称。用 $\boldsymbol\omega$ 表示其反偶向量，$\boldsymbol\Omega\vec x=\boldsymbol\omega\times\vec x$。于是体坐标相对惯性系的速度
$$\dot{\vec\rho}=\dot{\mathbf A}\vec\rho'=\boldsymbol\Omega\mathbf A\vec\rho'=\boldsymbol\omega\times\vec\rho.$$
欧拉角 $(\psi,\theta,\phi)$（进动、章动、自转）给出
$$\boldsymbol\omega=\dot\psi\hat k+\dot\theta\hat{\hat e}_{N}+\dot\phi\hat e_3,$$
展开为
$$\boldsymbol\omega=\begin{pmatrix}\dot\phi\sin\theta\sin\psi+\dot\theta\cos\psi\\ \dot\phi\sin\theta\cos\psi-\dot\theta\sin\psi\\ \dot\phi\cos\theta+\dot\psi\end{pmatrix}_{\text{惯性分量}}.$$

**物理意义**：角速度是姿态矩阵导数的反偶向量；欧拉角速度与角速度间存在表征姿态的 Jacobi 矩阵。

**特例**：小角度时 $\boldsymbol\omega\approx(\dot\phi,\dot\theta,\dot\psi)$。

**应用**：航天器姿态、陀螺仪、机器人末端姿态描述。

#### 知识点1.6 伽利略变换与相对运动

**定义**：惯性系 $O$ 与以 $\vec R(t)$ 平动、以 $\boldsymbol\omega(t)$ 转动的参考系 $O'$ 之间，任意向量 $\vec q$ 满足输运定理
$$\left(\frac{\mathrm d\vec q}{\mathrm dt}\right)_{O}=\left(\frac{\mathrm d\vec q}{\mathrm dt}\right)_{O'}+\boldsymbol\omega\times\vec q.$$

**假设**：$O$ 为惯性系，$O'$ 为任意运动参考系，$\vec q$ 为任意矢量。

**推导**：设 $\vec q=\sum_i q_i'\hat e_i'$，其中 $\hat e_i'$ 为 $O'$ 的基。求导
$$\dot{\vec q}=\sum_i\dot q_i'\hat e_i'+\sum_i q_i'\dot{\hat e}_i'.$$
对基向量有 $\dot{\hat e}_i'=\boldsymbol\omega\times\hat e_i'$（因 $|\hat e_i'|=1$ 且旋转生成），故
$$\dot{\vec q}=\left(\dot{\vec q}\right)_{O'}+\boldsymbol\omega\times\vec q.$$

**推论**：
$$\vec v_O=\vec v_{O'}+\vec R_{\!O O'}+\boldsymbol\omega\times\vec r',$$
$$\vec a_O=\vec a_{O'}+\vec A_{OO'}+2\boldsymbol\omega\times\vec v'+\boldsymbol\omega\times(\boldsymbol\omega\times\vec r')+\dot{\boldsymbol\omega}\times\vec r'.$$

**物理意义**：$2\boldsymbol\omega\times\vec v'$ 为科氏加速度，$\boldsymbol\omega\times(\boldsymbol\omega\times\vec r')$ 为向心加速度，$\dot{\boldsymbol\omega}\times\vec r'$ 为切向（欧拉）加速度。

**特例**：$O'$ 相对 $O$ 平动（$\boldsymbol\omega=0$）时得伽利略变换；$\vec q$ 静止于 $O'$ 时仅有牵连项。

**应用**：地球自转对弹道与傅科摆的影响、旋转机械振动。

### 第2章 牛顿质点动力学与中心力问题

#### 知识点2.1 牛顿三定律、动量定理与角动量定理

**定义**：动量 $\vec p=m\vec v$、角动量 $\vec L=\vec r\times\vec p$、力矩 $\vec M=\vec r\times\vec F$。

**假设**：质量 $m$ 常数，且所考察为惯性系。

**推导**：由第二定律
$$\vec F=\frac{\mathrm d\vec p}{\mathrm dt}=m\vec a.$$
角动量求导
$$\dot{\vec L}=\dot{\vec r}\times\vec p+\vec r\times\dot{\vec p}=\vec v\times(m\vec v)+\vec r\times\vec F=\vec M,$$
其中 $\vec v\times\vec v=0$。对系统求和并利用牛顿第三定律（内力成对且共线）
$$\sum_i\vec F_i^{(e)}=\frac{\mathrm d\vec P}{\mathrm dt},\qquad \vec P=\sum_i m_i\vec v_i,$$
$$\sum_i\vec M_i^{(e)}=\frac{\mathrm d\vec L_O}{\mathrm dt}.$$

**物理意义**：内力不改变总动量与总角动量；外力是系统动量与角动量变化的唯一来源。

**特例**：$\vec F=0$ 时动量守恒；$\vec M=0$ 时角动量守恒（开普勒第二定律的本质）。

**应用**：碰撞、火箭推进、卫星机动。

#### 知识点2.2 功–能定理与势能

**定义**：力 $\vec F$ 沿路径做功 $W=\int_{\vec r_1}^{\vec r_2}\vec F\cdot\mathrm d\vec r$，动能 $T=\frac12m v^2$。

**假设**：力可沿路径积分；若力为保守力，则存在势函数 $V(\vec r)$ 使 $\vec F=-\nabla V$。

**推导**：
$$\because \frac{\mathrm dT}{\mathrm dt}=\frac{\mathrm d}{\mathrm dt}\left(\frac12m\vec v\cdot\vec v\right)=m\vec v\cdot\dot{\vec v}=\vec F\cdot\vec v,$$
$$\therefore \int_{t_1}^{t_2}\frac{\mathrm dT}{\mathrm dt}\mathrm dt=T_2-T_1=\int\vec F\cdot\vec v\,\mathrm dt=\int\vec F\cdot\mathrm d\vec r=W.$$
若 $\vec F=-\nabla V$，则
$$W=-\int_{\vec r_1}^{\vec r_2}\nabla V\cdot\mathrm d\vec r=-(V_2-V_1),$$
故 $T_2+V_2=T_1+V_1$。

**意义**：保守系统中机械能守恒；非保守力使机械能改变。

**特例**：一维 $F(x)=-\mathrm dV/\mathrm dx$，平衡点 $F=0$ 对应势能极值。

**应用**：能量判据、振动稳定性、数值积分能量检验。

#### 知识点2.3 一维势能分析与平衡稳定性

**定义**：势能函数 $V(x)$，平衡位置满足 $V'(x_0)=0$。

**假设**：力保守且 $V$ 二次连续可微。

**推导**：在 $x_0$ 处 Taylor 展开
$$V(x_0+\xi)=V(x_0)+V'(x_0)\xi+\frac12V''(x_0)\xi^2+O(\xi^3).$$
因 $V'(x_0)=0$，
$$V(x_0+\xi)-V(x_0)\approx\frac12V''(x_0)\xi^2.$$
故 $V''(x_0)>0$ 时势能局部极小，为稳定平衡；$V''(x_0)<0$ 时不稳定；$V''(x_0)=0$ 需更高阶判断。运动方程 $m\ddot\xi=-V''(x_0)\xi$ 给出 $\omega=\sqrt{V''(x_0)/m}$。

**几何意义**：势能曲线的极小点为稳定平衡，极大点为不稳定平衡。

**特例**：单摆 $V=-mgl\cos\theta$，$\theta=0$ 稳定、$\theta=\pi$ 不稳定。

**应用**：结构屈曲、原子间势、分子振动频率。

#### 知识点2.4 中心力场与两体约化

**定义**：中心力 $\vec F=F(r)\hat e_r$，两体问题中约化质量 $\mu=m_1m_2/(m_1+m_2)$。

**假设**：力仅依赖距离，系统孤立。

**推导**：相对位矢 $\vec r=\vec r_1-\vec r_2$。由 $\vec F_{12}=-\vec F_{21}$，
$$m_1\ddot{\vec r}_1=\vec F_{12},\qquad m_2\ddot{\vec r}_2=-\vec F_{12}.$$
相减得到
$$\ddot{\vec r}=\ddot{\vec r}_1-\ddot{\vec r}_2=\left(\frac1{m_1}+\frac1{m_2}\right)\vec F_{12}=\frac{\vec F_{12}}{\mu}.$$
角动量守恒给出面积速度
$$\dot{\vec L}=\vec r\times\mu\ddot{\vec r}=\vec r\times\vec F=0\ \Rightarrow\ \vec L=\text{const},$$
平面运动且 $\frac12r^2\dot\theta=L/(2\mu)$ 恒定。

**意义**：两体问题化为一等效单体问题，中心力作用在质心系中决定相对轨道。

**特例**：$m_1\ll m_2$ 时 $\mu\to m_1$，回到行星绕太阳问题。

**应用**：双星、原子模型、卫星编队。

#### 知识点2.5 开普勒问题与圆锥曲线轨道

**定义**：引力 $F(r)=-GMm/r^2$，比角动量 $h=r^2\dot\theta$。

**假设**：牛顿引力、质量恒定、无其他扰动。

**推导**：由 $h=r^2\dot\theta$ 及径向方程
$$\ddot r-r\dot\theta^2=-\frac{GM}{r^2}.$$
作变换 $u=1/r$，利用 $\dot r=-\frac{h}{r^2}\frac{\mathrm dr}{\mathrm du}\frac{\mathrm du}{\mathrm d\theta}$ 得
$$\ddot r=-h^2u^2\frac{\mathrm d^2u}{\mathrm d\theta^2},\qquad r\dot\theta^2=h^2u^3.$$
代入径向方程
$$-h^2u^2\frac{\mathrm d^2u}{\mathrm d\theta^2}-h^2u^3=-\frac{GM}{r^2}=-GMu^2,$$
$$u''+u=\frac{GM}{h^2}.$$
通解
$$u=\frac{GM}{h^2}\bigl(1+e\cos(\theta-\theta_0)\bigr),$$
即
$$r=\frac{p}{1+e\cos(\theta-\theta_0)},\qquad p=\frac{h^2}{GM}.$$
由几何即得圆锥曲线：$e<1$ 椭圆、$e=1$ 抛物线、$e>1$ 双曲线。对椭圆，半长轴 $a=p/(1-e^2)$，周期由面积律
$$T=\frac{\pi ab}{h/2}=\frac{2\pi a^{3/2}}{\sqrt{GM}},$$
即开普勒第三定律 $T^2\propto a^3$。

**意义**：平方反比引力下轨道为圆锥曲线；能量 $E<0$ 对应椭圆。

**特例**：圆轨道 $e=0$，$v=\sqrt{GM/r}$。

**应用**：卫星轨道设计、霍曼转移、逃逸速度。

#### 知识点2.6 有心力场中的有效势与维里定理

**定义**：有效势 $V_{\text{eff}}(r)=V(r)+\dfrac{L^2}{2\mu r^2}$。

**假设**：中心力、平面运动、角动量守恒。

**推导**：径向动能 $T_r=\frac12\mu\dot r^2$，横向动能 $\frac12\mu r^2\dot\theta^2=\frac{L^2}{2\mu r^2}$，故
$$E=\frac12\mu\dot r^2+V_{\text{eff}}(r).$$
径向平衡 $V_{\text{eff}}'(r)=0$ 给出圆轨道的离心平衡条件。

**维里定理推导**：设束缚系统 $G=\sum_i\vec p_i\cdot\vec r_i$，求导并取时间平均
$$\frac{\mathrm dG}{\mathrm dt}=\sum_i\vec p_i\cdot\dot{\vec r}_i+\sum_i\dot{\vec p}_i\cdot\vec r_i=2T+\sum_i\vec F_i\cdot\vec r_i.$$
束缚运动时间平均 $\overline{\dot G}=0$，对平方反比力 $\sum_i\vec F_i\cdot\vec r_i=-V$，故
$$2\overline T+\overline V=0.$$
对 $n$ 次齐次势，$2\overline T=n\overline V$。

**意义**：连接动能与势能的统计平均，用于判断系统束缚性与稳定性。

**特例**：谐振子 $n=2$，$\overline T=\overline V$。

**应用**：星系动力学、分子团、束缚轨道估计。

### 第3章 振动理论与模态分析

#### 知识点3.1 单自由度无阻尼自由振动

**定义**：弹簧–质量系统 $m\ddot x+kx=0$，固有频率 $\omega_n=\sqrt{k/m}$。

**假设**：线性弹簧（胡克定律）、无阻尼、无外激励。

**推导**：牛顿第二定律
$$m\ddot x=-kx\ \Rightarrow\ \ddot x+\omega_n^2x=0.$$
设 $x=e^{\lambda t}$，得特征方程 $\lambda^2+\omega_n^2=0$，$\lambda=\pm i\omega_n$，故
$$x(t)=A\cos\omega_n t+B\sin\omega_n t=C\cos(\omega_n t-\varphi),$$
其中 $C=\sqrt{A^2+B^2}$，$\tan\varphi=B/A$。

**意义**：位移随时间做简谐振荡，周期 $T=2\pi/\omega_n$，频率不依赖振幅。

**特例**：初值 $x(0)=x_0,\ \dot x(0)=v_0$ 时 $A=x_0,\ B=v_0/\omega_n$，$C=\sqrt{x_0^2+v_0^2/\omega_n^2}$。

**应用**：结构基频估算、隔振设计、传感器标定。

#### 知识点3.2 单自由度阻尼自由振动

**定义**：黏性阻尼 $c$，方程 $m\ddot x+c\dot x+kx=0$，阻尼比 $\zeta=c/(2\sqrt{km})$。

**假设**：阻尼与速度成正比、线性。

**推导**：化为
$$\ddot x+2\zeta\omega_n\dot x+\omega_n^2x=0.$$
特征方程
$$\lambda^2+2\zeta\omega_n\lambda+\omega_n^2=0,\qquad \lambda_{1,2}=\omega_n\bigl(-\zeta\pm\sqrt{\zeta^2-1}\bigr).$$

- 欠阻尼 $\zeta<1$：$\lambda=-\zeta\omega_n\pm i\omega_d$，$\omega_d=\omega_n\sqrt{1-\zeta^2}$，
$$x=e^{-\zeta\omega_n t}\bigl(A\cos\omega_dt+B\sin\omega_dt\bigr).$$
- 临界阻尼 $\zeta=1$：$x=(A+Bt)e^{-\omega_nt}$。
- 过阻尼 $\zeta>1$：两个实指数组合。

对数衰减率
$$\delta=\ln\frac{x(t)}{x(t+T_d)}=\frac{2\pi\zeta}{\sqrt{1-\zeta^2}}.$$

**意义**：阻尼耗散能量，振幅按指数衰减；$\zeta$ 决定响应形态。

**特例**：$\zeta=0$ 退化为无阻尼；$\zeta\ll1$ 时 $\omega_d\approx\omega_n$，$\delta\approx2\pi\zeta$。

**应用**：振动测试识别阻尼、减振器设计。

#### 知识点3.3 单自由度受迫振动与频率响应

**定义**：外激励 $F_0\cos\Omega t$，方程 $m\ddot x+c\dot x+kx=F_0\cos\Omega t$。

**假设**：线性系统、稳态响应、激励谐波。

**推导**：设稳态 $x=X\cos(\Omega t-\phi)$，代入并用 $e^{i\Omega t}$ 复表示。令 $x=\Re\{\hat Xe^{i\Omega t}\}$，则
$$(-m\Omega^2+ic\Omega+k)\hat X=F_0,$$
$$\hat X=\frac{F_0}{k-m\Omega^2+ic\Omega}.$$
无量纲化（$\omega_n^2=k/m$，$\zeta=c/(2\sqrt{km})$，$r=\Omega/\omega_n$）
$$\frac{X}{F_0/k}=\frac1{\sqrt{(1-r^2)^2+(2\zeta r)^2}},\qquad \tan\phi=\frac{2\zeta r}{1-r^2}.$$
共振峰在 $r=\sqrt{1-2\zeta^2}$（$0<\zeta<1/\sqrt2$），峰值放大因子 $\approx 1/(2\zeta\sqrt{1-\zeta^2})$。

**意义**：响应幅值由静变形 $F_0/k$ 被动力放大因子放大，相位在共振处滞后 $90^\circ$。

**特例**：$\Omega\to0$ 时 $X\to F_0/k$；$\Omega\to\infty$ 时 $X\to0$；共振 $r=1$ 时 $X=F_0/(c\omega_n)$。

**应用**：旋转机械临界转速、地震响应谱、动力吸振器。

#### 知识点3.4 多自由度系统与广义特征值问题

**定义**：$n$ 自由度系统 $M\ddot{\vec q}+K\vec q=\vec0$，$M$ 对称正定，$K$ 对称。

**假设**：线性、无阻尼、保守、无陀螺项。

**推导**：设 $\vec q=\vec\phi e^{i\omega t}$，代入得
$$(K-\omega^2M)\vec\phi=0.$$
存在非零解要求
$$\det(K-\omega^2M)=0,$$
为关于 $\omega^2$ 的 $n$ 次特征方程，给出 $n$ 个特征值 $\omega_1^2\le\cdots\le\omega_n^2$ 与特征向量 $\vec\phi_i$。由 $M,K$ 对称正定，特征值实且正，特征向量满足
$$\vec\phi_i^\top M\vec\phi_j=m_i\delta_{ij},\qquad \vec\phi_i^\top K\vec\phi_j=k_i\delta_{ij}=m_i\omega_i^2\delta_{ij}.$$

**推导（正交性）**：对 $i\ne j$
$$\vec\phi_j^\top K\vec\phi_i=\omega_i^2\vec\phi_j^\top M\vec\phi_i,\qquad \vec\phi_i^\top K\vec\phi_j=\omega_j^2\vec\phi_i^\top M\vec\phi_j.$$
两式转置后比较，利用对称性得 $(\omega_i^2-\omega_j^2)\vec\phi_i^\top M\vec\phi_j=0$，故正交。

**意义**：模态是系统独立振动形态；响应用模态坐标解耦。

**特例**：$\omega=0$ 对应刚体模态（$K$ 奇异）。

**应用**：有限元模态分析、桥梁抗风、航天器频率规划。

#### 知识点3.5 模态叠加法与解耦

**定义**：物理坐标 $\vec q=\Phi\vec\eta$，$\Phi=[\vec\phi_1,\dots,\vec\phi_n]$ 为模态矩阵。

**假设**：模态完备（或取截断模态），比例阻尼。

**推导**：代入 $M\ddot{\vec q}+C\dot{\vec q}+K\vec q=\vec F$ 并左乘 $\Phi^\top$：
$$\Phi^\top M\Phi\ddot{\vec\eta}+\Phi^\top C\Phi\dot{\vec\eta}+\Phi^\top K\Phi\vec\eta=\Phi^\top\vec F.$$
利用正交性化为对角质量阵 $M_d$ 与刚度阵 $K_d$；若 $C$ 为比例阻尼 $C=\alpha M+\beta K$，则 $\Phi^\top C\Phi$ 亦对角。于是
$$\ddot\eta_i+2\zeta_i\omega_i\dot\eta_i+\omega_i^2\eta_i=\frac{f_i(t)}{m_i},\qquad f_i=\vec\phi_i^\top\vec F.$$
再按 SDOF 求解后叠加。

**意义**：将耦合大系统化为多个独立单自由度问题。

**特例**：仅激励某阶模态时仅该模态响应非零。

**应用**：结构动力学仿真、振动控制、模型降阶。

#### 知识点3.6 拉格朗日法建立振动方程

**定义**：保守系统 $L=T-V$，耗散力 $Q_i^{d}=-\partial\mathcal R/\partial\dot q_i$，瑞利耗散函数 $\mathcal R=\tfrac12\sum c_{ij}\dot q_i\dot q_j$。

**假设**：约束完整、力可写为势函数或广义力。

**推导**：对 $T=\frac12\dot{\vec q}^\top M\dot{\vec q}$，$V=\frac12\vec q^\top K\vec q$，
$$L=\frac12\dot{\vec q}^\top M\dot{\vec q}-\frac12\vec q^\top K\vec q.$$
代入
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_i}-\frac{\partial L}{\partial q_i}=Q_i,$$
得 $M\ddot{\vec q}+K\vec q=\vec Q$。加入阻尼
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial T}{\partial\dot q_i}+\frac{\partial\mathcal R}{\partial\dot q_i}+\frac{\partial V}{\partial q_i}=Q_i^{e},$$
即 $M\ddot{\vec q}+C\dot{\vec q}+K\vec q=\vec Q^{e}$。

**意义**：统一框架建立复杂系统振动方程，避免逐个隔离体分析。

**特例**：无阻尼保守系统 $C=0$。

**应用**：多关节机械臂、车辆悬架建模、转子–轴承系统。

### 第4章 虚功原理与达朗贝尔原理

#### 知识点4.1 约束、自由度与广义坐标

**定义**：完整约束 $f_\alpha(\vec r_1,\dots,\vec r_N,t)=0$ 或 $f_\alpha(q_1,\dots,q_n,t)=0$；广义坐标 $q_1,\dots,q_n$ 为独立参数，$n=3N-k$（$k$ 为独立约束数）。

**假设**：约束完整（可积）、理想（约束力不做虚功）、定常或非定常光滑。

**推导**：由约束方程微分
$$\mathrm df_\alpha=\sum_{i=1}^N\frac{\partial f_\alpha}{\partial\vec r_i}\cdot\mathrm d\vec r_i+\frac{\partial f_\alpha}{\partial t}\mathrm dt=0.$$
引入广义坐标后 $\vec r_i=\vec r_i(q,t)$，则
$$\delta\vec r_i=\sum_{j=1}^n\frac{\partial\vec r_i}{\partial q_j}\delta q_j,$$
对完整约束，$\delta q_j$ 彼此独立，为后续变分提供自由度。

**意义**：广义坐标消除约束方程，把 $3N$ 个坐标降为 $n$ 个独立坐标。

**特例**：完整定常约束下 $n$ 为常数；非完整约束（如滚动）仅能写为微分形式不可积。

**应用**：机械系统建模、机器人运动学、约束优化。

#### 知识点4.2 理想约束与虚功原理

**定义**：虚位移 $\delta\vec r_i$ 为满足约束的无穷小、与时间冻结的位移；理想约束要求约束力虚功之和为零：
$$\sum_i\vec R_i\cdot\delta\vec r_i=0.$$

**假设**：约束理想（无摩擦、刚性）、$\delta t=0$。

**推导**：静平衡时每个质点 $\vec F_i+\vec R_i=0$，取虚功并求和
$$\sum_i(\vec F_i+\vec R_i)\cdot\delta\vec r_i=0.$$
由理想约束
$$\sum_i\vec F_i\cdot\delta\vec r_i=0.$$
用广义力 $Q_j=\sum_i\vec F_i\cdot\partial\vec r_i/\partial q_j$，得
$$\sum_{j=1}^nQ_j\delta q_j=0.$$
因 $\delta q_j$ 独立，故
$$Q_j=0,\quad j=1,\dots,n.$$

**意义**：用能量（虚功）而非力平衡描述平衡，自动消去约束力。

**特例**：无约束时退化为 $n=3N$ 个力平衡条件。

**应用**：复杂机构静力平衡、结构内力分析、虚位移法求约束力。

#### 知识点4.3 达朗贝尔原理与动力学普遍方程

**定义**：惯性力 $-\vec p_i$，达朗贝尔原理将动力学问题转化为“动态平衡”。

**假设**：惯性系、理想约束。

**推导**：牛顿第二定律 $\vec F_i+\vec R_i=\dot{\vec p}_i$，移项
$$\vec F_i+\vec R_i-\dot{\vec p}_i=0,$$
视 $-\dot{\vec p}_i$ 为惯性力，应用虚功原理：
$$\sum_i(\vec F_i-\dot{\vec p}_i)\cdot\delta\vec r_i+\sum_i\vec R_i\cdot\delta\vec r_i=0.$$
由理想约束
$$\sum_i(\vec F_i-\dot{\vec p}_i)\cdot\delta\vec r_i=0.$$
用广义坐标改写
$$\sum_{j=1}^n\left(Q_j-\frac{\mathrm d}{\mathrm dt}\frac{\partial T}{\partial\dot q_j}+\frac{\partial T}{\partial q_j}\right)\delta q_j=0.$$

**意义**：达朗贝尔–拉格朗日方程是分析力学的出发点，统一处理约束与惯性。

**特例**：静力学时 $\dot{\vec p}_i=0$，退化为虚功原理。

**应用**：多体动力学、离心惯性力分析、动静法校核。

#### 知识点4.4 广义力与势能的引入

**定义**：广义力 $Q_j=\sum_i\vec F_i\cdot\dfrac{\partial\vec r_i}{\partial q_j}$；若存在势 $V$ 使 $\vec F_i=-\nabla_iV$，则 $Q_j=-\partial V/\partial q_j$。

**假设**：力保守或可分离为保守部分与耗散部分。

**推导**：
$$Q_j=\sum_i\vec F_i\cdot\frac{\partial\vec r_i}{\partial q_j}=-\sum_i\nabla_iV\cdot\frac{\partial\vec r_i}{\partial q_j}=-\frac{\partial V}{\partial q_j}.$$
故达朗贝尔–拉格朗日方程化为
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial T}{\partial\dot q_j}-\frac{\partial T}{\partial q_j}=-\frac{\partial V}{\partial q_j},$$
即
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial(T-V)}{\partial\dot q_j}-\frac{\partial(T-V)}{\partial q_j}=0.$$

**意义**：势能仅依赖坐标，进入方程右端负梯度。

**特例**：速度相关势（电磁场）$U(q,\dot q)$ 时令 $L=T-U$ 仍成立。

**应用**：电磁系统、弹性结构、重力场机构。

#### 知识点4.5 广义坐标下的动能结构

**定义**：系统动能
$$T=\frac12\sum_i m_i\dot{\vec r}_i\cdot\dot{\vec r}_i.$$

**假设**：$\vec r_i=\vec r_i(q,t)$ 二次可微。

**推导**：
$$\dot{\vec r}_i=\sum_j\frac{\partial\vec r_i}{\partial q_j}\dot q_j+\frac{\partial\vec r_i}{\partial t}.$$
代入
$$T=\frac12\sum_{j,k}m_{jk}\dot q_j\dot q_k+\sum_j a_j\dot q_j+a_0,$$
其中
$$m_{jk}=\sum_i m_i\frac{\partial\vec r_i}{\partial q_j}\cdot\frac{\partial\vec r_i}{\partial q_k},\quad a_j=\sum_i m_i\frac{\partial\vec r_i}{\partial q_j}\cdot\frac{\partial\vec r_i}{\partial t},\quad a_0=\frac12\sum_i m_i\left|\frac{\partial\vec r_i}{\partial t}\right|^2.$$
即 $T=T_2+T_1+T_0$，$T_2$ 为广义速度的二次型，$T_1$ 为线性项，$T_0$ 为独立项。

**意义**：动能对广义速度的结构决定了系统的惯性耦合与运动耦合。

**特例**：定常约束 $\partial\vec r_i/\partial t=0$ 时 $T_1=T_0=0$，$T=T_2$。

**应用**：非定常约束机构（如变长度摆）建模。

#### 知识点4.6 虚功率原理与 Jourdain 原理

**定义**：Jourdain 变分取位置与时间固定、仅改变速度的虚速度 $\delta\dot{\vec r}_i$；虚功率原理为
$$\sum_i(\vec F_i-\dot{\vec p}_i)\cdot\delta\dot{\vec r}_i=0.$$

**假设**：约束理想、$\delta\dot{\vec r}_i$ 与约束的速率形式相容；允许存在非完整约束
$$\sum_j a_{lj}(q,t)\dot q_j+b_l(q,t)=0.$$

**推导**：由达朗贝尔方程对速度变分，注意 $\dot{\vec p}_i=m_i\ddot{\vec r}_i$ 与 $\delta\dot{\vec r}_i$ 无关，而
$$\delta\dot{\vec r}_i=\sum_j\frac{\partial\vec r_i}{\partial q_j}\delta\dot q_j.$$
代入得
$$\sum_j\left(Q_j-\frac{\mathrm d}{\mathrm dt}\frac{\partial T}{\partial\dot q_j}+\frac{\partial T}{\partial q_j}\right)\delta\dot q_j=0.$$
非完整约束在速度层面给出 $\sum_j a_{lj}\delta\dot q_j=0$，引入乘子 $\mu_l$，由 $\delta\dot q_j$ 的独立性得
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial T}{\partial\dot q_j}-\frac{\partial T}{\partial q_j}=Q_j+\sum_l\mu_la_{lj}.$$

**意义**：虚功率原理对非完整约束同样有效，是分析力学对非完整系统的推广。

**特例**：无约束时 $\mu_l=0$，退化为完整系统方程。

**应用**：滚动球、冰刀、欠驱动机械系统的建模。

### 第5章 拉格朗日力学

#### 知识点5.1 拉格朗日方程的推导

**定义**：拉格朗日函数 $L=T-V$，作用量 $S=\int_{t_1}^{t_2}L(q,\dot q,t)\,\mathrm dt$。

**假设**：完整、理想约束；$L$ 二次连续可微；端点变分为零 $\delta q(t_1)=\delta q(t_2)=0$。

**推导（由达朗贝尔）**：已知
$$\sum_j\left(\frac{\mathrm d}{\mathrm dt}\frac{\partial T}{\partial\dot q_j}-\frac{\partial T}{\partial q_j}-Q_j\right)\delta q_j=0.$$
代入 $Q_j=-\partial V/\partial q_j$（$\partial V/\partial\dot q_j=0$）：
$$\sum_j\left[\frac{\mathrm d}{\mathrm dt}\frac{\partial(T-V)}{\partial\dot q_j}-\frac{\partial(T-V)}{\partial q_j}\right]\delta q_j=0.$$
由 $\delta q_j$ 独立得
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}=0,\qquad j=1,\dots,n.$$

**推导（由 Hamilton 原理）**：设 $\delta S=0$，
$$\delta S=\int_{t_1}^{t_2}\sum_j\left(\frac{\partial L}{\partial q_j}\delta q_j+\frac{\partial L}{\partial\dot q_j}\delta\dot q_j\right)\mathrm dt.$$
分部积分
$$\int_{t_1}^{t_2}\frac{\partial L}{\partial\dot q_j}\delta\dot q_j\,\mathrm dt=\left[\frac{\partial L}{\partial\dot q_j}\delta q_j\right]_{t_1}^{t_2}-\int_{t_1}^{t_2}\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}\delta q_j\,\mathrm dt.$$
边界项为零，故
$$\delta S=\int_{t_1}^{t_2}\sum_j\left(\frac{\partial L}{\partial q_j}-\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}\right)\delta q_j\,\mathrm dt=0.$$
由变分基本引理即得 Euler–Lagrange 方程。

**意义**：动力学由标量函数 $L$ 与变分原理完全确定，坐标无关。

**特例**：$L$ 不显含某坐标时该坐标守恒（循环积分）。

**应用**：复杂机械、场论、连续体离散化。

#### 知识点5.2 拉格朗日方程的能量积分与循环积分

**定义**：循环坐标 $q_k$ 满足 $\partial L/\partial q_k=0$；广义动量 $p_k=\partial L/\partial\dot q_k$。

**假设**：$L$ 不显含 $q_k$ 或显含时间与否。

**推导（循环积分）**：由 Euler–Lagrange
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_k}=\frac{\partial L}{\partial q_k}=0\ \Rightarrow\ p_k=\text{const}.$$

**推导（能量积分）**：计算
$$\frac{\mathrm dL}{\mathrm dt}=\sum_j\frac{\partial L}{\partial q_j}\dot q_j+\sum_j\frac{\partial L}{\partial\dot q_j}\ddot q_j+\frac{\partial L}{\partial t}.$$
代入 Euler–Lagrange $\partial L/\partial q_j=\dot p_j$：
$$\frac{\mathrm dL}{\mathrm dt}=\sum_j\left(\dot p_j\dot q_j+p_j\ddot q_j\right)+\frac{\partial L}{\partial t}=\frac{\mathrm d}{\mathrm dt}\sum_jp_j\dot q_j+\frac{\partial L}{\partial t}.$$
故
$$\frac{\mathrm d}{\mathrm dt}\left(\sum_jp_j\dot q_j-L\right)=-\frac{\partial L}{\partial t}.$$
定义哈密顿量 $H=\sum_jp_j\dot q_j-L$，则 $\partial L/\partial t=0$ 时 $H=$ 常数，即能量守恒。

**意义**：对称性对应守恒律（Noether 原理的初等形式）。

**特例**：定常约束下 $H=T+V$。

**应用**：保守系统能量检验、守恒量识别、数值积分。

#### 知识点5.3 拉格朗日乘子与完整约束

**定义**：约束 $f_\alpha(q,t)=0$ 以乘子 $\lambda_\alpha$ 引入。

**假设**：约束完整、乘子为待定的约束力系数。

**推导**：对受约束系统，$\delta q$ 不独立，需保留约束变分。由
$$\sum_j\left(\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}\right)\delta q_j=0,\qquad \delta f_\alpha=\sum_j\frac{\partial f_\alpha}{\partial q_j}\delta q_j=0.$$
引入乘子构造
$$\sum_j\left(\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}-\sum_\alpha\lambda_\alpha\frac{\partial f_\alpha}{\partial q_j}\right)\delta q_j=0.$$
取 $\lambda_\alpha$ 使系数全为零：
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}=\sum_\alpha\lambda_\alpha\frac{\partial f_\alpha}{\partial q_j}.$$
右端即广义约束力 $Q_j^{c}$。

**意义**：拉格朗日乘子给出约束力，可用于求约束反力。

**特例**：无约束时 $\lambda_\alpha=0$。

**应用**：多体系统、接触约束、机器人力控。

#### 知识点5.4 陀螺力、耗散力与非保守广义力

**定义**：广义力 $Q_j=Q_j^{c}+Q_j^{d}+Q_j^{e}$，其中 $Q_j^{d}=-\partial\mathcal R/\partial\dot q_j$。

**假设**：$Q_j$ 可分离为位势力、速度线性耗散力、其余广义力。

**推导**：完整拉格朗日方程为
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}=Q_j^{e}+Q_j^{d}.$$
若存在陀螺势 $V_g$ 使 $Q_j^{g}=\sum_k g_{jk}\dot q_k$（$g_{jk}=-g_{kj}$），则陀螺力不做功。耗散函数
$$\mathcal R=\frac12\sum_{j,k}c_{jk}\dot q_j\dot q_k,$$
耗散功率
$$P_d=-\sum_jQ_j^d\dot q_j=\sum_{j,k}c_{jk}\dot q_j\dot q_k=2\mathcal R\ge0.$$

**意义**：耗散力使能量单调下降，陀螺力改变方向但不改变能量。

**特例**：仅有黏性阻尼时 $\mathcal R=\frac12c\dot q^2$。

**应用**：转子陀螺效应、阻尼结构、端口哈密顿系统。

#### 知识点5.5 相对运动中的拉格朗日函数

**定义**：在转动参考系中位矢 $\vec r'$，牵连速度 $\boldsymbol\Omega\times\vec r'$。

**假设**：参考系以常角速度 $\boldsymbol\Omega$ 旋转（可推广到非定常）。

**推导**：惯性系速度
$$\vec v=\vec v'+\boldsymbol\Omega\times\vec r',\qquad \vec v'=\dot{\vec r}'.$$
动能
$$T=\frac12m(\vec v'+\boldsymbol\Omega\times\vec r')^2=\frac12m|\vec v'|^2+m\vec v'\cdot(\boldsymbol\Omega\times\vec r')+\frac12m|\boldsymbol\Omega\times\vec r'|^2.$$
势能 $V=V(\vec r')$，故
$$L=\frac12m|\dot{\vec r}'|^2+m\dot{\vec r}'\cdot(\boldsymbol\Omega\times\vec r')+\frac12m|\boldsymbol\Omega\times\vec r'|^2-V(\vec r').$$
展开并利用 $\dot{\vec r}'\cdot(\boldsymbol\Omega\times\vec r')=\boldsymbol\Omega\cdot(\vec r'\times\dot{\vec r}')$。

**意义**：拉格朗日框架自动产生科氏与离心项。

**特例**：地表局部坐标系中得地球自转修正项。

**应用**：地球物理流体、陀螺罗盘、旋转机械。

#### 知识点5.6 守恒律与 Noether 定理

**定义**：单参数变换 $q_j\to q_j+\epsilon\,\eta_j$ 使作用量不变，则存在守恒量。

**假设**：$L$ 在变换下至多相差时间全导数。

**推导**：设
$$\delta L=\sum_j\frac{\partial L}{\partial q_j}\eta_j+\sum_j\frac{\partial L}{\partial\dot q_j}\dot\eta_j.$$
利用 Euler–Lagrange
$$\delta L=\sum_j\left(\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}\right)\eta_j+\sum_j\frac{\partial L}{\partial\dot q_j}\dot\eta_j=\frac{\mathrm d}{\mathrm dt}\sum_j\frac{\partial L}{\partial\dot q_j}\eta_j.$$
若 $\delta L=\dfrac{\mathrm dF}{\mathrm dt}$，则
$$\frac{\mathrm d}{\mathrm dt}\left(\sum_j\frac{\partial L}{\partial\dot q_j}\eta_j-F\right)=0,$$
守恒量 $I=\sum_jp_j\eta_j-F$。

**意义**：每一连续对称性对应一个守恒律（时间平移→能量、空间平移→动量、转动→角动量）。

**特例**：空间平移对称给出动量守恒；转动对称给角动量守恒。

**应用**：守恒量识别、对称性约化、几何力学。

### 第6章 哈密顿力学与正则方程

#### 知识点6.1 Legendre 变换与哈密顿函数

**定义**：广义动量 $p_j=\partial L/\partial\dot q_j$，哈密顿量
$$H(q,p,t)=\sum_jp_j\dot q_j-L(q,\dot q,t).$$

**假设**：Hessian $\partial^2L/\partial\dot q_j\partial\dot q_k$ 非奇异（可逆）。

**推导**：由 $p_j=\partial L/\partial\dot q_j$ 反解 $\dot q_j=\dot q_j(q,p,t)$。取 $H$ 的微分
$$\mathrm dH=\sum_j\left(\dot q_j\mathrm dp_j+p_j\mathrm d\dot q_j-\frac{\partial L}{\partial q_j}\mathrm dq_j-\frac{\partial L}{\partial\dot q_j}\mathrm d\dot q_j\right)-\frac{\partial L}{\partial t}\mathrm dt.$$
利用 $p_j=\partial L/\partial\dot q_j$ 消去 $\mathrm d\dot q_j$ 项，并用 $\dot p_j=\partial L/\partial q_j$：
$$\mathrm dH=\sum_j\left(\dot q_j\mathrm dp_j-\dot p_j\mathrm dq_j\right)-\frac{\partial L}{\partial t}\mathrm dt.$$
比较得
$$\dot q_j=\frac{\partial H}{\partial p_j},\qquad \dot p_j=-\frac{\partial H}{\partial q_j},\qquad \frac{\partial H}{\partial t}=-\frac{\partial L}{\partial t}.$$

**意义**：将 $2n$ 个一阶方程取代 $n$ 个二阶方程，相空间几何更清晰。

**特例**：$H=T_2+V$（定常、保守）时即总能量。

**应用**：辛积分、相空间分析、量子化对应。

#### 知识点6.2 哈密顿正则方程与相空间流

**定义**：相空间坐标 $z=(q,p)$，正则方程 $\dot z=J\nabla H$，其中
$$J=\begin{pmatrix}0&I\\-I&0\end{pmatrix}.$$

**假设**：$H$ 二次连续可微，正则系统。

**推导**：将上节两式合并
$$\begin{pmatrix}\dot q\\\dot p\end{pmatrix}=\begin{pmatrix}0&I\\-I&0\end{pmatrix}\begin{pmatrix}\partial H/\partial q\\\partial H/\partial p\end{pmatrix}.$$
求能量的时间变化
$$\dot H=\sum_j\left(\frac{\partial H}{\partial q_j}\dot q_j+\frac{\partial H}{\partial p_j}\dot p_j\right)+\frac{\partial H}{\partial t}=\sum_j\left(\frac{\partial H}{\partial q_j}\frac{\partial H}{\partial p_j}-\frac{\partial H}{\partial p_j}\frac{\partial H}{\partial q_j}\right)+\frac{\partial H}{\partial t}=\frac{\partial H}{\partial t}.$$
故不显含时间时 $H$ 守恒。又 $\nabla\cdot\dot z=\nabla\cdot(J\nabla H)=0$，即相流保体积（Liouville 定理）。

**意义**：相空间体积守恒、辛结构保持。

**特例**：一维谐振子 $H=p^2/(2m)+\frac12m\omega^2q^2$，相轨迹为椭圆。

**应用**：辛算法、统计力学系综、稳定性相图。

#### 知识点6.3 泊松括号与运动积分

**定义**：泊松括号
$$\{f,g\}=\sum_j\left(\frac{\partial f}{\partial q_j}\frac{\partial g}{\partial p_j}-\frac{\partial f}{\partial p_j}\frac{\partial g}{\partial q_j}\right).$$

**假设**：$f,g$ 在相空间可微。

**推导**：计算
$$\dot f=\sum_j\left(\frac{\partial f}{\partial q_j}\dot q_j+\frac{\partial f}{\partial p_j}\dot p_j\right)+\frac{\partial f}{\partial t}=\sum_j\left(\frac{\partial f}{\partial q_j}\frac{\partial H}{\partial p_j}-\frac{\partial f}{\partial p_j}\frac{\partial H}{\partial q_j}\right)+\frac{\partial f}{\partial t}=\{f,H\}+\frac{\partial f}{\partial t}.$$
故 $f$ 为运动积分当且仅当 $\{f,H\}+\partial f/\partial t=0$。

**性质**：反对称性 $\{f,g\}=-\{g,f\}$；双线性；Leibniz 律；Jacobi 恒等式
$$\{f,\{g,h\}\}+\{g,\{h,f\}\}+\{h,\{f,g\}\}=0.$$
基本括号 $\{q_j,q_k\}=0,\ \{p_j,p_k\}=0,\ \{q_j,p_k\}=\delta_{jk}$。

**意义**：泊松括号是经典力学到量子力学的对应桥梁。

**特例**：$\{f,H\}=0$ 且不显含时间则守恒。

**应用**：可积系统、守恒量搜寻、哈密顿控制。

#### 知识点6.4 正则变换与生成函数

**定义**：变换 $(q,p)\to(Q,P)$ 保持正则方程形式者称正则变换。

**假设**：变换光滑且可逆，$2n$ 维相空间。

**推导**：要求新哈密顿 $K(Q,P,t)$ 满足
$$\sum_jp_j\mathrm dq_j-H\mathrm dt=\sum_jP_j\mathrm dQ_j-K\mathrm dt+\mathrm dF,$$
其中 $F$ 为生成函数。取 $F_1(q,Q,t)$，则
$$p_j=\frac{\partial F_1}{\partial q_j},\qquad P_j=-\frac{\partial F_1}{\partial Q_j},\qquad K=H+\frac{\partial F_1}{\partial t}.$$
若取 $F_2(q,P,t)=F_1+\sum_jQ_jP_j$，则
$$p_j=\frac{\partial F_2}{\partial q_j},\qquad Q_j=\frac{\partial F_2}{\partial P_j}.$$
生成函数的选择对应不同的正则化方案。

**意义**：通过坐标选择简化哈密顿量、求积运动。

**特例**：恒等变换 $F_2=\sum_jq_jP_j$。

**应用**：作用–角变量、Hamilton–Jacobi 方程、摄动理论。

#### 知识点6.5 Hamilton–Jacobi 方程

**定义**：寻找生成函数 $S(q,P,t)$ 使新哈密顿 $K=0$。

**假设**：$H$ 不显含时间或可通过分离变量处理。

**推导**：由 $F_2=S$，有 $p_j=\partial S/\partial q_j$，$Q_j=\partial S/\partial P_j$，$K=H+\partial S/\partial t$。令 $K=0$：
$$\frac{\partial S}{\partial t}+H\left(q,\frac{\partial S}{\partial q},t\right)=0.$$
此即 Hamilton–Jacobi 方程。若 $H$ 不显含 $t$，令 $S=W(q)-Et$，得
$$H\left(q,\frac{\partial W}{\partial q}\right)=E,$$
其中 $W$ 为特征函数。

**意义**：将 $2n$ 个常微分方程化为一个一阶偏微分方程；$S$ 为作用量。

**特例**：可分离系统 $W=\sum_jW_j(q_j)$ 时逐坐标积分。

**应用**：解析力学求积、WKB 近似、量子化条件。

#### 知识点6.6 作用–角变量与可积系统

**定义**：对周期运动，作用变量
$$J_k=\oint p_k\,\mathrm dq_k,\qquad \theta_k=\frac{\partial W}{\partial J_k}.$$

**假设**：系统完全可积（存在 $n$ 个对合守恒量）、运动有界。

**推导**：以 $W=\sum_kW_k(q_k;J)$ 代入 HJ 方程，得 $H=H(J)$。因 $H$ 不依赖 $\theta$，
$$\dot J_k=-\frac{\partial H}{\partial\theta_k}=0,\qquad \dot\theta_k=\frac{\partial H}{\partial J_k}=\nu_k(J)=\text{const}.$$
故 $\theta_k=\nu_kt+\theta_{k0}$，运动在环面上线性化。

**意义**：可积系统的运动被约化到不变环面上的线性流动。

**特例**：谐振子 $J=E/\nu$，$\nu=\omega$。

**应用**：天体力学、粒子加速器、慢变系统的绝热不变量。

### 第7章 刚体动力学

#### 知识点7.1 惯性张量与角动量

**定义**：对体坐标原点，惯性张量
$$\mathbf I=\int_V\rho\bigl[(\vec r\cdot\vec r)\mathbf 1-\vec r\vec r\bigr]\mathrm dV,$$
分量 $I_{ij}=\int\rho(r^2\delta_{ij}-x_ix_j)\mathrm dV$。

**假设**：刚体连续、密度可积、体坐标固定于刚体。

**推导（平行轴定理）**：设质点相对质心 $\vec r=\vec r_c+\vec\rho$，$\int\rho\vec\rho\,\mathrm dV=0$。则
$$\mathbf I_O=\int\rho\bigl[(\vec r_c+\vec\rho)^2\mathbf 1-(\vec r_c+\vec\rho)(\vec r_c+\vec\rho)\bigr]\mathrm dV.$$
展开并利用质心项为零：
$$\mathbf I_O=\mathbf I_c+m\bigl[(\vec r_c\cdot\vec r_c)\mathbf 1-\vec r_c\vec r_c\bigr]=\mathbf I_c+m(r_c^2\mathbf 1-\vec r_c\vec r_c).$$

**意义**：角动量 $\vec L=\mathbf I\boldsymbol\omega$；$\mathbf I$ 对称且正定。

**特例**：主轴系中 $\mathbf I=\mathrm{diag}(I_1,I_2,I_3)$。

**应用**：航天器姿态、飞轮储能、机械转子。

#### 知识点7.2 欧拉动力学方程推导

**定义**：体坐标下角动量 $\vec L=\mathbf I\boldsymbol\omega$，力矩 $\vec M$。

**假设**：体坐标原点取质心或固定点，主轴系。

**推导**：由角动量定理 $\dot{\vec L}_O=\vec M$，使用输运定理在体坐标中对时间求导：
$$\left(\frac{\mathrm d\vec L}{\mathrm dt}\right)_{\text{惯性}}=\left(\frac{\mathrm d\vec L}{\mathrm dt}\right)_{\text{体}}+\boldsymbol\omega\times\vec L=\vec M.$$
在主轴系 $\vec L=(I_1\omega_1,I_2\omega_2,I_3\omega_3)$，且 $\mathbf I$ 在体系中为常张量，故
$$I_1\dot\omega_1+(I_3-I_2)\omega_2\omega_3=M_1,$$
$$I_2\dot\omega_2+(I_1-I_3)\omega_3\omega_1=M_2,$$
$$I_3\dot\omega_3+(I_2-I_1)\omega_1\omega_2=M_3.$$

**意义**：陀螺力矩项 $\boldsymbol\omega\times\vec L$ 是转动坐标非惯性效应的体现。

**特例**：球对称 $I_1=I_2=I_3$ 时 $\boldsymbol\omega=$ 常数（无力矩）。

**应用**：陀螺仪、自旋稳定、姿态机动。

#### 知识点7.3 无力矩对称陀螺（自由刚体）

**定义**：$I_1=I_2=I_\perp$，$I_3=I_\parallel$，$\vec M=0$。

**假设**：无外力矩、主轴系。

**推导**：欧拉方程化为
$$I_\perp\dot\omega_1=(I_\perp-I_3)\omega_2\omega_3,\quad I_\perp\dot\omega_2=(I_3-I_\perp)\omega_3\omega_1,\quad I_3\dot\omega_3=0.$$
故 $\omega_3=\text{const}$。令 $\Omega=\dfrac{I_3-I_\perp}{I_\perp}\omega_3$，则
$$\dot\omega_1=-\Omega\omega_2,\qquad \dot\omega_2=\Omega\omega_1.$$
解为
$$\omega_1=A\cos\Omega t+B\sin\Omega t,\qquad \omega_2=A\sin\Omega t-B\cos\Omega t,$$
即横向角速度以频率 $\Omega$ 绕 $3$ 轴进动。

**意义**：角速度在体坐标系中绕自旋轴做圆周运动（体锥）；在惯性系中角动量守恒，自旋轴绕角动量进动（空间锥）。

**特例**：$I_3>I_\perp$（扁球）$\Omega>0$；$I_3<I_\perp$（长球）反向。

**应用**：地球自由章动、卫星自旋稳定、网球拍定理。

#### 知识点7.4 拉格朗日陀螺与欧拉角动力学

**定义**：拉格朗日陀螺为定点、$I_1=I_2=I_\perp$、$I_3$ 的刚体，重力矩作用。

**假设**：定点，重力沿惯性竖直方向，$L$ 不显含方位角。

**推导**：用欧拉角 $(\psi,\theta,\phi)$，
$$T=\frac12I_\perp(\dot\theta^2+\dot\psi^2\sin^2\theta)+\frac12I_3(\dot\phi+\dot\psi\cos\theta)^2,$$
$$V=mgl\cos\theta.$$
$L=T-V$ 不显含 $\psi,\phi$，故
$$p_\psi=I_\perp\dot\psi\sin^2\theta+I_3(\dot\phi+\dot\psi\cos\theta)\cos\theta=\text{const},$$
$$p_\phi=I_3(\dot\phi+\dot\psi\cos\theta)=\text{const}.$$
由能量守恒与两个循环积分完全可积。

**意义**：陀螺的进动、章动、自转耦合由三个守恒量描述。

**特例**：快陀螺近似下自转频率远大于章动频率，得规则进动。

**应用**：陀螺罗经、炮弹自旋稳定、行星自转。

#### 知识点7.5 刚体动能与广义惯性

**定义**：刚体动能 $T=\frac12m v_c^2+\frac12\boldsymbol\omega^\top\mathbf I_c\boldsymbol\omega$。

**假设**：质心速度为 $v_c$，角速度为 $\boldsymbol\omega$。

**推导（König 定理）**：将 $\vec v_i=\vec v_c+\boldsymbol\omega\times\vec\rho_i$，则
$$T=\frac12\sum_im_i(\vec v_c+\boldsymbol\omega\times\vec\rho_i)^2$$
$$=\frac12Mv_c^2+\vec v_c\cdot\left(\boldsymbol\omega\times\sum_im_i\vec\rho_i\right)+\frac12\sum_im_i|\boldsymbol\omega\times\vec\rho_i|^2.$$
因 $\sum_im_i\vec\rho_i=0$，交叉项为零。又
$$\sum_im_i|\boldsymbol\omega\times\vec\rho_i|^2=\sum_im_i\bigl[\omega^2\rho_i^2-(\boldsymbol\omega\cdot\vec\rho_i)^2\bigr]=\boldsymbol\omega^\top\mathbf I_c\boldsymbol\omega.$$

**意义**：平动与转动动能解耦（以质心为参考点）。

**特例**：绕固定轴转动 $T=\frac12I\omega^2$。

**应用**：滚动机构、飞轮、机器人连杆。

#### 知识点7.6 刚体微振动与自由刚体稳定性

**定义**：中间主惯量轴的永久转动稳定性。

**假设**：无力矩自由刚体，$I_1<I_2<I_3$。

**推导**：设绕 $1$ 轴永久转动 $\omega_1=\Omega+\delta_1,\ \omega_2=\delta_2,\ \omega_3=\delta_3$。线性化欧拉方程：
$$\dot\delta_1=0,\qquad \dot\delta_2=\frac{I_3-I_1}{I_2}\Omega\delta_3,\qquad \dot\delta_3=\frac{I_1-I_2}{I_3}\Omega\delta_2.$$
对 $2$、$3$ 求导
$$\ddot\delta_2=\frac{(I_3-I_1)(I_1-I_2)}{I_2I_3}\Omega^2\delta_2.$$
绕 $1$ 轴稳定要求系数为负，即 $(I_3-I_1)(I_1-I_2)<0$，对 $I_1<I_2<I_3$ 恒成立。绕中间轴 $2$ 时对应系数为正，故不稳定。

**意义**：自由刚体绕最大、最小惯量轴稳定，绕中间轴不稳定。

**特例**：网球拍翻转即绕中间轴不稳定运动。

**应用**：卫星姿态稳定性设计、空间碎片翻滚预测。

### 第8章 连续介质力学与应力应变

#### 知识点8.1 连续介质假设与应力张量

**定义**：应力张量 $\boldsymbol\sigma$ 由面力 $\vec t(\vec n)=\boldsymbol\sigma\vec n$ 定义。

**假设**：连续介质假设；Cauchy 应力原理（面力仅依赖法向 $\vec n$）；局部平衡。

**推导（Cauchy 四面体）**：取四面体，斜面法向 $\vec n$，面积为 $A$，其余三面面积 $A_i=n_iA$。体积力 $\vec f$，惯性 $\rho\vec a$。平衡
$$\vec t(\vec n)A+\sum_i\vec t(-\hat e_i)A_i+\rho\vec b V=\rho\vec a V.$$
令 $V/A\to0$（特征尺寸趋于零），面力项占优，利用 $\vec t(-\hat e_i)=-\vec t(\hat e_i)$，得
$$t_j(\vec n)=n_i\,\sigma_{ij},$$
即 $\vec t=\boldsymbol\sigma^\top\vec n$。对非极性介质 $\boldsymbol\sigma$ 对称，故 $\vec t=\boldsymbol\sigma\vec n$。

**意义**：一点的应力状态由张量完全描述。

**特例**：静水压 $\boldsymbol\sigma=-p\mathbf 1$，$\vec t=-p\vec n$。

**应用**：结构强度、土力学、流体黏性应力。

#### 知识点8.2 应力平衡方程与运动方程

**定义**：应力分量 $\sigma_{ij}$，体积力 $b_i$，加速度 $a_i$。

**假设**：连续、动量守恒、无体力偶。

**推导**：对任意体积 $V$（表面 $S$）动量定理
$$\frac{\mathrm d}{\mathrm dt}\int_V\rho v_i\,\mathrm dV=\int_V\rho b_i\,\mathrm dV+\oint_S\sigma_{ij}n_j\,\mathrm dS.$$
对面积分用散度定理
$$\oint_S\sigma_{ij}n_j\,\mathrm dS=\int_V\frac{\partial\sigma_{ij}}{\partial x_j}\,\mathrm dV.$$
对左端用 Reynolds 输运定理（$\rho\,\mathrm dV$ 守恒）：
$$\int_V\rho\frac{\mathrm dv_i}{\mathrm dt}\,\mathrm dV=\int_V\left(\rho b_i+\frac{\partial\sigma_{ij}}{\partial x_j}\right)\mathrm dV.$$
由体积任意性
$$\rho\frac{\mathrm dv_i}{\mathrm dt}=\rho b_i+\frac{\partial\sigma_{ij}}{\partial x_j},$$
矢量为 $\rho\dfrac{\mathrm d\vec v}{\mathrm dt}=\rho\vec b+\nabla\cdot\boldsymbol\sigma$。静力学时退化为 $\nabla\cdot\boldsymbol\sigma+\rho\vec b=0$。

**意义**：局部动量守恒，是固体与流体共同运动方程。

**特例**：流体中 $\sigma_{ij}=-p\delta_{ij}+\tau_{ij}$。

**应用**：有限元应力分析、流体 N–S 推导。

#### 知识点8.3 应变张量与几何方程

**定义**：位移场 $\vec u$，小应变张量
$$\varepsilon_{ij}=\frac12\!\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\right).$$

**假设**：小变形、小转动（$\|\nabla\vec u\|\ll1$）。

**推导**：取两邻近点 $\mathrm d\vec x$ 与 $\mathrm d\vec x+\mathrm d\vec u$，长度平方变化
$$\mathrm d\vec x'^2-\mathrm d\vec x^2=2\varepsilon_{ij}\mathrm dx_i\mathrm dx_j,$$
其中 $\varepsilon_{ij}=\frac12(\partial_ju_i+\partial_iu_j)$。正应变
$$\varepsilon_{xx}=\frac{\partial u_x}{\partial x},$$
工程剪应变
$$\gamma_{xy}=2\varepsilon_{xy}=\frac{\partial u_x}{\partial y}+\frac{\partial u_y}{\partial x}.$$

**意义**：对称化以消去刚体转动，纯变形由对称部分描述。

**特例**：刚体运动 $\vec u=\vec u_0+\boldsymbol\Omega\times\vec x$ 时 $\boldsymbol\varepsilon=0$。

**应用**：应变片测量、有限元插值、热应变。

#### 知识点8.4 相容方程

**定义**：应变场可积性条件（保证存在单值连续位移）。

**假设**：小应变、单连通区域。

**推导**：由定义
$$\varepsilon_{ij,k}+\varepsilon_{jk,i}-\varepsilon_{ik,j}-\varepsilon_{ji,k}=2\partial_i\partial_k u_j+2\partial_j\partial_i\cdots$$
消项后得 Saint-Venant 相容方程：
$$\varepsilon_{ij,kl}+\varepsilon_{kl,ij}-\varepsilon_{ik,jl}-\varepsilon_{jl,ik}=0.$$
二维情形可化为单一方程
$$\frac{\partial^2\varepsilon_{xx}}{\partial y^2}+\frac{\partial^2\varepsilon_{yy}}{\partial x^2}=2\frac{\partial^2\varepsilon_{xy}}{\partial x\partial y}.$$

**意义**：任意给定应变场不一定来自位移场，相容方程给出约束。

**特例**：单连通区域需全部相容方程；多连通需附加回路条件。

**应用**：应力函数方法、残余应力、变形协调。

#### 知识点8.5 主应力、应力不变量与莫尔圆

**定义**：主应力 $\sigma$ 满足 $(\boldsymbol\sigma-\sigma\mathbf 1)\vec n=0$，特征方程
$$\sigma^3-I_1\sigma^2+I_2\sigma-I_3=0,$$
$$I_1=\mathrm{tr}\boldsymbol\sigma,\quad I_2=\tfrac12\bigl[(\mathrm{tr}\boldsymbol\sigma)^2-\mathrm{tr}\boldsymbol\sigma^2\bigr],\quad I_3=\det\boldsymbol\sigma.$$

**假设**：$\boldsymbol\sigma$ 对称实张量。

**推导**：存在 $\vec n$ 使 $\boldsymbol\sigma\vec n=\sigma\vec n$，即 $(\boldsymbol\sigma-\sigma\mathbf 1)\vec n=0$。非零解要求行列式为零，展开得特征方程。对称矩阵特征值实、特征向量正交。主剪应力为
$$\tau_{\max}=\frac{\sigma_1-\sigma_3}{2}.$$
莫尔圆参数：圆心 $\frac{\sigma_1+\sigma_3}{2}$，半径 $\frac{\sigma_1-\sigma_3}{2}$。

**意义**：任意斜面应力轨迹为三个莫尔圆围成的区域。

**特例**：单向拉伸 $\sigma_1=\sigma$，$\sigma_2=\sigma_3=0$，$\tau_{\max}=\sigma/2$（45° 面）。

**应用**：强度准则、断裂面预测、压力容器设计。

#### 知识点8.6 大变形与变形梯度简介

**定义**：变形梯度 $\mathbf F=\partial\vec x/\partial\vec X$，右 Cauchy–Green 张量 $\mathbf C=\mathbf F^\top\mathbf F$，Green–Lagrange 应变 $\mathbf E=\frac12(\mathbf C-\mathbf 1)$。

**假设**：连续、可逆变形 $\det\mathbf F>0$。

**推导**：微元变换 $\mathrm d\vec x=\mathbf F\,\mathrm d\vec X$。长度平方
$$\mathrm d\vec x^2-\mathrm d\vec X^2=\mathrm d\vec X\cdot(\mathbf F^\top\mathbf F-\mathbf 1)\mathrm d\vec X=2\,\mathrm d\vec X\cdot\mathbf E\,\mathrm d\vec X.$$
在 $\|\nabla\vec u\|\ll1$ 时 $\mathbf F=\mathbf 1+\nabla\vec u$，
$$\mathbf E=\frac12(\nabla\vec u+\nabla\vec u^\top+\nabla\vec u^\top\nabla\vec u)\approx\boldsymbol\varepsilon.$$

**意义**：有限变形下应变的客观度量，剔除刚体转动。

**特例**：不可伸长体 $\det\mathbf F=1$（等容）。

**应用**：橡胶、软组织、大变形有限元。

### 第9章 弹性力学与能量方法

#### 知识点9.1 广义胡克定律与本构关系

**定义**：线弹性本构
$$\sigma_{ij}=C_{ijkl}\varepsilon_{kl}.$$

**假设**：小变形、线弹性、等温或绝热、材料均匀。

**推导（各向同性）**：由各向同性要求 $C_{ijkl}$ 为四阶各向同性张量，其一般形式
$$C_{ijkl}=\lambda\delta_{ij}\delta_{kl}+\mu(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}).$$
于是
$$\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij},$$
其中 $\lambda,\mu$ 为 Lamé 常数。求迹
$$\sigma_{kk}=(3\lambda+2\mu)\varepsilon_{kk}.$$
反解
$$\varepsilon_{ij}=\frac{1+\nu}{E}\sigma_{ij}-\frac{\nu}{E}\sigma_{kk}\delta_{ij},$$
其中
$$\mu=G=\frac{E}{2(1+\nu)},\qquad \lambda=\frac{E\nu}{(1+\nu)(1-2\nu)}.$$

**意义**：应力由应变唯一决定；$E$ 为杨氏模量、$\nu$ 为泊松比。

**特例**：不可压缩 $\nu\to1/2$ 时 $\lambda\to\infty$，需用体积模量表述。

**应用**：结构刚度、材料参数标定、热弹性。

#### 知识点9.2 Navier 位移方程

**定义**：以位移 $\vec u$ 为未知量的弹性静力学方程。

**假设**：各向同性线弹性、无体力（或含体力）。

**推导**：将本构代入平衡方程
$$\sigma_{ij,j}=\lambda\varepsilon_{kk,j}\delta_{ij}+2\mu\varepsilon_{ij,j}=0.$$
利用
$$\varepsilon_{kk,j}=\partial_j\partial_k u_k,\qquad \varepsilon_{ij,j}=\frac12(\partial_j\partial_j u_i+\partial_j\partial_i u_j)=\frac12\nabla^2u_i+\frac12\partial_i\partial_ju_j,$$
得
$$\lambda\partial_i(\nabla\cdot\vec u)+\mu\nabla^2u_i+\mu\partial_i(\nabla\cdot\vec u)+\rho b_i=0,$$
即
$$(\lambda+\mu)\nabla(\nabla\cdot\vec u)+\mu\nabla^2\vec u+\rho\vec b=0.$$

**意义**：将弹性问题化为位移边值问题求解。

**特例**：$\nabla\cdot\vec u=0$（等容）时退化为 $\mu\nabla^2\vec u=0$。

**应用**：有限元位移法、弹性波、断裂力学。

#### 知识点9.3 平面问题与 Airy 应力函数

**定义**：平面应力或平面应变，应力函数 $\phi(x,y)$：
$$\sigma_{xx}=\frac{\partial^2\phi}{\partial y^2},\qquad \sigma_{yy}=\frac{\partial^2\phi}{\partial x^2},\qquad \sigma_{xy}=-\frac{\partial^2\phi}{\partial x\partial y}.$$

**假设**：无体力、平面问题、单连通。

**推导**：无体力平衡方程
$$\frac{\partial\sigma_{xx}}{\partial x}+\frac{\partial\sigma_{xy}}{\partial y}=0,\qquad \frac{\partial\sigma_{xy}}{\partial x}+\frac{\partial\sigma_{yy}}{\partial y}=0.$$
第一个方程存在 $\phi$ 使 $\sigma_{xx}=\partial_y^2\phi$，$\sigma_{xy}=-\partial_x\partial_y\phi$；代入第二式得 $\sigma_{yy}=\partial_x^2\phi$。再由相容方程与各向同性本构可得
$$\nabla^4\phi=\left(\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}\right)^2\phi=0.$$

**意义**：弹性平面问题化为求解双调和方程。

**特例**：多项式解对应常应力场；极坐标下用于孔洞问题。

**应用**：孔边应力集中、裂纹尖端场、接触问题。

#### 知识点9.4 应变能与 Castigliano 定理

**定义**：应变能密度 $W=\frac12\sigma_{ij}\varepsilon_{ij}$，总应变能 $U=\int_VW\,\mathrm dV$。

**假设**：线弹性、保守（无耗散）、等温。

**推导**：线弹性下 $\sigma_{ij}=C_{ijkl}\varepsilon_{kl}$，
$$W=\frac12C_{ijkl}\varepsilon_{ij}\varepsilon_{kl},$$
且 $\sigma_{ij}=\partial W/\partial\varepsilon_{ij}$。对结构，外载荷做功 $W_e=\frac12\sum_iP_i\Delta_i$。由能量守恒 $U=W_e$。对载荷 $P_k$ 求偏导：
$$\Delta_k=\frac{\partial U}{\partial P_k}.$$
（Castigliano 第一定理）类似地，位移约束下 $P_k=\partial U^*/\partial\Delta_k$，$U^*$ 为余能。

**意义**：位移（载荷）可由能量对载荷（位移）的导数得到。

**特例**：静定桁架 $U=\sum\frac{N_i^2L_i}{2E_iA_i}$。

**应用**：结构位移计算、超静定、柔度矩阵。

#### 知识点9.5 最小势能原理与虚功

**定义**：总势能 $\Pi=U-W_e$。

**假设**：线弹性、保守载荷、位移边界满足。

**推导**：真实位移 $\vec u$，许可变分 $\delta\vec u$（满足位移边界）。虚功方程
$$\int_V\sigma_{ij}\delta\varepsilon_{ij}\mathrm dV=\int_Vb_i\delta u_i\,\mathrm dV+\int_{S_t}t_i\delta u_i\,\mathrm dS.$$
左端 $\delta U$，右端 $\delta W_e$，故 $\delta(U-W_e)=0$，即
$$\delta\Pi=0,\qquad \Pi=U-W_e.$$
其二阶变分 $\delta^2U=\int\delta\varepsilon^\top C\delta\varepsilon\,\mathrm dV>0$（$C$ 正定），故为极小值。

**意义**：弹性静力学等价于总势能极小；有限元由之离散。

**特例**：无载荷时 $U$ 最小。

**应用**：Ritz 法、Galerkin 法、有限元刚度阵。

#### 知识点9.6 最小余能原理与虚应力

**定义**：余能 $U^*=\int_VW^*\mathrm dV$，$W^*=\frac12\sigma_{ij}\varepsilon_{ij}$（线弹性下 $W^*=W$）。

**假设**：应力场满足平衡与力边界，变分 $\delta\sigma_{ij}$ 满足齐次平衡。

**推导**：由 $\delta\varepsilon_{ij}=\delta S_{ijkl}\sigma_{kl}$（$S$ 为柔度），
$$\int_V\varepsilon_{ij}\delta\sigma_{ij}\mathrm dV=\int_{S_u}u_i\delta t_i\,\mathrm dS,$$
故 $\delta U^*-\delta W_e^*=0$，其中 $W_e^*=\int_{S_u}u_i t_i\,\mathrm dS$。即 $\delta\Pi^*=0$，$\Pi^*=U^*-W_e^*$，且为极小。

**意义**：应力法求解弹性问题的能量基础。

**特例**：给定载荷的静定结构无余能约束。

**应用**：超静定结构分析、应力场近似、断裂能量释放率。

### 第10章 流体力学基础与 Navier–Stokes 方程

#### 知识点10.1 连续方程（质量守恒）

**定义**：密度 $\rho$，速度 $\vec v$。

**假设**：连续介质、无源无汇。

**推导**：对任意控制体
$$\frac{\mathrm d}{\mathrm dt}\int_V\rho\,\mathrm dV=0.$$
用 Reynolds 输运定理
$$\int_V\left(\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\vec v)\right)\mathrm dV=0.$$
由体积任意性
$$\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\vec v)=0.$$
展开并代入物质导数
$$\frac{\mathrm D\rho}{\mathrm D t}+\rho\nabla\cdot\vec v=0,\qquad \frac{\mathrm D}{\mathrm Dt}=\frac{\partial}{\partial t}+\vec v\cdot\nabla.$$

**意义**：密度与速度散度互相制约；不可压缩时 $\nabla\cdot\vec v=0$。

**特例**：定常、不可压缩流动 $\nabla\cdot(\rho\vec v)=0$，沿流管 $\rho Av=$ 常数。

**应用**：管道流量、激波关系、多相流。

#### 知识点10.2 欧拉方程与应力表述

**定义**：理想流体 $\sigma_{ij}=-p\delta_{ij}$。

**假设**：无黏、连续、牛顿第二定律。

**推导**：由连续介质动量方程
$$\rho\frac{\mathrm Dv_i}{\mathrm Dt}=\rho b_i+\frac{\partial\sigma_{ij}}{\partial x_j}.$$
代入 $\sigma_{ij}=-p\delta_{ij}$：
$$\rho\frac{\mathrm D\vec v}{\mathrm Dt}=-\nabla p+\rho\vec b,$$
展开
$$\rho\left(\frac{\partial\vec v}{\partial t}+(\vec v\cdot\nabla)\vec v\right)=-\nabla p+\rho\vec b.$$

**意义**：无黏流动的动量方程，非线性对流项是湍流与激波的根源。

**特例**：静力学 $\vec v=0$ 时 $\nabla p=\rho\vec b$（流体静压）。

**应用**：气动初步设计、大气运动、喷管流动。

#### 知识点10.3 Navier–Stokes 方程推导

**定义**：牛顿流体 $\tau_{ij}=\mu(\partial_ju_i+\partial_iu_j)+\lambda(\nabla\cdot\vec v)\delta_{ij}$。

**假设**：牛顿黏性、Stokes 假设 $\lambda=-\frac23\mu$、常黏度。

**推导**：总应力 $\sigma_{ij}=-p\delta_{ij}+\tau_{ij}$，代入动量方程
$$\rho\frac{\mathrm Dv_i}{\mathrm Dt}=-\partial_ip+\partial_j\tau_{ij}+\rho b_i.$$
计算
$$\partial_j\tau_{ij}=\mu\partial_j(\partial_ju_i+\partial_iu_j)+\lambda\partial_i(\partial_kv_k)=\mu\nabla^2v_i+(\mu+\lambda)\partial_i(\nabla\cdot\vec v).$$
故
$$\rho\left(\frac{\partial\vec v}{\partial t}+(\vec v\cdot\nabla)\vec v\right)=-\nabla p+\mu\nabla^2\vec v+\left(\mu+\lambda\right)\nabla(\nabla\cdot\vec v)+\rho\vec b.$$
不可压缩时
$$\rho\left(\frac{\partial\vec v}{\partial t}+(\vec v\cdot\nabla)\vec v\right)=-\nabla p+\mu\nabla^2\vec v+\rho\vec b.$$

**意义**：黏性耗散由 $\mu\nabla^2\vec v$ 体现；N–S 是湍流与边界层基础。

**特例**：$\mu=0$ 退化为欧拉方程；Stokes 流（$\mathrm{Re}\ll1$）忽略惯性。

**应用**：管道流动、边界层、CFD。

#### 知识点10.4 伯努利方程与沿流线积分

**定义**：定常、无黏、不可压缩、沿流线或势流。

**假设**：定常、正压、体积力有势、无黏。

**推导**：欧拉方程
$$\rho(\vec v\cdot\nabla)\vec v=-\nabla p-\rho\nabla\Phi.$$
用恒等式 $(\vec v\cdot\nabla)\vec v=\nabla(\tfrac12v^2)-\vec v\times(\nabla\times\vec v)$，点乘 $\mathrm d\vec r=\vec v\,\mathrm dt$：
$$\rho\left[\nabla(\tfrac12v^2)-\vec v\times\boldsymbol\omega\right]\cdot\mathrm d\vec r=-\nabla p\cdot\mathrm d\vec r-\rho\nabla\Phi\cdot\mathrm d\vec r.$$
沿流线 $\vec v\times\boldsymbol\omega\perp\vec v\parallel\mathrm d\vec r$，交叉项为零，得
$$\mathrm d\!\left(\frac12\rho v^2+p+\rho\Phi\right)=0,$$
即
$$\frac12\rho v^2+p+\rho\Phi=\text{const}.$$

**意义**：沿流线机械能守恒；压力能、动能、位能互换。

**特例**：无旋流场中常数全流场相同；不可压缩重力场 $p+\frac12\rho v^2+\rho gz=\text{const}$。

**应用**：文丘里流量计、机翼升力（配合环量）、驻点压力。

#### 知识点10.5 涡量与涡量方程

**定义**：涡量 $\boldsymbol\omega=\nabla\times\vec v$。

**假设**：正压、体积力有势、黏度常。

**推导**：对 N–S 取旋度。利用矢量恒等式
$$\nabla\times(\vec v\cdot\nabla\vec v)=\nabla\times(\boldsymbol\omega\times\vec v),$$
以及 $\nabla\times\nabla(\cdot)=0$，得
$$\frac{\partial\boldsymbol\omega}{\partial t}+\nabla\times(\boldsymbol\omega\times\vec v)=\frac{1}{\rho^2}\nabla\rho\times\nabla p+\nu\nabla^2\boldsymbol\omega.$$
展开 $\nabla\times(\boldsymbol\omega\times\vec v)=(\vec v\cdot\nabla)\boldsymbol\omega-(\boldsymbol\omega\cdot\nabla)\vec v+\boldsymbol\omega(\nabla\cdot\vec v)$，并利用连续方程化简得
$$\frac{\mathrm D\boldsymbol\omega}{\mathrm Dt}=(\boldsymbol\omega\cdot\nabla)\vec v-\boldsymbol\omega(\nabla\cdot\vec v)+\frac{\nabla\rho\times\nabla p}{\rho^2}+\nu\nabla^2\boldsymbol\omega.$$
不可压缩、正压、无黏时
$$\frac{\mathrm D\boldsymbol\omega}{\mathrm Dt}=(\boldsymbol\omega\cdot\nabla)\vec v.$$

**意义**：涡管因速度梯度被拉伸而强化（涡拉伸）；黏性与斜压项是涡量产生/耗散来源。

**特例**：平面流动 $\boldsymbol\omega\perp\nabla$，拉伸项为零。

**应用**：湍流、旋涡脱落、大气气旋。

#### 知识点10.6 势流、流函数与复势

**定义**：无旋流动 $\vec v=\nabla\phi$，不可压缩时 $\nabla^2\phi=0$；二维流函数 $\psi$ 满足 $u=\partial\psi/\partial y,\ v=-\partial\psi/\partial x$。

**假设**：无旋、不可压缩、二维（复势）。

**推导**：不可压缩连续方程 $\partial_xu+\partial_yv=0$ 允许引入 $\psi$。无旋条件 $\partial_xv-\partial_yu=0$ 给出 $\nabla^2\psi=0$。复势
$$w(z)=\phi+i\psi,\qquad z=x+iy,$$
则
$$\frac{\mathrm dw}{\mathrm dz}=u-iv.$$
Cauchy–Riemann 方程保证 $\phi,\psi$ 调和共轭。

**意义**：势流由 Laplace 方程控制，可叠加；流函数等值线为流线。

**特例**：均匀流 $w=Uz$；点源 $w=\frac{Q}{2\pi}\ln z$；偶极子、圆柱绕流叠加。

**应用**：机翼环量、绕流升力、渗流。

### 第11章 守恒律、涡动力学与变分原理

#### 知识点11.1 Reynolds 输运定理

**定义**：系统（物质体）物理量变化与控制体变化的关系。

**假设**：连续介质、可微速度场。

**推导**：对无系统源的积分量 $I=\int_{V_m}f\,\mathrm dV$，
$$\frac{\mathrm dI}{\mathrm dt}=\frac{\mathrm d}{\mathrm dt}\int_{V(t)}f\,\mathrm dV.$$
将 $V(t)$ 映射到参考体 $V_0$，$\mathrm dV=J\,\mathrm dV_0$，$\dot J=J\nabla\cdot\vec v$，故
$$\frac{\mathrm dI}{\mathrm dt}=\int_{V_0}\left(\dot fJ+f\dot J\right)\mathrm dV_0=\int_V\left(\frac{\mathrm Df}{\mathrm Dt}+f\nabla\cdot\vec v\right)\mathrm dV.$$
又由连续方程 $\mathrm D\rho/\mathrm Dt+\rho\nabla\cdot\vec v=0$，得
$$\frac{\mathrm d}{\mathrm dt}\int_V\rho b\,\mathrm dV=\int_V\rho\frac{\mathrm Db}{\mathrm Dt}\,\mathrm dV,$$
$$\frac{\mathrm d}{\mathrm dt}\int_V f\,\mathrm dV=\int_V\frac{\partial f}{\partial t}\mathrm dV+\oint_Sf\vec v\cdot\vec n\,\mathrm dS.$$

**意义**：联系物质导数与控制体形式，是守恒律的统一工具。

**特例**：定常流动中系统量的变化完全由通量决定。

**应用**：质量、动量、能量、熵守恒律推导。

#### 知识点11.2 Kelvin 环量定理与 Helmholtz 涡定理

**定义**：环量 $\Gamma=\oint_C\vec v\cdot\mathrm d\vec r$。

**假设**：无黏、正压、体积力有势、物质曲线。

**推导**：环量对时间求导，
$$\frac{\mathrm D\Gamma}{\mathrm Dt}=\oint_C\frac{\mathrm D\vec v}{\mathrm Dt}\cdot\mathrm d\vec r+\oint_C\vec v\cdot\frac{\mathrm D(\mathrm d\vec r)}{\mathrm Dt}.$$
第二项 $\oint\vec v\cdot\mathrm d\vec v=\oint\mathrm d(\tfrac12v^2)=0$（闭曲线）。用欧拉方程
$$\frac{\mathrm D\vec v}{\mathrm Dt}=-\frac{\nabla p}{\rho}-\nabla\Phi,$$
得
$$\frac{\mathrm D\Gamma}{\mathrm Dt}=-\oint_C\frac{\nabla p}{\rho}\cdot\mathrm d\vec r.$$
正压时 $\nabla p/\rho=\nabla\!\int\mathrm dp/\rho$ 为全微分，闭曲线积分为零，故 $\mathrm D\Gamma/\mathrm Dt=0$。

**Helmholtz 定理**：涡管强度守恒；涡管不能在流体内部终止；涡线随流体运动。

**意义**：无黏正压流动中涡量不产生、不消灭。

**特例**：从静止启动的流动保持无旋。

**应用**：机翼启动涡、环量升力、超流。

#### 知识点11.3 守恒律与 Noether 定理（场论形式）

**定义**：场 $\phi^a(x)$，拉格朗日密度 $\mathcal L(\phi^a,\partial_\mu\phi^a)$，作用量
$$S=\int\mathcal L\,\mathrm d^4x.$$

**假设**：局部场论、变分原理 $\delta S=0$。

**推导**：Euler–Lagrange 方程
$$\partial_\mu\frac{\partial\mathcal L}{\partial(\partial_\mu\phi^a)}-\frac{\partial\mathcal L}{\partial\phi^a}=0.$$
对连续对称变换 $\phi^a\to\phi^a+\epsilon\,\delta\phi^a$，若 $\delta\mathcal L=\partial_\mu(\epsilon K^\mu)$，定义
$$j^\mu=\frac{\partial\mathcal L}{\partial(\partial_\mu\phi^a)}\delta\phi^a-K^\mu.$$
则
$$\partial_\mu j^\mu=0,$$
即守恒流。对时间平移得能量–动量张量，对内部对称得电荷等守恒流。

**意义**：对称性与守恒律的普适对应。

**特例**：实标量场平移对称给出动量守恒；复场相位对称给出电荷守恒。

**应用**：量子场论、连续介质守恒律、变分断裂。

#### 知识点11.4 Hamilton 原理在连续体中的形式

**定义**：连续体位移 $\vec u(\vec x,t)$，作用量
$$S=\int_{t_1}^{t_2}\int_V(T-V)\,\mathrm dV\,\mathrm dt.$$

**假设**：超弹性、小变形或有限变形、保守。

**推导**：以弹性体为例，$T=\frac12\rho\dot u_i\dot u_i$，$V=\frac12C_{ijkl}\varepsilon_{ij}\varepsilon_{kl}$。令 $\delta S=0$，
$$\delta S=\int\!\!\int\left(\rho\dot u_i\delta\dot u_i-C_{ijkl}\varepsilon_{ij}\delta\varepsilon_{kl}\right)\mathrm dV\mathrm dt.$$
分部积分（时间与空间），利用边界条件，得
$$\rho\ddot u_i-\partial_j(C_{ijkl}\varepsilon_{kl})=0,$$
即弹性动力学方程。对流体可得 Euler–Lagrange 形式。

**意义**：连续介质动力学可由作用量泛函驻值统一描述。

**特例**：静力学时退化为最小势能原理。

**应用**：有限元变分基础、多物理场耦合、几何力学。

#### 知识点11.5 稳定性与能量方法（Lagrange–Dirichlet）

**定义**：平衡 $\vec x_e$；若任意邻域内初始扰动下运动始终有界且趋于平衡，则稳定/渐近稳定。

**假设**：系统保守或有阻尼，势能可微。

**推导（Lagrange–Dirichlet）**：设 $V$ 在 $\vec q_e$ 处取严格局部极小，$V(\vec q_e)=0$，$V(\vec q)>0$。取能量
$$E=T+V.$$
对小扰动 $E$ 有上界（$V$ 局部正定），由能量守恒 $E=$ 常数且 $T\ge0$ 得 $V\le E$，故运动保持在 $V$ 的等值线内，稳定。加阻尼后 $E$ 单调下降，渐近稳定。

**意义**：势能极小判据是力学稳定性最基础结果。

**特例**：$V$ 极大时不稳定（Lyapunov 不稳定性）。

**应用**：结构屈曲、陀螺稳定、航天器位形稳定。

#### 知识点11.6 能量方程、熵方程与耗散

**定义**：单位质量内能 $e$，总能 $E=\int_V\rho(e+\tfrac12v^2)\mathrm dV$；热流 $\vec q$，体热源 $r$。

**假设**：连续介质、局部热力学平衡、Clausius–Duhem 不等式成立。

**推导**：对控制体应用热力学第一定律
$$\frac{\mathrm d}{\mathrm dt}\int_V\rho\left(e+\tfrac12v^2\right)\mathrm dV=\int_V\rho(\vec b\cdot\vec v+r)\mathrm dV+\oint_S(\vec t\cdot\vec v-\vec q\cdot\vec n)\mathrm dS.$$
减去机械能方程（动量方程点乘 $\vec v$ 并积分）：
$$\rho\frac{\mathrm De}{\mathrm Dt}=\sigma_{ij}\frac{\partial v_i}{\partial x_j}-\nabla\cdot\vec q+\rho r.$$
对牛顿流体，$\sigma_{ij}\partial_jv_i=-p\nabla\cdot\vec v+\Phi$，其中 $\Phi=2\mu(\dot\varepsilon_{ij}\dot\varepsilon_{ij}-\tfrac13(\nabla\cdot\vec v)^2)\ge0$ 为黏性耗散。Fourier 定律 $\vec q=-k\nabla T$。第二定律给出
$$\rho\frac{\mathrm Ds}{\mathrm Dt}+\nabla\cdot\frac{\vec q}{T}-\frac{\rho r}{T}=\frac{\Phi}{T}+\frac{k|\nabla T|^2}{T^2}\ge0.$$

**意义**：机械能耗散为热，熵产生率由耗散与热传导不可逆性决定。

**特例**：理想流体 $\Phi=0$，绝热时可逆。

**应用**：可压缩流、燃烧、热弹性、耗散结构。

### 第12章 稳定性理论基础

#### 知识点12.1 Lyapunov 稳定性定义与直接法

**定义**：平衡 $\vec x=0$，$\dot{\vec x}=f(\vec x)$。若 $\forall\epsilon>0\ \exists\delta>0$ 使 $\|\vec x(0)\|<\delta\Rightarrow\|\vec x(t)\|<\epsilon\ \forall t\ge0$，则稳定；若进一步 $\vec x(t)\to0$，则渐近稳定。

**假设**：$f$ 局部 Lipschitz，原点为平衡。

**推导**：取 Lyapunov 函数 $V(\vec x)$ 连续可微，$V(0)=0$，$V(\vec x)>0\ (\vec x\ne0)$。计算
$$\dot V=\frac{\partial V}{\partial\vec x}\cdot f(\vec x).$$
若 $\dot V\le0$，则 $V(\vec x(t))\le V(\vec x(0))$，对任意 $\epsilon$ 取 $V$ 的水平集 $\{V\le c\}\subset B_\epsilon$，可选 $\delta$ 使 $B_\delta\subset\{V\le c\}$，得稳定。若 $\dot V<0\ (\vec x\ne0)$，则 $V$ 单调下降，结合 LaSalle 不变原理得渐近稳定。

**意义**：无需解方程即可判定稳定性。

**特例**：$V=\frac12\vec x^\top P\vec x$ 对线性系统给出 Lyapunov 方程。

**应用**：控制系统、非线性力学、电力系统。

#### 知识点12.2 线性化稳定性与特征值判据

**定义**：$\dot{\vec x}=A\vec x$，$A=\partial f/\partial\vec x|_{0}$。

**假设**：$f$ 在原点是 $C^1$。

**推导**：解 $\vec x(t)=e^{At}\vec x(0)$。若 $A$ 可对角化，$e^{At}=\sum_ke^{\lambda_kt}\vec v_k\vec w_k^\top$。故
$$\|\vec x(t)\|\le C e^{(\max\Re\lambda_k)t}\|\vec x(0)\|.$$
- 所有 $\Re\lambda_k<0$：渐近稳定；
- 存在 $\Re\lambda_k>0$：不稳定；
- $\Re\lambda_k\le0$ 且重根在虚轴：需 Jordan 块判断（可能多项式增长）。

**意义**：线性化谱决定非线性系统在双曲平衡点的稳定性（Hartman–Grobman）。

**特例**：纯虚特征值需非线性项判断（中心流形）。

**应用**：平衡点分岔、控制器极点配置。

#### 知识点12.3 Routh–Hurwitz 判据

**定义**：特征多项式 $p(s)=a_ns^n+\cdots+a_1s+a_0$。

**假设**：$a_n>0$，实系数。

**推导**：构造 Routh 表。首两行为
$$\begin{array}{c|cccc} s^n & a_n & a_{n-2} & a_{n-4} & \cdots\\ s^{n-1} & a_{n-1} & a_{n-3} & a_{n-5} & \cdots \end{array}$$
后续行由
$$b_1=\frac{a_{n-1}a_{n-2}-a_na_{n-3}}{a_{n-1}},\quad b_2=\frac{a_{n-1}a_{n-4}-a_na_{n-5}}{a_{n-1}},\ \dots$$
生成。判据：所有根实部为负的充要条件是 Routh 表首列所有元素同号（取正）。

**意义**：无需解根即可判稳，适合参数设计。

**特例**：二阶系统 $s^2+2\zeta\omega_ns+\omega_n^2$ 稳定要求 $\zeta>0,\omega_n>0$。

**应用**：PID 参数整定、飞机纵向稳定。

#### 知识点12.4 Nyquist 稳定判据

**定义**：开环传递函数 $L(s)$，闭环 $T=L/(1+L)$。

**假设**：$L(s)$ 解析、开环右半平面极点数 $P$ 已知。

**推导**：闭环特征方程 $1+L(s)=0$，其右半平面零点数 $Z$。由辐角原理
$$\oint\frac{(1+L)'}{1+L}\mathrm ds=2\pi i(Z-P).$$
令 $s$ 沿 Nyquist 围线（虚轴+右半大圆弧），则 $1+L$ 曲线绕原点圈数 $N=Z-P$。稳定要求 $Z=0$，故
$$N=-P.$$
即 $L(s)$ 曲线绕 $-1$ 点的圈数满足稳定条件。

**意义**：用开环频响判断闭环稳定，可含延迟。

**特例**：$P=0$ 时曲线不包围 $-1$ 点即稳定。

**应用**：频域控制器设计、鲁棒裕度（增益/相位裕度）。

#### 知识点12.5 结构稳定性与分岔初步

**定义**：参数 $\mu$，平衡 $\vec x_e(\mu)$ 随 $\mu$ 变化。

**假设**：平衡方程 $f(\vec x;\mu)=0$ 光滑。

**推导**：分岔条件为 Jacobian $D_xf(\vec x_e;\mu)$ 奇异。简单分岔（鞍结、跨临界、叉式）可由中心流形约化到一维。叉式分岔标准形 $\dot x=\mu x-x^3$ 有平衡 $x=0$ 与 $x=\pm\sqrt\mu$（$\mu>0$），线性化特征值 $\mu-3x^2$。

**意义**：稳定性的定性变化与参数相关，是失稳与突变的前兆。

**特例**：$\mu<0$ 仅有 $x=0$ 稳定；$\mu>0$ 时 $x=0$ 不稳定。

**应用**：压杆屈曲、颤振、化学反应振荡。

#### 知识点12.6 Lyapunov 方程与线性矩阵不等式

**定义**：连续 Lyapunov 方程 $A^\top P+PA=-Q$；线性矩阵不等式（LMI）形式 $A^\top P+PA\prec0,\ P\succ0$。

**假设**：$A$ 为定常矩阵、$P,Q$ 对称。

**推导**：取 $V=\vec x^\top P\vec x$，$P=P^\top\succ0$。沿 $\dot{\vec x}=A\vec x$
$$\dot V=\dot{\vec x}^\top P\vec x+\vec x^\top P\dot{\vec x}=\vec x^\top(A^\top P+PA)\vec x.$$
故 $A^\top P+PA=-Q\prec0$ 时 $\dot V<0$，渐近稳定。反之若 $A$ 稳定，取
$$P=\int_0^\infty e^{A^\top t}Qe^{At}\,\mathrm dt,$$
被积项有界且指数衰减。利用
$$\frac{\mathrm d}{\mathrm dt}\left(e^{A^\top t}Qe^{At}\right)=A^\top e^{A^\top t}Qe^{At}+e^{A^\top t}Qe^{At}A,$$
从 $0$ 到 $\infty$ 积分得
$$0-Q=A^\top P+PA.$$

**意义**：稳定性判据转化为凸优化（LMI）问题，便于多目标鲁棒设计。

**特例**：$Q=I$ 得 $P$ 唯一；离散系统为 $A^\top PA-P=-Q$。

**应用**：鲁棒控制、观测器设计、大系统稳定性。

### 第13章 控制理论基础（状态空间与频域）

#### 知识点13.1 线性时不变系统与状态空间

**定义**：状态方程
$$\dot{\vec x}=A\vec x+B\vec u,\qquad \vec y=C\vec x+D\vec u.$$

**假设**：线性时不变、有限维、连续时间。

**推导**：解的存在唯一性由 $A$ 局部 Lipschitz 保证。状态转移矩阵
$$\Phi(t)=e^{At}=\sum_{k=0}^\infty\frac{A^kt^k}{k!},$$
满足 $\dot\Phi=A\Phi$，$\Phi(0)=I$。零输入响应
$$\vec x(t)=e^{A(t-t_0)}\vec x(t_0),$$
全响应
$$\vec x(t)=e^{A(t-t_0)}\vec x(t_0)+\int_{t_0}^t e^{A(t-\tau)}B\vec u(\tau)\,\mathrm d\tau.$$

**意义**：状态空间统一描述多输入多输出系统与初值/输入作用。

**特例**：$A$ 稳定时 $e^{At}\to0$。

**应用**：现代控制、卡尔曼滤波、观测器。

#### 知识点13.2 Laplace 变换与传递函数

**定义**：传递函数 $G(s)=C(sI-A)^{-1}B+D$。

**假设**：零初值、单输入单输出或矩阵情形。

**推导**：对 $\dot{\vec x}=A\vec x+B\vec u$ 作 Laplace 变换：
$$sX(s)-\vec x(0)=AX(s)+BU(s).$$
零初值下 $X=(sI-A)^{-1}BU$，故
$$Y=CX+DU=\bigl[C(sI-A)^{-1}B+D\bigr]U=G(s)U.$$
由伴随矩阵 $(sI-A)^{-1}=\mathrm{adj}(sI-A)/\det(sI-A)$，极点为 $A$ 的特征值。

**意义**：输入–输出关系由有理函数描述，与状态实现相互转换。

**特例**：一阶系统 $G=K/(\tau s+1)$，二阶 $G=\omega_n^2/(s^2+2\zeta\omega_ns+\omega_n^2)$。

**应用**：频域设计、滤波器、建模辨识。

#### 知识点13.3 可控性与可观性判据

**定义**：可控性 Gram 矩阵 $W_c(t)=\int_0^te^{A\tau}BB^\top e^{A^\top\tau}\mathrm d\tau$；可控性矩阵 $\mathcal C=[B\ AB\ \cdots\ A^{n-1}B]$。

**假设**：线性时不变、有限维。

**推导**：由状态解
$$\vec x(t)=e^{At}\vec x_0+\int_0^te^{A(t-\tau)}B\vec u(\tau)\mathrm d\tau.$$
若可达目标 $\vec x_1$，需存在 $\vec u$。利用 $e^{A(t-\tau)}=\sum_k\alpha_k(\tau)A^k$，积分项落在 $\mathrm{span}\{B,AB,\dots,A^{n-1}B\}$。故可控当且仅当
$$\mathrm{rank}\,\mathcal C=n.$$
同理可观性要求
$$\mathrm{rank}\,[C^\top,\ A^\top C^\top,\ \dots,\ (A^\top)^{n-1}C^\top]=n.$$

**意义**：可达性与状态重构能力的代数刻画。

**特例**：单输入系统可控 iff $\det\mathcal C\ne0$。

**应用**：控制器/观测器存在性、模型降阶。

#### 知识点13.4 线性反馈与极点配置

**定义**：反馈律 $\vec u=-K\vec x+\vec r$，闭环 $\dot{\vec x}=(A-BK)\vec x+B\vec r$。

**假设**：$(A,B)$ 可控。

**推导**：可控标准形下，特征多项式系数与控制增益一一对应。对一般系统，可控性保证存在非奇异变换 $T$ 使
$$T^{-1}AT=\begin{pmatrix}0&1&\cdots&0\\ \vdots&&\ddots&\vdots\\ -a_0&-a_1&\cdots&-a_{n-1}\end{pmatrix},\qquad T^{-1}B=\begin{pmatrix}0\\\vdots\\1\end{pmatrix}.$$
闭环特征多项式 $\det(sI-A+BK)$ 的系数随 $K$ 线性可配，故可任意配置极点。Ackermann 公式
$$K=\begin{pmatrix}0&\cdots&0&1\end{pmatrix}\mathcal C^{-1}\phi(A),$$
$\phi$ 为期望特征多项式。

**意义**：可控系统可通过状态反馈任意配置动态。

**特例**：二阶系统可独立配置 $\omega_n,\zeta$。

**应用**：飞行控制、电机调速、振动抑制。

#### 知识点13.5 观测器与分离原理

**定义**：Luenberger 观测器
$$\dot{\hat{\vec x}}=A\hat{\vec x}+B\vec u+L(\vec y-C\hat{\vec x}).$$

**假设**：$(A,C)$ 可观。

**推导**：误差 $\vec e=\vec x-\hat{\vec x}$，
$$\dot{\vec e}=(A-LC)\vec e.$$
可观性保证 $A-LC$ 极点可任意配置，故 $L$ 使 $\vec e\to0$。整个系统
$$\begin{pmatrix}\dot{\vec x}\\\dot{\vec e}\end{pmatrix}=\begin{pmatrix}A-BK&BK\\0&A-LC\end{pmatrix}\begin{pmatrix}\vec x\\\vec e\end{pmatrix},$$
特征值集合为控制器极点与观测器极点之并（分离原理）。

**意义**：状态反馈与观测器可独立设计。

**特例**：全维观测器；降维观测器利用部分可测输出。

**应用**：Kalman 滤波器、无传感器控制、过程监测。

#### 知识点13.6 Kalman 滤波与最优状态估计

**定义**：随机系统
$$\dot{\vec x}=A\vec x+B\vec u+\vec w,\qquad \vec y=C\vec x+\vec v,$$
其中 $\vec w\sim N(0,Q),\ \vec v\sim N(0,R)$ 为白噪声。

**假设**：噪声零均值、不相关，初始状态 $\vec x_0\sim N(\hat{\vec x}_0,P_0)$。

**推导**：估计动态为 $\dot{\hat{\vec x}}=A\hat{\vec x}+B\vec u+L(\vec y-C\hat{\vec x})$。估计误差 $\vec e=\vec x-\hat{\vec x}$ 满足
$$\dot{\vec e}=(A-LC)\vec e+\vec w-L\vec v,$$
误差协方差 $P=\mathbb E[\vec e\vec e^\top]$ 满足
$$\dot P=(A-LC)P+P(A-LC)^\top+Q+LRL^\top.$$
展开
$$\dot P=AP+PA^\top+Q-LCP-PC^\top L^\top+LRL^\top.$$
对 $L$ 最小化 $\mathrm{tr}\,\dot P$：$\partial\,\mathrm{tr}\,\dot P/\partial L=0\Rightarrow -2PC^\top+2LR=0$，故
$$L=PC^\top R^{-1}.$$
代回得 Riccati 方程
$$\dot P=AP+PA^\top-PC^\top R^{-1}CP+Q.$$

**意义**：在噪声统计已知时为最小方差估计；与 LQR 构成对偶性。

**特例**：定常稳态时 $\dot P=0$，解代数 Riccati 方程。

**应用**：导航、目标跟踪、传感器融合。

### 第14章 经典频域控制与 PID

#### 知识点14.1 PID 控制器与结构

**定义**：
$$u(t)=K_pe(t)+K_i\int_0^te(\tau)\mathrm d\tau+K_d\frac{\mathrm de(t)}{\mathrm dt},$$
$$C(s)=K_p+\frac{K_i}{s}+K_ds.$$

**假设**：线性、误差 $e=r-y$；执行器动态被忽略或足够快。

**推导**：分析各环节作用。比例环节增益 $K_p$ 提供误差直接修正；积分环节使稳态误差趋于零（对阶跃输入系统型别提升）；微分环节提供超前相位，增加阻尼。

闭环一阶惯性被控对象 $P(s)=K/(\tau s+1)$ 采用 PI：
$$T(s)=\frac{(K_p+K_i/s)P}{1+(K_p+K_i/s)P}=\frac{K(K_ps+K_i)}{s(\tau s+1)+K(K_ps+K_i)},$$
稳态误差 $e_{ss}=\lim_{s\to0}\frac{sR(s)}{1+L(s)}=0$（I 作用）。

**意义**：三项分别处理稳态精度、响应速度、阻尼。

**特例**：PD 用于二阶系统增加 $\zeta$；PID 带积分饱和需抗饱和处理。

**应用**：过程控制、姿态控制、巡航。

#### 知识点14.2 Bode 图与稳定裕度

**定义**：$L(j\omega)$ 的幅频 $|L|$ 与相频 $\angle L$；增益穿越频率 $\omega_c$（$|L|=1$），相位穿越频率 $\omega_g$（$\angle L=-180^\circ$）。

**假设**：最小相位或已知非最小相位、开环稳定。

**推导**：闭环稳定裕度
$$\mathrm{GM}=\frac{1}{|L(j\omega_g)|}\ (\text{或 dB}),\qquad \mathrm{PM}=180^\circ+\angle L(j\omega_c).$$
由 Nyquist 判据，$P=0$ 时稳定要求 $L(j\omega)$ 不包围 $-1$。GM 衡量幅度还能增大多少倍才包围；PM 衡量相位还能滞后多少度。典型设计 $\mathrm{PM}\approx45^\circ\sim60^\circ$。

**意义**：频域裕度给出鲁棒性量化。

**特例**：一阶系统相位最多 $-90^\circ$，PM 总为正。

**应用**：环路整形、抗扰设计。

#### 知识点14.3 根轨迹法

**定义**：闭环极点随增益 $K$ 变化的轨迹，方程
$$1+KL(s)=0,\qquad L(s)=\frac{\prod(s-z_i)}{\prod(s-p_j)}.$$

**假设**：开环零极点已知。

**推导**：根轨迹满足
$$|KL(s)|=1,\qquad \angle L(s)=(2k+1)\pi.$$
由此得绘制规则：起点为开环极点（$K=0$），终点为零点（$K\to\infty$）；渐近线交角 $(2k+1)\pi/(n-m)$，交点
$$\sigma_a=\frac{\sum p_j-\sum z_i}{n-m};$$
实轴根轨迹在奇数个零极点左侧；分离点满足 $\mathrm dL/\mathrm ds=0$。

**意义**：直观显示增益对动态的影响。

**特例**：二阶系统根轨迹为从极点到零点或渐近线。

**应用**：增益选取、补偿器设计。

#### 知识点14.4 频域性能指标与带宽

**定义**：共振峰值 $M_r$、带宽 $\omega_b$、闭环带宽。

**假设**：二阶主导。

**推导**：二阶闭环
$$T(s)=\frac{\omega_n^2}{s^2+2\zeta\omega_ns+\omega_n^2},$$
$$|T(j\omega)|=\frac{\omega_n^2}{\sqrt{(\omega_n^2-\omega^2)^2+(2\zeta\omega_n\omega)^2}}.$$
峰值条件 $\mathrm d|T|^2/\mathrm d\omega=0$ 给出
$$\omega_r=\omega_n\sqrt{1-2\zeta^2},\qquad M_r=\frac{1}{2\zeta\sqrt{1-\zeta^2}}\quad(0<\zeta<1/\sqrt2).$$
带宽（$-3\,\mathrm{dB}$）
$$\omega_b=\omega_n\sqrt{(1-2\zeta^2)+\sqrt{(1-2\zeta^2)^2+1}}.$$

**意义**：$\zeta$ 决定超调与带宽，权衡快速性与平稳性。

**特例**：$\zeta=0.707$ 时 $\omega_b=\omega_n$，最平坦。

**应用**：伺服带宽、传感器选择、振动隔离。

#### 知识点14.5 稳态误差与系统型别

**定义**：开环传递函数 $L(s)$ 在原点极点个数 $\nu$ 称为系统型别；静态误差系数
$$K_p=\lim_{s\to0}L(s),\quad K_v=\lim_{s\to0}sL(s),\quad K_a=\lim_{s\to0}s^2L(s).$$

**假设**：闭环稳定、参考输入为阶跃/斜坡/抛物线。

**推导**：误差 $E(s)=\dfrac{R(s)}{1+L(s)}$。由终值定理
$$e_{ss}=\lim_{s\to0}sE(s)=\lim_{s\to0}\frac{sR(s)}{1+L(s)}.$$
- 阶跃 $R=1/s$，$\nu=0$：$e_{ss}=\dfrac{1}{1+K_p}$；$\nu\ge1$：$e_{ss}=0$。
- 斜坡 $R=1/s^2$，$\nu=1$：$e_{ss}=1/K_v$；$\nu\ge2$：$e_{ss}=0$。

**意义**：型别每提高一阶，可无静差跟踪高一阶输入，但稳定性裕度下降。

**特例**：一阶系统 $\nu=0$ 跟踪阶跃有静差。

**应用**：伺服系统精度指标、控制器选型。

#### 知识点14.6 超前–滞后校正器设计

**定义**：超前校正
$$C(s)=K\frac{Ts+1}{\alpha Ts+1},\qquad \alpha<1.$$

**假设**：最小相位对象、需提升相位裕度。

**推导**：令 $s=j\omega$，
$$C(j\omega)=K\frac{1+j\omega T}{1+j\alpha\omega T},$$
相角 $\phi(\omega)=\arctan(\omega T)-\arctan(\alpha\omega T)$。求极值：
$$\frac{\mathrm d\phi}{\mathrm d\omega}=\frac{T}{1+\omega^2T^2}-\frac{\alpha T}{1+\alpha^2\omega^2T^2}=0\ \Rightarrow\ \omega_m=\frac{1}{T\sqrt\alpha}.$$
代入得最大超前角
$$\sin\phi_m=\frac{1-\alpha}{1+\alpha}.$$
滞后校正 $\alpha>1$ 则产生相位滞后、提升低频增益，用于降低稳态误差而不显著牺牲带宽。

**意义**：利用零极点配置在指定频段整形环路。

**特例**：$\alpha\to0$ 时 $\phi_m\to90^\circ$（理想微分）。

**应用**：直流电机伺服、过程回路整定。

### 第15章 最优控制与变分方法

#### 知识点15.1 变分法与 Euler–Lagrange 方程

**定义**：泛函
$$J[y]=\int_{x_1}^{x_2}F(x,y,y')\mathrm dx.$$

**假设**：$F$ 二次可微、$y$ 端点固定。

**推导**：
$$\delta J=\int_{x_1}^{x_2}\left(\frac{\partial F}{\partial y}\delta y+\frac{\partial F}{\partial y'}\delta y'\right)\mathrm dx.$$
分部积分
$$\int_{x_1}^{x_2}\frac{\partial F}{\partial y'}\delta y'\,\mathrm dx=\left[\frac{\partial F}{\partial y'}\delta y\right]_{x_1}^{x_2}-\int_{x_1}^{x_2}\frac{\mathrm d}{\mathrm dx}\frac{\partial F}{\partial y'}\delta y\,\mathrm dx.$$
端点固定，故
$$\delta J=\int_{x_1}^{x_2}\left(\frac{\partial F}{\partial y}-\frac{\mathrm d}{\mathrm dx}\frac{\partial F}{\partial y'}\right)\delta y\,\mathrm dx=0,$$
得 Euler–Lagrange 方程
$$\frac{\partial F}{\partial y}-\frac{\mathrm d}{\mathrm dx}\frac{\partial F}{\partial y'}=0.$$

**意义**：泛函极值的必要条件；力学 Hamilton 原理的数学形式。

**特例**：Beltrami 恒等式 $\frac{\partial F}{\partial y'}y'-F=$ 常数（$F$ 不显含 $x$）。

**应用**：最短路径、悬链线、最优控制。

#### 知识点15.2 最优控制问题与 Pontryagin 极大值原理

**定义**：系统 $\dot{\vec x}=f(\vec x,\vec u,t)$，性能指标
$$J=\phi(\vec x(t_f),t_f)+\int_{t_0}^{t_f}L(\vec x,\vec u,t)\mathrm dt.$$

**假设**：$f,L,\phi$ 可微，控制 $\vec u\in\mathcal U$。

**推导**：引入协态 $\boldsymbol\lambda$，构造哈密顿函数
$$H(\vec x,\vec u,\boldsymbol\lambda,t)=L(\vec x,\vec u,t)+\boldsymbol\lambda^\top f(\vec x,\vec u,t).$$
由变分法可得协态方程
$$\dot{\boldsymbol\lambda}=-\frac{\partial H}{\partial\vec x},$$
状态方程
$$\dot{\vec x}=\frac{\partial H}{\partial\boldsymbol\lambda},$$
最优控制使
$$H(\vec x^*,\vec u^*,\boldsymbol\lambda,t)\le H(\vec x^*,\vec u,\boldsymbol\lambda,t)$$
（极大值原理，最小化问题时为最小）。横截条件
$$\boldsymbol\lambda(t_f)=\frac{\partial\phi}{\partial\vec x(t_f)},\qquad H(t_f)=-\frac{\partial\phi}{\partial t_f}.$$

**意义**：给出最优控制的必要条件，处理控制约束。

**特例**：无约束时 $\partial H/\partial\vec u=0$； bang–bang 控制。

**应用**：航天器轨道转移、最短时间控制、能源最优。

#### 知识点15.3 线性二次型调节器（LQR）

**定义**：线性系统 $\dot{\vec x}=A\vec x+B\vec u$，性能
$$J=\int_0^\infty(\vec x^\top Q\vec x+\vec u^\top R\vec u)\mathrm dt,\qquad Q\ge0,\ R>0.$$

**假设**：$(A,B)$ 可控、$(A,\sqrt Q)$ 可检测。

**推导**：设值函数 $V(\vec x)=\vec x^\top P\vec x$，$P$ 对称正定。Hamilton–Jacobi–Bellman 方程：
$$0=\min_{\vec u}\left[\vec x^\top Q\vec x+\vec u^\top R\vec u+\nabla V^\top(A\vec x+B\vec u)\right].$$
对 $\vec u$ 求极小：
$$\frac{\partial}{\partial\vec u}\left[\vec u^\top R\vec u+2\vec x^\top PB\vec u\right]=2R\vec u+2B^\top P\vec x=0,$$
$$\vec u^*=-R^{-1}B^\top P\vec x=-K\vec x.$$
代回得
$$\vec x^\top\left[A^\top P+PA-PBR^{-1}B^\top P+Q\right]\vec x=0,$$
对任意 $\vec x$ 成立，故代数 Riccati 方程
$$A^\top P+PA-PBR^{-1}B^\top P+Q=0.$$
闭环 $A-BK$ 稳定。

**意义**：最优状态反馈设计，兼顾状态误差与控制代价。

**特例**：$Q\to0$ 退化为最小能量；$R\to0$ 趋近快速控制。

**应用**：飞行控制、倒立摆、电力系统稳定。

#### 知识点15.4 动态规划与 Hamilton–Jacobi–Bellman 方程

**定义**：值函数
$$V(\vec x,t)=\min_{\vec u}\left[\phi(\vec x(t_f))+\int_t^{t_f}L\,\mathrm d\tau\right].$$

**假设**：最优性原理成立。

**推导**：由 Bellman 最优性原理，
$$V(\vec x,t)=\min_{\vec u}\left[L\Delta t+V(\vec x+f\Delta t,t+\Delta t)\right]+o(\Delta t).$$
Taylor 展开 $V(\vec x+f\Delta t,t+\Delta t)=V+\nabla_xV\cdot f\Delta t+\partial_tV\Delta t$，令 $\Delta t\to0$：
$$-\frac{\partial V}{\partial t}=\min_{\vec u}\left[L+\nabla_xV\cdot f\right].$$
即 HJB 方程，边界条件 $V(\vec x(t_f),t_f)=\phi(\vec x(t_f))$。

**意义**：最优控制的充分条件，全局最优。

**特例**：LQR 中 $V=\vec x^\top P\vec x$ 退化为 Riccati 方程。

**应用**：随机控制、自适应动态规划、强化学习。

#### 知识点15.5 两点边值问题与打靶法

**定义**：由极大值原理得
$$\dot{\vec x}=\frac{\partial H}{\partial\boldsymbol\lambda},\qquad \dot{\boldsymbol\lambda}=-\frac{\partial H}{\partial\vec x},$$
初值 $\vec x(t_0)$ 已知，横截条件 $\boldsymbol\lambda(t_f)$ 已知。

**假设**：边界条件充分、解存在。

**推导**：将状态与协态组合为 $2n$ 维系统 $\dot{\vec z}=g(\vec z)$，其中未知初值 $\boldsymbol\lambda(t_0)=\boldsymbol\lambda_0$ 通过求解得到 $\boldsymbol\lambda(t_f;\boldsymbol\lambda_0)=\partial\phi/\partial\vec x$。定义残差 $\vec r(\boldsymbol\lambda_0)=0$，用 Newton 法迭代
$$\boldsymbol\lambda_0^{k+1}=\boldsymbol\lambda_0^k-\left[\frac{\partial\vec r}{\partial\boldsymbol\lambda_0}\right]^{-1}\vec r(\boldsymbol\lambda_0^k),$$
其中 Jacobian 由变分方程积分得到。

**意义**：将最优控制化为数值边值问题。

**特例**：横截条件缺失时需附加约束。

**应用**：轨迹优化、航天任务设计、机器人规划。

#### 知识点15.6 模型预测控制初步

**定义**：在采样时刻 $k$ 求解有限时域优化
$$\min_{\{u_k,\dots,u_{k+N-1}\}}\sum_{i=0}^{N-1}\left(x_{k+i}^\top Qx_{k+i}+u_{k+i}^\top Ru_{k+i}\right)+x_{k+N}^\top P_fx_{k+N},$$
约束 $x_{k+i+1}=Ax_{k+i}+Bu_{k+i}$，$x_{k+i}\in\mathcal X,\ u_{k+i}\in\mathcal U$，仅实施 $u_k$，下一时刻滚动重解。

**假设**：模型精确或足够好、优化问题可行、计算实时。

**推导**：将预测状态展开为控制序列的线性函数
$$\vec X=\mathcal Ax_k+\mathcal B\vec U,$$
其中 $\vec U=(u_k,\dots,u_{k+N-1})$。性能指标化为二次型
$$J=\vec U^\top H\vec U+2x_k^\top F\vec U+x_k^\top Gx_k,\qquad H\succ0.$$
无约束时 $\vec U^*=-H^{-1}F^\top x_k$，即 $u_k=-Kx_k$。含约束时求解二次规划。终端代价 $P_f$ 满足离散 Riccati 方程可保证递推可行性。

**意义**：显式处理约束与多变量耦合，闭环为滚动优化。

**特例**：$N=1$、无约束时近似 LQR。

**应用**：化工过程、自动驾驶、航天器交会。

### 第16章 拉格朗日系统的控制与能量方法

#### 知识点16.1 拉格朗日系统方程与无源性

**定义**：机械系统
$$M(q)\ddot q+C(q,\dot q)\dot q+g(q)=\tau,$$
其中 $M$ 正定对称，$\dot M-2C$ 反对称。

**假设**：完整约束、无耗散或含阻尼。

**推导**：由拉格朗日方程，动能 $T=\frac12\dot q^\top M(q)\dot q$，得
$$M\ddot q+\dot M\dot q-\frac12\frac{\partial}{\partial q}(\dot q^\top M\dot q)=\tau,$$
定义 $C\dot q=\dot M\dot q-\frac12\frac{\partial}{\partial q}(\dot q^\top M\dot q)$。可验证 $\dot M-2C$ 反对称。系统能量
$$E=T+V,\qquad \dot E=\dot q^\top\tau-\dot q^\top D\dot q.$$
无外力且无阻尼时能量守恒；加入阻尼 $D>0$ 时 $\dot E=-\dot q^\top D\dot q\le0$，故 $(\dot q,q)$ 稳定（无源性）。

**意义**：机械系统固有能量耗散结构，是能量整形控制的基础。

**特例**：把 $M\ddot q+g=\tau$ 视为端口哈密顿系统。

**应用**：机器人柔顺控制、双臂协作、稳定抓取。

#### 知识点16.2 基于能量的 Lyapunov 控制

**定义**：总能量 $V(q,\dot q)=\frac12\dot q^\top M(q)\dot q+U(q)$。

**假设**：$U$ 有下界、$M$ 正定。

**推导**：取 Lyapunov 函数 $V$，计算
$$\dot V=\dot q^\top M\ddot q+\frac12\dot q^\top\dot M\dot q+\dot q^\top\frac{\partial U}{\partial q}.$$
由运动方程 $M\ddot q+C\dot q+\partial_qU=\tau$ 及反对称性，
$$\dot V=\dot q^\top\tau-\dot q^\top D\dot q.$$
选取控制
$$\tau=-K_d\dot q+g(q)+u,$$
则 $\dot V=-\dot q^\top(K_d+D)\dot q\le0$。由 LaSalle 原理，$\dot q\to0$，进而 $g(q)\to0$，收敛到期望位形。

**意义**：用系统能量构造 Lyapunov 函数实现稳定控制。

**特例**：PD+重力补偿控制使机械臂稳定到目标。

**应用**：机器人位置控制、多体系统镇定。

#### 知识点16.3 端口哈密顿系统与互联

**定义**：
$$\dot{\vec x}=(J-R)\frac{\partial H}{\partial\vec x}+B\vec u,$$
$J=-J^\top$（互联），$R=R^\top\ge0$（耗散）。

**假设**：无源端口结构。

**推导**：能量变化
$$\dot H=\frac{\partial H^\top}{\partial x}\dot x=\frac{\partial H^\top}{\partial x}(J-R)\frac{\partial H}{\partial x}+\frac{\partial H^\top}{\partial x}B u.$$
反对称项 $\nabla H^\top J\nabla H=0$，故
$$\dot H=-\nabla H^\top R\nabla H+\vec y^\top u,\qquad \vec y=B^\top\nabla H.$$
即端口无源：$\dot H\le\vec y^\top u$。互联两个无源端口仍无源（无源性定理）。

**意义**：统一描述机械、电气、液压等能量域；便于能量整形与互联控制。

**特例**：LC 电路、刚体、质量–弹簧系统均属此框架。

**应用**：能量整形、协同机器人、电力电子。

#### 知识点16.4 反馈线性化与计算力矩

**定义**：机械臂计算力矩控制
$$\tau=M(q)(\ddot q_d-K_pe-K_de)+C\dot q+g.$$

**假设**：$M$ 可逆、模型精确。

**推导**：代入 $M\ddot q+C\dot q+g=\tau$，消去非线性项：
$$M\ddot q=M(\ddot q_d-K_pe-K_de)\ \Rightarrow\ \ddot q=\ddot q_d-K_pe-K_de.$$
令 $e=q-q_d$，则
$$\ddot e+K_d\dot e+K_pe=0.$$
取 $K_p,K_d$ 正定可使 $e\to0$ 指数收敛。

**意义**：利用模型精确消去非线性、实现线性误差动态。

**特例**：$q_d$ 常数时退化为 PD+重力补偿。

**应用**：工业机器人轨迹跟踪、无人机姿态控制。

#### 知识点16.5 自适应控制与参数估计

**定义**：未知参数 $\boldsymbol\theta$，估计 $\hat{\boldsymbol\theta}$，自适应律 $\dot{\hat{\boldsymbol\theta}}$。

**假设**：参数定常或慢变、线性参数化。

**推导**：以 $M\ddot q+C\dot q+g=\tau$ 的线性参数化 $Y(q,\dot q,\ddot q)\boldsymbol\theta=\tau$ 为例。取控制律 $\tau=Y\hat{\boldsymbol\theta}-K_ds$，$s=\dot e+\Lambda e$。误差方程
$$M\dot s+Cs+K_ds=Y\tilde{\boldsymbol\theta},$$
其中 $\tilde{\boldsymbol\theta}=\boldsymbol\theta-\hat{\boldsymbol\theta}$。取 Lyapunov
$$V=\frac12s^\top Ms+\frac12\tilde{\boldsymbol\theta}^\top\Gamma^{-1}\tilde{\boldsymbol\theta}.$$
求导并选自适应律
$$\dot{\hat{\boldsymbol\theta}}=\Gamma Y^\top s,$$
得
$$\dot V=-s^\top K_ds\le0.$$
由 Barbalat 引理 $s\to0$，参数在持续激励下收敛。

**意义**：在线修正参数适应环境变化。

**特例**：无参数误差时退化为计算力矩。

**应用**：负载变化机器人、飞行器气动参数辨识。

#### 知识点16.6 滑模控制基础

**定义**：滑模面 $s(\vec x)=0$，控制 $\vec u=\vec u_{eq}-K\mathrm{sgn}(s)$。

**假设**：系统可写为 $\dot s=a(\vec x)+b(\vec x)u$，$b\ne0$。

**推导**：取 $V=\frac12s^2$，则
$$\dot V=s\dot s=s(a+bu).$$
选 $u=b^{-1}(-a-K\mathrm{sgn}(s))$，得
$$\dot V=-K|s|\le0,$$
且到达条件 $s\dot s\le-\eta|s|$（$\eta>0$）保证有限时间到达。到达后系统沿 $s=0$ 运动，对匹配不确定性不变。

**意义**：强鲁棒性，但需处理抖振（边界层、高阶滑模）。

**特例**：二阶系统滑模等价于 PD 约束。

**应用**：电机、飞行器、机器人抗扰。

### 第17章 典型工程系统贯通应用

#### 知识点17.1 航天器姿态动力学与控制

**定义**：刚体航天器姿态 $\boldsymbol\omega$，欧拉方程 $I\dot{\boldsymbol\omega}+\boldsymbol\omega\times(I\boldsymbol\omega)=\boldsymbol\tau$。

**假设**：刚体、力矩由推力器/动量轮提供。

**推导**：见第7章欧拉方程。加入反作用飞轮角动量 $\vec h_w$：
$$I\dot{\boldsymbol\omega}+\boldsymbol\omega\times(I\boldsymbol\omega+h_w)=-\dot h_w.$$
小角度姿态运动线性化（$\boldsymbol\omega\approx\dot{\boldsymbol\theta}$）：
$$I\ddot{\boldsymbol\theta}=\boldsymbol\tau.$$
以 $I_d=\mathrm{diag}(I_x,I_y,I_z)$，采用 PD 姿态控制
$$\boldsymbol\tau=-K_p\boldsymbol\theta-K_d\dot{\boldsymbol\theta},$$
得
$$I\ddot{\boldsymbol\theta}+K_d\dot{\boldsymbol\theta}+K_p\boldsymbol\theta=0,$$
选择 $K_p=I\omega_n^2,\ K_d=2I\zeta\omega_n$ 得指定带宽与阻尼。

**意义**：姿态动力学决定控制律设计；陀螺耦合需解耦。

**特例**：重力梯度稳定利用势能极小（俯仰轴指向地心）。

**应用**：卫星姿态捕获、对日定向、编队保持。

#### 知识点17.2 机械臂轨迹跟踪

**定义**：$n$ 关节机械臂 $M(q)\ddot q+C\dot q+g=\tau$，任务空间 $\vec x=k(q)$。

**假设**：刚性连杆、无摩擦或含摩擦、执行器理想。

**推导**：任务空间速度 $\dot{\vec x}=J(q)\dot q$，$J=\partial k/\partial q$。由计算力矩实现关节空间跟踪后，误差闭环 $\ddot e+K_d\dot e+K_pe=0$。若 $J$ 奇异则任务空间自由度退化，需阻尼最小二乘
$$\dot q=J^\top(JJ^\top+\lambda^2I)^{-1}\dot{\vec x}.$$

**意义**：关节空间与任务空间映射通过 Jacobian 描述。

**特例**：冗余机械臂有零空间运动 $N=I-J^+J$。

**应用**：焊接、装配、手术机器人。

#### 知识点17.3 气动弹性与颤振

**定义**：机翼弯扭耦合，升力由 Theodorsen 非定常理论给出；颤振为自激振动失稳。

**假设**：线性气动力、小变形、二维典型剖面。

**推导**：弯扭两自由度方程
$$\begin{pmatrix}m&S\\S&I_\alpha\end{pmatrix}\begin{pmatrix}\ddot h\\\ddot\alpha\end{pmatrix}+\begin{pmatrix}k_h&0\\0&k_\alpha\end{pmatrix}\begin{pmatrix}h\\\alpha\end{pmatrix}=\begin{pmatrix}L\\M_\alpha\end{pmatrix}.$$
气动力与 $\alpha,\dot h,\dot\alpha$ 相关，写成
$$M\ddot{\vec q}+C_a\dot{\vec q}+(K+K_a)\vec q=0.$$
稳定性由复特征值实部判断：当气动耦合使某模态阻尼由正变负时发生颤振。颤振速度 $V_F$ 由 $\Re\lambda(V)=0$ 确定。

**意义**：气动力引入的交叉刚度/阻尼导致自激振荡。

**特例**：经典弯扭颤振、失速颤振。

**应用**：机翼设计、飞行包线、桥梁风致振动。

#### 知识点17.4 车辆悬架与隔振

**定义**：1/4 车辆模型，簧上质量 $m_s$、簧下质量 $m_u$、悬架 $k_s,c_s$、轮胎 $k_t$。

**假设**：线性、路面激励 $z_r$。

**推导**：运动方程
$$m_s\ddot z_s=-k_s(z_s-z_u)-c_s(\dot z_s-\dot z_u),$$
$$m_u\ddot z_u=k_s(z_s-z_u)+c_s(\dot z_s-\dot z_u)-k_t(z_u-z_r).$$
写成状态空间
$$\dot{\vec x}=A\vec x+B\dot z_r+\text{disturbance}.$$
传递函数 $G_{z_s/z_r}(s)$ 给出乘坐舒适性；力传递率
$$T_f=\frac{|F_t|}{|k_tZ_r|}.$$
隔振要求低频隔离与共振峰抑制权衡。

**意义**：多体振动与被动/主动隔振设计。

**特例**：忽略簧下质量退化为基础隔振。

**应用**：汽车底盘、精密平台、建筑基础隔振。

#### 知识点17.5 飞行器纵向短周期近似

**定义**：纵向状态 $(\alpha,q)$，短周期近似忽略速度变化。

**假设**：定速飞行、小迎角、线性气动。

**推导**：纵向线性方程
$$\dot\alpha=\frac{Z_\alpha}{mV}\alpha+q+\frac{Z_\delta}{mV}\delta,$$
$$\dot q=M_\alpha\alpha+M_q q+M_\delta\delta.$$
短周期频率与阻尼
$$\omega_{sp}=\sqrt{-M_\alpha-\frac{Z_\alpha}{V}M_q},\qquad \zeta_{sp}=-\frac{M_q+\frac{Z_\alpha}{V}+\frac{M_\alpha\cdots}{}}{2\omega_{sp}}.$$
更完整地，特征方程
$$\lambda^2-\left(\frac{Z_\alpha}{V}+M_q\right)\lambda+\left(\frac{Z_\alpha}{V}M_q-M_\alpha\right)=0.$$
选择增稳系统增益使 $\zeta_{sp}$ 达到飞行品质要求。

**意义**：短周期模态决定操纵品质与阵风响应。

**特例**：$M_\alpha<0$ 静稳定时 $\omega_{sp}$ 增大。

**应用**：电传飞控、自动驾驶仪、阵风减缓。

#### 知识点17.6 贯通视角：从力学到控制的统一框架

**定义**：以 Lagrange–Hamilton 系统为被控对象，控制设计通过能量/无源/最优三条路径。

**假设**：系统可用 $L=T-V$ 描述、约束理想、模型已知或可辨识。

**推导**：
$$\text{力学模型}\xrightarrow{\ \text{虚功/变分}\ }L\xrightarrow{\ \text{Legendre}\ }H\xrightarrow{\ \text{无源性}\ }\dot H\le y^\top u.$$
对控制：
- PD+重力补偿：$\tau=-K_d\dot q-K_p(q-q_d)+g$，误差系统线性；
- 能量整形：$\tau=\nabla U_d-\nabla U+\text{damping}$，使闭环势能极小在期望位形；
- LQR：$J=\int(x^\top Qx+u^\top Ru)\mathrm dt$，Riccati 解给最优增益；
- 自适应/滑模：处理不确定性与扰动。
统一稳定性判据均为 Lyapunov 第二方法：
$$\dot V\le0\ \Rightarrow\ \text{稳定},\qquad \dot V<0\ \Rightarrow\ \text{渐近稳定}.$$

**意义**：动力学与控制不是割裂，而是同一能量/变分框架下的建模与干预。

**特例**：无源系统的能量整形控制可看作将闭环 Lyapunov 函数取为期望总能量。

**应用**：机器人、航天器、车辆、能源系统的一体化设计。

#### 知识点17.7 数值仿真与积分方案（力学–控制贯通）

**定义**：连续系统离散化 $\dot{\vec x}=f(\vec x,t)$，步长 $h$。

**假设**：$f$ 光滑、数值稳定性要求步长满足条件。

**推导**：显式 Euler
$$\vec x_{k+1}=\vec x_k+hf(\vec x_k,t_k),$$
局部截断误差 $O(h^2)$。经典 RK4
$$k_1=f(\vec x_k,t_k),\quad k_2=f(\vec x_k+\tfrac h2k_1,t_k+\tfrac h2),$$
$$k_3=f(\vec x_k+\tfrac h2k_2,t_k+\tfrac h2),\quad k_4=f(\vec x_k+hk_3,t_k+h),$$
$$\vec x_{k+1}=\vec x_k+\frac h6(k_1+2k_2+2k_3+k_4),$$
误差 $O(h^5)$。对哈密顿系统，辛积分器（Stormer–Verlet）
$$\vec p_{k+1/2}=\vec p_k-\frac h2\nabla V(\vec x_k),$$
$$\vec x_{k+1}=\vec x_k+hM^{-1}\vec p_{k+1/2},$$
$$\vec p_{k+1}=\vec p_{k+1/2}-\frac h2\nabla V(\vec x_{k+1})$$
保持辛结构与长期能量有界振荡。

**意义**：力学系统的几何性质要求结构保持算法。

**特例**：刚体转子长时间仿真需辛或 Lie 群积分器。

**应用**：轨道预报、多体动力学、实时控制仿真。

#### 知识点17.8 无量纲化与相似准则

**定义**：特征长度 $L$、特征速度 $U$、特征频率 $\omega_0$，无量纲变量 $\bar t=\omega_0t,\ \bar x=x/L$。

**假设**：可识别主导尺度、方程保持形式。

**推导**：以 N–S 方程为例，代入 $\vec v=U\vec v^*$，$x=Lx^*$，$t=(L/U)t^*$，$p=\rho U^2p^*$：
$$\frac{\partial\vec v^*}{\partial t^*}+(\vec v^*\cdot\nabla^*)\vec v^*=-\nabla^*p^*+\frac{\mu}{\rho UL}\nabla^{*2}\vec v^*.$$
定义 Reynolds 数 $\mathrm{Re}=\rho UL/\mu$，故
$$\frac{\mathrm D\vec v^*}{\mathrm Dt^*}=-\nabla^*p^*+\frac1{\mathrm{Re}}\nabla^{*2}\vec v^*.$$
类似地，Froude 数、Mach 数、Strouhal 数分别衡量重力、压缩性、非定常效应。

**意义**：相似准则决定模型实验与尺度律。

**特例**：$\mathrm{Re}\to\infty$ 退化为欧拉方程；$\mathrm{Re}\to0$ 为 Stokes 流。

**应用**：风洞实验、流固相似、热–流–固耦合缩放。

---
⬆ 上级领域：[[力学]]
