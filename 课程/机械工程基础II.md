---
课程名称: 机械工程基础II
学分: 3
总学时: 48
考核方式: 考试
课程属性: 必修
课程性质: 基础教育
开课单位: 机械与车辆学院
领域: 机械与工程基础
tags:
  - 课程
  - 领域/机械与工程基础
---

# 机械工程基础II


### 第1章 质点运动学

**知识点 1.1 参考系、质点模型与位置矢量**

定义：参考系是描述物体运动时选作基准的物体或坐标系；当物体尺寸远小于运动尺度时，可抽象为质点，其状态由位置矢量 $\vec r(t)$ 完全描述。

假设：空间为欧氏空间，可选惯性参考系；时间绝对、均匀流逝。

运动方程与轨迹：
$$
\vec r(t)=x(t)\hat i+y(t)\hat j+z(t)\hat k,\qquad
\begin{cases}x=x(t)\\y=y(t)\\z=z(t)\end{cases}
$$
消去时间 $t$ 得轨迹方程 $f(x,y,z)=0$。

推导速度定义：
$$
\vec v=\lim_{\Delta t\to0}\frac{\Delta\vec r}{\Delta t}=\frac{d\vec r}{dt}
=\dot x\hat i+\dot y\hat j+\dot z\hat k
$$
$$
\vec a=\frac{d\vec v}{dt}=\frac{d^2\vec r}{dt^2}
=\ddot x\hat i+\ddot y\hat j+\ddot z\hat k
$$

物理意义：$\vec v$ 沿轨迹切线，$\vec a$ 指向轨迹凹侧。极限情形：匀速直线运动 $\vec a=0$；匀加速运动 $\vec a=\text{const}$，积分得
$$
\vec v=\vec v_0+\vec a t,\qquad \vec r=\vec r_0+\vec v_0t+\tfrac12\vec a t^2
$$
应用：抛体运动、汽车加速、卫星轨道初步分析。

**知识点 1.2 自然坐标下的切向与法向加速度**

定义：以轨迹弧长 $s$ 为参数，取切向单位矢量 $\hat e_t$、主法向 $\hat e_n$、副法向 $\hat e_b$。

假设：轨迹二阶可导，曲率半径 $\rho\neq0$。

推导：由 $\vec v=v\hat e_t$ 求导，
$$
\vec a=\frac{d}{dt}(v\hat e_t)=\dot v\,\hat e_t+v\,\dot{\hat e}_t
$$
又因 $\hat e_t$ 随弧长转动，$\dfrac{d\hat e_t}{ds}=\dfrac{1}{\rho}\hat e_n$，故
$$
\dot{\hat e}_t=\frac{d\hat e_t}{ds}\frac{ds}{dt}=\frac{v}{\rho}\hat e_n
$$
代入得
$$
\vec a=a_t\hat e_t+a_n\hat e_n,\qquad
a_t=\frac{dv}{dt},\qquad a_n=\frac{v^2}{\rho}
$$

推广：三维时还有 $a_b=0$，且
$$
a=\sqrt{a_t^2+a_n^2}=\sqrt{\dot v^2+\frac{v^4}{\rho^2}}
$$

物理意义：$a_t$ 改变速率大小，$a_n$ 改变速度方向。极限：直线 $\rho\to\infty\Rightarrow a_n=0$；匀速率圆周 $a_t=0,\ a_n=v^2/R$。

应用：弯道限速、离心力估算、粒子在磁场中的回旋。

**知识点 1.3 极坐标与柱坐标下的速度、加速度**

定义：极坐标中位置由 $(r,\theta)$ 描述，$\hat e_r,\hat e_\theta$ 随 $\theta$ 转动。

推导：$\hat e_r=\cos\theta\,\hat i+\sin\theta\,\hat j$，$\hat e_\theta=-\sin\theta\,\hat i+\cos\theta\,\hat j$，则
$$
\dot{\hat e}_r=\dot\theta\,\hat e_\theta,\qquad \dot{\hat e}_\theta=-\dot\theta\,\hat e_r
$$
对 $\vec r=r\hat e_r$ 求导：
$$
\vec v=\dot r\hat e_r+r\dot{\hat e}_r
=\dot r\hat e_r+r\dot\theta\hat e_\theta
$$
再求导：
$$
\vec a=\frac{d}{dt}(\dot r\hat e_r)+\frac{d}{dt}(r\dot\theta\hat e_\theta)
$$
逐项：
$$
\frac{d}{dt}(\dot r\hat e_r)=\ddot r\hat e_r+\dot r\dot\theta\hat e_\theta
$$
$$
\frac{d}{dt}(r\dot\theta\hat e_\theta)
=(\dot r\dot\theta+r\ddot\theta)\hat e_\theta+r\dot\theta(-\dot\theta\hat e_r)
$$
合并
$$
\vec a=(\ddot r-r\dot\theta^2)\hat e_r+(r\ddot\theta+2\dot r\dot\theta)\hat e_\theta
$$

物理意义：$r\dot\theta^2$ 为离心项，$2\dot r\dot\theta$ 为科氏项。特例：匀速圆周 $r=R,\ \dot r=0,\ \ddot r=0$，则 $a_r=-R\dot\theta^2$。

柱坐标推广：$z$ 方向叠加 $\ddot z\hat k$。

应用：极坐标下中心力问题、机器人关节运动学。

**知识点 1.4 运动学基本问题的积分求解**

定义：第一类问题（已知运动求力/加速度）用求导；第二类问题（已知力求运动）用积分。

推导：由 $\ddot x=f(x,\dot x,t)$，逐次积分
$$
\dot x=\int_{0}^{t}f\,dt+C_1,\qquad x=\int_{0}^{t}\dot x\,dt+C_2
$$
初始条件 $x(0)=x_0,\ \dot x(0)=v_0$ 确定 $C_1=v_0,\ C_2=x_0$。

若力为速度的函数 $m\dot v=F(v)$，分离变量
$$
\int_{v_0}^{v}\frac{m\,dv}{F(v)}=t
$$

特例：重力场中 $m\ddot y=-mg$，积分得
$$
v=v_0-gt,\qquad y=y_0+v_0t-\tfrac12 gt^2
$$

应用：落体、刹车距离、阻尼运动。

**知识点 1.5 相对运动与平动参考系变换**

定义：设 $S'$ 相对 $S$ 平动，位矢 $\vec R_{OO'}$。

推导：由矢量三角形
$$
\vec r_A=\vec R_{OO'}+\vec r_{\,A'}
$$
对时间求导：
$$
\vec v_A=\vec v_{O'}+\vec v_{A'}\quad\Rightarrow\quad \vec v_{rel}=\vec v_A-\vec v_{O'}
$$
再求导：
$$
\vec a_A=\vec a_{O'}+\vec a_{rel}
$$

物理意义：平动参考系中加速度差别仅由牵连加速度引起，转动参考系还需补充科氏项。

应用：两车相遇、传送带、相对速度合成。

**知识点 1.6 转动参考系中的速度与加速度、科氏加速度**

定义：$S'$ 以角速度 $\vec\omega$ 绕轴转动，考察点相对 $S'$ 的位置 $\vec r'$。

推导：单位矢量转动满足
$$
\left(\frac{d\hat i'}{dt}\right)_S=\vec\omega\times\hat i',\quad\text{同理 } \hat j',\hat k'
$$
故任意相对矢量的绝对导数
$$
\left(\frac{d\vec A}{dt}\right)_S
=\left(\frac{d\vec A}{dt}\right)_{S'}+\vec\omega\times\vec A
$$
取 $\vec A=\vec r'$ 得速度合成：
$$
\vec v=\vec v_{O'}+\vec v'+ \vec\omega\times\vec r'
$$
再对速度求绝对导数并代用上式，得
$$
\vec a=\vec a_{O'}+\vec a'+2\vec\omega\times\vec v'
+\vec\alpha\times\vec r'
+\vec\omega\times(\vec\omega\times\vec r')
$$
其中 $2\vec\omega\times\vec v'$ 为科氏加速度，$\vec\omega\times(\vec\omega\times\vec r')$ 为向心加速度，$\vec\alpha\times\vec r'$ 为切向牵连。

物理意义：科氏项来自相对运动与参考系转动的耦合。特例：$\omega=0$ 退化为平动变换。

应用：地球自转对落体偏东、傅科摆、涡轮叶轮运动。

### 第2章 刚体运动学

**知识点 2.1 刚体、角速度矢量与点的速度**

定义：刚体是任意两点间距离不变的质点系。其运动可分解为随基点的平动与绕基点的转动。

推导：设基点 $O'$ 位矢 $\vec r_{O'}$，任一点 $P$ 相对 $O'$ 的矢径 $\vec\rho$ 长度不变，故
$$
\left(\frac{d\vec\rho}{dt}\right)_S=\vec\omega\times\vec\rho
$$
于是
$$
\vec v_P=\vec v_{O'}+\vec\omega\times\vec\rho
=\vec v_{O'}+\vec\omega\times(\vec r_P-\vec r_{O'})
$$

物理意义：$\vec\omega$ 是与基点选择无关的轴矢量。

特例：纯平动 $\vec\omega=0\Rightarrow\vec v_P=\vec v_{O'}$；定轴转动取轴上一点为基点，$\vec v=\vec\omega\times\vec r$。

应用：连杆机构速度分析、机构速度瞬心。

**知识点 2.2 速度瞬心与平面运动**

定义：平面运动中瞬时速度为零的点称速度瞬心 $C$。

推导：设平面内某点速度 $\vec v_C=0$，任一点速度
$$
\vec v_P=\vec v_C+\vec\omega\times(\vec r_P-\vec r_C)=\vec\omega\times(\vec r_P-\vec r_C)
$$
由 $|v_P|=\omega\,|PC|$，且方向垂直于 $PC$。为求瞬心：对两点 $A,B$ 作速度垂线，交点即瞬心。

物理意义：平面运动瞬时等价于绕瞬心转动。特例：轮子纯滚动时瞬心为接触点，$v_{center}=\omega R$。

应用：平面连杆、齿轮啮合、滚动轮速度场。

**知识点 2.3 刚体上点的加速度**

定义：加速度分解为基点加速度、切向与法向转动加速度。

推导：对 $\vec v_P=\vec v_{O'}+\vec\omega\times\vec\rho$ 求导：
$$
\vec a_P=\vec a_{O'}+\dot{\vec\omega}\times\vec\rho
+\vec\omega\times\frac{d\vec\rho}{dt}
$$
代入 $\dfrac{d\vec\rho}{dt}=\vec\omega\times\vec\rho$：
$$
\vec a_P=\vec a_{O'}+\vec\alpha\times\vec\rho
+\vec\omega\times(\vec\omega\times\vec\rho)
$$

物理意义：后两项分别为切向与向心（法向）加速度。对于平面运动且 $O'$ 取转动中心：
$$
\vec a_P=\vec a_{O'}+\vec\alpha\times\vec\rho-\omega^2\vec\rho
$$

特例：定轴匀速转动 $\vec\alpha=0$，只有向心项 $-\omega^2\vec\rho$。

应用：曲柄滑块加速度分析、凸轮机构动力学。

**知识点 2.4 欧拉角与角速度合成**

定义：用进动角 $\psi$、章动角 $\theta$、自转角 $\varphi$ 描述刚体姿态。

推导：Z-X-Z 欧拉角下，角速度由三次转动叠加：
$$
\vec\omega=\dot\psi\hat k+\dot\theta\hat e_{N}+\dot\varphi\hat e_{z'}
$$
将其写到本体坐标系 $(x',y',z')$：
$$
\begin{aligned}
\omega_{x'}&=\dot\psi\sin\theta\sin\varphi+\dot\theta\cos\varphi\\
\omega_{y'}&=\dot\psi\sin\theta\cos\varphi-\dot\theta\sin\varphi\\
\omega_{z'}&=\dot\psi\cos\theta+\dot\varphi
\end{aligned}
$$

物理意义：三个角速度仅两两组合，章动角为零时 $\psi,\varphi$ 不可分辨（欧拉角奇异）。

特例：定点转动且对称陀螺，$\theta=\text{const}$ 时 $\omega$ 在体轴上分量恒定。

应用：陀螺、飞行器姿态、卫星定向。

**知识点 2.5 角速度合成与瞬时轴**

定义：刚体同时参与多个转动时，角速度按矢量相加。

推导：设相对参考系 $S_1$ 的角速度 $\vec\omega_1$，$S_1$ 相对 $S$ 的角速度 $\vec\omega_e$，则绝对角速度
$$
\vec\omega=\vec\omega_1+\vec\omega_e
$$
由速度合成定理对任一点验证：
$$
\vec v=\vec v_e+\vec v_r
=\vec v_O+\vec\omega_e\times\vec r+\vec\omega_1\times\vec r
=\vec v_O+(\vec\omega_e+\vec\omega_1)\times\vec r
$$

特例：两转轴平行且反向、转速相等时 $\vec\omega=0$，刚体平动（如行星齿轮架）。

应用：齿轮系传动比计算、离心摆。

**知识点 2.6 点的复合运动在刚体分析中的应用**

定义：把刚体上任一点视为同时参与牵连运动与相对运动的动点。

推导：由科氏定理
$$
\vec a_P=\vec a_{O'}+\vec a_{r}+2\vec\omega\times\vec v_r
+\vec\omega\times(\vec\omega\times\vec r')
$$
其中相对加速度本身再由刚体转动给出：
$$
\vec a_{r}=\vec a_{O'}+\vec\alpha\times\vec\rho+\vec\omega\times(\vec\omega\times\vec\rho)
$$

物理意义：复合运动方法把复杂刚体运动拆成简单转动叠加，便于求约束点加速度。

特例：若动点相对体固定（$\vec v_r=0$），科氏项消失。

应用：滑块在转动导槽中运动、机构加速度图解法。

### 第3章 牛顿力学与质点动力学

**知识点 3.1 牛顿三定律与惯性系**

定义：第一定律（惯性律）定义惯性系；第二定律给出力与加速度的关系；第三定律给出作用力与反作用力。

假设：质量恒定，速度远小于光速。

牛顿第二定律：
$$
\vec F=\frac{d\vec p}{dt},\qquad \vec p=m\vec v
$$
当 $m$ 恒定时
$$
\vec F=m\vec a=m\ddot{\vec r}
$$

推导动量形式：$\vec F=\dfrac{d(m\vec v)}{dt}=\dot m\vec v+m\dot{\vec v}$，$m$ 恒定退化为 $m\vec a$。

物理意义：力是动量变化率，而非加速度本身。特例：$F=0\Rightarrow\vec p=\text{const}$。

应用：变质量火箭、碰撞问题。

**知识点 3.2 微分运动方程与求解方法**

定义：把牛顿第二定律投影到坐标系得到运动微分方程。

推导：直角坐标投影
$$
\begin{cases}
m\ddot x=F_x\\ m\ddot y=F_y\\ m\ddot z=F_z
\end{cases}
$$
自然坐标投影：
$$
\begin{cases}
m\dfrac{dv}{dt}=F_t\\[6pt]
m\dfrac{v^2}{\rho}=F_n\\[6pt]
0=F_b
\end{cases}
$$

物理意义：将矢量方程标量化，便于逐方向积分。特例：圆周运动仅法向方程决定约束反力。

应用：斜面滑动、圆锥摆、带电粒子运动。

**知识点 3.3 动量定理与动量守恒**

定义：力对时间的积分称冲量，动量变化等于冲量。

推导：对 $\vec F=\dfrac{d\vec p}{dt}$ 在 $[t_1,t_2]$ 积分：
$$
\int_{t_1}^{t_2}\vec F\,dt=\vec p_2-\vec p_1
$$
分量形式同理。

物理意义：冲量改变动量。特例：$\vec F=0\Rightarrow\vec p=\text{const}$，即动量守恒（单质点）。

应用：碰撞、冲击、爆炸。

**知识点 3.4 角动量定理（单质点）**

定义：对定点 $O$ 的角动量 $\vec L_O=\vec r\times\vec p$。

推导：求导
$$
\dot{\vec L}_O=\dot{\vec r}\times\vec p+\vec r\times\dot{\vec p}
=\vec v\times m\vec v+\vec r\times\vec F
=\vec r\times\vec F=\vec M_O
$$
其中用到 $\vec v\times\vec v=0$。

物理意义：力矩是角动量的变化率。特例：中心力 $\vec M_O=0\Rightarrow\vec L_O=\text{const}$，即开普勒第二定律。

应用：行星运动、圆周运动角动量守恒。

**知识点 3.5 保守力与势能**

定义：做功只与始末位置有关、与路径无关的力称保守力。

推导：保守力充要条件：
$$
\oint\vec F\cdot d\vec r=0\quad\Longleftrightarrow\quad\nabla\times\vec F=0\quad\Longleftrightarrow\quad\vec F=-\nabla V
$$
由势能定义
$$
V(\vec r)-V(\vec r_0)=-\int_{\vec r_0}^{\vec r}\vec F\cdot d\vec r
$$
两边微分：$dV=-\vec F\cdot d\vec r=-\nabla V\cdot d\vec r$，故 $\vec F=-\nabla V$。

特例：重力 $V=mgz$；弹簧 $V=\tfrac12kx^2$；万有引力 $V=-GMm/r$。

物理意义：势能是位置函数，零点是任意的，只有势能差有物理意义。

应用：势能曲线分析、平衡与稳定性。

**知识点 3.6 非惯性系中的惯性力**

定义：在加速度为 $\vec a_0$ 的非惯性系中，为使牛顿定律形式上成立而引入惯性力。

推导：由加速度合成
$$
\vec a=\vec a_0+\vec a'
$$
代入 $\vec F=m\vec a$：
$$
\vec F-m\vec a_0=m\vec a'=m\vec a_{rel}
$$
即
$$
\vec F_{真实}+\vec F_{惯性}=m\vec a_{rel},\qquad \vec F_{惯性}=-m\vec a_0
$$
转动参考系中还须加入离心力 $-m\vec\omega\times(\vec\omega\times\vec r')$ 与科氏力 $-2m\vec\omega\times\vec v'$。

物理意义：惯性力不是真实相互作用，而是参考系加速的效应。特例：自由落体电梯内 $\vec F_{惯性}=-m\vec g$，物体验不到重力。

应用：离心分离机、旋转流体面形、大气环流。

### 第4章 质点系动力学：动量与质心运动

**知识点 4.1 质点系、内力与外力、质心**

定义：质心位矢
$$
\vec r_C=\frac{\sum_i m_i\vec r_i}{\sum_i m_i}=\frac{1}{M}\int\vec r\,dm
$$

推导连续体情形：由微元求和取极限
$$
\vec r_C=\lim_{\Delta m_i\to0}\frac{\sum\vec r_i\Delta m_i}{\sum\Delta m_i}
=\frac{1}{M}\int_V\vec r\,\rho\,dV
$$
直角坐标分量：
$$
x_C=\frac{1}{M}\int x\,dm,\quad y_C=\frac{1}{M}\int y\,dm,\quad z_C=\frac{1}{M}\int z\,dm
$$

物理意义：质心是质量分布的一阶矩中心。特例：半圆板质心 $y_C=4R/(3\pi)$；三角形质心为三中线交点。

应用：重心、转动惯量、碰撞分析。

**知识点 4.2 质点系动量定理**

定义：质点系总动量 $\vec P=\sum m_i\vec v_i=M\vec v_C$。

推导：对每质点 $\vec F_i^{外}+\sum_j\vec F_{ij}^{内}=\dfrac{d}{dt}(m_i\vec v_i)$，求和：
$$
\sum_i\vec F_i^{外}+\sum_i\sum_j\vec F_{ij}^{内}=\frac{d\vec P}{dt}
$$
由牛顿第三定律内力成对相消：
$$
\sum_i\sum_j\vec F_{ij}^{内}=0
$$
故
$$
\frac{d\vec P}{dt}=\sum_i\vec F_i^{外}=\vec F^{外}
$$
又 $\vec P=M\vec v_C$，于是
$$
M\vec a_C=\vec F^{外}
$$

物理意义：质心运动如同全部质量集中于质心、受全部外力作用。特例：$\vec F^{外}=0\Rightarrow\vec P=\text{const}$，质心匀速直线。

应用：爆炸碎片、火箭、人船模型。

**知识点 4.3 变质量系统的运动：火箭方程**

定义：变质量系统指系统质量随时间变化，需考虑质量带走的动量。

假设：喷出质量相对主体的速度为定值 $\vec u$，无外力或外力已知。

推导：取 $t$ 时刻系统（质量 $m$，速度 $v$），$t+dt$ 时刻主体质量 $m+dm$（$dm<0$）、速度 $v+dv$，喷出质量 $-dm$、绝对速度 $v+u$。动量变化
$$
dp=[(m+dm)(v+dv)+(-dm)(v+u)]-mv
$$
展开并略去二阶小量 $dm\,dv$：
$$
dp=m\,dv+u\,dm
$$
由 $dp=F\,dt$：
$$
m\,dv+u\,dm=F\,dt
$$
无外力时
$$
m\,dv=-u\,dm\quad\Rightarrow\quad dv=-u\frac{dm}{m}
$$
积分
$$
\Delta v=u\ln\frac{m_0}{m}
$$

物理意义：速度增量与质量比的对数成正比。特例：$m\to0$ 时 $\Delta v\to\infty$（理想）。

应用：火箭推进、喷气发动机推力 $F_{推}=-u\dfrac{dm}{dt}$。

**知识点 4.4 质心参考系与二体问题**

定义：质心系中 $\vec P=0$，两体相对运动可化为单体问题。

推导：两体质量 $m_1,m_2$，相对矢径 $\vec r=\vec r_1-\vec r_2$，约化质量
$$
\mu=\frac{m_1m_2}{m_1+m_2}
$$
牛顿方程
$$
m_1\ddot{\vec r}_1=\vec F,\qquad m_2\ddot{\vec r}_2=-\vec F
$$
相减：
$$
\ddot{\vec r}_1-\ddot{\vec r}_2=\left(\frac{1}{m_1}+\frac{1}{m_2}\right)\vec F
=\frac{1}{\mu}\vec F
$$
即
$$
\mu\ddot{\vec r}=\vec F
$$

物理意义：二体相对运动等价于质量 $\mu$ 的单体在相同力场中运动。特例：$m_1\ll m_2$ 时 $\mu\approx m_1$。

应用：双星、氢原子、行星-卫星系统。

**知识点 4.5 碰撞与恢复系数**

定义：两体碰撞中法向相对速度之比称恢复系数 $e$。

推导：设碰前相对速度 $v_{rel}=v_1-v_2$，碰后 $v'_{rel}=v'_1-v'_2$，则
$$
e=-\frac{v'_{rel}}{v_{rel}},\qquad 0\le e\le1
$$
动量守恒：
$$
m_1v_1+m_2v_2=m_1v'_1+m_2v'_2
$$
联立解得
$$
v'_1=\frac{m_1-em_2}{m_1+m_2}v_1+\frac{(1+e)m_2}{m_1+m_2}v_2
$$
$$
v'_2=\frac{(1+e)m_1}{m_1+m_2}v_1+\frac{m_2-em_1}{m_1+m_2}v_2
$$

物理意义：$e=1$ 完全弹性，$e=0$ 完全非弹性。特例：等质量完全弹性正碰交换速度。

应用：台球、汽车碰撞、粒子散射。

**知识点 4.6 动能定理与柯尼希定理**

定义：质点系动能可分解为质心动能与相对质心动能。

推导：
$$
\vec v_i=\vec v_C+\vec v_i'
$$
则
$$
T=\sum\tfrac12 m_i v_i^2
=\sum\tfrac12 m_i(\vec v_C+\vec v_i')^2
$$
展开
$$
T=\tfrac12 Mv_C^2+\vec v_C\cdot\sum m_i\vec v_i'
+\sum\tfrac12 m_i v_i'^2
$$
由质心系定义 $\sum m_i\vec v_i'=0$，故
$$
T=\tfrac12 Mv_C^2+T'
$$

物理意义：系统动能等于质心平动动能加内部动能。特例：刚体平动 $T'=0$。

应用：碰撞能量损失、刚体滚动的动能分解。

### 第5章 角动量与刚体平面运动动力学

**知识点 5.1 质点系角动量定理**

定义：对定点 $O$ 的系统角动量 $\vec L_O=\sum_i\vec r_i\times m_i\vec v_i$。

推导：求导
$$
\dot{\vec L}_O=\sum_i\dot{\vec r}_i\times m_i\vec v_i+\sum_i\vec r_i\times m_i\dot{\vec v}_i
$$
第一项 $\vec v_i\times m_i\vec v_i=0$，第二项用牛顿方程：
$$
\dot{\vec L}_O=\sum_i\vec r_i\times\left(\vec F_i^{外}+\sum_j\vec F_{ij}^{内}\right)
$$
内力矩成对相消（强形式第三定律），故
$$
\dot{\vec L}_O=\sum_i\vec r_i\times\vec F_i^{外}=\vec M_O^{外}
$$

物理意义：只有外力矩改变系统角动量。特例：$\vec M_O^{外}=0\Rightarrow\vec L_O=\text{const}$。

应用：转台、花样滑冰收臂加速、直升机尾桨。

**知识点 5.2 转动惯量与平行轴定理**

定义：刚体绕轴 $z$ 的转动惯量
$$
I_z=\sum_i m_i\rho_i^2=\int\rho_m\rho^2\,dV
$$
其中 $\rho$ 为质量元到轴的距离。

推导平行轴定理：设质心轴 $z_C$ 与 $z$ 平行、相距 $d$，任一点 $\rho_i^2=|\vec r_i'+\vec d|^2$，则
$$
I_z=\sum m_i|\vec r_i'+\vec d|^2
=\sum m_i r_i'^2+2\vec d\cdot\sum m_i\vec r_i'+Md^2
$$
质心系中 $\sum m_i\vec r_i'=0$，故
$$
I_z=I_{z_C}+Md^2
$$

垂直轴定理（薄板）：
$$
I_z=I_x+I_y
$$
推导：$r_i^2=x_i^2+y_i^2$，故 $\sum m_i r_i^2=\sum m_i x_i^2+\sum m_i y_i^2$。

物理意义：远离质心轴转动惯量增大。特例：细杆绕端点 $I=mL^2/3$，绕质心 $I=mL^2/12$。

应用：飞轮储能、连杆等效惯量。

**知识点 5.3 刚体定轴转动微分方程**

定义：刚体绕固定轴 $z$ 转动，角加速度由合外力矩决定。

推导：将 $\vec L= I_z\dot\theta\hat k$ 代入角动量定理：
$$
\frac{d}{dt}(I_z\dot\theta)=M_z
$$
$I_z$ 恒定时
$$
I_z\ddot\theta=M_z
$$

物理意义：与牛顿第二定律对应，$I$ 是转动惯性度量。特例：$M_z=0\Rightarrow\dot\theta=\text{const}$。

应用：复摆
$$
I\ddot\theta=-mgd\sin\theta\approx-mgd\,\theta
\;\Rightarrow\;
T=2\pi\sqrt{\frac{I}{mgd}}
$$

**知识点 5.4 刚体平面运动微分方程**

定义：平面运动可分解为质心平动与绕质心转动。

推导：由质心运动定理与绕质心角动量定理：
$$
\begin{cases}
M\ddot x_C=F_x\\[4pt]
M\ddot y_C=F_y\\[4pt]
I_C\ddot\theta=M_C
\end{cases}
$$
其中 $M_C$ 为外力对质心的矩。

物理意义：三个方程解三个未知量 $\ddot x_C,\ddot y_C,\ddot\theta$。特例：滚动约束 $\ddot x_C=R\ddot\theta$。

应用：斜面滚轮、圆柱滚动竞赛
$$
a_C=\frac{g\sin\alpha}{1+k^2/R^2},\qquad k^2=I_C/M
$$

**知识点 5.5 欧拉动力学方程**

定义：在随体主惯性轴上写出刚体定点转动方程。

推导：取体轴为主轴，$\vec L=I_x\omega_x\hat i+I_y\omega_y\hat j+I_z\omega_z\hat k$。在惯性系求导需用
$$
\left(\frac{d\vec L}{dt}\right)_S=\left(\frac{d\vec L}{dt}\right)_{body}+\vec\omega\times\vec L
$$
得到
$$
\begin{cases}
I_x\dot\omega_x-(I_y-I_z)\omega_y\omega_z=M_x\\
I_y\dot\omega_y-(I_z-I_x)\omega_z\omega_x=M_y\\
I_z\dot\omega_z-(I_x-I_y)\omega_x\omega_y=M_z
\end{cases}
$$
无外力矩时称欧拉自由陀螺方程。

物理意义：交叉项体现角动量方向与角速度不一致（非主轴时更复杂）。特例：轴对称陀螺 $I_x=I_y$，绕 $z$ 匀速自转、$\omega_z$ 恒定。

应用：陀螺进动、网球拍定理（中间轴不稳）。

**知识点 5.6 陀螺进动与章动**

定义：快转陀螺在外力矩下自转轴绕外力矩方向缓慢转动称进动。

推导：设自转角动量 $L_s=I_z\omega_s$ 远大于横向，外力矩 $\vec M$ 垂直于 $\vec L_s$，则
$$
\vec M=\frac{d\vec L_s}{dt}=\vec\Omega_p\times\vec L_s
$$
解得进动角速度
$$
\Omega_p=\frac{M}{L_s\sin\gamma}
$$
仅重力矩时 $M=mgd\sin\gamma$，
$$
\Omega_p=\frac{mgd}{I_z\omega_s}
$$

物理意义：力矩垂直角动量时改变其方向而非大小。特例：$\omega_s\to\infty$ 时进动减缓。

应用：陀螺罗经、炮弹旋转稳定、自行车转向。

### 第6章 功、能与机械能守恒

**知识点 6.1 功与功率**

定义：力沿路径的功
$$
W=\int_{A}^{B}\vec F\cdot d\vec r
$$
功率
$$
P=\frac{dW}{dt}=\vec F\cdot\vec v
$$

推导：微元功 $dW=\vec F\cdot d\vec r=[F_x\,dx+F_y\,dy+F_z\,dz]$，沿曲线积分得总功。

物理意义：只有切向力分量做功。特例：向心力、法向约束力不做功。

应用：起重机做功、摩擦力耗散。

**知识点 6.2 动能定理（质点与质点系）**

定义：合外力做功等于动能增量。

推导（质点）：由 $\vec F=m\dfrac{d\vec v}{dt}$，两边点乘 $d\vec r=\vec v\,dt$：
$$
\vec F\cdot d\vec r=m\vec v\cdot d\vec v=\tfrac12 m\,d(v^2)
$$
积分
$$
W_{AB}=\tfrac12 mv_B^2-\tfrac12 mv_A^2
$$
质点系：对每质点求和，内力做功之和一般不为零。

物理意义：功能关系是牛顿定律的积分形式。特例：$W=0\Rightarrow T$ 不变。

应用：斜面、弹簧枪、机车牵引。

**知识点 6.3 势能与机械能守恒**

定义：系统机械能 $E=T+V$。

推导：保守力做功 $W=-[V_B-V_A]$，动能定理 $W=T_B-T_A$，故
$$
T_B-T_A=-(V_B-V_A)
\;\Rightarrow\;
T_A+V_A=T_B+V_B
$$
即 $E=\text{const}$。若存在非保守力 $W_{nc}$，则
$$
T_B+V_B-\left(T_A+V_A\right)=W_{nc}
$$

物理意义：能量在动能与势能间转化。特例：弹簧振子 $E=\tfrac12 m\dot x^2+\tfrac12 kx^2$，由 $dE/dt=0$ 得 $m\ddot x+kx=0$。

应用：摆、轨道运动、能量法求速度。

**知识点 6.4 广义力与广义坐标下的功**

定义：设位形由 $q_1,\dots,q_n$ 描述，$\vec r_i=\vec r_i(q,t)$，广义力
$$
Q_k=\sum_i\vec F_i\cdot\frac{\partial\vec r_i}{\partial q_k}
$$
推导：虚功
$$
\delta W=\sum_i\vec F_i\cdot\delta\vec r_i
=\sum_i\vec F_i\cdot\sum_k\frac{\partial\vec r_i}{\partial q_k}\delta q_k
=\sum_k Q_k\delta q_k
$$

物理意义：$Q_k$ 与广义坐标量纲匹配。特例：若力有势 $Q_k=-\partial V/\partial q_k$。

应用：拉格朗日方程右边项、机构广义力分析。

**知识点 6.5 势能曲线与平衡稳定性**

定义：一维系统 $E=T+V(x)$，平衡点满足 $V'(x_0)=0$。

推导：对平衡点小扰动，泰勒展开
$$
V(x)\approx V(x_0)+\tfrac12 V''(x_0)(x-x_0)^2
$$
运动方程
$$
m\ddot x=-V'(x)\approx-V''(x_0)(x-x_0)
$$
故 $V''(x_0)>0$ 时简谐振动，频率
$$
\omega=\sqrt{\frac{V''(x_0)}{m}}
$$
$V''(x_0)<0$ 时指数发散，不稳定。

物理意义：势能极小对应稳定平衡。特例：$V''=0$ 需更高阶判断（如 $x^3$ 处不稳定）。

应用：分子势、结构稳定性、双稳态。

**知识点 6.6 机械能守恒的判据与耗散函数**

定义：耗散力（如黏性阻尼）做负功，可用瑞利耗散函数 $D=\tfrac12 c\dot x^2$ 描述。

推导：阻尼力 $F_d=-c\dot x$，其虚功
$$
\delta W_d=-c\dot x\,\delta x=-\frac{\partial D}{\partial\dot x}\delta x
$$
故拉格朗日方程右边加 $-\dfrac{\partial D}{\partial\dot x}$。

物理意义：耗散函数为二次型，能量单调递减。

特例：无阻尼 $c=0$ 时机械能守恒。

应用：阻尼振动、减振器设计。

### 第7章 达朗贝尔原理、虚位移原理与拉格朗日方程

**知识点 7.1 达朗贝尔原理**

定义：将惯性力 $-m\vec a_i$ 视为对每质点的“外力”，则系统处于动力平衡。

推导：牛顿方程 $\vec F_i=m\vec a_i$ 移项：
$$
\vec F_i+(-m_i\vec a_i)=0
$$
对每质点成立，故系统在任意瞬时满足静力平衡形式。

物理意义：把动力学问题化为静力学问题。特例：静力学问题取 $\vec a=0$。

应用：动静法求解轴承反力、转子平衡。

**知识点 7.2 虚位移、理想约束与虚功原理**

定义：虚位移 $\delta\vec r_i$ 是满足约束的无穷小、任意、瞬时可能位移；理想约束指约束力总虚功为零：$\sum\vec N_i\cdot\delta\vec r_i=0$。

推导：达朗贝尔原理下，令 $\vec F_i$ 包含主动力 $\vec F_i^{a}$ 与约束力 $\vec N_i$：
$$
\sum_i(\vec F_i^{a}+\vec N_i-m_i\vec a_i)\cdot\delta\vec r_i=0
$$
理想约束项消失：
$$
\sum_i(\vec F_i^{a}-m_i\vec a_i)\cdot\delta\vec r_i=0
$$
静力学情形 $\ddot{\vec r}_i=0$ 得虚功原理
$$
\sum_i\vec F_i^{a}\cdot\delta\vec r_i=0
$$

物理意义：只需主动力即可判别平衡。特例：杠杆、滑轮组平衡。

应用：复杂机构平衡、约束力为零的条件。

**知识点 7.3 广义坐标与广义力表示的达朗贝尔原理**

定义：用独立坐标 $q_1,\dots,q_n$（自由度数）描述系统。

推导：把 $\delta\vec r_i=\sum_k\dfrac{\partial\vec r_i}{\partial q_k}\delta q_k$ 代入，并注意
$$
\sum_i m_i\vec a_i\cdot\delta\vec r_i
=\sum_k\left[\frac{d}{dt}\frac{\partial T}{\partial\dot q_k}-\frac{\partial T}{\partial q_k}\right]\delta q_k
$$
（此项推导见知识点7.4）。于是
$$
\sum_k\left[\frac{d}{dt}\frac{\partial T}{\partial\dot q_k}-\frac{\partial T}{\partial q_k}-Q_k\right]\delta q_k=0
$$
由 $\delta q_k$ 独立任意，得
$$
\frac{d}{dt}\frac{\partial T}{\partial\dot q_k}-\frac{\partial T}{\partial q_k}=Q_k
$$

物理意义：把约束力自动消去的动力学方程。特例：完整约束自由度等于独立坐标数。

应用：多体机构、机器人动力学。

**知识点 7.4 动能对广义速度偏导恒等式的证明**

定义：证明 $\sum_i m_i\vec a_i\cdot\dfrac{\partial\vec r_i}{\partial q_k}=\dfrac{d}{dt}\dfrac{\partial T}{\partial\dot q_k}-\dfrac{\partial T}{\partial q_k}$。

推导：先注意 $\dot{\vec r}_i=\sum_j\dfrac{\partial\vec r_i}{\partial q_j}\dot q_j+\dfrac{\partial\vec r_i}{\partial t}$，故
$$
\frac{\partial\dot{\vec r}_i}{\partial\dot q_k}=\frac{\partial\vec r_i}{\partial q_k}
$$
再计算
$$
\frac{d}{dt}\frac{\partial T}{\partial\dot q_k}
=\frac{d}{dt}\sum_i m_i\dot{\vec r}_i\cdot\frac{\partial\dot{\vec r}_i}{\partial\dot q_k}
=\frac{d}{dt}\sum_i m_i\dot{\vec r}_i\cdot\frac{\partial\vec r_i}{\partial q_k}
$$
$$
=\sum_i m_i\vec a_i\cdot\frac{\partial\vec r_i}{\partial q_k}
+\sum_i m_i\dot{\vec r}_i\cdot\frac{d}{dt}\frac{\partial\vec r_i}{\partial q_k}
$$
又因
$$
\frac{d}{dt}\frac{\partial\vec r_i}{\partial q_k}
=\sum_j\frac{\partial^2\vec r_i}{\partial q_j\partial q_k}\dot q_j+\frac{\partial^2\vec r_i}{\partial t\partial q_k}
=\frac{\partial\dot{\vec r}_i}{\partial q_k}
$$
所以第二项等于 $\sum_i m_i\dot{\vec r}_i\cdot\dfrac{\partial\dot{\vec r}_i}{\partial q_k}=\dfrac{\partial T}{\partial q_k}$。移项即证。

物理意义：该恒等式是拉格朗日方程成立的核心。

特例：$q_k$ 为直角坐标时退化为 $m\ddot x_k=F_k$。

应用：推导一般力学系统运动方程。

**知识点 7.5 拉格朗日方程（保守与含势）**

定义：保守系统 $Q_k=-\dfrac{\partial V}{\partial q_k}$，定义拉格朗日函数 $L=T-V$。

推导：从
$$
\frac{d}{dt}\frac{\partial T}{\partial\dot q_k}-\frac{\partial T}{\partial q_k}=Q_k
$$
代入 $Q_k$ 并注意 $V$ 不显含 $\dot q_k$：
$$
\frac{d}{dt}\frac{\partial(T-V)}{\partial\dot q_k}-\frac{\partial(T-V)}{\partial q_k}=0
$$
即
$$
\boxed{\frac{d}{dt}\frac{\partial L}{\partial\dot q_k}-\frac{\partial L}{\partial q_k}=0}
$$

物理意义：只需一个标量函数 $L$ 即可给出全部运动方程。特例：自由粒子 $L=\tfrac12 m\dot x^2$，得 $m\ddot x=0$。

应用：单摆、双摆、弹簧振子。

**知识点 7.6 含速度相关势与广义势**

定义：若广义力可写为
$$
Q_k=\frac{d}{dt}\frac{\partial U}{\partial\dot q_k}-\frac{\partial U}{\partial q_k}
$$
则称 $U$ 为广义势，$L=T-U$，拉格朗日方程形式不变。

推导：将 $Q_k$ 表达式代入基本拉格朗日方程
$$
\frac{d}{dt}\frac{\partial T}{\partial\dot q_k}-\frac{\partial T}{\partial q_k}
=\frac{d}{dt}\frac{\partial U}{\partial\dot q_k}-\frac{\partial U}{\partial q_k}
$$
移项即得含 $L=T-U$ 的方程。

物理意义：洛伦兹力、电磁场均可纳入广义势。特例：带电粒子在电磁场中
$$
U=q\phi-q\vec A\cdot\vec v,\qquad
L=\tfrac12 m\vec v^2-q\phi+q\vec A\cdot\vec v
$$

应用：电磁约束系统、旋转参考系中的有效势。

### 第8章 哈密顿原理与哈密顿力学

**知识点 8.1 泛函、变分与欧拉-拉格朗日方程**

定义：泛函 $J[y]=\int_{x_1}^{x_2}F(x,y,y')\,dx$，其变分 $\delta J=0$ 给出极值曲线。

推导：令 $y(x,\alpha)=y(x)+\alpha\eta(x)$，$\eta(x_1)=\eta(x_2)=0$，则
$$
\frac{dJ}{d\alpha}=\int_{x_1}^{x_2}\left(\frac{\partial F}{\partial y}\eta+\frac{\partial F}{\partial y'}\eta'\right)dx
$$
对第二项分部积分：
$$
\int\frac{\partial F}{\partial y'}\eta'\,dx
=\left[\frac{\partial F}{\partial y'}\eta\right]_{x_1}^{x_2}
-\int\frac{d}{dx}\frac{\partial F}{\partial y'}\eta\,dx
$$
边界项为零，故
$$
\delta J=\int\left[\frac{\partial F}{\partial y}-\frac{d}{dx}\frac{\partial F}{\partial y'}\right]\eta\,dx=0
$$
由 $\eta$ 任意，得欧拉-拉格朗日方程
$$
\frac{\partial F}{\partial y}-\frac{d}{dx}\frac{\partial F}{\partial y'}=0
$$

物理意义：泛函极值条件是微分方程。特例：最速降线、测地线。

应用：变分法、最优控制、结构优化。

**知识点 8.2 哈密顿原理（最小作用量原理）**

定义：系统真实运动使作用量
$$
S=\int_{t_1}^{t_2}L(q,\dot q,t)\,dt
$$
取驻值，且 $\delta q(t_1)=\delta q(t_2)=0$。

推导：对多自由度情形
$$
\delta S=\int_{t_1}^{t_2}\sum_k\left(\frac{\partial L}{\partial q_k}\delta q_k+\frac{\partial L}{\partial\dot q_k}\delta\dot q_k\right)dt
$$
其中 $\delta\dot q_k=\dfrac{d}{dt}\delta q_k$，分部积分
$$
\int\frac{\partial L}{\partial\dot q_k}\frac{d}{dt}\delta q_k\,dt
=\left[\frac{\partial L}{\partial\dot q_k}\delta q_k\right]_{t_1}^{t_2}
-\int\frac{d}{dt}\frac{\partial L}{\partial\dot q_k}\delta q_k\,dt
$$
边界项为零，故
$$
\delta S=\int\sum_k\left[\frac{\partial L}{\partial q_k}-\frac{d}{dt}\frac{\partial L}{\partial\dot q_k}\right]\delta q_k\,dt=0
$$
由 $\delta q_k$ 独立，得拉格朗日方程。

物理意义：真实路径是作用量泛函的驻值路径。特例：自由粒子匀速直线运动使 $S$ 最小。

应用：场论、连续介质、量子力学的经典极限。

**知识点 8.3 勒让德变换与哈密顿函数**

定义：共轭动量 $p_k=\dfrac{\partial L}{\partial\dot q_k}$，哈密顿量
$$
H=\sum_k p_k\dot q_k-L
$$

推导：对 $H$ 求全微分
$$
dH=\sum_k(p_k\,d\dot q_k+\dot q_k\,dp_k)-\sum_k\left(\frac{\partial L}{\partial q_k}dq_k+\frac{\partial L}{\partial\dot q_k}d\dot q_k\right)-\frac{\partial L}{\partial t}dt
$$
由 $p_k=\partial L/\partial\dot q_k$ 消去 $d\dot q_k$ 项，并用 $\dfrac{\partial L}{\partial q_k}=\dot p_k$：
$$
dH=\sum_k(\dot q_k\,dp_k-\dot p_k\,dq_k)-\frac{\partial L}{\partial t}dt
$$
故
$$
\dot q_k=\frac{\partial H}{\partial p_k},\qquad
\dot p_k=-\frac{\partial H}{\partial q_k}
$$

物理意义：相空间 $(q,p)$ 上的正则方程。特例：保守系统 $H=T+V$ 为总能量。

应用：天体力学、统计力学、混沌分析。

**知识点 8.4 泊松括号与守恒量**

定义：对相空间函数 $f(q,p,t),g(q,p,t)$
$$
\{f,g\}=\sum_k\left(\frac{\partial f}{\partial q_k}\frac{\partial g}{\partial p_k}-\frac{\partial f}{\partial p_k}\frac{\partial g}{\partial q_k}\right)
$$

推导：对 $f$ 沿时间求全导：
$$
\frac{df}{dt}=\sum_k\left(\frac{\partial f}{\partial q_k}\dot q_k+\frac{\partial f}{\partial p_k}\dot p_k\right)+\frac{\partial f}{\partial t}
$$
代入正则方程：
$$
\frac{df}{dt}=\{f,H\}+\frac{\partial f}{\partial t}
$$
故 $f$ 守恒 $\iff \{f,H\}+\partial_t f=0$；$f$ 不显含时间时守恒 $\iff\{f,H\}=0$。

物理意义：泊松括号给出守恒量的代数判据。特例：$\{q_i,q_j\}=0,\ \{p_i,p_j\}=0,\ \{q_i,p_j\}=\delta_{ij}$。

应用：角动量守恒、辛几何、量子化对应。

**知识点 8.5 诺特定理（对称性与守恒律）**

定义：若作用量在某一连续变换下不变，则存在对应守恒量。

推导：考虑单参数变换 $q_k\to q_k+\epsilon\,\Delta q_k$，$\delta L=\epsilon\dfrac{d}{dt}G$。由欧拉-拉格朗日方程
$$
\frac{d}{dt}\frac{\partial L}{\partial\dot q_k}=\frac{\partial L}{\partial q_k}
$$
故
$$
\frac{d}{dt}\left(\sum_k\frac{\partial L}{\partial\dot q_k}\Delta q_k-G\right)=0
$$
括号内即为守恒量。

物理意义：时间平移不变 $\to$ 能量守恒；空间平移不变 $\to$ 动量守恒；旋转不变 $\to$ 角动量守恒。

特例：中心力场旋转对称给出角动量守恒 $L_z=\text{const}$。

应用：场论中的守恒流、粒子物理。

**知识点 8.6 正则变换与哈密顿-雅可比方程**

定义：保持正则方程形式的相空间变换 $(q,p)\to(Q,P)$ 称正则变换。

推导：要求
$$
p\,dq-H\,dt=P\,dQ-K\,dt+dF
$$
其中 $F(q,Q,t)$ 为生成函数。比较微分系数：
$$
p=\frac{\partial F}{\partial q},\qquad P=-\frac{\partial F}{\partial Q},\qquad K=H+\frac{\partial F}{\partial t}
$$
取 $F$ 使 $K=0$，即
$$
H\!\left(q,\frac{\partial F}{\partial q},t\right)+\frac{\partial F}{\partial t}=0
$$
此为哈密顿-雅可比方程。分离变量 $F=S(q,t)$ 得完整积分后，$Q=\beta=\text{const}$。

物理意义：把动力学问题化为求解一阶偏微分方程。特例：可分离系统中 $S=\sum W_k(q_k)-Et$，得 $H=E$。

应用：天体轨道、可积系统、几何光学类比。

### 第9章 应力分析

**知识点 9.1 连续介质假设与应力矢量**

定义：连续介质假设认为物质连续充满空间，任意点有确定的密度、位移、温度等场量。面力集度称应力矢量
$$
\vec t^{(n)}=\lim_{\Delta A\to0}\frac{\Delta\vec F}{\Delta A}
$$

推导：取包含点 $P$、法向 $\vec n$ 的面元，合力 $\Delta\vec F$，则 $\vec t^{(n)}$ 既依赖位置又依赖法向。

物理意义：应力是接触面力，具有方向面与方向矢双标号性质。特例：均匀拉伸时 $\vec t=\sigma\vec n$。

应用：强度校核、有限元应力输出。

**知识点 9.2 应力张量与柯西四面体定理**

定义：存在二阶张量 $\boldsymbol\sigma$ 使得
$$
\vec t^{(n)}=\boldsymbol\sigma^T\vec n,\qquad t_i=\sigma_{ij}n_j
$$
其中 $\sigma_{ij}$ 为作用在 $j$ 面、方向 $i$ 的应力分量。

推导：取四面体，三个坐标面与斜面围成，设斜面积 $\Delta A$、法向 $\vec n=(n_1,n_2,n_3)$，则坐标面积 $\Delta A_i=n_i\Delta A$。受力平衡：
$$
\vec t^{(n)}\Delta A+\sum_i(-\vec t^{(i)})\Delta A_i+\vec f\Delta V=0
$$
令四面体收缩 $\Delta V/\Delta A\to0$：
$$
\vec t^{(n)}=\sum_i n_i\vec t^{(i)}
=\sum_i n_i\sigma_{ij}\hat e_j
$$
即 $t_i=\sigma_{ij}n_j$。

物理意义：任意斜面上应力由三个坐标面上的应力完全确定。特例：$\vec n=\hat e_x$ 时 $\vec t=\sigma_{xx}\hat e_x+\sigma_{xy}\hat e_y+\sigma_{xz}\hat e_z$。

应用：斜截面应力、应力莫尔圆基础。

**知识点 9.3 应力张量的对称性与平衡方程**

定义：$\sigma_{ij}=\sigma_{ji}$。

推导（对称性）：取微元立方体，对形心求力矩，体力矩与惯性力矩为高阶小量。绕 $z$ 轴：
$$
(\sigma_{xy}\,dy\,dz)\frac{dx}{2}+(\sigma_{xy}\,dy\,dz)\frac{dx}{2}
-(\sigma_{yx}\,dx\,dz)\frac{dy}{2}-(\sigma_{yx}\,dx\,dz)\frac{dy}{2}=0
$$
化简得 $\sigma_{xy}=\sigma_{yx}$，其余同理。

推导（平衡方程）：对 $x$ 方向力平衡：
$$
\left(\sigma_{xx}+\frac{\partial\sigma_{xx}}{\partial x}dx\right)dy\,dz-\sigma_{xx}dy\,dz
+\left(\sigma_{yx}+\frac{\partial\sigma_{yx}}{\partial y}dy\right)dx\,dz-\sigma_{yx}dx\,dz
+\left(\sigma_{zx}+\frac{\partial\sigma_{zx}}{\partial z}dz\right)dx\,dy-\sigma_{zx}dx\,dy
+f_x\,dx\,dy\,dz=0
$$
化简：
$$
\frac{\partial\sigma_{xx}}{\partial x}+\frac{\partial\sigma_{yx}}{\partial y}+\frac{\partial\sigma_{zx}}{\partial z}+f_x=0
$$
张量形式
$$
\sigma_{ji,j}+f_i=0\quad\Longleftrightarrow\quad\nabla\cdot\boldsymbol\sigma+\vec f=0
$$

物理意义：平衡方程是动量守恒在静止介质中的体现；动力学情形右边加 $\rho\ddot u_i$。特例：无体力 $\nabla\cdot\boldsymbol\sigma=0$。

应用：平面问题、轴对称问题求解。

**知识点 9.4 主应力、应力不变量与莫尔圆**

定义：主应力是使 $\vec t^{(n)}\parallel\vec n$ 的法向，即
$$
t_i=\sigma_{ij}n_j=\sigma n_i
\;\Rightarrow\;
(\sigma_{ij}-\sigma\delta_{ij})n_j=0
$$

推导：非零解条件为特征方程
$$
\det(\sigma_{ij}-\sigma\delta_{ij})=0
$$
展开得
$$
\sigma^3-I_1\sigma^2+I_2\sigma-I_3=0
$$
其中不变量
$$
I_1=\sigma_{ii}=\sigma_x+\sigma_y+\sigma_z
$$
$$
I_2=\tfrac12(\sigma_{ii}\sigma_{jj}-\sigma_{ij}\sigma_{ji})
=\sigma_x\sigma_y+\sigma_y\sigma_z+\sigma_z\sigma_x-\tau_{xy}^2-\tau_{yz}^2-\tau_{zx}^2
$$
$$
I_3=\det(\sigma_{ij})
$$
主应力满足 $\sigma_1\ge\sigma_2\ge\sigma_3$。

平面应力莫尔圆：对与 $x$ 轴成 $\theta$ 的斜截面
$$
\sigma_\theta=\frac{\sigma_x+\sigma_y}{2}+\frac{\sigma_x-\sigma_y}{2}\cos2\theta+\tau_{xy}\sin2\theta
$$
$$
\tau_\theta=-\frac{\sigma_x-\sigma_y}{2}\sin2\theta+\tau_{xy}\cos2\theta
$$
消去 $\theta$：
$$
\left(\sigma_\theta-\frac{\sigma_x+\sigma_y}{2}\right)^2+\tau_\theta^2
=\left(\frac{\sigma_x-\sigma_y}{2}\right)^2+\tau_{xy}^2
$$
即莫尔圆半径
$$
R=\sqrt{\left(\frac{\sigma_x-\sigma_y}{2}\right)^2+\tau_{xy}^2}
$$

物理意义：主应力面上剪应力为零；莫尔圆给出所有斜截面应力。特例：纯剪切 $\sigma_x=\sigma_y=0,\tau_{xy}=\tau$，主应力 $\pm\tau$。

应用：压力容器、梁的应力状态、脆性材料破坏判据。

**知识点 9.5 应力偏张量与球张量分解**

定义：$\sigma_{ij}=s_{ij}+p\delta_{ij}$，其中平均应力
$$
p=\tfrac13\sigma_{kk}=-p_{hydro},\qquad s_{ij}=\sigma_{ij}-\tfrac13\sigma_{kk}\delta_{ij}
$$

推导：由 $s_{ij}$ 定义，迹
$$
s_{ii}=\sigma_{ii}-\tfrac13\sigma_{kk}\delta_{ii}=\sigma_{kk}-\tfrac13\sigma_{kk}\cdot3=0
$$
故偏张量无迹。

物理意义：球张量描述体积变化，偏张量描述形状畸变与塑性。特例：静水压力只改变体积不改变形状，不引起屈服。

应用：塑性屈服准则、岩土力学。

**知识点 9.6 应力边界条件与坐标系变换**

定义：在表面 $\Gamma$ 上，给定面力 $\bar t_i$：
$$
\sigma_{ij}n_j=\bar t_i
$$

推导坐标变换：设新坐标轴方向余弦 $a_{ij}$，则
$$
\sigma'_{ij}=a_{ik}a_{jl}\sigma_{kl}
$$
验证：在新坐标下 $t'_i=\sigma'_{ij}n'_j$，代入原始关系并利用 $a_{ik}\sigma_{kl}n_l=t_k a_{ik}$，得证。

物理意义：应力分量随坐标系变化，但 $I_1,I_2,I_3$ 不变。

特例：二维旋转 $\theta$：
$$
\sigma'_{xx}=\sigma_{xx}\cos^2\theta+\sigma_{yy}\sin^2\theta+2\sigma_{xy}\sin\theta\cos\theta
$$

应用：应力应变花测量、各向异性分析。

### 第10章 应变与几何方程

**知识点 10.1 位移场与格林应变张量**

定义：位移 $\vec u(\vec r)=\vec r'-\vec r$。考察微线元 $d\vec r$ 变形后 $d\vec r'$。

推导：$d\vec r'=d\vec r+d\vec u$，其中
$$
d\vec u=\nabla\vec u\cdot d\vec r,\qquad du_i=\frac{\partial u_i}{\partial x_j}dx_j
$$
线元长度平方差
$$
|d\vec r'|^2-|d\vec r|^2
=\left(\delta_{ij}+\partial_i u_j+\partial_j u_i+\partial_i u_k\partial_j u_k-\delta_{ij}\right)dx_idx_j
$$
定义为
$$
2\varepsilon_{ij}dx_idx_j,\qquad
\varepsilon_{ij}=\tfrac12(\partial_i u_j+\partial_j u_i+\partial_k u_i\partial_k u_j)
$$
即格林（有限）应变张量。

物理意义：$\varepsilon_{ij}$ 度量线元长度与夹角的真实改变。特例：小变形忽略位移梯度乘积。

应用：大变形橡胶、薄膜、几何非线性分析。

**知识点 10.2 小变形应变张量与转动张量**

定义：小变形假设 $|\partial_i u_j|\ll1$，则
$$
\varepsilon_{ij}=\tfrac12(\partial_i u_j+\partial_j u_i)
$$

推导：把位移梯度张量分解为对称与反对称部分
$$
\partial_i u_j=\frac{1}{2}(\partial_i u_j+\partial_j u_i)+\frac{1}{2}(\partial_i u_j-\partial_j u_i)
=\varepsilon_{ij}+\omega_{ij}
$$
反对称部分
$$
\omega_{ij}=\tfrac12(\partial_i u_j-\partial_j u_i)
$$
对应转动矢量 $\omega_k=-\tfrac12\epsilon_{kij}\omega_{ij}$，即 $\vec\omega=\tfrac12\nabla\times\vec u$。

物理意义：$\varepsilon_{ij}$ 描述纯变形，$\omega_{ij}$ 描述刚体转动（不变形）。特例：刚体平移应变全零；刚体小转动 $\vec u=\vec\Omega\times\vec r$ 时 $\varepsilon=0$。

应用：应变片测量、有限元应变计算。

**知识点 10.3 正应变与剪应变的几何解释**

定义：沿方向 $\vec n$ 的正应变 $\varepsilon_n=\varepsilon_{ij}n_in_j$；两正交方向 $\vec n,\vec m$ 间剪应变 $\gamma_{nm}=2\varepsilon_{ij}n_im_j$。

推导：线元 $d\vec r=ds\,\vec n$，变形后长度
$$
ds'^2=ds^2+2\varepsilon_{ij}ds\,n_i\,ds\,n_j
=ds^2(1+2\varepsilon_n)
$$
故 $ds'\approx ds(1+\varepsilon_n)$。夹角改变由两线元内积变化得出：
$$
\vec n'\cdot\vec m'-\vec n\cdot\vec m=2\varepsilon_{ij}n_im_j=\gamma_{nm}
$$

物理意义：$\varepsilon_n$ 是相对伸长；$\gamma$ 是直角改变量。特例：$\varepsilon_x=\partial u/\partial x$；$\gamma_{xy}=\partial u/\partial y+\partial v/\partial x$。

应用：应变花、梁的弯曲应变。

**知识点 10.4 应变协调方程（相容方程）**

定义：6 个应变分量由 3 个位移分量确定，故不独立，须满足相容条件。

推导：对 $\varepsilon_{ij}$ 求二阶导，构造
$$
\varepsilon_{ij,kl}+\varepsilon_{kl,ij}-\varepsilon_{ik,jl}-\varepsilon_{jl,ik}=0
$$
验证：以 $\varepsilon_{xx}=u_{x,x}$，$\varepsilon_{xy}=\tfrac12(u_{x,y}+u_{y,x})$，$\varepsilon_{yy}=u_{y,y}$ 为例，二维有
$$
\frac{\partial^2\varepsilon_{xx}}{\partial y^2}+\frac{\partial^2\varepsilon_{yy}}{\partial x^2}
=2\frac{\partial^2\varepsilon_{xy}}{\partial x\partial y}
$$
这是因为
$$
\partial_y^2 u_{x,x}+\partial_x^2 u_{y,y}
=\partial_x\partial_y(u_{x,y}+u_{y,x})
=2\partial_x\partial_y\varepsilon_{xy}
$$

物理意义：保证由应变积分得到单值连续位移场。特例：平面问题仅一个独立相容方程。

应用：光弹性、逆解法、应力函数法。

**知识点 10.5 主应变与应变不变量**

定义：主应变由
$$
(\varepsilon_{ij}-\varepsilon\delta_{ij})n_j=0
$$
的特征值给出，满足
$$
\varepsilon^3-J_1\varepsilon^2+J_2\varepsilon-J_3=0
$$
其中
$$
J_1=\varepsilon_{ii},\quad
J_2=\tfrac12(\varepsilon_{ii}\varepsilon_{jj}-\varepsilon_{ij}\varepsilon_{ji}),\quad
J_3=\det(\varepsilon_{ij})
$$

推导：与应力特征方程完全类似，对称张量存在三个实特征值及正交特征方向。

物理意义：主应变方向相互正交且剪应变为零。特例：体积应变
$$
\theta=\varepsilon_{ii}=\varepsilon_1+\varepsilon_2+\varepsilon_3=\frac{\Delta V}{V_0}
$$

应用：塑性成形、体积模量、材料可压缩性。

**知识点 10.6 应变率张量与速度梯度分解**

定义：流体或大变形中，速度梯度
$$
\partial_i v_j=\tfrac12(\partial_i v_j+\partial_j v_i)+\tfrac12(\partial_i v_j-\partial_j v_i)
=D_{ij}+W_{ij}
$$

推导：对称部分（应变率张量）
$$
D_{ij}=\tfrac12(\partial_i v_j+\partial_j v_i)
$$
反对称部分对应涡量 $\vec\omega=\nabla\times\vec v$，且
$$
W_{ij}=\tfrac12(\partial_i v_j-\partial_j v_i)=-\tfrac12\epsilon_{ijk}\omega_k
$$
验证：$\epsilon_{ijk}W_{ij}$ 给出涡量分量。

物理意义：$D_{ij}$ 描述流体质点变形速率，$W_{ij}$ 描述刚体旋转速率。特例：刚体旋转 $\vec v=\vec\Omega\times\vec r$，$D=0$。

应用：黏性流体应力本构、湍流统计。

### 第11章 本构关系与弹性力学基本方程

**知识点 11.1 广义胡克定律与各向同性弹性**

定义：线弹性本构 $\sigma_{ij}=C_{ijkl}\varepsilon_{kl}$；各向同性时
$$
\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij}
$$

推导：各向同性四阶张量只能由 $\delta_{ij}$ 组合：
$$
C_{ijkl}=\lambda\delta_{ij}\delta_{kl}+\mu(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk})+\nu\delta_{ij}\delta_{kl}\text{ 组合}
$$
代入 $\sigma_{ij}=C_{ijkl}\varepsilon_{kl}$ 并用对称性，得上述拉梅形式，$\lambda,\mu$ 为拉梅常数。

逆形式：
$$
\varepsilon_{ij}=\frac{1+\nu}{E}\sigma_{ij}-\frac{\nu}{E}\sigma_{kk}\delta_{ij}
$$
其中 $E$ 杨氏模量、$\nu$ 泊松比，关系
$$
\lambda=\frac{E\nu}{(1+\nu)(1-2\nu)},\qquad
\mu=G=\frac{E}{2(1+\nu)}
$$

物理意义：$\lambda,\mu$ 分别反映体积与剪切弹性。特例：单轴拉伸得 $\sigma=E\varepsilon$，横向收缩 $\nu$。

应用：结构刚度、地震波速（$v_p,v_s$）。

**知识点 11.2 弹性力学基本方程与纳维方程**

定义：15 个未知量（6 应力+6 应变+3 位移）由 15 个方程封闭：平衡 3、几何 6、本构 6。

推导：将本构代入平衡，并用几何表示应变：
$$
\sigma_{ij}=\lambda u_{k,k}\delta_{ij}+\mu(u_{i,j}+u_{j,i})
$$
代入 $\sigma_{ij,j}+f_i=0$：
$$
\lambda u_{k,ki}+\mu(u_{i,jj}+u_{j,ij})+f_i=0
$$
合并
$$
(\lambda+\mu)u_{j,ji}+\mu u_{i,jj}+f_i=0
$$
矢量形式
$$
(\lambda+\mu)\nabla(\nabla\cdot\vec u)+\mu\nabla^2\vec u+\vec f=0
$$
动力学情形右边加 $\rho\ddot{\vec u}$，称纳维方程。

物理意义：位移法的控制方程，是弹性力学有限元基础。特例：无体力时
$$
(\lambda+\mu)\nabla(\nabla\cdot\vec u)+\mu\nabla^2\vec u=0
$$

应用：弹性波、地基沉降、接触问题。

**知识点 11.3 弹性常数关系与工程常数**

定义：拉伸试验给出 $E,\nu$；体积模量 $K=-\dfrac{p}{\Delta V/V}$；剪切模量 $G$。

推导：体积应变
$$
\theta=\varepsilon_{kk}=\frac{1-2\nu}{E}\sigma_{kk}
$$
故
$$
K=\frac{E}{3(1-2\nu)}=\lambda+\tfrac23\mu
$$
又 $\mu=G=\dfrac{E}{2(1+\nu)}$。

物理意义：$K$ 描述抗体积压缩能力。特例：$\nu=0.5$ 时 $K\to\infty$，材料不可压（橡胶、流体）。

应用：高压实验、声速测量：
$$
v_p=\sqrt{\frac{K+4G/3}{\rho}},\qquad v_s=\sqrt{\frac{G}{\rho}}
$$

**知识点 11.4 平面应力与平面应变**

定义：平面应力 $\sigma_z=\tau_{xz}=\tau_{yz}=0$；平面应变 $\varepsilon_z=\gamma_{xz}=\gamma_{yz}=0$。

推导平面应力：由本构
$$
\varepsilon_x=\frac1E(\sigma_x-\nu\sigma_y),\quad
\varepsilon_y=\frac1E(\sigma_y-\nu\sigma_x),\quad
\gamma_{xy}=\frac{\tau_{xy}}{G}
$$
平面应变的应力：
$$
\sigma_x=\frac{E}{(1+\nu)(1-2\nu)}[(1-\nu)\varepsilon_x+\nu\varepsilon_y]
$$
$$
\sigma_y=\frac{E}{(1+\nu)(1-2\nu)}[\nu\varepsilon_x+(1-\nu)\varepsilon_y]
$$
$$
\tau_{xy}=G\gamma_{xy}
$$

物理意义：平面应力适用于薄板，平面应变适用于长柱体。特例：平面应力下 $\varepsilon_z=-\dfrac{\nu}{E}(\sigma_x+\sigma_y)$。

应用：薄壁结构、坝体、隧道。

**知识点 11.5 艾里应力函数与双调和方程**

定义：无体力平面问题中令
$$
\sigma_x=\frac{\partial^2\Phi}{\partial y^2},\quad
\sigma_y=\frac{\partial^2\Phi}{\partial x^2},\quad
\tau_{xy}=-\frac{\partial^2\Phi}{\partial x\partial y}
$$

推导：平衡方程
$$
\frac{\partial\sigma_x}{\partial x}+\frac{\partial\tau_{xy}}{\partial y}=0,\quad
\frac{\partial\tau_{xy}}{\partial x}+\frac{\partial\sigma_y}{\partial y}=0
$$
代入 $\Phi$ 自动满足。相容方程
$$
\left(\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}\right)(\sigma_x+\sigma_y)=0
$$
而 $\sigma_x+\sigma_y=\nabla^2\Phi$，故
$$
\nabla^4\Phi=0
$$
即双调和方程。

物理意义：把求解应力化为求双调和函数。特例：多项式 $\Phi=ax^2+bxy+cy^2$ 给均匀应力场。

应用：孔边应力集中、裂纹尖端场、梁的解。

**知识点 11.6 强度理论与屈服准则**

定义：复杂应力下材料失效判据。

推导（最大剪应力 Treca）：
$$
\tau_{max}=\frac{\sigma_1-\sigma_3}{2}\le\frac{\sigma_s}{2}
\;\Rightarrow\;
\sigma_1-\sigma_3\le\sigma_s
$$
推导（von Mises）：畸变能密度
$$
v_d=\frac{1+\nu}{6E}\left[(\sigma_1-\sigma_2)^2+(\sigma_2-\sigma_3)^2+(\sigma_3-\sigma_1)^2\right]
$$
单轴屈服时 $\sigma_1=\sigma_s,\sigma_2=\sigma_3=0$，$v_d^0=\dfrac{1+\nu}{3E}\sigma_s^2$。令 $v_d=v_d^0$：
$$
\sqrt{\tfrac12[(\sigma_1-\sigma_2)^2+(\sigma_2-\sigma_3)^2+(\sigma_3-\sigma_1)^2]}=\sigma_s
$$
等效应力即 von Mises 应力。

物理意义：Tresca 偏保守；von Mises 与多数金属吻合。特例：纯剪切屈服 $\tau_s=\sigma_s/\sqrt3$。

应用：压力容器、管道、金属成形。

**知识点 11.7 热弹性本构**

定义：温度变化引起热应变 $\varepsilon^{th}_{ij}=\alpha\Delta T\,\delta_{ij}$。

推导：总应变 = 弹性应变 + 热应变：
$$
\varepsilon_{ij}=\frac{1+\nu}{E}\sigma_{ij}-\frac{\nu}{E}\sigma_{kk}\delta_{ij}+\alpha\Delta T\,\delta_{ij}
$$
逆解
$$
\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij}-(3\lambda+2\mu)\alpha\Delta T\,\delta_{ij}
$$

物理意义：温升产生膨胀，受约束时产生热应力。特例：全约束杆升温 $\sigma=-E\alpha\Delta T/(1-2\nu)$（三维）。

应用：热障涂层、焊接残余应力、航天热结构。

### 第12章 能量原理与变分法

**知识点 12.1 应变能与余应变能**

定义：弹性体单位体积应变能
$$
w=\int_0^{\varepsilon_{ij}}\sigma_{ij}\,d\varepsilon_{ij}
=\tfrac12\sigma_{ij}\varepsilon_{ij}
$$
（线弹性）。总应变能 $U=\int_V w\,dV$。

推导：对单轴弹簧，$dU=\sigma A\,d(\varepsilon L)=\sigma\varepsilon\,V$，积分得 $\tfrac12\sigma\varepsilon V$。三维由功共轭推广，线弹性时二次型。

余应变能
$$
w^*=\sigma_{ij}\varepsilon_{ij}-w=\tfrac12\sigma_{ij}\varepsilon_{ij}=w
$$
线弹性时二者相等。

物理意义：应变能储存于弹性体，可释放做功。特例：受拉杆 $U=\dfrac{P^2L}{2EA}=\dfrac12P\delta$。

应用：卡氏定理、有限元刚度矩阵。

**知识点 12.2 卡氏第一定理与第二定理**

定义：若应变能表为位移的函数 $U(\delta_i)$，则 $\dfrac{\partial U}{\partial\delta_i}=P_i$（第一定理）；若表为载荷的函数 $U(P_i)$，则 $\dfrac{\partial U}{\partial P_i}=\delta_i$（第二定理）。

推导（第二定理）：总势能 $\Pi=U-\sum P_i\delta_i$。对于线弹性，$U=\tfrac12\sum P_i\delta_i$。由平衡条件 $\partial\Pi/\partial\delta_i=0$ 及 $P_i=\partial U/\partial\delta_i$，得
$$
\delta_i=\frac{\partial U}{\partial P_i}
$$
更一般地，对 $U^*=U$，可证明。

物理意义：位移由能量对载荷偏导给出。特例：悬臂梁端部力 $P$，$U=\dfrac{P^2L^3}{6EI}$，则 $\delta=\dfrac{\partial U}{\partial P}=\dfrac{PL^3}{3EI}$。

应用：静定/静不定结构位移计算、单位载荷法。

**知识点 12.3 虚功原理与虚余功原理**

定义：对变形体，外力虚功等于内力虚功。

推导：平衡微元上 $\sigma_{ij,j}+f_i=0$，乘虚位移 $\delta u_i$ 并积分：
$$
\int_V(\sigma_{ij,j}+f_i)\delta u_i\,dV=0
$$
分部积分
$$
\int_V\sigma_{ij,j}\delta u_i\,dV
=\int_S\sigma_{ij}n_j\delta u_i\,dS-\int_V\sigma_{ij}\delta\varepsilon_{ij}\,dV
$$
代入得
$$
\int_V\sigma_{ij}\delta\varepsilon_{ij}\,dV
=\int_S t_i\delta u_i\,dS+\int_V f_i\delta u_i\,dV
$$
即虚功方程。

物理意义：任意满足位移边界条件的虚位移均成立。特例：刚体虚位移给 $\int t\,dS+\int f\,dV=0$（整体平衡）。

应用：有限元、能量法、极限分析。

**知识点 12.4 最小势能原理**

定义：总势能
$$
\Pi=U-W_{ext}=U-\int_S\bar t_iu_i\,dS-\int_V f_iu_i\,dV
$$
真实位移使 $\Pi$ 取极小。

推导：令 $\delta\Pi=0$：
$$
\delta U-\int\bar t_i\delta u_i\,dS-\int f_i\delta u_i\,dV=0
$$
而 $\delta U=\int\sigma_{ij}\delta\varepsilon_{ij}dV$，由虚功原理恰为该式，故真实解满足 $\delta\Pi=0$。二阶变分 $\delta^2\Pi=\int\delta\varepsilon_{ij}C_{ijkl}\delta\varepsilon_{kl}dV>0$（$C$ 正定），故为极小。

物理意义：弹性平衡是能量极小态。特例：线弹性、小变形、保守载荷。

应用：瑞利-里兹法、有限元变分基础。

**知识点 12.5 瑞利-里兹法与伽辽金法**

定义：假设近似位移 $u_i\approx\sum_{m}a_m\varphi^{(m)}_i(\vec r)$，$\varphi^{(m)}$ 满足位移边界条件。

推导（Ritz）：代入 $\Pi$ 并对系数求极值：
$$
\frac{\partial\Pi}{\partial a_m}=0,\qquad m=1,2,\dots
$$
得线性方程组 $K_{mn}a_n=f_m$，其中
$$
K_{mn}=\int_V \varepsilon_{ij}(\varphi^{(m)})C_{ijkl}\varepsilon_{kl}(\varphi^{(n)})dV
$$
（Galerkin）：控制方程残差 $R=L(u)-f$，要求与试函数正交：
$$
\int_V R\,\psi_m\,dV=0
$$

物理意义：将无限维问题降为有限维代数问题。特例：取三角函数作试函数可得梁固有频率近似。

应用：板壳分析、特征值问题、有限元。

**知识点 12.6 能量法求临界载荷与稳定性**

定义：用势能二阶变分判断平衡稳定性。

推导：总势能
$$
\Pi=U_{bend}-P\Delta
$$
柱屈曲取挠曲线 $w(x)=a\phi(x)$，则
$$
\Pi(a)=\frac{EI}{2}\int_0^L(\phi'')^2dx\,a^2-\frac{P}{2}\int_0^L(\phi')^2dx\,a^2
$$
平衡 $\partial\Pi/\partial a=0$ 给 $a=0$（直线）或
$$
P_{cr}=\frac{EI\int(\phi'')^2dx}{\int(\phi')^2dx}
$$
即瑞利商。由 $\phi$ 的最优选择得最小临界载荷。

物理意义：屈曲是势能由正定变为非正定的分岔点。特例：取 $\phi=\sin(\pi x/L)$ 得
$$
P_{cr}=\frac{\pi^2EI}{L^2}
$$
与应用精确解一致。

应用：压杆、薄板、壳体稳定性。

### 第13章 杆件强度、刚度与稳定性

**知识点 13.1 轴向拉压与剪切**

定义：杆件轴向受力时横截面上正应力
$$
\sigma=\frac{N}{A}
$$

推导：取截面法，轴力 $N$ 平衡外力；由圣维南原理，远离加载端的截面上应力均匀：
$$
N=\int_A\sigma\,dA=\sigma A
$$
变形
$$
\Delta L=\int_0^L\frac{N(x)}{EA}dx,\qquad \varepsilon=\frac{\Delta L}{L}=\frac{\sigma}{E}
$$
泊松效应
$$
\varepsilon_{lat}=-\nu\varepsilon_{axial}
$$

物理意义：$EA$ 为拉压刚度。特例：阶梯杆总伸长 $\sum N_iL_i/(E_iA_i)$。

应用：桁架、螺栓、拉杆。

**知识点 13.2 圆轴扭转**

定义：圆轴受扭矩 $T$，剪应力
$$
\tau=\frac{T\rho}{I_p},\qquad I_p=\frac{\pi d^4}{32}
$$

推导：平面假设（横截面保持平面且无翘曲）、相邻截面绕轴相对转动 $\varphi$。距轴 $\rho$ 处剪应变
$$
\gamma=\rho\frac{d\varphi}{dx}
$$
线弹性 $\tau=G\gamma=G\rho\dfrac{d\varphi}{dx}$。扭矩平衡
$$
T=\int_A\tau\rho\,dA=G\frac{d\varphi}{dx}\int_A\rho^2dA=GI_p\frac{d\varphi}{dx}
$$
故
$$
\frac{d\varphi}{dx}=\frac{T}{GI_p},\qquad
\tau=\frac{T\rho}{I_p},\qquad
\tau_{max}=\frac{T}{W_p},\ W_p=\frac{\pi d^3}{16}
$$
扭转角
$$
\varphi=\int_0^L\frac{T(x)}{GI_p}dx
$$

物理意义：$GI_p$ 为扭转刚度。特例：空心圆轴 $I_p=\dfrac{\pi(D^4-d^4)}{32}$。

应用：传动轴、弹簧、扭转试验。

**知识点 13.3 梁的弯曲正应力与挠曲线**

定义：纯弯曲正应力
$$
\sigma=\frac{My}{I_z}
$$

推导：平面假设与单向受力假设，应变 $\varepsilon=\dfrac{y}{\rho}$，由胡克定律 $\sigma=E\dfrac{y}{\rho}$。轴力与弯矩平衡：
$$
N=\int_A\sigma\,dA=\frac{E}{\rho}\int_A y\,dA=0\ \Rightarrow\ \text{中性轴过形心}
$$
$$
M=\int_A\sigma y\,dA=\frac{E}{\rho}\int_A y^2dA=\frac{EI_z}{\rho}
$$
故曲率
$$
\frac{1}{\rho}=\frac{M}{EI_z},\qquad
\sigma=\frac{My}{I_z},\qquad
\sigma_{max}=\frac{M}{W_z}
$$
小变形挠曲线
$$
\frac{d^2w}{dx^2}=\frac{M(x)}{EI_z}
$$

物理意义：$EI_z$ 为弯曲刚度，$W_z$ 为抗弯截面系数。特例：矩形 $I_z=bh^3/12$。

应用：梁的设计、挠度校核、超静定梁。

**知识点 13.4 梁的剪应力与组合变形**

定义：横截面剪应力（矩形截面）
$$
\tau=\frac{FS_z^*}{I_zb}
$$

推导：取微段 $dx$，距中性轴 $y$ 以上部分受力差
$$
dT=\frac{dM}{dx}\frac{S_z^*}{I_z}dx=\frac{F S_z^*}{I_z}dx
$$
该力由剪应力平衡 $\tau b\,dx$，故
$$
\tau=\frac{FS_z^*}{I_zb}
$$
其中 $S_z^*$ 为面积矩。

组合变形：拉弯
$$
\sigma=\frac{N}{A}+\frac{M y}{I}
$$
弯扭
$$
\sigma=\frac{M}{W},\qquad\tau=\frac{T}{W_p}
$$
按第三强度理论
$$
\sigma_{eq4}=\sqrt{\sigma^2+3\tau^2}
$$

物理意义：剪应力沿高度抛物线分布，中性轴处最大。特例：圆形截面 $\tau_{max}=\dfrac{4F}{3A}$。

应用：吊车梁、传动轴、桥梁。

**知识点 13.5 压杆稳定的欧拉公式**

定义：细长压杆临界压力 $P_{cr}$。

推导：设两端铰支，轴力 $P$ 在挠曲后产生附加弯矩 $Pw$，微段平衡：
$$
EI\frac{d^2w}{dx^2}=-Pw
$$
令 $k^2=P/(EI)$：
$$
w''+k^2w=0
$$
通解
$$
w=A\sin kx+B\cos kx
$$
边界条件 $w(0)=0\Rightarrow B=0$；$w(L)=0\Rightarrow A\sin kL=0$。非平凡解 $\sin kL=0$：
$$
kL=n\pi\ \Rightarrow\ P=\frac{n^2\pi^2EI}{L^2}
$$
最小临界力（$n=1$）：
$$
P_{cr}=\frac{\pi^2EI}{L^2}
$$
一般约束用有效长度
$$
P_{cr}=\frac{\pi^2EI}{(\mu L)^2},\qquad \mu L=L_{eff}
$$
两端铰支 $\mu=1$、两端固定 $\mu=0.5$、一端固定一端自由 $\mu=2$、一端固定一端铰支 $\mu=0.7$。

物理意义：屈曲在压应力远低于屈服时发生。特例：柔度 $\lambda=\mu L/i$，大柔度用欧拉公式。

应用：立柱、桁架压杆、框架稳定。

**知识点 13.6 稳定性判据与能量法**

定义：对任意挠曲线，瑞利商给出临界载荷上界。

推导：弯曲应变能
$$
\Delta U=\frac12\int_0^L EI(w'')^2dx
$$
载荷势能变化
$$
\Delta V=-\frac12\int_0^L P(w')^2dx
$$
临界条件 $\Delta U+\Delta V=0$：
$$
P_{cr}=\frac{\int_0^L EI(w'')^2dx}{\int_0^L (w')^2dx}
$$
对精确 $w$ 得精确 $P_{cr}$，对试函数得上界。

物理意义：能量法避免求解变系数微分方程。特例：悬臂杆取 $w=\delta(1-\cos\frac{\pi x}{2L})$ 得 $P_{cr}\approx\dfrac{\pi^2EI}{4L^2}$。

应用：变截面杆、组合杆、复合材料杆稳定。

**知识点 13.7 疲劳与应力集中**

定义：交变应力下裂纹萌生与扩展导致失效。应力集中系数
$$
K_t=\frac{\sigma_{max}}{\sigma_{nom}}
$$

推导：对椭圆孔，Inglis 解给出孔边最大应力
$$
\sigma_{max}=\sigma_0\left(1+2\frac{a}{b}\right)
$$
裂纹（$b\to0$）时 $K_t\to\infty$。断裂力学用应力强度因子
$$
K_I=Y\sigma\sqrt{\pi a}
$$
裂纹扩展 Paris 律
$$
\frac{da}{dN}=C(\Delta K)^m
$$

物理意义：缺口处局部应力远高于名义应力。特例：圆孔 $\sigma_{max}=3\sigma_0$。

应用：飞机蒙皮、焊接接头、压力容器寿命评估。

### 第14章 连续介质与流体运动学

**知识点 14.1 连续介质假设与描述方法**

定义：流体作为连续介质，用场量 $\rho(\vec r,t),\vec v(\vec r,t),p(\vec r,t)$ 描述。拉格朗日法跟踪质点，欧拉法观察空间点。

推导：两种描述的换算。设拉格朗日轨迹 $\vec r=\vec r(\vec r_0,t)$，物理量 $f$。质点导数（物质导数）
$$
\frac{Df}{Dt}=\left.\frac{df}{dt}\right|_{trajectory}
=\frac{\partial f}{\partial t}+\frac{\partial f}{\partial x_i}\frac{dx_i}{dt}
=\frac{\partial f}{\partial t}+(\vec v\cdot\nabla)f
$$

物理意义：$D/Dt$ 是随流体质点的变化率。特例：定常场 $\partial_t=0$，但质点仍可因对流改变 $f$。

应用：气象、叶轮机械、河流。

**知识点 14.2 连续性方程**

定义：质量守恒。

推导：取控制体 $V$，质量变化率等于净流入：
$$
\frac{d}{dt}\int_V\rho\,dV=-\oint_S\rho\vec v\cdot\vec n\,dS
$$
对左侧用欧拉法把 $V$ 固定：
$$
\int_V\frac{\partial\rho}{\partial t}dV=-\int_V\nabla\cdot(\rho\vec v)dV
$$
由体积 $V$ 任意：
$$
\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\vec v)=0
$$
展开
$$
\frac{\partial\rho}{\partial t}+\vec v\cdot\nabla\rho+\rho\nabla\cdot\vec v=0
\;\Rightarrow\;
\frac{D\rho}{Dt}+\rho\nabla\cdot\vec v=0
$$
不可压时 $D\rho/Dt=0$，得
$$
\nabla\cdot\vec v=0
$$

物理意义：连续性方程是质量守恒的微分表述。特例：定常一维管流 $\rho_1v_1A_1=\rho_2v_2A_2$。

应用：管道、喷嘴、气象、血液流。

**知识点 14.3 流线与迹线、散度与旋度**

定义：流线满足 $d\vec r\times\vec v=0$，即
$$
\frac{dx}{v_x}=\frac{dy}{v_y}=\frac{dz}{v_z}
$$
迹线是质点轨迹。定常流中二者重合。

推导：$dx_i/ds=v_i$ 积分得流线。涡量
$$
\vec\omega=\nabla\times\vec v
$$
无旋流 $\vec\omega=0$ 时存在速度势 $\vec v=\nabla\phi$。

物理意义：流线直观显示流速方向，涡量描述局部旋转。特例：点涡 $\vec v=\dfrac{\Gamma}{2\pi r}\hat e_\theta$，$\nabla\times\vec v=0$（除原点），环量 $\Gamma=\oint\vec v\cdot d\vec l$。

应用：机翼绕流、旋涡、流场可视化。

**知识点 14.4 雷诺输运定理**

定义：系统（物质体）物理量变化率与控制体关系。

推导：设 $B=\int_{sys}\beta\,dm$（$\beta$ 为单位质量物理量），则
$$
\frac{dB}{dt}=\frac{d}{dt}\int_{CV}\rho\beta\,dV+\oint_{CS}\rho\beta\vec v_r\cdot\vec n\,dS
$$
更常用的固定控制体形式：
$$
\frac{dB_{sys}}{dt}
=\int_{CV}\frac{\partial(\rho\beta)}{\partial t}dV+\oint_{CS}\rho\beta\vec v\cdot\vec n\,dS
$$
推导思路：物质体在 $t$ 与 $t+dt$ 的重叠，令 $dt\to0$，体积分的局部变化加上通过表面的通量。

物理意义：把随体导数转化为欧拉场积分。特例：$\beta=1$ 给连续性方程；$\beta=\vec v$ 给动量方程。

应用：动量方程、能量方程、涡轮机分析。

**知识点 14.5 变形率与涡量方程（运动学）**

定义：流体微团速度梯度分解
$$
\partial_i v_j=D_{ij}+W_{ij}
$$
应变率张量 $D_{ij}$ 与自旋张量 $W_{ij}$。

推导：流体微团相对运动
$$
dv_i=\partial_j v_i\,dx_j
$$
分解后，对称部分导致线元伸缩与剪切，反对称部分导致刚体旋转 $\vec\omega/2$。

物理意义：$D_{ij}$ 决定黏性应力；$\vec\omega$ 决定涡动力学。特例：剪切流 $v_x=\dot\gamma y$ 有 $D_{xy}=\dot\gamma/2$、$\omega_z=-\dot\gamma/2$。

应用：混合、湍流、流变学。

**知识点 14.6 环量与开尔文定理**

定义：速度环量 $\Gamma=\oint_C\vec v\cdot d\vec l$。

推导：对时间求导
$$
\frac{d\Gamma}{dt}=\oint_C\frac{D\vec v}{Dt}\cdot d\vec l+\oint_C\vec v\cdot\frac{d}{dt}(d\vec l)
$$
第二项可证为零（闭合线元导数积分为零），故
$$
\frac{d\Gamma}{dt}=\oint_C\frac{D\vec v}{Dt}\cdot d\vec l
$$
对理想正压流体、有势力，$\dfrac{D\vec v}{Dt}=-\nabla(p/\rho+\Pi)$ 为梯度，闭合积分为零：
$$
\frac{d\Gamma}{dt}=0
$$
即开尔文环量守恒定理。

物理意义：无黏、正压、有势外力下涡量不能产生或消灭。特例：从静止启动的理想流始终无旋。

应用：机翼起动涡、涡环、超流。

### 第15章 流体动力学与纳维-斯托克斯方程

**知识点 15.1 欧拉方程与伯努利方程**

定义：理想流体（无黏）运动方程。

推导：对流体微元用牛顿第二定律，表面力只有压力：
$$
\rho\frac{D\vec v}{Dt}=-\nabla p+\rho\vec g
$$
写成分量
$$
\rho\left(\frac{\partial v_i}{\partial t}+v_j\partial_j v_i\right)
=-\partial_i p+\rho g_i
$$
沿流线积分（定常、不可压、有势体力）：
$$
\rho\,d\!\left(\frac{v^2}{2}\right)=-dp-\rho\,d\Pi
$$
积分得
$$
\frac{p}{\rho}+\frac{v^2}{2}+\Pi=\text{const}
$$
重力场 $\Pi=gz$：
$$
\frac{p}{\rho}+\frac{v^2}{2}+gz=\text{const}
$$

物理意义：压力能、动能、位能之和沿流线守恒。特例：托里拆利出流 $v=\sqrt{2gh}$。

应用：皮托管、文丘里管、机翼升力。

**知识点 15.2 牛顿流体本构与纳维-斯托克斯方程**

定义：牛顿流体应力
$$
\sigma_{ij}=-p\delta_{ij}+\tau_{ij},\qquad
\tau_{ij}=2\mu D_{ij}+\lambda_v D_{kk}\delta_{ij}
$$
通常取体积黏度 $\lambda_v=-\tfrac23\mu$ 使 $\tau_{kk}=0$，不可压时 $\tau_{ij}=2\mu D_{ij}$。

推导：把 $\tau_{ij}$ 代入柯西动量方程
$$
\rho\frac{Dv_i}{Dt}=\partial_j\sigma_{ij}+\rho f_i
$$
不可压（$\nabla\cdot\vec v=0$）：
$$
\partial_j\tau_{ij}=\mu\partial_j(\partial_i v_j+\partial_j v_i)
=\mu(\partial_i\partial_jv_j+\partial_j^2v_i)=\mu\nabla^2v_i
$$
故
$$
\rho\left(\frac{\partial\vec v}{\partial t}+(\vec v\cdot\nabla)\vec v\right)
=-\nabla p+\mu\nabla^2\vec v+\rho\vec f
$$
即可压缩情形
$$
\rho\frac{D\vec v}{Dt}=-\nabla p+\mu\nabla^2\vec v+\left(\lambda_v+\mu\right)\nabla(\nabla\cdot\vec v)+\rho\vec f
$$

物理意义：黏性使动量扩散，非线性对流项导致湍流。特例：$\mu=0$ 退化为欧拉方程。

应用：CFD、润滑、管道流动、生物流体。

**知识点 15.3 库埃特流与泊肃叶流**

定义：两平行板间与管道内黏性流动。

推导（平面库埃特）：不可压、定常、全发展，$v=(u(y),0,0)$，连续方程自动满足。N-S 的 $x$ 分量
$$
0=\mu\frac{d^2u}{dy^2}
\;\Rightarrow\;
u=C_1y+C_2
$$
边界 $u(0)=0,u(h)=U$：
$$
u(y)=\frac{U}{h}y
$$
剪切应力 $\tau=\mu U/h$。

推导（泊肃叶，圆管）：轴向速度 $u(r)$，N-S 的 $z$ 分量
$$
0=-\frac{dp}{dz}+\mu\frac{1}{r}\frac{d}{dr}\left(r\frac{du}{dr}\right)
$$
积分一次
$$
r\frac{du}{dr}=\frac{1}{2\mu}\frac{dp}{dz}r^2+C_1
$$
轴心有限 $C_1=0$；壁面 $u(R)=0$：
$$
u(r)=\frac{1}{4\mu}\left(-\frac{dp}{dz}\right)(R^2-r^2)
$$
流量
$$
Q=\int_0^R u\,2\pi r\,dr
=\frac{\pi R^4}{8\mu}\left(-\frac{dp}{dz}\right)
$$
平均速度 $v_m=Q/(\pi R^2)=\dfrac{R^2}{8\mu}\dfrac{\Delta p}{L}$。

物理意义：抛物型速度剖面。特例：达西-魏斯巴赫摩擦系数 $f=64/Re$。

应用：润滑、毛细血管、液压管路。

**知识点 15.4 涡量方程与亥姆霍兹定理**

定义：对 N-S 取旋度得涡量动力学方程。

推导：利用矢量恒等式
$$
(\vec v\cdot\nabla)\vec v=\nabla\left(\frac{v^2}{2}\right)-\vec v\times(\nabla\times\vec v)
$$
对动量方程取旋度，压力梯度旋度为零（正压），得
$$
\frac{\partial\vec\omega}{\partial t}+\nabla\times(\vec\omega\times\vec v)
=\frac{\mu}{\rho}\nabla^2\vec\omega
$$
再用
$$
\nabla\times(\vec\omega\times\vec v)
=(\vec v\cdot\nabla)\vec\omega-(\vec\omega\cdot\nabla)\vec v+\vec\omega(\nabla\cdot\vec v)
$$
不可压时
$$
\frac{D\vec\omega}{Dt}
=(\vec\omega\cdot\nabla)\vec v+\nu\nabla^2\vec\omega
$$
其中涡拉伸项 $(\vec\omega\cdot\nabla)\vec v$ 可放大涡量（湍流能量级串），$\nu\nabla^2\vec\omega$ 为黏性扩散。

物理意义：涡量因拉伸与黏性变化。特例：二维流 $\vec\omega\perp\nabla\vec v$，涡拉伸为零，仅扩散。

应用：湍流、旋涡分离、涡方法。

**知识点 15.5 无量纲化与相似准则**

定义：用特征量 $L,U,\rho,\mu$ 无量纲化，得到雷诺数等。

推导：设 $\vec v=U\vec v^*$，$\vec r=L\vec r^*$，$t=(L/U)t^*$，$p=\rho U^2 p^*$。代入 N-S：
$$
\rho\frac{U^2}{L}\frac{D\vec v^*}{Dt^*}
=-\frac{\rho U^2}{L}\nabla^*p^*
+\frac{\mu U}{L^2}\nabla^{*2}\vec v^*
$$
除以 $\rho U^2/L$：
$$
\frac{D\vec v^*}{Dt^*}=-\nabla^*p^*+\frac{1}{Re}\nabla^{*2}\vec v^*
$$
其中
$$
Re=\frac{\rho UL}{\mu}=\frac{UL}{\nu}
$$
类似可得马赫数 $Ma=U/c$、弗劳德数 $Fr=U/\sqrt{gL}$、斯特劳哈尔数 $St=fL/U$、普朗特数 $Pr=\nu/\alpha$。

物理意义：相似流动具有相同无量纲数。$Re$ 表示惯性/黏性比。特例：$Re\ll1$ 蠕动流，$Re\gg1$ 边界层/湍流。

应用：风洞缩比、船舶模型、血流模拟。

**知识点 15.6 量纲分析与白金汉 π 定理**

定义：若物理过程涉及 $n$ 个变量、$k$ 个基本量纲，则存在 $n-k$ 个独立无量纲数。

推导：设变量 $q_1,\dots,q_n$ 满足 $f(q_1,\dots,q_n)=0$，量纲矩阵秩为 $k$。取 $k$ 个重复变量，构造幂乘积使量纲抵消：
$$
\pi_j=q_{j+k}\prod_{i=1}^{k}q_i^{a_{ij}}
$$
解线性方程组得 $a_{ij}$，从而 $F(\pi_1,\dots,\pi_{n-k})=0$。

物理意义：减少实验变量数目，指导相似试验。特例：单摆周期 $T=\sqrt{L/g}\,f(\theta_0)$。由 $n=4$（$T,L,g,\theta_0$）、$k=2$（长度、时间）得两个 $\pi$：$T\sqrt{g/L}$ 与 $\theta_0$。

应用：阻力-升力相似、传热传质、化学工程。

### 第16章 边界层、湍流与相似理论

**知识点 16.1 普朗特边界层方程**

定义：高 $Re$ 流动中黏性效应集中在物面附近薄层，厚度 $\delta\ll L$。

推导：对二维定常不可压 N-S 做量级分析，$u\sim U$，$x\sim L$，$y\sim\delta$，由连续方程
$$
\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0
\;\Rightarrow\;
v\sim U\frac{\delta}{L}
$$
$y$ 动量方程各项量级
$$
u\frac{\partial v}{\partial x}\sim U^2\frac{\delta}{L^2},\qquad
v\frac{\partial v}{\partial y}\sim U^2\frac{\delta}{L^2},\qquad
\nu\frac{\partial^2 v}{\partial y^2}\sim \frac{\nu U\delta}{L\delta^2}
$$
对流项远大于黏性项（当 $Re_L\gg1$），故 $\partial p/\partial y\approx0$，压力由外部势流给定。$x$ 方程保留黏性项（$u_{,yy}\sim U/\delta^2$ 主导），得
$$
u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}
=U\frac{dU}{dx}+\nu\frac{\partial^2 u}{\partial y^2}
$$
$$
\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0
$$
边界条件 $y=0:u=v=0$；$y\to\infty:u\to U(x)$。

物理意义：边界层把黏性效应压缩到薄层，外部用无黏流。特例：平板 $dU/dx=0$，得布拉修斯方程。

应用：机翼摩阻、热边界层、分离预测。

**知识点 16.2 布拉修斯相似解**

定义：平板边界层速度剖面自相似。

推导：引入流函数 $\psi$，$u=\partial\psi/\partial y,\ v=-\partial\psi/\partial x$。令
$$
\eta=y\sqrt{\frac{U}{\nu x}},\qquad
\psi=\sqrt{\nu U x}\,f(\eta)
$$
则
$$
u=U f'(\eta),\qquad
v=\frac{1}{2}\sqrt{\frac{\nu U}{x}}(\eta f'-f)
$$
代入边界层方程，各项化为
$$
u\frac{\partial u}{\partial x}=-\frac{U^2}{2x}\eta f'f''
$$
$$
v\frac{\partial u}{\partial y}=\frac{U^2}{2x}(f f''-\eta f'f'')
$$
$$
\nu\frac{\partial^2 u}{\partial y^2}=\frac{U^2}{2x}f'''
$$
相加得布拉修斯方程
$$
f'''+\tfrac12 f f''=0
$$
边界条件
$$
f(0)=f'(0)=0,\qquad f'(\infty)=1
$$
数值解给出壁面剪切
$$
\tau_w=\mu\left.\frac{\partial u}{\partial y}\right|_{y=0}
=0.332\rho U^2\sqrt{\frac{\nu}{Ux}}
\;\Rightarrow\;
C_f=\frac{0.664}{\sqrt{Re_x}}
$$
边界层厚度 $\delta_{99}\approx5.0\sqrt{\nu x/U}$。

物理意义：相似解把偏微分方程降为常微分方程。特例：$Re_x$ 增大，$\delta$ 增大而 $C_f$ 减小。

应用：平板摩阻、传热类比的出发点。

**知识点 16.3 冯·卡门动量积分方程**

定义：对边界层做动量守恒积分，得工程近似法。

推导：对边界层方程沿 $y$ 从 $0$ 到 $\delta$ 积分。连续方程给
$$
v=-\int_0^y\frac{\partial u}{\partial x}dy
$$
动量方程积分：
$$
\int_0^\delta u\frac{\partial u}{\partial x}dy+\int_0^\delta v\frac{\partial u}{\partial y}dy
=U\frac{dU}{dx}\delta+\nu\left[\frac{\partial u}{\partial y}\right]_0^\delta
$$
利用分部积分与边界条件，化为
$$
\frac{d}{dx}\int_0^\delta u(U-u)\,dy+\frac{dU}{dx}\int_0^\delta(U-u)\,dy
=\frac{\tau_w}{\rho}
$$
定义排移厚度与动量损失厚度：
$$
\delta^*=\int_0^\delta\left(1-\frac{u}{U}\right)dy,\qquad
\theta=\int_0^\delta\frac{u}{U}\left(1-\frac{u}{U}\right)dy
$$
则
$$
\frac{d}{dx}(U^2\theta)+U\frac{dU}{dx}\delta^*=\frac{\tau_w}{\rho}
$$
平板 $U$ 常数时
$$
\frac{d\theta}{dx}=\frac{C_f}{2}
$$

物理意义：避开细节，用壁面摩擦与剖面参数关系求解。特例：取 $u/U=2\eta-\eta^2$，
$$
\theta=\frac{2}{15}\delta,\quad
\frac{\tau_w}{\rho}=\frac{2\nu U}{\delta}
\;\Rightarrow\;
\delta=5.48\sqrt{\frac{\nu x}{U}},\quad C_f=\frac{0.730}{\sqrt{Re_x}}
$$

应用：湍流边界层、分离判据、工程摩阻估算。

**知识点 16.4 边界层分离与压差阻力**

定义：逆压梯度 $\dfrac{dp}{dx}>0$ 使近壁流体减速，最终回流，边界层脱离物面称分离。

推导：壁面处动量方程
$$
\mu\left.\frac{\partial^2 u}{\partial y^2}\right|_{y=0}=\frac{dp}{dx}
$$
顺压梯度 $\dfrac{dp}{dx}<0$ 时壁面速度剖面凸向上，稳定；逆压梯度时壁面处 $u_{yy}>0$，剖面出现拐点，随后壁面附近 $u<0$，分离点由
$$
\left.\frac{\partial u}{\partial y}\right|_{y=0}=0
$$
确定。

物理意义：分离导致尾流与压差阻力，升力下降（失速）。特例：圆柱绕流 $Re\approx2\times10^5$ 边界层转捩使分离点后移，阻力骤降（阻力危机）。

应用：机翼失速、扩压管、汽车减阻。

**知识点 16.5 湍流与雷诺平均方程**

定义：湍流是随机、多尺度、有旋的流动。瞬时量分解为均值与脉动：
$$
u_i=\bar u_i+u_i',\qquad \bar{u_i'}=0
$$

推导：将分解代入 N-S 并对时间平均。不可压连续性
$$
\partial_i\bar u_i=0,\qquad \partial_i u_i'=0
$$
动量方程平均，非线性对流项
$$
\overline{u_j\partial_j u_i}
=\bar u_j\partial_j\bar u_i
+\overline{u_j'\partial_j u_i'}
$$
把脉动项移到右边，得到雷诺平均 N-S 方程（RANS）：
$$
\rho\left(\frac{\partial\bar u_i}{\partial t}+\bar u_j\partial_j\bar u_i\right)
=-\partial_i\bar p+\mu\partial_j^2\bar u_i
+\partial_j\left(-\rho\overline{u_i'u_j'}\right)+\rho f_i
$$
其中雷诺应力
$$
\tau_{ij}^{t}=-\rho\overline{u_i'u_j'}
$$
需要湍流模型（$k$-$\varepsilon$、$k$-$\omega$、SST）封闭。

物理意义：湍流脉动对平均流施加减应力，增强动量交换。特例：均匀剪切湍流有 $\overline{u'v'}\neq0$。

应用：CFD 工程湍流模拟、大气边界层、管道流。

**知识点 16.6 雷诺数与湍流转捩**

定义：层流向湍流转变称转捩，判据为临界雷诺数。

推导：对扰动做线性稳定性分析，扰动方程
$$
\partial_t u'+\bar U\partial_x u'=-(u'\cdot\nabla)\bar U-\nabla p'+\frac{1}{Re}\nabla^2u'
$$
设模态 $u'=\hat u(y)e^{i(kx-\omega t)}$，得奥尔-索末菲方程。当 $Re>Re_{cr}$ 时存在增长模态，流动失稳。圆管 $Re_{cr}\approx2300$，平板边界层 $Re_{x,cr}\approx5\times10^5$。

物理意义：惯性超过黏性耗散时扰动放大。特例：低 $Re$ 层流稳定，高 $Re$ 湍流充分发展。

应用：管流设计、边界层抽吸、减阻。

**知识点 16.7 升力、阻力与相似律**

定义：物体受流体合力分解为升力（垂直来流）与阻力（沿来流）。

推导：由动量定理与伯努利，机翼环量 $\Gamma$ 产生升力（库塔-儒可夫斯基定理）
$$
L'=\rho U_\infty\Gamma
$$
推导思路：绕翼型的环量使上表面速度高、压力低，积分压力差得升力；或由涡的动量反作用。无量纲化
$$
C_L=\frac{L}{\tfrac12\rho U^2S},\qquad
C_D=\frac{D}{\tfrac12\rho U^2S}
$$
升力线斜率（薄翼理论）
$$
C_L=2\pi\alpha
$$
对于有限翼展，诱导阻力
$$
C_{D,i}=\frac{C_L^2}{\pi AR\,e},\qquad AR=\frac{b^2}{S}
$$

物理意义：环量是升力来源，诱导阻力源于翼尖涡。特例：对称翼型零迎角无升力；$Re$ 与 $Ma$ 相似保证风洞数据可用。

应用：飞行器气动设计、风力机、汽车气动。

---
⬆ 上级领域：[[机械与工程基础]]
