---
课程名称: 动力学与控制核心贯通课II
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

# 动力学与控制核心贯通课II


### 第1章 质点动力学与牛顿运动定律

**知识点1.1　参考系、质点模型与运动学描述**

定义：质点是忽略几何尺度、只保留质量的理想点；参考系是用于度量位置与时间的坐标系与钟的组合；惯性参考系是使惯性定律成立的参考系。

假设：在所研究问题中物体尺度远小于运动特征尺度，且转动与变形对平动无本质影响时可用质点模型；时间与空间绝对，即 $\Delta t$ 与参考系无关。

位置矢量与速度、加速度：
$$\mathbf r=\mathbf r(t),\qquad \mathbf v=\frac{\mathrm d\mathbf r}{\mathrm dt},\qquad \mathbf a=\frac{\mathrm d\mathbf v}{\mathrm dt}=\frac{\mathrm d^2\mathbf r}{\mathrm dt^2}.$$

推导自然坐标系分解。设弧长 $s$，切向单位矢量 $\boldsymbol\tau=\mathrm d\mathbf r/\mathrm ds$，则
$$\mathbf v=\frac{\mathrm d\mathbf r}{\mathrm dt}=\frac{\mathrm d\mathbf r}{\mathrm ds}\frac{\mathrm ds}{\mathrm dt}=v\boldsymbol\tau.$$
对时间求导：
$$\mathbf a=\frac{\mathrm dv}{\mathrm dt}\boldsymbol\tau+v\frac{\mathrm d\boldsymbol\tau}{\mathrm dt}.$$
因 $\boldsymbol\tau\cdot\boldsymbol\tau=1\Rightarrow \boldsymbol\tau\cdot\mathrm d\boldsymbol\tau/\mathrm dt=0$，故 $\mathrm d\boldsymbol\tau/\mathrm dt$ 垂直于 $\boldsymbol\tau$，方向指向曲率中心。设曲率半径 $\rho$，$\mathrm d\boldsymbol\tau/\mathrm ds=\mathbf n/\rho$，于是
$$\mathbf a=\dot v\,\boldsymbol\tau+\frac{v^2}{\rho}\mathbf n.$$

物理意义：切向分量改变速率快慢，法向分量 $v^2/\rho$ 改变速度方向。极限特例：直线运动 $\rho\to\infty$ 时 $a_n=0$；匀速圆周 $\dot v=0$ 时只剩向心加速度。

典型应用：飞行器航迹的过载分析，法向加速度决定飞行员承受的过载 $n=v^2/(\rho g)$。

**知识点1.2　牛顿三定律与惯性系**

定义：牛顿第一定律（惯性定律）指出自由质点保持静止或匀速直线运动；第二定律给出力与动量变化率的关系；第三定律指出作用力与反作用力等大反向共线。

假设：质量为常量、速度远小于光速（非相对论）、惯性系内成立。

第二定律的原始形式与常用形式：
$$\mathbf F=\frac{\mathrm d\mathbf p}{\mathrm dt},\qquad \mathbf p=m\mathbf v,$$
当 $m$ 为常量时 $\mathbf F=m\mathbf a$。

推导动量形式对质量变化情形的推广（变质量系统）。设 $t$ 时刻系统质量为 $m$、速度为 $\mathbf v$；$t+\mathrm dt$ 时刻主体质量 $m+\mathrm dm$、速度 $\mathbf v+\mathrm d\mathbf v$，同时有质量 $-\mathrm dm>0$ 被抛出、相对主体速度为 $\mathbf u$（绝对速度为 $\mathbf v+\mathbf u$）。系统总动量守恒（无外力）：
$$m\mathbf v=(m+\mathrm dm)(\mathbf v+\mathrm d\mathbf v)+(-\mathrm dm)(\mathbf v+\mathbf u).$$
展开并忽略二阶小量 $\mathrm dm\,\mathrm d\mathbf v$：
$$m\mathbf v=m\mathbf v+m\,\mathrm d\mathbf v+\mathbf v\,\mathrm dm-\mathbf v\,\mathrm dm-\mathbf u\,\mathrm dm,$$
$$\Rightarrow m\,\mathrm d\mathbf v=-\mathbf u\,\mathrm dm\ \Rightarrow\ m\frac{\mathrm d\mathbf v}{\mathrm dt}=-\mathbf u\frac{\mathrm dm}{\mathrm dt}.$$
有外力时 $\mathbf F_{\text{外}}+m\dot{\mathbf v}=-\mathbf u\dot m$，即
$$m\frac{\mathrm d\mathbf v}{\mathrm dt}=\mathbf F_{\text{外}}+\mathbf F_{\text{推力}},\qquad \mathbf F_{\text{推力}}=-\mathbf u\frac{\mathrm dm}{\mathrm dt}.$$

物理意义：反作用力体现动量流的输运。极限特例：$\dot m=0$ 退化为常规第二定律。

典型应用：火箭方程（下一章推导）。

**知识点1.3　质点运动微分方程的积分与三类问题**

定义：质点运动微分方程在直角坐标为
$$m\ddot x=F_x,\quad m\ddot y=F_y,\quad m\ddot z=F_z.$$

假设：力已知（或为位置、速度、时间的函数）。

第一类问题（已知运动求力）：由 $\mathbf r(t)$ 两次求导得 $\mathbf a$，代入即可。

第二类问题（已知力求运动）：解二阶常微分方程。一维情形 $m\ddot x=F(x)$，用变换
$$\ddot x=\frac{\mathrm dv}{\mathrm dt}=v\frac{\mathrm dv}{\mathrm dx},$$
推导能量积分：
$$mv\frac{\mathrm dv}{\mathrm dx}=F(x)\ \Rightarrow\ mv\,\mathrm dv=F(x)\,\mathrm dx,$$
$$\frac12 mv^2-\frac12 mv_0^2=\int_{x_0}^{x}F(x')\,\mathrm dx'.$$
即动能定理的一维形式。

推导与时间无关的力 $F(x)$ 的运动解：
$$\int \frac{\mathrm dv}{F(x)/m}=\int \mathrm dt\ \Rightarrow\ \int_{v_0}^{v}\frac{\mathrm dv'}{F(x)/m}=t-t_0.$$
再由 $v=\mathrm dx/\mathrm dt$ 分离变量积分得 $x(t)$。

物理意义：一次积分降低一阶，力只依赖位置时可得到能量积分。极限特例：$F=$ 常量给出 $x=x_0+v_0t+\frac12 at^2$。

典型应用：弹簧振子 $m\ddot x=-kx$ 的解为简谐运动 $x=A\cos(\omega t+\varphi)$，$\omega=\sqrt{k/m}$。

**知识点1.4　动量定理与冲量定理**

定义：动量 $\mathbf p=m\mathbf v$；冲量 $\mathbf I=\int_{t_1}^{t_2}\mathbf F\,\mathrm dt$。

假设：合力可在时间区间内积分。

推导动量定理：由第二定律
$$\mathbf F=\frac{\mathrm d\mathbf p}{\mathrm dt}\ \Rightarrow\ \int_{t_1}^{t_2}\mathbf F\,\mathrm dt=\int_{t_1}^{t_2}\mathrm d\mathbf p=\mathbf p_2-\mathbf p_1.$$
即 $\mathbf I=\Delta\mathbf p$。

物理意义：力对时间的累积等于动量增量；碰撞等短时强相互作用中，平均力 $\bar{\mathbf F}=\Delta\mathbf p/\Delta t$。

极限特例：弹性碰撞中动量守恒；非弹性完全碰撞中动能损失但动量仍守恒。

典型应用：缓冲公式 $F_{\max}\approx\Delta p/\Delta t$，安全带延长碰撞时间以减小峰值力。

**知识点1.5　角动量定理与面积速度**

定义：对定点 $O$ 的角动量 $\mathbf L=\mathbf r\times\mathbf p$；力矩 $\mathbf M=\mathbf r\times\mathbf F$。

假设：定点为惯性系中固定点。

推导：对 $\mathbf L$ 求导
$$\frac{\mathrm d\mathbf L}{\mathrm dt}=\frac{\mathrm d\mathbf r}{\mathrm dt}\times\mathbf p+\mathbf r\times\frac{\mathrm d\mathbf p}{\mathrm dt}.$$
因 $\dot{\mathbf r}\times\mathbf p=\mathbf v\times m\mathbf v=\mathbf 0$，又 $\dot{\mathbf p}=\mathbf F$，故
$$\frac{\mathrm d\mathbf L}{\mathrm dt}=\mathbf r\times\mathbf F=\mathbf M.$$
中心力场中 $\mathbf M=\mathbf r\times F(r)\hat{\mathbf r}=\mathbf 0$，故 $\mathbf L=$ 常量，运动为平面运动。

推导面积速度定理（开普勒第二定律）：
$$\mathrm d\mathbf A=\frac12|\mathbf r\times\mathrm d\mathbf r|\ \Rightarrow\ \frac{\mathrm d\mathbf A}{\mathrm dt}=\frac12|\mathbf r\times\mathbf v|=\frac{L}{2m}=\text{常量}.$$

物理意义：角动量守恒等价于面积速度恒定；轨道近日点速度快、远日点速度慢。

极限特例：直线运动过原点时 $\mathbf L=0$。

典型应用：行星运动、卫星轨道面不变。

**知识点1.6　动能定理、功与功率**

定义：功 $W=\int\mathbf F\cdot\mathrm d\mathbf r$；动能 $T=\frac12 mv^2$；功率 $P=\mathbf F\cdot\mathbf v$。

推导动能定理：
$$m\frac{\mathrm d\mathbf v}{\mathrm dt}=\mathbf F\ \Rightarrow\ m\mathbf v\cdot\frac{\mathrm d\mathbf v}{\mathrm dt}=\mathbf F\cdot\mathbf v.$$
因 $\frac{\mathrm d}{\mathrm dt}\left(\frac12 mv^2\right)=m\mathbf v\cdot\dot{\mathbf v}$，对时间积分：
$$\int_{t_1}^{t_2}\frac{\mathrm d}{\mathrm dt}\left(\frac12 mv^2\right)\mathrm dt=\int_{t_1}^{t_2}\mathbf F\cdot\mathbf v\,\mathrm dt,$$
$$T_2-T_1=\int_{t_1}^{t_2}\mathbf F\cdot\frac{\mathrm d\mathbf r}{\mathrm dt}\mathrm dt=\int_{\mathbf r_1}^{\mathbf r_2}\mathbf F\cdot\mathrm d\mathbf r=W.$$

物理意义：合外力做功等于动能增量，做功只依赖路径积分，摩擦力做功与路径有关。

极限特例：匀速圆周运动合外力不做功（$\mathbf F\perp\mathbf v$），动能不变。

典型应用：能量法求解复杂约束问题，例如滑块-弹簧系统。

**知识点1.7　保守力与势能**

定义：若力场满足 $\oint\mathbf F\cdot\mathrm d\mathbf r=0$（等价于路径无关），则称保守力，可定义势能 $V$ 使
$$\mathbf F=-\nabla V.$$

推导保守力与旋度关系（斯托克斯定理）：
$$\oint_{\partial S}\mathbf F\cdot\mathrm d\mathbf r=\iint_S(\nabla\times\mathbf F)\cdot\mathrm d\mathbf S.$$
路径无关 $\Leftrightarrow\oint=0$ 对任意回路 $\Leftrightarrow \nabla\times\mathbf F=\mathbf 0$。

推导一元情形：$F(x)=-\mathrm dV/\mathrm dx$。由动能定理
$$T_2-T_1=\int_{x_1}^{x_2}F\,\mathrm dx=-\int_{x_1}^{x_2}\frac{\mathrm dV}{\mathrm dx}\mathrm dx=-(V_2-V_1),$$
$$\Rightarrow T_1+V_1=T_2+V_2=E=\text{常量}.$$

物理意义：机械能守恒源于力的管量场性质；势能零点可任意选取。

极限特例：重力 $V=mgz$；弹簧 $V=\frac12 kx^2$；万有引力 $V=-GMm/r$。

典型应用：有心力场中的有效势 $V_{\text{eff}}=V+L^2/(2mr^2)$。

**知识点1.8　约束与广义坐标初步**

定义：约束是限制质点位置或速度的几何/运动学条件。完整约束 $f(\mathbf r_i,t)=0$；非完整约束不可积分为位置关系。

假设：约束方程可微。

推导：$N$ 个质点、$s$ 个完整约束时自由度 $n=3N-s$。引入广义坐标 $q_1,\dots,q_n$：
$$\mathbf r_i=\mathbf r_i(q_1,\dots,q_n,t).$$
速度关系
$$\mathbf v_i=\sum_{j=1}^{n}\frac{\partial\mathbf r_i}{\partial q_j}\dot q_j+\frac{\partial\mathbf r_i}{\partial t}.$$
虚位移（$t$ 冻结）
$$\delta\mathbf r_i=\sum_{j=1}^{n}\frac{\partial\mathbf r_i}{\partial q_j}\delta q_j.$$

物理意义：广义坐标自动满足完整约束，将问题维数从 $3N$ 降为 $n$。

极限特例：无约束时广义坐标即直角坐标。

典型应用：单摆用 $\theta$ 一个广义坐标描述；双摆用 $\theta_1,\theta_2$。

### 第2章 质点系动力学与守恒定律

**知识点2.1　质心运动定理**

定义：质点系质心
$$\mathbf r_C=\frac{1}{m}\sum_i m_i\mathbf r_i,\qquad m=\sum_i m_i.$$

推导：对每个质点 $m_i\ddot{\mathbf r}_i=\mathbf F_i^{\text{外}}+\sum_{j\neq i}\mathbf F_{ij}$，求和：
$$\sum_i m_i\ddot{\mathbf r}_i=\sum_i\mathbf F_i^{\text{外}}+\sum_i\sum_{j\neq i}\mathbf F_{ij}.$$
内力成对出现，由牛顿第三定律 $\mathbf F_{ij}=-\mathbf F_{ji}$ 相互抵消，故
$$m\ddot{\mathbf r}_C=\sum_i\mathbf F_i^{\text{外}}.$$

物理意义：质心像集中了全部质量、受全部外力的质点一样运动，内力不影响质心运动。

极限特例：合外力为零时质心匀速直线运动。

典型应用：爆炸碎片、火箭抛壳的质心轨迹。

**知识点2.2　动量守恒定律**

定义：系统总动量 $\mathbf P=\sum_i m_i\mathbf v_i=m\mathbf v_C$。

推导：由质心运动定理
$$\frac{\mathrm d\mathbf P}{\mathrm dt}=\mathbf F^{\text{外}}.$$
若 $\mathbf F^{\text{外}}=\mathbf 0$，则 $\mathbf P=$ 常量，即
$$\sum_i m_i\mathbf v_i=\sum_i m_i\mathbf v_i'.$$

物理意义：守恒是分量式的，某一方向受力为零则该方向动量守恒。

极限特例：完全非弹性碰撞后共速 $\mathbf v'=\frac{\sum m_i\mathbf v_i}{\sum m_i}$。

典型应用：反冲运动、两体碰撞。

**知识点2.3　角动量定理与角动量守恒**

定义：对定点 $O$ 的系统角动量 $\mathbf L_O=\sum_i\mathbf r_i\times m_i\mathbf v_i$；总外力矩 $\mathbf M_O^{\text{外}}=\sum_i\mathbf r_i\times\mathbf F_i^{\text{外}}$。

推导：
$$\frac{\mathrm d\mathbf L_O}{\mathrm dt}=\sum_i\dot{\mathbf r}_i\times m_i\mathbf v_i+\sum_i\mathbf r_i\times m_i\dot{\mathbf v}_i.$$
第一项 $\sum_i \mathbf v_i\times m_i\mathbf v_i=0$。第二项
$$\sum_i\mathbf r_i\times\mathbf F_i=\sum_i\mathbf r_i\times\mathbf F_i^{\text{外}}+\sum_i\sum_{j\neq i}\mathbf r_i\times\mathbf F_{ij}.$$
内力矩之和
$$\mathbf r_i\times\mathbf F_{ij}+\mathbf r_j\times\mathbf F_{ji}=(\mathbf r_i-\mathbf r_j)\times\mathbf F_{ij}=\mathbf 0$$
（内力沿连线时；一般中心力成立），故
$$\frac{\mathrm d\mathbf L_O}{\mathrm dt}=\mathbf M_O^{\text{外}}.$$
当 $\mathbf M_O^{\text{外}}=\mathbf 0$ 时 $\mathbf L_O$ 守恒。

物理意义：转动惯性；对质心亦成立 $\dot{\mathbf L}_C=\mathbf M_C^{\text{外}}$。

极限特例：中心力问题角动量守恒即开普勒第二定律。

典型应用：花样滑冰收臂加速、直升机尾桨。

**知识点2.4　机械能守恒与内力功**

定义：系统动能 $T=\sum_i\frac12 m_iv_i^2$，势能 $V$，机械能 $E=T+V$。

推导：对系统动能定理
$$\mathrm dT=\sum_i\mathbf F_i^{\text{外}}\cdot\mathrm d\mathbf r_i+\sum_i\sum_{j\neq i}\mathbf F_{ij}\cdot\mathrm d\mathbf r_i.$$
当内力保守且外力保守时，总功等于势能减少 $-\mathrm dV$，故 $\mathrm d(T+V)=0$。注意内力可以做功（如爆炸、摩擦生热），但保守内力功已并入势能。

由质心系分解（柯尼希定理）：
$$T=\frac12 m v_C^2+T_C,\qquad T_C=\sum_i\frac12 m_i v_i'^{\,2}.$$
推导：$\mathbf v_i=\mathbf v_C+\mathbf v_i'$，
$$T=\sum_i\frac12 m_i(\mathbf v_C+\mathbf v_i')^2=\frac12 m v_C^2+\mathbf v_C\cdot\frac{\mathrm d}{\mathrm dt}\sum_i m_i\mathbf r_i'+\sum_i\frac12 m_iv_i'^{\,2}.$$
质心系中 $\sum_i m_i\mathbf r_i'=\mathbf 0$，故交叉项为零，得证。

物理意义：机械能 = 质心平动动能 + 内动能。

极限特例：刚性碰撞中 $T_C$ 转换为热。

典型应用：刚体平面运动 $T=\frac12 mv_C^2+\frac12 I_C\omega^2$。

**知识点2.5　二体问题与约化质量**

定义：两体在中心力 $\mathbf F_{12}=-\nabla_1 V(|\mathbf r_1-\mathbf r_2|)$ 下运动。

推导：$\mathbf r=\mathbf r_1-\mathbf r_2$，两方程
$$m_1\ddot{\mathbf r}_1=\mathbf F,\qquad m_2\ddot{\mathbf r}_2=-\mathbf F.$$
相除消元：$\ddot{\mathbf r}_1-\ddot{\mathbf r}_2=\left(\frac1{m_1}+\frac1{m_2}\right)\mathbf F$，令
$$\mu=\frac{m_1m_2}{m_1+m_2},\qquad \mu\ddot{\mathbf r}=\mathbf F(\mathbf r).$$
质心 $\mathbf R=(m_1\mathbf r_1+m_2\mathbf r_2)/(m_1+m_2)$ 满足 $\ddot{\mathbf R}=\mathbf 0$。

物理意义：相对运动等价于质量为 $\mu$ 的单粒子在中心势中运动。

极限特例：$m_2\gg m_1$ 时 $\mu\approx m_1$，退化为固定中心问题。

典型应用：开普勒轨道能量 $E=\frac12\mu v^2-GMm\mu/r$。

**知识点2.6　碰撞与恢复系数**

定义：恢复系数 $e=\dfrac{v_2'-v_1'}{v_1-v_2}$（一维，分离速度与接近速度之比）。

推导：动量守恒
$$m_1v_1+m_2v_2=m_1v_1'+m_2v_2'.$$
结合恢复系数定义 $v_2'-v_1'=e(v_1-v_2)$，解得
$$v_1'=\frac{m_1v_1+m_2v_2-m_2e(v_1-v_2)}{m_1+m_2},\quad v_2'=\frac{m_1v_1+m_2v_2+m_1e(v_1-v_2)}{m_1+m_2}.$$
相对动能损失
$$\Delta T=\frac12\mu(1-e^2)(v_1-v_2)^2,\quad \mu=\frac{m_1m_2}{m_1+m_2}.$$

物理意义：$e=1$ 弹性、$e=0$ 完全非弹性。

极限特例：等质量弹性碰撞交换速度。

典型应用：碰撞摆测弹速、航天器引力助推。

**知识点2.7　火箭方程（齐奥尔科夫斯基公式）**

定义：无外力、火箭以相对速率 $u$ 喷气。

推导：由知识点1.2，
$$m\frac{\mathrm dv}{\mathrm dt}=-u\frac{\mathrm dm}{\mathrm dt}\ \Rightarrow\ \mathrm dv=-u\frac{\mathrm dm}{m}.$$
积分：
$$\int_0^t\mathrm dv=-u\int_{m_0}^{m}\frac{\mathrm dm}{m}\ \Rightarrow\ v-v_0=u\ln\frac{m_0}{m}.$$
$$\Delta v=u\ln\frac{m_0}{m_f}=u\ln\frac{m_0}{m_0-m_p}.$$

物理意义：$\Delta v$ 由比冲与质量比的对数决定；引力场中需加损失项。

极限特例：$m_p\to0$ 时 $\Delta v\to0$。

典型应用：多级火箭、霍曼转移总 $\Delta v$ 预算。

**知识点2.8　有心力场与有效势**

定义：有心力 $\mathbf F=F(r)\hat{\mathbf r}$，势能 $V(r)$。

推导：角动量守恒 $L=mr^2\dot\theta=$ 常量。能量
$$E=\frac12 m(\dot r^2+r^2\dot\theta^2)+V(r).$$
用 $\dot\theta=L/(mr^2)$ 消去：
$$E=\frac12 m\dot r^2+\frac{L^2}{2mr^2}+V(r)=\frac12 m\dot r^2+V_{\text{eff}}(r),$$
$$V_{\text{eff}}(r)=V(r)+\frac{L^2}{2mr^2}.$$
径向方程 $m\ddot r=-V_{\text{eff}}'(r)$。

物理意义：离心项作为排斥势，圆轨道条件 $V_{\text{eff}}'(r_0)=0$ 给出 $F(r_0)=-L^2/(mr_0^3)$。

极限特例：$L=0$ 时 $V_{\text{eff}}=V$，径向自由落体。

典型应用：轨道稳定性判据 $V_{\text{eff}}''(r_0)>0$。

### 第3章 拉格朗日力学

**知识点3.1　约束分类与自由度**

定义：完整约束可由 $f_\alpha(\mathbf r_i,t)=0$ 表示；非完整约束形如 $\sum_i\mathbf A_i\cdot\mathrm d\mathbf r_i+B\,\mathrm dt=0$ 且不可积分。

假设：约束力不做虚功（理想约束）。

推导自由度：$N$ 个质点有 $3N$ 个坐标，$s$ 个独立完整约束使独立坐标减为 $n=3N-s$，可取 $q_1,\dots,q_n$。

物理意义：理想约束包括光滑面、刚性杆、滚动无滑动等。

极限特例：非完整约束（冰刀、滚动圆盘）不能用 $\delta q$ 全独立处理。

典型应用：滚动圆盘 5 个坐标 1 个约束、1 个非完整约束。

**知识点3.2　虚功原理**

定义：虚位移是满足约束且时间冻结的无穷小位移；理想约束下约束力虚功为零。

推导：静力学平衡 $\mathbf F_i+\mathbf N_i=\mathbf 0$，乘以 $\delta\mathbf r_i$ 求和：
$$\sum_i\mathbf F_i\cdot\delta\mathbf r_i+\sum_i\mathbf N_i\cdot\delta\mathbf r_i=0.$$
由理想约束第二项为零，得虚功原理
$$\delta W=\sum_i\mathbf F_i\cdot\delta\mathbf r_i=0.$$

物理意义：平衡的充要条件（对理想约束）是主动力虚功为零，是全部静力学的基础。

极限特例：单约束时退化为力平衡。

典型应用：杠杆、滑轮系统平衡分析。

**知识点3.3　达朗贝尔原理与拉格朗日方程推导**

定义：达朗贝尔惯性力 $-\dot{\mathbf p}_i$。

推导：动力学方程 $\mathbf F_i-\dot{\mathbf p}_i+\mathbf N_i=0$，取虚位移并利用理想约束：
$$\sum_i(\mathbf F_i-\dot{\mathbf p}_i)\cdot\delta\mathbf r_i=0.$$
以 $q_j$ 表示，$\delta\mathbf r_i=\sum_j\frac{\partial\mathbf r_i}{\partial q_j}\delta q_j$，$\delta q_j$ 独立，故
$$\sum_i(\mathbf F_i-\dot{\mathbf p}_i)\cdot\frac{\partial\mathbf r_i}{\partial q_j}=0,\quad j=1,\dots,n.$$
定义广义力 $Q_j=\sum_i\mathbf F_i\cdot\partial\mathbf r_i/\partial q_j$。处理惯性项：
$$\sum_i\dot{\mathbf p}_i\cdot\frac{\partial\mathbf r_i}{\partial q_j}=\sum_i m_i\ddot{\mathbf r}_i\cdot\frac{\partial\mathbf r_i}{\partial q_j}.$$
利用两恒等式
$$\frac{\partial\mathbf r_i}{\partial q_j}=\frac{\partial\dot{\mathbf r}_i}{\partial\dot q_j},\qquad \frac{\mathrm d}{\mathrm dt}\frac{\partial\mathbf r_i}{\partial q_j}=\frac{\partial\dot{\mathbf r}_i}{\partial q_j},$$
得
$$\sum_i m_i\ddot{\mathbf r}_i\cdot\frac{\partial\mathbf r_i}{\partial q_j}=\frac{\mathrm d}{\mathrm dt}\frac{\partial T}{\partial\dot q_j}-\frac{\partial T}{\partial q_j},$$
其中 $T=\sum_i\frac12 m_i\dot{\mathbf r}_i^2$。于是
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial T}{\partial\dot q_j}-\frac{\partial T}{\partial q_j}=Q_j.$$

物理意义：动力学方程从矢量变为 $n$ 个标量方程，约束力自动消去。

极限特例：$Q_j=0$ 时为自由运动。

典型应用：机械臂、机构动力学建模。

**知识点3.4　拉格朗日函数与最小作用量形式**

定义：对有势力 $Q_j=-\partial V/\partial q_j$，定义 $L=T-V$。

推导：由 $\partial V/\partial\dot q_j=0$，
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial(T-V)}{\partial\dot q_j}-\frac{\partial(T-V)}{\partial q_j}=0,$$
即欧拉-拉格朗日方程
$$\boxed{\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}=0.}$$

物理意义：$L$ 是标量、与坐标选择无关；方程在广义坐标变换下形式不变（协变）。

极限特例：$L=T$ 且无势时为自由粒子 $q_j$ 匀速。

典型应用：$L=\frac12 m\dot x^2-\frac12 kx^2$ 给出简谐振子。

**知识点3.5　广义势与洛伦兹力**

定义：广义势 $U(q,\dot q,t)$ 使 $Q_j=\frac{\mathrm d}{\mathrm dt}\frac{\partial U}{\partial\dot q_j}-\frac{\partial U}{\partial q_j}$。

推导：拉格朗日方程写为 $L=T-U$。对电磁场
$$U=q_e\phi-q_e\mathbf A\cdot\mathbf v.$$
验证洛伦兹力：
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial U}{\partial\dot x}=\frac{\mathrm d}{\mathrm dt}(-q_eA_x)=-q_e\left(\frac{\partial A_x}{\partial t}+\dot x\frac{\partial A_x}{\partial x}+\dot y\frac{\partial A_x}{\partial y}+\dot z\frac{\partial A_x}{\partial z}\right),$$
$$-\frac{\partial U}{\partial x}=-q_e\frac{\partial\phi}{\partial x}+q_e\left(\dot x\frac{\partial A_x}{\partial x}+\dot y\frac{\partial A_y}{\partial x}+\dot z\frac{\partial A_z}{\partial x}\right).$$
相加得
$$Q_x=q_e\left(-\frac{\partial\phi}{\partial x}-\frac{\partial A_x}{\partial t}+\dot y\left(\frac{\partial A_y}{\partial x}-\frac{\partial A_x}{\partial y}\right)+\dot z\left(\frac{\partial A_z}{\partial x}-\frac{\partial A_x}{\partial z}\right)\right)$$
$$=q_e\left(E_x+(\mathbf v\times\mathbf B)_x\right).$$

物理意义：速度相关力也可纳入拉格朗日框架。

极限特例：无磁场时退化为电场力。

典型应用：带电粒子在电磁场中的运动。

**知识点3.6　循环坐标与守恒量**

定义：不出现在 $L$ 中的坐标 $q_k$ 称为循环坐标。

推导：若 $\partial L/\partial q_k=0$，则欧拉-拉格朗日方程给出
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_k}=0\ \Rightarrow\ p_k=\frac{\partial L}{\partial\dot q_k}=\text{常量}.$$
（诺特定理特例。）

物理意义：对称性对应守恒量；平移对称→动量，旋转对称→角动量，时间平移→能量。

极限特例：中心力问题中 $\theta$ 循环，$p_\theta=mr^2\dot\theta=L$。

典型应用：陀螺运动中的进动守恒量分析。

**知识点3.7　能量函数与广义能量积分**

定义：能量函数 $h=\sum_j\dot q_j\frac{\partial L}{\partial\dot q_j}-L$。

推导：
$$\frac{\mathrm dh}{\mathrm dt}=\sum_j\ddot q_j\frac{\partial L}{\partial\dot q_j}+\sum_j\dot q_j\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}-\frac{\mathrm dL}{\mathrm dt}.$$
而
$$\frac{\mathrm dL}{\mathrm dt}=\sum_j\frac{\partial L}{\partial q_j}\dot q_j+\sum_j\frac{\partial L}{\partial\dot q_j}\ddot q_j+\frac{\partial L}{\partial t}.$$
由拉格朗日方程 $\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}=\frac{\partial L}{\partial q_j}$，代入相消得
$$\frac{\mathrm dh}{\mathrm dt}=-\frac{\partial L}{\partial t}.$$

物理意义：若 $L$ 不显含时间，则 $h=$ 常量（广义能量守恒）。当势能与速度无关且约束定常时 $h=T+V$。

极限特例：非惯性势（含时约束）时 $h\neq T+V$。

典型应用：旋转参考系中雅可比积分。

**知识点3.8　非惯性系中的拉格朗日方法**

定义：以平动加速 + 转动参考系。

推导：惯性系坐标 $\mathbf r$ 与动系 $\mathbf r'$ 关系
$$\mathbf r=\mathbf R+\mathbf r'.$$
代入 $L=\frac12 m\dot{\mathbf r}^2-V$，动系速度
$$\mathbf v=\dot{\mathbf R}+\dot{\mathbf r}'+\boldsymbol\omega\times\mathbf r'.$$
展开动能得离心势 $-\frac12 m(\boldsymbol\omega\times\mathbf r')^2$ 与科氏势 $-m\mathbf v'\cdot(\boldsymbol\omega\times\mathbf r')$ 等项。拉格朗日方程自动给出惯性力：
$$\mathbf F_{\text{离心}}=-m\boldsymbol\omega\times(\boldsymbol\omega\times\mathbf r'),\qquad \mathbf F_{\text{科氏}}=-2m\boldsymbol\omega\times\mathbf v'.$$

物理意义：非惯性效应通过动力学项的势能化处理。

极限特例：$\omega=0$ 退化为平动坐标系。

典型应用：地球自转对抛体、傅科摆的影响。

### 第4章 变分原理与哈密顿原理

**知识点4.1　泛函与泛函极值**

定义：泛函 $J[y]=\int_{x_0}^{x_1}F(x,y,y')\,\mathrm dx$，其变量是函数 $y(x)$。

推导变分定义：设 $y(x;\alpha)=y(x)+αη(x)$，$\eta(x_0)=\eta(x_1)=0$。极值条件
$$\left.\frac{\partial J}{\partial\alpha}\right|_{\alpha=0}=0.$$
计算
$$\frac{\partial J}{\partial\alpha}=\int_{x_0}^{x_1}\left(\frac{\partial F}{\partial y}\eta+\frac{\partial F}{\partial y'}\eta'\right)\mathrm dx.$$
分部积分第二项
$$\int_{x_0}^{x_1}\frac{\partial F}{\partial y'}\eta'\,\mathrm dx=\left[\frac{\partial F}{\partial y'}\eta\right]_{x_0}^{x_1}-\int_{x_0}^{x_1}\frac{\mathrm d}{\mathrm dx}\frac{\partial F}{\partial y'}\eta\,\mathrm dx.$$
边界项因 $\eta=0$ 为零，故
$$\delta J=\int_{x_0}^{x_1}\left(\frac{\partial F}{\partial y}-\frac{\mathrm d}{\mathrm dx}\frac{\partial F}{\partial y'}\right)\eta\,\mathrm dx=0.$$
由 $\eta$ 任意得欧拉-拉格朗日方程
$$\frac{\partial F}{\partial y}-\frac{\mathrm d}{\mathrm dx}\frac{\partial F}{\partial y'}=0.$$

物理意义：泛函驻值条件给出微分方程；变分法把全局极值问题化为局域方程。

极限特例：$F=w y\sqrt{1+y'^2}$ 时给出最速降线的摆线。

典型应用：最速降线、最小旋转面。

**知识点4.2　最速降线的完整推导与应用**

定义：在重力场中求质点从 $O$ 到 $P$ 时间最短的曲线。

推导：设 $y$ 向下为正，速度由能量守恒
$$\frac12 mv^2=mgy\ \Rightarrow\ v=\sqrt{2gy}.$$
总时间
$$J[y]=\int\frac{\mathrm ds}{v}=\int_0^{x_1}\frac{\sqrt{1+y'^2}}{\sqrt{2gy}}\mathrm dx.$$
被积函数 $F=\sqrt{(1+y'^2)/(2gy)}$ 不显含 $x$，用贝尔特拉米恒等式
$$F-y'\frac{\partial F}{\partial y'}=C.$$
计算
$$\frac{\partial F}{\partial y'}=\frac{y'}{\sqrt{2gy}\sqrt{1+y'^2}},$$
$$F-y'\frac{\partial F}{\partial y'}=\frac{1}{\sqrt{2gy}\sqrt{1+y'^2}}=C.$$
即 $y(1+y'^2)=1/(2gC^2)=2a$。解
$$y=2a\sin^2\frac{\theta}{2}=a(1-\cos\theta),\quad x=a(\theta-\sin\theta),$$
是摆线。

物理意义：最快路径不是直线而是摆线，说明速度大处应走更陡的路。

极限特例：重力趋于零时退化为直线（此时无时间差）。

典型应用：滑雪道、光学费马原理类问题。

**知识点4.3　哈密顿原理**

定义：真实运动使作用量
$$S=\int_{t_1}^{t_2}L(q,\dot q,t)\,\mathrm dt$$
取驻值，且端点 $q(t_1),q(t_2)$ 固定。

推导：由知识点4.1 对每个坐标 $q_j$ 变分，直接得欧拉-拉格朗日方程。即
$$\delta S=0\ \Leftrightarrow\ \frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}=0.$$
更严格地，$\delta S=\int_{t_1}^{t_2}\sum_j\left(\frac{\partial L}{\partial q_j}\delta q_j+\frac{\partial L}{\partial\dot q_j}\delta\dot q_j\right)\mathrm dt$，其中 $\delta\dot q_j=\mathrm d(\delta q_j)/\mathrm dt$，分部积分并用端点固定。

物理意义：哈密顿原理是积分形式的最小作用量原理，对保守与非保守（广义力）系统均适用。

极限特例：自由粒子 $L=\frac12 m\dot x^2$ 给出匀速直线运动。

典型应用：连续体、场论中的拉格朗日密度方法。

**知识点4.4　诺特定理（对称性与守恒律）**

定义：单参数连续变换 $q_j\to q_j+\epsilon\,\phi_j(q)$ 保持作用量不变。

推导：作用量变分
$$\delta S=\epsilon\int_{t_1}^{t_2}\sum_j\left(\frac{\partial L}{\partial q_j}\phi_j+\frac{\partial L}{\partial\dot q_j}\dot\phi_j\right)\mathrm dt.$$
由欧拉-拉格朗日方程 $\partial L/\partial q_j=\frac{\mathrm d}{\mathrm dt}\partial L/\partial\dot q_j$，
$$\delta S=\epsilon\int_{t_1}^{t_2}\frac{\mathrm d}{\mathrm dt}\sum_j\frac{\partial L}{\partial\dot q_j}\phi_j\,\mathrm dt=\epsilon\left[\sum_j p_j\phi_j\right]_{t_1}^{t_2}.$$
若 $\delta S=0$ 对任意区间，则被积函数时间导数为零：
$$Q=\sum_j p_j\phi_j=\text{常量}.$$

物理意义：每一种连续对称性对应一个守恒量。

极限特例：时间平移对称 $\phi_j=\dot q_j$ 给出能量守恒；空间平移给出动量；旋转给出角动量。

典型应用：场论中能量-动量张量、电荷守恒。

**知识点4.5　哈密顿原理与非保守力**

定义：非保守广义力 $Q_j^{\text{nc}}$。

推导：对非保守系统，
$$\delta S+\int_{t_1}^{t_2}\sum_j Q_j^{\text{nc}}\delta q_j\,\mathrm dt=0,$$
得
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}=Q_j^{\text{nc}}.$$
若存在耗散函数 $\mathcal F=\frac12\sum_{jk}c_{jk}\dot q_j\dot q_k$，则 $Q_j^{\text{nc}}=-\partial\mathcal F/\partial\dot q_j$。

物理意义：耗散函数使能量单调减少。

极限特例：瑞利阻尼 $\mathcal F=\frac12 c\dot q^2$ 给出阻尼力 $-c\dot q$。

典型应用：结构阻尼振动建模。

**知识点4.6　雅可比度量原理与最短路径类比**

定义：对能量守恒系统，作用量 $W=\int\sum_j p_j\,\mathrm dq_j$（雅可比作用量）。

推导：由 $T=\frac12\sum_{jk}a_{jk}\dot q_j\dot q_k$，能量 $E=T+V$ 得 $\mathrm dt=\mathrm ds/\sqrt{2(E-V)}$，其中
$$\mathrm ds^2=2(E-V)\sum_{jk}a_{jk}\,\mathrm dq_j\mathrm dq_k.$$
雅可比作用量
$$W=\int\sqrt{2(E-V)\sum_{jk}a_{jk}\,\mathrm dq_j\mathrm dq_k},$$
其变分等价于黎曼流形上的测地线问题。

物理意义：保守力学系统轨道等价于度规为 $g_{jk}=2(E-V)a_{jk}$ 的空间中的测地线。

极限特例：$V=0$ 时即自由粒子沿直线。

典型应用：定能量轨道几何化、天体力学。

### 第5章 哈密顿力学与正则方程

**知识点5.1　勒让德变换与哈密顿函数**

定义：广义动量 $p_j=\partial L/\partial\dot q_j$；哈密顿函数
$$H(q,p,t)=\sum_j p_j\dot q_j-L(q,\dot q,t),$$
其中 $\dot q$ 用 $p,q$ 表示。

推导：设关系 $p_j=p_j(q,\dot q,t)$ 可逆解出 $\dot q_j=\dot q_j(q,p,t)$。取微分
$$\mathrm dH=\sum_j\dot q_j\,\mathrm dp_j+\sum_j p_j\,\mathrm d\dot q_j-\sum_j\frac{\partial L}{\partial q_j}\mathrm dq_j-\sum_j\frac{\partial L}{\partial\dot q_j}\mathrm d\dot q_j-\frac{\partial L}{\partial t}\mathrm dt.$$
因 $p_j=\partial L/\partial\dot q_j$，含 $\mathrm d\dot q_j$ 的项相消：
$$\mathrm dH=\sum_j\dot q_j\,\mathrm dp_j-\sum_j\frac{\partial L}{\partial q_j}\mathrm dq_j-\frac{\partial L}{\partial t}\mathrm dt.$$
比较得
$$\dot q_j=\frac{\partial H}{\partial p_j},\qquad \frac{\partial L}{\partial q_j}=-\frac{\partial H}{\partial q_j},\qquad \frac{\partial H}{\partial t}=-\frac{\partial L}{\partial t}.$$
由拉格朗日方程 $\dot p_j=\partial L/\partial q_j=-\partial H/\partial q_j$，于是得到哈密顿正则方程。

物理意义：$q,p$ 对称，相空间运动由单一标量 $H$ 生成。

极限特例：保守系统 $H=T+V$ 且 $H=E=$ 常量。

典型应用：数值积分辛格式、天体轨道与量子化的历史来源。

**知识点5.2　哈密顿正则方程及其性质**

定义：
$$\dot q_j=\frac{\partial H}{\partial p_j},\qquad \dot p_j=-\frac{\partial H}{\partial q_j},\qquad j=1,\dots,n.$$

推导：见 5.1。能量守恒：
$$\frac{\mathrm dH}{\mathrm dt}=\sum_j\frac{\partial H}{\partial q_j}\dot q_j+\sum_j\frac{\partial H}{\partial p_j}\dot p_j+\frac{\partial H}{\partial t}=\sum_j\frac{\partial H}{\partial q_j}\frac{\partial H}{\partial p_j}-\sum_j\frac{\partial H}{\partial p_j}\frac{\partial H}{\partial q_j}+\frac{\partial H}{\partial t}=\frac{\partial H}{\partial t}.$$
故 $H$ 不显含 $t$ 时守恒。

物理意义：相空间速度场是无散度的（刘维尔定理来源）。

极限特例：一维 $H=p^2/(2m)+V(q)$ 给出 $\dot p=-V'(q)$。

典型应用：辛积分器（Störmer-Verlet）保持相体积。

**知识点5.3　泊松括号**

定义：
$$\{f,g\}=\sum_j\left(\frac{\partial f}{\partial q_j}\frac{\partial g}{\partial p_j}-\frac{\partial f}{\partial p_j}\frac{\partial g}{\partial q_j}\right).$$

推导任意力学量演化：
$$\frac{\mathrm df}{\mathrm dt}=\sum_j\left(\frac{\partial f}{\partial q_j}\dot q_j+\frac{\partial f}{\partial p_j}\dot p_j\right)+\frac{\partial f}{\partial t}=\sum_j\left(\frac{\partial f}{\partial q_j}\frac{\partial H}{\partial p_j}-\frac{\partial f}{\partial p_j}\frac{\partial H}{\partial q_j}\right)+\frac{\partial f}{\partial t},$$
即
$$\frac{\mathrm df}{\mathrm dt}=\{f,H\}+\frac{\partial f}{\partial t}.$$
基本括号 $\{q_i,q_j\}=0,\ \{p_i,p_j\}=0,\ \{q_i,p_j\}=\delta_{ij}$。

物理意义：守恒量满足 $\{f,H\}=0$。

极限特例：雅可比恒等式 $\{f,\{g,h\}\}+\{g,\{h,f\}\}+\{h,\{f,g\}\}=0$。

典型应用：判断守恒量、量子力学对应。

**知识点5.4　正则变换**

定义：变换 $Q_j=Q_j(q,p,t),P_j=P_j(q,p,t)$ 保持正则方程形式。

推导：要求 $\sum_j p_j\mathrm dq_j-H\mathrm dt=\sum_j P_j\mathrm dQ_j-K\mathrm dt+\mathrm dF$。所谓母函数生成。若取 $F=F_1(q,Q,t)$，则
$$p_j=\frac{\partial F_1}{\partial q_j},\qquad P_j=-\frac{\partial F_1}{\partial Q_j},\qquad K=H+\frac{\partial F_1}{\partial t}.$$
若取 $F=F_2(q,P,t)=F_1+\sum_j Q_jP_j$，则
$$p_j=\frac{\partial F_2}{\partial q_j},\qquad Q_j=\frac{\partial F_2}{\partial P_j},\qquad K=H+\frac{\partial F_2}{\partial t}.$$

物理意义：正则变换是保持辛结构 $\{Q_i,P_j\}=\delta_{ij}$ 的相空间坐标变换。

极限特例：恒等变换 $F_2=\sum_j q_jP_j$。

典型应用：谐振子对角化、哈密顿-雅可比方程。

**知识点5.5　哈密顿-雅可比方程**

定义：寻找母函数 $S(q,P,t)$ 使 $K=0$。

推导：由 $F_2=S$，$p_j=\partial S/\partial q_j$，$Q_j=\partial S/\partial P_j$，$K=H+\partial S/\partial t=0$，即
$$H\left(q_1,\dots,q_n,\frac{\partial S}{\partial q_1},\dots,\frac{\partial S}{\partial q_n},t\right)+\frac{\partial S}{\partial t}=0.$$
若 $H$ 不显含 $t$，分离 $S=W(q,P)-Et$，得定态方程
$$H\left(q,\frac{\partial W}{\partial q}\right)=E.$$

物理意义：一阶偏微分方程，其解给出轨道族；特征线与正则方程一致。

极限特例：一维 $W'=\pm\sqrt{2m(E-V)}$，积分给出 $\int\sqrt{2m(E-V)}\,\mathrm dq=Et$。

典型应用：可分离系统的求解、几何光学类比。

**知识点5.6　刘维尔定理**

定义：相空间密度 $\rho(q,p,t)$。

推导：相空间流不可压缩，$\nabla\cdot\mathbf v_{\text{相}}=0$，
$$\sum_j\left(\frac{\partial\dot q_j}{\partial q_j}+\frac{\partial\dot p_j}{\partial p_j}\right)=\sum_j\left(\frac{\partial^2H}{\partial q_j\partial p_j}-\frac{\partial^2H}{\partial p_j\partial q_j}\right)=0.$$
连续性方程
$$\frac{\partial\rho}{\partial t}+\sum_j\left(\frac{\partial(\rho\dot q_j)}{\partial q_j}+\frac{\partial(\rho\dot p_j)}{\partial p_j}\right)=0$$
化为
$$\frac{\partial\rho}{\partial t}+\{\rho,H\}=0,\qquad\text{即}\quad\frac{\mathrm d\rho}{\mathrm dt}=0.$$

物理意义：相空间体积在演化中保持，是统计力学的几何基础。

极限特例：平衡态时 $\{\rho,H\}=0$。

典型应用：微正则系综、辛积分器稳定性。

**知识点5.7　作用量-角变量**

定义：对周期运动，角变量 $w_j=\partial W/\partial J_j$，作用变量 $J_j=\oint p_j\,\mathrm dq_j$。

推导：对一维周期系统
$$J=\oint p\,\mathrm dq=\oint\sqrt{2m(E-V)}\,\mathrm dq.$$
因 $H=H(J)$ 与 $w$ 无关，正则方程
$$\dot w=\frac{\partial H}{\partial J}=\nu(J)=\text{常量},\qquad \dot J=0.$$
故 $w=\nu t+w_0$ 线性增长。

物理意义：作用变量在绝热不变量意义上近似守恒；角变量是频率乘以时间。

极限特例：谐振子 $J=E/\omega$，$\nu=\omega/2\pi$。

典型应用：摄动论、玻尔-索末菲量子化条件 $\oint p\,\mathrm dq=nh$。

### 第6章 刚体动力学

**知识点6.1　刚体运动学与欧拉角**

定义：刚体上任意两点距离不变；用质心位置和姿态确定位形。

推导欧拉角：取固连坐标系，从参考系经三次转动 $R_z(\psi)R_x(\theta)R_z(\phi)$（$z$-$x$-$z$ 约定）得到。方向余弦矩阵 $\mathbf R\in SO(3)$，$\mathbf R^{\mathsf T}\mathbf R=\mathbf I,\ \det\mathbf R=1$。

物理意义：刚体有 6 个自由度（质心 3 + 姿态 3）。

极限特例：万向节死锁发生在 $\theta=0,\pi$。

典型应用：飞行器姿态、陀螺导航。

**知识点6.2　角速度与欧拉运动学方程**

定义：角速度 $\boldsymbol\omega$ 满足对固连矢量 $\dot{\mathbf e}=\boldsymbol\omega\times\mathbf e$。

推导欧拉角速率与角速度分量（$z$-$x$-$z$）：
$$\dot{\mathbf R}\mathbf R^{\mathsf T}=\tilde{\boldsymbol\omega}=\begin{pmatrix}0&-\omega_z&\omega_y\\ \omega_z&0&-\omega_x\\ -\omega_y&\omega_x&0\end{pmatrix}.$$
在体轴系下
$$\omega_x=\dot\phi\sin\theta\sin\psi+\dot\theta\cos\psi,$$
$$\omega_y=\dot\phi\sin\theta\cos\psi-\dot\theta\sin\psi,$$
$$\omega_z=\dot\phi\cos\theta+\dot\psi.$$
（或对称地写为 $\boldsymbol\omega=\dot\psi\hat{\mathbf z}+\dot\theta\hat{\mathbf n}+\dot\phi\hat{\mathbf e}_3$。）

物理意义：角速度是伪矢量，不满足平行四边形简单交换。

极限特例：平面运动 $\theta=0$，$\omega_z=\dot\phi+\dot\psi$。

典型应用：惯导姿态解算。

**知识点6.3　惯量张量与平行轴定理**

定义：对质心的惯量张量
$$\mathbf I_C=\int(\mathbf r'^2\mathbf I_3-\mathbf r'\otimes\mathbf r')\,\mathrm dm.$$

推导平行轴定理：设参考点 $O$，$\mathbf r=\mathbf R_C+\mathbf r'$，则
$$\mathbf I_O=\int[(\mathbf R_C+\mathbf r')^2\mathbf I_3-(\mathbf R_C+\mathbf r')\otimes(\mathbf R_C+\mathbf r')]\mathrm dm.$$
展开后一次项含 $\int\mathbf r'\,\mathrm dm=\mathbf 0$（质心系），得
$$\mathbf I_O=\mathbf I_C+m(R_C^2\mathbf I_3-\mathbf R_C\otimes\mathbf R_C).$$
标量形式对轴 $a$：$I_a=I_{C,a}+md^2$。

物理意义：惯量张量对称正定，可用主轴对角化。

极限特例：细杆绕端点 $I=\frac13 mL^2$，绕中心 $\frac1{12}mL^2$。

典型应用：机械臂惯量计算。

**知识点6.4　角动量与惯量张量关系**

定义：$\mathbf L=\mathbf I\boldsymbol\omega$（对质心或定点）。

推导：$\mathbf L=\int\mathbf r\times(\boldsymbol\omega\times\mathbf r)\,\mathrm dm=\int[r^2\boldsymbol\omega-(\mathbf r\cdot\boldsymbol\omega)\mathbf r]\mathrm dm$，写成矩阵即 $\mathbf L=\mathbf I\boldsymbol\omega$。

物理意义：$\mathbf L$ 与 $\boldsymbol\omega$ 一般不平行，仅主轴方向平行。

极限特例：球对称刚体 $I$ 各向同性，$\mathbf L=I\boldsymbol\omega$。

典型应用：陀螺力矩分析。

**知识点6.5　欧拉动力学方程**

定义：对体轴系的角动量定理。

推导：在惯性系 $\dot{\mathbf L}= \mathbf M$。转到体轴系，矢量求导需加 $\boldsymbol\omega\times$：
$$\left(\frac{\mathrm d\mathbf L}{\mathrm dt}\right)_{\text{惯性}}=\left(\frac{\mathrm d\mathbf L}{\mathrm dt}\right)_{\text{体}}+\boldsymbol\omega\times\mathbf L.$$
取主轴系 $\mathbf L=(I_1\omega_1,I_2\omega_2,I_3\omega_3)$，得
$$I_1\dot\omega_1+(I_3-I_2)\omega_2\omega_3=M_1,$$
$$I_2\dot\omega_2+(I_1-I_3)\omega_3\omega_1=M_2,$$
$$I_3\dot\omega_3+(I_2-I_1)\omega_1\omega_2=M_3.$$

物理意义：陀螺项 $(I_3-I_2)\omega_2\omega_3$ 是角动量方向变化产生的耦合。

极限特例：自由欧拉陀螺 $M_i=0$ 有对 $\mathbf L$ 的进动。

典型应用：卫星自旋稳定、网球拍定理。

**知识点6.6　自由对称陀螺与欧拉角解**

定义：$I_1=I_2=I_\perp\neq I_3$，无力矩。

推导：由欧拉方程第三式
$$I_3\dot\omega_3=0\ \Rightarrow\ \omega_3=\text{常量}.$$
第一式 $I_\perp\dot\omega_1=(I_\perp-I_3)\omega_2\omega_3$，第二式 $I_\perp\dot\omega_2=(I_3-I_\perp)\omega_3\omega_1$。令 $\Omega=\frac{I_3-I_\perp}{I_\perp}\omega_3$，则
$$\dot\omega_1=-\Omega\,\omega_2,\qquad \dot\omega_2=\Omega\,\omega_1,$$
解得横向角速度以 $\Omega$ 旋转：
$$\omega_1=A\cos\Omega t,\quad\omega_2=A\sin\Omega t.$$
再由欧拉角运动学得进动率 $\dot\psi$ 与自旋 $\dot\phi$ 比。

物理意义：自旋轴绕总角动量做圆锥进动。

极限特例：球陀螺 $I_3=I_\perp$ 时 $\Omega=0$，$\boldsymbol\omega$ 恒定。

典型应用：子弹自旋稳定、陀螺罗盘。

**知识点6.7　拉格朗日陀螺与进动**

定义：重力场中定点陀螺，拉格朗日函数
$$L=\frac12 I_\perp(\dot\theta^2+\dot\phi^2\sin^2\theta)+\frac12 I_3(\dot\phi\cos\theta+\dot\psi)^2-Mgl\cos\theta.$$

推导：$\psi,\phi$ 为循环坐标，得守恒量
$$p_\psi=I_3(\dot\psi+\dot\phi\cos\theta)=L_3=\text{常量},$$
$$p_\phi=I_\perp\dot\phi\sin^2\theta+I_3(\dot\psi+\dot\phi\cos\theta)\cos\theta=L_z=\text{常量}.$$
能量守恒给出 $E=\frac12 I_\perp\dot\theta^2+V_{\text{eff}}(\theta)$，其中
$$V_{\text{eff}}(\theta)=\frac{(L_z-L_3\cos\theta)^2}{2I_\perp\sin^2\theta}+\frac{L_3^2}{2I_3}+Mgl\cos\theta.$$
由 $V_{\text{eff}}'(\theta_0)=0$ 得稳定进动条件。

物理意义：陀螺的规则进动来自两个守恒量的组合。

极限特例：$L_z=L_3\cos\theta$ 时 $\dot\phi=0$，纯自旋。

典型应用：陀螺仪、弹道姿态。

**知识点6.8　刚体能量与柯尼希定理（三维）**

定义：刚体动能 $T=\frac12 m v_C^2+\frac12\boldsymbol\omega\cdot\mathbf I_C\boldsymbol\omega$。

推导：由第二章柯尼希定理，体动能为质心系动能
$$T_C=\frac12\int(\boldsymbol\omega\times\mathbf r')^2\mathrm dm.$$
展开 $(\boldsymbol\omega\times\mathbf r')^2=(\boldsymbol\omega\times\mathbf r')\cdot(\boldsymbol\omega\times\mathbf r')=\boldsymbol\omega\cdot[\mathbf r'\times(\boldsymbol\omega\times\mathbf r')]=\boldsymbol\omega\cdot[(\mathbf r'^2\mathbf I_3-\mathbf r'\otimes\mathbf r')\boldsymbol\omega]$，故
$$T_C=\frac12\boldsymbol\omega\cdot\mathbf I_C\boldsymbol\omega.$$

物理意义：三维转动动能是角速度的二次型，由惯量张量控制。

极限特例：主轴下 $T_C=\frac12(I_1\omega_1^2+I_2\omega_2^2+I_3\omega_3^2)$。

典型应用：多刚体系统拉格朗日建模。

### 第7章 连续介质运动学与应力分析

**知识点7.1　连续介质假设与两种描述方法**

定义：连续介质假设认为物质在宏观上充满空间，物理量（密度、速度、应力）是空间与时间的连续可微函数；忽略分子尺度涨落。

拉格朗日描述：跟踪物质点，独立变量 $(\mathbf X,t)$，$\mathbf x=\boldsymbol\chi(\mathbf X,t)$。
欧拉描述：在空间点观察，独立变量 $(\mathbf x,t)$。

推导物质导数：设 $f=f(\mathbf x,t)$，沿质点轨迹 $\mathbf x=\boldsymbol\chi(\mathbf X,t)$，
$$\frac{\mathrm Df}{\mathrm Dt}=\left.\frac{\partial f}{\partial t}\right|_{\mathbf X}=\frac{\partial f}{\partial t}+\sum_i\frac{\partial f}{\partial x_i}\frac{\mathrm dx_i}{\mathrm dt}=\frac{\partial f}{\partial t}+\mathbf v\cdot\nabla f.$$

物理意义：物质导数 = 局地变化 + 对流变化。

极限特例：定常场 $\partial f/\partial t=0$ 时随体变化全由对流项产生。

典型应用：流体加速、传热。

**知识点7.2　变形梯度与应变张量**

定义：变形梯度 $\mathbf F=\partial\mathbf x/\partial\mathbf X$，$F_{iJ}=\partial x_i/\partial X_J$（大写指标材料坐标，小写指标空间坐标）。

推导：两物质点微元 $\mathrm d\mathbf X$ 变形后
$$\mathrm dx_i=\sum_J F_{iJ}\,\mathrm dX_J.$$
长度平方差
$$|\mathrm d\mathbf x|^2-|\mathrm d\mathbf X|^2=\mathrm d\mathbf X^{\mathsf T}(\mathbf F^{\mathsf T}\mathbf F-\mathbf I)\mathrm d\mathbf X=2\,\mathrm d\mathbf X^{\mathsf T}\mathbf E\,\mathrm d\mathbf X,$$
其中格林-拉格朗日应变
$$\mathbf E=\frac12(\mathbf F^{\mathsf T}\mathbf F-\mathbf I).$$
小变形下 $\mathbf F=\mathbf I+\nabla_{\mathbf X}\mathbf u$，得小应变张量
$$\boldsymbol\varepsilon=\frac12(\nabla\mathbf u+\nabla\mathbf u^{\mathsf T}),\qquad \varepsilon_{ij}=\frac12\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\right).$$

物理意义：应变张量度量长度变化的对称部分；反对称部分是刚体转动。

极限特例：纯平移 $\mathbf u=$ 常量时 $\boldsymbol\varepsilon=\mathbf 0$。

典型应用：材料应变测量、有限元。

**知识点7.3　小应变张量的几何意义与协调方程**

定义：$\varepsilon_{ii}$ 为线应变，$\varepsilon_{ij}(i\neq j)$ 为工程剪应变之半。

推导线应变几何：取材料线元沿 $\hat{\mathbf e}_1$，长度变化率
$$\frac{\mathrm d\ell-\mathrm dL}{\mathrm dL}\approx\varepsilon_{11}=\frac{\partial u_1}{\partial x_1}.$$
剪应变：两正交线元夹角变化 $\gamma_{12}=2\varepsilon_{12}=\partial u_1/\partial x_2+\partial u_2/\partial x_1$。

推导协调方程：位移单值要求应变满足
$$\frac{\partial^2\varepsilon_{11}}{\partial x_2^2}+\frac{\partial^2\varepsilon_{22}}{\partial x_1^2}=2\frac{\partial^2\varepsilon_{12}}{\partial x_1\partial x_2},$$
及其循环置换。

物理意义：应变场不能任意，必须由位移场导出。

极限特例：均匀应变场自动满足协调。

典型应用：应变花测量、光弹性。

**知识点7.4　速度梯度、应变率与旋度**

定义：速度梯度 $\mathbf L=\nabla\mathbf v$，$L_{ij}=\partial v_i/\partial x_j$。

推导分解：
$$\mathbf L=\mathbf D+\mathbf W,\quad \mathbf D=\frac12(\mathbf L+\mathbf L^{\mathsf T}),\quad \mathbf W=\frac12(\mathbf L-\mathbf L^{\mathsf T}).$$
$\mathbf D$ 为应变率张量（对称），$\mathbf W$ 为旋转率张量（反对称）。角速度矢量 $\boldsymbol\omega$ 满足 $W_{ij}=-\epsilon_{ijk}\omega_k$，即
$$\boldsymbol\omega=\frac12\nabla\times\mathbf v.$$

物理意义：流体微团运动 = 平移 + 线变形 + 剪切 + 旋转。

极限特例：刚性旋转 $\mathbf D=0$，$\nabla\times\mathbf v=2\boldsymbol\omega$。

典型应用：涡量动力学、湍流。

**知识点7.5　雷诺输运定理**

定义：系统（物质体）物理量随体变化与固定控制体的关系。

推导：设系统在 $t$ 占据 $V(t)$，$\Phi=\int_{V(t)}\phi\,\mathrm dV$，
$$\frac{\mathrm D\Phi}{\mathrm Dt}=\frac{\mathrm d}{\mathrm dt}\int_{V(t)}\phi\,\mathrm dV.$$
对含时体积求导含体积变化项。体积微元变化率 $\frac{1}{\mathrm dV}\frac{\mathrm D(\mathrm dV)}{\mathrm Dt}=\nabla\cdot\mathbf v$，故
$$\frac{\mathrm D}{\mathrm Dt}\int_V\phi\,\mathrm dV=\int_V\left(\frac{\partial\phi}{\partial t}+\nabla\cdot(\phi\mathbf v)\right)\mathrm dV=\int_V\frac{\partial\phi}{\partial t}\mathrm dV+\oint_S\phi\mathbf v\cdot\mathrm d\mathbf S.$$

物理意义：系统量的变化 = 控制体内局地变化 + 通过边界的通量。

极限特例：不可压 $\nabla\cdot\mathbf v=0$ 时体积元守恒。

典型应用：连续性方程、动量方程的积分形式。

**知识点7.6　应力张量与柯西公式**

定义：应力矢量 $\mathbf t^{(\mathbf n)}=\lim_{\Delta S\to0}\Delta\mathbf F/\Delta S$，表示法向为 $\mathbf n$ 的面上单位面积作用力。

推导柯西公式：取四面体微元，斜面积 $\Delta S$，三个坐标面投影 $\Delta S_i=\Delta S\,n_i$。动量平衡在趋于零时惯性项与体积力为高阶小量，
$$\mathbf t^{(\mathbf n)}\Delta S=\mathbf t^{(1)}\Delta S_1+\mathbf t^{(2)}\Delta S_2+\mathbf t^{(3)}\Delta S_3,$$
$$\Rightarrow \mathbf t^{(\mathbf n)}=\sum_i n_i\mathbf t^{(i)}=\boldsymbol\sigma\mathbf n,\qquad t_j=\sum_i\sigma_{ij}n_i.$$
即柯西应力定理。

物理意义：应力是二阶张量，只需三个坐标面上的应力矢量即可确定任意面应力。

极限特例：流体静止时 $\sigma_{ij}=-p\delta_{ij}$。

典型应用：结构强度分析。

**知识点7.7　应力张量对称性与应力平衡方程**

定义：$\boldsymbol\sigma$ 的对称性由角动量守恒证明。

推导力矩平衡：对微元取矩，
$$\int_V(\mathbf r\times\rho\mathbf b)\mathrm dV+\oint_S\mathbf r\times(\boldsymbol\sigma\mathbf n)\mathrm dS=\frac{\mathrm D}{\mathrm Dt}\int_V\mathbf r\times\rho\mathbf v\,\mathrm dV.$$
取体积趋于零，内力矩项给出 $\sigma_{ij}=\sigma_{ji}$（应力张量对称）。

推导平衡方程：对微元用高斯散度定理，
$$\oint_S\boldsymbol\sigma\mathbf n\,\mathrm dS=\int_V\nabla\cdot\boldsymbol\sigma\,\mathrm dV,$$
代入动量方程得
$$\nabla\cdot\boldsymbol\sigma+\rho\mathbf b=\rho\frac{\mathrm D\mathbf v}{\mathrm Dt},\qquad \frac{\partial\sigma_{ij}}{\partial x_j}+\rho b_i=\rho\dot v_i.$$
静力平衡时右端为零：
$$\frac{\partial\sigma_{ij}}{\partial x_j}+\rho b_i=0.$$

物理意义：应力散度是单位体积内力，平衡方程是连续介质运动定律。

极限特例：静水压力 $\sigma_{ij}=-p\delta_{ij}$ 给出 $\nabla p=\rho\mathbf b$。

典型应用：弹性静力学、流体静力学。

**知识点7.8　主应力与应力不变量**

定义：主应力是使 $\mathbf t^{(\mathbf n)}\parallel\mathbf n$ 的 $\sigma$。

推导：$(\boldsymbol\sigma-\sigma\mathbf I)\mathbf n=\mathbf 0$ 有非零解要求
$$\det(\boldsymbol\sigma-\sigma\mathbf I)=0.$$
展开得特征方程
$$\sigma^3-I_1\sigma^2+I_2\sigma-I_3=0,$$
$$I_1=\sigma_{kk},\quad I_2=\frac12(\sigma_{ii}\sigma_{jj}-\sigma_{ij}\sigma_{ji}),\quad I_3=\det\boldsymbol\sigma.$$
$I_1,I_2,I_3$ 为坐标不变量。

物理意义：主应力与主方向表征一点应力状态本质；最大剪应力 $\tau_{\max}=( \sigma_1-\sigma_3)/2$。

极限特例：$\sigma_{ij}=-p\delta_{ij}$ 时三主应力相等。

典型应用：屈服准则（von Mises、Tresca）。

### 第8章 守恒定律与Navier-Stokes方程

**知识点8.1　质量守恒与连续性方程**

定义：系统质量守恒 $\frac{\mathrm D}{\mathrm Dt}\int_V\rho\,\mathrm dV=0$。

推导：由雷诺输运定理
$$\frac{\mathrm D}{\mathrm Dt}\int_V\rho\,\mathrm dV=\int_V\left(\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf v)\right)\mathrm dV=0.$$
对任意体积成立，故
$$\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf v)=0,\qquad \frac{\partial\rho}{\partial t}+\frac{\partial(\rho v_i)}{\partial x_i}=0.$$
利用物质导数展开 $\nabla\cdot(\rho\mathbf v)=\mathbf v\cdot\nabla\rho+\rho\nabla\cdot\mathbf v$，得
$$\frac{\mathrm D\rho}{\mathrm Dt}+\rho\nabla\cdot\mathbf v=0.$$

物理意义：密度变化与速度散度相消；不可压缩时 $\nabla\cdot\mathbf v=0$。

极限特例：定常一维 $(\rho v)'=0$，$\rho v=$ 常量。

典型应用：管道流动、激波关系（R-H条件）。

**知识点8.2　动量守恒与动量方程推导**

定义：动量变化率等于体积力与面力之和。

推导：对物质体
$$\frac{\mathrm D}{\mathrm Dt}\int_V\rho v_i\,\mathrm dV=\int_V\rho b_i\,\mathrm dV+\oint_S\sigma_{ij}n_j\,\mathrm dS.$$
由雷诺输运定理与散度定理，
$$\int_V\left(\frac{\partial(\rho v_i)}{\partial t}+\frac{\partial(\rho v_i v_j)}{\partial x_j}\right)\mathrm dV=\int_V\left(\rho b_i+\frac{\partial\sigma_{ij}}{\partial x_j}\right)\mathrm dV.$$
故
$$\frac{\partial(\rho v_i)}{\partial t}+\frac{\partial(\rho v_iv_j)}{\partial x_j}=\rho b_i+\frac{\partial\sigma_{ij}}{\partial x_j}.$$
利用连续性方程化简左端：
$$\frac{\partial(\rho v_i)}{\partial t}+\frac{\partial(\rho v_iv_j)}{\partial x_j}=\rho\frac{\mathrm Dv_i}{\mathrm Dt},$$
得柯西运动方程
$$\rho\frac{\mathrm Dv_i}{\mathrm Dt}=\frac{\partial\sigma_{ij}}{\partial x_j}+\rho b_i.$$

物理意义：单位体积惯性力 = 体积力 + 应力散度。

极限特例：无粘静止得欧拉静力方程。

典型应用：流体动力学、固体动力学统一形式。

**知识点8.3　角动量守恒与应力对称（三维）**

定义：角动量守恒要求应力的力矩为零。

推导：对物质体角动量定理
$$\frac{\mathrm D}{\mathrm Dt}\int_V\mathbf r\times\rho\mathbf v\,\mathrm dV=\int_V\mathbf r\times\rho\mathbf b\,\mathrm dV+\oint_S\mathbf r\times(\boldsymbol\sigma\mathbf n)\mathrm dS.$$
利用运动方程，左端可化为 $\int_V\mathbf r\times\rho\frac{\mathrm D\mathbf v}{\mathrm Dt}\mathrm dV$ 与角速度无关项。面力项经散度定理展开，角动量方程最终给出 $\epsilon_{ijk}\sigma_{jk}=0$，即 $\sigma_{ij}=\sigma_{ji}$。

物理意义：三维中应力张量对称由角动量守恒保证，无体偶矩时成立。

极限特例：有体偶矩（极性介质）时不再对称。

典型应用：液晶、颗粒介质（微极理论）。

**知识点8.4　能量守恒与能量方程**

定义：总能量 = 动能 + 内能；外力做功 + 热流 = 能量变化。

推导：动能方程由动量方程点乘 $\mathbf v$：
$$\rho\frac{\mathrm D}{\mathrm Dt}\left(\frac12 v^2\right)=\frac{\partial(v_i\sigma_{ij})}{\partial x_j}-\sigma_{ij}\frac{\partial v_i}{\partial x_j}+\rho b_iv_i.$$
总能量方程
$$\rho\frac{\mathrm D}{\mathrm Dt}\left(e+\frac12v^2\right)=\frac{\partial(v_i\sigma_{ij})}{\partial x_j}+\rho b_iv_i-\frac{\partial q_j}{\partial x_j}+\rho r,$$
其中 $e$ 比内能，$q_j$ 热流，$r$ 体积热源。两式相减得内能方程
$$\rho\frac{\mathrm De}{\mathrm Dt}=\sigma_{ij}D_{ij}-\frac{\partial q_j}{\partial x_j}+\rho r,\qquad D_{ij}=\frac12\left(\frac{\partial v_i}{\partial x_j}+\frac{\partial v_j}{\partial x_i}\right).$$

物理意义：应力在变形率上做功（应力功率）转化为内能，热流散度与热源是其余项。

极限特例：无粘无热导时可逆绝热。

典型应用：传热与流体耦合、燃烧。

**知识点8.5　本构方程与牛顿流体**

定义：本构关系联系应力与运动学量；斯托克斯流体假设：应力与应变率线性、各向同性。

推导牛顿流体应力：
$$\sigma_{ij}=-p\delta_{ij}+\tau_{ij},\qquad \tau_{ij}=\lambda D_{kk}\delta_{ij}+2\mu D_{ij}.$$
由各向同性张量表示定理，二阶对称张量只能是 $\delta_{ij}$ 与 $D_{ij}$ 的组合。斯托克斯假设 $\lambda=-\frac23\mu$，则
$$\tau_{ij}=2\mu\left(D_{ij}-\frac13D_{kk}\delta_{ij}\right),\qquad \sigma_{ij}=-p\delta_{ij}+2\mu D_{ij}-\frac23\mu D_{kk}\delta_{ij}.$$

物理意义：切应力正比于应变率，压力为各向同性部分。

极限特例：不可压 $D_{kk}=0$ 时 $\sigma_{ij}=-p\delta_{ij}+2\mu D_{ij}$。

典型应用：水、空气等牛顿流体。

**知识点8.6　Navier-Stokes方程完整推导**

定义：牛顿流体动量方程。

推导：将本构代入柯西运动方程。以不可压为例，
$$\frac{\partial\sigma_{ij}}{\partial x_j}=-\frac{\partial p}{\partial x_i}+\mu\frac{\partial}{\partial x_j}\left(\frac{\partial v_i}{\partial x_j}+\frac{\partial v_j}{\partial x_i}\right).$$
因不可压 $\partial v_j/\partial x_j=0$ 且
$$\frac{\partial}{\partial x_j}\frac{\partial v_j}{\partial x_i}=\frac{\partial}{\partial x_i}(\nabla\cdot\mathbf v)=0,$$
故
$$\rho\left(\frac{\partial\mathbf v}{\partial t}+(\mathbf v\cdot\nabla)\mathbf v\right)=-\nabla p+\mu\nabla^2\mathbf v+\rho\mathbf b.$$
可压形式：
$$\rho\frac{\mathrm D\mathbf v}{\mathrm Dt}=-\nabla p+\nabla(\lambda\nabla\cdot\mathbf v)+\nabla\cdot(2\mu\mathbf D)+\rho\mathbf b.$$

物理意义：惯性、压力、粘性、体积力平衡；粘性项为拉普拉斯型扩散。

极限特例：$\mu=0$ 得欧拉方程；$\mathrm{Re}\ll1$ 得斯托克斯流。

典型应用：CFD、边界层、湍流。

**知识点8.7　无量纲化与雷诺数**

定义：选取特征长度 $L$、速度 $U$、密度 $\rho$，无量纲量 $\mathbf x^*=\mathbf x/L$ 等。

推导：不可压 N-S 代入
$$\rho\frac{U^2}{L}\frac{\mathrm D\mathbf v^*}{\mathrm Dt^*}=-\frac{1}{L}\nabla^*p+\mu\frac{U}{L^2}\nabla^{*2}\mathbf v^*,$$
两边乘 $L/(\rho U^2)$：
$$\frac{\mathrm D\mathbf v^*}{\mathrm Dt^*}=-\nabla^*p^*+\frac{\mu}{\rho UL}\nabla^{*2}\mathbf v^*,\qquad \mathrm{Re}=\frac{\rho UL}{\mu}.$$

物理意义：$\mathrm{Re}$ 度量惯性力与粘性力之比，决定流态。

极限特例：$\mathrm{Re}\to0$ 粘性主导（蠕动流）；$\mathrm{Re}\to\infty$ 惯性主导（边界层、湍流）。

典型应用：风洞相似准则、管道流动转捩。

**知识点8.8　伯努利方程与涡量方程**

定义：无粘、定常、不可压、沿流线。

推导伯努利：欧拉方程
$$\rho(\mathbf v\cdot\nabla)\mathbf v=-\nabla p+\rho\mathbf b.$$
用 $(\mathbf v\cdot\nabla)\mathbf v=\nabla(\frac12 v^2)-\mathbf v\times(\nabla\times\mathbf v)$，沿流线 $\mathrm d\mathbf r\parallel\mathbf v$ 点乘，$\mathbf v\times\boldsymbol\omega\cdot\mathbf v=0$，重力 $\mathbf b=-\nabla(gz)$：
$$\nabla\left(\frac12\rho v^2+p+\rho gz\right)=0\ \Rightarrow\ \frac12\rho v^2+p+\rho gz=\text{常量}.$$

推导涡量方程：取动量方程旋度，无体积力、
$$\frac{\partial\boldsymbol\omega}{\partial t}+\nabla\times(\boldsymbol\omega\times\mathbf v)=\nu\nabla^2\boldsymbol\omega,$$
用恒等式 $\nabla\times(\boldsymbol\omega\times\mathbf v)=(\mathbf v\cdot\nabla)\boldsymbol\omega-(\boldsymbol\omega\cdot\nabla)\mathbf v+\boldsymbol\omega(\nabla\cdot\mathbf v)$，不可压得
$$\frac{\mathrm D\boldsymbol\omega}{\mathrm Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf v+\nu\nabla^2\boldsymbol\omega.$$

物理意义：涡线随流体拉伸增强（涡拉伸），粘性使其扩散耗散。

极限特例：二维不可压无粘时涡量物质守恒。

典型应用：机翼尾涡、龙卷风。

### 第9章 弹性本构与弹性力学基本方程

**知识点9.1　应变能密度与超弹性**

定义：应变能密度 $W(\boldsymbol\varepsilon)$ 使 $\sigma_{ij}=\partial W/\partial\varepsilon_{ij}$。

推导：对弹性体，外力功以应变能形式储存，可逆。由热力学，
$$\mathrm dW=\sigma_{ij}\,\mathrm d\varepsilon_{ij}.$$
若 $W$ 存在，则
$$\sigma_{ij}=\frac{\partial W}{\partial\varepsilon_{ij}}.$$
对线性弹性 $W=\frac12 C_{ijkl}\varepsilon_{ij}\varepsilon_{kl}$，得 $\sigma_{ij}=C_{ijkl}\varepsilon_{kl}$。

物理意义：超弹性材料存在势函数，保证路径无关与 Maxwell 互等。

极限特例：$W=\frac12 E\varepsilon^2$ 一维弹簧。

典型应用：橡胶、软组织本构。

**知识点9.2　广义胡克定律与弹性常数**

定义：线性各向同性弹性
$$\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij}.$$

推导：由各向同性张量表示，$C_{ijkl}=\lambda\delta_{ij}\delta_{kl}+\mu(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk})$。代入得上式。反解：
$$\varepsilon_{ij}=\frac{1+\nu}{E}\sigma_{ij}-\frac{\nu}{E}\sigma_{kk}\delta_{ij}.$$
Lamé 常数与工程常数关系：
$$\lambda=\frac{E\nu}{(1+\nu)(1-2\nu)},\qquad \mu=G=\frac{E}{2(1+\nu)}.$$
推导：由 $\sigma_{kk}=(3\lambda+2\mu)\varepsilon_{kk}$，及单轴拉伸 $\sigma_{11}=\sigma$，
$$\varepsilon_{11}=\frac{\lambda+\mu}{\mu(3\lambda+2\mu)}\sigma=\frac{1}{E}\sigma,$$
$$\varepsilon_{22}=\varepsilon_{33}=-\frac{\lambda}{2\mu(3\lambda+2\mu)}\sigma=-\frac{\nu}{E}\sigma.$$
联立解出 $E,\nu$ 与 $\lambda,\mu$ 的关系。

物理意义：两个独立弹性常数描述各向同性线弹性。

极限特例：不可压 $\nu\to1/2$ 时 $\lambda\to\infty$。

典型应用：材料参数标定、有限元刚度阵。

**知识点9.3　弹性力学基本方程（15方程）**

定义：未知量 $u_i,\varepsilon_{ij},\sigma_{ij}$ 共 15 个。

推导：由
- 平衡 $\partial_j\sigma_{ij}+\rho b_i=0$：3 个；
- 几何 $\varepsilon_{ij}=\frac12(\partial_iu_j+\partial_ju_i)$：6 个；
- 本构 $\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij}$：6 个。
合 15 个方程。

位移法（Lamé-Navier）推导：将本构与几何代入平衡，
$$\partial_j[\lambda\partial_ku_k\delta_{ij}+\mu(\partial_iu_j+\partial_ju_i)]+\rho b_i=0,$$
$$\Rightarrow(\lambda+\mu)\partial_i(\partial_ku_k)+\mu\partial_j\partial_ju_i+\rho b_i=0,$$
即
$$(\lambda+\mu)\nabla(\nabla\cdot\mathbf u)+\mu\nabla^2\mathbf u+\rho\mathbf b=\mathbf 0.$$
矢量形式 $(\lambda+2\mu)\nabla(\nabla\cdot\mathbf u)-\mu\nabla\times(\nabla\times\mathbf u)+\rho\mathbf b=\mathbf 0$。

物理意义：弹性静力学是椭圆型边值问题。

极限特例：无体力时位移满足双调和方程（Beltrami-Michell）。

典型应用：有限元、边界元。

**知识点9.4　平面应力与平面应变**

定义：平面应力 $\sigma_{33}=\sigma_{13}=\sigma_{23}=0$；平面应变 $\varepsilon_{33}=\varepsilon_{13}=\varepsilon_{23}=0$。

推导平面应力本构：
$$\varepsilon_{11}=\frac1E(\sigma_{11}-\nu\sigma_{22}),\quad\varepsilon_{22}=\frac1E(\sigma_{22}-\nu\sigma_{11}),\quad\varepsilon_{12}=\frac{1+\nu}{E}\sigma_{12}.$$
平面应变：由 $\varepsilon_{33}=0$ 得 $\sigma_{33}=\nu(\sigma_{11}+\sigma_{22})$，本构为
$$\sigma_{11}=\frac{E}{(1+\nu)(1-2\nu)}[(1-\nu)\varepsilon_{11}+\nu\varepsilon_{22}],$$
$$\sigma_{22}=\frac{E}{(1+\nu)(1-2\nu)}[\nu\varepsilon_{11}+(1-\nu)\varepsilon_{22}],\quad\sigma_{12}=2\mu\varepsilon_{12}.$$

物理意义：薄板近似平面应力，长柱近似平面应变。

极限特例：$\nu=0$ 时两者相同。

典型应用：板壳、隧道。

**知识点9.5　应力函数与艾里解**

定义：艾里应力函数 $\phi(x,y)$ 使
$$\sigma_{11}=\frac{\partial^2\phi}{\partial y^2},\quad\sigma_{22}=\frac{\partial^2\phi}{\partial x^2},\quad\sigma_{12}=-\frac{\partial^2\phi}{\partial x\partial y}.$$

推导：无体力平面问题平衡 $\partial_1\sigma_{11}+\partial_2\sigma_{12}=0,\ \partial_1\sigma_{12}+\partial_2\sigma_{22}=0$ 自动满足。代入协调方程
$$\nabla^4\phi=\frac{\partial^4\phi}{\partial x^4}+2\frac{\partial^4\phi}{\partial x^2\partial y^2}+\frac{\partial^4\phi}{\partial y^4}=0.$$
即双调和方程。

物理意义：把 3 个未知应力化为 1 个标量函数，自动满足平衡。

极限特例：多项式解给出均匀应力。

典型应用：孔边应力集中（Kirsch 解）。

**知识点9.6　弹性动力学与波动方程**

定义：弹性动力学方程
$$\rho\ddot u_i=\partial_j\sigma_{ij}+\rho b_i.$$

推导：代入本构与几何，无体力得
$$\rho\ddot{\mathbf u}=(\lambda+\mu)\nabla(\nabla\cdot\mathbf u)+\mu\nabla^2\mathbf u.$$
用 Helmholtz 分解 $\mathbf u=\nabla\varphi+\nabla\times\boldsymbol\psi$，$\nabla\cdot\boldsymbol\psi=0$。纵波（膨胀波）：
$$\nabla\left[\rho\ddot\varphi-(\lambda+2\mu)\nabla^2\varphi\right]=0\ \Rightarrow\ \ddot\varphi=c_p^2\nabla^2\varphi,$$
$$c_p=\sqrt{\frac{\lambda+2\mu}{\rho}}.$$
横波（剪切波）：
$$\ddot{\boldsymbol\psi}=c_s^2\nabla^2\boldsymbol\psi,\qquad c_s=\sqrt{\frac{\mu}{\rho}}.$$

物理意义：弹性体支持两种体波，$c_p>c_s$。

极限特例：$\mu=0$（理想流体）时 $c_s=0$。

典型应用：地震波、超声无损检测。

**知识点9.7　能量原理与唯一性定理**

定义：总势能 $\Pi=\int_V W\,\mathrm dV-\int_V\mathbf b\cdot\mathbf u\,\mathrm dV-\int_{S_t}\bar{\mathbf t}\cdot\mathbf u\,\mathrm dS$。

推导最小势能原理：取容许位移变分 $\delta\mathbf u$，在 $S_u$ 上为零。一阶变分
$$\delta\Pi=\int_V\sigma_{ij}\delta\varepsilon_{ij}\mathrm dV-\int_Vb_i\delta u_i\mathrm dV-\int_{S_t}\bar t_i\delta u_i\mathrm dS.$$
由高斯散度定理
$$\int_V\sigma_{ij}\delta\varepsilon_{ij}\mathrm dV=\int_V\sigma_{ij}\partial_j(\delta u_i)\mathrm dV=\int_{S_t}\sigma_{ij}n_j\delta u_i\mathrm dS-\int_V\partial_j\sigma_{ij}\delta u_i\mathrm dV.$$
代入平衡与面力条件得 $\delta\Pi=0$。

物理意义：真实位移使总势能取驻值（稳定平衡时取极小）。

极限特例：线性弹性时 $\Pi$ 为凸，极小唯一。

典型应用：瑞利-里兹法、有限元。

**知识点9.8　功的互等定理（Betti-Maxwell）**

定义：线弹性体两种载荷状态 $(1),(2)$。

推导：由 $\sigma_{ij}^{(2)}\varepsilon_{ij}^{(1)}=\sigma_{ij}^{(1)}\varepsilon_{ij}^{(2)}$（弹性张量对称），积分并利用平衡得
$$\int_V b_i^{(1)}u_i^{(2)}\mathrm dV+\int_{S}\bar t_i^{(1)}u_i^{(2)}\mathrm dS=\int_V b_i^{(2)}u_i^{(1)}\mathrm dV+\int_S\bar t_i^{(2)}u_i^{(1)}\mathrm dS.$$
即第一组力在第二组位移上做功等于第二组力在第一组位移上做功。

物理意义：柔度矩阵对称 $K_{ij}=K_{ji}$，是结构分析基础。

极限特例：单点力给出位移互等 $u_{ij}=u_{ji}$。

典型应用：影响线、位移法。

### 第10章 能量原理与变分法

**知识点10.1　虚功原理与虚位移**

定义：对弹性体，虚位移 $\delta\mathbf u$ 满足位移边界条件，虚功原理
$$\int_V\sigma_{ij}\delta\varepsilon_{ij}\,\mathrm dV=\int_Vb_i\delta u_i\,\mathrm dV+\int_{S_t}\bar t_i\delta u_i\,\mathrm dS.$$

推导：由平衡方程 $-\partial_j\sigma_{ij}=\rho b_i$（静力），乘 $\delta u_i$ 积分，
$$\int_V(\partial_j\sigma_{ij}+b_i)\delta u_i\,\mathrm dV=0.$$
分部积分
$$\int_V\partial_j\sigma_{ij}\delta u_i\mathrm dV=\int_S\sigma_{ij}n_j\delta u_i\mathrm dS-\int_V\sigma_{ij}\partial_j\delta u_i\mathrm dV.$$
在 $S_u$ 上 $\delta u_i=0$，在 $S_t$ 上 $\sigma_{ij}n_j=\bar t_i$，并由对称性 $\sigma_{ij}\partial_j\delta u_i=\sigma_{ij}\delta\varepsilon_{ij}$，得虚功原理。

物理意义：平衡的弱形式，是有限元离散的出发点。

极限特例：刚体虚位移退化为刚体平衡。

典型应用：有限元单元刚度方程。

**知识点10.2　最小势能原理与稳定平衡**

定义：总势能 $\Pi[\mathbf u]$。

推导：见 9.7。进一步，二阶变分
$$\delta^2\Pi=\int_V C_{ijkl}\delta\varepsilon_{ij}\delta\varepsilon_{kl}\,\mathrm dV=\int_V 2\delta W\,\mathrm dV>0$$
（正定弹性张量），故驻点为极小，平衡稳定。

物理意义：能量极小对应稳定平衡；非凸能量可能出现失稳与相变。

极限特例：$W$ 非凸时极小值多解。

典型应用：屈曲分析、相场法。

**知识点10.3　余能原理与卡氏第一定理**

定义：余能 $\Pi^*=\int_V W^*\,\mathrm dV-\int_{S_u}\bar{\mathbf u}\cdot\mathbf t\,\mathrm dS$，$W^*=\frac12 C^{-1}_{ijkl}\sigma_{ij}\sigma_{kl}$（线弹性 $W^*=W$）。

推导：对静力容许应力场变分 $\delta\boldsymbol\sigma$ 满足平衡，$\delta\Pi^*=0$ 给出协调方程。即最小余能原理等价于几何协调。

推导卡氏第一定理：设 $\Pi$ 为给定广义位移 $\Delta_i$ 的函数，则广义力
$$P_i=\frac{\partial\Pi}{\partial\Delta_i}.$$
证明：$\Pi=U(\Delta)-\sum_i P_i\Delta_i$，在平衡时 $\partial U/\partial\Delta_i=P_i$，而 $\partial\Pi/\partial\Delta_i=0$ 并非本式；更准确地，应变能 $U$ 对位移导数等于力。

物理意义：能量对位移的导数给出力；对力的导数给出位移。

极限特例：线性弹簧 $U=\frac12k\Delta^2$，$P=k\Delta$。

典型应用：桁架位移求解。

**知识点10.4　卡氏第二定理与单位载荷法**

定义：应变能 $U(P_1,\dots,P_n)$ 对广义力 $P_i$ 的偏导等于对应位移
$$\Delta_i=\frac{\partial U}{\partial P_i}.$$

推导：由余能 $\Pi^*=\int A(\sigma)\mathrm dV-\sum\Delta_iP_i$，$\partial U^*/\partial P_i=\Delta_i$；线弹性时 $U^*=U$，得证。

物理意义：求任意点位移的通用方法。

极限特例：$P=0$ 处需施加虚拟单位力。

典型应用：莫尔积分、桁架/刚架位移。

**知识点10.5　瑞利-里兹法**

定义：取试解 $\mathbf u\approx\sum_{m=1}^{N}c_m\boldsymbol\phi_m$，代入泛函 $\Pi[c_m]$，令 $\partial\Pi/\partial c_m=0$。

推导：$\delta\Pi=\sum_m\frac{\partial\Pi}{\partial c_m}\delta c_m=0$，$\delta c_m$ 独立，得
$$\frac{\partial\Pi}{\partial c_m}=0,\quad m=1,\dots,N.$$
线性弹性时得代数方程组 $\mathbf K\mathbf c=\mathbf f$，其中
$$K_{mn}=\int_V C_{ijkl}\partial_j\phi_m^i\,\partial_l\phi_n^k\,\mathrm dV,\qquad f_m=\int_V b_i\phi_m^i\mathrm dV+\int_{S_t}\bar t_i\phi_m^i\mathrm dS.$$

物理意义：泛函驻值条件化为有限维代数问题，是有限元与数值变分的基础。

极限特例：$N\to\infty$ 且 $\phi_m$ 完备时收敛到精确解；$N=1$ 给近似下界（势能）。

典型应用：梁振动频率估计。

**知识点10.6　伽辽金法与加权残值**

定义：取试解 $\mathbf u_N=\sum c_m\boldsymbol\phi_m$，要求残值 $\mathbf R=\mathcal L\mathbf u_N-\mathbf f$ 与试探函数正交：
$$\int_V\mathbf R\cdot\boldsymbol\phi_m\,\mathrm dV=0.$$

推导：弱形式（虚功）用同一组 $\boldsymbol\phi_m$ 作权函数即得伽辽金法。对线性问题，伽辽金法与瑞利-里兹法在本构对称时一致；对非自伴问题伽辽金法仍可用。

物理意义：投影到试探子空间，误差正交意义下最小。

极限特例：配点法、最小二乘法是不同权函数选择。

典型应用：流体有限元、谱方法。

**知识点10.7　哈密顿原理用于连续体**

定义：拉格朗日密度 $\mathcal L=T_{\text{密}}-V_{\text{密}}$。

推导：作用量
$$S=\int_{t_1}^{t_2}\int_V\mathcal L\,\mathrm dV\,\mathrm dt.$$
对场变量 $\phi(x,t)$ 变分，$\delta S=0$ 给出
$$\frac{\partial}{\partial t}\frac{\partial\mathcal L}{\partial\dot\phi}+\frac{\partial}{\partial x_i}\frac{\partial\mathcal L}{\partial(\partial_i\phi)}-\frac{\partial\mathcal L}{\partial\phi}=0.$$
以弹性杆为例 $\mathcal L=\frac12\rho\dot u^2-\frac12 E(u')^2$，得波动方程
$$\rho\ddot u-Eu''=0.$$

物理意义：连续系统动力学由场论型欧拉-拉格朗日方程支配。

极限特例：线性弹性退化为波动方程。

典型应用：梁板壳动力学、场论。

### 第11章 线性振动系统

**知识点11.1　单自由度自由振动**

定义：$m\ddot x+kx=0$。

推导：令 $\omega_n=\sqrt{k/m}$，特征方程 $\lambda^2+\omega_n^2=0$，$\lambda=\pm\mathrm i\omega_n$，通解
$$x(t)=A\cos\omega_n t+B\sin\omega_n t=C\cos(\omega_n t-\varphi),$$
$$C=\sqrt{A^2+B^2},\quad\tan\varphi=B/A.$$
由初值 $x(0)=x_0,\dot x(0)=v_0$：$A=x_0,B=v_0/\omega_n$。

物理意义：能量在动能与势能间周期交换，总能量 $E=\frac12kC^2$ 守恒。

极限特例：$k\to0$ 时 $\omega_n\to0$，自由粒子。

典型应用：简谐振子、单摆小角近似。

**知识点11.2　阻尼自由振动**

定义：$m\ddot x+c\dot x+kx=0$。

推导：设 $\zeta=\frac{c}{2\sqrt{km}}$，$\omega_n=\sqrt{k/m}$，特征方程
$$\lambda_{1,2}=-\zeta\omega_n\pm\omega_n\sqrt{\zeta^2-1}.$$
欠阻尼 $\zeta<1$：$\omega_d=\omega_n\sqrt{1-\zeta^2}$，
$$x=e^{-\zeta\omega_n t}(A\cos\omega_d t+B\sin\omega_d t).$$
对数衰减率
$$\delta=\ln\frac{x(t)}{x(t+T_d)}=2\pi\zeta/\sqrt{1-\zeta^2}.$$
过阻尼 $\zeta>1$：两个实负根，无振荡。临界 $\zeta=1$：重根 $x=(A+Bt)e^{-\omega_n t}$。

物理意义：阻尼耗散能量，振幅指数衰减。

极限特例：$\zeta=1$ 最快无振荡趋于平衡。

典型应用：车辆悬架、隔振设计。

**知识点11.3　受迫振动与共振**

定义：$m\ddot x+c\dot x+kx=F_0\cos\omega t$。

推导：稳态解 $x_p=X\cos(\omega t-\phi)$，代入得
$$X=\frac{F_0/k}{\sqrt{(1-r^2)^2+(2\zeta r)^2}},\qquad r=\frac{\omega}{\omega_n},$$
$$\tan\phi=\frac{2\zeta r}{1-r^2}.$$
位移放大因子 $\beta=Xk/F_0$，最大响应满足
$$r_{\text{res}}=\sqrt{1-2\zeta^2}\ (\zeta<1/\sqrt2).$$

物理意义：共振附近相位滞后 $90^\circ$；阻尼限制共振峰高。

极限特例：$\zeta\to0$ 且 $r=1$ 时振幅发散（能量持续注入）。

典型应用：地震响应谱、旋转机械临界转速。

**知识点11.4　多自由度系统与特征值问题**

定义：$\mathbf M\ddot{\mathbf q}+\mathbf K\mathbf q=\mathbf 0$。

推导：设 $\mathbf q=\boldsymbol\phi e^{\mathrm i\omega t}$，得
$$(\mathbf K-\omega^2\mathbf M)\boldsymbol\phi=\mathbf 0.$$
非零解要求
$$\det(\mathbf K-\omega^2\mathbf M)=0,$$
得到 $n$ 个特征值 $\omega_1^2\le\dots\le\omega_n^2$。$\mathbf M,\mathbf K$ 对称正定时特征值实正，特征向量关于 $\mathbf M$ 正交：
$$\boldsymbol\phi_i^{\mathsf T}\mathbf M\boldsymbol\phi_j=\delta_{ij}m_i,\qquad \boldsymbol\phi_i^{\mathsf T}\mathbf K\boldsymbol\phi_j=\delta_{ij}k_i.$$
瑞利商
$$R(\mathbf v)=\frac{\mathbf v^{\mathsf T}\mathbf K\mathbf v}{\mathbf v^{\mathsf T}\mathbf M\mathbf v},\qquad \omega_1^2=\min_{\mathbf v}R(\mathbf v).$$

物理意义：系统有 $n$ 个固有频率与主振型，正交性解耦方程。

极限特例：对角质量矩阵与刚度矩阵时频率 $=\sqrt{k_i/m_i}$。

典型应用：建筑模态分析、转子动力学。

**知识点11.5　模态叠加法**

定义：$\mathbf q(t)=\sum_i\boldsymbol\phi_i\,\eta_i(t)$。

推导：代入运动方程并左乘 $\boldsymbol\phi_j^{\mathsf T}$，利用正交性得解耦模态方程
$$m_j\ddot\eta_j+k_j\eta_j=\boldsymbol\phi_j^{\mathsf T}\mathbf F(t).$$
每个模态是单自由度系统，解出后叠加。

物理意义：任意响应是各阶主振型的线性组合。

极限特例：初始条件只激发某阶时系统做纯模态振动。

典型应用：振动测试、模型降阶。

**知识点11.6　连续系统振动（弦、梁）**

定义：弦横向位移 $y(x,t)$，张力 $T$，线密度 $\rho$。

推导弦方程：微元受力竖直分量
$$T\sin\theta_2-T\sin\theta_1\approx T\left.\frac{\partial y}{\partial x}\right|_{x_2}-T\left.\frac{\partial y}{\partial x}\right|_{x_1}\approx T\frac{\partial^2y}{\partial x^2}\mathrm dx,$$
由牛顿第二定律
$$\rho\,\mathrm dx\frac{\partial^2y}{\partial t^2}=T\frac{\partial^2y}{\partial x^2}\mathrm dx\ \Rightarrow\ \frac{\partial^2y}{\partial t^2}=c^2\frac{\partial^2y}{\partial x^2},\quad c=\sqrt{T/\rho}.$$
分离变量 $y=X(x)U(t)$，$X''/X=U''/(c^2U)=-\kappa^2$，得模态 $X_n=\sin(n\pi x/L)$，$\omega_n=n\pi c/L$。

梁的欧拉-伯努利方程：
$$\rho A\frac{\partial^2w}{\partial t^2}+EI\frac{\partial^4w}{\partial x^4}=0,$$
特征频率 $\omega_n=(\beta_n L)^2\sqrt{EI/(\rho A L^4)}$。

物理意义：连续系统有无穷多模态，频率与边界条件有关。

极限特例：固定-自由梁 $\beta_nL$ 为超越方程根。

典型应用：弦乐器、桥梁风致振动。

**知识点11.7　隔振与动力吸振器**

定义：基础激励 $y(t)$，质量块运动 $x(t)$，传递率 $T_r$。

推导：隔振系统 $m\ddot x+c(\dot x-\dot y)+k(x-y)=0$，设 $y=Y e^{\mathrm i\omega t}$，稳态位移传递率
$$T_r=\frac{|X|}{|Y|}=\sqrt{\frac{1+(2\zeta r)^2}{(1-r^2)^2+(2\zeta r)^2}}.$$
隔振条件 $r>\sqrt2$。

动力吸振器：主质量 $m_1$、附加 $m_2$，调谐 $\omega_2=\sqrt{k_2/m_2}=\omega_1$，使主质量在 $\omega_1$ 处响应为零。

物理意义：附加共振吸收能量。

极限特例：无阻尼吸振器在主质量共振处幅值为零。

典型应用：TMD 高塔减振。

### 第12章 稳定性理论基础

**知识点12.1　平衡点与李雅普诺夫稳定性定义**

定义：系统 $\dot{\mathbf x}=\mathbf f(\mathbf x)$，平衡点 $\mathbf f(\mathbf x_e)=\mathbf 0$。

- 稳定：$\forall\varepsilon>0,\exists\delta>0$，$\|\mathbf x(0)-\mathbf x_e\|<\delta\Rightarrow\|\mathbf x(t)-\mathbf x_e\|<\varepsilon,\forall t\ge0$。
- 渐近稳定：稳定且 $\mathbf x(t)\to\mathbf x_e$。
- 不稳定：否则。

推导：稳定性是解对初扰动的连续依赖。

物理意义：局部性质，与吸引域有关。

极限特例：线性系统全局稳定与局部一致。

典型应用：飞行器平衡状态、控制回路。

**知识点12.2　李雅普诺夫第二方法**

定义：设 $V(\mathbf x)$ 在 $\mathbf x_e$ 邻域正定，$\dot V=\nabla V\cdot\mathbf f(\mathbf x)$。

定理（充分条件）：若 $\dot V\le0$ 则稳定；若 $\dot V<0$（负定）则渐近稳定；若 $\dot V>0$ 则不稳定。

推导：对任意轨线
$$\frac{\mathrm d}{\mathrm dt}V(\mathbf x(t))=\nabla V\cdot\dot{\mathbf x}=\nabla V\cdot\mathbf f(\mathbf x)=\dot V.$$
若 $\dot V\le0$，则 $V(\mathbf x(t))\le V(\mathbf x(0))$。由 $V$ 正定连续，存在比较函数 $\alpha$ 使 $\alpha(\|\mathbf x-\mathbf x_e\|)\le V$，故 $\|\mathbf x-\mathbf x_e\|$ 有界，稳定。负定时 $V$ 单调减且有下界，极限存在，轨线趋于 $\dot V=0$ 集合；结合径向无界得渐近稳定。

物理意义：用广义能量函数判断稳定性，无需解方程。

极限特例：$\dot V=0$ 仅说明稳定，可能极限环。

典型应用：非线性控制、拉塞尔不变集。

**知识点12.3　线性化稳定性与特征值判据**

定义：$\mathbf A=\partial\mathbf f/\partial\mathbf x|_{\mathbf x_e}$。

推导：设 $\mathbf x=\mathbf x_e+\delta\mathbf x$，泰勒展开
$$\dot{\delta\mathbf x}=\mathbf A\,\delta\mathbf x+O(\|\delta\mathbf x\|^2).$$
线性系统解 $\delta\mathbf x(t)=e^{\mathbf A t}\delta\mathbf x(0)$。若 $\mathbf A$ 所有特征值实部 $<0$，则渐近稳定；若存在实部 $>0$ 则不稳定；若临界（实部为零）需非线性项判断。

物理意义：双曲平衡点的稳定性由线性化决定（Hartman-Grobman）。

极限特例：中心 $\pm\mathrm i\omega$ 的线性化不能判断（可能稳定或不稳定）。

典型应用：极点配置、结构稳定。

**知识点12.4　劳斯-赫尔维茨判据**

定义：特征多项式 $p(s)=a_0s^n+a_1s^{n-1}+\dots+a_n$，$a_0>0$。

推导：劳斯表递推
$$b_1=\frac{a_1a_2-a_0a_3}{a_1},\quad b_2=\frac{a_1a_4-a_0a_5}{a_1},\dots$$
所有第一列元素 $>0$ 时系统稳定（根全在左半平面）。二阶条件 $a_1a_2>a_0a_3$；三阶 $a_1a_2>a_0a_3$。

物理意义：无需解根即可判稳。

极限特例：$n=2$ 时 $a_1>0,a_2>0$ 即稳定。

典型应用：控制器参数范围设计。

**知识点12.5　中心流形与分岔初步**

定义：平衡点处特征值分实部为零（临界）与负实部。

推导：将系统分块为临界子空间 $\mathbf y$、稳定子空间 $\mathbf z$：
$$\dot{\mathbf y}=\mathbf A\mathbf y+\mathbf g(\mathbf y,\mathbf z),\quad \dot{\mathbf z}=\mathbf B\mathbf z+\mathbf h(\mathbf y,\mathbf z),$$
$\mathbf A$ 特征值零实部，$\mathbf B$ 负实部。存在中心流形 $\mathbf z=\boldsymbol\pi(\mathbf y)$ 满足
$$\frac{\partial\boldsymbol\pi}{\partial\mathbf y}[\mathbf A\mathbf y+\mathbf g(\mathbf y,\boldsymbol\pi)]=\mathbf B\boldsymbol\pi+\mathbf h(\mathbf y,\boldsymbol\pi).$$
稳定性由降维系统 $\dot{\mathbf y}=\mathbf A\mathbf y+\mathbf g(\mathbf y,\boldsymbol\pi(\mathbf y))$ 决定。

物理意义：临界模态主导慢时间尺度动力学。

极限特例：无稳定子空间时系统直接是中心流形。

典型应用：Hopf 分岔、机翼颤振。

**知识点12.6　结构稳定性与李雅普诺夫方程**

定义：线性系统 $\dot{\mathbf x}=\mathbf A\mathbf x$。

推导：选 $V=\mathbf x^{\mathsf T}\mathbf P\mathbf x$，$\mathbf P=\mathbf P^{\mathsf T}>0$，
$$\dot V=\mathbf x^{\mathsf T}(\mathbf A^{\mathsf T}\mathbf P+\mathbf P\mathbf A)\mathbf x=-\mathbf x^{\mathsf T}\mathbf Q\mathbf x.$$
给定 $\mathbf Q>0$，若李雅普诺夫方程
$$\mathbf A^{\mathsf T}\mathbf P+\mathbf P\mathbf A=-\mathbf Q$$
有唯一正定解 $\mathbf P$，则系统渐近稳定。

物理意义：把稳定性判定转化为线性矩阵方程求解。

极限特例：$\mathbf Q=\mathbf I$ 给 $\mathbf A^{\mathsf T}\mathbf P+\mathbf P\mathbf A=-\mathbf I$。

典型应用：LQR、鲁棒稳定性。

### 第13章 经典控制理论

**知识点13.1　拉普拉斯变换与传递函数**

定义：单边拉普拉斯变换 $\mathcal L\{f(t)\}=\int_0^\infty f(t)e^{-st}\mathrm dt=F(s)$。

推导：线性性质与微分性质
$$\mathcal L\{\dot f\}=sF(s)-f(0^-).$$
对线性定常系统
$$a_n y^{(n)}+\dots+a_0y=b_m u^{(m)}+\dots+b_0u,$$
零初值下取变换
$$(a_ns^n+\dots+a_0)Y(s)=(b_ms^m+\dots+b_0)U(s),$$
得传递函数
$$G(s)=\frac{Y(s)}{U(s)}=\frac{b_ms^m+\dots+b_0}{a_ns^n+\dots+a_0}.$$

物理意义：零初始条件下输出与输入拉氏变换之比，描述系统固有动态。

极限特例：$u=\delta$ 时 $Y(s)=G(s)$，即脉冲响应 $g(t)=\mathcal L^{-1}\{G\}$，且 $y=g*u$。

典型应用：频域分析与控制器设计。

**知识点13.2　典型环节与系统响应**

定义：比例 $K$、积分 $1/s$、微分 $s$、一阶 $K/(Ts+1)$、二阶 $\omega_n^2/(s^2+2\zeta\omega_ns+\omega_n^2)$。

推导一阶阶跃响应：$G=K/(Ts+1)$，$U=1/s$，
$$Y=\frac{K}{s(Ts+1)}=\frac{K}{s}-\frac{KT}{Ts+1},$$
$$y(t)=K(1-e^{-t/T}),\quad t\ge0.$$
二阶阶跃响应（欠阻尼）：
$$y(t)=1-\frac{e^{-\zeta\omega_nt}}{\sqrt{1-\zeta^2}}\sin\!\left(\omega_dt+\arccos\zeta\right).$$
上升时间、超调量 $M_p=e^{-\pi\zeta/\sqrt{1-\zeta^2}}$。

物理意义：时域指标由 $\zeta,\omega_n$ 决定。

极限特例：$\zeta=0$ 无阻尼等幅振荡。

典型应用：伺服系统整定。

**知识点13.3　系统稳定性与罗斯判据应用**

定义：闭环特征方程 $1+G(s)H(s)=0$。

推导：将开环传递函数代入得闭环多项式，用第12章劳斯判据求使系统稳定的参数范围（如增益 $K$ 上限）。

物理意义：反馈可能失稳，增益过大会引入右半平面极点。

极限特例：二阶系统恒稳定（正系数）。

典型应用：PID 增益稳定边界。

**知识点13.4　根轨迹法**

定义：闭环极点随增益 $K$ 变化的轨迹，$1+KG(s)=0$。

推导：幅值条件与相角条件
$$|KG(s)|=1,\qquad \angle G(s)=(2k+1)\pi.$$
起点为开环极点 $K=0$，终点为开环零点 $K\to\infty$。渐近线角度
$$\theta_a=\frac{(2k+1)\pi}{n-m},\qquad \sigma_a=\frac{\sum p_i-\sum z_j}{n-m}.$$

物理意义：直观显示参数对闭环极点位置的影响。

极限特例：比例反馈只能落在根轨迹上。

典型应用：增益与补偿器设计。

**知识点13.5　频率响应与Bode图**

定义：$s=\mathrm i\omega$，$G(\mathrm i\omega)$ 的幅频与相频。

推导：对正弦输入 $u=A\sin\omega t$，线性定常系统稳态输出
$$y=|G(\mathrm i\omega)|A\sin(\omega t+\angle G(\mathrm i\omega)).$$
幅值以分贝 $\mathrm{dB}=20\log_{10}|G|$，相位以度表示。

物理意义：Bode 图用对数坐标显示增益裕度、相位裕度。

极限特例：一阶环节转折频率 $\omega=1/T$ 处相位 $-45^\circ$。

典型应用：频域整定、鲁棒性分析。

**知识点13.6　奈奎斯特稳定判据**

定义：$Z=N+P$，$N$ 为开环绕 $(-1,\mathrm j0)$ 的净圈数，$P$ 为右半平面开环极点数，$Z$ 为右半平面闭环极点数。

推导：由幅角原理，对闭合围线 $\Gamma$ 有
$$N=\frac{1}{2\pi}\Delta\arg F(s),\qquad F(s)=1+G(s)H(s).$$
$\Gamma$ 包住整个右半平面（含虚轴绕行），$F$ 的零点即闭环极点，极点即开环极点，故 $Z=N+P$。稳定要求 $Z=0$，即 $N=-P$。

物理意义：用开环频响判断闭环稳定，适用于非最小相位系统。

极限特例：开环稳定 $P=0$ 时要求不包围临界点。

典型应用：时滞系统稳定分析。

**知识点13.7　PID控制器**

定义：$u(t)=K_pe+K_i\int e\,\mathrm dt+K_d\dot e$，$e=r-y$。

推导：频域
$$G_c(s)=K_p+\frac{K_i}{s}+K_ds=K_p\left(1+\frac{1}{T_is}+T_ds\right).$$
闭环二阶近似下，$K_p$ 提高带宽、减小稳态误差；$K_i$ 消除常值误差但降低相位裕度；$K_d$ 提供超前相位改善阻尼。

物理意义：P、I、D 分别调节刚度、消除偏差、增加阻尼。

极限特例：纯 P 对一阶系统有稳态误差，PI 消除之。

典型应用：飞行控制、工业过程控制。

**知识点13.8　稳态误差与系统型别**

定义：误差传递函数 $E(s)=R(s)/(1+G(s)H(s))$。

推导：终值定理
$$e_{ss}=\lim_{t\to\infty}e(t)=\lim_{s\to0}sE(s).$$
对单位阶跃、斜坡、加速度输入，稳态误差由开环传递函数的型别 $N$（原点极点数）决定：
- $N=0$：阶跃误差 $1/(1+K_p)$；
- $N=1$：阶跃零误差、斜坡误差 $1/K_v$；
- $N=2$：斜坡零误差、加速度误差 $1/K_a$。

物理意义：型别决定系统跟踪能力。

极限特例：$N\ge$ 输入阶数时无稳态误差。

典型应用：跟踪系统设计。

### 第14章 状态空间与现代控制

**知识点14.1　状态空间模型与解**

定义：$\dot{\mathbf x}=\mathbf A\mathbf x+\mathbf B\mathbf u,\ \mathbf y=\mathbf C\mathbf x+\mathbf D\mathbf u$。

推导：齐次解
$$\mathbf x(t)=e^{\mathbf A t}\mathbf x(0),\qquad e^{\mathbf A t}=\sum_{k=0}^\infty\frac{(\mathbf A t)^k}{k!}.$$
非齐次解（常数变易）
$$\mathbf x(t)=e^{\mathbf A t}\mathbf x(0)+\int_0^t e^{\mathbf A(t-\tau)}\mathbf B\mathbf u(\tau)\,\mathrm d\tau.$$
验证：求导得
$$\dot{\mathbf x}=\mathbf A e^{\mathbf A t}\mathbf x(0)+\mathbf B\mathbf u(t)+\int_0^t\mathbf A e^{\mathbf A(t-\tau)}\mathbf B\mathbf u\,\mathrm d\tau=\mathbf A\mathbf x+\mathbf B\mathbf u.$$

物理意义：状态包含系统全部记忆。

极限特例：$\mathbf A$ 可对角化时 $e^{\mathbf A t}=\mathbf V e^{\boldsymbol\Lambda t}\mathbf V^{-1}$。

典型应用：数字仿真、观测器。

**知识点14.2　能控性**

定义：存在控制使任意初态到任意终态。能控性矩阵
$$\mathcal C=\begin{pmatrix}\mathbf B&\mathbf A\mathbf B&\cdots&\mathbf A^{n-1}\mathbf B\end{pmatrix}.$$

推导（充分性）：若 $\mathrm{rank}\,\mathcal C=n$，对目标 $\mathbf x(T)=\mathbf 0$，取
$$\mathbf u(t)=-\mathbf B^{\mathsf T}e^{\mathbf A^{\mathsf T}(T-t)}\mathbf W^{-1}\mathbf x(0),$$
其中能控性格拉姆矩阵
$$\mathbf W=\int_0^T e^{\mathbf A\tau}\mathbf B\mathbf B^{\mathsf T}e^{\mathbf A^{\mathsf T}\tau}\mathrm d\tau.$$
代入状态解可验证到达零。必要性由凯莱-哈密顿定理与反证给出。

物理意义：$\mathrm{rank}\,\mathcal C=n$ 即所有状态模式可由输入激励。

极限特例：单输入系统 $\mathcal C$ 为方阵，可逆即能控。

典型应用：极点配置前提。

**知识点14.3　能观性**

定义：由输出与输入确定初态。能观性矩阵
$$\mathcal O=\begin{pmatrix}\mathbf C\\ \mathbf C\mathbf A\\ \vdots\\ \mathbf C\mathbf A^{n-1}\end{pmatrix}.$$

推导：$\mathbf y(t)=\mathbf C e^{\mathbf A t}\mathbf x(0)$（$\mathbf u=0$）。求各阶导数：
$$\mathbf y^{(k)}(0)=\mathbf C\mathbf A^k\mathbf x(0),$$
联立得 $\mathcal O\mathbf x(0)=[\mathbf y(0);\dot{\mathbf y}(0);\dots]$。有唯一解要求 $\mathrm{rank}\,\mathcal O=n$。

物理意义：能观性对偶于能控性（$\mathbf A\leftrightarrow\mathbf A^{\mathsf T},\mathbf B\leftrightarrow\mathbf C^{\mathsf T}$）。

极限特例：输出直接等于状态时能观。

典型应用：状态观测器存在性。

**知识点14.4　状态反馈与极点配置**

定义：$\mathbf u=-\mathbf K\mathbf x+\mathbf r$，闭环
$$\dot{\mathbf x}=(\mathbf A-\mathbf B\mathbf K)\mathbf x+\mathbf B\mathbf r.$$

推导：闭环特征多项式 $\det(s\mathbf I-\mathbf A+\mathbf B\mathbf K)$。当 $(\mathbf A,\mathbf B)$ 能控时，可通过选择 $\mathbf K$ 任意配置极点。单输入时，将系统化为能控标准型，比较系数即可解出 $\mathbf K$。

物理意义：反馈改变系统固有频率与阻尼。

极限特例：全状态反馈可任意配置；输出反馈受限。

典型应用：自动驾驶、姿态控制。

**知识点14.5　状态观测器**

定义：$\dot{\hat{\mathbf x}}=\mathbf A\hat{\mathbf x}+\mathbf B\mathbf u+\mathbf L(\mathbf y-\mathbf C\hat{\mathbf x})$。

推导：误差 $\mathbf e=\mathbf x-\hat{\mathbf x}$ 满足
$$\dot{\mathbf e}=(\mathbf A-\mathbf L\mathbf C)\mathbf e.$$
故当 $(\mathbf C,\mathbf A)$ 能观时可选 $\mathbf L$ 使 $\mathbf A-\mathbf L\mathbf C$ 特征值任意配置，$\mathbf e\to\mathbf 0$。分离原理：观测器与控制器极点可独立配置。

物理意义：用输出偏差校正估计。

极限特例：$\mathbf L=\mathbf 0$ 为开环估计，误差不衰减。

典型应用：无传感器控制、卡尔曼滤波基础。

**知识点14.6　线性二次型最优控制（LQR）**

定义：最小化
$$J=\int_0^\infty(\mathbf x^{\mathsf T}\mathbf Q\mathbf x+\mathbf u^{\mathsf T}\mathbf R\mathbf u)\,\mathrm dt,\quad \mathbf Q\ge0,\mathbf R>0.$$

推导：设最优值 $V=\mathbf x^{\mathsf T}\mathbf P\mathbf x$。哈密顿-雅可比-贝尔曼方程（或直接极小化哈密顿量）给出代数黎卡提方程
$$\mathbf A^{\mathsf T}\mathbf P+\mathbf P\mathbf A-\mathbf P\mathbf B\mathbf R^{-1}\mathbf B^{\mathsf T}\mathbf P+\mathbf Q=\mathbf 0,$$
最优控制
$$\mathbf u^*=-\mathbf R^{-1}\mathbf B^{\mathsf T}\mathbf P\mathbf x=-\mathbf K\mathbf x.$$
验证：$\dot V=\mathbf x^{\mathsf T}(\mathbf A^{\mathsf T}\mathbf P+\mathbf P\mathbf A)\mathbf x+2\mathbf x^{\mathsf T}\mathbf P\mathbf B\mathbf u+\mathbf u^{\mathsf T}\mathbf R\mathbf u$，对 $\mathbf u$ 求极值 $\partial/\partial\mathbf u=2\mathbf B^{\mathsf T}\mathbf P\mathbf x+2\mathbf R\mathbf u=0$，得 $\mathbf u^*$，代回并令剩余项 $+\mathbf x^{\mathsf T}\mathbf Q\mathbf x=0$ 得黎卡提方程。

物理意义：折中状态偏差与控制代价，保证闭环稳定。

极限特例：$\mathbf R\to\infty$ 时控制增益趋于零。

典型应用：飞行控制、倒立摆。

**知识点14.7　卡尔曼滤波（连续形式）**

定义：对随机系统
$$\dot{\mathbf x}=\mathbf A\mathbf x+\mathbf w,\quad \mathbf y=\mathbf C\mathbf x+\mathbf v,$$
$\mathbf w\sim(0,\mathbf Q_f)$，$\mathbf v\sim(0,\mathbf R_f)$。

推导：估计误差协方差 $\mathbf P=E[\mathbf e\mathbf e^{\mathsf T}]$。最小方差估计的黎卡提微分方程
$$\dot{\mathbf P}=\mathbf A\mathbf P+\mathbf P\mathbf A^{\mathsf T}-\mathbf P\mathbf C^{\mathsf T}\mathbf R_f^{-1}\mathbf C\mathbf P+\mathbf Q_f,$$
卡尔曼增益 $\mathbf K_f=\mathbf P\mathbf C^{\mathsf T}\mathbf R_f^{-1}$，估计更新
$$\dot{\hat{\mathbf x}}=\mathbf A\hat{\mathbf x}+\mathbf K_f(\mathbf y-\mathbf C\hat{\mathbf x}).$$

物理意义：兼顾模型预测与测量信息，是线性高斯系统的最优估计。

极限特例：稳态时退化为代数黎卡提方程。

典型应用：惯导组合、目标跟踪。

### 第15章 最优控制

**知识点15.1　最优控制问题与变分法**

定义：系统 $\dot{\mathbf x}=\mathbf f(\mathbf x,\mathbf u,t)$，性能指标
$$J=\Phi(\mathbf x(t_f),t_f)+\int_{t_0}^{t_f}L(\mathbf x,\mathbf u,t)\,\mathrm dt.$$

推导：引入拉格朗日乘子 $\boldsymbol\lambda$，构造增广泛函
$$\bar J=\Phi+\int_{t_0}^{t_f}\left[L+\boldsymbol\lambda^{\mathsf T}(\mathbf f-\dot{\mathbf x})\right]\mathrm dt.$$
分部积分 $\int\boldsymbol\lambda^{\mathsf T}\dot{\mathbf x}\,\mathrm dt=[\boldsymbol\lambda^{\mathsf T}\mathbf x]-\int\dot{\boldsymbol\lambda}^{\mathsf T}\mathbf x\,\mathrm dt$，对 $\mathbf x,\mathbf u,\boldsymbol\lambda$ 分别变分，得
$$\dot{\mathbf x}=\frac{\partial H}{\partial\boldsymbol\lambda},\quad \dot{\boldsymbol\lambda}=-\frac{\partial H}{\partial\mathbf x},\quad \frac{\partial H}{\partial\mathbf u}=0,$$
其中哈密顿量
$$H=L+\boldsymbol\lambda^{\mathsf T}\mathbf f.$$
横截条件
$$\left.\left(\frac{\partial\Phi}{\partial\mathbf x}-\boldsymbol\lambda\right)\right|_{t_f}=0.$$

物理意义：最优控制满足正则方程与驻值条件。

极限特例：终端固定时无横截条件。

典型应用：轨迹优化。

**知识点15.2　庞特里亚金极大值原理**

定义：$H(\mathbf x,\mathbf u,\boldsymbol\lambda,t)$。

定理：最优控制使哈密顿量逐点全局最小（或最大，取决于符号约定）：
$$\mathbf u^*=\arg\min_{\mathbf u\in U}H(\mathbf x^*,\mathbf u,\boldsymbol\lambda^*,t).$$
协态方程
$$\dot{\boldsymbol\lambda}=-\frac{\partial H}{\partial\mathbf x},\qquad \dot{\mathbf x}=\frac{\partial H}{\partial\boldsymbol\lambda}.$$

推导：由变分法，对 $\mathbf u$ 的变分 $\delta J=\int(\partial H/\partial\mathbf u)\delta\mathbf u\,\mathrm dt$。当 $\mathbf u$ 无约束时 $\partial H/\partial\mathbf u=0$；当 $\mathbf u$ 有界约束时，驻值条件改为全局极小。

物理意义：把无穷维优化化为逐点代数优化。

极限特例：无约束时退化为变分法条件。

典型应用：时间最优（bang-bang）、航天器燃料最优。

**知识点15.3　哈密顿-雅可比-贝尔曼方程**

定义：最优代价函数 $V(\mathbf x,t)=\min_{\mathbf u}\left[\int_t^{t_f}L\,\mathrm d\tau+\Phi\right]$。

推导：由动态规划原理
$$V(\mathbf x,t)=\min_{\mathbf u}\left\{L(\mathbf x,\mathbf u,t)\Delta t+V(\mathbf x+\mathbf f\Delta t,t+\Delta t)\right\}.$$
泰勒展开
$$V(\mathbf x+\mathbf f\Delta t,t+\Delta t)\approx V+\nabla_{\mathbf x}V\cdot\mathbf f\Delta t+\frac{\partial V}{\partial t}\Delta t.$$
令 $\Delta t\to0$：
$$-\frac{\partial V}{\partial t}=\min_{\mathbf u}\left\{L+\nabla_{\mathbf x}V\cdot\mathbf f\right\}=H^*.$$
对 $\mathbf u$ 求极值给出最优控制。终端条件 $V(\mathbf x,t_f)=\Phi$。

物理意义：连续时间动态规划，给出全局最优反馈。

极限特例：LQR 中 $V=\mathbf x^{\mathsf T}\mathbf P\mathbf x$ 给出黎卡提方程。

典型应用：非线性最优反馈。

**知识点15.4　线性二次型与两点边值问题**

定义：线性系统 $\dot{\mathbf x}=\mathbf A\mathbf x+\mathbf B\mathbf u$，二次型指标。

推导：哈密顿量
$$H=\frac12\mathbf x^{\mathsf T}\mathbf Q\mathbf x+\frac12\mathbf u^{\mathsf T}\mathbf R\mathbf u+\boldsymbol\lambda^{\mathsf T}(\mathbf A\mathbf x+\mathbf B\mathbf u).$$
$\partial H/\partial\mathbf u=\mathbf R\mathbf u+\mathbf B^{\mathsf T}\boldsymbol\lambda=0$ 得 $\mathbf u=-\mathbf R^{-1}\mathbf B^{\mathsf T}\boldsymbol\lambda$。协态方程
$$\dot{\boldsymbol\lambda}=-\mathbf Q\mathbf x-\mathbf A^{\mathsf T}\boldsymbol\lambda.$$
与状态方程构成两点边值问题：
$$\begin{pmatrix}\dot{\mathbf x}\\ \dot{\boldsymbol\lambda}\end{pmatrix}=\begin{pmatrix}\mathbf A&-\mathbf B\mathbf R^{-1}\mathbf B^{\mathsf T}\\ -\mathbf Q&-\mathbf A^{\mathsf T}\end{pmatrix}\begin{pmatrix}\mathbf x\\ \boldsymbol\lambda\end{pmatrix}.$$
设 $\boldsymbol\lambda=\mathbf P(t)\mathbf x$，代入得黎卡提微分方程
$$-\dot{\mathbf P}=\mathbf P\mathbf A+\mathbf A^{\mathsf T}\mathbf P-\mathbf P\mathbf B\mathbf R^{-1}\mathbf B^{\mathsf T}\mathbf P+\mathbf Q.$$

物理意义：最优控制是状态的线性反馈，时变增益由黎卡提方程给出。

极限特例：无限时域稳态给出代数黎卡提方程。

典型应用：有限时间交会、制导。

**知识点15.5　时间最优控制与bang-bang**

定义：$\dot{\mathbf x}=\mathbf A\mathbf x+\mathbf B u$，$|u|\le1$，最小化 $J=\int_0^{t_f}\mathrm dt$。

推导：$H=1+\boldsymbol\lambda^{\mathsf T}(\mathbf A\mathbf x+\mathbf Bu)$ 对 $u$ 线性。极小化要求
$$u^*=-\mathrm{sgn}\left(\boldsymbol\lambda^{\mathsf T}\mathbf B\right).$$
控制幅值取边界，切换次数由切换函数零点决定。

物理意义：最速控制总是取最大或最小控制值。

极限特例：双积分系统的时间最优为有限次切换。

典型应用：电机、机械臂快速定位。

**知识点15.6　奇异最优控制**

定义：当 $\partial H/\partial\mathbf u=0$ 对某区间成立且不能解出 $\mathbf u$，称奇异弧。

推导：对 $\frac{\mathrm d}{\mathrm dt}\frac{\partial H}{\partial\mathbf u}=0$ 反复求导直到出现 $\mathbf u$，得奇异控制。其最优性需奇异弧的二阶条件（Legendre-Clebsch 推广）判断。

物理意义：奇异弧出现在约束与状态耦合退化情形。

极限特例：线性系统哈密顿量对 $u$ 线性时，$\lambda^{\mathsf T}B=0$ 给出奇异情形。

典型应用：燃油最优、奇异弧最优轨道。

### 第16章 非线性动力学与分岔

**知识点16.1　相平面与平衡点分类**

定义：二维自治系统 $\dot x=f(x,y),\dot y=g(x,y)$。

推导：消去时间 $\frac{\mathrm dy}{\mathrm dx}=\frac{g}{f}$。线性化矩阵
$$\mathbf J=\begin{pmatrix}f_x&f_y\\ g_x&g_y\end{pmatrix}.$$
特征值 $\lambda_{1,2}$ 决定类型：结点、鞍点、焦点、中心。迹 $\tau=\mathrm{tr}\mathbf J$、行列式 $\Delta=\det\mathbf J$，$\lambda=\frac{\tau\pm\sqrt{\tau^2-4\Delta}}{2}$。

物理意义：相图全局结构由平衡点与连接轨道组成。

极限特例：$\Delta<0$ 必为鞍点。

典型应用：单摆、捕食模型。

**知识点16.2　极限环**

定义：孤立周期轨道。

推导（Poincaré-Bendixson定理）：有界平面区域内若无平衡点，则存在极限环。判断稳定性用庞加莱映射 $P:\Sigma\to\Sigma$ 的导数。

物理意义：自激振荡，与保守系统的连续周期轨道本质区别。

极限特例：范德坡振子
$$\ddot x-\mu(1-x^2)\dot x+x=0$$
存在稳定极限环。

典型应用：心脏起搏、振荡电路。

**知识点16.3　分岔（鞍结、跨临界、叉式、Hopf）**

定义：参数变化时相图拓扑改变。

推导鞍结分岔：一维 $\dot x=\mu-x^2$，平衡 $x=\pm\sqrt\mu$，$\mu<0$ 无、$\mu>0$ 两个，$\mu=0$ 合并。

推导 Hopf 分岔：当 $\tau=0,\Delta>0$ 且特征值穿过虚轴，产生小极限环。范式（极坐标）
$$\dot r=\mu r-ar^3,\qquad \dot\theta=\omega+br^2.$$
平衡径向方程 $\mu r-ar^3=0$ 得 $r=\sqrt{\mu/a}$（超临界）或 $r=\sqrt{-\mu/a}$（亚临界）。

物理意义：极限环幅值随参数连续增长。

极限特例：$\mu<0$ 稳定焦点。

典型应用：颤振、失速。

**知识点16.4　李雅普诺夫指数与混沌**

定义：相邻轨线分离 $\delta(t)\sim\delta_0 e^{\lambda t}$。

推导：
$$\lambda=\lim_{t\to\infty}\frac1t\ln\frac{\|\delta(t)\|}{\|\delta_0\|}.$$
$\lambda>0$ 表示敏感依赖初值（混沌）。对 Logistic 映射 $x_{n+1}=rx_n(1-x_n)$，Lyapunov 指数
$$\lambda=\lim_{N\to\infty}\frac1N\sum_{n=0}^{N-1}\ln|r-2rx_n|.$$
$r>r_\infty\approx3.5699$ 时有正指数与周期倍化级联，Feigenbaum 常数
$$\delta_F=\lim_{n\to\infty}\frac{r_n-r_{n-1}}{r_{n+1}-r_n}=4.6692\dots$$

物理意义：确定性系统可有不可预测随机性。

极限特例：$\lambda=0$ 对应临界，$\lambda<0$ 收敛。

典型应用：湍流、保密通信。

**知识点16.5　周期倍化分岔**

定义：周期轨道稳定性的翻转。

推导：映射 $x_{n+1}=f(x_n)$，不动点 $x^*=f(x^*)$，乘子
$$\Lambda=f'(x^*).$$
$|\Lambda|<1$ 稳定。当 $\Lambda=-1$ 时，二次迭代 $f(f(x))$ 出现新不动点，周期 $T$ 变 $2T$。

物理意义：Feigenbaum 级联是通向混沌的典型路径。

极限特例：$\Lambda=1$ 对应鞍结或跨临界。

典型应用：Logistic 映射、Duffing 振子。

**知识点16.6　庞加莱截面与分形维数**

定义：庞加莱截面把连续流降为离散映射。

推导：对 $n$ 维流，截面 $n-1$ 维，返回映射 $P$。混沌吸引子的关联维
$$D_2=\lim_{\epsilon\to0}\frac{\ln C(\epsilon)}{\ln\epsilon},\qquad C(\epsilon)=\frac{1}{N^2}\sum_{i,j}\Theta(\epsilon-\|\mathbf x_i-\mathbf x_j\|).$$
Kaplan-Yorke 猜想
$$D_{KY}=k+\frac{\sum_{i=1}^{k}\lambda_i}{|\lambda_{k+1}|},$$
其中 $\sum_{i=1}^k\lambda_i\ge0>\sum_{i=1}^{k+1}\lambda_i$。

物理意义：混沌吸引子具有分形结构，维数非整数。

极限特例：周期轨道关联维为整数。

典型应用：奇怪吸引子分析（Lorenz、Rössler）。

---
⬆ 上级领域：[[力学]]
