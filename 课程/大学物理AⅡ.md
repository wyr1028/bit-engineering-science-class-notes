---
课程名称: 大学物理AⅡ
学分: 4
总学时: 64
考核方式: 考试
课程属性: 必修
课程性质: 专业基础课程
开课单位: 物理学院
领域: 物质科学（物理·化学·生命）
tags:
  - 课程
  - 领域/物质科学（物理·化学·生命）
---

# 大学物理AⅡ


### 第1章 质点运动学

**1. 位矢、速度与加速度**

定义：质点在参考系中的位置用位置矢量 $\mathbf r(t)$ 描述；速度与加速度定义为
$$\mathbf v=\frac{\mathrm d\mathbf r}{\mathrm dt},\qquad \mathbf a=\frac{\mathrm d\mathbf v}{\mathrm dt}=\frac{\mathrm d^2\mathbf r}{\mathrm dt^2}.$$
假设：存在理想时钟与刚性量尺，可测 $\mathbf r(t)$ 且时间可微。
推导：由平均速度 $\bar{\mathbf v}=\Delta\mathbf r/\Delta t$ 取极限得瞬时速度；同理对 $\mathbf v$ 取极限得 $\mathbf a$。在直角坐标 $\mathbf r=x\mathbf i+y\mathbf j+z\mathbf k$ 中，基矢为常量，故
$$\mathbf v=\dot x\,\mathbf i+\dot y\,\mathbf j+\dot z\,\mathbf k,\qquad \mathbf a=\ddot x\,\mathbf i+\ddot y\,\mathbf j+\ddot z\,\mathbf k.$$
意义：$\mathbf v$ 沿轨迹切线，$\mathbf a$ 指向速度变化方向（不一定指向轨迹内侧）。
特例：匀速直线运动 $\mathbf a=0$；匀变速 $\mathbf v=\mathbf v_0+\mathbf a t$。
应用：卫星轨道参数由 $\mathbf r(t)$ 与 $\mathbf v(t)$ 递推确定。

**2. 自然坐标系中的切向与法向加速度**

定义：以轨迹切线 $\boldsymbol\tau$ 与主法线 $\mathbf n$ 为基，沿轨迹用弧长 $s$ 描述。
假设：轨迹光滑且曲率 $\kappa=1/\rho$ 有限。
公式：
$$\mathbf a=a_\tau\boldsymbol\tau+a_n\mathbf n,\qquad a_\tau=\frac{\mathrm dv}{\mathrm dt}=\frac{\mathrm d^2s}{\mathrm dt^2},\qquad a_n=\frac{v^2}{\rho}=v^2\kappa.$$
推导：因 $\mathrm d\boldsymbol\tau/\mathrm ds=\kappa\mathbf n$，有
$$\mathbf v=\frac{\mathrm ds}{\mathrm dt}\boldsymbol\tau,\quad \mathbf a=\frac{\mathrm d}{\mathrm dt}\!\left(\frac{\mathrm ds}{\mathrm dt}\boldsymbol\tau\right)=\frac{\mathrm d^2s}{\mathrm dt^2}\boldsymbol\tau+\frac{\mathrm ds}{\mathrm dt}\frac{\mathrm d\boldsymbol\tau}{\mathrm dt}.$$
而 $\dfrac{\mathrm d\boldsymbol\tau}{\mathrm dt}=\dfrac{\mathrm d\boldsymbol\tau}{\mathrm ds}\dfrac{\mathrm ds}{\mathrm dt}=\kappa v\,\mathbf n$，故 $a_n=\kappa v^2$。
意义：$a_\tau$ 改变速率大小，$a_n$ 改变速度方向。
特例：匀速圆周运动 $a_\tau=0,\ a_n=v^2/R$；直线运动 $\rho\to\infty,\ a_n=0$。
应用：公路弯道超高与限速设计。

**3. 圆周运动的角量与线量**

定义：角位置 $\theta$，角速度 $\omega=\dot\theta$，角加速度 $\alpha=\dot\omega$。
公式：
$$s=R\theta,\quad v=R\omega,\quad a_\tau=R\alpha,\quad a_n=R\omega^2.$$
推导：对 $s=R\theta$ 逐次求导即得。
意义：角量对所有半径点相同，线量随 $R$ 线性放大。
特例：匀角加速 $\omega=\omega_0+\alpha t,\ \theta=\theta_0+\omega_0t+\tfrac12\alpha t^2$。
应用：飞轮储能、旋转机械齿轮传动。

**4. 平面极坐标中的速度与加速度**

定义：$\mathbf r=r\mathbf e_r$，$\mathbf e_r,\mathbf e_\theta$ 随 $\theta$ 变化。
假设：平面运动，$\theta$ 连续可微。
推导：因 $\dot{\mathbf e}_r=\dot\theta\mathbf e_\theta,\ \dot{\mathbf e}_\theta=-\dot\theta\mathbf e_r$，
$$\mathbf v=\frac{\mathrm d(r\mathbf e_r)}{\mathrm dt}=\dot r\mathbf e_r+r\dot{\mathbf e}_r=\dot r\mathbf e_r+r\dot\theta\mathbf e_\theta.$$
再求导：
$$\mathbf a=\ddot r\mathbf e_r+\dot r\dot\theta\mathbf e_\theta+\dot r\dot\theta\mathbf e_\theta+r\ddot\theta\mathbf e_\theta+r\dot\theta(-\dot\theta\mathbf e_r),$$
$$\Rightarrow\ \mathbf a=(\ddot r-r\dot\theta^2)\mathbf e_r+(r\ddot\theta+2\dot r\dot\theta)\mathbf e_\theta.$$
意义：$-\!r\dot\theta^2$ 为向心项，$2\dot r\dot\theta$ 为科里奥利项。
特例：匀速圆周 $r$ 常数、$\dot\theta$ 常数时 $a_r=-r\omega^2$。
应用：行星运动、离心机构分析。

**5. 相对运动与伽利略变换**

定义：设 $S'$ 相对 $S$ 平动，位矢 $\mathbf R$。
公式：
$$\mathbf r=\mathbf r'+\mathbf R,\quad \mathbf v=\mathbf v'+\mathbf V,\quad \mathbf a=\mathbf a'+\mathbf A.$$
推导：对位矢关系逐次求时间导数，$\mathbf R$ 为两系原点相对位矢。
假设：$t'=t$（绝对时空）。
意义：加速度变换仅差牵连加速度，是惯性力分析的起点。
特例：$\mathbf A=0$ 时两系加速度相同，牛顿第二定律形式不变。
应用：风相对飞机、雨中行人持伞角度。

**6. 抛体运动与阻力模型**

定义：重力场中初速 $\mathbf v_0$ 的抛体。
假设：无阻力时 $\mathbf a=-g\mathbf j$。
推导：由 $\ddot x=0,\ \ddot y=-g$ 两次积分：
$$x=v_{0x}t,\qquad y=v_{0y}t-\tfrac12gt^2.$$
消去 $t$ 得轨道
$$y=x\tan\theta_0-\frac{gx^2}{2v_0^2\cos^2\theta_0}.$$
射程 $R=\dfrac{v_0^2\sin2\theta_0}{g}$。
意义：抛物线轨迹，$45^\circ$ 射程最大。
特例：考虑线性阻力 $\mathbf f=-b\mathbf v$ 时速度呈指数趋近，轨迹偏离抛物线且末端趋近竖直。
应用：弹道计算、喷雾液滴沉积。

**7. 单摆的非线性运动方程**

定义：长为 $l$ 的刚性摆，摆角 $\theta$。
推导：切向受力 $-mg\sin\theta$，弧长 $s=l\theta$ 故 $a_\tau=l\ddot\theta$：
$$ml\ddot\theta=-mg\sin\theta\ \Rightarrow\ \ddot\theta+\frac{g}{l}\sin\theta=0.$$
小角近似 $\sin\theta\approx\theta$：$\ddot\theta+\omega_0^2\theta=0,\ \omega_0=\sqrt{g/l}$。
意义：非线性恢复力使周期随振幅增大。
特例：大角周期 $T=4\sqrt{l/g}\,K(\sin(\theta_0/2))$，$K$ 为第一类完全椭圆积分。
应用：时钟、地震仪、能量采集。

### 第2章 牛顿运动定律与惯性系

**1. 牛顿三定律与惯性系**

定义：惯性系中不受力物体保持静止或匀速直线运动；$\mathbf F=\mathrm d\mathbf p/\mathrm dt$；作用与反作用等大反向。
假设：存在质量、力与惯性参考系。
推导：对常质量 $m$，$\mathbf F=\dfrac{\mathrm d(m\mathbf v)}{\mathrm dt}=m\mathbf a$。
意义：力是动量变化率，是自然界最基本的动力学公设。
特例：变质量时须用 $\mathbf F=\dot m\mathbf v+m\dot{\mathbf v}$。
应用：火箭、雨滴增长。

**2. 动量定理与冲量**

定义：冲量 $\mathbf J=\int_{t_1}^{t_2}\mathbf F\,\mathrm dt$。
推导：$\mathbf F=\dfrac{\mathrm d\mathbf p}{\mathrm dt}\Rightarrow \dfrac{\mathrm d\mathbf p}{\mathrm dt}\,\mathrm dt=\mathrm d\mathbf p$，积分得
$$\mathbf J=\int_{t_1}^{t_2}\mathbf F\,\mathrm dt=\mathbf p_2-\mathbf p_1.$$
意义：力的时间累积改变动量。
特例：碰撞中冲击力大而作用时间短，用平均力 $\bar F=J/\Delta t$。
应用：安全气囊、缓冲器设计。

**3. 变质量系统与火箭方程**

定义：系统质量随时间变化，排出物相对系统速度为 $\mathbf u$。
假设：无外力、排气流相对速度恒定。
推导：设时刻 $t$ 系统质量 $m$、速度 $v$；$t+\mathrm dt$ 排出质量 $-\mathrm dm>0$、速度 $v-u$。动量守恒：
$$mv=m(v+\mathrm dv)+(-\mathrm dm)(v-u),$$
展开：$mv=mv+m\,\mathrm dv-u\,\mathrm dm+\mathrm dm\,\mathrm dv$，略去二阶量：
$$m\,\mathrm dv=u\,\mathrm dm\ \Rightarrow\ \mathrm dv=-u\frac{\mathrm dm}{m}.$$
积分得齐奥尔科夫斯基公式
$$v-v_0=u\ln\frac{m_0}{m}.$$
意义：速度增量由喷速与质量比对数决定。
特例：有外力时 $m\dot{\mathbf v}=\mathbf F_{\rm ext}+\mathbf u\,\dot m$（注意 $\dot m<0$）。
应用：多级火箭、变质量雨滴下落。

**4. 非惯性系与惯性力**

定义：相对惯性系加速的参考系称为非惯性系。
推导：伽利略变换 $\mathbf a=\mathbf a'+\mathbf A$，代入 $\mathbf F=m\mathbf a$：
$$m\mathbf a'=\mathbf F-m\mathbf A.$$
故在非惯性系中须引入惯性力 $\mathbf F_{\rm in}=-m\mathbf A$。
意义：惯性力无施力物体，只在非惯性系中作为形式力出现。
特例：自由落体电梯中 $\mathbf F_{\rm in}$ 抵消重力，失重。
应用：加速车辆中乘客倾覆、飞机机动载荷。

**5. 转动参考系中的离心力与科里奥利力**

定义：以角速度 $\boldsymbol\omega$ 转动的参考系。
公式：
$$\mathbf F_{\rm cen}=m\omega^2\mathbf r_\perp,\qquad \mathbf F_{\rm Cor}=-2m\,\boldsymbol\omega\times\mathbf v'.$$
推导：转动系中速度变换 $\mathbf v=\mathbf v'+\boldsymbol\omega\times\mathbf r'$，再对时间求导并注意基矢旋转 $\dot{\mathbf e}=\boldsymbol\omega\times\mathbf e$，整理得
$$\mathbf a=\mathbf a'+2\boldsymbol\omega\times\mathbf v'+\boldsymbol\omega\times(\boldsymbol\omega\times\mathbf r')+\dot{\boldsymbol\omega}\times\mathbf r'.$$
移项即得惯性力。
意义：离心力使旋转系物体外推，科里奥利力使运动轨迹偏转。
特例：$\boldsymbol\omega=0$ 退化为平动非惯性系。
应用：傅科摆、信风偏转、离心分离机、涡轮叶轮载荷。

**6. 万有引力与重力随纬度变化**

定义：$\mathbf F=-G\dfrac{Mm}{r^2}\hat{\mathbf r}$。
推导：地表随地球自转的参考系中含离心加速度 $r\omega^2\cos\lambda$，有效重力
$$g(\lambda)=g_0-\omega^2R\cos^2\lambda\ (\text{近似}),$$
更精确
$$g(\lambda)=\frac{GM}{R^2}-\omega^2R\cos^2\lambda.$$
意义：赤道重力最小、两极最大。
特例：忽略自转 $\omega=0$ 时 $g=GM/R^2$。
应用：重力测量、卫星轨道摄动。

### 第3章 功、能与守恒定律

**1. 功与动能定理**

定义：力沿路径做功 $W=\int_a^b\mathbf F\cdot\mathrm d\mathbf r$。
推导：对质点用 $\mathbf F=m\dfrac{\mathrm d\mathbf v}{\mathrm dt}$，
$$\mathbf F\cdot\mathrm d\mathbf r=m\frac{\mathrm d\mathbf v}{\mathrm dt}\cdot\mathbf v\,\mathrm dt=m\mathbf v\cdot\mathrm d\mathbf v=\mathrm d\!\left(\tfrac12mv^2\right).$$
积分得
$$W=\int_a^b\mathbf F\cdot\mathrm d\mathbf r=\tfrac12mv_b^2-\tfrac12mv_a^2=\Delta K.$$
意义：合外力功等于动能增量。
特例：$\mathbf F\perp\mathbf v$ 时 $W=0$（洛伦兹力、向心力）。
应用：能量法求解变力位移。

**2. 保守力与势能**

定义：做功只依赖端点的力为保守力，$\oint\mathbf F\cdot\mathrm d\mathbf r=0$。
推导：路径无关 $\Rightarrow$ 存在标量 $U$ 使 $\mathbf F\cdot\mathrm d\mathbf r=-\mathrm dU$，于是 $\mathbf F=-\nabla U$。反之若 $\mathbf F=-\nabla U$，则 $\oint\mathbf F\cdot\mathrm d\mathbf r=-\oint\mathrm dU=0$。
意义：势能是位置的函数，零点可任选。
特例：重力 $U=mgh$；弹力 $U=\tfrac12kx^2$；万有引力 $U=-GMm/r$。
应用：保守系统用势能曲线判断运动范围。

**3. 机械能守恒**

定义：$E=K+U$。
推导：由动能定理 $W_{\rm all}=\Delta K$。将力分为保守与非保守：$W_c+W_{nc}=\Delta K$，而 $W_c=-\Delta U$，故
$$W_{nc}=\Delta(K+U)=\Delta E.$$
若 $W_{nc}=0$ 则 $E$ 恒为常数。
意义：耗散为零时总机械能守恒。
特例：有摩擦时 $\Delta E=-W_f<0$。
应用：摆、弹簧振子、滑轨滑块。

**4. 势能曲线与平衡稳定性**

定义：一维势能 $U(x)$，力 $F=-U'(x)$。
推导：平衡点 $F=0\Rightarrow U'(x_0)=0$。在 $x_0$ 附近泰勒展开
$$U(x)\approx U(x_0)+\tfrac12U''(x_0)(x-x_0)^2.$$
$U''>0$ 时能量极小、稳定；$U''<0$ 时极大、不稳定。
意义：曲率决定微小振动频率 $\omega=\sqrt{U''/m}$。
特例：$\tfrac12kx^2$ 势给出简谐振动。
应用：分子键、双稳电路。

**5. 角动量与力矩**

定义：$\mathbf L=\mathbf r\times\mathbf p$，$\boldsymbol\tau=\mathbf r\times\mathbf F$。
推导：
$$\frac{\mathrm d\mathbf L}{\mathrm dt}=\frac{\mathrm d\mathbf r}{\mathrm dt}\times\mathbf p+\mathbf r\times\frac{\mathrm d\mathbf p}{\mathrm dt}=\mathbf v\times(m\mathbf v)+\mathbf r\times\mathbf F=\boldsymbol\tau,$$
因 $\mathbf v\times\mathbf v=0$。
意义：$\boldsymbol\tau=0$ 时角动量守恒。
特例：中心力场 $\mathbf r\parallel\mathbf F$ 故 $\boldsymbol\tau=0$，$\mathbf L$ 守恒。
应用：开普勒第二定律、滑冰收臂加速。

**6. 中心力场与面积速度**

定义：中心力 $\mathbf F=f(r)\hat{\mathbf r}$。
推导：$\mathbf L=m\mathbf r\times\mathbf v$ 守恒，方向恒定说明轨道在过力心的平面内；大小
$$L=mr^2\dot\theta.$$
扇形面积 $\mathrm dA=\tfrac12r^2\mathrm d\theta$，故
$$\frac{\mathrm dA}{\mathrm dt}=\tfrac12r^2\dot\theta=\frac{L}{2m}=\text{const}.$$
意义：面积速度恒为常数（开普勒第二定律）。
特例：圆轨道 $r$ 常数时 $\dot\theta$ 常数。
应用：行星与卫星运动、有心力散射。

**7. 两体约化与有效势**

定义：两体质量 $m_1,m_2$，相对位矢 $\mathbf r=\mathbf r_1-\mathbf r_2$。
推导：质心 $\mathbf R$ 匀速运动；相对运动
$$\mu\ddot{\mathbf r}=f(r)\hat{\mathbf r},\qquad \mu=\frac{m_1m_2}{m_1+m_2}.$$
径向能量
$$E=\tfrac12\mu\dot r^2+\frac{L^2}{2\mu r^2}+U(r)=K_r+U_{\rm eff}(r),$$
$$U_{\rm eff}(r)=U(r)+\frac{L^2}{2\mu r^2}.$$
意义：离心势垒阻止落入力心。
特例：$L=0$ 时纯径向运动。
应用：双星、氢原子、行星轨道稳定性。

### 第4章 质点系动力学与刚体

**1. 质心运动定理**

定义：$\mathbf R_c=\dfrac{1}{M}\sum_im_i\mathbf r_i,\ M=\sum m_i$。
推导：内力成对抵消：
$$M\ddot{\mathbf R}_c=\sum_i\mathbf F_i^{(\rm ext)}.$$
意义：质心像受全部外力的单质点一样运动。
特例：合外力为零时质心匀速直线运动。
应用：爆炸碎片、双星系统。

**2. 变质量与碰撞守恒**

定义：完全弹性、完全非弹性、部分弹性碰撞。
推导：动量守恒 $m_1\mathbf v_1+m_2\mathbf v_2=m_1\mathbf v_1'+m_2\mathbf v_2'$。弹性碰撞联立动能守恒
$$\tfrac12m_1v_1^2+\tfrac12m_2v_2^2=\tfrac12m_1v_1'^2+\tfrac12m_2v_2'^2,$$
解得
$$v_1'=\frac{m_1-m_2}{m_1+m_2}v_1+\frac{2m_2}{m_1+m_2}v_2.$$
意义：恢复系数 $e=(v_2'-v_1')/(v_1-v_2)$，$e=1$ 弹性、$e=0$ 完全非弹性。
特例：等质量弹性碰撞交换速度。
应用：台球、粒子对撞。

**3. 转动惯量**

定义：$I=\sum_im_ir_i^2$，连续体 $I=\int r^2\,\mathrm dm$。
推导：动能
$$K=\sum_i\tfrac12m_iv_i^2=\sum_i\tfrac12m_i(\omega r_i)^2=\tfrac12\!\left(\sum_im_ir_i^2\right)\omega^2=\tfrac12I\omega^2.$$
意义：转动惯量是转动惯性的度量，依赖质量分布与转轴。
特例：细杆绕中心 $I=\tfrac1{12}ml^2$；圆盘 $I=\tfrac12mR^2$；球 $I=\tfrac25mR^2$。
应用：飞轮、机械臂、陀螺。

**4. 平行轴定理**

定义：绕通过质心的轴转动惯量 $I_c$，平行移轴距离 $d$。
推导：
$$I=\sum m_i|\mathbf r_i'+\mathbf d|^2=\sum m_ir_i'^2+2\mathbf d\cdot\sum m_i\mathbf r_i'+Md^2.$$
质心系 $\sum m_i\mathbf r_i'=0$，故
$$I=I_c+Md^2.$$
意义：绕质心轴转动惯量最小。
特例：细杆绕端点 $I=\tfrac1{12}ml^2+m(l/2)^2=\tfrac13ml^2$。
应用：结构件惯性设计。

**5. 刚体定轴转动方程**

定义：力矩 $M_z=\sum r_iF_{i\perp}$。
推导：$L_z=I\omega$，故
$$M_z=\frac{\mathrm dL_z}{\mathrm dt}=\frac{\mathrm d(I\omega)}{\mathrm dt}=I\alpha\ (\text{刚体}).$$
意义：转动与平动 $\mathbf F=m\mathbf a$ 对应。
特例：$M_z=0$ 时角动量守恒。
应用：滑轮系统、陀螺稳定。

**6. 刚体平面运动与滚动**

定义：平面运动 = 质心平动 + 绕质心转动。
推导：总动能
$$K=\tfrac12Mv_c^2+\tfrac12I_c\omega^2.$$
纯滚动条件 $v_c=R\omega$。
意义：滚动分解为平动与转动。
特例：纯滚动无摩擦损耗，静摩擦不做功。
应用：车轮、圆柱沿斜面滚下，$a=\dfrac{g\sin\theta}{1+I_c/(MR^2)}$。

**7. 欧拉动力学方程**

定义：在主轴坐标系中角动量 $\mathbf L=I_1\omega_1\mathbf e_1+I_2\omega_2\mathbf e_2+I_3\omega_3\mathbf e_3$。
推导：在惯性系 $\dfrac{\mathrm d\mathbf L}{\mathrm dt}=\mathbf N$；转动系中
$$\frac{\mathrm d\mathbf L}{\mathrm dt}=\left(\frac{\mathrm d\mathbf L}{\mathrm dt}\right)_{\rm body}+\boldsymbol\omega\times\mathbf L.$$
分量展开：
$$I_1\dot\omega_1-(I_2-I_3)\omega_2\omega_3=N_1,$$
$$I_2\dot\omega_2-(I_3-I_1)\omega_3\omega_1=N_2,$$
$$I_3\dot\omega_3-(I_1-I_2)\omega_1\omega_2=N_3.$$
意义：描述刚体定点转动。
特例：自由陀螺 $N=0$，欧拉情形。网球拍定理（中间轴不稳定）由此可证。
应用：航天器姿态动力学、陀螺仪。

**8. 拉格朗日陀螺与进动**

定义：轴对称自由陀螺，自转角速度 $\dot\psi$、进动 $\dot\phi$、章动 $\dot\theta$。
推导：拉格朗日函数
$$L=\tfrac12I_1(\dot\theta^2+\dot\phi^2\sin^2\theta)+\tfrac12I_3(\dot\psi+\dot\phi\cos\theta)^2-Mgl\cos\theta.$$
$p_\psi,p_\phi$ 为循环动量，得
$$I_3(\dot\psi+\dot\phi\cos\theta)=\text{const},\quad I_1\dot\phi\sin^2\theta+I_3(\dot\psi+\dot\phi\cos\theta)\cos\theta=\text{const}.$$
意义：重力矩使陀螺进动而不倒。
特例：快转近似进动角速度 $\Omega\approx Mgl/L$。
应用：导弹/卫星姿态稳定、罗盘。

### 第5章 拉格朗日力学

**1. 广义坐标与约束**

定义：用 $q_1,\dots,q_n$ 描述系统位形，$n=3N-k$（$k$ 完整约束数）。
假设：约束为完整、理想（约束力不做虚功）。
推导：由虚功原理 $\sum_a(\mathbf F_a-m_a\ddot{\mathbf r}_a)\cdot\delta\mathbf r_a=0$，约束力项消失。
意义：自动消去未知约束力，减少方程数。
特例：单摆用 $\theta$ 一个广义坐标。
应用：机构运动学、机器人。

**2. 达朗贝尔原理到拉格朗日方程**

定义：$\displaystyle\sum_a(\mathbf F_a-m_a\ddot{\mathbf r}_a)\cdot\delta\mathbf r_a=0$。
推导：由 $\mathbf r_a=\mathbf r_a(q,t)$，
$$\delta\mathbf r_a=\sum_i\frac{\partial\mathbf r_a}{\partial q_i}\delta q_i,\qquad
\mathbf v_a=\sum_i\frac{\partial\mathbf r_a}{\partial q_i}\dot q_i+\frac{\partial\mathbf r_a}{\partial t}.$$
关键是两个恒等式
$$\frac{\partial\mathbf v_a}{\partial\dot q_i}=\frac{\partial\mathbf r_a}{\partial q_i},\qquad
\frac{\mathrm d}{\mathrm dt}\frac{\partial\mathbf r_a}{\partial q_i}=\frac{\partial\mathbf v_a}{\partial q_i}.$$
于是虚功中惯性项
$$\sum_am_a\ddot{\mathbf r}_a\cdot\frac{\partial\mathbf r_a}{\partial q_i}
=\frac{\mathrm d}{\mathrm dt}\frac{\partial T}{\partial\dot q_i}-\frac{\partial T}{\partial q_i},$$
其中 $T=\tfrac12\sum m_a\dot{\mathbf r}_a^2$。定义广义力 $Q_i=\sum_a\mathbf F_a\cdot\dfrac{\partial\mathbf r_a}{\partial q_i}$，得
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial T}{\partial\dot q_i}-\frac{\partial T}{\partial q_i}=Q_i.$$
若力来自势能 $Q_i=-\partial V/\partial q_i$，令 $L=T-V$，则
$$\boxed{\ \frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_i}-\frac{\partial L}{\partial q_i}=0\ }$$
意义：坐标无关、形式统一，便于处理复杂约束。
特例：直角坐标下退化为牛顿方程。
应用：多体系统、机械臂动力学。

**3. 循环坐标与广义动量守恒**

定义：不显含于 $L$ 的坐标 $q_j$ 称循环坐标。
推导：拉格朗日方程中 $\partial L/\partial q_j=0\Rightarrow\dfrac{\mathrm d}{\mathrm dt}\dfrac{\partial L}{\partial\dot q_j}=0$，故
$$p_j=\frac{\partial L}{\partial\dot q_j}=\text{const}.$$
意义：对称性对应守恒律（诺特）。
特例：中心力场中 $\theta$ 循环，$L_z$ 守恒。
应用：角动量、平移动量守恒。

**4. 单摆与双摆的拉格朗日建立**

定义：单摆 $L=\tfrac12ml^2\dot\theta^2+mgl\cos\theta$。
推导：$\dfrac{\partial L}{\partial\dot\theta}=ml^2\dot\theta,\ \dfrac{\partial L}{\partial\theta}=-mgl\sin\theta$，得 $ml^2\ddot\theta+mgl\sin\theta=0$。
双摆：设角 $\theta_1,\theta_2$，动能
$$T=\tfrac12m_1l_1^2\dot\theta_1^2+\tfrac12m_2\!\left(l_1^2\dot\theta_1^2+l_2^2\dot\theta_2^2+2l_1l_2\dot\theta_1\dot\theta_2\cos(\theta_1-\theta_2)\right),$$
$$V=-(m_1+m_2)gl_1\cos\theta_1-m_2gl_2\cos\theta_2.$$
意义：耦合非线性方程，展示混沌。
特例：小振动线性化得简正模。
应用：桁架摆、机械臂。

**5. 约束力与拉格朗日乘子**

定义：非完整或需显式约束力时引入乘子 $\lambda_k$。
推导：约束 $f_k(\mathbf r,t)=0$，虚功中约束力项 $\sum_k\lambda_k\nabla f_k\cdot\delta\mathbf r$，得
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_i}-\frac{\partial L}{\partial q_i}=\sum_k\lambda_k\frac{\partial f_k}{\partial q_i}.$$
意义：乘子正比于约束力。
特例：斜面滑块可求法向力。
应用：接触力、关节反力计算。

**6. 广义势与电磁场中的拉格朗日函数**

定义：速度相关势 $U(q,\dot q,t)$，力 $Q_i=-\dfrac{\partial U}{\partial q_i}+\dfrac{\mathrm d}{\mathrm dt}\dfrac{\partial U}{\partial\dot q_i}$。
推导：带电粒子洛伦兹力 $\mathbf F=q(\mathbf E+\mathbf v\times\mathbf B)$ 可由
$$U=q\varphi-q\mathbf A\cdot\mathbf v$$
导出，故
$$L=\tfrac12m\dot{\mathbf r}^2-q\varphi+q\mathbf A\cdot\dot{\mathbf r}.$$
意义：洛伦兹力自然嵌入拉格朗日框架。
特例：无磁场时退化为 $L=T-q\varphi$。
应用：加速器、等离子体。

### 第6章 哈密顿力学

**1. 勒让德变换与哈密顿函数**

定义：$H(q,p,t)=\sum_i p_i\dot q_i-L$，其中 $p_i=\partial L/\partial\dot q_i$。
推导：对 $L$ 微分
$$\mathrm dL=\sum_i\frac{\partial L}{\partial q_i}\mathrm dq_i+\sum_i\frac{\partial L}{\partial\dot q_i}\mathrm d\dot q_i+\frac{\partial L}{\partial t}\mathrm dt
=\sum_i\dot p_i\,\mathrm dq_i+\sum_ip_i\,\mathrm d\dot q_i+\frac{\partial L}{\partial t}\mathrm dt.$$
而 $\mathrm dH=\sum_i\dot q_i\,\mathrm dp_i+\sum_ip_i\,\mathrm d\dot q_i-\mathrm dL$
$=\sum_i\dot q_i\,\mathrm dp_i-\sum_i\dot p_i\,\mathrm dq_i-\dfrac{\partial L}{\partial t}\mathrm dt$。
比较系数得正则方程。
意义：$H$ 通常为总能量。
特例：稳定保守系统 $H=T+V=E$。
应用：相空间分析、统计力学。

**2. 哈密顿正则方程**

定义：$\dot q_i=\dfrac{\partial H}{\partial p_i},\ \dot p_i=-\dfrac{\partial H}{\partial q_i}$。
推导：由勒让德变换比较系数直接得到；或
$$\dot q_i=\frac{\partial H}{\partial p_i},\qquad -\frac{\partial H}{\partial q_i}=\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_i}=\dot p_i.$$
意义：一阶对称方程，便于正则变换与量子化。
特例：谐振子 $H=\dfrac{p^2}{2m}+\dfrac12m\omega^2q^2$。
应用：辛积分、分子动力学。

**3. 泊松括号与守恒量**

定义：$\{A,B\}=\sum_i\!\left(\dfrac{\partial A}{\partial q_i}\dfrac{\partial B}{\partial p_i}-\dfrac{\partial A}{\partial p_i}\dfrac{\partial B}{\partial q_i}\right)$。
推导：对 $A(q,p,t)$ 求全导
$$\frac{\mathrm dA}{\mathrm dt}=\frac{\partial A}{\partial t}+\sum_i\!\left(\frac{\partial A}{\partial q_i}\dot q_i+\frac{\partial A}{\partial p_i}\dot p_i\right)
=\frac{\partial A}{\partial t}+\{A,H\}.$$
故 $A$ 守恒 $\Leftrightarrow\partial A/\partial t+\{A,H\}=0$。
意义：守恒量代数结构对应李代数。
特例：$\{q_i,p_j\}=\delta_{ij}$，$\{L_x,L_y\}=L_z$。
应用：可积系统、量子对应 $[\hat A,\hat B]=i\hbar\{A,B\}$。

**4. 正则变换与生成函数**

定义：使正则方程保持形式的变换 $(q,p)\to(Q,P)$。
推导：要求 $\sum p\,\mathrm dq-H\,\mathrm dt=\sum P\,\mathrm dQ-K\,\mathrm dt+\mathrm dF$，其中 $F$ 为生成函数。取 $F_1(q,Q,t)$ 则
$$p=\frac{\partial F_1}{\partial q},\quad P=-\frac{\partial F_1}{\partial Q},\quad K=H+\frac{\partial F_1}{\partial t}.$$
意义：通过变换简化哈密顿量。
特例：恒等变换 $F_1=qQ$。
应用：简谐振子化、可积系统。

**5. 哈密顿-雅可比方程**

定义：寻找使 $K=0$ 的生成函数 $S(q,P,t)$（哈密顿主函数）。
推导：由 $F_2=S(q,P,t)$ 有 $p=\partial S/\partial q$，$Q=\partial S/\partial P$。要求 $K=H+\partial S/\partial t=0$：
$$H\!\left(q,\frac{\partial S}{\partial q},t\right)+\frac{\partial S}{\partial t}=0,$$
$$S=\int L\,\mathrm dt+\text{const}.$$
意义：将运动化为求偏微分方程的完全积分。
特例：稳定 $H$ 时 $S=W(q,P)-Et$，得 $H(q,\partial W/\partial q)=E$（哈密顿特征函数）。
应用：可积系统、量子力学 WKB。

**6. 刘维尔定理**

定义：相空间体积元 $\mathrm d\Gamma=\prod_i\mathrm dq_i\mathrm dp_i$。
推导：相流速度 $\dot{\mathbf z}=(\dot q,\dot p)$，散度
$$\nabla\!\cdot\!\dot{\mathbf z}=\sum_i\!\left(\frac{\partial\dot q_i}{\partial q_i}+\frac{\partial\dot p_i}{\partial p_i}\right)
=\sum_i\!\left(\frac{\partial^2H}{\partial q_i\partial p_i}-\frac{\partial^2H}{\partial p_i\partial q_i}\right)=0.$$
故相体积随流不变：$\dfrac{\mathrm d\Gamma}{\mathrm dt}=0$。
意义：哈密顿流不可压缩，是统计力学基础。
特例：耗散系统相体积收缩。
应用：系综理论、辛几何。

### 第7章 变分原理与对称性

**1. 泛函与变分**

定义：泛函 $J[y]=\int_{x_1}^{x_2}f(x,y,y')\mathrm dx$。
推导：令 $y\to y+\epsilon\eta$，$\eta(x_1)=\eta(x_2)=0$，
$$\delta J=\int\left(\frac{\partial f}{\partial y}\eta+\frac{\partial f}{\partial y'}\eta'\right)\mathrm dx
=\int\!\left(\frac{\partial f}{\partial y}-\frac{\mathrm d}{\mathrm dx}\frac{\partial f}{\partial y'}\right)\eta\,\mathrm dx,$$
末步用分部积分 $[\ \partial f/\partial y'\,\eta\ ]_{x_1}^{x_2}=0$。
意义：$\delta J=0$ 对应泛函极值。
特例：$f=\sqrt{1+y'^2}$ 给出直线。
应用：最速降线、测地线。

**2. 欧拉-拉格朗日方程**

定义：令 $\delta J=0$ 得
$$\frac{\partial f}{\partial y}-\frac{\mathrm d}{\mathrm dx}\frac{\partial f}{\partial y'}=0.$$
推导：由基本引理，若 $\int g(x)\eta(x)\mathrm dx=0$ 对任意 $\eta$ 成立，则 $g(x)=0$。
意义：变分问题化为微分方程。
特例：多变量与高阶导推广。
应用：场论、最优控制。

**3. 哈密顿原理**

定义：真实路径使作用量 $S=\int_{t_1}^{t_2}L(q,\dot q,t)\mathrm dt$ 取驻值。
推导：$\delta S=0$ 逐分量给出
$$\delta S=\int\sum_i\!\left(\frac{\partial L}{\partial q_i}\delta q_i+\frac{\partial L}{\partial\dot q_i}\delta\dot q_i\right)\mathrm dt
=\int\sum_i\!\left(\frac{\partial L}{\partial q_i}-\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_i}\right)\delta q_i\,\mathrm dt=0,$$
即拉格朗日方程。
意义：整个力学统一为单一变分原理。
特例：牛顿、光学费马原理皆为其特例。
应用：场论、路径积分。

**4. 诺特定理**

定义：作用量在连续变换下的不变性对应守恒流。
推导：设变换 $q_i\to q_i+\epsilon K_i$ 使 $\delta L=\epsilon\dfrac{\mathrm dF}{\mathrm dt}$。由 $\delta L=\sum\dfrac{\partial L}{\partial q_i}\delta q_i+\sum p_i\delta\dot q_i$ 及欧拉-拉格朗日，得
$$\frac{\mathrm d}{\mathrm dt}\!\left(\sum_ip_iK_i-F\right)=0.$$
守恒量
$$Q=\sum_ip_iK_i-F.$$
意义：时间平移→能量，空间平移→动量，转动→角动量。
特例：中心力场转动对称给出 $L_z$。
应用：粒子物理守恒律、场论。

**5. 测地线与最速降线**

定义：弯曲空间最短路径泛函 $\displaystyle J=\int\sqrt{g_{ij}\dot x^i\dot x^j}\,\mathrm ds$。
推导：取 $f=\sqrt{g_{ij}\dot x^i\dot x^j}$，欧拉-拉格朗日给
$$\ddot x^k+\Gamma^k_{ij}\dot x^i\dot x^j=0,$$
$$\Gamma^k_{ij}=\tfrac12g^{kl}(\partial_ig_{lj}+\partial_jg_{il}-\partial_lg_{ij}).$$
意义：测地线方程是广义相对论中自由粒子运动方程。
特例：平直空间退化为直线。
应用：广义相对论、导航。

**6. 守恒量与循环坐标的变分解释**

定义：循环坐标对应平移对称。
推导：$L$ 不含 $q_i$ 时 $p_i=\partial L/\partial\dot q_i$ 守恒；这正是诺特定理中平移 $K_i=1$ 的情形，$Q=p_i$。
意义：对称性与守恒律一一对应。
特例：时间平移不变 $L$ 不显含 $t$ 时 $H$ 守恒。
应用：解析力学、场论诺特流。

### 第8章 中心力场与天体力学

**1. 两体问题约化**

定义：约化质量 $\mu=\dfrac{m_1m_2}{m_1+m_2}$，相对坐标 $\mathbf r$。
推导：总动能分解
$$T=\tfrac12M\dot{\mathbf R}^2+\tfrac12\mu\dot{\mathbf r}^2,$$
势仅依赖 $|\mathbf r|$，故质心运动与相对运动解耦。
意义：两体问题化为单体问题。
特例：$m_2\to\infty$ 时 $\mu\to m_1$，固定力心。
应用：双星、地月系统。

**2. 有效势与运动分类**

定义：$U_{\rm eff}=U(r)+\dfrac{L^2}{2\mu r^2}$。
推导：径向能量 $E=\tfrac12\mu\dot r^2+U_{\rm eff}(r)$，由 $\dot r^2\ge0$ 得 $E\ge U_{\rm eff}$。
意义：$U_{\rm eff}$ 的极小与零点决定轨道类型。
特例：引力场 $U=-k/r$：$E<0$ 椭圆，$E=0$ 抛物，$E>0$ 双曲。
应用：轨道能量判据。

**3. 比内公式与轨道方程**

定义：$u=1/r$，$u(\theta)$。
推导：由 $L=\mu r^2\dot\theta$，
$$\dot r=-\frac{L}{\mu}\frac{\mathrm du}{\mathrm d\theta},\qquad
\ddot r=-\frac{L^2u^2}{\mu^2}\frac{\mathrm d^2u}{\mathrm d\theta^2}.$$
径向方程 $\mu(\ddot r-r\dot\theta^2)=f(1/u)$ 化为比内公式
$$f\!\left(\frac1u\right)=-\frac{L^2u^2}{\mu}\!\left(\frac{\mathrm d^2u}{\mathrm d\theta^2}+u\right).$$
意义：由力律求轨道或反之。
特例：$f=-k/r^2$ 则 $u''+u=\mu k/L^2$。
应用：轨道反演问题。

**4. 平方反比引力与圆锥曲线轨道**

定义：$f(r)=-k/r^2$，$k=GM\mu$ 或 $GMm$。
推导：解 $u''+u=\mu k/L^2$：
$$u=\frac{\mu k}{L^2}(1+e\cos(\theta-\theta_0)).$$
故
$$r=\frac{L^2/(\mu k)}{1+e\cos(\theta-\theta_0)}
=\frac{p}{1+e\cos(\theta-\theta_0)},\qquad p=\frac{L^2}{\mu k}.$$
意义：轨道为圆锥曲线，$e$ 为离心率。
特例：$e=0$ 圆、$0<e<1$ 椭圆、$e=1$ 抛物、$e>1$ 双曲。
应用：开普勒轨道、卫星设计。

**5. 开普勒三定律**

定义：行星轨道、面积速度、周期。
推导：第一定律即上一式椭圆。第二定律由 $L=\mu r^2\dot\theta$ 常数。第三定律：椭圆面积 $A=\pi ab$，
$$T=\frac{A}{\mathrm dA/\mathrm dt}=\frac{2\pi ab\mu}{L}.$$
由 $b=a\sqrt{1-e^2}$、$p=a(1-e^2)$、$L^2=\mu kp=\mu k\,a(1-e^2)$ 得
$$T^2=\frac{4\pi^2\mu^2a^2b^2}{L^2}=\frac{4\pi^2a^3}{G(m_1+m_2)}.$$
意义：$T^2\propto a^3$。
特例：圆轨道 $a=R$，$T=2\pi\sqrt{R^3/GM}$。
应用：测定天体质量、GPS 轨道。

**6. 引力势能与逃逸速度**

定义：$U=-GMm/r$。
推导：由 $\displaystyle U(r)-U(\infty)=-\int_\infty^r\!\left(-\frac{GMm}{r'^2}\right)\mathrm dr'=-\frac{GMm}{r}$。
逃逸条件 $\tfrac12mv_e^2+U(R)=0$：
$$v_e=\sqrt{\frac{2GM}{R}}=\sqrt{2gR}.$$
意义：总能量非负即可逃逸。
特例：地球 $v_e\approx11.2\,\mathrm{km/s}$。
应用：火箭发射、黑洞视界。

**7. 卢瑟福散射**

定义：带电粒子在库仑场散射，$k=Z_1Z_2e^2/4\pi\epsilon_0$。
推导：轨道 $r=p/(1+e\cos\theta)$，渐近角由 $r\to\infty$ 给 $1+e\cos(\theta_\infty-\theta_0)=0$。散射角
$$\cot\!\left(\frac{\vartheta}{2}\right)=\frac{\mu v_0^2b}{k},$$
其中 $b$ 为瞄准距。微分截面
$$\frac{\mathrm d\sigma}{\mathrm d\Omega}=\left(\frac{k}{4E\sin^2(\vartheta/2)}\right)^2.$$
意义：验证原子核式结构。
特例：同种电荷排斥 $e>1$ 双曲。
应用：卢瑟福实验、离子束分析。

### 第9章 小振动与简正模

**1. 平衡与线性化**

定义：平衡 $\partial V/\partial q_i=0$。
推导：在平衡点展开 $V(q)=\tfrac12\sum K_{ij}\xi_i\xi_j$，$K_{ij}=\partial^2V/\partial q_i\partial q_j$；动能 $T=\tfrac12\sum M_{ij}\dot\xi_i\dot\xi_j$。
意义：稳定平衡 $K$ 正定。
特例：单自由度 $K=k$，$\omega=\sqrt{k/M}$。
应用：分子振动、结构模态。

**2. 拉格朗日方程与特征值问题**

定义：$L=\tfrac12\dot{\boldsymbol\xi}^{\rm T}M\dot{\boldsymbol\xi}-\tfrac12\boldsymbol\xi^{\rm T}K\boldsymbol\xi$。
推导：拉格朗日方程
$$M\ddot{\boldsymbol\xi}+K\boldsymbol\xi=0.$$
设 $\boldsymbol\xi=\mathbf a e^{i\omega t}$：
$$(K-\omega^2M)\mathbf a=0\ \Rightarrow\ \det(K-\omega^2M)=0.$$
意义：$n$ 个实特征频率（简正频率）。
特例：简并时需正交化。
应用：桥梁、机翼颤振分析。

**3. 简正坐标与解耦**

定义：由广义特征向量 $\mathbf a^{(r)}$ 构造简正坐标 $\eta_r$。
推导：$\boldsymbol\xi=\sum_r\eta_r\mathbf a^{(r)}$，代入得
$$\ddot\eta_r+\omega_r^2\eta_r=0.$$
意义：任意小振动是独立简谐模的叠加。
特例：两个耦合振子出现拍频。
应用：声学、分子光谱。

**4. 两耦合振子与拍**

定义：质量 $m$，弹簧 $k$，耦合弹簧 $k_c$。
推导：方程
$$m\ddot x_1=-kx_1-k_c(x_1-x_2),\quad m\ddot x_2=-kx_2+k_c(x_1-x_2).$$
对称/反对称组合 $x_\pm=x_1\pm x_2$：
$$\ddot x_+=-\frac{k}{m}x_+,\qquad \ddot x_-=-\frac{k+2k_c}{m}x_-.$$
意义：两种简正模，同相与反相。
特例：$k_c\ll k$ 时拍频 $\Delta\omega=\omega_+-\omega_-\approx k_c/\sqrt{mk}$。
应用：无线电拍频、声拍。

**5. 阻尼与受迫振动**

定义：阻尼 $\gamma$，驱动力 $F_0\cos\Omega t$。
推导：
$$\ddot x+2\gamma\dot x+\omega_0^2x=\frac{F_0}{m}\cos\Omega t.$$
稳态解 $x=A\cos(\Omega t-\delta)$，
$$A=\frac{F_0/m}{\sqrt{(\omega_0^2-\Omega^2)^2+4\gamma^2\Omega^2}},\quad
\tan\delta=\frac{2\gamma\Omega}{\omega_0^2-\Omega^2}.$$
意义：$\Omega\approx\omega_0$ 共振，幅值受 $\gamma$ 限制。
特例：无阻尼 $\gamma=0$ 共振发散。
应用：振动控制、共振放大。

**6. 参数共振**

定义：参数随时间周期变化，如 $\omega^2(t)=\omega_0^2(1+h\cos\gamma t)$。
推导：Mathieu 方程
$$\ddot x+\omega_0^2(1+h\cos\gamma t)x=0.$$
在 $\gamma\approx2\omega_0$ 时出现不稳定（参数共振）。
意义：即使无周期外力，参数调制也可致振动放大。
特例：$h=0$ 退化为简谐。
应用：秋千起振、开关电容、参量放大器。

**7. 瑞利耗散函数**

定义：$F=\tfrac12\sum b_{ij}\dot q_i\dot q_j$，广义力 $Q_i=-\partial F/\partial\dot q_i$。
推导：拉格朗日方程
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_i}-\frac{\partial L}{\partial q_i}+\frac{\partial F}{\partial\dot q_i}=0.$$
能量耗散率 $\dot E=-2F$。
意义：统一处理线性阻尼。
特例：$F=\tfrac12b\dot q^2$。
应用：结构阻尼、电路 RLC。

### 第10章 非线性动力学与稳定性

**1. 相空间与轨线**

定义：一阶系统 $\dot{\mathbf x}=\mathbf f(\mathbf x)$，$\mathbf x\in\mathbb R^n$。
推导：轨线由时间消去 $\dfrac{\mathrm dx_2}{\mathrm dx_1}=\dfrac{f_2}{f_1}$ 得；不自交（除平衡点）。
意义：几何化研究动力学。
特例：保守系统相体积守恒。
应用：摆、激光、神经元模型。

**2. 平衡点与线性稳定性**

定义：$\mathbf f(\mathbf x_0)=0$，雅可比 $J=\partial\mathbf f/\partial\mathbf x|_{\mathbf x_0}$。
推导：$\delta\dot{\mathbf x}=J\delta\mathbf x$，解 $\delta\mathbf x=\sum c_ke^{\lambda_kt}$。$\mathrm{Re}\,\lambda_k<0$ 渐近稳定，存在 $>0$ 不稳定，纯虚为中心。
意义：局部稳定性由特征值决定。
特例：二维 $\lambda_{1,2}$，$\det J<0$ 鞍点。
应用：生态、电路、控制。

**3. 李雅普诺夫稳定性**

定义：正定函数 $V(\mathbf x)$ 沿轨线 $\dot V\le0$ 则稳定；$\dot V<0$ 渐近稳定。
推导：由链式法则 $\dot V=\nabla V\cdot\mathbf f(\mathbf x)$；能量函数常作 $V$。
意义：无需解方程即可判稳。
特例：有阻尼摆 $V=\tfrac12\dot\theta^2+(g/l)(1-\cos\theta)$，$\dot V=-b\dot\theta^2$。
应用：电力系统暂态稳定。

**4. 分岔**

定义：参数变化使相图定性改变。
推导：鞍结分岔 $\dot x=\mu-x^2$，$\mu>0$ 有两个平衡点 $x=\pm\sqrt\mu$ 相撞消失；叉式分岔 $\dot x=\mu x-x^3$ 在 $\mu=0$ 由一稳变两稳。
意义：结构不稳定性。
特例：Hopf 分岔产生极限环。
应用：激光阈值、种群暴发。

**5. 极限环与庞加莱-本迪克松定理**

定义：孤立周期闭轨。
推导：二维系统中若存在有界闭不变区域且无平衡点，则内部存在闭轨。
意义：自持振荡的数学刻画。
特例：范德波尔振子 $\ddot x-\mu(1-x^2)\dot x+x=0$。
应用：心脏起搏、电子振荡器。

**6. 混沌与洛伦兹系统**

定义：对初值敏感、有界、非周期的确定性运动。
推导：洛伦兹方程
$$\dot x=\sigma(y-x),\quad \dot y=x(\rho-z)-y,\quad \dot z=xy-\beta z,$$
$\rho>24.74$ 时出现奇怪吸引子，李氏指数 $\lambda_1>0$。
意义：确定性非周期与不可预测性。
特例：Rössler、逻辑斯蒂映射 $x_{n+1}=rx_n(1-x_n)$，$r>3.57$ 混沌。
应用：天气预报、保密通信。

### 第11章 连续介质运动学与守恒方程

**1. 连续介质假设与场描述**

定义：物质性质用连续场 $\rho(\mathbf x,t),\mathbf v(\mathbf x,t)$ 描述。
假设：存在特征长度远大于分子尺度（Knudsen 数 $\ll1$）。
推导：由微观量平均定义宏观场，如 $\rho=\lim_{\Delta V\to\Delta V^*}\dfrac{\Delta m}{\Delta V}$。
意义：将多粒子系统化为场。
特例：稀薄气体需动理学。
应用：空气动力学、固体力学。

**2. 物质导数**

定义：随流体微团的时间变化率
$$\frac{\mathrm D}{\mathrm Dt}=\frac{\partial}{\partial t}+\mathbf v\cdot\nabla.$$
推导：对标量场 $f(\mathbf x,t)$ 沿轨线 $\mathbf x(t)$ 求导
$$\frac{\mathrm df}{\mathrm dt}=\frac{\partial f}{\partial t}+\frac{\partial f}{\partial x_i}\frac{\mathrm dx_i}{\mathrm dt}
=\frac{\partial f}{\partial t}+v_i\frac{\partial f}{\partial x_i}.$$
意义：局部变化 + 对流变化。
特例：定常场 $\partial/\partial t=0$。
应用：输运、污染物扩散。

**3. 雷诺输运定理**

定义：系统在时刻 $t$ 占据控制体 $V(t)$。
推导：对广延量 $B=\int_{V(t)}\rho b\,\mathrm dV$，
$$\frac{\mathrm DB}{\mathrm Dt}=\frac{\mathrm d}{\mathrm dt}\int_{V(t)}\rho b\,\mathrm dV
=\int_V\frac{\partial(\rho b)}{\partial t}\mathrm dV+\oint_{\partial V}\rho b\,\mathbf v\cdot\mathbf n\,\mathrm dA.$$
用高斯定理合并：
$$\frac{\mathrm DB}{\mathrm Dt}=\int_V\!\left[\frac{\partial(\rho b)}{\partial t}+\nabla\!\cdot\!(\rho b\mathbf v)\right]\mathrm dV.$$
意义：系统与控制体桥梁。
特例：$b=1$ 给质量守恒。
应用：动量、能量方程推导。

**4. 连续性方程**

定义：质量守恒。
推导：取 $b=1$，$B=\int\rho\,\mathrm dV$ 守恒，$\mathrm DB/\mathrm Dt=0$：
$$\int_V\!\left[\frac{\partial\rho}{\partial t}+\nabla\!\cdot\!(\rho\mathbf v)\right]\mathrm dV=0\ \Rightarrow\ \frac{\partial\rho}{\partial t}+\nabla\!\cdot\!(\rho\mathbf v)=0.$$
不可压 $\rho$ 常数化为 $\nabla\cdot\mathbf v=0$。
意义：密度变化由通量决定。
特例：定常 $\nabla\cdot(\rho\mathbf v)=0$。
应用：管道流量、地下水渗流。

**5. 动量方程（柯西方程）**

定义：微元体受面力 $\mathbf t=\boldsymbol\sigma\cdot\mathbf n$ 与体力 $\mathbf f$。
推导：对动量 $\int\rho\mathbf v\,\mathrm dV$ 用输运定理并令其等于外力
$$\int_V\!\left[\frac{\partial(\rho\mathbf v)}{\partial t}+\nabla\!\cdot\!(\rho\mathbf v\mathbf v)\right]\mathrm dV
=\oint_{\partial V}\boldsymbol\sigma\cdot\mathbf n\,\mathrm dA+\int_V\mathbf f\,\mathrm dV.$$
用高斯定理 $\oint\boldsymbol\sigma\cdot\mathbf n\,\mathrm dA=\int\nabla\cdot\boldsymbol\sigma\,\mathrm dV$，并利用连续性方程得
$$\rho\frac{\mathrm D\mathbf v}{\mathrm Dt}=\nabla\!\cdot\!\boldsymbol\sigma+\mathbf f.$$
意义：连续体动量平衡。
特例：无应力即牛顿第二定律。
应用：固体、流体统一框架。

**6. 角动量方程与应力对称**

定义：角动量守恒 $\dfrac{\mathrm D}{\mathrm Dt}\int\mathbf x\times\rho\mathbf v\,\mathrm dV=\int\mathbf x\times\mathbf f+\oint\mathbf x\times(\boldsymbol\sigma\cdot\mathbf n)$。
推导：展开后仅在无体力偶假设下得
$$\boldsymbol\sigma=\boldsymbol\sigma^{\rm T},$$
即 $\sigma_{ij}=\sigma_{ji}$。
意义：应力张量对称，独立分量由 9 降为 6。
特例：有极化体力偶的液晶等情况不成立。
应用：有限元本构。

**7. 能量方程**

定义：总能量 $\int\rho(e+\tfrac12v^2)\mathrm dV$，$e$ 比内能。
推导：利用输运定理与动量方程，得
$$\rho\frac{\mathrm De}{\mathrm Dt}=\boldsymbol\sigma:\nabla\mathbf v-\nabla\!\cdot\!\mathbf q+r,$$
其中 $\mathbf q$ 热流，$r$ 体热源。黏性耗散 $\boldsymbol\sigma:\nabla\mathbf v$。
意义：内能变化 = 变形功 + 热传导。
特例：固体小变形线性化得热弹性耦合。
应用：燃烧、气动加热。

### 第12章 应力、应变与本构关系

**1. 应力张量**

定义：面力 $\mathbf t^{(\mathbf n)}=\boldsymbol\sigma\cdot\mathbf n$。
推导：四面体微元力平衡（柯西公式）：设斜面上 $\mathbf n=n_i\mathbf e_i$，由动量平衡忽略惯性
$$\mathbf t^{(\mathbf n)}\,\mathrm dA=\sum_i\mathbf t^{(\mathbf e_i)}\,\mathrm dA_i,\quad
\mathrm dA_i=n_i\mathrm dA,$$
故 $t_j=\sigma_{ij}n_i$。
意义：任意截面应力由张量决定。
特例：静水压 $\boldsymbol\sigma=-p\mathbf I$。
应用：应力分析、有限元。

**2. 应力变换与主应力**

定义：坐标变换 $\boldsymbol\sigma'=Q\boldsymbol\sigma Q^{\rm T}$。
推导：由张量性质，主应力满足特征方程
$$\det(\boldsymbol\sigma-\sigma\mathbf I)=0\ \Rightarrow\ \sigma^3-I_1\sigma^2+I_2\sigma-I_3=0,$$
$$I_1=\mathrm{tr}\,\boldsymbol\sigma,\quad I_2=\tfrac12[(\mathrm{tr}\boldsymbol\sigma)^2-\mathrm{tr}\boldsymbol\sigma^2],\quad I_3=\det\boldsymbol\sigma.$$
意义：不变量与坐标无关。
特例：莫尔圆几何求解。
应用：强度理论、断裂。

**3. 小应变张量**

定义：位移场 $\mathbf u(\mathbf x)$，位移梯度 $\partial u_i/\partial x_j$。
推导：两点距离平方变化
$$\mathrm ds'^2-\mathrm ds^2=2\varepsilon_{ij}\mathrm dx_i\mathrm dx_j,$$
$$\varepsilon_{ij}=\tfrac12(\partial_iu_j+\partial_ju_i).$$
反对称部分 $\omega_{ij}=\tfrac12(\partial_iu_j-\partial_ju_i)$ 为刚体转动。
意义：对称应变描述纯变形。
特例：均匀拉伸 $\varepsilon_{11}=\partial u_1/\partial x_1$。
应用：应变花测量。

**4. 广义胡克定律**

定义：线弹性 $\sigma_{ij}=C_{ijkl}\varepsilon_{kl}$。
推导：各向同性要求 $C_{ijkl}=\lambda\delta_{ij}\delta_{kl}+\mu(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk})$，故
$$\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij},$$
$$\varepsilon_{ij}=\frac{1}{2\mu}\!\left(\sigma_{ij}-\frac{\lambda}{3\lambda+2\mu}\sigma_{kk}\delta_{ij}\right).$$
意义：拉梅常数 $\lambda,\mu$；$\mu$ 为剪切模量。
特例：单轴拉伸得 $E=\dfrac{\mu(3\lambda+2\mu)}{\lambda+\mu}$，$\nu=\dfrac{\lambda}{2(\lambda+\mu)}$。
应用：结构强度、地震波。

**5. 弹性常数关系**

定义：$E$ 杨氏模量、$\nu$ 泊松比、$K$ 体积模量、$G=\mu$ 剪切模量。
推导：由单轴拉伸 $\sigma_{11}=\sigma$，$\varepsilon_{11}=\sigma/E$，$\varepsilon_{22}=\varepsilon_{33}=-\nu\sigma/E$，代入广义胡克定律比较得
$$G=\frac{E}{2(1+\nu)},\qquad K=\frac{E}{3(1-2\nu)},\qquad
\lambda=\frac{E\nu}{(1+\nu)(1-2\nu)}.$$
意义：四个常数仅两个独立。
特例：不可压 $\nu=1/2$，$K\to\infty$。
应用：材料参数互换、波速 $c_p=\sqrt{(K+4G/3)/\rho}$。

**6. 平面应力与平面应变**

定义：平面应力 $\sigma_{33}=0$；平面应变 $\varepsilon_{33}=0$。
推导：平面应力下 $\varepsilon_{33}=-\dfrac{\nu}{E}(\sigma_{11}+\sigma_{22})$；平面应变下 $\sigma_{33}=\nu(\sigma_{11}+\sigma_{22})$。
统一记为
$$\varepsilon_\alpha=\frac{1}{E'}(\sigma_\alpha-\nu'\sigma_\beta),$$
平面应力 $E'=E,\nu'=\nu$；平面应变 $E'=E/(1-\nu^2),\nu'=\nu/(1-\nu)$。
意义：二维问题两大类。
特例：薄板与长坝。
应用：板壳、挡土墙。

**7. 应变能与余能**

定义：应变能密度 $w=\tfrac12\sigma_{ij}\varepsilon_{ij}$。
推导：对弹性加载路径积分
$$w=\int_0^{\varepsilon_{ij}}\sigma_{ij}\,\mathrm d\varepsilon_{ij}
=\tfrac12\lambda(\varepsilon_{kk})^2+\mu\varepsilon_{ij}\varepsilon_{ij}.$$
意义：储能与材料刚度正定。
特例：各向同性正定条件 $\mu>0,\ 3\lambda+2\mu>0$。
应用：能量法、断裂能释放率。

**8. 黏弹性本构**

定义：应力依赖应变历史。
推导：麦克斯韦模型 $\dot\varepsilon=\dfrac{\dot\sigma}{E}+\dfrac{\sigma}{\eta}$；开尔文-沃伊特模型 $\sigma=E\varepsilon+\eta\dot\varepsilon$。
意义：滞后与松弛。
特例：阶跃应变给出应力松弛 $\sigma=\sigma_0e^{-Et/\eta}$。
应用：高分子、生物组织、减振。

### 第13章 弹性力学与能量法

**1. 平衡方程与纳维方程**

定义：静力平衡 $\nabla\cdot\boldsymbol\sigma+\mathbf f=0$。
推导：将本构代入应变-位移关系：
$$\sigma_{ij}=\lambda\delta_{ij}\partial_ku_k+\mu(\partial_iu_j+\partial_ju_i),$$
$$\partial_j\sigma_{ij}=\lambda\partial_i(\partial_ku_k)+\mu\partial_j\partial_iu_j+\mu\partial_j\partial_ju_i
=(\lambda+\mu)\partial_i(\nabla\!\cdot\!\mathbf u)+\mu\nabla^2u_i.$$
故
$$(\lambda+\mu)\nabla(\nabla\!\cdot\!\mathbf u)+\mu\nabla^2\mathbf u+\mathbf f=0.$$
意义：弹性静力学基本方程（拉梅-纳维）。
特例：无体力求通解。
应用：有限元、边界元。

**2. 弹性波与拉梅解**

定义：无体力动力方程 $\rho\ddot{\mathbf u}=(\lambda+\mu)\nabla(\nabla\!\cdot\!\mathbf u)+\mu\nabla^2\mathbf u$。
推导：亥姆霍兹分解 $\mathbf u=\nabla\phi+\nabla\times\boldsymbol\psi$，$\nabla\cdot\boldsymbol\psi=0$，代入得
$$\nabla\!\left[(\lambda+2\mu)\nabla^2\phi-\rho\ddot\phi\right]+\nabla\times\!\left[\mu\nabla^2\boldsymbol\psi-\rho\ddot{\boldsymbol\psi}\right]=0,$$
故
$$\nabla^2\phi=\frac{1}{c_p^2}\ddot\phi,\qquad \nabla^2\boldsymbol\psi=\frac{1}{c_s^2}\ddot{\boldsymbol\psi},$$
$$c_p=\sqrt{\frac{\lambda+2\mu}{\rho}},\qquad c_s=\sqrt{\frac{\mu}{\rho}}.$$
意义：纵波（P）与横波（S）分离。
特例：$c_p>c_s$，地震先见 P 波。
应用：地震学、超声检测。

**3. 梁的弯曲**

定义：欧拉-伯努利梁，挠度 $w(x)$。
推导：平截面假设下应变 $\varepsilon=-\dfrac{y}{\rho}$，$\rho$ 为曲率半径；$\sigma=E\varepsilon$。弯矩
$$M=\int_A\sigma y\,\mathrm dA=-\frac{E}{\rho}\int_Ay^2\mathrm dA=-\frac{EI}{\rho}.$$
小挠度 $\dfrac1\rho\approx w''$，得
$$EI\,w''=M(x),\qquad \frac{\mathrm d^2}{\mathrm dx^2}\!\left(EI w''\right)=q(x).$$
意义：梁弯曲微分方程。
特例：简支梁集中力解为三次曲线。
应用：桥梁、机翼梁。

**4. 最小势能原理**

定义：总势能 $\Pi=\int_V w\,\mathrm dV-\int_V\mathbf f\cdot\mathbf u\,\mathrm dV-\int_{\partial V}\mathbf t\cdot\mathbf u\,\mathrm dA$。
推导：对 $\mathbf u\to\mathbf u+\delta\mathbf u$ 取一阶变分
$$\delta\Pi=\int_V\boldsymbol\sigma:\delta\boldsymbol\varepsilon\,\mathrm dV-\int\mathbf f\cdot\delta\mathbf u-\int\mathbf t\cdot\delta\mathbf u=0,$$
分部积分即得平衡方程与力边界条件。
意义：稳定平衡使势能最小。
特例：线性弹性下唯一极小。
应用：里兹法、有限元。

**5. 卡氏定理与单位载荷法**

定义：位移 $u_i=\partial U/\partial P_i$，$U$ 为应变能。
推导：对载荷 $P_i$ 给微小增量，能量增量 $\delta U=u_i\,\delta P_i$，而由互等定理与线性关系 $\delta U=\dfrac{\partial U}{\partial P_i}\delta P_i$，比较得。
意义：求指定点位移的便捷工具。
特例：莫尔积分 $u=\int\dfrac{M M_0}{EI}\mathrm dx$。
应用：桁架挠度、框架分析。

**6. 虚功原理**

定义：任意虚位移下 $\delta W_{\rm ext}=\delta U$。
推导：由能量守恒一阶变分 $\delta(\Pi)=0$ 展开即得；对刚体 $\delta U=0$。
意义：静力学统一表述。
特例：拉格朗日-达朗贝尔原理的静力版本。
应用：机构力分析、极限分析。

**7. 应力集中与断裂能**

定义：孔边应力集中系数 $K_t=\sigma_{\max}/\sigma_{\rm nom}$。
推导：圆孔在单向拉伸下由 Airy 应力函数求解得
$$\sigma_{\theta\theta}\big|_{\theta=\pm\pi/2}=3\sigma_\infty\ \Rightarrow\ K_t=3.$$
格里菲斯能量判据：裂纹扩展释放弹性能 $\mathrm dU_e$ 与表面能 $2\gamma\,\mathrm da$ 平衡
$$\frac{\partial}{\partial a}(U_e-2\gamma a)=0\ \Rightarrow\ G=\frac{\partial U_e}{\partial a}=2\gamma.$$
意义：应力集中与断裂能量门槛。
特例：无限大板中心裂纹 $K_I=\sigma\sqrt{\pi a}$。
应用：结构完整性、损伤容限。

### 第14章 流体力学与纳维-斯托克斯方程

**1. 理想流体与欧拉方程**

定义：无黏、无热传导流体。
推导：柯西方程中 $\boldsymbol\sigma=-p\mathbf I$，故
$$\rho\frac{\mathrm D\mathbf v}{\mathrm Dt}=-\nabla p+\mathbf f.$$
意义：忽略黏性的动量方程。
特例：静流体 $\mathbf v=0$ 得 $\nabla p=\mathbf f$。
应用：机翼升力（远场）、大气动力学。

**2. 伯努利方程**

定义：定常、无黏、不可压、沿流线。
推导：欧拉方程沿流线投影
$$\mathbf v\cdot\nabla\!\left(\tfrac12v^2\right)=-\frac{1}{\rho}\mathbf v\cdot\nabla p+\frac{1}{\rho}\mathbf v\cdot\mathbf f.$$
若有势体力 $\mathbf f=-\rho\nabla\Phi$，则
$$\mathbf v\cdot\nabla\!\left(\tfrac12v^2+\frac{p}{\rho}+\Phi\right)=0,$$
沿流线
$$\tfrac12v^2+\frac{p}{\rho}+\Phi=\text{const}.$$
意义：速度大处压强小。
特例：重力场 $\Phi=gz$：$\tfrac12v^2+p/\rho+gz=\text{const}$。
应用：文丘里流量计、皮托管。

**3. 连续性方程与流函数**

定义：不可压二维流 $\partial_xu+\partial_yv=0$。
推导：令 $u=\dfrac{\partial\psi}{\partial y},\ v=-\dfrac{\partial\psi}{\partial x}$，则守恒自动满足。
意义：流函数等值线为流线。
特例：势流 $\nabla^2\psi=0$。
应用：机翼绕流、地下水。

**4. 势流与复势**

定义：无旋流 $\mathbf v=\nabla\phi$，不可压 $\nabla^2\phi=0$。
推导：引入复势 $F(z)=\phi+i\psi$，$z=x+iy$。柯西-黎曼方程给出
$$u-v\,\mathrm i=\frac{\mathrm dF}{\mathrm dz}.$$
意义：保角变换可解复杂边界。
特例：圆柱绕流、偶极子、点涡。
应用：气动外形、翼型设计（儒可夫斯基）。

**5. 牛顿黏性定律与纳维-斯托克斯方程**

定义：牛顿流体 $\sigma_{ij}=-p\delta_{ij}+2\mu\varepsilon_{ij}-\tfrac23\mu\varepsilon_{kk}\delta_{ij}$。
推导：将本构代入柯西动量方程：
$$\partial_j\sigma_{ij}=-\partial_ip+\mu\partial_j(\partial_iu_j+\partial_ju_i)-\tfrac23\mu\partial_i(\partial_ku_k).$$
不可压时最后一项为零，得
$$\boxed{\ \rho\!\left(\frac{\partial\mathbf v}{\partial t}+\mathbf v\!\cdot\!\nabla\mathbf v\right)=-\nabla p+\mu\nabla^2\mathbf v+\mathbf f\ }$$
意义：流体运动的基本方程。
特例：$\mu=0$ 退化为欧拉方程；$\mathbf v=0$ 退化为静力平衡。
应用：湍流、润滑、血流、气候。

**6. 量纲分析与雷诺数**

定义：惯性力/黏性力 $Re=\rho UL/\mu=UL/\nu$。
推导：N-S 方程无量纲化，令 $\mathbf v=U\mathbf v^*,\ \mathbf x=L\mathbf x^*,\ t=(L/U)t^*$：
$$\frac{\partial\mathbf v^*}{\partial t^*}+\mathbf v^*\!\cdot\!\nabla^*\mathbf v^*=-\nabla^*p^*+\frac{1}{Re}\nabla^{*2}\mathbf v^*.$$
意义：$Re$ 决定流态（层流/湍流）。
特例：$Re\ll1$ 斯托克斯流，$Re\gg1$ 高雷诺数。
应用：风洞缩比、管道设计。

**7. 泊肃叶流动**

定义：圆管层流，定常、轴对称。
推导：轴向导 $z$，方程化为
$$\frac{1}{r}\frac{\mathrm d}{\mathrm dr}\!\left(r\frac{\mathrm dv_z}{\mathrm dr}\right)=\frac{1}{\mu}\frac{\mathrm dp}{\mathrm dz}.$$
积分并施以 $v_z(R)=0$ 与轴对称为有限：
$$v_z(r)=\frac{1}{4\mu}\!\left(-\frac{\mathrm dp}{\mathrm dz}\right)(R^2-r^2).$$
流量
$$Q=\int_0^R v_z\,2\pi r\,\mathrm dr=\frac{\pi R^4}{8\mu}\!\left(-\frac{\mathrm dp}{\mathrm dz}\right).$$
意义：流量正比于半径四次方。
特例：哈根-泊肃叶定律。
应用：毛细管黏度计、血管流量。

**8. 边界层与普朗特方程**

定义：壁面附近黏性显著薄层，厚度 $\delta(x)$。
推导：对 N-S 做量级分析：$u\sim U,\ x\sim L,\ y\sim\delta$，连续性得 $v\sim U\delta/L$。当 $Re\gg1$ 时
$$u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=-\frac{1}{\rho}\frac{\mathrm dp}{\mathrm dx}+\nu\frac{\partial^2u}{\partial y^2}.$$
意义：黏性效应局限薄层，外流近似无黏。
特例：平板布拉修斯解。
应用：摩阻、分离、失速。

**9. 斯托克斯阻力**

定义：$Re\ll1$ 小球缓慢运动。
推导：斯托克斯流解给出
$$\mathbf F_d=-6\pi\mu R\mathbf v.$$
与重力、浮力平衡得沉降末速
$$v_t=\frac{2R^2(\rho_s-\rho_f)g}{9\mu}.$$
意义：低雷诺数阻力线性于速度与半径。
特例：米利肯油滴实验。
应用：沉降、气溶胶、细胞分离。

### 第15章 涡动力学与波动

**1. 涡量与环量**

定义：$\boldsymbol\omega=\nabla\times\mathbf v$，环量 $\Gamma=\oint_C\mathbf v\cdot\mathrm d\mathbf l$。
推导：斯托克斯定理
$$\Gamma=\oint_C\mathbf v\cdot\mathrm d\mathbf l=\int_S(\nabla\times\mathbf v)\cdot\mathbf n\,\mathrm dA=\int_S\boldsymbol\omega\cdot\mathbf n\,\mathrm dA.$$
意义：环量为涡通量。
特例：无旋流 $\boldsymbol\omega=0$，$\Gamma=0$。
应用：升力环量、旋涡识别。

**2. 涡量方程**

定义：对 N-S 取旋度。
推导：利用恒等式
$$\nabla\times(\mathbf v\cdot\nabla\mathbf v)=(\mathbf v\cdot\nabla)\boldsymbol\omega-(\boldsymbol\omega\cdot\nabla)\mathbf v+(\nabla\!\cdot\!\mathbf v)\boldsymbol\omega,$$
对不可压：
$$\frac{\partial\boldsymbol\omega}{\partial t}+\mathbf v\!\cdot\!\nabla\boldsymbol\omega
=(\boldsymbol\omega\!\cdot\!\nabla)\mathbf v+\nu\nabla^2\boldsymbol\omega.$$
意义：涡量被对流、拉伸、黏性扩散。
特例：二维 $\boldsymbol\omega\cdot\nabla\mathbf v=0$，仅有扩散。
应用：涡环、湍流级联。

**3. 开尔文环量定理**

定义：理想正压流体，体力有势。
推导：由涡量方程取 $\nu=0$ 得
$$\frac{\mathrm D\boldsymbol\omega}{\mathrm Dt}=(\boldsymbol\omega\!\cdot\!\nabla)\mathbf v.$$
故
$$\frac{\mathrm D\Gamma}{\mathrm Dt}=\int_S\!\left(\frac{\mathrm D\boldsymbol\omega}{\mathrm Dt}-\boldsymbol\omega\!\cdot\!\nabla\mathbf v\right)\!\cdot\!\mathrm d\mathbf A=0.$$
意义：理想流体中环量与涡管强度守恒。
特例：有黏性时环量衰减。
应用：翼型起动涡、超流。

**4. 亥姆霍兹涡定理**

定义：涡管强度守恒、涡线随流体运动。
推导：由 $\nabla\cdot\boldsymbol\omega=0$ 与开尔文定理，涡管任一截面通量相同；涡线由物质点组成。
意义：涡管不能中断，只能成环或延伸至边界。
特例：涡环稳定迁移。
应用：涡环推进、龙卷风。

**5. 声波方程**

定义：小扰动 $\rho=\rho_0+\rho',\ p=p_0+p',\ \mathbf v=\mathbf v'$。
推导：线性化连续性 $\partial_t\rho'+\rho_0\nabla\!\cdot\!\mathbf v=0$，欧拉方程 $\rho_0\partial_t\mathbf v=-\nabla p'$。由等熵 $p'=c^2\rho'$，取散度与时间导数消元：
$$\nabla^2p'=\frac{1}{c^2}\frac{\partial^2p'}{\partial t^2},\qquad c=\sqrt{\frac{\partial p}{\partial\rho}\Big|_s}.$$
意义：声速由介质可压缩性决定。
特例：理想气体 $c=\sqrt{\gamma k_BT/m}$。
应用：声纳、超声、气动噪声。

**6. 弹性波与色散**

定义：波速依赖波长称色散。
推导：试解 $u=e^{i(kx-\omega t)}$ 代入波动方程得 $\omega=ck$；对梁弯曲方程 $EIw''''+\rho Aw_{tt}=0$ 得
$$\omega=\sqrt{\frac{EI}{\rho A}}k^2.$$
意义：色散使波形变化、群速不同。
特例：相速 $v_p=\omega/k$，群速 $v_g=\mathrm d\omega/\mathrm dk$；无色散 $v_g=v_p$。
应用：地震波、光纤、等离子体。

**7. 波的能量与能流**

定义：能量密度 $e$，能流密度（强度） $\mathbf I=p'\mathbf v'$。
推导：对声波
$$e=\frac{1}{2}\rho_0v'^2+\frac{1}{2}\frac{p'^2}{\rho_0c^2},\qquad \mathbf I=\overline{p'\mathbf v'}=\frac{p_{\max}^2}{2\rho_0c}\hat{\mathbf k}.$$
意义：强度与幅值平方成正比。
特例：球面波强度 $\propto1/r^2$。
应用：声级、地震震级。

**8. 多普勒效应**

定义：源或观察者运动引起频率变化。
推导：设源速 $v_s$、观察者速 $v_o$，波速 $c$，沿连线：
$$f'=\frac{c\pm v_o}{c\mp v_s}f.$$
意义：相对运动压缩/拉伸波长。
特例：马赫锥当 $v_s>c$。
应用：雷达测速、红移。

### 第16章 相对论力学

**1. 洛伦兹变换**

定义：惯性系间保持光速不变与间隔不变的线性变换。
推导：设 $x'=\gamma(x-vt),\ t'=\gamma(t-vx/c^2)$，由 $x^2-c^2t^2$ 不变可定
$$\gamma=\frac{1}{\sqrt{1-v^2/c^2}}.$$
意义：时空混合，同时性相对。
特例：$v\ll c$ 退化为伽利略变换。
应用：粒子寿命、GPS 修正。

**2. 相对论动量与能量**

定义：四维动量 $p^\mu=(E/c,\mathbf p)$。
推导：四维速度 $u^\mu=\gamma(c,\mathbf v)$，不变质量 $m$ 给出
$$\mathbf p=\gamma m\mathbf v,\qquad E=\gamma mc^2.$$
由 $p^\mu p_\mu=m^2c^2$：
$$E^2-p^2c^2=m^2c^4.$$
意义：质量与能量等价。
特例：$v=0$ 时 $E=mc^2$；$v\to c$ 需无穷能量。
应用：核反应、加速器。

**3. 相对论动能与力**

定义：$K=E-mc^2$。
推导：泰勒展开
$$\gamma mc^2=mc^2+\tfrac12mv^2+\frac{3}{8}\frac{mv^4}{c^2}+\cdots,$$
故 $K=(\gamma-1)mc^2$。力
$$\mathbf F=\frac{\mathrm d\mathbf p}{\mathrm dt}=\gamma m\mathbf a+\gamma^3m\frac{\mathbf v\!\cdot\!\mathbf a}{c^2}\mathbf v.$$
意义：惯性随速度增大且方向相关。
特例：$v\ll c$ 时 $K\to\tfrac12mv^2$。
应用：电子加速、核结合能。

**4. 相对论拉格朗日函数**

定义：自由粒子作用量 $S=-\displaystyle\int mc^2\,\mathrm d\tau$。
推导：$\mathrm d\tau=\mathrm dt/\gamma=\mathrm dt\sqrt{1-v^2/c^2}$，故
$$L=-mc^2\sqrt{1-\frac{v^2}{c^2}}.$$
广义动量
$$p=\frac{\partial L}{\partial v}=\frac{mv}{\sqrt{1-v^2/c^2}},$$
哈密顿
$$H=pv-L=\sqrt{p^2c^2+m^2c^4}.$$
意义：相对论动力学可由最小作用量导出。
特例：$v\ll c$ 得经典 $L=\tfrac12mv^2$（差常数 $mc^2$）。
应用：场论、高能物理。

**5. 四维形式与协变动力学**

定义：$x^\mu=(ct,\mathbf x)$，$\partial_\mu=(\partial_{ct},\nabla)$，度规 $\eta=\mathrm{diag}(1,-1,-1,-1)$。
推导：洛伦兹变换 $x'^\mu=\Lambda^\mu{}_\nu x^\nu$，$\eta_{\mu\nu}\Lambda^\mu{}_\alpha\Lambda^\nu{}_\beta=\eta_{\alpha\beta}$。电磁场张量 $F^{\mu\nu}=\partial^\mu A^\nu-\partial^\nu A^\mu$，运动方程
$$m\frac{\mathrm du^\mu}{\mathrm d\tau}=qF^{\mu\nu}u_\nu.$$
意义：物理定律洛伦兹协变。
特例：$\mu=0$ 给能量方程，$\mu=i$ 给洛伦兹力。
应用：电动力学、粒子物理。

**6. 质能关系与结合能**

定义：静止能量 $E_0=mc^2$，质量亏损 $\Delta m$。
推导：核反应前后四动量守恒，静能差
$$Q=(m_{\rm initial}-m_{\rm final})c^2.$$
意义：结合能来源。
特例：太阳核聚变每秒约 $4\times10^{26}\,\mathrm J$。
应用：核能、恒星演化。

**7. 相对论多普勒效应**

定义：沿视线相对速度 $v$。
推导：四维波矢 $k^\mu=(\omega/c,\mathbf k)$ 在洛伦兹变换下给出
$$\omega'=\gamma\omega(1-v\cos\theta/c).$$
意义：红移/蓝移与横向多普勒。
特例：纵向 $\theta=0$ 得 $\omega'=\omega\sqrt{(1-\beta)/(1+\beta)}$；横向 $\theta=\pi/2$ 得时间膨胀红移 $\omega'=\omega/\gamma$。
应用：宇宙膨胀、激光冷却。

---
⬆ 上级领域：[[物质科学（物理·化学·生命）]]
