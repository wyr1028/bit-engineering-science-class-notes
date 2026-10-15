---
课程名称: 大学物理AⅠ
学分: 4
总学时: 64
考核方式: 考试
课程属性: 必修
课程性质: 公共基础课程
开课单位: 物理学院
领域: 物质科学（物理·化学·生命）
tags:
  - 课程
  - 领域/物质科学（物理·化学·生命）
---

# 大学物理AⅠ


### 第1章 质点运动学

**知识点1：参考系、质点模型与位矢**
定义：研究物体机械运动时，被选作标准（视为静止）的物体叫参考系；当物体线度远小于所研究运动的特征尺度时，可抽象为质点。
假设：空间为欧氏空间，时间绝对、均匀流逝；可建立直角坐标系 $Oxyz$。
位置矢量：$\vec r(t)=x(t)\hat i+y(t)\hat j+z(t)\hat k$。
推导轨迹：由参数式 $x=x(t),\,y=y(t),\,z=z(t)$ 消去 $t$，得轨迹方程 $F(x,y,z)=0$。
位移与路程：$\Delta\vec r=\vec r(t+\Delta t)-\vec r(t)$；$\Delta s=\int_{t}^{t+\Delta t}|\mathrm d\vec r|$。
意义：位矢相对参考系原点定义；位移是矢量，路程是标量。
极限：$\Delta t\to 0$ 时 $|\Delta\vec r|\to \Delta s$，即 $\mathrm ds=|\mathrm d\vec r|$。
应用：行星、导弹质心运动均先按质点处理。

**知识点2：速度、加速度及其坐标表示**
定义：$\vec v=\dfrac{\mathrm d\vec r}{\mathrm dt}$，$\vec a=\dfrac{\mathrm d\vec v}{\mathrm dt}=\dfrac{\mathrm d^2\vec r}{\mathrm dt^2}$。
推导：对 $\vec r=x\hat i+y\hat j+z\hat k$ 求导，单位矢量不变，
$$\vec v=\dot x\hat i+\dot y\hat j+\dot z\hat k=v_x\hat i+v_y\hat j+v_z\hat k,\qquad \vec a=\ddot x\hat i+\ddot y\hat j+\ddot z\hat k.$$
速率 $v=|\vec v|=\sqrt{v_x^2+v_y^2+v_z^2}$。
由 $\vec v=\dfrac{\mathrm ds}{\mathrm dt}\hat\tau$，再求导
$$\vec a=\frac{\mathrm dv}{\mathrm dt}\hat\tau+v\frac{\mathrm d\hat\tau}{\mathrm dt}.$$
因 $\hat\tau\cdot\hat\tau=1\Rightarrow \dfrac{\mathrm d\hat\tau}{\mathrm dt}\perp\hat\tau$，定义曲率 $\kappa=\left|\dfrac{\mathrm d\hat\tau}{\mathrm ds}\right|$，得
$$\vec a=\frac{\mathrm dv}{\mathrm dt}\hat\tau+\frac{v^2}{R}\hat n,\qquad R=\frac{1}{\kappa}.$$
意义：切向加速度改变速率大小，法向加速度改变速度方向。
特例：匀速直线 $a=0$；匀速率圆周运动 $a_t=0,\ a_n=v^2/R$。
应用：已知 $x(t)$ 求 $v,a$，或由 $a$ 积分求运动（需初条件）。

**知识点3：匀变速直线运动与微积分求解**
定义：$a=$ 常矢量且与 $v$ 共线。
推导：由 $a=\dfrac{\mathrm dv}{\mathrm dt}$，$\int_{v_0}^{v}\mathrm dv=\int_0^t a\,\mathrm dt$，得 $v=v_0+at$；
再由 $v=\dfrac{\mathrm dx}{\mathrm dt}$，$\int_{x_0}^{x}\mathrm dx=\int_0^t(v_0+at)\,\mathrm dt$，得
$$x=x_0+v_0t+\tfrac12 at^2.$$
消 $t$：$t=(v-v_0)/a$，代入得 $v^2-v_0^2=2a(x-x_0)$。
意义：三式为运动学基本关系。
极限：$a\to0$ 退化为匀速直线。
应用：自由落体、斜坡滑动。

**知识点4：抛体运动**
假设：忽略空气阻力，$a=-g\hat j$，初速 $v_0$ 与水平成 $\theta$。
推导：
$$\ddot x=0\Rightarrow x=(v_0\cos\theta)t,\qquad \ddot y=-g\Rightarrow y=(v_0\sin\theta)t-\tfrac12 gt^2.$$
消 $t=\dfrac{x}{v_0\cos\theta}$ 得轨迹
$$y=x\tan\theta-\frac{gx^2}{2v_0^2\cos^2\theta}.$$
射程：令 $y=0$，$x=R=\dfrac{v_0^2\sin 2\theta}{g}$；最大高度 $H=\dfrac{v_0^2\sin^2\theta}{2g}$。
意义：轨迹为抛物线（重力场中）。
特例：$\theta=45^\circ$ 射程最大；$\theta$ 与 $90^\circ-\theta$ 射程相同。
应用：无阻力弹道计算。

**知识点5：圆周运动的角量描述**
定义：角位置 $\theta$，角速度 $\omega=\dfrac{\mathrm d\theta}{\mathrm dt}$，角加速度 $\beta=\dfrac{\mathrm d\omega}{\mathrm dt}$。
推导：弧长 $s=R\theta$，两边对 $t$ 求导 $v=R\omega$；再求导
$$a_t=\frac{\mathrm dv}{\mathrm dt}=R\frac{\mathrm d\omega}{\mathrm dt}=R\beta,\qquad a_n=\frac{v^2}{R}=R\omega^2.$$
匀角加速：$\omega=\omega_0+\beta t,\ \theta=\theta_0+\omega_0t+\tfrac12\beta t^2,\ \omega^2=\omega_0^2+2\beta(\theta-\theta_0)$。
意义：角量与线量由半径联系。
特例：$\beta=0$ 匀速圆周运动。
应用：飞轮、齿轮传动。

**知识点6：相对运动与伽利略变换**
定义：设 $S'$ 相对 $S$ 以 $\vec R_0(t)$ 平动，则
$$\vec r=\vec r'+\vec R_0,\qquad \vec v=\vec v'+\vec V_0,\qquad \vec a=\vec a'+\vec A_0.$$
推导：对位矢关系逐次求导，$\vec V_0=\dfrac{\mathrm d\vec R_0}{\mathrm dt}$，$\vec A_0=\dfrac{\mathrm d\vec V_0}{\mathrm dt}$。
意义：绝对量 = 相对量 + 牵连量。
特例：$S'$ 匀速平动（$\vec A_0=0$）时 $\vec a=\vec a'$，牛顿定律形式不变（伽利略相对性），空间间隔与时间间隔不变：$\Delta\vec r=\Delta\vec r',\ \Delta t=\Delta t'$。
应用：船上抛球、雨滴相对车速。

**知识点7：极坐标下的速度与加速度**
定义：$\vec r=r\hat e_r$，$\hat e_r=\cos\theta\,\hat i+\sin\theta\,\hat j,\ \hat e_\theta=-\sin\theta\,\hat i+\cos\theta\,\hat j$。
推导：$\dfrac{\mathrm d\hat e_r}{\mathrm dt}=\dfrac{\mathrm d\hat e_r}{\mathrm d\theta}\dfrac{\mathrm d\theta}{\mathrm dt}=\dot\theta\hat e_\theta$，同理 $\dfrac{\mathrm d\hat e_\theta}{\mathrm dt}=-\dot\theta\hat e_r$。
$$\vec v=\frac{\mathrm d}{\mathrm dt}(r\hat e_r)=\dot r\hat e_r+r\dot\theta\hat e_\theta.$$
$$\vec a=(\ddot r-r\dot\theta^2)\hat e_r+(r\ddot\theta+2\dot r\dot\theta)\hat e_\theta.$$
意义：径向加速度含离心项 $-r\dot\theta^2$；横向含科里奥利型 $2\dot r\dot\theta$。
特例：圆周运动 $r=R$ 时 $a_r=-R\omega^2,\ a_\theta=R\beta$，与知识点5一致。
应用：有心力问题、行星轨道。

### 第2章 牛顿运动定律与动力学

**知识点1：牛顿三定律与惯性系**
内容：第一定律（惯性定律）存在惯性系，质点不受力时保持静止或匀速直线运动；第二定律 $\vec F=\dfrac{\mathrm d\vec p}{\mathrm dt}=m\vec a$（$m$ 不变）；第三定律 $\vec F_{12}=-\vec F_{21}$。
说明：第二定律是矢量方程，可分解为三个分量方程；适用于惯性系、宏观低速。
推导弹簧振子：取 $F=-kx$，则
$$m\ddot x=-kx\Rightarrow \ddot x+\omega^2x=0,\quad \omega=\sqrt{k/m}.$$
其解 $x=A\cos(\omega t+\varphi)$，代入验证 $-\omega^2A\cos+\omega^2A\cos=0$。
意义：力是改变运动状态的原因，不是维持运动的原因。
特例：$\vec F=0$ 退化为惯性定律。
应用：连接体、斜面、绳张力问题。

**知识点2：变质量方程与火箭方程**
假设：系统受外力 $\vec F_{\text{ext}}$，质量连续变化且抛射物相对主体速度为 $\vec u$。
推导：取 $t$ 时系统（主体 $m$，速度 $\vec v$），$t+\mathrm dt$ 时主体 $m+\mathrm dm$（$\mathrm dm<0$）速度 $\vec v+\mathrm d\vec v$，喷出质量 $-\mathrm dm$ 速度 $\vec v+\vec u$。
动量守恒（无外力）：
$$m\vec v=(m+\mathrm dm)(\vec v+\mathrm d\vec v)+(-\mathrm dm)(\vec v+\vec u).$$
展开忽略二阶小量：$0=m\,\mathrm d\vec v-\vec u\,\mathrm dm\Rightarrow m\dfrac{\mathrm d\vec v}{\mathrm dt}=\vec u\dfrac{\mathrm dm}{\mathrm dt}$。
有外力时 $m\dfrac{\mathrm d\vec v}{\mathrm dt}=\vec F_{\text{ext}}+\vec u\dfrac{\mathrm dm}{\mathrm dt}$。
取 $\vec u=-v_e\hat v$（向后喷气），一维：$m\dfrac{\mathrm dv}{\mathrm dt}=F_{\text{ext}}-v_e\dfrac{\mathrm dm}{\mathrm dt}$；无外力时 $\mathrm dv=-v_e\dfrac{\mathrm dm}{m}$。
积分：$\Delta v=v_e\ln\dfrac{m_0}{m_f}$（齐奥尔科夫斯基公式）。
意义：推力来自质量喷出，$\vec u\,\mathrm dm/\mathrm dt$ 称反推力。
特例：$m_f\to0$ 时 $\Delta v\to\infty$（理想）。
应用：火箭、喷气推进。

**知识点3：微分方程型动力学问题**
方法：由 $\vec F=m\ddot{\vec r}$ 写出二阶方程，分离变量积分两次，用初条件定常数。
推导例（空气阻力正比速度，一维下落）：$m\dot v=mg-kv$。
分离：$\dfrac{\mathrm dv}{mg-kv}=\dfrac{\mathrm dt}{m}$，积分
$$-\frac1k\ln(mg-kv)\Big|_{v_0}^{v}=\frac tm\Rightarrow v=\frac{mg}{k}+\left(v_0-\frac{mg}{k}\right)e^{-kt/m}.$$
终端速度：$t\to\infty,\ v\to mg/k$。
意义：阻尼使速度趋近极限。
特例：$k\to0$ 退化为 $v=gt$（取 $v_0=0$）。
应用：雨滴、降落伞。

**知识点4：约束力与绳/绳张力**
定义：约束力是限制物体运动方向的正压力、张力、摩擦力。
光滑面约束：正压力 $N\perp$ 面，不做功。
推导圆锥摆：绳长 $l$，张角 $\theta$，竖直平衡 $T\cos\theta=mg$，水平向心 $T\sin\theta=m\dfrac{v^2}{R}$，$R=l\sin\theta$。
联立：$\tan\theta=\dfrac{v^2}{gl\sin\theta}\Rightarrow v=\sqrt{gl\sin\theta\tan\theta}$，$T=\dfrac{mg}{\cos\theta}$。
意义：约束力由运动状态决定。
特例：$\theta\to0$ 时 $v\to0$。
应用：火车转弯、圆锥摆。

**知识点5：动量定理**
定义：冲量 $\vec I=\int_{t_1}^{t_2}\vec F\,\mathrm dt$，动量 $\vec p=m\vec v$。
推导：由 $\vec F=\dfrac{\mathrm d\vec p}{\mathrm dt}$，$\int_{t_1}^{t_2}\vec F\,\mathrm dt=\vec p_2-\vec p_1$。
平均力 $\bar{\vec F}=\dfrac{\Delta\vec p}{\Delta t}$。
意义：冲量等于动量增量。
特例：碰撞中 $\Delta t$ 小、$\bar F$ 大。
应用：碰撞、缓冲。

**知识点6：质心与质心运动定理**
定义：$\vec r_C=\dfrac{\sum m_i\vec r_i}{\sum m_i}$，连续体 $\vec r_C=\dfrac1M\int\vec r\,\mathrm dm$。
推导：$M\vec r_C=\sum m_i\vec r_i$，求导 $M\vec v_C=\sum m_i\vec v_i=\vec p$；再求导 $M\vec a_C=\sum\vec F_i^{\text{ext}}=\vec F^{\text{ext}}$（内力成对抵消）。
意义：系统平动可用质心等效，内力不影响质心运动。
特例：$F^{\text{ext}}=0$ 时质心匀速。
应用：爆炸、人船模型。

**知识点7：功、动能定理**
定义：$\mathrm dW=\vec F\cdot\mathrm d\vec r$，$W=\int\vec F\cdot\mathrm d\vec r$。
推导：牛顿第二定律两边点乘 $\mathrm d\vec r=\vec v\,\mathrm dt$：
$$m\frac{\mathrm d\vec v}{\mathrm dt}\cdot\vec v\,\mathrm dt=m\vec v\cdot\mathrm d\vec v=\mathrm d\left(\tfrac12 mv^2\right).$$
积分得 $W=\tfrac12mv_2^2-\tfrac12mv_1^2$。
意义：合力功等于动能增量。
特例：$W=0$ 时速率不变（如向心力）。
应用：变力做功、能量法解题。

**知识点8：保守力、势能与机械能守恒**
定义：做功与路径无关、沿闭合路径做功为零的力为保守力；$E_p$ 满足 $W_{\text{保}}=-\Delta E_p$。
推导：一维 $F=-\dfrac{\mathrm dE_p}{\mathrm dx}$；三维 $\vec F=-\nabla E_p$。
由动能定理，合外力的功中保守力部分 $W_c=-\Delta E_p$，非保守功 $W_{nc}=\Delta E_k$，故
$$W_{nc}=\Delta E_k+\Delta E_p=\Delta E.$$
当 $W_{nc}=0$ 时 $E_k+E_p=\text{const}$。
意义：能量转化与守恒。
特例：重力 $E_p=mgh$；弹簧 $E_p=\frac12kx^2$；引力 $E_p=-GMm/r$。
应用：能量法求速率、势能曲线。

**知识点9：势能曲线与平衡稳定性**
定义：由 $E_p(x)$ 分析运动范围与平衡点。
推导：平衡条件 $F=-\dfrac{\mathrm dE_p}{\mathrm dx}=0$；稳定性由二阶导判断：$\dfrac{\mathrm d^2E_p}{\mathrm dx^2}>0$ 稳定，$<0$ 不稳定，$=0$ 需更高阶。
小振动频率：在极小点展开 $E_p\approx E_p(x_0)+\tfrac12E_p''(x_0)(x-x_0)^2$，$k=E_p''$，$\omega=\sqrt{E_p''/m}$。
意义：势能极小对应稳定平衡。
特例：谐振子 $E_p\propto x^2$。
应用：分子键、轨道稳定性。

### 第3章 非惯性系与惯性力

**知识点1：平动非惯性系**
定义：相对惯性系加速平动的参考系。
推导：伽利略变换 $\vec a=\vec a'+\vec A$，代入 $\vec F=m\vec a$ 得
$$m\vec a'=\vec F-m\vec A=\vec F+\vec F_{\text{惯}},\qquad \vec F_{\text{惯}}=-m\vec A.$$
意义：在非惯性系中引入惯性力可沿用牛顿第二定律形式。
特例：$\vec A=0$ 退化为惯性系。
应用：加速车厢内单摆、电梯超重失重。

**知识点2：转动参考系中的惯性力**
推导：设 $S'$ 相对 $S$ 以 $\vec\omega$ 转动。位矢在 $S$ 与 $S'$ 的导数关系
$$\left(\frac{\mathrm d\vec A}{\mathrm dt}\right)_S=\left(\frac{\mathrm d\vec A}{\mathrm dt}\right)_{S'}+\vec\omega\times\vec A.$$
坐标速度 $\vec v=\vec v'+\vec V_0+\vec\omega\times\vec r'$。
对速度再用一次上式求加速度：
$$\vec a=\vec a'+2\vec\omega\times\vec v'+\vec\omega\times(\vec\omega\times\vec r')+\dot{\vec\omega}\times\vec r'.$$
代入得
$$\vec F+m\vec a'=\vec F-2m\vec\omega\times\vec v'-m\vec\omega\times(\vec\omega\times\vec r')-m\dot{\vec\omega}\times\vec r'.$$
意义：科里奥利力 $-2m\vec\omega\times\vec v'$、离心力 $-m\vec\omega\times(\vec\omega\times\vec r')$。
特例：$\vec\omega=0$ 退化为平动非惯性系。
应用：傅科摆、信风、地球自转偏向。

**知识点3：离心力与地球自转效应**
推导：$\vec\omega\perp\vec r'$ 分量时离心加速度 $\omega^2r'$ 指向外，离心力 $m\omega^2r'$。
重力修正：$\vec g_{\text{eff}}=\vec g-\vec\omega\times(\vec\omega\times\vec r)$，纬度 $\lambda$ 处
$$g(\lambda)\approx g_0\left(1-\frac{\omega^2R}{g_0}\cos^2\lambda\right).$$
意义：赤道有效重力略小、地球略扁。
特例：两极 $\lambda=90^\circ$ 无离心项。
应用：重力测量、地球形状。

**知识点4：惯性离心势能与旋转液面**
推导：距转轴 $\rho$ 处离心力 $F_c=m\omega^2\rho$，可写为 $F_c=-\dfrac{\partial U_c}{\partial\rho}$，积分得 $U_c=-\tfrac12m\omega^2\rho^2$。旋转液体表面为等压面，$\mathrm dp=\rho_g(-g\,\mathrm dz+\omega^2\rho\,\mathrm d\rho)=0$，积分
$$z=\frac{\omega^2\rho^2}{2g}+\text{const}.$$
意义：表面为旋转抛物面。特例：$\omega\to0$ 时为水平面。应用：离心分离、天文镜面。

**知识点5：科里奥利力与傅科摆进动**
推导：小摆水平位移 $(x,y)$，铅直分量 $\Omega_z=\Omega\sin\lambda$，运动方程
$$\ddot x-2\Omega_z\dot y+\omega_0^2x=0,\qquad \ddot y+2\Omega_z\dot x+\omega_0^2y=0.$$
令 $u=x+\mathrm iy$，合并为 $\ddot u+2\mathrm i\Omega_z\dot u+\omega_0^2u=0$，其解含因子 $e^{-\mathrm i\Omega_z t}$，故振动面以角速度 $\Omega_z$ 进动，进动周期 $T=2\pi/(\Omega\sin\lambda)$。
意义：直接验证地球自转。特例：赤道 $\lambda=0$ 不进动；两极最快。应用：傅科摆、陀螺罗盘。

**知识点6：转动系的拉格朗日量与有效势**
推导：粒子在角速度 $\vec\Omega$ 的转动系中，惯性系动能
$$T=\tfrac12m\left[\dot r^2+r^2(\dot\theta+\Omega)^2\right].$$
展开后含耦合项 $mr^2\dot\theta\Omega$，对应科里奥利力。中心力下定义有效势
$$V_{\text{eff}}(r)=V(r)+\frac{L_z^2}{2mr^2}.$$
由 $V_{\text{eff}}'(r)=0$ 求圆轨道，由 $V_{\text{eff}}''(r)>0$ 判稳。
意义：把离心效应并入势能。特例：$L_z=0$ 为径向运动。应用：轨道稳定性、旋转星系。

### 第4章 动量、角动量与守恒定律

**知识点1：角动量与力矩**
定义：$\vec L=\vec r\times\vec p$，$\vec M=\vec r\times\vec F$。
推导：
$$\frac{\mathrm d\vec L}{\mathrm dt}=\frac{\mathrm d\vec r}{\mathrm dt}\times\vec p+\vec r\times\frac{\mathrm d\vec p}{\mathrm dt}=\vec v\times m\vec v+\vec r\times\vec F=\vec M.$$
意义：力矩是角动量变化率。
特例：$\vec M=0$ 时 $\vec L$ 守恒。
应用：开普勒第二定律、花样滑冰。

**知识点2：角动量守恒与开普勒第二定律**
推导：有心力 $\vec F\parallel\vec r$，故 $\vec M=0$，$\vec L=\text{const}$。
面积速度：$\mathrm dA=\tfrac12|\vec r\times\mathrm d\vec r|=\tfrac12|\vec r\times\vec v|\mathrm dt=\dfrac{L}{2m}\mathrm dt$，
$$\frac{\mathrm dA}{\mathrm dt}=\frac{L}{2m}=\text{const}.$$
意义：行星近日点快、远日点慢。
特例：圆轨道面积为周期相关。
应用：行星运动、卫星。

**知识点3：冲量-角动量定理**
推导：$\int\vec M\,\mathrm dt=\vec L_2-\vec L_1$（冲量矩等于角动量增量）。
意义：力矩的时间累积。
应用：碰撞中的转动。

**知识点4：质点系角动量定理**
推导：内力矩成对相消，故 $\dfrac{\mathrm d\vec L}{\mathrm dt}=\vec M^{\text{ext}}$。
特例：无外力矩时系统角动量守恒。
应用：旋转平台、陀螺。

**知识点5：碰撞与恢复系数**
定义：恢复系数 $e=\dfrac{v_2-v_1}{v_{10}-v_{20}}$（分离速度／接近速度）。
推导：一维两球碰撞，动量守恒 $m_1v_{10}+m_2v_{20}=m_1v_1+m_2v_2$，与恢复系数定义联立解得
$$v_1=\frac{(m_1-em_2)v_{10}+(1+e)m_2v_{20}}{m_1+m_2},\qquad v_2=\frac{(1+e)m_1v_{10}+(m_2-em_1)v_{20}}{m_1+m_2}.$$
动能损失 $\Delta E_k=\dfrac{m_1m_2}{2(m_1+m_2)}(1-e^2)(v_{10}-v_{20})^2$。
意义：$e=1$ 完全弹性，$e=0$ 完全非弹性。特例：等质量完全弹性碰撞交换速度。应用：台球、碰撞安全。

**知识点6：守恒律的层次与适用条件**
内容：动量守恒（合外力为零）、角动量守恒（合外力矩为零）、机械能守恒（非保守力不做功）。推导：三者分别对应时间、空间平移与空间旋转对称性（诺特定理），即系统拉格朗日量在这些变换下不变。意义：对称性决定守恒律。特例：碰撞中动量恒守恒，但机械能仅在 $e=1$ 时守恒。应用：多过程问题联立求解。

### 第5章 刚体动力学

**知识点1：转动惯量与平行轴定理**
定义：$I=\sum m_ir_i^2=\int r^2\,\mathrm dm$。
平行轴定理推导：设质心 $C$，轴过 $C$ 与过 $O$ 平行，$O$ 相对 $C$ 位矢 $\vec d$，则
$$I_O=\int|\vec r|^2\mathrm dm=\int|\vec r_C+\vec d|^2\mathrm dm=I_C+2\vec d\cdot\int\vec r_C\,\mathrm dm+Md^2.$$
因 $\int\vec r_C\,\mathrm dm=0$，故 $I_O=I_C+Md^2$。
意义：轴平移后惯量增 $Md^2$。
特例：细杆绕端点 $I=\frac13ml^2$（由 $I_C=\frac1{12}ml^2$）。
应用：复摆、飞轮。

**知识点2：转动定律**
推导：刚体视为质点系，$\mathrm d\vec L=\vec r\times\mathrm d\vec p$，绕定轴 $z$：
$$M_z=\int r_\perp\,\mathrm dF_\perp.$$
由 $\vec M=\dfrac{\mathrm d\vec L}{\mathrm dt}$，定轴时 $\vec L=I\vec\omega$，$I$ 不变，故
$$M=I\beta=I\frac{\mathrm d\omega}{\mathrm dt}.$$
意义：转动定律对应 $F=ma$。
特例：$M=0$ 匀速转动。
应用：滑轮组、绞车。

**知识点3：刚体动能与功**
推导：$E_k=\sum\tfrac12m_iv_i^2=\tfrac12\sum m_i(\omega r_i)^2=\tfrac12I\omega^2$。
力矩功：$\mathrm dW=M\,\mathrm d\theta$，$W=\int M\,\mathrm d\theta$。
功能定理：$W=\tfrac12I\omega_2^2-\tfrac12I\omega_1^2$。
意义：转动动能与平动动能对应。
特例：滚动物体动能 $=\tfrac12mv_C^2+\tfrac12I_C\omega^2$。
应用：滚下斜面速率 $v=\sqrt{\dfrac{2gh}{1+I_C/(mR^2)}}$。

**知识点4：复摆**
推导：刚体绕水平轴微摆，$\theta$ 小，重力矩 $M=-mgd\sin\theta\approx-mgd\,\theta$，$d$ 为质心到轴距。
$$I\ddot\theta=-mgd\,\theta\Rightarrow \ddot\theta+\frac{mgd}{I}\theta=0,\quad T=2\pi\sqrt{\frac{I}{mgd}}.$$
意义：由周期测质心与惯量。
特例：单摆 $I=md^2$，$T=2\pi\sqrt{d/g}$。
应用：物理摆、可倒摆。

**知识点5：刚体角动量守恒**
推导：$\vec M^{\text{ext}}=0\Rightarrow I\omega=\text{const}$。
变惯量：$I_1\omega_1=I_2\omega_2$。
意义：伸臂减速、收臂加速。
特例：质心过轴、人偶系统。
应用：跳水、体操。

**知识点5：刚体的滚动与瞬心**
推导：纯滚动条件 $v_C=R\omega$，动能
$$E_k=\tfrac12mv_C^2+\tfrac12I_C\omega^2.$$
由机械能守恒滚下高度 $h$：
$$v_C=\sqrt{\frac{2gh}{1+I_C/(mR^2)}}.$$
瞬心法：接触点瞬时速度为零，$v_C=\omega R$，得力矩与动量关系。
意义：纯滚动摩擦力不做功。特例：实心球 $I_C=\tfrac25mR^2$，$v_C=\sqrt{10gh/7}$。应用：轮、圆柱、球的下滚。

### 第6章 万有引力与天体运动

**知识点1：万有引力定律**
定义：$\vec F=-G\dfrac{m_1m_2}{r^2}\hat r$。
推导：由开普勒第三定律 $T^2\propto a^3$ 与圆周运动 $F=m\dfrac{4\pi^2}{T^2}r$，得 $F\propto\dfrac{m}{r^2}$，再由第三定律对称性引入 $M$。
意义：引力与质量乘积成正比、与距离平方成反比。
特例：$r\to0$ 时牛顿定律失效（需广义相对论）。
应用：卫星、潮汐。

**知识点2：引力势能与宇宙速度**
推导：$E_p(r)=-\int_\infty^r\left(-\dfrac{GMm}{r'^2}\right)\mathrm dr'=-\dfrac{GMm}{r}$。
逃逸速度：$\tfrac12mv^2-\dfrac{GMm}{R}=0\Rightarrow v=\sqrt{2GM/R}$。
第一宇宙速度（近地圆轨道）：$mg=\dfrac{mv^2}{R}\Rightarrow v=\sqrt{gR}\approx7.9\ \text{km/s}$。
意义：能量判据。
特例：$v<v_e$ 为椭圆/圆轨道。
应用：航天器发射。

**知识点3：开普勒三定律推导**
第一：轨道为椭圆，太阳在一焦点。
第二：面积速度守恒（见第4章）。
第三：$T^2=\dfrac{4\pi^2}{GM}a^3$。
推导：圆轨道 $GM/r^2=v^2/r$ 与 $T=2\pi r/v$ 联立即得 $T^2=\dfrac{4\pi^2}{GM}r^3$；椭圆由 $r=a$ 推广。
意义：统一描述行星运动。
特例：$a=r$ 圆轨道。
应用：由周期测质量。

**知识点4：潮汐力**
推导：距中心 $r$ 处 E 点引力加速度 $g_E=-GM/r^2$，地心 C 的 $g_C=-GM/R^2$，潮汐加速度
$$a_t=g_E-g_C.$$
设 E 沿径向偏离 $R$ 距离 $x$，$r=R+x$，展开
$$a_t=\frac{GM}{R^2}-\frac{GM}{(R+x)^2}\approx\frac{2GM}{R^3}x.$$
意义：潮汐力正比 $1/R^3$、正比偏离量。特例：近地侧与远地侧都相对地心向外隆起。应用：潮汐、洛希极限。

**知识点5：二体问题与约化质量**
推导：设相对位矢 $\vec r=\vec r_1-\vec r_2$，两体方程相减
$$\mu\ddot{\vec r}=-\frac{GMm}{r^2}\hat r,\qquad \mu=\frac{m_1m_2}{m_1+m_2}.$$
意义：二体问题约化为单体问题，力心取在质心。特例：$m_2\gg m_1$ 时 $\mu\to m_1$。应用：双星、行星-卫星。

**知识点6：轨道能量与圆锥曲线**
推导：能量 $E=\tfrac12\mu v^2-\dfrac{GM\mu}{r}$，角动量 $L$，偏心率
$$e=\sqrt{1+\frac{2EL^2}{\mu(GM)^2}}.$$
$E<0$ 椭圆，$E=0$ 抛物线，$E>0$ 双曲线。
意义：能量与角动量决定轨道类型。特例：圆轨道 $e=0$，$E=-\dfrac{GM\mu}{2a}$。应用：卫星轨道设计。

### 第7章 振动

**知识点1：简谐运动方程与能量**
定义：回复力 $F=-kx$，$m\ddot x=-kx$。
推导：令 $\omega^2=k/m$，方程 $\ddot x+\omega^2x=0$；试探解 $x=A\cos(\omega t+\varphi)$：
$$\dot x=-A\omega\sin(\omega t+\varphi),\quad \ddot x=-A\omega^2\cos=-\omega^2x.$$
能量：$E_k=\tfrac12m\dot x^2=\tfrac12mA^2\omega^2\sin^2$，$E_p=\tfrac12kx^2=\tfrac12kA^2\cos^2$，
$$E=E_k+E_p=\tfrac12kA^2=\text{const}.$$
意义：简谐运动机械能守恒。
特例：相位 $\pm\omega t$ 决定初状态。
应用：弹簧振子、LC 回路。

**知识点2：单摆与复摆**
单摆：切向 $ml\ddot\theta=-mg\sin\theta$，小角 $\sin\theta\approx\theta$，$\ddot\theta+\dfrac gl\theta=0$，$T=2\pi\sqrt{l/g}$。
复摆见第5章。
意义：小振幅等时性。
特例：大角度非线性，周期随振幅增大。
应用：计时。

**知识点3：阻尼振动**
假设：阻力 $F_f=-\gamma v$。
推导：$m\ddot x=-kx-\gamma\dot x\Rightarrow \ddot x+2\beta\dot x+\omega_0^2x=0$，$\beta=\dfrac{\gamma}{2m}$。
特征方程 $\lambda^2+2\beta\lambda+\omega_0^2=0$，$\lambda=-\beta\pm\sqrt{\beta^2-\omega_0^2}$。
欠阻尼 $\beta<\omega_0$：$\omega=\sqrt{\omega_0^2-\beta^2}$，$x=Ae^{-\beta t}\cos(\omega t+\varphi)$。
过阻尼、临界阻尼分别讨论。
意义：振幅指数衰减。
特例：$\beta\to0$ 退化为简谐。
应用：减震器、仪表阻尼。

**知识点4：受迫振动与共振**
推导：$m\ddot x=-kx-\gamma\dot x+F_0\cos\omega t$，写成
$$\ddot x+2\beta\dot x+\omega_0^2x=\frac{F_0}{m}\cos\omega t.$$
稳态解 $x=A\cos(\omega t-\delta)$，代入得
$$A=\frac{F_0/m}{\sqrt{(\omega_0^2-\omega^2)^2+4\beta^2\omega^2}},\quad \tan\delta=\frac{2\beta\omega}{\omega_0^2-\omega^2}.$$
共振：$\omega_r=\sqrt{\omega_0^2-2\beta^2}$ 为振幅极大点。
意义：驱动频率接近固有频率时振幅剧增。
特例：$\beta\to0$，$\omega\to\omega_0$ 共振发散。
应用：共振破坏、隔振。

**知识点5：同方向简谐振动合成**
推导：$x_1=A_1\cos(\omega t+\varphi_1)$，$x_2=A_2\cos(\omega t+\varphi_2)$，由旋转矢量相加，
$$A=\sqrt{A_1^2+A_2^2+2A_1A_2\cos\Delta\varphi},\quad \tan\varphi=\frac{A_1\sin\varphi_1+A_2\sin\varphi_2}{A_1\cos\varphi_1+A_2\cos\varphi_2}.$$
意义：相位差决定增强或减弱。
特例：$\Delta\varphi=0$ 加强，$\pi$ 减弱。
应用：拍、干涉。

**知识点6：非线性振动与相图**
定义：$\ddot x+f(x)=0$，能量 $E=\tfrac12\dot x^2+V(x)$ 守恒。
推导：由能量式解出 $\dot x=\pm\sqrt{2(E-V(x))}$，相轨迹为闭曲线时对应周期运动；单摆能量 $\tfrac12\dot\theta^2-\tfrac gl\cos\theta=E$，小振动为椭圆，大振幅为闭曲线，分界线（separatrix）满足 $E=g/l$。
意义：相图几何可判稳定性。特例：小角近似为简谐。应用：单摆、等离子体振荡。

### 第8章 波动

**知识点1：波动方程与波速**
定义：一维绳波满足 $\dfrac{\partial^2 y}{\partial x^2}=\dfrac1{v^2}\dfrac{\partial^2 y}{\partial t^2}$。
推导（微元）：绳张力 $T$，线密度 $\mu$，取微元 $\mathrm dx$，竖直合力
$$\Delta F_y=T\left(\frac{\partial y}{\partial x}\Big|_{x+\mathrm dx}-\frac{\partial y}{\partial x}\Big|_x\right)=T\frac{\partial^2y}{\partial x^2}\mathrm dx.$$
由牛顿第二定律 $T\dfrac{\partial^2y}{\partial x^2}\mathrm dx=\mu\,\mathrm dx\dfrac{\partial^2y}{\partial t^2}$，得 $v=\sqrt{T/\mu}$。
意义：一切满足该方程的扰动以 $v$ 传播。
特例：行波解 $y=f(x\mp vt)$。
应用：声波、地震波。

**知识点2：平面简谐波**
定义：$y=A\cos[\omega(t-x/v)+\varphi_0]=A\cos(kx-\omega t+\varphi_0)$，$k=2\pi/\lambda$。
推导：将 $t'=t-x/v$ 代入简谐振动式。
波数关系：$\lambda=vT=2\pi v/\omega$。
意义：空间周期性 $\lambda$，时间周期性 $T$。
特例：沿 $+x$ 传播；$x$ 增大相位滞后。
应用：声波、光波。

**知识点3：波的能量与能流密度**
推导：微元动能与势能
$$\mathrm dE_k=\tfrac12(\mu\,\mathrm dx)\left(\frac{\partial y}{\partial t}\right)^2,\quad \mathrm dE_p=\tfrac12T\left(\frac{\partial y}{\partial x}\right)^2\mathrm dx.$$
由行波关系 $\partial y/\partial t=-v\,\partial y/\partial x$，得 $\mathrm dE_k=\mathrm dE_p$。
平均能流密度（强度）$I=\tfrac12\rho v\omega^2A^2$。
意义：能流与振幅平方、频率平方成正比。
特例：球面波 $A\propto1/r$，$I\propto1/r^2$。
应用：声强、能量衰减。

**知识点4：惠更斯原理与波的叠加**
内容：波前上每一点都是子波源，子波包络为新的波前。
推导：两同频波叠加 $y=y_1+y_2$，合振幅 $A=\sqrt{A_1^2+A_2^2+2A_1A_2\cos\Delta\varphi}$，$\Delta\varphi=\dfrac{2\pi}{\lambda}\Delta r$。
意义：说明衍射、反射、折射。
特例：相干条件频率相同、相位差恒定。
应用：干涉、衍射。

**知识点5：驻波**
推导：两反向行波
$$y_1=A\cos(kx-\omega t),\quad y_2=A\cos(kx+\omega t),$$
相加用和差化积：$y=2A\cos kx\cos\omega t$。
节点 $\cos kx=0\Rightarrow x=(2n+1)\lambda/4$；腹点 $x=n\lambda/2$。
意义：波形不移动、能量不传播。
特例：弦两端固定 $L=n\lambda/2$，基频 $f_1=\dfrac{v}{2L}$。
应用：乐器、激光腔。

**知识点6：多普勒效应**
推导：设波源、观察者沿连线速度 $v_s,v_o$，波速 $u$。
源动：波长变 $\lambda'=(u-v_s)T=u/f-v_s/f$；观察者接收频率
$$f'=\frac{u+v_o}{\lambda'}=\frac{u+v_o}{u-v_s}f.$$
意义：相对运动改变接收频率。
特例：$v_o=v_s=0$ 时 $f'=f$。
应用：雷达测速、天体红移。

### 第9章 狭义相对论力学基础

**知识点1：洛伦兹变换**
假设：①相对性原理 ②光速不变原理。
推导：设 $S'$ 相对 $S$ 以 $u$ 沿 $x$ 运动，线性变换 $x'=\gamma(x-ut)$，$x=\gamma(x'+ut')$。
代入 $x=ct$ 时 $x'=ct'$：
$$ct'=\gamma(ct-ut),\quad ct=\gamma(ct'+ut').$$
两式相乘消去，得 $\gamma^2(1-u^2/c^2)=1\Rightarrow \gamma=\dfrac1{\sqrt{1-u^2/c^2}}$。
完整：
$$x'=\gamma(x-ut),\ y'=y,\ z'=z,\ t'=\gamma\left(t-\frac{ux}{c^2}\right).$$
意义：时间与空间耦合。
特例：$u\ll c$，$\gamma\to1$，退化为伽利略变换。
应用：高能粒子。

**知识点2：时间膨胀与长度收缩**
推导：同地两事件 $\Delta x'=0$，$\Delta t'=\Delta\tau$，则 $\Delta t=\gamma\Delta\tau$（固有时最短）。
测长：同时测两端 $\Delta t=0$，$l=|x_2'-x_1'|=\gamma^{-1}l_0$，即动尺缩短。
意义：运动时钟变慢、运动尺子缩短。
特例：$u\to c$，$\gamma\to\infty$。
应用：μ 子寿命、GPS 修正。

**知识点3：相对论质量、动量与能量**
定义：$\vec p=\gamma m_0\vec v$，总能量 $E=\gamma m_0c^2$，静能 $E_0=m_0c^2$。
推导：由 $E^2=p^2c^2+m_0^2c^4$，及恒等式 $\gamma^2-\gamma^2v^2/c^2=1$，
$$\Rightarrow E^2-c^2p^2=m_0^2c^4.$$
动能 $E_k=E-E_0=(\gamma-1)m_0c^2$。
意义：质能等价。
特例：$v\ll c$，$E_k\approx\tfrac12m_0v^2$。
应用：核能、粒子物理。

**知识点4：同时性的相对性与因果律**
推导：由洛伦兹变换
$$\Delta t'=\gamma\left(\Delta t-\frac{u\Delta x}{c^2}\right).$$
若在 $S$ 中同时（$\Delta t=0$），则 $\Delta t'=-\gamma u\Delta x/c^2\neq0$，即同时性相对。
意义：不同惯性系对“同时”判断不同。特例：同地事件（$\Delta x=0$）同时性绝对。应用：类时间隔的因果顺序不变。

**知识点5：相对论速度变换**
推导：由 $v_x'=\dfrac{\mathrm dx'}{\mathrm dt'}$，代入洛伦兹变换微分式：
$$v_x'=\frac{v_x-u}{1-uv_x/c^2},\quad v_y'=\frac{v_y}{\gamma(1-uv_x/c^2)},\quad v_z'=\frac{v_z}{\gamma(1-uv_x/c^2)}.$$
意义：速度合成不超过 $c$。特例：若 $v_x=c$，则 $v_x'=c$，与光速不变一致。应用：相对论喷流。

**知识点6：相对论动力学方程与质能应用**
推导：由 $\vec F=\dfrac{\mathrm d}{\mathrm dt}(\gamma m_0\vec v)$；质能关系 $E=mc^2$，$\Delta E=\Delta m\,c^2$。
意义：质量与能量等价。特例：光子 $m_0=0$、$E=pc$。应用：核裂变与聚变的能量释放。

### 第10章 分析力学：拉格朗日方程

**知识点1：约束与广义坐标**
定义：完整约束减少自由度数；用广义坐标 $q_1,\dots,q_s$ 描述，$s=3N-k$。
推导：位矢 $\vec r_i=\vec r_i(q_1,\dots,q_s,t)$，故
$$\dot{\vec r}_i=\sum_j\frac{\partial\vec r_i}{\partial q_j}\dot q_j+\frac{\partial\vec r_i}{\partial t}.$$
意义：自动消去约束力。
特例：无约束时取直角坐标。
应用：多体系统降维。

**知识点2：虚功原理与达朗贝尔原理**
定义：理想约束下 $\sum_i(\vec F_i-m_i\vec a_i)\cdot\delta\vec r_i=0$。
推导：平衡时 $\sum_i\vec F_i\cdot\delta\vec r_i=0$（约束力虚功为零）；动力学加惯性力 $-\vec a_i$ 即可。
意义：虚位移原理统一静力学与动力学。
特例：静止时退化为虚功原理。
应用：结构静力。

**知识点3：拉格朗日方程推导**
推导：由达朗贝尔原理，
$$\sum_i(\vec F_i-m_i\ddot{\vec r}_i)\cdot\delta\vec r_i=0,$$
代入 $\delta\vec r_i=\sum_j\dfrac{\partial\vec r_i}{\partial q_j}\delta q_j$，得
$$\sum_j\left[\sum_i\vec F_i\cdot\frac{\partial\vec r_i}{\partial q_j}-\sum_i m_i\ddot{\vec r}_i\cdot\frac{\partial\vec r_i}{\partial q_j}\right]\delta q_j=0.$$
利用 $\dfrac{\partial\vec r_i}{\partial q_j}=\dfrac{\partial\dot{\vec r}_i}{\partial \dot q_j}$ 与分部求导
$$\sum_i m_i\ddot{\vec r}_i\cdot\frac{\partial\vec r_i}{\partial q_j}=\frac{\mathrm d}{\mathrm dt}\frac{\partial T}{\partial\dot q_j}-\frac{\partial T}{\partial q_j},$$
定义 $Q_j=\sum_i\vec F_i\cdot\partial\vec r_i/\partial q_j$，若 $\vec F_i=-\nabla_iV$ 则 $Q_j=-\partial V/\partial q_j$。
由于 $\delta q_j$ 任意：
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}=0,\qquad L=T-V.$$
意义：牛顿力学的推广，便于处理约束。
特例：$L$ 不显含 $q_j$ 时 $\partial L/\partial\dot q_j$ 守恒（循环坐标）。
应用：双摆、中心力。

**知识点4：循环坐标与守恒量**
推导：若 $\dfrac{\partial L}{\partial q_j}=0$，则 $p_j=\dfrac{\partial L}{\partial\dot q_j}$ 守恒。
意义：对称性对应守恒量（诺特）。
特例：中心力 $L$ 与 $\theta$ 无关，$r^2\dot\theta$ 守恒。
应用：行星运动。

**知识点5：拉格朗日方程应用示例**
中心力：$L=\tfrac12m(\dot r^2+r^2\dot\theta^2)-V(r)$。
对 $r$：$m\ddot r-mr\dot\theta^2=-V'(r)$；对 $\theta$：$mr^2\dot\theta=\text{const}$。
意义：可导出开普勒问题。
特例：$V=-GMm/r$。
应用：轨道积分。

**知识点6：拉格朗日乘子与约束力**
推导：完整约束 $f_\alpha(q,t)=0$，增广泛函
$$S^*=\int\left(L+\sum_\alpha\lambda_\alpha f_\alpha\right)\mathrm dt,$$
变分得
$$\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}=\sum_\alpha\lambda_\alpha\frac{\partial f_\alpha}{\partial q_j}.$$
意义：$\lambda_\alpha$ 即约束力分量。特例：刚性杆的轴力。应用：多刚体系统。

### 第11章 分析力学：哈密顿力学

**知识点1：勒让德变换与哈密顿量**
定义：$H=\sum_jp_j\dot q_j-L$，$p_j=\dfrac{\partial L}{\partial\dot q_j}$。
推导：由 $L$ 的微元
$$\mathrm dL=\sum_j\frac{\partial L}{\partial q_j}\mathrm dq_j+\sum_j\frac{\partial L}{\partial\dot q_j}\mathrm d\dot q_j=\sum_j\dot p_j\,\mathrm dq_j+\sum_jp_j\,\mathrm d\dot q_j.$$
而 $H=\sum_jp_j\dot q_j-L$，$\mathrm dH=\sum_j(\dot q_j\,\mathrm dp_j-\dot p_j\,\mathrm dq_j)$，故
$$\dot q_j=\frac{\partial H}{\partial p_j},\qquad \dot p_j=-\frac{\partial H}{\partial q_j}.$$
意义：一阶正则方程组，相空间流。
特例：$L=\tfrac12m\dot x^2-V$，$H=\dfrac{p^2}{2m}+V=E$。
应用：统计力学、量子化。

**知识点2：泊松括号**
定义：$\{A,B\}=\sum_j\left(\dfrac{\partial A}{\partial q_j}\dfrac{\partial B}{\partial p_j}-\dfrac{\partial A}{\partial p_j}\dfrac{\partial B}{\partial q_j}\right)$。
推导：$\dfrac{\mathrm dA}{\mathrm dt}=\{A,H\}+\dfrac{\partial A}{\partial t}$。
意义：$A$ 守恒 $\Leftrightarrow\{A,H\}=0$。
特例：$\{q_i,p_j\}=\delta_{ij}$。
应用：对称性与守恒量。

**知识点3：正则变换**
定义：变换 $(q,p)\to(Q,P)$ 保持正则方程形式。
推导：要求 $\sum p\,\mathrm dq-H\,\mathrm dt=\sum P\,\mathrm dQ-K\,\mathrm dt+\mathrm dF$，$F$ 为生成函数，得变换关系。
意义：选取可积坐标（作用量-角变量）。
特例：恒等变换。
应用：天体力学。

**知识点4：哈密顿量守恒与能量关系**
推导：由正则方程
$$\frac{\mathrm dH}{\mathrm dt}=\frac{\partial H}{\partial t}=-\frac{\partial L}{\partial t}.$$
若 $L$ 不显含 $t$，则 $H$ 守恒。
意义：时间平移不变对应能量守恒。特例：当动能 $T$ 为 $\dot q$ 的二次齐次函数时 $H=T+V=E$。应用：判别耗散系统。

**知识点5：哈密顿-雅可比方程**
推导：取生成函数 $S(q,t)$，$p=\partial S/\partial q$，代入 $H(q,p,t)+\partial S/\partial t=0$ 得
$$H\left(q,\frac{\partial S}{\partial q},t\right)+\frac{\partial S}{\partial t}=0.$$
意义：将运动化为求解一阶偏微分方程。特例：可分离变量时得作用量积分。应用：天体力学、WKB 近似。

**知识点6：相空间与刘维尔定理**
推导：相流速度 $(\dot q,\dot p)=(\partial H/\partial p,-\partial H/\partial q)$，其散度
$$\nabla\cdot\vec V=\frac{\partial\dot q}{\partial q}+\frac{\partial\dot p}{\partial p}=\frac{\partial^2H}{\partial q\partial p}-\frac{\partial^2H}{\partial p\partial q}=0.$$
由连续性方程 $\partial_t\rho+\nabla\cdot(\rho\vec V)=0$ 得 $\mathrm d\rho/\mathrm dt=0$。
意义：相空间体积守恒。特例：无耗散保守系统。应用：统计力学基础。

### 第12章 变分原理与对称性

**知识点1：最小作用量原理**
定义：$\delta S=\delta\displaystyle\int_{t_1}^{t_2}L(q,\dot q,t)\,\mathrm dt=0$。
推导：设 $q(t)+\delta q(t)$，$\delta q(t_1)=\delta q(t_2)=0$，
$$\delta S=\int\left(\frac{\partial L}{\partial q}\delta q+\frac{\partial L}{\partial\dot q}\delta\dot q\right)\mathrm dt.$$
分部积分第二项：$\int\dfrac{\partial L}{\partial\dot q}\dfrac{\mathrm d}{\mathrm dt}\delta q\,\mathrm dt=\left[\dfrac{\partial L}{\partial\dot q}\delta q\right]-\int\dfrac{\mathrm d}{\mathrm dt}\dfrac{\partial L}{\partial\dot q}\delta q\,\mathrm dt$。
边界项为零，故
$$\delta S=\int\left(\frac{\partial L}{\partial q}-\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q}\right)\delta q\,\mathrm dt=0\Rightarrow \frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q}-\frac{\partial L}{\partial q}=0.$$
意义：真实路径使作用量取极值。
特例：最短路径、费马原理同构。
应用：场论、最优控制。

**知识点2：诺特定理**
内容：每个连续对称性对应一个守恒量。
推导：设变换 $q_i\to q_i+\epsilon\,\eta_i$ 使 $\delta L=\epsilon\dfrac{\mathrm dF}{\mathrm dt}$，则
$$J=\sum_i\frac{\partial L}{\partial\dot q_i}\eta_i-F=\text{const}.$$
意义：时间平移→能量，空间平移→动量，转动→角动量。
特例：循环坐标 $F=0$，$p_j$ 守恒。
应用：守恒律统一解释。

**知识点3：最速降线与变分法应用**
定义：最短时间泛函 $T=\displaystyle\int\sqrt{\frac{1+y'^2}{2gy}}\,\mathrm dx$。
推导：被积函数不显含 $x$，由欧拉-拉格朗日方程的首次积分 $F-y'\partial F/\partial y'=C$，解得摆线
$$x=a(\theta-\sin\theta),\qquad y=a(1-\cos\theta).$$
意义：变分法经典范例。特例：直线并非最短时间路径。应用：最优控制、光学。

**知识点4：哈密顿原理与场的推广**
定义：$\delta\displaystyle\int_{t_1}^{t_2}L\,\mathrm dt=0$（位形空间路径），与最小作用量原理等价。
推导：对被积函数 $L(q,\dot q,t)$ 变分并分部积分，边界 $\delta q=0$，即得拉格朗日方程。连续场推广：拉格朗日密度 $\mathcal L(\phi,\partial_\mu\phi)$，
$$\partial_\mu\frac{\partial\mathcal L}{\partial(\partial_\mu\phi)}-\frac{\partial\mathcal L}{\partial\phi}=0.$$
意义：从有限自由度到场的统一描述。特例：平衡态最小势能。应用：弹性动力学、电磁场。

**知识点5：诺特定理与三大守恒律**
推导：分别取时间平移 $t\to t+\epsilon$、空间平移 $\vec r\to\vec r+\epsilon\hat n$、转动 $\delta\vec r=\delta\vec\theta\times\vec r$，代入
$$J=\sum_i\frac{\partial L}{\partial\dot q_i}\eta_i-F=\text{const}$$
得能量、动量、角动量守恒。
意义：对称性对应守恒流。特例：循环坐标给出动量守恒。应用：粒子物理中的守恒荷。

**知识点6：场论初步的守恒荷**
定义：若四维流满足 $\partial_\mu j^\mu=0$，则荷
$$Q=\int j^0\,\mathrm d^3x$$
守恒。推导：对时空体积积分并利用高斯定理，边界流为零。意义：连续对称性产生守恒荷。特例：电磁四维流。应用：量子场论。

### 第13章 连续介质力学基础

**知识点1：应力张量**
定义：面元 $\mathrm d\vec A=\hat n\,\mathrm dA$ 上的力 $\mathrm d\vec F=\boldsymbol\sigma\cdot\mathrm d\vec A$。
推导：由微元四面体平衡（柯西定理），
$$T_i^{(n)}=\sigma_{ij}n_j.$$
意义：应力张量对称 $\sigma_{ij}=\sigma_{ji}$（无体力矩）。
特例：静水压力 $\sigma_{ij}=-p\delta_{ij}$。
应用：固体、流体。

**知识点2：应变张量与几何方程**
定义：位移 $\vec u$，小应变张量
$$\varepsilon_{ij}=\frac12\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\right).$$
推导：两邻近点位移差 $\mathrm du_i=\partial_j u_i\,\mathrm dx_j$，分解为对称（应变）与反对称（转动）部分。
意义：$\varepsilon_{ii}$ 为体应变。
特例：纯刚体运动时 $\varepsilon_{ij}=0$。
应用：变形测量。

**知识点3：广义胡克定律**
假设：线弹性、小变形。
公式：$\sigma_{ij}=C_{ijkl}\varepsilon_{kl}$，各向同性时
$$\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij}.$$
推导：由各向同性张量结构 $C_{ijkl}=\lambda\delta_{ij}\delta_{kl}+\mu(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk})$ 代入即得。
逆式：$\varepsilon_{ij}=\dfrac{1+\nu}{E}\sigma_{ij}-\dfrac{\nu}{E}\sigma_{kk}\delta_{ij}$。
意义：应力应变线性关系。
特例：单轴拉伸得 $E,\nu$ 与 $\lambda,\mu$ 关系 $\lambda=\dfrac{E\nu}{(1+\nu)(1-2\nu)},\ \mu=\dfrac{E}{2(1+\nu)}$。
应用：结构强度。

**知识点4：平衡微分方程**
推导：微元受力平衡，面力与体力
$$\int_{\partial V}\boldsymbol\sigma\cdot\hat n\,\mathrm dA+\int_V\vec f\,\mathrm dV=0.$$
用高斯定理化为体积分 $\int_V(\nabla\cdot\boldsymbol\sigma+\vec f)\mathrm dV=0$，得
$$\frac{\partial\sigma_{ij}}{\partial x_j}+f_i=0.$$
意义：应力散度平衡体力。
特例：无体力 $\nabla\cdot\boldsymbol\sigma=0$。
应用：梁、板的应力分布。

**知识点5：应变协调方程与平面问题**
推导：由几何方程 $\varepsilon_{ij}=\tfrac12(\partial_j u_i+\partial_i u_j)$ 消去位移 $u_i$，得圣维南协调方程
$$\partial_k\partial_l\varepsilon_{ij}+\partial_i\partial_j\varepsilon_{kl}-\partial_i\partial_l\varepsilon_{jk}-\partial_j\partial_k\varepsilon_{il}=0.$$
意义：保证应变场可由单值位移积分得到。特例：平面应力 $\sigma_z=0$，平面应变 $\varepsilon_z=0$。应用：板、平面弹性问题。

**知识点6：应变能与卡氏定理**
推导：线弹性应变能密度 $u=\tfrac12\sigma_{ij}\varepsilon_{ij}$，总应变能 $U=\int u\,\mathrm dV$。卡氏第二定理
$$\delta_i=\frac{\partial U}{\partial F_i}.$$
意义：能量法求位移与超静定。特例：轴力杆 $U=\dfrac{F^2l}{2EA}$。应用：桁架、框架。

### 第14章 流体力学与 N-S 方程

**知识点1：连续性方程**
假设：质量守恒。
推导：控制体质量变化率等于净流入：
$$\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\vec v)=0.$$
一维定常不可压：$A_1v_1=A_2v_2$。
意义：质量守恒的场形式。
特例：不可压 $\nabla\cdot\vec v=0$。
应用：管道流量。

**知识点2：欧拉方程**
推导：对理想流体（无粘）微元用牛顿第二定律，
$$\rho\frac{\mathrm D\vec v}{\mathrm Dt}=-\nabla p+\vec f,\qquad \frac{\mathrm D}{\mathrm Dt}=\frac{\partial}{\partial t}+\vec v\cdot\nabla.$$
意义：理想流体动量方程。
特例：静流体 $\nabla p=\rho\vec g$。
应用：气动。

**知识点3：纳维-斯托克斯方程推导**
假设：牛顿粘性流体，应力 $\sigma_{ij}=-p\delta_{ij}+\tau_{ij}$，$\tau_{ij}=2\mu\varepsilon_{ij}-\tfrac23\mu\varepsilon_{kk}\delta_{ij}$。
推导：代入 Cauchy 方程 $\rho\dfrac{\mathrm Dv_i}{\mathrm Dt}=\partial_j\sigma_{ij}+f_i$，
$$\partial_j\tau_{ij}=\mu\partial_j(\partial_jv_i+\partial_iv_j)-\tfrac23\mu\partial_i(\partial_jv_j)=\mu\nabla^2v_i+\tfrac13\mu\partial_i(\nabla\cdot\vec v).$$
不可压 $\nabla\cdot\vec v=0$ 时得
$$\rho\left(\frac{\partial\vec v}{\partial t}+\vec v\cdot\nabla\vec v\right)=-\nabla p+\mu\nabla^2\vec v+\vec f.$$
意义：粘性流体运动基本方程。
特例：$\mu=0$ 退化为欧拉方程；$\vec v=0$ 退化为流体静力学。
应用：湍流、边界层。

**知识点4：伯努利方程**
推导：沿流线对定常、不可压、无粘欧拉方程积分：
$$\vec v\cdot\nabla\left(\tfrac12v^2\right)=-\frac1\rho\vec v\cdot\nabla p.$$
沿流线 $\mathrm d\vec r\parallel\vec v$，得
$$\mathrm d\left(\tfrac12v^2\right)+\frac{\mathrm dp}{\rho}=0\Rightarrow \tfrac12\rho v^2+p+\rho gz=\text{const}.$$
意义：流速大处压强小。
特例：静止流体退化为 $p+\rho gz=\text{const}$。
应用：机翼升力、文丘里管。

**知识点5：雷诺数与粘性流态**
定义：$Re=\dfrac{\rho vL}{\mu}$。
推导：惯性项 $\rho v^2/L$ 与粘性项 $\mu v/L^2$ 之比。
意义：$Re$ 大惯性主导易湍流，小粘性主导层流。
特例：圆管临界 $Re\approx2300$。
应用：管道、飞行器。

**知识点6：泊肃叶流与库埃特流**
推导：圆管定常层流 N-S 方程化为
$$0=-\frac{\mathrm dp}{\mathrm dx}+\frac{\mu}{r}\frac{\mathrm d}{\mathrm dr}\left(r\frac{\mathrm dv}{\mathrm dr}\right),$$
积分并用无滑移边界 $v(R)=0$，得抛物线
$$v(r)=\frac{\Delta p}{4\mu L}(R^2-r^2),$$
流量 $Q=\dfrac{\pi R^4\Delta p}{8\mu L}$（泊肃叶定律）。库埃特流 $v(y)=Uy/h$。
意义：粘性解析解。特例：流量正比 $R^4$。应用：血流、润滑。

### 第15章 弹性杆件与梁的弯曲

**知识点1：拉压与剪切**
定义：正应力 $\sigma=F/A$，正应变 $\varepsilon=\Delta l/l$，$\sigma=E\varepsilon$。
推导：胡克定律一维形式；剪切 $\tau=G\gamma$，$G=\dfrac{E}{2(1+\nu)}$。
意义：杆件轴向刚度 $k=EA/l$。
特例：超静定由位移协调解。
应用：桁架。

**知识点2：欧拉-伯努利梁方程**
假设：平面截面保持平面；小变形；忽略剪切变形。
推导：曲率与应变 $\varepsilon_x=-\dfrac{y}{R}\approx-y\dfrac{\mathrm d^2w}{\mathrm dx^2}$，应力 $\sigma_x=E\varepsilon_x$。
弯矩 $M=\int_A(-y)\sigma_x\,\mathrm dA=E\dfrac{\mathrm d^2w}{\mathrm dx^2}\int_A y^2\mathrm dA=EI\dfrac{\mathrm d^2w}{\mathrm dx^2}$。
由微元平衡 $\dfrac{\mathrm dM}{\mathrm dx}=Q,\ \dfrac{\mathrm dQ}{\mathrm dx}=-q$，得
$$EI\frac{\mathrm d^4w}{\mathrm dx^4}=q.$$
意义：梁挠度控制方程。
特例：简支梁集中力得 $w_{\max}=\dfrac{FL^3}{48EI}$。
应用：桥梁、机翼。

**知识点3：梁的应力与强度**
推导：$\sigma_x=-Ey\dfrac{\mathrm d^2w}{\mathrm dx^2}=-\dfrac{My}{I}$，最大在边缘。
意义：中性轴应力为零。
特例：矩形截面 $I=bh^3/12$。
应用：选材、截面优化。

**知识点4：圆轴扭转**
推导：扭转变形 $\gamma=\rho\dfrac{\mathrm d\phi}{\mathrm dx}$，$\tau=G\gamma=G\rho\dfrac{\mathrm d\phi}{\mathrm dx}$，扭矩
$$T=\int\rho\tau\,\mathrm dA=GI_p\frac{\mathrm d\phi}{\mathrm dx},\quad I_p=\frac{\pi d^4}{32},$$
最大切应力 $\tau_{\max}=\dfrac{TR}{I_p}$。
意义：圆轴扭转刚度 $GI_p/L$。特例：圆截面极惯矩 $\propto d^4$。应用：传动轴。

**知识点5：瑞利-里兹能量法**
推导：设试函数 $w(x)=\sum a_i\phi_i(x)$，总势能泛函
$$\Pi=\tfrac12\int EI(w'')^2\mathrm dx-\int qw\,\mathrm dx,$$
由 $\partial\Pi/\partial a_i=0$ 求系数。
意义：变分近似法。特例：简支梁取正弦试函数得精确挠度与屈曲载荷。应用：复杂结构强度分析。

**知识点6：强度理论与组合变形**
推导：平面应力主应力
$$\sigma_{1,2}=\frac{\sigma_x+\sigma_y}{2}\pm\sqrt{\left(\frac{\sigma_x-\sigma_y}{2}\right)^2+\tau_{xy}^2}.$$
von Mises 等效应力 $\sigma_e=\sqrt{\sigma_x^2+3\tau^2}$（平面应力）。
意义：多轴应力失效判据。特例：纯剪切屈服 $\tau_Y=\sigma_Y/\sqrt3$。应用：压力容器、轴。

### 第16章 稳定性与分岔基础

**知识点1：平衡稳定性判据**
定义：受扰后能返回原平衡为稳定。
推导：势能极小 $\delta^2E_p>0$ 稳定；对动力学系统，线性化 $\delta\ddot x=-\dfrac{E_p''}{m}\delta x$。
意义：特征值实部判稳。
特例：$\omega^2>0$ 振荡稳定；$\omega^2<0$ 指数发散。
应用：结构、轨道。

**知识点2：欧拉屈曲**
推导：压杆 $EIw''''+Pw''=0$，设 $w=A\sin kx$，$k^2=P/EI$，边界 $w(0)=w(L)=0$ 要求 $kL=n\pi$。
$$P_{cr}=\frac{n^2\pi^2EI}{L^2},\quad n=1\text{ 为临界力}.$$
意义：细长杆失稳载荷。
特例：两端铰支；其他边界用有效长度 $\mu L$。
应用：柱、桁架。

**知识点3：分岔与李雅普诺夫稳定性**
定义：参数变化使平衡解数目/稳定性改变称分岔。
推导：一维 $\dot x=f(x,\mu)$，平衡 $f=0$；在平衡处线性化特征值 $\lambda=\partial f/\partial x$；$\lambda=0$ 处发生分岔（如叉式 $\dot x=\mu x-x^3$，$\mu>0$ 出现新平衡 $x=\pm\sqrt\mu$）。
意义：稳定性交换、结构突变。
特例：$\mu<0$ 唯一稳定零点。
应用：屈曲、相变、生态。

**知识点4：李雅普诺夫第二方法**
定义：若存在正定函数 $V(\vec x)$ 使沿轨线 $\dot V\le0$，则原点稳定；$\dot V<0$ 渐近稳定。
推导：沿轨线 $\dot V=\nabla V\cdot f(\vec x)$，结合 $V$ 正定即得稳定性结论。
意义：无需解方程即可判稳。特例：取机械能作 $V$。应用：控制理论。

**知识点5：后屈曲与缺陷敏感性**
推导：含初始缺陷 $\epsilon$ 时载荷-挠度关系在分枝点附近近似
$$P\approx P_{cr}\left(1-\frac{\epsilon}{w}\right),$$
缺陷使实际临界载荷低于理论值。
意义：实际结构稳定性对缺陷敏感。特例：$\epsilon\to0$ 回到分岔结果。应用：薄壳、薄板。

**知识点6：极限环与非线性稳定性**
定义：相平面中孤立闭轨对应自激振荡。
推导：范德波尔方程 $\ddot x-\mu(1-x^2)\dot x+x=0$；当 $\mu>0$ 时存在稳定极限环（能量耗散与注入平衡）。
意义：非线性系统可自持周期运动。特例：$\mu=0$ 时退化为简谐（连续闭轨，非孤立）。应用：电子振荡器、心跳节律。

---
⬆ 上级领域：[[物质科学（物理·化学·生命）]]
