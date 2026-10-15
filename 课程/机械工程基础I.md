---
课程名称: 机械工程基础I
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

# 机械工程基础I


### 第1章 质点运动学与质点动力学

#### 1.1 位矢、速度与加速度及其坐标表示
质点的运动由位置矢量 $\mathbf r(t)=x(t)\mathbf i+y(t)\mathbf j+z(t)\mathbf k$ 描述。
定义速度与加速度
$$
\mathbf v=\frac{\mathrm d\mathbf r}{\mathrm dt}=\dot{\mathbf r},\qquad \mathbf a=\frac{\mathrm d\mathbf v}{\mathrm dt}=\ddot{\mathbf r}.
$$
对定义直接求导即得
$$
\mathbf v=\dot x\,\mathbf i+\dot y\,\mathbf j+\dot z\,\mathbf k,\qquad
\mathbf a=\ddot x\,\mathbf i+\ddot y\,\mathbf j+\ddot z\,\mathbf k.
$$
物理意义：速度是位矢对时间的一阶变化率，加速度是二阶变化率；二者均为矢量，方向分别沿轨迹切向与速度变化方向。特例：$\mathbf a=\mathbf 0$ 对应匀速直线运动；$\mathbf v=\text{const}$ 亦然。应用：任何后续动力学方程都可化为对 $\mathbf r(t)$ 的常微分方程。

#### 1.2 自然坐标系的切向与法向加速度
设轨迹弧长 $s$，切向单位矢量 $\mathbf e_t=\mathrm d\mathbf r/\mathrm ds$，曲率 $\kappa=|\mathrm d\mathbf e_t/\mathrm ds|$，曲率半径 $\rho=1/\kappa$。速度
$$
\mathbf v=\frac{\mathrm ds}{\mathrm dt}\mathbf e_t=v\,\mathbf e_t.
$$
求导得
$$
\mathbf a=\frac{\mathrm d}{\mathrm dt}(v\mathbf e_t)
=\dot v\,\mathbf e_t+v\frac{\mathrm d\mathbf e_t}{\mathrm dt}
=\dot v\,\mathbf e_t+v\frac{\mathrm d\mathbf e_t}{\mathrm ds}\frac{\mathrm ds}{\mathrm dt}
=\dot v\,\mathbf e_t+v^2\kappa\,\mathbf e_n
=\dot v\,\mathbf e_t+\frac{v^2}{\rho}\mathbf e_n,
$$
其中 $\mathbf e_n$ 为主法向单位矢量。物理意义：切向项改变速率大小，法向项（向心加速度 $v^2/\rho$）改变速度方向。特例：直线运动 $\rho\to\infty$，$a_n=0$。应用：圆周运动 $a_n=\omega^2 R$。

#### 1.3 极坐标与柱坐标下的速度、加速度
平面极坐标中 $\mathbf r=r\,\mathbf e_r$，$\mathbf e_r=\cos\theta\,\mathbf i+\sin\theta\,\mathbf j$，$\mathbf e_\theta=-\sin\theta\,\mathbf i+\cos\theta\,\mathbf j$。由
$$
\dot{\mathbf e}_r=\dot\theta\,\mathbf e_\theta,\qquad \dot{\mathbf e}_\theta=-\dot\theta\,\mathbf e_r,
$$
得
$$
\mathbf v=\dot r\,\mathbf e_r+r\dot\theta\,\mathbf e_\theta,
$$
$$
\mathbf a=(\ddot r-r\dot\theta^2)\mathbf e_r+(r\ddot\theta+2\dot r\dot\theta)\mathbf e_\theta.
$$
推导：对 $\mathbf v=\dot r\mathbf e_r+r\dot\theta\mathbf e_\theta$ 再求导，代入 $\dot{\mathbf e}_r,\dot{\mathbf e}_\theta$ 并合并同类项即得。物理意义：$-r\dot\theta^2$ 为离心项，$2\dot r\dot\theta$ 为科里奥利项。应用：中心力问题、机器人臂运动学。

#### 1.4 牛顿运动定律与动量定理
第二定律原始形式
$$
\mathbf F=\frac{\mathrm d\mathbf p}{\mathrm dt},\qquad \mathbf p=m\mathbf v.
$$
当质量不变时 $\mathbf F=m\mathbf a$。对时间积分得动量定理
$$
\int_{t_1}^{t_2}\mathbf F\,\mathrm dt=m\mathbf v_2-m\mathbf v_1=\Delta\mathbf p,
$$
即冲量等于动量增量。物理意义：力是动量的变化率；内力不改变系统总动量。特例：$\mathbf F=0$ 时 $\mathbf p=\text{const}$（动量守恒）。应用：碰撞、变质量系统（火箭）方程
$$
m\frac{\mathrm d\mathbf v}{\mathrm dt}=\mathbf F_{\text{ext}}+\mathbf u\frac{\mathrm dm}{\mathrm dt}.
$$

#### 1.5 功、动能定理与保守力
力所做的功
$$
W=\int_{\mathbf r_1}^{\mathbf r_2}\mathbf F\cdot\mathrm d\mathbf r.
$$
对 $m\mathbf a=\mathbf F$ 两边点乘 $\mathrm d\mathbf r$ 并积分：
$$
\because\ m\frac{\mathrm d\mathbf v}{\mathrm dt}\cdot\mathrm d\mathbf r=m\mathbf v\cdot\mathrm d\mathbf v
=\mathrm d\!\left(\tfrac12 mv^2\right)
$$
$$
\therefore\ W=\int \mathrm d\!\left(\tfrac12 mv^2\right)=\tfrac12 mv_2^2-\tfrac12 mv_1^2=\Delta T.
$$
保守力定义为做功与路径无关，等价于 $\oint\mathbf F\cdot\mathrm d\mathbf r=0$，从而可引入势能 $\mathbf F=-\nabla V$。物理意义：动能定理是牛顿第二定律的能量积分。应用：功能原理求速度。

#### 1.6 机械能守恒与势能
若仅保守力做功，由动能定理与 $W=-\Delta V$ 得
$$
\Delta T=-\Delta V\ \Rightarrow\ T+V=E=\text{const}.
$$
推导：$W=\Delta T$ 且 $W=\int-\nabla V\cdot\mathrm d\mathbf r=-\Delta V$，两式相等。物理意义：机械能守恒是时间平移对称性的结果。特例：含摩擦时 $E$ 减少，$\Delta E=W_{\text{nc}}$。应用：重力场 $V=mgh$，弹簧 $V=\tfrac12 kx^2$，引力 $V=-GMm/r$。

#### 1.7 角动量、力矩与角动量定理
对固定点 $O$：$\mathbf L=\mathbf r\times\mathbf p$，$\mathbf M=\mathbf r\times\mathbf F$。求导
$$
\dot{\mathbf L}=\dot{\mathbf r}\times\mathbf p+\mathbf r\times\dot{\mathbf p}
=\mathbf v\times m\mathbf v+\mathbf r\times\mathbf F
=\mathbf 0+\mathbf M=\mathbf M.
$$
物理意义：力矩是角动量的变化率；中心力（$\mathbf M=0$）下角动量守恒，即开普勒第二定律（面积速度守恒）。特例：$\mathbf M=0\Rightarrow \mathbf L=\text{const}$。应用：行星运动、刚体转动。

#### 1.8 万有引力与开普勒问题（Binet 方程）
中心力 $\mathbf F=F(r)\mathbf e_r$。由角动量守恒 $L=mr^2\dot\theta=\text{const}$。令 $u=1/r$，可导出轨道方程
$$
\frac{\mathrm d^2u}{\mathrm d\theta^2}+u=-\frac{F(1/u)}{mh^2u^2},\qquad h=\frac{L}{m}=r^2\dot\theta.
$$
推导：$\dot r=\frac{\mathrm dr}{\mathrm d\theta}\dot\theta=-\frac1{u^2}\frac{\mathrm du}{\mathrm d\theta}\frac{L u^2}{m}=-\frac{L}{m}\frac{\mathrm du}{\mathrm d\theta}$，再对时间求导得 $\ddot r$，代入径向方程 $\ddot r-r\dot\theta^2=F/m$ 化简便得。对引力 $F=-GMm/r^2$：
$$
\frac{\mathrm d^2u}{\mathrm d\theta^2}+u=\frac{GM}{h^2}
\ \Rightarrow\ u=\frac{GM}{h^2}+A\cos(\theta-\theta_0),
$$
即圆锥曲线轨道。应用：卫星轨道、椭圆/抛物线/双曲线判别。

### 第2章 刚体运动学

#### 2.1 刚体自由度与欧拉定理
刚体上任两点距离不变，$N$ 个质点的 $3N$ 个坐标受 $\binom{N}{2}$ 个距离约束，净剩 6 个自由度：3 个平动 + 3 个转动。欧拉定理：定点运动刚体的任意位移等价于绕过该点某轴的一次转动。物理意义：刚体构型由一点位置与一个方向（姿态）确定。应用：机构分析、飞行器姿态。

#### 2.2 方向余弦矩阵与欧拉角
体坐标系基 $\{\mathbf e_1,\mathbf e_2,\mathbf e_3\}$ 与空间基 $\{\mathbf E_1,\mathbf E_2,\mathbf E_3\}$ 的关系
$$
\mathbf e_i=\sum_{j}R_{ij}\mathbf E_j,\qquad R=\{R_{ij}\}\in SO(3),\quad R^{\mathrm T}R=RR^{\mathrm T}=I,\ \det R=1.
$$
欧拉角 $(\psi,\theta,\phi)$ 下 $R=R_3(\psi)R_1(\theta)R_3(\phi)$（3-1-3 约定）。$SO(3)$ 条件来自保持内积与右手性：$\mathbf e_i\cdot\mathbf e_j=\delta_{ij}\Rightarrow R^{\mathrm T}R=I$。应用：刚体姿态描述、机器人 DH 参数。

#### 2.3 角速度矢量与泊松方程
体基随刚体转动，存在角速度 $\boldsymbol\omega$ 使
$$
\dot{\mathbf e}_i=\boldsymbol\omega\times\mathbf e_i.
$$
推导：由 $\mathbf e_i$ 为单位正交基，$\dot{\mathbf e}_i$ 与 $\mathbf e_i$ 正交，可写作 $\dot{\mathbf e}_i=\boldsymbol\Omega\,\mathbf e_i$，反对称矩阵 $\boldsymbol\Omega$ 对应矢量 $\boldsymbol\omega$，即 $\boldsymbol\Omega\mathbf x=\boldsymbol\omega\times\mathbf x$。物理意义：$\boldsymbol\omega$ 描述刚体瞬时转动，一般与体轴方向不同。特例：定轴转动 $\boldsymbol\omega=\dot\theta\mathbf k$。

#### 2.4 速度合成定理（刚体速度场）
对刚体上任意点 $P$，相对基点 $O$ 有 $\mathbf r_P=\mathbf r_O+\boldsymbol\rho$，$\boldsymbol\rho$ 在体基中分量不变，则
$$
\mathbf v_P=\dot{\mathbf r}_O+\dot{\boldsymbol\rho}
=\mathbf v_O+\boldsymbol\omega\times\boldsymbol\rho.
$$
物理意义：刚体速度场由平动部分与绕基点的转动部分叠加；$\boldsymbol\omega$ 与基点选择无关。应用：平面运动速度分析、瞬心法（$\mathbf v_P=\mathbf 0$ 之点即瞬心）。

#### 2.5 加速度合成与科里奥利加速度
在固定基与转动基（角速度 $\boldsymbol\omega$）之间，对 $\mathbf v=\mathbf v_O+\boldsymbol\omega\times\boldsymbol\rho+\mathbf v_{\text{rel}}$ 求导得
$$
\mathbf a=\mathbf a_O+\boldsymbol\alpha\times\boldsymbol\rho
+\boldsymbol\omega\times(\boldsymbol\omega\times\boldsymbol\rho)
+2\boldsymbol\omega\times\mathbf v_{\text{rel}}+\mathbf a_{\text{rel}}.
$$
推导：分别对牵连速度、相对速度求导，其中相对矢量的方向变化贡献 $2\boldsymbol\omega\times\mathbf v_{\text{rel}}$（科氏项），向心项来自 $\boldsymbol\omega\times(\boldsymbol\omega\times\boldsymbol\rho)$。物理意义：科氏加速度在旋转参考系中使运动物体偏转。应用：地球自转效应、涡轮机械。

#### 2.6 定点转动与瞬时转轴
刚体绕固定点 $O$ 转动时，其上与 $\boldsymbol\omega$ 平行的点瞬时速度为零，这些点构成瞬时转轴。若已知体上两点速度 $\mathbf v_A,\mathbf v_B$，则
$$
\mathbf v_A-\mathbf v_B=\boldsymbol\omega\times(\mathbf r_A-\mathbf r_B),
$$
两边叉乘 $\mathbf r_A-\mathbf r_B$ 可解出 $\boldsymbol\omega$ 的轴向。推导：由速度场 $\mathbf v=\boldsymbol\omega\times\mathbf r$ 相减即得，再与 $\mathbf r_A-\mathbf r_B$ 做适当组合消去平行分量。物理意义：定点运动任意时刻等价于绕瞬时轴纯转动，转轴方向不断变化。特例：欧拉盘、陀螺的瞬时轴。应用：刚体姿态运动学。

### 第3章 刚体动力学与欧拉方程

#### 3.1 惯性张量与角动量
对固定点 $O$，连续体角动量
$$
\mathbf L=\int \mathbf r\times(\boldsymbol\omega\times\mathbf r)\,\mathrm dm
=I\,\boldsymbol\omega,
$$
其中惯性张量
$$
I=\int\big(r^2\delta_{ij}-r_ir_j\big)\,\mathrm dm
=\begin{pmatrix}
I_{xx}&-I_{xy}&-I_{xz}\\
-I_{xy}&I_{yy}&-I_{yz}\\
-I_{xz}&-I_{yz}&I_{zz}
\end{pmatrix}.
$$
推导：用矢量恒等式 $\mathbf r\times(\boldsymbol\omega\times\mathbf r)=r^2\boldsymbol\omega-(\mathbf r\cdot\boldsymbol\omega)\mathbf r$，写成矩阵形式即得。$I$ 对称正定。物理意义：$I$ 反映质量相对转轴的分布。

#### 3.2 平行轴定理与主惯性轴
质心 $C$ 与任意点 $O$：$\mathbf r=\mathbf R+\boldsymbol\rho$，则
$$
I_O=I_C+m\big(R^2\delta_{ij}-R_iR_j\big).
$$
推导：将 $\int(\mathbf R+\boldsymbol\rho)\times[\,\cdot\,]\mathrm dm$ 展开，交叉项含 $\int\boldsymbol\rho\,\mathrm dm=m\boldsymbol\rho_C=\mathbf 0$。惯性张量可对角化，特征方向即主惯性轴，特征值即主惯性矩。物理意义：绕质心主轴转动最“简单”。应用：陀螺、卫星姿态动力学。

#### 3.3 欧拉动力学方程（完整推导）
角动量定理在体基中：$\big(\dot{\mathbf H}\big)_{\text{space}}=\big(\dot{\mathbf H}\big)_{\text{body}}+\boldsymbol\omega\times\mathbf H=\mathbf M$。取主轴 $\mathbf H=(I_1\omega_1,I_2\omega_2,I_3\omega_3)$，得
$$
\begin{aligned}
I_1\dot\omega_1-(I_2-I_3)\omega_2\omega_3&=M_1,\\
I_2\dot\omega_2-(I_3-I_1)\omega_3\omega_1&=M_2,\\
I_3\dot\omega_3-(I_1-I_2)\omega_1\omega_2&=M_3.
\end{aligned}
$$
推导：$\boldsymbol\omega\times\mathbf H=(\omega_2H_3-\omega_3H_2,\dots)$，代入 $H_i=I_i\omega_i$ 并整理。物理意义：陀螺力矩项 $-(I_j-I_k)\omega_j\omega_k$ 使转动耦合。特例：球对称 $I_1=I_2=I_3$ 得 $\dot{\boldsymbol\omega}=\mathbf 0$。

#### 3.4 自由欧拉陀螺与稳定性
无力矩自由刚体：$I_1\dot\omega_1=(I_2-I_3)\omega_2\omega_3$ 等。以中间轴为例做线性化扰动 $\omega=\boldsymbol\omega_0+\delta$，得
$$
\ddot\delta\sim-\frac{(I_1-I_2)(I_2-I_3)}{I_1I_3}\omega_2^2\delta.
$$
物理意义：绕最大、最小主惯量轴转动稳定，绕中间轴不稳定（网球拍定理）。应用：航天器姿态稳定、德雅尼效应。

#### 3.5 陀螺近似理论
高速自转 $\omega_3\gg\omega_1,\omega_2$ 时，$I_3\dot\omega_3\approx0$，其余两式给出进动
$$
I_3\omega_3\,\boldsymbol\Omega\approx\mathbf M,\qquad \boldsymbol\Omega=\frac{\mathbf M}{I_3\omega_3}.
$$
推导：忽略 $\dot\omega_{1,2}$ 与二阶小量，主轴方程化为 $\mathbf M=\boldsymbol\omega_3\mathbf e_3\times\mathbf H$ 的进动关系。物理意义：陀螺力矩抵抗外矩，产生进动而非倾倒。应用：陀螺仪、自行车稳定性。

#### 3.6 惯量椭球与主惯性矩性质
以 $O$ 为中心作曲面
$$
\mathbf r^{\mathrm T}I\,\mathbf r=I_1x_1^2+I_2x_2^2+I_3x_3^2=1,
$$
其半轴长 $1/\sqrt{I_i}$，称为惯量椭球。绕单位矢量 $\mathbf n$ 的转动惯量 $I(\mathbf n)=\mathbf n^{\mathrm T}I\mathbf n$。推导：将 $I$ 对角化后沿主轴写出二次型，曲面方程即为 $I_1x_1^2+I_2x_2^2+I_3x_3^2=1$。物理意义：惯量椭球直观呈现各方向转动惯量分布；绕主轴分别取得极值。特例：球对称时退化为球面。应用：陀螺动力学、刚体定向。

### 第4章 拉格朗日力学

#### 4.1 广义坐标与约束
系统用 $n$ 个广义坐标 $q_j$ 描述，完整约束 $f(\mathbf r_i,t)=0$ 可消去自由度。位矢表示为 $\mathbf r_i=\mathbf r_i(q_1,\dots,q_n,t)$，速度
$$
\mathbf v_i=\sum_{j=1}^n\frac{\partial\mathbf r_i}{\partial q_j}\dot q_j+\frac{\partial\mathbf r_i}{\partial t}.
$$
物理意义：广义坐标把约束自动编码，避免约束力。特例：$n=3N-k$，$k$ 为完整约束数。

#### 4.2 虚功原理与达朗贝尔原理
对系统施加虚位移 $\delta\mathbf r_i$（满足约束），理想约束力做功为零，则
$$
\sum_i(\mathbf F_i-m_i\mathbf a_i)\cdot\delta\mathbf r_i=0.
$$
定义广义力 $Q_j=\sum_i\mathbf F_i\cdot\dfrac{\partial\mathbf r_i}{\partial q_j}$，有 $\sum_i\mathbf F_i\cdot\delta\mathbf r_i=\sum_jQ_j\delta q_j$。物理意义：动力学问题化为“主动力—惯性力”的静力平衡（达朗贝尔）。应用：复杂机构动力学。

#### 4.3 拉格朗日方程的完整推导
由 $\mathbf v_i=\sum_j\frac{\partial\mathbf r_i}{\partial q_j}\dot q_j+\frac{\partial\mathbf r_i}{\partial t}$，利用
$$
\frac{\partial\mathbf v_i}{\partial\dot q_j}=\frac{\partial\mathbf r_i}{\partial q_j},\qquad
\frac{\mathrm d}{\mathrm dt}\frac{\partial\mathbf r_i}{\partial q_j}=\frac{\partial\mathbf v_i}{\partial q_j},
$$
惯性项化为
$$
\sum_i m_i\mathbf a_i\cdot\frac{\partial\mathbf r_i}{\partial q_j}
=\frac{\mathrm d}{\mathrm dt}\frac{\partial T}{\partial\dot q_j}-\frac{\partial T}{\partial q_j},
\qquad T=\frac12\sum_i m_iv_i^2.
$$
代入达朗贝尔原理得
$$
\frac{\mathrm d}{\mathrm dt}\frac{\partial T}{\partial\dot q_j}-\frac{\partial T}{\partial q_j}=Q_j.
$$
物理意义：用动能与广义力即得运动方程，无需约束力。特例：保守系统 $Q_j=-\partial V/\partial q_j$，令 $L=T-V$ 得欧拉—拉格朗日方程
$$
\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}=0.
$$

#### 4.4 循环坐标与广义动量守恒
若 $\partial L/\partial q_j=0$，则 $q_j$ 为循环坐标，$p_j=\partial L/\partial\dot q_j$ 守恒：
$$
\frac{\mathrm d}{\mathrm dt}p_j=\frac{\partial L}{\partial q_j}=0\ \Rightarrow\ p_j=\text{const}.
$$
物理意义：坐标不显含于拉氏量对应守恒量（诺特思想）。特例：自由粒子 $L=\tfrac12 m\dot x^2$，$p_x=m\dot x$ 守恒。应用：中心力问题角动量守恒。

#### 4.5 拉格朗日乘子与约束力
含非完整/未消约束 $f_\alpha(q,t)=0$ 时，引进乘子 $\lambda_\alpha$：
$$
\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}
=\sum_\alpha\lambda_\alpha\frac{\partial f_\alpha}{\partial q_j}.
$$
推导：在达朗贝尔原理中约束力做功不为零，对每个约束以 $\lambda_\alpha\delta f_\alpha$ 计入虚功，重整即得。物理意义：$\lambda_\alpha$ 正比于约束力。应用：斜面滚动、绳张力求解。

#### 4.6 耗散函数与广义势
有黏性阻力时引入瑞利耗散函数
$$
\mathcal D=\frac12\sum_{i,j}c_{ij}\dot q_i\dot q_j,\qquad Q_j^{\text{diss}}=-\frac{\partial\mathcal D}{\partial\dot q_j}.
$$
对广义势 $V(q,\dot q,t)$ 定义 $Q_j=-\partial V/\partial q_j+\mathrm d(\partial V/\partial\dot q_j)/\mathrm dt$，则
$$
\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}+\frac{\partial\mathcal D}{\partial\dot q_j}=0,\qquad L=T-V.
$$
推导：达朗贝尔原理中把耗散力与广义势对应的广义力分别计入右端，阻力线性时即为 $\mathcal D$ 的偏导。物理意义：$2\mathcal D$ 为功率耗散。特例：$c_{ij}=c\,\delta_{ij}$ 时线性阻尼。应用：阻尼结构动力学、电路类比。

### 第5章 哈密顿力学

#### 5.1 勒让德变换与哈密顿函数
定义广义动量 $p_j=\partial L/\partial\dot q_j$，若 $\det(\partial^2L/\partial\dot q_i\partial\dot q_j)\neq0$ 可反解 $\dot q_j=\dot q_j(q,p,t)$。定义
$$
H(q,p,t)=\sum_j p_j\dot q_j-L.
$$
取全微分
$$
\mathrm dH=\sum_j\Big(\dot q_j\,\mathrm dp_j-\frac{\partial L}{\partial q_j}\mathrm dq_j\Big)-\frac{\partial L}{\partial t}\mathrm dt,
$$
与 $H$ 作为 $(q,p,t)$ 函数的微分比较即得正则方程。物理意义：$H$ 通常为总能量。

#### 5.2 哈密顿正则方程
$$
\dot q_j=\frac{\partial H}{\partial p_j},\qquad
\dot p_j=-\frac{\partial H}{\partial q_j},\qquad
\frac{\partial H}{\partial t}=-\frac{\partial L}{\partial t}.
$$
推导见 5.1 的比较系数。物理意义：$2n$ 个一阶方程，相空间中的流。特例：$H=p^2/2m+V(q)$ 给出 $\dot q=p/m$，$\dot p=-V'(q)$，即牛顿方程。应用：数值积分（辛算法）、统计力学。

#### 5.3 泊松括号与守恒量
$$
\{f,g\}=\sum_j\Big(\frac{\partial f}{\partial q_j}\frac{\partial g}{\partial p_j}
-\frac{\partial f}{\partial p_j}\frac{\partial g}{\partial q_j}\Big).
$$
运动方程
$$
\dot f=\{f,H\}+\frac{\partial f}{\partial t}.
$$
推导：对 $f(q,p,t)$ 求全导数并代入正则方程。物理意义：$f$ 守恒 $\iff\{f,H\}=0$。特例：$\{q_i,p_j\}=\delta_{ij}$。应用：判断守恒量、量子化对应。

#### 5.4 正则变换与生成函数
变换 $(q,p)\to(Q,P)$ 保持正则方程，当且仅当存在生成函数 $F(q,P,t)$ 使
$$
p_j=\frac{\partial F}{\partial q_j},\qquad Q_j=\frac{\partial F}{\partial P_j},\qquad
H'=H+\frac{\partial F}{\partial t}.
$$
推导：要求 $\sum p\,\mathrm dq-H\mathrm dt=\sum P\,\mathrm dQ-H'\mathrm dt+\mathrm dF$，比较微分即得。物理意义：正则变换保持相空间结构（辛结构）。应用：哈密顿—雅可比理论、作用量—角变量。

#### 5.5 哈密顿—雅可比方程
寻找生成函数 $S(q,t)$ 使新哈密顿 $H'=0$，则 $P_j=\text{const}$，$Q_j=\text{const}$，由
$$
p_j=\frac{\partial S}{\partial q_j},\qquad
H'=H\Big(q,\frac{\partial S}{\partial q},t\Big)+\frac{\partial S}{\partial t}=0
$$
得
$$
\frac{\partial S}{\partial t}+H\Big(q,\frac{\partial S}{\partial q},t\Big)=0.
$$
物理意义：一阶偏微分方程，$S$ 为哈密顿主函数，其等值面为波阵面。特例：$H$ 不显含时间时 $S=W(q)-Et$，$W$ 为特征函数。应用：可积系统、几何光学类比。

#### 5.6 刘维尔定理
相空间体积在哈密顿流下不变：
$$
\frac{\mathrm d}{\mathrm dt}\int_\Omega\mathrm dq\,\mathrm dp
=\int_\Omega\nabla_{\text{ph}}\cdot(\dot{\mathbf x})\,\mathrm dq\,\mathrm dp=0.
$$
推导：由正则方程 $\nabla\cdot\dot{\mathbf x}=\sum_j(\partial\dot q_j/\partial q_j+\partial\dot p_j/\partial p_j)=\sum_j(\partial^2H/\partial q_j\partial p_j-\partial^2H/\partial p_j\partial q_j)=0$。物理意义：相空间流不可压缩。应用：统计力学系综。

### 第6章 变分原理与作用量

#### 6.1 泛函与变分
泛函 $J[y]=\int_{x_1}^{x_2}F(x,y,y')\,\mathrm dx$。变分
$$
\delta J=\int\Big(\frac{\partial F}{\partial y}\delta y+\frac{\partial F}{\partial y'}\delta y'\Big)\mathrm dx.
$$
物理意义：变分是函数空间中的“微分”，$\delta y$ 为固定端点的任意小扰动。特例：$F=\sqrt{1+y'^2}$ 给出弧长，泛函极小对应最短路径。

#### 6.2 欧拉—拉格朗日方程（变分法基本引理推导）
对 $\delta J$ 分部积分：
$$
\because\ \int\frac{\partial F}{\partial y'}\frac{\mathrm d(\delta y)}{\mathrm dx}\mathrm dx
=\Big[\frac{\partial F}{\partial y'}\delta y\Big]_{x_1}^{x_2}
-\int\frac{\mathrm d}{\mathrm dx}\Big(\frac{\partial F}{\partial y'}\Big)\delta y\,\mathrm dx,
$$
固定端点 $\delta y(x_1)=\delta y(x_2)=0$，故
$$
\delta J=\int\Big[\frac{\partial F}{\partial y}-\frac{\mathrm d}{\mathrm dx}\Big(\frac{\partial F}{\partial y'}\Big)\Big]\delta y\,\mathrm dx=0.
$$
由变分法基本引理（$\delta y$ 任意且连续），被积函数为零：
$$
\frac{\partial F}{\partial y}-\frac{\mathrm d}{\mathrm dx}\Big(\frac{\partial F}{\partial y'}\Big)=0.
$$
物理意义：泛函极值条件化为微分方程。特例：$F$ 不显含 $x$ 时得贝尔特拉米恒等式 $F-y'\partial F/\partial y'=\text{const}$。

#### 6.3 哈密顿最小作用量原理
定义作用量 $S=\int_{t_1}^{t_2}L(q,\dot q,t)\,\mathrm dt$。真实运动使 $\delta S=0$。由 6.2 取 $x=t$，$y=q$，$F=L$，即得拉格朗日方程
$$
\frac{\mathrm d}{\mathrm dt}\frac{\partial L}{\partial\dot q_j}-\frac{\partial L}{\partial q_j}=0.
$$
物理意义：全部力学由一标量 $L$ 与一个变分原理推出，体现自然界的极值/平稳性。应用：场论、最优控制（变分法）。

#### 6.4 测地线与最短路径
曲面上弧长泛函 $J=\int\sqrt{g_{ij}\dot q^i\dot q^j}\,\mathrm dt$（取 $t=s$ 得 $\int\sqrt{g_{ij}\mathrm dq^i\mathrm dq^j}$）。若取 $F=\tfrac12 g_{ij}\dot q^i\dot q^j$，欧拉—拉格朗日方程给出测地线方程
$$
\ddot q^k+\Gamma^k_{ij}\dot q^i\dot q^j=0,\qquad
\Gamma^k_{ij}=\tfrac12 g^{kl}(\partial_i g_{jl}+\partial_j g_{il}-\partial_l g_{ij}).
$$
推导：将 $F$ 代入欧拉—拉格朗日方程并对 $g$ 求导，与克氏符号定义对照即得。物理意义：自由粒子在弯曲空间走测地线。特例：平面 $g_{ij}=\delta_{ij}$ 得直线。

#### 6.5 诺特定理（对称性与守恒律）
若作用量在单参数变换 $q\to q+\epsilon\,\phi(q,t)$ 下不变，则
$$
\frac{\mathrm d}{\mathrm dt}\Big(\frac{\partial L}{\partial\dot q_j}\phi_j\Big)=0
\quad\Rightarrow\quad Q=\frac{\partial L}{\partial\dot q_j}\phi_j=\text{const}.
$$
推导：令 $q_j(\epsilon)$ 满足变换，$\delta L=\epsilon\,(\partial L/\partial q_j\,\phi_j+\partial L/\partial\dot q_j\,\dot\phi_j)$，由欧拉—拉格朗日方程替换 $\partial L/\partial q_j$ 得全导数
$$
\delta L=\epsilon\frac{\mathrm d}{\mathrm dt}\Big(\frac{\partial L}{\partial\dot q_j}\phi_j\Big),
$$
不变则守恒。物理意义：时间平移$\to$能量，空间平移$\to$动量，转动$\to$角动量。应用：场论守恒流。

#### 6.6 自然边界条件与横截条件
若端点自由（$\delta y$ 在端点任意），分部积分边界项要求
$$
\Big[\frac{\partial F}{\partial y'}\Big]_{x_1}^{x_2}=0\ \Rightarrow\ \frac{\partial F}{\partial y'}=0\ (x=x_1,x_2),
$$
称自然边界条件。若端点 $x_2$ 可变，则横截条件
$$
\Big[F-y'\frac{\partial F}{\partial y'}\Big]_{x=x_2}=0.
$$
推导：保留分部积分中的边界项 $\big[\partial F/\partial y'\,\delta y\big]$，当 $\delta y$ 端点任意时令其为零得自然边界；端点可动时还需对端点变分做一阶展开并令总变分为零。物理意义：边界条件由变分原理自动给出，区别于本质边界条件。应用：最优控制、梁自由端边界。

### 第7章 应力张量与应变张量

#### 7.1 应力矢量与柯西公式
在法向 $\mathbf n$ 的面上，应力矢量 $\mathbf t(\mathbf n)$ 分量 $t_i=\sigma_{ji}n_j$（张量记法）或 $\mathbf t=\boldsymbol\sigma^{\mathrm T}\mathbf n$；由角动量守恒 $\boldsymbol\sigma=\boldsymbol\sigma^{\mathrm T}$，故 $\mathbf t=\boldsymbol\sigma\mathbf n$。
推导（四面体）：对微小四面体受力平衡，面积 $S_i=n_iS$（$S$ 为斜面面积），体积力为高阶小量，
$$
\mathbf t(\mathbf n)S+\sum_i(-S_i)\mathbf t(-\mathbf e_i)=\mathbf 0
\ \Rightarrow\ t_j=\sigma_{ij}n_i.
$$
物理意义：应力是二阶张量，只需一个面上的应力即可求出任意面。特例：静水压 $\boldsymbol\sigma=-p\mathbf I$ 得 $\mathbf t=-p\mathbf n$。

#### 7.2 应力张量对称性推导
对微元立方体求力矩平衡（对角动量方程取体积 $\to0$）：
$$
M_i=\epsilon_{ijk}\int x_j(\sigma_{kl}n_l)\mathrm dS+\int(\mathbf x\times\mathbf f)_i\mathrm dV.
$$
面积分化为体积分后涉及 $\sigma_{jk}-\sigma_{kj}$ 乘 $x$，令体积 $\to0$ 得 $\sigma_{ij}=\sigma_{ji}$。物理意义：无体力矩时应力张量对称，独立分量由 9 减为 6。

#### 7.3 主应力与应力不变量
特征值问题 $\boldsymbol\sigma\mathbf n=\sigma\mathbf n$，特征方程
$$
\det(\boldsymbol\sigma-\sigma\mathbf I)=-\sigma^3+I_1\sigma^2-I_2\sigma+I_3=0,
$$
$$
I_1=\operatorname{tr}\boldsymbol\sigma,\quad
I_2=\tfrac12\big[(\operatorname{tr}\boldsymbol\sigma)^2-\operatorname{tr}\boldsymbol\sigma^2\big],\quad
I_3=\det\boldsymbol\sigma.
$$
物理意义：主应力为极值正应力，对应面无剪应力；$I_k$ 与坐标无关。特例：平面应力得莫尔圆。

#### 7.4 小应变张量的定义与几何意义
位移场 $\mathbf u(\mathbf x)$，定义
$$
\varepsilon_{ij}=\frac12\Big(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\Big).
$$
推导：微元线段 $\mathrm d\mathbf x$ 变形后长变
$$
|\mathrm d\mathbf x+\mathrm d\mathbf u|^2-|\mathrm d\mathbf x|^2
=2\varepsilon_{ij}\mathrm dx_i\mathrm dx_j+O(|\nabla\mathbf u|^2).
$$
物理意义：$\varepsilon_{ii}$ 为沿 $x_i$ 的相对伸长，$\varepsilon_{ij}(i\neq j)$ 为剪切角的一半。特例：刚体转动 $\mathbf u=\boldsymbol\omega\times\mathbf x$ 给出 $\varepsilon_{ij}=0$（反对称部分不产生应变）。应用：应变片、数字图像相关。

#### 7.5 应变协调方程
六个应变分量由三个位移决定，须满足协调条件
$$
\varepsilon_{ij,kl}+\varepsilon_{kl,ij}-\varepsilon_{ik,jl}-\varepsilon_{jl,ik}=0.
$$
推导：对 $\varepsilon_{ij}=\tfrac12(u_{i,j}+u_{j,i})$ 求二阶导并组合，$\partial_k\partial_l\partial_i u_j$ 的对称性使含 $u$ 项抵消。物理意义：保证应变场可积分出单值连续位移，无裂缝/重叠。特例：平面问题化为一个方程
$$
\frac{\partial^2\varepsilon_{xx}}{\partial y^2}+\frac{\partial^2\varepsilon_{yy}}{\partial x^2}=2\frac{\partial^2\varepsilon_{xy}}{\partial x\partial y}.
$$

#### 7.6 体积应变、应变偏张量与应变率张量
体积应变
$$
\varepsilon_V=\varepsilon_{kk}=\frac{\Delta V}{V}\approx\nabla\cdot\mathbf u.
$$
应变偏张量 $e_{ij}=\varepsilon_{ij}-\tfrac13\varepsilon_{kk}\delta_{ij}$，满足 $e_{kk}=0$。由速度场定义应变率张量
$$
\dot\varepsilon_{ij}=\frac12\Big(\frac{\partial v_i}{\partial x_j}+\frac{\partial v_j}{\partial x_i}\Big).
$$
推导：体积 $V=\int\det(\mathbf I+\nabla\mathbf u)\,\mathrm dV_0\approx\int(1+\operatorname{tr}\nabla\mathbf u)\mathrm dV_0$，故 $\Delta V/V=\varepsilon_{kk}$。物理意义：$\varepsilon_V$ 描述体积改变，偏量描述形状畸变；应变率张量用于黏性流体。特例：不可压 $\varepsilon_V=0\Leftrightarrow\nabla\cdot\mathbf u=0$。应用：塑性体积不变、黏性本构。

### 第8章 弹性本构关系与弹性力学基本方程

#### 8.1 广义胡克定律
小变形线弹性材料中应力与应变成线性关系：
$$
\sigma_{ij}=C_{ijkl}\varepsilon_{kl},\qquad \varepsilon_{ij}=S_{ijkl}\sigma_{kl},
$$
四阶刚度张量 $C$ 与柔度张量 $S$ 互逆。由 $\sigma_{ij}=\sigma_{ji}$、$\varepsilon_{kl}=\varepsilon_{lk}$ 得对称性
$$
C_{ijkl}=C_{jikl}=C_{ijlk}=C_{jilk},
$$
独立分量由 81 减为 36；由应变能存在（$C_{ijkl}=C_{klij}$）再减为 21。推导：应变能 $U=\tfrac12\sigma_{ij}\varepsilon_{ij}$ 的势性 $\partial U/\partial\varepsilon_{ij}=\sigma_{ij}$ 迫使 $\partial^2U/\partial\varepsilon_{ij}\partial\varepsilon_{kl}$ 对称。物理意义：本构方程不依赖运动定律，是材料属性。

#### 8.2 各向同性材料与弹性常数
各向同性时 $C_{ijkl}=\lambda\delta_{ij}\delta_{kl}+\mu(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk})$，即
$$
\sigma_{ij}=\lambda\varepsilon_{kk}\delta_{ij}+2\mu\varepsilon_{ij},
$$
$\lambda,\mu$ 为拉梅常数（$\mu$ 即剪切模量 $G$）。反解
$$
\varepsilon_{ij}=\frac{1}{2\mu}\Big(\sigma_{ij}-\frac{\lambda}{3\lambda+2\mu}\sigma_{kk}\delta_{ij}\Big).
$$
推导：对各向同性张量按 $\delta_{ij}\delta_{kl}$ 与 $(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk})$ 两个基展开，代入对称性确定系数，再用逆张量公式。物理意义：两个独立弹性常数刻画各向同性材料。特例：$\sigma_{kk}=3\lambda\varepsilon_{kk}+2\mu\varepsilon_{kk}=(3\lambda+2\mu)\varepsilon_{kk}=3K\varepsilon_{kk}$，体积模量 $K=\lambda+\tfrac23\mu$。

#### 8.3 平衡微分方程
取微元体，应力面积分化为体积分（高斯定理），令体积力 $\mathbf f$：
$$
\sigma_{ji,j}+f_i=0.
$$
推导：对微元受力平衡
$$
\oint_{\partial V}\boldsymbol\sigma\mathbf n\,\mathrm dS+\int_V\mathbf f\,\mathrm dV=\mathbf 0,
$$
利用 $\int_{\partial V}\sigma_{ji}n_j\mathrm dS=\int_V\sigma_{ji,j}\mathrm dV$，因 $V$ 任意得方程。物理意义：静力平衡的局部形式，是椭圆型方程（静力边值问题）。

#### 8.4 位移法：纳维方程
将本构代入平衡并用几何方程 $\varepsilon_{ij}=\tfrac12(u_{i,j}+u_{j,i})$：
$$
\sigma_{ji,j}=\lambda(\varepsilon_{kk})_{,i}+2\mu\varepsilon_{ji,j}
=\lambda u_{k,ki}+\mu(u_{j,ji}+u_{i,jj}).
$$
于是
$$
(\lambda+\mu)u_{k,ki}+\mu u_{i,jj}+f_i=0,
\quad\text{即}\quad (\lambda+\mu)\nabla(\nabla\cdot\mathbf u)+\mu\nabla^2\mathbf u+\mathbf f=\mathbf 0.
$$
物理意义：以位移为未知量的弹性静力基本方程，自动满足协调。特例：无体力时变为 $(\lambda+\mu)\nabla\nabla\cdot\mathbf u+\mu\nabla^2\mathbf u=\mathbf 0$。应用：有限元位移法。

#### 8.5 贝尔特拉米—米歇尔协调方程
应力法中以应力为未知量，须补充协调方程
$$
\nabla^2\sigma_{ij}+\frac{1}{1+\nu}\Theta_{,ij}=0,\qquad \Theta=\sigma_{kk},
$$
无体力时取迹得 $\nabla^2\Theta=0$。推导：将协和方程用胡克定律换为应力，再用平衡方程消去体力项与部分导数，整理即得；取 $i=j$ 缩并得到 $\nabla^2\Theta=0$。物理意义：应力场必须先满足此相容条件才可能由单值位移产生。特例：平面问题退化为艾里应力函数的双调和方程。

#### 8.6 弹性应变能与应变能密度
$$
u=\frac12\sigma_{ij}\varepsilon_{ij}
=\frac12\lambda\varepsilon_{kk}^2+\mu\varepsilon_{ij}\varepsilon_{ij}.
$$
推导：对线弹性体从零应力到终态加载，每增量 $\delta u=\sigma_{ij}\delta\varepsilon_{ij}$，沿线性路径积分得 $\tfrac12\sigma_{ij}\varepsilon_{ij}$。物理意义：弹性变形储存的势能（格林弹性材料正定要求 $\mu>0$，$3\lambda+2\mu>0$）。特例：单轴拉伸 $u=\tfrac12\sigma\varepsilon=\sigma^2/(2E)$。应用：能量法、断裂力学能量释放率。

#### 8.7 弹性常数关系
$$
E=2\mu(1+\nu)=\frac{\mu(3\lambda+2\mu)}{\lambda+\mu},\quad
\nu=\frac{\lambda}{2(\lambda+\mu)},\quad
G=\mu=\frac{E}{2(1+\nu)},
$$
$$
\lambda=\frac{E\nu}{(1+\nu)(1-2\nu)},\quad K=\frac{E}{3(1-2\nu)}.
$$
推导：由单轴应力状态 $\sigma_{11}=\sigma$，其余为零，代入 8.2 反解得 $\varepsilon_{11}=\sigma/E$，$\varepsilon_{22}=\varepsilon_{33}=-\nu\sigma/E$，与定义对照即得。物理意义：工程常数 $(E,\nu)$ 与理论常数 $(\lambda,\mu)$ 的换算。特例：$\nu=0.5$ 不可压，$K\to\infty$。

#### 8.8 边界条件与弹性力学边值问题
位移边界 $\mathbf u=\bar{\mathbf u}$，应力边界 $\boldsymbol\sigma\mathbf n=\bar{\mathbf t}$。完整定解问题由 3 个平衡方程、6 个几何方程、6 个本构方程共 15 个方程与相应边界条件组成（15 个未知量：6 应力 + 6 应变 + 3 位移）。物理意义：弹性静力学是适定的椭圆边值问题。特例：圣维南原理允许用等效静力边界替代局部复杂分布。

### 第9章 能量原理与变分法

#### 9.1 虚功原理与虚位移原理
对弹性体，虚位移 $\delta\mathbf u$ 满足位移边界，则平衡必要充分条件为
$$
\int_V\sigma_{ij}\delta\varepsilon_{ij}\,\mathrm dV=\int_V f_i\delta u_i\,\mathrm dV+\int_{\partial V}\bar t_i\delta u_i\,\mathrm dS.
$$
推导：由平衡 $\sigma_{ji,j}+f_i=0$ 点乘 $\delta u_i$ 积分，分部积分并用高斯定理将应力散度项化为边界与应变项。物理意义：功的互等性与能量格式的出发点。

#### 9.2 最小势能原理
定义总势能
$$
\Pi[\mathbf u]=\int_V\Big(\tfrac12\sigma_{ij}\varepsilon_{ij}-f_iu_i\Big)\mathrm dV-\int_{\partial V}\bar t_iu_i\,\mathrm dS.
$$
真实位移使 $\delta\Pi=0$，且为极小（稳定平衡）。推导：对 $\Pi$ 变分，第一变分给出 9.1 的虚功平衡条件；第二变分 $\delta^2\Pi=\int\delta\varepsilon_{ij}C_{ijkl}\delta\varepsilon_{kl}\,\mathrm dV>0$（应变能正定）保证极小。物理意义：弹性平衡是势能泛函的极小点。应用：瑞利—里兹法、有限元。

#### 9.3 最小余能原理
定义余能
$$
\Pi_c[\boldsymbol\sigma]=\int_V\tfrac12\sigma_{ij}\varepsilon_{ij}\,\mathrm dV-\int_{\partial V_u}\mathbf t\cdot\bar{\mathbf u}\,\mathrm dS,
$$
满足平衡与应力边界之应力场使 $\delta\Pi_c=0$ 且为极小。推导：对 $\Pi_c$ 变分并用协调方程与虚位移原理对偶性即得。物理意义：与最小势能原理对偶，适合应力法。

#### 9.4 卡氏定理
线弹性体应变能对广义力的偏导等于相应广义位移：
$$
U_i=\frac{\partial U}{\partial F_i},\qquad
\delta_i=\frac{\partial U}{\partial M_i}.
$$
推导：由 $U=\tfrac12\sum F_i\delta_i$（线性），对 $F_i$ 求偏导得
$$
\frac{\partial U}{\partial F_i}=\delta_i+\sum_j F_j\frac{\partial\delta_j}{\partial F_i}.
$$
利用柔度矩阵对称性 $\partial\delta_j/\partial F_i=\partial\delta_i/\partial F_j$ 与 $\sum_jF_j\partial\delta_j/\partial F_i=\partial U/\partial F_i$ 自洽，即得。物理意义：求位移的简便能量法。应用：桁架、梁挠度。

#### 9.5 瑞利—里兹法
设近似解 $\mathbf u\approx\sum_{k=1}^N c_k\boldsymbol\phi_k$，代入 $\Pi$ 得代数泛函 $\Pi(c_1,\dots,c_N)$，由
$$
\frac{\partial\Pi}{\partial c_k}=0,\qquad k=1,\dots,N
$$
求系数。物理意义：把无限维变分问题离散为有限维（有限元雏形），得到的是上界逼近。特例：$N=1$ 给出工程近似解。应用：板壳、振动频率估计。

#### 9.6 伽辽金法
取近似解使残差与试探函数正交：
$$
\int_V(\mathcal L\mathbf u-\mathbf f)\cdot\boldsymbol\phi_k\,\mathrm dV=0,
$$
$\mathcal L$ 为平衡算子。物理意义：加权残差法之特例，与虚功原理等价。应用：非自伴问题、流体与结构耦合。

#### 9.7 互等定理（贝蒂—麦克斯韦）
线弹性体上，第一组力在第二组力位移上做功等于第二组力在第一组力位移上做功：
$$
\sum_iF_i^{(1)}\delta_i^{(2)}=\sum_iF_i^{(2)}\delta_i^{(1)}.
$$
推导：由 $U=\tfrac12\mathbf F^{\mathrm T}\mathbf d$ 的双线性与柔度矩阵对称性 $\delta_i=\sum_jC_{ij}F_j$，$C_{ij}=C_{ji}$。物理意义：影响系数对称。应用：求复杂载荷下位移（单位载荷法）。

### 第10章 杆件强度、刚度与组合变形

#### 10.1 轴向拉伸与压缩
等直杆横截面正应力
$$
\sigma=\frac{N}{A},\qquad \varepsilon=\frac{\sigma}{E}=\frac{\Delta l}{l},\qquad \Delta l=\frac{Nl}{EA}.
$$
推导：由平截面假设，轴力 $N=\int_A\sigma\,\mathrm dA$，均匀分布得 $\sigma=N/A$；胡克定律积分得伸长。物理意义：$EA$ 为抗拉刚度。特例：变轴力 $\Delta l=\int_0^l\frac{N(x)}{EA}\mathrm dx$。

#### 10.2 圆轴扭转
圆轴横截面剪应力
$$
\tau=\frac{T\rho}{I_p},\qquad I_p=\frac{\pi d^4}{32},
$$
扭转角
$$
\varphi=\frac{Tl}{GI_p},\qquad \tau_{\max}=\frac{T}{W_p},\ W_p=\frac{\pi d^3}{16}.
$$
推导：设截面刚性转动 $\varphi$，剪应变 $\gamma=\rho\,\mathrm d\varphi/\mathrm dx$，胡克定律 $\tau=G\rho\theta$（$\theta$ 单位长度扭转角），扭矩 $T=\int_A\rho\tau\,\mathrm dA=G\theta I_p$，故 $\theta=T/(GI_p)$，回代得应力。物理意义：$GI_p$ 抗扭刚度。特例：空心圆轴 $I_p=\tfrac{\pi}{32}(D^4-d^4)$。

#### 10.3 梁的弯曲正应力
纯弯曲梁正应力
$$
\sigma=\frac{My}{I_z},\qquad \sigma_{\max}=\frac{M}{W_z},\ W_z=\frac{I_z}{y_{\max}}.
$$
推导：平截面假设 + 中性层假设：纵向应变 $\varepsilon=\frac{y}{\rho}$（$\rho$ 曲率半径），胡克定律 $\sigma=E\frac{y}{\rho}$；轴力为零给中性轴过形心 $\int_A y\mathrm dA=0$；弯矩 $M=\int_A\sigma y\,\mathrm dA=\frac{E}{\rho}I_z$，故 $1/\rho=M/(EI_z)$，回代得 $\sigma=My/I_z$。物理意义：$EI_z$ 抗弯刚度，应力线性分布。特例：矩形 $I_z=bh^3/12$。

#### 10.4 弯曲剪应力与挠曲线
矩形截面剪应力（儒拉夫斯基公式）
$$
\tau=\frac{QS_z^*}{I_zb},
$$
$S_z^*$ 为所求点以外面积对中性轴的静矩。挠曲线微分方程
$$
EI_z\frac{\mathrm d^2w}{\mathrm dx^2}=M(x)\ \Rightarrow\ 
EI_zw''''=q(x).
$$
推导：由微元纵向力平衡导出 $\tau=\frac{QS^*}{I_zb}$；由平面曲线曲率 $1/\rho\approx w''$ 与 $M/EI$ 关系得挠曲线方程，对 $x$ 求两次导并用 $M''=q$ 得四阶方程。物理意义：$q$ 为分布载荷。特例：悬臂/简支边界条件。

#### 10.5 组合变形与强度理论
拉弯组合 $\sigma=\frac{N}{A}+\frac{M_z y}{I_z}+\frac{M_y z}{I_y}$；弯扭组合按第三、第四强度理论：
$$
\sigma_{\text{eq3}}=\sqrt{\sigma^2+4\tau^2},\qquad
\sigma_{\text{eq4}}=\sqrt{\sigma^2+3\tau^2}.
$$
推导：第三强度理论以最大剪应力判据，$\tau_{\max}=\sqrt{(\sigma/2)^2+\tau^2}$，等效应力 $\sigma_{\text{eq3}}=2\tau_{\max}$；第四强度理论以畸变能密度判据，用一般应力状态畸变能表达式代入即得。物理意义：把复杂应力化为单向拉伸等效。应用：轴类零件强度校核。

#### 10.6 平面应力莫尔圆
主应力
$$
\sigma_{1,2}=\frac{\sigma_x+\sigma_y}{2}\pm\sqrt{\Big(\frac{\sigma_x-\sigma_y}{2}\Big)^2+\tau_{xy}^2},
$$
主方向 $\tan2\theta_p=\frac{2\tau_{xy}}{\sigma_x-\sigma_y}$。最大剪应力
$$
\tau_{\max}=\sqrt{\Big(\frac{\sigma_x-\sigma_y}{2}\Big)^2+\tau_{xy}^2}=\frac{\sigma_1-\sigma_2}{2}.
$$
推导：对斜截面 $\sigma_\theta=\sigma_x\cos^2\theta+\sigma_y\sin^2\theta+2\tau_{xy}\sin\theta\cos\theta$，用倍角公式化为圆方程 $\big(\sigma-\frac{\sigma_x+\sigma_y}{2}\big)^2+\tau^2=R^2$。物理意义：应力状态用圆几何表示。应用：实验应力分析。

#### 10.7 压杆稳定的欧拉临界力（见第11章）
见第11章，此处给出细长压杆柔度与临界应力关系
$$
\sigma_{cr}=\frac{\pi^2E}{\lambda^2},\qquad \lambda=\frac{\mu l}{i},\ i=\sqrt{\frac{I}{A}},
$$
$\mu$ 为长度系数。物理意义：长细比越大越易失稳。应用：压杆设计。

### 第11章 稳定性理论与屈曲

#### 11.1 平衡稳定性与李亚普诺夫定义
平衡点 $\mathbf x_0$：对任意 $\epsilon>0$ 存在 $\delta>0$，当 $\|\mathbf x(0)-\mathbf x_0\|<\delta$ 时 $\|\mathbf x(t)-\mathbf x_0\|<\epsilon$，则稳定；若还 $\mathbf x(t)\to\mathbf x_0$ 则渐近稳定。物理意义：小扰动下不无限偏离平衡。特例：简摆最低点稳定、最高点不稳定。

#### 11.2 李亚普诺夫第二方法
若存在正定函数 $V(\mathbf x)$（$V(\mathbf 0)=0$，$V>0$）使 $\dot V\le0$，则原点稳定；$\dot V<0$ 渐近稳定。
$$
\dot V=\nabla V\cdot\dot{\mathbf x}.
$$
推导：$\dot V\le0$ 使 $V$ 沿轨线不增，从而轨线被 $V$ 的等值线限制在原点邻域内。物理意义：能量观点判定稳定性，无需解方程。应用：非线性系统、控制系统。

#### 11.3 欧拉压杆屈曲的完整推导
等直细长杆受轴向压力 $P$，设弯曲挠度 $y(x)$，弯矩 $M=-Py$。由挠曲线方程
$$
EI y''=-Py\ \Rightarrow\ y''+k^2y=0,\qquad k^2=\frac{P}{EI}.
$$
通解 $y=A\sin kx+B\cos kx$。两端铰支边界 $y(0)=0,\ y(l)=0$ 给 $B=0$，$A\sin kl=0$。非平凡解要求 $\sin kl=0$，即 $kl=n\pi$，于是
$$
P_{cr}=\frac{n^2\pi^2EI}{l^2},\qquad n=1,2,\dots
$$
最小临界力 $n=1$：
$$
P_{cr}=\frac{\pi^2EI}{l^2}.
$$
物理意义：分岔点，载荷超过 $P_{cr}$ 直线平衡失稳。特例：其他边界以有效长度 $\mu l$ 替代 $l$。应用：压杆设计、结构稳定。

#### 11.4 不同边界条件下的压杆
一般
$$
P_{cr}=\frac{\pi^2EI}{(\mu l)^2},
$$
两端铰支 $\mu=1$；一端固定一端自由 $\mu=2$；两端固定 $\mu=0.5$；一端固定一端铰支 $\mu\approx0.7$。推导：将各边界条件代入 $y=A\sin kx+B\cos kx+Cx+D$ 的弯矩/挠度条件，得到相应的特征方程 $\tan kl=kl$ 等超越方程，取最低根换算。物理意义：约束越强有效长度越短、临界力越高。

#### 11.5 能量法求临界力
用势能二阶变分：压杆总势能
$$
\Pi=\frac12\int EI(y'')^2\mathrm dx-\frac{P}{2}\int(y')^2\mathrm dx.
$$
临界力由中性稳定 $\delta^2\Pi=0$：
$$
P_{cr}=\min\frac{\int EI(y'')^2\mathrm dx}{\int(y')^2\mathrm dx}.
$$
推导：轴向缩短量 $\Delta=\tfrac12\int(y')^2\mathrm dx$，压力做功 $P\Delta$ 变为势能降低项；令第二变分为零即得瑞利商。物理意义：瑞利商给出上界。应用：复杂边界/变截面压杆。

#### 11.6 初始缺陷与后屈曲
考虑初始挠度 $y_0$，方程
$$
EI y''+P(y+y_0)=0.
$$
令 $y_0=a\sin\frac{\pi x}{l}$，得
$$
y=\frac{P/P_{cr}}{1-P/P_{cr}}\,a\sin\frac{\pi x}{l}.
$$
物理意义：$P\to P_{cr}$ 时挠度发散（渐近），实际结构因缺陷在临界附近发生大变形。特例：$P\ll P_{cr}$ 挠度可忽略。应用：缺陷敏感结构（薄壳）。

### 第12章 流体运动学与连续性方程

#### 12.1 描述方法：拉格朗日与欧拉
拉格朗日法跟踪质点 $\mathbf x=\mathbf x(\mathbf X,t)$；欧拉法给出场 $\mathbf v(\mathbf x,t)$。物理意义：前者便于固体，后者便于流动。特例：稳态流场 $\partial/\partial t=0$。

#### 12.2 物质导数
质点物理量 $\phi$ 的变化率
$$
\frac{\mathrm D\phi}{\mathrm Dt}=\frac{\partial\phi}{\partial t}+(\mathbf v\cdot\nabla)\phi.
$$
推导：由 $\phi=\phi(\mathbf x(t),t)$，链式法则
$$
\frac{\mathrm d\phi}{\mathrm dt}=\frac{\partial\phi}{\partial t}+\frac{\partial\phi}{\partial x_i}\frac{\mathrm dx_i}{\mathrm dt}
=\partial_t\phi+v_i\partial_i\phi.
$$
物理意义：局部变化 + 对流变化。特例：定常场下物质导数 = 对流项。应用：加速度 $\mathbf a=\partial_t\mathbf v+(\mathbf v\cdot\nabla)\mathbf v$。

#### 12.3 雷诺输运定理
对系统（随流物质）广延量 $\Phi=\int_{V_s}\phi\,\mathrm dV$：
$$
\frac{\mathrm D}{\mathrm Dt}\int_{V(t)}\phi\,\mathrm dV
=\int_V\frac{\partial\phi}{\partial t}\mathrm dV+\oint_{\partial V}\phi\mathbf v\cdot\mathbf n\,\mathrm dS.
$$
推导：对 $\Phi(t)$ 在运动控制体上求导，用 $\frac{\mathrm d}{\mathrm dt}\int_{V(t)}\phi\,\mathrm dV=\int_V\dot\phi\,\mathrm dV+\int_V\phi\nabla\cdot\mathbf v\,\mathrm dV$，再由高斯定理把 $\phi\nabla\cdot\mathbf v$ 项转成面积分并合并。物理意义：系统与控制体之间的转换桥梁。

#### 12.4 连续性方程
质量守恒（$\phi=\rho$，系统质量不变）：
$$
\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf v)=0.
$$
推导：由雷诺输运定理，$\frac{\mathrm D}{\mathrm Dt}\int\rho\,\mathrm dV=0$，对任意控制体被积函数为零
$$
\partial_t\rho+\nabla\cdot(\rho\mathbf v)=0.
$$
不可压时 $\rho=$ const：
$$
\nabla\cdot\mathbf v=\frac{\partial v_x}{\partial x}+\frac{\partial v_y}{\partial y}+\frac{\partial v_z}{\partial z}=0.
$$
物理意义：体积守恒（不可压）。应用：管流 $A_1v_1=A_2v_2$。

#### 12.5 流线、迹线与脉线
流线满足 $\mathrm d\mathbf x\times\mathbf v=0$，即
$$
\frac{\mathrm dx}{v_x}=\frac{\mathrm dy}{v_y}=\frac{\mathrm dz}{v_z}.
$$
物理意义：流线为瞬时速度场切线，定常时与迹线重合。特例：均匀流流线为直线。应用：流动显示。

#### 12.6 无旋流动与速度势
涡量 $\boldsymbol\omega=\nabla\times\mathbf v$。若 $\boldsymbol\omega=\mathbf 0$ 则存在速度势 $\phi$，$\mathbf v=\nabla\phi$，不可压时
$$
\nabla\cdot\mathbf v=\nabla^2\phi=0.
$$
推导：由矢量恒等式 $\nabla\times(\nabla\phi)=0$，且单连通域内无旋场必为势场（庞加莱引理）。物理意义：流体运动由拉普拉斯方程支配。特例：均匀流 $\phi=Ux$，点源 $\phi=\frac{Q}{2\pi}\ln r$。应用：势流理论、绕流。

#### 12.7 涡量与开尔文定理
开尔文环量定理：正压理想流体在有势力作用下，沿封闭物质曲线的环量守恒：
$$
\frac{\mathrm D\Gamma}{\mathrm Dt}=0,\qquad \Gamma=\oint\mathbf v\cdot\mathrm d\mathbf l.
$$
推导：$\frac{\mathrm D}{\mathrm Dt}\oint\mathbf v\cdot\mathrm d\mathbf l=\oint\frac{\mathrm D\mathbf v}{\mathrm Dt}\cdot\mathrm d\mathbf l$（其余项抵消），代入欧拉方程 $-\frac1\rho\nabla p-\nabla\Phi$，正压与有势条件下环路积分为零。物理意义：无黏正压有势流动中涡量不生不灭。应用：机翼升力、涡动力学。

### 第13章 流体动力学：欧拉方程、伯努利方程与守恒定律

#### 13.1 欧拉方程（无黏流体运动方程）
对微元流体应用牛顿第二定律 $\rho\frac{\mathrm D\mathbf v}{\mathrm Dt}=-\nabla p+\mathbf f$，展开
$$
\frac{\partial\mathbf v}{\partial t}+(\mathbf v\cdot\nabla)\mathbf v=-\frac1\rho\nabla p+\frac{\mathbf f}{\rho}.
$$
推导：控制体受力来自压强与体积力，$\oint_{\partial V}-p\mathbf n\,\mathrm dS+\int_V\mathbf f\,\mathrm dV=-\int_V\nabla p\,\mathrm dV+\int_V\mathbf f\,\mathrm dV$，用雷诺输运定理取动量守恒并令体积趋于零。物理意义：无黏理想流体动量方程。特例：静流体 $\mathbf v=\mathbf0$ 得 $\nabla p=\mathbf f$（流体静力学）。

#### 13.2 伯努利方程（沿流线推导）
沿流线取弧长微元，将欧拉方程投影到流线方向：
$$
\rho\Big(\frac{\partial v}{\partial t}+v\frac{\partial v}{\partial s}\Big)=-\frac{\partial p}{\partial s}-\rho\frac{\partial\Phi}{\partial s}.
$$
定常、正压、有势、无黏沿流线积分：
$$
\frac{v^2}{2}+\int\frac{\mathrm dp}{\rho}+\Phi=C\quad(\text{沿流线}).
$$
不可压且重力场 $\Phi=gz$：
$$
\frac{v^2}{2g}+\frac{p}{\rho g}+z=H=\text{const}.
$$
物理意义：动能、压能、位能之和不变。特例：文丘里管、托里拆利公式 $v=\sqrt{2gh}$。应用：流量测量、机翼升力（配合环量）。

#### 13.3 动量方程与动量矩方程
对控制体 CV，定常不可压：
$$
\sum\mathbf F=\int_{CS}\rho\mathbf v(\mathbf v\cdot\mathbf n)\,\mathrm dS,
$$
$$
\sum\mathbf M=\int_{CS}\rho(\mathbf r\times\mathbf v)(\mathbf v\cdot\mathbf n)\,\mathrm dS.
$$
推导：对控制体应用雷诺输运定理于动量 $\int_{CV}\rho\mathbf v\,\mathrm dV$ 与角动量，定常时左边为零。物理意义：净出口动量流等于合外力/力矩。特例：管流弯头受力、射流冲击平板 $F=\rho Qv$。应用：涡轮机械、推进。

#### 13.4 能量方程
对控制体，热力学第一定律：
$$
\frac{\mathrm d}{\mathrm dt}\int_{CV}\rho\Big(e+\frac{v^2}{2}\Big)\mathrm dV
=-\oint_{CS}\rho\Big(h+\frac{v^2}{2}\Big)\mathbf v\cdot\mathbf n\,\mathrm dS+\dot Q-\dot W_s,
$$
$h=e+p/\rho$ 为比焓。推导：对总能量应用雷诺输运定理，压力做功与黏性功、导热通量分别计入。物理意义：能量守恒（含热）。特例：绝能流 $h+\frac{v^2}{2}=\text{const}$。应用：喷管、扩压器。

#### 13.5 旋转参考系中的运动方程
在角速度 $\boldsymbol\Omega$ 的旋转系中，绝对加速度
$$
\mathbf a_{\text{abs}}=\mathbf a_{\text{rel}}+2\boldsymbol\Omega\times\mathbf v_{\text{rel}}+\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf r)+\dot{\boldsymbol\Omega}\times\mathbf r.
$$
物理意义：科氏力与离心力为惯性力。特例：地转风平衡（科氏力与气压梯度力平衡）：
$$
-f\mathbf k\times\mathbf v=-\frac1\rho\nabla_hp,
$$
$f=2\Omega\sin\varphi$ 为科氏参数。应用：大气、海洋环流。

#### 13.6 涡量输运方程与伯努利推广
对 $\mathbf v$ 用恒等式 $(\mathbf v\cdot\nabla)\mathbf v=\nabla(\tfrac12v^2)-\mathbf v\times\boldsymbol\omega$，无黏正压有势力时欧拉方程变为
$$
\frac{\partial\mathbf v}{\partial t}+\nabla\Big(\frac{v^2}{2}+\int\frac{\mathrm dp}{\rho}+\Phi\Big)=\mathbf v\times\boldsymbol\omega.
$$
两边取旋度并利用 $\nabla\times(\mathbf v\times\boldsymbol\omega)$ 展开，得
$$
\frac{\mathrm D\boldsymbol\omega}{\mathrm Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf v-\boldsymbol\omega(\nabla\cdot\mathbf v).
$$
对不可压流
$$
\frac{\mathrm D\boldsymbol\omega}{\mathrm Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf v.
$$
物理意义：涡线随流体拉伸、倾斜（涡管强度守恒），故伯努利能量函数在无旋时全局为常数、有旋时仅沿流线为常数。特例：二维无黏流中 $(\boldsymbol\omega\cdot\nabla)\mathbf v=0$，涡量沿迹线守恒。应用：涡环、升力线、涡量法 CFD。

### 第14章 纳维—斯托克斯方程与边界层理论

#### 14.1 牛顿流体的黏性应力
牛顿流体
$$
\tau_{ij}=\mu\Big(\frac{\partial v_i}{\partial x_j}+\frac{\partial v_j}{\partial x_i}\Big)+\lambda(\nabla\cdot\mathbf v)\delta_{ij},
$$
$\mu$ 动力黏度，$\lambda$ 第二黏度，斯托克斯假设 $\lambda=-\tfrac23\mu$。物理意义：应力与应变率线性。特例：不可压 $\nabla\cdot\mathbf v=0$ 时 $\tau_{ij}=\mu(v_{i,j}+v_{j,i})$。

#### 14.2 纳维—斯托克斯方程的完整推导
将应力写成 $\boldsymbol\sigma=-p\mathbf I+\boldsymbol\tau$，代入动量守恒
$$
\rho\frac{\mathrm Dv_i}{\mathrm Dt}=\frac{\partial\sigma_{ji}}{\partial x_j}+f_i.
$$
对牛顿流体且 $\mu$ 为常数：
$$
\frac{\partial\tau_{ji}}{\partial x_j}
=\mu\frac{\partial}{\partial x_j}\Big(\frac{\partial v_j}{\partial x_i}+\frac{\partial v_i}{\partial x_j}\Big)+(\lambda+\mu)\frac{\partial}{\partial x_i}(\nabla\cdot\mathbf v)
=\mu\nabla^2v_i+(\mu+\lambda)\frac{\partial}{\partial x_i}(\nabla\cdot\mathbf v).
$$
于是
$$
\rho\frac{\mathrm D\mathbf v}{\mathrm Dt}=-\nabla p+\mu\nabla^2\mathbf v+(\mu+\lambda)\nabla(\nabla\cdot\mathbf v)+\mathbf f.
$$
不可压时
$$
\rho\Big(\frac{\partial\mathbf v}{\partial t}+(\mathbf v\cdot\nabla)\mathbf v\Big)=-\nabla p+\mu\nabla^2\mathbf v+\mathbf f,
$$
常记运动黏度 $\nu=\mu/\rho$。物理意义：动量输运的完整方程。特例：$\mu=0$ 退化为欧拉方程。应用：CFD、湍流、润滑。

#### 14.3 无量纲化与雷诺数
取特征长度 $L$、速度 $U$，无量纲化 N-S 方程：
$$
\frac{\partial\tilde{\mathbf v}}{\partial\tilde t}+(\tilde{\mathbf v}\cdot\tilde\nabla)\tilde{\mathbf v}
=-\tilde\nabla\tilde p+\frac{1}{Re}\tilde\nabla^2\tilde{\mathbf v},\qquad
Re=\frac{\rho UL}{\mu}.
$$
推导：令 $\mathbf v=U\tilde{\mathbf v}$，$\mathbf x=L\tilde{\mathbf x}$，$t=(L/U)\tilde t$，$p=\rho U^2\tilde p$，代入 N-S 方程除以 $\rho U^2/L$ 即得。物理意义：惯性力与黏性力之比。特例：$Re\to0$ 斯托克斯流（蠕动流），$Re\gg1$ 惯性主导（可形成边界层）。应用：流动相似、模型试验。

#### 14.4 边界层方程
普朗特量级分析：$\delta\ll L$。令 $v_x\sim U$，$v_y\sim U\delta/L$，$\partial_x\sim1/L$，$\partial_y\sim1/\delta$，展开得
$$
\frac{\partial v_x}{\partial x}+\frac{\partial v_y}{\partial y}=0,
$$
$$
v_x\frac{\partial v_x}{\partial x}+v_y\frac{\partial v_x}{\partial y}
=-\frac1\rho\frac{\partial p}{\partial x}+\nu\frac{\partial^2v_x}{\partial y^2},
$$
且 $\partial p/\partial y\approx0$。推导：比较各项量级保留 $O(1)$ 项，黏性项由 $\partial^2/\partial x^2$ 降为 $\partial^2/\partial y^2$，压力由外部势流给定。物理意义：薄层内黏性不可忽略。应用：摩擦阻力、转捩。

#### 14.5 布拉修斯相似解
对平板层流边界层，引入流函数 $\psi=\sqrt{\nu Ux}\,f(\eta)$，$\eta=y\sqrt{U/(\nu x)}$，得布拉修斯方程
$$
f'''+\tfrac12 ff''=0,\quad f(0)=f'(0)=0,\ f'(\infty)=1.
$$
推导：由速度 $v_x=\partial\psi/\partial y=Uf'(\eta)$，$v_y=-\partial\psi/\partial x$，代入边界层方程，各项化归为 $\eta$ 的常微分方程。物理意义：$\delta\sim\sqrt{\nu x/U}$，壁面剪应力
$$
\tau_w=\mu\frac{\partial v_x}{\partial y}\Big|_{y=0}=\rho U^2\frac{0.332}{\sqrt{Re_x}}.
$$
应用：摩擦阻力估算。

#### 14.6 泊肃叶流动与库埃特流动
圆管层流泊肃叶流：轴向动量方程
$$
0=-\frac{\mathrm dp}{\mathrm dx}+\mu\frac1r\frac{\mathrm d}{\mathrm dr}\Big(r\frac{\mathrm dv}{\mathrm dr}\Big).
$$
积分两次并用壁面与轴对称条件得
$$
v(r)=\frac{1}{4\mu}\Big(-\frac{\mathrm dp}{\mathrm dx}\Big)(R^2-r^2),\qquad
Q=\frac{\pi R^4}{8\mu}\Big(-\frac{\mathrm dp}{\mathrm dx}\Big).
$$
平行板库埃特流：移动板带动，$v(y)=U y/h$（无压差），叠加压差得
$$
v(y)=\frac{U y}{h}-\frac{1}{2\mu}\frac{\mathrm dp}{\mathrm dx}\big(hy-y^2\big).
$$
推导：将 N-S 方程取定常、一维、不可压简化，直接积分。物理意义：$Q\propto R^4$（哈根—泊肃叶定律）。应用：黏度测量、微流控。

#### 14.7 斯托克斯流与绕球阻力
$Re\ll1$ 时惯性项可忽略：
$$
\nabla p=\mu\nabla^2\mathbf v,\qquad \nabla\cdot\mathbf v=0.
$$
球体绕流解给出斯托克斯阻力
$$
F_d=6\pi\mu a U.
$$
推导：以球坐标流函数求解 $E^4\psi=0$（$E^2$ 为斯托克斯算子），施加远场与球面无滑移条件，积分应力得阻力。物理意义：黏性主导，$F_d\propto U$（线性）。特例：加入惯性修正的奥森解 $C_d=\frac{24}{Re}(1+\tfrac{3}{16}Re)$。应用：沉降、微颗粒运动。

### 第15章 振动理论

#### 15.1 单自由度无阻尼自由振动
弹簧—质量系统，牛顿第二定律
$$
m\ddot x=-kx\ \Rightarrow\ \ddot x+\omega_n^2x=0,\qquad \omega_n=\sqrt{\frac km}.
$$
通解
$$
x(t)=A\cos(\omega_nt-\varphi),\qquad A=\sqrt{x_0^2+(\dot x_0/\omega_n)^2},\qquad \tan\varphi=\frac{\dot x_0}{\omega_nx_0}.
$$
推导：特征方程 $r^2+\omega_n^2=0$，$r=\pm i\omega_n$，取实部即得。物理意义：等时简谐运动，周期 $T=2\pi/\omega_n$。特例：能量 $E=\tfrac12kA^2$ 守恒。应用：固有频率估算。

#### 15.2 单自由度有阻尼自由振动
$$
m\ddot x+c\dot x+kx=0\ \Rightarrow\ \ddot x+2\zeta\omega_n\dot x+\omega_n^2x=0,
$$
$\zeta=\frac{c}{2\sqrt{km}}$。特解 $x=e^{-\zeta\omega_nt}(A\cos\omega_dt+B\sin\omega_dt)$，$\omega_d=\omega_n\sqrt{1-\zeta^2}$。推导：特征方程 $r^2+2\zeta\omega_nr+\omega_n^2=0$，$r=-\zeta\omega_n\pm i\omega_n\sqrt{1-\zeta^2}$。物理意义：对数衰减率 $\delta=\ln\frac{x(t)}{x(t+T_d)}=\frac{2\pi\zeta}{\sqrt{1-\zeta^2}}$。特例：$\zeta=1$ 临界阻尼，$\zeta>1$ 过阻尼。应用：减振、地震仪。

#### 15.3 受迫振动与共振
简谐激励 $F_0\cos\Omega t$：
$$
\ddot x+2\zeta\omega_n\dot x+\omega_n^2x=\frac{F_0}{m}\cos\Omega t.
$$
稳态解 $x=X\cos(\Omega t-\phi)$，其中
$$
X=\frac{F_0/k}{\sqrt{(1-r^2)^2+(2\zeta r)^2}},\qquad
\tan\phi=\frac{2\zeta r}{1-r^2},\qquad r=\frac{\Omega}{\omega_n}.
$$
推导：设特解 $x=Xe^{i(\Omega t-\phi)}$ 代入复方程，得 AM 放大因子 $H(\Omega)=\frac{1}{1-r^2+i2\zeta r}$，取模与辐角。物理意义：$\Omega\approx\omega_n$ 共振，$\zeta=0$ 时振幅发散；$r=\sqrt{1-2\zeta^2}$ 处峰值。特例：基础激励 $X/Y=\sqrt{\frac{1+(2\zeta r)^2}{(1-r^2)^2+(2\zeta r)^2}}$（隔振）。应用：动力吸振器、隔振设计。

#### 15.4 多自由度系统的振动
以拉格朗日方程建立
$$
M\ddot{\mathbf q}+K\mathbf q=\mathbf 0,
$$
$M$ 对称正定，$K$ 对称。设 $\mathbf q=\mathbf u e^{i\omega t}$ 得广义特征值问题
$$
(K-\omega^2M)\mathbf u=\mathbf 0\ \Rightarrow\ \det(K-\omega^2M)=0.
$$
物理意义：$n$ 个固有频率与模态。特例：两自由度频率
$$
\omega^2=\frac12\Big(\frac{k_1+k_2}{m_1}+\frac{k_2}{m_2}\Big)\pm\sqrt{\dots}.
$$
应用：结构模态分析。

#### 15.5 模态叠加与正交性
模态满足正交关系
$$
\mathbf u_i^{\mathrm T}M\mathbf u_j=m_i\delta_{ij},\qquad \mathbf u_i^{\mathrm T}K\mathbf u_j=k_i\delta_{ij}.
$$
推导：对 $(K-\omega_i^2M)\mathbf u_i=\mathbf 0$ 左乘 $\mathbf u_j^{\mathrm T}$，利用 $M,K$ 对称，交换 $i,j$ 相减得 $(\omega_i^2-\omega_j^2)\mathbf u_j^{\mathrm T}M\mathbf u_i=0$。物理意义：可用模态坐标解耦。特例：任意响应 $\mathbf q=\sum_i\eta_i(t)\mathbf u_i$，各 $\eta_i$ 独立单自由度方程。应用：振型叠加法、有限元。

#### 15.6 连续系统振动（杆的纵振）
均质杆纵振波动方程
$$
\rho A\frac{\partial^2u}{\partial t^2}=EA\frac{\partial^2u}{\partial x^2}\ \Rightarrow\ \frac{\partial^2u}{\partial t^2}=c^2\frac{\partial^2u}{\partial x^2},\quad c=\sqrt{E/\rho}.
$$
推导：取微元受力 $N(x+\mathrm dx)-N(x)$，$N=EA\partial u/\partial x$，牛顿第二定律即得。分离变量 $u=X(x)T(t)$ 给
$$
\ddot T+\omega^2T=0,\qquad X''+\frac{\omega^2}{c^2}X=0.
$$
两端自由杆频率 $\omega_n=n\pi c/l$。物理意义：无限自由度连续模态。应用：超声、结构声学。

#### 15.7 梁的弯曲振动
欧拉—伯努利梁
$$
\rho A\frac{\partial^2w}{\partial t^2}+EI\frac{\partial^4w}{\partial x^4}=0.
$$
设 $w=X(x)e^{i\omega t}$ 得
$$
X''''-\beta^4X=0,\quad \beta^4=\frac{\rho A\omega^2}{EI},
$$
通解 $X=C_1\cos\beta x+C_2\sin\beta x+C_3\cosh\beta x+C_4\sinh\beta x$，由边界条件得频率方程；简支梁 $\omega_n=(n\pi)^2\sqrt{EI/(\rho Al^4)}$。推导：将 $w$ 代入并分离变量，特征方程 $r^4=\beta^4$。物理意义：频率与 $\sqrt{EI/\rho A}/l^2$ 成正比。应用：桥梁、叶片振动。

#### 15.8 参数振动与稳定性（马蒂厄方程）
当刚度或长度周期变化，运动方程化为马蒂厄型
$$
\ddot x+\big(\delta+\varepsilon\cos t\big)x=0.
$$
在参数共振区（如 $\Omega\approx2\omega_n$）振幅指数增长。推导：将周期系数展开并作摄动/弗洛凯分析，得不稳定区（ tongues）。物理意义：周期参数激励可在非共振频率引发失稳。特例：$\varepsilon\to0$ 主不稳定区收缩至 $\Omega=2\omega_n$。应用：转子、山体滑坡、秋千驱动。

---
⬆ 上级领域：[[机械与工程基础]]
