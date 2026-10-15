---
课程名称: 流体力学核心贯通课I
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

# 流体力学核心贯通课I


### 第1章 连续介质假设与流体物理属性

**1. 流体的连续介质假设**

定义：把由大量离散分子组成的流体视为连续填充空间、其物理量（密度、速度、压强、温度）为空间与时间的连续可微函数。设分子平均自由程为 $\lambda$，特征宏观尺度为 $L$，则连续介质假设成立的条件为
$$ \lambda \ll L . $$
假设：在特征体积 $\Delta V$（称为质点体积）内分子数 $N$ 足够大，但 $\Delta V$ 相对宏观尺度足够小。密度定义为
$$ \rho(\mathbf{x},t)=\lim_{\Delta V\to \Delta V_0}\frac{\Delta m}{\Delta V}, $$
其中 $\Delta V_0$ 是使 $\Delta m/\Delta V$ 趋于稳定的最小体积。

推导：设体积 $\Delta V$ 内平均分子数为 $\bar N=\rho \Delta V/m_m$（$m_m$ 为分子质量）。密度统计涨落相对量为 $\sim 1/\sqrt{\bar N}$。要求涨落 $\Delta\rho/\rho \ll 1$，即
$$ \bar N=\frac{\rho\Delta V}{m_m}\gg 1 \Rightarrow \Delta V\gg \frac{m_m}{\rho}, $$
同时局部性要求 $\Delta V^{1/3}\ll L$。两式相容要求 $\lambda\ll L$。$\square$

物理意义：连续介质假设把分子运动“抹平”，使微积分工具可用；它把流体力学与分子动理论联系起来。极限：稀薄气体 $\lambda\sim L$ 时假设失效，须用玻尔兹曼方程或 DSMC。典型应用：常规空气、水在宏观尺度下均满足；高空稀薄大气、微通道稀薄流动需修正。

**2. 密度、比容与重度**

定义：密度 $\rho=m/V$；比容 $v=1/\rho$；重度 $\gamma=\rho g$。对气体由状态方程控制。

推导：由质量守恒，同一种流体均质时 $m=\rho V$，微分 $dm=0\Rightarrow d\rho/\rho+dV/V=0$，即体积压缩与密度增加相反。

意义：$\rho$ 是欧拉方程、N-S方程中的惯性标度。极限：不可压缩流体 $\rho=\text{const}$。应用：水 $\rho\approx 1000\,\mathrm{kg/m^3}$，海平面空气 $\rho\approx 1.225\,\mathrm{kg/m^3}$。

**3. 压缩性与体积模量**

定义：体积模量
$$ K=-\frac{dp}{dV/V}=\rho\frac{dp}{d\rho}. $$
推导：由 $v=1/\rho\Rightarrow dv=-\ d\rho/\rho^2$，故
$$ \frac{dV}{V}=\frac{d\rho/\rho^2}{1/\rho}\cdot(-1)^{-1}\ \Rightarrow\ \frac{dV}{V}=-\frac{d\rho}{\rho}, $$
所以 $K=dp/(d\rho/\rho)=\rho dp/d\rho$。$\square$

意义：$K$ 越大越难压缩；$dp=\rho c^2 d\rho/\rho$ 引入声速 $c=\sqrt{\partial p/\partial\rho}$。极限：$K\to\infty$ 为不可压缩。应用：水锤、声学、气动。

**4. 热膨胀系数**

定义：
$$ \beta=-\frac{1}{\rho}\left(\frac{\partial\rho}{\partial T}\right)_p=\frac{1}{V}\left(\frac{\partial V}{\partial T}\right)_p. $$
推导：$\rho=m/V$，固定 $m$ 对 $T$ 求导：
$$ \left(\frac{\partial\rho}{\partial T}\right)_p=-\frac{m}{V^2}\left(\frac{\partial V}{\partial T}\right)_p=-\rho\frac{1}{V}\left(\frac{\partial V}{\partial T}\right)_p, $$
移项即得。$\square$

意义：浮力对流、热对流的驱动源（瑞利-贝纳德）。方程 $\rho=\rho_0[1-\beta(T-T_0)]$ 用于布辛涅斯克近似。应用：自然对流、大气稳定度。

**5. 粘性与牛顿内摩擦定律**

定义：流体抵抗剪切变形的性质。牛顿内摩擦定律（一维）：
$$ \tau=\mu\frac{du}{dy}, $$
$\mu$ 为动力粘度，运动粘度 $\nu=\mu/\rho$。

推导（动量输运）：两层流体间，分子随机运动穿越界面，把上层大动量输运到下层。设分子平均速度热运动为 $\bar c$，平均自由程 $\lambda$，速度梯度 $du/dy$，则净动量通量（切应力）
$$ \tau\sim \rho \bar c\lambda\frac{du}{dy}. $$
与 $\tau=\mu du/dy$ 比较得
$$ \mu\sim\rho\bar c\lambda. \quad\square $$
由气体动理论更严格有 $\mu=\tfrac13\rho\bar c\lambda$。

意义：粘性源于分子动量交换（气体）与分子间引力（液体）。极限：$\mu=0$ 理想流体。应用：润滑、边界层、管流摩阻。

**6. 牛顿流体与非牛顿流体**

定义：凡满足 $\tau=\mu\dot\gamma$（$\dot\gamma=du/dy$ 为剪切率）且 $\mu$ 与 $\dot\gamma$ 无关者为牛顿流体。更一般本构 $\tau=f(\dot\gamma)$。常见模型：
- 幂律流体 $\tau=K\dot\gamma^n$（$n<1$ 剪切变稀，$n>1$ 剪切增稠）；
- 宾汉塑性体 $\tau=\tau_y+\mu_p\dot\gamma$（$\tau<\tau_y$ 不流动）。

推导（表观粘度）：对幂律流体定义 $\mu_{\rm app}=\tau/\dot\gamma=K\dot\gamma^{n-1}$，故 $n<1$ 时 $\mu_{\rm app}$ 随 $\dot\gamma$ 减小。应用：血液、聚合物熔体、钻井泥浆。

**7. 表面张力**

定义：单位长度上的表面收缩力 $\sigma$（N/m），或单位面积表面自由能。杨-拉普拉斯方程（球面液滴）：
$$ \Delta p=p_{\rm in}-p_{\rm out}=\frac{2\sigma}{R}. $$
推导：对半径 $R$ 的液滴，内压做功使半径增大 $dR$，表面能增加 $dE=\sigma\,dA=\sigma\,8\pi R\,dR$，压力做功 $dW=\Delta p\,dV=\Delta p\,4\pi R^2 dR$，平衡 $dW=dE$：
$$ \Delta p\,4\pi R^2 dR=8\pi\sigma R\,dR\Rightarrow \Delta p=\frac{2\sigma}{R}. \quad\square $$
一般曲面 $\Delta p=\sigma\left(\frac1{R_1}+\frac1{R_2}\right)$。极限：$R\to\infty$ 平界面无压差。应用：毛细、液滴、气泡。

**8. 理想气体状态方程与声速**

定义：
$$ p=\rho R_g T,\qquad R_g=\frac{R_u}{M}. $$
声速定义为等熵小扰动传播速度
$$ c=\sqrt{\left(\frac{\partial p}{\partial\rho}\right)_s}. $$
推导：等熵关系 $p/\rho^\gamma=\text{const}$，取对数微分 $dp/p=\gamma d\rho/\rho$，故
$$ \left(\frac{\partial p}{\partial\rho}\right)_s=\gamma\frac p\rho=\gamma R_g T\Rightarrow c=\sqrt{\gamma R_g T}. \quad\square $$
意义：马赫数 $Ma=u/c$ 度量可压缩性。极限：$Ma<0.3$ 近似不可压。应用：喷管、激波、气动设计。

**9. 输运性质：导热与扩散**

定义：傅里叶导热定律 $\mathbf q=-k\nabla T$；菲克扩散定律 $\mathbf J=-D\nabla C$；牛顿粘性定律张量 $\boldsymbol\tau=\mu(\nabla\mathbf u+\nabla\mathbf u^{\rm T})+\lambda(\nabla\cdot\mathbf u)\mathbf I$。

推导（类比）：三者均为梯度驱动的线性输运，通量 $\propto$ 梯度，输运系数分别为 $k,D,\mu$。由分子动理论 $k\sim\tfrac13\rho c_v\bar c\lambda$，$D\sim\tfrac13\bar c\lambda$，故 $\mu/(\rho D)\sim 1$，$k/(\mu c_v)\sim1$（普朗特数 $Pr=\mu c_p/k$ 量级 1）。应用：热边界层、传质、气动加热。

### 第2章 流体静力学

**1. 静压强与欧拉平衡方程**

定义：静止流体中一点压强各向同性，即应力张量为 $\boldsymbol\sigma=-p\mathbf I$。平衡条件为表面力与体积力平衡。

推导：取四面体微元，斜面上的应力为 $\mathbf t^{(n)}$，三个坐标面上应力分别为 $\mathbf t^{(x)},\mathbf t^{(y)},\mathbf t^{(z)}$。动量平衡（惯性为零）：
$$ \mathbf t^{(n)}dA_n+\mathbf t^{(x)}dA_x+\mathbf t^{(y)}dA_y+\mathbf t^{(z)}dA_z+\mathbf f\rho\,dV=0. $$
令微元尺度 $\to0$，体积力项为高阶小量。由 $dA_x=n_x dA_n$ 等，得
$$ \mathbf t^{(n)}=n_x\mathbf t^{(x)}+n_y\mathbf t^{(y)}+n_z\mathbf t^{(z)}\Rightarrow \mathbf t^{(n)}=\boldsymbol\sigma\cdot\mathbf n. $$
各向同性条件下 $\mathbf t^{(n)}=-p\mathbf n$。再对流体微元六面体列力平衡：$x$ 方向
$$ p\,dy\,dz-\left(p+\frac{\partial p}{\partial x}dx\right)dy\,dz+f_x\rho\,dx\,dy\,dz=0, $$
$$ \Rightarrow -\frac{\partial p}{\partial x}+f_x\rho=0. $$
三方向合并：
$$ \nabla p=\rho\mathbf f. \quad\square $$

物理意义：静止流体中压强梯度与单位体积体积力平衡。极限：$\mathbf f=0\Rightarrow p=$ 常数（无重力自由球）。应用：重力场、离心分离。

**2. 等压面**

定义：$p=$ 常数之面。由 $\nabla p=\rho\mathbf f$ 得等压面法向平行于 $\mathbf f$。

推导：沿等压面取切向 $\delta\mathbf l$，则 $dp=\nabla p\cdot\delta\mathbf l=0\Rightarrow \rho\mathbf f\cdot\delta\mathbf l=0$，即 $\mathbf f\perp\delta\mathbf l$。故等压面垂直于体积力，且两流体分界面为等压面。应用：自由面、连通器。

**3. 重力场中静压强分布**

由 $\nabla p=\rho\mathbf g$，取 $z$ 向上，$\mathbf g=-g\mathbf e_z$：
$$ \frac{dp}{dz}=-\rho g\Rightarrow p(z)=p_0+\int_z^{z_0}\rho g\,dz'. $$
不可压缩流体 $\rho=$ 常数：
$$ p=p_0+\rho g(z_0-z),\qquad p+\rho g z=\text{const}. $$
推导：直接积分 $dp=-\rho g\,dz$。意义：测压管高度 $h=p/(\rho g)$。可压缩等温大气：$p=\rho R_g T_0$，$dp/p=-g dz/(R_gT_0)$，积分
$$ p=p_0\exp\!\left[-\frac{g(z-z_0)}{R_gT_0}\right]. $$
应用：气压高度公式、水坝压力。

**4. 绝对压强、表压与真空度**

定义：绝对压强 $p_{\rm abs}$ 以真空为零；表压 $p_{\rm gage}=p_{\rm abs}-p_{\rm atm}$；真空度 $p_{\rm vac}=p_{\rm atm}-p_{\rm abs}$。推导：三者由参考点定义直接给出。注意液体汽化当 $p_{\rm abs}\le p_v$，产生空化。

**5. 测压原理**

U 形管压差计：两测点压差
$$ p_1-p_2=(\rho_m-\rho)g h, $$
$\rho_m$ 为指示液密度，$h$ 为液柱差。推导：从点1沿测压管内流体至点2列等压面：在指示液等压面上
$$ p_1+\rho g z_1=p_2+\rho g z_2+(\rho_m-\rho)g h \Rightarrow p_1-p_2=\rho g(z_2-z_1)+(\rho_m-\rho)gh, $$
对等高测点即得上式。微压计（斜管）放大倍数 $1/\sin\theta$。应用：风洞测压、皮托管。

**6. 平面上的静压力与压力中心**

总压力大小
$$ F=\int_A p\,dA=\left(p_0+\rho g h_C\right)A, $$
作用点（压力中心）深度
$$ y_{cp}=\frac{\int_A y\,p\,dA}{\int_A p\,dA}=y_C+\frac{I_C}{y_CA}, $$
其中 $I_C$ 为对通过形心水平轴（平行于交线）的惯性矩，$y_C$ 为形心沿斜面的坐标。

推导：设平板与自由面夹角 $\theta$，深度 $h=y\sin\theta$。取 $p_0=0$，$p=\rho g y\sin\theta$。则
$$ F=\int_A \rho g y\sin\theta\,dA=\rho g\sin\theta\,\bar y A=\rho g h_C A. $$
对 $O$（延长面与自由面交线）取矩：
$$ F\,y_{cp}=\int_A y\,p\,dA=\rho g\sin\theta\int_A y^2\,dA. $$
由平行轴定理 $\int y^2 dA=I_O=I_C+A y_C^2$，故
$$ y_{cp}=\frac{\rho g\sin\theta(I_C+Ay_C^2)}{\rho g\sin\theta\,\bar y A}=y_C+\frac{I_C}{y_CA}. \quad\square $$
意义：压力中心总在形心之下。特例：矩形宽 $b$ 高 $h$，$I_C=bh^3/12$，$y_{cp}=y_C+h^2/(12y_C)$。应用：闸门、船体。

**7. 曲面总压力与浮力**

曲面上的水平分力等于曲面在铅垂投影面上的压力；铅垂分力等于曲面所围曲边柱体内液重（或“压力体”重量）：
$$ F_x=\rho g h_{C}A_x,\qquad F_z=\rho g V_{\rm body}. $$
浮力（阿基米德原理）：
$$ F_B=\rho_f g V_{\rm disp}. $$
推导：对浸没物体表面积分 $-\oint p\,\mathbf n\,dA$，用高斯定理 $\oint p\,\mathbf n\,dA=\int_V\nabla p\,dV=\int_V\rho_f\mathbf g\,dV$，故
$$ \mathbf F_B=-\int_V\rho_f\mathbf g\,dV=-\rho_f V\mathbf g\ (\text{向上}). \quad\square $$
意义：浮力为压力梯度的体积分；与物体自身密度无关。应用：船、气球、浮力稳定性。

**8. 相对平衡**

匀速直线运动与静止等效。匀加速直线运动，设加速度 $\mathbf a$，则
$$ \nabla p=\rho(\mathbf g-\mathbf a), $$
自由面垂直于有效重力 $\mathbf g-\mathbf a$。绕铅垂轴以 $\omega$ 旋转：
$$ dp=\rho\left[\omega^2 r\,dr-g\,dz\right]\Rightarrow p=\frac12\rho\omega^2r^2-\rho gz+C, $$
自由面为旋转抛物面 $z=\dfrac{\omega^2r^2}{2g}$。

推导：非惯性系中体积力 $\mathbf f=\mathbf g-\mathbf a$（含离心惯性力 $\omega^2r\,\mathbf e_r$）。代入 $\nabla p$ 积分即得。应用：离心分离机、旋转铸件、液面测转速。

### 第3章 流体运动学

**1. 拉格朗日与欧拉描述**

定义：拉格朗日以质点标记 $\mathbf x=\mathbf x(\mathbf a,t)$ 追踪；欧拉以场 $\mathbf u(\mathbf x,t)$ 描述。质点位置满足
$$ \frac{d\mathbf x}{dt}=\mathbf u(\mathbf x,t),\qquad \mathbf x(\mathbf a,0)=\mathbf a. $$
推导：迹线微分方程 $\dfrac{dx}{u}=\dfrac{dy}{v}=\dfrac{dz}{w}=dt$。

**2. 物质导数**

定义：
$$ \frac{D\phi}{Dt}=\frac{\partial\phi}{\partial t}+\mathbf u\cdot\nabla\phi. $$
推导：$\phi=\phi(\mathbf x(t),t)$，链式法则
$$ \frac{d\phi}{dt}=\frac{\partial\phi}{\partial t}+\frac{\partial\phi}{\partial x}\frac{dx}{dt}+\cdots=\frac{\partial\phi}{\partial t}+\mathbf u\cdot\nabla\phi. \quad\square $$
速度的物质导数即加速度
$$ \mathbf a=\frac{D\mathbf u}{Dt}=\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u. $$
意义：局部导数（非定常）$+$ 对流导数（非线性）。应用：N-S 方程的惯性项。

**3. 迹线、流线与脉线**

定义：流线处处与速度相切，$\dfrac{d\mathbf x}{u}=\dfrac{d\mathbf y}{v}=\dfrac{d\mathbf z}{w}$；迹线为质点轨迹；脉线为同一源点连续质点的连线。定常流三者重合。推导：流线相切条件 $d\mathbf x\parallel\mathbf u$ 展开即上式。应用：烟线、染色实验。

**4. 亥姆霍兹速度分解**

定理：任一局部速度场可分解为平移、旋转、变形：
$$ \mathbf u(\mathbf x+\delta\mathbf x)=\mathbf u(\mathbf x)+\boldsymbol\Omega\times\delta\mathbf x+\mathbf S\cdot\delta\mathbf x, $$
其中 $\boldsymbol\Omega=\tfrac12\nabla\times\mathbf u$ 为角速度矢量，$\mathbf S=\tfrac12(\nabla\mathbf u+\nabla\mathbf u^{\rm T})$ 为变形率张量。

推导：泰勒展开 $\mathbf u(\mathbf x+\delta\mathbf x)=\mathbf u(\mathbf x)+\nabla\mathbf u\cdot\delta\mathbf x$。把 $\nabla\mathbf u=\mathbf S+\mathbf A$，$\mathbf A=\tfrac12(\nabla\mathbf u-\nabla\mathbf u^{\rm T})$ 反对称。反对称张量作用等价于叉乘：$\mathbf A\cdot\delta\mathbf x=\boldsymbol\Omega\times\delta\mathbf x$，且 $\mathbf A$ 分量 $A_{ij}=\tfrac12(\partial_i u_j-\partial_j u_i)=\tfrac12\epsilon_{ijk}\omega_k$。故得分解。$\square$ 意义：变形由对称部分贡献，旋转由反对称部分贡献。

**5. 应变率与体积膨胀率**

定义：线应变率 $\varepsilon_{xx}=\partial u/\partial x$，剪切率 $\varepsilon_{xy}=\tfrac12(\partial u/\partial y+\partial v/\partial x)$。体积膨胀率（速度散度）
$$ \nabla\cdot\mathbf u=\lim_{\delta V\to0}\frac1{\delta V}\frac{D(\delta V)}{Dt}. $$
推导：取微元长方体 $\delta V=\delta x\delta y\delta z$，$\frac{1}{\delta V}\frac{D\delta V}{Dt}=\frac{1}{dx}\frac{Ddx}{Dt}+\cdots$。而 $\frac{D(dx)}{Dt}=du=(\partial u/\partial x)dx$（相对伸长），故
$$ \frac1{\delta V}\frac{D\delta V}{Dt}=\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}+\frac{\partial w}{\partial z}=\nabla\cdot\mathbf u. \quad\square $$
意义：不可压缩条件即 $\nabla\cdot\mathbf u=0$。

**6. 涡量与有旋性**

定义：涡量 $\boldsymbol\omega=\nabla\times\mathbf u$。有旋流 $\boldsymbol\omega\neq0$，无旋流 $\boldsymbol\omega=0$。旋度与角速度关系 $\boldsymbol\Omega=\boldsymbol\omega/2$。推导：由分解定理直接得出。意义：涡量是流动旋转强度的度量，涡线为 $\boldsymbol\omega$ 的场线。

**7. 速度势与流函数**

定义：无旋流 $\nabla\times\mathbf u=0$ 时存在速度势 $\phi$，$\mathbf u=\nabla\phi$。不可压缩二维时存在流函数 $\psi$，$u=\partial\psi/\partial y,\ v=-\partial\psi/\partial x$。

推导：无旋 $\nabla\times\nabla\phi\equiv0$，故构造的 $\mathbf u=\nabla\phi$ 自动无旋。二维不可压连续性 $\partial u/\partial x+\partial v/\partial y=0$，令 $u=\psi_y,\ v=-\psi_x$，则 $u_x+v_y=\psi_{yx}-\psi_{xy}=0$。$\square$ 意义：$\psi=$ 常数为流线；无旋不可压时二者满足柯西-黎曼关系，$\phi,\psi$ 为共轭调和函数，$\nabla^2\phi=0,\ \nabla^2\psi=0$。

**8. 连续微分方程（拉格朗日推导）**

质量守恒：
$$ \frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0, $$
不可压缩 $\nabla\cdot\mathbf u=0$。

推导（拉格朗日/物质体积法）：物质体积 $\delta V$ 质量守恒，$\frac{D(\rho\delta V)}{Dt}=0$：
$$ \delta V\frac{D\rho}{Dt}+\rho\frac{D(\delta V)}{Dt}=0. $$
由体积膨胀率 $\frac1{\delta V}\frac{D\delta V}{Dt}=\nabla\cdot\mathbf u$，代入
$$ \frac{D\rho}{Dt}+\rho\nabla\cdot\mathbf u=0. $$
展开 $\frac{D\rho}{Dt}=\partial_t\rho+\mathbf u\cdot\nabla\rho$，并注意 $\nabla\cdot(\rho\mathbf u)=\rho\nabla\cdot\mathbf u+\mathbf u\cdot\nabla\rho$，即得守恒形式。$\square$ 应用：所有守恒方程质量侧。

**9. 柱坐标与球坐标运动学公式**

柱坐标 $(r,\theta,z)$：
$$ \nabla\cdot\mathbf u=\frac1r\frac{\partial(ru_r)}{\partial r}+\frac1r\frac{\partial u_\theta}{\partial\theta}+\frac{\partial u_z}{\partial z}, $$
$$ (\mathbf u\cdot\nabla)\mathbf u\ \text{的 }r\text{分量}=\frac{\partial u_r}{\partial t}+u_r\frac{\partial u_r}{\partial r}+\frac{u_\theta}{r}\frac{\partial u_r}{\partial\theta}+u_z\frac{\partial u_r}{\partial z}-\frac{u_\theta^2}{r}. $$
推导：由坐标变换 $\mathbf u=u_r\mathbf e_r+u_\theta\mathbf e_\theta+u_z\mathbf e_z$，注意 $\partial\mathbf e_r/\partial\theta=\mathbf e_\theta$ 等，代入物质导数导数单位矢量变化项即得离心项 $-u_\theta^2/r$。应用：圆管、涡、喷流。

**10. 雷诺输运定理（运动学形式）**

设 $\Phi=\int_V \rho\phi\,dV$ 为广延量，则
$$ \frac{D}{Dt}\int_{V(t)}\rho\phi\,dV=\int_V\frac{\partial(\rho\phi)}{\partial t}dV+\oint_S \rho\phi\,\mathbf u\cdot\mathbf n\,dA. $$
推导：先对固定控制体与运动物质体积之差积分，或用雅可比行列式 $J=D(\mathbf x)/D(\mathbf a)$：$\frac{DJ}{Dt}=J\nabla\cdot\mathbf u$，则
$$ \frac{D}{Dt}\int_{V_0}\rho\phi J\,dV_0=\int_{V_0}\left[\frac{D(\rho\phi)}{Dt}+\rho\phi\nabla\cdot\mathbf u\right]J\,dV_0. $$
再由连续性 $\frac{D(\rho\phi)}{Dt}+\rho\phi\nabla\cdot\mathbf u=\frac{\partial(\rho\phi)}{\partial t}+\nabla\cdot(\rho\phi\mathbf u)$，用高斯定理得通量形式。$\square$ 应用：守恒方程积分形式的统一工具。

### 第4章 积分形式守恒方程与伯努利方程

**1. 质量守恒积分形式**

$$ \frac{d}{dt}\int_V\rho\,dV+\oint_S\rho\mathbf u\cdot\mathbf n\,dA=0. $$
定常不可压：出入流质量流量相等 $\sum_{\rm out}Q=\sum_{\rm in}Q$。推导：在雷诺输运中取 $\phi=1$ 得。应用：管流、激波管。

**2. 动量定理积分形式**

$$ \sum\mathbf F=\frac{d}{dt}\int_V\rho\mathbf u\,dV+\oint_S\rho\mathbf u(\mathbf u\cdot\mathbf n)\,dA. $$
推导：由雷诺输运 $\phi=\mathbf u$，$\frac{D}{Dt}\int\rho\mathbf u\,dV=\sum\mathbf F$（牛顿第二定律对物质体），展开得。定常时
$$ \sum\mathbf F=\oint_S\rho\mathbf u(\mathbf u\cdot\mathbf n)\,dA. $$
应用：弯管受力、射流冲击、推力。

**3. 动量矩方程**

$$ \sum\mathbf M_O=\frac{d}{dt}\int_V\rho(\mathbf r\times\mathbf u)dV+\oint_S\rho(\mathbf r\times\mathbf u)(\mathbf u\cdot\mathbf n)dA. $$
推导：对物质体角动量求导，雷诺输运取 $\phi=\mathbf r\times\mathbf u$，左侧为外力矩。应用：涡轮机欧拉方程、水轮机。

**4. 能量方程积分形式与伯努利方程**

对控制体（含轴功与传热）：
$$ \dot Q-\dot W_s=\frac{d}{dt}\int_V\rho e\,dV+\oint_S\left(p+\rho e\right)\mathbf u\cdot\mathbf n\,dA, $$
其中 $e=u_{\rm int}+\tfrac12|\mathbf u|^2+gz$。

沿流线、定常、无粘、不可压、无轴功与传热的伯努利方程：
$$ \frac{p}{\rho}+\frac{|\mathbf u|^2}{2}+gz=\text{const}. $$
推导：定常欧拉方程沿流线 $\mathbf u\cdot\nabla\mathbf u=-\nabla p/\rho-\nabla\Phi$，其中 $\mathbf u\cdot\nabla\mathbf u=\nabla(|\mathbf u|^2/2)-\mathbf u\times\boldsymbol\omega$。沿流线 $\mathbf u\cdot(\mathbf u\times\boldsymbol\omega)=0$，故
$$ \mathbf u\cdot\nabla\left(\frac{|\mathbf u|^2}{2}\right)=-\mathbf u\cdot\nabla\left(\frac p\rho+\Phi\right), $$
对不可压 $\rho$ 常数、$\Phi=gz$ 积分得 $\frac{|\mathbf u|^2}{2}+\frac p\rho+gz=\text{const}$。$\square$ 意义：机械能守恒（有旋流亦沿流线成立）。可压缩等熵推广：$\frac{|\mathbf u|^2}{2}+\frac{\gamma}{\gamma-1}\frac p\rho=\text{const}$。

**5. 伯努利方程的应用**

皮托管（测总压）：$p_0=p+\tfrac12\rho u^2\Rightarrow u=\sqrt{2(p_0-p)/\rho}$。文丘里管：
$$ Q=A_1\sqrt{\frac{2(p_1-p_2)}{\rho\!\left[1-(A_1/A_2)^2\right]}}. $$
推导：对文丘里管进、喉列伯努利与连续方程 $A_1u_1=A_2u_2$，消去 $u_2$ 解得。应用：流量计、空速管。

**6. 射流冲击与反作用力**

平板（法向）受力 $F=\rho Q u$；斜板 $F_n=\rho Q u\sin\alpha$；弯管动量差。推导：定常动量定理，入口动量 $\rho Qu$，出口 $\mathbf u_2$，$F=\rho Q(u_1-u_2)$ 向量式，投影得。应用：水轮机叶片、喷水推进、火箭推力
$$ F=\dot m_{\rm e}u_{\rm e}+(p_{\rm e}-p_{\rm a})A_{\rm e}. $$

**7. 非惯性系动量方程**

$$ \rho\frac{D\mathbf u_{\rm rel}}{Dt}=-\nabla p+\nabla\cdot\boldsymbol\tau+\rho\mathbf g-\rho\mathbf a_0-\rho(2\boldsymbol\Omega\times\mathbf u_{\rm rel})-\rho[\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf r)], $$
含科氏力与离心力。推导：由绝对加速度与相对加速度关系 $\mathbf a_{\rm abs}=\mathbf a_{\rm rel}+\mathbf a_0+2\boldsymbol\Omega\times\mathbf u_{\rm rel}+\boldsymbol\Omega\times(\boldsymbol\Omega\times\mathbf r)+\dot{\boldsymbol\Omega}\times\mathbf r$，代入牛顿第二定律。应用：旋转机械、地球物理流（地转风）。

**8. 火箭方程（齐奥尔科夫斯基）**

无外力、定常喷流：
$$ m\frac{d\mathbf v}{dt}=-\dot m u_{\rm e}\Rightarrow \Delta v=u_{\rm e}\ln\frac{m_0}{m_f}. $$
推导：对系统（箭体+推进剂）用动量守恒，$d(mv)=0$，箭体动量变化 $m\,dv$，喷出质量 $-dm$ 带走速度 $v-u_{\rm e}$，得 $m\,dv+(-dm)(v-u_{\rm e})-mv=0\Rightarrow m\,dv=u_{\rm e}\,dm$，积分且 $m=dm/dt<0$ 得 $\Delta v=u_{\rm e}\ln(m_0/m_f)$。$\square$ 应用：航天推进、姿态控制。

### 第5章 微分形式守恒方程与N-S方程

**1. 连续性方程（守恒与对流形式）**

$$ \frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0,\qquad \frac{D\rho}{Dt}+\rho\nabla\cdot\mathbf u=0. $$
推导见第3章第8条。不可压缩 $\nabla\cdot\mathbf u=0$；均质定常流亦然。应用：所有CFD求解。

**2. 欧拉方程（无粘动量方程）**

$$ \frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=-\frac1\rho\nabla p+\mathbf f. $$
推导：对无粘流体微元用牛顿第二定律，表面力只有压强，$\sum\mathbf F=-\nabla p\,dV+\rho\mathbf f\,dV$，单位质量即得。应用：气动、势流、激波（积分形式）。

**3. 应力张量与柯西动量方程**

定义应力张量 $\boldsymbol\sigma$，$\mathbf t^{(n)}=\boldsymbol\sigma\cdot\mathbf n$。动量方程
$$ \rho\frac{D\mathbf u}{Dt}=\nabla\cdot\boldsymbol\sigma+\rho\mathbf f. $$
推导：微元六面体各面力为 $\boldsymbol\sigma$ 的面积分，用高斯定理转体积分：
$$ \oint_S\boldsymbol\sigma\cdot\mathbf n\,dA=\int_V\nabla\cdot\boldsymbol\sigma\,dV. $$
与惯性项 $\int\rho\mathbf a\,dV$ 平衡，被积函数相等即得。$\square$ 意义：普遍形式，尚未用本构。应用：粘性、弹性介质。

**4. 本构方程与斯托克斯假设**

牛顿流体本构：
$$ \boldsymbol\sigma=-p\mathbf I+\boldsymbol\tau,\qquad \boldsymbol\tau=\mu(\nabla\mathbf u+\nabla\mathbf u^{\rm T})+\lambda(\nabla\cdot\mathbf u)\mathbf I. $$
斯托克斯假设 $\lambda=-\tfrac23\mu$，使体积粘性系数 $\kappa=\lambda+\tfrac23\mu=0$。

推导（斯托克斯假设）：设应力球量与体积应变率线性：$\tfrac13{\rm tr}(\boldsymbol\sigma)=-p+\kappa\nabla\cdot\mathbf u$，其中 $\kappa=\lambda+\tfrac23\mu$。要求平衡压强 $p$ 与 $\nabla\cdot\mathbf u$ 无关，则 $\kappa=0$，即 $\lambda=-\tfrac23\mu$。$\square$ 意义：牛顿流体在静止时退化为各向同性压强。

**5. Navier-Stokes方程推导**

把本构代入柯西方程：
$$ \rho\frac{D\mathbf u}{Dt}=-\nabla p+\nabla\cdot\left[\mu(\nabla\mathbf u+\nabla\mathbf u^{\rm T})+\lambda(\nabla\cdot\mathbf u)\mathbf I\right]+\rho\mathbf f. $$
对常粘度不可压：
$$ \nabla\cdot(\nabla\mathbf u+\nabla\mathbf u^{\rm T})=\nabla^2\mathbf u+\nabla(\nabla\cdot\mathbf u)=\nabla^2\mathbf u, $$
于是
$$ \boxed{\ \frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u=-\frac1\rho\nabla p+\nu\nabla^2\mathbf u+\mathbf f\ },\quad \nabla\cdot\mathbf u=0. $$
推导细节：$\nabla\cdot(\nabla\mathbf u)=\nabla^2\mathbf u$（矢量的拉普拉斯），$\nabla\cdot(\nabla\mathbf u^{\rm T})=\nabla(\nabla\cdot\mathbf u)=0$。$\square$ 意义：非线性对流项与粘性扩散项的竞争决定层流/湍流。极限：$\nu=0$ 欧拉；$Re\gg1$ 边界层；$Re\ll1$ 斯托克斯流。

**6. 能量方程微分形式**

内能形式：
$$ \rho\frac{D u_{\rm int}}{Dt}=-p\nabla\cdot\mathbf u+\Phi+\nabla\cdot(k\nabla T), $$
其中耗散函数
$$ \Phi=\mu\left[\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\right)-\frac23(\nabla\cdot\mathbf u)\delta_{ij}\right]\frac{\partial u_i}{\partial x_j}\ \ge 0. $$
总焓形式（无体积力、无源）：
$$ \rho\frac{D}{Dt}\!\left(h+\frac{|\mathbf u|^2}{2}\right)=\frac{\partial p}{\partial t}+\nabla\cdot(\boldsymbol\tau\cdot\mathbf u)+\nabla\cdot(k\nabla T). $$
推导：对物质能量守恒，机械能方程与总能方程相减得内能方程；机械能由动量方程点乘 $\mathbf u$ 得：
$$ \rho\mathbf u\cdot\frac{D\mathbf u}{Dt}=\rho\frac{D}{Dt}\!\left(\frac{|\mathbf u|^2}{2}\right)=\mathbf u\cdot(\nabla\cdot\boldsymbol\sigma)+\rho\mathbf u\cdot\mathbf f, $$
与总能方程相减消去机械能，得内能方程及耗散 $\Phi=\boldsymbol\tau:\nabla\mathbf u$。$\square$ 应用：气动加热、燃烧、传热。

**7. 状态方程与热力学关系**

$$ p=\rho R_gT,\qquad h=c_pT,\qquad u_{\rm int}=c_vT,\qquad \gamma=\frac{c_p}{c_v},\qquad R_g=c_p-c_v. $$
推导（理想气体）：迈耶公式由焓定义 $h=u_{\rm int}+p/\rho=u_{\rm int}+R_gT$，对 $T$ 求导 $c_p=c_v+R_g$。等熵关系由 $ds=c_v\,dT/T-R_g\,dp/p=0$ 积分得 $p/\rho^\gamma=\text{const}$。应用：可压缩流、声速。

**8. 初始条件与边界条件**

物理边界：壁面无滑移 $\mathbf u=\mathbf u_{\rm wall}$，不可穿透 $\mathbf u\cdot\mathbf n=\mathbf u_{\rm wall}\cdot\mathbf n$；无温跃变 $T=T_{\rm wall}$ 或热流条件。自由面运动学 $\mathbf u\cdot\nabla F+F_t=0$；动力条件 $[\![p]\!]=2\sigma H$。远场、进出口、周期性条件。推导：无滑移源于分子与壁面动量平衡。应用：CFD 定解。

**9. N-S 无量纲化与雷诺数**

取 $L,U,\rho U^2$ 为非定常时间、长度、压强与粘性应力尺度：
$$ \frac{\partial\mathbf u^*}{\partial t^*}+(\mathbf u^*\cdot\nabla^*)\mathbf u^*=-\nabla^*p^*+\frac1{Re}\nabla^{*2}\mathbf u^*+\frac1{Fr^2}\mathbf f^*, $$
$$ Re=\frac{\rho UL}{\mu}=\frac{UL}{\nu},\qquad Fr=\frac{U}{\sqrt{gL}}. $$
推导：令 $\mathbf u=U\mathbf u^*, \mathbf x=L\mathbf x^*, t=(L/U)t^*, p=\rho U^2p^*$，代入 N-S，对流项系数为 1，粘性项系数 $1/Re$。$\square$ 意义：$Re$ 度量惯性/粘性；$Re\ll1$ 斯托克斯流，$Re\gg1$ 边界层/湍流。应用：相似、模型试验。

### 第6章 涡动力学

**1. 涡线、涡管与涡通量**

定义：涡线为处处切于 $\boldsymbol\omega$ 的曲线；涡管为通过闭合曲线的一族涡线组成的管。涡通量 $J=\int_A\boldsymbol\omega\cdot\mathbf n\,dA$。由 $\nabla\cdot\boldsymbol\omega=\nabla\cdot(\nabla\times\mathbf u)=0$，涡场无源。推导：$\nabla\cdot(\nabla\times\mathbf u)=\partial_i\epsilon_{ijk}\partial_ju_k=0$（指标对称反对称缩并）。应用：涡管强度守恒。

**2. 斯托克斯定理与环量**

环量 $\Gamma=\oint_C\mathbf u\cdot d\mathbf l=\int_A\boldsymbol\omega\cdot\mathbf n\,dA$。推导：由斯托克斯定理 $\oint_C\mathbf u\cdot d\mathbf l=\int_A(\nabla\times\mathbf u)\cdot\mathbf n\,dA$。应用：测环量、升力。

**3. 涡管强度守恒（海姆霍兹第一定理）**

同一涡管内各横截面涡通量相等：
$$ \int_{A_1}\boldsymbol\omega\cdot\mathbf n_1\,dA_1=\int_{A_2}\boldsymbol\omega\cdot\mathbf n_2\,dA_2. $$
推导：对涡管侧面与两截面围成的闭曲面用 $\oint\boldsymbol\omega\cdot\mathbf n\,dA=\int_V\nabla\cdot\boldsymbol\omega\,dV=0$，侧面 $\boldsymbol\omega\perp\mathbf n$ 无通量，故两截面通量相等。$\square$ 意义：涡管不能在流体内部终止，只能成环或延伸至边界/自由面。

**4. 涡量输运方程**

对不可压牛顿流体：
$$ \frac{D\boldsymbol\omega}{Dt}=(\boldsymbol\omega\cdot\nabla)\mathbf u+\nu\nabla^2\boldsymbol\omega+\frac1{\rho^2}\nabla\rho\times\nabla p+\nabla\times\mathbf f. $$
推导：对动量方程取旋度。$\nabla\times[(\mathbf u\cdot\nabla)\mathbf u]=\nabla\times(\boldsymbol\Omega\times\mathbf u)$，利用矢量恒等式 $(\mathbf u\cdot\nabla)\mathbf u=\nabla(|\mathbf u|^2/2)-\mathbf u\times\boldsymbol\omega$，其旋度为 $\nabla\times(\boldsymbol\omega\times\mathbf u)$。再用恒等式
$$ \nabla\times(\boldsymbol\omega\times\mathbf u)=(\mathbf u\cdot\nabla)\boldsymbol\omega-(\boldsymbol\omega\cdot\nabla)\mathbf u+\boldsymbol\omega(\nabla\cdot\mathbf u)-\mathbf u(\nabla\cdot\boldsymbol\omega). $$
不可压、$\nabla\cdot\boldsymbol\omega=0$ 时化为 $(\mathbf u\cdot\nabla)\boldsymbol\omega-(\boldsymbol\omega\cdot\nabla)\mathbf u$。粘性项 $\nabla\times(\nu\nabla^2\mathbf u)=\nu\nabla^2\boldsymbol\omega$（常粘度）。斜压项来自 $-\nabla p\times\nabla(1/\rho)$。$\square$ 意义：涡拉伸项 $(\boldsymbol\omega\cdot\nabla)\mathbf u$ 使湍流涡增强（三维），是湍流能量串级机制；斜压项在斜压流体中生成涡（如海陆风）。

**5. 开尔文环量定理**

无粘、正压（$\rho=\rho(p)$）、有势体积力下，沿任一物质闭合曲线的环量守恒：
$$ \frac{D\Gamma}{Dt}=0. $$
推导：$\Gamma=\oint_C\mathbf u\cdot d\mathbf l$，
$$ \frac{D\Gamma}{Dt}=\oint_C\frac{D\mathbf u}{Dt}\cdot d\mathbf l+\oint_C\mathbf u\cdot\frac{D(d\mathbf l)}{Dt}=\oint_C\frac{D\mathbf u}{Dt}\cdot d\mathbf l+\oint_C\mathbf u\cdot d\mathbf u. $$
第二项 $\oint d(|\mathbf u|^2/2)=0$。代入欧拉方程 $\frac{D\mathbf u}{Dt}=-\nabla p/\rho+\mathbf f$，正压时 $-\nabla p/\rho=-\nabla\!\int dp/\rho$，有势力 $\mathbf f=-\nabla\Phi$，均为全微分，闭合积分零。故 $D\Gamma/Dt=0$。$\square$ 意义：无粘正压流中涡不能生成/消灭；解释机翼启动涡与附着环量。特例：无粘、初始无旋则永远无旋。

**6. 比奥-萨伐尔定理**

由涡量场求诱导速度：
$$ \mathbf u(\mathbf x)=-\frac1{4\pi}\int_V\frac{(\mathbf x-\mathbf x')\times\boldsymbol\omega(\mathbf x')}{|\mathbf x-\mathbf x'|^3}dV'+\mathbf u_\infty. $$
对涡线 $\boldsymbol\omega\,dV=\Gamma\,d\mathbf l$：
$$ \mathbf u=\frac{\Gamma}{4\pi}\int_C\frac{d\mathbf l\times(\mathbf x-\mathbf x')}{|\mathbf x-\mathbf x'|^3}. $$
推导：由 $\boldsymbol\omega=\nabla\times\mathbf u$、$\nabla\cdot\mathbf u=0$，引入矢量势 $\mathbf u=\nabla\times\mathbf A$，得泊松方程 $\nabla^2\mathbf A=-\boldsymbol\omega$，解 $\mathbf A=\frac1{4\pi}\int\frac{\boldsymbol\omega}{|\mathbf x-\mathbf x'|}dV'$，取旋度并用 $\nabla\times(\boldsymbol\omega/r)=\nabla(1/r)\times\boldsymbol\omega$ 得比奥-萨伐尔式。$\square$ 应用：涡方法、下洗、诱导阻力。

**7. 直线涡与兰肯涡**

无限长直涡诱导环向速度：
$$ u_\theta=\frac{\Gamma}{2\pi r}. $$
推导：取半径 $r$ 的圆，$\Gamma=\oint\mathbf u\cdot d\mathbf l=2\pi r u_\theta\Rightarrow u_\theta=\Gamma/(2\pi r)$。兰肯涡：核内刚体旋转 $u_\theta=\omega r$（$r\le a$），外部 $u_\theta=\Gamma/(2\pi r)$（$r>a$），$a$ 为核半径。核内涡量 $\omega=2\omega_0$，外部无旋。压强由径向平衡 $\frac{dp}{dr}=\rho u_\theta^2/r$ 积分：外部 $p=p_\infty-\rho\Gamma^2/(8\pi^2r^2)$，核心中心低压 $\Delta p=\rho\omega_0^2a^2/2$。应用：龙卷风、翼尖涡、旋流器。

**8. 涡对、涡环与自诱导**

两个等强度反号点涡（涡对）以 $U=\Gamma/(2\pi d)$ 平移（$d$ 为间距）。同号涡绕公共中心旋转。涡环以
$$ U_{\rm ring}\approx\frac{\Gamma}{4\pi R}\left(\ln\frac{8R}{a}-\frac14\right) $$
自诱导平移。推导：涡对的每个涡由另一个涡诱导速度 $U=\Gamma/(2\pi d)$，方向相同故整体平移；同号涡诱导方向垂直连线故旋转。涡环用比奥-萨伐尔积分得含曲率的对数修正。应用：飞机尾涡、烟圈、真菌孢子传播。

**9. 涡的生成与破裂**

涡在壁面由无滑移与压强梯度生成；涡拉伸使涡管变细、涡量增强（角动量守恒）。涡破裂在强旋流中出现滞止泡或螺旋破裂，由瑞利准则 $V_\theta/r$ 或斯沃茨准则判定。推导（涡拉伸）：涡管体积守恒，轴向拉伸使截面缩小，由涡管强度守恒 $\omega A=$ 常数，故 $\omega$ 增大。应用：旋流燃烧、前体涡、龙卷风。

### 第7章 势流理论

**1. 速度势与拉普拉斯方程**

不可压无旋：$\mathbf u=\nabla\phi$，连续性给出
$$ \nabla^2\phi=0. $$
推导：$\nabla\cdot\mathbf u=\nabla\cdot\nabla\phi=\nabla^2\phi=0$。$\square$ 应用：绕流、地下水、电渗类比。

**2. 基本解与叠加**

均匀流：$\phi=Ux,\ \psi=Uy$。点源（强度 $q=2\pi m$）：$\phi=\frac{m}{2\pi}\ln r,\ u_r=\frac{m}{r}$。点涡：$\phi=\frac{\Gamma}{2\pi}\theta,\ u_\theta=\frac{\Gamma}{2\pi r}$。偶极子（强度 $\mu$）：$\phi=-\frac{\mu}{2\pi}\frac{x}{x^2+y^2}$。

推导：对 $\nabla^2\phi=0$ 求轴对称解 $\phi=\phi(r)$：$\frac1r\frac{d}{dr}(r\phi')=0\Rightarrow\phi=A\ln r+B$，得源汇；点涡由环量归一。叠加原理源于方程线性。应用：绕流构造。

**3. 圆柱绕流与表面压力**

叠加均匀流与偶极子（$\mu=2\pi Ua^2$）：
$$ \phi=U\!\left(r+\frac{a^2}{r}\right)\cos\theta,\qquad \psi=U\!\left(r-\frac{a^2}{r}\right)\sin\theta. $$
速度
$$ u_r=U\!\left(1-\frac{a^2}{r^2}\right)\cos\theta,\quad u_\theta=-U\!\left(1+\frac{a^2}{r^2}\right)\sin\theta. $$
表面 $r=a$：$u_r=0,\ u_\theta=-2U\sin\theta$。压强由伯努利
$$ C_p=\frac{p-p_\infty}{\tfrac12\rho U^2}=1-4\sin^2\theta. $$
推导：均匀流 $\phi_1=Ux=Ur\cos\theta$ 与偶极子 $\phi_2=-\frac{\mu}{2\pi}\frac{\cos\theta}{r}$ 叠加，令 $r=a$ 处 $u_r=0$ 定 $\mu=2\pi Ua^2$。$\square$ 意义：前后对称无阻力（达朗贝尔佯谬）。极限：$\theta=90^\circ$ 处 $C_p=-3$（最大吸力）。

**4. 带环量圆柱与马格努斯效应**

表面 $u_\theta=-2U\sin\theta+\Gamma/(2\pi a)$，升力（单位长度）
$$ L'=\rho U\Gamma. $$
推导：压强 $p=p_\infty+\tfrac12\rho U^2-\tfrac12\rho u_\theta^2$，对 $\theta$ 积分 $F_y=-\int_0^{2\pi}p\sin\theta\,a\,d\theta$。$u_\theta^2$ 中的 $(-2U\sin\theta+\Gamma/2\pi a)^2$ 交叉项 $-\tfrac12\rho\cdot2(-2U\sin\theta)(\Gamma/2\pi a)$ 贡献非零积分，其余对称项积分为零，得 $F_y=\rho U\Gamma$。$\square$ 应用：马格努斯效应、旋转球、菲涅尔旋筒。

**5. 镜像法与地面效应**

平壁附近源/涡用镜像源/涡处理，边界上法向速度为零。半无限流中涡距壁 $h$，镜像涡在 $-h$ 反号，使壁面法向速度为零；实际壁面对涡有诱导速度。应用：机翼地面效应、涡-壁相互作用。

**6. 保角变换与儒可夫斯基翼型**

解析变换 $w=f(z)$ 保角且保持 $\phi,\psi$。儒可夫斯基变换
$$ \zeta=z+\frac{c^2}{z}, $$
把圆心偏离原点的圆映射为翼型。绕流在 $\zeta$ 平面升力
$$ L=\rho U_\infty\Gamma,\qquad \Gamma=4\pi U_\infty c\sin(\alpha+\beta). $$
推导：变换保持环量与源汇强度；由库塔条件取后缘有限速度（后缘为变换的尖点），确定环量。$\square$ 应用：翼型设计、升力线理论。

**7. 库塔-儒可夫斯基定理**

$$ L'=\rho U_\infty\Gamma,\qquad M'=\rho U_\infty\Gamma x_{\rm cp}. $$
推导：由上述带环量圆柱积分，或由动量与环量关系。物理意义：升力源于环量与来流耦合。应用：机翼升力、启动涡。

**8. 达朗贝尔佯谬与粘性作用**

无粘不可压势流中，任意封闭物体总阻力为零。推导：绕物体压力分布前后对称，或由动量积分与涡通量。实际物体有粘性产生边界层与尾迹、分离，阻力不为零。意义：揭示理想流体模型的局限；处理需边界层/湍流。

**9. 非定常势流与附加质量**

物体在流体中加速时，周围流体动能变化等效于附加质量 $m_a$。球体 $m_a=\tfrac12\rho V_d$。推导：势流总动能
$$ T=\frac12\int\rho|\nabla\phi|^2dV=\frac12 m_a U^2, $$
对球体用偶极子势计算积分得 $m_a=\tfrac12\rho\cdot\tfrac43\pi a^3$。应用：船舶、颗粒加速、水下运动。

### 第8章 量纲分析与相似理论

**1. 量纲与量纲齐次性**

基本量纲：$M,L,T,\Theta$（质量、长度、时间、温度）。任何物理量 $q$ 的量纲 $[q]=M^aL^bT^c\Theta^d$。物理方程量纲齐次：各项量纲相同。推导：方程两边的量纲必须一致（否则单位制改变会破坏方程），这约束经验公式形式。

**2. 白金汉 π 定理**

若问题涉及 $n$ 个物理量、$k$ 个基本量纲，则可组成 $n-k$ 个独立无量纲数 $\pi_1,\dots,\pi_{n-k}$，且物理关系等价于 $f(\pi_1,\dots,\pi_{n-k})=0$。

证明：设量 $q_1,\dots,q_n$，量纲矩阵 $A$（$k\times n$），秩 $r=k$。齐次线性方程组 $A\mathbf x=0$ 的解空间维数 $n-k$。对每对解向量 $\mathbf x^{(j)}$ 构造乘积 $\pi_j=\prod q_i^{x_i^{(j)}}$，即无量纲。任一无量纲量皆可表为 $\pi_j$ 的函数；任一物理关系 $g(q_1,\dots,q_n)=0$ 可写成 $F(\pi_1,\dots,\pi_{n-k})=0$。$\square$ 应用：实验设计、关联数据、相似准则。

**3. 相似准则与物理意义**

$$ Re=\frac{UL}{\nu},\quad Fr=\frac{U}{\sqrt{gL}},\quad Eu=\frac{p}{\rho U^2},\quad Ma=\frac{U}{c},\quad We=\frac{\rho U^2L}{\sigma},\quad St=\frac{fL}{U}. $$
推导：由无量纲 N-S 得 $Re,Fr$；由自由面与表面张力条件得 $We$；马赫数由可压缩 Euler 方程得。应用：风洞、水槽、血流动、喷雾。

**4. N-S 方程无量纲化与相似条件**

见第5章第9条。完全相似需所有无量纲数相等，通常不可能同时满足，故按主导力选择（部分相似）。应用：缩比模型、雷诺数效应修正。

**5. 常用经验关联**

管流摩擦系数 $f=f(Re,\varepsilon/D)$（穆迪图）；阻力系数 $C_D=C_D(Re,Ma,\text{形状})$；传热 $Nu=f(Re,Pr)$。这些由 π 定理组织，指数由实验/量纲匹配确定。应用：工程估算。

**6. 量纲分析实例：球体阻力**

设为 $D=f(\rho,U,d,\mu)$，取 $\rho,U,d$ 为基，得
$$ C_D=\frac{D}{\tfrac12\rho U^2 d^2}=f\!\left(\frac{\rho Ud}{\mu}\right)=f(Re). $$
推导：$\pi_1=D/(\rho U^2d^2)$，$\pi_2=\mu/(\rho Ud)=1/Re$。斯托克斯区 $C_D=24/Re$；牛顿区 $C_D\approx0.44$。应用：沉降、气溶胶。

### 第9章 边界层理论

**1. 边界层概念与普朗特量级分析**

定义：大雷诺数下贴近壁面、粘性主导的薄层，厚度 $\delta\ll L$。量级：$x\sim L,\ y\sim\delta,\ u\sim U$。由连续性 $v\sim U\delta/L$。推导：$x$ 方向 N-S 惯性 $\sim U^2/L$，粘性 $\sim\nu U/\delta^2$。平衡得
$$ \delta\sim\sqrt{\frac{\nu L}{U}}=\frac{L}{\sqrt{Re_L}}. \quad\square $$

**2. 边界层方程推导**

从 N-S 出发，量级比较：$\partial^2 u/\partial x^2\sim U/L^2\ll \partial^2u/\partial y^2\sim U/\delta^2$。$y$ 方向动量方程各项 $\sim\delta/L$，得 $\partial p/\partial y\approx0$，即
$$ p=p_e(x). $$
得到普朗特边界层方程
$$ \frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}=-\frac1\rho\frac{dp_e}{dx}+\nu\frac{\partial^2u}{\partial y^2},\qquad \frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}=0. $$
边界条件：$y=0:u=v=0$；$y\to\infty:u\to U_e(x)$。

**3. 边界层厚度定义**

位移厚度
$$ \delta^*=\int_0^\infty\left(1-\frac u{U_e}\right)dy, $$
动量厚度
$$ \theta=\int_0^\infty\frac u{U_e}\left(1-\frac u{U_e}\right)dy, $$
形状因子 $H=\delta^*/\theta$。推导（位移厚度）：因边界层内速度亏损，外部流线外移 $\delta^*$，由质量流量等价
$$ \int_0^\infty(U_e-u)dy=U_e\delta^*. \quad\square $$
意义：断面有效位移。湍流平板 $H\approx1.4$，层流 $\approx2.6$。

**4. 卡门动量积分方程**

$$ \frac{\tau_w}{\rho}=\frac{d}{dx}\!\left(U_e^2\theta\right)+U_e\frac{dU_e}{dx}\delta^*. $$
推导：对边界层方程沿 $y$ 积分。先积分连续性得 $v=-\int_0^y \partial_xu\,dy$，代入动量方程，利用 $y\to\infty$ 处 $u=U_e$ 及伯努利 $U_edU_e/dx=-\frac1\rho dp_e/dx$，整理得。$\square$ 应用：近似求解各种边界层。

**5. 平板层流边界层的布拉修斯解**

引入相似变量 $\eta=y\sqrt{U/\nu x}$，流函数 $\psi=\sqrt{\nu U x}f(\eta)$，则 $u=Uf'(\eta)$。布拉修斯方程
$$ 2f'''+ff''=0,\quad f(0)=f'(0)=0,\ f'(\infty)=1. $$
推导：由 $u=\psi_y=Uf'$，$v=\sqrt{\nu U/x}(\eta f'-f)$，逐项代入边界层方程化简得。数值：$f''(0)=0.3321$，
$$ \delta=5.0\sqrt{\frac{\nu x}{U}},\quad \tau_w=\frac{\rho U^2}{\sqrt{Re_x}}f''(0),\quad C_f=\frac{1.328}{\sqrt{Re_L}}. $$
应用：平板摩阻、传热类比。

**6. 平板湍流边界层**

$1/7$ 次幂速度剖面 $u/U_e=(y/\delta)^{1/7}$。卡门积分给出
$$ \delta\approx0.37x\,Re_x^{-1/5},\qquad \frac{\tau_w}{\rho}=0.0225U_e^2Re_\delta^{-1/4}. $$
平板总摩阻 $C_f=0.074Re_L^{-1/5}$（$Re<10^7$），普朗特-施利希廷 $C_f=0.455/(\ln Re_L)^{2.58}$。推导：把幂律剖面代入卡门方程积分，用经验壁面律 $\tau_w=0.0225\rho U_e^2(U_e\delta/\nu)^{-1/4}$ 求解 $\delta(x)$。应用：船舶、机翼摩阻估算。

**7. 边界层分离**

逆压梯度 $dp/dx>0$ 使近壁流体减速，壁面剪切 $\tau_w=0$ 处分离，出现回流。分离判据 $(\partial u/\partial y)_{y=0}=0$。推导：在 $y=0$ 处边界层方程给出 $\mu(\partial^2u/\partial y^2)_0=dp_e/dx$。逆压时壁面速度剖面为负曲率，当 $\tau_w=0$ 及壁面附近出现 $u<0$，即分离。应用：失速、钝体尾流、扩散器。

**8. 波尔豪森近似（积分法）**

假设多项式速度剖面
$$ \frac u{U_e}=2\eta-2\eta^3+\eta^4,\quad \eta=y/\delta, $$
代入卡门方程求 $\delta(x)$。对平板层流得 $\delta=5.5\sqrt{\nu x/U}$，$C_f=1.46/\sqrt{Re_L}$，与布拉修斯接近。推导：由多项式满足边界条件（壁面 $u=0,u_{yy}$ 由压力、外缘 $u_y=0,u_{yy}=0$）构造，计算 $\delta^*,\theta,\tau_w$ 代入积分方程，解一阶 ODE。应用：快速工程估算。

### 第10章 层流管流与湍流基础

**1. 雷诺实验与转捩**

层流规则，湍流随机脉动。圆管转捩临界 $Re_c=\rho UD/\mu\approx2300$。推导：量纲分析得唯一参数 $Re$；实验确定临界值。应用：管流设计、湍流控制。

**2. 哈根-泊肃叶流（层流圆管）**

充分发展管流速度抛物分布
$$ u(r)=\frac{1}{4\mu}\left(-\frac{dp}{dx}\right)(R^2-r^2),\qquad u_{\max}=\frac{R^2}{8\mu}\left(-\frac{dp}{dx}\right). $$
流量 $Q=\dfrac{\pi R^4}{8\mu}\left(-\dfrac{dp}{dx}\right)$，平均速度 $\bar u=Q/(\pi R^2)=u_{\max}/2$，摩阻系数 $f=64/Re$。

推导：轴向 N-S 充分发展时 $\partial u/\partial x=0$，惯性为零，$\frac1r\frac{d}{dr}(r\frac{du}{dr})=\frac1\mu\frac{dp}{dx}$。积分两次，用 $u(R)=0$ 与 $du/dr(0)=0$ 得。$\square$ 应用：粘性计量、微流控、润滑。

**3. 平板间库埃特-泊肃叶流**

固定下板、上板速度 $U$：
$$ u(y)=\frac{U y}{h}-\frac{1}{2\mu}\frac{dp}{dx}y(h-y). $$
推导：$\frac{d^2u}{dy^2}=\frac1\mu\frac{dp}{dx}$，两次积分，用 $u(0)=0,u(h)=U$ 定常数。应用：轴承、涂布。

**4. 湍流统计描述**

雷诺分解 $u_i=\bar u_i+u_i'$，$\overline{u_i'}=0$。均方根 $u_{\rm rms}=\sqrt{\overline{u'^2}}$，湍流强度 $Tu=u_{\rm rms}/\bar u$。相关函数
$$ R_{ij}(\tau)=\overline{u_i'(t)u_j'(t+\tau)},\qquad \Lambda=\int_0^\infty R_{11}(\tau)\,d\tau. $$
推导：由随机过程定义。应用：风工程、湍流测量。

**5. 雷诺平均 N-S 方程（RANS）**

$$ \frac{\partial\bar u_i}{\partial t}+\bar u_j\frac{\partial\bar u_i}{\partial x_j}=-\frac1\rho\frac{\partial\bar p}{\partial x_i}+\nu\frac{\partial^2\bar u_i}{\partial x_j\partial x_j}-\frac{\partial}{\partial x_j}\overline{u_i'u_j'}. $$
雷诺应力 $\tau_{ij}^{(t)}=-\rho\overline{u_i'u_j'}$。推导：对 N-S 取平均，利用 $\overline{\bar u}=0$ 与 $\overline{u_i'}=0$。非线性项的分解给出
$$ \overline{u_j\partial_j u_i}=\bar u_j\partial_j\bar u_i+\overline{u_j'\partial_j u_i'}=\bar u_j\partial_j\bar u_i+\partial_j\overline{u_i'u_j'}. $$
移项得额外应力散度。$\square$ 意义：湍流封闭问题（未知量多于方程）。

**6. 布辛涅斯克涡粘模型与混合长度**

$$ -\rho\overline{u_i'u_j'}=\mu_t\left(\frac{\partial\bar u_i}{\partial x_j}+\frac{\partial\bar u_j}{\partial x_i}\right)-\frac23\rho k\delta_{ij}, $$
$$ \mu_t=\rho l_m^2\left|\frac{d\bar u}{dy}\right|,\qquad k=\frac12\overline{u_i'u_i'}. $$
推导（混合长度）：普朗特假设流体团在混合长度 $l_m$ 内保持原动量，掺混后速度脉动 $u'\sim l_m\,d\bar u/dy$，故 $\overline{u'v'}\sim -l_m^2|d\bar u/dy|d\bar u/dy$。$\square$ 应用：壁面律、自由剪切流。

**7. 湍流能谱与柯尔莫哥洛夫理论**

能谱 $E(\kappa)$ 在惯性子区
$$ E(\kappa)=C_K\varepsilon^{2/3}\kappa^{-5/3}. $$
柯尔莫哥洛夫尺度 $\eta=(\nu^3/\varepsilon)^{1/4}$，含能尺度 $L\sim k^{3/2}/\varepsilon$，$Re_\lambda=(u'\lambda/\nu)$。推导：量纲分析，惯性区仅依赖 $\varepsilon$ 与 $\kappa$，$[E]=L^3T^{-2}$，$\varepsilon=[L^2T^{-3}]$，得 $E=C\varepsilon^{2/3}\kappa^{-5/3}$。$\square$ 应用：湍流建模、大气与工程湍流。

**8. 壁面律与摩阻系数**

粘性底层 $u^+=y^+$，对数律
$$ u^+=\frac1\kappa\ln y^++B,\quad \kappa\approx0.41,\ B\approx5.0, $$
$u^+=u/u_\tau,\ y^+=yu_\tau/\nu,\ u_\tau=\sqrt{\tau_w/\rho}$。圆管摩阻：层流 $f=64/Re$，湍流布拉修斯 $f=0.316Re^{-1/4}$，科尔布鲁克公式
$$ \frac1{\sqrt f}=-2\ln\!\left(\frac{\varepsilon/D}{3.7}+\frac{2.51}{Re\sqrt f}\right). $$
推导：壁面律由量纲分析＋实验（或混合长度积分）；科尔布鲁克式由对数律与粗糙度叠加。应用：管网、穆迪图。

**9. 局部损失与管网**

局部损失 $h_m=\zeta\dfrac{\bar u^2}{2g}$，总损失 $h_f=\lambda\dfrac{L}{D}\dfrac{\bar u^2}{2g}$；串联流量相同损失叠加，并联损失相同流量叠加。推导：由伯努利与量纲分析，局部损失与动压成正比。应用：管路计算。

### 第11章 可压缩一维流动

**1. 热力学关系与滞止量**

$$ h_0=h+\frac{u^2}{2},\quad T_0=T+\frac{u^2}{2c_p},\quad p_0=p\left(1+\frac{\gamma-1}{2}Ma^2\right)^{\frac{\gamma}{\gamma-1}}, $$
$$ \frac{T_0}{T}=1+\frac{\gamma-1}{2}Ma^2,\quad \frac{\rho_0}{\rho}=\left(1+\frac{\gamma-1}{2}Ma^2\right)^{\frac1{\gamma-1}}. $$
推导：由能量守恒 $h_0=$ 常数，$h=c_pT$，$Ma=u/c$，$c=\sqrt{\gamma R_gT}$，代入等熵关系 $p_0/p=(T_0/T)^{\gamma/(\gamma-1)}$。$\square$ 应用：气动测总压总温。

**2. 等熵关系**

$$ \frac p{\rho^\gamma}=\text{const},\quad \frac T{\rho^{\gamma-1}}=\text{const},\quad \frac p{T^{\gamma/(\gamma-1)}}=\text{const}. $$
推导：$ds=c_v\,dT/T-R_g\,dp/p=0$，用 $p=\rho R_gT$ 消元得。应用：喷管、扩压器。

**3. 一维等熵变截面流**

连续性 $\rho uA=$ 常数，动量 $dp+\rho u\,du=0$，得面积-速度关系
$$ \frac{dA}{A}=(Ma^2-1)\frac{du}{u}. $$
推导：对数微分连续性 $d\rho/\rho+du/u+dA/A=0$；等熵 $dp=\gamma p/\rho\,d\rho$；动量 $dp=-\rho u\,du$，故 $d\rho/\rho=-Ma^2\,du/u$。代入消去 $d\rho$ 得。$\square$ 意义：$Ma<1$ 收缩加速、扩张减速；$Ma>1$ 收缩减速、扩张加速；声速点 $Ma=1$ 处 $dA=0$（喉部）。应用：拉瓦尔喷管。

**4. 拉瓦尔喷管**

喉部 $Ma=1$，面积比
$$ \frac A{A^*}=\frac1{Ma}\left[\frac{2}{\gamma+1}\left(1+\frac{\gamma-1}{2}Ma^2\right)\right]^{\frac{\gamma+1}{2(\gamma-1)}}. $$
质量流量
$$ \dot m=\frac{p_0A^*}{\sqrt{T_0}}\sqrt{\frac\gamma{R_g}}\left(\frac2{\gamma+1}\right)^{\frac{\gamma+1}{2(\gamma-1)}}. $$
推导：在喉部 $Ma=1$，$T^*=T_0\frac2{\gamma+1}$，$p^*=p_0(\frac2{\gamma+1})^{\gamma/(\gamma-1)}$，$\rho^*=\rho_0(\frac2{\gamma+1})^{1/(\gamma-1)}$，代入 $\dot m=\rho^*u^*A^*$，$u^*=c^*=\sqrt{\gamma R_gT^*}$，化简得。面积比由等熵关系与 $u=Ma\,c$ 得。应用：火箭、喷气发动机、风洞。

**5. 正激波与兰金-雨贡纽关系**

连续、动量、能量：
$$ \rho_1u_1=\rho_2u_2,\quad p_1+\rho_1u_1^2=p_2+\rho_2u_2^2,\quad h_1+\tfrac12u_1^2=h_2+\tfrac12u_2^2. $$
得兰金-雨贡纽
$$ \frac{p_2}{p_1}=\frac{2\gamma Ma_1^2-(\gamma-1)}{\gamma+1},\qquad \frac{\rho_2}{\rho_1}=\frac{(\gamma+1)Ma_1^2}{2+(\gamma-1)Ma_1^2}, $$
$$ Ma_2^2=\frac{(\gamma-1)Ma_1^2+2}{2\gamma Ma_1^2-(\gamma-1)}. $$
推导：由动量与连续消去 $u$：$\rho_1u_1=\rho_2u_2=m$，$p_2-p_1=m^2(1/\rho_1-1/\rho_2)$，结合能量 $h_2-h_1=\frac12(u_1^2-u_2^2)$ 及 $h=\frac{\gamma}{\gamma-1}\frac p\rho$，解出压强比与密度比。

性质：$Ma_1>1$ 必 $Ma_2<1$；总压损失
$$ \frac{p_{02}}{p_{01}}=\left[\frac{(\gamma+1)Ma_1^2}{2+(\gamma-1)Ma_1^2}\right]^{\frac{\gamma}{\gamma-1}}\left[\frac{\gamma+1}{2\gamma Ma_1^2-(\gamma-1)}\right]^{\frac1{\gamma-1}}<1. $$
熵增 $s_2-s_1=c_p\ln(T_2/T_1)-R_g\ln(p_2/p_1)>0$，激波不可逆。应用：超音速进气道、激波风洞。

**6. 斜激波与激波极线**

斜激波角 $\beta$、偏转角 $\theta$，用相对法向分量 $Ma_{1n}=Ma_1\sin\beta$ 代入正激波关系，切向分量守恒。$\theta$-$\beta$-$Ma$ 关系
$$ \tan\theta=2\cot\beta\,\frac{Ma_1^2\sin^2\beta-1}{Ma_1^2(\gamma+\cos2\beta)+2}. $$
推导：把速度分解为激波法向/切向，法向用正激波关系，切向不变，几何合成得。$\square$ 应用：超音速机翼、进气道。

**7. 普朗特-迈耶膨胀波**

凸角膨胀等熵，偏转角
$$ \nu(Ma)=\sqrt{\frac{\gamma+1}{\gamma-1}}\arctan\sqrt{\frac{\gamma-1}{\gamma+1}(Ma^2-1)}-\arctan\sqrt{Ma^2-1}. $$
推导：对等熵动量方程 $d\theta=-\frac{\sqrt{Ma^2-1}}{1+\frac{\gamma-1}{2}Ma^2}\frac{dMa}{Ma}$ 积分，得普朗特-迈耶函数。应用：超音速膨胀、喷流。

**8. 范诺流（等截面摩擦管流）**

有摩擦无传热，总温不变、总压下降、熵增。微分关系
$$ \frac{dp}{p}=-\frac{\gamma Ma^2}{1-Ma^2}\frac{4f\,dx}{Dh}\left(\text{一种形式}\right), $$
$$ \frac{4fL^*}{Dh}=\frac{1-Ma^2}{\gamma Ma^2}+\frac{\gamma+1}{2\gamma}\ln\frac{(\gamma+1)Ma^2}{2+(\gamma-1)Ma^2}. $$
推导：连续性、动量（含壁面摩擦 $-\tau_w P\,dx$）、能量 $T_0=$ 常数联立，消元得 $Ma$ 沿程 ODE，积分得临界长度 $L^*$。应用：绝热管道、喷管摩擦修正。

**9. 瑞利流（等截面加热管流）**

无摩擦有传热。关系
$$ \frac{T_0}{T_0^*}=\frac{2(\gamma+1)Ma^2(1+\frac{\gamma-1}{2}Ma^2)}{(1+\gamma Ma^2)^2}, $$
$$ \frac{p}{p^*}=\frac{1+\gamma}{1+\gamma Ma^2},\qquad \frac{T}{T^*}=\left[\frac{(\gamma+1)Ma}{1+\gamma Ma^2}\right]^2. $$
推导：连续、动量（无摩擦）、能量（加热）与状态方程联立，积分得。亚音速加热加速趋近 $Ma=1$，超音速加热减速趋近 $Ma=1$（热壅塞）。应用：燃烧室、加力燃烧。

### 第12章 流体力学变分原理与稳定性

**1. 泛函与变分、欧拉-拉格朗日方程**

设泛函 $J[y]=\int_{x_1}^{x_2}F(x,y,y')\,dx$。极值条件
$$ \frac{\partial F}{\partial y}-\frac{d}{dx}\!\left(\frac{\partial F}{\partial y'}\right)=0. $$
推导：令 $y\to y+\varepsilon\eta$，$\delta J=\int\left(F_y\eta+F_{y'}\eta'\right)dx$。对第二项分部积分，
$$ \delta J=\int\left[F_y-\frac{d}{dx}F_{y'}\right]\eta\,dx+\left[F_{y'}\eta\right]_{x_1}^{x_2}. $$
端点固定 $\eta(x_1)=\eta(x_2)=0$，由 $\delta J=0$ 对任意 $\eta$ 得 EL 方程。$\square$ 应用：最小作用量、流体哈密顿表述。

**2. 哈密顿原理与流体作用量**

定义作用量
$$ S=\int_{t_1}^{t_2}L\,dt,\qquad L=T-V. $$
真实运动使 $\delta S=0$。连续体拉格朗日密度 $\mathcal L=\tfrac12\rho|\mathbf u|^2-\rho e_{\rm int}-\rho\Phi$，作用量
$$ S=\int\int_V\mathcal L\,dV\,dt. $$
推导：对拉格朗日描述 $\mathbf x(\mathbf a,t)$ 变分，得连续性、动量与能量方程（通过诺特/变分）。$\square$ 应用：守恒律、辛结构、数值格式。

**3. 拉格朗日与哈密顿描述**

拉格朗日量 $L(q,\dot q)$，动量 $p=\partial L/\partial\dot q$，哈密顿量
$$ H=p\cdot\dot q-L,\quad \dot q=\frac{\partial H}{\partial p},\quad \dot p=-\frac{\partial H}{\partial q}. $$
推导：对 $H$ 求微分并利用 EL 方程得正则方程。流体中用速度势/克莱布施表示，哈密顿结构给出哈密顿方程。应用：稳定性、守恒格式。

**4. 诺特定理与守恒律**

作用量在时间平移、空间平移、旋转下不变分别给出能量、动量、角动量守恒。推导：设 $L$ 不显含时间，则
$$ \frac{d}{dt}\!\left(\sum\dot q_i\frac{\partial L}{\partial\dot q_i}-L\right)=0\Rightarrow E=\text{const}. $$
空间平移不变给出 $\sum\partial L/\partial\dot q_i$ 守恒（动量）。旋转不变给出角动量。$\square$ 应用：守恒律来源、对称性分析。

**5. 流体平衡稳定性（瑞利-泰勒）**

密度大的流体位于轻流体之上不稳定，扰动增长率
$$ \sigma=\sqrt{Ag\,k-\frac{\sigma k^3}{\rho_1+\rho_2}}, $$
$A=(\rho_2-\rho_1)/(\rho_2+\rho_1)$ 为阿特伍德数。推导：对界面扰动 $\eta=\eta_0e^{ikx+\sigma t}$，用线性化欧拉方程与界面动力学/运动学条件，匹配速度势，得色散关系。短波由表面张力稳定。应用：惯性约束聚变、云、超新星。

**6. 开尔文-亥姆霍兹不稳定性**

两平行流有速度差时界面不稳定，增长率
$$ \sigma=k\frac{\rho_1U_1+\rho_2U_2}{\rho_1+\rho_2}\pm\cdots $$
无表面张力、无重力时对任意 $k$ 不稳定；表面张力稳定短波。推导：线性化势流，界面运动学与动力学条件给出色散关系，虚部即增长率。应用：风浪、云街、剪切层。

**7. 瑞利-贝纳德对流**

加热流体层从下方，临界瑞利数
$$ Ra_c=\frac{g\beta\Delta T d^3}{\nu\kappa}\approx1708, $$
对应临界波数 $k_c\approx3.117$。推导：布辛涅斯克近似线性稳定性分析，扰动方程
$$ (\partial_t-\nu\nabla^2)(\partial_t-\kappa\nabla^2)\nabla^2 w = Ra\,\kappa\nabla_H^2 w, $$
设正规模 $w\propto e^{ikx+\sigma t}\sin(n\pi z)$，临界条件得 $Ra_c=27\pi^4/4\approx1708$（自由-自由边界）。$\square$ 应用：对流胞、太阳对流、地幔。

**8. 平行流的稳定性（奥尔-索末菲）**

二维平行流扰动满足
$$ (U-c)(\phi''-k^2\phi)-U''\phi=-\frac{i}{kRe}(\phi''''-2k^2\phi''+k^4\phi). $$
推导：线性化 N-S，设流函数扰动 $\psi=\phi(y)e^{ik(x-ct)}$，舍非线性项得。应用：边界层转捩、稳定性。

**9. 旋转流稳定性（瑞利准则）**

无粘旋流稳定条件：角动量 $r^2\Omega(r)$ 向外增大不减小，即
$$ \frac{d}{dr}(r^2\Omega)^2\ge0\ \text{或}\ \frac{d}{dr}(r u_\theta)^2\ge0. $$
推导：对轴对称扰动做能量/变分分析，离心力与压力梯度平衡；违反则惯性不稳定。应用：泰勒-库埃特流、涡破碎。

### 第13章 界面现象与多相流基础

**1. 杨-拉普拉斯方程与曲率**

一般曲面：
$$ \Delta p=\sigma\left(\frac1{R_1}+\frac1{R_2}\right)=2\sigma H. $$
推导见第1章第7条。对柱面 $R_2=\infty$，$\Delta p=\sigma/R$。应用：毛细管液面、气泡。

**2. 毛细上升与接触角**

毛细管内上升高度
$$ h=\frac{2\sigma\cos\theta}{\rho g R}. $$
推导：由杨-拉普拉斯得弯月面压差 $\Delta p=2\sigma\cos\theta/R$（近似），与静压 $\rho g h$ 平衡。$\square$ 应用：多孔介质、微流控。

**3. 接触角与润湿**

杨方程 $\sigma_{sv}=\sigma_{sl}+\sigma_{lv}\cos\theta$。$\theta<90^\circ$ 亲液，$\theta>90^\circ$ 疏液。推导：接触线处三界面张力水平分力平衡。应用：涂层、超疏水。

**4. 液滴与气泡动力学**

瑞利-普莱塞特气泡自然振荡频率
$$ f=\frac1{2\pi R_0}\sqrt{\frac{3\gamma p_0}{\rho_l}}. $$
推导：气泡半径扰动 $R=R_0(1+\epsilon)$，由瑞利-普莱塞特方程
$$ R\ddot R+\tfrac32\dot R^2=\frac{p_{gas}-p_\infty}{\rho_l}-\frac{2\sigma}{\rho_lR} $$
线性化得谐振，$\omega^2=3\gamma p_0/(\rho_lR_0^2)$。应用：超声、空化。

**5. 颗粒沉降与斯托克斯阻力**

低雷诺数球体阻力 $F_D=6\pi\mu aU$，终端速度
$$ U_t=\frac{2(\rho_p-\rho_f)ga^2}{9\mu}. $$
推导：绕球斯托克斯解得阻力 $F_D=6\pi\mu aU$，重力-浮力=阻力，解出 $U_t$。$\square$ 应用：沉降、离心分离、血液。

**6. 两相流流型与基本参数**

流型：泡状、弹状、环状、雾状。含气率 $\alpha=A_g/A$，表观速度 $j_g=Q_g/A$，滑移比 $S=u_g/u_l$。均相模型混合物密度
$$ \rho_m=\alpha\rho_g+(1-\alpha)\rho_l. $$
推导：体积平均。应用：管道输送、锅炉、油气。

### 第14章 湍流统计理论与数值方法基础

**1. 湍流的基本特征**

随机性、多尺度、耗散、扩散、非线性。雷诺数足够大时层流失稳。湍流脉动能量从大尺度向小尺度传递（能级串），最终粘性耗散。

**2. 雷诺应力输运方程与封闭问题**

$$ \frac{D\overline{u_i'u_j'}}{Dt}=P_{ij}+\Pi_{ij}+D_{ij}-\varepsilon_{ij}, $$
生成项 $P_{ij}=-\overline{u_i'u_k'}\partial_k\bar u_j-\overline{u_j'u_k'}\partial_k\bar u_i$。未知量与方程数不匹配，需湍流模型（$k$-$\varepsilon$、$k$-$\omega$、雷诺应力模型）。推导：对脉动动量方程乘脉动速度取平均。应用：CFD湍流建模。

**3. 大涡模拟（LES）**

对 N-S 做空间滤波：
$$ \frac{\partial\bar u_i}{\partial t}+\bar u_j\frac{\partial\bar u_i}{\partial x_j}=-\frac1\rho\frac{\partial\bar p}{\partial x_i}+\nu\nabla^2\bar u_i-\frac{\partial\tau_{ij}^{\rm SGS}}{\partial x_j}, $$
亚格子应力 $\tau_{ij}^{\rm SGS}=\overline{u_iu_j}-\bar u_i\bar u_j$。斯马戈林斯基模型
$$ \tau_{ij}^{\rm SGS}-\tfrac13\tau_{kk}\delta_{ij}=-2\nu_t\bar S_{ij},\quad \nu_t=(C_s\Delta)^2|\bar S|. $$
推导：滤波 N-S 得滤波项与亚格子项，采用涡粘假设封闭。应用：分离流、气动噪声。

**4. 直接数值模拟（DNS）**

直接求解所有尺度，网格需满足 $\Delta x\sim\eta$（柯尔莫哥洛夫尺度），网格数 $\sim Re^{9/4}$。推导：由 $\eta=L Re^{-3/4}$，三维网格 $(L/\eta)^3=Re^{9/4}$。应用：湍流机理研究。

**5. 有限体积法离散**

对控制方程在控制体积上积分：
$$ \int_{V}\frac{\partial(\rho\phi)}{\partial t}dV+\oint_S\rho\phi\mathbf u\cdot\mathbf n\,dA=\oint_S\Gamma\nabla\phi\cdot\mathbf n\,dA+\int_V S_\phi\,dV. $$
离散得代数方程组 $a_P\phi_P=\sum a_{nb}\phi_{nb}+b$。推导：高斯定理转面积分，对流通量用迎风/中心格式，扩散用中心差分。应用：CFD 主流软件。

**6. 有限差分与稳定性**

时间推进显式格式受 CFL 条件
$$ C=\frac{|u|\Delta t}{\Delta x}\le1. $$
推导：冯·诺依曼稳定性分析，放大因子 $|G|\le1$。应用：数值模拟步长选择。

**7. 涡方法与格子玻尔兹曼简介**

涡方法用离散涡元追踪涡量输运（比奥-萨伐尔求速度）。格子玻尔兹曼解离散速度分布 $f_i$：
$$ f_i(\mathbf x+\mathbf e_i\Delta t,t+\Delta t)=f_i(\mathbf x,t)+\Omega_i(f). $$
推导：玻尔兹曼方程 BGK 碰撞近似离散化，查普曼-恩斯科格展开恢复 N-S。应用：多孔、多相。

**8. 湍流边界层的对数律与摩阻预测**

见第10章第8条。工程摩阻关联与 $k$-$\varepsilon$ 壁面函数法应用。应用：飞行器阻力、管网能耗。

**9. 数值误差与验证**

离散误差、迭代误差、舍入误差；网格收敛性（GCI）。验证（与解析/实验）与确认（物理模型）。应用：CFD质量保证。

**10. 流体力学的守恒律数值格式（黎曼求解器）**

一维 Euler 方程守恒形式
$$ \frac{\partial\mathbf U}{\partial t}+\frac{\partial\mathbf F}{\partial x}=0,\quad \mathbf U=[\rho,\rho u,E]^{\rm T}. $$
Godunov 格式在界面解黎曼问题，通量 $\mathbf F_{i+1/2}=\mathbf F(\mathbf U^*)$。迎风/通量分裂（Roe、HLLC）处理激波。应用：可压缩CFD、激波捕获。

---
⬆ 上级领域：[[力学]]
