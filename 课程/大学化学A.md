---
课程名称: 大学化学A
学分: 2
总学时: 32
考核方式: 考试
课程属性: 选修
课程性质: 公共基础课程
开课单位: 化学与化工学院
领域: 物质科学（物理·化学·生命）
tags:
  - 课程
  - 领域/物质科学（物理·化学·生命）
---

# 大学化学A


### 第1章 气体的状态方程与分子运动论

**知识点1：理想气体状态方程**

- 定义：理想气体是忽略分子本身体积与分子间作用力、且分子间碰撞完全弹性的假想气体。描述其状态参量 $(p,V,n,T)$ 之间关系的方程称为状态方程。
- 假设：分子无体积；分子间除弹性碰撞外无相互作用；分子数足够多以统计平均；温度以热力学温标计。
- 公式：$pV=nRT$，其中 $R=8.314\ \mathrm{J\,mol^{-1}K^{-1}}$。
- 推导：由三条经验定律综合。
  Boyle 定律（$T,n$ 不变）：$V\propto \dfrac{1}{p}$；
  Charles 定律（$p,n$ 不变）：$V\propto T$；
  Avogadro 定律（$p,T$ 不变）：$V\propto n$。
  三者相乘得
$$
\because\quad V\propto n\cdot T\cdot \frac{1}{p}
\ \Rightarrow\ V=R\,\frac{nT}{p}
\ \therefore\quad pV=nRT.
$$
  量纲验证：$[pV]=(\mathrm{Pa}\,\mathrm{m^3})=\mathrm{J}$，$[nRT]=\mathrm{mol}\cdot\mathrm{J\,mol^{-1}K^{-1}}\cdot\mathrm{K}=\mathrm{J}$，一致。
- 意义：$R$ 是联系宏观量 $(p,V,T)$ 与微观量（分子数）的普适常数；可改写为 $p=\dfrac{N}{V}k_BT$，$k_B=R/N_A$。
- 特例：$n=1$、标准状况（$273.15\ \mathrm{K},101.325\ \mathrm{kPa}$）下摩尔体积 $V_m=22.414\ \mathrm{L\,mol^{-1}}$。
- 应用：气体密度 $\rho=\dfrac{pM}{RT}$，可用于由密度测摩尔质量 $M=\dfrac{\rho RT}{p}$。

**知识点2：道尔顿分压定律与阿马加分体积定律**

- 定义：混合气体中某组分单独占据总体积时产生的压力称为分压 $p_i$；某组分单独处于总压时占有的体积称为分体积 $V_i$。
- 假设：各组分互不反应，均服从理想气体状态方程。
- 公式：$p_i=x_ip_{\text{总}}$，$V_i=x_iV_{\text{总}}$，$x_i=\dfrac{n_i}{\sum n_j}$。
- 推导：
$$
\because\quad p_iV=n_iRT,\qquad p_{\text{总}}V=\Big(\sum_j n_j\Big)RT
$$
$$
\therefore\quad \frac{p_i}{p_{\text{总}}}=\frac{n_i}{\sum_j n_j}=x_i
\ \Rightarrow\ p_i=x_ip_{\text{总}},\qquad p_{\text{总}}=\sum_i p_i.
$$
  同理对分体积，
$$
\because\quad p_{\text{总}}V_i=n_iRT,\qquad p_{\text{总}}V=\Big(\sum_j n_j\Big)RT
\ \therefore\quad V_i=x_iV.
$$
- 意义：总压等于各组分分压之和；分体积之和等于总体积。
- 特例：只对理想气体严格成立；真实气体在低压高温下近似成立。
- 应用：气体分析中由体积分数求摩尔分数；呼吸气体分压计算。

**知识点3：气体分子运动论的基本方程（压强公式）**

- 定义：从分子微观运动出发，求大量分子对器壁碰撞产生的宏观压强。
- 假设：分子数为 $N$，体积 $V$，分子质量 $m$；分子速度各向同性（$\overline{v_x^2}=\overline{v_y^2}=\overline{v_z^2}=\frac{1}{3}\overline{v^2}$）；分子间碰撞与器壁碰撞均为完全弹性；分子间距离远大于分子线度。
- 公式：$p=\dfrac{1}{3}\dfrac{N}{V}m\overline{v^2}=\dfrac{2}{3}\dfrac{N}{V}\overline{\varepsilon}_k$，其中 $\overline{\varepsilon}_k=\frac12m\overline{v^2}$。
- 推导：设沿 $x$ 方向运动、速度为 $v_x$ 的分子碰壁后以 $-v_x$ 反弹，单次动量变化 $\Delta p_x=2mv_x$。在时间 $\Delta t$ 内能碰到面积 $A$ 器壁的分子位于底面积 $A$、高 $v_x\Delta t$ 的柱体内，数目为 $\frac{N}{V}A v_x\Delta t$。合力
$$
F=\frac{N}{V}A v_x\Delta t\cdot\frac{2mv_x}{\Delta t}=2\frac{N}{V}Amv_x^2.
$$
  对速度分布取平均并对三个方向求和（只计 $v_x>0$，系数已含），
$$
F=\frac{N}{V}Am\overline{v^2}\cdot\frac{1}{3}\cdot 2\cdot\frac12
=\frac{1}{3}\frac{N}{V}A m\overline{v^2}
$$
  故
$$
p=\frac{F}{A}=\frac{1}{3}\frac{N}{V}m\overline{v^2}.
$$
- 意义：宏观压强是分子对器壁碰撞动量传递的统计结果，与分子数密度和平均平动动能成正比。
- 特例：若各方向速度不等，则应写 $p=\frac{N}{V}m\,\overline{v_x^2}$。
- 应用：导出玻意耳定律与温度微观解释。

**知识点4：温度的微观本质与能量均分定理**

- 定义：温度是分子平均平动动能的量度。
- 公式：$\overline{\varepsilon}_k=\dfrac{3}{2}k_BT$；能量均分定理每一可积二次自由度贡献 $\frac12k_BT$。
- 推导：比较 $pV=nRT=Nk_BT$ 与 $pV=\dfrac{1}{3}Nm\overline{v^2}$：
$$
\because\quad pV=\frac{1}{3}Nm\overline{v^2}=Nk_BT
\ \Rightarrow\ \frac{1}{2}m\overline{v^2}=\frac{3}{2}k_BT.
$$
  对每个二次型自由度 $q$，哈密顿量含 $\frac12 a q^2$，由 Boltzmann 分布
$$
\overline{\frac12 a q^2}=\frac{\displaystyle\int \frac12 aq^2e^{-\beta aq^2/2}dq}{\displaystyle\int e^{-\beta aq^2/2}dq}=\frac{1}{2\beta}=\frac12k_BT,\quad \beta=\frac{1}{k_BT}.
$$
- 意义：温度只与平均平动动能有关，与分子种类无关。
- 特例：单原子理想气体 $U=\frac32nRT$；双原子（常温转动）$U=\frac52nRT$。
- 应用：计算理想气体内能与热容 $C_{V,m}=\frac32R$（单原子）。

**知识点5：麦克斯韦速率分布律**

- 定义：$f(v)\,dv$ 表示速率在 $[v,v+dv]$ 内分子数占总数的比例。
- 假设：平衡态、经典统计、分子间无作用（或作用可忽略）。
- 公式：$f(v)=4\pi\Big(\dfrac{m}{2\pi k_BT}\Big)^{3/2}v^2e^{-mv^2/(2k_BT)}$。
- 推导：由 Boltzmann 分布，速度处于 $d^3v$ 的分子数
$$
dN=N\Big(\frac{m}{2\pi k_BT}\Big)^{3/2}e^{-m(v_x^2+v_y^2+v_z^2)/(2k_BT)}dv_xdv_y dv_z.
$$
  在速度空间取半径 $v$、厚度 $dv$ 的球壳，体积 $d^3v=4\pi v^2dv$，得
$$
f(v)dv=\Big(\frac{m}{2\pi k_BT}\Big)^{3/2}e^{-mv^2/(2k_BT)}4\pi v^2dv.
$$
- 意义：分布存在最概然速率，速率大或小的分子数均少。
- 特例：最概然速率 $v_p=\sqrt{2k_BT/m}$（令 $df/dv=0$）；平均速率 $\bar v=\sqrt{8k_BT/(\pi m)}$；方均根速率 $v_{rms}=\sqrt{3k_BT/m}$，且 $v_p<\bar v<v_{rms}$。
- 应用：计算分子碰撞频率、平均自由程。

**知识点6：范德华方程**

- 定义：考虑分子体积与分子间引力的真实气体状态方程。
- 假设：分子本身占体积 $b$；分子间存在吸引使实测压强减小，修正量正比于浓度平方。
- 公式：$\Big(p+\dfrac{an^2}{V^2}\Big)(V-nb)=nRT$。
- 推导：体积修正：分子可自由活动空间由 $V$ 减为 $V-nb$，故理想气体方程中 $V\to V-nb$。压强修正：内层分子受周围分子净引力指向气体内部，使撞击器壁的速度减小，实测压强 $p$ 比"理想压强" $p_{id}$ 小。每个分子受引力正比于密度 $n/V$，碰壁分子数也正比于 $n/V$，故
$$
p_{id}-p\propto \Big(\frac{n}{V}\Big)^2\ \Rightarrow\ p_{id}=p+\frac{an^2}{V^2}.
$$
  代入理想方程 $p_{id}(V-nb)=nRT$，即得
$$
\Big(p+\frac{an^2}{V^2}\Big)(V-nb)=nRT.
$$
- 意义：$a$ 反映分子间引力，$b$ 反映分子体积；低压时 $an^2/V^2\to0$、$nb\ll V$，退化为理想气体。
- 特例：$a=b=0$ 得 $pV=nRT$；等温线在低温下出现极大极小与拐点（临界点），临界条件 $\partial p/\partial V=0$、$\partial^2 p/\partial V^2=0$。
- 应用：估算临界参数 $V_c=3nb$，$p_c=a/(27b^2)$，$T_c=8a/(27Rb)$。

**知识点7：压缩因子与对应状态原理**

- 定义：压缩因子 $Z=\dfrac{pV_m}{RT}$ 衡量真实气体对理想的偏差。
- 公式：$Z=1+\dfrac{B}{V_m}+\dfrac{C}{V_m^2}+\cdots$（维里方程）。
- 推导：由范德华方程 $p=\dfrac{RT}{V_m-b}-\dfrac{a}{V_m^2}$，
$$
Z=\frac{pV_m}{RT}=\frac{V_m}{V_m-b}-\frac{a}{V_mRT}
=\frac{1}{1-b/V_m}-\frac{a}{V_mRT}.
$$
  当 $b/V_m\ll1$，用 $(1-x)^{-1}=1+x+x^2+\cdots$，
$$
Z=1+\Big(b-\frac{a}{RT}\Big)\frac{1}{V_m}+\frac{b^2}{V_m^2}+\cdots
$$
  故 $B=b-a/(RT)$，$C=b^2$。令 $Z$ 用临界量的约化量表示，$p_r=p/p_c$，$V_r=V_m/V_c$，$T_r=T/T_c$，代入可消去 $a,b$，得 $Z=f(p_r,T_r)$。
- 意义：所有气体在同一对比状态 $(p_r,T_r)$ 下 $Z$ 近似相同——对应状态原理。
- 特例：$Z=1$ 理想气体；$Z<1$ 引力占优；$Z>1$ 体积排斥占优。
- 应用：无 $a,b$ 数据时用 $Z$ 图查算真实气体 $p,V,T$。

### 第2章 热力学第一定律与热化学

**知识点1：系统、状态函数与功热**

- 定义：被研究的对象为系统，与其相互作用的外部为环境。状态函数只取决于状态而与路径无关（$U,H,S,G,A$）；功与热是过程量。
- 假设：系统宏观静止、无整体动能与位能变化（热力学约定）。
- 公式：体积功 $\delta W=-p_{\text{外}}\,dV$。
- 推导：活塞面积 $A$，外压 $p_{\text{外}}$ 推动活塞位移 $dx$，对外做功
$$
\delta W=F_{\text{外}}dx=(p_{\text{外}}A)\frac{dV}{A}=p_{\text{外}}dV
$$
  按热力学惯例系统对外做功取负，故
$$
\delta W=-p_{\text{外}}dV.
$$
  有限过程 $W=-\displaystyle\int_{V_1}^{V_2}p_{\text{外}}\,dV$。
- 意义：功的正负取决于 $dV$ 与 $p_{\text{外}}$；状态函数全微分可写 $dU$。
- 特例：恒外压 $W=-p_{\text{外}}(V_2-V_1)$；自由膨胀 $p_{\text{外}}=0$，$W=0$。
- 应用：判断过程是否可逆（可逆过程 $p_{\text{外}}=p\pm dp$，做功最大）。

**知识点2：热力学第一定律与可逆过程**

- 定义：能量守恒在热现象中的表述：$\Delta U=Q+W$（$Q$ 吸热为正，$W$ 对外做功为负）。
- 假设：孤立系统内能守恒。
- 公式：$dU=\delta Q+\delta W=\delta Q-p_{\text{外}}dV$。
- 推导：对循环过程，系统回到原态 $\oint dU=0$，故 $\oint\delta Q=-\oint\delta W$，即净吸热等于对外净做功。对一般过程，由状态函数 $U$ 的性质直接有 $\Delta U=Q+W$。
- 意义：第一类永动机不可能。
- 特例：绝热 $\delta Q=0\Rightarrow \Delta U=W$；恒容 $dV=0\Rightarrow \Delta U=Q_V$。
- 应用：计算反应内能变与体积功。

**知识点3：焓与等压热效应**

- 定义：$H\equiv U+pV$。
- 公式：$\Delta H=Q_p$（恒压、只做体积功）。
- 推导：
$$
\because\quad H=U+pV\ \Rightarrow\ \Delta H=\Delta U+\Delta(pV)
$$
  恒压且只做体积功时 $\Delta(pV)=p\Delta V$，而
$$
\Delta U=Q_p-p\Delta V
\ \Rightarrow\ \Delta U+p\Delta V=Q_p
\ \therefore\quad \Delta H=Q_p.
$$
- 意义：化学反应多在恒压下进行，$\Delta H$ 直接可测（量热）。
- 特例：理想气体 $\Delta H=\Delta U+\Delta n_gRT$。
- 应用：由 $\Delta H$、$\Delta U$ 相互换算：$\Delta H=\Delta U+\Delta n_gRT$。

**知识点4：热容及 $C_p$ 与 $C_V$ 关系**

- 定义：$C_V=\Big(\dfrac{\partial U}{\partial T}\Big)_V$，$C_p=\Big(\dfrac{\partial H}{\partial T}\Big)_p$。
- 公式：理想气体 $C_p-C_V=nR$。
- 推导：理想气体 $U=U(T)$，$H=U+pV=U+nRT$，故
$$
C_p=\Big(\frac{\partial H}{\partial T}\Big)_p=\Big(\frac{\partial U}{\partial T}\Big)_p+nR=C_V+nR.
$$
- 意义：等压加热时系统还需对外做膨胀功，故 $C_p>C_V$。
- 特例：单原子理想气体 $C_{V,m}=\frac32R$，$C_{p,m}=\frac52R$，$\gamma=C_p/C_V=5/3$。
- 应用：绝热可逆方程推导。

**知识点5：绝热可逆过程方程**

- 定义：绝热且可逆（准静态）且理想气体的过程。
- 公式：$pV^\gamma=\text{const}$，$TV^{\gamma-1}=\text{const}$，$T^\gamma p^{1-\gamma}=\text{const}$，$\gamma=C_p/C_V$。
- 推导：绝热 $dU=\delta W$，即 $nC_{V,m}dT=-p\,dV$。由 $pV=nRT$ 微分 $p\,dV+V\,dp=nR\,dT$，
$$
\because\quad nC_{V,m}dT=-p\,dV,\qquad nR\,dT=p\,dV+V\,dp
$$
$$
\therefore\quad \frac{C_{V,m}}{R}(p\,dV+V\,dp)=-p\,dV
\ \Rightarrow\ C_{V,m}(p\,dV+V\,dp)=-Rp\,dV.
$$
  代入 $R=C_p-C_V$，
$$
C_V p\,dV+C_V V\,dp=-(C_p-C_V)p\,dV
\ \Rightarrow\ C_V V\,dp=-C_p p\,dV
$$
$$
\Rightarrow\quad \frac{dp}{p}=-\gamma\frac{dV}{V}
\ \therefore\quad \ln p+\gamma\ln V=\text{const}\ \Rightarrow\ pV^\gamma=\text{const}.
$$
- 意义：绝热线比等温线陡（斜率 $-p/(\gamma V)$ vs $-p/V$）。
- 特例：$\gamma=1$（等温）退化为 Boyle 定律。
- 应用：内燃机循环、声速计算 $c=\sqrt{\gamma p/\rho}$。

**知识点6：盖斯定律**

- 定义：反应热只取决于始末态，与途径无关。
- 假设：恒压或恒容，只做体积功。
- 公式：$\Delta H_{\text{总}}=\sum_i\Delta H_i$。
- 推导：$H$ 是状态函数，$\Delta H=\displaystyle\int_{\text{始}}^{\text{末}}dH$ 只与端点有关。若把反应 $A\to B$ 分成两步 $A\to C\to B$，
$$
\because\quad \oint dH=0\ \Rightarrow\ \Delta H_{A\to B}=\Delta H_{A\to C}+\Delta H_{C\to B}.
$$
- 意义：可用已知反应焓组合求难测反应焓。
- 特例：燃烧焓、生成焓均可作为"标准反应"。
- 应用：$\Delta_rH^\ominus=\sum\nu_i\Delta_fH_i^\ominus$。

**知识点7：基尔霍夫定律**

- 定义：反应焓随温度变化的关系。
- 公式：$\dfrac{d\Delta_rH^\ominus}{dT}=\Delta_rC_p^\ominus$。
- 推导：对 $p$ 恒定，$H=H(T)$，
$$
\frac{d}{dT}\Big[H(\text{产物})-H(\text{反应物})\Big]
=\sum_i\nu_i\frac{dH_i}{dT}
=\sum_i\nu_iC_{p,i}=\Delta_rC_p.
$$
  积分：
$$
\Delta_rH^\ominus(T_2)=\Delta_rH^\ominus(T_1)+\int_{T_1}^{T_2}\Delta_rC_p^\ominus\,dT.
$$
- 意义：温度区间不大时 $\Delta_rC_p$ 近似常数，$\Delta_rH(T_2)=\Delta_rH(T_1)+\Delta_rC_p\Delta T$。
- 特例：$\Delta_rC_p=0$ 时反应焓与温度无关。
- 应用：由 298 K 数据求高温反应焓。

**知识点8：生成焓与键焓**

- 定义：标准摩尔生成焓 $\Delta_fH_m^\ominus$ 为由参考态单质生成 1 mol 化合物的焓变；键焓为断 1 mol 键的焓变。
- 公式：$\Delta_rH^\ominus=\sum\nu_i\Delta_fH_i^\ominus=-\sum(\text{键焓})$。
- 推导：构造循环：反应物 $\to$ 组成单质 $\to$ 产物。由盖斯定律
$$
\Delta_rH^\ominus=\sum_{\text{产物}}\nu\Delta_fH^\ominus-\sum_{\text{反应物}}\nu\Delta_fH^\ominus.
$$
  断开反应物全部键吸热，形成产物键放热，故反应焓等于断键吸收减去成键放出。
- 意义：可用于估算反应焓。
- 特例：参考态单质 $\Delta_fH^\ominus=0$。
- 应用：估算未知化合物生成焓。

### 第3章 热力学第二定律与熵

**知识点1：卡诺循环与热机效率**

- 定义：由两个等温可逆与两个绝热可逆过程组成的可逆循环。
- 公式：$\eta=\dfrac{|W|}{Q_1}=1-\dfrac{T_2}{T_1}$。
- 推导：等温 $T_1$ 膨胀吸热 $Q_1=nRT_1\ln(V_2/V_1)$；绝热 $T_1\to T_2$ 有 $T_1V_2^{\gamma-1}=T_2V_3^{\gamma-1}$；等温 $T_2$ 压缩放热 $|Q_2|=nRT_2\ln(V_3/V_4)$；绝热回程 $T_2V_4^{\gamma-1}=T_1V_1^{\gamma-1}$。由两绝热式相除
$$
\frac{V_2}{V_1}=\frac{V_3}{V_4}\ \Rightarrow\ \frac{Q_1}{T_1}=\frac{|Q_2|}{T_2}.
$$
  故
$$
\eta=\frac{Q_1-|Q_2|}{Q_1}=1-\frac{|Q_2|}{Q_1}=1-\frac{T_2}{T_1}.
$$
- 意义：效率只由两个热源温度决定，与工作物质无关。
- 特例：$T_2=0\ \mathrm{K}$ 时 $\eta=1$，但绝对零度不可达。
- 应用：所有热机效率上限。

**知识点2：卡诺定理与克劳修斯不等式**

- 定义：在相同两热源间工作的一切热机，可逆者效率最大。
- 公式：$\eta_{\text{不可逆}}\le\eta_{\text{可逆}}=1-T_2/T_1$，即 $\dfrac{Q_1}{T_1}+\dfrac{Q_2}{T_2}\le0$。
- 推导：用反证法。若某热机效率高于卡诺热机，则联合其驱动逆卡诺热机，可使热从低温流向高温而不留其他变化，违反第二定律。故 $\eta\le\eta_C$。对任意循环，把循环用无穷多条绝热线分割成无穷多个卡诺微循环，每个微循环 $\dfrac{\delta Q_1}{T_1}+\dfrac{\delta Q_2}{T_2}\le0$，求和（取极限）
$$
\oint\frac{\delta Q}{T}\le0.
$$
- 意义：可逆过程 $\oint \delta Q/T=0$，定义熵；不可逆过程 $\oint \delta Q/T<0$。
- 特例：可逆绝热 $\delta Q=0\Rightarrow dS=0$。
- 应用：判断过程方向。

**知识点3：熵的定义与克劳修斯不等式**

- 定义：可逆过程中 $dS\equiv \dfrac{\delta Q_{\text{rev}}}{T}$；$S$ 是状态函数。
- 公式：$dS\ge\dfrac{\delta Q}{T}$，等号可逆。
- 推导：设从态 1 到态 2 有可逆路径 $a$ 与不可逆路径 $b$，组成循环 $1\xrightarrow{a}2\xrightarrow{b}1$。由克劳修斯不等式
$$
\oint\frac{\delta Q}{T}=\int_1^2\frac{\delta Q_{\text{rev}}}{T}-\int_1^2\frac{\delta Q_{\text{irr}}}{T}\le0
$$
$$
\therefore\quad \int_1^2 dS=\int_1^2\frac{\delta Q_{\text{rev}}}{T}\ge\int_1^2\frac{\delta Q_{\text{irr}}}{T}
\ \Rightarrow\ dS\ge\frac{\delta Q}{T}.
$$
- 意义：熵变可沿任一可逆路径计算，与真实路径无关。
- 特例：孤立系统 $\delta Q=0\Rightarrow dS\ge0$（熵增原理）。
- 应用：计算任意过程的熵变。

**知识点4：熵变的计算**

- 定义：$S$ 是状态函数，可设计可逆路径。
- 公式：理想气体
$$
\Delta S=nC_{V,m}\ln\frac{T_2}{T_1}+nR\ln\frac{V_2}{V_1}
=nC_{p,m}\ln\frac{T_2}{T_1}-nR\ln\frac{p_2}{p_1}.
$$
- 推导：由 $dU=\delta Q_{\text{rev}}-p\,dV$，$\delta Q_{\text{rev}}=nC_{V,m}dT+p\,dV$，
$$
dS=\frac{nC_{V,m}dT}{T}+\frac{p}{T}dV
=\frac{nC_{V,m}dT}{T}+nR\frac{dV}{V}.
$$
  积分得第一式。再用 $pV=nRT$ 换元得第二式。
- 意义：理想气体熵依赖 $T$ 与 $V$（或 $p$）。
- 特例：等温 $\Delta S=nR\ln(V_2/V_1)$；等容 $\Delta S=nC_V\ln(T_2/T_1)$；相变 $\Delta S=\Delta H_{\text{相变}}/T_{\text{相变}}$。
- 应用：设计等温+等容两步替代实际过程。

**知识点5：玻尔兹曼熵与统计意义**

- 定义：$S=k_B\ln\Omega$，$\Omega$ 为宏观态对应的微观态数。
- 公式：$S=k_B\ln\Omega$。
- 推导：设两独立系统，$\Omega=\Omega_1\Omega_2$，熵应相加 $S=S_1+S_2$。唯一满足 $f(xy)=f(x)+f(y)$ 的函数为对数，故 $S=k_B\ln\Omega$。由等温自由膨胀：$N$ 个分子体积加倍，每分子微观选择加倍，$\Omega=2^N$，
$$
\Delta S=k_B\ln 2^N=Nk_B\ln2=nR\ln2,
$$
  与热力学 $\Delta S=nR\ln(V_2/V_1)=nR\ln2$ 一致。
- 意义：熵是系统混乱度的量度。
- 特例：完美晶体 $T=0$ 时 $\Omega=1\Rightarrow S=0$（热力学第三定律）。
- 应用：解释自发过程向概率大的方向进行。

**知识点6：吉布斯自由能与亥姆霍兹自由能**

- 定义：$G\equiv H-TS=U+pV-TS$；$A\equiv U-TS$。
- 公式：恒温恒压 $dG\le0$，恒温恒容 $dA\le0$。
- 推导：对任意过程 $dS\ge\delta Q/T$，即 $\delta Q\le TdS$。由第一定律 $\delta W\le dU-TdS$。
  恒温恒容（$dT=0,dV=0$，只做体积功 $W=0$）：
$$
\because\quad dU-\delta Q\le TdS\ \Rightarrow\ d(U-TS)\le0
\ \therefore\quad dA\le0.
$$
  恒温恒压（$p$ 恒定，只做体积功 $\delta W=-p\,dV$）：
$$
\because\quad dU+p\,dV-TdS\le0\ \Rightarrow\ d(H-TS)\le0
\ \therefore\quad dG\le0.
$$
- 意义：$G$ 是恒温恒压自发性的判据；$A$ 是恒温恒容判据。
- 特例：$\Delta G=0$ 平衡；$\Delta G<0$ 自发；$\Delta G>0$ 非自发。
- 应用：相平衡、化学平衡的方向判断。

**知识点7：热力学基本方程与麦克斯韦关系**

- 定义：联系状态函数微分的方程。
- 公式：
$$
dU=TdS-p\,dV,\quad dH=TdS+Vdp,\quad dA=-SdT-p\,dV,\quad dG=-SdT+Vdp.
$$
- 推导：由第一、第二定律 $\delta Q_{\text{rev}}=TdS$，
$$
dU=TdS-p\,dV.
$$
  对 $H=U+pV$，$dH=dU+pdV+Vdp=TdS+Vdp$。对 $A=U-TS$，$dA=dU-TdS-SdT=-SdT-p\,dV$。对 $G=H-TS$，$dG=dH-TdS-SdT=-SdT+Vdp$。
  因 $dG=-SdT+Vdp$ 为全微分，混合偏导相等：
$$
-\Big(\frac{\partial S}{\partial p}\Big)_T=\Big(\frac{\partial V}{\partial T}\Big)_p,
$$
  同理
$$
\Big(\frac{\partial T}{\partial V}\Big)_S=-\Big(\frac{\partial p}{\partial S}\Big)_V,\quad
\Big(\frac{\partial T}{\partial p}\Big)_S=\Big(\frac{\partial V}{\partial S}\Big)_p,\quad
\Big(\frac{\partial S}{\partial V}\Big)_T=\Big(\frac{\partial p}{\partial T}\Big)_V.
$$
- 意义：可由可测量 $(p,V,T)$ 求不可测的量（如熵随体积、压强变化）。
- 特例：理想气体 $(\partial U/\partial V)_T=0$ 可由此证明。
- 应用：推导内能、焓随 $p,V$ 的关系。

### 第4章 化学平衡

**知识点1：标准平衡常数与反应商**

- 定义：对反应 $0=\sum\nu_i B_i$，标准平衡常数
$$
K^\ominus=\prod_i \Big(\frac{p_i/p^\ominus}{}\Big)^{\nu_i}\ (\text{气相})
=\prod_i a_i^{\nu_i}.
$$
  反应商 $Q$ 用任意时刻活度写相同形式。
- 假设：理想气体或理想溶液，恒温。
- 公式：$\Delta_rG_m=\Delta_rG_m^\ominus+RT\ln Q$。
- 推导：对气体组分，$\mu_i=\mu_i^\ominus+RT\ln(p_i/p^\ominus)$。反应的 $G$
$$
\Delta_rG_m=\sum_i\nu_i\mu_i=\sum_i\nu_i\mu_i^\ominus+RT\sum_i\nu_i\ln(p_i/p^\ominus)
=\Delta_rG_m^\ominus+RT\ln\prod_i(p_i/p^\ominus)^{\nu_i}=\Delta_rG_m^\ominus+RT\ln Q.
$$
- 意义：$Q$ 与 $K^\ominus$ 的相对大小决定方向。
- 特例：平衡时 $\Delta_rG_m=0,\ Q=K^\ominus$。
- 应用：判断反应方向与限度。

**知识点2：$\Delta_rG_m^\ominus=-RT\ln K^\ominus$**

- 定义：标准平衡常数与标准反应吉布斯能的关系。
- 公式：$\Delta_rG_m^\ominus=-RT\ln K^\ominus$，即 $K^\ominus=e^{-\Delta_rG_m^\ominus/RT}$。
- 推导：平衡时 $\Delta_rG_m=0$ 且 $Q=K^\ominus$，代入 $\Delta_rG_m=\Delta_rG_m^\ominus+RT\ln Q$：
$$
0=\Delta_rG_m^\ominus+RT\ln K^\ominus
\ \therefore\quad \Delta_rG_m^\ominus=-RT\ln K^\ominus.
$$
- 意义：$\Delta_rG_m^\ominus<0\Rightarrow K^\ominus>1$，反应在标准态下有利。
- 特例：$\Delta_rG_m^\ominus=0\Rightarrow K^\ominus=1$。
- 应用：由热力学数据算 $K^\ominus$。

**知识点3：范特霍夫方程**

- 定义：平衡常数随温度的变化。
- 公式：$\dfrac{d\ln K^\ominus}{dT}=\dfrac{\Delta_rH_m^\ominus}{RT^2}$。
- 推导：由 $\ln K^\ominus=-\Delta_rG_m^\ominus/(RT)$，故
$$
\frac{d\ln K^\ominus}{dT}=-\frac{1}{R}\frac{d}{dT}\Big(\frac{\Delta_rG_m^\ominus}{T}\Big).
$$
  利用 Gibbs–Helmholtz 公式 $\dfrac{d(G/T)}{dT}=-\dfrac{H}{T^2}$，得
$$
\frac{d\ln K^\ominus}{dT}=-\frac{1}{R}\Big(-\frac{\Delta_rH_m^\ominus}{T^2}\Big)
=\frac{\Delta_rH_m^\ominus}{RT^2}.
$$
  定积分：$\ln\dfrac{K_2^\ominus}{K_1^\ominus}=\dfrac{\Delta_rH_m^\ominus}{R}\Big(\dfrac{1}{T_1}-\dfrac{1}{T_2}\Big)$。
- 意义：吸热反应升温 $K$ 增大，放热反应升温 $K$ 减小。
- 特例：$\Delta_rH_m^\ominus$ 为常数时上式严格。
- 应用：由两个温度下的 $K$ 求反应焓。

**知识点4：勒夏特列原理**

- 定义：平衡系统受外界扰动时，平衡向削弱该扰动的方向移动。
- 公式：对浓度扰动 $Q$ 与 $K$ 比较；对压强扰动由 $\Delta n_g$ 判断。
- 推导：设压强增大 $\lambda$ 倍，$Q'=\lambda^{\Delta\nu_g}Q$。若 $\Delta\nu_g<0$，则 $Q'<Q=K$，系统向正向移动以消除压强增大；若 $\Delta\nu_g>0$，则 $Q'>K$，向逆向移动。温度扰动由 van't Hoff 方程判断。
- 意义：统一的定性判据。
- 特例：$\Delta\nu_g=0$ 时压强不影响平衡。
- 应用：工业合成氨条件选择。

### 第5章 化学动力学

**知识点1：反应速率的定义与速率方程**

- 定义：对反应 $aA+bB\to cC+dD$，反应速率
$$
v=-\frac{1}{a}\frac{d[A]}{dt}=-\frac{1}{b}\frac{d[B]}{dt}=\frac{1}{c}\frac{d[C]}{dt}=\frac{1}{d}\frac{d[D]}{dt}.
$$
  经验速率方程 $v=k[A]^m[B]^n$，$m,n$ 为级数。
- 假设：恒容、均相。
- 公式：$v=k\prod_i[B_i]^{n_i}$。
- 推导：各物种生成速率之比等于化学计量系数之比，故选 $-d[A]/(a\,dt)$ 消除计量系数，使 $v$ 唯一定义。
- 意义：反应级数由实验确定，不一定等于计量系数。
- 特例：零级 $v=k$；一级 $v=k[A]$；二级 $v=k[A]^2$。
- 应用：由实验数据确定级数。

**知识点2：零级、一级、二级反应积分速率方程**

- 定义：将微分速率方程积分，得浓度随时间关系。
- 公式：
$$
[A]=[A]_0-kt\ (\text{零级}),\quad [A]=[A]_0e^{-kt}\ (\text{一级}),\quad \frac{1}{[A]}=\frac{1}{[A]_0}+kt\ (\text{二级}).
$$
- 推导：零级 $-\dfrac{d[A]}{dt}=k$，积分 $[A]=[A]_0-kt$。一级 $-\dfrac{d[A]}{dt}=k[A]$，
$$
\int_{[A]_0}^{[A]}\frac{d[A]}{[A]}=-k\int_0^t dt\ \Rightarrow\ \ln\frac{[A]}{[A]_0}=-kt\ \Rightarrow\ [A]=[A]_0e^{-kt}.
$$
  二级 $-\dfrac{d[A]}{dt}=k[A]^2$，
$$
\int \frac{d[A]}{[A]^2}=-k\int dt\ \Rightarrow\ \frac{1}{[A]}-\frac{1}{[A]_0}=kt.
$$
- 意义：半衰期 $t_{1/2}$ 特征：零级 $t_{1/2}=[A]_0/(2k)$；一级 $t_{1/2}=\ln2/k$ 与初始浓度无关；二级 $t_{1/2}=1/(k[A]_0)$。
- 特例：一级反应是所有核衰变与多数分解的基础。
- 应用：由线性化图判断级数。

**知识点3：阿伦尼乌斯方程与活化能**

- 定义：$k=A e^{-E_a/(RT)}$，$A$ 指前因子，$E_a$ 活化能。
- 公式：$\ln k=\ln A-\dfrac{E_a}{RT}$。
- 推导：由过渡态理论或碰撞理论，反应速率正比于能量超过阈值的分子比例。由 Boltzmann 分布，能量 $\ge E_a$ 的分子比例 $\propto e^{-E_a/(RT)}$，故
$$
k=A e^{-E_a/(RT)}.
$$
  对两温度作差：
$$
\ln\frac{k_2}{k_1}=\frac{E_a}{R}\Big(\frac{1}{T_1}-\frac{1}{T_2}\Big).
$$
- 意义：$E_a$ 越大温度敏感性越强。
- 特例：$E_a=0$ 时 $k=A$ 与温度无关。
- 应用：由 $\ln k$–$1/T$ 直线斜率求 $E_a$。

**知识点4：稳态近似与速率决定步**

- 定义：对高活性中间体，设其浓度近似不随时间变化。
- 公式：$\dfrac{d[I]}{dt}\approx0$。
- 推导：以 $A\xrightarrow{k_1}I\xrightarrow{k_2}B$ 为例，
$$
\frac{d[I]}{dt}=k_1[A]-k_2[I]\approx0\ \Rightarrow\ [I]=\frac{k_1}{k_2}[A].
$$
  故 $v=k_2[I]=k_1[A]$，与一级一致。
- 意义：把复杂机理化简为可测速率方程。
- 特例：$k_2\gg k_1$ 时稳态成立。
- 应用：推导酶催化、链反应速率方程。

**知识点5：催化作用**

- 定义：催化剂改变反应途径、降低活化能而不改变平衡。
- 公式：$k=A e^{-E_a/(RT)}$，$K^\ominus$ 不变。
- 推导：催化剂提供新路径，使 $E_a$ 由 $E_{a,1}$ 降为 $E_{a,2}$，速率比
$$
\frac{k_2}{k_1}=e^{(E_{a,1}-E_{a,2})/(RT)}.
$$
  因催化剂同等降低正逆反应活化能，$\Delta_rH^\ominus$ 不变，故 $K^\ominus$ 不变。
- 意义：只加速达到平衡，不改变平衡位置。
- 特例：负催化剂（阻化剂）抑制反应。
- 应用：工业催化、生物酶。

### 第6章 电化学

**知识点1：原电池与电极电势**

- 定义：将氧化还原反应拆为两个半反应，电子经外电路传递，产生电流的装置。
- 公式：$E_{\text{cell}}=E_{\text{正极}}-E_{\text{负极}}$。
- 推导：正极发生还原，负极发生氧化，电子从负极流向正极；外电路开路时两端电势差即电池电动势。
- 意义：把化学能转化为电能。
- 特例：标准氢电极 $E^\ominus=0$。
- 应用：由标准电极电势表计算 $E^\ominus$。

**知识点2：能斯特方程**

- 定义：浓度（分压）对电极电势的影响。
- 公式：对 $aA+bB\to cC+dD$，
$$
E=E^\ominus-\frac{RT}{zF}\ln Q=E^\ominus-\frac{0.0592}{z}\lg Q\ (298.15\ \mathrm{K}).
$$
- 推导：由 $\Delta_rG_m=\Delta_rG_m^\ominus+RT\ln Q$ 及 $\Delta_rG_m=-zFE$：
$$
-zFE=-zFE^\ominus+RT\ln Q
\ \Rightarrow\ E=E^\ominus-\frac{RT}{zF}\ln Q.
$$
  298.15 K 时 $\dfrac{RT}{F}\ln10=0.0592\ \mathrm{V}$。
- 意义：电池电动势随反应商变化。
- 特例：平衡时 $E=0$，$Q=K^\ominus$，得 $\lg K^\ominus=zFE^\ominus/(2.303RT)$。
- 应用：计算非标准态电池电动势。

**知识点3：电动势与吉布斯能**

- 定义：$\Delta_rG_m=-zFE$。
- 公式：$\Delta_rG_m^\ominus=-zFE^\ominus$。
- 推导：外电路可逆放电时电功 $W_{\text{电}}=-zFE$（$z$ 为转移电子数）。恒温恒压可逆过程最大非体积功等于 $-\Delta G$：
$$
\Delta_rG_m=-W_{\text{电,最大}}=-zFE.
$$
- 意义：连接电化学与热力学。
- 特例：$E>0$ 即 $\Delta G<0$，反应自发。
- 应用：由 $E^\ominus$ 求平衡常数。

**知识点4：法拉第电解定律**

- 定义：电解时电解析出物质的质量与通过电量成正比。
- 公式：$m=\dfrac{QM}{zF}=\dfrac{ItM}{zF}$。
- 推导：转移 $z$ mol 电子析出 1 mol 物质，故析出 $n$ mol 需电量 $Q=znF=z\dfrac{m}{M}F$，整理得
$$
m=\frac{MQ}{zF}.
$$
- 意义：电量与物质量严格对应。
- 特例：$1\ \mathrm{F}=96485\ \mathrm{C\,mol^{-1}}$。
- 应用：电镀、电解计算电流效率。

**知识点5：金属腐蚀与防护**

- 定义：金属被氧化侵蚀的电化学过程。
- 公式：吸氧腐蚀 $O_2+2H_2O+4e^-\to4OH^-$；析氢腐蚀 $2H^++2e^-\to H_2$。
- 推导：金属表面形成微电池，阳极 $M\to M^{z+}+ze^-$，阴极发生氧还原或氢还原。两极接通使金属溶解。
- 意义：腐蚀是自发的电化学电池反应。
- 特例：锌块阴极保护（牺牲阳极）。
- 应用：镀层、缓蚀剂、外加电流保护。

### 第7章 原子结构与量子力学基础

**知识点1：波粒二象性与德布罗意关系**

- 定义：微观粒子既有粒子性又具有波动性。
- 公式：$\lambda=\dfrac{h}{p}=\dfrac{h}{mv}$。
- 推导：由光子能量 $E=h\nu$ 与动量 $p=E/c$，得 $p=h/\lambda$。de Broglie 假设实物粒子也满足该关系，即
$$
\lambda=\frac{h}{p}.
$$
- 意义：电子、中子等可发生衍射。
- 特例：宏观物体 $m$ 大，$\lambda$ 小到不可观测。
- 应用：电子显微镜分辨率 $\delta\sim\lambda$。

**知识点2：海森堡不确定原理**

- 定义：位置与动量不能同时精确确定。
- 公式：$\Delta x\,\Delta p_x\ge\dfrac{\hbar}{2}$，$\Delta E\,\Delta t\ge\dfrac{\hbar}{2}$。
- 推导：一维波包由波数分布 $\Delta k$ 叠加，位置不确定度 $\Delta x\sim1/\Delta k$。由 $p=\hbar k$，$\Delta p=\hbar\Delta k$，故 $\Delta x\Delta p\sim\hbar$。严格用 Fourier 变换的均方偏差得 $\Delta x\Delta p_x\ge\hbar/2$。
- 意义：不是仪器误差，而是波的本性。
- 特例：束缚态基态有零点能。
- 应用：解释自然线宽、隧穿。

**知识点3：从经典力学到薛定谔方程**

- 定义：量子态由波函数 $\psi(x,t)$ 描述，其模方为概率密度。
- 公式：$i\hbar\dfrac{\partial\psi}{\partial t}=\hat H\psi$。
- 推导：从经典哈密顿量 $H=\dfrac{p^2}{2m}+V(x)$（这正是牛顿第二定律的能量形式）。引入算符
$$
\hat p=-i\hbar\nabla,\qquad \hat E=i\hbar\frac{\partial}{\partial t},\qquad \hat H=-\frac{\hbar^2}{2m}\nabla^2+V.
$$
  将 $E=H$ 视为算符作用于 $\psi$ 的本征方程
$$
i\hbar\frac{\partial\psi}{\partial t}=\Big(-\frac{\hbar^2}{2m}\nabla^2+V\Big)\psi.
$$
  定态时 $\psi=\phi(\mathbf r)e^{-iEt/\hbar}$，代入得定态方程
$$
-\frac{\hbar^2}{2m}\nabla^2\phi+V\phi=E\phi.
$$
  这也与拉格朗日量 $L=T-V$ 通过勒让德变换 $H=\sum p_i\dot q_i-L$ 得到的经典哈密顿量一致，量子化即把泊松括号 $\{,\}$ 换为对易子 $\frac{1}{i\hbar}[\ ,\ ]$。
- 意义：量子力学的基本运动方程。
- 特例：自由粒子 $V=0$ 得平面波解。
- 应用：求解一切束缚态问题。

**知识点4：一维势箱的完整求解**

- 定义：粒子在 $0<x<L$ 内 $V=0$，箱外 $V=\infty$。
- 假设：波函数连续、单值、有限，边界处为零。
- 公式：$E_n=\dfrac{n^2h^2}{8mL^2}$，$\psi_n=\sqrt{\dfrac{2}{L}}\sin\dfrac{n\pi x}{L}$，$n=1,2,3,\dots$。
- 推导：箱内方程
$$
-\frac{\hbar^2}{2m}\frac{d^2\psi}{dx^2}=E\psi
\ \Rightarrow\ \frac{d^2\psi}{dx^2}+k^2\psi=0,\quad k^2=\frac{2mE}{\hbar^2}.
$$
  通解 $\psi=A\sin kx+B\cos kx$。边界 $\psi(0)=0\Rightarrow B=0$；$\psi(L)=0\Rightarrow A\sin kL=0$。非平凡解要求 $kL=n\pi$，故
$$
k=\frac{n\pi}{L}\ \Rightarrow\ E_n=\frac{\hbar^2k^2}{2m}=\frac{n^2h^2}{8mL^2}.
$$
  归一化：$\int_0^L|A|^2\sin^2(n\pi x/L)dx=1\Rightarrow |A|^2=2/L$。
- 意义：能量量子化是边界条件的结果；零点能 $E_1=h^2/(8mL^2)>0$ 与不确定原理自洽。
- 特例：$n$ 增大节点数 $n-1$，能级间隔 $\Delta E=\dfrac{(2n+1)h^2}{8mL^2}$。
- 应用：共轭分子光谱估算、量子点。

**知识点5：氢原子薛定谔方程与量子数**

- 定义：电子在核库仑势 $V=-e^2/(4\pi\varepsilon_0 r)$ 中的定态问题。
- 公式：$E_n=-\dfrac{Z^2}{n^2}\cdot13.6\ \mathrm{eV}$，$n=1,2,\dots$。
- 推导：球坐标下方程可分离 $\psi=R(r)Y(\theta,\phi)$。角向方程给出 $L^2=l(l+1)\hbar^2$，$l=0,1,\dots,n-1$；磁量子数 $m=-l,\dots,l$。径向方程
$$
-\frac{\hbar^2}{2m}\frac{1}{r^2}\frac{d}{dr}\Big(r^2\frac{dR}{dr}\Big)+\Big[V+\frac{l(l+1)\hbar^2}{2mr^2}\Big]R=ER.
$$
  要求在 $r\to0,\infty$ 有界，得到 $n=1,2,\dots$ 及
$$
E_n=-\frac{m e^4}{2(4\pi\varepsilon_0)^2\hbar^2}\frac{Z^2}{n^2}=-\frac{13.6Z^2}{n^2}\ \mathrm{eV}.
$$
- 意义：能级只由主量子数决定；轨道形状由 $l$，取向由 $m$，自旋由 $m_s$ 描述。
- 特例：基态 $n=1$，$E=-13.6\ \mathrm{eV}$，玻尔半径 $a_0=\dfrac{4\pi\varepsilon_0\hbar^2}{me^2}=0.529\ \text{Å}$。
- 应用：原子光谱、周期律。

**知识点6：多电子原子与屏蔽效应**

- 定义：电子间排斥使实际核电荷被部分屏蔽。
- 公式：$E_n=-13.6\dfrac{(Z-\sigma)^2}{n^2}\ \mathrm{eV}$，$\sigma$ 为屏蔽常数。
- 推导：近似把内层电子看作均匀负电云，对外层电子屏蔽 $\sigma$ 个核电荷，故有效核电荷 $Z^*=Z-\sigma$。代入类氢能级公式即得。
- 意义：能级随 $n$ 与 $l$ 分裂，发生能级交错。
- 特例：$n$ 相同 $l$ 越小越穿透，能量越低（$ns<np<nd$）。
- 应用：解释电子排布与周期表分区。

### 第8章 化学键与分子结构

**知识点1：价键理论与杂化轨道**

- 定义：共价键由两原子自旋相反的单电子配对、轨道重叠形成。
- 假设：电子定域于两核之间；最大重叠原理。
- 公式：键能来源于交换积分；sp、sp²、sp³ 杂化夹角 $180^\circ,120^\circ,109.5^\circ$。
- 推导：以 sp³ 为例，一个 2s 与三个 2p 线性组合
$$
\phi_i=a_is+b_i p_x+c_i p_y+d_i p_z,\quad i=1,\dots,4.
$$
  由正交归一与等价性条件 $\int\phi_i\phi_j d\tau=\delta_{ij}$，可解得四个等价方向互相成 $109.5^\circ$，其矢量和为零
$$
\sum_{i=1}^4 \mathbf u_i=0.
$$
- 意义：解释甲烷正四面体结构。
- 特例：sp 杂化直线形；sp² 平面三角形。
- 应用：预测分子几何与键角。

**知识点2：分子轨道理论与 $H_2^+$**

- 定义：分子轨道由原子轨道线性组合（LCAO）产生，电子填充分子轨道。
- 公式：$\psi_\pm=N(\phi_A\pm\phi_S)$，$E_\pm=\dfrac{H_{AA}\pm H_{AB}}{1\pm S_{AB}}$。
- 推导：变分法。试探波函数 $\psi=c_A\phi_A+c_B\phi_B$，能量
$$
E=\frac{\int\psi^*\hat H\psi\,d\tau}{\int\psi^*\psi\,d\tau}
=\frac{c_A^2H_{AA}+2c_Ac_BH_{AB}+c_B^2H_{BB}}{c_A^2+2c_Ac_BS+c_B^2},\quad S=\int\phi_A\phi_B d\tau.
$$
  令 $\partial E/\partial c_A=\partial E/\partial c_B=0$ 得久期方程。对称 $H_{AA}=H_{BB}$，解得
$$
E_\pm=\frac{H_{AA}\pm H_{AB}}{1\pm S}.
$$
- 意义：成键轨道能量低于原子轨道，反键轨道高于原子轨道。
- 特例：$H_2^+$ 一个电子填入 $\sigma_{1s}$，键级 $=1/2$。
- 应用：解释第二周期分子磁性、键级。

**知识点3：键级与键能**

- 定义：键级 $=\dfrac{1}{2}(\text{成键电子数}-\text{反键电子数})$。
- 公式：键能近似与键级正相关，键长与键级负相关。
- 推导：由分子轨道能级图，净成键电子数决定稳定化能
$$
E_{\text{稳定化}}=\sum_{i\in\text{occ}}n_iE_i.
$$
  键级是净成键电子数的一半。
- 意义：预测键长、键能、磁性。
- 特例：$O_2$ 键级 2 且有 2 个未成对电子，解释顺磁性（价键理论无法解释）。
- 应用：判断分子稳定性。

**知识点4：价层电子对互斥（VSEPR）**

- 定义：价层电子对尽可能远离。
- 公式：电子对数 $m=\dfrac{1}{2}(\text{中心原子价电子数}+\text{配体提供电子数}\pm\text{电荷})$。
- 推导：由库仑排斥最小化，等价点电荷在球面分布：$m=2$ 直线、$m=3$ 平面三角、$m=4$ 四面体、$m=5$ 三角双锥、$m=6$ 八面体。孤对占据较大角空间，压缩键角。
- 意义：预测分子形状（只算配体位置）。
- 特例：$H_2O$（4 对，2 孤对）键角 $104.5^\circ<109.5^\circ$。
- 应用：判断分子几何。

**知识点5：分子光谱与转动振动**

- 定义：分子吸收光子发生转动、振动、电子跃迁。
- 公式：转动能级 $E_J=\dfrac{\hbar^2}{2I}J(J+1)$；振动 $E_v=(v+\frac12)h\nu_0$。
- 推导：刚性转子哈密顿量 $\hat H=\dfrac{\hat L^2}{2I}$，本征值 $E_J=\dfrac{\hbar^2}{2I}J(J+1)$。谐振子势 $V=\frac12kx^2$，本征值 $E_v=(v+\frac12)\hbar\omega$，$\omega=\sqrt{k/\mu}$。
- 意义：红外光谱给出官能团信息。
- 特例：转动跃迁选律 $\Delta J=\pm1$，谱线间隔 $2B$，$B=\hbar^2/(2I)$。
- 应用：红外、拉曼、微波谱。

### 第9章 分子间作用力与晶体结构

**知识点1：范德华力与氢键**

- 定义：分子间定向力、诱导力、色散力统称范德华力；氢键为 X–H···Y 的强偶极作用。
- 公式：色散能 $E\approx-\dfrac{3}{4}\dfrac{\alpha_1\alpha_2I_1I_2}{(I_1+I_2)r^6}$。
- 推导：两瞬时偶极由涨落耦合，二阶微扰能量
$$
E^{(2)}=\sum_{n\ne0}\frac{|\langle 0|\hat V|n\rangle|^2}{E_0-E_n}\propto-\frac{\alpha_1\alpha_2}{r^6}.
$$
  故吸引势 $\propto-1/r^6$，配合近距排斥 $\propto1/r^{12}$ 得 Lennard-Jones 势
$$
V(r)=4\varepsilon\Big[\Big(\frac{\sigma}{r}\Big)^{12}-\Big(\frac{\sigma}{r}\Big)^{6}\Big].
$$
- 意义：决定沸点、溶解度、生物大分子构象。
- 特例：$\partial V/\partial r=0$ 得平衡距离 $r_0=2^{1/6}\sigma$。
- 应用：解释 $\mathrm{H_2O}$、HF、$\mathrm{NH_3}$ 异常高沸点。

**知识点2：离子键与晶格能（Born–Landé）**

- 定义：离子晶体中正负离子静电吸引形成离子键，晶格能为拆散 1 mol 晶体所需能量。
- 公式：$U=-\dfrac{N_A M z^+z^-e^2}{4\pi\varepsilon_0 r_0}\Big(1-\dfrac{1}{n}\Big)$。
- 推导：一对离子势能 $-\dfrac{z^+z^-e^2}{4\pi\varepsilon_0 r}$，Madelung 常数 $M=\sum_j \pm1/r_j$ 汇总所有离子对。排斥能 $B/r^n$，每个离子
$$
u(r)=-\frac{M z^+z^-e^2}{4\pi\varepsilon_0 r}+\frac{B}{r^n}.
$$
  平衡时 $\dfrac{du}{dr}=0$：
$$
\frac{M z^+z^-e^2}{4\pi\varepsilon_0 r_0^2}=\frac{nB}{r_0^{n+1}}
\ \Rightarrow\ B=\frac{M z^+z^-e^2r_0^{n-1}}{4\pi\varepsilon_0 n}.
$$
  代回得单对能量，再乘 $N_A$ 并入 1/2 修正计数，
$$
U=-\frac{N_A M z^+z^-e^2}{4\pi\varepsilon_0 r_0}\Big(1-\frac{1}{n}\Big).
$$
- 意义：解释离子晶体高熔点、高硬度。
- 特例：NaCl 型 $M=1.748$。
- 应用：由晶格能估算溶解度、稳定性。

**知识点3：晶体类型与离子极化**

- 定义：离子、原子、分子、金属晶体；离子极化指阳离子使阴离子电子云变形。
- 公式：极化力 $\propto\dfrac{z}{r^2}$；变形性 $\propto r^3$。
- 推导：外电场中原子感生偶极矩 $\mu=\alpha E$，$\alpha\approx4\pi\varepsilon_0 r^3$。阳离子电荷密度高则极化力强。
- 意义：极化增大使键由离子性向共价性过渡。
- 特例：AgI、ZnS 因极化显著而溶解度小。
- 应用：解释过渡金属化合物颜色与溶解度。

### 第10章 溶液与依数性

**知识点1：拉乌尔定律**

- 定义：稀溶液中溶剂蒸气压等于纯溶剂蒸气压乘以溶剂摩尔分数。
- 公式：$p_A=p_A^*x_A$。
- 推导：液面上溶剂分子蒸发速率正比于其表面分数，近似为摩尔分数 $x_A$，故蒸气压 $p_A\propto x_A$，比例系数由纯溶剂确定即 $p_A^*$。
- 意义：非挥发性溶质使蒸气压下降 $\Delta p=p_A^*x_B$。
- 特例：理想溶液全浓度范围成立。
- 应用：测摩尔质量、蒸馏。

**知识点2：亨利定律**

- 定义：稀溶液中挥发性溶质蒸气压正比于其摩尔分数。
- 公式：$p_B=k_x x_B$。
- 推导：溶质分子在液面附近的逃逸速率正比于其浓度，故 $p_B=kx_B$。与拉乌尔定律形式相似但比例系数为亨利常数。
- 意义：气体溶解度随分压增大。
- 特例：$k_x$ 随温度升高而增大（溶解度降低）。
- 应用：潜水减压病、碳酸饮料。

**知识点3：沸点升高与凝固点下降**

- 定义：非挥发性溶质使沸点升高、凝固点下降。
- 公式：$\Delta T_b=K_b b_B$，$\Delta T_f=K_f b_B$。
- 推导：由 Clausius–Clapeyron 方程 $\dfrac{d\ln p}{dT}=\dfrac{\Delta H_{\text{相变}}}{RT^2}$。沸点时溶液蒸气压等于外压，较纯溶剂需更高温度补偿蒸气压下降 $\Delta p=p_A^*x_B$。对稀溶液 $\ln x_A\approx-x_B$，积分得
$$
\Delta T_b=\frac{RT_b^{*2}}{\Delta_{\text{vap}}H_m}b_B\equiv K_b b_B.
$$
  凝固点类似：
$$
\Delta T_f=\frac{RT_f^{*2}}{\Delta_{\text{fus}}H_m}b_B\equiv K_f b_B.
$$
- 意义：依数性只与粒子数有关，与溶质种类无关。
- 特例：电解质需乘 van't Hoff 因子 $i$：$\Delta T=iKb$。
- 应用：测摩尔质量、防冻液。

**知识点4：渗透压**

- 定义：半透膜两侧溶液与纯溶剂间的压力差。
- 公式：$\Pi=cRT$（van't Hoff 公式）。
- 推导：把溶质视为理想气体在溶液体积中产生的"渗透压"。平衡时化学势相等，$\mu_A(\text{溶液})=\mu_A^*(\text{纯})+RT\ln x_A$；加外压使纯溶剂化学势升高 $V_{m,A}\Pi$。令化学势相等
$$
RT\ln x_A+V_{m,A}\Pi=0
\ \Rightarrow\ \Pi=-\frac{RT}{V_{m,A}}\ln x_A.
$$
  稀溶液 $\ln x_A=\ln(1-x_B)\approx-x_B\approx-n_B/n_A$，且 $n_AV_{m,A}=V$，故
$$
\Pi=\frac{n_BRT}{V}=c_BRT.
$$
- 意义：依数性中灵敏度最高。
- 特例：电解质 $\Pi=icRT$。
- 应用：测高聚物摩尔质量、反渗透。

### 第11章 酸碱平衡与沉淀溶解平衡

**知识点1：酸碱质子理论与水的电离**

- 定义：酸是质子供体，碱是质子受体；共轭酸碱对相差一个质子。
- 公式：$K_w=[H^+][OH^-]=1.0\times10^{-14}$（298 K），$\mathrm{pH}=-\lg[H^+]$。
- 推导：水自偶电离 $\mathrm{2H_2O\rightleftharpoons H_3O^++OH^-}$，平衡常数 $K^\ominus=a(\mathrm{H_3O^+})a(\mathrm{OH^-})$。稀溶液中活度近似浓度，故 $K_w=[H^+][OH^-]$。
- 意义：$\mathrm{pH+pOH=14}$。
- 特例：温度升高 $K_w$ 增大，中性 pH 下降。
- 应用：pH 计算。

**知识点2：弱酸弱碱电离平衡**

- 定义：弱电解质部分电离。
- 公式：$K_a=\dfrac{[H^+][A^-]}{[HA]}$；$[H^+]=\sqrt{K_ac}$（最简式）。
- 推导：设电离度 $\alpha$，初始浓度 $c$，
$$
HA\rightleftharpoons H^++A^-,\quad K_a=\frac{(c\alpha)^2}{c(1-\alpha)}=\frac{c\alpha^2}{1-\alpha}.
$$
  当 $\alpha\ll1$，$1-\alpha\approx1$，得 $\alpha=\sqrt{K_a/c}$，
$$
[H^+]=c\alpha=\sqrt{K_ac}.
$$
- 意义：稀释使电离度增大（越稀越电离）。
- 特例：$c/K_a>500$ 时最简式误差 $<5\%$。
- 应用：计算弱酸 pH。

**知识点3：缓冲溶液与 Henderson–Hasselbalch 方程**

- 定义：能抵抗外加少量酸碱或稀释而 pH 基本不变的溶液。
- 公式：$\mathrm{pH}=\mathrm{p}K_a+\lg\dfrac{[A^-]}{[HA]}$。
- 推导：由 $K_a=\dfrac{[H^+][A^-]}{[HA]}$，
$$
[H^+]=K_a\frac{[HA]}{[A^-]}
\ \Rightarrow\ -\lg[H^+]=-\lg K_a-\lg\frac{[HA]}{[A^-]}
\ \Rightarrow\ \mathrm{pH}=\mathrm{p}K_a+\lg\frac{[A^-]}{[HA]}.
$$
- 意义：缓冲容量在 $\mathrm{pH}=\mathrm{p}K_a$ 附近最大。
- 特例：$[A^-]=[HA]$ 时 $\mathrm{pH}=\mathrm{p}K_a$。
- 应用：生物体液、化学实验 pH 控制。

**知识点4：溶度积与溶解度**

- 定义：难溶电解质沉淀溶解平衡的平衡常数 $K_{sp}$。
- 公式：$\mathrm{A_mB_n}$：$K_{sp}=[A^{n+}]^m[B^{m-}]^n$。
- 推导：$\mathrm{A_mB_n(s)\rightleftharpoons mA^{n+}+nB^{m-}}$，平衡时固相活度为 1，故
$$
K_{sp}=a(A^{n+})^m a(B^{m-})^n\approx[A^{n+}]^m[B^{m-}]^n.
$$
- 意义：判断沉淀生成与溶解。
- 特例：同离子效应使溶解度降低。
- 应用：分步沉淀、沉淀转化。

**知识点5：沉淀的生成与转化**

- 定义：溶液离子积 $Q_{sp}$ 与 $K_{sp}$ 比较判断。
- 公式：$Q_{sp}>K_{sp}$ 生成沉淀；$Q_{sp}<K_{sp}$ 溶解。
- 推导：由反应商判据，$Q_{sp}=K_{sp}$ 平衡；增大 $Q_{sp}$ 使平衡向左移动析出沉淀。
- 意义：定量控制沉淀条件。
- 特例：沉淀转化由 $K_{sp}$ 与配位/酸碱副反应决定。
- 应用：分离金属离子。

### 第12章 配位化合物与晶体场理论

**知识点1：配位化合物的组成与命名**

- 定义：中心离子（原子）与配体通过配位键结合形成的复杂离子或分子。
- 公式：配位数、内界 $[ML_n]^{z\pm}$、外界离子。
- 推导：配体孤对电子填入中心离子空轨道形成 $\sigma$ 配位键，故配位数等于配体给予配位原子的电子对数。
- 意义：解释配合物性质。
- 特例：$\mathrm{[Cu(NH_3)_4]^{2+}}$ 配位数 4。
- 应用：分析化学、催化、生物金属酶。

**知识点2：晶体场理论与分裂能**

- 定义：配体静电场使中心离子 d 轨道能级分裂。
- 公式：八面体场中 $E_{e_g}-E_{t_{2g}}=\Delta_o=10Dq$；$E_{t_{2g}}=-4Dq$，$E_{e_g}=+6Dq$。
- 推导：d 轨道在八面体配位场中，$e_g$（$d_{z^2},d_{x^2-y^2}$）正对配体，受排斥大、能量升高；$t_{2g}$（$d_{xy},d_{yz},d_{xz}$）指向配体间隙，能量降低。分裂能守恒（重心不变）：
$$
2E_{e_g}+3E_{t_{2g}}=0,\qquad E_{e_g}-E_{t_{2g}}=10Dq
$$
  解得
$$
E_{e_g}=6Dq,\qquad E_{t_{2g}}=-4Dq.
$$
- 意义：解释配合物颜色与磁性。
- 特例：四面体场分裂能 $\Delta_t=\dfrac{4}{9}\Delta_o$，且 $e<t_2$ 顺序颠倒。
- 应用：判断高自旋/低自旋。

**知识点3：配位平衡**

- 定义：配合物生成与解离的平衡。
- 公式：$K_f=\dfrac{[ML_n]}{[M][L]^n}$（稳定常数），$K_d=1/K_f$。
- 推导：$\mathrm{M+nL\rightleftharpoons ML_n}$，平衡常数由质量作用定律
$$
K_f=\frac{[ML_n]}{[M][L]^n}.
$$
- 意义：稳定常数越大配合物越稳定。
- 特例：逐级稳定常数乘积等于总稳定常数 $\beta_n=\prod K_i$。
- 应用：配位滴定、掩蔽。

### 第13章 统计热力学基础

**知识点1：玻尔兹曼分布**

- 定义：平衡态下独立粒子在各能级上的最概然分布。
- 公式：$N_i=\dfrac{N}{q}g_ie^{-\varepsilon_i/(k_BT)}$，$q=\sum_i g_ie^{-\varepsilon_i/(k_BT)}$。
- 推导：在总粒子数 $\sum N_i=N$ 与总能量 $\sum N_i\varepsilon_i=U$ 约束下，使微观态数 $W=\dfrac{N!}{\prod N_i!}$ 最大。取对数并用 Stirling 公式 $\ln N!\approx N\ln N-N$，
$$
\ln W=N\ln N-\sum_i N_i\ln N_i.
$$
  用拉格朗日乘子 $\alpha,\beta$ 构建
$$
\mathcal L=\ln W-\alpha\sum_iN_i-\beta\sum_iN_i\varepsilon_i,
$$
  令 $\partial\mathcal L/\partial N_i=0$，得 $\ln(N_i/g_i)=-(1+\alpha)-\beta\varepsilon_i$，即 $N_i=g_ie^{-1-\alpha}e^{-\beta\varepsilon_i}$。由归一化消去 $e^{-1-\alpha}$，并由热力学恒等式识别 $\beta=1/(k_BT)$，得
$$
N_i=\frac{N g_ie^{-\varepsilon_i/(k_BT)}}{\sum_j g_je^{-\varepsilon_j/(k_BT)}}.
$$
- 意义：连接微观能级与宏观热力学。
- 特例：能级连续时过渡为麦克斯韦分布。
- 应用：计算内能、熵、平衡常数。

**知识点2：配分函数与热力学函数**

- 定义：$q=\sum_i g_ie^{-\beta\varepsilon_i}$，是状态和。
- 公式：$U=Nk_BT^2\Big(\dfrac{\partial\ln q}{\partial T}\Big)_V$；$S=Nk_B\ln q+\dfrac{U}{T}$。
- 推导：内能
$$
U=\sum_iN_i\varepsilon_i=\frac{N}{q}\sum_i\varepsilon_ig_ie^{-\beta\varepsilon_i}
=-\frac{N}{q}\frac{\partial q}{\partial\beta}
=-N\frac{\partial\ln q}{\partial\beta}
=Nk_BT^2\frac{\partial\ln q}{\partial T}.
$$
  由 $A=-k_BT\ln(Q/N!)$，用 Stirling 得 $A=-Nk_BT\ln(eq/N)$，再由 $S=(U-A)/T$：
$$
S=Nk_B\ln q+\frac{U}{T}+Nk_B.
$$
  更严格含 $N!$ 修正时 $S=Nk_B\Big(\ln\dfrac qN+1\Big)+\dfrac{U}{T}$。
- 意义：由 $q$ 计算全部热力学量。
- 特例：$q$ 可分解为平动、转动、振动、电子配分函数之积。
- 应用：由光谱数据计算热力学函数。

**知识点3：经典极限与能量均分**

- 定义：高温或大质量使能级密集，量子统计趋于经典。
- 公式：平动配分函数 $q_t=\Big(\dfrac{2\pi mk_BT}{h^2}\Big)^{3/2}V$。
- 推导：平动为三维势箱，能级 $\varepsilon_n=\dfrac{h^2}{8mL^2}(n_x^2+n_y^2+n_z^2)$。当能级密集求和对能级连续求和转为相空间积分：
$$
q_t=\sum_{n_x,n_y,n_z}e^{-\varepsilon/(k_BT)}
\approx\frac{1}{h^3}\int e^{-(p_x^2+p_y^2+p_z^2)/(2mk_BT)}d^3x\,d^3p.
$$
  位置积分给 $V$，动量高斯积分给 $(2\pi mk_BT)^{3/2}$，故
$$
q_t=\frac{V}{h^3}(2\pi mk_BT)^{3/2}=\Big(\frac{2\pi mk_BT}{h^2}\Big)^{3/2}V.
$$
  由此得平动内能 $U_t=\dfrac32Nk_BT$，即能量均分。
- 意义：量子统计在高温极限回到经典结果。
- 特例：低温或大能级间隔必须用量子求和。
- 应用：气体热容、化学反应平动配分函数。

### 第14章 表面现象与胶体化学

**知识点1：表面张力与表面能**

- 定义：使液体表面收缩的力为表面张力 $\gamma$。
- 公式：$\gamma=\Big(\dfrac{\partial G}{\partial A}\Big)_{T,p}$。
- 推导：恒温恒压下可逆扩大表面积 $dA$，需非体积功 $\delta W_{\text{表面}}=\gamma\,dA$，由 $dG\le\delta W_{\text{非体积}}$ 取可逆等号
$$
dG=\gamma\,dA\ \Rightarrow\ \gamma=\Big(\frac{\partial G}{\partial A}\Big)_{T,p}.
$$
- 意义：表面能是单位面积吉布斯能。
- 特例：小液滴由于附加压强更高，蒸气压更大（Kelvin 方程）：
$$
\ln\frac{p_r}{p_\infty}=\frac{2\gamma V_m}{rRT}.
$$
- 应用：润湿、毛细现象、乳化。

**知识点2：朗缪尔吸附等温式**

- 定义：单分子层吸附，吸附质在表面达到动态平衡。
- 假设：每个吸附位吸附一个分子；吸附热与覆盖度无关；分子间无作用。
- 公式：$\theta=\dfrac{bp}{1+bp}$。
- 推导：设吸附速率 $v_a=k_a p(1-\theta)$（空位比例 $1-\theta$），脱附速率 $v_d=k_d\theta$。平衡 $v_a=v_d$：
$$
k_ap(1-\theta)=k_d\theta\ \Rightarrow\ \theta=\frac{k_ap}{k_d+k_ap}=\frac{bp}{1+bp},\quad b=\frac{k_a}{k_d}.
$$
- 意义：描述气固吸附。
- 特例：低压 $\theta\approx bp$（正比）；高压 $\theta\to1$（饱和）。
- 应用：比表面测定、多相催化。

**知识点3：胶体的性质与结构**

- 定义：分散相粒子直径 1–100 nm 的分散体系。
- 公式：扩散、布朗运动、电泳；胶团结构。
- 推导：胶体粒子受流体分子不平衡碰撞做布朗运动，其均方位移由 Einstein 关系
$$
\overline{x^2}=2Dt,\qquad D=\frac{RT}{6\pi\eta N_A r}.
$$
- 意义：胶体粒子介于真溶液与粗分散系之间，具有巨大比表面。
- 特例：加入电解质使胶体聚沉，聚沉能力与离子价数有关（Schulze–Hardy 规则）。
- 应用：血液、涂料、纳米材料。

**知识点4：临界胶束浓度与表面活性剂**

- 定义：表面活性剂在溶液中形成胶束的最低浓度。
- 公式：$\mathrm{CMC}$ 附近表面张力—浓度曲线出现拐点。
- 推导：表面活性剂分子由亲水头与疏水尾组成，倾向富集于液面降低表面能；浓度增大后转为自组装成胶束以隐藏疏水尾，使体系自由能降低。胶束形成后表面张力不再显著下降。
- 意义：解释洗涤、乳化、增溶。
- 特例：临界胶束浓度随温度、电解质变化。
- 应用：去污、药物载体、微乳液。

### 第15章 传递现象与连续介质（力学衔接）

**知识点1：菲克第一定律与扩散**

- 定义：扩散通量正比于浓度梯度。
- 公式：$\mathbf J=-D\nabla c$。
- 推导：从唯象线性响应假设出发，或由随机行走 $D=\dfrac{1}{2}\dfrac{\langle\Delta x^2\rangle}{\Delta t}$。在体积元内对守恒方程
$$
\frac{\partial c}{\partial t}+\nabla\cdot\mathbf J=0
$$
  代入 $\mathbf J=-D\nabla c$，得菲克第二定律
$$
\frac{\partial c}{\partial t}=D\nabla^2c.
$$
- 意义：菲克定律是线性不可逆过程的一个典型。
- 特例：稳态 $\nabla^2c=0$（拉普拉斯方程）。
- 应用：材料渗碳、反应器传质。

**知识点2：傅里叶定律与热传导**

- 定义：热流密度正比于温度梯度。
- 公式：$\mathbf q=-k\nabla T$。
- 推导：与扩散类比，能量守恒 $\rho c_p\dfrac{\partial T}{\partial t}=-\nabla\cdot\mathbf q$，代入得
$$
\frac{\partial T}{\partial t}=\alpha\nabla^2T,\qquad \alpha=\frac{k}{\rho c_p}.
$$
- 意义：导热系数 $k$ 为物性参数。
- 特例：稳态一维平板 $\dfrac{d^2T}{dx^2}=0\Rightarrow T$ 线性分布，$q=k\Delta T/L$。
- 应用：热防护、传热。

**知识点3：牛顿黏性定律与纳维–斯托克斯方程**

- 定义：流体层间剪应力正比于速度梯度。
- 公式：$\tau=\mu\dfrac{du}{dy}$；不可压缩 N–S：
$$
\rho\Big(\frac{\partial\mathbf v}{\partial t}+\mathbf v\cdot\nabla\mathbf v\Big)
=-\nabla p+\mu\nabla^2\mathbf v+\mathbf f.
$$
- 推导：对流体微元应用牛顿第二定律。加速度
$$
\mathbf a=\frac{D\mathbf v}{Dt}=\frac{\partial\mathbf v}{\partial t}+\mathbf v\cdot\nabla\mathbf v.
$$
  表面力含压强梯度 $-\nabla p$ 与黏性应力散度 $\nabla\cdot\boldsymbol\tau$。对牛顿流体本构 $\boldsymbol\tau=\mu(\nabla\mathbf v+\nabla\mathbf v^T)-\frac23\mu(\nabla\cdot\mathbf v)\mathbf I$，不可压缩 $\nabla\cdot\mathbf v=0$ 时 $\nabla\cdot\boldsymbol\tau=\mu\nabla^2\mathbf v$。于是
$$
\rho\frac{D\mathbf v}{Dt}=-\nabla p+\mu\nabla^2\mathbf v+\mathbf f.
$$
- 意义：与连续性方程 $\dfrac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf v)=0$ 联立封闭。
- 特例：低雷诺数忽略惯性项得 Stokes 方程 $-\nabla p+\mu\nabla^2\mathbf v=0$；理想流体 $\mu=0$ 得欧拉方程。
- 应用：气动、化工流体、微流控；同时是化学场中反应-扩散-流动耦合的基础。

**知识点4：热力学稳定性与凸性判据**

- 定义：稳定平衡要求自由能对扰动取极小。
- 公式：$\Big(\dfrac{\partial^2 G}{\partial p^2}\Big)_T=-V\kappa_T<0$，即压缩系数正定。
- 推导：由 $dG=-SdT+Vdp$，二阶导
$$
\Big(\frac{\partial^2 G}{\partial p^2}\Big)_T=\Big(\frac{\partial V}{\partial p}\Big)_T=-V\kappa_T.
$$
  稳定态要求 $G$ 的 Hessian 负定（对 $p$ 凹），故 $\kappa_T>0$；否则微小压强扰动会自发放大，系统失稳。同理定容热容要求 $C_V>0$。
- 意义：把力学中势能极小/稳定平衡的判据推广到热力学。
- 特例：气液共存区 $(\partial p/\partial V)_T=0$ 处出现拐点，稳定性失效。
- 应用：相分离、临界现象、自催化反应器的热稳定性分析。

---
⬆ 上级领域：[[物质科学（物理·化学·生命）]]
