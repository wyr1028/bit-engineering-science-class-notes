---
课程名称: 电工和电子技术AⅡ
学分: 4
总学时: 64
考核方式: 考试
课程属性: 必修
课程性质: 专业基础课程
开课单位: 集成电路与电子学院
领域: 电子·控制·系统
tags:
  - 课程
  - 领域/电子·控制·系统
---

# 电工和电子技术AⅡ


### 第1章 电路模型、参考方向与网络基本定律

- **电流、电压与功率的参考方向**
  定义：电流 $i$ 为单位时间通过导体横截面的电荷量，$i=\dfrac{dq}{dt}$；电压 $u$ 为单位正电荷从一点移到另一点电场力所做的功 $u=\dfrac{dW}{dq}$。
  假设：集总参数假设，即电路尺寸远小于电磁波波长，可忽略电磁场沿导线的传播效应。
  公式：关联参考方向下吸收功率 $p=ui$；能量 $W=\int_{t_0}^{t}u i\,dt$。
  推导：由定义 $i=dq/dt$，$dW=u\,dq$，故
  $$\because \frac{dW}{dt}=u\frac{dq}{dt}\quad \therefore p=\frac{dW}{dt}=ui .$$
  物理意义：$p>0$ 表示元件吸收功率，$p<0$ 表示发出功率。
  极限：无源符号约定对线性与非线性元件均适用。
  应用：判断电源与负载、计算效率。

- **基尔霍夫电流定律（KCL）**
  定义：任一集总节点，流出电流的代数和为零。
  假设：电荷守恒，节点无电荷积累。
  公式：$\sum_{k=1}^{n} i_k=0$。
  推导：设节点内电荷 $q$，流入电流使 $q$ 增加，则
  $$\because \frac{dq}{dt}=\sum i_{\text{in}}-\sum i_{\text{out}},\qquad q=\text{const}\ \Rightarrow \frac{dq}{dt}=0 ,$$
  $$\therefore \sum i_{\text{in}}=\sum i_{\text{out}}\ \Leftrightarrow \sum_{k} i_k=0 .$$
  物理意义：电流连续性，电荷不灭。
  特例：广义节点（封闭面）KCL 同样成立。
  应用：节点电压法、电路节点分析。

- **基尔霍夫电压定律（KVL）**
  定义：任一集总回路，各支路电压代数和为零。
  假设：集总参数电路，电场沿闭合路径的环流为零（准静态）。
  公式：$\sum_{k=1}^{n} u_k=0$。
  推导：单位正电荷沿回路一周回到起点，总功为零：
  $$\because \oint_{\Gamma}\vec{E}\cdot d\vec{l}=-\frac{d\Phi_B}{dt}\approx 0\ (\text{准静态}),$$
  $$u_k=\int_{\text{支路}k}\vec{E}\cdot d\vec{l},\quad \therefore \sum u_k=\oint_{\Gamma}\vec{E}\cdot d\vec{l}=0 .$$
  物理意义：能量守恒，回路电位单值。
  特例：含时变磁通时 KVL 需引入感应电动势项。
  应用：回路电流法、桥式电路分析。

- **电阻串联、并联与分压分流**
  定义：串联时电流相同，并联时电压相同。
  公式：$R_{eq}=\sum R_k$（串）；$\dfrac{1}{R_{eq}}=\sum\dfrac{1}{R_k}$（并）。
  推导：串联由 KVL，$u=\sum R_k i=R_{eq}i$，故 $R_{eq}=\sum R_k$。
  并联由 KCL，$i=\sum u/R_k=u\sum 1/R_k$，故 $1/R_{eq}=\sum 1/R_k$。
  分压：$u_k=\dfrac{R_k}{\sum R}R$（总电压 $u$ 下），由 $i=u/\sum R_k$。分流：$i_k=\dfrac{G_k}{\sum G}i$。
  物理意义：串联等效电阻增大，并联减小。
  极限：两电阻并联 $R_{eq}=\dfrac{R_1R_2}{R_1+R_2}$。
  应用：电压表、电流表扩程。

- **电压源、电流源及其等效变换**
  定义：理想电压源出口电压恒定，理想电流源输出电流恒定，内阻分别为 $0$、$\infty$。
  公式：实电压源 $u=u_s-R_s i$；实电流源 $i=i_s-\dfrac{u}{R_s}$。
  推导：两式对应同一外特性 $u=u_s-R_s i$ 与 $i=i_s-u/R_s$，令 $i_s=u_s/R_s$ 即等效：
  $$\therefore u=u_s-R_s i=u_s-R_s i .$$
  物理意义：端口伏安特性相同即对外等效。
  特例：理想电压源短路电流无穷大；理想电流源开路电压无穷大，均不可直接等效。
  应用：化简含源一端口。

- **受控源模型**
  定义：受控电压源 $u=\mu u_x$ 或受控电流源 $i=\beta i_x$、$g_m u_x$、$r_m i_x$。
  公式：VCCS $i=g_m u$，CCVS $u=r_m i$ 等。
  推导：由线性四端网络，输出量与输入量线性相关，泰勒展开一阶项即得。
  物理意义：描述晶体管、运放等器件的线性受控行为。
  应用：小信号等效电路。

- **线性、叠加定理与齐次性**
  定义：线性电路响应与激励成线性关系。
  公式：若激励 $x_1\to y_1$，$x_2\to y_2$，则 $a x_1+b x_2\to a y_1+b y_2$。
  推导：线性方程组 $A y=x$，$y=A^{-1}x$，由 $A^{-1}$ 的线性性直接得 $\therefore y(a x_1+b x_2)=a y_1+b y_2$。
  物理意义：各独立源单独作用的响应可叠加。
  特例：功率非线性，不能叠加 $p\ne p_1+p_2$。
  应用：多电源电路分析。

- **戴维南与诺顿定理**
  定义：任意含源线性一端口可等效为电压源 $u_{oc}$ 与 $R_{eq}$ 串联，或电流源 $i_{sc}$ 与 $R_{eq}$ 并联，且 $i_{sc}=u_{oc}/R_{eq}$。
  推导：由叠加定理，端口电压
  $$u=u_{oc}-R_{eq}i .$$
  取 $i=0$ 得 $u=u_{oc}$；取 $u=0$ 得 $i_{sc}=u_{oc}/R_{eq}$。
  $$\therefore R_{eq}=\frac{u_{oc}}{i_{sc}},\qquad \text{诺顿源} i_{sc}=\frac{u_{oc}}{R_{eq}} .$$
  物理意义：一端口对外仅由两个参数表征。
  应用：求负载最大功率、信号源等效。

- **最大功率传输**
  定义：负载从给定一端口获得最大功率的条件。
  公式：$R_L=R_{eq}$，$P_{max}=\dfrac{u_{oc}^2}{4R_{eq}}$。
  推导：负载功率
  $$P_L=i^2R_L=\left(\frac{u_{oc}}{R_{eq}+R_L}\right)^2R_L .$$
  令 $\dfrac{dP_L}{dR_L}=0$：
  $$\frac{d}{dR_L}\frac{R_L}{(R_{eq}+R_L)^2}=\frac{(R_{eq}+R_L)^2-2R_L(R_{eq}+R_L)}{(R_{eq}+R_L)^4}=0,$$
  $$\therefore R_{eq}+R_L-2R_L=0\ \Rightarrow R_L=R_{eq},\quad P_{max}=\frac{u_{oc}^2}{4R_{eq}} .$$
  物理意义：阻抗匹配。
  特例：此时效率仅 $50\%$。
  应用：射频、传感器接口。

### 第2章 直流电阻电路的方程法分析

- **支路电流法与方程独立性**
  定义：以支路电流为未知量列 KCL、KVL 求解。
  公式：$b$ 条支路、$n$ 个节点，独立 KCL 方程 $n-1$ 个，独立 KVL 方程 $b-n+1$ 个。
  推导：节点电流方程组系数矩阵行和为 $0$，秩为 $n-1$；网孔（基本回路）数由图论 Euler 公式 $b-n+1$ 给出。
  物理意义：未知量与独立方程数相等，解唯一。
  应用：小规模电路手算。

- **节点电压法**
  定义：以节点电压为未知量，对独立节点列 KCL。
  公式：标准形式 $G_{ii}v_i-\sum_{j\ne i}G_{ij}v_j=i_{si}$。
  推导：支路电流 $i_{ij}=G_{ij}(v_i-v_j)$，代入 KCL $\sum i=0$，整理得
  $$\sum_j G_{ij}(v_i-v_j)=i_{si}\ \Rightarrow\ G_{ii}v_i-\sum_{j\ne i}G_{ij}v_j=i_{si}.$$
  物理意义：自导纳乘本节点电压减互导纳乘邻节点电压等于注入电流。
  特例：含电压源时用超节点。
  应用：大规模网络、SPICE 内核。

- **网孔电流法**
  定义：以网孔电流为未知量列 KVL。
  公式：$R_{ii}i_i-\sum R_{ij}i_j=u_{si}$。
  推导：网孔电流自动满足 KCL，对每个网孔列 KVL，支路电压用网孔电流差表示即得。
  物理意义：回路自由度 $b-n+1$。
  应用：平面电路。

- **互易定理**
  定义：线性无源网络激励与响应互换位置后比值不变。
  公式：$\dfrac{u_2}{i_1}\bigg|_{1\to2}=\dfrac{u_1}{i_2}\bigg|_{2\to1}$。
  推导：由线性网络 $Z$ 矩阵对称性 $Z_{12}=Z_{21}$（互易网络），$\therefore$ 传递函数对称。
  应用：测量与等效。

- **Tellegen 定理**
  定义：网络各支路电压电流乘积之和为零，$\sum u_k i_k=0$。
  推导：由 KCL、KVL，$u^T i=(A^T v_n)^T i=v_n^T A i=0$。
  物理意义：功率守恒的拓扑表达。
  应用：网络灵敏度、稳定性证明。

- **Δ–Y 变换**
  定义：三端无源网络两种接法对外等效。
  公式：$R_Y=\dfrac{R_\Delta\text{相邻两电阻之积}}{\sum R_\Delta}$。
  推导：令两种接法在任两端间等效电阻相等，联立三式：
  $$R_{ab}=R_a+R_b=\frac{R_{ab}^\Delta(R_{ac}^\Delta+R_{bc}^\Delta)}{R_{ab}^\Delta+R_{bc}^\Delta+R_{ca}^\Delta},$$
  解方程组即得 $R_a=\dfrac{R_{ab}R_{ca}}{R_{ab}+R_{bc}+R_{ca}}$ 等。
  应用：桥式电路化简。

### 第3章 正弦稳态交流电路的相量分析

- **正弦量的三要素与有效值**
  定义：$u(t)=U_m\cos(\omega t+\varphi)$，三要素为幅值、角频率、初相。
  公式：有效值 $U=\dfrac{U_m}{\sqrt2}$。
  推导：均方根
  $$U=\sqrt{\frac1T\int_0^T U_m^2\cos^2(\omega t+\varphi)dt}
  =\sqrt{\frac{U_m^2}{2}}=\frac{U_m}{\sqrt2},$$
  其中用 $\int_0^T\cos^2=T/2$。
  物理意义：有效值等于产生相同热效应的直流值。
  应用：220 V 市电为有效值。

- **相量表示与相量法**
  定义：$\dot U=U\angle\varphi$，对应 $u(t)=\sqrt2 U\cos(\omega t+\varphi)$。
  公式：$\dfrac{d}{dt}\leftrightarrow j\omega$，$\int dt\leftrightarrow \dfrac{1}{j\omega}$。
  推导：由 Euler 公式 $u=\mathrm{Re}[\sqrt2\dot U e^{j\omega t}]$，求导
  $$\frac{du}{dt}=\mathrm{Re}[\sqrt2(j\omega\dot U)e^{j\omega t}],\ \therefore \frac{d}{dt}\to j\omega .$$
  物理意义：微分方程化为代数方程。
  应用：交流网络分析。

- **阻抗与导纳**
  定义：$\dot Z=U/I=Z\angle\varphi_Z$，$Z=\sqrt{R^2+X^2}$，$\varphi_Z=\arctan(X/R)$。
  公式：$Z_R=R$，$Z_L=j\omega L$，$Z_C=\dfrac{1}{j\omega C}$。
  推导：R：$u=Ri\Rightarrow Z_R=R$；L：$u=L di/dt\Rightarrow Z_L=j\omega L$；C：$i=Cdu/dt\Rightarrow \dot U=\dot I/(j\omega C)$。
  物理意义：幅值比值与相位差。
  特例：串联谐振 $X=0$，并联谐振 $B=0$。
  应用：滤波器设计。

- **阻抗串并联与分压分流**
  公式：$\dot Z_{eq}=\sum\dot Z_k$（串），$1/\dot Z_{eq}=\sum1/\dot Z_k$（并）。
  推导：同 KCL/KVL 相量形式，代数结构同直流。
  应用：交流电桥。

- **瞬时功率、有功、无功与视在功率**
  定义：$p=ui$；有功 $P=\dfrac1T\int p\,dt$；无功 $Q$；视在 $S=UI$。
  公式：$p=UI\cos\varphi(1+\cos2\omega t)+UI\sin\varphi\sin2\omega t$；
  $$P=UI\cos\varphi,\quad Q=UI\sin\varphi,\quad S=UI,\quad S^2=P^2+Q^2 .$$
  推导：设 $u=\sqrt2U\cos\omega t$，$i=\sqrt2I\cos(\omega t-\varphi)$，用积化和差
  $$p=UI[\cos\varphi+\cos(2\omega t-\varphi)]=UI\cos\varphi(1+\cos2\omega t)+UI\sin\varphi\sin2\omega t .$$
  取平均 $\overline{\cos2\omega t}=\overline{\sin2\omega t}=0$，得 $P,Q$。
  物理意义：$Q$ 表征储能元件与电源往返交换。
  特例：纯电阻 $\varphi=0$，$Q=0$；纯电感 $\varphi=90^\circ$，$P=0$。
  应用：功率因数补偿。

- **复功率与功率守恒**
  定义：$\tilde S=\dot U\dot I^*=P+jQ$。
  推导：$\dot U=U\angle\varphi_u,\dot I=I\angle\varphi_i$，$\dot U\dot I^*=UI\angle(\varphi_u-\varphi_i)=UI\cos\varphi+jUI\sin\varphi$。
  物理意义：复功率满足守恒 $\sum\tilde S_k=0$（由 KCL/KVL）。
  应用：交流网络功率计算。

- **RLC 串联谐振**
  定义：$X_L=X_C$ 时电压电流同相。
  公式：$\omega_0=\dfrac{1}{\sqrt{LC}}$，$Q_0=\dfrac{\omega_0 L}{R}=\dfrac{1}{R}\sqrt{\dfrac LC}$，带宽 $BW=\omega_0/Q_0$。
  推导：$Z=R+j(\omega L-1/\omega C)$，令虚部为零
  $$\omega L-\frac1{\omega C}=0\ \Rightarrow\ \omega_0=\frac1{\sqrt{LC}} .$$
  品质因数 $Q_0=\dfrac{\omega_0L}{R}=\dfrac{U_L}{U}=\dfrac{U_C}{U}$。
  物理意义：谐振时电压最高、电流最大。
  应用：选频、无线电接收。

- **功率因数提高**
  定义：用并联电容提高 $\cos\varphi$。
  公式：$C=\dfrac{P}{\omega U^2}(\tan\varphi_1-\tan\varphi_2)$。
  推导：补偿前后有功不变 $P=UI_1\cos\varphi_1=UI_2\cos\varphi_2$，无功减少量
  $$\Delta Q=P(\tan\varphi_1-\tan\varphi_2)=\omega C U^2\ \Rightarrow\ C=\frac{P}{\omega U^2}(\tan\varphi_1-\tan\varphi_2).$$
  应用：工厂无功补偿。

### 第4章 三相交流电路

- **对称三相电源**
  定义：三个同幅值、同频率、相位互差 $120^\circ$ 的电压。
  公式：$\dot U_A=U\angle0^\circ$，$\dot U_B=U\angle-120^\circ$，$\dot U_C=U\angle120^\circ$。
  推导：$\dot U_A+\dot U_B+\dot U_C=U(1+e^{-j120^\circ}+e^{j120^\circ})=U(1-1)=0$。
  物理意义：对称三相瞬时和为零。
  应用：三相输电。

- **星形（Y）接法线电压与相电压**
  公式：线电压 $U_l=\sqrt3 U_p$，相位超前 $30^\circ$。
  推导：$\dot U_{AB}=\dot U_A-\dot U_B=\dot U_A(1-e^{-j120^\circ})=\sqrt3\dot U_A\angle30^\circ$。
  物理意义：可获得两种电压。
  特例：中线使负载不对称时相电压仍对称。
  应用：低压配电。

- **三角形（Δ）接法**
  公式：$U_l=U_p$，$I_l=\sqrt3 I_p$。
  推导：$\dot I_A=\dot I_{AB}-\dot I_{CA}$，对称时 $=\sqrt3\dot I_{AB}\angle-30^\circ$。
  应用：大功率负载。

- **三相功率**
  公式：$P=\sqrt3 U_l I_l\cos\varphi$，$Q=\sqrt3 U_l I_l\sin\varphi$，$S=\sqrt3 U_l I_l$。
  推导：三相总有功 $P=3U_pI_p\cos\varphi$，代入 $U_l=\sqrt3U_p$ 得 $P=\sqrt3U_lI_l\cos\varphi$。
  物理意义：与接法无关，仅用线量表示。
  应用：电动机铭牌。

- **不对称三相与中线作用**
  定义：负载阻抗不等即不对称。
  公式：中性点位移 $\dot U_{N'N}=\dfrac{\sum \dot U_k Y_k}{\sum Y_k}$。
  推导：对中性点列 KCL $\sum(\dot U_k-\dot U_{N'N})Y_k=0$，解出上式。
  物理意义：中线抑制中性点位移。
  应用：照明负载须用中线。

- **对称分量法初步**
  定义：任意三相量分解为正序、负序、零序。
  公式：$\dot U_A=\dot U_{A1}+\dot U_{A2}+\dot U_{A0}$，算子 $a=e^{j120^\circ}$。
  推导：由 $\dot U_B=a^2\dot U_{A1}+a\dot U_{A2}+\dot U_{A0}$ 等三式联立解得各序分量。
  应用：故障分析。

### 第5章 非正弦周期电流电路与频率响应

- **周期函数的傅里叶级数**
  定义：满足狄里赫利条件的周期函数可分解为各次谐波。
  公式：$f(t)=a_0+\sum_{n=1}^\infty(a_n\cos n\omega t+b_n\sin n\omega t)$。
  推导：利用三角函数正交性
  $$\int_0^T\cos m\omega t\cos n\omega t\,dt=\frac T2\delta_{mn},$$
  两端乘 $\cos n\omega t$ 并积分得 $a_n=\dfrac2T\int_0^T f\cos n\omega t\,dt$，同理 $b_n$。
  物理意义：频域分解。
  应用：整流波形分析。

- **非正弦周期量的有效值与功率**
  公式：$U=\sqrt{U_0^2+\sum U_n^2}$，$P=P_0+\sum P_n$。
  推导：由正交性，$\overline{f^2}=a_0^2+\sum\frac{a_n^2+b_n^2}{2}=U_0^2+\sum U_n^2$。
  物理意义：各次谐波功率（无效功率）独立叠加。
  应用：THD 计算。

- **网络函数与频率响应**
  定义：$H(j\omega)=\dfrac{\dot Y}{\dot X}$。
  公式：幅频 $|H|$，相频 $\angle H$。
  推导：对线性定常系统，正弦输入经相量运算得输出比，即 $H$。
  物理意义：稳态正弦传递特性。
  应用：滤波器、放大器带宽。

- **一阶 RC 低通与高通**
  公式：低通 $H=\dfrac{1}{1+j\omega RC}$，$f_c=\dfrac1{2\pi RC}$。
  推导：分压 $\dot U_o=\dot U_i\dfrac{1/j\omega C}{R+1/j\omega C}$，整理得 $H$。
  物理意义：$f_c$ 处幅值降至 $1/\sqrt2$（$-3\,\mathrm{dB}$）。
  应用：抗混叠滤波。

- **RLC 串联带通与 $Q$ 值**
  公式：$H=\dfrac{1}{1+jQ_0(\frac{\omega}{\omega_0}-\frac{\omega_0}{\omega})}$。
  推导：$H=\dfrac{R}{R+j(\omega L-1/\omega C)}$，分子分母同除 $R$ 并代入 $Q_0,\omega_0$。
  物理意义：$Q_0$ 越高选择性越好。
  应用：中频放大器。

- **波特图渐近线**
  定义：以 $\mathrm{dB}$ 表示幅频、对数频率表示横轴。
  公式：$20\lg|H|$，一阶极点斜率 $-20\,\mathrm{dB/dec}$。
  推导：$|H|=\dfrac1{\sqrt{1+(\omega/\omega_c)^2}}$，当 $\omega\gg\omega_c$，$20\lg|H|\approx-20\lg(\omega/\omega_c)$。
  应用：系统辨识。

### 第6章 动态电路的时域分析

- **电容与电感及其储能**
  公式：$i_C=C\dfrac{du_C}{dt}$，$u_L=L\dfrac{di_L}{dt}$；$W_C=\frac12Cu_C^2$，$W_L=\frac12Li_L^2$。
  推导：$W_C=\int u\,i\,dt=\int C u\,du=\frac12Cu^2$；$W_L=\int L i\,di=\frac12Li^2$。
  物理意义：储能元件状态连续，$u_C,i_L$ 不能跃变。
  应用：换路定则。

- **换路定则与初始值**
  公式：$u_C(0^+)=u_C(0^-)$，$i_L(0^+)=i_L(0^-)$。
  推导：能量不能突变（有限电流电压下），$\therefore W$ 连续。
  应用：求初始条件。

- **一阶 RC 电路零输入响应**
  公式：$u_C(t)=U_0e^{-t/\tau}$，$\tau=RC$。
  推导：$RC\dfrac{du_C}{dt}+u_C=0$，分离变量 $\dfrac{du_C}{u_C}=-\dfrac{dt}{RC}$，积分
  $$\ln u_C=-\frac t{RC}+K\ \Rightarrow u_C=U_0e^{-t/\tau}.$$
  物理意义：按指数衰减，$\tau$ 为时间常数。
  应用：放电电路。

- **一阶 RC 零状态与全响应**
  公式：$u_C=U_s(1-e^{-t/\tau})$；全响应 $=U_s+(U_0-U_s)e^{-t/\tau}$。
  推导：$RC\dfrac{du_C}{dt}+u_C=U_s$，通解 $u_C=U_s+Ae^{-t/\tau}$，代入初值定 $A=U_0-U_s$。
  物理意义：全响应 = 稳态分量 + 暂态分量 = 零输入 + 零状态。
  应用：数字电路上升时间。

- **一阶 RL 电路**
  公式：$i_L(t)=I_\infty+(I_0-I_\infty)e^{-t/\tau}$，$\tau=L/R$。
  推导：$L\dfrac{di_L}{dt}+Ri_L=U_s$，同法求解。
  应用：继电器、电感充放电。

- **二阶 RLC 串联电路的微分方程**
  公式：$LC\dfrac{d^2u_C}{dt^2}+RC\dfrac{du_C}{dt}+u_C=U_s$。
  推导：KVL $u_L+u_R+u_C=U_s$，$u_R=RC\dfrac{du_C}{dt}$，$u_L=LC\dfrac{d^2u_C}{dt^2}$，代入得式。
  物理意义：特征方程 $s^2+\dfrac RLs+\dfrac1{LC}=0$。

- **二阶电路三种阻尼**
  定义：$\alpha=\dfrac R{2L}$，$\omega_0=\dfrac1{\sqrt{LC}}$。
  公式：$\Delta=\alpha^2-\omega_0^2$；过阻尼 $\Delta>0$，临界 $\Delta=0$，欠阻尼 $\Delta<0$。
  推导：特征根 $s_{1,2}=-\alpha\pm\sqrt{\alpha^2-\omega_0^2}$；欠阻尼时 $s=-\alpha\pm j\omega_d$，$\omega_d=\sqrt{\omega_0^2-\alpha^2}$，解
  $$u_C=e^{-\alpha t}(A\cos\omega_d t+B\sin\omega_d t).$$
  物理意义：阻尼振荡。
  应用：RLC 滤波、振荡器。

### 第7章 拉普拉斯变换与网络函数

- **拉普拉斯变换定义与性质**
  定义：$F(s)=\mathcal L\{f(t)\}=\displaystyle\int_0^\infty f(t)e^{-st}dt$，$s=\sigma+j\omega$。
  性质：$\mathcal L\{f'\}=sF(s)-f(0^-)$；$\mathcal L\{\int f\}=F(s)/s$；卷积 $\mathcal L\{f*g\}=F(s)G(s)$。
  推导：分部积分
  $$\int_0^\infty f'e^{-st}dt=[fe^{-st}]_0^\infty+s\int_0^\infty fe^{-st}dt=sF(s)-f(0^-).$$
  应用：解微分方程。

- **常用变换对**
  公式：$\delta(t)\leftrightarrow1$，$1\leftrightarrow1/s$，$e^{-at}\leftrightarrow1/(s+a)$，$\cos\omega t\leftrightarrow s/(s^2+\omega^2)$。
  推导：$e^{-at}$ 直接积分 $\int_0^\infty e^{-(s+a)t}dt=\dfrac1{s+a}$。
  应用：逆变换查表。

- **运算等效电路**
  公式：电容 $U_C(s)=\dfrac{I_C(s)}{sC}+\dfrac{u_C(0^-)}s$；电感 $U_L(s)=sLI_L(s)-Li_L(0^-)$。
  推导：对 $i_C=Cdu_C/dt$ 取拉氏变换得 $I_C=sCU_C-Cu_C(0^-)$，整理。
  应用：直接列代数方程。

- **网络函数与冲激响应**
  定义：$H(s)=Y(s)/X(s)$，零状态条件下。
  公式：$h(t)=\mathcal L^{-1}\{H(s)\}$。
  推导：$Y=HX$，逆变换 $y=h*x=\int_0^t h(\tau)x(t-\tau)d\tau$，$\therefore h$ 为冲激响应。
  应用：系统辨识。

- **极点、零点与稳定性**
  定义：$H(s)=N(s)/D(s)$，$D$ 根为极点，$N$ 根为零点。
  判据：全部极点实部 $<0$ 则稳定。
  推导：极点 $s_k$ 对应 $e^{s_k t}$，若 $\mathrm{Re}(s_k)>0$ 则发散，$\therefore$ 须左半平面。
  应用：控制系统。

- **初值与终值定理**
  公式：$f(0^+)=\lim_{s\to\infty}sF(s)$，$f(\infty)=\lim_{s\to0}sF(s)$。
  推导：由微分性质取极限即得。
  应用：快速求稳态。

### 第8章 磁路与变压器

- **磁场基本量与安培环路定律**
  定义：磁感应强度 $\vec B$、磁场强度 $\vec H$、磁通 $\Phi=\int\vec B\cdot d\vec S$。
  公式：$\oint\vec H\cdot d\vec l=\sum Ni=NI$（全电流定律）。
  推导：由安培定律 $\oint\vec B\cdot d\vec l=\mu_0 I_{enc}$，含介质时 $\vec B=\mu\vec H$，得 $\oint\vec H\cdot d\vec l=NI$。
  应用：铁心磁路计算。

- **磁路欧姆定律**
  定义：磁动势 $F=NI$，磁阻 $R_m$。
  公式：$\Phi=\dfrac{F}{R_m}$，$R_m=\dfrac{l}{\mu A}$。
  推导：均匀磁路 $H l=NI$，$\Phi=BA=\mu HA=\dfrac{\mu A}{l}NI=\dfrac{F}{R_m}$。
  物理意义：类比电路欧姆定律。
  特例：气隙磁阻 $R_\delta=\dfrac{\delta}{\mu_0 A}$。
  应用：电机、变压器铁心。

- **铁磁材料的磁化与磁滞**
  定义：磁化曲线 $B$–$H$ 非线性，磁滞回线。
  公式：$B=\mu H$，$\mu$ 非常数；剩磁 $B_r$，矫顽力 $H_c$。
  物理意义：磁畴转向与钉扎。
  应用：永磁、硅钢片。

- **涡流与铁心损耗**
  公式：涡流损耗 $P_e\propto f^2B_m^2$，磁滞损耗 $P_h\propto fB_m^n$。
  推导：感应电动势 $E\propto fB_m$，涡流 $I\propto E$，$\therefore P_e\propto f^2B_m^2$。
  应用：叠片铁心与高频铁氧体。

- **理想变压器**
  定义：无损耗、无漏磁、$\mu\to\infty$。
  公式：$\dfrac{u_1}{u_2}=\dfrac{N_1}{N_2}=n$，$\dfrac{i_1}{i_2}=-\dfrac1n$，$Z_{in}=n^2Z_L$。
  推导：由 $\Phi$ 相同，$u=N\frac{d\Phi}{dt}$，比值即匝比；由磁动势平衡 $N_1i_1+N_2i_2=0$（理想），得电流关系；$u_1/u_2=n$，$i_1=-i_2/n$，$\therefore Z_{in}=u_1/i_1=n^2Z_L$。
  应用：变压、阻抗匹配。

- **变压器感应电动势公式**
  公式：$E_1=4.44fN_1\Phi_m$。
  推导：$\Phi=\Phi_m\sin\omega t$，$e_1=N_1\dfrac{d\Phi}{dt}=N_1\omega\Phi_m\cos\omega t$，有效值 $E_1=\dfrac{N_1\omega\Phi_m}{\sqrt2}=4.44fN_1\Phi_m$。
  应用：铁心设计。

- **实际变压器等效电路**
  定义：漏磁、铜损、铁损。
  公式：$Z_{eq}=R_1+jX_1+n^2(R_2'+jX_2')$，励磁支路 $R_m\parallel jX_m$。
  推导：将副边折算到原边，$E_2'=nE_2$，$\therefore$ 得 T 型等效电路。
  应用：电压调整率、效率 $\eta=\dfrac{P_2}{P_1}$。

### 第9章 半导体物理基础与 PN 结

- **本征与杂质半导体**
  定义：本征载流子 $n_i$；N 型多子电子、P 型多子空穴。
  公式：$n_i^2=BT^3e^{-E_g/(kT)}$；$n_0p_0=n_i^2$（质量作用定律）。
  推导：热平衡下电子空穴产生复合率相等，费米分布积分得 $n_i$，质量作用由费米能级统一化学势推出。
  应用：器件温度特性。

- **漂移与扩散**
  公式：漂移 $J_{drift}=q(n\mu_n+p\mu_p)E$；扩散 $J_{diff}=q(D_n\nabla n-D_p\nabla p)$。
  推导：玻尔兹曼分布下扩散系数与迁移率满足爱因斯坦关系 $D/\mu=kT/q$。
  应用：PN 结、BJT 电流。

- **PN 结形成与内建电势**
  定义：P、N 接触后扩散形成空间电荷区与内建电场。
  公式：$V_{bi}=\dfrac{kT}{q}\ln\dfrac{N_AN_D}{n_i^2}$。
  推导：平衡时净电流为零，$J_n=q\mu_nnE+qD_n\dfrac{dn}{dx}=0$，用 Einstein 关系
  $$\frac{1}{n}dn=-\frac q{kT}dV\ \Rightarrow\ V_{bi}=\frac{kT}{q}\ln\frac{n_{n0}}{n_{p0}}=\frac{kT}{q}\ln\frac{N_AN_D}{n_i^2}.$$
  物理意义：扩散与漂移动态平衡。
  应用：结电压、温度传感器。

- **耗尽层与势垒电容**
  公式：$W=\sqrt{\dfrac{2\varepsilon}{q}\left(\dfrac1{N_A}+\dfrac1{N_D}\right)(V_{bi}-V)}$；$C_j=\dfrac{\varepsilon A}{W}$。
  推导：泊松方程 $\dfrac{d^2V}{dx^2}=-\dfrac{\rho}{\varepsilon}$，在突变结中分段积分，电中性 $N_Ax_p=N_Dx_n$，解出 $W$。
  应用：变容二极管。

- **PN 结伏安特性（Shockley 方程）**
  公式：$I=I_s\left(e^{V/V_T}-1\right)$，$V_T=kT/q$。
  推导：少数载流子扩散，准中性区连续性方程 $\dfrac{d^2(\Delta p)}{dx^2}=\dfrac{\Delta p}{L_p^2}$，解得指数分布，在耗尽层边缘
  $$\Delta p(x_n)=\Delta p_n e^{-x/L_p},\quad I_p=qAD_p\frac{\Delta p_{n0}}{L_p}(e^{V/V_T}-1),\ \therefore I=I_s(e^{V/V_T}-1).$$
  物理意义：单向导电性。
  特例：$V\gg V_T$ 时 $I\approx I_se^{V/V_T}$；$V<0$ 时 $I\to-I_s$。
  应用：整流、检波。

- **PN 结击穿**
  定义：齐纳击穿（隧穿）与雪崩击穿（碰撞电离）。
  公式：击穿电压 $V_{BR}$ 随掺杂浓度升高而降低。
  应用：稳压二极管。

### 第10章 二极管及其应用电路

- **理想与折线模型**
  定义：理想二极管导通压降 $0$、反向电流 $0$；折线模型导通压降 $V_{on}$、正向电阻 $r_D$。
  公式：导通 $i=\dfrac{u-V_{on}}{R+r_D}$。
  推导：用分段线性近似真实指数曲线。
  应用：快速手算。

- **整流电路**
  定义：半波、全波（中心抽头）、桥式。
  公式：半波 $U_o=0.45U_2$；全波/桥式 $U_o=0.9U_2$。
  推导：设 $u_2=\sqrt2U_2\sin\omega t$，半波平均
  $$U_o=\frac1T\int_0^{T/2}\sqrt2U_2\sin\omega t\,dt=\frac{\sqrt2U_2}{\pi}=0.45U_2 .$$
  全波为半波两倍 $\therefore 0.9U_2$。
  应用：直流电源。

- **电容滤波**
  公式：$U_o\approx U_2$（全波），纹波 $\Delta U\approx\dfrac{I_L}{2fC}$。
  推导：放电期间 $i_L=C\dfrac{du}{dt}$，近似线性，$\Delta U=\dfrac{I_L}{C}\cdot\dfrac{1}{2f}$。
  应用：电源滤波。

- **稳压二极管电路**
  公式：$U_o=U_Z$，限流 $R=\dfrac{U_i-U_Z}{I_Z+I_L}$。
  应用：基准电压、钳位。

- **限幅与钳位电路**
  定义：利用二极管导通电压限制信号幅度。
  公式：上限 $U_Z+0.7$，下限 $-0.7$。
  应用：波形整形、保护。

- **二极管小信号模型**
  公式：$r_d=\dfrac{V_T}{I_D}$。
  推导：$g_d=\dfrac{dI}{dV}\big|_{I_D}=\dfrac{I_s}{V_T}e^{V/V_T}\approx\dfrac{I_D}{V_T}$，$\therefore r_d=1/g_d=V_T/I_D$。
  应用：小信号检波。

### 第11章 双极型晶体管及基本放大电路

- **BJT 结构与电流放大**
  定义：发射区、基区、集电区，NPN/PNP。
  公式：$I_E=I_B+I_C$，$\bar\beta=I_C/I_B$，$\alpha=I_C/I_E$。
  推导：由 KCL，$\therefore I_E=I_B+I_C$；定义电流放大系数。
  关系：$\beta=\dfrac{\alpha}{1-\alpha}$，$\alpha=\dfrac{\beta}{1+\beta}$。
  应用：信号放大。

- **放大区电流方程**
  公式：$I_C=I_S e^{V_{BE}/V_T}$，$I_B=I_C/\beta$。
  推导：发射结正偏注入电子，基区输运与复合比决定 $\beta$，集电结反偏收集，$\therefore I_C$ 受 $V_{BE}$ 指数控制。
  物理意义：压控电流源。
  应用：偏置设计。

- **三种组态与特性曲线**
  定义：共射（CE）、共集（CC）、共基（CB）。
  公式：CE 电压增益大且反相，CC 电压增益 $\approx1$，CB 电流增益 $\approx1$。
  应用：按需求选组态。

- **固定偏置与分压式偏置**
  公式：$I_B=\dfrac{V_{CC}-V_{BE}}{R_B}$；分压偏置 $V_B=\dfrac{R_2}{R_1+R_2}V_{CC}$，$I_C\approx\dfrac{V_B-0.7}{R_E}$。
  推导：分压偏置由戴维南等效，基极电阻小则 $V_B$ 稳定，$\therefore I_C$ 与 $\beta$ 弱相关。
  应用：稳定静态工作点。

- **跨导小信号模型**
  公式：$g_m=\dfrac{\partial I_C}{\partial V_{BE}}=\dfrac{I_C}{V_T}$，$r_\pi=\dfrac{\beta V_T}{I_C}$，$r_o\approx\dfrac{V_A}{I_C}$。
  推导：$g_m=\dfrac{d}{dV_{BE}}I_Se^{V_{BE}/V_T}=\dfrac{I_C}{V_T}$；$r_\pi=\dfrac{\partial V_{BE}}{\partial I_B}=\dfrac{\beta}{g_m}$。
  应用：交流等效电路。

- **共射放大电路增益**
  公式：$A_u=-g_m(R_C\parallel R_L)$，$R_i=r_\pi\parallel R_B$，$R_o=R_C$。
  推导：输出 $u_o=-g_mu_{be}(R_C\parallel R_L)$，输入 $u_i=u_{be}$，$\therefore A_u=u_o/u_i$。
  物理意义：反相放大。
  应用：电压放大。

- **共集（射极跟随器）**
  公式：$A_u=\dfrac{(1+\beta)(R_E\parallel R_L)}{r_\pi+(1+\beta)(R_E\parallel R_L)}\approx1$，$R_o\approx\dfrac{R_s\parallel r_\pi}{1+\beta}$。
  推导：由小信号等效电路列 KVL 求解。
  应用：缓冲、阻抗变换。

- **频率响应与密勒效应**
  公式：$f_H$ 受 $C_\pi(1+|A_u|)$ 影响（密勒电容）。
  推导：跨接电容折算为输入等效 $C_{in}=C_\pi(1+A)$，增大时间常数降低带宽。
  应用：宽带放大器。

### 第12章 场效应晶体管放大电路

- **JFET/MOSFET 结构与工作原理**
  定义：电压控制型器件，栅极电压控制沟道导电。
  公式：$I_D=f(V_{GS},V_{DS})$。
  物理意义：高输入阻抗。
  应用：高频、低功耗。

- **转移特性与跨导**
  公式：增强型 NMOS 饱和区 $I_D=\dfrac12k_n'(W/L)(V_{GS}-V_{TH})^2$；$g_m=\dfrac{\partial I_D}{\partial V_{GS}}=k_n'(W/L)(V_{GS}-V_{TH})=\sqrt{2k_n'(W/L)I_D}$。
  推导：对转移特性求导。沟道电荷与漂移电流积分（肖克利模型）得平方律。
  应用：增益设计。

- **MOSFET 偏置**
  公式：分压自偏置 $V_G=\dfrac{R_2}{R_1+R_2}V_{DD}$，$V_{GS}=V_G-I_DR_S$。
  推导：栅极电流为零，列输入输出回路方程联立。
  应用：Q 点设置。

- **共源放大电路**
  公式：$A_u=-g_m(R_D\parallel R_L)$，$R_i=R_G$，$R_o=R_D$。
  推导：与小信号模型同。
  应用：通用放大。

- **共漏（源极跟随器）**
  公式：$A_u=\dfrac{g_m(R_S\parallel R_L)}{1+g_m(R_S\parallel R_L)}\approx1$。
  应用：源极缓冲。

- **CMOS 反相器/放大基础**
  定义：NMOS 与 PMOS 互补。
  公式：阈值损失、开关功耗 $P=\alpha CV^2f$。
  应用：数字集成电路。

### 第13章 集成运算放大器

- **理想运放模型**
  定义：$A_{od}\to\infty$，$R_{id}\to\infty$，$R_o\to0$，带宽 $\to\infty$。
  公式：$u_o=A_{od}(u_+-u_-)$。
  推导：差分放大，理想时 $u_+-u_-\to0$ 即“虚短”。
  应用：线性运算电路。

- **虚短与虚断**
  定义：深度负反馈下 $u_+\approx u_-$，$i_+\approx i_-\approx0$。
  推导：$u_o=A_{od}(u_+-u_-)$ 有限，$A_{od}\to\infty\Rightarrow u_+-u_-\to0$；输入电阻无穷大 $\Rightarrow i=0$。
  应用：快速分析。

- **反相放大电路**
  公式：$A_u=-\dfrac{R_f}{R_1}$，$R_i=R_1$。
  推导：节点 N 虚地 $u_-\approx0$，KCL $\dfrac{u_i}{R_1}+\dfrac{u_o}{R_f}=0\Rightarrow u_o=-\dfrac{R_f}{R_1}u_i$。
  应用：信号反相放大。

- **同相放大电路**
  公式：$A_u=1+\dfrac{R_f}{R_1}$。
  推导：$u_-=\dfrac{R_1}{R_1+R_f}u_o=u_i\Rightarrow A_u=\dfrac{R_1+R_f}{R_1}$。
  应用：高输入阻抗放大。

- **求和、减法、积分与微分**
  公式：加法 $u_o=-R_f\sum\dfrac{u_i}{R_i}$；积分 $u_o=-\dfrac1{RC}\int u_i dt$；微分 $u_o=-RC\dfrac{du_i}{dt}$。
  推导：积分器由 KCL $C\dfrac{d}{dt}(0-u_o)=\dfrac{u_i}R$，积分得式。
  应用：模拟计算。

- **电压比较器**
  公式：$u_o=\begin{cases}U_{sat+},&u_+>u_-\\U_{sat-},&u_+<u_-\end{cases}$。
  应用：阈值检测、方波发生。

- **实际运放的失调与共模抑制**
  公式：$u_o=A_{od}u_d+A_{cm}u_{cm}$，CMRR $=20\lg\dfrac{A_{od}}{A_{cm}}$。
  应用：精密测量。

### 第14章 反馈放大电路与信号产生

- **反馈的四种类型**
  定义：电压/电流反馈、串联/并联反馈。
  公式：闭环增益 $A_f=\dfrac{A}{1+AF}$。
  推导：由 $X_o=AX_i'$，$X_i'=X_i-FX_o$，消去中间量
  $$X_o=A(X_i-FX_o)\ \Rightarrow\ X_o(1+AF)=AX_i\ \Rightarrow\ A_f=\frac{A}{1+AF}.$$
  物理意义：负反馈稳定增益。
  应用：放大器设计。

- **负反馈对输入输出电阻的影响**
  公式：串联反馈 $R_{if}=R_i(1+AF)$；并联反馈 $R_{if}=\dfrac{R_i}{1+AF}$；电压反馈 $R_{of}=\dfrac{R_o}{1+AF}$；电流反馈 $R_{of}=R_o(1+AF)$。
  推导：由环路方程求等效端口电阻。
  应用：阻抗匹配。

- **负反馈对带宽与失真的改善**
  公式：$BW_f=BW(1+AF)$，增益-带宽积恒定。
  推导：$A_f=A/(1+AF)$，$\therefore A_f\cdot BW_f=A\cdot BW$。
  应用：宽带放大。

- **反馈稳定性与自激**
  定义：环路增益 $AF$ 在相移 $180^\circ$ 时幅值 $\ge1$ 即自激。
  判据：相位裕度 $>45^\circ$，增益裕度 $>10\,\mathrm{dB}$。
  推导：闭环极点由 $1+AF=0$ 决定，右半平面极点则发散。
  应用：补偿网络。

- **巴克豪森振荡条件**
  公式：$|AF|=1$，$\angle AF=0^\circ$。
  推导：由自激时 $X_i=0$ 但 $X_o\ne0$，$X_o=AFX_o\Rightarrow AF=1$。
  应用：正弦振荡器。

- **RC 桥式振荡器（文氏电桥）**
  公式：$f_0=\dfrac1{2\pi RC}$，$R_f/R_1=2$。
  推导：文氏网络 $H=\dfrac{1}{3+j(\omega RC-1/\omega RC)}$，谐振 $\omega_0=1/RC$ 时 $H=1/3$，$\therefore$ 放大倍数须为 $3$。
  应用：音频信号源。

- **LC 振荡器**
  公式：考毕兹 $f_0=\dfrac1{2\pi\sqrt{LC}}$。
  应用：射频振荡。

### 第15章 功率放大电路与直流稳压电源

- **功率放大器的分类**
  定义：甲类（$360^\circ$）、乙类（$180^\circ$）、甲乙类（$180^\circ\sim360^\circ$）、丁类（开关）。
  物理意义：导通角决定效率。
  应用：音频、射频功放。

- **甲类放大器效率**
  公式：$\eta_{max}=25\%$（电阻负载）。
  推导：最大输出功率 $P_o=U_{CC}^2/(8R_L)$，电源功率 $P_V=U_{CC}^2/(2R_L)$，$\therefore \eta=P_o/P_V=1/4$。
  应用：小功率高保真。

- **乙类互补对称功放效率**
  公式：$\eta_{max}=\dfrac\pi4\approx78.5\%$，$P_{om}=\dfrac{(U_{CC}-U_{CES})^2}{2R_L}$。
  推导：每管导通半周，输出 $P_o=\dfrac{1}{2}\dfrac{U_{om}^2}{R_L}$，电源 $P_V=\dfrac{2}{T}\int_0^{T/2}\dfrac{U_{om}}{R_L}\sin\omega t\cdot U_{CC}dt=\dfrac{2U_{CC}U_{om}}{\pi R_L}$，
  $$\therefore \eta=\frac{\pi}{4}\frac{U_{om}}{U_{CC}}\xrightarrow{U_{om}\to U_{CC}}\frac\pi4 .$$
  物理意义：效率高于甲类。
  特例：交越失真，需给偏置成甲乙类。
  应用：音频功放。

- **串联型稳压电源**
  定义：调整管、取样、比较放大、基准。
  公式：$U_o=\dfrac{R_1+R_2}{R_2}U_{REF}$。
  推导：负反馈使取样电压等于基准，即 $\dfrac{R_2}{R_1+R_2}U_o=U_{REF}$。
  应用：线性稳压。

- **稳压系数与输出电阻**
  公式：$S_r=\dfrac{\Delta U_o/U_o}{\Delta U_i/U_i}$，$R_o=\dfrac{\Delta U_o}{\Delta I_o}$。
  应用：电源指标。

- **开关电源原理**
  定义：PWM 控制占空比。
  公式：Buck $U_o=DU_i$。
  推导：电感伏秒平衡 $D(U_i-U_o)=(1-D)U_o$，$\therefore U_o=DU_i$。
  应用：高效电源。

### 第16章 数字逻辑基础与门电路

- **数制与码制**
  定义：二进制、十六进制、BCD。
  公式：$(N)_2=\sum a_i2^i$。
  推导：按权展开。
  应用：数字系统表示。

- **逻辑代数基本定律**
  公式：$A+AB=A$，$A+A'B=A+B$，摩根 $\overline{A+B}=A'B'$。
  推导：摩根由真值表穷举验证；吸收律 $A+AB=A(1+B)=A$。
  应用：逻辑化简。

- **逻辑门与真值表**
  定义：与、或、非、与非、或非、异或。
  公式：$Y=A\cdot B$，$Y=A+B$，$Y=\overline A$，$Y=A\oplus B=A'B+AB'$。
  应用：组合电路。

- **最小项与卡诺图化简**
  定义：$n$ 变量最小项共 $2^n$ 个。
  方法：相邻格合并消去变量，圈尽可能大。
  推导：$AB+AB'=A(B+B')=A$。
  应用：两级门优化。

- **TTL 与非门电路**
  定义：输入级、倒相级、输出级。
  公式：$V_{OH}\ge2.4\,\mathrm{V}$，$V_{OL}\le0.4\,\mathrm{V}$。
  应用：传统逻辑。

- **CMOS 门电路**
  定义：NMOS 下拉网络、PMOS 上拉网络。
  公式：静态功耗近似为零。
  应用：低功耗 IC。

- **竞争冒险**
  定义：信号延迟导致毛刺。
  消除：加冗余项、输出加滤波电容。
  应用：时序设计。

### 第17章 组合逻辑电路

- **组合逻辑分析与设计**
  定义：输出仅由当前输入决定。
  步骤：写真值表 → 化简 → 画电路。
  应用：数字逻辑设计。

- **加法器**
  公式：半加 $S=A\oplus B$，$C=AB$；全加 $S_i=A_i\oplus B_i\oplus C_i$。
  推导：由真值表并结合异或定义。
  应用：ALU。

- **编码器与优先编码器**
  公式：$2^n$ 输入对应 $n$ 位输出。
  应用：键盘编码。

- **译码器**
  公式：$n$ 位输入对应 $2^n$ 输出，每个输出为最小项。
  应用：地址译码、显示。

- **数据选择器（MUX）**
  公式：$Y=\sum m_i D_i$。
  推导：由使能端与地址线选择某一路数据输出。
  应用：逻辑函数实现。

- **数值比较器**
  公式：$A>B$、$A=B$、$A<B$ 输出。
  应用：大小判断。

### 第18章 时序逻辑电路、定时器与数据转换

- **锁存器与触发器**
  定义：SR、D、JK、T 触发器。
  公式：$Q^{n+1}=D$；$Q^{n+1}=JQ'^n+K'Q^n$。
  推导：由特征表化简得特征方程。
  应用：存储、计数。

- **寄存器和移位寄存器**
  定义：并行/串行输入输出。
  应用：数据缓存。

- **计数器**
  公式：模 $N$ 计数器需 $\lceil\log_2N\rceil$ 个触发器。
  推导：$n$ 个触发器有 $2^n$ 状态。
  应用：分频、计时。

- **555 定时器**
  公式：单稳态 $t_w=1.1RC$；多谐振荡 $T=0.693(R_1+2R_2)C$。
  推导：电容充放电由阈值 $1/3V_{CC}$、$2/3V_{CC}$ 决定，利用 RC 指数
  $$t=R C\ln\frac{V_\infty-V_0}{V_\infty-V_{th}} .$$
  计算即得。
  应用：脉冲产生。

- **D/A 转换器**
  公式：$u_o=-\dfrac{U_{REF}}{2^n}\sum_{i=0}^{n-1}D_i2^i$。
  推导：权电阻/ R-2R 网络由虚地求和。
  应用：波形发生。

- **A/D 转换器**
  公式：量化误差 $\pm\dfrac12\mathrm{LSB}$，LSB $=\dfrac{U_{REF}}{2^n}$。
  推导：逐次逼近、双积分原理。
  应用：数据采集。

---
⬆ 上级领域：[[电子·控制·系统]]
