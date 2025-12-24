# Coulomb Distorted Exact WKB Method for Narrow Resonances

## 用于极窄共振（如α衰变）的库仑畸变精确WKB方法

---

## 目录

1. [引言与动机](#1-引言与动机)
2. [理论框架](#2-理论框架)
3. [库仑势的WKB分析](#3-库仑势的wkb分析)
4. [Stokes图结构](#4-stokes图结构)
5. [库仑修正的连接公式](#5-库仑修正的连接公式)
6. [S矩阵的解析表达式](#6-s矩阵的解析表达式)
7. [与有效程展开的结合](#7-与有效程展开的结合)
8. [数值实现](#8-数值实现)
9. [参考文献](#9-参考文献)

---

## 1. 引言与动机

### 1.1 问题背景

对于α衰变这样的极窄共振，共振宽度可以小到 Γ ~ 10⁻²⁰ eV，对应的半衰期从微秒到数十亿年不等。传统的散射相移方法在这种情况下面临严重困难：

- 相移在共振附近的变化 Δδ ~ π 发生在极窄的能量范围内
- 直接数值计算无法分辨如此精细的结构
- 需要解析延拓到复能量/复动量平面

### 1.2 两种互补方法

**方法A：S矩阵极点法（您的笔记）**
- 在复k平面搜索 1/S(k) = 0
- 用有效程展开(ERF)参数化低能散射
- 数值求解薛定谔方程匹配库仑波函数

**方法B：Exact WKB（arXiv:2508.09211）**
- 用Borel重求和处理WKB级数
- 通过Stokes现象分析共振态
- 给出S矩阵的解析表达式

**本文目标**：将两种方法结合，发展**库仑畸变的Exact WKB方法**。

---

## 2. 理论框架

### 2.1 薛定谔方程

考虑α粒子与子核系统，约化质量为μ，相对坐标为r：

```
[-ℏ²/(2μ) d²/dr² + V_eff(r) - E] u(r) = 0
```

有效势能：

```
V_eff(r) = V_N(r) + V_C(r) + V_cent(r)
```

其中：
- **核势** V_N(r)：短程吸引势，如Woods-Saxon
  ```
  V_N(r) = -V₀ / [1 + exp((r-R)/a)]
  ```
- **库仑势** V_C(r)：长程排斥势
  ```
  V_C(r) = Z₁Z₂e²/r = 2ηℏ²k/(μr)    (r > R_C)
  ```
- **离心势** V_cent(r)：
  ```
  V_cent(r) = ℓ(ℓ+1)ℏ²/(2μr²)
  ```

### 2.2 渐近行为

**短程势（标准WKB）**：
```
u(r → ∞) ~ sin(kr - ℓπ/2 + δ_ℓ)
         ~ (1/2i)[e^{i(kr-ℓπ/2+δ)} - e^{-i(kr-ℓπ/2+δ)}]
```

**库仑势（Coulomb distorted）**：
```
u(r → ∞) ~ sin(kr - ηln(2kr) - ℓπ/2 + σ_ℓ + δ_ℓ)
```

或用库仑波函数表示：
```
u(r → ∞) ~ A·F_ℓ(η,kr) + B·G_ℓ(η,kr)
```

其中：
- η = Z₁Z₂e²μ/(ℏ²k) 是Sommerfeld参数
- σ_ℓ = arg Γ(ℓ+1+iη) 是库仑相移
- F_ℓ, G_ℓ 是正则和非正则库仑波函数

### 2.3 S矩阵定义

**库仑修正的S矩阵**：
```
S_ℓ = e^{2i(σ_ℓ + δ_ℓ)}
```

用入射波和出射波表示：
```
u(r → ∞) ~ (1/S_ℓ)·φ_c⁻(kr) + φ_c⁺(kr)
```

其中：
```
φ_c⁺(kr) = G_ℓ(η,kr) + iF_ℓ(η,kr)  (出射波)
φ_c⁻(kr) = G_ℓ(η,kr) - iF_ℓ(η,kr)  (入射波)
```

**共振条件**：S矩阵在复k平面的极点，即 1/S_ℓ(k) = 0

---

## 3. 库仑势的WKB分析

### 3.1 标准WKB解

定义局部动量：
```
p(r) = √[2μ(E - V_eff(r))]
     = ℏ√[k² - 2ηk/r - ℓ(ℓ+1)/r² - 2μV_N(r)/ℏ²]
```

标准WKB波函数：
```
ψ_WKB(r) = 1/√p(r) · exp(±(i/ℏ)∫p(r')dr')
```

### 3.2 Langer修正

**关键替换**：对于库仑型势，需要做Langer修正
```
ℓ(ℓ+1) → (ℓ + 1/2)²
```

修正后的有效动量：
```
k_eff(r) = √[k² - 2ηk/r - (ℓ+1/2)²/r² - 2μV_N(r)/ℏ²]
```

**物理意义**：
- 消除WKB在r=0处的奇异性
- 对纯库仑势给出精确能谱
- 保证正确的渐近行为

### 3.3 转折点分析

转折点由 k_eff(r_t) = 0 定义：
```
k² - 2ηk/r_t - (ℓ+1/2)²/r_t² - 2μV_N(r_t)/ℏ² = 0
```

对于α衰变典型情况（库仑势垒 + 核势阱）：

**三个转折点**：
- r₁：内转折点（核势阱边缘）
- r₂：中间转折点（势垒内侧）
- r₃：外转折点（势垒外侧）

```
        V(r)
          │
    V_B ──┼───────╮
          │       │╲
          │       │ ╲
    E ────┼───────┼──╲────────────
          │   ┌───┤   ╲
          │   │   │    ╲__________
          │   │   │
    0 ────┼───┴───┴───────────────→ r
          │   r₁  r₂  r₃
          │
   V_N ───┘
```

---

## 4. Stokes图结构

### 4.1 Exact WKB中的Stokes图

**定义**：Stokes线是复平面上满足以下条件的曲线
```
Im ∫_{r_t}^{r} k_eff(r')dr' = 0
```

**物理意义**：
- 在Stokes线上，WKB解的主导项发生转换
- 穿越Stokes线时需要应用连接公式
- 共振态对应特定的Stokes图拓扑

### 4.2 纯库仑势的Stokes图

对于纯库仑势（无核势），有效动量：
```
k_C(r) = √[k² - 2ηk/r - (ℓ+1/2)²/r²]
```

**转折点位置**（令 k_C = 0）：
```
r_± = (η ± √(η² + (ℓ+1/2)²)) / k
```

对于 E > 0 的散射态：
- r₊ > 0 是物理转折点（经典反射点）
- r₋ < 0 是非物理的

**复平面结构**：
```
                Im r
                  │
                  │    Stokes线
                  │   ╱
                  │  ╱
                  │ ╱
    ──────────────┼────────────→ Re r
           r₋    │0    r₊
                  │ ╲
                  │  ╲
                  │   ╲
                  │
```

### 4.3 库仑+核势的Stokes图

加入短程核势后，结构变得更复杂：

**低能情况（E < V_barrier）**：
```
                Im r
                  │
           (-)   │   (+)
                  │╲    ╱
                  │ ╲  ╱  Stokes线
    ──────────────┼──╳────●────●────→ Re r
                  │ ╱r₁╲ r₂   r₃
           (+)   │╱    ╲
                  │   (-)
                  │
```

**共振能量**：当能量接近准束缚态时，Stokes图发生拓扑变化。

### 4.4 Siegert边界条件

**物理要求**：共振态只有出射波，无入射波

在Exact WKB框架中，这对应于：
```
ψ⁻ 在 Im r → -∞ 处占主导（Re r > 0）
```

其中 ψ⁻ ∝ exp(-∫S_odd dr) 是衰减解。

**与S矩阵极点的关系**：
```
Siegert条件 ⟺ 1/S(k) = 0 ⟺ 入射波系数为零
```

---

## 5. 库仑修正的连接公式

### 5.1 标准WKB连接公式（Airy函数）

在转折点 r_t 附近，势能近似线性：
```
V_eff(r) ≈ V_eff(r_t) + V'_eff(r_t)·(r - r_t)
```

标准连接公式（从经典禁戒区到允许区）：
```
1/√|p| exp(-1/ℏ ∫|p|dr) ↔ 2/√p cos(1/ℏ ∫p dr - π/4)
```

### 5.2 库仑修正的连接公式

对于库仑势，需要用**库仑波函数**而非Airy函数作为匹配函数。

**Froman-Froman方法**：

定义相位积分：
```
w(r) = ∫_{r_t}^{r} k_eff(r')dr'
```

在库仑场中的连接公式：
```
经典禁戒区 (r < r_t):
    ψ ~ 1/√|k_eff| · exp(-|w|)

经典允许区 (r > r_t):
    ψ ~ 1/√k_eff · [A·sin(w + φ_C) + B·cos(w + φ_C)]
```

其中 φ_C 是库仑修正相位：
```
φ_C = σ_ℓ - ηln(2kr_t) + π/4 + O(1/η)
```

### 5.3 多转折点情况

对于α衰变的三转折点问题，需要依次穿越：

**区域划分**：
- 区域 I：r < r₁（核势阱内，振荡）
- 区域 II：r₁ < r < r₂（势垒下方，指数衰减）
- 区域 III：r₂ < r < r₃（势垒内，指数）
- 区域 IV：r > r₃（渐近区，振荡）

**穿透因子**：
```
T = exp(-2γ)
```

其中Gamow因子：
```
γ = (1/ℏ) ∫_{r₂}^{r₃} |p(r)|dr
  = ∫_{r₂}^{r₃} √[2μ(V_eff(r) - E)/ℏ²] dr
```

---

## 6. S矩阵的解析表达式

### 6.1 WKB近似下的S矩阵

**单转折点情况**（高能散射）：
```
S_ℓ = exp(2iδ_ℓ)
```

其中WKB相移：
```
δ_ℓ^WKB = lim_{r→∞} [∫_{r_t}^{r} k_eff(r')dr' - kr + ηln(2kr)] - (ℓ+1/2)π/2 + σ_ℓ
```

### 6.2 势垒穿透的S矩阵

对于有势垒的情况，S矩阵包含穿透和反射贡献：
```
S_ℓ = S_ℓ^(refl) + S_ℓ^(trans)
```

**反射项**（在外转折点r₃反射）：
```
S_ℓ^(refl) = -exp(2i[σ_ℓ + δ_ℓ^(out)])
```

**穿透项**（隧穿整个势垒）：
```
S_ℓ^(trans) = exp(2i[σ_ℓ + δ_ℓ^(in)]) · T
```

其中 T = exp(-2γ) 是穿透因子。

### 6.3 共振处的S矩阵

在共振能量 E_R 附近，S矩阵具有极点结构：
```
S_ℓ(E) ≈ S_ℓ^(bg) · (E - E_R + iΓ/2) / (E - E_R - iΓ/2)
```

**共振参数的WKB表达式**：

共振能量 E_R 由Bohr-Sommerfeld条件确定：
```
∫_{r₁}^{r₂} k_eff(r)dr = (n + 1/2)π    (n = 0, 1, 2, ...)
```

共振宽度 Γ 由穿透因子给出：
```
Γ = (ℏv/2R) · T = (ℏv/2R) · exp(-2γ)
```

其中：
- v = ℏk/μ 是粒子速度
- R 是核半径
- γ 是Gamow因子

### 6.4 Exact WKB的S矩阵表达式

参考arXiv:2508.09211的方法，对于可解析求解的势能：
```
1/T = 1/√A - √A
```

其中 A 是由超几何函数构成的表达式。

**推广到库仑+核势**：

需要用**库仑超几何函数**（Confluent hypergeometric function）：
```
F_ℓ(η,ρ) = C_ℓ(η) · ρ^{ℓ+1} · e^{-iρ} · ₁F₁(ℓ+1+iη, 2ℓ+2, 2iρ)
```

S矩阵可以表示为：
```
S_ℓ = e^{2iσ_ℓ} · Γ(ℓ+1+iη)/Γ(ℓ+1-iη) · [匹配系数比]
```

---

## 7. 与有效程展开的结合

### 7.1 有效程展开(ERF)回顾

**您笔记中的方法**：

库仑修正的有效程公式：
```
C_ℓ²(η) · k^{2ℓ+1} · cot(δ_ℓ) = -1/a_ℓ + (1/2)r_ℓ k² + ...
```

其中 C_ℓ²(η) 是库仑穿透因子：
```
C₀²(η) = 2πη / (e^{2πη} - 1)
C_ℓ²(η) = C₀²(η) · ∏_{s=1}^{ℓ} (s² + η²) / s²
```

### 7.2 WKB参数的ERF解释

**散射长度 a_ℓ 的WKB表达式**：
```
1/a_ℓ ≈ k · cot(δ_ℓ^WKB)|_{k→0}
```

在低能极限：
```
δ_ℓ^WKB → ∫_{r_t}^{∞} [k_eff(r) - k]dr - ηln(2kr_t) + O(k)
```

**有效程 r_ℓ 的WKB表达式**：
```
r_ℓ ≈ 2 · d(k·cot δ_ℓ)/dk²|_{k→0}
```

### 7.3 S矩阵极点的解析延拓

**ERF方法**（您的笔记）：
```
找 k_R 使得：k_R · C_ℓ²(η_R) · cot(δ_ℓ) = i

即：-1/a_ℓ + (1/2)r_ℓ k_R² + ... = i · k_R · C_ℓ²(η_R)
```

**Exact WKB方法**：
```
找 k_R 使得：A(k_R) = 1  (见公式 1/T = 1/√A - √A)

即透射系数 T 发散 → S矩阵有极点
```

**两种方法的等价性**：

在低能极限，ERF参数可以从WKB相移展开系数提取：
```
WKB相移 δ_ℓ(k) → 展开得到 a_ℓ, r_ℓ, ... → ERF方程 → S极点 k_R
```

### 7.4 结合策略

**步骤1**：用Exact WKB计算相移 δ_ℓ(k) 的解析表达式

**步骤2**：从 δ_ℓ(k) 的低能展开提取ERF参数
```
a_ℓ = -lim_{k→0} tan(δ_ℓ)/k
r_ℓ = 2 · d[k·cot(δ_ℓ)]/dk²|_{k→0}
```

**步骤3**：用ERF公式解析延拓到复k平面
```
k_R{...}·i = -1/a_ℓ + (1/2)r_ℓ k_R² + ...
```

**步骤4**：求解得到共振位置 k_R = k_r - i·k_i
```
E_R = ℏ²k_r²/(2μ) - ℏ²k_i²/(2μ)
Γ = 2ℏ²k_r·k_i/μ
```

---

## 8. 数值实现

### 8.1 算法流程

```python
def compute_alpha_decay_resonance(Z1, Z2, A, V0, R0, a, ell=0):
    """
    计算α衰变共振参数

    参数:
        Z1, Z2: 电荷数
        A: 质量数
        V0, R0, a: Woods-Saxon势参数
        ell: 角动量量子数

    返回:
        E_R: 共振能量 (MeV)
        Gamma: 共振宽度 (MeV)
        T_half: 半衰期 (s)
    """

    # 1. 定义势能
    mu = reduced_mass(A)
    V_eff = lambda r: woods_saxon(r, V0, R0, a) + \
                      coulomb(r, Z1, Z2) + \
                      centrifugal(r, ell, mu)

    # 2. 找转折点
    r1, r2, r3 = find_turning_points(V_eff, E_guess)

    # 3. 计算WKB相位积分
    def phase_integral(E):
        k = sqrt(2*mu*E)/hbar
        eta = sommerfeld_parameter(Z1, Z2, mu, k)

        # 内部区域积分
        w_in = integrate(k_eff, r1, r2)

        # 势垒区域积分 (Gamow因子)
        gamma = integrate(abs_k_eff, r2, r3)

        return w_in, gamma

    # 4. Bohr-Sommerfeld量子化 → E_R
    def quantization_condition(E):
        w_in, _ = phase_integral(E)
        return w_in - (n + 0.5) * pi

    E_R = root_finder(quantization_condition, E_guess)

    # 5. 计算穿透因子 → Γ
    _, gamma = phase_integral(E_R)
    T = exp(-2 * gamma)

    k_R = sqrt(2*mu*E_R)/hbar
    v = hbar * k_R / mu
    Gamma = hbar * v / (2 * R0) * T

    # 6. 半衰期
    T_half = hbar * log(2) / Gamma

    return E_R, Gamma, T_half
```

### 8.2 库仑修正的数值细节

**Sommerfeld参数**：
```python
def sommerfeld_parameter(Z1, Z2, mu, k):
    """η = Z1*Z2*e²*μ/(ℏ²k)"""
    alpha = 1/137.036  # 精细结构常数
    return Z1 * Z2 * alpha * mu * c / (hbar * k)
```

**库仑波函数**（使用mpmath或scipy）：
```python
from scipy.special import coulomb_f, coulomb_g

def coulomb_waves(ell, eta, rho):
    """计算正则和非正则库仑波函数"""
    F = coulomb_f(ell, eta, rho)
    G = coulomb_g(ell, eta, rho)
    return F, G
```

**Gamow因子的数值积分**：
```python
def gamow_factor(E, V_eff, r2, r3, mu):
    """计算Gamow穿透因子"""
    def integrand(r):
        return sqrt(2*mu*(V_eff(r) - E)) / hbar

    gamma, _ = quad(integrand, r2, r3)
    return gamma
```

### 8.3 测试用例

**²¹²Po → ²⁰⁸Pb + α**：
```
输入参数:
  Z1 = 2 (α)
  Z2 = 82 (Pb)
  A_daughter = 208
  Q = 8.95 MeV

期望输出:
  T_1/2 ≈ 0.3 μs
```

**²³⁸U → ²³⁴Th + α**：
```
输入参数:
  Z1 = 2 (α)
  Z2 = 90 (Th)
  A_daughter = 234
  Q = 4.27 MeV

期望输出:
  T_1/2 ≈ 4.5 × 10⁹ years
```

---

## 9. 参考文献

### 核心文献

1. **Exact WKB与共振**
   - O. Morikawa and S. Ogawa, arXiv:2508.09211 [quant-ph] (2025)
   - O. Morikawa and S. Ogawa, arXiv:2503.18741 [hep-th] (2025)
   - O. Morikawa and S. Ogawa, arXiv:2505.02301 [hep-th] (2025)

2. **S矩阵极点与α衰变**
   - arXiv:1612.04135 - "An estimate of Alpha decay half-life from the poles of S-matrix"
   - Phys. Rev. C 96, 044602 (2017) - α-induced scattering and resonance-pole-decay model

3. **WKB方法**
   - N. Fröman and P.O. Fröman, "JWKB Approximation: Contributions to the Theory" (North Holland, 1965)
   - Phys. Rev. A 36, 5523 (1987) - Generalized WKB for screened Coulomb

4. **有效程展开**
   - Phys. Lett. B 833, 137273 (2022) - Effective range expansion for narrow resonances
   - arXiv:2110.07484 - ERF for near-threshold resonances

5. **库仑波函数**
   - Abramowitz & Stegun, Chapter 14
   - DLMF (Digital Library of Mathematical Functions), Chapter 33

### 教科书

- L.D. Landau and E.M. Lifshitz, "Quantum Mechanics" (Pergamon, 1977), §134-135
- J.R. Taylor, "Scattering Theory" (Wiley, 1972), Chapter 12
- A. Messiah, "Quantum Mechanics" (Dover, 1999), Chapter XIX

---

## 附录A：常用公式

### A.1 库仑波函数的渐近形式

大ρ展开：
```
F_ℓ(η,ρ) → sin(ρ - ηln(2ρ) - ℓπ/2 + σ_ℓ)
G_ℓ(η,ρ) → cos(ρ - ηln(2ρ) - ℓπ/2 + σ_ℓ)
```

小ρ展开：
```
F_ℓ(η,ρ) → C_ℓ(η) · ρ^{ℓ+1}
G_ℓ(η,ρ) → ρ^{-ℓ} / [(2ℓ+1)C_ℓ(η)]
```

### A.2 库仑相移

```
σ_ℓ = arg Γ(ℓ + 1 + iη)
    = Σ_{n=1}^{ℓ} arctan(η/n) - γη + η[ln|η| - 1] + O(1/η)
```

其中 γ ≈ 0.5772 是Euler常数。

### A.3 Gamow公式

α衰变半衰期的Gamow近似：
```
log₁₀(T_1/2) = a·Z/√E + b
```

其中 a ≈ 1.5, b 是与核结构相关的常数。

更精确的Geiger-Nuttall关系：
```
log₁₀(λ) = A - B·Z/√Q
```

---

## 附录B：符号表

| 符号 | 含义 | 单位 |
|-----|------|-----|
| E | 能量 | MeV |
| k | 波数 | fm⁻¹ |
| η | Sommerfeld参数 | 无量纲 |
| σ_ℓ | 库仑相移 | rad |
| δ_ℓ | 核相移 | rad |
| Γ | 共振宽度 | MeV |
| γ | Gamow因子 | 无量纲 |
| T | 穿透因子 | 无量纲 |
| a_ℓ | 散射长度 | fm |
| r_ℓ | 有效程 | fm |

---

*文档版本: 1.0*
*最后更新: 2025-12-22*
